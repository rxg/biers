#lang racket

;(require test-engine/racket-tests) ;; for check-expect
(require rackunit)
(require math/distributions)

;; stub for unimplemented code
(define ... (λ args (error "not yet implemented!")))

;; Used for testing multiple-value functions
(module+ test
  (define-syntax list-args
    (syntax-rules ()
      [(_ e) (call-with-values (λ () e) list)])))

;;
;; modelizer.rkt - Some tools for approximately fitting
;;                 Bayesian statistical models
;;
;; Author: Ron Garcia <rxg@cs.ubc.ca>
;; Date: February 22, 2019
;;


;; RG: Introduce notions of "fit environment" and "compatibility environment"
;;  currently called "data" and "env"

;;
;; Quilt: Quadratic Approximation Intermediate Language with Types
;; Model Specification:
;;
;; This is an explicitly typed stat-model language.  Later we'll want to
;; infer types from the combination of an implicitly-typed model specification
;; and the data provided to fit with it.
;;

;; Syntax:
;; v ∈ Variable
;; tdecl ∈ TypeDeclaration
;; vdef ∈ VariableDefinition
;; i ∈ IndexMetavariable
;; k ∈ Index
;; r ∈ VariablePattern
;; q ∈ VariableReference
;; n ∈ Number
;; m ∈ Model
;; e ∈ Expression
;; d ∈ DataType
;; t ∈ TypeSpec
;; s ∈ Symbol
;; p ∈ PriorDistribution

;; m ::=
;; (model
;;   [type-decls tdecl ...]
;;   [var-decls [r d] ...]
;;   [var-defs vdef ...]
;;   [data ...]                ;; data used to fit the model, if any
;;   [fit ...])                ;; posterior, if the model has been fit to data

;; tdecl ::= [i Row]                 ;; Row index metavariable
;;         | [i d (Enum s ...)]      ;; Enum index, type name, and elements
;; Alternate uniform decomposition:
;;    tdecl ::= [i . t]  where t ::= (Row) | (d (Enum s ...))    

;; vdef ::= [r . ~ . p]              ;; prior distribution
;;        | [r . ~ . irrelevant]     ;; omitted prior (for predictors)
;;        | [r . = . e]              ;; definition 

;; q ::= v           ;; atomic variable
;;     | (v e)       ;; indexed variable
;;     | #(q q ...)  ;; vector of variables

;; r ::= v           ;; atomic variable
;;     | (v i)       ;; indexed family of variables
;;     | (v e)       ;; indexed variable
;;     | #(r r ...)  ;; vector of variables

;; p ::= (normal e e)
;;     | (log-normal e e)
;;     | (multivariate-normal #(e e ...)
;;                            #(#(e e ...)))
;;     | (cauchy e e)
;;     | (laplace e)
;;     | (discrete (e n) ...)
;;     | (exponential e)
;;     | (binomial e e)
;;     | (uniform e e)

;; e ::= n
;;     | s
;;     | v
;;     | (v e)
;;     | (+ e e)
;;     | (* e e)
;;     | (expt e e)
;;     | (exp e)
;;     | #(e e ...)

;; k ::= s | n  ;; actual indices are either symbols or numbers

;; Notes:
;; - dependency order for var-defs is last-to-first, to match stat models.
;;   Think of it like let* in reverse. No re-ordering to match dependencies.
;; - data type, index variables, and variable are disjoint
;;   by convention, data types are capitalized, index 

;; Terminological notes:
;; - I will distinguish between several notions of variable:
;;   1) variable name (e.g. μ, α, x): just the unique root name
;;   2) variable scheme (e.g. (μ i), (α s), x): with relevant index metavariable
;;   3) variable reference (e.g. (μ 7), (α 'male), x):  a particular instance.


;; Chapter 3: Gaussian model of Kalamari forager height
(define ex1
  `(model
    [type-decls [i Row]]
    [var-decls  [(h i) Number] [μ Number] [σ Number]]
    [var-defs   [(h i) . ~ . (normal μ σ)]
                [μ     . ~ . (normal 178 20)]
                [σ     . ~ . (uniform 0 50)]]))

;; Chapter 3: Schematic of a linear regression (force 0 intercept)
(define ex2
  `(model
    [type-decls [i Row]]
    [var-decls  [(y i) Number] [(μ i) Number] [β Number] [σ Number]
                [(x i) Number]]
    [var-defs   [(y i) . ~ . (normal (μ i) σ)]
                [(μ i) . = . (* β (x i))]
                [β     . ~ . (normal 0 10)]
                [σ     . ~ . (exponential 1)]
                [(x i) . ~ . (normal 0 1)]]))

;; Chapter 5:
;; Indicator variable for male (0 = not male, 1 = male)
;; This is essentially Student's t-test encoded as a linear model,
;; without the decision (just the estimation)
(define ex3
  `(model
    [type-decls [i Row]]
    [var-decls  [(h i) Number] [(μ i) Number] [α Number] [βm Number]
                [σ Number] [(m i) Number]]
    [var-defs   [(h i) . ~ . (normal μ σ)]
                [(μ i) . = . (+ α (* βm (m i)))]
                [α     . ~ . (normal 178 20)]
                [βm    . ~ . (normal 0 1)]
                [σ     . ~ . (uniform 0 50)]
                [(m i) . ~ . (discrete (0 0.5) (1 0.5))]]))

;; Chapter 5: Index variable, which is less hacky than an indicator variable IMO
;; This is not exactly the same as above: models (α 'male) and (α 'female)
;; instead of α and (+ α βm) (yuk!)
(define ex4
  `(model
    [type-decls [i Row] [j Sex (Enum 'male 'female)]]
    [var-decls  [(h i) Number] [(μ i) Number] [(α j) Number] [(sex i) Sex]
                [σ Number] [(m i) Number]]
    [var-defs   [(h i) . ~ . (normal μ σ)]
                [(μ i) . = . (α (sex i))]
                [(α j) . ~ . (normal 178 20)]
                [σ     . ~ . (uniform 0 50)]
                [(sex i) . ~ . (discrete ('female 0.5) ('male 0.5))]]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Construct a compatibility function for target parameters by
;; analyzing the model with respect to some data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; To compute:
;; going from last-to-first, compute log-compatibilities
;; - target variable priors just turn into log-probabilities given
;;   the provided target value. Contribute no new bindings.
;; - non-target parameters (i.e., equations) just turn into environment
;;   bindings for use in "earlier" priors.  May form a "vector" or otherwise
;;   indexed set of data. Contribute nothing to log probability
;; - observed variable priors also turn into log-probabilities, though
;;   splayed across the Row observations (so there's an implicit map)
;;   contribute no new bindings.
 


;; RG: Env-related data types are defined toward the end of the file

;; Data is Env
;; first elements are a "data table": equal-length vectors of data. More
;; entries may be added by the analyzer after the data table.

;; Model Data -> (Targets -> Real)
;; Given a model and some data to fit it to, produce a log-compatibility
;; function for its target parameters, suitable for quadratic approximation
(define (make-log-compatibility-fit-function model data)
  ;; instr does the real work of computing the log-compatibility
  (define instr (analyze-model model data))
  ;; targets determines the order of the inputs to our function
  (define targets (get-target-variables model data))

  ;; Targets -> Real
  ;; compute the log compatibility of the given target parameter values
  (λ (tgt-values)
    (define env0 (make-env targets tgt-values))
    (define lcsf0 0)
    (let-values ([(env lcsf) (instr env0 lcsf0)])
      ;; discard the environment once we know the final log-compatibility
      lcsf)))

(module+ test
  (check-equal?
   ((make-log-compatibility-fit-function
     `(model
       [type-decls [i Row]]
       [var-decls  [(h i) Number] [μ Number] [σ Number]]
       [var-defs   [(h i) . ~ . (normal μ σ)]
                   [μ     . ~ . (normal 178 20)]
                   [σ     . ~ . (uniform 0 50)]])
     (make-env '(h) '(#(0 1 2))))
    '(160 4))
   (+ (pdf (normal-dist 160 4) 0 true)
      (pdf (normal-dist 160 4) 1 true)
      (pdf (normal-dist 160 4) 2 true)
      (pdf (normal-dist 178 20) 160 true)
      (pdf (uniform-dist 0 50) 4 true)))

  (check-equal?
   ((make-log-compatibility-fit-function
     `(model
       [type-decls [i Row]]
       [var-decls  [(y i) Number] [(μ i) Number] [β Number] [σ Number]
                   [(x i) Number]]
       [var-defs   [(y i) . ~ . (normal (μ i) σ)]
                   [(μ i) . = . (* β (x i))]
                   [β     . ~ . (normal 0 10)]
                   [σ     . ~ . (exponential 1)]
                   [(x i) . ~ . irrelevant]])
     (make-env '(y x) '(#(1 2 3) #(4 5 6))))
    '(9 12))
   (+ (pdf (normal-dist (* 9 4) 12) 1 true)
      (pdf (normal-dist (* 9 5) 12) 2 true)
      (pdf (normal-dist (* 9 6) 12) 3 true)
      (pdf (normal-dist 0 10) 9 true)
      (pdf (exponential-dist 1) 12 true))
   )
  )

;; Instr is Env LCSF -> Env LCSF
;; our higher-order internal "bytecode" representation for
;; incremental log-compatibility calculations.  An Instr takes an environment
;; of derived parameter bindings and a 
;; log-compatibility-so-far (lcsf) accumulator and updates both according to
;; the encoded instructions. Each Instr can bind new variables and/or
;; contribute to the log-compatibility.

;; LCSF is Real
;; log-compatibility so-far accumulator

;; Model Data -> Instr
;; produce an Instr object that fits the given model to the given data
(define (analyze-model model data)
  (define vdefs (get-variable-defs model))
  (define tdecls (get-tdecls model))
  (analyze-vdefs vdefs data tdecls))

(module+ test
  (check-equal?
   (list-args
    ((analyze-model
      `(model
        [type-decls [i Row]]
        [var-decls  [(h i) Number] [μ Number] [σ Number]]
        [var-defs   [(h i) . ~ . (normal μ σ)]
                    [μ     . ~ . (normal 178 20)]
                    [σ     . ~ . (uniform 0 50)]])
      (make-env '(h) '(#(0 1 2))))
     (make-env '(μ σ) '(160 4))
     0))
   (list (make-env '(μ σ) '(160 4))
         (+ (pdf (normal-dist 160 4) 0 true)
            (pdf (normal-dist 160 4) 1 true)
            (pdf (normal-dist 160 4) 2 true)
            (pdf (normal-dist 178 20) 160 true)
            (pdf (uniform-dist 0 50) 4 true))))
  )

;; VDefs Data TDecls -> Instr
;; turn a list of vdefs (in reverse-binding order) into a sequence Instr
(define (analyze-vdefs vdefs data tdecls)
  (let ([instr* (for/list ([vdef vdefs])
                  (analyze-vdef vdef data tdecls))])
    ;; vdefs are interpreted last-to-first
    (sequence-instrs (reverse instr*))))

(module+ test
  (let () 
    (check-equal?
     (list-args ((analyze-vdefs
                  (list
                   '(σ . ~ . (normal μ 1))
                   '(μ . = . 0))
                  (extend-env empty-env 'σ 2)
                  '([i Row]))
                 empty-env 99))
     `(,(extend-env empty-env 'μ 0)
       ,(+ 99 (pdf (normal-dist 0 1) 2 true))))))

;; Turn a vdef and data table into an inner ll fn which
;; threads internal defs and log-likelihood-so-far.

;; Vdef Data TDecls -> Instr
(define (analyze-vdef vdef data tdecls)
  (match vdef
    [`((,v ,i) . = . ,e)    (analyze-multi-binding v i e data tdecls)]
    [`(,v . = . ,e)         (analyze-single-binding v e data)]
    [`(,r . ~ . irrelevant) (λ (env lcsf) (values env lcsf))] ;; no-op
    [`((,v ,i) . ~ . ,p)    (analyze-multi-prior v i p data tdecls)]
    ;; RG : Need a case for multi-binding via multivariate normal
    ;; RG : Do I want to support simultaneous binding with a vector too?
    ;; [`#(v ...) . ~ . ,p) (analyze-vector-prior v i p data tdecls)]
    [`(,v . ~ . ,p)         (analyze-single-prior v p data)]
    [`,otherwise
     (error 'get-derived-variables "Bad variable definition: ~a"
            otherwise)]))

(module+ test
  (check-equal? (list-args ((analyze-vdef '(σ . = . 9) empty-env '([i Row]))
                            empty-env 99))
                (list-args ((analyze-single-binding 'σ 9 empty-env)
                            empty-env 99)))
  (let ([data (extend-env empty-env 'h #(1 2 3))])
    (check-equal? (list-args ((analyze-vdef '((σ i) . = . 9)
                                            data
                                            '([i Row]))
                              empty-env 99))
                  (list-args ((analyze-multi-binding 'σ 'i 9 data '([i Row]))
                              empty-env 99))))

  (check-equal? (list-args ((analyze-vdef
                             '(σ . ~ . (normal 0 1)) empty-env '([i Row]))
                            (extend-env empty-env 'σ 2) 99))
                `(,(extend-env empty-env 'σ 2)
                  ,(+ 99 (pdf (normal-dist 0 1) 2 true))))
  (let ([data (extend-env empty-env 'σ #(1 2 3))])
    (check-equal? (list-args ((analyze-vdef '((σ i) . ~ . (normal 0 1))
                                            data
                                            '([i Row]))
                              empty-env 99))
                  (list-args ((analyze-multi-prior
                               'σ 'i '(normal 0 1) data '([i Row]))
                              empty-env 99))))
  )

;; (q . = . e)
;; Evaluate e and bind it to q in the dynamic environment
;; RG: Introduces notion of ExprInstr : Env -> Value
(define (analyze-single-binding q e data)
  (define einstr (analyze-expr e data))
  (λ (env lcsf)
    (values (extend-env env q (einstr env))
            lcsf)))

(module+ test
  (let ([old-lcsf 0.3])
    (define data (make-env '(μ σ) '(7 12)))
    (check-equal?
     (let-values ([(env lcsf)
                   ((analyze-single-binding 'ρ '(+ μ σ) data)
                    empty-env old-lcsf)])
       (list env lcsf))
     (list (extend-env empty-env 'ρ 19) old-lcsf))))


;; ((v i) . = . e)
(define (analyze-multi-binding v i e data tdecl)
  (define indices (get-indices-for-metavariable i (get-row-indices data) tdecl))
  (define instrs
    (for/list ([idx indices])
      (analyze-single-binding `(,v ,idx)  e (extend-env data i idx))))
  (sequence-instrs instrs))

(module+ test
  (let ([data (extend-env empty-env 'h #(1 2 3))])
    (check-equal? (list-args
                   ((analyze-multi-binding 'σ 'i 9 data '([i Row]))
                    empty-env 99))
                  '((((σ 2) 9) ((σ 1) 9) ((σ 0) 9))   99))
    (check-equal? (list-args
                   ((analyze-multi-binding 'σ 'i 'i data '([i Row]))
                    empty-env 99))
                  '((((σ 2) 2) ((σ 1) 1) ((σ 0) 0))   99))
    ))


;; (q . ~ . p) 
(define (analyze-single-prior q p data)
  (define pinstr (analyze-prior q p data))
  (λ (env lcsf) (values env
                        (+ lcsf (pinstr env)))))

(module+ test
  (let ([old-lcsf 0.3])
    (define data (make-env '(μ σ) '(7 12)))
    (check-equal?
     (let-values ([(env lcsf)
                   ((analyze-single-prior 'μ '(normal 0 1)
                                          data) empty-env old-lcsf)])
       (list env lcsf))
     (list empty-env (+ old-lcsf (pdf (normal-dist 0 1) 7 true))))))

;;
;; A probably-broken attempt to make a log-normal distribution by hand
;;
(define (log-normal-pdf npdf)
  (λ (in . args)
    (apply npdf (cons (log in) args))))

(define (log-normal-sample nsmp)
  (λ args
    (cond
      [(empty? args) (exp (apply nsmp args))]
      [else (map exp (apply nsmp args))])))

(define (log-normal-dist mean stddev)
  (define nd (normal-dist mean stddev))
  (distribution (log-normal-pdf (distribution-pdf nd))
                (log-normal-sample (distribution-sample nd))))

  


;; Expr Data -> PriorInstr
(define (analyze-prior q p data)
  (define q-inst (analyze-ref q data)) ;; look up q's value
  ;; RG - looks like the dist-constructors could be factored out and the args
  ;; map-analyzed, except for multivariate-normal
  (define dist-fn
    (match p
      [`(normal ,e1 ,e2)
       (let ([i1 (analyze-expr e1 data)]
             [i2 (analyze-expr e2 data)])
         (λ (env) (normal-dist (i1 env) (i2 env))))]
      [`(log-normal ,e1 ,e2)
       (let ([i1 (analyze-expr e1 data)]
             [i2 (analyze-expr e2 data)])
         ;; RG: Implementation of log-normal is sketchy!
         (λ (env) (log-normal-dist (i1 env) (i2 env))))]
      [`(cauchy ,e1 ,e2)
       (let ([i1 (analyze-expr e1 data)]
             [i2 (analyze-expr e2 data)])
         (λ (env) (cauchy-dist (i1 env) (i2 env))))]
      [`(laplace ,e)
       (let ([i (analyze-expr e data)])
         (λ (env) (... (i env))))]
      ;; RG : pairs of discrete-value and discrete-prob
      [`(discrete (,e1* ,n*) ...)
       (λ (env) (... #f))] 
      [`(exponential ,e)
       (let ([i (analyze-expr e data)])
         (λ (env) (exponential-dist (i env))))]
      [`(binomial ,e1 ,e2)
       (let ([i1 (analyze-expr e1 data)]
             [i2 (analyze-expr e2 data)])
         (λ (env) (binomial-dist (i1 env) (i2 env))))]
      [`(uniform ,e1 ,e2)
       (let ([i1 (analyze-expr e1 data)]
             [i2 (analyze-expr e2 data)])
         (λ (env) (uniform-dist (i1 env) (i2 env))))]
      ;; RG : this one's somewhat special...does it go in this function?
      [`(multivariate-normal #(,e1 ,e1* ...) #(#(,e2* ,e2** ...)))
       (λ (env) #f)]))
  ;; calculate log-distribution
  (λ (env) (pdf (dist-fn env) (q-inst env) true)))

(module+ test
  (let ()
    (define data (make-env '(μ σ) '(7 12)))
    (check-equal? ((analyze-prior 'μ '(normal 0 1)
                                  data) empty-env)
                  (pdf (normal-dist 0 1) 7 true))
    (check-equal? ((analyze-prior 'μ '(normal (* 0 9) (* σ 1))
                                  data) empty-env)
                  (pdf (normal-dist 0 12) 7 true))

    (check-equal? ((analyze-prior 'μ '(cauchy (* 0 9) (* σ 1))
                                  data) empty-env)
                  (pdf (cauchy-dist 0 12) 7 true))
    (check-equal? ((analyze-prior 'μ '(uniform (* 0 9) (* σ 1))
                                  data) empty-env)
                  (pdf (uniform-dist 0 12) 7 true))
    (check-equal? ((analyze-prior 'μ '(binomial (* 0 9) (* σ 1))
                                  data) empty-env)
                  (pdf (binomial-dist 0 12) 7 true))
    (check-equal? ((analyze-prior 'μ '(exponential  (* σ 1))
                                  data) empty-env)
                  (pdf (exponential-dist 12) 7 true))    
    ))


;; RG : merge this with analyze-multi-binding?
;; ((v i) . ~ . p)
(define (analyze-multi-prior v i p data tdecl)
  (define indices (get-indices-for-metavariable i (get-row-indices data) tdecl))
  (define instrs
    (for/list ([idx indices])
      (analyze-single-prior `(,v ,idx)  p (extend-env data i idx))))
  (sequence-instrs instrs))

(module+ test
  (check-equal?
   (list-args ((analyze-multi-prior 'μ 'i '(normal i 1)
                                    (make-env '(x) '(#(0 1 2)))
                                    '([i Row]))
               (make-env '(μ) '(#(9 8 7))) 12))
   (list (make-env '(μ) '(#(9 8 7)))
         (+ 12
            (pdf (normal-dist 0 1) 9 true)
            (pdf (normal-dist 1 1) 8 true)
            (pdf (normal-dist 2 1) 7 true))))                          
  )

;; Expr Data -> ExprInstr
;; translate the given expression (in the context of fit data) into an instr
(define (analyze-expr e data)
  (match e
    [`,n #:when (number? n) (λ (env) n)]
    [`(quote ,s) #:when (symbol? s) (λ (env) s)]  
    [`(+ ,e1 ,e2)
     (let ([i1 (analyze-expr e1 data)]
           [i2 (analyze-expr e2 data)])
       (λ (env) (+ (i1 env) (i2 env))))]
    [`(* ,e1 ,e2)
     (let ([i1 (analyze-expr e1 data)]
           [i2 (analyze-expr e2 data)])
       (λ (env) (* (i1 env) (i2 env))))]
    [`(expt ,e1 ,e2)
     (let ([i1 (analyze-expr e1 data)]
           [i2 (analyze-expr e2 data)])
       (λ (env) (expt (i1 env) (i2 env))))]
    [`(exp ,e)
     (let ([i (analyze-expr e data)])
       (λ (env) (exp (i env))))]
    [`,v #:when (symbol? v) (analyze-ref v data)]
    [`(,v ,e)
     #:when (symbol? v)
     (let ([i (analyze-expr e data)]       
           [idx-i (analyze-indexed-ref v data)])
       (λ (env) ;; evaluate the index, then do the lookup
         (idx-i env (i env))))]
    [`#(,e ,e* ...)
     (let ([i (analyze-expr e data)]
           [i* (for/list ([e^ e*]) (analyze-expr e^ data))])
       (λ (env) `#(,(i env) ,@(for/list ([i^ i*]) (i^ env)))))]))

(module+ test
  (let ()
    (define (ee expr) ((analyze-expr expr empty-env) empty-env))
    (check-equal?  (ee 5) 5)
    (check-equal? (ee '(+ 5 6)) 11)
    (check-equal? (ee '(* 5 6)) 30)
    (check-equal? (ee '(expt 5 6)) 15625)
    (check-equal? (ee '(exp 5)) (exp 5))
    (check-equal? (ee '(* (+ 5 6) (+ 5 6))) (* 11 11))
    (check-equal? (ee '(+ (* 5 6) (expt 5 6))) 15655)
    (check-equal?
     ((analyze-expr '(μ (+ 1 1)) (make-env '(μ) '(#(7 8 9)))) empty-env)
     9)
    (check-equal?
     ((analyze-expr '(μ (+ 1 1)) empty-env) (make-env '((μ 2)) '(9)))
     9)
    (let ([data (make-env '(μ) '(#(7 8 9)))]
          [env (make-env '((σ 9)) '(22))])
      (check-equal?
       ((analyze-expr '(σ (μ (+ 1 1))) data) env)
       22))
    (check-equal?
     ((analyze-expr '(μ 'female) empty-env) (make-env '((μ female)) '(9)))
     9)
    (check-exn exn:fail?
               (λ () ((analyze-expr '(μ female) empty-env)
                      (make-env '((μ female)) '(9)))))
    ))


;; Variable Data -> (Env Idx -> Value)
;; prepare to look up an indexed variable value in data table or environment
;; Note: Results in a number-indexed instr!
(define (analyze-indexed-ref v data)
  (with-handlers ([exn:fail? ;; Not in data table: look up in env
                   (λ (exn)                              
                     (λ (env idx) (lookup-value env `(,v ,idx))))])
    (let ([v* (lookup-family data v)])
      (λ (env idx)
        ;; found in data table, just return value at appropriate index
        (vector-ref v* idx)))))

(module+ test
  (let ([data (make-env '(μ) '(#(7 8 9)))])
    (check-equal? ((analyze-indexed-ref 'μ data) empty-env 0) 7)
    (check-equal? ((analyze-indexed-ref 'σ empty-env)
                   (extend-env empty-env '(σ female) 2) 'female)
                  2)))


;; look up the variable's value either in the data table or in the environment
(define (analyze-ref q data)
  (with-handlers ([exn:fail? ;; Not in data table: look up in env
                   (λ (exn) (λ (env) (lookup-value env q)))])
    (let ([v (lookup-value data q)])
      ;; found in data table, just return its value
      (λ (env) v))))

(module+ test
  (check-equal? ((analyze-ref 'a (extend-env empty-env 'a 10)) empty-env) 10)
  (check-equal? ((analyze-ref 'a empty-env) (extend-env empty-env 'a 5)) 5))


;; Linearize a bunch of closure-represented instructions
(define (sequence-instrs instr*)
  (let loop ([instr* instr*])
    (cond
      [(empty? instr*) (λ (env lcsf) (values env lcsf))]
      [else
       (let ([first-instr (first instr*)]
             [rest-instr* (loop (rest instr*))])
         ;; plumb the functions together
         (λ (env lcsf)
           (call-with-values (λ () (first-instr env lcsf))
                             rest-instr*)))])))

(module+ test
  (let ()
    (define i1 (λ (env lcsf) (values (extend-env env 'one 1) (+ lcsf 1))))
    (define i2 (λ (env lcsf) (values (extend-env env 'two 2) (+ lcsf 2))))
    (define i3 (λ (env lcsf) (values (extend-env env 'three 3) (+ lcsf 3))))
    (define i* (sequence-instrs (list i1 i2 i3)))
    (check-equal?
     (let-values ([(env lcsf) (i* '() 0)])
       (list (for/list ([v '(three two one)]) (list v (lookup-value env v)))
             lcsf))
     (list '((three 3) (two 2) (one 1)) 6))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Model Operations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Model -> (listOf VariableScheme)
;; get the schemes of the variables in the model
(define (get-tdecls model)
  (match model
    [`(model
       [type-decls ,tdecls ...]
       [var-decls  [,r* ,dr*] ...]
       [var-defs   ,vdef* ...])
     tdecls]))

(module+ test
  (check-equal? (get-tdecls ex1) '([i Row]))
  (check-equal? (get-tdecls ex4) '([i Row] [j Sex (Enum 'male 'female)])))

;; Model -> (listOf VariableScheme)
;; get the schemes of the variables in the model
(define (get-variable-schemes model)
  (match model
    [`(model
       [type-decls [,i* . ,t*] ...]
       [var-decls  [,r* ,dr*] ...]
       [var-defs   ,vdef* ...])
     r*]))

(module+ test
  (check-equal? (get-variable-schemes ex1) '((h i) μ σ)))


;; Model -> (listof VariableName)
;; get the names of variables in the model
(define (get-variable-names model)
  (map scheme->name (get-variable-schemes model)))

(module+ test
  (check-equal? (get-variable-names ex1) '(h μ σ)))


;; Data -> (listof Natural)
;; Row indices are determined by the number of rows in the fit data
(define (get-row-indices data)
  ;; row indices are just the length of every data row (so pick the first)
  (range (vector-length (second (first data)))))


;; 
(define (get-indices-for-metavariable index row-indices tdecls)
  (match (assq index tdecls)
    [`[,i Row] row-indices]
    [`(,i ,S (Enum ,s* ...)) `(,s* ...)]))


;; Model -> (listof VarDef)
;; get the variable definitions
(define (get-variable-defs model)
  (match model
    [`(model
       [type-decls [,i* . ,t*] ...]
       [var-decls  [,r* ,dr*] ...]
       [var-defs   ,vdef* ...])
     vdef*]))

(module+ test
  (check-equal? (get-variable-defs ex1)
                '([(h i) . ~ . (normal μ σ)]
                  [μ     . ~ . (normal 178 20)]
                  [σ     . ~ . (uniform 0 50)])))


;; VariableScheme -> VariableName
;; get the variable name from a variable scheme
(define (scheme->name r [fn 'scheme->name])
  (match r
    [`,v #:when (symbol? v) v]
    [`(,v ,i) #:when (symbol? v) v]
    [`,otherwise (error fn "Bad variable scheme ~a" otherwise)]))

(module+ test
  (check-equal? (scheme->name 'x) 'x)
  (check-equal? (scheme->name '(x i) 'parent-fn) 'x)
  (check-exn exn:fail? (λ () (scheme->name '(foo bar baz)))))


;; Model -> (listof VariablePattern)
;; get the variables defined using "="
(define (get-derived-variables model)
  (define vdef* (get-variable-defs model))
  (filter-map (λ (vdef)
                (match vdef
                  [`(,r . = . ,e) (scheme->name r)]
                  [`(,r . ~ . ,e) #f]
                  [`,otherwise
                   (error 'get-derived-variables "Bad variable definition: ~a"
                          otherwise)]))
              vdef*))

(module+ test
  (check-equal? (get-derived-variables ex1) '())
  (check-equal? (get-derived-variables ex2) '(μ)))


;; Model -> (listof VariableName)
;; the variables defined using "~"
(define (get-original-variables model)
  (define vdef* (get-variable-defs model))
  (filter-map (λ (vdef)
                (match vdef
                  [`(,r . = . ,e) #f]
                  [`(,r . ~ . ,e) (scheme->name r)]
                  [`,otherwise
                   (error 'get-original-variables "Bad variable definition: ~n"
                          otherwise)]))
              vdef*))

(module+ test
  (check-equal? (get-original-variables ex1) '(h μ σ))
  (check-equal? (get-original-variables ex2) '(y β σ x)))


;; Model Data -> (listof VariableName)
;; get the names of the observed variables that will be used to
;; condition the rest of the model. Ignore any extra data
(define (get-observed-variables model data)
  (define data-vars (env-refs data))
  (define model-vars (get-original-variables model))
  (filter (λ (dv) (member dv model-vars)) data-vars))

(module+ test
  (check-equal? (get-observed-variables ex1
                                        (make-env (list 'h)
                                                  (list (vector 1 2 3))))
                '(h)))


;; get the names of the unobserved variables (i.e., the rest)
(define (get-unobserved-variables model data)
  (define all-vars (get-variable-names model))
  (define observed-vars (get-observed-variables model data))
  (filter (λ (v) (not (member v observed-vars))) all-vars))

(module+ test
  (check-equal? (get-unobserved-variables ex1
                                          (make-env (list 'h)
                                                    (list (vector 1 2 3))))
                '(μ σ))

  (check-equal? (get-unobserved-variables ex2
                                          (make-env (list 'y 'x)
                                                    (list (vector 1 2 3)
                                                          (vector 4 5 6))))
                '(μ β σ)))


;; Model Data -> VariableName
;; get the names of the target variables of inference
(define (get-target-variables model data)
  (define unobserved-vars (get-unobserved-variables model data))
  (define derived-vars (get-derived-variables model))
  (filter (λ (v) (not (member v derived-vars))) unobserved-vars))

(module+ test
  (check-equal? (get-target-variables ex1
                                      (make-env (list 'h)
                                                (list (vector 1 2 3))))
                '(μ σ))

  (check-equal? (get-target-variables ex2
                                      (make-env (list 'y 'x)
                                                (list (vector 1 2 3)
                                                      (vector 4 5 6))))
                '(β σ)))


;;
;;  Environments (internal and data-table)
;;


;; Ref is one of:
;; - Symbol               (atomic variable ref)
;; - Symbol               (family variable name)
;; - (list Symbol Index)  (family variable ref)

;; Binding is one of:
;; - Value                (atomic binding)
;; - (vectorOf Value)     (family binding)
;; For now we'll be opaque about values, but so far they're numbers and symbols


;; Env is (listof (list Ref Binding))
;; And environment will either map a "real reference" Ref to a Value
;; or a "synthetic" vector of references to a vector.
;; e.g.  if σ ↦ #(5 6 7), that's like (σ 0) ↦ 5, (σ 1) ↦ 6, (σ 2) ↦ 7
(module+ test
  (define env0
  `([σ #(5 6 7)] [(α male) 2.0] [(α female) 3.9]
                 [μ 9] [s #('male 'female 'male)])))

(define empty-env empty)

;; given a list of targets and values, produce an environment.
;; Will need to deal with indexed target parameters somehow.
;; (listof Symbol) (listof Value) -> Env
(define (make-env names bindings)
  (map list names bindings))

(module+ test
  (check-equal? (make-env '(a b c) (list #(1 2 3) 4 5))
                '((a #(1 2 3)) (b 4) (c 5))))

;; Env Ref Binding -> Env
(define (extend-env env ref binding)
  (cons (list ref binding) env))


;; RG - revisit the next one:  what do I *really* want?
;
;; Env -> (listof Ref)
;; produce the reference names of variables bound by the given environment
(define (env-refs env)
  (map first env))

;; Env Ref -> (vectorOf Value)
;; produce the binding associated with a variable
(define (lookup-family env var [context-fn 'lookup-family])
  (unless (symbol? var)
    (error context-fn "Bad variable family: ~a" var))
  (cond
    [(assoc var env) => (λ (b)
                          (unless (vector? (second b))
                            (error context-fn "Bad family binding: ~a" b))
                          (second b))]
    [else (error context-fn "Bad lookup: ~a" var)]))


(module+ test
  (check-equal? (lookup-family env0 'σ) #(5 6 7))
  (check-equal? (lookup-family env0 's) #('male 'female 'male))
  (check-exn exn:fail? (λ () (lookup-family env0 0 'α))))


;; Quilt indices
(define (index? x)
  (or (symbol? x)
      (natural? x)))


;; Quilt Values
(define (value? x)
  (or (symbol? x)
      (number? x)))


;; Env Ref -> Value
;; produce the value associated with a fully-resolved binding
(define (lookup-value env ref [context 'lookup-value])
  (define (bad-ref) (error context "Bad lookup: ~a" ref))
  (match ref
    [`,v #:when (symbol? v)
         (let ([b (lookup-env env v context)])
           (if (value? b)
               b
               (error context "Reference ~a not bound to non-value ~a" v b)))]
    [`(,v ,i) #:when (and (symbol? v) (index? i))
              (let ([b (with-handlers ([exn:fail? (λ args #f)])
                         (lookup-env env v context))])
                (cond
                  ;; reference to a vector family
                  [(and (vector? b) (natural? i)
                        (< i (vector-length b)))
                   (vector-ref b i)]
                  ;; direct reference
                  [else (lookup-env env `(,v ,i))]))]
    [`,else (bad-ref)]))

(module+ test
  (check-equal? (lookup-value env0 '(σ 1)) 6)
  (check-equal? (lookup-value env0 '(σ 0)) 5)
  (check-equal? (lookup-value env0 'μ) 9)
  (check-equal? (lookup-value env0 '(α female)) 3.9)
  
  (check-exn exn:fail? (λ () (lookup-value env0 '(σ 3))))
  (check-exn exn:fail? (λ () (lookup-value env0 'σ)))
  (check-exn exn:fail? (λ () (lookup-value env0 'α)))
  (check-exn exn:fail? (λ () (lookup-value env0 'x))))
    

;; Env Ref -> Value or (vectorOf Value)
;; produce the binding associated with a variable, otherwise raise an error
(define (lookup-env env ref [context 'lookup-env])
  (cond
    [(assoc ref env) => second]
    [else (error context "Bad lookup: ~a" ref)]))
