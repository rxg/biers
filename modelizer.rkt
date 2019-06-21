#lang racket

;(require test-engine/racket-tests) ;; for check-expect
(require rackunit)

;;
;; modelizer.rkt - Some tools for approximately fitting
;;                 Bayesian statistical models
;;
;; Author: Ron Garcia <rxg@cs.ubc.ca>
;; Date: February 22, 2019
;;


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
;;        | [r . = . e]              ;; definition 

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

;; Examples:


;; 1) Linear Regression
;; Inference will expect to get height and weight inputs.
;; α, β, and σ are subject to inference
(define example1
  '(model
    [type-decls [i Row]]
    [var-decls  [(h i) Number] [(μ i) Number] [α Number] [β Number]
                [(w i) Number] [σ Number]]
    [var-defs   [(h i) . ~ . (normal (μ i) σ)]
                [(μ i) . = . (+ α (* β (w i)))]
                [α     . ~ . (normal 0 1)]
                [β     . ~ . (normal 0 1)]
                [(w i) . ~ . (normal 80 5)]
                [σ     . ~ . (log-normal 100 100)]]))


;; 2) Stratified Linear Regression Model
;; Inference will expect to get height and weight inputs.
;; α, β, and σ are subject to inference
(define example2
  '(model
    [type-decls [i Row] [m  Sex (Enum 'male 'female)]]
    [var-decls [(h i) Number] [(μ i) Number] [(α m) Number] [β Number]
               [(w i) Number] [(g i) Sex] [σ Number]]
    [var-defs   [(h i) . ~ . (normal (μ i) σ)]
                [(μ i) . = . (+ (α (g i)) (* β (w i)))]
                [(α m) . ~ . (normal 0 1)]
                [β     . ~ . (normal 0 1)]
                [(w i) . ~ . (normal 80 5)]
                [(g i) . ~ . (discrete ('male 1) ('female 1))]
                [σ     . ~ . (log-normal 100 100)]]))


;; Operations:
;; (draw-samples model [#:number n] variable-name ...)
;; variables should be in the same "tier" of variable
;; 

;; (fit-model model data) ;; what form should "data" take?


;; Some examples taken from McElreath's quap implementation:

;;;;;;;;;;
;; flist0 <- list(
;;     dist ~ dnorm( mean=a+b*speed , sd=sigma )
;; )
;; linear model integrated right into the likelihood
(define flist0a
  '(model
    [types     [i Row]]
    [var-decls [(dist i) Number] [α Number] [β Number]
               [(speed i) Number] [σ Number]]
    [var-defs  [(dist i)  . ~ . (normal (+ α (* β (speed i)) σ))]
               [α         . ~ . (normal 0 1)]
               [β         . ~ . (normal 0 1)]
               [(speed i) . ~ . (normal 80 5)]
               [σ         . ~ . (log-normal 100 100)]]))


;;;;;;;;;;
;;flist0 <- list(
;;    dist ~ dnorm( mean=mu , sd=sigma ) ,
;;    mu ~ a+b*speed
;;)
;; linear model defined separately
(define flist0b
  '(model
    [types     [i Row]]
    [var-decls [(dist i) Number] [α Number] [β Number]
               [(speed i) Number] [σ Number]]
    [var-defs  [(dist i)  . ~ . (normal (μ i) σ)]
               [(μ i)     . = . (+ α (* β (speed i)))]
               [α         . ~ . (normal 0 1)]
               [β         . ~ . (normal 0 1)]
               [(speed i) . ~ . (normal 80 5)]
               [σ         . ~ . (log-normal 100 100)]]))

;;;;;;;;;;
;;flist1 <- list(
;;    dist ~ dnorm( mean=a+b*speed , sd=sigma ) ,
;;    b ~ dnorm(0,1) ,
;;    sigma ~ dcauchy(0,1)
;;)
;;
(define flist1a
  '(model
    [types     [i Row]]
    [var-decls [(dist i) Number] [α Number] [β Number]
               [(speed i) Number] [σ Number]]
    [var-defs  [(dist i)  . ~ . (normal (+ α (* β (speed i))) σ)]
               [α         . ~ . (normal 0 1)]
               [β         . ~ . (normal 0 1)]
               [(speed i) . ~ . (normal 80 5)]
               [σ         . ~ . (cauchy 0 1)]]))


;;;;;;;;;;
;;flist1 <- list(
;;    dist ~ dnorm( mean=mu , sd=sigma ) ,
;;    mu ~ a+b*speed ,
;;    b ~ dnorm(0,1) ,
;;    sigma ~ dcauchy(0,1)
;;)
;;
(define flist1b
  '(model
    [types     [i Row]]
    [var-decls [(dist i) Number] [α Number] [β Number]
               [(speed i) Number] [σ Number]]
    [var-defs  [(dist i)  . ~ . (normal (μ i) σ)]
               [(μ i)     . = . (+ α (* β (speed i)))]
               [α         . ~ . (normal 0 1)]
               [β         . ~ . (normal 0 1)]
               [(speed i) . ~ . (normal 80 5)]
               [σ         . ~ . (cauchy 0 1)]]))


;;;;;;;;;;
;;flist2 <- list(
;;    dist ~ dnorm( mean=a+b*speed , sd=sigma ) ,
;;    c(a,b) ~ dnorm(0,10) , 
;;    sigma ~ dcauchy(0,1)
;;)
;; use vector to assign a and b the same distribution (hurm...)
(define flist2
  '(model
    [types     [i Row]]
    [var-decls [(dist i) Number] [α Number] [β Number]
               [(speed i) Number] [σ Number]]
    [var-defs  [(dist i)  . ~ . (normal (+ α (* β (speed i))) σ)]
               [α         . ~ . (normal 0 10)]
               [β         . ~ . (normal 0 10)]
               [(speed i) . ~ . (normal 80 5)]
               [σ         . ~ . (cauchy 0 1)]]))


;;;;;;;;;;
;;flist3 <- list(
;;    dist ~ dnorm( mean=a+b*speed , sd=sigma ) ,
;;    b ~ dlaplace(1) , 
;;    sigma ~ dcauchy(0,1)
;;)
;;
(define flist3
  '(model
    [types     [i Row]]
    [var-decls [(dist i) Number] [α Number] [β Number]
               [(speed i) Number] [σ Number]]
    [var-defs  [(dist i)  . ~ . (normal (+ α (* β (speed i))) σ)]
               [α         . ~ . (normal 0 1)]
               [β         . ~ . (laplace 1)]
               [(speed i) . ~ . (normal 80 5)]
               [σ         . ~ . (cauchy 0 1)]]))

;; Example of fitting
;;fit <- map( flist1,
;;            start=list(a=40,b=0.1,sigma=20),
;;            data=cars,
;;            debug=FALSE )
;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Everything below here need not be supported, since (if I understand
;; correctly) quadratic approximation is not good for any of it.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;#########
;;
;;library(rethinking)
;;data(chimpanzees)
;;

;;;;;;;;;;
;;flist1 <- list(
;;    pulled.left ~ dbinom( prob=logistic( a + b*prosoc.left ) , size=1 ),
;;    c(a,b) ~ dnorm(0,1)
;;)
;;

;;;;;;;;;;
;;flist2 <- list(
;;    pulled.left ~ dbinom( size=1 , prob=logistic( a + b*prosoc.left ) ),
;;    b ~ dnorm(0,1)
;;)
;;

;;;;;;;;;;
;;flist4 <- alist(
;;    pulled.left ~ dbinom( prob=p , size=1 ),
;;    logit(p) <- a + b*prosoc.left ,
;;    c(a,b) ~ dnorm(0,1)
;;)
;;

;;;;;;;;;;
;;fit2 <- map( flist4 , data=chimpanzees , start=list(a=0,b=0) , debug=FALSE )
;;


;;########
;;# regularized logistic regression example
;;y <- c( rep(0,10) , rep(1,10) )
;;x <- c( rep(-1,9) , rep(1,11) )
;;

;;;;;;;;;;
;;flist0 <- list(
;;    y ~ dbinom( prob=logistic( a + b*x ) , size=1 )
;;)
;;

;;;;;;;;;;
;;flist1 <- list(
;;    y ~ dbinom( prob=logistic( a + b*x ) , size=1 ),
;;    c(a,b) ~ dnorm(0,10)
;;)
;;

;;;;;;;;;;
;;fit3a <- map( flist0 , data=list(y=y,x=x) , start=list(a=0,b=0) )
;;

;;;;;;;;;;
;;fit3b <- map( flist1 , data=list(y=y,x=x) , start=list(a=0,b=0) )
;;
;;plot( y ~ x )
;;p <- sample.naive.posterior(fit3b)
;;xseq <- seq(-1,1,length.out=20)
;;pi.mu <- sapply( xseq , function(x) mean(logistic(p$a+p$b*x)) )
;;pi.ci <- sapply( xseq , function(x) PCI(logistic(p$a+p$b*x)) )
;;lines( xseq , pi.mu )
;;shade( pi.ci , xseq )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Let's code!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;
;; make-log-compat-fn:
;;   Given a model and data to fit it tp, produce a function
;;   that maps parameter values to relative log compatibility values

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

;; 
#;(make-log-compat-fn example-model)

;; !!!
;; Given a model and some data to fit to it, create a log-compatible function
;; suitable for quadratic approximation
(define (make-log-compat-fn model data)
  (define variables (get-variable-schemes model))
  (define index (get-row-index model))
  (define tdecls (get-tdecls model))
  (define var-defs (get-variable-defs model))
  (define observed (get-observed-variables model data))
  (define targets (get-target-variables model data))
  ;; targets gives us an ordering for values

  ;; create a function from targets to Real
  ;; given a set of targets, compute the log-compatibility
  ;; of all of the targets, given the observations
  ;; note that the observations are indexed on Row
  ;; (that's how we know how many observation variables there are)

  ;; To compute:
  ;; going from last-to-first, compute log-likelihoods
  ;; - target variable priors just turn into log-probabilities given
  ;;   the provided target value. Contribute no new bindings.
  ;; - non-target parameters (i.e., equations) just turn into environment
  ;;   bindings for use in "earlier" priors.  May form a "vector" or otherwise
  ;;   indexed set of data. Contribute nothing to log probability
  ;; - observed variable priors also turn into log-probabilities, though
  ;;   splayed across the Row observations (so there's an implicit map)
  ;;   contribute no new bindings.
  ;;
  ;;   So evaluation is essentially a "foldr" that produces a pair of an
  ;;   environment-and-cumulative-log-joint-compatibility value (summed up).

  ;;  Since we interpret vdefs last-to-first, we can foldr and connect them
  ;;  backwards!
  (define inner-lc-fn
    ;; accumulator llsf is the log-likelihood so far, starts at 0
    ;; env is any new derived parameter values, like a linear model
    ;; - at construction time we can determine what needs to be looked up
    ;;   in the env (anything that's not in data)
    (let loop ([vdef* var-defs])
      (cond
        [(empty? vdef*) (λ (env llsf) (values env llsf))]
        [else
         (let ([this-fn (analyze-vdef (first vdef*) data)]
               [rest-fn (loop (rest vdef*))])
           ;; plumb the functions together
           ;; (note:  if I foldl this, I can call directly w/o values
           (λ (env llsf)
             (call-with-values (λ () (rest-fn env llsf))
                               this-fn)))])))
  
  (define (lc-fn tgt-values)
    (define env0 (make-env targets tgt-values))
    (define likelihood0 0)
    #f)
  lc-fn)


;; !!!
;; Turn a vdef and data table into an inner ll fn which
;; threads internal defs and log-likelihood-so-far.
(define (analyze-vdef vdef data)
  (match vdef
    ;; RG : implement these!
    [`(,r . = . ,e) (λ (env llsf) (values env llsf))]
    [`(,r . ~ . ,e) (λ (env llsf) (values env llsf))]
    [`,otherwise
     (error 'get-derived-variables "Bad variable definition: ~a"
            otherwise)]))

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

;; Model -> Symbol
;; get the name of the row index variable.  There must be one
(define (get-row-index model)
  (match model
    [`(model
       [type-decls [,i* . ,t*] ...]
       [var-decls  [,r* ,dr*] ...]
       [var-defs   ,vdef* ...])
     ;; Find the first declaration with "type" '(Row)
     (let ([result (assoc '(Row) (map cons t* i*))])
       (if result
           (cdr result)
           (error 'get-row-index "No Row index specified.")))]))

(module+ test
  (check-equal? (get-row-index ex1) 'i))


;; Helper:
;; (listof X) -> (listof X)
;; return a list of the unique elements (in no certain order)
(define (uniquify-list lst)
  (set->list (list->set lst)))


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
(define env0
  `([σ #(5 6 7)] [(α male) 2.0] [(α female) 3.9]
                 [μ 9] [s #('male 'female 'male)]))


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


;; Indices in Quilt
(define (index? x)
  (or (symbol? x)
      (natural? x)))


;; Values in Quilt
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
                ;; reference to a vector family
                (cond
                  [(and (vector? b) (natural? i)
                        (< i (vector-length b)))
                   (vector-ref b i)]
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
  


#;
(define env0
  `([σ #(5 6 7)] [(α male) 2.0] [(α female) 3.9]
                 [μ 9] [s #('male 'female 'male)]))
    

;; Env Ref -> Value or (vectorOf Value)
;; produce the binding associated with a variable
(define (lookup-env env ref [context 'lookup-env])
  (cond
    [(assoc ref env) => second]
    [else (error context "Bad lookup: ~a" ref)]))
