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
;; i ∈ Index
;; r ∈ VariablePattern
;; n ∈ Number
;; m ∈ Model
;; e ∈ Expression
;; d ∈ DataType
;; t ∈ TypeSpec
;; s ∈ Symbol
;; p ∈ Prior

;; t ::= (Index)             ;; metavariable index
;;     | (d (Enum s ...))    ;; type name and structure

;; m ::=
;; (new-model
;;   [type-decls [i . t] ...]
;;   [var-decls [r d] ...]
;;   [var-defs vdef ...]
;;   [data ...])   ;; really want data to be optional

;; vdef ::= [r . ~ . p]              ;; prior distribution
;;        | [r . = . e]              ;; definition 

;; r ::= v           ;; atomic variable
;;     | (v d)       ;; indexed family of variables
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


;; Notes:
;; - dependency order for var-defs is last-to-first, to match stat models.
;;   think of it like let* in reverse
;;   For now, we will not attempt to re-order to match dependencies.
;; - data type and variable names must be disjoint

;; Examples:

;; 0) Normal Model (0-ary Linear Regression)

;; (new-model
;;   [type-decls [i Index]]
;;   [var-decls  [(h i) Number] [μ Number] [σ Number]]
;;   [var-defs   [(h i) . ~ . (normal (μ σ)]
;;               [μ     . ~ . (normal 178 20)]
;;               [σ     . ~ . (uniform 0 50)]])

;; 1) Linear Regression

;; (new-model
;;   [type-decls [i Index]]
;;   [var-decls  [(h i) Number] [(μ i) Number] [α Number] [β Number]
;;               [(w i) Number] [σ Number]]
;;   [var-defs   [(h i) . ~ . (normal (μ i) σ)]
;;               [(μ i) . = . (+ α (* β (w i))]
;;               [α     . ~ . (normal 0 1)]
;;               [β     . ~ . (normal 0 1)]
;;               [(w i) . ~ . (normal 80 5)]]
;;               [σ     . ~ . (log-normal 100 100)]])
;;
;; Inference will expect to get height and weight inputs.
;; α, β, and σ are subject to inference


;; 2) Stratified Linear Regression Model

;; (new-model
;;   [type-decls [i Index] [m  Sex (Enum 'male 'female)]]
;;   [var-decls [(h i) Number] [(μ i) Number] [(α m) Number] [β Number]
;;              [(w i) Number] [(g i) Sex] [σ Number]]
;;   [var-defs   [(h i) . ~ . (normal (μ i) σ)]
;;               [(μ i) . = . (+ (α (g i)) (* β (w i))]
;;               [(α m) . ~ . (normal 0 1)]
;;               [β     . ~ . (normal 0 1)]
;;               [(w i) . ~ . (normal 80 5)]]
;;               [(g i) . ~ . (discrete ('male 1) ('female 1))]
;;               [σ     . ~ . (log-normal 100 100)]])
;;
;; Inference will expect to get height and weight inputs.
;; α, β, and σ are subject to inference


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

;; (define flist0
;;   (new-model
;;     [types     [i Index]]
;;     [var-decls [(dist i) Number] [α Number] [β Number]
;;                [(speed i) Number] [σ Number]]
;;     [var-defs  [(dist i)  . ~ . (normal (+ α (* β (speed i)) σ)]
;;                [α         . ~ . (normal 0 1)]
;;                [β         . ~ . (normal 0 1)]
;;                [(speed i) . ~ . (normal 80 5)]]
;;                [σ         . ~ . (log-normal 100 100)]]))


;;
;;;;;;;;;;
;;flist0 <- list(
;;    dist ~ dnorm( mean=mu , sd=sigma ) ,
;;    mu ~ a+b*speed
;;)
;;

;; (define flist0
;;   (new-model
;;     [types     [i Index]]
;;     [var-decls [(dist i) Number] [α Number] [β Number]
;;                [(speed i) Number] [σ Number]]
;;     [var-defs  [(dist i)  . ~ . (normal (μ i) σ)]
;;                [(μ i)     . = . (+ α (* β (speed i))]
;;                [α         . ~ . (normal 0 1)]
;;                [β         . ~ . (normal 0 1)]
;;                [(speed i) . ~ . (normal 80 5)]]
;;                [σ         . ~ . (log-normal 100 100)]]))




;;;;;;;;;;
;;flist1 <- list(
;;    dist ~ dnorm( mean=a+b*speed , sd=sigma ) ,
;;    b ~ dnorm(0,1) ,
;;    sigma ~ dcauchy(0,1)
;;)
;;

;; (define flist1
;;   (new-model
;;     [types     [i Index]]
;;     [var-decls [(dist i) Number] [α Number] [β Number]
;;                [(speed i) Number] [σ Number]]
;;     [var-defs  [(dist i)  . ~ . (normal (+ α (* β (speed i)) σ)]
;;                [α         . ~ . (normal 0 1)]
;;                [β         . ~ . (normal 0 1)]
;;                [(speed i) . ~ . (normal 80 5)]]
;;                [σ         . ~ . (cauchy 0 1)]]))


;;;;;;;;;;
;;flist1 <- list(
;;    dist ~ dnorm( mean=mu , sd=sigma ) ,
;;    mu ~ a+b*speed ,
;;    b ~ dnorm(0,1) ,
;;    sigma ~ dcauchy(0,1)
;;)
;;

;; (define flist1
;;   (new-model
;;     [types     [i Index]]
;;     [var-decls [(dist i) Number] [α Number] [β Number]
;;                [(speed i) Number] [σ Number]]
;;     [var-defs  [(dist i)  . ~ . (normal (μ i) σ)]
;;                [(μ i)     . = . (+ α (* β (speed i))]
;;                [α         . ~ . (normal 0 1)]
;;                [β         . ~ . (normal 0 1)]
;;                [(speed i) . ~ . (normal 80 5)]]
;;                [σ         . ~ . (cauchy 0 1)]]))


;;;;;;;;;;
;;flist2 <- list(
;;    dist ~ dnorm( mean=a+b*speed , sd=sigma ) ,
;;    c(a,b) ~ dnorm(0,10) , 
;;    sigma ~ dcauchy(0,1)
;;)
;;



;;;;;;;;;;
;;flist3 <- list(
;;    dist ~ dnorm( mean=a+b*speed , sd=sigma ) ,
;;    b ~ dlaplace(1) , 
;;    sigma ~ dcauchy(0,1)
;;)
;;


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
;;   Given a model and data to fit to it, produce a function
;;   that maps parameter values to relative log compatibility values

;; Chapter 3: Gaussian model of Kalamari forager height
(define ex1
  `(new-model
    [type-decls [i Index]]
    [var-decls  [(h i) Number] [μ Number] [σ Number]]
    [var-defs   [(h i) . ~ . (normal μ σ)]
                [μ     . ~ . (normal 178 20)]
                [σ     . ~ . (uniform 0 50)]]))

;; Chapter 3: Schematic of a linear regression (with a prior on the predictor)
(define ex2
  `(new-model
    [type-decls [i Index]]
    [var-decls  [(y i) Number] [(μ i) Number] [β Number] [σ Number]
                [(x i) Number]]
    [var-defs   [(y i) . ~ . (normal μ σ)]
                [(μ i) . = . (* β (x i))]
                [β     . ~ . (normal 0 10)]
                [σ     . ~ . (exponential 1)]
                [(x i) . ~ . (normal 0 1)]]))

;; Indicator variable for male (0 = not male, 1 = male)
;; This is essentially Student's t-test encoded as a linear model,
;; without the decision (just the estimation)
(define ex3
  `(new-model
    [type-decls [i Index]]
    [var-decls  [(h i) Number] [(μ i) Number] [α Number] [βm Number]
                [σ Number] [(m i) Number]]
    [var-defs   [(h i) . ~ . (normal μ σ)]
                [(μ i) . = . (+ α (* βm (m i)))]
                [α     . ~ . (normal 178 20)]
                [βm    . ~ . (normal 0 1)]
                [σ     . ~ . (uniform 0 50)]
                [(m i) . ~ . (discrete (0 0.5) (1 0.5))]]))

;; Chapter ??: Index variable
(define ex4
  `(new-model
    [type-decls [i Index] [j Sex (Enum 'male 'female)]]
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
  (define variables (get-variable-names model))
  (define index (get-index model))
  (define var-defs (get-variable-defs model))
  (define observed (get-observed-variables model data))
  (define targets (get-target-variables model data))
  ;; targets gives us an ordering for values

  ;; create a function from targets to Real
  ;; given a set of targets, compute the log-compatibility
  ;; of all of the targets, given the observations
  ;; note that the observations are indexed on Index
  ;; (that's how we know how many observation variables there are)

  ;; To compute:
  ;; going from last-to-first, compute log-likelihoods
  ;; - target variable priors just turn into log-probabilities given
  ;;   the provided target value. Contribute no new bindings.
  ;; - non-target parameters (i.e., equations) just turn into environment
  ;;   bindings for use in "earlier" priors.  May form a "vector" or otherwise
  ;;   indexed set of data. Contribute nothing to log probability
  ;; - observed variable priors also turn into log-probabilities, though
  ;;   splayed across the Index observations (so there's an implicit map)
  ;;   contribute no new bindings.
  ;;
  ;;   So evaluation is essentially a "foldr" that produces a pair of an
  ;;   environment-and-cumulative-log-joint-compatibility value (summed up).

  (define (lc-fn tgt-values)
    (define env0 (make-env targets tgt-values))
    (define likelihood0 0)
    #f)
  lc-fn)

;; m ::=
;; (new-model
;;   [type-decls [i . t] ...]
;;   [var-decls [r d] ...]
;;   [var-defs vdef ...])


;; Model -> (listOf Symbol)
;; get the names of the parameters
(define (get-variable-names model)
  (match model
    [`(new-model
       [type-decls [,i* . ,t*] ...]
       [var-decls  [,r* ,dr*] ...]
       [var-defs   ,vdef* ...])
     r*]))

(module+ test
  (check-equal? (get-variable-names ex1) '((h i) μ σ)))

;; Model -> Symbol
;; get the name of the index variable.  There must be one
(define (get-index model)
  (match model
    [`(new-model
       [type-decls [,i* . ,t*] ...]
       [var-decls  [,r* ,dr*] ...]
       [var-defs   ,vdef* ...])
     ;; Find the first type declaration with type 'Index
     (let ([result (assoc '(Index) (map cons t* i*))])
       (if result
           (cdr result)
           (error 'get-index "No Index specified.")))]))

(module+ test
  (check-equal? (get-index ex1) 'i))


;; Model -> (listof VarDef)
;; get the variable definitions
(define (get-variable-defs model)
  (match model
    [`(new-model
       [type-decls [,i* . ,t*] ...]
       [var-decls  [,r* ,dr*] ...]
       [var-defs   ,vdef* ...])
     vdef*]))

(module+ test
  (check-equal? (get-variable-defs ex1)
                '([(h i) . ~ . (normal μ σ)]
                  [μ     . ~ . (normal 178 20)]
                  [σ     . ~ . (uniform 0 50)])))


;; Model -> (listof VariablePattern)
;; the variables defined using "="
(define (get-derived-variables model)
  (define vdef* (get-variable-defs model))
  (filter-map (λ (vdef)
                (match vdef
                  [`(,r . = . ,e) r]
                  [`(,r . ~ . ,e) #f]
                  [`,otherwise
                   (error 'get-derived-variables "Bad variable definition: ~n"
                          otherwise)]))
              vdef*))

(module+ test
  (check-equal? (get-derived-variables ex1) '())
  (check-equal? (get-derived-variables ex2) '((μ i))))


;; Model -> (listof VariablePattern)
;; the variables defined using "="
(define (get-original-variables model)
  (define vdef* (get-variable-defs model))
  (filter-map (λ (vdef)
                (match vdef
                  [`(,r . = . ,e) #f]
                  [`(,r . ~ . ,e) r]
                  [`,otherwise
                   (error 'get-original-variables "Bad variable definition: ~n"
                          otherwise)]))
              vdef*))

(module+ test
  (check-equal? (get-original-variables ex1) '((h i) μ σ))
  (check-equal? (get-original-variables ex2) '((y i) β σ (x i))))


;; !!!
;; get the names of the observed variables that will be used to
;; condition the rest of the model.
(define (get-observed-variables model data)
  ;; the observed variables are named in "data" though without the
  ;; index variable.  We should confirm that the data and model are in sync
  ;; in that regard.  (any original variable named "(v i)" where "i"
  ;; is the indexx should appear as a vector "v" in the data table (I think)
  empty)


;; !!!
;; get the names of the unobserved variables (i.e., the rest)
(define (get-unobserved-variables model data)
  ;; (get-variable-names model data)
  ;; minus
  ;; (get-observed-variables)
  empty)


;; TERMINOLOGY: for lack of a better word, I will use the name
;; "target variables" for the direct targets of inference, i.e., those
;; unobserved variables that will be directly conditioned with respect to the
;; data. This omits derived variables e.g., linear models.

;; !!!
;; get the names of the target variables of inference
(define (get-target-variables model data)
  ;; target variables are the unobserved variables that have associated
  ;; distributions (i.e., are NOT defined using "=")
  ;; (get-unobserved-variables model data)
  ;; minus
  ;; (get-derived-variables model)
  empty)


;; !!!
;; given a list of targets and values, produce an environment.
;; Will need to deal with indexed target parameters somehow.
(define (make-env names values) empty)