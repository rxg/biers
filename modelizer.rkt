#lang racket


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
;; vdef ∈ VariableDefinition
;; r ∈ Variable Pattern
;; n ∈ Number
;; m ∈ Model
;; e ∈ Expression
;; d ∈ DataType
;; t ∈ TypeSpec
;; s ∈ Symbol
;; p ∈ Prior

;; t ::= Index | Number | (Enum s ...) 

;; m ::=
;; (new-model
;;   [type-decls [d t] ...]
;;   [var-decls [r d] ...]
;;   [var-defs vdef ...]
;;   [data ...])

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
;;   [type-decls (i Index) (m (Enum 'male 'female))]
;;   [var-decls [(h i) Number] [(μ i) Number] [(α m) Number] [β Number]
;;              [(w i) Number] [(g i) m] [σ Number]]
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
;; make-log-plaus-fn:
;;   Given a model and data to fit to it, produce a function
;;   that maps parameter values to relative log plausibility values

(make-log-plaus-fn
 `(new-model
   [type-decls [i Index]]
   [var-decls  [(h i) Number] [μ Number] [σ Number]]
   [var-defs   [(h i) . ~ . (normal μ σ)]
               [μ     . ~ . (normal 178 20)]
               [σ     . ~ . (uniform 0 50)]]))

(define (make-log-plaus-fn model data)
  (define variables (get-variables model))
  (define index (get-index model))
  (define observed (get-observed-variables model data))
  (define targets (get-target-variables model data))
  #f)


;; Model -> (listOf Symbol)
;; get the names of the parameters
(define (get-variables model) empty)

;; Model -> Symbol
;; get the type name of the index variable.  There must be one
(define (get-index model) 'i)

;; get the names of the observed variables that will be used to
;; condition the rest of the model.
(define (get-observed-variables model data) empty)

;; get the names of the unobserved variables (i.e., the rest)
(define (get-unobserved-variables model data) empty)

;; TERMINOLOGY: for lack of a better word, I will use the name
;; "target variables" for the direct targets of inference, i.e., those
;; unobserved variables that will be directly conditioned with respect to the
;; data. This omits derived variables e.g., linear models.

;; get the names of the target variables of inference
(define (get-target-variables model data) empty)
