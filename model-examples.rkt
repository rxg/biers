#lang racket

;;
;; model-examples.rkt - some example Modelizer models, taken from McElreath
;;

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

