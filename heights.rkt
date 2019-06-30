#lang racket

(require plot)
(require math/base)
(require math/distributions)
(require "Racket_NumericalMethods/Numerical_Integration.rkt")
(require racket/vector) ;; for vector-take
(require (only-in "utils.rkt" with-arity-of))
;;
;; heights.rkt - Inference over height data
;;  Based on McElreath, Statistical Rethinking, Chapter 3
;;
;; Author: Ron Garcia <rxg@cs.ubc.ca>
;; Date: November, 2018
;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Stochastic simulator of individual heights
;;
(define (model)
  (let ([sigma  (sample (uniform-dist 0 50))]
        [mu (sample (normal-dist 178 20))])
    (let ([h_i (sample (normal-dist mu sigma))])
      h_i)))

;; Query the generative model to visualize the distribution

;; math/plot's density function produces and plots a continuous
;; density function that was estimated from a list of samples. See
;;   https://en.wikipedia.org/wiki/Density_estimation
;; for the general idea.
;;
;; I do not know the basis for this particular implementation, but based on
;; the docs, it appears to use Gaussian kernel smoothing:
;;   https://en.wikipedia.org/wiki/Kernel_density_estimation
;;
;; For comparison, See R's density function
;;   https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/density
;; as well as it's options for choosing a bandwidth:
;; https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/bandwidth
;;
;; The R manual page for density recommends calculating bandwidth using
;;     the method "SJ" of:
;;     Sheather, S. J. and Jones, M. C. (1991). A reliable data-based bandwidth
;;         selection method for kernel density estimation.
;;         Journal of the Royal Statistical Society series B, 53, 683--690.
;;         http://www.jstor.org/stable/2345597.
;;  though R's default, for historical purposes, is "Silverman's rule of thumb":
;;  Silverman, B. W. (1986). Density Estimation. London: Chapman and Hall.
;;
;; Racket appears to use a method related to R's bw.nrd.
#;
(plot (density (for/list ([i (in-range 1000)]) (model)) 2))





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Inference Engine

;; Random variables:
;; h* = (list h_i ...) : individual height samples (in cm)
;;      (assume independent, identically distributed (i.i.d.))
;; μ : mean height (in cm)
;; σ : standard deviation (in cm)

;; Statistical Model
;; I = sample size (not random)
;; { h_i ∼ Normal(μ,σ) } for i ∈ I (implies i.i.d.)
;;                     (mathematically models "regression to the mean"?
;;                      Regession to the mean is an *empirical* phenomenon that
;;                      models should account for to enable good inferences.)
;; μ ∼ Normal(178,20) (McElreath is 178cm tall
;;                     96% weight on the mean being in 178 ± 3·20 (so [118,238])
;;                     Normal prior on mean height weighs against extreme means
;;                     Data must be overwhelming to move the mean to extremes.)
;; σ ∼ Uniform(0,50)  (standard deviation of heights is SOMEWHERE in [0,50])
;;                     so "negligible" (but non-zero) probability density
;;                     on negative heights (i.e. impossible heights)
;;                     Now that we have computers, can we force it to zero?
;;                     Maybe! see Truncated Gaussian and Rectified Gaussian on
;;                     Wikipedia. How to use responsibly?)

;; We perform inference on all of the data, the individual values h_i.
;; This is going to be analogous to the urn model: we can do it all at once
;; since h_i are assumed to be i.i.d.

;; parameterized "merge densities" operation
;; when combining probabilities, use multiplication;
;; when combining log-probabilities, use addition
(define (merge-d [log? #f])
  (if log? + *))

;; prior (log-)density
;; prior is 2-dimensional: P(μ,σ|Γ) (Γ embodies initial knowledge)
;; Assume μ,σ independent: P(μ,σ|Γ) = P(μ|Γ)P(σ|Γ)
(define (prior-pdf μ σ [log? #f])
  (define cmb (if log? + *))
  (let ([σ-pdf (distribution-pdf (uniform-dist 0 50))]
        [μ-pdf (distribution-pdf (normal-dist 178 20))])
    (cmb (μ-pdf μ log?) (σ-pdf σ log?))))

;; Let's graph the prior pdf in 3D!
;; Include 4 std-dev for μ; extend σ past the support values
#;
(plot3d (surface3d prior-pdf 98 258 -5 55))


;; multi-observation (log-)likelihood (Γ makes the relativity of prior explicit)
;; P(h*|μ,σ,Γ) = Π_i P(h_i|μ,σ,Γ)
(define (aggregate-likelihood h* μ σ [log? #f])
  (define dist (normal-dist μ σ))
  ;; P(h_i,μ,σ|Γ) - single-observation likelihood
  (define (likelihood h_i) (pdf dist h_i log?))
  (define cmb (merge-d log?))
  (for/fold ([p-acc (sequence-ref h* 0)])
            ([h_i (sequence-tail h* 1)])
    (cmb (likelihood h_i) p-acc)))


;; joint (log-)probability density (modulo initial prior knowledge Γ)
;; P(h*μ,σ|Γ) = P(h*|μ,σ,Γ)P(μ,σ|Γ)
(define (nnjoint prior h* μ σ [log? #f])
  (define cmb (merge-d log?))
  (cmb (aggregate-likelihood h* μ σ log?) (prior μ σ log?)))

;; curried version of joint
(define (make-joint h* [log? #f])
  (λ (μ σ)
    (nnjoint prior-pdf h* μ σ log?)))
  

;;
;; Terrifyingly inefficient numerical integration-based posterior
;;  I can't even get a single value back from posterior... :(
;;  At least it heats up the room a little ;)
;;


;; double-integration?  No, not so much.
(define (marginal h*)
  (adaptive-integrate
   (λ (μ)
     (adaptive-integrate 
      (λ (σ)
        (nnjoint h* μ σ))
      0 50 .1)) ;; range of the uniform distribution
   120 240 .1)) ;; roughly 3σ range around the mean

(define (posterior μ σ h*)
  (/ (nnjoint h* μ σ)
     (marginal h*)))


;;
;; Grid Sampling
;;


;; Build a 1D grid of values
(define (seq lo hi size)
  (range lo hi (/ (- hi lo) size)))

;; grid log posterior as per R code 4.14 of McElreath
(define (fit-heights h* μ-lo μ-hi σ-lo σ-hi)
  (local
    ;; Build a 200 x 200 grid of μ,σ values
    [(define μ*-grid (seq μ-lo μ-hi 200))
     (define σ*-grid (seq σ-lo σ-hi 200))

     ;; I think it slows things down to build this up progressively:
     ;; Do it as much in aggregate as possible

     ;; list of `#(,μ ,σ) vector entries
     #;(define μσ* (for*/list ([μ μ*][σ σ*])
                     (vector μ σ)))
     ;; 
     #;
     (define ll* (for*/list ([μ μ*][σ σ*])
                   (let ([dist (distribution-pdf (normal-dist μ σ))])
                     (vector μ σ (for/sum ([d d*])
                                   (log (dist d)))))))

     (define prod*
       (for*/list ([μ (in-list μ*-grid)]
                   [σ (in-list σ*-grid)])
         (let ([dist (distribution-pdf (normal-dist μ σ))])
           (vector μ σ
                   (+ (for/sum ([d h*])
                        (log (dist d)))
                      (log (prior-pdf μ σ)))))))

     ;; Normalize so that the max is 1
     ;; This is not the "correct" normalization, but nice for visualization
     ;; 1 "should" be the volume under the curve, not the max value
     (define prob*
       (let ([max-prod (for/fold ([mx (vector-ref (first prod*) 2)])
                                 ([p (in-list (rest prod*))])
                         (max mx (vector-ref p 2)))])
         (for/list ([p (in-list prod*)])
           (vector (vector-ref p 0)
                   (vector-ref p 1)
                   (exp (- (vector-ref p 2) max-prod))))))]
    prob*))

(define (plot-μσ-grid-posterior prob*)
  (plot3d
   (points3d prob* #:size 1)
   #:title "Posterior Normal Distribution of Heights in cm (Howell Data)"
   #:x-label "Mean (μ)"
   #:y-label "Standard Deviation (σ)"
   #:z-label "Relative Plausibility"))

;; Example Fit
#;
(begin
  (require "howell.rkt")
  ;; grab the first 20 adult heights, just to get started
  (define first-20-heights (in-vector adult-heights 0 20))
  #;(define random-20-heights (take (shuffle (vector->list adult-heights)) 20))
  (define post-20 (fit-heights first-20-heights 150 170 4 20))
  (define post-adults (fit-heights adult-heights 153 156 6 9))
  (map plot-μσ-grid-posterior (list post-20 post-adults))
  ;; Sample from the posterior
  (define-values (μσ-20* w-20)
    (for/lists (uss w) ([p post-adults])
      (values (vector-take p 2) (vector-ref p 2))))
  (define post-dist-20 (discrete-dist μσ-20* w-20))
  (define s-20* (sample post-dist-20 10000))
  (void))

;; Example of drawing samples from the posterior
#;
(begin
  (require (submod #;"." "heights.rkt" howell))
  (define post (fit-heights adult-heights 153 156 6 9))
  (define-values (μσ* w)
    (for/lists (uss w) ([p post]) (values (vector-take p 2) (vector-ref p 2))))
  (define post-dist (discrete-dist μσ* w))
  (define s* (sample post-dist 10000))
  ;; "heat map" of samples
  (plot (points s* #:alpha 0.4 #:sym 'fullcircle1 #:color 'blue))
  ;; density plot of mu
  (plot (density (for/list ([s s*]) (vector-ref s 0))))
  ;; density plot of sigma
  (plot (density (for/list ([s s*]) (vector-ref s 1))))
  (void))

;;  Repeating the same, but with only 20 samples so as to inspect sigma
#;
(begin
  (require "howell.rkt")
  (define random-20-heights (take (shuffle (vector->list adult-heights)) 20))
  (define post (fit-heights random-20-heights 150 170 4 20))
  (define-values (μσ* w)
    (for/lists (uss w) ([p post]) (values (vector-take p 2) (vector-ref p 2))))
  (define post-dist (discrete-dist μσ* w))
  (define s* (sample post-dist 10000))
  (list
   ;; "heat map" of samples
   (plot (points s* #:alpha 0.4 #:sym 'fullcircle1 #:color 'blue))
   ;; density plot of mu
   (plot (density (for/list ([s s*]) (vector-ref s 0))))
   ;; density plot of sigma
   (plot (density (for/list ([s s*]) (vector-ref s 1)))))
  (void))


;;
;; Laplace Approximation
;;

#;
(begin
  (require "howell.rkt")           ; Howell Data
  (require "laplace-approx.rkt")   ; Laplace Approximation Engine
  (require "modelizer.rkt")        ; Interpret model specifications 
  (plot-height 300) (plot-width 300)
  ;;
  ;; Fit the first 20 heights
  ;;
  
  ;; grab the first 20 adult heights, just to get started, and help w/ testing
  (define first-20-heights (in-vector adult-heights 0 20))
  ;; for variety, select a random 20 adult heights
  #;(define random-20-heights (take (shuffle (vector->list adult-heights)) 20))
  ;; non-normalized log posterior (compositionally constructed)
  (define lnpf20 (make-joint first-20-heights 'log))

  (define mlnpf20 
    (make-log-compatibility-function
     `(model
       [type-decls [i Row]]
       [var-decls  [(h i) Number] [μ Number] [σ Number]]
       [var-defs   [(h i) . ~ . (normal μ σ)]
                   [μ     . ~ . (normal 178 20)]
                   [σ     . ~ . (uniform 0 50)]])
     (extend-env empty-env 'h first-20-heights)
     ))
    
  (define laf20 (laplace-approx lnpf20 '((140 160) (4 20))))
  (define laf20 (laplace-approx mlnpf20 '((140 160) (4 20))))
  (define qpdff20 (lapprox-pdf laf20))
  (define qp3df20
    (plot3d (surface3d qpdff20 140 160 4 20)
            #:title
            "Approximate Posterior Normal Distributions of Heights (n=20)"
            #:x-label "Mean (μ) in cm"
            #:y-label "Standard Deviation (σ) in cm"
            #:z-label "Probability"))

  ;; given non-normalized log posterior and mode,
  ;; produce unit-normalized posterior
  (define (make-updf nnlp la)
    (let ([mx (apply nnlp (lapprox-means la))])
      (with-arity-of nnlp
        (λ x*
          (exp (- (apply nnlp x*) mx))))))
    
  ;; for comparison: unit-normalized posterior
  (define unpostf20 (make-updf lnpf20 laf20))
  (define unp3df20
    (plot3d (surface3d unpostf20 140 160 4 20)
            #:title
            "Relative Posterior Normal Distributions of Heights (n=20)"
            #:x-label "Mean (μ) in cm"
            #:y-label "Standard Deviation (σ) in cm"
            #:z-label "Relative Plausibility"))

  #;(list unp3df20 qp3df20) ;; to see them juxtaposed

  ;;
  ;; Fit the adult heights
  ;;
  
  ;; non-normalized log posterior (compositionally constructed)
  (define lnp (make-joint adult-heights 'log))
  (define la (laplace-approx lnp '((140 160) (4 20)))) 

  (define qpdf (lapprox-pdf la))
  (define qp3d
    (plot3d (surface3d qpdf 153 157 7 11)
            #:title
            "Approximate Posterior Normal Distributions of Heights (n=352)"
            #:x-label "Mean (μ) in cm"
            #:y-label "Standard Deviation (σ) in cm"
            #:z-label "Probability")) ;; much tighter!

  ;; for comparison: non-normalized posterior
  (define unpost (make-updf lnp la))
  (define unp3d
    (plot3d (surface3d unpost 153 157 7 11)
            #:title
            "Relative Posterior Normal Distributions of Heights (n=352)"
            #:x-label "Mean (μ) in cm"
            #:y-label "Standard Deviation (σ) in cm"
            #:z-label "Relative Plausibility"))

  #;(list unp3d qp3d) ;; to see them juxtaposed
  (void))

#;
(begin
  (require "howell.rkt")           ; Howell Data
  (require "laplace-approx.rkt")   ; Laplace Approximation Engine
  (require "modelizer.rkt")        ; Interpret model specifications 
  (plot-height 300) (plot-width 300)
  ;;
  ;; Fit the first 20 heights
  ;;
  
  ;; grab the first 20 adult heights, just to get started, and help w/ testing
  (define first-20-heights (in-vector adult-heights 0 20))
  ;; for variety, select a random 20 adult heights
  #;(define random-20-heights (take (shuffle (vector->list adult-heights)) 20))
  ;; non-normalized log posterior (compositionally constructed)
  (define lnpf20 (make-joint first-20-heights 'log))

  (define mlnpf20 
    (make-log-compatibility-function
     `(model
       [type-decls [i Row]]
       [var-decls  [(h i) Number] [μ Number] [σ Number]]
       [var-defs   [(h i) . ~ . (normal μ σ)]
                   [μ     . ~ . (normal 178 20)]
                   [σ     . ~ . (uniform 0 50)]])
     (extend-env empty-env 'h (for/vector ([f first-20-heights]) f))
     ))
    
  (define laf20 (time (laplace-approx lnpf20 '((140 160) (4 20)))))
  (define mlaf20 (time (laplace-approx (λ (μ σ) (mlnpf20 `(,μ ,σ)))
                                 '((140 160) (4 20)))))

  (define qpdff20 (lapprox-pdf laf20))
  (define mqpdff20 (lapprox-pdf mlaf20))

  (define qp3df20
    (plot3d (surface3d qpdff20 140 160 4 20)
            #:title
            "Approximate Posterior Normal Distributions of Heights (n=20)"
            #:x-label "Mean (μ) in cm"
            #:y-label "Standard Deviation (σ) in cm"
            #:z-label "Probability"))

  (define mqp3df20
    (plot3d (surface3d mqpdff20 140 160 4 20)
            #:title
            "Approximate Posterior Normal Distributions of Heights (n=20)(*)"
            #:x-label "Mean (μ) in cm"
            #:y-label "Standard Deviation (σ) in cm"
            #:z-label "Probability"))

  ;; given non-normalized log posterior and mode,
  ;; produce unit-normalized posterior
  (define (make-updf nnlp la)
    (let ([mx (apply nnlp (lapprox-means la))])
      (with-arity-of nnlp
        (λ x*
          (exp (- (apply nnlp x*) mx))))))
    
  ;; for comparison: unit-normalized posterior
  (define unpostf20 (make-updf lnpf20 laf20))
  (define unp3df20
    (plot3d (surface3d unpostf20 140 160 4 20)
            #:title
            "Relative Posterior Normal Distributions of Heights (n=20)"
            #:x-label "Mean (μ) in cm"
            #:y-label "Standard Deviation (σ) in cm"
            #:z-label "Relative Plausibility"))

  (list unp3df20 qp3df20 mqp3df20) ;; to see them juxtaposed

  )
