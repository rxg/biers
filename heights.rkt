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
(plot (density (for/list ([i (range 1000)]) (model)) 2))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Real data (from Nancy Howell's data set)

;;
;; Experimenting with csv-reading and data-frame  Racket packages
;; (unfortunately the data-frame package's csv facilities are pretty weak)
;;
;; to load from the REPL or another module:
;; (require (submod "heights.rkt" howell))
;; or to load within this file::
;; (require (submod "." howell))

(module howell racket
  
  (require csv-reading) ;; for reading McElreath's csv file
  ;; McElreath uses data frames a lot in R: are Racket's data frames helpful?
  ;; Possibly not, since R may either not have a tradition of using
  ;; rich data abstractions, or possibly just not have them at all?
  (require data-frame)  
  (provide (all-defined-out))
  
  ;;
  ;; Helper Functions
  ;;
  
  ;; transform a vector of row vectors into a vector of column vectors
  (define (rows->cols data)
    (let ([count (vector-length (vector-ref data 0))])
      (for/list ([idx (range count)])
        (for/vector ([row data]) (vector-ref row idx)))))

  ;; create a data frame packed with the provided series data
  (define (make-populated-data-frame series-names data-rows)
    (let ([df (make-data-frame)]
          [data-vecs (rows->cols data-rows)])
      (begin ;; imperative update, then return
        (for ([name series-names]
              [data data-vecs])
          (df-add-series df (make-series name #:data data)))
        df)))

  ;; low-level interface:
  ;; create a new data-frame by filtering an old one 
  (define (df-filter df filter-fn series-names)
    (let ([data-rows (apply df-select* df series-names  #:filter filter-fn)])
      (make-populated-data-frame series-names data-rows)))
  
  ;; lambda-like syntax for writing data frame filters
  ;; parameterized on series names.  E.g.: (λf (age) (>= age 18))
  (define-syntax λf
    (syntax-rules ()
      [(_ (nm* ...) exp)
       (λ (series-names)
         (λ (row)
           (let ([nm*
                  (vector-ref row
                              (index-of series-names (symbol->string 'nm*)))]
                 ...)
             exp)))]))

  ;; high-level interface:
  ;; create a new data-frame by filtering an old one 
  (define (filter-df filter df)
    (let* ([series-names (df-series-names df)]
           [fn (filter series-names)])
      (df-filter df fn series-names)))
  
  ;;
  ;; Data frames for Howell data
  ;;
  
  (define csv-reader
    (make-csv-reader (open-input-file "Howell1.csv")
                     '((separator-chars #\;))))
  
  (define csv-data (csv->list csv-reader))
  
  ;; get column headers interpreted as strings
  (define headers (first csv-data))
  
  ;; get row data interpreted as exact numbers (rest skips the header)
  ;; yields a vector of row vectors (compatible with df-filter helpers below)
  (define data-rows
    (for/vector ([row (rest csv-data)])
      (for/vector ([datum row])
        (string->number datum 10 'read 'decimal-as-exact))))

  (define df (make-populated-data-frame headers data-rows))
  
  (define heights (df-select df "height"))

  ;; Let's get the adults:
  (define adults-df (filter-df (λf (age) (>= age 18)) df))

  (define adult-heights (df-select adults-df "height"))
  
  #;
  (begin
    (define coords
      (for/list ([w (in-vector (df-select adults-df "weight"))]
                 [h (in-vector (df-select adults-df "height"))])
        (vector w h)))
    ;; Scatter-plot of adult weight vs. height (see the correlation)
    (define psp (plot (points coords)))
    ;; density plot of adult height (see that it looks normal)
    (define pd (plot (density (df-select adults-df "height") 2)))
    (list psp pd))  
  (void)) ; module howell



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

;; prior is 2-dimensional: P(μ,σ|Γ) (using Γ to expose open-ended assumptions)
;; Assume μ,σ independent: P(μ,σ|Γ) = P(μ|Γ)P(σ|Γ)
;; can now produce log distribution too!
(define (prior-pdf μ σ [log? #f])
  (define cmb (if log? + *))
  (let ([σ-pdf (distribution-pdf (uniform-dist 0 50))]
        [μ-pdf (distribution-pdf (normal-dist 178 20))])
    (cmb (μ-pdf μ log?) (σ-pdf σ log?))))

;; Let's graph the prior pdf in 3D!
;; Going out 2 std-dev on μ
;; Going past 50 on σ 
#;
(plot3d (surface3d prior-pdf 138 218 0 60))


;; single-joint
;; P(h_i,μ,σ|Γ) = P(h_i|μ,σ,Γ)P(μ,σ|Γ)
;; conditionally independent likelihood assumption: = P(h_i|μ,σ)P(μ,σ|Γ)
#;(define (joint-i cprob prior h_i μ σ)
    (* (cprob h_i μ σ) (prior μ σ)))

;; multi-joint
;; P(h*,μ,σ|Γ) = Π_i P(h_i,μ,σ|Γ)
#;(define (joint cprob prior h* μ σ)
    (for/fold ([p 1.0])
              ([h_i h*])
      (* (joint-i cprob prior h_i μ σ) p)))


;; when combining probabilities, use multiplication;
;; when combining log-probabilities, use multiplication
(define (merge-pdf [log? #f])
  (if log? + *))

;; A bit more streamlined
;; and can now produce log distribution
#;(define (joint cprob prior h* μ σ [log? #f])
  (define unit (if log? 0.0 1.0))
  (define cmb (if log? + *))
  (for/fold ([p unit])
            ([h_i h*])
    (cmb (cprob h_i μ σ log?) (prior μ σ log?) p)))

;; TAKE 3:  don't factor in the prior for EACH DATA POINT!  They are
;; conditionally independent wrt μ and σ (by assumption)
(define (joint-likelihood cprob h* μ σ [log? #f])
  (define cmb (merge-pdf log?))
  (for/fold ([p-acc (sequence-ref h* 0)])
            ([h_i (sequence-tail h* 1)])
    (cmb (cprob h_i μ σ log?) p-acc)))


(define (correct-joint cprob prior h* μ σ [log? #f])
  (define cmb (merge-pdf log?))
  (cmb (joint-likelihood cprob h* μ σ log?) (prior μ σ log?)))

(define (make-correct-nnpost h* [log? #f])
  (define (cprob h_i μ σ log?)
    (pdf (normal-dist μ σ) h_i log?))
  (λ (μ σ)
    (correct-joint cprob prior-pdf h* μ σ log?)))

;; RG: I THINK THAT THE FOLLOWING IS WRONG!

;; Likelihood function for all the heights
;; heights are independent and identically distributed
;; P(h* | μ,σ) = Π_i P(h_i|μ,σ)
(define (likelihood h* μ σ)
  (let ([h-pdf (distribution-pdf (normal-dist μ σ))])
    (apply * (for/list ([h h*]) (h-pdf h)))))

;; Joint Distribution
;; P(h*,μ,σ) =  P(h*|μ,σ)P(μ,σ) = (Π_i P(h_i|μ,σ))P(μ)P(σ)

;; non-normalized posterior
#;(define (make-post h* [log? #f])
  (λ (μ σ)
    (joint
     (λ (h_i μ σ log?)
       (pdf (normal-dist μ σ) h_i log?))
     prior-pdf h* μ σ log?)))
  

;;
;; Terrifyingly inefficient numerical integration-based posterior
;;  I can't even get a single value back from posterior... :(
;;  At least it heats up the room a little ;)
;;


;; double-integration?  No, not so much.
(define (avg-likelihood h*)
  (adaptive-integrate
   (λ (μ)
     (adaptive-integrate 
      (λ (σ)
        (likelihood h* μ σ))
      0 50 .1)) ;; range of the uniform distribution
   120 240 .1)) ;; roughly 3σ range around the mean

(define (posterior μ σ h*)
  (/ (likelihood h* μ σ)
     (avg-likelihood h*)))




  
;; From here I can do grid sampling...

;;
;; Grid Sampling
;;


;; Build a 1D grid of values
(define (seq lo hi size)
  (range lo hi (/ (- hi lo) size)))


(define (fit-heights h* μ-lo μ-hi σ-lo σ-hi)
  (local
    ;; Build a 200 x 200 grid of σ,μ values
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
       (for*/list ([μ (in-list μ*-grid)][σ (in-list σ*-grid)])
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

(define (plot-μσ-posterior prob*)
  (plot3d
   (points3d prob* #:size 1)
   #:title "Posterior Normal Distribution of Heights in cm (Howell Data)"
   #:x-label "Mean (μ)"
   #:y-label "Standard Deviation (σ)"
   #:z-label "Relative Plausibility"))

;; Example Fit
#;
(begin
  (require (submod "." howell))
  ;; grab the first 20 adult heights, just to get started
  (define first-20-heights (in-vector adult-heights 0 20))
  (define random-20-heights (take (shuffle (vector->list adult-heights)) 20))
  (define post-20 (fit-heights first-20-heights 150 170 4 20))
  (define post (fit-heights adult-heights 153 156 6 9))
  (map plot-μσ-posterior (list post-20 post))
  ;; Sample from the posterior
  (define-values (μσ* w)
    (for/lists (uss w) ([p post]) (values (vector-take p 2) (vector-ref p 2))))
  (define post-dist (discrete-dist μσ* w))
  (define s* (sample post-dist 10000))
  (void))

;; Example of drawing samples from the posterior
#;
(begin
  (require (submod "." howell))
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
  (require (submod "." howell))
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
;; WARNING:  nlopt still a bit broken, possible memory safety error!
;;


(begin
  (require (submod "." #;"heights.rkt" howell))
  (require "laplace-approx.rkt")

  ;;
  ;; Fit the first 20 heights
  ;;
  
  ;; grab the first 20 adult heights, just to get started, and help w/ testing
  (define first-20-heights (in-vector adult-heights 0 20))
  ;; for variety, select a random 20 adult heights
  ;(define random-20-heights (take (shuffle (vector->list adult-heights)) 20))


  ;(define post-20 (fit-heights first-20-heights 150 170 4 20))
  ;; non-normalized log posterior (compositionally constructed)
  (define lnpf20 (make-correct-nnpost first-20-heights 'log))
  (define laf20 (laplace-approx lnpf20 '((150 170) (4 20))))
  (define modef20 (lapprox-means laf20))
  (define qpdff20 (lapprox-pdf laf20))
  (define qp3df20 (plot3d (surface3d qpdff20 150 170 4 20)))

  ;; given non-normalized log posterior and mode,
  ;; produce unit-normalized posterior
  (define (make-updf nnlp mode)
    (let ([mx (apply nnlp mode)])
      (with-arity-of nnlp
        (λ x*
        (exp (- (apply nnlp x*) mx))))))
    
  ;; for comparison: unit-normalized posterior
  (define unpostf20 (make-updf lnpf20 modef20))
  (define unp3df20 (plot3d (surface3d unpostf20 150 170 4 20)))

  #;(list qp3df20 unp3df20) ;; to see them juxtaposed

  ;;
  ;; Fit the adult heights
  ;;
  
  ;; non-normalized log posterior (compositionally constructed)
  (define lnp (make-correct-nnpost adult-heights 'log))
  (define la (laplace-approx lnp '((150 170) (4 20)))) 
  (define mode (lapprox-means la))

  (define qpdf (lapprox-pdf la))
  (define qp3d (plot3d (surface3d qpdf 156 160 7 11))) ;; much tighter!

  ;; for comparison: non-normalized posterior
  (define unpost (make-updf lnp mode))
  (define unp3d (plot3d (surface3d unpost 156 160 7 11)))

  #;(list qp3d unp3d) ;; to see them juxtaposed
  (void))
