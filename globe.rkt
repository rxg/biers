#lang racket

(provide (all-defined-out))
(require plot)
(require (only-in plot/utils linear-seq))
(require math/number-theory)
(require math/distributions)
(require math/statistics)
(require "utils.rkt")
;; Integration, from https://github.com/mkierzenka/Racket_NumericalMethods
(require "Racket_NumericalMethods/Numerical_Integration.rkt") 

;;
;; globe.rkt - Globe inference model from Chapter 2 of
;;   McElreath, Statistical Rethinking
;;   Plus sampling-based techniques from Chapter 3

;; What proportion of the Earth's surface is water?

;; McElreath offers an experiment where you toss a model of the globe in
;; the air, and track where your index finger lands when you catch it:
;; on water, or on land?  You can apply this data (outcomes) to a statistical
;; inference engine to acquire statistical information about plausible
;; proportions of land and water.


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Stochastic Simulator - 'cause I'm too lazy to toss an inflatable globe
;;   simulator acts as "ground truth".  In practice, this simulator is TOO
;;   closely tied to the statistical model and not real life:  in practice
;;   there are likely to be correlations between globe tosses because land is
;;   not evenly distributed on the earth and people suck at tossing globes.
;;   In short: models "are not even wrong."

;; UNESCO says the proportion of water covering the surface of Earth is
;; between 71% and 72% (it changes with tides, etc.)
(define WATER-PROPORTION 71)
(define LAND-PROPORTION (- 100 WATER-PROPORTION))

;; Obs is 'water or 'land

;; (Samples n) is { l ∈ (listof Obs) | (length l) = n }
;; - a list of n samples

;; Summary is (list Natural Natural)
;; (list w n) is the number of water outcomes w
;;            out of a total number of trials n

;; n:Natural -> (Samples n)
;; draw n samples from the distribution of land vs. water
(define (draw-globe-samples n)
  (let ([dist (discrete-dist (list 'land 'water)
                             (list LAND-PROPORTION
                                   WATER-PROPORTION))])
    (sample dist n)))

;; Samples(n) -> Summary
;; Given samples, report the summary
(define (summarize-samples samples)
  (list (length (filter (λ (s) (equal? s 'water)) samples))
        (length samples)))

;; Samples(n) -> (listof Summary)
;; produce a list of running summaries from a list of samples including 0 case
(define (running-summary samples)
  (for/list ([n (in-range (add1 (length samples)))])
    (summarize-samples (take samples n))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Inference Engine
;;


;; Stochastic Variables
;; n - Number of globe tosses overall
;; p - Estimate of the proportion of water covering the globe
;; w - Number of water observations from tossing the globe


;; likelihood function - binomial distribution 
;; Pr(w | n,p)
;; the count of water observations w is distributed binomially, with
;; probability p of water on each toss, and n tosses total.
(define (exact-likelihood w n p)
  (*
   (/ (factorial n)
      (* (factorial w) (factorial (- n w))))
   (expt p w)
   (expt (- 1 p) (- n w))))

;; ...or build one using the binomial distribution that comes with Racket
(define (likelihood w n p)
  (let ([d (binomial-dist n p)])
    (pdf d w)))

;;
;; Candidate Prior distributions on p:
;;

;; uniform distribution (it's 1 everywhere!)
(define (flat-prior p)
  (let ([d (uniform-dist 0 1)])
    ;; uniform-dist produces the inexact number 1.0: fix it
    (inexact->exact (pdf d p))))


;; DEFINITELY more than half of the globe is water
(define (step-prior p)
  (if (< p 1/2)
      0
      2))


;; funky spikey prior
(define (peaked-prior p)
  (exp (* -5 (abs (- p 0.5)))))


;;
;; Joint Probability: Likelihood * Prior
;;
(define (joint w n p prior)
  (* (likelihood w n p)
     (prior p)))

;;
;; Bayesian Inference Using Numerical Integration
;;

;; Calculate the "average likelihood" a.k.a. marginal via numerical integration
(define (marginal w n prior)
  (adaptive-integrate (λ (p) (joint w n p prior)) 0 1 0.001))

;; Posterior pdf generator (convenient for plotting and queries)
(define (mk-posterior w n prior)
  (let ([avg (marginal w n prior)]) ;; compute the marginal once
    (λ (p)
      (/ (joint w n p prior)
         avg))))

;; plot the posterior for the given summary data
;; use #f for y-max to automate the choice
(define (plot-posterior w n prior y-max)
  (plot (function (mk-posterior w n prior) 0 1)
        #:x-label #f #:y-label #f #:y-min 0 #:y-max y-max))


;;
;; Bayesian Inference Using Grid Approximation
;;

;; McElreath's approximate marginal is faux-discrete: kinda wonky.
;; Instead we approximate the integral with boxes a la Gelman et al.

;; a dead simple integral approximation. 
;; doesn't try to be smart at the boundaries
(define (box-integrate ys step)
  #;(sum (for/list ([y ys]) (* y step)))
  ;; exploit distributivity of multiplication over sums
  (* (sum ys) step))

;; compute joint probability at a set of points
(define (joints w n prior points)
  (for/list ([p points]) (joint w n p prior)))

;; approximate the posteriors at the given points
(define (grid-points-posterior w n prior points)
  (let* ([joint* (joints w n prior points)]
         [step (grid-step-size points)] 
         [approx-marginal (box-integrate joint* step)])
    (for/list ([joint joint*]) (/ joint approx-marginal))))

;; Grid(m,n,count) is (listOf Rational)
;; grid of count evenly distributed points from m to n.

;; create a grid of count points from m to n
(define (make-grid m n count) (linear-seq m n count))

;; recover the grid step-size (we know it's uniform)
(define (grid-step-size points)
  (- (second points) (first points)))

;; Produce a grid of "count" posterior points
(define (grid-count-posterior w n prior count)
  (let ([points (make-grid 0 1 count)])
    (grid-points-posterior w n prior points)))

;; collate points and posteriors into plotting coordinates
(define (line-coords points posteriors)
  (map vector points posteriors))

;; given w n and count, generate a list of posterior coords
(define (grid-approx-coords w n prior count)
  (let* ([points (make-grid 0 1 count)]
         [posteriors (grid-points-posterior w n prior points)])
    (line-coords points posteriors)))

;; plot coordinates as a line and maybe label the points
(define (plot-coords coords labels?)
  (let ([maybe-points (if (not labels?)
                          empty
                          (map point-label coords))])
    (plot (cons (lines coords) maybe-points)
          #:y-min 0 #:x-label #f #:y-label #f)))

;; plot "count" grid posterior points for the summary (w,n)
(define (plot-grid-approx w n prior count labels?)
  (let* ([coords (grid-approx-coords w n prior count)])
    (plot-coords coords labels?)))

;;
;; More Visualization Tools
;;

;; interleave the contents of two lists
;; to show each sample and its effect on the posterior
(define (interleave lsta lstb)
  (cond
    [(empty? lsta) lstb]
    [else (cons (first lsta)
                (interleave lstb (rest lsta)))]))

;; given a list of samples, trace sequential Bayesian inference.
;; parameterized on inference/plotting engine
(define (simulate-bayes obs* infer-then-plot)
  (let ([sm* (running-summary obs*)])
    (let ([plots
           (for/list ([sm sm*])
             (match-let ([(list w n) sm])
               (infer-then-plot w n)))])
      (interleave plots obs*))))

;; wrapper to simulate inference using grid approximation
(define (simulate-grid-bayes prior size obs*)
  (simulate-bayes obs*
                  (λ (w n) (plot-grid-approx w n prior size #f))))

;; wrapper to simulate inference using numerical integration
(define (simulate-integral-bayes prior obs* y-max)
  (simulate-bayes obs*
                  (λ (w n) (plot-posterior w n prior y-max))))

;; Better visualization using both!
;; a) use grid approximation to find the max posterior
;; b) use integration to plot simulation, now with common plot axes
(define (plot-sequential-inference prior obs*)
  (let ([sm* (running-summary obs*)]
        [grid-size 20])
    ;; pre-compute approximations to get y-axis right
    (let ([grids (for/list ([sm sm*])
                   (match-let ([(list w n) sm])
                     (grid-count-posterior w n prior grid-size)))])
      (let ([y-max (apply max
                          (apply append grids))]
            [y-scale-factor 1.1]) ;; add some breathing room at top of plot
        (simulate-integral-bayes prior obs* (* y-max y-scale-factor))))))

;; Globally set plot dimensions
(plot-width 120)
(plot-height 120)

;; Example use
#;
(begin  
  (define s* (draw-globe-samples 10))
  s*
  (plot-sequential-inference flat-prior s*))


;; Example from McElreath Chapter 2:
;; comparing Quadratic Approximation to the ideal (Figure 2.8, p. 43)
#;
(plot
 (list
  (function (distribution-pdf (beta-dist 7 4)) 0 1 #:color "blue")
  (function (distribution-pdf (normal-dist 0.67 0.16)) 0 1 #:color "black")))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Chapter 3: Sampling-based analysis/summarization of the posterior

;;
;; Tools for working with grid posteriors
;;
;; grid Distribution object
(define-struct gd (grid density dist))

;; produce a count-sized grid distribution for the posterior of the parameters
(define (gd-posterior w n prior count)
  (define p-grid (make-grid 0 1 count))
  (define posterior (grid-points-posterior w n prior p-grid))
  (define dist (discrete-dist p-grid posterior))
  (make-gd p-grid posterior dist))

;; Render the distribution of gd for 2D plotting
(define (render-gd gd)
  (lines (map vector (gd-grid gd) (gd-density gd))))

;; Draw n samples from the given grid posterior
(define (sample-gd gd n)
  (sample (gd-dist gd) n))

(define (render-samples samples)
  (points (map vector (range (length samples)) samples)))
;; The following example uses grid approximation to produce a posterior,
;; samples from it, and then plots the samples several ways, comparing
;; the resulting density to the posterior itself.
#;
(begin
  ;; grid distribution
  (plot-width 360) (plot-height 360)
  (define p1e3 (gd-posterior 6 9 flat-prior 1000))
  ;; Draw samples from the grid approximation (Chapter 3)
  (define samples (sample-gd p1e3 10000))
  ;; points ordered by draw (comment due to Suzanna)
  (plot (render-samples samples)
        #:y-min 0 #:y-max 1)
  ;; From Suzanna, put the points in order by value
  (plot (render-samples (sort samples <))
        #:y-min 0 #:y-max 1)
  ;; density plot of the samples
  (plot (density samples 2)
        #:x-min 0 #:x-max 1)
  ;; for comparison, the source of the samples
  (plot (render-gd p1e3))
  (void))

;;
;; Summary Statistics, computed either from the posterior or from samples
;;

;; Intervals of defined boundaries, based on the grid posterior
;; compute the mass of probability that satisfies pred
(define (pr-posterior pred p-grid posterior)
  (let* ([coords (map list p-grid posterior)]
         [relevant-coords (filter (λ (p) (pred (first p))) coords)]
         [relevant-mass (map second relevant-coords)])
    (/ (sum relevant-mass)
       (sum posterior))))
    
#;
(begin
  (pr-posterior (λ (p) (< p 0.5)) p-grid posterior)
  (void))

;; Intervals of defined boundaries, based on the samples
(define (pr-samples pred samples)
  (/ (count pred samples)
     (length samples)))

#;
(begin
  (pr-samples (λ (p) (< p 0.5)) samples)
  (pr-samples (λ (p) (< 0.5 p 0.75)) samples)
  (void))


;; lower and upper boundaries of the *middle* q percent of probability mass
(define (percentile-interval samples q)
  (let* ([split (/ q 2)]
         [lo (- 1/2 split)]
         [hi (+ 1/2 split)])
    (for/list ([p (list lo hi)])
      (list p (quantile p < samples)))))

;; Highest Posterior Density (HPD) Interval
;; calculate the minimum interval that contains q percent of probability mass
(define (hpdi samples q)
  ;; WHY WHY use multiple return values?!?
  (define-values (lo hi) (real-hpd-interval q samples))
  (list lo hi))


;; calculate the loss (wrt loss function loss-fn) for point p-guess against
;; the grid posterior (p-grid,posterior)
(define (calculate-loss p-guess p-grid posterior loss-fn)
  (/ (sum (for/list ([x p-grid]
                     [y posterior])
            (* (loss-fn p-guess x) y)))
     (sum posterior)))

;; approximate a point estimate with respect to the given loss-function
;; for the approximate posterior (pgrid,posterior)
(define (point-estimate loss-fn p-grid posterior)
  (argmin (λ (pg) (calculate-loss pg p-grid posterior loss-fn))
          p-grid))

;; Example loss functions

;; absolute loss: corresponds to the mean (expected value) of the posterior
(define (absolute-loss pg p) (abs (- pg p)))

;; quadratic loss: corresponds to the median of the posterior
(define (quadratic-loss pg p) (sqr (- pg p)))
         
;; 0-1 loss: corresponds to *a* mode of the posterior (not necessarily unique)
(define (0-1-loss pg p)  (if (= pg p) 0 1))

;;
;; Laplace Approximation:
;;   Compute the posterior mode, then approximate the posterior using a
;;   Gaussian distribution (by approximating the log-posterior using a
;;   quadratic function).
;;

#;
(begin
  (require "laplace-approx.rkt")
  (define post (mk-posterior 6 9 flat-prior))
  (define lpost (compose log post))
  (define mode (maximize-fn lpost '((0 1))))
  (define la (laplace-approx lpost (list mode)))
  (define la-post (lapprox-pdf la))
  ;; We can draw samples from the Laplace approximation!
  #;(define samples (lapprox-sample la 10))
  (plot-width 360)
  (plot-height 360)
  (define cmp
    (plot
     (list (function post 0 1
                     #:label "Numerical Integration Posterior" #:color 'black)
           (function la-post 0 1
                     #:label "Laplace Approximation Posterior" #:color 'cyan))))
  cmp)


;;
;; Sampling to simulate prediction
;;


;; Let's make a likelihood sampler.
;; produces a list of w's for given n and p
(define (likelihood-sampler n p count)
  (let ([d (binomial-dist n p)])
    (sample d count)))


;; draw likelihood samples proportional to samples from a given distribution
(define (predictive-dist-samples dist-samples n)
  (for/fold ([samples empty])
            ([p dist-samples])
    (append (likelihood-sampler n p 10000) samples)))

;;
;; The following are previous alternative implementations with space/time issues
;;

;; THIS ONE IS A SPACE HOG, especially as the numbers get big! 
(define (predictive-dist-samples2 dist-samples n)
  (apply append
         (for/list ([p dist-samples])
           (likelihood-sampler n p 10000))))

;; THIS ONE IS SLOW: maybe loops can be tuned, but  
(define (predictive-dist-samples3 dist-samples n)
  (let* ([size (* (length dist-samples) 10000)]
         [vec (make-vector size)])
    (for ([i (in-range 0 size 10000)]
          [p (in-list dist-samples)])
      (for ([s (in-list (likelihood-sampler n p 10000))]
            [j (in-range i (+ i 10000))])
        (vector-set! vec j s)))
    vec))



;; Timing Code
#;
(time (begin (predictive-dist-samples0 post*) (void)))


#;
(begin
  (likelihood-sampler 9 (/ WATER-PROPORTION 100) 5)
  (count-samples (likelihood-sampler 2 0.7 100000))
  (plot (render-hist (likelihood-sampler 9 0.7 100000)))
  (void))

;; Build a grid prior for sampling
(define (gd-prior prior count)
  (define p-grid (make-grid 0 1 count))
  (define grid-prior (map prior p-grid))
  (define dist (discrete-dist p-grid grid-prior))
  (make-gd p-grid grid-prior dist))


;; Example of constructing a posterior predictive distribution 
;; and prior predictive distributions via sampling
#;
(begin
  (plot-width 360) (plot-height 360)
  ;; grid posterior
  (define post (gd-posterior 6 9 flat-prior 10000))

  ;; grid priors 
  (define flat-gd (gd-prior flat-prior 1000))  
  (define step-gd (gd-prior step-prior 1000))   
  (define peaked-gd (gd-prior peaked-prior 1000))

  ;; Draw samples from the grid posterior approximation (Chapter 3)
  (define post* (sample-gd post 10000))
  ;; predictive posterior samples
  (define predict-post* (predictive-dist-samples post* 9))
  ;; plot the posterior predictive distribution
  (plot (render-hist predict-post*))
  (void))