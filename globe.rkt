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
(define (plot-posterior w n prior [y-max #f])
  (plot (function (mk-posterior w n prior) 0 1)
        #:x-label #f #:y-label #f #:y-min 0 #:y-max y-max
        #:width 120 #:height 120))


;;
;; Bayesian Inference Using Grid Approximation
;;

;; Grid(m,n,count) is (listOf Rational)
;; grid of count evenly distributed points from m to n.

;; create a grid of count points from lo to hi
(define (make-grid lo hi count) (linear-seq lo hi count))

;; lo value of grid
(define (grid-lo g) (first g))

;; hi value of grid
(define (grid-hi g) (last g))

; number of grid points
(define (grid-count g) (length g))

;; grid step-size (well-defined because our grids are uniform)
(define (grid-step-size g)
  (- (second g) (first g)))

;; McElreath's approximate marginal is faux-discrete: kinda wonky.
;; Instead we approximate the integral with boxes a la Gelman et al. BDA3

;; a dead simple integral approximation. 
;; doesn't try to be smart at the boundaries
(define (box-integrate ys step)
  #;(sum (for/list ([y ys]) (* y step)))
  ;; exploit distributivity of multiplication over sums
  (* (sum ys) step))

;; compute joint probability at a set of points
(define (joints w n prior p*)
  (for/list ([p p*]) (joint w n p prior)))

;; approximate the posteriors at the given points
(define (grid-points-posterior w n prior grid)
  (let* ([joint* (joints w n prior grid)]
         [step (grid-step-size grid)] 
         [approx-marginal (box-integrate joint* step)])
    (for/list ([joint joint*]) (/ joint approx-marginal))))

;;
;; Tools for working with grid posteriors
;;
;; grid distribution object
(define-struct gd (grid density dist))

;; forwarding grid selectors
(define (gd-lo gd) (grid-lo (gd-grid gd)))
(define (gd-hi gd) (grid-hi (gd-grid gd)))
(define (gd-count gd) (grid-count (gd-grid gd)))
(define (gd-step-size gd) (grid-step-size (gd-grid gd)))


;; produce a count-sized grid distribution for the posterior of the parameters
(define (gd-posterior w n prior count)
  (define p-grid (make-grid 0 1 count))
  (define posterior (grid-points-posterior w n prior p-grid))
  (define dist (discrete-dist p-grid posterior))
  (make-gd p-grid posterior dist))

(define (gd-coords gd)
  (map vector (gd-grid gd) (gd-density gd)))

(define (gd-zero-coords gd)
  (define n (grid-count (gd-grid gd)))
  (map vector  (gd-grid gd)
       (for/list ([i n]) 0)))
  
;; Render the distribution of gd for 2D plotting and maybe label the points
(define (render-gd gd [labels? #f])
  (define coords (gd-coords gd))
  (define maybe-points (if (not labels?)                          
                           empty
                           (map point-label coords)))
  (cons (lines (gd-coords gd)) maybe-points))

(define (plot-gd gd [labels? #f])
  (plot (render-gd gd labels?)))


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
                  (λ (w n) (plot-gd (gd-posterior w n prior size) #f))))

;; wrapper to simulate inference using numerical integration
(define (simulate-integral-bayes prior obs* y-max)
  (simulate-bayes obs*
                  (λ (w n) (plot-posterior w n prior y-max))))

;; Better visualization using both!
;; a) use grid approximation to find the max posterior
;; b) use integration to sequentially plot each inference step,
;; using the same axis scales throughout 
(define (plot-sequential-inference prior obs*)
  (let ([sm* (running-summary obs*)]
        [grid-size 20])
    ;; pre-compute approximations to get y-axis right
    (let ([grids (for/list ([sm sm*])
                   (match-let ([(list w n) sm])
                     (gd-density (gd-posterior w n prior grid-size))))])
      (let ([y-max (apply max
                          (apply append grids))]
            [y-scale-factor 1.1]) ;; add some breathing room at top of plot
        (simulate-integral-bayes prior obs* (* y-max y-scale-factor))))))


;; Example use
#;
(begin  
  (define s* (draw-globe-samples 10))
  s*
  (plot-sequential-inference flat-prior s*))


;; Example from McElreath Chapter 2:
;; Comparing Quadratic Approximation to the ideal (Figure 2.8, p. 43)
;; Bonus: including integral approximation
#;
(begin
  (require "laplace-approx.rkt")
  (define la-plots
    (for/list ([w '(6 12 18)]
               [n '(9 18 27)])
      (define (log-joint p) (log (joint w n p flat-prior)))
      (define lapprox-d (lapprox-pdf
                         (laplace-approx log-joint '((0 1)))))
      (define analytic-d (distribution-pdf (beta-dist (add1 w) (add1 (- n w)))))
      (define infer-d (mk-posterior w n flat-prior))
      (plot (list (function lapprox-d  0 1 #:color "blue")
                  (function analytic-d 0 1 #:color "pink")
                  (function  infer-d 0 1 #:color "black"
                             #:style 'dot #:width 2))
            #:title (format "n = ~a" n)
            #:y-label "Density" #:x-label "proportion water")))
  (void))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Chapter 3: Sampling-based analysis/summarization of the posterior


;; draw n samples from the given grid posterior
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
  (define p1e3 (gd-posterior 6 9 flat-prior 1000))
  ;; draw samples from the grid approximation (Chapter 3)
  (define samples (sample-gd p1e3 10000))
  ;; points ordered by draw (comment due to Suzanna)
  (define sample-plot1 (plot (render-samples samples)
                             #:y-min 0 #:y-max 1))
  ;; from Suzanna, points ordered by value
  (define sample-plot2 (plot (render-samples (sort samples <))
                             #:y-min 0 #:y-max 1))
  ;; density plot of the samples
  (define sample-plot3 (plot (density samples 2)
                             #:x-min 0 #:x-max 1))
  ;; for comparison, the source of the samples
  (define sample-plot4 (plot-gd p1e3))
  #;(list sample-plot1 sample-plot2 sample-plot3 sample-plot4)
  (void))

;;
;; Summary Statistics, computed either from the posterior or from samples
;;

;; Intervals of defined boundaries, based on the grid posterior
;; compute the mass of probability that satisfies pred
(define (pr-posterior pred? p-grid posterior)
  (define relevant-mass
    (for/list ([p p-grid]
               [d posterior]
               #:when (pred? p))
      d))
  (/ (sum relevant-mass)
     (sum posterior)))

;; (p -> Boolean) -> Grid -> Prob
;; how much probability lies in the posterior region satisfying pred?
;; typically pred defines a bounded region
(define (gd-mass-in pred gd)
  (pr-posterior pred (gd-grid gd) (gd-density gd)))

#;
(begin
  (define post-69 (gd-posterior 6 9 flat-prior 1000))
  (define p<0.5 (gd-mass-in (λ (p) (< p 0.5)) post-69))
  (void))


;; Intervals of defined boundaries, based on the samples
(define (pr-samples pred samples)
  (/ (count pred samples)
     (length samples)))

;; unified: sampling-based variant of gd-mass-in
(define (gd-mass-in-samples pred gd)
  (define samples (sample-gd gd 10000))
  (pr-samples pred samples))

#;
(begin
  (define post-69 (gd-posterior 6 9 flat-prior 1000))
  (define p<0.5 (gd-mass-in-samples (λ (p) (< p 0.5)) post-69))
  (define p0.5-0.75 (gd-mass-in-samples (λ (p) (< 0.5 p 0.75)) post-69))
  (void))


;; lower and upper boundaries of the *middle* q percent of probability mass
(define (percentile-interval samples q)
  (let* ([split (/ q 2)]
         [lo (- 1/2 split)]
         [hi (+ 1/2 split)])
    (for/list ([p (list lo hi)])
      (quantile p < samples)
      #;(list p (quantile p < samples)))))

;; grid-posterior wrapper for percentile-interval
(define (gd-compatibility-interval gd q)
  (percentile-interval (sample-gd gd 10000) q))

;; Highest Posterior Density (HPD) Interval
;; calculate the minimum interval that contains q percent of probability mass
(define (hpdi samples q)
  ;; WHY WHY use multiple return values?!?
  (define-values (lo hi) (real-hpd-interval q samples))
  (list lo hi))

;; grid-posterior wrapper for hpdi
(define (gd-hpdi gd q)
  (hpdi (sample-gd gd 10000) q))

;;
;; exploit credibility intervals in plots
;;

;; fill the lo-to-hi region of gd-post posterior plot
(define (plot-with-interval gd-post lo hi)
  (define zeros (gd-zero-coords gd-post))
  (define post (gd-coords gd-post))
  (plot (list
         (lines-interval zeros post
                         #:x-min lo #:x-max hi
                         #:line1-style 'transparent #:line2-style 'transparent)
         (lines post))))

;; hpdi version 
(define (plot-with-hpdi gd-post q)
  (match-define `(,lo ,hi) (gd-hpdi gd-post q))
  (plot-with-interval gd-post lo hi))

;; compatibility interval version
(define (plot-with-compatibility-interval gd-post q)
  ;; use samples to determine credibility interval
  (match-define `(,lo ,hi) (gd-compatibility-interval gd-post q))
  (plot-with-interval gd-post lo hi))

#;
(begin
  (define gd (gd-posterior 6 9 flat-prior 1000))
  (define plot-ci (plot-with-compatibility-interval gd 0.5))
  (define plot-hpdi (plot-with-hpdi gd 0.5))
  (void))
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