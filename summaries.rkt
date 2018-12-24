#lang racket

(provide (all-defined-out))
(provide median mean stddev variance)

;;
;; summary-stats.rkt - some techniques for computing summary statistics
;;    
(require plot)
(require math/base)
(require math/number-theory)
(require math/distributions)
(require math/statistics)

;;
;; Summary Statistics, computed either from the posterior or from samples
;;

;; Intervals of defined boundaries, based on the grid posterior
;; compute the mass of probability that satisfies pred
(define (pr-mass-posterior pred p-grid posterior)
  (let* ([coords (map list p-grid posterior)]
         [relevant-coords (filter (λ (p) (pred (first p))) coords)]
         [relevant-mass (map second relevant-coords)])
    (/ (sum relevant-mass)
       (sum posterior))))
    
#;
(pr-mass-posterior (λ (p) (< p 0.5)) p-grid posterior)

;; Intervals of defined boundaries, based on the samples
;; compute the mass of probability that satisfies pred
(define (pr-mass-samples pred samples)
  (/ (count pred samples)
     (length samples)))

#;(pr-mass-samples (λ (p) (< p 0.5)) samples)
#;(pr-mass-samples (λ (p) (< 0.5 p 0.75)) samples)

;; DEPRECATED!: quantile is already defined in math/statistics.
;; produce the UPPER boundary of the lower q percent of probability mass
;; (the LOWER boundary is 0).
(define (my-quantile q < samples)
  (let ([pos (inexact->exact (ceiling (* q (length samples))))])
    (list-ref (sort samples <) pos)))

;; lower and upper boundaries of the *middle* q percent of probability mass
(define (percentile-interval samples q)
  (let* ([split (/ q 2)]
         [lo (- 1/2 split)]
         [hi (+ 1/2 split)])
    (map list
         (list lo hi)
         (map (λ (p) (my-quantile p < samples))
              (list lo hi)))))

;; Highest Posterior Density (HPD) Interval
;; calculate the minimum interval that contains q percent of probability mass
(define (hpdi samples q)
  ;; WHY WHY use multiple return values?!?
  (define-values (lo hi) (real-hpd-interval q samples))
 (list lo hi))


;; calculate the loss (wrt loss function loss-fn) for point p-guess against
;; the grid posterior (p-grid,posterior)
(define (calculate-loss p-guess p-grid posterior loss-fn)
  (/ (sum (map (λ (x y) (* (loss-fn p-guess x) y)) p-grid posterior))
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
;; Maximum a posterior (for a posterior function)
;;


;; Using nlopt (wrapper for NLopt library) to find map.

(module maximize racket
  (require nlopt)
  (provide maximize-posterior)
  (define (maximize-posterior post)
    (define-values
      (fst snd)
      (optimize/args post 1 
                   #:maximize #t
                   #:bounds '((0.0 . 1.0))
                   #:method 'GN_DIRECT_L))
    ;; no idea what the first value is about!
    snd)
) ; module maximize
(require (submod "." maximize))


;;
;; Sampling to simulate prediction
;;

;; code was very specialized...


;; Custom histogram plotting ('cause the plot density plotter
;; is sometimes weird)
(define (mk-hist-coords samples)
  (let-values ([(x* y*) (count-samples samples)])
    (let ([pre-coords (map vector x* y*)])
      (sort pre-coords < #:key (λ (v) (vector-ref v 0))))))

(define (plot-simple-hist samples)
  (plot (discrete-histogram (mk-hist-coords samples))))
