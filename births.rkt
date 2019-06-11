#lang racket

;; Chapter 3 Hard exercises, using the birth data set
;; The data for first birth and second birth

;; globe's binomial model is reusable here.
;; code "boy" as "water"
(require "globe.rkt") 

;; first birth: index indicates family
(define birth1
  (vector 'boy 'girl 'girl 'girl 'boy 'boy 'girl 'boy 'girl 'boy
          'girl 'girl 'boy 'boy 'girl 'boy 'boy 'girl 'girl 'girl 'boy
          'girl 'girl 'girl 'boy 'girl 'girl 'girl 'girl 'boy 'boy 'boy
          'girl 'boy 'girl 'boy 'boy 'boy 'girl 'boy 'girl 'boy 'boy 'girl
          'boy 'girl 'girl 'boy 'boy 'girl 'boy 'girl 'girl 'girl 'girl
          'girl 'girl 'girl 'boy 'boy 'girl 'boy 'girl 'girl 'boy 'girl
          'girl 'girl 'boy 'girl 'girl 'boy 'boy 'boy 'boy 'girl 'boy 'girl
          'boy 'boy 'boy 'boy 'boy 'girl 'girl 'boy 'girl 'boy 'boy 'girl
          'boy 'girl 'boy 'boy 'boy 'girl 'boy 'boy 'boy 'boy))

;; second birth: index indicates family
(define birth2
  (vector 'girl 'boy 'girl 'boy 'girl 'boy 'boy 'boy 'girl 'girl
          'boy 'boy 'boy 'boy 'boy 'girl 'girl 'boy 'boy 'boy 'girl 'girl
          'boy 'boy 'boy 'girl 'boy 'boy 'boy 'girl 'boy 'boy 'boy 'girl
          'boy 'girl 'girl 'boy 'boy 'boy 'boy 'girl 'girl 'boy 'girl 'boy
          'boy 'boy 'boy 'boy 'boy 'boy 'boy 'boy 'boy 'boy 'boy 'boy 'boy
          'boy 'boy 'girl 'boy 'boy 'girl 'boy 'boy 'girl 'boy 'boy 'boy
          'girl 'girl 'girl 'girl 'girl 'girl 'boy 'girl 'girl 'girl 'boy
          'boy 'girl 'girl 'boy 'girl 'girl 'boy 'boy 'girl 'girl 'girl
          'boy 'boy 'boy 'girl 'girl 'girl 'girl))

(define (birth? x) (or (boy? x) (girl? x)))
(define (boy? b) (eq? b 'boy))
(define (girl? b) (eq? b 'girl))

(define num-births (+ (vector-length birth1)
                      (vector-length birth2)))
  
(define num-boys (+ (vector-count boy? birth1)
                    (vector-count boy? birth2)))

(define num-girls (+ (vector-count girl? birth1)
                     (vector-count girl? birth2)))

(define post (gd-posterior num-boys num-births flat-prior 10000))

;; which parameter value maximizes the posterior probability?
(define post-max
  (let ([max-coord (argmax (λ (coord) (vector-ref coord 1))
                           (gd-coords post))])
    (vector-ref max-coord 0)))

;; Draw 10,000 samples from the posterior

(define post-samples (gd-sample post 10000))
(define h (hpdi post-samples 0.5))

;; this capability was already packaged
(define h2 (gd-hpdi post 0.5))

(define h2plot1 (plot-with-hpdi post 0.5))

(define h2plot2 (plot-with-hpdi post 0.89))

(define h2plot3 (plot-with-hpdi post 0.97))

;; predictive posterior distribution
(define pred (predictive-dist-from-p-samples post-samples 200))
(require plot)
(require "utils.rkt")
(define ppost-plot (plot (render-hist pred)))

;; Now consider samples of 100 boys (to compare against first-born batch)
(define num-first-boys (vector-count boy? birth1))
(define pred100 (predictive-dist-from-p-samples post-samples 100))
(define ppost-plot100 (plot (density pred100)))

;; what about kids that followed girls
(define after-girl
  (for/list ([b1 birth1]
             [b2 birth2]
             #:when (equal? b1 'girl))
    b2))

(define num-girl-first-borns (length after-girl))
(define num-boys-after-girls (count boy? after-girl))
(define pred49 (predictive-dist-from-p-samples post-samples 49))
(define ppost-plot49 (plot (density pred49)))


