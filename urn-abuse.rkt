#lang racket

;; Transcripts from Suzanna and I toying with our Bayesian engine's mind

;; What if the "empirical phenomenon" is essentially coin flips, i.e., your
;; model of plausible phenomena is woefully wrong?

;; to be used with urn.rkt

(define (fake-sample)
  (let ([index (random 2)])
    (if (zero? index) 'white 'blue)))

(define (draw-fake-sample n)
  (map (λ (_) (fake-sample)) (range n)))

;; Contrary to Ron's expectations, the inference engine quickly
;; develops confidence as the mean rate of blue balls moves away from .5!
;; at .5, though, it IS split perfectly between 2 and 3

;; 25 white balls followed by 25 blue balls
(define sequential (append (map (λ (n) 'white) (range 25))
                           (map (λ (n) 'blue) (range 25))))

(define (ratio samples)
  (/ (length (filter (λ (c) (equal? c 'blue))
                     samples))
     (length samples)))

;; alternating 25 white and 25 blue balls
(define alternate
  (map (λ (n) (if (even? n) 'white 'blue)) (range 50)))

;; (confirm that order doesn't affect inference)
(equal? (last (simulate-bayes sequential))
        (last (simulate-bayes alternate)))