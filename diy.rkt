#lang racket

;;
;; diy.rkt - hand-crafted implementations of some components that already exist
;;           in Racket.  Still useful for learning!
;;

(module+ test
  (require rackunit))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Urn stochastic sampler -- originally from urn.rkt
;; Color is 'blue | 'white


(define URN-SIZE 5)
;; UrnBlues is [0,URN-SIZE]
;; - number of blue balls in an urn of 5 balls

;; UrnBlues [0,URN-SIZE) -> Color
;; What is the color of the idxth ball (assuming blue ones are first)?
(define (get-ball-color urn-blues idx)
  (if (< idx urn-blues)
      'blue
      'white))

(module+ test
  (check-equal? (get-ball-color 0 0) 'white)
  (check-equal? (get-ball-color 3 0) 'blue)
  (check-equal? (get-ball-color 3 1) 'blue)
  (check-equal? (get-ball-color 3 2) 'blue)
  (check-equal? (get-ball-color 3 3) 'white)
  (check-equal? (get-ball-color 3 4) 'white)
  (check-equal? (get-ball-color 5 4) 'blue))


;; UrnBlues -> Color
;; sampling with replacement from the urn
(define (draw-sample urn-blues)
  (let ([index (random URN-SIZE)])
    (get-ball-color urn-blues index)))

;; UrnBlues n:Natural -> { l ∈ (listof Color) | (length l) = n }
(define (draw-samples urn-blues n)
  (for/list ([_ (in-range n)]) (draw-sample urn-blues)))

(module+ test
  ;; Urn Natural -> Rational
  ;; produce the blue rate of sampling from an urn count times
  ;; Use this function to sanity-check that drawing samples is well-behaved
  (define (rate urn count)
    (/ (length (filter (λ (c) (equal? c 'blue))
                       (draw-samples urn count)))
       count))

  (let* ([URN-BLUES 3]
         [TOTAL 1000]
         [TOLERANCE 1/10])
    (check-within (/ URN-BLUES URN-SIZE) (rate URN-BLUES 1000) TOLERANCE)))
                                              
                                             
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Formerly from globe.rkt

;; DEPRECATED!: quantile is already defined in math/statistics.
;; produce the UPPER boundary of the lower q percent of probability mass
;; (the LOWER boundary is 0).
(define (my-quantile q < samples)
  (let ([pos (inexact->exact (ceiling (* q (length samples))))])
    (list-ref (sort samples <) pos)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;