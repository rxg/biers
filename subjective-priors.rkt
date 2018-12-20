#lang racket

;;
;; subjective-priors.rkt
;;   Based on a scenario in the appendix of The Theory That Would Not Die
;;   In the absence of much data, the (subjective) prior dominates the data.
;;

(require plot)


;; F = Coin is False (two Heads) versus legitimate (not False)
;; D = HHH (three flips yield heads)


;; P(F|D) where prior P(F) ∈ [0,1] is the input
;; (pfd 1/2) = 8/9  (pretty convinced!)
;; (pfd 0.11) = 0.497
;; (pfd 0.12) = 0.522
(define (pfd pf)
  (let (;; P(not F)
        [pnf (- 1 pf)]
        ;; P(D|F)
        [pdf 1]
        ;; P(D|not F)
        [pdnf 1/8])
    (/ (* pdf pf)
       (+ (* pdf pf) (* pdnf pnf)))))

(plot (function pfd 0 1) #:x-label "P(F)" #:y-label "P(F|D)")
