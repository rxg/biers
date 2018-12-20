#lang racket


;;
;; dice.rkt - Simulation of The Dice Problem
;;  Adapted from "Pascal and the Invention of Probability Theory" by Oystein Ore
;;

;; How many times must you throw two dice before you have a better-than 50% odds
;; of rolling two 6's.  Answer is solution to:
;; min{ n ∈ N | p(n) > 1/2 }
;; where
;; p(n) = 1 - (35/36)^n
;; turns out the answer is 25, where p(25) = .5055

;; Let's simulate playing this out:

(define (p n)
  (- 1 (expt 35/36 n)))

(define (iterate-dice)
  (let loop ([n 1]) ;; n - minimum number of rolls it might take this time
    (let ([d1 (add1 (random 6))]
          [d2 (add1 (random 6))])
      (if (and (= d1 6) (= d2 6))
          n
          (loop (add1 n))))))

;; Interesting!  The expected value of "number of rolls before double-six"
;; does not seem to be 25, but rather around 35 or 36!
(define (average-number-rolls i)
  (let ([rolls (map (λ (_) (iterate-dice)) (range i))])
    (/ (apply + rolls) i)))

;; But let's instead turn this into win/lose with a MAX of N.  That seems
;; more like the right "game"
(define (play rounds)
  (let loop ([n 0])
    (if (= n rounds)
        #f
        (let ([d1 (add1 (random 6))]
              [d2 (add1 (random 6))])
          (if (and (= d1 6) (= d2 6))
              #t
              (loop (add1 n)))))))

;; Ahh, ok given 25 rounds minimum to win, seems to hit about .5
(define (average-wins rounds i)
  (let ([wins (for/list ([_ (in-range i)]) (if (play rounds) 1 0))])
    (/ (foldr + 0 wins) i)))

        