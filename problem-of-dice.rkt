#lang racket

(require rackunit)

;; The Problem of Dice
;; From correspondence between Fermat and Pascal

;; Q: What is "fairness" in gambling
;; - what is fair proportion of stakes to take to forego a play?
;;   (Problem of Dice)
;; - what is fair proportion of stakes to take to leave a game unfinished?
;;   (Problem of Points)

;; Game: try in 8 dice rolls to get one 6

;; Probability
;; probability of winning in any round
(define P-WIN-ROUND 1/6)

;; Probability
;; probability of losing a round (= probability of not winning)
(define P-LOSE-ROUND (- 1 P-WIN-ROUND))

;; (m ∈ Natural) -> Probability
;; INVARIANT: 0 < m
;; probability of losing at rounds 1 through m 
(define (P-lose-at m) (expt P-LOSE-ROUND m))

(check-equal? (P-lose-at 1) 5/6)
(check-equal? (P-lose-at 2) (* 5/6 5/6))

;; (m ∈ Natural) -> Probability
;; INVARIANT: 0 < m
;; probability of winning on the mth round after losing rounds 1 to (m-1)
(define (P-win-at m)
  (cond [(= m 1)  P-WIN-ROUND]
        [else (* (P-lose-at (sub1 m)) P-WIN-ROUND)]))

(check-equal? (P-win-at 1) 1/6)
(check-equal? (P-win-at 2) (* 5/6 1/6))
(check-equal? (P-win-at 3) (* 5/6 5/6 1/6))

;; (n ∈ Natural) -> Probability
;; INVARIANT: 0 < n
;; probability of winning in a game within n rounds
(define (P-win-in n)
  (cond [(= n 1) (P-win-at 1)]
        [else (+ (P-win-at n) (P-win-in (sub1 n)))]))

(check-equal? (P-win-in 1) P-WIN-ROUND)
(check-equal? (P-win-in 2) (+ (P-win-at 1) (P-win-at 2)))
(check-equal? (P-win-in 3) (+ (P-win-at 1) (P-win-at 2) (P-win-at 3))) 

;; Claim: for all n>0, (P-win-in n) + (P-lose-at n) = 1
;; try it 10 times
(for ([repetitions (range 10)])
  (let ([n (add1 (random 20))])
    (check-equal? (+ (P-win-in n) (P-lose-at n)) 1)))