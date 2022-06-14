#lang racket


;; -> -1 or +1
;; flip a coin, +1 for heads, -1 for tails
(define (flip) (if (zero? (random 2)) +1 -1))

;; Natural -> (listof Natural)
;; produce a running list of differences between heads and tails
(define (flips n0)
  (for/list ([n (in-range n0)]) (flip)))
        
(define (running-sum seq0)
  (let ([sum 0])
    (cons sum 
          (for/list ([n seq0])
            (set! sum (+ sum n))
            sum))))

(define (absy seq) (for/list ([n seq]) (abs n)))


(define (pts seq) (for/list ([x (in-range (length seq))]
                             [y seq])
                    (list x y)))

(define (do-winnings n) (pts (running-sum (flips n))))

(define (winnings-per-round seq)
  (for/list ([p seq])
    (match-let ([`(,x ,y) p]) `(,x ,(/ y (add1 x))))))


(define (absit n) (pts (absy (running-sum (flips n)))))


(require plot)

#;
(let ()
    (define winnings (do-winnings 100000))
    (define winnings/round (winnings-per-round winnings))
    (values
     (plot (points winnings) #:x-label "Round" #:y-label "Winnings")
     (plot (points winnings/round)
           #:x-label "Round" #:y-label "Winnings per round")))