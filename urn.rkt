#lang racket

;; urn.rkt:
;;  5-ball urn statistical process:
;;  - A Sampling Engine 
;;  - An Inference Engine
;;  - Some Viz/Rendering tools (under construction)
(require math/distributions)
(require pict)
(require plot)
(require "utils.rkt")

(module+ test
  (require rackunit))

;; Color is 'blue | 'white

;;
;; Stocastic Simulator
;;

(define URN-SIZE 5)
;; UrnBlues is [0,URN-SIZE]
;; - number of blue balls in an urn of 5 balls

;; UrnBlues n:Natural -> { l ∈ (listof Color) | (length l) = n }
;; draw n uniformly random samples, with replacement, from the given urn
(define (draw-samples urn-blues n)
  (let ([urn-whites (- URN-SIZE urn-blues)])
    (let ([dist (discrete-dist (list 'blue 'white)
                               (list urn-blues urn-whites))])
      (sample dist n))))

(define (draw-samples-vector urn-blues n)
  ; using list->vector is MUCH faster than using build-vector and single samples
  (list->vector (draw-samples urn-blues n)))

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


;;
;; Inference Engine -- Given a list of urn samples, predict how many
;;   blue balls are in the urn
;;

;; Random Variables:
;; C=c - a ball is 'white or 'blue
;; B=i - the urn has i blue balls out of 5
;; Cs='(c1 ...) - list of sampled ball colors in order are '(c1 ...)
;;    (well...technically this is a Nat-indexed set of random variables
;;    Cs_i, each representing the results of the first *i* samples)

;; The inference engine incrementally updates the model's "beliefs"
;; P(B=i|Cs='(c1 ...))


;; Likelihoods -- these remain constant

;; P(C='blue|B=i) 
;; Factor for updating your prior on seeing a blue ball
;; i.e., likelihood that you draw a blue ball if the urn has i blue balls;
(define PBI*
  '(0 1/5 2/5 3/5 4/5 5/5))

;; P(C='white|B=i) 
;; Factor for updating your prior on seeing a white ball 
;; i.e., likelihood that you draw a white ball if the urn has i blue balls;
(define PWI*
  '(5/5 4/5 3/5 2/5 1/5 0))

(define (likelihood c i)
  (case c
    [(white) (list-ref PWI* i)]
    [(blue) (list-ref PBI* i)]))

;; Initial prior -- deduction must start from some assumption.

;; P(B=i|Cs='())
;; Given no data, we assume no preconceived inclination toward any
;; urn structure, modeled with uniform initial probabilities.
(define uniform-prior '(1/6 1/6 1/6 1/6 1/6 1/6))

;; Sequential Bayesian Inference: see
;; https://stats.stackexchange.com/questions/244396/\
;; bayesian-updating-coin-tossing-example
;; as well as
;; http://www.stats.ox.ac.uk/~steffen/teaching/bs2HT9/kalman.pdf

;; P(B=i|Cs='(c1 ... . cn)) =
;; P(C=cn|B=i)P(B=i|Cs='(c1 ...)) /  (Σ_i P(C=cn|B=i)P(B=i|Cs='(c1 ...)))

;; In particular, P(C=cn|B=i,Cs='(c1 ...)) = P(C=cn|B=i) because
;; C ⊥ Cs | B  (conditional independence)
;; otherwise we'd have to deal with P(C=cn|B=i,Cs='(c1 ...))
;; https://en.wikipedia.org/wiki/Conditional_independence

;; Note: we're being explicit about accumulating a sequence of
;; observations at each step,rather than saying that:
;; "the posterior P(B=i|C=c) becomes the new prior P(B=i) in the next step."
;; :: shudder! ::  

;; average likelihood
;; P(C=c|Cs='(c1 ...)) = (Σ_i P(C=c|B=i)P(B=i|Cs='(c1 ...)))
;; C and Cs are NOT independent, only *conditionally* independent
(define (pcc* pci* pi*)
  (sum (map * pci* pi*)))

;; Initial P(C='white|Cs='()): 1/2  (changes as your observations change)
#;(define pw0 (pcc* PWI* uniform-prior))

;; Initial P(C='blue|Cs='()): 1/2   (changes as your observations change
#;(define pb0 (pcc* PBI* uniform-prior))

;; perform an inference step given P(C=c|B=i), P(B=i|Cs='(c1 ...))
(define (posterior pci* pi*)
  (let ([the-pcc* (pcc* pci* pi*)])
    (map (λ (pcb pb) (/ (* pcb pb) the-pcc*))
         pci* pi*)))

;; likelihood distribution given a ball color
(define (pci* c)
  (cond
    [(equal? c 'blue) PBI*]
    [(equal? c 'white) PWI*]))

;; update belief based on an observation
;; c:Color pi:P(B=i|Cs='(ci ...)) -> P(B=i|Cs='(ci ... . c))
(define (fit c pi*)
  (posterior (pci* c) pi*))

;; (listof Color) -> (listof P(B=i))
;; given a list of samples, trace sequential Bayesian inference.
(define (simulate-bayes c*0)
  (let loop ([c* c*0] [pi* uniform-prior])
    (cond
      [(empty? c*) (list pi*)]
      [else
       (cons pi* (loop (rest c*)
                       (fit (first c*) pi*)))])))


;; (listof Color) -> Urn
;; given a list of samples, predict the contents of the urn
;; maximum a posterior
(define (predict c*)
  (let ([stats (last (simulate-bayes c*))])
    (argmax (λ (i) (list-ref stats i)) (range URN-SIZE))))

(define (wifey-demands-units c*)
  (let ([n (predict c*)])
    (printf "~a blue balls" n)))


;;
;; Rendering tools - Makes the wifey happy with visualizations
;; using Racket Plot Library
;;

;; Color -> Pict
;; represent the given ball visually
(define (render-ball color)
  (disk 20 #:color (symbol->string color)))

;; (listof Color) -> Pict
;; lay out a list of balls horizontally
(define (render-balls color*)
  (apply ht-append (map render-ball color*)))


(plot-width 120)
(plot-height 120)
;; Render probability mass function 
(define (plot-probs p*)
  (let* ([x-lbls (range 6)]
         [entries (map vector x-lbls p*)])
    (plot (discrete-histogram entries)
          #:y-max 1 #:x-label #f #:y-label #f)))
  

;; (listof Color) -> Pict
;; render a sequence of Bayesian inference steps
(define (render-bayes c*0)
  (let ([b* (map render-ball c*0)])
    (let ([chart* (map plot-probs (simulate-bayes c*0))])
      (interleave chart* b*))))

;;
;; Example Use
;;
#;
(begin
  (define blues 3) ;; 3 blue and 2 white balls in urn
  (define samples (draw-samples blues 10)) ;; draw 10 samples with replacement
  (render-balls samples) ;; visualize the samples
  (render-bayes samples) ;; visualize sequential Bayesian updating
  (void))