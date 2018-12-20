#lang racket

;;
;; utils.rkt - Various utility functions
;; Ronald Garcia <rxg@cs.ubc.ca>
;; Date: December 12, 2018
;;

(provide interleave
         sample-hist-coords
         render-hist
         sum)

(require (only-in math/statistics count-samples)
         (only-in plot discrete-histogram))

;; interleave the contents of two lists
(define (interleave lsta lstb)
  (cond
    [(empty? lsta) lstb]
    [else (cons (first lsta)
                (interleave lstb (rest lsta)))]))

;; Hand-crafted discrete histogram tools
(define (sample-hist-coords samples)
  (let-values ([(x* y*) (count-samples samples)])
    (sort (map vector x* y*)
          (λ (a b) (< (vector-ref a 0) (vector-ref b 0))))))

(define (render-hist samples)
  (discrete-histogram (sample-hist-coords samples)))

(define (sum ls*)
  (apply + ls*))