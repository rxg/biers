#lang racket

;;
;; utils.rkt - Various utility functions
;; Ronald Garcia <rxg@cs.ubc.ca>
;; Date: December 12, 2018
;;
(require (only-in math/statistics count-samples)
         (only-in plot discrete-histogram))

(provide interleave
         sample-hist-coords
         render-hist
         sum
         with-arity-of)

(module+ test
  (require rackunit)
  (define-syntax check-error
    (syntax-rules ()
      [(_ e) (check-exn exn:fail? (λ () e))])))


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

;; produce a copy of f that inherits the arity of f0 
(define (with-arity-of f0 f)
  (let ([arity (procedure-arity f0)])
    (procedure-reduce-arity f arity)))

(module+ test
  (let ([f0 (λ () 0)]
        [f1 (λ (x) 0)]
        [f2 (λ (x y) 0)]
        [f* (λ x* 0)])
    (check-equal? (procedure-arity (with-arity-of f0 f*))
                  (procedure-arity f0))
    (check-equal? (procedure-arity (with-arity-of f1 f*))
                  (procedure-arity f1))
    (check-equal? (procedure-arity (with-arity-of f2 f*))
                  (procedure-arity f2))
    (check-equal? (procedure-arity (with-arity-of f* f*))
                  (procedure-arity f*))))
