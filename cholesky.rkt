#lang racket

;;
;; cholesky.rkt - Cholesky decomposition
;; From https://rosettacode.org/wiki/Cholesky_decomposition
;; Copied and Adapted by: Ronald Garcia <rxg@cs.ubc.ca>
;; Date: December 10 2018  (added RackUnit Tests)
;;
;; Change Log:
;; - Modified cholesky to return a Racket matrix.
;; - Added RackUnit Tests, based on the Rosetta Code examples.
;;

(require math)
(provide cholesky)

(module+ test
  (require rackunit))

;; Matrix -> Matrix
;; produce the Cholesky Factor of A
(define (cholesky A)
  (define mref matrix-ref)
  (define n (matrix-num-rows A))
  (define L (for/vector ([_ n]) (for/vector ([_ n]) 0)))
  (define (set L i j x) (vector-set! (vector-ref L i) j x))
  (define (ref L i j) (vector-ref (vector-ref L i) j))
  (for* ([i n] [k n])
    (set L i k
         (cond 
           [(= i k) 
            (sqrt (- (mref A i i) (for/sum ([j k]) (sqr (ref L k j)))))]
           [(> i k) 
            (/ (- (mref A i k) (for/sum ([j k]) (* (ref L i j) (ref L k j))))
               (ref L k k))]
           [else 0])))
  (build-matrix n n (λ (i j) (ref L i j))))

(module+ test
  (define EPSILON 0.001)
  (check-within
   (cholesky (matrix [[25 15 -5]
                      [15 18  0]
                      [-5  0 11]]))
   (matrix [[5 0 0]
            [3 3 0]
            [-1 1 3]])
   EPSILON)

  (check-within
   (cholesky (matrix [[18 22  54 42]
                      [22 70  86 62]
                      [54 86 174 134]
                      [42 62 134 106]]))
   (matrix [[4.242640687119285 0 0 0]
            [5.185449728701349 6.565905201197403 0 0]
            [12.727922061357857 3.0460384954008553 1.6497422479090704 0]
            [9.899494936611665 1.6245538642137891
                               1.849711005231382 1.3926212476455924]])
   EPSILON)
  (void))