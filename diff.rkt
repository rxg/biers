#lang racket

;;
;; diff.rkt - Numerical differentiation tools
;; Author: Ron Garcia <rxg@cs.ubc.ca>
;; Date: 26 Nov 2018
;;

(provide gradient scalar-gradient hessian hessian->covariance)

(module+ test
  (require rackunit)
  (define-syntax check-error
    (syntax-rules ()
      [(_ e) (check-exn exn:fail? (λ () e))])))

;;
;; Helpers
;;

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
  
;; "infinitesimal" perturbation for differentiation
;; This is the value that TI-8X calculators used.
;; (see citations at Wikipedia:Numerical_differentiation)
(define DIFF-EPSILON (make-parameter 0.001)) 


;; coerce a 1-element "vector" to a scalar
(define (scalar-demote x*)
  (unless (and (list? x*) (= (length x*) 1) (number? (first x*)))
    (error 'scalar-demote "Expected 1 element list result: got ~a" x*))
  (first x*))

(module+ test
  (check-equal? (scalar-demote '(7)) 7)
  (check-error (scalar-demote 7))
  (check-error (scalar-demote '(7 6))))

  
;; coerce any scalar to a 1-element "vector"
(define (scalar-promote x)
  (unless (number? x)
    (error 'scalar-promote "Expected a number: got ~a" x))
  (list x))

(module+ test
  (check-equal? (scalar-promote 7) '(7))
  (check-error (scalar-promote '(7))))


;; coerce a scalar-producing function to a "vector"-producing function
;; while preserving the arity of the original (compose does not do this)
(define (scalar-promote-fn f)
  (with-arity-of f (compose scalar-promote f)))

(module+ test
  (check-equal? ((scalar-promote-fn (λ () 7))) '(7))
  (check-error ((scalar-promote-fn (λ () '(7))))))
                                  

;; R^n {i ∈ N | i < n} R -> R^n
;; shift the ith position of a "vector" of numbers by a delta
(define (shift x* i δ)
  (unless (and (exact-nonnegative-integer? i) (< i (length x*)))
    (error 'shift "Bad index: ~a (vector length ~a)" i (length x*)))
  (for/list ([x x*] [j (in-range 0 (length x*))])
    (+ x (if (= i j) δ 0))))

(module+ test
  (check-equal? (shift (list 0 0 0 0) 3 2) (list 0 0 0 2))
  (check-error (shift (list 0 0 0 0) -3 -2))
  (check-error (shift (list 0 0 0 0) 0.5 -2)))


;;
;; Derivatives and Gradients
;; Adapted from differentiation example code in the Racket documentation
;; https://docs.racket-lang.org/htdp-langs/intermediate.html#\
;;  %28form._%28%28lib._lang%2Fhtdp-intermediate..rkt%29._check-range%29%29
;;


;; (R^n -> R^m) -> R^n -> R^m
;; m-ary ith partial derivative operator
(define (∂ i f)
  (let ([n (procedure-arity f)])
    (unless (exact-nonnegative-integer? n)
      (error '∂ "Bad procedure arity: ~a" n))
    (unless (and (exact-nonnegative-integer? i) (< i n))
      (error '∂ "Bad coordinate ~a for arity ~a\n" i n))
    (with-arity-of f
      (λ x*       
        (let* ([ε (DIFF-EPSILON)]
               [x-h (shift x* i (- ε))]
               [x+h (shift x* i (+ ε))])
          (for/list ([fl (apply f x-h)]
                     [fr (apply f x+h)])
            (/ (- fr fl)
               2 ε)))))))

(module+ test
  (define (sc n) (list (sin n) (cos n))) 
  (check-equal? (sc 0) '(0 1))
  (define d0sc (∂ 0 sc))
  (check-within (d0sc 0) '(1 0) (DIFF-EPSILON))
  (check-within (d0sc pi) '(-1 0) (DIFF-EPSILON))
  (check-within (d0sc (/ pi 2)) '(0 -1) (DIFF-EPSILON)))


;; (R^n -> R^1) -> (R^n -> R^n)
;; gradient operator
(define (gradient f)
  (define n (procedure-arity f))
  (with-arity-of f
    (λ x*
      (for/list ([i (in-range 0 n)])
        (scalar-demote
         (apply (∂ i f) x*))))))

(module+ test
  (local 
    [(define (fun x y) (list (+ (expt x 2) (* 3 y))))
     (define gf (gradient fun))]
    (check-within (gf 0 0) '(0.0 3.0) (DIFF-EPSILON))
    (check-within (gf 1 1) '(2.0 3.0) (DIFF-EPSILON))
    (check-within
     (parameterize ([DIFF-EPSILON 0.0001]) (gf 1 1))
     '(2.0 3.0)
     (DIFF-EPSILON))))


;; (R^n -> R) -> (R^n -> R^n)
;; operator that computes the numerical gradient of the given function
(define (scalar-gradient f)
  (gradient (scalar-promote-fn f)))

(module+ test
  (check-within ((scalar-gradient sin) 0) '(1) (DIFF-EPSILON))
  (check-within ((scalar-gradient cos) 0) '(0) (DIFF-EPSILON))
  (local [(define (fn x) (+ (* 5 x) 7))]
    (check-within ((scalar-gradient fn) 9) '(5) (DIFF-EPSILON))))

;; (R -> R) -> (R -> R)
;; operator that produces the derivative of f
(define (differentiate f)
  (with-arity-of f
    (compose scalar-demote (scalar-gradient f))))

(module+ test
  (check-within ((differentiate sin) 0) 1 (DIFF-EPSILON))
  (check-within ((differentiate cos) 0) 0 (DIFF-EPSILON))
  (check-within ((differentiate
                  (λ (x) (+ (* 5 x) 7))) 9)
                5
                (DIFF-EPSILON)))

;;
;; Example: Use numerical differentiation to understand the standard Gaussian
;;
(module* gaussian #f
  (require plot)
  ;; standard gaussian
  (define (sg x) (exp (* (- 1/2) (expt x 2))))
  (define dsg (differentiate sg))
  (define ddsg (differentiate dsg))
  (plot (list
         (axes)
         (function sg -3 3 #:label "Standard Gaussian" #:color 'black)
         (function (compose log sg) -3 3 #:label "Log Gaussian" #:style 'dot) 
         (function dsg -3 3 #:label "First derivative" #:color 'cyan)
         (function ddsg -3 3 #:label "Second derivative" #:color 'magenta))
        #:y-max 3
        #:title "Standard Gaussian, its Log, and its Derivatives")
  ;; TODO:  Multivariate Normal Distribution and sampling (fix Racket-ML)
  )

#;(require (submod "diff.rkt" gaussian))


;;
;; Jacobians and Hessians
;;


;; (R^n -> R^m) -> (R^n -> R^{nxm})
;; operator that produces transposed Jacobian values for f (good for Hessian)
(define (jacobian-trans f)
  (let ([n (procedure-arity f)])
    (with-arity-of f 
      (λ x*
        (for/list ([i (in-range 0 n)])
          (apply (∂ i f) x*))))))

(module+ test
  ;; stolen from Wikipedia:Jacobian Matrix and Determinant
  (define (f x y) (list (* (expt x 2) y)
                        (+ (* 5 x) (sin y))))
  (define Jf (jacobian-trans f))
  (define (analytic-Jf x y) (list (list (* 2 x y) 5)
                                  (list (expt x 2) (cos y))))
  (check-within (Jf 0 0) (analytic-Jf 0 0) (DIFF-EPSILON))
  (check-within (Jf 1 1) (analytic-Jf 1 1) (DIFF-EPSILON))
  (check-within (Jf 3 3) (analytic-Jf 3 3) (DIFF-EPSILON))
  (check-within (Jf 10 30) (analytic-Jf 10 30) (DIFF-EPSILON)))

  
;; (R^n -> R) -> (R^n -> R^{nxn})
;; approximate a Hessian (second derivative) operator for f
(define (wobbly-hessian f)
  (jacobian-trans (scalar-gradient f)))

(module+ test
  (local []
    ;; example from Khan Academy
    ;; https://www.khanacademy.org/math/multivariable-calculus/\
    ;; applications-of-multivariable-derivatives/modal/v/\
    ;; quadratic-approximation-example
    (define (f x y) (* (exp (/ x 2)) (sin y)))
    (define Hf (wobbly-hessian f))
    (check-within (Hf 0 (/ pi 2)) '[[0.25 0] [0 -1]] (DIFF-EPSILON))
    ;'((0.2500000208516262 0.0) (0.0 -0.9999996666842925))
    )

  (local []
    ;; example from Wikipedia
    ;; https://en.wikipedia.org/wiki/\
    ;; Taylor_series#Taylor_series_in_several_variables
    (define (g x y) (* (exp x) (log (add1 y))))
    (define Hg (wobbly-hessian g))
    (check-within (Hg 0 0) '[[0 1] [1 -1]] (DIFF-EPSILON))
    ; '((0 1.0000005000002088) (1.0000005000002392 -1.0000020000053844))
    )
  )


;; produce the transpose of a square matrix
(define (transpose m**)
  (let ([n (length m**)])
    (let loop ([m** m**])
      (cond [(empty? (rest m**)) (map list (first m**))]
            [else
             (let ([col* (loop (rest m**))])
               (for/list ([c (first m**)]
                          [col col*])
                 (cons c col)))])))
  )

(module+ test
  (check-equal? (transpose '((0 1) (2 3))) '((0 2) (1 3)))
  (check-equal? (transpose '((0 1 2) (2 3 4) (3 4 5)))
                '((0 2 3) (1 3 4) (2 4 5))))


;; transform ith matrix row to make "lower-triangle-dominant" symmetric matrix
(define (make-symmetric-row row col i)
  (let* ([n (length row)]
         [row-chunk (take row (add1 i))]
         [col-chunk (drop col (add1 i))])
    (append row-chunk col-chunk)))

(module+ test
  (check-equal? (make-symmetric-row '(0 1 2 3) '(0 2 4 6) 0) '(0 2 4 6))
  (check-equal? (make-symmetric-row '(1 2 3 4) '(0 2 4 6) 1) '(1 2 4 6))
  (check-equal? (make-symmetric-row '(2 3 4 5) '(0 2 4 6) 2) '(2 3 4 6))
  (check-equal? (make-symmetric-row '(3 4 5 6) '(0 2 4 6) 3) '(3 4 5 6)))


;; make symmetric matrix using mx's lower triangle
(define (symmetrify mx)
  (let ([tmx (transpose mx)])
    (for/list ([row mx] [col tmx] [i (in-range (length mx))])
      (make-symmetric-row row col i))))

(module+ test
  (check-equal? (symmetrify '((0 1) (3 4))) '((0 3) (3 4)))
  (check-equal? (symmetrify '((0 1 2) (3 4 5) (6 7 8)))
                '((0 3 6) (3 4 7) (6 7 8))))
  

(define (hessian f)
  (let ([Hf (wobbly-hessian f)])
    (with-arity-of Hf
      (compose symmetrify Hf))))

(module+ test
  (local []
    ;; example from Khan Academy
    ;; https://www.khanacademy.org/math/multivariable-calculus/\
    ;; applications-of-multivariable-derivatives/modal/v/\
    ;; quadratic-approximation-example
    (define (f x y) (* (exp (/ x 2)) (sin y)))
    (define Hf (hessian f))
    (define mx (Hf 0 (/ pi 2)))
    (check-within mx '[[0.25 0] [0 -1]] (DIFF-EPSILON))
    ;'((0.2500000208516262 0.0) (0.0 -0.9999996666842925))
    (check-equal? mx (transpose mx)))

  (local []
    ;; example from Wikipedia
    ;; https://en.wikipedia.org/wiki/\
    ;; Taylor_series#Taylor_series_in_several_variables
    (define (g x y) (* (exp x) (log (add1 y))))
    (define Hg (hessian g))
    (define mx (Hg 0 0))
    (check-within mx '[[0 1] [1 -1]] (DIFF-EPSILON))
    ; '((0 1.0000005000002088) (1.0000005000002392 -1.0000020000053844))
    (check-equal? mx (transpose mx)))
  )


;;
;; Covariance Matrix (for Laplacian Approximation)
;;

(require math/matrix)

;; (listof (listof Number)) -> Matrix
;; convert a list-based Hessian matrix to a Racket covariance matrix
(define (hessian->covariance mx)
  (let* ([hmx (list*->matrix mx)]
         [nhmx (matrix-scale hmx -1)]
         [covmx (matrix-inverse nhmx)]
         )
    covmx))



