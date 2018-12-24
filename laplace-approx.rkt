#lang racket

(provide lapprox lapprox-dimension lapprox-means lapprox-covariance lapprox-dist
         laplace-approx laplace-approx-at maximize-fn
         lapprox-sample lapprox-pdf ar-compose)

(require nlopt/highlevel) ;; for optimize/args
(require nlopt/unsafe) ;; low-level interface
(require math/flonum)
(require "diff.rkt") ;; function differentiation
(require "racket-ml/multivariate-normal.rkt") ;; normal distribution
(require math/distributions) ;; for uniform-dist
(require math/matrix)
(require ffi/unsafe)
;;
;; laplace-approx.rkt - Laplacian Approximation of posterior distributions
;;
;; Author: Ron Garcia <rxg@cs.ubc.ca>
;; Date: December 15, 2018
;;

;; Normal curves are useful for approximating posterior distributions resulting
;; from Bayesian Inference. Inspired by McElreath, Statistical Rethinking.
;;
;; References:
;; Gelman et al. Bayesian Data Analysis , 3rd Edition
;; A few web pages:
;; http://www.stat.cmu.edu/~ryurko/post/bayesian-baby-steps-normal-next-steps/
;; http://www.sumsar.net/blog/2013/11/easy-laplace-approximation/


;; Helper Function
;; arity-respecting compose function (because compose does not...)
(define (ar-compose f g)
  (procedure-reduce-arity 
   (compose f g)
   (procedure-arity g)))

(define-struct lapprox (dimension means covariance dist))
;; Lapprox is (make-lapprox n means cov dist)
;; Laplacian approximation of a posterior distribution
;; n is the dimensionality of the distribution
;; means is the list of means
;; cov is the covariance matrix
;; dist is the multivariate normal distribution object
;; technically all we need to store is means and cov, but whatever...


;; given a log-posterior bounds, and an initial search parameter, produce
;; a Laplacian approximation of the posterior
(define (laplace-approx fn bounds [maybe-x0 #f] #:method [method 'GN_DIRECT_L])
  (define mode (maximize-fn fn bounds maybe-x0 #:method method))
  (laplace-approx-at fn mode))

;; (R^n -> R) R^n -> Lapprox
;; produce Laplace approximation for log-posterior fn at given mode coordinates
(define (laplace-approx-at fn mode)
  (define cov (hessian->covariance (apply (hessian fn) mode)))
  (make-lapprox (length mode)
                mode
                cov
                (multivariate-normal-dist (->col-matrix mode) cov)))


;; (R^n -> R) (listof (list R R)) -> R^n
;; RG: nlopt interface appears to be buggy.  Check FFI.
;; Find an input to fn that maximizes it, starting from x0 or bounds center
(define (maximize-fn fn bounds [maybe-x0 #f] #:method [method 'GN_DIRECT_L])
  (define n (procedure-arity fn))
  (unless (number? n)
    (error 'maximize-fn "function argument has non-numeric arity ~a." n))
  (unless (= (sequence-length bounds) n)
    (error 'maximize-fn "function has arity ~a but only ~a bounds provided"
           n (length bounds)))

  ;; Use the provided initial coordinate, else the center of the bounds
  (define x0
    (if maybe-x0
        (begin
          (and (= (length maybe-x0) (length bounds))
               (for/list ([x maybe-x0] [bound bounds])
                 (unless (<= (first bound) x (second bound))
                   (error 'maximize-fn
                          "bound violation: ~a <= ~a <= ~a"
                          (first bound) x (second bound)))))
          maybe-x0)
        (for/list ([bound bounds])
          (let ([x (* 0.5 (+ (first bound) (second bound)))])
            (when (member x (list +inf.0 -inf.0))
              (error 'maximize-fn
                     "Cannot deduce initial guess: infinite bound detected."))
            x))))
  #;(printf "x0 = ~a\n" x0)
  (define bounds-pairs (map (λ (p) (cons (first p) (second p))) bounds))
  ;; globally search for the ballpark
  (define-values
    (gvalue gx-max)
    (optimize-fn #;optimize/args fn x0
                   #:maximize #t
                   #:bounds bounds-pairs
                   #:method method))
  gx-max)

;; RG - turn this into a test case
#;
(begin
  (maximize-fn (λ (a b) (+ a b)) `((0 ,pi) (1  5)))
  (void))


;; The following is adapted from code in nlopt/safe.rkt
;; (Flonum ... -> Flonum) ->
;; _int (_ptr i _double) (_ptr o _double) _pointer -> Flonum
(define (wrap-as-raw-fn f)
  (lambda (n raw-x raw-grad data)
    (define x-flv (make-flvector n))
    ;; don't use gradient-based methods, and don't expect data
    (when raw-grad
      ;; RG - I wish I could call force-stop here, but I have no opt!
      ;; not hard to fix: pass opt to wrap-as-raw-fn...
      (error 'wrappy "Gradient expected"))
    (memcpy (flvector->cpointer x-flv) raw-x n _double)
    (define x (flvector->list x-flv))
    (real->double-flonum (apply f x))))


;;
;; Attempt to use the low-level interface to replace optimize/args
;; interface made to match optimize/args
;;
;; Seeming bugs in nlopt/highlighlevel and nlopt/unsafe led to me needing
;; max-wrap to transport results past memory errors
;; (though they were still memory errors)
;; pull requests filed with nlopt.

#;
(define (max-wrap f)
  (define arity (procedure-arity f))
  (define vmax (box -inf.0))
  (define x (box #f))
  (define (wrap-fn . x*)
    (define v (apply f x*))
    (when (< (unbox vmax) v)
      (set-box! vmax v)
      (set-box! x x*))
    v)
  (define (report) (values (unbox vmax) (unbox x)))
  (values wrap-fn report))


;; debugging diagnostics macro
(define-syntax dprintf
  (syntax-rules ()
    #;[(_ e ...) (printf e ...)]
    [(_ e ...) (void)]
    ))


;; Custom replacement for nlopt's optimize/arg. Instrumented for debugging, and
;; does not compute gradients, so runs faster.  
(define (optimize-fn fn x0-list
                     #:maximize maximize #:bounds bounds #:method method)
  (define arity (procedure-arity fn))
  (define opt (create method arity)) ;; CREATE
  (define lb* (list->flvector (map car bounds))) 
  (define ub* (list->flvector (map cdr bounds)))
  (define x0 (list->flvector x0-list))
  #;(define-values (max-fn report) (max-wrap fn))
  (define wrap-fn (wrap-as-raw-fn fn #;max-fn))
  (set-lower-bounds opt (flvector->cpointer lb*)) ;; SET-LOWER-BOUNDS
  (set-upper-bounds opt (flvector->cpointer ub*)) ;; SET-UPPER-BOUNDS
  (set-max-objective opt wrap-fn #f) ;; SET-MAX-OBJECTIVE
  (set-ftol-rel opt (real->double-flonum 1e-8)) ;; SET-FTOL-REL
  #;(set-xtol-rel opt (real->double-flonum 1e-8)) ;; SET-XTOL-REL
  ;; print out debug info
  (dprintf "arity: ~a\n" arity)
  (let ()
    (define lbs (make-flvector (get-dimension opt)))
    (define res (get-lower-bounds opt (flvector->cpointer lbs)))
    (dprintf "lower-bounds: ~a ~a\n" lbs res))
  (let ()
    (define ubs (make-flvector (get-dimension opt)))
    (define res (get-upper-bounds opt (flvector->cpointer ubs)))
    (dprintf "upper-bounds: ~a ~a\n" ubs res))
  (dprintf "x0: ~a\n" x0)
  (dprintf "method: ~a\n" method)
  (dprintf "stopval: ~a\n" (get-stopval opt))
  (dprintf "ftol-rel: ~a\n" (get-ftol-rel opt))
  (dprintf "ftol-abs: ~a\n" (get-ftol-abs opt))
  (dprintf "xtol-rel: ~a\n" (get-xtol-rel opt))
  (let ()
    (define xtols (make-flvector (get-dimension opt)))
    (define res (get-xtol-abs opt (flvector->cpointer xtols)))
    (dprintf "xtol-abs: ~a ~a\n" xtols res))
  (dprintf "maxeval: ~a\n" (get-maxeval opt))
  (dprintf "maxtime: ~a\n" (get-maxtime opt))
  (define x0-carray (malloc _double arity 'atomic-interior))
  (memcpy x0-carray (flvector->cpointer x0) arity _double)
  (define-values (flag value) (optimize opt x0-carray)) ;; OPTIMIZE
  (define x-opt (make-flvector arity))
  (memcpy (flvector->cpointer x-opt) x0-carray arity _double)
  ;; print out result information
  (dprintf "RESULTS\n")
  (dprintf "result: ~a\n" flag)
  (dprintf "value: ~a\n" value)
  (dprintf "x-opt: ~a\n" x-opt)
  (dprintf "num-evals: ~a\n" (get-numevals opt))
  #;(define-values (max-value max-x*) (report))
  #;(dprintf "max-value: ~a\n" max-value)
  #;(dprintf "max-x*: ~a\n" max-x*)
  (values value (flvector->list x-opt)))



;; print some stats on the distribution
;; TODO: IMPLEMENT!
(define (precis la)
  (void))

;; Lapprox -> (listof R^n) where n = (lapprox-dimension la)
;; ToDo: scalarize the one-dimensional case?
(define (lapprox-sample la n)
  (map matrix->list (sample (lapprox-dist la) n)))

;; Lapprox -> (R^n -> [0,1])  where n = (lapprox-dimension la)
;; produce a probability distribution function
(define (lapprox-pdf la)
  (λ x*
    (pdf (lapprox-dist la) (->col-matrix x*))))