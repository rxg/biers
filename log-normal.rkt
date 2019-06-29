#lang racket

;;
;; log-normal.rkt - implementing the log normal distribution on top of
;;   the normal distribution
;;
(require math/distributions)
(provide log-normal-dist)

;;
;; A probably-broken attempt to make a log-normal distribution by hand
;;
(define (log-normal-pdf npdf)
  (λ (in . args)
    (apply npdf (cons (log in) args))))

(define (log-normal-sample nsmp)
  (λ args
    (cond
      [(empty? args) (exp (apply nsmp args))]
      [else (map exp (apply nsmp args))])))

(define (log-normal-dist mean stddev)
  (define nd (normal-dist mean stddev))
  (distribution (log-normal-pdf (distribution-pdf nd))
                (log-normal-sample (distribution-sample nd))))

  
