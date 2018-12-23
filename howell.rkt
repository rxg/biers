#lang racket

(provide all-defined-out)
;;
;; howell.rkt - load and examine Nancy Howell's height data
;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Real data (from Nancy Howell's data set)

;;
;; Experimenting with csv-reading and data-frame  Racket packages
;; (unfortunately the data-frame package's csv facilities are pretty weak)
;;
;; to load from the REPL or another module:
;; (require (submod "heights.rkt" howell))
;; or to load within this file::
;; (require (submod "." howell))

(module howell racket
  (require plot)
  (require csv-reading) ;; for reading McElreath's csv file
  ;; McElreath uses data frames a lot in R: are Racket's data frames helpful?
  ;; Possibly not, since R may either not have a tradition of using
  ;; rich data abstractions, or possibly just not have them at all?
  (require data-frame)  
  (provide (all-defined-out))
  
  ;;
  ;; Helper Functions
  ;;
  
  ;; transform a vector of row vectors into a vector of column vectors
  (define (rows->cols data)
    (let ([count (vector-length (vector-ref data 0))])
      (for/list ([idx (range count)])
        (for/vector ([row data]) (vector-ref row idx)))))

  ;; create a data frame packed with the provided series data
  (define (make-populated-data-frame series-names data-rows)
    (let ([df (make-data-frame)]
          [data-vecs (rows->cols data-rows)])
      (begin ;; imperative update, then return
        (for ([name series-names]
              [data data-vecs])
          (df-add-series df (make-series name #:data data)))
        df)))

  ;; low-level interface:
  ;; create a new data-frame by filtering an old one 
  (define (df-filter df filter-fn series-names)
    (let ([data-rows (apply df-select* df series-names  #:filter filter-fn)])
      (make-populated-data-frame series-names data-rows)))
  
  ;; lambda-like syntax for writing data frame filters
  ;; parameterized on series names.  E.g.: (λf (age) (>= age 18))
  (define-syntax λf
    (syntax-rules ()
      [(_ (nm* ...) exp)
       (λ (series-names)
         (λ (row)
           (let ([nm*
                  (vector-ref row
                              (index-of series-names (symbol->string 'nm*)))]
                 ...)
             exp)))]))

  ;; high-level interface:
  ;; create a new data-frame by filtering an old one 
  (define (filter-df filter df)
    (let* ([series-names (df-series-names df)]
           [fn (filter series-names)])
      (df-filter df fn series-names)))
  
  ;;
  ;; Data frames for Howell data
  ;;
  
  (define csv-reader
    (make-csv-reader (open-input-file "Howell1.csv")
                     '((separator-chars #\;))))
  
  (define csv-data (csv->list csv-reader))
  
  ;; get column headers interpreted as strings
  (define headers (first csv-data))
  
  ;; get row data interpreted as exact numbers (rest skips the header)
  ;; yields a vector of row vectors (compatible with df-filter helpers below)
  (define data-rows
    (for/vector ([row (rest csv-data)])
      (for/vector ([datum row])
        (string->number datum 10 'read 'decimal-as-exact))))

  (define df (make-populated-data-frame headers data-rows))
  
  (define heights (df-select df "height"))

  ;; Let's get the adults:
  (define adults-df (filter-df (λf (age) (>= age 18)) df))

  (define adult-heights (df-select adults-df "height"))
  
  (define adult-weight-height
    (for/list ([w (in-vector (df-select adults-df "weight"))]
               [h (in-vector (df-select adults-df "height"))])
      (vector w h)))
  ;; Scatter-plot of adult weight vs. height (see the correlation)
  (define plot-adult-weight-height (plot (points adult-weight-height)))
  ;; density plot of adult height (see that it looks normal)
  (define plot-adult-height-density
    (plot (density (df-select adults-df "height") 2)))  
  (void)) ; module howell
