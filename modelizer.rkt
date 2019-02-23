#lang racket


;;
;; modelizer.rkt - Some tools for approximating Bayesian statistical models
;;
;; Author: Ron Garcia <rxg@cs.ubc.ca>
;; Date: February 22, 2019
;;


;; Model Specification:

;; v ∈ Variable
;; n ∈ Number
;; m ∈ Model
;; e ∈ Expression
;; d ∈ DataType
;; t ∈ TypeSpec
;; s ∈ Symbol
;; t ::= Index | Number | (Enum s ...) | (Range n n) | (Dist t)
;; m ::=
;; (new-model
;;   [types (d t) ...]
;;   [variables [v e]                    ;; single variable
;;              [(v d) d.e]              ;; indexed variable 
;;              [(v v ...) e) ...]]      ;; multivariate
;;   [data ...])
;; e ::= (normal e e)
;;     | (uniform e e)
;;     | (multivariate-normal #(e e ...)
;;                            #(#(e e ...)))
;;     | n
;;     | s
;;     | v
;;     | (v e)
;;     | (+ e e)
;;     | (* e e)
;;     | (expt e e)
;;     | (exp e)


;; Examples

;; Linear Regression

;; (new-model
;;   [types (i Index)]
;;   [variables [(h i) (normal (μ i) σ)]
;;              [(μ i) (+ α (* β (w i))]
;;              [α     (normal 0 1)]
;;              [β     (normal 0 1)]
;;              [(w i) (normal 80 5)]]
;;              [σ     (log-normal 100 100)]])
;;
;; Inference will expect to get height and weight inputs.
;; α, β, and σ are subject to inference

;; Stratified Linear Regression Model

;; (new-model
;;   [types (i Index) (m (enum 'male 'female))]
;;   [variables [(h i) (normal (μ i) σ)]
;;              [(μ i) (+ (α (g i)) (* β (w i))]
;;              [(α m) (normal 0 1)]
;;              [β     (normal 0 1)]
;;              [(w i) (normal 80 5)]]
;;              [(g i) (discrete ('male 1) ('female 1))]
;;              [σ     (log-normal 100 100)]])
;;
;; Inference will expect to get height and weight inputs.
;; α, β, and σ are subject to inference


;; Operations:
;; (draw-samples model [#:number n] variable-name ...)
;; variables should be in the same "tier" of variable
;; 
