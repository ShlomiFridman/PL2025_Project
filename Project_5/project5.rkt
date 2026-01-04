#lang racket
(require flomat)

(define epsilon 1e-6)
(define scale_pow 16)

(define mat (matrix '((1 2) (3 4)) ))
(define n 2)

(define (calc-vec-norm v)
  (sqrt
   (dot v v)))


(define find-bk-rec (lambda (A b_curr) (
    let*
        (
            [b_next_vec (times A b_curr)]
            [b_next_norm (/ 1 (calc-vec-norm b_next_vec))]
            [b_next (times b_next_vec b_next_norm)]
            [b_diff (minus b_curr b_next)]
            [b_diff_norm_abs (abs (calc-vec-norm b_diff))]
        )
        (
            if (< b_diff_norm_abs epsilon)
            b_curr
            (find-bk-rec A b_next)
        )
)))

(define find-bk (lambda (A)
    (find-bk-rec A (ones n 1))
))

(find-bk mat)

