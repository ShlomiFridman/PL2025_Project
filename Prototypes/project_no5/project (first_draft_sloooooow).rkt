#lang racket

(define pretty-print (λ (lst)
    (for ([row (take-submatrix mat 5)])
    (display row)
    (newline))
))
(define take-submatrix (λ (mat n)
  (take (map (lambda (row) (take row n)) mat) n)
))
(define last-element (λ (lst)
  (list-ref lst (- (length lst) 1))
))
(define timing (λ(F . Arg)
                 (let* [
                        [t_start (current-inexact-milliseconds)]
                        [ret (apply F Arg)]
                        [t_end (current-inexact-milliseconds)]
                        [t (- t_end t_start)]
                       ]
                       [list ret t]
)))

(define mat-transpose (λ (mat)
    (apply map list mat)
))

(define dot-prod (λ (row col)
    (apply + (map * row col))
))

(define mat-mul (λ (A B)
    (let*
        [(Bt (mat-transpose B))]
        [map (λ (a-row)
            (map (λ (b-col)
                    (dot-prod a-row b-col))
                Bt))
            A]
    )
))

(define mat-add (λ (A B)
    (map (λ (row-a row-b)
         (map + row-a row-b))
       A B)
))

(define mat-scale (λ (k mat)
  (map (λ (row)
         (map (λ (x) (* k x)) row))
       mat)
))

(define gen-i-mat (λ (n)
    (build-list n (λ (i)
        (build-list n (λ (j) (if (= i j) 1 0)))
    ))
))

(define mat-max-norm (λ (mat)
  (apply max (map (λ (row) 
                (apply max (map abs row))) 
              mat))
))

(define mat-exp (λ (A epsilon)
    (let ([n (length A)])
        (let loop
            (
                [k 1]
                [current-term (gen-i-mat n)]
                [total-sum (gen-i-mat n)]
            )
            (let*
                (
                    [next-term-raw (mat-mul current-term A)]
                    [next-term (mat-scale (/ 1.0 k) next-term-raw)]
                    [term-magnitude (mat-max-norm next-term)]
                )
                (if (< term-magnitude epsilon)
                    (mat-add total-sum next-term)
                    (loop
                        (+ k 1)
                        next-term
                        (mat-add total-sum next-term)
                    )
                )
            )
        )
    )
))

(define read-csv (λ (filename)
  (call-with-input-file filename
    (lambda (in)
      (for/list ([line (in-lines in)])
        (map string->number
             (string-split line ","))
        )))
))

(display "Loading data...") (newline)
(define mat (read-csv "exp_data.csv"))
(display "Data loaded.") (newline)

(display "Displaying the 5x5 sub-matrix") (newline)
(pretty-print (take-submatrix mat 5))

;; Calculate with high precision
(display "Begins calculations...") (newline)
(define result (timing mat-exp mat 0.001))

(display "Results are in, the time it took in sec:") (newline)
(display (last-element result))