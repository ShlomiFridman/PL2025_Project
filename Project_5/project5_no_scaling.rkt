#lang racket
(require flomat)

(define timing (λ(F . Arg)
                 (let* [
                        [t_start (current-inexact-milliseconds)]
                        [ret (apply F Arg)]
                        [t_end (current-inexact-milliseconds)]
                        [t (- t_end t_start)]
                       ]
                       [list ret t]
)))

(define path "./exp_data.csv")
(define delim #rx"([ ]*(,)[ ]*)|([ ]+)")

(define next-line-it (λ(file) (
    let ((line (read-line file 'any)))
        (if (eof-object? line)
            (list)
            (cons line (next-line-it file))
        )
)))

(define read-data (λ(path delim) (
    let*
        (
            [file (open-input-file path)]
            [lines (call-with-input-file path next-line-it)]
            [numbers (map (λ(x) (string-split x delim)) lines)]
        )
        (map (λ(x)(map string->number x)) numbers)
)))


(define epsilon 1e-6)

(define find-bk-rec (λ(A b_curr) (
    let*
        (
            [b_next_vec (times A b_curr)]
            [b_next_norm (/ 1 (norm b_next_vec))]
            [b_next (times b_next_vec b_next_norm)]
            [b_diff (minus b_curr b_next)]
            [b_diff_norm (norm b_diff)]
        )
        (
            if (< b_diff_norm epsilon)
                b_curr
                (find-bk-rec A b_next)
        )
)))

(define find-bk (λ(A n)
    (find-bk-rec A (ones n 1))
))

(define calc-abs-mu (λ(A bk) (
    let*
    (
        [bk-t (transpose bk)]
        [m1 (times bk-t A)]
        [m2 (times m1 bk)]
        [m3 (ref m2 0 0)]
        [bk-norm (dot bk bk)]
        [norm (/ m3 bk-norm)]
    )
    (abs norm)
)))

(define find-m-rec (λ(mu-org mu-acc m-fact m) (
    if (< (/ mu-acc m-fact) epsilon)
        m  
        (let*
            (
                [mu-acc-next (* mu-acc mu-org)]
                [m-next (+ m 1)]
                [m-fact-next (* m-fact m-next)]
            )
            (find-m-rec mu-org mu-acc-next m-fact-next m-next)
        )
)))

(define find-m (λ(A bk) (
    let
        ([mu-norm (calc-abs-mu A bk)])
        (find-m-rec mu-norm mu-norm 1 1)
)))

(define calc-exp-mat-rec (λ(A-org A-acc k-fact res-mat loop-ind end-ind) (
    if (= loop-ind end-ind)
    res-mat
    (
        let*
            (
                [A-acc-next (times A-acc A-org)]
                [k-fact-next (* k-fact loop-ind)]
                [A-to-add (times A-acc-next (/ 1 k-fact-next))]
                [res-mat-next (plus res-mat A-to-add)]
                [loop-ind-next (+ loop-ind 1)]
            )
            (calc-exp-mat-rec A-org A-acc-next k-fact-next res-mat-next loop-ind-next end-ind)
    )
)))

(define calc-exp-mat (λ(A n m) (
    let
        (
            [eye-n (eye n)]
            [m1 (+ m 1)]
        )
        (calc-exp-mat-rec A eye-n 1 eye-n 1 m1)
)))

(define run-prog (λ(A n) (
    let*
        (
            [bk (find-bk A n)]
            [m (find-m A bk)]
        )
        (calc-exp-mat A n m)
)))

(display "Loading matrix from file...   ")
(define data (read-data path delim))
(define mat (lists->flomat data))
(define mat-n (nrows mat))
(display "Done\n\n")

(define res (timing run-prog mat mat-n))
(define res-mat (car res))
(define res-time (car (cdr res)))

(display "Input sub-matrix (5x5):\n")
(display (sub mat 0 0 5 5))
(newline)
(newline)

(display "Total time (sec): ")
(display (/ res-time 1000))
(newline)
(newline)


(display "Result sub-matrix (5x5):\n")
(display (sub res-mat 0 0 5 5))
(newline)
