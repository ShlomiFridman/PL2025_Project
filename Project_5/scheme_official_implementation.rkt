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

(display "Loading matrix from file...   ")
(define data (read-data path delim))
(define mat (lists->flomat data))
(display "Done\n\n")

(define res (timing expm mat))
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
