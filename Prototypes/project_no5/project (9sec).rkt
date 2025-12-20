#lang racket

(require flomat)

(define pretty-print (λ (lst)
    (for ([row (take-submatrix lst 5)])
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

(define path "/home/shlomi/Documents/PL_project/exp_data.txt")
(define delim #rx"([ ]*(,)[ ]*)|([ ]+)")

(define (next-line-it file)
  (let ((line (read-line file 'any)))
    (if (eof-object? line)
        (list)
        (cons line (next-line-it file)))))

(define read_data (λ(path delim)
      (define file (open-input-file path))
      (define lines (call-with-input-file path next-line-it))
      (define numbers (map (λ(x) (string-split x delim)) lines))
      (define mat (map (λ(x)(map string->number x)) numbers))
      (close-input-port file)
       mat))


(define (identity n)
  (build-flomat n n (λ (i j) (if (= i j) 1.0 0.0))))

(define (series-exp mat eps)
  (define n (nrows mat))
  (define I (identity n))

  (let loop ([S I]
             [T I]
             [k 1])
    (if (< (norm T 'max) eps)
        S
        (let*
            (
                [T-next (times (times T mat) (/ 1.0 k))]
                [S-next (plus S T-next)]
            )
            (loop S-next T-next (+ k 1))))))


;; Load the data
(display "Loading matrix...") (newline)
(define data (read_data path delim))
(display "Matrix loaded...") (newline)
(display "Formatting matrix...") (newline)
(define mat (lists->flomat data))
(display "Formatting done") (newline)
(define epsilon 0.0001)

;; Calculate with high precision
(display "Begins calculations...") (newline)
(define result (timing series-exp mat epsilon))

(display "Results are in:") (newline)
(display result)