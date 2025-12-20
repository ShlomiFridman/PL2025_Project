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

;; -------------------------
;; Factorial helper
;; -------------------------
(define (factorial n)
  (if (zero? n) 1
      (* n (factorial (sub1 n)))))

;; -------------------------
;; Matrix exponential via power series
;; -------------------------
(define (matrix-exp-series A epsilon max-terms)
  (define N (nrows A))
  (define I (eye N N 0))
  (define result I)
  (define term I)
  
  (let loop ([k 1] [term I] [result I])
    (if (> k max-terms)
        result
        (let ([new-term (times (times term A) (/ 1 k))]
              [new-result (plus result (times (times term A) (/ 1 k)))])
          (if (< (norm new-term 'max) epsilon)
              new-result
              (loop (add1 k) new-term new-result))))))



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
(define result (timing matrix-exp-series mat epsilon 100000000000))

(display "Results are in:") (newline)
(display result)