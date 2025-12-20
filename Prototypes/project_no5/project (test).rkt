#lang racket

(require flomat)

(define m1 (eye 4 4 0)) ; Prints the matrix
(define m2 (eye 4 4 0)) ; Prints the matrix

(diag '(1 2))
(eig (diag '(1 2)))