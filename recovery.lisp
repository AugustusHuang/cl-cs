;;;; Algorithms based on compressive sensing
;;;; Author: Augustus Huang
;;;; Date: June 9, 2015

(in-package :cl-cs)

(defun hard-threshold (vec thres)
  "Hard threshold helper function."
  (declare (type vector vec))
  )

(defun soft-threshold (vec thres)
  "Soft threshold helper function."
  (declare (type vector vec))
  )

(defun threshold (vec thres)
  "Threshold helper function."
  (declare (type vector vec))
  )

(defun stagewise-omp (smatrix residual &key (stages 10) (err 1e-5) (sens 0.5))
  "Stagewise orthogonal matching pursuit, in Sparse Solution of Underdetermined
Linear Equations by Stagewise Orthogonal Matching Pursuit by D. L. Donoho,
Y. Tsaig, I. Drori and J-L. Starck."
  (declare (type sparse-matrix smatrix)
	   (type vector residual))
  (let ((row (array-dimension (sparse-matrix-values smatrix) 0))
	(col (array-dimension (sparse-matrix-values smatrix) 1))
	)))

(defun lasso (smatrix vec &key (stages 0) (err 1e-5) (lambda-stop 0) (residual-stop 0))
  "Lasso algorithm, in Least Angle Regression by B. Efron, T. Hastie
I. Johnstone and R. Tibshirani."
  (declare (type sparse-matrix smatrix)
	   (type vector vec))
  )