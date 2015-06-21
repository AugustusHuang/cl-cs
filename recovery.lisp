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

(defun least-square (mat vec &key (stages -1))
  "Iterative Least Square method for solving linear equations."
  )

(defun stagewise-omp (smat svec &key (stages 10) (err 1e-5) (sens 0.5))
  "Stagewise orthogonal matching pursuit, in Sparse Solution of Underdetermined
Linear Equations by Stagewise Orthogonal Matching Pursuit by D. L. Donoho,
Y. Tsaig, I. Drori and J-L. Starck."
  (declare (type sparse-matrix smat)
	   (type sparse-vector svec))
  (let ((len (sparse-vector-len vec))
	(stage 1)
	(residual svec)
	(norm (vector-norm svec))
	(active ())
	(x-i ())
	(corr 0))
    (loop for i from stage to stages do
	 (setf corr (sparse-matrix*vec (sparse-matrix*const (sqrt len)
							    (transpose-sparse-matrix smat))
				       (normalize-vec residual)))
	 )))

(defun lasso (smat vec &key (stages 0) (err 1e-5) (lambda-stop 0) (residual-stop 0))
  "Lasso algorithm, in Least Angle Regression by B. Efron, T. Hastie
I. Johnstone and R. Tibshirani."
  (declare (type sparse-matrix smat)
	   (type vector vec))
  )