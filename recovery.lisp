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

(defun stagewise-omp (smat vec &key (stages 10) (err 1e-5) (sens 0.5))
  "Stagewise orthogonal matching pursuit, in Sparse Solution of Underdetermined
Linear Equations by Stagewise Orthogonal Matching Pursuit by D. L. Donoho,
Y. Tsaig, I. Drori and J-L. Starck."
  (declare (type sparse-matrix smat)
	   (type vector vec))
  (let ((len (length vec))
	(stage 1)
	(residual vec)
	(norm (vector-norm vec))
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