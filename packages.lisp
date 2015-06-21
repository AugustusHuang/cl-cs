;;;; General compressive sensing packages
;;;; Author: Augustus Huang
;;;; Date: June 9, 2015

(defpackage :cl-cs
  (:nicknames :cs)
  (:use :cl :cl-cs-sparse :cl-cs-util)
  (:export :least-square
	   :stagewise-omp
	   :lasso
	   ))

(defpackage :cl-cs-sparse
  (:nicknames :cs-sparse)
  (:use :cl :cl-cs-util)
  (:export :make-sparse-matrix-with-matrix
	   :make-sparse-vector-with-vector
	   :aref-sparse-matrix
	   :aref-sparse-vector
	   :transpose-sparse-matrix
	   :negative-sparse-matrix
	   :negative-sparse-vector
	   :sparse-m-*
	   :sparse-m-+
	   :sparse-m--
	   :sparse-v-*
	   :sparse-v-+
	   :sparse-v--
	   :sparse-matrix-*-vector
	   :sparse-matrix-*-sparse-vector
	   :matrix-*-sparse-vector
	   :sparse-matrix-*-const
	   :sparse-inner-product))

(defpackage :cl-cs-util
  (:nicknames :cs-util)
  (:use :cl)
  (:export :dovec
	   :doseq
	   :random-array
	   :random-matrix
	   :with-gensyms
	   :list-dimensions
	   :list-to-array
	   :1d-array-to-list
	   :ignore-trailing-zero
	   :matrix-invert
	   :matrix-determinant
	   :inner-product
	   :vector-norm
	   :sparse-vector-norm
	   :vector-abs
	   :sparse-vector-abs))
