;;;; Sparse matrix representation
;;;; Author: Augustus Huang
;;;; Date: June 7, 2015

(in-package :cl-cs)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defstruct sparse-matrix
    (values (make-array '(0 0)) :type matrix)
    (present-p (make-array '(0 0) :element-type '(mod 2)))
    (sparsity 0 :type fixnum)
    (capacity 0 :type fixnum)))

;;; Don't use make-sparse-matrix, use this instead!
(declaim (inline make-sparse-matrix-with-element))
(defun make-sparse-matrix-with-element (row column &key (element-type t))
  (make-sparse-matrix
   :values (make-array `(,row ,column) :element-type element-type)
   :present-p (make-array `(,row ,column) :element-type (quote (mod 2)) :initial-element 0)
   :capacity (* row column)))

(defun aref-sparse-matrix (smat row column)
  "Returns the (row, column)th item in the sparse-matrix."
  (declare (type sparse-matrix smat)
	   (type fixnum row column))
  (aref (sparse-matrix-values smat) row column))

(defun (setf aref-sparse-matrix) (smat row column value)
  "Sets the (row, column)th item in the sparse-matrix with 'value'."
  (declare (type sparse-matrix smat)
	   (type fixnum row column))
  (progn
    (setf (aref (sparse-matrix-values smat) row column) value
	  (aref (sparse-matrix-present-p smat) row column) 1)
    (incf (sparse-matrix-sparsity smat))))

(defun resize-sparse-matrix (smat row column)
  "Removes the (row, column)th item in the sparse-matrix."
  (declare (type sparse-matrix smat))
  (progn
    (setf (aref (sparse-matrix-values smat) row column) 0
	  (aref (sparse-matrix-present-p smat) index) 0)
    (decf (sparse-matrix-sparsity smat))))

(defun transpose-sparse-matrix (smat)
  "Transposes a sparse-matrix."
  (declare (type sparse-matrix smat))
  (let* ((row (array-dimension (sparse-matrix-values smat) 0))
	 (col (array-dimension (sparse-matrix-values smat) 1))
	 (element-type (type-of (aref (sparse-matrix-values smat) 0 0)))
	 (sout (make-sparse-matrix-with-element `(,col ,row) :element-type element-type)))
    (loop for i from 0 to (1- row) do
	 (loop for j from 0 to (1- col) do
	      (setf (aref-sparse-matrix sout i j) (aref-sparse-matrix smat j i))))))

(defun sparse-+-2 (smat1 smat2)
  "Helper function of general sparse-+."
  (declare (type sparse-matrix smat1 smat2))
  (let ((m1-row (array-dimension (sparse-matrix-values smat1) 0))
	(m1-col (array-dimension (sparse-matrix-values smat1) 1))
	(m2-row (array-dimension (sparse-matrix-values smat2) 0))
	(m2-col (array-dimension (sparse-matrix-values smat2) 1)))
    (if (or (/= m1-col m2-col)
	    (/= m1-row m2-row))
	(error "size mismatch")
	(let ((sout (make-sparse-matrix-with-element m1-row m1-col)))
	  (loop for i from 0 to (1- m1-row) do
	       (loop for j from 0 to (1- m1-col) do
		    (if (boole boole-ior (aref (sparse-matrix-present-p smat1) i j) (aref (sparse-matrix-present-p smat2) i j))
			(setf (aref-sparse-matrix sout i j) (+ (aref-sparse-matrix smat1 i j)
							       (aref-sparse-matrix smat2 i j))))))
	  sout))))

(defun sparse---2 (smat1 smat2)
  "Helper function of general sparse--."
  (declare (type sparse-matrix smat1 smat2))
  (let ((m1-row (array-dimension (sparse-matrix-values smat1) 0))
	(m1-col (array-dimension (sparse-matrix-values smat1) 1))
	(m2-row (array-dimension (sparse-matrix-values smat2) 0))
	(m2-col (array-dimension (sparse-matrix-values smat2) 1)))
    (if (or (/= m1-row m2-row)
	    (/= m1-col m2-col))
	(error "size mismatch")
	(let ((sout (make-sparse-matrix-with-element m1-row m1-col)))
	  (loop for i from 0 to (1- m1-row) do
	       (loop for j from 0 to (1- m1-col) do
		    (if (boole boole-and (aref (sparse-matrix-present-p smat1) i j) (aref (sparse-matrix-present-p smat2) i j))
			(setf (aref-sparse-matrix sout i j) (- (aref-sparse-matrix smat1 i j)
							       (aref-sparse-matrix smat2 i j))))))
	  sout))))

(defun sparse-* (smat &rest more)
  "Returns product of sparse matrices, from left to right."
  (declare (type sparse-matrix smat))
  (reduce #'sparse-*-2 (cons smat more)))

(defun sparse-+ (smat &rest more)
  "Returns sum of sparse matrices, from left to right."
  (declare (type sparse-matrix smat))
  (reduce #'sparse-+-2 (cons smat more)))

(defun sparse-- (smat &rest more)
  "Returns difference of sparse matrices, from left to right."
  (declare (type sparse-matrix smat))
  (reduce #'sparse---2 (cons smat more)))

(defun sparse-matrix*vec (smat vec)
  "Returns the multiplication of a sparse-matrix and a vector."
  (declare (type sparse-matrix smat)
	   (type vector vec))
  )

(defun sparse-matrix-invert (smat)
  "Returns the inverse matrix of a sparse-matrix."
  (declare (type sparse-matrix smat))
  )

(defun sparse-matrix-determinant (smat)
  "Returns the deteminant of a sparse-matrix."
  (declare (type sparse-matrix smat))
  )

(defun sparse-matrix-rank (smat)
  "Returns the rank of a sparse-matrix."
  (declare (type sparse-matrix smat))
  )