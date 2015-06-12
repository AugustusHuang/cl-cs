;;;; Sparse matrix representation
;;;; Author: Augustus Huang
;;;; Date: June 7, 2015

(in-package :cl-cs)

;;; I decide to use CSR format instead... Should be space compact.
(eval-when (:compile-toplevel :load-toplevel :execute)
  (defstruct sparse-matrix
    (values (make-array 0) :type simple-array)
    (col-index (make-array 0 :element-type 'fixnum) :type simple-array)
    (row-ptr (make-array 0 :element-type 'fixnum) :type simple-array))

  (defstruct adjustable-sparse-matrix
    (values (make-array 0 :adjustable t))
    (col-index (make-array 0 :element-type 'fixnum :adjustable t))
    (row-ptr (make-array 0 :element-type 'fixnum :adjustable t))))

(defun list-dimensions (list depth)
  "Counts the dimension of a list."
  (loop repeat depth
       collect (length list)
       do (setf list (car list))))

(defun list-to-array (list depth)
  "Makes an array from a given list."
  (make-array (list-dimensions list depth) :initial-contents list))

(defun list-to-array-adjustable (list depth)
  "Makes an adjustable array from a given list."
  (make-array (list-dimensions list depth) :initial-contents list :adjustable t))

;;; Transform a matrix into a sparse matrix.
(defun make-sparse-matrix-with-matrix (matrix)
  "Make function of a sparse matrix with a generic matrix."
  (declare (type matrix matrix))
  (let ((row (array-dimension matrix 0))
	(col (array-dimension matrix 1))
	(v-length 0)
	(values-list ())
	(index-list ())
	(ptr-list ()))
    (loop for i from 0 to (1- row) do
	 (setf v-length (length values-list))
	 (loop for j from 0 to (1- col) do
	      (let ((val (aref matrix i j)))
		(if (/= 0 val)
		    (progn
		      (push val values-list)
		      (push j index-list)))))
	 (push v-length ptr-list))
    (push (length index-list) ptr-list)
    (make-sparse-matrix :values (list-to-array (reverse values-list) 1)
			:col-index (list-to-array (reverse index-list) 1)
			:row-ptr (list-to-array (reverse ptr-list) 1))))

(defun make-adjustable-sparse-matrix-with-matrix (matrix)
  "Make function of an adjustable sparse matrix with a generic matrix."
  (declare (type matrix matrix))
  (let ((row (array-dimension matrix 0))
	(col (array-dimension matrix 1))
	(v-length 0)
	(values-list ())
	(index-list ())
	(ptr-list ()))
    (loop for i from 0 to (1- row) do
	 (setf v-length (length values-list))
	 (loop for j from 0 to (1- col) do
	      (let ((val (aref matrix i j)))
		(if (/= 0 val)
		    (progn
		      (push val values-list)
		      (push j index-list)))))
	 (push v-length ptr-list))
    (push (length index-list) ptr-list)
    (make-adjustable-sparse-matrix :values (list-to-array-adjustable (reverse values-list) 1)
			:col-index (list-to-array-adjustable (reverse index-list) 1)
			:row-ptr (list-to-array-adjustable (reverse ptr-list) 1))))

(defun aref-sparse-matrix (smat row column)
  "Aref function of a sparse matrix."
  (declare (type sparse-matrix smat)
	   (type fixnum row column))
  (let ((row-start (aref (sparse-matrix-row-ptr smat) row))
	(row-end (aref (sparse-matrix-row-ptr smat) (1+ row))))
    (loop for i from row-start to (1- row-end) do
	 (if (= column (aref (sparse-matrix-col-index smat) i))
	     (return (aref (sparse-matrix-values smat) i))))
    0))

(defun aref-adjustable-sparse-matrix (smat row column)
  "Aref function of an adjustable sparse matrix."
  (declare (type adjustable-sparse-matrix smat)
	   (type fixnum row column))
  (let ((row-start (aref (adjustable-sparse-matrix-row-ptr smat) row))
	(row-end (aref (adjustable-sparse-matrix-row-ptr smat) (1+ row))))
    (loop for i from row-start to (1- row-end) do
	 (if (= column (aref (adjustable-sparse-matrix-col-index smat) i))
	     (return (aref (adjustable-sparse-matrix-values smat) i))))
    0))

(defun (setf aref-adjustable-sparse-matrix) (smat row column value)
  "(setf aref) function of an adjustable sparse matrix."
  (declare (type adjustable-sparse-matrix smat)
	   (type fixnum row column)))

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