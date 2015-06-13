;;;; Sparse matrix representation
;;;; Author: Augustus Huang
;;;; Date: June 7, 2015

(in-package :cl-cs)

;;; I decide to use CSR format instead... Should be space compact.
(eval-when (:compile-toplevel :load-toplevel :execute)
  (defstruct sparse-matrix
    (values (make-array 0) :type simple-array)
    (col-index (make-array 0 :element-type 'fixnum) :type simple-array)
    (row-ptr (make-array 0 :element-type 'fixnum) :type simple-array)
    (cols 0 :type fixnum)))

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
			:row-ptr (list-to-array (reverse ptr-list) 1)
			:cols col)))

(defun aref-sparse-matrix (smat row column)
  "Aref function of a sparse matrix."
  (declare (type sparse-matrix smat)
	   (type fixnum row column))
  (let ((row-start (aref (sparse-matrix-row-ptr smat) row))
	(row-end (aref (sparse-matrix-row-ptr smat) (1+ row))))
    (loop for i from row-start to (1- row-end) do
	 (if (= column (aref (sparse-matrix-col-index smat) i))
	     (return-from aref-sparse-matrix (aref (sparse-matrix-values smat) i))))
    0))

(defun (setf aref-sparse-matrix) (value smat row column)
  "(setf aref) function of a sparse matrix."
  (declare (type sparse-matrix smat)
	   (type fixnum row column))
  (let ((row-start (aref (sparse-matrix-row-ptr smat) row))
	(row-end (aref (sparse-matrix-row-ptr smat) (1+ row))))
    (loop for i from row-start to (1- row-end) do
	 (let ((col (aref (sparse-matrix-col-index smat) i)))
	   (if (<= column col)
	       (if (= column col)
		   ;; Directly use the new value instead.
		   (progn
		     (setf (aref (sparse-matrix-values smat) i) value)
		     (return-from aref-sparse-matrix))
		   ;; The first time we meet a column index larger than
		   ;; where we want to insert a new value,
		   ;; adjust these three arrays.
		   (let ((helper1 ())
			 (helper2 ())
			 (values-list (1d-array-to-list (sparse-matrix-values smat)))
			 (index-list (1d-array-to-list (sparse-matrix-col-index smat))))
		     ;; Firstly we push the items before the new item into
		     ;; helper list, so are the indices.
		     (loop for j from 0 to (- i 2) do
			  (push (pop values-list) helper1)
			  (push (pop index-list) helper2))
		     
		     ;; And then push the item we wanna change from 0.
		     (push value helper1)
		     (push column helper2)

		     ;; And push back, turn them into arrays and change the
		     ;; row information.
		     (loop for j from 0 to (1- i) do
			  (push (pop helper1) values-list)
			  (push (pop helper2) index-list))
		     (setf (sparse-matrix-values smat) (list-to-array values-list 1)
			   (sparse-matrix-col-index smat) (list-to-array index-list 1)))))))
    (loop for j from (1+ row) to (1- (length (sparse-matrix-row-ptr smat))) do
	 (incf (aref (sparse-matrix-row-ptr smat) j)))))

(defun transpose-sparse-matrix (smat)
  "Transpose function of a sparse matrix."
  (declare (type sparse-matrix smat))
  (let (())))

(defun sparse-*-2 (smat1 smat2)
  "Helper function of general sparse-*, size mismatch won't be checked."
  (declare (type sparse-matrix smat1 smat2))
  (let* ((m1-row (1- (array-dimension (sparse-matrix-row-ptr smat1) 0)))
	 (m1-col (sparse-matrix-cols smat1))
	 (m2-col (sparse-matrix-cols smat2))
	 ;; work array will have the same length with m1-col.
	 (work (make-array m1-col :initial-element 0))
	 ;; And iout will have the same length with m1-row.
	 (iout (make-array m1-row :initial-element 0))
	 ;; They will be casted into arrays at last.
	 (jout ())
	 (vout ())
	 (len 0))
    (loop for i from 0 to (1- m1-row) do
	 (let ((row-start (aref (sparse-matrix-row-ptr smat1) i))
	       (row-end (aref (sparse-matrix-row-ptr smat2) (1+ i))))
	   (loop for j from row-start to (1- row-end) do
		(let (())))))))

(defun sparse-+-2 (smat1 smat2)
  "Helper function of general sparse-+."
  (declare (type sparse-matrix smat1 smat2))
  (let* ((m1-row (1- (array-dimension (sparse-matrix-row-ptr smat1) 0)))
	 (m1-col (sparse-matrix-cols smat1))
	 (m2-col (sparse-matrix-cols smat2))
	 (work (make-array m1-col :initial-element 0))
	 (iout (make-array m1-row :initial-element 0))
	 (jout ())
	 (vout ())
	 (len 0))
    (loop for i from 0 to (1- m1-row) do
	 (let ((row-start (aref (sparse-matrix-row-ptr smat1) i))
	       (row-end (aref (sparse-matrix-row-ptr smat1) (1+ i))))
	   (loop for j from row-start to (1- row-end) do
		(let ((col (aref (sparse-matrix-col-index smat1) j)))
		  (incf len)
		  (push col jout)
		  (push (aref (sparse-matrix-values smat1) j) vout)
		  (setf (aref work j) len)))
	   (loop for j from row-start to (1- row-end) do
		(let* ((col (aref (sparse-matrix-col-index smat2) j))
		       (pos (aref work col)))
		  (if (= pos 0)
		      (progn
			(incf len)
			(push col jout)
			(push (aref (sparse-matrix-values smat2) j) vout)
			(setf (aref work j) len))
		      (incf (nth (- (length vout) pos) vout) (aref (sparse-matrix-values smat2) j)))))
	   (loop for k from row-start to len do
		(setf (aref work (nth (- (length jout) k) jout)) 0))
	   (setf (aref iout (1+ i)) (1+ len))))
    (make-sparse-matrix :values (list-to-array (reverse vout) 1)
			:col-index (list-to-array (reverse jout) 1)
			:row-ptr iout
			:cols m2-col)))

(defun negative-sparse (smat)
  "Returns the negative matrix of a sparse matrix."
  (declare (type sparse-matrix smat))
  (loop for i from 0 to (1- (length (sparse-matrix-values smat))) do
       (setf (aref (sparse-matrix-values smat) i) (- (aref (sparse-matrix-values smat) i)))))

(defun sparse---2 (smat1 smat2)
  "Helper function of general sparse--."
  (declare (type sparse-matrix smat1 smat2))
  (sparse-+-2 smat1 (negative-sparse smat2)))

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
  "Returns the multiplication of a sparse matrix and a vector."
  (declare (type sparse-matrix smat)
	   (type vector vec))
  )

(defun sparse-matrix-rank (smat)
  "Returns the rank of a sparse matrix."
  (declare (type sparse-matrix smat))
  )