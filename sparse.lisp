;;;; Sparse matrix representation
;;;; Author: Augustus Huang
;;;; Date: June 7, 2015

(in-package :cl-cs-sparse)

;;; I decide to use CSR format instead... Should be space compact.
(eval-when (:compile-toplevel :load-toplevel :execute)
  (defstruct sparse-matrix
    (values (make-array 0) :type simple-array)
    (col-index (make-array 0 :element-type 'fixnum) :type simple-array)
    (row-ptr (make-array 0 :element-type 'fixnum) :type simple-array)
    (cols 0 :type fixnum))

  (defstruct sparse-vector
    (values (make-array 0) :type simple-array)
    (index (make-array 0 :element-type 'fixnum) :type simple-array)
    (len 0 :type fixnum)))

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

;;; Transform a vector into a sparse vector.
(defun make-sparse-vector-with-vector (vector)
  "Make function of a sparse vector with a generic vector."
  (declare (type vector vector))
  (let ((len (length vector))
	(values ())
	(index ()))
    (loop for i from 0 to (1- len) do
	 (let ((value (aref vector i)))
	   (if (/= 0 value)
	       (progn
		 (push value values)
		 (push i index)))))
    (make-sparse-vector :values (list-to-array (reverse values) 1)
			:index (list-to-array (reverse index) 1)
			:len len)))

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

(defun aref-sparse-vector (svec index)
  "Aref function of a sparse vector."
  (declare (type sparse-vector svec)
	   (type fixnum index))
  (loop for i from 0 to (1- (length (sparse-vector-index svec))) do
       (if (= index (aref (sparse-vector-index svec) i))
	   (return-from aref-sparse-vector (aref (sparse-vector-values svec) i))))
  0)

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

(defun (setf aref-sparse-vector) (value svec index)
  "(setf aref) function of a sparse vector."
  (declare (type sparse-vector svec)
	   (type fixnum index))
  (loop for i from 0 to (1- (length (sparse-vector-index svec))) do
       (let ((ind (aref (sparse-vector-index svec) i)))
	 (if (<= index ind)
	     ;; Just substitute it...
	     (if (= index ind)
		 (progn
		   (setf (aref (sparse-vector-values svec) index) value)
		   (return-from aref-sparse-vector))
		 ;; The first time we meet a index greater than our goal.
		 (let ((helper1 ())
		       (helper2 ())
		       (values-list (1d-array-to-list (sparse-vector-values svec)))
		       (index-list (1d-array-to-list (sparse-vector-index svec))))
		   (loop for j from 0 to (- i 2) do
			(push (pop values-list) helper1)
			(push (pop index-list) helper2))

		   (push value helper1)
		   (push index helper2)

		   (loop for j from 0 to (1- i) do
			(push (pop helper1) values-list)
			(push (pop helper2) index-list))
		   (setf (sparse-vector-values svec) (list-to-array values-list 1)
			 (sparse-vector-index svec) (list-to-array index-list 1))))))))

(defun transpose-sparse-matrix (smat)
  "Transpose function of a sparse matrix."
  (declare (type sparse-matrix smat))
  (let* ((row (1- (length (sparse-matrix-row-ptr smat))))
	 ;; Since we know transposition won't change slots...
	 (vout (make-array (length sparse-matrix-values smat) :initial-element 0))
	 (jout (make-array (length vout) :initial-element 0))
	 (temp1 0)
	 (temp2 0)
	 (iout (make-array (1+ row) :initial-element 0)))
    (loop for i from 0 to (1- row) do
	 (let ((row-start (aref (sparse-matrix-row-ptr smat) i))
	       (row-end (aref (sparse-matrix-row-ptr smat) (1+ i))))
	   (loop for j from row-start to (1- row-end) do
		(setf temp1 (incf (aref (sparse-matrix-col-index smat) j)))
		(incf (aref iout temp)))))
    (loop for i from 0 to (1- row) do
	 (incf (aref iout (1+ i)) (aref iout i)))
    (loop for i from 0 to (1- row) do
	 (let ((row-start (aref (sparse-matrix-row-ptr smat) i))
	       (row-end (aref (sparse-matrix-row-ptr smat) (1+ i))))
	   (loop for j from row-start to (1- row-end) do
		(setf temp1 (aref (sparse-matrix-col-index smat))
		      temp2 (aref iout temp1)
		      (aref vout temp2) (sparse-matrix-values j)
		      (aref jout temp2) i
		      (aref iout temp1) (incf temp2)))))
    (loop for i from (1- row) downto 0 do
	 (setf (aref iout (1+ i) (aref iout i))))
    (make-sparse-matrix :values vout
			:col-index jout
			:row-ptr iout
			:cols row)))

(defun sparse-m-*-2 (smat1 smat2)
  "Helper function of general sparse matrix multiplication, size mismatch won't be checked."
  (declare (type sparse-matrix smat1 smat2))
  (let* ((m1-row (1- (length (sparse-matrix-row-ptr smat1))))
	 (m1-col (sparse-matrix-cols smat1))
	 (m2-col (sparse-matrix-cols smat2))
	 ;; work array will have the same length with m1-col.
	 (work (make-array m1-col :initial-element 0))
	 ;; And iout will have the same length with m1-row.
	 (iout (make-array m1-row :initial-element 0))
	 ;; They will be cropped into simple arrays at last.
	 (jout (make-array (+ m1-row m2-col) :initial-element 0))
	 (vout (make-array (+ m1-row m2-col) :initial-element 0))
	 (len 0))
    (loop for i from 0 to (1- m1-row) do
	 (let ((row-start (aref (sparse-matrix-row-ptr smat1) i))
	       (row-end (aref (sparse-matrix-row-ptr smat1) (1+ i))))
	   (loop for j from row-start to (1- row-end) do
		(let* ((temp1 (aref (sparse-matrix-col-index smat1) j))
		       (temp2 (aref (sparse-matrix-values smat1) j))
		       (row2-start (aref (sparse-matrix-row-ptr smat2) temp2))
		       (row2-end (aref (sparse-matrix-row-ptr smat2) (1+ temp2))))
		  (loop for k from row2-start to (1- row2-end) do
		       (let* ((temp3 (aref (sparse-matrix-col-index smat2) k))
			      (temp4 (aref work temp3)))
			 (if (= 0 temp4)
			     (progn
			       (incf len)
			       (setf (aref work temp3) len)
			       (if (>= (1- len) (length jout))
				   (setf jout (adjust-array jout (* 2 (length jout)))
					 vout (adjust-array vout (* 2 (length vout)))
					 (aref jout (1- len)) temp3
					 (aref vout (1- len)) (* (aref (sparse-matrix-values smat2) k)
								 temp1))
				   (setf (aref jout (1- len)) temp3
					 (aref vout (1- len)) (* (aref (sparse-matrix-values smat2) k)
								 temp1))))
			     (if (>= temp4 (length vout))
				 (progn
				   (setf vout (adjust-array vout (* 2 (length vout))))
				   (incf (aref vout temp4) (* (aref (sparse-matrix-values smat2) k)
							      temp1)))
				 (incf (aref vout temp4) (* (aref (sparse-matrix-values smat2) k)
							    temp1))))))))
	   (loop for j from (aref iout i) to (1- len) do
		(setf (aref work (aref jout j)) 0))))
    (setf (aref iout (1+ i)) (1+ len))
    ;; Since there will be no 0s in vout,
    ;; we are safe to remove them totally...
    ;; And then jout should have the same length.
    (make-sparse-matrix :values (setf vout (ignore-trailing-zero vout))
			:col-index (adjust-array jout (length vout))
			:row-ptr iout
			:cols m2-col)))

(defun sparse-m-+-2 (smat1 smat2)
  "Helper function of general sparse matrix addition."
  (declare (type sparse-matrix smat1 smat2))
  (let* ((m1-row (1- (length (sparse-matrix-row-ptr smat1))))
	 (m1-col (sparse-matrix-cols smat1))
	 (m2-col (sparse-matrix-cols smat2))
	 (work (make-array m1-col :initial-element 0))
	 (iout (make-array m1-row :initial-element 0))
	 (jout (make-array (+ m1-row m2-col) :initial-element 0))
	 (vout (make-array (+ m1-row m2-col) :initial-element 0))
	 (len 0))
    (loop for i from 0 to (1- m1-row) do
	 (let ((row-start (aref (sparse-matrix-row-ptr smat1) i))
	       (row-end (aref (sparse-matrix-row-ptr smat1) (1+ i))))
	   (loop for j from row-start to (1- row-end) do
		(let ((col (aref (sparse-matrix-col-index smat1) j)))
		  (incf len)
		  (if (> len (length jout))
		      (setf jout (adjust-array jout (* 2 (length jout)))
			    vout (adjust-array vout (* 2 (length vout)))
			    (aref jout (1- len)) col
			    (aref vout (1- len)) (aref (sparse-matrix-values smat1) j))
		      (setf (aref jout (1- len)) col
			    (aref vout (1- len)) (aref (sparse-matrix-values smat1) j)))
		  (setf (aref work col) (1- len))))
	   (loop for j from row-start to (1- row-end) do
		(let* ((col (aref (sparse-matrix-col-index smat2) j))
		       (pos (aref work col)))
		  (if (= pos 0)
		      (progn
			(incf len)
			(if (> len (length jout))
			    (setf jout (adjust-array jout (* 2 (length jout)))
				  vout (adjust-array vout (* 2 (length vout)))
				  (aref jout (1- len)) col
				  (aref vout (1- len)) (aref (sparse-matrix-values smat2) j)
				  (aref work col) (1- len))))
		      (if (>= pos (length vout))
			  (progn
			    (setf vout (adjust-array vout (* 2 (length vout))))
			    (incf (aref vout pos) (aref (sparse-matrix-values smat2) j)))
			  (incf (aref vout pos) (aref (sparse-matrix-values smat2) j))
	   (loop for k from row-start to len do
		(setf (aref work (nth (- (length jout) k) jout)) 0))
	   (setf (aref iout (1+ i)) (1+ len))))
    (make-sparse-matrix :values (setf vout (ignore-trailing-zero vout))
			:col-index (adjust-array jout (length vout))
			:row-ptr iout
			:cols m2-col)))))))

(defun negative-sparse-matrix (smat)
  "Negating function of a sparse matrix."
  (declare (type sparse-matrix smat))
  (loop for i from 0 to (1- (length (sparse-matrix-values smat))) do
       (setf (aref (sparse-matrix-values smat) i) (- (aref (sparse-matrix-values smat) i)))))

(defun negative-sparse-vector (svec)
  "Negating function of a sparse vector."
  (declare (type sparse-vector svec))
  (loop for i from 0 to (1- (length (sparse-vector-values svec))) do
       (setf (aref (sparse-vector-values svec) i) (- (aref (sparse-vector-values svec) i)))))

(defun sparse-m---2 (smat1 smat2)
  "Helper function of general sparse matrix subtraction."
  (declare (type sparse-matrix smat1 smat2))
  (sparse-m-+-2 smat1 (negative-sparse-matrix smat2)))

(defun sparse-m-* (smat &rest more)
  "Product function of sparse matrices, from left to right."
  (declare (type sparse-matrix smat))
  (reduce #'sparse-m-*-2 (cons smat more)))

(defun sparse-m-+ (smat &rest more)
  "Addition function of sparse matrices, from left to right."
  (declare (type sparse-matrix smat))
  (reduce #'sparse-m-+-2 (cons smat more)))

(defun sparse-m-- (smat &rest more)
  "Subtraction function of sparse matrices, from left to right."
  (declare (type sparse-matrix smat))
  (reduce #'sparse-m---2 (cons smat more)))

(defun sparse-v-+-2 (svec1 svec2)
  "Helper function of general sparse vector addition."
  (declare (type sparse-vector svec1 svec2))
  (let ((len (sparse-vector-len svec1))
	(ilen (length (sparse-vector-index svec2)))
	(value-index (mapcar #'list
			     (1d-array-to-list (sparse-vector-values svec1))
			     (1d-array-to-list (sparse-vector-index svec2)))))
    (flet ((list-second (lst) (mapcar #'second lst))
	   (second-< (lst1 lst2)
	     (if (< (second lst1) (second lst2))
		 t
		 nil))
	   (find-nth (item lst)
	     (let ((i 1))
	       (dolist (obj lst)
		 (if (/= item obj)
		     (incf i)))
	       i)))
      ;; If we have found a slot are nonzero in both vectors,
      ;; add the value of svec2 to the value of svec1,
      ;; if there's a new slot in svec2, push it to the back.
      (loop for i from 0 to (1- ilen) do
	   (let ((item (aref (sparse-vector-index svec2) i))
		 (val (aref (sparse-vector-values svec2) i))
		 (nth 0))
	     ;; Here in order to evaluate T use n+1 instead of n...
	     (if (setf nth (find-nth item (list-second value-index)))
		 (incf (first (nth (1- nth) value-index)) val)
		 ;; Not found, push it to the back.
		 (push (list val item) value-index))))
      (setf value-index (sort value-index #'second-<)))
    (make-sparse-vector :values (list-to-array (mapcar #'first value-index) 1)
			:index (list-to-array (mapcar #'second value-index) 1)
			:len len)))

(defun sparse-v---2 (svec1 svec2)
  "Helper function of general sparse vector subtraction."
  (declare (type sparse-vector svec1 svec2))
  (sparse-v-+-2 svec1 (negative-sparse-vector svec2)))

(defun sparse-v-+ (svec &rest more)
  "Addition function of sparse vectors, from left to right."
  (declare (type sparse-vector svec))
  (reduce #'sparse-v-+-2 (cons svec more)))

(defun sparse-v-- (svec &rest more)
  "Subtraction function of sparse vectors, from left to right."
  (declare (type sparse-vector svec))
  (reduce #'sparse-v---2 (cons svec more)))

;;; All functions won't check the size mismatch. Use them carefully.
(defun sparse-matrix-*-vector (smat vec)
  "Multiplication routine of a sparse matrix and a vector."
  (declare (type sparse-matrix smat)
	   (type vector vec))
  (let* ((row (1- (array-dimension (sparse-matrix-row-ptr smat) 0)))
	 (len (length vec))
	 (out (make-array row :initial-element 0)))
    (loop for i from 0 to (1- row) do
	 (let ((row-start (aref (sparse-matrix-row-ptr smat) i))
	       (row-end (aref (sparse-matrix-row-ptr smat) (1+ i))))
	   (loop for j from row-start to (1- row-end) do
		(incf (aref out i) (* (aref (sparse-matrix-values smat) j)
				      (aref vec (aref (sparse-matrix-col-index smat) j)))))))
    out))

(defun sparse-matrix-*-sparse-vector (smat svec)
  "Mutiplication routine of a sparse matrix and a sparse vector."
  (declare (type sparse-matrix smat)
	   (type sparse-vector svec))
  ;; Since only the columns with index corresponding to nonzero slots in the
  ;; sparse vector are contributing, ignore other columns.
  (let ((slen (sparse-vector-len svec))
	(srow (1- (length (sparse-matrix-row-ptr smat))))
	(scol (sparse-matrix-cols smat))
	(out (make-array srow :initial-element 0)))
    (flet (())
      (loop for i from 0 to (1- srow) do
	   (let ((row-start (aref (sparse-matrix-row-ptr smat) i))
		 (row-end (aref (sparse-matrix-row-ptr smat) (1+ i))))
	     (loop for j from row-start to (1- row-end) do
		  (incf (aref out i) (* (aref (sparse-matrix-values smat) j)
					))))))))

(defun matrix-*-sparse-vector (mat svec)
  "Multiplication routine of a general matrix and a sparse vector."
  (declare (type matrix mat)
	   (type sparse-vector svec))
  (let ((slen (sparse-vector-len svec))
	(mlen (array-dimension mat 0))
	(out (make-array mlen :initial-element 0)))
    (loop for i from 0 to (1- mlen) do
	 (loop for j from 0 to (1- slen) do
	      (incf (aref out i) (* (aref mat i (aref (sparse-vector-index svec) j))
				    (aref (sparse-vector-values svec) j)))))
    out))

(defun sparse-matrix-*-const (smat const)
  "Multiplication routine of a sparse matrix and a constant."
  (declare (type sparse-matrix smat)
	   (type fixnum const))
  (let ((v-length (length (sparse-matrix-values smat))))
    (loop for i from 0 to (1- v-length) do
	 (setf (aref (sparse-matrix-values smat) i)
	       (* (aref (sparse-matrix-values smat) i) const)))))

(defun sparse-inner-product (svec1 svec2)
  "Inner product of two sparse vectors."
  (declare (type sparse-vector svec1 svec2))
  (let* ((ilen1 (length (sparse-vector-index svec1)))
	 (ilen2 (length (sparse-vector-index svec2)))
	 (vout (make-array ilen1 :initial-element 0))
	 (iout (make-array ilen1 :initial-element 0))
	 (iptr 0))
    ;; Loop through svec1, if there's an item of svec2 in slot i nonzero,
    ;; multiply them and add it to the new vector, or do nothing.
    (loop for i from 0 to (1- ilen1) do
	 (loop for j from 0 to (1- ilen2) do
	      (if (= (aref (sparse-vector-index svec2) j)
		     (aref (sparse-vector-index svec1) i))
		  (setf (aref vout iptr)
			(* (aref (sparse-vector-values svec1) i)
			   (aref (sparse-vector-values svec2) j))
			(aref iout iptr) (aref (sparse-vector-index svec1) i)
			(iptr (1+ iptr))))))
    (make-sparse-vector :values (ignore-trailing-zero vout)
			:index (ignore-trailing-zero iout)
			:len (sparse-vector-len svec1))))
