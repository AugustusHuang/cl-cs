;;;; The MIT License (MIT)

;;;; Copyright (c) 2015 Huang Xuxing

;;;; Permission is hereby granted, free of charge, to any person obtaining
;;;; a copy of this software and associated documentation files
;;;; (the "Software"), to deal in the Software without restriction,
;;;; including without limitation the rights to use, copy, modify, merge,
;;;; publish, distribute, sublicense, and/or sell copies of the Software,
;;;; and to permit persons to whom the Software is furnished to do so,
;;;; subject to the following conditions:

;;;; The above copyright notice and this permission notice shall be included
;;;; in all copies or substantial portions of the Software.

;;;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
;;;; IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;;;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
;;;; THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;;;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
;;;; ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
;;;; OTHER DEALINGS IN THE SOFTWARE.

;;;; Algorithms based on compressive sensing.

(in-package :cl-cs)

(defun sort-index (seq predicate)
  "Sort function with sorted index stored and returned."
  (declare (type sequence seq)
	   (type function predicate))
  (typecase seq
    (list
     (let ((seq-index (mapcar #'list seq (make-1-to-n-list (length seq)))))
       (flet ((first-predicate (lst1 lst2)
		(if (funcall predicate (first lst1) (first lst2))
		    t
		    nil)))
	 (mapcar #'second (sort seq-index #'first-predicate)))))
    (vector
     (let* ((seq-list (:texonomy-util::1d-array-to-list seq))
	    (seq-index (mapcar #'list seq-list (make-1-to-n-list (length seq)))))
       (flet ((first-predicate (lst1 lst2)
		(if (funcall predicate (first lst1) (first lst2))
		    t
		    nil)))
	 (mapcar #'second (sort seq-index #'first-predicate)))))))

(defun erf (vec)
  "Error function, iterating through the vector."
  (declare (type vector vec))
  (let* ((len (length vec))
	 (a1 0.254829592)
	 (a2 -0.284496736)
	 (a3 1.421413741)
	 (a4 -1.453152027)
	 (a5 1.061405429)
	 (p 0.3275911)
	 (out (make-array len :initial-element 0)))
    (loop for i from 0 to (1- len) do
	 (let* ((sign (if (>= (aref vec i) 0) 1 -1))
		(abs-value (abs (aref vec i)))
		(temp (/ 1.0 (+ 1.0 (* p abs-value)))))
	   (setf (aref out i) (* sign
				 (- 1.0 (* temp
					   (exp (* abs-value (- abs-value)))
					   (+ a1 (* temp (+ a2 (* temp (+ a3 (* temp (+ a4 (* temp a5))))))))))))))
    out))

(defun erfc (vec)
  "Compliment of error function, iterate over the vector."
  (declare (type vector vec))
  (let* ((len (length vec))
	 (out (make-array len :initial-element 0))
	 (erf-vec (erf vec)))
    (loop for i from 0 to (1- len) do
	 (setf (aref out i) (1- (aref erf-vec i))))
    out))
(defun fdrthresh (vec param)
  "Get the fdr threshold of VEC with PARAM."
  (declare (type vector vec)
	   (type real param))
  (let* ((abs-vec (vector-abs vec))
	 (len (length vec))
	 (sorted-vec (sort vec #'<))
	 (sort-index (sort-index vec #'<))
	 (pobs (erfc (./ sorted-vec (sqrt 2)))))
    (let* ((n (make-1-to-n-vector len))
	   (pnull (./ n len))
	   (maximum 0))
      (loop for i from 0 to (1- (length pobs)) do
	 ;; pobs has the same length as sorted-vec,
	 ;; and pnull has the same length as vec.
	   (let ((good (<= (aref probs (- (1- n) i)) (* param (aref pnull i)))))
	     (if (and good (> (aref n i) maximum))
		 (setf maximum (aref n i)))))
      (if (/= maximum 0)
	  ;; We've found some GOOD non-nil.
	  (aref abs-vec (nth (- (1+ len) maximum) sort-index))
	  ;; All GOODs are nil, return trivial value.
	  (+ 0.01 (vector-max abs-vec))))))
(defun hardthresh (vec param)
  "Apply the hard threshold PARAM to VEC."
  (declare (type vector vec)
	   (type real param))
  (let* ((len (length vec))
	 (out (make-array len :initial-element 0)))
    (loop for i from 0 to (1- len) do
	 (let ((ith (aref vec i)))
	   (if (> ith param)
	       (setf (aref out i) ith))))
    out))
(defun softthresh (vec param)
  "Apply the soft threshold PARAM to VEC."
  (declare (type vector vec)
	   (type real param))
  (let* ((len (length vec))
	 (out (make-array len :initial-element 0)))
    (loop for i from 0 to (1- len) do
	 (let* ((sign (if (>= (aref vec i) 0) 1 -1))
		(abs-value (abs (aref vec i)))
		(temp (- abs-value param)))
	   (set (aref out i) (* sign (/ (+ temp (abs temp)) 2)))))
    out))

(defun least-square (mat vec &key (stages -1))
  "Iterative Least Square method for solving linear equations."
  )

(defun stagewise-omp (matrix vec &key (thresh #'fdrthresh) (param 0.5) (iter 10) (err 1e-5))
  "Stagewise Orthogonal Matching Pursuit algorithm, get an approximating solution to L_1 minimization problem."
  (declare (type matrix matrix)
	   ;; VEC will be a vector, it's not necessary be sparse.
	   (type vector vector))
  (assert (= (array-dimension matrix 0)
	     (length vec))
	  (matrix vec thresh param iter err)
	  "Size mismatch, matrix of size ~D-by-~D but vector of length ~D."
	  (array-dimension matrix 0)
	  (array-dimension matrix 1)
	  (length vec))
  (let* ((len (length vec))
	 (col (array-dimension matrix 1))
	 ;; Initially the residual will be VEC.
	 (residual vec)
	 (vnorm (norm vec))
	 (i-full (make-1-to-n-vector col))
	 (i-now ())
	 (active ())
	 (j-active ())
	 (x-i (make-array col :initial-element 0))
	 ;; Now the output sparse vector is zero-sparse vector.
	 (out (make-sparse-vector :len col)))
    ;; Iterate ITER times and output the result in sparse vector form.
    (loop for i from 0 to (1- iter) do
	 (let* ((corr (./ (matrix-*-vector (.* (transpose matrix) (sqrt n))
					   residual)
			  (norm residual)))
		(thr (funcall thresh corr param)))
	   (setf i-now (1d-array-to-list (hardthresh (vector-abs corr) thr))
		 ;; UNION apply on lists, change them...
		 j-active (union active i-now))
	   (if (= (length j-active) (length active))
	       ;; Maybe we shall use some more gentle way?
	       (go done))
	   (setf active j-active
		 x-i (qr-solve (mask-matrix matrix active) vec)
		 residual (m- vec (matrix-*-vector (mask-matrix matrix active)
						   x-i)))
	   (if (<= (norm residual) (* err vnorm))
	       (go done))))
    ;; No way to assign a list to a vector... FIXME
    (setf (sparse-vector-values out) x-i
	  (sparse-vector-index out) active)
    out))

(defun lasso (smat vec &key (stages 0) (err 1e-5) (lambda-stop 0) (residual-stop 0))
  "Lasso algorithm, in Least Angle Regression by B. Efron, T. Hastie
I. Johnstone and R. Tibshirani."
  (declare (type sparse-matrix smat)
	   (type vector vec))
  )
