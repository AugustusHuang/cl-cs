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

;;;; General compressive sensing packages

(defpackage :cl-cs
  (:nicknames :cs)
  (:use :cl :cl-cs-sparse :cl-cs-util)
  (:export
   :least-square
   :stagewise-omp
   :lasso))

(defpackage :cl-cs-sparse
  (:nicknames :cs-sparse)
  (:use :cl :cl-cs-util)
  (:export
   :make-sparse-matrix-with-matrix
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
  (:export
   :make-1-to-n-list
   :make-1-to-n-vector
   :dovec
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
