;;;; General compressive sensing packages
;;;; Author: Augustus Huang
;;;; Date: June 9, 2015

(defpackage :cl-cs
  (:nicknames :cs)
  (:use :cl :cl-cs-util)
  (:export))

(defpackage :cl-cs-util
  (:nicknames :cs-util)
  (:use :cl)
  (:export :dovec
	   :doseq
	   :random-array
	   :with-gensyms))