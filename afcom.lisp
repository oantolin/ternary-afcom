(ql:quickload :cl-permutation)
(ql:quickload :alexandria)

(defpackage ternary-afcom
  (:use #:cl)
  (:nicknames #:afcom)
  (:local-nicknames (#:alex #:alexandria)))

(in-package :ternary-afcom)

;;; Compute affinely commuting triples

(defun afcom (G)
  "Return the list of affinely commuting triples in G."
  (uiop:while-collecting (relate)
    (alexandria:map-combinations
     (lambda (triple)
       (destructuring-bind (x y z) triple
         (let* ((x-1 (perm:perm-inverse x))
                (a (perm:perm-compose x-1 y))
                (b (perm:perm-compose x-1 z)))
           (when (perm:perm= (perm:perm-compose a b) (perm:perm-compose b a))
             (relate triple)))))
     (uiop:while-collecting (collect)
       (perm:do-group-elements (x G)
         (collect x)))
     :length 3)))

;;; Some constructions of groups

(defun dihedral (n)
  "The dihedral group of order 2n."
  (perm:group-from
   (list
    (cons n (alex:iota (1- n) :start 1))
    (append (alex:iota (1- n) :start (1- n) :step -1) (list n)))))

(defun semidirect-product-cyclics (a n)
  "The group Z/n x| <a>, where a is an invertible element of Z/n."
  (assert (= (gcd a n) 1))
  (perm:group-from
   (list
    (cons n (alex:iota (1- n) :start 1))
    (loop for k below n collect (1+ (mod (* a k) n))))))

(defun symmetric (n)
  "Symmetric group on n letters."
  (perm:group-from-cycles
   (list (list (perm:make-cycle 1 2))
         (list (apply #'perm:make-cycle (alex:iota n :start 1))))
   n))

(defun alternating (n)
  "Alternating group on n letters."
  (perm:group-from-cycles
   (loop for k from 3 to n collect (list (perm:make-cycle 1 2 k)))
   n))

;; TODO: semidihedral, generalized quaternion

;;; Isomorphism of small hypergraphs

(defun lex< (xs ys &key (eq #'=) (lt #'<))
  "Lexicographical ordering realtive to EQ and LT."
  (loop for xt on xs and yt on ys
        while (and xt yt (funcall eq (car xt) (car yt)))
        finally (return (if (null xt)
                            (not (null yt))
                            (and yt (funcall lt (car xt) (car yt)))))))

(defun map-vertices (fn hypergraph)
  "Apply FN to each vertex of HYPERGRAPH."
  (sort (mapcar (lambda (edge) (sort (mapcar fn edge) #'<)) hypergraph) #'lex<))

(defun normalize (hypergraph &key (test #'equal))
  "Return an isomorphic copy of HYPERGRAPH with vertices {1, 2, ..., n}."
  (let ((vertices (reduce (lambda (x y) (union x y :test test))
                          hypergraph :initial-value nil)))
    (map-vertices
     (lambda (v) (1+ (position v vertices :test test))) hypergraph)))

(defun number-of-vertices (hypergraph)
  "Return number of vertices of HYPERGRAPH.
This assumes the vertices are {1, 2, ..., n} for some n."
  (reduce #'max (reduce #'append hypergraph) :initial-value 0))

(defun canonicalize (hypergraph)
  "Return a canonical representative of the isomorphism class of HYPERGRAPH.
This assumes the vertices are {1, 2, ..., n} for some n. You can
compare two hypergraphs for isomorphism by comparing their canonical
representatives for equality."
  (when hypergraph
    (let ((canon hypergraph))
      (perm:doperms (p (number-of-vertices hypergraph))
        (let ((h (map-vertices (lambda (v) (perm:perm-eval p v)) hypergraph)))
          (when (lex< h canon :eq #'equal :lt #'lex<) (setf canon h))))
      canon)))

(defun map-subsets (fn set)
  "Apply FN to all subsets of SET."
  (if (null set)
      (funcall fn nil)
      (let ((x (car set)))
        (map-subsets (lambda (s) (funcall fn (cons x s)) (funcall fn s))
                     (cdr set)))))

(defun collect-hypergraphs (process &key (canonicalize t))
  "Return a deduplicated collection of hypergraphs produced by PROCESS.
PROCESS should be a function of one argument and will be called with a
function that should be called on each hypergraph to be collected. If
CANONICALIZE is non-nil, each hypergraph will be canonicalized before
collection."
  (let ((hypergraphs (make-hash-table :test #'equal)))
    (funcall
     process
     (lambda (h)
       (setf (gethash (if canonicalize (canonicalize h) h) hypergraphs) t)))
    (alex:hash-table-keys hypergraphs)))

(defun complete-uniform-hypergraph (n k)
  "Return the complete K-uniform hypergraph on N vertices."
  (uiop:while-collecting (save-edge)
    (alex:map-combinations
     #'save-edge
     (alex:iota n :start 1)
     :length k :copy t)))

(defun uniform-hypergraphs (n k)
  "Return a list of K-uniform hypergraphs on N vertices."
  (collect-hypergraphs
   (lambda (save-hypergraph)
     (map-subsets save-hypergraph (complete-uniform-hypergraph n k)))
   :canonicalize t))

(defun induced-subhypergraph (hypergraph subset)
  "Return the induced subhypergraph of HYPERGRAPH on a SUBSET of vertices.
The induced subhypergraph is canonicalized."
  (canonicalize
   (normalize
    (remove-if-not (lambda (edge) (subsetp edge subset)) hypergraph)
    :test #'eql)))

(defun map-induced-subhypergraphs (fn hypergraph size)
  "Apply FN to each induced subhypergraph of HYPERGRAPH of given SIZE."
  (let ((n (number-of-vertices hypergraph)))
    (when (>= n size)
     (alex:map-combinations
     (lambda (subset)
       (funcall fn (induced-subhypergraph hypergraph subset)))
     (alex:iota n :start 1)
     :length size))))

(defun all-induced-subhypergraphs (hypergraph size)
  "Return list of all induced subhypergraphs of HYPERGRAPH of given SIZE."
  (collect-hypergraphs
   (lambda (save-hypergraph)
     (map-induced-subhypergraphs save-hypergraph hypergraph size))
   :canonicalize nil))

(defun has-induced-subhypergraph-p (hypergraph subhypergraph)
  "Does the HYPERGRAPH have the SUBHYPERGRAPH as an induced subhypergraph?"
  (map-induced-subhypergraphs
   (lambda (h)
     (when (equal h subhypergraph) (return-from has-induced-subhypergraph-p t)))
   hypergraph (number-of-vertices subhypergraph))
  nil)
