(ql:quickload :cl-permutation)
(ql:quickload :alexandria)

(defpackage ternary-afcom
  (:use #:cl)
  (:nicknames #:afcom)
  (:local-nicknames (#:alex #:alexandria)))

(in-package :ternary-afcom)

;;; Compute affinely commuting triples

(defun commute-affinely-p (x y z)
  "Does the permutations X, Y & Z commute affinely?"
  (let* ((x-1 (perm:perm-inverse x))
         (a (perm:perm-compose x-1 y))
         (b (perm:perm-compose x-1 z)))
    (perm:perm= (perm:perm-compose a b) (perm:perm-compose b a))))

(defun afcom (G)
  "Return the list of affinely commuting triples in G."
  (uiop:while-collecting (relate)
    (alex:map-combinations
     (lambda (triple)
       (when (apply #'commute-affinely-p triple)
         (relate triple)))
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

(defun quaternion (n)
  "Generalized quaternion group of order 4n.
This is the subgroup of the unit quaternions generated by a j and a
2n-th root of unity."
  ;; Let w be a 2n-th root of unity, then the number k corresponds to
  ;; w^k for k in [1,2n], and to w^k j for k in [2n+1,4n]
  (perm:group-from-cycles
   (list
    (list (apply #'perm:make-cycle (alex:iota (* 2 n) :start 1))
          (apply #'perm:make-cycle (alex:iota (* 2 n) :start (1+ (* 2 n)))))
    (loop for k from 1 to n
          collect (perm:make-cycle k  (+ (* 3 n) k)  (+ n k) (+ (* 2 n) k))))
   (* 4 n)))

;; TODO: semidihedral, what else?

;;; Isomorphism of small hypergraphs

(defun lex< (xs ys &key (eq #'=) (lt #'<))
  "Lexicographical ordering relative to EQ and LT."
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
  (if (< n k)
      '(())
      (collect-hypergraphs
       (lambda (save-hypergraph)
         (let ((attach (complete-uniform-hypergraph (1- n) (1- k))))
           (mapc (lambda (h)
                   (map-subsets
                    (lambda (s)
                      (funcall
                       save-hypergraph
                       (append h (mapcar (lambda (e) (append e (list n))) s))))
                    attach))
                 (uniform-hypergraphs (1- n) k))))
       :canonicalize t)))

;;; Induced subhypergraphs

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

;;; Searching for hypergraphs inside AfCom(G)

(defun missing-uniform-hypergraphs (n hypergraphs &rest groups)
  "List N-vertex 3-uniform HYPERGRAPHS not found in afcom of any of the GROUPS.
The HYPERGRAPHS parameter can be the symbol :all meaning all 3-uniform
hypergraphs on N vertices which do not have an induced copy of
((1 2 3) (1 2 4) (1 3 4)), since that cannot occur in an (afcom G)."
  (reduce
   (lambda (hypergraphs group)
     (set-difference hypergraphs
                     (all-induced-subhypergraphs (normalize (afcom group)) n)
                     :test #'equal))
   groups
   :initial-value
   (if (eq hypergraphs :all)
       (remove-if
        (lambda (hypergraph)
          (has-induced-subhypergraph-p hypergraph '((1 2 3) (1 2 4) (1 3 4))))
        (uniform-hypergraphs n 3))
       hypergraphs)))

(defun prepare (hypergraph)
  "Parse HYPERGRAPH into arrays of affinely commuting conditions.
Returns the number of vertices of the HYPERGRAPH and two arrays, in
the first the k-th entry lists pairs (i j) such that (i j k) should be
an affinely commuting triple. The second array lists the remaining
pairs, those for which the triple should not be affinely commutative."
  (let* ((v (number-of-vertices hypergraph))
         (yes (make-array (1+ v) :initial-element nil))
         (no (make-array (1+ v) :initial-element nil)))
    (dolist (h hypergraph)
      (push (butlast h) (aref yes (car (last h)))))
    (dotimes (u (1+ v))
      (setf (aref no u)
            (loop for i from 1 below u
                  append (loop for j from (1+ i) below u
                               for e = (list i j)
                               unless (member e (aref yes u) :test #'equal)
                                 collect e))))
    (values v yes no)))

(defun realize (hypergraph n)
  "Attempt to find HYPERGRAPH as an induced subhypergraph of AfCom(S_N)."
  (multiple-value-bind (v yes no) (prepare hypergraph)
    (when (< v 3) (return-from realize t))
    (let ((r (make-array (1+ v) :initial-element nil)))
      (labels
          ((try (u)
             (if (> u v)
                 r
                 (perm:doperms (p n)
                   (and
                    (loop for i from 1 below u never (perm:perm= p (aref r i)))
                    (loop for (i j) in (aref yes u)
                          always (commute-affinely-p (aref r i) (aref r j) p))
                    (loop for (i j) in (aref no u)
                          never (commute-affinely-p (aref r i) (aref r j) p))
                    (setf (aref r u) p)
                    (when (try (1+ u)) (return-from try r)))))))
        (setf (aref r 1) (perm:perm-identity n))
        (try 2)))))

(defun size-5-not-in-sym-6 ()
  "List of 3-uniform HYPERGRAPHS with 5 vertices not found in AfCom(S_6)."
  (remove-if
   (lambda (h)
     (or (has-induced-subhypergraph-p h '((1 2 3) (1 2 4) (1 3 4)))
         ;; (equal h '((1 2 3) (1 2 4) (1 3 5) (2 4 5) (3 4 5)))
         (realize h 6)))
   (uniform-hypergraphs 5 3)))

;;; Some hypergraphs

(defparameter snub-disphenoid
  '((1 2 3) (1 2 6) (1 3 4) (1 4 5) (1 5 6)
    (2 3 8) (2 7 8) (3 4 8) (4 5 8) (5 7 8)
    (5 6 7)))
