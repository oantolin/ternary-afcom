(ql:quickload :cl-permutation)
(ql:quickload :alexandria)

;;; Compute affinely commuting triples


(defun afcom (G)
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
    (cons n (alexandria:iota (1- n) :start 1))
    (append (alexandria:iota (1- n) :start (1- n) :step -1) (list n)))))

;; TODO: semidihedral, generalized quaternion, C_m x| C_n, etc.

;;; Isomorphism of small hypergraphs

(defun lex< (xs ys &key (eq #'=) (lt #'<))
  (loop for xt on xs and yt on ys
        while (and xt yt (funcall eq (car xt) (car yt)))
        finally (return (if (null xt)
                            (not (null yt))
                            (and yt (funcall lt (car xt) (car yt)))))))

(defun map-vertices (fn hypergraph)
  (sort (mapcar (lambda (edge) (sort (mapcar fn edge) #'<)) hypergraph) #'lex<))

(defun normalize (hypergraph &key (test #'perm:perm=))
  (let ((vertices (reduce (lambda (x y) (union x y :test test))
                          hypergraph)))
    (map-vertices
     (lambda (v) (1+ (position v vertices :test test))) hypergraph)))

(defun number-of-vertices (hypergraph)
  (reduce #'max (reduce #'append hypergraph) :initial-value 0))

(defun canonicalize (hypergraph)
  (when hypergraph
    (let ((canon hypergraph))
      (perm:doperms (p (number-of-vertices hypergraph))
        (let ((h (map-vertices (lambda (v) (perm:perm-eval p v)) hypergraph)))
          (when (lex< h canon :eq #'equal :lt #'lex<) (setf canon h))))
      canon)))

(defun map-subsets (fn set)
  (if (null set)
      (funcall fn nil)
      (let ((x (car set)))
        (map-subsets (lambda (s) (funcall fn (cons x s)) (funcall fn s))
                     (cdr set)))))

(defun collect-hypergraphs (process &key (canonicalize t))
  (let ((hypergraphs (make-hash-table :test #'equal)))
    (funcall
     process
     (lambda (h)
       (setf (gethash (if canonicalize (canonicalize h) h) hypergraphs) t)))
    (alexandria:hash-table-keys hypergraphs)))

(defun complete-uniform-hypergraph (n k)
  (uiop:while-collecting (save-edge)
    (alexandria:map-combinations
     #'save-edge
     (alexandria:iota n :start 1)
     :length k :copy t)))

(defun uniform-hypergraphs (n k)
  (collect-hypergraphs
   (lambda (save-hypergraph)
     (map-subsets save-hypergraph (complete-uniform-hypergraph n k)))
   :canonicalize t))

(defun induced-subhypergraph (hypergraph subset)
  (canonicalize
   (remove-if-not (lambda (edge) (subsetp edge subset)) hypergraph)))

(defun map-induced-subhypergraphs (fn hypergraph size)
  (alexandria:map-combinations
   (lambda (subset)
     (funcall fn (induced-subhypergraph hypergraph subset)))
   (alexandria:iota (number-of-vertices hypergraph) :start 1)
   :length size))

(defun all-induced-hypergraphs (hypergraph size)
  (collect-hypergraphs
   (lambda (save-hypergraph)
     (map-induced-subhypergraphs save-hypergraph hypergraph size))
   :canonicalize nil))

(defun has-induced-subhypergraphp (hypergraph subhypergraph)
  (map-induced-subhypergraphs
   (lambda (h)
     (when (equal h subhypergraph) (return-from has-induced-subhypergraphp t)))
   hypergraph (number-of-vertices subhypergraph))
  nil)
