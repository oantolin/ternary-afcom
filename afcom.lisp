(ql:quickload :cl-permutation)
(ql:quickload :alexandria)

;;; Compute affinely commuting triples

(defun afcom (G)
  (uiop:while-collecting (relate)
    (alexandria:map-combinations
     (lambda (triple)
       (destructuring-bind ((i . x) (j . y) (k . z)) triple
         (let* ((x-1 (perm:perm-inverse x))
                (a (perm:perm-compose x-1 y))
                (b (perm:perm-compose x-1 z)))
           (when (perm:perm= (perm:perm-compose a b) (perm:perm-compose b a))
             (relate (list i j k))))))
     (mapcar (let ((i 0)) (lambda (x) (cons (incf i) x)))
             (uiop:while-collecting (collect)
               (perm:do-group-elements (x G)
                 (collect x))))
     :length 3)))

;;; Some constructions of groups

(defun dihedral (n)
  "The dihedral group of order 2n."
  (perm:group-from
   (list
    (cons n (alexandria:iota (1- n) :start 1))
    (append (alexandria:iota (1- n) :start (1- n) :step -1) (list n)))))

;; TODO: semidihedral, generalized quaternion, C_m x| C_n, etc.

;;; Isomorphism of small 3-uniform hypergraphs

;; TODO!
