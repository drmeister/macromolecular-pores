;;;
;;; Initialize everything
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Work with the structures
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(progn
  (setq *default-pathname-defaults*
	(pathname "/src/Dropbox/Development/molecular-designs/conrad-macrocycles/"))
  (defparameter *cd*
    (with-open-file
        (fin (probe-file "2h2h2macrocycle.cdxml") :direction :input)
      (chem:make-chem-draw fin)))
  (defparameter *agg* (chem:as-aggregate *cd*)))
(print "Done")

(progn
  (defparameter *stereocenters*
    (sort (cando:gather-stereocenters *agg*) #'string< :key #'chem:get-name))
  (cando:set-stereoisomer-func *stereocenters* (constantly :S) :show t)
  (let ((quat-matcher (core:make-cxx-object 'chem:chem-info)))
    (chem:compile-smarts quat-matcher "[C&H0&D4]")
    (chem:map-atoms nil (lambda (a) (when (chem:matches quat-matcher a)
				      (chem:set-configuration a :S)
				      (format t "Set atom ~a to :S~%" (chem:get-name a))))
		    *agg*)))

(defparameter *fix-atoms*
  (sort (select:atoms-with-property *agg* :fix) #'string<
	:key (lambda (a) (string (getf (chem:properties a) :fix)))))

(progn
  (defparameter *fix-points* (anchor:circle-points 40 (length *fix-atoms*)))
;;; Anchor the :fix atoms to *fixed-points*
  (anchor:on-points *fix-atoms* *fix-points*))


(defparameter *ff* (energy:setup-amber))
(cando:jostle *agg* 40)
(energy:minimize *agg* :force-field *ff* :restraints-on nil)

(save-mol2 *agg* "mc-test.mol2")




(cando:chimera *agg*)
(energy:minimize *agg* :force-field *ff* :restraints-on nil
		 :cg-tolerance 0.1
		 :max-tn-steps 100)

(cando:save-mol2 *agg* "tricycle.mol2")

