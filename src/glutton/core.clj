; Clojure implementation of the Giddings Lab GFS software
; Authors: Dennis Crenshaw, Chris Maier, Brian Risk
; Date: December 27, 2010

(ns glutton.core
  (:use (clojure.contrib [duck-streams :only [write-lines]])
        (glutton [glutton :only [digest]]
                 [genomic-utils :only [single-fasta aa-mass water-mass]])))

(defn peptides-from-fasta
  "Given a FASTA file-str, this generates peptides from the genome sequence(s)"
  [file-str]
  (flatten ;; Bucket of peptides, no need to keep the seq/frame structure
   (for [seq (single-fasta file-str)
         frame (:frames seq)]
     (digest frame))))
