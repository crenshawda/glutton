; Clojure implementation of the Giddings Lab GFS software
; Authors: Dennis Crenshaw, Chris Maier, Brian Risk
; Date: January 19, 2011

(ns glutton.core
  (:use (clojure.contrib [duck-streams :only [write-lines]])
        (glutton [glutton :only [digest]]
                 ;; [enzyme-free :only [digest]]
                 [genomic-utils :only [single-fasta aa-mass water-mass]])))

(defn peptides-from-fasta
  "Given a FASTA file-str and an optional organism source (eg. ecoli), this generates peptides from the genome sequence(s)"
  ([file-str]
     (peptides-from-fasta file-str "unknown"))
  ([file-str source]
     (flatten ;; Bucket of peptides, no need to keep the seq/frame structure
      (for [seq (single-fasta file-str)
            frame (:frames seq)]
        (digest frame :source source)))))
