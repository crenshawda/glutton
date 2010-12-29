; Clojure implementation of the Giddings Lab GFS software
; Authors: Dennis Crenshaw, Chris Maier, Brian Risk
; Date: December 27, 2010

(ns glutton.core
  [:require [glutton
             [glutton :as glt]
             [genomic-utils :as gu]
             [file-utils :as file]]])

(defn peptides-from-fasta
  "Given a FASTA file-str, this generates peptides from the genome sequence(s)"
  [file-str]
  (flatten ;; Bucket of peptides, no need to keep the seq/frame structure
   (for [seq (gu/single-fasta file-str)]
     ;; Alternatively: (pmap glt/digest (:frames seq))
     (for [frame (:frames seq)]
       (glt/digest frame)))))

(defn peptide-count
  "Given a FASTA file-str, this counts peptides generated from the genome sequence"
  [file-str]
  (reduce + 0 (map count (peptides-from-fasta file-str))))
