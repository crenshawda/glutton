; Clojure implementation of the Giddings Lab GFS software
; Authors: Dennis Crenshaw, Chris Maier, Brian Risk
; Date: December 27, 2010

(ns glutton.core
  [:require [glutton
             [glutton :as glt]
             [genomic-utils :as gu]
             [file-utils :as file]]])

;; to-amino-acids should be done after digestion as an output option?
(defn peptides-from-fasta
  "Given a FASTA file-str, this generates peptides from the genome sequence"
  [file-str]
  (pmap glt/digest
        (gu/to-amino-acids
         (:frames
          ;; This first is fragile-- it should be able to comprehend better than this!
          (first (gu/single-fasta file-str))))))

(defn peptide-count
  "Given a FASTA file-str, this counts peptides generated from the genome sequence"
  [file-str]
  (reduce + 0 (map count (peptides-from-fasta file-str))))
