; Clojure implementation of the Giddings Lab GFS software
; Authors: Dennis Crenshaw, Chris Maier, Brian Risk
; Date: September 28, 2010

(ns glutton.core
    (:require (glutton [glutton :as glt]
                       [genomic-utils :as gu]
                       [file-utils :as file])))

(defn peptides-from-fasta-file
  [file-string]
  (pmap glt/digest (pmap gu/frame->amino-acids (:frames (file/single-fasta file-string)))))

(defn peptide-count
  "Given a file string, takes a fasta file and count generated peptides from the genome sequence"
  [file-string]
  (reduce + 0 (pmap count (peptides-from-fasta-file file-string))))
