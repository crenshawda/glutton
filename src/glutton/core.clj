; Clojure implementation of the Giddings Lab GFS software
; Authors: Dennis Crenshaw, Chris Maier, Brian Risk
; Date: January 27, 2011

(ns glutton.core
  (:use (glutton [genomic-utils :only [fasta->clj
                                       nucleotides->frames]]
                 [glutton :only [digest]])))

(defn ->peptides
  "Given a FASTA file-str and an optional organism source (eg. ecoli), this generates peptides from the genome sequence(s)"
  ([file-str]
     (->peptides file-str "unknown"))
  ([file-str source]
     (flatten ;; Bucket of peptides, no need to keep the seq/frame structure
      (for [frames (map nucleotides->frames (fasta->clj file-str))]
        (map #(digest % :source source) frames)))))
