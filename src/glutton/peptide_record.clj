(ns glutton.peptide-record
  (:use (glutton [genomic-utils :as gu])))

(defprotocol Extend
  "Makes stuff longer"
  (extend-with [this item] "Adds item to a thing"))

(defrecord Peptide [sequence nucleotide-start mass breaks source digestion]
  Extend
  (extend-with [this aa]
    (-> this
        (update-in [:sequence] conj aa)
        (update-in [:mass] + (gu/aa-mass aa)))))

(defn initiate-peptide
  ([aa position]
     (initiate-peptide [aa] position "" ""))
  ([aa position source digestion]
     (Peptide. [aa] position (gu/aa-mass aa) 0 source digestion)))
