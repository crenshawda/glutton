(ns glutton.peptide-record
  (:use (glutton (genomic-utils :only [aa-mass water-mass]))))

(defprotocol Extend
  "Makes stuff longer"
  (extend-with [this item] "Adds item to a thing"))



(defrecord Peptide [sequence nucleotide-start mass breaks source digestion]
  Extend
  (extend-with [this aa]
    (-> this
        (update-in [:sequence] conj aa)
        (update-in [:mass] + (aa-mass aa)))))

(defn initiate-peptide
  ([aa position]
     (initiate-peptide [aa] position "" ""))
  ([aa position breaks source digestion]
     (Peptide. [aa] position
               (+ (aa-mass aa)
                  (water-mass))
               breaks source digestion)))
