(ns glutton.peptide-record
  (:use (glutton [lexicon :as lex]
                 [genomic-utils :as gu])))

(defprotocol Extend
  "Makes stuff longer"
  (extend-with [this item] "Adds item to a thing"))

(defrecord Peptide [sequence nucleotide-start mass breaks]
  Extend
  (extend-with [this aa]
    (-> this
        (update-in [:sequence] conj aa)
        (update-in [:mass] + (gu/aa-mass aa)))))

(defn initiate-peptide [aa position]
  (Peptide. [aa] position (gu/aa-mass aa) 0))
