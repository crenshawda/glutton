(ns glutton.peptide-record
  (:use [glutton [lexicon :as lex]]))

; This needs to be in genomic-utils as a public multi-method
(defn- aa-mass
  ([amino-acid]
     (aa-mass amino-acid :monoisotopic-mass))
  ([amino-acid mass-type]
     (mass-type (lex/*amino-acid-dictionary* amino-acid))))

(defprotocol Extend
  "Makes stuff longer"
  (extend-with [this item] "Adds item to a thing"))

(defrecord Peptide [sequence nucleotide-start mass breaks]
  Extend
  (extend-with [this aa]
    (-> this
        (update-in [:sequence] conj aa)
        (update-in [:mass] + (aa-mass aa)))))
