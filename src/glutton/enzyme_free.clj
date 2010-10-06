(ns glutton.enzyme-free
  (:use [clojure.contrib.seq-utils :only [indexed]])
  (:use [glutton [lexicon :as lex]]))

(defn- aa-mass
  ([amino-acid]
     (aa-mass amino-acid :monoisotopic-mass))
  ([amino-acid mass-type]
     (mass-type (lex/*amino-acid-dictionary* amino-acid))))

(defprotocol Extend
  "Makes stuff longer"
  (extend-with [this item] "Adds item to a thing"))

(defrecord Peptide [sequence nucleotide-start mass]
  Extend
  (extend-with [this aa]
    (-> this
        (update-in [:sequence] conj aa)
        (update-in [:mass] + (aa-mass aa)))))

(defn create-peptide [aa position]
  (Peptide. [aa] position (aa-mass aa)))

(defn- all-sub-peptides [[[position first-aa] & other-aas]]
  (if (seq? other-aas)
    (lazy-cat [(reductions (fn [p [pos aa]] (extend-with p aa))
                           (create-peptide first-aa (* position 3)) other-aas)]
              (all-sub-peptides other-aas))
    [[(create-peptide first-aa position)]]))

(defn digest
  "Return all sub-peptides of a translation of a nucleotide sequence whose masses fall within
  a specified tolerance of a target mass."
  [genome mass-target mass-tolerance]
  (flatten
   (for [sub-peptides (all-sub-peptides (indexed genome))]
     (take-while #(<= (:mass %) (+ mass-target mass-tolerance))
                 (drop-while #(< (:mass %) (- mass-target mass-tolerance))
                             sub-peptides)))))
