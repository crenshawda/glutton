(ns glutton.enzyme-free
  (:use [clojure.contrib.seq-utils :only [indexed]]
        [glutton [peptide-record]])
  (:import [glutton.peptide-record Peptide]))

(defn create-peptide [aa position]
  "Constructor creates peptide with 0 breaks, as they aren't important in this algorithm"
  (Peptide. [aa] position (aa-mass aa) 0))

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
