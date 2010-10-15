(ns glutton.enzyme-free
  (:use [clojure.contrib.seq-utils :only [indexed]]
        [glutton [peptide-record]])
  (:import [glutton.peptide-record Peptide]))

(defn- all-sub-peptides [[[position first-aa] & other-aas]]
  (if (seq? other-aas)
    (lazy-cat [(reductions (fn [p [pos aa]] (extend-with p aa))
                           (create-peptide first-aa (* position 3)) other-aas)]
              (all-sub-peptides other-aas))
    [[(initiate-peptide first-aa position)]]))

(defn digest
  "Return all sub-peptides of a translation of a nucleotide sequence whose masses fall within
  a specified tolerance of a target mass."
  [genome mass-target mass-tolerance]
  (flatten
   (for [sub-peptides (all-sub-peptides (indexed genome))]
     (take-while #(<= (:mass %) (+ mass-target mass-tolerance))
                 (drop-while #(< (:mass %) (- mass-target mass-tolerance))
                             sub-peptides)))))
