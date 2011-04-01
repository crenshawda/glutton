(ns glutton.enzyme-free
  (:use [clojure.contrib
         [def :only [defnk]]]
        (glutton [peptide-record]
                 (util :only [indexed]))))

(defn- all-sub-peptides [[[position first-aa] & other-aas] config]
  (let [{:keys [source digestion]} config]
    (if (seq? other-aas)
      (lazy-cat [(reductions (fn [p [pos aa]] (extend-with p aa))
                           (initiate-peptide first-aa (* position 3) source digestion)
                           other-aas)]
                (all-sub-peptides other-aas config))
      [[(initiate-peptide first-aa position source digestion)]])))

(defn- digest*
  [genome config]
  (let [{:keys [mass-target mass-tolerance]} config]
    (flatten
     (for [sub-peptides (all-sub-peptides genome config)]
       (take-while #(<= (:mass %) (+ mass-target mass-tolerance))
                   (drop-while #(< (:mass %) (- mass-target mass-tolerance))
                               sub-peptides))))))

(defnk digest
  "Return all sub-peptides of a translation of a nucleotide sequence whose masses fall within
  a target mass-tolerance.

  Optional Keyword Parameters:
  :mass-target    -> Mass target center, used to calculate the range of aceptable masses.  Defaults to 500 daltons.
  :mass-threshold -> Mass threshold is a ± value that is used to calculate the mass range.  Deaults to 100 daltons.
  :source         -> Organism that the genome came from.  Defaults to unknown."
  [genome :mass-target 500 :mass-tolerance 100 :source ""]
  (let [config {:mass-target mass-target
                :mass-tolerance mass-tolerance
                :source source
                :digestion "enzyme-free"}]
    (digest* (indexed genome) config)))
