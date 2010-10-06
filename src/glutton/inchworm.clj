(ns glutton.inchworm
  (:use [clojure.contrib.seq-utils :only [indexed]])
  (:use [glutton [lexicon :as lex]]))

;; (def mass-target 500)
;; (def mass-tolerance 100)

(defn- in-range? [mass target tolerance]
  (and (>= mass (- target tolerance))
       (<= mass (+ target tolerance))))

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
    (lazy-cat [(reductions (fn [p [pos aa]] (extend-with p aa ))
                           (create-peptide first-aa position) other-aas)]
              (all-sub-peptides other-aas))
    [[(create-peptide first-aa position)]]))

(defn- filtered-by [genome mass-target mass-tolerance]
  (flatten
   (remove nil?
           (for [sub-peptides (all-sub-peptides (indexed genome))]
             (let [partitions (partition-by #(in-range? (:mass %) mass-target mass-tolerance)
                                            sub-peptides)]
               (if (in-range? (:mass (first (first partitions))) mass-target mass-tolerance)
                 (first partitions)
                 (second partitions)))))))
