(ns glutton.peptide
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

(defn peptide
  "Simple constructor function for a Peptide record."
  [sequence nucleotide-start mass breaks source digestion]
  (Peptide. sequence nucleotide-start mass breaks source digestion))

(defn initiate-peptide
  ([aa position source digestion]
     (initiate-peptide aa position "n/a" source digestion))
  ([aa position breaks source digestion]
     (Peptide. [aa] position
               (+ (aa-mass aa)
                  (water-mass))
               breaks source digestion)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;






(defprotocol PeptideCandidate
  (extend-peptide [this amino-acid])
  (finish-candidate [this]))

(def nucleotide-candidate-impl
  {:extend-peptide (fn [this amino-acid]
                     (-> this
                      (update-in [:peptide] conj amino-acid)
                      (update-in [:mass] + (aa-mass amino-acid (:mass-type this)))))})


(defrecord DNAPeptideCandidate [peptide mass mass-type strand start stop breaks sequence-id sequence-length])
(defrecord PeptideFromDNA [peptide sequence strand start stop mass internal-breaks])

(extend DNAPeptideCandidate
  PeptideCandidate
  (merge nucleotide-candidate-impl
         {:finish-candidate (fn [this]
                              (let [peptide (:peptide this)
                                    strand (:strand this)
                                    start (:start this)
                                    stop (:stop this)
                                    breaks (:breaks this)
                                    length-in-nt (* 3 (count peptide))]
                                (PeptideFromDNA. (apply str (map name peptide))
                                                 (:sequence-id this)
                                                 strand
                                                 (if (= strand :+)
                                                   start
                                                   (- stop length-in-nt))
                                                 (if (= strand :+)
                                                   (+ start length-in-nt)
                                                   stop)
                                                 (:mass this)
                                                 (if (zero? breaks)
                                                   breaks
                                                   (if (= (last peptide) :.)
                                                     breaks
                                                     (dec breaks))))))}))

(defn dna-peptide-candidate
  "Constructor function for a DNA-based peptide candidate"
  [initial-amino-acid strand strand-start mass-type sequence-id overall-sequence-length & [initial-breaks]]
  (let [strand-dependent-start (if (= strand :+) strand-start)
        strand-dependent-stop (if (= strand :-) (- overall-sequence-length strand-start))]
    (DNAPeptideCandidate. [initial-amino-acid]
                          (+ (aa-mass initial-amino-acid mass-type)
                             (water-mass))
                          mass-type
                          strand
                          strand-dependent-start
                          strand-dependent-stop
                          (or initial-breaks 0)
                          sequence-id
                          overall-sequence-length)))



(defrecord RNAPeptideCandidate [peptide mass mass-type start stop sequence-id breaks])
(defrecord PeptideFromRNA [peptide sequence start stop mass internal-breaks])

(defn rna-peptide-candidate [initial-amino-acid start mass-type sequence-id & [initial-breaks]]
  (RNAPeptideCandidate. [initial-amino-acid]
                        (+ (aa-mass initial-amino-acid mass-type)
                           (water-mass))
                        mass-type
                        start
                        nil
                        sequence-id
                        (or initial-breaks 0)))

(extend RNAPeptideCandidate
  PeptideCandidate
  (merge nucleotide-candidate-impl
         {:finish-candidate (fn [this]
                              (let [peptide (:peptide this)
                                    start (:start this)
                                    breaks (:breaks this)]
                                (PeptideFromRNA. (apply str (map name peptide))
                                                 (:sequence-id this)
                                                 start
                                                 (+ start (* 3 (count peptide)))
                                                 (:mass this)
                                                 (if (zero? breaks)
                                                   breaks
                                                   (if (= (last peptide) :.)
                                                     breaks
                                                     (dec breaks))))))}))

(defrecord ProteinPeptideCandidate [peptide mass mass-type start stop sequence breaks])
(defrecord PeptideFromProtein [peptide sequence start stop mass internal-breaks])

(extend ProteinPeptideCandidate
  PeptideCandidate
  (merge nucleotide-candidate-impl
         {:finish-candidate (fn [this]
                              (let [peptide (:peptide this)
                                    start (:start this)
                                    breaks (:breaks this)]
                                (PeptideFromProtein. (apply str (map name peptide))
                                                     (:sequence this)
                                                     start
                                                     (+ start (count peptide))
                                                     (:mass this)
                                                     (if (zero? breaks)
                                                       breaks
                                                       (if (= (last peptide) :.)
                                                         breaks
                                                         (dec breaks))))))}))

(defn protein-peptide-candidate
  "Constructor function for a protein-based peptide candidate"
  [initial-amino-acid start mass-type sequence-id  & [initial-breaks]]
  (ProteinPeptideCandidate. [initial-amino-acid]
                            (+ (aa-mass initial-amino-acid mass-type)
                               (water-mass))
                            mass-type
                            start
                            nil
                            sequence-id
                            (or initial-breaks 0)))
