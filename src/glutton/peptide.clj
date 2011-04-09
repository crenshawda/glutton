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


(defrecord PeptideFromDNA [peptide sequence strand start stop mass internal-breaks])
(defrecord PeptideFromRNA [peptide sequence start stop mass internal-breaks])
(defrecord PeptideFromProtein [peptide sequence start stop mass internal-breaks])

(defprotocol PeptideCandidate
  (extend-peptide [this amino-acid])
  (finish-candidate [this]))

(defrecord DNAPeptideCandidate [peptide mass mass-type strand start stop breaks sequence-id sequence-length]
  PeptideCandidate
  (extend-peptide [this amino-acid]
                  (-> this
                      (update-in [:peptide] conj amino-acid)
                      (update-in [:mass] + (aa-mass amino-acid mass-type))))
  (finish-candidate [this]
                    (let [length-in-nt (* 3 (count peptide))]
                      (PeptideFromDNA. (apply str (map name peptide))
                                       sequence-id
                                       strand
                                       (if (= strand :+)
                                         start
                                         (- stop length-in-nt))
                                       (if (= strand :+)
                                         (+ start length-in-nt)
                                         stop)
                                       mass
                                       (if (zero? breaks)
                                         breaks
                                         (if (= (last peptide) :.)
                                           breaks
                                           (dec breaks)))))))

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



(defrecord RNAPeptideCandidate [peptide mass mass-type start stop breaks])
(defrecord ProteinPeptideCandidate [peptide mass mass-type start stop breaks])
