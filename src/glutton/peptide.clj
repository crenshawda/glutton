(ns glutton.peptide
  (:use (glutton (mass :only [water-mass aa-mass]))))

;; TODO: Finished peptides have an `internal-breaks` field, but this is kind of meaningless without
;;       an indication as to what the protease is

(defprotocol PeptideCandidate
  "Behaviors common to all PeptideCandidates, whether they come from DNA, RNA, or protein sequences."
  (extend-peptide [this amino-acid]
    "Appends `amino-acid` to the C-terminus of the peptide candidate")
  (finish-candidate [this]
    "Finalizes peptide coordinates and generates a finished peptide"))

(def nucleotide-candidate-impl
  {:extend-peptide (fn [this amino-acid]
                     (-> this
                      (update-in [:peptide] conj amino-acid)
                      (update-in [:mass] + (aa-mass amino-acid (:mass-type this)))))})

;; ## DNA Peptide Candidates

;;TODO sequence-length is not actually used in the record... just its constructor
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

;; TODO initial-breaks should be a boolean, not a number
;; TODO Extract the initial amino-acid plus water into a function

(defn dna-peptide-candidate
  "Constructor function for a DNA-based peptide candidate.

  `initial-amino-acid` will be the first (N-terminal) amino acid of the peptide.
  `strand` is either :+ or :-, and influences what the final start and stop positions will be.
  `strand-start` is the zero-based position from the 3' end of the strand that this position starts at.
  `mass-type` is a keyword indicating which masses to calculate the peptide's weight (e.g., :monoisotopic-mass)
  `sequence-id` is an identifier for the sequence being digested (e.g., \"chr1\")
  `overall-sequence-length` is the length of the entire original sequence being digested.  This is required
      in order for accurate positions on the reverse strand to be computed.
  `initial-breaks` is an optional parameter indicating whether the first amino acid is also a cleavage site
      for the protease in question.  If it is missing, the first amino acid is NOT considered a break."
  [initial-amino-acid strand strand-start mass-type sequence-id overall-sequence-length & [initial-breaks]]
  (let [strand-dependent-start (if (= strand :+) strand-start)
        strand-dependent-stop (if (= strand :-) (- overall-sequence-length strand-start))]
    (DNAPeptideCandidate. [initial-amino-acid]
                          (+ (aa-mass initial-amino-acid mass-type)
                             (water-mass mass-type))
                          mass-type
                          strand
                          strand-dependent-start
                          strand-dependent-stop
                          (or initial-breaks 0)
                          sequence-id
                          overall-sequence-length)))

;; ## RNA Peptide Candidates

(defrecord RNAPeptideCandidate [peptide mass mass-type start stop sequence-id breaks])
(defrecord PeptideFromRNA [peptide sequence start stop mass internal-breaks])

(defn rna-peptide-candidate
  "Constructor function for an RNA-based peptide candidate.

  It differs from `dna-peptide-candidate`, in that there is no strand to keep track of.

  `initial-amino-acid` will be the first (N-terminal) amino acid of the peptide.
  `start` is the zero-based position from the 3' end of the RNA at which this peptide starts.
  `mass-type`  is a keyword indicating which masses to calculate the peptide's weight (e.g., :monoisotopic-mass)
  `sequence-id`  is an identifier for the sequence being digested (e.g., \"my RNA molecule\")
  `initial-breaks` is an optional parameter indicating whether the first amino acid is also a cleavage site
      for the protease in question.  If it is missing, the first amino acid is NOT considered a break."
  [initial-amino-acid start mass-type sequence-id & [initial-breaks]]
  (RNAPeptideCandidate. [initial-amino-acid]
                        (+ (aa-mass initial-amino-acid mass-type)
                           (water-mass mass-type))
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

;; ## Protein Peptide Candidates

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
  "Constructor function for a protein-based peptide candidate.

  `initial-amino-acid` will be the first (N-terminal) amino acid of the peptide.
  `start` is the zero-based position from the N-terminus of the protein at which this peptide starts.
  `mass-type`  is a keyword indicating which masses to calculate the peptide's weight (e.g., :monoisotopic-mass)
  `sequence-id`  is an identifier for the sequence being digested (e.g., \"GSK3-Beta\")
  `initial-breaks` is an optional parameter indicating whether the first amino acid is also a cleavage site
      for the protease in question.  If it is missing, the first amino acid is NOT considered a break."
  [initial-amino-acid start mass-type sequence-id  & [initial-breaks]]
  (ProteinPeptideCandidate. [initial-amino-acid]
                            (+ (aa-mass initial-amino-acid mass-type)
                               (water-mass mass-type))
                            mass-type
                            start
                            nil
                            sequence-id
                            (or initial-breaks 0)))
