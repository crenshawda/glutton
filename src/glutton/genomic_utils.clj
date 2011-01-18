(ns glutton.genomic-utils
  (:use (clojure [string :only [upper-case join]])
        (glutton [lexicon :only [amino-acid-dictionary
                                 codon-translation-matrix
                                 nucleotide-base-pair-dictionary
                                 water-constants]]
                 [file-utils :only [file->records]])))

(defn- parse-dna-string
  "This is an accessory method to help parse the nucleotide sequence from a fasta file, returns a seq of codon keywords."
  [nucleotide-sequence]
  (for [codon (partition 3 nucleotide-sequence)]
    (keyword (upper-case (join codon)))))

(defn- compliment-base-pairs
  "Compliments nucleotide base pairs (mostly for reverse frame readings)"
  [nucleotide-sequence]
    (replace nucleotide-base-pair-dictionary nucleotide-sequence))

(defn- to-amino-acids
  "Translates a nucleotide frame or sequence of frames into amino-acids"
  [frame-seq]
  (if (seq? (first frame-seq))
    (pmap to-amino-acids frame-seq)
    (replace codon-translation-matrix frame-seq)))

(defn aa-mass
  ([amino-acid]
     (aa-mass amino-acid :monoisotopic-mass))
  ([amino-acid mass-type]
     (mass-type (amino-acid-dictionary amino-acid))))


(defn water-mass
  ([]
     (water-mass :monoisotopic-mass))
  ([mass-type]
     (mass-type water-constants)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; FASTA Utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- parse-fasta
  [header sequence]
  (let [sequence (join sequence)
        ;; reverse isn't lazy! is this killing me?
        reverse-compliment-sequence (compliment-base-pairs (reverse sequence))
        ->aas #(for [n (range 3)]
                 (to-amino-acids (parse-dna-string (drop n %))))]
    {:header header
     ;; Reading Frame Indexes 0-2 Forward, 3-5 Reverse-compliment
     :frames (lazy-cat (->aas sequence)
                       (->aas reverse-compliment-sequence))}))

(defn single-fasta
  [file-str]
  (for [fasta-seq (file->records file-str)]
    (let [[header sequence] fasta-seq]
      (parse-fasta header sequence))))
