(ns glutton.genomic-utils
  [:require
   [clojure.string :as string]
   [glutton
    [lexicon :as lex]
    [file-utils :as file]]])

(defn- parse-dna-string
  "This is an accessory method to help parse the nucleotide sequence from a fasta file, returns a seq of codon keywords."
  [nucleotide-sequence]
  (for [codon (partition 3 nucleotide-sequence)]
    (keyword (string/upper-case (string/join codon)))))

(defn- compliment-base-pairs
  "Compliments nucleotide base pairs (mostly for reverse frame readings)"
  [nucleotide-sequence]
    (replace lex/nucleotide-base-pair-dictionary nucleotide-sequence))

(defn- to-amino-acids
  "Translates a nucleotide frame or sequence of frames into amino-acids"
  [frame-seq]
  (if (seq? (first frame-seq))
    (pmap to-amino-acids frame-seq)
    (replace lex/codon-translation-matrix frame-seq)))

(defn aa-mass
  ([amino-acid]
     (aa-mass amino-acid :monoisotopic-mass))
  ([amino-acid mass-type]
     (mass-type (lex/amino-acid-dictionary amino-acid))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; FASTA Utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- parse-fasta
  [header sequence]
  (let [sequence (string/join sequence)
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
  (for [fasta-seq (file/file->records file-str)]
    (let [[header sequence] fasta-seq]
      (parse-fasta header sequence))))
