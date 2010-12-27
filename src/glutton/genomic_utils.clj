(ns glutton.genomic-utils
  [:require
   [glutton
    [lexicon :as lex]
    [file-utils :as file]]])

(defn normalize-genomic-data
  "Normalizes a genomic sequence by making a single, upper-cased string"
  [genomic-seq]
  (->> genomic-seq
      (apply str)
      (.toUpperCase)))

(defn parse-dna-string
  "This is an accessory method to help parse the necleotide sequence from a fasta file, returns a seq of codon keywords."
  [nucleotide-sequence]
  (for [codon (partition 3 nucleotide-sequence)]
     (keyword (apply str codon))))

(defn compliment-base-pairs
  "Compliments nucleotide base pairs (mostly for reverse frame readings)"
  [nucleotide-sequence]
    (replace lex/nucleotide-base-pair-dictionary nucleotide-sequence))

(defn to-amino-acids
  "Translates a nucleotide frame or a sequence of frame into amino-acids"
  [frame-seq]
  (if (seq? (first frame-seq))
    (pmap to-amino-acids frame-seq)
    (replace lex/codon-translation-matrix frame-seq)))

(defn aa-mass
  ([amino-acid]
     (aa-mass amino-acid :monoisotopic-mass))
  ([amino-acid mass-type]
     (mass-type (lex/amino-acid-dictionary amino-acid))))
