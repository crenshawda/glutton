(ns glutton.genomic-utils
  (:require (glutton [lexicon :as lex])
            (clojure.contrib [seq-utils :as su])))

(defn parse-dna-string
  "This is an accessory method to help parse the necleotide sequence from a fasta file, returns a seq of codon keywords."
  [nucleotide-sequence]
  (vec
   (for [codon (partition 3 nucleotide-sequence)]
     (keyword (apply str codon))))) ; Maybe I don't need vec here?
; This was the previous implementation, it seemed a good bit faster based on completely inadequate testing.
(comment
  (loop [codons []
	 sequence nucleotide-sequence]
    (if (empty? sequence)
      codons
      (recur
       (conj codons (keyword (apply str (take 3 sequence))))
       (drop 3 sequence)))))

(defn compliment-base-pairs
  "Compliments base pairs (ostensibly for reverse frame readings)"
  [nucleotide-sequence]
    (replace lex/*nucleotide-base-pair-dictionary* nucleotide-sequence))

(defn normalize-genomic-data
  "Normalizes a genomic data to help ensure parser compatability"
  [genomic-data-seq]
  (let [g-data (apply str genomic-data-seq)]
    (-> g-data
        (.toUpperCase))))

(defn frame->amino-acids
  "Pass this a reading frame of codon keywords and it will translate it into amino acids"
  [frame]
  (replace lex/*codon-translation-matrix* frame))

(defn translate-frames
  "Give this the list of the frames you want translated and it will process them"
  [frames]
  (pmap frame->amino-acids frames))

; I think this needs to be moved/refactored into file-utils
(defn parse-fasta
  "Give this a FASTA header and sequence and it'll parse it out into the important bits"
  [header sequence]
  (let [original-sequence (vec (normalize-genomic-data sequence))
        reverse-compliment-sequence (compliment-base-pairs (rseq original-sequence))]
    {:header header
     :frames [(parse-dna-string original-sequence)  ; Reading Frame Indexes 0-2 Forward, 3-5 Reverse-compliment
	      (parse-dna-string (drop 1 original-sequence))
	      (parse-dna-string (drop 2 original-sequence))
	      (parse-dna-string reverse-compliment-sequence)
	      (parse-dna-string (drop 1 reverse-compliment-sequence))
	      (parse-dna-string (drop 2 reverse-compliment-sequence))]}))

(defn aa-mass
  ([amino-acid]
     (aa-mass amino-acid :monoisotopic-mass))
  ([amino-acid mass-type]
     (mass-type (lex/*amino-acid-dictionary* amino-acid))))
