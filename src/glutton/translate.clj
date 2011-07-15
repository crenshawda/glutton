(ns glutton.translate
  "DNA sequence translation functions.  These have been tested both manually and
  against the [Bioinformatics.org Sequence Manipulation Tools Suite][1]

  [1]: http://www.bioinformatics.org/sms/index.html"
  (:require (clojure (string :as string))))

(def standard-genetic-code
  "Based on the [DNA codon table](http://en.wikipedia.org/wiki/DNA_codon_table).
  Keys are DNA codon strings, and values are single-letter amino acid codes. Stop
  codons are represented by \"*\"."
  {"GCT" "A" "GCC" "A" "GCA" "A" "GCG" "A" ; Alanine
   "TGT" "C" "TGC" "C"                     ; Cysteine
   "GAT" "D" "GAC" "D"                     ; Aspartic Acid
   "GAA" "E" "GAG" "E"                     ; Aspartate
   "TTT" "F" "TTC" "F"                     ; Phenylalanine
   "GGT" "G" "GGC" "G" "GGA" "G" "GGG" "G" ; Glycine
   "CAT" "H" "CAC" "H"                     ; Histidine
   "ATT" "I" "ATC" "I" "ATA" "I"           ; Isoleucine
   "AAA" "K" "AAG" "K"                     ; Lysine
   "TTA" "L" "TTG" "L" "CTT" "L" "CTC" "L" "CTA" "L" "CTG" "L" ; Leucine
   "ATG" "M"                               ; Methionine
   "AAT" "N" "AAC" "N"                     ; Asparagine
   "CCT" "P" "CCC" "P" "CCA" "P" "CCG" "P" ; Proline
   "CAA" "Q" "CAG" "Q"                     ; Glutamine
   "CGT" "R" "CGC" "R" "CGA" "R" "CGG" "R" "AGA" "R" "AGG" "R" ; Arginine
   "TCT" "S" "TCC" "S" "TCA" "S" "TCG" "S" "AGT" "S" "AGC" "S" ; Serine
   "ACT" "T" "ACC" "T" "ACA" "T" "ACG" "T" ; Threonine
   "GTT" "V" "GTC" "V" "GTA" "V" "GTG" "V" ; Valine
   "TGG" "W"                               ; Tryptophan
   "TAT" "Y" "TAC" "Y"                     ; Tyrosine
   "TAA" "*" "TGA" "*" "TAG" "*"           ; STOP
   })

(defn to-rna-code
  "Given a genetic code map (codons to amino acids), convert all codons to use
  RNA bases (i.e., tranform all thymines to uracils).  If you pass in an RNA
  code map, it will remain unchanged."
  [genetic-code]
  (into {}
        (for [[codon acid] genetic-code]
          [(apply str (replace {\T \U} codon)) acid])))

(defn translate
  "Translate a string of nucleotides to a string of amino acids.  Accepts optional
  keyword arguments for the genetic `:code` (a map of codon strings to amino acid
  single-letter codes) and the `:frame` in which to do the translation.  If no
  keyword arguments are given, defaults to translation by the Standard Genetic Code
  in frame 0.

  Currently, ambiguity codes are not handled intelligently with respect to degeneracy,
  such that \"GCN\" is translated to \"A\" (Alanine); any ambiguity code in a codon
  means the codon is translated to \"X\", the standard code for \"any amino acid\"."
  [^String nucleotides & {:keys [code frame]
                          :or {code standard-genetic-code
                               frame 0}}]
  {:pre [(contains? #{0 1 2} frame)]}
  (let [sb (new StringBuffer)
        max-pos (- (.length nucleotides) 3)
        upcased-nucleotides (string/upper-case nucleotides)]
    (loop [position frame]
      (if (<= position max-pos)
        (let [codon-end (+ position 3)]
          (.append sb (get code
                           (.substring upcased-nucleotides position codon-end)
                           "X"))
          (recur codon-end))))
    (.toString sb)))


(defn complement-nucleotides
  "Converts each nucleotide character to its complement.  All letters in the nucletoide
  string are upper-cased prior to complementation, and all resulting letters are upper-cased."
  [^String nucleotides]
  (apply str (replace {\A \T
                       \T \A
                       \C \G
                       \G \C} (string/upper-case nucleotides))))

(defn reverse-complement
  "Reverses and complements a nucleotide string.  All resulting letters are upper-cased."
  [^String nucleotides]
  (complement-nucleotides (string/reverse nucleotides)))

(defn translate-fasta
  "Given a FASTA record map, create a new FASTA record map with the translated peptide
  sequence, and \" Frame `X`\" appended to the header (where `X` is 0, 1, or 2)."
  [{:keys [header sequence] :as fasta-record}
   & {:keys [code frame]
      :or {code standard-genetic-code
           frame 0}
      :as params}]
  {:header (str header " Frame " frame)
   :sequence (apply translate sequence (flatten (seq params)))})

(defn reverse-complement-fasta
  "Given a FASTA record, return a new FASTA record map with the sequence reverse-complemented,
who  and \" Reverse Complement\" appended to the header."
  [{:keys [header sequence] :as fasta}]
  {:header (str header " Reverse Complement")
   :sequence (reverse-complement sequence)})

(defn six-frame-translation
  "For a given FASTA record, returns a sequence of 6 FASTA records, one for each of the six
  frames of translation (frames 0, 1, and 2, forward and reverse complement)"
  [fasta]
  (concat
   (map translate-fasta
        (repeat fasta)
        (repeat :frame)
        (range 3))
   (map translate-fasta
        (repeat (reverse-complement-fasta fasta))
        (repeat :frame)
        (range 3))))
