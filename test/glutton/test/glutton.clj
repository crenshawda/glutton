(ns glutton.test.glutton
  (:use glutton.glutton :reload)
  (:use clojure.test)
  (:use (clojure.java (io :only [file resource]))
        (glutton (genomic-utils :only [fasta->clj]))))

(deftest test-sequence-type
  (let [dna "ACGATCGATGGCTACGTTAGTGCTACGTGATCGTAGCGTCGATGACGATCGATGGCTACGTTAGTGCTACGTGATCGTAGCGTCGATG"
        rna "CGAUUGACUGUGACGUAGCUAUGCUGCGAUGUAUGCUAUGCUGAUGCUGGAUGCCGAUUGACUGUGACGUAGCUAUGCUGCGAUGUAUGCUAUGCUGAUGCUGGAUGC"
        protein "MDEPPGKPLSCEEKEKLKEKLAFLKREYSKTLARLQRAQRAEKIKHSIKKTVEEQDCLSQQDLSPQLKHSEPKNKICVYDKLHIKTHLDEETGEKTSITLDVGPESFNPGDGPGGLPIQR"]
    (are [sequence type] (= (sequence-type sequence) type)
         dna :dna
         rna :rna
         protein :protein)))

(let [ecoli-fasta (file (resource "ecoli.fasta"))
      ecoli-dna (:sequence (first (fasta->clj ecoli-fasta)))]
  (deftest ^{:integration true}
    ecoli-peptide-count
    (let [peps (loop-digest "ECOLI" ecoli-dna {:break-after [:K :R]
                                               :start-with [:M]
                                               :mass-threshold 500
                                               :missed-cleavages 2})]
      (is (= 3234797
             (count peps))))))
