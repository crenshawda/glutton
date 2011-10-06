(ns glutton.test.glutton
  (:use glutton.glutton :reload)
  (:use clojure.test)
  (:use (clojure.java (io :only [file resource]))
        (glutton (fasta :only [parse-fasta]))))

(deftest test-sequence-type
  (let [dna "ACGATCGATGGCTACGTTAGTGCTACGTGATCGTAGCGTCGATGACGATCGATGGCTACGTTAGTGCTACGTGATCGTAGCGTCGATG"
        rna "CGAUUGACUGUGACGUAGCUAUGCUGCGAUGUAUGCUAUGCUGAUGCUGGAUGCCGAUUGACUGUGACGUAGCUAUGCUGCGAUGUAUGCUAUGCUGAUGCUGGAUGC"
        protein "MDEPPGKPLSCEEKEKLKEKLAFLKREYSKTLARLQRAQRAEKIKHSIKKTVEEQDCLSQQDLSPQLKHSEPKNKICVYDKLHIKTHLDEETGEKTSITLDVGPESFNPGDGPGGLPIQR"]
    (are [sequence type] (= (sequence-type sequence) type)
         dna :dna
         rna :rna
         protein :protein)))

(let [ecoli-fasta (file (resource "ecoli.fasta"))
      ecoli-dna (:sequence (first (parse-fasta ecoli-fasta)))]
  (deftest ^{:integration true}
    ecoli-peptide-count
    (let [peps (digest "ECOLI" ecoli-dna
                       :break-after ["K" "R"]
                       :start-with ["M"]
                       :mass-threshold 500
                       :missed-cleavages 2)]
      (is (= 3234797
             (count peps))))))



(let [hmg20b-fasta (file (resource "hmg20b.fasta"))
      hmg20b-rna (:sequence (first (parse-fasta hmg20b-fasta)))]
  (deftest ^{:integration true}
    hmg20b-peptide-count
    (let [peps (digest "HMG20B mRNA" hmg20b-rna
                       :break-after ["K" "R"]
                       :start-with ["M"]
                       :mass-threshold 500
                       :missed-cleavages 2)]
      (is (= 602
             (count peps))))))


(let [beta-catenin-fasta (file (resource "beta-catenin.fasta"))
      beta-catenin-aa (:sequence (first (parse-fasta beta-catenin-fasta)))]
  (deftest ^{:integration true}
    beta-catenin-peptide-count
    (let [peps (digest "beta catenin" beta-catenin-aa
                       :break-after ["K" "R"]
                       :start-with ["M"]
                       :mass-threshold 500
                       :missed-cleavages 2)]
      (is ;; this is almost certainly wrong
       (= 304
             (count peps))))))
