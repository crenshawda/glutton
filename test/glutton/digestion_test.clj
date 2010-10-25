(ns glutton.digestion-test
  [:use [lazytest.describe :only [describe it]]
        [glutton.peptide-record :as pr] :reload-all]
  [:import [glutton.peptide-record Peptide]]
  [:require [glutton
             [glutton :as gu]
             [enzyme-free :as ef]]])

(def break-acids #{:M :K :R :.})

(def test-cases
     {:single-breaks (for [break-acid break-acids]
                       [:S :F break-acid :S])
      :double-breaks (for [break-acid break-acids]
                       [:S :F break-acid break-acid :S])})

; Define small test cases and compile a cannonical amino acid vector
(def test-seq-roll (vec (flatten (vals test-cases))))

; Glutton defaults-- :missed-cleavages 2 :break-after [:K :R] :start-with [:M] :mass-threshold 500
(def test-seq-roll-gu-expected
     (list
      (pr/Peptide. [:S :F :. :S :S :F :R] 0 711.33402 1)
      (pr/Peptide. [:S :F :. :S :S :F :R :S :S :F :K] 0 1160.56145 2)
      (pr/Peptide. [:S :F :. :S :S :F :R :S :S :F :K :S :S :F :M :S :S :F :. :. :S :S :F :R] 0 2411.1004600000006 3)
      (pr/Peptide. [:S :S :F :K :S :S :F :M :S :S :F :. :. :S :S :F :R] 21 1699.7664400000003 2)
      (pr/Peptide. [:S :S :F :M :S :S :F :. :. :S :S :F :R] 33 1250.53901 1)
      (pr/Peptide. [:M :S :S :F :. :. :S :S :F :R] 42 929.40654 1)
      (pr/Peptide. [:S :S :F :. :. :S :S :F :R] 45 798.3660500000001 1)
      (pr/Peptide. [:S :S :F :K :S :S :F :M :S :S :F :. :. :S :S :F :R :R] 21 1855.8675500000004 3)
      (pr/Peptide. [:S :S :F :M :S :S :F :. :. :S :S :F :R :R] 33 1406.64012 2)
      (pr/Peptide. [:M :S :S :F :. :. :S :S :F :R :R] 42 1085.50765 2)
      (pr/Peptide. [:S :S :F :. :. :S :S :F :R :R] 45 954.4671600000001 2)
      (pr/Peptide. [:S :S :F :M :S :S :F :. :. :S :S :F :R :R :S :S :F :K] 33 1855.8675500000002 3)
      (pr/Peptide. [:M :S :S :F :. :. :S :S :F :R :R :S :S :F :K] 42 1534.7350800000002 3)
      (pr/Peptide. [:S :S :F :. :. :S :S :F :R :R :S :S :F :K] 45 1403.6945900000003 3)
      (pr/Peptide. [:R :S :S :F :K] 72 605.3285400000001 1)
      (pr/Peptide. [:R :S :S :F :K :K] 72 733.4235000000001 2)
      (pr/Peptide. [:S :S :F :K :K] 75 577.32239 2)
      (pr/Peptide. [:R :S :S :F :K :K :S :S :F :M :M :S] 72 1403.6689800000004 2)
      (pr/Peptide. [:S :S :F :K :K :S :S :F :M :M :S] 75 1247.56787 2)
      (pr/Peptide. [:K :S :S :F :M :M :S] 87 798.34044 0)
      (pr/Peptide. [:S :S :F :M :M :S] 90 670.2454799999999 0)))

(def test-seq-roll-ef-expected
     (list
      (pr/Peptide. [:S :S :F :. :. :S :S] 45 495.19653000000005 0)
      (pr/Peptide. [:S :F :K :K] 78 490.29035999999996 0)
      (pr/Peptide. [:F :K :K :S] 81 490.29036 0)
      (pr/Peptide. [:S :F :M :M] 93 496.18142 0)
      (pr/Peptide. [:F :M :M :S] 96 496.18142000000006 0)))

; Use take on this to create arbitrary-length random amino-acid sequences for stress-testing and benchmarking
(def amino-acid-generator
     (repeatedly #(let [aas (distinct (vals glutton.lexicon/codon-translation-matrix))
                        aa-count (count aas)]
                    (nth aas (rand aa-count)))))

(describe gu/digest
  (it "takes an aa-seq and using a trypsin digest, returns mass-filtered peptides with < N internal cleavages"
    (= test-seq-roll-gu-expected (gu/digest test-seq-roll))))

(describe ef/digest
  (it "takes an aa-seq, a mass target and tolerance and returns a mass-range targeted enzyme free digestion of a genome"
    (= test-seq-roll-ef-expected (ef/digest test-seq-roll 500 10))))
