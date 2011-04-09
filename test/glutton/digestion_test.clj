(ns glutton.digestion-test
  (:use (glutton (peptide :only [peptide]))
        clojure.test)
  (:require (glutton (glutton :as gu)
                     (enzyme-free :as ef))))

(def test-seq-roll [:S :F :. :S
                    :S :F :R :S
                    :S :F :K :S
                    :S :F :M :S
                    :S :F :. :. :S
                    :S :F :R :R :S
                    :S :F :K :K :S
                    :S :F :M :M :S
                    ])

; Glutton defaults-- :missed-cleavages 2 :break-after [:K :R] :start-with [:M] :mass-threshold 500
(def test-seq-roll-gu-expected
     [(peptide [:S :S :F :R :S :S :F :K] 9 944.47157471 2 "" "glutton")
      (peptide [:S :S :F :R :S :S :F :K :S :S :F :M :S :S :F :.] 9 1717.7770047100005 2 "" "glutton")
      (peptide [:S :S :F :K :S :S :F :M :S :S :F :.] 21 1240.54342471 1 "" "glutton")
      (peptide [:S :S :F :M :S :S :F :.] 33 791.3159947099999 0 "" "glutton")
      (peptide [:S :S :F :R :R] 60 651.3452547100001 2 "" "glutton")
      (peptide [:S :S :F :R :R :S :S :F :K] 60 1100.57268471 3 "" "glutton")
      (peptide [:R :S :S :F :K] 72 623.33910471 2 "" "glutton")
      (peptide [:R :S :S :F :K :K] 72 751.43406471 3 "" "glutton")
      (peptide [:S :S :F :K :K] 75 595.3329547100001 2 "" "glutton")
      (peptide [:S :S :F :K :K :S :S :F :M :M :S] 75 1265.5784347100002 2 "" "glutton")
      (peptide [:K :S :S :F :M :M :S] 87 816.35100471 1 "" "glutton")
      (peptide [:S :S :F :M :M :S] 90 688.25604471 0 "" "glutton")])

(def test-seq-roll-ef-expected
     [(peptide [:. :S :S :F :R] 6 495.24414471000006 "n/a" "" "enzyme-free")
      (peptide [:S :S :F :R] 9 495.24414471000006 "n/a" "" "enzyme-free")
      (peptide  [:S :F :R :S] 12 495.24414471000006 "n/a" "" "enzyme-free")
      (peptide  [:F :R :S :S] 15 495.24414471000006 "n/a" "" "enzyme-free")
      (peptide  [:R :S :S :F] 18 495.24414471 "n/a" "" "enzyme-free")
      (peptide  [:. :. :S :S :F :R] 54 495.24414471000006 "n/a" "" "enzyme-free")
      (peptide  [:. :S :S :F :R] 57 495.24414471000006 "n/a" "" "enzyme-free")
      (peptide  [:S :S :F :R] 60 495.24414471000006 "n/a" "" "enzyme-free")
      (peptide  [:R :R :S :S] 69 504.2768447100001 "n/a" "" "enzyme-free")
      (peptide  [:R :S :S :F] 72 495.24414471 "n/a" "" "enzyme-free")
      (peptide  [:S :F :K :K] 78 508.30092471 "n/a" "" "enzyme-free")
      (peptide  [:F :K :K :S] 81 508.30092471000006 "n/a" "" "enzyme-free")])

; Use take on this to create arbitrary-length random amino-acid sequences for stress-testing and benchmarking
(def amino-acid-generator
     (repeatedly #(let [aas (distinct (vals glutton.lexicon/codon-translation-matrix))
                        aa-count (count aas)]
                    (nth aas (rand aa-count)))))

(deftest glutton-digest
  (is (= test-seq-roll-gu-expected (gu/digest test-seq-roll))))

(deftest enzyme-free-digest
  (is (= test-seq-roll-ef-expected (ef/digest test-seq-roll :mass-target 500 :mass-tolerance 10))))
