(ns glutton.test.peptide
  (:use glutton.peptide)
  (:use clojure.test)
  (:import (glutton.peptide DNAPeptideCandidate
                            PeptideFromDNA
                            RNAPeptideCandidate
                            PeptideFromRNA
                            ProteinPeptideCandidate
                            PeptideFromProtein)))

(deftest test-dna-peptide-candidate
  (let [candidate (dna-peptide-candidate "M" :+ 123 :monoisotopic-mass "TEST-SEQUENCE" 1000000)]
    (is (= candidate
           (DNAPeptideCandidate. ["M"] 149.05105471000002 :monoisotopic-mass :+ 123 nil 0 "TEST-SEQUENCE" 1000000)))
    (is (= (extend-peptide candidate "W")
           (DNAPeptideCandidate. ["M" "W"] 335.13036471 :monoisotopic-mass :+ 123 nil 0 "TEST-SEQUENCE" 1000000)))
    (is (= (finish-candidate (DNAPeptideCandidate. ["M" "W"] 335.13036471 :monoisotopic-mass :+ 123 nil 0 "TEST-SEQUENCE" 1000000)))
        (PeptideFromDNA. "MW" "TEST-SEQUENCE" :+ 123 129 335.13036471 0)))

  (let [candidate (dna-peptide-candidate "M" :- 1230 :monoisotopic-mass "TEST-SEQUENCE" 1000000)]
    (is (= candidate
           (DNAPeptideCandidate. ["M"] 149.05105471000002 :monoisotopic-mass :- nil 998770 0 "TEST-SEQUENCE" 1000000)))
    (is (= (extend-peptide candidate "W")
           (DNAPeptideCandidate. ["M" "W"] 335.13036471 :monoisotopic-mass :- nil 998770 0 "TEST-SEQUENCE" 1000000)))
    (is (= (finish-candidate (DNAPeptideCandidate. ["M" "W"] 335.13036471 :monoisotopic-mass :- nil 998770 0 "TEST-SEQUENCE" 1000000)))
        (PeptideFromDNA. "MW" "TEST-SEQUENCE" :- 998764 998770 335.13036471 0))))


(deftest test-rna-peptide-candidate
  (let [candidate (rna-peptide-candidate "M" 123 :monoisotopic-mass "TEST-SEQUENCE")]
    (is (= candidate
           (RNAPeptideCandidate. ["M"] 149.05105471000002 :monoisotopic-mass 123 nil "TEST-SEQUENCE" 0)))
    (is (= (extend-peptide candidate "W")
           (RNAPeptideCandidate. ["M" "W"] 335.13036471 :monoisotopic-mass 123 nil "TEST-SEQUENCE" 0)))
    (is (= (finish-candidate (RNAPeptideCandidate. ["M" "W"] 335.13036471 :monoisotopic-mass 123 nil "TEST-SEQUENCE" 0)))
        (PeptideFromRNA. "MW" "TEST-SEQUENCE" 123 129 335.13036471 0))))


(deftest test-protein-peptide-candidate
  (let [candidate (protein-peptide-candidate "M" 123 :monoisotopic-mass "TEST-SEQUENCE")]
    (is (= candidate
           (ProteinPeptideCandidate. ["M"] 149.05105471000002 :monoisotopic-mass 123 nil "TEST-SEQUENCE" 0)))
    (is (= (extend-peptide candidate "W")
           (ProteinPeptideCandidate. ["M" "W"] 335.13036471 :monoisotopic-mass 123 nil "TEST-SEQUENCE" 0)))
    (is (= (finish-candidate (ProteinPeptideCandidate. ["M" "W"] 335.13036471 :monoisotopic-mass 123 nil "TEST-SEQUENCE" 0)))
        (PeptideFromProtein. "MW" "TEST-SEQUENCE" 123 129 335.13036471 0))))
