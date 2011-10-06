(ns glutton.test.translate
  (:use glutton.translate :reload)
  (:use clojure.test))

(let [test-sequence "ACTGCATACGGCATACGACT"]
  (deftest translate-dna
    (are [frame protein] (= (translate test-sequence :frame frame)
                            protein)
         0 "TAYGIR"
         1 "LHTAYD"
         2 "CIRHTT")
    (let [nucleotides "ACTNNNGCATACGGCATACGACT"]
      (are [frame protein] (= (translate nucleotides :frame frame)
                              protein)
           0 "TXAYGIR"
           1 "XXHTAYD"
           2 "XXIRHTT")))

  (deftest test-translate-fasta
    (let [fasta {:header ">Testing"
                 :sequence test-sequence}]
      (are [frame translated] (= (translate-fasta fasta :frame frame)
                                 translated)
           0 {:header ">Testing Frame 0"
              :sequence "TAYGIR"}
           1 {:header ">Testing Frame 1"
              :sequence "LHTAYD"}
           2 {:header ">Testing Frame 2"
              :sequence "CIRHTT"}))))

(deftest test-complement
  (let [normal-output "TGACGTATGCCGTATGCTGA"
        output-with-ns "TGACGTNNNATGCCGTATGCTGA"]
    (are [input output] (= (complement-nucleotides input)
                           output)
         "ACTGCATACGGCATACGACT"
         normal-output

         "actgcatacggcatacgact"
         normal-output

         "ACTgcatacGGCATACGact"
         normal-output

         ;; with Ns
         "ACTGCANNNTACGGCATACGACT"
         output-with-ns

         "actgcannntacggcatacgact"
         output-with-ns

         "ACTGcanNNTACGGCAtacgACT"
         output-with-ns))

  (are [input] (= (complement-nucleotides))))

(deftest test-reverse-complement
  (let [normal-output "ACTGCATACGGCATACGACT"
        output-with-ns "ACTGCATACGGCATACGNNNACT"]
    (are [input output] (= (reverse-complement input)
                           output)
         "AGTCGTATGCCGTATGCAGT"
         normal-output

         "agtcgtatgccgtatgcagt"
         normal-output

         "AGTCgtatgcCGTATGcagt"
         normal-output

         "AGTNNNCGTATGCCGTATGCAGT"
         output-with-ns

         "agtnnncgtatgccgtatgcagt"
         output-with-ns

         "AGTNnncgtaTGCCGTAtgcaGT"
         output-with-ns)))

(deftest test-reverse-complement-fasta
  (is (= (reverse-complement-fasta {:header ">Testing"
                                    :sequence "ACTGCATACGGCATACGACT"})
         {:header ">Testing Reverse Complement"
          :sequence "AGTCGTATGCCGTATGCAGT"})))

(deftest test-six-frame-translation
  (let [normal-output (seq [{:header ">Testing Frame 0"
                             :sequence "TAYGIR"}
                            {:header ">Testing Frame 1"
                             :sequence "LHTAYD"}
                            {:header ">Testing Frame 2"
                             :sequence "CIRHTT"}
                            {:header ">Testing Reverse Complement Frame 0"
                             :sequence "SRMPYA"}
                            {:header ">Testing Reverse Complement Frame 1"
                             :sequence "VVCRMQ"}
                            {:header ">Testing Reverse Complement Frame 2"
                             :sequence "SYAVCS"}])
        output-with-ns (seq [{:header ">Testing Frame 0"
                             :sequence "TXAYGIR"}
                            {:header ">Testing Frame 1"
                             :sequence "XXHTAYD"}
                            {:header ">Testing Frame 2"
                             :sequence "XXIRHTT"}
                            {:header ">Testing Reverse Complement Frame 0"
                             :sequence "SRMPYXX"}
                            {:header ">Testing Reverse Complement Frame 1"
                             :sequence "VVCRMXX"}
                            {:header ">Testing Reverse Complement Frame 2"
                             :sequence "SYAVCXS"}])]
    (are [input output] (= (six-frame-translation input)
                           output)
         {:header ">Testing"
          :sequence "ACTGCATACGGCATACGACT"}
         normal-output

         {:header ">Testing"
          :sequence "actgcatacggcatacgact"}
         normal-output

         {:header ">Testing"
          :sequence "actgCATAcggcATACGact"}
         normal-output

         {:header ">Testing"
          :sequence "ACTNNNGCATACGGCATACGACT"}
         output-with-ns

         {:header ">Testing"
          :sequence "actnnngcatacggcatacgact"}
         output-with-ns

         {:header ">Testing"
          :sequence "ACtnnngCATACGgcataCGACT"}
         output-with-ns)))

(deftest test-rna-code-conversion
  (let [rna-standard-code
        {"AUG" "M" "UAA" "*" "UGA" "*" "UAG" "*" "GCU" "A" "GCC" "A" "GCA" "A" "GCG" "A"
         "UUA" "L" "UUG" "L" "CUU" "L" "CUC" "L" "CUA" "L" "CUG" "L" "CGU" "R" "CGC" "R"
         "CGA" "R" "CGG" "R" "AGA" "R" "AGG" "R" "AAA" "K" "AAG" "K" "AAU" "N" "AAC" "N"
         "GAU" "D" "GAC" "D" "UUU" "F" "UUC" "F" "UGU" "C" "UGC" "C" "CCU" "P" "CCC" "P"
         "CCA" "P" "CCG" "P" "CAA" "Q" "CAG" "Q" "UCU" "S" "UCC" "S" "UCA" "S" "UCG" "S"
         "AGU" "S" "AGC" "S" "GAA" "E" "GAG" "E" "ACU" "T" "ACC" "T" "ACA" "T" "ACG" "T"
         "GGU" "G" "GGC" "G" "GGA" "G" "GGG" "G" "UGG" "W" "CAU" "H" "CAC" "H" "UAU" "Y"
         "UAC" "Y" "AUU" "I" "AUC" "I" "AUA" "I" "GUU" "V" "GUC" "V" "GUA" "V" "GUG" "V"}]
    (is (= (to-rna-code standard-genetic-code)
           rna-standard-code)
        "DNA to RNA code conversion works properly")
    (is (= (to-rna-code rna-standard-code)
           rna-standard-code)
        "Converting an RNA code to an RNA code should be idemptoent")))
