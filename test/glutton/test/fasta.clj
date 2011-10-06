(ns glutton.test.fasta
  (:use glutton.fasta :reload)
  (:use clojure.test)
  (:use (clojure.java (io :only [resource file reader]))))

(let [fasta (file (resource "drosophila-2L-small.fa"))
      fasta-record {:header ">chr2L"
              :sequence "CgacaatgcacgacagaggaagcagaacagatatttagattgcctctcattttctctcccatattatagggagaaatatgatcgcgtatgcgagagtagtgccaacatattgtgctctttgattttttggcaacccaaaatggtggcggatgaaCGAGATGATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCATTGCAACGTTAAATACAGCACAATATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTAATGAGTGCCTCTCGTTCTCTGTCTTATATTACCGCAAACCCAAAAAgacaatacacgacagagagagagagcagcggagatatttagattgcctattaaatatgatcgcgtatgcgagagtagtgccaacatattgtgctctCTATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAAACCCAAatcgacaatgcacgacagaggaagcagaacagatatttagattgcctctcattttctctcccatattatagggagaaatatgatcgcgtatgcg"}]
  (deftest read-fasta-file
    (is (= (parse-fasta fasta)
           (seq [fasta-record]))))
  (deftest fasta-lines
    (is (= (fasta->lines fasta-record)
           '(">chr2L"
             "Cgacaatgcacgacagaggaagcagaacagatatttagattgcctctcat"
             "tttctctcccatattatagggagaaatatgatcgcgtatgcgagagtagt"
             "gccaacatattgtgctctttgattttttggcaacccaaaatggtggcgga"
             "tgaaCGAGATGATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCA"
             "TTGCAACGTTAAATACAGCACAATATATGATCGCGTATGCGAGAGTAGTG"
             "CCAACATATTGTGCTAATGAGTGCCTCTCGTTCTCTGTCTTATATTACCG"
             "CAAACCCAAAAAgacaatacacgacagagagagagagcagcggagatatt"
             "tagattgcctattaaatatgatcgcgtatgcgagagtagtgccaacatat"
             "tgtgctctCTATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAA"
             "ACCCAAatcgacaatgcacgacagaggaagcagaacagatatttagattg"
             "cctctcattttctctcccatattatagggagaaatatgatcgcgtatgcg"))))

  (deftest write-fasta-file
    (are [fastas expected-file] (let [temp (java.io.File/createTempFile "testing" "fa")]
                                  (write-fasta temp fastas)
                                  (with-open [original (reader expected-file)
                                              new (reader temp)]
                                    (is (= (line-seq original)
                                           (line-seq new)))))
         [fasta-record] fasta

         [{:header ">Testing Frame 0"
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
           :sequence "SYAVCS"}]
         (file (resource "multiple-test.fa")))))
