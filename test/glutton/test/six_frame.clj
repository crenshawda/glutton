(ns glutton.test.six-frame
  (:use glutton.six-frame :reload)
  (:use clojure.test)
  (:use (clojure.java (io :only [file reader resource]))))

(deftest test-file-translation
  (let [fasta (file (resource "drosophila-2L-small.fa"))
        output (java.io.File/createTempFile "testing" "fa")]
    (-main fasta output)
    (with-open [expected (reader (file (resource "drosophila-2L-small_six_frame_translation.fa")))
                actual (reader output)]
      (is (= (line-seq expected)
             (line-seq actual))))))
