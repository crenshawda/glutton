(ns glutton.test.genomic-utils
  (:use glutton.genomic-utils :reload)
  (:use clojure.test))

(let [sequence "ACTGCATAGGGG"]
  (deftest test-complement
    (is (= "TGACGTATCCCC"
           (complement-nucleotides sequence))))

  (deftest test-reverse-complement
    (is (= "CCCCTATGCAGT"
           (reverse-complement sequence)))))
