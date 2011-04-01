(ns glutton.file-utils
  (:use (clojure.java (io :only [copy file reader writer]))))

(defn file->records [file]
  (partition 2
             (partition-by #(.startsWith % ">")
                           (line-seq (reader file)))))

(defn write-lines
  "Takes a sequence of Strings and prints them to file `f`, separated by newlines.

  (This is a re-implementation of the old clojure.contrib.duck-streams/write-lines)"
  [f lines]
  (copy (file f) (apply str (interpose "\n" lines))))

;; Alternate version with JUST peptide and start postion
(defn ->tab-delimited
  "Outputs peptide records into a simple tab-delimited format"
  ([peptide-records]
     (->tab-delimited peptide-records "output.txt"))
  ([peptide-records filename]
     (write-lines filename
                  (for [rec peptide-records]
                    (str
                     (apply str (map name (:sequence rec)))
                     "\t"
                     (:nucleotide-start rec)))))) ;; (interpose "\t" (next (vals rec)))
