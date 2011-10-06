(ns glutton.fasta
  "Parsing and writing FASTA files"
  (:use (clojure.java (io :only [file reader writer]))
        (clojure (string :only [join]))))

(defn parse-fasta
  "Parse a FASTA file into a sequence of maps, one for each record present in the file.
  Each map contains the complete record header (including the \">\") under the `:header` key,
  and the concatenation of all the sequence lines under the `:sequence` key."
  [file]
  (for [[header sequence] (partition 2
                                     (partition-by #(.startsWith % ">")
                                                   (line-seq (reader file))))]
    {:header (first header)
     :sequence (join sequence)}))

(defn fasta->lines
  "Given a FASTA record map, creates a sequence of strings for outputting to a file.  The sequence
  is split into 50-character lines."
  [record]
  (list* (:header record)
         (map (partial apply str)
              (partition-all 50 (:sequence record)))))

(defn write-lines
  "Takes a sequence of Strings and prints them to file `f`, separated by newlines.

  (This is a re-implementation of the old `clojure.contrib.duck-streams/write-lines`)"
  [f lines]
  (with-open [w (writer f)]
    (loop [line (first lines) remaining (next lines)]
      (doto w
        (.write line)
        (.newLine))
      (if (seq? remaining)
        (recur (first remaining) (next remaining))))))

(defn write-fasta
  "Write a sequence of FASTA `records` to the file `filename`"
  [filename records]
  (write-lines filename
               (flatten (map fasta->lines records))))
