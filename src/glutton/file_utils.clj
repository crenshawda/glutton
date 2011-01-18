(ns glutton.file-utils
  (:use (clojure.contrib [duck-streams :only [write-lines]])
        (clojure.java [io :only [writer reader]])))

;; THIS: Is stolen from `encode-mongo.import.varsplice
(defn file->records [file]
  (partition 2
             (partition-by #(.startsWith % ">")
                           (line-seq (reader file)))))

;; TODO: Figgur out a way to do this without duck-streams
(defn ->tab-delimited
  "Outputs peptide records into a simple tab-delimited format"
  ([peptide-records]
     (->tab-delimited peptide-records "output.txt"))
  ([peptide-records filename]
     (write-lines filename
                  (map #(apply str (interpose "\t" (vals %)))
                       peptide-records))))
