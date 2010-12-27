(ns glutton.file-utils
  [:use
   [clojure.java.io :only [writer]]])

;; THIS: Is stolen from `encode-mongo.import.varsplice
(defn file->records [file]
  (partition 2
             (partition-by #(.startsWith % ">")
                           (line-seq (reader file)))))

