(ns glutton.file-utils
  [:use
   [clojure.java.io :only [writer]]])

;; THIS: Is stolen from `encode-mongo.import.varsplice
(defn file->records [file]
  (partition 2
             (partition-by #(.startsWith % ">")
                           (line-seq (reader file)))))

(defn ->tab-delimited
  "Outputs peptide records into a simple tab-delimited format"
  ([peptide-records]
     (->tab-delimited peptide-records "output.txt"))
  ([peptide-records filename]
     ;; Doesn't check to see if the file already exists before appending
     (with-open [w (writer filename :append true)]
       (doseq [frame peptide-records]
         (doseq [record frame]
           (.write w
                 (str (apply str (interpose "\t" (vals record))) "\n")))))))
