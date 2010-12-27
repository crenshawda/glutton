(ns glutton.file-utils
  [:use
   [clojure.java.io :only [writer]]])

;; THIS: Is stolen from `encode-mongo.import.varsplice
(defn file->records [file]
  (partition 2
             (partition-by #(.startsWith % ">")
                           (line-seq (reader file)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                        ; Fasta file functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn partition-multiseq-fasta
  "Give this a sequence of fasta files and it will parse out the important stuff for each"
  [fasta-lines]
  (partition-by #(.startsWith % ">") fasta-lines))


; TODO: Refactor ALL this to the right places
(defn single-fasta
  [fasta-file-str]
  ; TODO: (first because multiseq gives back a list of 1?
  (let [[head seq] (partition-multiseq-fasta (flatten (ds/read-lines fasta-file-str)))]
    (gu/parse-fasta head seq)))
