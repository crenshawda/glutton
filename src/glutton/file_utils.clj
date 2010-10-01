(ns glutton.file-utils
  (:require (clojure.contrib [duck-streams :as ds])
            (glutton [genomic-utils :as gu])))

(defn file-roller
  "Give it a directory and it will give a list of all the files within"
  [directory]
  (file-seq (new java.io.File directory)))

(defn suffix-filter
  "Give this a list of files and a suffix it will filter on it"
  [file-list suffix]
  (filter #(.endsWith (.getName %) suffix) file-list))

(defn get-fasta-files
  "Give this a directory string and it will grab all the fasta files out of that directory. This is a convenience composition of the file-roller and suffix-filter methods."
  [dir-string]
  (suffix-filter (file-roller dir-string) ".fasta"))

(defn lines-from-files
  "Give this a list of files and it will return a list of each files lines"
  [files]
  (pmap ds/read-lines files))


; TODO: Refactor ALL this to the right places
(defn single-fasta
  [fasta-file-str]
  ; TODO: (first because multiseq gives back a list of 1?
  (let [[head seq] (partition-multiseq-fasta (flatten (ds/read-lines fasta-file-str)))]
    (gu/parse-fasta head seq)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                        ; Fasta file functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn partition-multiseq-fasta
  "Give this a sequence of fasta files and it will parse out the important stuff for each"
  [fasta-lines]
  (partition-by #(.startsWith % ">") fasta-lines))
