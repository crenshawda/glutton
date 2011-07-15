(ns glutton.six-frame
  (:use (clojure.java (io :only [file]))
        (clojure (string :only [replace-first]))
        (glutton (translate :only [six-frame-translation])
                 (fasta :only [parse-fasta write-fasta]))))

(defn file-extension
  "Return the file extension of a given file `f`; that is, everything between the
  final `.` and the end of the filename."
  [f]
  (re-find #"(?<=\.).[^\.]*$" (.getName (file f))))

(defn filename-without-extension
  "Given a \"fileable\" thing `f` (such as a `java.io.File` or a filename: \"/tmp/foo/bar.txt\",
  returns the base name of the file (no leading path, and no trailing extension, e.g., \"bar\")"
  [f]
  (let [f (file f)
        [name ext] ((juxt #(.getName %) file-extension) f)]
    (if ext
      (replace-first name
                     (re-pattern (str "." ext "$"))
                     "")
      name)))

(defn -main
  "Given a FASTA file `f` of DNA sequences, generate an output file containing a six-frame
  translation of each sequence in the original file.

  If no `output` file name is given, the results will be sent to a file with the same name as
  the original, with \"_six_frame_translation\" appended to the base name of `f`, with the same
  extension as `f` (e.g. \"testing.fa\" becomes \"testing_six_frame_translation.fa\")."
  ([f]
     (-main f (str (filename-without-extension f)
                   "_six_frame_translation."
                   (file-extension f))))
  ([f output]
     (->> (file f)
          parse-fasta
          (mapcat six-frame-translation)
          (write-fasta output))))
