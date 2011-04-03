(ns glutton.genomic-utils
  (:use (clojure (string :only [upper-case join]
                         :as string))
        (glutton [lexicon :only [amino-acid-dictionary
                                 codon-translation-matrix
                                 water-constants]]
                 [file-utils :only [file->records]])))

;;;;;;;;;;;;;;;;
;; Mass Utils ;;
;;;;;;;;;;;;;;;;

(defn aa-mass
  ([amino-acid]
     (aa-mass amino-acid :monoisotopic-mass))
  ([amino-acid mass-type]
     (mass-type (amino-acid-dictionary amino-acid))))

(defn water-mass
  ([]
     (water-mass :monoisotopic-mass))
  ([mass-type]
     (mass-type water-constants)))

;;;;;;;;;;;;;;;;;
;; FASTA Utils ;;
;;;;;;;;;;;;;;;;;

(defrecord FASTA [header sequence])

(defn fasta->clj
  "Returns a sequence of fasta-records, works on multi-seq FASTA files"
  [file-str]
  (for [fasta-seq (file->records file-str)]
    (let [[header sequence] fasta-seq
          header (first header)
          sequence (join sequence)]
      (FASTA. header sequence))))

;;;;;;;;;;;;;;;;;;;;;;
;; Nucleotide Utils ;;
;;;;;;;;;;;;;;;;;;;;;;

(defn- reverse-char-seq
  "So frame can have a completely lazy genome seq both ways"
  [s]
  (letfn [(rec [i]
               (lazy-seq
                (cons (.charAt s i)
                      (if (pos? i)
                        (rec (dec i))))))]
    (rec (dec (.length s)))))

(defn- frame [num dir s]
  {:pre [(contains? #{:F :R} dir)
         (contains? #{0 1 2} num)]}
  (partition 3 (drop num (if (= :F dir)
                           s
                           (reverse-char-seq s)))))

(defn complement-nucleotides [^String nucleotides]
  (apply str (replace {\A \T \T \A \C \G \G \C} nucleotides)))

(defn reverse-complement [^String nucleotides]
  (complement-nucleotides (string/reverse nucleotides)))

(defn- as-codon [s] (keyword (upper-case (join s))))
(defn- ->aa [codon] (codon codon-translation-matrix))

;; Killing the JVM STILL unclogs a few 100K AAs after choking... something's up here...
(defn nucleotides->frames
  "Returns a sequence of amino-acid frames translated from a sequence of nucleotides.
   NOTE: Reading Frame Indexes 0-2 Forward, 3-5 Reverse-compliment"
  [fasta-record]
  (for [dir #{:F :R}
        num (range 3)]
    (let [frame (frame num dir (:sequence fasta-record))]
      (if (= dir :F)
        (->> frame
             (map as-codon)
             (map ->aa))
        (->> frame
             (map complement-nucleotides)
             (map as-codon)
             (map ->aa))))))
