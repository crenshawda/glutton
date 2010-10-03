(ns glutton.glutton
  (:use [clojure.contrib.seq-utils :only [indexed]]
        [clojure.contrib.def :only [defnk]]))

(defrecord Peptide [sequence breaks nucleotide-start])

(defn- start? [starts aa] (contains? starts aa))
(defn- break? [breaks aa] (contains? breaks aa))

(defn- process [candidates current-aa prev-aa loc {:keys [break-after start-with]}]
  (lazy-cat (for [c candidates]
              (let [c (update-in c [:sequence] conj current-aa)]
                (if (break? break-after current-aa)
                  (update-in c [:breaks] inc)
                  c)))
            (if (or (break? break-after prev-aa)    ;after you cleave, you need to begin a new one
                    (start? start-with current-aa)  ;obviously
                    (= :M prev-aa))                 ;N-terminal methionines often get removed
                                                    ; (TODO: does this happen for all organisms?)
              [(Peptide. [current-aa] 0 (* 3 loc))])))

(defn- digest*
  [[[loc [prev-aa current-aa]] & other-aas :as aas] candidates config]
  (let [new-candidates (process candidates current-aa prev-aa loc config)
        {:keys [missed-cleavages break-after]} config]
    (if (seq other-aas)
      (if-not (break? break-after current-aa)
        (recur other-aas new-candidates config)
        (lazy-cat new-candidates (digest* other-aas
                                          (remove #(> (:breaks %) missed-cleavages)
                                                  new-candidates)
                                          config)))
      new-candidates)))

(defnk digest
  "Perform a synthetic enzymatic digest of a peptide sequence.

   aas = sequence of amino acids, given as single-letter code keywords, e.g., [:G :L :Y :T: V].

   Optional Keyword Parameters:
   :missed-cleavages -> Maximum number of internal missed cleavage sites allowed for a candidate peptide.
   Defaults to 2
   :break-after      -> Amino acids that signal a break.  Candidate peptides will end with one of these
   amino acids.  Defaults to [:K :R] (i.e., trypsin digestion).
   :start-with       -> Amino acids that signal the start of a new candidate peptide.
   Defaults to [:M].  Note that the start of the digested peptide sequence always begins a new candidate,
   whether it is in this list or not.

   Returns a lazy sequence of Peptides."
  [aas :missed-cleavages 2 :break-after [:K :R] :start-with [:M]]
  (let [break-after (set (conj break-after nil))
        start-with (set start-with)]
    (digest* (indexed (partition 2 1 (cons nil aas)))
             []
             {:missed-cleavages missed-cleavages
              :break-after break-after
              :start-with start-with})))
