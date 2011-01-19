(ns glutton.glutton
  ( :use (clojure.contrib [def :only [defnk]]
                          [seq-utils :only [indexed]])
         (glutton [peptide-record :only [extend-with initiate-peptide]])))

(defn- above-threshold
  "Return only those peptides whose mass is greater than or equal to 'mass-threshold'"
  [mass-threshold peptides]
  (filter #(>= (:mass %)
               mass-threshold)
          peptides))

(defn- start? [starts aa]
  (contains? starts aa))

(defn- break? [breaks aa]
  (contains? breaks aa))

(defn- stop-codon? [aa]
  (= aa :.))

(defn- process [candidates current-aa prev-aa loc {:keys [break-after start-with source digestion]}]
  (lazy-cat (for [c candidates]
              (let [c (extend-with c current-aa)]
                (if (break? break-after current-aa)
                  (update-in c [:breaks] inc) ; Should this be a protocol as well?
                  c)))
            (if (or (break? break-after prev-aa) ;after you cleave, you need to begin a new one
                    (start? start-with current-aa) ;obviously
                    (= :M prev-aa) ;N-terminal methionines often get removed                                         ; (TODO: does this happen for all organisms?)
                    (stop-codon? prev-aa))
              [(initiate-peptide current-aa (* 3 loc)
                                 (if (break? break-after current-aa)
                                   1
                                   0)
                                 source digestion)])))

(defn- digest*
  [[[loc [prev-aa current-aa]] & other-aas :as aas] candidates config]
  (let [new-candidates (process candidates current-aa prev-aa loc config)
        {:keys [missed-cleavages break-after mass-threshold]} config]
    (if (seq other-aas)
      (cond (break? break-after current-aa)  (lazy-cat (above-threshold mass-threshold new-candidates)
                                                       (digest* other-aas
                                                                (remove #(> (:breaks %)
                                                                            missed-cleavages)
                                                                        new-candidates)
                                                                config))
            (and (stop-codon? current-aa)
                 (break? break-after prev-aa)) (lazy-cat [] (digest* other-aas [] config))
            (stop-codon? current-aa) (lazy-cat (above-threshold mass-threshold new-candidates)
                                               (digest* other-aas [] config))

            :default (recur other-aas new-candidates config))
      (above-threshold mass-threshold new-candidates))))

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
   :mass-threshold   -> All candidate peptides with mass lower than this should be removed.  Defaults to 500
    daltons.
   :source         -> Organism that the genome came from.  Defaults to unknown.

   Returns a lazy sequence of Peptides."
  [aas :missed-cleavages 2 :break-after [:K :R] :start-with [:M] :mass-threshold 500 :source ""]
  (let [break-after (set (conj break-after nil))
        start-with (set start-with)
        config {:missed-cleavages missed-cleavages
                :break-after break-after
                :start-with start-with
                :mass-threshold mass-threshold
                :source source
                :digestion "glutton"}]
    (digest* (indexed (partition 2 1 (cons nil aas)))
             []
             config)))
