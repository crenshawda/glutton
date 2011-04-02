(ns glutton.glutton
  ( :use (clojure.contrib [def :only [defnk]])
         (glutton [peptide-record :only [extend-with initiate-peptide]]
                  [util :only [indexed]])))

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

(def my-codon-translation-matrix
     {"ATG" :M ; synonymous with :START codon
      "TAA" :.
      "TGA" :.
      "TAG" :. ; synonymous with :STOP codon
      "GCT" :A
      "GCC" :A
      "GCA" :A
      "GCG" :A
      "TTA" :L
      "TTG" :L
      "CTT" :L
      "CTC" :L
      "CTA" :L
      "CTG" :L
      "CGT" :R
      "CGC" :R
      "CGA" :R
      "CGG" :R
      "AGA" :R
      "AGG" :R
      "AAA" :K
      "AAG" :K
      "AAT" :N
      "AAC" :N
      "GAT" :D
      "GAC" :D
      "TTT" :F
      "TTC" :F
      "TGT" :C
      "TGC" :C
      "CCT" :P
      "CCC" :P
      "CCA" :P
      "CCG" :P
      "CAA" :Q
      "CAG" :Q
      "TCT" :S
      "TCC" :S
      "TCA" :S
      "TCG" :S
      "AGT" :S
      "AGC" :S
      "GAA" :E
      "GAG" :E
      "ACT" :T
      "ACC" :T
      "ACA" :T
      "ACG" :T
      "GGT" :G
      "GGC" :G
      "GGA" :G
      "GGG" :G
      "TGG" :W
      "CAT" :H
      "CAC" :H
      "TAT" :Y
      "TAC" :Y
      "ATT" :I
      "ATC" :I
      "ATA" :I
      "GTT" :V
      "GTC" :V
      "GTA" :V
      "GTG" :V})


(defn- amino-acid [codon]
  (my-codon-translation-matrix (String. codon)))









(defn- loop-digest
  [^String nucleotides {:keys [mass-threshold
                               missed-cleavages
                               break-after
                               start-with]}]
  (let [all-peptides (atom [])
        candidates (atom [])

        length (.length nucleotides)

        break? (set break-after)

        start? (set start-with)

        extend-candidates (fn [candidates start aa last-aa] ;;TODO start might be different on reverse
                            (let [len (count candidates)]
                              (loop [i 0
                                     cs candidates]
                                (if (not= i len)
                                  (recur (inc i) (assoc cs
                                                         i
                                                         (-> (get cs i)
                                                             (extend-with aa)
                                                             (#(if (break? aa)
                                                                 (update-in % [:breaks] inc)
                                                                 %)))))
                                   (if (or (break? last-aa)
                                             (start? aa)
                                             (= :M last-aa)
                                             (= :. last-aa))
                                       (conj cs (initiate-peptide aa start
                                                                  (if (break? aa) 1 0)
                                                                  "SOURCE" "GLUTTON"))
                                       cs)))))


        add-filtered-candidates (fn [all-peptides]
                                  (into all-peptides (above-threshold mass-threshold @candidates)))
        remove-breaks (fn [peptide-candidates]
                        (vec (remove #(> (:breaks %)
                                          missed-cleavages)
                                      peptide-candidates)))
        empty-candidates (fn [_] [])
        ]

    (loop [position 0 codon (char-array 3) codon-index 0 last-aa :-]
      (if (not (= position length))
        (if (= codon-index 3)
          ;; translate codon to amino acid; add that to list
          (let [aa (amino-acid codon)]
            (swap! candidates extend-candidates (- position 3) aa last-aa)

            (cond (break? aa)
                  (do
                     ;; copy all candidate peptides that satisfy the mass filter to peptides collection
                    (swap! all-peptides add-filtered-candidates)
                    ;; Then, remove any candidates that have all their breaks
                    (swap! candidates remove-breaks))

                  (and (= :. aa) (break? last-aa))
                  ;; clear out all candidates, flushing nothing
                  (swap! candidates empty-candidates)

                  (= :. aa) ;; if it's just a stop, then copy all proper mass candidates out
                  ;; begin anew with no candidates

                  (do (swap! all-peptides add-filtered-candidates)
                      (swap! candidates empty-candidates)))
            (recur (inc position)
                   (char-array 3 [(.charAt nucleotides position) \- \-])
                   1
                   aa))
          (do (aset-char codon codon-index (.charAt nucleotides position))
              (recur (inc position)
                     codon
                     (inc codon-index)
                     last-aa)))
        (do (swap! all-peptides add-filtered-candidates)
            @all-peptides)))))



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
