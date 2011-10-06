;; glutton v0.2

(ns glutton.glutton
  ( :use (clojure (set :only [superset?]))
         (glutton [peptide :only [dna-peptide-candidate
                                  extend-peptide finish-candidate
                                  rna-peptide-candidate
                                  protein-peptide-candidate]]
                  (translate :only [reverse-complement standard-genetic-code to-rna-code STOP]))))

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
  (= aa STOP))

(defn sequence-type
  "Based on sampling the first 100 monomers of the given biopolymer sequence,
  determines if the sequence is DNA, RNA, or protein."
  [^String sequence]
  (let [dna-bases #{\A \C \G \T}
        rna-bases #{\A \C \G \U}
        amino-acids #{\A \R \N \D \C
                      \E \Q \G \H \I
                      \L \K \M \F \P
                      \S \T \W \Y \V}
        monomers (set (take 100 sequence))]
    (condp superset? monomers
      dna-bases :dna
      rna-bases :rna
      amino-acids :protein)))

(defn digest
  "Perform a synthetic trypsin digest of a sequence.

   If the sequence is DNA, it is translated to protein in all six reading frames before digestion.
   If the sequence is RNA, it is translated to protein in all three reading frames before digestion.
   If the sequence is protein, it is just digested.

   * `sequence-id` is an identifier for the sequence.
   * `sequence` is the actual biopolymer sequence being digested, as a string.  Standard single-letter
   codes should be used

   Keyword parameters

   * `mass-threshold`: only peptides that have a mass greater than or equal to this (in Daltons) will
     be returned.  Useful for ignoring peptides that are too small to be of interest.  Defaults to `500`.
   * `missed-cleavages`: proteases don't cleave with 100% efficiency.  This parameter simulates this by
     allowing a peptide to contain some number of internal cleavage sites.  Defaults to 2.
   * `break-after`: a vector of single-letter amino acid codes after which a peptide should be cleaved.
     Defaults to `[\"K\" \"R\"] for trypsin.
   * `start-with`: a vector of single-letter amino acid codes at which a new peptide should be started.
     Defaults to `[\"M\"]` (methionine).  The first amino acid of the (possibly translated) sequence will
     begin a new peptide, whether it is in this vector or not.

  Returns a non-lazy sequence of peptides (see `glutton.peptide`)."
  [sequence-id ^String primary-sequence & {:keys [mass-threshold
                                                  missed-cleavages
                                                  break-after
                                                  start-with]
                                           :or {mass-threshold 500
                                                missed-cleavages 2
                                                break-after ["K" "R"]
                                                start-with ["M"]}}]
  (let [all-peptides (atom [])

        length (.length primary-sequence)

        sequence-type (sequence-type primary-sequence)
        _ (println "Sequence Type:" sequence-type)

        step-size (if (contains? #{:dna :rna} sequence-type)
                    3 1)

        to-aa (condp = sequence-type
                  :dna standard-genetic-code
                  :rna (to-rna-code standard-genetic-code)
                  :protein identity)

        break? (set break-after)
        start? (set start-with)

        extend-candidates (fn [candidates start aa last-aa strand]
                            (loop [i 0
                                   cs candidates]
                              (if (not= i (count candidates))
                                (recur (inc i) (assoc cs
                                                 i
                                                 (-> (get cs i)
                                                     (extend-peptide aa)
                                                     (#(if (break? aa)
                                                         (update-in % [:breaks] inc)
                                                         %)))))
                                (if (or (break? last-aa)
                                        (start? aa)
                                        (= "M" last-aa) ;; N-term Ms are often cleaved
                                        (= STOP last-aa))
                                  (conj cs
                                        (condp = sequence-type
                                            :dna (dna-peptide-candidate aa
                                                                        strand
                                                                        start
                                                                        :monoisotopic-mass
                                                                        sequence-id
                                                                        length
                                                                        (if (break? aa) 1 0))
                                            :rna (rna-peptide-candidate aa
                                                                        start
                                                                        :monoisotopic-mass
                                                                        sequence-id
                                                                        (if (break? aa) 1 0))
                                            :protein (protein-peptide-candidate aa
                                                                                start
                                                                                :monoisotopic-mass
                                                                                sequence-id
                                                                                (if (break? aa) 1 0))))
                                  cs))))

        above-threshold (fn [peptide-candidates]
                          (filter #(>= (:mass %) mass-threshold)
                                  peptide-candidates))

        add-filtered-candidates (fn [all-peptides candidates]
                                  (into all-peptides (map finish-candidate (above-threshold candidates))))

        remove-breaks (fn [peptide-candidates]
                        (vec (remove #(> (:breaks %) missed-cleavages)
                                     peptide-candidates)))


        ;; needs next-chunk, which is a function that depends on the sequence
        the-main-event (fn [frame next-chunk strand]
                         (time
                          (do
                            (with-local-vars [candidates []] ;; these are thread-local, non-interned vars
                              (loop [position frame
                                     last-aa :-]
                                (if (not (>= position (- length step-size))) ;; TODO: check this math

                                  ;; translate codon to amino acid; add that to list
                                  (let [aa (to-aa (next-chunk position))]

                                    ;; STRAND!! That's only for DNA!
                                    ;; consider passing a map of extra information instead of individual params
                                    (var-set candidates (extend-candidates @candidates position aa last-aa strand))
                                    (cond (break? aa)
                                          (do
                                            ;; copy all candidate peptides that satisfy the mass filter to peptides collection
                                            (swap! all-peptides add-filtered-candidates @candidates)
                                            ;; Then, remove any candidates that have all their breaks
                                            (var-set candidates
                                                     (remove-breaks @candidates)))

                                          (and (= STOP aa) (break? last-aa))
                                          ;; clear out all candidates, flushing nothing
                                          (var-set candidates [])

                                          (= STOP aa) ;; if it's just a stop, then copy all proper mass candidates out
                                          ;; begin anew with no candidates

                                          (do (swap! all-peptides add-filtered-candidates @candidates)
                                              (var-set candidates [])))
                                    (recur (+ position step-size) aa))
                                  (swap! all-peptides add-filtered-candidates @candidates)))))))]
    (condp = sequence-type
        :dna (doseq [strand [:+ :-]]
               (let [nt (if (= strand :+) primary-sequence (reverse-complement primary-sequence))
                     next-chunk (fn [position] (.substring nt position (+ position step-size)))]
                 (doseq [frame [0 1 2]]
                   (the-main-event frame next-chunk strand))))
        :rna (doseq [strand [:+]]
               (let [nt primary-sequence
                     next-chunk (fn [position] (.substring nt position (+ position step-size)))]
                 (doseq [frame [0 1 2]]
                   (println "Digesting Frame" frame "of" sequence-id)
                   (the-main-event frame next-chunk strand))))
        :protein (let [nt primary-sequence
                       next-chunk (fn [position] (.substring nt position (+ position step-size)))]
                   (println "Digesting " sequence-id)
                   (the-main-event 0 next-chunk nil)))
    @all-peptides))
