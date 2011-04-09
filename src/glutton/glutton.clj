(ns glutton.glutton
  ( :use (clojure.contrib [def :only [defnk]])
         (glutton [peptide :only [extend-with initiate-peptide dna-peptide-candidate extend-peptide finish-candidate]]
                  [util :only [indexed]]
                  (genomic-utils :only [reverse-complement]))))

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
     {"ATG" :M
      "TAA" :.
      "TGA" :.
      "TAG" :.
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

(defn loop-digest
  [sequence-id ^String nucleotides {:keys [mass-threshold
                                           missed-cleavages
                                           break-after
                                           start-with]}]
  (let [all-peptides (atom [])
        length (.length nucleotides)
        ^String reverse-complement (reverse-complement nucleotides)
        break? (set break-after)
        start? (set start-with)
        extend-candidates (fn [candidates start aa last-aa strand]
                            (loop [i 0 cs candidates]
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
                                        (= :M last-aa) ;; N-term Ms are often cleaved
                                        (= :. last-aa))
                                  (conj cs
                                        (dna-peptide-candidate aa
                                                               strand
                                                               start
                                                               :monoisotopic-mass
                                                               sequence-id
                                                               length
                                                               (if (break? aa) 1 0)))
                                  cs))))
        above-threshold (fn [peptide-candidates]
                          (filter #(>= (:mass %) mass-threshold)
                                  peptide-candidates))
        add-filtered-candidates (fn [all-peptides candidates]
                                  (into all-peptides (map finish-candidate (above-threshold candidates))))
        remove-breaks (fn [peptide-candidates]
                        (vec (remove #(> (:breaks %) missed-cleavages)
                                     peptide-candidates)))]
    ;; TODO: coordinates will be screwed up on the reverse
    (doseq [strand [:+ :-]]
      (let [nt (if (= strand :+) nucleotides reverse-complement)
            get-codon (fn [position] (.substring nt position (+ position 3)))]
        (doseq [frame [0 1 2]]
          (time
           (do
             (println "Digesting frame" frame "of strand" (name strand))
             (with-local-vars [candidates []] ;; these are thread-local, non-interned vars
               (loop [position frame last-aa :-]
                 (if (not (>= position (- length 3))) ;; TODO: check this math

                   ;; translate codon to amino acid; add that to list

                   (let [aa (my-codon-translation-matrix (get-codon position))]
                     (var-set candidates (extend-candidates @candidates position aa last-aa strand))
                     (cond (break? aa)
                           (do
                             ;; copy all candidate peptides that satisfy the mass filter to peptides collection
                             (swap! all-peptides add-filtered-candidates @candidates)
                             ;; Then, remove any candidates that have all their breaks
                             (var-set candidates
                                      (remove-breaks @candidates)))

                           (and (= :. aa) (break? last-aa))
                           ;; clear out all candidates, flushing nothing
                           (var-set candidates [])

                           (= :. aa) ;; if it's just a stop, then copy all proper mass candidates out
                           ;; begin anew with no candidates

                           (do (swap! all-peptides add-filtered-candidates @candidates)
                               (var-set candidates [])))
                     (recur (+ position 3) aa))
                   (swap! all-peptides add-filtered-candidates @candidates)))))))))
    @all-peptides))



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
