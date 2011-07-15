(ns glutton.glutton
  ( :use (clojure.contrib [def :only [defnk]])
         (clojure (set :only [superset?]))
         (glutton [peptide :only [dna-peptide-candidate
                                  extend-peptide finish-candidate
                                  rna-peptide-candidate
                                  protein-peptide-candidate]]
                  (translate :only [reverse-complement]))))

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

(def my-codon-translation-matrix
     {"ATG" :M "TAA" :. "TGA" :. "TAG" :. "GCT" :A "GCC" :A "GCA" :A "GCG" :A
      "TTA" :L "TTG" :L "CTT" :L "CTC" :L "CTA" :L "CTG" :L "CGT" :R "CGC" :R
      "CGA" :R "CGG" :R "AGA" :R "AGG" :R "AAA" :K "AAG" :K "AAT" :N "AAC" :N
      "GAT" :D "GAC" :D "TTT" :F "TTC" :F "TGT" :C "TGC" :C "CCT" :P "CCC" :P
      "CCA" :P "CCG" :P "CAA" :Q "CAG" :Q "TCT" :S "TCC" :S "TCA" :S "TCG" :S
      "AGT" :S "AGC" :S "GAA" :E "GAG" :E "ACT" :T "ACC" :T "ACA" :T "ACG" :T
      "GGT" :G "GGC" :G "GGA" :G "GGG" :G "TGG" :W "CAT" :H "CAC" :H "TAT" :Y
      "TAC" :Y "ATT" :I "ATC" :I "ATA" :I "GTT" :V "GTC" :V "GTA" :V "GTG" :V})


(def my-rna-codon-translation-matrix
     {"AUG" :M "UAA" :. "UGA" :. "UAG" :. "GCU" :A "GCC" :A "GCA" :A "GCG" :A
      "UUA" :L "UUG" :L "CUU" :L "CUC" :L "CUA" :L "CUG" :L "CGU" :R "CGC" :R
      "CGA" :R "CGG" :R "AGA" :R "AGG" :R "AAA" :K "AAG" :K "AAU" :N "AAC" :N
      "GAU" :D "GAC" :D "UUU" :F "UUC" :F "UGU" :C "UGC" :C "CCU" :P "CCC" :P
      "CCA" :P "CCG" :P "CAA" :Q "CAG" :Q "UCU" :S "UCC" :S "UCA" :S "UCG" :S
      "AGU" :S "AGC" :S "GAA" :E "GAG" :E "ACU" :T "ACC" :T "ACA" :T "ACG" :T
      "GGU" :G "GGC" :G "GGA" :G "GGG" :G "UGG" :W "CAU" :H "CAC" :H "UAU" :Y
      "UAC" :Y "AUU" :I "AUC" :I "AUA" :I "GUU" :V "GUC" :V "GUA" :V "GUG" :V})

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

(defn loop-digest
  [sequence-id ^String primary-sequence {:keys [mass-threshold
                                           missed-cleavages
                                           break-after
                                           start-with]}]
  (let [all-peptides (atom [])

        length (.length primary-sequence)

        sequence-type (sequence-type primary-sequence)
        _ (println "Sequence Type:" sequence-type)

        nucleotide? #{:dna :rna}

        step-size (if (nucleotide? sequence-type)
                    3 1)

        to-aa (condp = sequence-type
                  :dna my-codon-translation-matrix
                  :rna my-rna-codon-translation-matrix
                  :protein keyword
                    )

        ;; Only if doing DNA
        ^String reverse-complement (reverse-complement primary-sequence)

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
                              (loop [position frame last-aa :-]
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

                                          (and (= :. aa) (break? last-aa))
                                          ;; clear out all candidates, flushing nothing
                                          (var-set candidates [])

                                          (= :. aa) ;; if it's just a stop, then copy all proper mass candidates out
                                          ;; begin anew with no candidates

                                          (do (swap! all-peptides add-filtered-candidates @candidates)
                                              (var-set candidates [])))
                                    (recur (+ position step-size) aa))
                                  (swap! all-peptides add-filtered-candidates @candidates)))))))


        ]


    (condp = sequence-type
        :dna (doseq [strand [:+ :-]]
               (let [nt (if (= strand :+) primary-sequence reverse-complement)
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



;; (defnk digest
;;   "Perform a synthetic enzymatic digest of a peptide sequence.

;;    aas = sequence of amino acids, given as single-letter code keywords, e.g., [:G :L :Y :T: V].

;;    Optional Keyword Parameters:
;;    :missed-cleavages -> Maximum number of internal missed cleavage sites allowed for a candidate peptide.
;;    Defaults to 2
;;    :break-after      -> Amino acids that signal a break.  Candidate peptides will end with one of these
;;    amino acids.  Defaults to [:K :R] (i.e., trypsin digestion).
;;    :start-with       -> Amino acids that signal the start of a new candidate peptide.
;;    Defaults to [:M].  Note that the start of the digested peptide sequence always begins a new candidate,
;;    whether it is in this list or not.
;;    :mass-threshold   -> All candidate peptides with mass lower than this should be removed.  Defaults to 500
;;     daltons.
;;    :source         -> Organism that the genome came from.  Defaults to unknown.

;;    Returns a lazy sequence of Peptides."
;;   [aas :missed-cleavages 2 :break-after [:K :R] :start-with [:M] :mass-threshold 500 :source ""]
;;   (let [break-after (set (conj break-after nil))
;;         start-with (set start-with)
;;         config {:missed-cleavages missed-cleavages
;;                 :break-after break-after
;;                 :start-with start-with
;;                 :mass-threshold mass-threshold
;;                 :source source
;;                 :digestion "glutton"}]
;;     (digest* (indexed (partition 2 1 (cons nil aas)))
;;              []
;;              config)))
