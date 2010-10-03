(ns glutton.glutton
  (:use [clojure.contrib.seq-utils :only [indexed]]
        [clojure.contrib.def :only [defnk]]))

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
              [{:sequence [current-aa] :breaks 0 :nucleotide-start (* 3 loc)}])))

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

(defnk digest [aas :missed-cleavages 2 :break-after [:K :R] :start-with [:M]]
  (let [break-after (set (conj break-after nil))
        start-with (set start-with)]
    (digest* (indexed (partition 2 1 (cons nil aas)))
             []
             {:missed-cleavages missed-cleavages
              :break-after break-after
              :start-with start-with})))
