(ns glutton.glutton
  (:use [clojure.contrib.seq-utils :only [indexed]]))

(def breaks #{:K :R :. nil}) ; Consider Proline
(def starts #{:M})
(def max-breaks 2)

(defn- start? [aa] (contains? starts aa))
(defn- break? [aa] (contains? breaks aa))

(defn- process [candidates current-aa prev-aa loc]
  (lazy-cat (for [c candidates]
              (let [c (update-in c [:sequence] conj current-aa)]
                (if (break? current-aa)
                  (update-in c [:breaks] inc)
                  c)))
            (if (or (break? prev-aa)    ;after you cleave, you need to begin a new one
                    (start? current-aa) ;obviously
                    (= :M prev-aa))     ;N-terminal methionines often get removed
                                        ; (TODO: does this happen for all organisms?)
                     [{:sequence [current-aa] :breaks 0 :nucleotide-start (* 3 loc)}])))

(defn- digest*
  [[[loc [prev-aa current-aa]] & other-aas :as aas] candidates]
  (let [new-candidates (process candidates current-aa prev-aa loc)]
    (if (seq other-aas)
      (if-not (break? current-aa)
        (recur other-aas new-candidates)
        (lazy-cat new-candidates (digest* other-aas (remove #(> (:breaks %) max-breaks)
                                                        new-candidates))))
      new-candidates)))

(defn digest [aas]
  (digest* (indexed (partition 2 1 (cons nil aas))) []))
