(ns glutton.glutton
  (:use [clojure.contrib.seq-utils :only [indexed]])
  (:use [clojure.contrib.duck-streams :only [read-lines]]))

(def begin-flag :START)
(def breaks #{:K :R :. begin-flag}) ; Consider Proline
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
            (if (or (break? prev-aa)
                    (start? current-aa))
              [{:sequence [current-aa] :breaks 0 :nucleotide-start (* 3 loc)}])))

(defn- digest*
  [[[loc [prev-aa current-aa]] & other-aas :as aas] candidates]
  (if (seq other-aas)
    (let [candidates (process candidates current-aa prev-aa loc)]
      (if-not (break? current-aa)
        (recur other-aas candidates)
        (lazy-cat candidates (digest* other-aas (remove #(> (:breaks %) max-breaks)
                                                        candidates)))))
    (process candidates current-aa prev-aa loc)))

(defn digest
  [aas]
  (digest* (indexed (partition 2 1 (cons begin-flag aas))) []))
