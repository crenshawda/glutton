(ns glutton.glutton
  (:use [clojure.contrib.seq-utils :only [indexed]]
        [glutton [peptide-record]])
  (:import [glutton.peptide-record Peptide]))

(defn create-peptide [aa position]
  (Peptide. [aa] position 0 0))

(def breaks #{:K :R :. nil}) ; Consider Proline
(def starts #{:M})
(def max-breaks 2)


(defn- start? [aa] (contains? starts aa))
(defn- break? [aa] (contains? breaks aa))

(defn- process [candidates current-aa prev-aa loc]
  (lazy-cat (for [c candidates]
              (let [c (extend-with c current-aa)]
                (if (break? current-aa)
                  (update-in c [:breaks] inc) ; Should this be a protocol as well?
                  c)))
            (if (or (break? prev-aa)    ;after you cleave, you need to begin a new one
                    (start? current-aa) ;obviously
                    (= :M prev-aa))     ;N-terminal methionines often get removed
                                        ; (TODO: does this happen for all organisms?)
              [(create-peptide [current-aa] (* 3 loc))])))

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
