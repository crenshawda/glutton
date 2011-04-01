(ns glutton.util)

(defn indexed [s]
  (map vector (iterate inc 0) s))
