(defproject glutton "0.0.1"
  :description "Fast, efficient, extensible genome digestion algorithm"
  :dependencies [[org.clojure/clojure "1.2.0"]
		 [org.clojure/clojure-contrib "1.2.0"]]
  :dev-dependencies [[swank-clojure "1.2.1"]
                     [org.clojars.ninjudd/lazytest "1.1.3-SNAPSHOT"]]
  :jvm-opts ["-server"
             "-Xmx1g"])
