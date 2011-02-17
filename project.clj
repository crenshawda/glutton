(defproject glutton "0.0.1"
  :description "Fast, efficient, extensible genome digestion algorithm"
  :dependencies [[org.clojure/clojure "1.3.0-alpha4"]
                 [org.clojure.contrib/standalone "1.3.0-alpha4"]]
  :dev-dependencies [[swank-clojure "1.3.0-SNAPSHOT"]
                     [org.clojars.ninjudd/lazytest "1.1.3-SNAPSHOT"]]
  :jvm-opts ["-server"
             "-Xmx1g"])
