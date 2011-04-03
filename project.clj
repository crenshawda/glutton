(defproject glutton "0.0.1"
  :description "Fast, efficient, extensible genome digestion algorithm"
  :dependencies [[org.clojure/clojure "1.3.0-alpha4"]
                 [org.clojure.contrib/standalone "1.3.0-alpha4"]
                 [robert/hooke "1.1.0"]]
  :dev-dependencies [[swank-clojure "1.3.0-SNAPSHOT"]]
  :test-selectors {:default (complement :integration)
                   :integration :integration
                   :all (constantly true)}
  :jvm-opts ["-server"
             "-Xmx2g"])
