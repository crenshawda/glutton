(defproject glutton "0.0.1"
  :description "Fast, efficient, extensible genome digestion algorithm"
  :dependencies [[org.clojure/clojure "1.3.0-alpha4"]
                 [robert/hooke "1.1.0"]]
  :dev-dependencies [[lein-difftest "1.3.1"]
                     [lein-retest "1.0.1"]
                     [lein-test-out "0.1.0"]
                     [swank-clojure "1.3.0-SNAPSHOT"]]
  :hooks [leiningen.hooks.difftest
          leiningen.hooks.retest]
  :test-selectors {:default (complement :integration)
                   :integration :integration
                   :all (constantly true)}
  :jvm-opts ["-server"
             "-Xmx2g"])
