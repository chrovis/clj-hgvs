(defproject clj-hgvs "0.1.0-SNAPSHOT"
  :description "HGVS library for Clojure(Script)"
  :url "http://example.com/FIXME"
  :dependencies [[org.clojure/clojure "1.8.0"]]
  :plugins [[lein-cljsbuild "1.1.5"]
            [lein-doo "0.1.7"]]
  :profiles {:dev {:dependencies [[org.clojure/clojurescript "1.9.495"]]}}
  :cljsbuild {:builds {:test {:source-paths ["src" "test"]
                              :compiler {:output-to "target/testable.js"
                                         :output-dir "target"
                                         :main clj-hgvs.runner
                                         :optimizations :simple}}}}
  :doo {:build "test"})
