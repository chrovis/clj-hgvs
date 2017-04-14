(defproject clj-hgvs "0.1.0-SNAPSHOT"
  :description "Clojure(Script) library for handling HGVS"
  :url "https://github.com/chrovis/clj-hgvs"
  :license {:name "Apache License, Version 2.0"
            :url "http://www.apache.org/licenses/LICENSE-2.0"}
  :dependencies [[org.clojure/clojure "1.8.0" :scope "provided"]
                 [org.clojure/clojurescript "1.9.518" :scope "provided"]]
  :plugins [[lein-cljsbuild "1.1.5"]
            [lein-codox "0.10.3"]
            [lein-doo "0.1.7"]]
  :profiles {:1.7 {:dependencies [[org.clojure/clojure "1.7.0"]]}
             :1.9 {:dependencies [[org.clojure/clojure "1.9.0-alpha15"]]}}
  :cljsbuild {:builds {:test {:source-paths ["src" "test"]
                              :compiler {:output-to "target/testable.js"
                                         :output-dir "target"
                                         :main clj-hgvs.runner
                                         :optimizations :simple}}}}
  :codox {:namespaces [#"^clj-hgvs\.(?!internal)"]
          :output-path "docs"}
  :doo {:build "test"})
