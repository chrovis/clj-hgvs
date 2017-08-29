(defproject clj-hgvs "0.2.2"
  :description "Clojure(Script) library for handling HGVS"
  :url "https://github.com/chrovis/clj-hgvs"
  :license {:name "Apache License, Version 2.0"
            :url "http://www.apache.org/licenses/LICENSE-2.0"}
  :dependencies [[org.clojure/clojure "1.8.0" :scope "provided"]
                 [org.clojure/clojurescript "1.9.908" :scope "provided"]]
  :plugins [[lein-cljsbuild "1.1.7"]
            [lein-cloverage "1.0.9"]
            [lein-codox "0.10.3"]
            [lein-doo "0.1.7"]]
  :profiles {:1.7 {:dependencies [[org.clojure/clojure "1.7.0"]]}
             :1.9 {:dependencies [[org.clojure/clojure "1.9.0-alpha16"]]}}
  :deploy-repositories [["snapshots" {:url "https://clojars.org/repo/"
                                      :username [:env/clojars_username :gpg]
                                      :password [:env/clojars_password :gpg]}]]
  :cljsbuild {:builds {:test {:source-paths ["src" "test"]
                              :compiler {:output-to "target/testable.js"
                                         :output-dir "target"
                                         :main clj-hgvs.runner
                                         :optimizations :simple
                                         :process-shim false}}}} ; workaround (cf. bensu/doo#141)
  :codox {:namespaces [#"^clj-hgvs\.(?!internal)"]
          :output-path "docs"
          :source-uri "https://github.com/chrovis/clj-hgvs/blob/{version}/{filepath}#L{line}"}
  :doo {:build "test"})
