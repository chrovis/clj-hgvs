(defproject clj-hgvs "0.4.5"
  :description "Clojure(Script) library for handling HGVS"
  :url "https://github.com/chrovis/clj-hgvs"
  :license {:name "Apache License, Version 2.0"
            :url "http://www.apache.org/licenses/LICENSE-2.0"}
  :dependencies [[org.clojure/clojure "1.9.0" :scope "provided"]
                 [org.clojure/clojurescript "1.10.520" :scope "provided"]]
  :plugins [[lein-cljsbuild "1.1.7"]
            [lein-cloverage "1.1.2"]
            [lein-codox "0.10.7"]
            [lein-doo "0.1.11"]]
  :profiles {:dev {:dependencies [[org.clojure/test.check "0.10.0"]
                                  [codox-theme-rdash "0.1.2"]]}
             :1.8 {:dependencies [[org.clojure/clojure "1.8.0"]
                                  [clojure-future-spec "1.9.0"]]}
             :1.9 {:dependencies [[org.clojure/clojure "1.9.0"]]}
             :1.10 {:dependencies [[org.clojure/clojure "1.10.3"]]}}
  :deploy-repositories [["snapshots" {:url "https://clojars.org/repo/"
                                      :username [:env/clojars_username :gpg]
                                      :password [:env/clojars_password :gpg]}]]
  :cljsbuild {:builds {:test {:source-paths ["src" "test"]
                              :compiler {:output-to "target/testable.js"
                                         :output-dir "target"
                                         :main clj-hgvs.runner
                                         :optimizations :none}}}}
  :codox {:project {:name "clj-hgvs"}
          :themes [:rdash]
          :namespaces [#"^clj-hgvs\.(?!internal)"]
          :output-path "docs"
          :source-uri "https://github.com/chrovis/clj-hgvs/blob/{version}/{filepath}#L{line}"}
  :doo {:build "test"}
  :signing {:gpg-key "developer@xcoo.jp"})
