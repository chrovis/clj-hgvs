(ns clj-hgvs.test-common
  (:require [clojure.spec.test.alpha :as st]
            clj-hgvs.coordinate
            clj-hgvs.core
            clj-hgvs.mutation))

(st/instrument #?(:clj  (st/enumerate-namespace '[clj-hgvs.coordinate
                                                  clj-hgvs.core
                                                  clj-hgvs.mutation])
                  :cljs '[clj-hgvs.coordinate
                          clj-hgvs.core
                          clj-hgvs.mutation]))
