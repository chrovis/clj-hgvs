(ns clj-hgvs.runner
  (:require [doo.runner :refer-macros [doo-tests]]
            [clj-hgvs.coordinate-test]
            [clj-hgvs.core-test]
            [clj-hgvs.mutation-test]))

(doo-tests 'clj-hgvs.coordinate-test
           'clj-hgvs.core-test
           'clj-hgvs.mutation-test)
