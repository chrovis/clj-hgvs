(ns clj-hgvs.runner
  (:require [doo.runner :refer-macros [doo-tests]]
            [clj-hgvs.core-test]))

(doo-tests 'clj-hgvs.core-test)
