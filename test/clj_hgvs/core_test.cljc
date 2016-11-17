(ns clj-hgvs.core-test
  (:require #?(:clj [clojure.test :refer :all]
               :cljs [cljs.test :refer-macros [deftest is testing]])
            [clj-hgvs.core :as hgvs]))

(deftest parse-test
  (testing
      (is (= (hgvs/parse "NM_005228.3:c.2361G>A") {:transcript "NM_005228.3"
                                                   :kind :coding-dna
                                                   :mutations '({:numbering "2361"
                                                                 :type :substitution
                                                                 :ref "G"
                                                                 :alt "A"})}))))
