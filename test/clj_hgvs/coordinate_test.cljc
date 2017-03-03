(ns clj-hgvs.coordinate-test
  (:require #?(:clj [clojure.test :refer :all]
               :cljs [cljs.test :refer-macros [deftest are testing]])
            [clj-hgvs.coordinate :as coord]))

(deftest format-cdna-coordinate-test
  (testing "returns a string expression of a CDNA coordinate"
    (are [m s] (= (coord/format m) s)
      (coord/->CDNACoordinate 3 nil nil) "3"
      (coord/->CDNACoordinate 3 :up nil) "-3"
      (coord/->CDNACoordinate 3 :down nil) "*3"
      (coord/->CDNACoordinate 87 nil 3) "87+3"
      (coord/->CDNACoordinate 88 nil -1) "88-1"
      (coord/->CDNACoordinate 88 nil 0) "88"
      (coord/->CDNACoordinate 85 :up 3) "-85+3"
      (coord/->CDNACoordinate 37 :down 3) "*37+3")))
