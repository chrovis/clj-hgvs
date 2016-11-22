(ns clj-hgvs.mutation-test
  (:require #?(:clj [clojure.test :refer :all]
               :cljs [cljs.test :refer-macros [deftest is testing]])
            [clj-hgvs.mutation :as mut]))

(deftest ->long-amino-acid-test
  (testing "converts a short amino acid to a long one"
    (is (= (mut/->long-amino-acid "S") "Ser")))
  (testing "returns itself when a long amino acid is passed"
    (is (= (mut/->long-amino-acid "Ser") "Ser")))
  (testing "returns nil when an illegal string is passed"
    (is (nil? (mut/->long-amino-acid "Foo")))
    (is (nil? (mut/->long-amino-acid "")))
    (is (nil? (mut/->long-amino-acid nil)))))

(deftest ->short-amino-acid-test
  (testing "converts a long amino acid to a short one"
    (is (= (mut/->short-amino-acid "Ser") "S")))
  (testing "returns itself when a short amino acid is passed"
    (is (= (mut/->short-amino-acid "S") "S")))
  (testing "returns nil when an illegal string is passed"
    (is (nil? (mut/->short-amino-acid "Z")))
    (is (nil? (mut/->short-amino-acid "")))
    (is (nil? (mut/->short-amino-acid nil)))))
