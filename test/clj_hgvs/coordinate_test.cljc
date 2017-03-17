(ns clj-hgvs.coordinate-test
  (:require #?(:clj [clojure.test :refer :all]
               :cljs [cljs.test :refer-macros [deftest are is testing]])
            [clj-hgvs.coordinate :as coord]))

(deftest genomic-coordinate-test
  (testing "validates an input and returns GenomicCoordinate"
    (is (= (coord/genomic-coordinate 3) (coord/->GenomicCoordinate 3))))
  (testing "throws an error if an input is illegal"
    (are [p] (thrown? #?(:clj Error, :cljs js/Error) (coord/genomic-coordinate p))
      0
      3.5
      "3"
      nil)))

(deftest mitochondrial-coordinate-test
  (testing "validates an input and returns MitochondrialCoordinate"
    (is (= (coord/mitochondrial-coordinate 3) (coord/->MitochondrialCoordinate 3))))
  (testing "throws an error if an input is illegal"
    (are [p] (thrown? #?(:clj Error, :cljs js/Error) (coord/mitochondrial-coordinate p))
      0
      3.5
      "3"
      nil)))

(deftest cdna-coordinate-test
  (testing "validates inputs and returns CDNACoordinate"
    (are [p o r] (= (coord/cdna-coordinate p o r) (coord/->CDNACoordinate p o r))
      3 0 nil
      3 0 :upstream
      3 0 :downstream
      87 3 nil
      88 -1 nil
      85 3 :upstream
      37 3 :downstream))
  (testing "throws an error if inputs are illegal"
    (are [p o r] (thrown? #?(:clj Error, :cljs js/Error) (coord/cdna-coordinate p o r))
      0 0 nil
      3.5 0 nil
      "3" 0 nil
      nil 0 nil
      3 0 :illegal
      3 8.5 nil
      3 "3" nil)))

(deftest format-cdna-coordinate-test
  (testing "returns a string expression of a CDNA coordinate"
    (are [m s] (= (coord/format m) s)
      (coord/cdna-coordinate 3 0 nil) "3"
      (coord/cdna-coordinate 3 0 :upstream) "-3"
      (coord/cdna-coordinate 3 0 :downstream) "*3"
      (coord/cdna-coordinate 87 3 nil) "87+3"
      (coord/cdna-coordinate 88 -1 nil) "88-1"
      (coord/cdna-coordinate 88 0 nil) "88"
      (coord/cdna-coordinate 85 3 :upstream) "-85+3"
      (coord/cdna-coordinate 37 3 :downstream) "*37+3")))

(deftest ncdna-coordinate-test
  (testing "validates an input and returns NCDNACoordinate"
    (is (= (coord/ncdna-coordinate 3) (coord/->NCDNACoordinate 3))))
  (testing "throws an error if an input is illegal"
    (are [p] (thrown? #?(:clj Error, :cljs js/Error) (coord/ncdna-coordinate p))
      0
      3.5
      "3"
      nil)))

(deftest rna-coordinate-test
  (testing "validates inputs and returns RNACoordinate"
    (are [p s i] (= (coord/rna-coordinate p s i) (coord/->RNACoordinate p s i))
      3 nil nil
      3 nil :upstream
      3 nil :downstream
      87 3 nil
      88 -1 nil
      85 3 :upstream
      37 3 :downstream))
  (testing "throws an error if inputs are illegal"
    (are [p s i] (thrown? #?(:clj Error, :cljs js/Error) (coord/rna-coordinate p s i))
      0 nil nil
      3.5 nil nil
      "3" nil nil
      nil nil nil
      3 nil :illegal
      3 8.5 nil
      3 "3" nil)))

(deftest protein-coordinate-test
  (testing "validates an input and returns ProteinCoordinate"
    (is (= (coord/protein-coordinate 3) (coord/->ProteinCoordinate 3))))
  (testing "throws an error if an input is illegal"
    (are [p] (thrown? #?(:clj Error, :cljs js/Error) (coord/protein-coordinate p))
      0
      3.5
      "3"
      nil)))
