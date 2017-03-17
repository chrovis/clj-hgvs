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
    (are [p s i] (= (coord/cdna-coordinate p s i) (coord/->CDNACoordinate p s i))
      3 nil nil
      3 :up nil
      3 :down nil
      87 nil 3
      88 nil -1
      85 :up 3
      37 :down 3))
  (testing "throws an error if inputs are illegal"
    (are [p s i] (thrown? #?(:clj Error, :cljs js/Error) (coord/cdna-coordinate p s i))
      0 nil nil
      3.5 nil nil
      "3" nil nil
      nil nil nil
      3 :illegal nil
      3 nil 8.5
      3 nil "3")))

(deftest format-cdna-coordinate-test
  (testing "returns a string expression of a CDNA coordinate"
    (are [m s] (= (coord/format m) s)
      (coord/cdna-coordinate 3 nil nil) "3"
      (coord/cdna-coordinate 3 :up nil) "-3"
      (coord/cdna-coordinate 3 :down nil) "*3"
      (coord/cdna-coordinate 87 nil 3) "87+3"
      (coord/cdna-coordinate 88 nil -1) "88-1"
      (coord/cdna-coordinate 88 nil 0) "88"
      (coord/cdna-coordinate 85 :up 3) "-85+3"
      (coord/cdna-coordinate 37 :down 3) "*37+3")))

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
      3 :up nil
      3 :down nil
      87 nil 3
      88 nil -1
      85 :up 3
      37 :down 3))
  (testing "throws an error if inputs are illegal"
    (are [p s i] (thrown? #?(:clj Error, :cljs js/Error) (coord/rna-coordinate p s i))
      0 nil nil
      3.5 nil nil
      "3" nil nil
      nil nil nil
      3 :illegal nil
      3 nil 8.5
      3 nil "3")))

(deftest protein-coordinate-test
  (testing "validates an input and returns ProteinCoordinate"
    (is (= (coord/protein-coordinate 3) (coord/->ProteinCoordinate 3))))
  (testing "throws an error if an input is illegal"
    (are [p] (thrown? #?(:clj Error, :cljs js/Error) (coord/protein-coordinate p))
      0
      3.5
      "3"
      nil)))
