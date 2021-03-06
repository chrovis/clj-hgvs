(ns clj-hgvs.coordinate-test
  (:require [clojure.test :refer [are deftest is testing]]
            [clj-hgvs.coordinate :as coord]
            clj-hgvs.test-common))

;;; uncertain coordinate

(deftest uncertain-coordinate-test
  (testing "validates an input and returns UncertainCoordinate"
    (are [s e] (= (coord/uncertain-coordinate s e) (coord/->UncertainCoordinate s e))
      (coord/genomic-coordinate 123456) (coord/genomic-coordinate 234567)
      (coord/unknown-coordinate) (coord/genomic-coordinate 1)))
  (testing "throws an error if an input is illegal"
    (are [s e] (thrown? #?(:clj Throwable, :cljs js/Error)
                        (coord/uncertain-coordinate s e))
      "123456" "234567"
      (coord/genomic-coordinate 123456) (coord/protein-coordinate 234567))))

(deftest parse-uncertain-coordinate-test
  (testing "parses input string, returning UncertainCoordinate"
    (are [s t e] (= (coord/parse-uncertain-coordinate s t) e)
      "(123456_234567)" :genomic
      (coord/uncertain-coordinate (coord/genomic-coordinate 123456)
                                  (coord/genomic-coordinate 234567))

      "(?_1)" :genomic
      (coord/uncertain-coordinate (coord/unknown-coordinate)
                                  (coord/genomic-coordinate 1))

      "4072-?" :coding-dna
      (coord/uncertain-coordinate (coord/unknown-coordinate)
                                  (coord/coding-dna-coordinate 4072 -1 nil))

      "5154+?" :coding-dna
      (coord/uncertain-coordinate (coord/coding-dna-coordinate 5154 1 nil)
                                  (coord/unknown-coordinate))))
  (testing "throws an error if any input is illegal"
    (are [s t] (thrown? #?(:clj Throwable, :cljs js/Error)
                        (coord/parse-uncertain-coordinate s t))
      "(123456_234567" :genomic
      "4072-1" :coding-dna
      "(123456_234567)" :illegal)))

(deftest format-uncertain-coordinate-test
  (testing "returns a string expression of a uncertain coordinate"
    (are [m s] (= (coord/format m) s)
      (coord/uncertain-coordinate (coord/genomic-coordinate 123456)
                                  (coord/genomic-coordinate 234567))
      "(123456_234567)"

      (coord/uncertain-coordinate (coord/unknown-coordinate)
                                  (coord/genomic-coordinate 1))
      "(?_1)"

      (coord/uncertain-coordinate (coord/unknown-coordinate)
                                  (coord/coding-dna-coordinate 4072 -1 nil))
      "4072-?"

      (coord/uncertain-coordinate (coord/coding-dna-coordinate 5154 1 nil)
                                  (coord/unknown-coordinate))
      "5154+?")))

(deftest uncertain-coordinate-calc-test
  (is (thrown? #?(:clj IllegalArgumentException, :cljs js/Error)
               (coord/plus (coord/uncertain-coordinate (coord/genomic-coordinate 2)
                                                       (coord/genomic-coordinate 7))
                           2)))
  (is (thrown? #?(:clj IllegalArgumentException, :cljs js/Error)
               (coord/minus (coord/uncertain-coordinate (coord/genomic-coordinate 2)
                                                        (coord/genomic-coordinate 7))
                            2))))

(deftest plain-uncertain-coordinate-test
  (testing "returns a plain map representing UncertainCoordinate"
    (is (= (coord/plain (coord/uncertain-coordinate (coord/unknown-coordinate)
                                                    (coord/genomic-coordinate 1)))
           {:coordinate "uncertain"
            :start {:coordinate "unknown"}
            :end {:coordinate "genomic"
                  :position 1}}))))

(deftest restore-uncertain-coordinate-test
  (testing "restores a plain map to UncertainCoordinate"
    (is (= (coord/restore {:coordinate "uncertain"
                           :start {:coordinate "unknown"}
                           :end {:coordinate "genomic"
                                 :position 1}})
           (coord/uncertain-coordinate (coord/unknown-coordinate)
                                       (coord/genomic-coordinate 1))))))

;;; unknown coordinate

(deftest format-unknown-coordinate-test
  (testing "returns a string expression of UnknownCoordinate"
    (is (= (coord/format (coord/unknown-coordinate)) "?"))))

(deftest unknown-coordinate-calc-test
  (is (thrown? #?(:clj IllegalArgumentException, :cljs js/Error)
               (coord/plus (coord/unknown-coordinate) 2)))
  (is (thrown? #?(:clj IllegalArgumentException, :cljs js/Error)
               (coord/minus (coord/unknown-coordinate) 2))))

(deftest plain-unknown-coordinate-test
  (testing "returns a plain map representing UnknownCoordinate"
    (is (= (coord/plain (coord/unknown-coordinate)) {:coordinate "unknown"}))))

(deftest restore-unknown-coordinate-test
  (testing "restores a plain map to UnknownCoordinate"
    (is (= (coord/restore {:coordinate "unknown"}) (coord/unknown-coordinate)))))

;;; genomic coordinate

(deftest genomic-coordinate-test
  (testing "validates an input and returns GenomicCoordinate"
    (is (= (coord/genomic-coordinate 3) (coord/->GenomicCoordinate 3))))
  (testing "throws an error if an input is illegal"
    (are [p] (thrown? #?(:clj Throwable, :cljs js/Error) (coord/genomic-coordinate p))
      0
      3.5
      "3"
      nil)))

(deftest genomic-compare-test
  (testing "compares actual positions"
    (are [c1 c2 e] (= (compare c1 c2) e)
      (coord/genomic-coordinate 2) (coord/genomic-coordinate 2) 0
      (coord/genomic-coordinate 3) (coord/genomic-coordinate 2) 1
      (coord/genomic-coordinate 2) (coord/genomic-coordinate 3) -1)))

(deftest parse-genomic-coordinate-test
  (testing "parses input string, returning GenomicCoordinate"
    (is (= (coord/parse-genomic-coordinate "3") (coord/genomic-coordinate 3))))
  (testing "returns UnknownCoordinate if input is \"?\""
    (is (= (coord/parse-genomic-coordinate "?") (coord/unknown-coordinate)))))

(deftest format-genomic-coordinate-test
  (testing "returns a string expression of GenomicCoordinate"
    (is (= (coord/format (coord/genomic-coordinate 3)) "3"))))

(deftest genomic-coordinate-calc-test
  (are [c n e] (= (coord/plus c n) (coord/minus c (- n)) e)
    (coord/genomic-coordinate 2) 1  (coord/genomic-coordinate 3)
    (coord/genomic-coordinate 2) -1 (coord/genomic-coordinate 1)
    (coord/genomic-coordinate 2) 0  (coord/genomic-coordinate 2))
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/plus c n))
    (coord/genomic-coordinate 2) -2
    (coord/genomic-coordinate 2) -3)
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/minus c n))
    (coord/genomic-coordinate 2) 2
    (coord/genomic-coordinate 2) 3))

(deftest plain-genomic-coordinate-test
  (testing "returns a plain map representing GenomicCoordinate"
    (is (= (coord/plain (coord/genomic-coordinate 3))
           {:coordinate "genomic", :position 3}))))

(deftest restore-genomic-coordinate-test
  (testing "restores a plain map to GenomicCoordinate"
    (is (= (coord/restore {:coordinate "genomic", :position 3})
           (coord/genomic-coordinate 3)))))

;;; mitochondrial coordinate

(deftest mitochondrial-coordinate-test
  (testing "validates an input and returns MitochondrialCoordinate"
    (is (= (coord/mitochondrial-coordinate 3) (coord/->MitochondrialCoordinate 3))))
  (testing "throws an error if an input is illegal"
    (are [p] (thrown? #?(:clj Throwable, :cljs js/Error) (coord/mitochondrial-coordinate p))
      0
      3.5
      "3"
      nil)))

(deftest parse-mitochondrial-coordinate-test
  (testing "parses input string, returning MitochondrialCoordinate"
    (is (= (coord/parse-mitochondrial-coordinate "3") (coord/mitochondrial-coordinate 3))))
  (testing "returns UnknownCoordinate if input is \"?\""
    (is (= (coord/parse-mitochondrial-coordinate "?") (coord/unknown-coordinate)))))

(deftest format-mitochondrial-coordinate-test
  (testing "returns a string expression of MitochondrialCoordinate"
    (is (= (coord/format (coord/mitochondrial-coordinate 3)) "3"))))

(deftest mitochondrial-coordinate-calc-test
  (are [c n e] (= (coord/plus c n) (coord/minus c (- n)) e)
    (coord/mitochondrial-coordinate 2) 1  (coord/mitochondrial-coordinate 3)
    (coord/mitochondrial-coordinate 2) -1 (coord/mitochondrial-coordinate 1)
    (coord/mitochondrial-coordinate 2) 0  (coord/mitochondrial-coordinate 2))
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/plus c n))
    (coord/mitochondrial-coordinate 2) -2
    (coord/mitochondrial-coordinate 2) -3)
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/minus c n))
    (coord/mitochondrial-coordinate 2) 2
    (coord/mitochondrial-coordinate 2) 3))

(deftest plain-mitochondrial-coordinate-test
  (testing "returns a plain map representing MitochondrialCoordinate"
    (is (= (coord/plain (coord/mitochondrial-coordinate 3))
           {:coordinate "mitochondrial", :position 3}))))

(deftest restore-mitochondrial-coordinate-test
  (testing "restores a plain map to MitochondrialCoordinate"
    (is (= (coord/restore {:coordinate "mitochondrial", :position 3})
           (coord/mitochondrial-coordinate 3)))))

;;; coding DNA coordinate

(deftest coding-dna-coordinate-test
  (testing "validates inputs and returns CodingDNACoordinate"
    (are [p o r] (= (coord/coding-dna-coordinate p o r) (coord/->CodingDNACoordinate p o r))
      3 0 nil
      3 0 :upstream
      3 0 :downstream
      87 3 nil
      88 -1 nil
      85 3 :upstream
      37 3 :downstream))
  (testing "throws an error if inputs are illegal"
    (are [p o r] (thrown? #?(:clj Throwable, :cljs js/Error) (coord/coding-dna-coordinate p o r))
      0 0 nil
      3.5 0 nil
      "3" 0 nil
      nil 0 nil
      3 0 :illegal
      3 8.5 nil
      3 "3" nil)))

(deftest codina-dna-compare-test
  (testing "compares actual positions"
    (are [c1 c2 e] (= (compare c1 c2) e)
      (coord/coding-dna-coordinate 2) (coord/coding-dna-coordinate 2) 0
      (coord/coding-dna-coordinate 3) (coord/coding-dna-coordinate 2) 1
      (coord/coding-dna-coordinate 2) (coord/coding-dna-coordinate 3) -1
      (coord/coding-dna-coordinate 2 0 :upstream) (coord/coding-dna-coordinate 3 0 :upstream) 1
      (coord/coding-dna-coordinate 3 0 :upstream) (coord/coding-dna-coordinate 2 0 :upstream) -1
      (coord/coding-dna-coordinate 2 3 nil) (coord/coding-dna-coordinate 2 1 nil) 1
      (coord/coding-dna-coordinate 2 1 nil) (coord/coding-dna-coordinate 2 3 nil) -1
      (coord/coding-dna-coordinate 2 -1 nil) (coord/coding-dna-coordinate 2 -3 nil) 1
      (coord/coding-dna-coordinate 2 -3 nil) (coord/coding-dna-coordinate 2 -1 nil) -1
      (coord/coding-dna-coordinate 2 3 nil) (coord/coding-dna-coordinate 3 -1 nil) -1
      (coord/coding-dna-coordinate 3 0 :upstream) (coord/coding-dna-coordinate 3 0 nil) -1
      (coord/coding-dna-coordinate 3 0 :downstream) (coord/coding-dna-coordinate 3 0 nil) 1)))

(deftest parse-coding-dna-coordinate-test
  (testing "parses input string, returning CodingDNACoordinate"
    (are [s e] (= (coord/parse-coding-dna-coordinate s) e)
      "3" (coord/coding-dna-coordinate 3 0 nil)
      "-3" (coord/coding-dna-coordinate 3 0 :upstream)
      "*3" (coord/coding-dna-coordinate 3 0 :downstream)
      "87+3" (coord/coding-dna-coordinate 87 3 nil)
      "88-1" (coord/coding-dna-coordinate 88 -1 nil)
      "88" (coord/coding-dna-coordinate 88 0 nil)
      "-85+3" (coord/coding-dna-coordinate 85 3 :upstream)
      "*37+3" (coord/coding-dna-coordinate 37 3 :downstream)))
  (testing "returns UnknownCoordinate if input is \"?\""
    (is (= (coord/parse-coding-dna-coordinate "?") (coord/unknown-coordinate)))))

(deftest format-coding-dna-coordinate-test
  (testing "returns a string expression of a coding DNA coordinate"
    (are [m s] (= (coord/format m) s)
      (coord/coding-dna-coordinate 3 0 nil) "3"
      (coord/coding-dna-coordinate 3 0 :upstream) "-3"
      (coord/coding-dna-coordinate 3 0 :downstream) "*3"
      (coord/coding-dna-coordinate 87 3 nil) "87+3"
      (coord/coding-dna-coordinate 88 -1 nil) "88-1"
      (coord/coding-dna-coordinate 88 0 nil) "88"
      (coord/coding-dna-coordinate 85 3 :upstream) "-85+3"
      (coord/coding-dna-coordinate 37 3 :downstream) "*37+3")))

(deftest coding-dna-in-exon?-test
  (testing "detects a coordinate is in exon ranges"
    (is (true? (coord/in-exon? (coord/coding-dna-coordinate 3 0 nil))))
    (is (false? (coord/in-exon? (coord/coding-dna-coordinate 87 3 nil))))))

(deftest coding-dna-coordinate-calc-test
  (are [c n e] (= (coord/plus c n) (coord/minus c (- n)) e)
    (coord/coding-dna-coordinate 3 0 nil) 0  (coord/coding-dna-coordinate 3 0 nil)
    (coord/coding-dna-coordinate 3 0 nil) 2  (coord/coding-dna-coordinate 5 0 nil)
    (coord/coding-dna-coordinate 3 0 nil) -2 (coord/coding-dna-coordinate 1 0 nil)
    (coord/coding-dna-coordinate 3 0 nil) -3 (coord/coding-dna-coordinate 1 0 :upstream)
    (coord/coding-dna-coordinate 3 0 nil) -5 (coord/coding-dna-coordinate 3 0 :upstream)
    ;; upstream
    (coord/coding-dna-coordinate 3 0 :upstream) 0  (coord/coding-dna-coordinate 3 0 :upstream)
    (coord/coding-dna-coordinate 3 0 :upstream) 2  (coord/coding-dna-coordinate 1 0 :upstream)
    (coord/coding-dna-coordinate 3 0 :upstream) -2 (coord/coding-dna-coordinate 5 0 :upstream)
    (coord/coding-dna-coordinate 3 0 :upstream) 3  (coord/coding-dna-coordinate 1 0 nil)
    (coord/coding-dna-coordinate 3 0 :upstream) 5  (coord/coding-dna-coordinate 3 0 nil)
    ;; downstream
    (coord/coding-dna-coordinate 3 0 :downstream) 0  (coord/coding-dna-coordinate 3 0 :downstream)
    (coord/coding-dna-coordinate 3 0 :downstream) 2  (coord/coding-dna-coordinate 5 0 :downstream)
    (coord/coding-dna-coordinate 3 0 :downstream) -2 (coord/coding-dna-coordinate 1 0 :downstream)
    ;; positive offset
    (coord/coding-dna-coordinate 3 2 nil) 0  (coord/coding-dna-coordinate 3 2 nil)
    (coord/coding-dna-coordinate 3 2 nil) 1  (coord/coding-dna-coordinate 3 3 nil)
    (coord/coding-dna-coordinate 3 2 nil) -1 (coord/coding-dna-coordinate 3 1 nil)
    (coord/coding-dna-coordinate 3 2 nil) -2 (coord/coding-dna-coordinate 3 0 nil)
    (coord/coding-dna-coordinate 3 2 nil) -4 (coord/coding-dna-coordinate 1 0 nil)
    (coord/coding-dna-coordinate 3 2 nil) -5 (coord/coding-dna-coordinate 1 0 :upstream)
    ;; positive offset and upstream
    (coord/coding-dna-coordinate 3 2 :upstream) 2  (coord/coding-dna-coordinate 3 4 :upstream)
    (coord/coding-dna-coordinate 3 2 :upstream) -2 (coord/coding-dna-coordinate 3 0 :upstream)
    (coord/coding-dna-coordinate 3 2 :upstream) -4 (coord/coding-dna-coordinate 5 0 :upstream)
    ;; positive offset and downstream
    (coord/coding-dna-coordinate 3 2 :downstream) 2  (coord/coding-dna-coordinate 3 4 :downstream)
    (coord/coding-dna-coordinate 3 2 :downstream) -2 (coord/coding-dna-coordinate 3 0 :downstream)
    (coord/coding-dna-coordinate 3 2 :downstream) -4 (coord/coding-dna-coordinate 1 0 :downstream)
    ;; negative offset
    (coord/coding-dna-coordinate 3 -2 nil) 0  (coord/coding-dna-coordinate 3 -2 nil)
    (coord/coding-dna-coordinate 3 -2 nil) 1  (coord/coding-dna-coordinate 3 -1 nil)
    (coord/coding-dna-coordinate 3 -2 nil) -1 (coord/coding-dna-coordinate 3 -3 nil)
    (coord/coding-dna-coordinate 3 -2 nil) 2  (coord/coding-dna-coordinate 3 0 nil)
    (coord/coding-dna-coordinate 3 -2 nil) 4  (coord/coding-dna-coordinate 5 0 nil)
    ;; negative offset and upstream
    (coord/coding-dna-coordinate 3 -2 :upstream) 2  (coord/coding-dna-coordinate 3 0 :upstream)
    (coord/coding-dna-coordinate 3 -2 :upstream) -2 (coord/coding-dna-coordinate 3 -4 :upstream)
    (coord/coding-dna-coordinate 3 -2 :upstream) 4  (coord/coding-dna-coordinate 1 0 :upstream)
    ;; negative offset and downstream
    (coord/coding-dna-coordinate 3 -2 :downstream) 2  (coord/coding-dna-coordinate 3 0 :downstream)
    (coord/coding-dna-coordinate 3 -2 :downstream) -2 (coord/coding-dna-coordinate 3 -4 :downstream)
    (coord/coding-dna-coordinate 3 -2 :downstream) 4  (coord/coding-dna-coordinate 5 0 :downstream))
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/plus c n))
    (coord/coding-dna-coordinate 3 0 :downstream) -3
    (coord/coding-dna-coordinate 3 0 :downstream) -5)
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/minus c n))
    (coord/coding-dna-coordinate 3 0 :downstream) 3
    (coord/coding-dna-coordinate 3 0 :downstream) 5))

(deftest plain-coding-dna-coordinate-test
  (testing "returns a plain map representing CodingDNACoordinate"
    (is (= (coord/plain (coord/coding-dna-coordinate 3 0 nil))
           {:coordinate "coding-dna", :position 3, :offset 0, :region nil}))))

(deftest restore-coding-dna-coordinate-test
  (testing "restores a plain map to CodingDNACoordinate"
    (is (= (coord/restore {:coordinate "coding-dna", :position 3, :offset 0, :region nil})
           (coord/coding-dna-coordinate 3 0 nil)))))

;;; non-coding DNA coordinate

(deftest non-coding-dna-coordinate-test
  (testing "validates an input and returns NonCodingDNACoordinate"
    (is (= (coord/non-coding-dna-coordinate 3) (coord/->NonCodingDNACoordinate 3))))
  (testing "throws an error if an input is illegal"
    (are [p] (thrown? #?(:clj Throwable, :cljs js/Error) (coord/non-coding-dna-coordinate p))
      0
      3.5
      "3"
      nil)))

(deftest non-coding-dna-compare-test
  (testing "compares actual positions"
    (are [c1 c2 e] (= (compare c1 c2) e)
      (coord/non-coding-dna-coordinate 2) (coord/non-coding-dna-coordinate 2) 0
      (coord/non-coding-dna-coordinate 3) (coord/non-coding-dna-coordinate 2) 1
      (coord/non-coding-dna-coordinate 2) (coord/non-coding-dna-coordinate 3) -1)))

(deftest parse-non-coding-dna-coordinate-test
  (testing "parses input string, returning NonCodingDNACoordinate"
    (is (= (coord/parse-non-coding-dna-coordinate "3") (coord/non-coding-dna-coordinate 3))))
  (testing "returns UnknownCoordinate if input is \"?\""
    (is (= (coord/parse-non-coding-dna-coordinate "?") (coord/unknown-coordinate)))))

(deftest format-non-coding-dna-coordinate-test
  (testing "returns a string expression of NonCodingDNACoordinate"
    (is (= (coord/format (coord/non-coding-dna-coordinate 3)) "3"))))

(deftest non-coding-dna-coordinate-calc-test
  (are [c n e] (= (coord/plus c n) (coord/minus c (- n)) e)
    (coord/non-coding-dna-coordinate 2) 1  (coord/non-coding-dna-coordinate 3)
    (coord/non-coding-dna-coordinate 2) -1 (coord/non-coding-dna-coordinate 1)
    (coord/non-coding-dna-coordinate 2) 0  (coord/non-coding-dna-coordinate 2))
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/plus c n))
    (coord/non-coding-dna-coordinate 2) -2
    (coord/non-coding-dna-coordinate 2) -3)
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/minus c n))
    (coord/non-coding-dna-coordinate 2) 2
    (coord/non-coding-dna-coordinate 2) 3))

(deftest plain-non-coding-dna-coordinate-test
  (testing "returns a plain map representing NonCodingDNACoordinate"
    (is (= (coord/plain (coord/non-coding-dna-coordinate 3))
           {:coordinate "non-coding-dna", :position 3}))))

(deftest restore-non-coding-dna-coordinate-test
  (testing "restores a plain map to NonCodingDNACoordinate"
    (is (= (coord/restore {:coordinate "non-coding-dna", :position 3})
           (coord/non-coding-dna-coordinate 3)))))

;;; circular DNA coordinate

(deftest circular-dna-coordinate-test
  (testing "validates an input and returns CircularDNACoordinate"
    (is (= (coord/circular-dna-coordinate 3) (coord/->CircularDNACoordinate 3))))
  (testing "throws an error if an input is illegal"
    (are [p] (thrown? #?(:clj Throwable, :cljs js/Error) (coord/circular-dna-coordinate p))
      0
      3.5
      "3"
      nil)))

(deftest parse-circular-dna-coordinate-test
  (testing "parses input string, returning CircularDNACoordinate"
    (is (= (coord/parse-circular-dna-coordinate "3") (coord/circular-dna-coordinate 3))))
  (testing "returns UnknownCoordinate if input is \"?\""
    (is (= (coord/parse-circular-dna-coordinate "?") (coord/unknown-coordinate)))))

(deftest format-circular-dna-coordinate-test
  (testing "returns a string expression of CircularDNACoordinate"
    (is (= (coord/format (coord/circular-dna-coordinate 3)) "3"))))

(deftest circular-dna-coordinate-calc-test
  (are [c n e] (= (coord/plus c n) (coord/minus c (- n)) e)
    (coord/circular-dna-coordinate 2) 1  (coord/circular-dna-coordinate 3)
    (coord/circular-dna-coordinate 2) -1 (coord/circular-dna-coordinate 1)
    (coord/circular-dna-coordinate 2) 0  (coord/circular-dna-coordinate 2))
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/plus c n))
    (coord/circular-dna-coordinate 2) -2
    (coord/circular-dna-coordinate 2) -3)
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/minus c n))
    (coord/circular-dna-coordinate 2) 2
    (coord/circular-dna-coordinate 2) 3))

(deftest plain-circular-dna-coordinate-test
  (testing "returns a plain map representing CircularDNACoordinate"
    (is (= (coord/plain (coord/circular-dna-coordinate 3))
           {:coordinate "circular-dna", :position 3}))))

(deftest restore-circular-dna-coordinate-test
  (testing "restores a plain map to CircularDNACoordinate"
    (is (= (coord/restore {:coordinate "circular-dna", :position 3})
           (coord/circular-dna-coordinate 3)))))

;;; RNA coordinate

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
    (are [p s i] (thrown? #?(:clj Throwable, :cljs js/Error) (coord/rna-coordinate p s i))
      0 nil nil
      3.5 nil nil
      "3" nil nil
      nil nil nil
      3 nil :illegal
      3 8.5 nil
      3 "3" nil)))

(deftest rna-compare-test
  (testing "compares actual positions"
    (are [c1 c2 e] (= (compare c1 c2) e)
      (coord/rna-coordinate 2 nil nil) (coord/rna-coordinate 2 nil nil) 0
      (coord/rna-coordinate 3 nil nil) (coord/rna-coordinate 2 nil nil) 1
      (coord/rna-coordinate 2 nil nil) (coord/rna-coordinate 3 nil nil) -1
      (coord/rna-coordinate 2 0 :upstream) (coord/rna-coordinate 3 0 :upstream) 1
      (coord/rna-coordinate 3 0 :upstream) (coord/rna-coordinate 2 0 :upstream) -1
      (coord/rna-coordinate 2 3 nil) (coord/rna-coordinate 2 1 nil) 1
      (coord/rna-coordinate 2 1 nil) (coord/rna-coordinate 2 3 nil) -1
      (coord/rna-coordinate 2 -1 nil) (coord/rna-coordinate 2 -3 nil) 1
      (coord/rna-coordinate 2 -3 nil) (coord/rna-coordinate 2 -1 nil) -1
      (coord/rna-coordinate 2 3 nil) (coord/rna-coordinate 3 -1 nil) -1
      (coord/rna-coordinate 3 0 :upstream) (coord/rna-coordinate 3 0 nil) -1
      (coord/rna-coordinate 3 0 :downstream) (coord/rna-coordinate 3 0 nil) 1)))

(deftest parse-rna-coordinate-test
  (testing "parses input string, returning RNACoordinate"
    (are [s e] (= (coord/parse-rna-coordinate s) e)
      "3" (coord/rna-coordinate 3 nil nil)
      "-3" (coord/rna-coordinate 3 0 :upstream)
      "*3" (coord/rna-coordinate 3 0 :downstream)
      "87+3" (coord/rna-coordinate 87 3 nil)
      "88-1" (coord/rna-coordinate 88 -1 nil)
      "88" (coord/rna-coordinate 88 nil nil)
      "-85+3" (coord/rna-coordinate 85 3 :upstream)
      "*37+3" (coord/rna-coordinate 37 3 :downstream)))
  (testing "returns UnknownCoordinate if input is \"?\""
    (is (= (coord/parse-rna-coordinate "?") (coord/unknown-coordinate)))))

(deftest format-rna-coordinate-test
  (testing "returns a string expression of a RNA coordinate"
    (are [m s] (= (coord/format m) s)
      (coord/rna-coordinate 3 0 nil) "3"
      (coord/rna-coordinate 3 0 :upstream) "-3"
      (coord/rna-coordinate 3 0 :downstream) "*3"
      (coord/rna-coordinate 87 3 nil) "87+3"
      (coord/rna-coordinate 88 -1 nil) "88-1"
      (coord/rna-coordinate 88 0 nil) "88"
      (coord/rna-coordinate 85 3 :upstream) "-85+3"
      (coord/rna-coordinate 37 3 :downstream) "*37+3")))

(deftest rna-in-exon?-test
  (testing "detects a coordinate is in exon ranges"
    (is (true? (coord/in-exon? (coord/rna-coordinate 3 0 nil))))
    (is (true? (coord/in-exon? (coord/rna-coordinate 3 nil nil))))
    (is (false? (coord/in-exon? (coord/rna-coordinate 87 3 nil))))))

(deftest rna-coordinate-calc-test
  (are [c n e] (= (coord/plus c n) (coord/minus c (- n)) e)
    (coord/rna-coordinate 3 0 nil) 0  (coord/rna-coordinate 3 0 nil)
    (coord/rna-coordinate 3 0 nil) 2  (coord/rna-coordinate 5 0 nil)
    (coord/rna-coordinate 3 0 nil) -2 (coord/rna-coordinate 1 0 nil)
    (coord/rna-coordinate 3 0 nil) -3 (coord/rna-coordinate 1 0 :upstream)
    (coord/rna-coordinate 3 0 nil) -5 (coord/rna-coordinate 3 0 :upstream)
    ;; upstream
    (coord/rna-coordinate 3 0 :upstream) 0  (coord/rna-coordinate 3 0 :upstream)
    (coord/rna-coordinate 3 0 :upstream) 2  (coord/rna-coordinate 1 0 :upstream)
    (coord/rna-coordinate 3 0 :upstream) -2 (coord/rna-coordinate 5 0 :upstream)
    (coord/rna-coordinate 3 0 :upstream) 3  (coord/rna-coordinate 1 0 nil)
    (coord/rna-coordinate 3 0 :upstream) 5  (coord/rna-coordinate 3 0 nil)
    ;; downstream
    (coord/rna-coordinate 3 0 :downstream) 0  (coord/rna-coordinate 3 0 :downstream)
    (coord/rna-coordinate 3 0 :downstream) 2  (coord/rna-coordinate 5 0 :downstream)
    (coord/rna-coordinate 3 0 :downstream) -2 (coord/rna-coordinate 1 0 :downstream)
    ;; positive offset
    (coord/rna-coordinate 3 2 nil) 0  (coord/rna-coordinate 3 2 nil)
    (coord/rna-coordinate 3 2 nil) 1  (coord/rna-coordinate 3 3 nil)
    (coord/rna-coordinate 3 2 nil) -1 (coord/rna-coordinate 3 1 nil)
    (coord/rna-coordinate 3 2 nil) -2 (coord/rna-coordinate 3 0 nil)
    (coord/rna-coordinate 3 2 nil) -4 (coord/rna-coordinate 1 0 nil)
    (coord/rna-coordinate 3 2 nil) -5 (coord/rna-coordinate 1 0 :upstream)
    ;; positive offset and upstream
    (coord/rna-coordinate 3 2 :upstream) 2  (coord/rna-coordinate 3 4 :upstream)
    (coord/rna-coordinate 3 2 :upstream) -2 (coord/rna-coordinate 3 0 :upstream)
    (coord/rna-coordinate 3 2 :upstream) -4 (coord/rna-coordinate 5 0 :upstream)
    ;; positive offset and downstream
    (coord/rna-coordinate 3 2 :downstream) 2  (coord/rna-coordinate 3 4 :downstream)
    (coord/rna-coordinate 3 2 :downstream) -2 (coord/rna-coordinate 3 0 :downstream)
    (coord/rna-coordinate 3 2 :downstream) -4 (coord/rna-coordinate 1 0 :downstream)
    ;; negative offset
    (coord/rna-coordinate 3 -2 nil) 0  (coord/rna-coordinate 3 -2 nil)
    (coord/rna-coordinate 3 -2 nil) 1  (coord/rna-coordinate 3 -1 nil)
    (coord/rna-coordinate 3 -2 nil) -1 (coord/rna-coordinate 3 -3 nil)
    (coord/rna-coordinate 3 -2 nil) 2  (coord/rna-coordinate 3 0 nil)
    (coord/rna-coordinate 3 -2 nil) 4  (coord/rna-coordinate 5 0 nil)
    ;; negative offset and upstream
    (coord/rna-coordinate 3 -2 :upstream) 2  (coord/rna-coordinate 3 0 :upstream)
    (coord/rna-coordinate 3 -2 :upstream) -2 (coord/rna-coordinate 3 -4 :upstream)
    (coord/rna-coordinate 3 -2 :upstream) 4  (coord/rna-coordinate 1 0 :upstream)
    ;; negative offset and downstream
    (coord/rna-coordinate 3 -2 :downstream) 2  (coord/rna-coordinate 3 0 :downstream)
    (coord/rna-coordinate 3 -2 :downstream) -2 (coord/rna-coordinate 3 -4 :downstream)
    (coord/rna-coordinate 3 -2 :downstream) 4  (coord/rna-coordinate 5 0 :downstream))
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/plus c n))
    (coord/rna-coordinate 3 0 :downstream) -3
    (coord/rna-coordinate 3 0 :downstream) -5)
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/minus c n))
    (coord/rna-coordinate 3 0 :downstream) 3
    (coord/rna-coordinate 3 0 :downstream) 5))

(deftest plain-rna-coordinate-test
  (testing "returns a plain map representing RNACoordinate"
    (is (= (coord/plain (coord/rna-coordinate 3 0 nil))
           {:coordinate "rna", :position 3, :offset 0, :region nil}))))

(deftest restore-rna-coordinate-test
  (testing "restores a plain map to RNACoordinate"
    (is (= (coord/restore {:coordinate "rna", :position 3, :offset 0, :region nil})
           (coord/rna-coordinate 3 0 nil)))))

;;; protein coordinate

(deftest protein-coordinate-test
  (testing "validates an input and returns ProteinCoordinate"
    (is (= (coord/protein-coordinate 3) (coord/->ProteinCoordinate 3))))
  (testing "throws an error if an input is illegal"
    (are [p] (thrown? #?(:clj Throwable, :cljs js/Error) (coord/protein-coordinate p))
      0
      3.5
      "3"
      nil)))

(deftest protein-compare-test
  (testing "compares actual positions"
    (are [c1 c2 e] (= (compare c1 c2) e)
      (coord/protein-coordinate 2) (coord/protein-coordinate 2) 0
      (coord/protein-coordinate 3) (coord/protein-coordinate 2) 1
      (coord/protein-coordinate 2) (coord/protein-coordinate 3) -1)))

(deftest parse-protein-coordinate-test
  (testing "parses input string, returning ProteinCoordinate"
    (is (= (coord/parse-protein-coordinate "3") (coord/protein-coordinate 3))))
  (testing "returns UnknownCoordinate if input is \"?\""
    (is (= (coord/parse-protein-coordinate "?") (coord/unknown-coordinate)))))

(deftest format-protein-coordinate-test
  (testing "returns a string expression of ProteinCoordinate"
    (is (= (coord/format (coord/protein-coordinate 3)) "3"))))

(deftest protein-coordinate-calc-test
  (are [c n e] (= (coord/plus c n) (coord/minus c (- n)) e)
    (coord/protein-coordinate 2) 1  (coord/protein-coordinate 3)
    (coord/protein-coordinate 2) -1 (coord/protein-coordinate 1)
    (coord/protein-coordinate 2) 0  (coord/protein-coordinate 2))
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/plus c n))
    (coord/protein-coordinate 2) -2
    (coord/protein-coordinate 2) -3)
  (are [c n] (thrown? #?(:clj ArithmeticException, :cljs js/RangeError)
                      (coord/minus c n))
    (coord/protein-coordinate 2) 2
    (coord/protein-coordinate 2) 3))

(deftest plain-protein-coordinate-test
  (testing "returns a plain map representing ProteinCoordinate"
    (is (= (coord/plain (coord/protein-coordinate 3))
           {:coordinate "protein", :position 3}))))

(deftest restore-protein-coordinate-test
  (testing "restores a plain map to ProteinCoordinate"
    (is (= (coord/restore {:coordinate "protein", :position 3})
           (coord/protein-coordinate 3)))))
