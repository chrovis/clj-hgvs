(ns clj-hgvs.coordinate-test
  (:require #?(:clj [clojure.test :refer :all]
               :cljs [cljs.test :refer-macros [deftest are is testing]])
            [clj-hgvs.coordinate :as coord]))

;;; uncertain coordinate

(deftest uncertain-coordinate-test
  (testing "validates an input and returns UncertainCoordinate"
    (are [s e] (= (coord/uncertain-coordinate s e) (coord/->UncertainCoordinate s e))
      (coord/genomic-coordinate 123456) (coord/genomic-coordinate 234567)
      (coord/unknown-coordinate) (coord/genomic-coordinate 1)))
  (testing "throws an error if an input is illegal"
    (are [s e] (thrown? #?(:clj Error, :cljs js/Error)
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
                                  (coord/genomic-coordinate 1))))
  (testing "throws an error if any input is illegal"
    (are [s t] (thrown? #?(:clj Throwable, :cljs js/Error)
                        (coord/parse-uncertain-coordinate s t))
      "(123456_234567" :genomic
      "(123456_234567)" :illegal)))

(deftest format-uncertain-coordinate-test
  (testing "returns a string expression of a uncertain coordinate"
    (are [m s] (= (coord/format m) s)
      (coord/uncertain-coordinate (coord/genomic-coordinate 123456)
                                  (coord/genomic-coordinate 234567))
      "(123456_234567)"

      (coord/uncertain-coordinate (coord/unknown-coordinate)
                                  (coord/genomic-coordinate 1))
      "(?_1)")))

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
    (are [p] (thrown? #?(:clj Error, :cljs js/Error) (coord/genomic-coordinate p))
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
    (are [p] (thrown? #?(:clj Error, :cljs js/Error) (coord/mitochondrial-coordinate p))
      0
      3.5
      "3"
      nil)))

(deftest mitochondrial-compare-test
  (testing "compares actual positions"
    (are [c1 c2 e] (= (compare c1 c2) e)
      (coord/mitochondrial-coordinate 2) (coord/mitochondrial-coordinate 2) 0
      (coord/mitochondrial-coordinate 3) (coord/mitochondrial-coordinate 2) 1
      (coord/mitochondrial-coordinate 2) (coord/mitochondrial-coordinate 3) -1)))

(deftest parse-mitochondrial-coordinate-test
  (testing "parses input string, returning MitochondrialCoordinate"
    (is (= (coord/parse-mitochondrial-coordinate "3") (coord/mitochondrial-coordinate 3))))
  (testing "returns UnknownCoordinate if input is \"?\""
    (is (= (coord/parse-mitochondrial-coordinate "?") (coord/unknown-coordinate)))))

(deftest format-mitochondrial-coordinate-test
  (testing "returns a string expression of MitochondrialCoordinate"
    (is (= (coord/format (coord/mitochondrial-coordinate 3)) "3"))))

(deftest plain-mitochondrial-coordinate-test
  (testing "returns a plain map representing MitochondrialCoordinate"
    (is (= (coord/plain (coord/mitochondrial-coordinate 3))
           {:coordinate "mitochondrial", :position 3}))))

(deftest restore-mitochondrial-coordinate-test
  (testing "restores a plain map to MitochondrialCoordinate"
    (is (= (coord/restore {:coordinate "mitochondrial", :position 3})
           (coord/mitochondrial-coordinate 3)))))

;;; coding DNA coordinate

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

(deftest cdna-compare-test
  (testing "compares actual positions"
    (are [c1 c2 e] (= (compare c1 c2) e)
      (coord/cdna-coordinate 2) (coord/cdna-coordinate 2) 0
      (coord/cdna-coordinate 3) (coord/cdna-coordinate 2) 1
      (coord/cdna-coordinate 2) (coord/cdna-coordinate 3) -1
      (coord/cdna-coordinate 2 0 :upstream) (coord/cdna-coordinate 3 0 :upstream) 1
      (coord/cdna-coordinate 3 0 :upstream) (coord/cdna-coordinate 2 0 :upstream) -1
      (coord/cdna-coordinate 2 3 nil) (coord/cdna-coordinate 2 1 nil) 1
      (coord/cdna-coordinate 2 1 nil) (coord/cdna-coordinate 2 3 nil) -1
      (coord/cdna-coordinate 2 -1 nil) (coord/cdna-coordinate 2 -3 nil) 1
      (coord/cdna-coordinate 2 -3 nil) (coord/cdna-coordinate 2 -1 nil) -1
      (coord/cdna-coordinate 2 3 nil) (coord/cdna-coordinate 3 -1 nil) -1
      (coord/cdna-coordinate 3 0 :upstream) (coord/cdna-coordinate 3 0 nil) -1
      (coord/cdna-coordinate 3 0 :downstream) (coord/cdna-coordinate 3 0 nil) 1)))

(deftest parse-cdna-coordinate-test
  (testing "parses input string, returning CDNACoordinate"
    (are [s e] (= (coord/parse-cdna-coordinate s) e)
      "3" (coord/cdna-coordinate 3 0 nil)
      "-3" (coord/cdna-coordinate 3 0 :upstream)
      "*3" (coord/cdna-coordinate 3 0 :downstream)
      "87+3" (coord/cdna-coordinate 87 3 nil)
      "88-1" (coord/cdna-coordinate 88 -1 nil)
      "88" (coord/cdna-coordinate 88 0 nil)
      "-85+3" (coord/cdna-coordinate 85 3 :upstream)
      "*37+3" (coord/cdna-coordinate 37 3 :downstream)))
  (testing "returns UnknownCoordinate if input is \"?\""
    (is (= (coord/parse-cdna-coordinate "?") (coord/unknown-coordinate)))))

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

(deftest cdna-in-exon?-test
  (testing "detects a coordinate is in exon ranges"
    (is (true? (coord/in-exon? (coord/cdna-coordinate 3 0 nil))))
    (is (false? (coord/in-exon? (coord/cdna-coordinate 87 3 nil))))))

(deftest plain-cdna-coordinate-test
  (testing "returns a plain map representing CDNACoordinate"
    (is (= (coord/plain (coord/cdna-coordinate 3 0 nil))
           {:coordinate "cdna", :position 3, :offset 0, :region nil}))))

(deftest restore-cdna-coordinate-test
  (testing "restores a plain map to CDNACoordinate"
    (is (= (coord/restore {:coordinate "cdna", :position 3, :offset 0, :region nil})
           (coord/cdna-coordinate 3 0 nil)))))

;;; non-coding DNA coordinate

(deftest ncdna-coordinate-test
  (testing "validates an input and returns NCDNACoordinate"
    (is (= (coord/ncdna-coordinate 3) (coord/->NCDNACoordinate 3))))
  (testing "throws an error if an input is illegal"
    (are [p] (thrown? #?(:clj Error, :cljs js/Error) (coord/ncdna-coordinate p))
      0
      3.5
      "3"
      nil)))

(deftest ncdna-compare-test
  (testing "compares actual positions"
    (are [c1 c2 e] (= (compare c1 c2) e)
      (coord/ncdna-coordinate 2) (coord/ncdna-coordinate 2) 0
      (coord/ncdna-coordinate 3) (coord/ncdna-coordinate 2) 1
      (coord/ncdna-coordinate 2) (coord/ncdna-coordinate 3) -1)))

(deftest parse-ncdna-coordinate-test
  (testing "parses input string, returning NCDNACoordinate"
    (is (= (coord/parse-ncdna-coordinate "3") (coord/ncdna-coordinate 3))))
  (testing "returns UnknownCoordinate if input is \"?\""
    (is (= (coord/parse-ncdna-coordinate "?") (coord/unknown-coordinate)))))

(deftest format-ncdna-coordinate-test
  (testing "returns a string expression of NCDNACoordinate"
    (is (= (coord/format (coord/ncdna-coordinate 3)) "3"))))

(deftest plain-ncdna-coordinate-test
  (testing "returns a plain map representing NCDNACoordinate"
    (is (= (coord/plain (coord/ncdna-coordinate 3))
           {:coordinate "ncdna", :position 3}))))

(deftest restore-ncdna-coordinate-test
  (testing "restores a plain map to NCDNACoordinate"
    (is (= (coord/restore {:coordinate "ncdna", :position 3})
           (coord/ncdna-coordinate 3)))))

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
    (are [p s i] (thrown? #?(:clj Error, :cljs js/Error) (coord/rna-coordinate p s i))
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
    (are [p] (thrown? #?(:clj Error, :cljs js/Error) (coord/protein-coordinate p))
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

(deftest plain-protein-coordinate-test
  (testing "returns a plain map representing ProteinCoordinate"
    (is (= (coord/plain (coord/protein-coordinate 3))
           {:coordinate "protein", :position 3}))))

(deftest restore-protein-coordinate-test
  (testing "restores a plain map to ProteinCoordinate"
    (is (= (coord/restore {:coordinate "protein", :position 3})
           (coord/protein-coordinate 3)))))
