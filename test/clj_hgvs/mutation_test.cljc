(ns clj-hgvs.mutation-test
  (:require [clojure.test :refer [are deftest is testing]]
            [clj-hgvs.coordinate :as coord]
            [clj-hgvs.mutation :as mut]
            clj-hgvs.test-common))

(deftest ->long-amino-acid-test
  (testing "converts a short amino acid to a long one"
    (is (= (mut/->long-amino-acid "S") "Ser"))
    (is (= (mut/->long-amino-acid \S) "Ser")))
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
    (is (= (mut/->short-amino-acid "S") "S"))
    (is (= (mut/->short-amino-acid \S) "S")))
  (testing "returns nil when an illegal string is passed"
    (is (nil? (mut/->short-amino-acid "Z")))
    (is (nil? (mut/->short-amino-acid "")))
    (is (nil? (mut/->short-amino-acid nil)))))

;;; Uncertain mutation

(def uncertain-mutation1s "(?)")
(def uncertain-mutation1k :rna)
(def uncertain-mutation1 (mut/uncertain-mutation (mut/rna-unknown-mutation)))
(def uncertain-mutation1m {:mutation "uncertain-mutation"
                           :content-mutation {:mutation "rna-unknown"}})

(def uncertain-mutation2s "(306g>u)")
(def uncertain-mutation2k :rna)
(def uncertain-mutation2 (mut/uncertain-mutation
                          (mut/rna-substitution (coord/rna-coordinate 306 nil nil) "g" "u")))

(def uncertain-mutation3s "(Cys123Gly)")
(def uncertain-mutation3ss "(C123G)")
(def uncertain-mutation3k :protein)
(def uncertain-mutation3 (mut/uncertain-mutation
                          (mut/protein-substitution "Cys"
                                                    (coord/protein-coordinate 123)
                                                    "Gly")))

(deftest format-uncertain-mutation-test
  (testing "returns a string expression of an uncertain mutation"
    (are [m o s] (= (mut/format m o) s)
      uncertain-mutation1 nil uncertain-mutation1s
      uncertain-mutation2 nil uncertain-mutation2s
      uncertain-mutation3 nil uncertain-mutation3s
      uncertain-mutation3 {:amino-acid-format :short} uncertain-mutation3ss)))

(deftest parse-uncertain-mutation-test
  (testing "returns a correct UncertainMutation"
    (are [s k m] (= (mut/parse-uncertain-mutation s k) m)
      uncertain-mutation1s uncertain-mutation1k uncertain-mutation1
      uncertain-mutation2s uncertain-mutation2k uncertain-mutation2
      uncertain-mutation3s uncertain-mutation3k uncertain-mutation3
      uncertain-mutation3ss uncertain-mutation3k uncertain-mutation3)))

(deftest plain-uncertain-mutation-test
  (testing "returns a plain map representing UncertainMutation"
    (is (= (mut/plain uncertain-mutation1) uncertain-mutation1m))))

(deftest restore-uncertain-mutation-test
  (testing "restores a plain map to UncertainMutation"
    (is (= (mut/restore uncertain-mutation1m) uncertain-mutation1))))

;;; DNA mutations

;;; DNA - substitution

(def dna-substitution1s "45576A>C")
(def dna-substitution1k :genome)
(def dna-substitution1 (mut/dna-substitution (coord/genomic-coordinate 45576) "A" ">" "C"))
(def dna-substitution1m {:mutation "dna-substitution"
                         :coord (coord/plain (coord/genomic-coordinate 45576))
                         :ref "A"
                         :type ">"
                         :alt "C"})

(def dna-substitution2s "88+1G>T")
(def dna-substitution2k :coding-dna)
(def dna-substitution2 (mut/dna-substitution (coord/coding-dna-coordinate 88 1 nil) "G" ">" "T"))

(def dna-substitution3s "123G=")
(def dna-substitution3k :coding-dna)
(def dna-substitution3 (mut/dna-substitution (coord/coding-dna-coordinate 123) "G" "="))

(def dna-substitution4s "85C=/>T")
(def dna-substitution4k :coding-dna)
(def dna-substitution4 (mut/dna-substitution (coord/coding-dna-coordinate 85) "C" "=/>" "T"))

(def dna-substitution5s "85C=//>T")
(def dna-substitution5k :coding-dna)
(def dna-substitution5 (mut/dna-substitution (coord/coding-dna-coordinate 85) "C" "=//>" "T"))

(deftest format-dna-substitution-test
  (testing "returns a string expression of a DNA substitution"
    (are [m s] (= (mut/format m nil) s)
      dna-substitution1 dna-substitution1s
      dna-substitution2 dna-substitution2s
      dna-substitution3 dna-substitution3s
      dna-substitution4 dna-substitution4s
      dna-substitution5 dna-substitution5s))
  (testing "alt is omitted if type is ="
    (is (= (mut/format (mut/dna-substitution (coord/coding-dna-coordinate 123) "G" "=" "G"))
           "123G="))))

(deftest parse-dna-substitution-test
  (testing "returns a correct DNASubstitution"
    (are [s k m] (= (mut/parse-dna-substitution s k) m)
      dna-substitution1s dna-substitution1k dna-substitution1
      dna-substitution2s dna-substitution2k dna-substitution2
      dna-substitution3s dna-substitution3k dna-substitution3
      dna-substitution4s dna-substitution4k dna-substitution4
      dna-substitution5s dna-substitution5k dna-substitution5)))

(deftest plain-dna-substitution-test
  (testing "returns a plain map representing DNASubstitution"
    (is (= (mut/plain dna-substitution1) dna-substitution1m))))

(deftest restore-dna-substitution-test
  (testing "restores a plain map to DNASubstitution"
    (is (= (mut/restore dna-substitution1m) dna-substitution1))))

;;; DNA - deletion

(def dna-deletion1s "7del")
(def dna-deletion1k :genome)
(def dna-deletion1 (mut/dna-deletion (coord/genomic-coordinate 7) nil))
(def dna-deletion1m {:mutation "dna-deletion"
                     :coord-start (coord/plain (coord/genomic-coordinate 7))
                     :coord-end nil
                     :ref nil})

(def dna-deletion2s "6_8del")
(def dna-deletion2k :genome)
(def dna-deletion2 (mut/dna-deletion (coord/genomic-coordinate 6)
                                     (coord/genomic-coordinate 8)))

(def dna-deletion3s "6_8delTGC")
(def dna-deletion3ss "6_8del")
(def dna-deletion3k :genome)
(def dna-deletion3 (mut/dna-deletion (coord/genomic-coordinate 6)
                                     (coord/genomic-coordinate 8)
                                     "TGC"))

(def dna-deletion4s "120_123+48del")
(def dna-deletion4k :coding-dna)
(def dna-deletion4 (mut/dna-deletion (coord/coding-dna-coordinate 120)
                                     (coord/coding-dna-coordinate 123 48 nil)))

(def dna-deletion5s "(4071+1_4072-1)_(5145+1_5146-1)del")
(def dna-deletion5k :coding-dna)
(def dna-deletion5 (mut/dna-deletion (coord/uncertain-coordinate
                                      (coord/coding-dna-coordinate 4071 1 nil)
                                      (coord/coding-dna-coordinate 4072 -1 nil))
                                     (coord/uncertain-coordinate
                                      (coord/coding-dna-coordinate 5145 1 nil)
                                      (coord/coding-dna-coordinate 5146 -1 nil))))

(def dna-deletion6s "(?_-30)_(12+1_13-1)del")
(def dna-deletion6k :coding-dna)
(def dna-deletion6 (mut/dna-deletion (coord/uncertain-coordinate
                                      (coord/unknown-coordinate)
                                      (coord/coding-dna-coordinate 30 0 :upstream))
                                     (coord/uncertain-coordinate
                                      (coord/coding-dna-coordinate 12 1 nil)
                                      (coord/coding-dna-coordinate 13 -1 nil))))

(def dna-deletion7s "(?_-1)_(*1_?)del")
(def dna-deletion7k :coding-dna)
(def dna-deletion7 (mut/dna-deletion (coord/uncertain-coordinate
                                      (coord/unknown-coordinate)
                                      (coord/coding-dna-coordinate 1 0 :upstream))
                                     (coord/uncertain-coordinate
                                      (coord/coding-dna-coordinate 1 0 :downstream)
                                      (coord/unknown-coordinate))))

(deftest format-dna-deletion-test
  (testing "returns a string expression of a DNA deletion"
    (are [m o s] (= (mut/format m o) s)
      dna-deletion1 nil dna-deletion1s
      dna-deletion2 nil dna-deletion2s
      dna-deletion3 nil dna-deletion3ss
      dna-deletion3 {:show-bases? true} dna-deletion3s
      dna-deletion4 nil dna-deletion4s
      dna-deletion5 nil dna-deletion5s
      dna-deletion6 nil dna-deletion6s
      dna-deletion7 nil dna-deletion7s)))

(deftest parse-dna-deletion-test
  (testing "returns a correct DNADeletion"
    (are [s k m] (= (mut/parse-dna-deletion s k) m)
      dna-deletion1s dna-deletion1k dna-deletion1
      dna-deletion2s dna-deletion2k dna-deletion2
      dna-deletion3s dna-deletion3k dna-deletion3
      dna-deletion4s dna-deletion4k dna-deletion4
      dna-deletion5s dna-deletion5k dna-deletion5
      dna-deletion6s dna-deletion6k dna-deletion6
      dna-deletion7s dna-deletion7k dna-deletion7))
  (testing "invalid DNA deletion"
    (are [s k] (thrown? #?(:clj Throwable, :cljs js/Error)
                        (mut/parse-dna-deletion s k))
      "8_6del" :genome)))

(deftest plain-dna-deletion-test
  (testing "returns a plain map representing DNADeletion"
    (is (= (mut/plain dna-deletion1) dna-deletion1m))))

(deftest restore-dna-deletion-test
  (testing "restores a plain map to DNADeletion"
    (is (= (mut/restore dna-deletion1m) dna-deletion1))))

;;; DNA - duplication

(def dna-duplication1s "7dup")
(def dna-duplication1k :genome)
(def dna-duplication1 (mut/dna-duplication (coord/genomic-coordinate 7) nil))
(def dna-duplication1m {:mutation "dna-duplication"
                        :coord-start (coord/plain (coord/genomic-coordinate 7))
                        :coord-end nil
                        :ref nil})

(def dna-duplication2s "6_8dup")
(def dna-duplication2k :genome)
(def dna-duplication2 (mut/dna-duplication (coord/genomic-coordinate 6)
                                           (coord/genomic-coordinate 8)))

(def dna-duplication3s "6_8dupTGC")
(def dna-duplication3ss "6_8dup")
(def dna-duplication3k :genome)
(def dna-duplication3 (mut/dna-duplication (coord/genomic-coordinate 6)
                                           (coord/genomic-coordinate 8)
                                           "TGC"))

(def dna-duplication4s "120_123+48dup")
(def dna-duplication4k :coding-dna)
(def dna-duplication4 (mut/dna-duplication (coord/coding-dna-coordinate 120)
                                           (coord/coding-dna-coordinate 123 48 nil)))

(def dna-duplication5s "(4071+1_4072-1)_(5145+1_5146-1)dup")
(def dna-duplication5k :coding-dna)
(def dna-duplication5 (mut/dna-duplication (coord/uncertain-coordinate
                                            (coord/coding-dna-coordinate 4071 1 nil)
                                            (coord/coding-dna-coordinate 4072 -1 nil))
                                           (coord/uncertain-coordinate
                                            (coord/coding-dna-coordinate 5145 1 nil)
                                            (coord/coding-dna-coordinate 5146 -1 nil))))

(def dna-duplication6s "(?_-30)_(12+1_13-1)dup")
(def dna-duplication6k :coding-dna)
(def dna-duplication6 (mut/dna-duplication (coord/uncertain-coordinate
                                            (coord/unknown-coordinate)
                                            (coord/coding-dna-coordinate 30 0 :upstream))
                                           (coord/uncertain-coordinate
                                            (coord/coding-dna-coordinate 12 1 nil)
                                            (coord/coding-dna-coordinate 13 -1 nil))))

(def dna-duplication7s "(?_-1)_(*1_?)dup")
(def dna-duplication7k :coding-dna)
(def dna-duplication7 (mut/dna-duplication (coord/uncertain-coordinate
                                            (coord/unknown-coordinate)
                                            (coord/coding-dna-coordinate 1 0 :upstream))
                                           (coord/uncertain-coordinate
                                            (coord/coding-dna-coordinate 1 0 :downstream)
                                            (coord/unknown-coordinate))))

(deftest format-dna-duplication-test
  (testing "returns a string expression of a DNA duplication"
    (are [m o s] (= (mut/format m o) s)
      dna-duplication1 nil dna-duplication1s
      dna-duplication2 nil dna-duplication2s
      dna-duplication3 nil dna-duplication3ss
      dna-duplication3 {:show-bases? true} dna-duplication3s
      dna-duplication4 nil dna-duplication4s
      dna-duplication5 nil dna-duplication5s
      dna-duplication6 nil dna-duplication6s
      dna-duplication7 nil dna-duplication7s)))

(deftest parse-dna-duplication-test
  (testing "returns a correct DNADuplication"
    (are [s k m] (= (mut/parse-dna-duplication s k) m)
      dna-duplication1s dna-duplication1k dna-duplication1
      dna-duplication2s dna-duplication2k dna-duplication2
      dna-duplication3s dna-duplication3k dna-duplication3
      dna-duplication4s dna-duplication4k dna-duplication4
      dna-duplication5s dna-duplication5k dna-duplication5
      dna-duplication6s dna-duplication6k dna-duplication6
      dna-duplication7s dna-duplication7k dna-duplication7))
  (testing "invalid DNA duplication"
    (are [s k] (thrown? #?(:clj Throwable, :cljs js/Error)
                        (mut/parse-dna-duplication s k))
      "8_6dup" :genome)))

(deftest plain-dna-duplication-test
  (testing "returns a plain map representing DNADuplication"
    (is (= (mut/plain dna-duplication1) dna-duplication1m))))

(deftest restore-dna-duplication-test
  (testing "restores a plain map to DNADuplication"
    (is (= (mut/restore dna-duplication1m) dna-duplication1))))

;;; DNA - insertion

(def dna-insertion1s "5756_5757insAGG")
(def dna-insertion1k :genome)
(def dna-insertion1 (mut/dna-insertion (coord/genomic-coordinate 5756)
                                       (coord/genomic-coordinate 5757)
                                       "AGG"))
(def dna-insertion1m {:mutation "dna-insertion"
                      :coord-start (coord/plain (coord/genomic-coordinate 5756))
                      :coord-end (coord/plain (coord/genomic-coordinate 5757))
                      :alt "AGG"})

(def dna-insertion2s "123_124insL37425.1:23_361")
(def dna-insertion2k :genome)
(def dna-insertion2 (mut/dna-insertion (coord/genomic-coordinate 123)
                                       (coord/genomic-coordinate 124)
                                       {:transcript "L37425.1"
                                        :coord-start (coord/genomic-coordinate 23)
                                        :coord-end (coord/genomic-coordinate 361)}))

(def dna-insertion3s "122_123ins123_234inv")
(def dna-insertion3k :genome)
(def dna-insertion3 "TODO")

(def dna-insertion4s "122_123ins213_234invinsAins123_211inv")
(def dna-insertion4k :genome)
(def dna-insertion4 "TODO")

(def dna-insertion5s "549_550insN")
(def dna-insertion5k :genome)
(def dna-insertion5 (mut/dna-insertion (coord/genomic-coordinate 549)
                                       (coord/genomic-coordinate 550)
                                       "N"))

(def dna-insertion6sb "1134_1135insNNNNNNNNNN")
(def dna-insertion6sc "1134_1135insN[10]")
(def dna-insertion6k :genome)
(def dna-insertion6 (mut/dna-insertion (coord/genomic-coordinate 1134)
                                       (coord/genomic-coordinate 1135)
                                       "NNNNNNNNNN"))

(def dna-insertion7s "?_?insNC_000023.10:(12345_23456)_(34567_45678)")
(def dna-insertion7k :genome)
(def dna-insertion7 (mut/dna-insertion (coord/unknown-coordinate)
                                       (coord/unknown-coordinate)
                                       {:transcript "NC_000023.10"
                                        :coord-start (coord/uncertain-coordinate
                                                      (coord/genomic-coordinate 12345)
                                                      (coord/genomic-coordinate 23456))
                                        :coord-end (coord/uncertain-coordinate
                                                    (coord/genomic-coordinate 34567)
                                                    (coord/genomic-coordinate 45678))}))

(deftest format-dna-insertion-test
  (testing "returns a string expression of a DNA insertion"
    (are [m o s] (= (mut/format m o) s)
      dna-insertion1 nil dna-insertion1s
      dna-insertion2 nil dna-insertion2s
      ;; dna-insertion3 nil dna-insertion3s ; TODO
      ;; dna-insertion4 nil dna-insertion4s ; TODO
      dna-insertion5 nil dna-insertion5s
      dna-insertion6 {:ins-format :auto} dna-insertion6sc
      dna-insertion6 {:ins-format :bases} dna-insertion6sb
      dna-insertion6 {:ins-format :count} dna-insertion6sc
      dna-insertion7 nil dna-insertion7s)))

(deftest parse-dna-insertion-test
  (testing "returns a correct DNAInsertion"
    (are [s k m] (= (mut/parse-dna-insertion s k) m)
      dna-insertion1s dna-insertion1k dna-insertion1
      dna-insertion2s dna-insertion2k dna-insertion2
      ;; dna-insertion3s dna-insertion3k dna-insertion3 ; TODO
      ;; dna-insertion4s dna-insertion4k dna-insertion4 ; TODO
      dna-insertion5s dna-insertion5k dna-insertion5
      dna-insertion6sb dna-insertion6k dna-insertion6
      dna-insertion6sc dna-insertion6k dna-insertion6
      dna-insertion7s dna-insertion7k dna-insertion7))
  (testing "invalid DNA insertion"
    (are [s k] (thrown? #?(:clj Throwable, :cljs js/Error)
                        (mut/parse-dna-insertion s k))
      "5756insAGG" :genome
      "5757_5756insAGG" :genome)))

(deftest plain-dna-insertion-test
  (testing "returns a plain map representing DNAInsertion"
    (is (= (mut/plain dna-insertion1) dna-insertion1m))))

(deftest restore-dna-insertion-test
  (testing "restores a plain map to DNAInsertion"
    (is (= (mut/restore dna-insertion1m) dna-insertion1))))

;;; DNA - inversion

(def dna-inversion1s "1077_1080inv")
(def dna-inversion1k :genome)
(def dna-inversion1 (mut/dna-inversion (coord/genomic-coordinate 1077)
                                       (coord/genomic-coordinate 1080)))

(def dna-inversion2s "77_80inv")
(def dna-inversion2k :coding-dna)
(def dna-inversion2 (mut/dna-inversion (coord/coding-dna-coordinate 77)
                                       (coord/coding-dna-coordinate 80)))

(def dna-inversion3s "77-?_80+?inv")
(def dna-inversion3k :coding-dna)
(def dna-inversion3 (mut/dna-inversion (coord/uncertain-coordinate
                                        (coord/unknown-coordinate)
                                        (coord/coding-dna-coordinate 77 -1 nil))
                                       (coord/uncertain-coordinate
                                        (coord/coding-dna-coordinate 80 1 nil)
                                        (coord/unknown-coordinate))))

(deftest format-dna-inversion-test
  (testing "returns a string expression of a DNA inversion"
    (are [m s] (= (mut/format m nil) s)
      dna-inversion1 dna-inversion1s
      dna-inversion2 dna-inversion2s
      dna-inversion3 dna-inversion3s)))

(deftest parse-dna-inversion-test
  (testing "returns a correct DNAInversion"
    (are [s k m] (= (mut/parse-dna-inversion s k) m)
      dna-inversion1s dna-inversion1k dna-inversion1
      dna-inversion2s dna-inversion2k dna-inversion2
      dna-inversion3s dna-inversion3k dna-inversion3))
  (testing "invalid DNA inversion"
    (are [s k] (thrown? #?(:clj Throwable, :cljs js/Error)
                        (mut/parse-dna-inversion s k))
      "1077inv" :genome
      "1080_1077inv" :genome)))

;;; DNA - conversion

(def dna-conversion1s "333_590con1844_2101")
(def dna-conversion1k :genome)
(def dna-conversion1 (mut/dna-conversion (coord/genomic-coordinate 333)
                                         (coord/genomic-coordinate 590)
                                         {:transcript nil
                                          :kind nil
                                          :coord-start (coord/genomic-coordinate 1844)
                                          :coord-end (coord/genomic-coordinate 2101)}))
(def dna-conversion1m {:mutation "dna-conversion"
                       :coord-start (coord/plain (coord/genomic-coordinate 333))
                       :coord-end (coord/plain (coord/genomic-coordinate 590))
                       :alt {:transcript nil
                             :kind nil
                             :coord-start (coord/plain (coord/genomic-coordinate 1844))
                             :coord-end (coord/plain (coord/genomic-coordinate 2101))}})

(def dna-conversion2s "415_1655conAC096506.5:g.409_1683")
(def dna-conversion2k :genome)
(def dna-conversion2 (mut/dna-conversion (coord/genomic-coordinate 415)
                                         (coord/genomic-coordinate 1655)
                                         {:transcript "AC096506.5"
                                          :kind :genome
                                          :coord-start (coord/genomic-coordinate 409)
                                          :coord-end (coord/genomic-coordinate 1683)}))

(def dna-conversion3s "15_355conNM_004006.1:20_360")
(def dna-conversion3k :coding-dna)
(def dna-conversion3 (mut/dna-conversion (coord/coding-dna-coordinate 15)
                                         (coord/coding-dna-coordinate 355)
                                         {:transcript "NM_004006.1"
                                          :kind nil
                                          :coord-start (coord/coding-dna-coordinate 20)
                                          :coord-end (coord/coding-dna-coordinate 360)}))

(deftest format-dna-conversion-test
  (testing "returns a string expression of a DNA conversion"
    (are [m s] (= (mut/format m nil) s)
      dna-conversion1 dna-conversion1s
      dna-conversion2 dna-conversion2s
      dna-conversion3 dna-conversion3s)))

(deftest parse-dna-conversion-test
  (testing "returns a correct DNAConversion"
    (are [s k m] (= (mut/parse-dna-conversion s k) m)
      dna-conversion1s dna-conversion1k dna-conversion1
      dna-conversion2s dna-conversion2k dna-conversion2
      dna-conversion3s dna-conversion3k dna-conversion3))
  (testing "invalid DNA conversion"
    (are [s k] (thrown? #?(:clj Throwable, :cljs js/Error)
                        (mut/parse-dna-conversion s k))
      "333con1844_2101" :genome
      "590_333con1844_2101" :genome)))

(deftest plain-dna-conversion-test
  (testing "returns a plain map representing DNAConversion"
    (is (= (mut/plain dna-conversion1) dna-conversion1m))))

(deftest restore-dna-conversion-test
  (testing "restores a plain map to DNAConversion"
    (is (= (mut/restore dna-conversion1m) dna-conversion1))))

;;; DNA - indel

(def dna-indel1s "6775delinsGA")
(def dna-indel1k :genome)
(def dna-indel1 (mut/dna-indel (coord/genomic-coordinate 6775) nil nil "GA"))
(def dna-indel1m {:mutation "dna-indel"
                  :coord-start (coord/plain (coord/genomic-coordinate 6775))
                  :coord-end nil
                  :ref nil
                  :alt "GA"})

(def dna-indel2s "6775delTinsGA")
(def dna-indel2ss "6775delinsGA")
(def dna-indel2k :genome)
(def dna-indel2 (mut/dna-indel (coord/genomic-coordinate 6775) nil "T" "GA"))

(def dna-indel3s "145_147delinsTGG")
(def dna-indel3k :coding-dna)
(def dna-indel3 (mut/dna-indel (coord/coding-dna-coordinate 145)
                               (coord/coding-dna-coordinate 147)
                               nil
                               "TGG"))

(def dna-indel4sb "1134_1138delinsNNNNNNNNNN")
(def dna-indel4sc "1134_1138delinsN[10]")
(def dna-indel4k :genome)
(def dna-indel4 (mut/dna-indel (coord/genomic-coordinate 1134)
                               (coord/genomic-coordinate 1138)
                               nil
                               "NNNNNNNNNN"))

(deftest format-dna-indel-test
  (testing "returns a string expression of a DNA indel"
    (are [m o s] (= (mut/format m o) s)
      dna-indel1 nil dna-indel1s
      dna-indel2 nil dna-indel2ss
      dna-indel2 {:show-bases? true} dna-indel2s
      dna-indel3 nil dna-indel3s
      dna-indel4 {:ins-format :auto} dna-indel4sc
      dna-indel4 {:ins-format :bases} dna-indel4sb
      dna-indel4 {:ins-format :count} dna-indel4sc)))

(deftest parse-dna-indel-test
  (testing "returns a correct DNAIndel"
    (are [s k m] (= (mut/parse-dna-indel s k) m)
      dna-indel1s dna-indel1k dna-indel1
      dna-indel2s dna-indel2k dna-indel2
      dna-indel3s dna-indel3k dna-indel3
      dna-indel4sb dna-indel4k dna-indel4
      dna-indel4sc dna-indel4k dna-indel4))
  (testing "invalid DNA indel"
    (are [s k] (thrown? #?(:clj Throwable, :cljs js/Error)
                        (mut/parse-dna-indel s k))
      "147_145delinsTGG" :coding-dna)))

(deftest plain-dna-indel-test
  (testing "returns a plain map representing DNAIndel"
    (is (= (mut/plain dna-indel1) dna-indel1m))))

(deftest restore-dna-indel-test
  (testing "restores a plain map to DNAIndel"
    (is (= (mut/restore dna-indel1m) dna-indel1))))

;;; DNA - alleles

(def dna-alleles1s "[123G>A;345del]")
(def dna-alleles1k :genome)
(def dna-alleles1 (mut/dna-alleles [(mut/dna-substitution (coord/genomic-coordinate 123) "G" ">" "A")
                                    (mut/dna-deletion (coord/genomic-coordinate 345) nil)]
                                   nil))

(def dna-alleles2s "[123G>A];[345del]")
(def dna-alleles2k :genome)
(def dna-alleles2 (mut/dna-alleles [(mut/dna-substitution (coord/genomic-coordinate 123) "G" ">" "A")]
                                   [(mut/dna-deletion (coord/genomic-coordinate 345) nil)]))

(def dna-alleles3s "2376[G>C];[G>C]")
(def dna-alleles3k :coding-dna)
(def dna-alleles3 (mut/dna-alleles [(mut/dna-substitution (coord/coding-dna-coordinate 2376)
                                                          "G" ">" "C")]
                                   [(mut/dna-substitution (coord/coding-dna-coordinate 2376)
                                                          "G" ">" "C")]))

(def dna-alleles4s "123_124[14];[18]")
(def dna-alleles4k :genome)
(def dna-alleles4 (mut/dna-alleles [(mut/dna-repeated-seqs (coord/genomic-coordinate 123)
                                                           (coord/genomic-coordinate 124)
                                                           nil 14)]
                                   [(mut/dna-repeated-seqs (coord/genomic-coordinate 123)
                                                           (coord/genomic-coordinate 124)
                                                           nil 18)]))

(deftest format-dna-alleles-test
  (testing "returns a string expression of a DNA alleles"
    (are [m s] (= (mut/format m) s)
      dna-alleles1 dna-alleles1s
      dna-alleles2 dna-alleles2s
      dna-alleles3 dna-alleles3s
      dna-alleles4 dna-alleles4s)))

(deftest parse-dna-alleles-test
  (testing "returns a correct DNAAlleles"
    (are [s k m] (= (mut/parse-dna-alleles s k) m)
      dna-alleles1s dna-alleles1k dna-alleles1
      dna-alleles2s dna-alleles2k dna-alleles2
      dna-alleles3s dna-alleles3k dna-alleles3
      dna-alleles4s dna-alleles4k dna-alleles4)))

;;; DNA - repeated sequences

(def dna-repeated-seqs1s-c "123_124[14]")
(def dna-repeated-seqs1s-b "123TG[14]")
(def dna-repeated-seqs1k :genome)
(def dna-repeated-seqs1-c (mut/dna-repeated-seqs (coord/genomic-coordinate 123)
                                                 (coord/genomic-coordinate 124)
                                                 nil
                                                 14))
(def dna-repeated-seqs1-b (mut/dna-repeated-seqs (coord/genomic-coordinate 123)
                                                 nil
                                                 "TG"
                                                 14))
(def dna-repeated-seqs1-a (mut/dna-repeated-seqs (coord/genomic-coordinate 123)
                                                 (coord/genomic-coordinate 124)
                                                 "TG"
                                                 14))
(def dna-repeated-seqs1m-c {:mutation "dna-repeated-seqs"
                            :coord-start (coord/plain (coord/genomic-coordinate 123))
                            :coord-end (coord/plain (coord/genomic-coordinate 124))
                            :ref nil
                            :ncopy 14})

(def dna-repeated-seqs2s "-128_-126[(600_800)]")
(def dna-repeated-seqs2k :coding-dna)
(def dna-repeated-seqs2 (mut/dna-repeated-seqs (coord/coding-dna-coordinate 128 0 :upstream)
                                               (coord/coding-dna-coordinate 126 0 :upstream)
                                               nil
                                               [600 800]))

(deftest format-dna-repeated-seqs-test
  (testing "returns a string expression of a DNA repeated sequences"
    (are [m o s] (= (mut/format m o) s)
      dna-repeated-seqs1-c nil dna-repeated-seqs1s-c
      dna-repeated-seqs1-b nil dna-repeated-seqs1s-b
      dna-repeated-seqs1-a nil dna-repeated-seqs1s-b
      dna-repeated-seqs1-a {:range-format :bases} dna-repeated-seqs1s-b
      dna-repeated-seqs1-a {:range-format :coord} dna-repeated-seqs1s-c
      dna-repeated-seqs2 nil dna-repeated-seqs2s))
  (testing "throws exception"
    (are [m o] (thrown? #?(:clj Exception, :cljs js/Error) (mut/format m o))
      dna-repeated-seqs1-c {:range-format :bases}
      dna-repeated-seqs1-b {:range-format :coord})))

(deftest parse-dna-repeated-seqs-test
  (testing "returns a correct DNARepeatedSeqs"
    (are [s k m] (= (mut/parse-dna-repeated-seqs s k) m)
      dna-repeated-seqs1s-c dna-repeated-seqs1k dna-repeated-seqs1-c
      dna-repeated-seqs1s-b dna-repeated-seqs1k dna-repeated-seqs1-b
      dna-repeated-seqs2s dna-repeated-seqs2k dna-repeated-seqs2))
  (testing "invalid DNA repeated sequences"
    (are [s k] (thrown? #?(:clj Throwable, :cljs js/Error)
                        (mut/parse-dna-repeated-seqs s k))
      "124_123[14]" :genome)))

(deftest plain-dna-repeated-seqs-test
  (testing "returns a plain map representing DNARepeatedSeqs"
    (is (= (mut/plain dna-repeated-seqs1-c) dna-repeated-seqs1m-c))))

(deftest restore-dna-repeated-seqs-test
  (testing "restores a plain map to DNARepeatedSeqs"
    (is (= (mut/restore dna-repeated-seqs1m-c) dna-repeated-seqs1-c))))

;;; DNA - general

(deftest parse-dna-test
  (are [s c] (instance? c (mut/parse-dna s :genome))
    "45576A>C"            clj_hgvs.mutation.DNASubstitution
    "6_8delTGC"           clj_hgvs.mutation.DNADeletion
    "6_8dupTGC"           clj_hgvs.mutation.DNADuplication
    "5756_5757insAGG"     clj_hgvs.mutation.DNAInsertion
    "1077_1080inv"        clj_hgvs.mutation.DNAInversion
    "333_590con1844_2101" clj_hgvs.mutation.DNAConversion
    "6775delinsGA"        clj_hgvs.mutation.DNAIndel
    "[123G>A;345del]"     clj_hgvs.mutation.DNAAlleles
    "123_124[14]"         clj_hgvs.mutation.DNARepeatedSeqs))

;;; RNA mutations

;;; RNA - substitution

(def rna-substitution1s "76a>c")
(def rna-substitution1 (mut/rna-substitution (coord/rna-coordinate 76 nil nil) "a" "c"))
(def rna-substitution1m {:mutation "rna-substitution"
                         :coord (coord/plain (coord/rna-coordinate 76 nil nil))
                         :ref "a"
                         :alt "c"})

(def rna-substitution2s "-14g>c")
(def rna-substitution2 (mut/rna-substitution (coord/rna-coordinate 14 0 :upstream) "g" "c"))

(def rna-substitution3s "*46u>a")
(def rna-substitution3 (mut/rna-substitution (coord/rna-coordinate 46 0 :downstream) "u" "a"))

(deftest format-rna-substitution-test
  (testing "returns a string expression of a RNA substitution"
    (are [m s] (= (mut/format m nil) s)
      rna-substitution1 rna-substitution1s
      rna-substitution2 rna-substitution2s
      rna-substitution3 rna-substitution3s)))

(deftest parse-rna-substitution-test
  (testing "returns a correct RNASubstitution"
    (are [s m] (= (mut/parse-rna-substitution s) m)
      rna-substitution1s rna-substitution1
      rna-substitution2s rna-substitution2
      rna-substitution3s rna-substitution3)))

(deftest plain-rna-substitution-test
  (testing "returns a plain map representing RNASubstitution"
    (is (= (mut/plain rna-substitution1) rna-substitution1m))))

(deftest restore-rna-substitution-test
  (testing "restores a plain map to RNASubstitution"
    (is (= (mut/restore rna-substitution1m) rna-substitution1))))

;;; RNA - deletion

(def rna-deletion1s "7del")
(def rna-deletion1 (mut/rna-deletion (coord/rna-coordinate 7 nil nil) nil))
(def rna-deletion1m {:mutation "rna-deletion"
                     :coord-start (coord/plain (coord/rna-coordinate 7 nil nil))
                     :coord-end nil
                     :ref nil})

(def rna-deletion2s "7delu")
(def rna-deletion2ss "7del")
(def rna-deletion2 (mut/rna-deletion (coord/rna-coordinate 7 nil nil) nil "u"))

(def rna-deletion3s "6_8del")
(def rna-deletion3 (mut/rna-deletion (coord/rna-coordinate 6 nil nil)
                                     (coord/rna-coordinate 8 nil nil)))

(def rna-deletion4s "(4072_5145)del")
(def rna-deletion4 "TODO")

(deftest format-rna-deletion-test
  (testing "returns a string expression of a RNA deletion"
    (are [m o s] (= (mut/format m o) s)
      rna-deletion1 nil rna-deletion1s
      rna-deletion2 nil rna-deletion2ss
      rna-deletion2 {:show-bases? true} rna-deletion2s
      rna-deletion3 nil rna-deletion3s
      ;; rna-deletion4 nil rna-deletion4s ; TODO
      )))

(deftest parse-rna-deletion-test
  (testing "returns a correct RNADeletion"
    (are [s m] (= (mut/parse-rna-deletion s) m)
      rna-deletion1s rna-deletion1
      rna-deletion2s rna-deletion2
      rna-deletion3s rna-deletion3
      ;; rna-deletion4s rna-deletion4 ; TODO
      ))
  (testing "invalid RNA deletion"
    (are [s] (thrown? #?(:clj Throwable, :cljs js/Error)
                      (mut/parse-rna-deletion s))
      "8_6del")))

(deftest plain-rna-deletion-test
  (testing "returns a plain map representing RNADeletion"
    (is (= (mut/plain rna-deletion1) rna-deletion1m))))

(deftest restore-rna-deletion-test
  (testing "restores a plain map to RNADeletion"
    (is (= (mut/restore rna-deletion1m) rna-deletion1))))

;;; RNA - duplication

(def rna-duplication1s "7dup")
(def rna-duplication1 (mut/rna-duplication (coord/rna-coordinate 7 nil nil) nil))
(def rna-duplication1m {:mutation "rna-duplication"
                        :coord-start (coord/plain (coord/rna-coordinate 7 nil nil))
                        :coord-end nil
                        :ref nil})

(def rna-duplication2s "7dupu")
(def rna-duplication2ss "7dup")
(def rna-duplication2 (mut/rna-duplication (coord/rna-coordinate 7 nil nil) nil "u"))

(def rna-duplication3s "6_8dup")
(def rna-duplication3 (mut/rna-duplication (coord/rna-coordinate 6 nil nil)
                                           (coord/rna-coordinate 8 nil nil)))

(deftest format-rna-duplication-test
  (testing "returns a string expression of a RNA duplication"
    (are [m o s] (= (mut/format m o) s)
      rna-duplication1 nil rna-duplication1s
      rna-duplication2 nil rna-duplication2ss
      rna-duplication2 {:show-bases? true} rna-duplication2s
      rna-duplication3 nil rna-duplication3s)))

(deftest parse-rna-duplication-test
  (testing "returns a correct RNADuplication"
    (are [s m] (= (mut/parse-rna-duplication s) m)
      rna-duplication1s rna-duplication1
      rna-duplication2s rna-duplication2
      rna-duplication3s rna-duplication3))
  (testing "invalid RNA duplication"
    (are [s] (thrown? #?(:clj Throwable, :cljs js/Error)
                      (mut/parse-rna-duplication s))
      "8_6dup")))

(deftest plain-rna-duplication-test
  (testing "returns a plain map representing RNADuplication"
    (is (= (mut/plain rna-duplication1) rna-duplication1m))))

(deftest restore-rna-duplication-test
  (testing "restores a plain map to RNADuplication"
    (is (= (mut/restore rna-duplication1m) rna-duplication1))))

;;; RNA - insertion

(def rna-insertion1s "756_757insacu")
(def rna-insertion1 (mut/rna-insertion (coord/rna-coordinate 756 nil nil)
                                       (coord/rna-coordinate 757 nil nil)
                                       "acu"))
(def rna-insertion1m {:mutation "rna-insertion"
                      :coord-start (coord/plain (coord/rna-coordinate 756 nil nil))
                      :coord-end (coord/plain (coord/rna-coordinate 757 nil nil))
                      :alt "acu"})

(def rna-insertion2s "431_432insn[5]")
(def rna-insertion2 (mut/rna-insertion (coord/rna-coordinate 431 nil nil)
                                       (coord/rna-coordinate 432 nil nil)
                                       "nnnnn"))

(def rna-insertion3s "123_124insL37425.1:23_361")
(def rna-insertion3 (mut/rna-insertion (coord/rna-coordinate 123 nil nil)
                                       (coord/rna-coordinate 124 nil nil)
                                       {:genbank "L37425.1"
                                        :coord-start 23
                                        :coord-end 361}))

(deftest format-rna-insertion-test
  (testing "returns a string expression of a RNA insertion"
    (are [m s] (= (mut/format m nil) s)
      rna-insertion1 rna-insertion1s
      rna-insertion2 rna-insertion2s
      rna-insertion3 rna-insertion3s)))

(deftest parse-rna-insertion-test
  (testing "returns a correct RNAInsertion"
    (are [s m] (= (mut/parse-rna-insertion s) m)
      rna-insertion1s rna-insertion1
      rna-insertion2s rna-insertion2
      rna-insertion3s rna-insertion3))
  (testing "invalid RNA insertion"
    (are [s] (thrown? #?(:clj Throwable, :cljs js/Error)
                      (mut/parse-rna-insertion s))
      "756insacu"
      "757_756insacu")))

(deftest plain-rna-insertion-test
  (testing "returns a plain map representing RNAInsertion"
    (is (= (mut/plain rna-insertion1) rna-insertion1m))))

(deftest restore-rna-insertion-test
  (testing "restores a plain map to RNAInsertion"
    (is (= (mut/restore rna-insertion1m) rna-insertion1))))

;;; RNA - inversion

(def rna-inversion1s "177_180inv")
(def rna-inversion1 (mut/rna-inversion (coord/rna-coordinate 177 nil nil)
                                       (coord/rna-coordinate 180 nil nil)))
(def rna-inversion1m {:mutation "rna-inversion"
                      :coord-start (coord/plain (coord/rna-coordinate 177 nil nil))
                      :coord-end (coord/plain (coord/rna-coordinate 180 nil nil))})

(deftest format-rna-inversion-test
  (testing "returns a string expression of a RNA inversion"
    (are [m s] (= (mut/format m nil) s)
      rna-inversion1 rna-inversion1s)))

(deftest parse-rna-inversion-test
  (testing "returns a correct RNAInversion"
    (are [s m] (= (mut/parse-rna-inversion s) m)
      rna-inversion1s rna-inversion1))
  (testing "invalid RNA inversion"
    (are [s] (thrown? #?(:clj Throwable, :cljs js/Error)
                      (mut/parse-rna-inversion s))
      "177inv"
      "180_177inv")))

(deftest plain-rna-inversion-test
  (testing "returns a plain map representing RNAInversion"
    (is (= (mut/plain rna-inversion1) rna-inversion1m))))

(deftest restore-rna-inversion-test
  (testing "restores a plain map to RNAInversion"
    (is (= (mut/restore rna-inversion1m) rna-inversion1))))

;;; RNA - conversion

(def rna-conversion1s "123_345con888_1110")
(def rna-conversion1 (mut/rna-conversion (coord/rna-coordinate 123 nil nil)
                                         (coord/rna-coordinate 345 nil nil)
                                         {:transcript nil
                                          :coord-start (coord/rna-coordinate 888 nil nil)
                                          :coord-end (coord/rna-coordinate 1110 nil nil)}))
(def rna-conversion1m {:mutation "rna-conversion"
                       :coord-start (coord/plain (coord/rna-coordinate 123 nil nil))
                       :coord-end (coord/plain (coord/rna-coordinate 345 nil nil))
                       :alt {:transcript nil
                             :coord-start (coord/plain (coord/rna-coordinate 888 nil nil))
                             :coord-end (coord/plain (coord/rna-coordinate 1110 nil nil))}})

(def rna-conversion2s "415_1655conAC096506.5:409_1649")
(def rna-conversion2 (mut/rna-conversion (coord/rna-coordinate 415 nil nil)
                                         (coord/rna-coordinate 1655 nil nil)
                                         {:transcript "AC096506.5"
                                          :coord-start (coord/rna-coordinate 409 nil nil)
                                          :coord-end (coord/rna-coordinate 1649 nil nil)}))

(deftest format-rna-conversion-test
  (testing "returns a string expression of a RNA conversion"
    (are [m s] (= (mut/format m nil) s)
      rna-conversion1 rna-conversion1s
      rna-conversion2 rna-conversion2s)))

(deftest parse-rna-conversion-test
  (testing "returns a correct RNAConversion"
    (are [s m] (= (mut/parse-rna-conversion s) m)
      rna-conversion1s rna-conversion1
      rna-conversion2s rna-conversion2))
  (testing "invalid RNA conversion"
    (are [s] (thrown? #?(:clj Throwable, :cljs js/Error)
                      (mut/parse-rna-conversion s))
      "123con888_1110"
      "345_123con888_1110")))

(deftest plain-rna-conversion-test
  (testing "returns a plain map representing RNAConversion"
    (is (= (mut/plain rna-conversion1) rna-conversion1m))))

(deftest restore-rna-conversion-test
  (testing "restores a plain map to RNAConversion"
    (is (= (mut/restore rna-conversion1m) rna-conversion1))))

;;; RNA - indel

(def rna-indel1s "775delinsga")
(def rna-indel1 (mut/rna-indel (coord/rna-coordinate 775 nil nil) nil nil "ga"))
(def rna-indel1m {:mutation "rna-indel"
                  :coord-start (coord/plain (coord/rna-coordinate 775 nil nil))
                  :coord-end nil
                  :ref nil
                  :alt "ga"})

(def rna-indel2s "775deluinsga")
(def rna-indel2ss "775delinsga")
(def rna-indel2 (mut/rna-indel (coord/rna-coordinate 775 nil nil) nil "u" "ga"))

(def rna-indel3s "775_777delinsc")
(def rna-indel3 (mut/rna-indel (coord/rna-coordinate 775 nil nil)
                               (coord/rna-coordinate 777 nil nil)
                               nil
                               "c"))

(deftest format-rna-indel-test
  (testing "returns a string expression of a RNA indel"
    (are [m o s] (= (mut/format m o) s)
      rna-indel1 nil rna-indel1s
      rna-indel2 nil rna-indel2ss
      rna-indel2 {:show-bases? true} rna-indel2s
      rna-indel3 nil rna-indel3s)))

(deftest parse-rna-indel-test
  (testing "returns a correct RNAIndel"
    (are [s m] (= (mut/parse-rna-indel s) m)
      rna-indel1s rna-indel1
      rna-indel2s rna-indel2
      rna-indel3s rna-indel3))
  (testing "invalid RNA indel"
    (are [s] (thrown? #?(:clj Throwable, :cljs js/Error)
                      (mut/parse-rna-indel s))
      "777_775delinsc")))

(deftest plain-rna-indel-test
  (testing "returns a plain map representing RNAIndel"
    (is (= (mut/plain rna-indel1) rna-indel1m))))

(deftest restore-rna-indel-test
  (testing "restores a plain map to RNAIndel"
    (is (= (mut/restore rna-indel1m) rna-indel1))))

;;; RNA - alleles

(def rna-alleles1s "[76a>u;103del]")
(def rna-alleles1 (mut/rna-alleles [(mut/rna-substitution (coord/rna-coordinate 76 nil nil) "a" "u")
                                    (mut/rna-deletion (coord/rna-coordinate 103 nil nil) nil)]
                                   nil))

(def rna-alleles2s "[76a>u];[103del]")
(def rna-alleles2 (mut/rna-alleles [(mut/rna-substitution (coord/rna-coordinate 76 nil nil) "a" "u")]
                                   [(mut/rna-deletion (coord/rna-coordinate 103 nil nil) nil)]))

(def rna-alleles3s "76[a>u];[a>u]")
(def rna-alleles3 (mut/rna-alleles
                   [(mut/rna-substitution (coord/rna-coordinate 76 nil nil)
                                          "a" "u")]
                   [(mut/rna-substitution (coord/rna-coordinate 76 nil nil)
                                          "a" "u")]))

(def rna-alleles4s "-124_-123[14];[18]")
(def rna-alleles4 (mut/rna-alleles
                   [(mut/rna-repeated-seqs (coord/rna-coordinate 124 0 :upstream)
                                           (coord/rna-coordinate 123 0 :upstream)
                                           nil
                                           14)]
                   [(mut/rna-repeated-seqs (coord/rna-coordinate 124 0 :upstream)
                                           (coord/rna-coordinate 123 0 :upstream)
                                           nil
                                           18)]))

(deftest format-rna-alleles-test
  (testing "returns a string expression of a RNA alleles"
    (are [m s] (= (mut/format m) s)
      rna-alleles1 rna-alleles1s
      rna-alleles2 rna-alleles2s
      rna-alleles3 rna-alleles3s
      rna-alleles4 rna-alleles4s)))

(deftest parse-rna-alleles-test
  (testing "returns a correct RNAAlleles"
    (are [s m] (= (mut/parse-rna-alleles s) m)
      rna-alleles1s rna-alleles1
      rna-alleles2s rna-alleles2
      rna-alleles3s rna-alleles3
      rna-alleles4s rna-alleles4)))

;;; RNA - repeated sequences

(def rna-repeated-seqs1s-c "-124_-123[14]")
(def rna-repeated-seqs1s-b "-124ug[14]")
(def rna-repeated-seqs1-c (mut/rna-repeated-seqs (coord/rna-coordinate 124 0 :upstream)
                                                 (coord/rna-coordinate 123 0 :upstream)
                                                 nil
                                                 14))
(def rna-repeated-seqs1-b (mut/rna-repeated-seqs (coord/rna-coordinate 124 0 :upstream)
                                                 nil
                                                 "ug"
                                                 14))
(def rna-repeated-seqs1-a (mut/rna-repeated-seqs (coord/rna-coordinate 124 0 :upstream)
                                                 (coord/rna-coordinate 123 0 :upstream)
                                                 "ug"
                                                 14))
(def rna-repeated-seqs1m-c {:mutation "rna-repeated-seqs"
                            :coord-start (coord/plain (coord/rna-coordinate 124 0 :upstream))
                            :coord-end (coord/plain (coord/rna-coordinate 123 0 :upstream))
                            :ref nil
                            :ncopy 14})

(def rna-repeated-seqs2s "-128_-126[(600_800)]")
(def rna-repeated-seqs2 (mut/rna-repeated-seqs (coord/rna-coordinate 128 0 :upstream)
                                               (coord/rna-coordinate 126 0 :upstream)
                                               nil
                                               [600 800]))

(deftest format-rna-repeated-seqs-test
  (testing "returns a string expression of a RNA repeated-seqs"
    (are [m o s] (= (mut/format m o) s)
      rna-repeated-seqs1-c nil rna-repeated-seqs1s-c
      rna-repeated-seqs1-b nil rna-repeated-seqs1s-b
      rna-repeated-seqs1-a nil rna-repeated-seqs1s-b
      rna-repeated-seqs1-a {:range-format :bases} rna-repeated-seqs1s-b
      rna-repeated-seqs1-a {:range-format :coord} rna-repeated-seqs1s-c
      rna-repeated-seqs2 nil rna-repeated-seqs2s))
  (testing "throws exception"
    (are [m o] (thrown? #?(:clj Exception, :cljs js/Error) (mut/format m o))
      dna-repeated-seqs1-c {:range-format :bases}
      dna-repeated-seqs1-b {:range-format :coord})))

(deftest parse-rna-repeated-seqs-test
  (testing "returns a correct RNARepeatedSeqs"
    (are [s m] (= (mut/parse-rna-repeated-seqs s) m)
      rna-repeated-seqs1s-c rna-repeated-seqs1-c
      rna-repeated-seqs1s-b rna-repeated-seqs1-b
      rna-repeated-seqs2s rna-repeated-seqs2))
  (testing "invalid RNA repeated-seqs"
    (are [s] (thrown? #?(:clj Throwable, :cljs js/Error)
                      (mut/parse-rna-repeated-seqs s))
      "-123_-124[14]")))

(deftest plain-rna-repeated-seqs-test
  (testing "returns a plain map representing RNARepeatedSeqs"
    (is (= (mut/plain rna-repeated-seqs1-c) rna-repeated-seqs1m-c))))

(deftest restore-rna-repeated-seqs-test
  (testing "restores a plain map to RNARepeatedSeqs"
    (is (= (mut/restore rna-repeated-seqs1m-c) rna-repeated-seqs1-c))))

;;; RNA - no RNA detected

(def no-rna-s "0")
(def no-rna (mut/no-rna))
(def no-rna-m {:mutation "no-rna"})

(deftest format-no-rna-test
  (is (= (mut/format no-rna) no-rna-s)))

(deftest plain-no-rna-test
  (is (= (mut/plain no-rna) no-rna-m)))

(deftest restore-no-rna-test
  (is (= (mut/restore no-rna-m) no-rna)))

;;; RNA - unknown mutation

(def rna-unknown-s "?")
(def rna-unknown (mut/rna-unknown-mutation))
(def rna-unknown-m {:mutation "rna-unknown"})

(deftest format-rna-unknown-mutation-test
  (is (= (mut/format rna-unknown) rna-unknown-s)))

(deftest plain-rna-unknown-mutation-test
  (is (= (mut/plain rna-unknown) rna-unknown-m)))

(deftest restore-rna-unknown-mutation-test
  (is (= (mut/restore rna-unknown-m) rna-unknown)))

;;; RNA - splice affected

(def rna-splice-affected-s "spl")
(def rna-splice-affected (mut/rna-splice-affected))
(def rna-splice-affected-m {:mutation "rna-splice-affected"})

(deftest format-rna-splice-affected-test
  (is (= (mut/format rna-splice-affected) rna-splice-affected-s)))

(deftest plain-rna-splice-affected-test
  (is (= (mut/plain rna-splice-affected) rna-splice-affected-m)))

(deftest restore-rna-splice-affected-test
  (is (= (mut/restore rna-splice-affected-m) rna-splice-affected)))

;;; RNA - general

(deftest parse-rna-test
  (are [s c] (instance? c (mut/parse-rna s))
    "76a>c"              clj_hgvs.mutation.RNASubstitution
    "6_8del"             clj_hgvs.mutation.RNADeletion
    "6_8dup"             clj_hgvs.mutation.RNADuplication
    "756_757insacu"      clj_hgvs.mutation.RNAInsertion
    "177_180inv"         clj_hgvs.mutation.RNAInversion
    "123_345con888_1110" clj_hgvs.mutation.RNAConversion
    "775_777delinsc"     clj_hgvs.mutation.RNAIndel
    "[76a>u;103del]"     clj_hgvs.mutation.RNAAlleles
    "-124_-123[14]"      clj_hgvs.mutation.RNARepeatedSeqs
    "(?)"                clj_hgvs.mutation.UncertainMutation
    "0"                  clj_hgvs.mutation.NoRNA
    "?"                  clj_hgvs.mutation.RNAUnknownMutation
    "="                  clj_hgvs.mutation.RNANoEffect
    "spl"                clj_hgvs.mutation.RNASpliceAffected))

;;; Protein mutations

;;; Protein - substitution

(def protein-substitution1s "Arg54Ser")
(def protein-substitution1ss "R54S")
(def protein-substitution1 (mut/protein-substitution "Arg"
                                                     (coord/protein-coordinate 54)
                                                     "Ser"))
(def protein-substitution1m {:mutation "protein-substitution"
                             :ref "Arg"
                             :coord (coord/plain (coord/protein-coordinate 54))
                             :alt "Ser"})

(def protein-substitution2s "Cys123=")
(def protein-substitution2ss "C123=")
(def protein-substitution2 (mut/protein-substitution "Cys"
                                                     (coord/protein-coordinate 123)
                                                     "Cys"))

(def protein-substitution3s "Ter123=")
(def protein-substitution3ss "*123=")
(def protein-substitution3 (mut/protein-substitution "Ter"
                                                     (coord/protein-coordinate 123)
                                                     "Ter"))

(deftest format-protein-substitution-test
  (testing "returns a string expression of a protein substitution"
    (are [m o s] (= (mut/format m o) s)
      protein-substitution1 nil protein-substitution1s
      protein-substitution1 {:amino-acid-format :short} protein-substitution1ss
      protein-substitution2 nil protein-substitution2s
      protein-substitution2 {:amino-acid-format :short} protein-substitution2ss
      protein-substitution3 nil protein-substitution3s
      protein-substitution3 {:amino-acid-format :short} protein-substitution3ss)))

(deftest parse-protein-substitution-test
  (testing "returns a correct ProteinSubstitution"
    (are [s m] (= (mut/parse-protein-substitution s) m)
      protein-substitution1s protein-substitution1
      protein-substitution1ss protein-substitution1
      protein-substitution2s protein-substitution2
      protein-substitution2ss protein-substitution2
      protein-substitution3s protein-substitution3
      protein-substitution3ss protein-substitution3)))

(deftest plain-protein-substitution-test
  (testing "returns a plain map representing ProteinSubstitution"
    (is (= (mut/plain protein-substitution1) protein-substitution1m))))

(deftest restore-protein-substitution-test
  (testing "restores a plain map to ProteinSubstitution"
    (is (= (mut/restore protein-substitution1m) protein-substitution1))))

;;; Protein - deletion

(def protein-deletion1s "Ala3del")
(def protein-deletion1ss "A3del")
(def protein-deletion1 (mut/protein-deletion "Ala" (coord/protein-coordinate 3)))
(def protein-deletion1m {:mutation "protein-deletion"
                         :ref-start "Ala"
                         :coord-start (coord/plain (coord/protein-coordinate 3))
                         :ref-end nil
                         :coord-end nil})

(def protein-deletion2s "Cys76_Glu79del")
(def protein-deletion2ss "C76_E79del")
(def protein-deletion2 (mut/protein-deletion "Cys" (coord/protein-coordinate 76)
                                             "Glu" (coord/protein-coordinate 79)))
(deftest format-protein-deletion-test
  (testing "returns a string expression of a protein deletion"
    (are [m o s] (= (mut/format m o) s)
      protein-deletion1 nil protein-deletion1s
      protein-deletion1 {:amino-acid-format :short} protein-deletion1ss
      protein-deletion2 nil protein-deletion2s
      protein-deletion2 {:amino-acid-format :short} protein-deletion2ss)))

(deftest parse-protein-deletion-test
  (testing "returns a correct ProteinDeletion"
    (are [s m] (= (mut/parse-protein-deletion s) m)
      protein-deletion1s protein-deletion1
      protein-deletion1ss protein-deletion1
      protein-deletion2s protein-deletion2
      protein-deletion2ss protein-deletion2))
  (testing "invalid protein deletion"
    (are [s] (thrown? #?(:clj Throwable, :cljs js/Error)
                      (mut/parse-protein-deletion s))
      "Glu79_Cys76del")))

(deftest plain-protein-deletion-test
  (testing "returns a plain map representing ProteinDeletion"
    (is (= (mut/plain protein-deletion1) protein-deletion1m))))

(deftest restore-protein-deletion-test
  (testing "restores a plain map to ProteinDeletion"
    (is (= (mut/restore protein-deletion1m) protein-deletion1))))

;;; Protein - duplication

(def protein-duplication1s "Ala3dup")
(def protein-duplication1ss "A3dup")
(def protein-duplication1 (mut/protein-duplication "Ala" (coord/protein-coordinate 3)))
(def protein-duplication1m {:mutation "protein-duplication"
                            :ref-start "Ala"
                            :coord-start (coord/plain (coord/protein-coordinate 3))
                            :ref-end nil
                            :coord-end nil})

(def protein-duplication2s "Ala3_Ser5dup")
(def protein-duplication2ss "A3_S5dup")
(def protein-duplication2 (mut/protein-duplication "Ala" (coord/protein-coordinate 3)
                                                   "Ser" (coord/protein-coordinate 5)))

(deftest format-protein-duplication-test
  (testing "returns a string expression of a protein duplication"
    (are [m o s] (= (mut/format m o) s)
      protein-duplication1 nil protein-duplication1s
      protein-duplication1 {:amino-acid-format :short} protein-duplication1ss
      protein-duplication2 nil protein-duplication2s
      protein-duplication2 {:amino-acid-format :short} protein-duplication2ss)))

(deftest parse-protein-duplication-test
  (testing "returns a correct ProteinDuplication"
    (are [s m] (= (mut/parse-protein-duplication s) m)
      protein-duplication1s protein-duplication1
      protein-duplication1ss protein-duplication1
      protein-duplication2s protein-duplication2
      protein-duplication2ss protein-duplication2))
  (testing "invalid protein duplication"
    (are [s] (thrown? #?(:clj Throwable, :cljs js/Error)
                      (mut/parse-protein-duplication s))
      "Ser5_Ala3dup")))

(deftest plain-protein-duplication-test
  (testing "returns a plain map representing ProteinDuplication"
    (is (= (mut/plain protein-duplication1) protein-duplication1m))))

(deftest restore-protein-duplication-test
  (testing "restores a plain map to ProteinDuplication"
    (is (= (mut/restore protein-duplication1m) protein-duplication1))))

;;; Protein - insertion

(def protein-insertion1s "Lys23_Leu24insArgSerGln")
(def protein-insertion1ss "K23_L24insRSQ")
(def protein-insertion1 (mut/protein-insertion "Lys" (coord/protein-coordinate 23)
                                               "Leu" (coord/protein-coordinate 24)
                                               ["Arg" "Ser" "Gln"]))
(def protein-insertion1m {:mutation "protein-insertion"
                          :ref-start "Lys"
                          :coord-start (coord/plain (coord/protein-coordinate 23))
                          :ref-end "Leu"
                          :coord-end (coord/plain (coord/protein-coordinate 24))
                          :alts ["Arg" "Ser" "Gln"]})

(def protein-insertion2s "Lys23_Leu24insArgSerTer")
(def protein-insertion2ss "K23_L24insRS*")
(def protein-insertion2 (mut/protein-insertion "Lys" (coord/protein-coordinate 23)
                                               "Leu" (coord/protein-coordinate 24)
                                               ["Arg" "Ser" "Ter"]))

(def protein-insertion3s "Arg78_Gly79ins5")
(def protein-insertion3 (mut/protein-insertion "Arg" (coord/protein-coordinate 78)
                                               "Gly" (coord/protein-coordinate 79)
                                               ["Xaa" "Xaa" "Xaa" "Xaa" "Xaa"]))

(deftest format-protein-insertion-test
  (testing "returns a string expression of a protein insertion"
    (are [m o s] (= (mut/format m o) s)
      protein-insertion1 nil protein-insertion1s
      protein-insertion1 {:amino-acid-format :short} protein-insertion1ss
      protein-insertion2 nil protein-insertion2s
      protein-insertion2 {:amino-acid-format :short} protein-insertion2ss
      protein-insertion3 nil protein-insertion3s)))

(deftest parse-protein-insertion-test
  (testing "returns a correct ProteinInsertion"
    (are [s m] (= (mut/parse-protein-insertion s) m)
      protein-insertion1s protein-insertion1
      protein-insertion1ss protein-insertion1
      protein-insertion2s protein-insertion2
      protein-insertion2ss protein-insertion2
      protein-insertion3s protein-insertion3))
  (testing "invalid protein insertion"
    (are [s] (thrown? #?(:clj Throwable, :cljs js/Error)
                      (mut/parse-protein-insertion s))
      "Lys23insArgSerGln"
      "Leu24_Lys23insArgSerGln")))

(deftest plain-protein-insertion-test
  (testing "returns a plain map representing ProteinInsertion"
    (is (= (mut/plain protein-insertion1) protein-insertion1m))))

(deftest restore-protein-insertion-test
  (testing "restores a plain map to ProteinInsertion"
    (is (= (mut/restore protein-insertion1m) protein-insertion1))))

;;; Protein - indel

(def protein-indel1s "Cys28delinsTrpVal")
(def protein-indel1ss "C28delinsWV")
(def protein-indel1 (mut/protein-indel "Cys" (coord/protein-coordinate 28)
                                       nil nil
                                       ["Trp" "Val"]))
(def protein-indel1m {:mutation "protein-indel"
                      :ref-start "Cys"
                      :coord-start (coord/plain (coord/protein-coordinate 28))
                      :ref-end nil
                      :coord-end nil
                      :alts ["Trp" "Val"]})

(def protein-indel2s "Cys28_Lys29delinsTrp")
(def protein-indel2ss "C28_K29delinsW")
(def protein-indel2 (mut/protein-indel "Cys" (coord/protein-coordinate 28)
                                       "Lys" (coord/protein-coordinate 29)
                                       ["Trp"]))

(def protein-indel3s "Cys28_Lys29delinsTer")
(def protein-indel3ss "C28_K29delins*")
(def protein-indel3 (mut/protein-indel "Cys" (coord/protein-coordinate 28)
                                       "Lys" (coord/protein-coordinate 29)
                                       ["Ter"]))

(deftest format-protein-indel-test
  (testing "returns a string expression of a protein indel"
    (are [m o s] (= (mut/format m o) s)
      protein-indel1 nil protein-indel1s
      protein-indel1 {:amino-acid-format :short} protein-indel1ss
      protein-indel2 nil protein-indel2s
      protein-indel2 {:amino-acid-format :short} protein-indel2ss
      protein-indel3 nil protein-indel3s
      protein-indel3 {:amino-acid-format :short} protein-indel3ss)))

(deftest parse-protein-indel-test
  (testing "returns a correct ProteinIndel"
    (are [s m] (= (mut/parse-protein-indel s) m)
      protein-indel1s protein-indel1
      protein-indel1ss protein-indel1
      protein-indel2s protein-indel2
      protein-indel2ss protein-indel2
      protein-indel3s protein-indel3
      protein-indel3ss protein-indel3))
  (testing "invalid protein indel"
    (are [s] (thrown? #?(:clj Throwable, :cljs js/Error)
                      (mut/parse-protein-indel s))
      "Lys29_Cys28delinsTrp")))

(deftest plain-protein-indel-test
  (testing "returns a plain map representing ProteinIndel"
    (is (= (mut/plain protein-indel1) protein-indel1m))))

(deftest restore-protein-indel-test
  (testing "restores a plain map to ProteinIndel"
    (is (= (mut/restore protein-indel1m) protein-indel1))))

;;; Protein - alleles

(def protein-alleles1s "[Ser73Arg;Asn603del]")
(def protein-alleles1 (mut/protein-alleles
                       [(mut/protein-substitution "Ser" (coord/protein-coordinate 73) "Arg")
                        (mut/protein-deletion "Asn" (coord/protein-coordinate 603))]
                       nil))

(def protein-alleles2s "[Ser73Arg];[Asn603del]")
(def protein-alleles2 (mut/protein-alleles
                       [(mut/protein-substitution "Ser" (coord/protein-coordinate 73) "Arg")]
                       [(mut/protein-deletion "Asn" (coord/protein-coordinate 603))]))

(def protein-alleles3s "[Ser73Arg];[Ser73=]")
(def protein-alleles3 (mut/protein-alleles
                       [(mut/protein-substitution "Ser" (coord/protein-coordinate 73) "Arg")]
                       [(mut/protein-substitution "Ser" (coord/protein-coordinate 73) "Ser")]))

(def protein-alleles4s "Ala2[10];[11]")
(def protein-alleles4 (mut/protein-alleles
                       [(mut/protein-repeated-seqs "Ala" (coord/protein-coordinate 2)
                                                   nil nil
                                                   10)]
                       [(mut/protein-repeated-seqs "Ala" (coord/protein-coordinate 2)
                                                   nil nil
                                                   11)]))

(deftest format-protein-alleles-test
  (testing "returns a string expression of a protein alleles"
    (are [m s] (= (mut/format m) s)
      protein-alleles1 protein-alleles1s
      protein-alleles2 protein-alleles2s
      protein-alleles3 protein-alleles3s
      protein-alleles4 protein-alleles4s)))

(deftest parse-protein-alleles-test
  (testing "returns a correct ProteinAlleles"
    (are [s m] (= (mut/parse-protein-alleles s) m)
      protein-alleles1s protein-alleles1
      protein-alleles2s protein-alleles2
      protein-alleles3s protein-alleles3
      protein-alleles4s protein-alleles4)))

;;; Protein - repeated sequences

(def protein-repeated-seqs1s "Ala2[10]")
(def protein-repeated-seqs1ss "A2[10]")
(def protein-repeated-seqs1 (mut/protein-repeated-seqs "Ala" (coord/protein-coordinate 2)
                                                       nil nil
                                                       10))
(def protein-repeated-seqs1m {:mutation "protein-repeated-seqs"
                              :ref-start "Ala"
                              :coord-start (coord/plain (coord/protein-coordinate 2))
                              :ref-end nil
                              :coord-end nil
                              :ncopy 10})

(def protein-repeated-seqs2s "Arg65_Ser67[12]")
(def protein-repeated-seqs2ss "R65_S67[12]")
(def protein-repeated-seqs2 (mut/protein-repeated-seqs "Arg" (coord/protein-coordinate 65)
                                                       "Ser" (coord/protein-coordinate 67)
                                                       12))

(def protein-repeated-seqs3s "Gln18[(70_80)]")
(def protein-repeated-seqs3 (mut/protein-repeated-seqs "Gln" (coord/protein-coordinate 18)
                                                       nil nil
                                                       [70 80]))

(deftest format-protein-repeated-seqs-test
  (testing "returns a string expression of a protein repeated-seqs"
    (are [m o s] (= (mut/format m o) s)
      protein-repeated-seqs1 nil protein-repeated-seqs1s
      protein-repeated-seqs1 {:amino-acid-format :short} protein-repeated-seqs1ss
      protein-repeated-seqs2 nil protein-repeated-seqs2s
      protein-repeated-seqs2 {:amino-acid-format :short} protein-repeated-seqs2ss
      protein-repeated-seqs3 nil protein-repeated-seqs3s)))

(deftest parse-protein-repeated-seqs-test
  (testing "returns a correct ProteinRepeatedSeqs"
    (are [s m] (= (mut/parse-protein-repeated-seqs s) m)
      protein-repeated-seqs1s protein-repeated-seqs1
      protein-repeated-seqs1ss protein-repeated-seqs1
      protein-repeated-seqs2s protein-repeated-seqs2
      protein-repeated-seqs2ss protein-repeated-seqs2
      protein-repeated-seqs3s protein-repeated-seqs3))
  (testing "invalid protein repeated-seqs"
    (are [s] (thrown? #?(:clj Throwable, :cljs js/Error)
                      (mut/parse-protein-repeated-seqs s))
      "Ser67_Arg65[12]")))

(deftest plain-protein-repeated-seqs-test
  (testing "returns a plain map representing ProteinRepeatedSeqs"
    (is (= (mut/plain protein-repeated-seqs1) protein-repeated-seqs1m))))

(deftest restore-protein-repeated-seqs-test
  (testing "restores a plain map to ProteinRepeatedSeqs"
    (is (= (mut/restore protein-repeated-seqs1m) protein-repeated-seqs1))))

;;; Protein - frame shift

(def protein-frame-shift1s "Arg97ProfsTer23")
(def protein-frame-shift1sn "Arg97Profs")
(def protein-frame-shift1 (mut/protein-frame-shift "Arg"
                                                   (coord/protein-coordinate 97)
                                                   "Pro"
                                                   (coord/protein-coordinate 23)))
(def protein-frame-shift1m {:mutation "protein-frame-shift"
                            :ref "Arg"
                            :coord (coord/plain (coord/protein-coordinate 97))
                            :alt "Pro"
                            :new-ter-site (coord/plain (coord/protein-coordinate 23))})

(def protein-frame-shift2s "Arg97fs")
(def protein-frame-shift2ss "R97fs")
(def protein-frame-shift2 (mut/protein-frame-shift "Arg"
                                                   (coord/protein-coordinate 97)
                                                   nil
                                                   nil))

(def protein-frame-shift3s "Pro661=fs")
(def protein-frame-shift3 (mut/protein-frame-shift "Pro"
                                                   (coord/protein-coordinate 661)
                                                   "Pro"
                                                   nil))

(def protein-frame-shift4s "Ile327Argfs*?")
(def protein-frame-shift4 (mut/protein-frame-shift "Ile"
                                                   (coord/protein-coordinate 327)
                                                   "Arg"
                                                   (coord/unknown-coordinate)))

(def protein-frame-shift5s "Gln151Thrfs*9")
(def protein-frame-shift5 (mut/protein-frame-shift "Gln"
                                                   (coord/protein-coordinate 151)
                                                   "Thr"
                                                   (coord/protein-coordinate 9)))
(deftest format-protein-frame-shift-test
  (testing "returns a string expression of a protein frame-shift"
    (are [m o s] (= (mut/format m o) s)
      protein-frame-shift1 nil protein-frame-shift1sn
      protein-frame-shift1 {:show-ter-site? true} protein-frame-shift1s
      protein-frame-shift2 nil protein-frame-shift2s
      protein-frame-shift2 {:amino-acid-format :short} protein-frame-shift2ss
      protein-frame-shift3 nil protein-frame-shift3s
      protein-frame-shift4 {:show-ter-site? true, :ter-format :short} protein-frame-shift4s
      protein-frame-shift5 {:show-ter-site? true, :ter-format :short} protein-frame-shift5s)))

(deftest parse-protein-frame-shift-test
  (testing "returns a correct ProteinFrameShift"
    (are [s m] (= (mut/parse-protein-frame-shift s) m)
      protein-frame-shift1s protein-frame-shift1
      protein-frame-shift2s protein-frame-shift2
      protein-frame-shift2ss protein-frame-shift2
      protein-frame-shift3s protein-frame-shift3
      protein-frame-shift4s protein-frame-shift4
      protein-frame-shift5s protein-frame-shift5)))

(deftest plain-protein-frame-shift-test
  (testing "returns a plain map representing ProteinFrameShift"
    (is (= (mut/plain protein-frame-shift1) protein-frame-shift1m))))

(deftest restore-protein-frame-shift-test
  (testing "restores a plain map to ProteinFrameShift"
    (is (= (mut/restore protein-frame-shift1m) protein-frame-shift1))))

;;; Protein - extension

(def protein-extension1s "Met1ext-5")
(def protein-extension1ss "M1ext-5")
(def protein-extension1 (mut/protein-extension "Met"
                                               (coord/protein-coordinate 1)
                                               nil
                                               :upstream
                                               (coord/protein-coordinate 5)))
(def protein-extension1m {:mutation "protein-extension"
                          :ref "Met"
                          :coord (coord/plain (coord/protein-coordinate 1))
                          :alt nil
                          :region "upstream"
                          :new-site (coord/plain (coord/protein-coordinate 5))})

(def protein-extension2s "Met1Valext-12")
(def protein-extension2ss "M1Vext-12")
(def protein-extension2 (mut/protein-extension "Met"
                                               (coord/protein-coordinate 1)
                                               "Val"
                                               :upstream
                                               (coord/protein-coordinate 12)))

(def protein-extension3s "Ter110Glnext*17")
(def protein-extension3ss "*110Glnext*17")
(def protein-extension3 (mut/protein-extension "Ter"
                                               (coord/protein-coordinate 110)
                                               "Gln"
                                               :downstream
                                               (coord/protein-coordinate 17)))

(def protein-extension4s "Ter327Argext*?")
(def protein-extension4 (mut/protein-extension "Ter"
                                               (coord/protein-coordinate 327)
                                               "Arg"
                                               :downstream
                                               (coord/unknown-coordinate)))

(deftest format-protein-extension-test
  (testing "returns a string expression of a protein extension"
    (are [m o s] (= (mut/format m o) s)
      protein-extension1 nil protein-extension1s
      protein-extension1 {:amino-acid-format :short} protein-extension1ss
      protein-extension2 nil protein-extension2s
      protein-extension2 {:amino-acid-format :short} protein-extension2ss
      protein-extension3 nil protein-extension3s
      protein-extension3 {:ter-format :short} protein-extension3ss
      protein-extension4 nil protein-extension4s)))

(deftest parse-protein-extension-test
  (testing "returns a correct ProteinExtension"
    (are [s m] (= (mut/parse-protein-extension s) m)
      protein-extension1s protein-extension1
      protein-extension1ss protein-extension1
      protein-extension2s protein-extension2
      protein-extension2ss protein-extension2
      protein-extension3s protein-extension3
      protein-extension3ss protein-extension3
      protein-extension4s protein-extension4)))

(deftest plain-protein-extension-test
  (testing "returns a plain map representing ProteinExtension"
    (is (= (mut/plain protein-extension1) protein-extension1m))))

(deftest restore-protein-extension-test
  (testing "restores a plain map to ProteinExtension"
    (is (= (mut/restore protein-extension1m) protein-extension1))))

;;; Protein - general

(deftest parse-protein-test
  (are [s c] (instance? c (mut/parse-protein s))
    "Arg54Ser"                clj_hgvs.mutation.ProteinSubstitution
    "Cys76_Glu79del"          clj_hgvs.mutation.ProteinDeletion
    "Ala3_Ser5dup"            clj_hgvs.mutation.ProteinDuplication
    "Lys23_Leu24insArgSerGln" clj_hgvs.mutation.ProteinInsertion
    "Cys28_Lys29delinsTrp"    clj_hgvs.mutation.ProteinIndel
    "[Ser73Arg;Asn603del]"    clj_hgvs.mutation.ProteinAlleles
    "Arg65_Ser67[12]"         clj_hgvs.mutation.ProteinRepeatedSeqs
    "Arg97ProfsTer23"         clj_hgvs.mutation.ProteinFrameShift
    "Met1Valext-12"           clj_hgvs.mutation.ProteinExtension
    "(Arg2371Ser)"            clj_hgvs.mutation.UncertainMutation
    "0"                       clj_hgvs.mutation.NoProtein
    "?"                       clj_hgvs.mutation.ProteinUnknownMutation
    "Met1?"                   clj_hgvs.mutation.ProteinUnknownMutation
    "="                       clj_hgvs.mutation.ProteinNoEffect))

;;; others

(deftest equiv-test
  (testing "General"
    (are [mut1 mut2] (true? (mut/equiv mut1 mut2))
      (mut/parse-dna "2361G>A" :coding-dna) (mut/parse-dna "2361G>A" :coding-dna))
    (are [mut1 mut2] (false? (mut/equiv mut1 mut2))
      (mut/parse-dna "2361G>A" :coding-dna) (mut/parse-dna "2371G>A" :coding-dna)
      (mut/parse-dna "2361G>A" :coding-dna) (mut/parse-dna "2361del" :coding-dna)))

  (testing "DNA - deletion"
    (are [s1 s2] (true? (mut/equiv (mut/parse-dna-deletion s1 :genome)
                                   (mut/parse-dna-deletion s2 :genome)))
      "6_8del"    "6_8del"
      "6_8delTGC" "6_8delTGC"
      "6_8delTGC" "6_8del")
    (are [s1 s2] (false? (mut/equiv (mut/parse-dna-deletion s1 :genome)
                                    (mut/parse-dna-deletion s2 :genome)))
      "6_8del"    "6_9del"
      "6_8delTGC" "6_8delTGA"))

  (testing "DNA - duplication"
    (are [s1 s2] (true? (mut/equiv (mut/parse-dna-duplication s1 :genome)
                                   (mut/parse-dna-duplication s2 :genome)))
      "6_8dup"    "6_8dup"
      "6_8dupTGC" "6_8dupTGC"
      "6_8dupTGC" "6_8dup")
    (are [s1 s2] (false? (mut/equiv (mut/parse-dna-duplication s1 :genome)
                                    (mut/parse-dna-duplication s2 :genome)))
      "6_8dup"    "6_9dup"
      "6_8dupTGC" "6_8dupTGA"))

  (testing "DNA - indel"
    (are [s1 s2] (true? (mut/equiv (mut/parse-dna-indel s1 :genome)
                                   (mut/parse-dna-indel s2 :genome)))
      "6775delinsGA"  "6775delinsGA"
      "6775delTinsGA" "6775delTinsGA"
      "6775delTinsGA" "6775delinsGA")
    (are [s1 s2] (false? (mut/equiv (mut/parse-dna-indel s1 :genome)
                                    (mut/parse-dna-indel s2 :genome)))
      "6775delinsGA"  "6776delinsGA"
      "6775delTinsGA" "6775delCinsGA"))

  (testing "DNA - repeated sequences"
    (are [s1 s2] (true? (mut/equiv (mut/parse-dna-repeated-seqs s1 :genome)
                                   (mut/parse-dna-repeated-seqs s2 :genome)))
      "123_124[14]" "123_124[14]"
      "123TG[14]"   "123TG[14]"
      "123_124[14]" "123TG[14]"
      "123TG[14]"   "123_124[14]")
    (are [s1 s2] (false? (mut/equiv (mut/parse-dna-repeated-seqs s1 :genome)
                                    (mut/parse-dna-repeated-seqs s2 :genome)))
      "123_124[14]" "123_125[14]"
      "123TG[14]"   "123TC[14]"
      "123_124[14]" "123_124[15]"
      "123_124[14]" "123TGC[14]"
      "123TGC[14]"  "123_124[14]"))

  (testing "RNA - deletion"
    (are [s1 s2] (true? (mut/equiv (mut/parse-rna-deletion s1)
                                   (mut/parse-rna-deletion s2)))
      "6_8del"    "6_8del"
      "6_8delugc" "6_8delugc"
      "6_8delugc" "6_8del")
    (are [s1 s2] (false? (mut/equiv (mut/parse-rna-deletion s1)
                                    (mut/parse-rna-deletion s2)))
      "6_8del"    "6_9del"
      "6_8delugc" "6_8deluga"))

  (testing "RNA - duplication"
    (are [s1 s2] (true? (mut/equiv (mut/parse-rna-duplication s1)
                                   (mut/parse-rna-duplication s2)))
      "6_8dup"    "6_8dup"
      "6_8dupugc" "6_8dupugc"
      "6_8dupugc" "6_8dup")
    (are [s1 s2] (false? (mut/equiv (mut/parse-rna-duplication s1)
                                    (mut/parse-rna-duplication s2)))
      "6_8dup"    "6_9dup"
      "6_8dupugc" "6_8dupuga"))

  (testing "RNA - indel"
    (are [s1 s2] (true? (mut/equiv (mut/parse-rna-indel s1)
                                   (mut/parse-rna-indel s2)))
      "6775delinsga"  "6775delinsga"
      "6775deluinsga" "6775deluinsga"
      "6775deluinsga" "6775delinsga")
    (are [s1 s2] (false? (mut/equiv (mut/parse-rna-indel s1)
                                    (mut/parse-rna-indel s2)))
      "6775delinsga"  "6776delinsga"
      "6775deluinsga" "6775delcinsga"))

  (testing "RNA - repeated sequences"
    (are [s1 s2] (true? (mut/equiv (mut/parse-rna-repeated-seqs s1)
                                   (mut/parse-rna-repeated-seqs s2)))
      "-124_-123[14]" "-124_-123[14]"
      "-124ug[14]"    "-124ug[14]"
      "-124_-123[14]" "-124ug[14]"
      "-124ug[14]"    "-124_-123[14]")
    (are [s1 s2] (false? (mut/equiv (mut/parse-rna-repeated-seqs s1)
                                    (mut/parse-rna-repeated-seqs s2)))
      "-124_-123[14]" "-124_-122[14]"
      "-124ug[14]"    "-124uc[14]"
      "-124_-123[14]" "-124_-123[15]"
      "-124_-123[14]" "-124ugc[14]"
      "-124ugc[14]"   "-124_-123[14]"))

  (testing "Protein - frame shift"
    (are [mut1 mut2] (true? (mut/equiv mut1 mut2))
      (mut/parse-protein "K53Afs*9") (mut/parse-protein "K53Afs*9")
      (mut/parse-protein "K53Afs*9") (mut/parse-protein "K53fs*9")
      (mut/parse-protein "K53Afs*9") (mut/parse-protein "K53Afs")
      (mut/parse-protein "K53Afs*9") (mut/parse-protein "K53fs"))
    (are [mut1 mut2] (false? (mut/equiv mut1 mut2))
      (mut/parse-protein "K53Afs*9") (mut/parse-protein "K63Afs*9")
      (mut/parse-protein "K53Afs*9") (mut/parse-protein "K53Vfs*9")
      (mut/parse-protein "K53Afs*9") (mut/parse-protein "K53Afs*8")))

  (testing "Protein - extension"
    (are [mut1 mut2] (true? (mut/equiv mut1 mut2))
      (mut/parse-protein "Ter110Glnext*17") (mut/parse-protein "Ter110Glnext*17")
      (mut/parse-protein "Ter110Glnext*17") (mut/parse-protein "Ter110ext*17")
      (mut/parse-protein "Ter110Glnext*17") (mut/parse-protein "Ter110Glnext*?"))
    (are [mut1 mut2] (false? (mut/equiv mut1 mut2))
      (mut/parse-protein "Ter110Glnext*17") (mut/parse-protein "Ter120Glnext*17")
      (mut/parse-protein "Ter110Glnext*17") (mut/parse-protein "Ter110Gluext*17")
      (mut/parse-protein "Ter110Glnext*17") (mut/parse-protein "Ter110Glnext*18")))

  (testing "Uncertain mutation"
    (are [mut1 mut2] (true? (mut/equiv mut1 mut2))
      (mut/parse-rna "(?)") (mut/parse-rna "(?)")
      (mut/parse-protein "(K53Afs*9)") (mut/parse-protein "(K53fs)"))
    (are [mut1 mut2] (false? (mut/equiv mut1 mut2))
      (mut/parse-rna "(?)") (mut/parse-rna "(306g>u)")
      (mut/parse-protein "(K53Afs*9)") (mut/parse-protein "(K53Afs*8)"))))
