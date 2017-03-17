(ns clj-hgvs.mutation-test
  (:require #?(:clj [clojure.test :refer :all]
               :cljs [cljs.test :refer-macros [deftest is are testing]])
            [clj-hgvs.coordinate :as coord]
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

;;; DNA mutations

;;; DNA - substitution

(def dna-substitution1s "45576A>C")
(def dna-substitution1k :genome)
(def dna-substitution1 (mut/map->DNASubstitution {:coord (coord/genomic-coordinate 45576)
                                                  :ref "A"
                                                  :type ">"
                                                  :alt "C"}))

(def dna-substitution2s "88+1G>T")
(def dna-substitution2k :cdna)
(def dna-substitution2 (mut/map->DNASubstitution {:coord (coord/cdna-coordinate 88 1 nil)
                                                  :ref "G"
                                                  :type ">"
                                                  :alt "T"}))

(def dna-substitution3s "123G=")
(def dna-substitution3k :cdna)
(def dna-substitution3 (mut/map->DNASubstitution {:coord (coord/cdna-coordinate 123)
                                                  :ref "G"
                                                  :type "="
                                                  :alt nil}))

(def dna-substitution4s "85C=/>T")
(def dna-substitution4k :cdna)
(def dna-substitution4 (mut/map->DNASubstitution {:coord (coord/cdna-coordinate 85)
                                                  :ref "C"
                                                  :type "=/>"
                                                  :alt "T"}))

(def dna-substitution5s "85C=//>T")
(def dna-substitution5k :cdna)
(def dna-substitution5 (mut/map->DNASubstitution {:coord (coord/cdna-coordinate 85)
                                                  :ref "C"
                                                  :type "=//>"
                                                  :alt "T"}))

(deftest format-dna-substitution-test
  (testing "returns a string expression of a DNA substitution"
    (are [m s] (= (mut/format m nil) s)
      dna-substitution1 dna-substitution1s
      dna-substitution2 dna-substitution2s
      dna-substitution3 dna-substitution3s
      dna-substitution4 dna-substitution4s
      dna-substitution5 dna-substitution5s)))

(deftest parse-dna-substitution-test
  (testing "returns a correct DNASubstitution"
    (are [s k m] (= (mut/parse-dna-substitution s k) m)
      dna-substitution1s dna-substitution1k dna-substitution1
      dna-substitution2s dna-substitution2k dna-substitution2
      dna-substitution3s dna-substitution3k dna-substitution3
      dna-substitution4s dna-substitution4k dna-substitution4
      dna-substitution5s dna-substitution5k dna-substitution5)))

;;; DNA - deletion

(def dna-deletion1s "7del")
(def dna-deletion1k :genome)
(def dna-deletion1 (mut/map->DNADeletion {:coord-start (coord/genomic-coordinate 7)
                                          :coord-end nil
                                          :ref nil}))

(def dna-deletion2s "6_8del")
(def dna-deletion2k :genome)
(def dna-deletion2 (mut/map->DNADeletion {:coord-start (coord/genomic-coordinate 6)
                                          :coord-end (coord/genomic-coordinate 8)
                                          :ref nil}))

(def dna-deletion3s "6_8delTGC")
(def dna-deletion3ss "6_8del")
(def dna-deletion3k :genome)
(def dna-deletion3 (mut/map->DNADeletion {:coord-start (coord/genomic-coordinate 6)
                                          :coord-end (coord/genomic-coordinate 8)
                                          :ref "TGC"}))

(def dna-deletion4s "120_123+48del")
(def dna-deletion4k :cdna)
(def dna-deletion4 (mut/map->DNADeletion {:coord-start (coord/cdna-coordinate 120)
                                          :coord-end (coord/cdna-coordinate 123 48 nil)
                                          :ref nil}))

(def dna-deletion5s "(4071+1_4072-1)_(5145+1_5146-1)del")
(def dna-deletion5k :cdna)
(def dna-deletion5 (mut/map->DNADeletion {:coord-start [(coord/cdna-coordinate 4071 1 nil)
                                                        (coord/cdna-coordinate 4072 -1 nil)]
                                          :coord-end [(coord/cdna-coordinate 5145 1 nil)
                                                      (coord/cdna-coordinate 5146 -1 nil)]
                                          :ref nil}))

(def dna-deletion6s "(?_-30)_(12+1_13-1)del")
(def dna-deletion6k :cdna)
(def dna-deletion6 (mut/map->DNADeletion {:coord-start [(coord/unknown-coordinate)
                                                        (coord/cdna-coordinate 30 0 :upstream)]
                                          :coord-end [(coord/cdna-coordinate 12 1 nil)
                                                      (coord/cdna-coordinate 13 -1 nil)]
                                          :ref nil}))

(def dna-deletion7s "(?_-1)_(*1_?)del")
(def dna-deletion7k :cdna)
(def dna-deletion7 (mut/map->DNADeletion {:coord-start [(coord/unknown-coordinate)
                                                        (coord/cdna-coordinate 1 0 :upstream)]
                                          :coord-end [(coord/cdna-coordinate 1 0 :downstream)
                                                      (coord/unknown-coordinate)]
                                          :ref nil}))

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
      dna-deletion7s dna-deletion7k dna-deletion7)))

;;; DNA - duplication

(def dna-duplication1s "7dup")
(def dna-duplication1k :genome)
(def dna-duplication1 (mut/map->DNADuplication {:coord-start (coord/genomic-coordinate 7)
                                                :coord-end nil
                                                :ref nil}))

(def dna-duplication2s "6_8dup")
(def dna-duplication2k :genome)
(def dna-duplication2 (mut/map->DNADuplication {:coord-start (coord/genomic-coordinate 6)
                                                :coord-end (coord/genomic-coordinate 8)
                                                :ref nil}))

(def dna-duplication3s "6_8dupTGC")
(def dna-duplication3ss "6_8dup")
(def dna-duplication3k :genome)
(def dna-duplication3 (mut/map->DNADuplication {:coord-start (coord/genomic-coordinate 6)
                                                :coord-end (coord/genomic-coordinate 8)
                                                :ref "TGC"}))


(def dna-duplication4s "120_123+48dup")
(def dna-duplication4k :cdna)
(def dna-duplication4 (mut/map->DNADuplication {:coord-start (coord/cdna-coordinate 120)
                                                :coord-end (coord/cdna-coordinate 123 48 nil)
                                                :ref nil}))

(def dna-duplication5s "(4071+1_4072-1)_(5145+1_5146-1)dup")
(def dna-duplication5k :cdna)
(def dna-duplication5 (mut/map->DNADuplication {:coord-start [(coord/cdna-coordinate 4071 1 nil)
                                                              (coord/cdna-coordinate 4072 -1 nil)]
                                                :coord-end [(coord/cdna-coordinate 5145 1 nil)
                                                            (coord/cdna-coordinate 5146 -1 nil)]
                                                :ref nil}))

(def dna-duplication6s "(?_-30)_(12+1_13-1)dup")
(def dna-duplication6k :cdna)
(def dna-duplication6 (mut/map->DNADuplication {:coord-start [(coord/unknown-coordinate)
                                                              (coord/cdna-coordinate 30 0 :upstream)]
                                                :coord-end [(coord/cdna-coordinate 12 1 nil)
                                                            (coord/cdna-coordinate 13 -1 nil)]
                                                :ref nil}))

(def dna-duplication7s "(?_-1)_(*1_?)dup")
(def dna-duplication7k :cdna)
(def dna-duplication7 (mut/map->DNADuplication {:coord-start [(coord/unknown-coordinate)
                                                              (coord/cdna-coordinate 1 0 :upstream)]
                                                :coord-end [(coord/cdna-coordinate 1 0 :downstream)
                                                            (coord/unknown-coordinate)]
                                                :ref nil}))

(deftest format-dna-duplication-test
  (testing "returns a string expression of a DNA duplication"
    (are [m o s] (= (mut/format m o ) s)
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
      dna-duplication7s dna-duplication7k dna-duplication7)))

;;; DNA - insertion

(def dna-insertion1s "5756_5757insAGG")
(def dna-insertion1k :genome)
(def dna-insertion1 (mut/map->DNAInsertion {:coord-start (coord/genomic-coordinate 5756)
                                            :coord-end (coord/genomic-coordinate 5757)
                                            :alt "AGG"}))

(def dna-insertion2s "123_124insL37425.1:23_361")
(def dna-insertion2k :genome)
(def dna-insertion2 (mut/map->DNAInsertion {:coord-start (coord/genomic-coordinate 123)
                                            :coord-end (coord/genomic-coordinate 124)
                                            :alt {:transcript "L37425.1"
                                                  :coord-start (coord/genomic-coordinate 23)
                                                  :coord-end (coord/genomic-coordinate 361)}}))

(def dna-insertion3s "122_123ins123_234inv")
(def dna-insertion3k :genome)
(def dna-insertion3 "TODO")

(def dna-insertion4s "122_123ins213_234invinsAins123_211inv")
(def dna-insertion4k :genome)
(def dna-insertion4 "TODO")

(deftest format-dna-insertion-test
  (testing "returns a string expression of a DNA insertion"
    (are [m s] (= (mut/format m nil) s)
      dna-insertion1 dna-insertion1s
      dna-insertion2 dna-insertion2s
      ;; dna-insertion3 dna-insertion3s ; TODO
      ;; dna-insertion4 dna-insertion4s ; TODO
      )))

(deftest parse-dna-insertion-test
  (testing "returns a correct DNAInsertion"
    (are [s k m] (= (mut/parse-dna-insertion s k) m)
      dna-insertion1s dna-insertion1k dna-insertion1
      dna-insertion2s dna-insertion2k dna-insertion2
      ;; dna-insertion3s dna-insertion3k dna-insertion3 ; TODO
      ;; dna-insertion4s dna-insertion4k dna-insertion4 ; TODO
      )))

;;; DNA - inversion

(def dna-inversion1s "1077_1080inv")
(def dna-inversion1k :genome)
(def dna-inversion1 (mut/map->DNAInversion {:coord-start (coord/genomic-coordinate 1077)
                                            :coord-end (coord/genomic-coordinate 1080)}))

(def dna-inversion2s "77_80inv")
(def dna-inversion2k :cdna)
(def dna-inversion2 (mut/map->DNAInversion {:coord-start (coord/cdna-coordinate 77)
                                            :coord-end (coord/cdna-coordinate 80)}))

(deftest format-dna-inversion-test
  (testing "returns a string expression of a DNA inversion"
    (are [m s] (= (mut/format m nil) s)
      dna-inversion1 dna-inversion1s
      dna-inversion2 dna-inversion2s)))

(deftest parse-dna-inversion-test
  (testing "returns a correct DNAInversion"
    (are [s k m] (= (mut/parse-dna-inversion s k) m)
      dna-inversion1s dna-inversion1k dna-inversion1
      dna-inversion2s dna-inversion2k dna-inversion2)))


;;; DNA - conversion

(def dna-conversion1s "333_590con1844_2101")
(def dna-conversion1k :genome)
(def dna-conversion1 (mut/map->DNAConversion {:coord-start (coord/genomic-coordinate 333)
                                              :coord-end (coord/genomic-coordinate 590)
                                              :alt {:transcript nil
                                                    :kind nil
                                                    :coord-start (coord/genomic-coordinate 1844)
                                                    :coord-end (coord/genomic-coordinate 2101)}}))

(def dna-conversion2s "415_1655conAC096506.5:g.409_1683")
(def dna-conversion2k :genome)
(def dna-conversion2 (mut/map->DNAConversion {:coord-start (coord/genomic-coordinate 415)
                                              :coord-end (coord/genomic-coordinate 1655)
                                              :alt {:transcript "AC096506.5"
                                                    :kind :genome
                                                    :coord-start (coord/genomic-coordinate 409)
                                                    :coord-end (coord/genomic-coordinate 1683)}}))

(def dna-conversion3s "15_355conNM_004006.1:20_360")
(def dna-conversion3k :cdna)
(def dna-conversion3 (mut/map->DNAConversion {:coord-start (coord/cdna-coordinate 15)
                                              :coord-end (coord/cdna-coordinate 355)
                                              :alt {:transcript "NM_004006.1"
                                                    :kind nil
                                                    :coord-start (coord/cdna-coordinate 20)
                                                    :coord-end (coord/cdna-coordinate 360)}}))

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
      dna-conversion3s dna-conversion3k dna-conversion3)))

;;; DNA - indel

(def dna-indel1s "6775delinsGA")
(def dna-indel1k :genome)
(def dna-indel1 (mut/map->DNAIndel {:coord-start (coord/genomic-coordinate 6775)
                                    :coord-end nil
                                    :alt "GA"}))

(def dna-indel2s "145_147delinsTGG")
(def dna-indel2k :cdna)
(def dna-indel2 (mut/map->DNAIndel {:coord-start (coord/cdna-coordinate 145)
                                    :coord-end (coord/cdna-coordinate 147)
                                    :alt "TGG"}))

(deftest format-dna-indel-test
  (testing "returns a string expression of a DNA indel"
    (are [m s] (= (mut/format m nil) s)
      dna-indel1 dna-indel1s
      dna-indel2 dna-indel2s)))

(deftest parse-dna-indel-test
  (testing "returns a correct DNAIndel"
    (are [s k m] (= (mut/parse-dna-indel s k) m)
      dna-indel1s dna-indel1k dna-indel1
      dna-indel2s dna-indel2k dna-indel2)))

;;; DNA - repeated sequences

(def dna-repeated-seqss "123_124[14]")
(def dna-repeated-seqsk :genome)
(def dna-repeated-seqs (mut/map->DNARepeatedSeqs {:coord-start (coord/genomic-coordinate 123)
                                                  :coord-end (coord/genomic-coordinate 124)
                                                  :ncopy 14}))

(deftest format-dna-repeated-seqs-test
  (testing "returns a string expression of a DNA repeated sequences"
    (are [m s] (= (mut/format m nil) s)
      dna-repeated-seqs dna-repeated-seqss)))

(deftest parse-dna-repeated-seqs-test
  (testing "returns a correct DNARepeatedSeqs"
    (are [s k m] (= (mut/parse-dna-repeated-seqs s k) m)
      dna-repeated-seqss dna-repeated-seqsk dna-repeated-seqs)))

;;; RNA mutations

;;; RNA - substitution

(def rna-substitution1s "76a>c")
(def rna-substitution1 (mut/map->RNASubstitution {:coord (coord/rna-coordinate 76 nil nil)
                                                  :ref "a"
                                                  :alt "c"}))

(def rna-substitution2s "-14g>c")
(def rna-substitution2 (mut/map->RNASubstitution {:coord (coord/rna-coordinate 14 0 :upstream)
                                                  :ref "g"
                                                  :alt "c"}))

(def rna-substitution3s "*46u>a")
(def rna-substitution3 (mut/map->RNASubstitution {:coord (coord/rna-coordinate 46 0 :downstream)
                                                  :ref "u"
                                                  :alt "a"}))

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

;;; RNA - deletion

(def rna-deletion1s "7del")
(def rna-deletion1 (mut/map->RNADeletion {:coord-start (coord/rna-coordinate 7 nil nil)
                                          :coord-end nil
                                          :ref nil}))

(def rna-deletion2s "7delu")
(def rna-deletion2ss "7del")
(def rna-deletion2 (mut/map->RNADeletion {:coord-start (coord/rna-coordinate 7 nil nil)
                                          :coord-end nil
                                          :ref "u"}))

(def rna-deletion3s "6_8del")
(def rna-deletion3 (mut/map->RNADeletion {:coord-start (coord/rna-coordinate 6 nil nil)
                                          :coord-end (coord/rna-coordinate 8 nil nil)
                                          :ref nil}))

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
      )))

;;; RNA - duplication

(def rna-duplication1s "7dup")
(def rna-duplication1 (mut/map->RNADuplication {:coord-start (coord/rna-coordinate 7 nil nil)
                                                :coord-end nil
                                                :ref nil}))

(def rna-duplication2s "7dupu")
(def rna-duplication2ss "7dup")
(def rna-duplication2 (mut/map->RNADuplication {:coord-start (coord/rna-coordinate 7 nil nil)
                                                :coord-end nil
                                                :ref "u"}))

(def rna-duplication3s "6_8dup")
(def rna-duplication3 (mut/map->RNADuplication {:coord-start (coord/rna-coordinate 6 nil nil)
                                                :coord-end (coord/rna-coordinate 8 nil nil)
                                                :ref nil}))

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
      rna-duplication3s rna-duplication3)))

;;; RNA - insertion

(def rna-insertion1s "756_757insacu")
(def rna-insertion1 (mut/map->RNAInsertion {:coord-start (coord/rna-coordinate 756 nil nil)
                                            :coord-end (coord/rna-coordinate 757 nil nil)
                                            :alt "acu"}))

(def rna-insertion2s "431_432ins(5)")
(def rna-insertion2 (mut/map->RNAInsertion {:coord-start (coord/rna-coordinate 431 nil nil)
                                            :coord-end (coord/rna-coordinate 432 nil nil)
                                            :alt "nnnnn"}))

(def rna-insertion3s "123_124insL37425.1:23_361")
(def rna-insertion3 (mut/map->RNAInsertion {:coord-start (coord/rna-coordinate 123 nil nil)
                                            :coord-end (coord/rna-coordinate 124 nil nil)
                                            :alt {:genbank "L37425.1"
                                                  :coord-start 23
                                                  :coord-end 361}}))

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
      rna-insertion3s rna-insertion3)))

;;; RNA - inversion

(def rna-inversion1s "177_180inv")
(def rna-inversion1 (mut/map->RNAInversion {:coord-start (coord/rna-coordinate 177 nil nil)
                                            :coord-end (coord/rna-coordinate 180 nil nil)}))

(deftest format-rna-inversion-test
  (testing "returns a string expression of a RNA inversion"
    (are [m s] (= (mut/format m nil) s)
      rna-inversion1 rna-inversion1s)))

(deftest parse-rna-inversion-test
  (testing "returns a correct RNAInversion"
    (are [s m] (= (mut/parse-rna-inversion s) m)
      rna-inversion1s rna-inversion1)))

;;; RNA - conversion

(def rna-conversion1s "123_345con888_1110")
(def rna-conversion1 (mut/map->RNAConversion {:coord-start (coord/rna-coordinate 123 nil nil)
                                              :coord-end (coord/rna-coordinate 345 nil nil)
                                              :alt {:transcript nil
                                                    :coord-start (coord/rna-coordinate 888 nil nil)
                                                    :coord-end (coord/rna-coordinate 1110 nil nil)}}))

(def rna-conversion2s "415_1655conAC096506.5:409_1649")
(def rna-conversion2 (mut/map->RNAConversion {:coord-start (coord/rna-coordinate 415 nil nil)
                                              :coord-end (coord/rna-coordinate 1655 nil nil)
                                              :alt {:transcript "AC096506.5"
                                                    :coord-start (coord/rna-coordinate 409 nil nil)
                                                    :coord-end (coord/rna-coordinate 1649 nil nil)}}))

(deftest format-rna-conversion-test
  (testing "returns a string expression of a RNA conversion"
    (are [m s] (= (mut/format m nil) s)
      rna-conversion1 rna-conversion1s
      rna-conversion2 rna-conversion2s)))

(deftest parse-rna-conversion-test
  (testing "returns a correct RNAConversion"
    (are [s m] (= (mut/parse-rna-conversion s) m)
      rna-conversion1s rna-conversion1
      rna-conversion2s rna-conversion2)))

;;; RNA - indel

(def rna-indel1s "775delinsga")
(def rna-indel1 (mut/map->RNAIndel {:coord-start (coord/rna-coordinate 775 nil nil)
                                    :coord-end nil
                                    :alt "ga"}))

(def rna-indel2s "775_777delinsc")
(def rna-indel2 (mut/map->RNAIndel {:coord-start (coord/rna-coordinate 775 nil nil)
                                    :coord-end (coord/rna-coordinate 777 nil nil)
                                    :alt "c"}))

(deftest format-rna-indel-test
  (testing "returns a string expression of a RNA indel"
    (are [m s] (= (mut/format m nil) s)
      rna-indel1 rna-indel1s
      rna-indel2 rna-indel2s)))

(deftest parse-rna-indel-test
  (testing "returns a correct RNAIndel"
    (are [s m] (= (mut/parse-rna-indel s) m)
      rna-indel1s rna-indel1
      rna-indel2s rna-indel2)))

;;; RNA - repeated sequences

(def rna-repeated-seqs1s "-124_-123[14]")
(def rna-repeated-seqs1 (mut/map->RNARepeatedSeqs {:coord-start (coord/rna-coordinate 124 0 :upstream)
                                                   :coord-end (coord/rna-coordinate 123 0 :upstream)
                                                   :ref nil
                                                   :ncopy 14
                                                   :ncopy-other nil}))

(def rna-repeated-seqs2s "-124ug[14]")
(def rna-repeated-seqs2 (mut/map->RNARepeatedSeqs {:coord-start (coord/rna-coordinate 124 0 :upstream)
                                                   :coord-end nil
                                                   :ref "ug"
                                                   :ncopy 14
                                                   :ncopy-other nil}))

(def rna-repeated-seqs3s "-124_-123[14];[18]")
(def rna-repeated-seqs3 (mut/map->RNARepeatedSeqs {:coord-start (coord/rna-coordinate 124 0 :upstream)
                                                   :coord-end (coord/rna-coordinate 123 0 :upstream)
                                                   :ref nil
                                                   :ncopy 14
                                                   :ncopy-other 18}))

(deftest format-rna-repeated-seqs-test
  (testing "returns a string expression of a RNA repeated-seqs"
    (are [m s] (= (mut/format m nil) s)
      rna-repeated-seqs1 rna-repeated-seqs1s
      rna-repeated-seqs2 rna-repeated-seqs2s
      rna-repeated-seqs3 rna-repeated-seqs3s)))

(deftest parse-rna-repeated-seqs-test
  (testing "returns a correct RNARepeatedSeqs"
    (are [s m] (= (mut/parse-rna-repeated-seqs s) m)
      rna-repeated-seqs1s rna-repeated-seqs1
      rna-repeated-seqs2s rna-repeated-seqs2
      rna-repeated-seqs3s rna-repeated-seqs3)))

;;; Protein mutations

;;; Protein - substitution

(def protein-substitution1s "Arg54Ser")
(def protein-substitution1ss "R54S")
(def protein-substitution1 (mut/map->ProteinSubstitution {:ref "Arg"
                                                          :coord (coord/protein-coordinate 54)
                                                          :alt "Ser"}))

(def protein-substitution2s "Cys123=")
(def protein-substitution2ss "C123=")
(def protein-substitution2 (mut/map->ProteinSubstitution {:ref "Cys"
                                                          :coord (coord/protein-coordinate 123)
                                                          :alt "Cys"}))

(deftest format-protein-substitution-test
  (testing "returns a string expression of a protein substitution"
    (are [m o s] (= (mut/format m o) s)
      protein-substitution1 nil protein-substitution1s
      protein-substitution1 {:amino-acid-format :short} protein-substitution1ss
      protein-substitution2 nil protein-substitution2s
      protein-substitution2 {:amino-acid-format :short} protein-substitution2ss)))

(deftest parse-protein-substitution-test
  (testing "returns a correct ProteinSubstitution"
    (are [s m] (= (mut/parse-protein-substitution s) m)
      protein-substitution1s protein-substitution1
      protein-substitution1ss protein-substitution1
      protein-substitution2s protein-substitution2
      protein-substitution2ss protein-substitution2)))

;;; Protein - deletion

(def protein-deletion1s "Ala3del")
(def protein-deletion1ss "A3del")
(def protein-deletion1 (mut/map->ProteinDeletion {:ref-start "Ala"
                                                  :coord-start (coord/protein-coordinate 3)
                                                  :ref-end nil
                                                  :coord-end nil}))

(def protein-deletion2s "Cys76_Glu79del")
(def protein-deletion2ss "C76_E79del")
(def protein-deletion2 (mut/map->ProteinDeletion {:ref-start "Cys"
                                                  :coord-start (coord/protein-coordinate 76)
                                                  :ref-end "Glu"
                                                  :coord-end (coord/protein-coordinate 79)}))
(deftest format-protein-deletion-test
  (testing "returns a string expression of a protein deletion"
    (are [m o s] (= (mut/format m o) s)
      protein-deletion1 nil protein-deletion1s
      protein-deletion1 {:amino-acid-format :short} protein-deletion1ss
      protein-deletion2 nil protein-deletion2s
      protein-deletion2 {:amino-acid-format :short} protein-deletion2ss))
  (testing "not show last amino acid if range size is 1."
    (is (= (mut/format (mut/map->ProteinDeletion
                        {:ref-start "Ala"
                         :coord-start (coord/protein-coordinate 3)
                         :ref-end "Ala"
                         :coord-end (coord/protein-coordinate 3)}) nil)
           "Ala3del"))))

(deftest parse-protein-deletion-test
  (testing "returns a correct ProteinDeletion"
    (are [s m] (= (mut/parse-protein-deletion s) m)
      protein-deletion1s protein-deletion1
      protein-deletion1ss protein-deletion1
      protein-deletion2s protein-deletion2
      protein-deletion2ss protein-deletion2)))

;;; Protein - duplication

(def protein-duplication1s "Ala3dup")
(def protein-duplication1ss "A3dup")
(def protein-duplication1 (mut/map->ProteinDuplication {:ref-start "Ala"
                                                        :coord-start (coord/protein-coordinate 3)
                                                        :ref-end nil
                                                        :coord-end nil}))

(def protein-duplication2s "Ala3_Ser5dup")
(def protein-duplication2ss "A3_S5dup")
(def protein-duplication2 (mut/map->ProteinDuplication {:ref-start "Ala"
                                                        :coord-start (coord/protein-coordinate 3)
                                                        :ref-end "Ser"
                                                        :coord-end (coord/protein-coordinate 5)}))

(deftest format-protein-duplication-test
  (testing "returns a string expression of a protein duplication"
    (are [m o s] (= (mut/format m o) s)
      protein-duplication1 nil protein-duplication1s
      protein-duplication1 {:amino-acid-format :short} protein-duplication1ss
      protein-duplication2 nil protein-duplication2s
      protein-duplication2 {:amino-acid-format :short} protein-duplication2ss))
  (testing "not show last amino acid if range size is 1."
    (is (= (mut/format (mut/map->ProteinDuplication
                        {:ref-start "Ala"
                         :coord-start (coord/protein-coordinate 3)
                         :ref-end "Ala"
                         :coord-end (coord/protein-coordinate 3)}) nil)
           "Ala3dup"))))

(deftest parse-protein-duplication-test
  (testing "returns a correct ProteinDuplication"
    (are [s m] (= (mut/parse-protein-duplication s) m)
      protein-duplication1s protein-duplication1
      protein-duplication1ss protein-duplication1
      protein-duplication2s protein-duplication2
      protein-duplication2ss protein-duplication2)))

;;; Protein - insertion

(def protein-insertions "Lys23_Leu24insArgSerGln")
(def protein-insertionss "K23_L24insRSQ")
(def protein-insertion (mut/map->ProteinInsertion {:ref-start "Lys"
                                                   :coord-start (coord/protein-coordinate 23)
                                                   :ref-end "Leu"
                                                   :coord-end (coord/protein-coordinate 24)
                                                   :alts ["Arg" "Ser" "Gln"]}))

(deftest format-protein-insertion-test
  (testing "returns a string expression of a protein insertion"
    (are [m o s] (= (mut/format m o) s)
      protein-insertion nil protein-insertions
      protein-insertion {:amino-acid-format :short} protein-insertionss)))

(deftest parse-protein-insertion-test
  (testing "returns a correct ProteinInsertion"
    (are [s m] (= (mut/parse-protein-insertion s) m)
      protein-insertions protein-insertion
      protein-insertionss protein-insertion)))

;;; Protein - indel

(def protein-indel1s "Cys28delinsTrpVal")
(def protein-indel1ss "C28delinsWV")
(def protein-indel1 (mut/map->ProteinIndel {:ref-start "Cys"
                                            :coord-start (coord/protein-coordinate 28)
                                            :ref-end nil
                                            :coord-end nil
                                            :alts ["Trp" "Val"]}))

(def protein-indel2s "Cys28_Lys29delinsTrp")
(def protein-indel2ss "C28_K29delinsW")
(def protein-indel2 (mut/map->ProteinIndel {:ref-start "Cys"
                                            :coord-start (coord/protein-coordinate 28)
                                            :ref-end "Lys"
                                            :coord-end (coord/protein-coordinate 29)
                                            :alts ["Trp"]}))

(deftest format-protein-indel-test
  (testing "returns a string expression of a protein indel"
    (are [m o s] (= (mut/format m o) s)
      protein-indel1 nil protein-indel1s
      protein-indel1 {:amino-acid-format :short} protein-indel1ss
      protein-indel2 nil protein-indel2s
      protein-indel2 {:amino-acid-format :short} protein-indel2ss)))

(deftest parse-protein-indel-test
  (testing "returns a correct ProteinIndel"
    (are [s m] (= (mut/parse-protein-indel s) m)
      protein-indel1s protein-indel1
      protein-indel1ss protein-indel1
      protein-indel2s protein-indel2
      protein-indel2ss protein-indel2)))

;;; Protein - repeated sequences

(def protein-repeated-seqs1s "Ala2[10]")
(def protein-repeated-seqs1ss "A2[10]")
(def protein-repeated-seqs1 (mut/map->ProteinRepeatedSeqs {:ref-start "Ala"
                                                           :coord-start (coord/protein-coordinate 2)
                                                           :ref-end nil
                                                           :coord-end nil
                                                           :ncopy 10
                                                           :ncopy-other nil}))

(def protein-repeated-seqs2s "Ala2[10];[11]")
(def protein-repeated-seqs2ss "A2[10];[11]")
(def protein-repeated-seqs2 (mut/map->ProteinRepeatedSeqs {:ref-start "Ala"
                                                           :coord-start (coord/protein-coordinate 2)
                                                           :ref-end nil
                                                           :coord-end nil
                                                           :ncopy 10
                                                           :ncopy-other 11}))

(def protein-repeated-seqs3s "Arg65_Ser67[12]")
(def protein-repeated-seqs3ss "R65_S67[12]")
(def protein-repeated-seqs3 (mut/map->ProteinRepeatedSeqs {:ref-start "Arg"
                                                           :coord-start (coord/protein-coordinate 65)
                                                           :ref-end "Ser"
                                                           :coord-end (coord/protein-coordinate 67)
                                                           :ncopy 12
                                                           :ncopy-other nil}))

(deftest format-protein-repeated-seqs-test
  (testing "returns a string expression of a protein repeated-seqs"
    (are [m o s] (= (mut/format m o) s)
      protein-repeated-seqs1 nil protein-repeated-seqs1s
      protein-repeated-seqs1 {:amino-acid-format :short} protein-repeated-seqs1ss
      protein-repeated-seqs2 nil protein-repeated-seqs2s
      protein-repeated-seqs2 {:amino-acid-format :short} protein-repeated-seqs2ss
      protein-repeated-seqs3 nil protein-repeated-seqs3s
      protein-repeated-seqs3 {:amino-acid-format :short} protein-repeated-seqs3ss)))

(deftest parse-protein-repeated-seqs-test
  (testing "returns a correct ProteinRepeatedSeqs"
    (are [s m] (= (mut/parse-protein-repeated-seqs s) m)
      protein-repeated-seqs1s protein-repeated-seqs1
      protein-repeated-seqs1ss protein-repeated-seqs1
      protein-repeated-seqs2s protein-repeated-seqs2
      protein-repeated-seqs2ss protein-repeated-seqs2
      protein-repeated-seqs3s protein-repeated-seqs3
      protein-repeated-seqs3ss protein-repeated-seqs3)))

;;; Protein - frame shift

(def protein-frame-shift1s "Arg97ProfsTer23")
(def protein-frame-shift1 (mut/map->ProteinFrameShift {:ref "Arg"
                                                       :coord (coord/protein-coordinate 97)
                                                       :alt "Pro"
                                                       :new-ter-site (coord/protein-coordinate 23)}))

(def protein-frame-shift2s "Arg97fs")
(def protein-frame-shift2ss "R97fs")
(def protein-frame-shift2 (mut/map->ProteinFrameShift {:ref "Arg"
                                                       :coord (coord/protein-coordinate 97)
                                                       :alt nil
                                                       :new-ter-site nil}))

(def protein-frame-shift3s "Ile327Argfs*?")
(def protein-frame-shift3 (mut/map->ProteinFrameShift {:ref "Ile"
                                                       :coord (coord/protein-coordinate 327)
                                                       :alt "Arg"
                                                       :new-ter-site (coord/unknown-coordinate)}))

(def protein-frame-shift4s "Gln151Thrfs*9")
(def protein-frame-shift4 (mut/map->ProteinFrameShift {:ref "Gln"
                                                       :coord (coord/protein-coordinate 151)
                                                       :alt "Thr"
                                                       :new-ter-site (coord/protein-coordinate 9)}))
(deftest format-protein-frame-shift-test
  (testing "returns a string expression of a protein frame-shift"
    (are [m o s] (= (mut/format m o) s)
      protein-frame-shift1 nil protein-frame-shift1s
      protein-frame-shift2 nil protein-frame-shift2s
      protein-frame-shift2 {:amino-acid-format :short} protein-frame-shift2ss
      protein-frame-shift3 {:ter-format :short} protein-frame-shift3s
      protein-frame-shift4 {:ter-format :short} protein-frame-shift4s)))

(deftest parse-protein-frame-shift-test
  (testing "returns a correct ProteinFrameShift"
    (are [s m] (= (mut/parse-protein-frame-shift s) m)
      protein-frame-shift1s protein-frame-shift1
      protein-frame-shift2s protein-frame-shift2
      protein-frame-shift2ss protein-frame-shift2
      protein-frame-shift3s protein-frame-shift3
      protein-frame-shift4s protein-frame-shift4)))

;;; Protein - extension

(def protein-extension1s "Met1ext-5")
(def protein-extension1ss "M1ext-5")
(def protein-extension1 (mut/map->ProteinExtension {:ref "Met"
                                                    :coord (coord/protein-coordinate 1)
                                                    :alt nil
                                                    :new-site "-5"}))

(def protein-extension2s "Met1Valext-12")
(def protein-extension2ss "M1Vext-12")
(def protein-extension2 (mut/map->ProteinExtension {:ref "Met"
                                                    :coord (coord/protein-coordinate 1)
                                                    :alt "Val"
                                                    :new-site "-12"}))

(def protein-extension3s "Ter110GlnextTer17")
(def protein-extension3 (mut/map->ProteinExtension {:ref "Ter"
                                                    :coord (coord/protein-coordinate 110)
                                                    :alt "Gln"
                                                    :new-site "Ter17"}))

(deftest format-protein-extension-test
  (testing "returns a string expression of a protein extension"
    (are [m o s] (= (mut/format m o) s)
      protein-extension1 nil protein-extension1s
      protein-extension1 {:amino-acid-format :short} protein-extension1ss
      protein-extension2 nil protein-extension2s
      protein-extension2 {:amino-acid-format :short} protein-extension2ss
      protein-extension3 nil protein-extension3s)))

(deftest parse-protein-extension-test
  (testing "returns a correct ProteinExtension"
    (are [s m] (= (mut/parse-protein-extension s) m)
      protein-extension1s protein-extension1
      protein-extension1ss protein-extension1
      protein-extension2s protein-extension2
      protein-extension2ss protein-extension2
      protein-extension3s protein-extension3)))
