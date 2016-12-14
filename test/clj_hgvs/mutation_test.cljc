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
(def dna-substitution1 (mut/map->DNASubstitution {:coord-start (coord/->GenomicCoordinate 45576)
                                                  :coord-end nil
                                                  :ref "A"
                                                  :type ">"
                                                  :alt "C"}))

(def dna-substitution2s "88+1G>T")
(def dna-substitution2k :cdna)
(def dna-substitution2 (mut/map->DNASubstitution {:coord-start (coord/->CDNACoordinate 88 nil 1)
                                                  :coord-end nil
                                                  :ref "G"
                                                  :type ">"
                                                  :alt "T"}))

(def dna-substitution3s "123G=")
(def dna-substitution3k :cdna)
(def dna-substitution3 (mut/map->DNASubstitution {:coord-start (coord/->CDNACoordinate 123 nil nil)
                                                  :coord-end nil
                                                  :ref "G"
                                                  :type "="
                                                  :alt nil}))

(def dna-substitution4s "85C=/>T")
(def dna-substitution4k :cdna)
(def dna-substitution4 (mut/map->DNASubstitution {:coord-start (coord/->CDNACoordinate 85 nil nil)
                                                  :coord-end nil
                                                  :ref "C"
                                                  :type "=/>"
                                                  :alt "T"}))

(def dna-substitution5s "85C=//>T")
(def dna-substitution5k :cdna)
(def dna-substitution5 (mut/map->DNASubstitution {:coord-start (coord/->CDNACoordinate 85 nil nil)
                                                  :coord-end nil
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
(def dna-deletion1 (mut/map->DNADeletion {:coord-start (coord/->GenomicCoordinate 7)
                                          :coord-end nil
                                          :ref nil}))

(def dna-deletion2s "6_8del")
(def dna-deletion2k :genome)
(def dna-deletion2 (mut/map->DNADeletion {:coord-start (coord/->GenomicCoordinate 6)
                                          :coord-end (coord/->GenomicCoordinate 8)
                                          :ref nil}))

(def dna-deletion3s "120_123+48del")
(def dna-deletion3k :cdna)
(def dna-deletion3 (mut/map->DNADeletion {:coord-start (coord/->CDNACoordinate 120 nil nil)
                                          :coord-end (coord/->CDNACoordinate 123 nil 48)
                                          :ref nil}))

(def dna-deletion4s "(4071+1_4072-1)_(5145+1_5146-1)del")
(def dna-deletion4k :cdna)
(def dna-deletion4 (mut/map->DNADeletion {:coord-start [(coord/->CDNACoordinate 4071 nil 1)
                                                        (coord/->CDNACoordinate 4072 nil -1)]
                                          :coord-end [(coord/->CDNACoordinate 5145 nil 1)
                                                      (coord/->CDNACoordinate 5146 nil -1)]
                                          :ref nil}))

(def dna-deletion5s "(?_-30)_(12+1_13-1)del")
(def dna-deletion5k :cdna)
(def dna-deletion5 (mut/map->DNADeletion {:coord-start [(coord/->UnknownCoordinate)
                                                        (coord/->CDNACoordinate 30 :up nil)]
                                          :coord-end [(coord/->CDNACoordinate 12 nil 1)
                                                      (coord/->CDNACoordinate 13 nil -1)]
                                          :ref nil}))

(def dna-deletion6s "(?_-1)_(*1_?)del")
(def dna-deletion6k :cdna)
(def dna-deletion6 (mut/map->DNADeletion {:coord-start [(coord/->UnknownCoordinate)
                                                        (coord/->CDNACoordinate 1 :up nil)]
                                          :coord-end [(coord/->CDNACoordinate 1 :down nil)
                                                      (coord/->UnknownCoordinate)]
                                          :ref nil}))

(deftest format-dna-deletion-test
  (testing "returns a string expression of a DNA deletion"
    (are [m s] (= (mut/format m nil) s)
      dna-deletion1 dna-deletion1s
      dna-deletion2 dna-deletion2s
      dna-deletion3 dna-deletion3s
      dna-deletion4 dna-deletion4s
      dna-deletion5 dna-deletion5s
      dna-deletion6 dna-deletion6s)))

(deftest parse-dna-deletion-test
  (testing "returns a correct DNADeletion"
    (are [s k m] (= (mut/parse-dna-deletion s k) m)
      dna-deletion1s dna-deletion1k dna-deletion1
      dna-deletion2s dna-deletion2k dna-deletion2
      dna-deletion3s dna-deletion3k dna-deletion3
      dna-deletion4s dna-deletion4k dna-deletion4
      dna-deletion5s dna-deletion5k dna-deletion5
      dna-deletion6s dna-deletion6k dna-deletion6)))

;;; DNA - duplication

(def dna-duplication1s "7dup")
(def dna-duplication1k :genome)
(def dna-duplication1 (mut/map->DNADuplication {:coord-start (coord/->GenomicCoordinate 7)
                                                :coord-end nil
                                                :ref nil}))

(def dna-duplication2s "6_8dup")
(def dna-duplication2k :genome)
(def dna-duplication2 (mut/map->DNADuplication {:coord-start (coord/->GenomicCoordinate 6)
                                                :coord-end (coord/->GenomicCoordinate 8)
                                                :ref nil}))

(def dna-duplication3s "120_123+48dup")
(def dna-duplication3k :cdna)
(def dna-duplication3 (mut/map->DNADuplication {:coord-start (coord/->CDNACoordinate 120 nil nil)
                                                :coord-end (coord/->CDNACoordinate 123 nil 48)
                                                :ref nil}))

(def dna-duplication4s "(4071+1_4072-1)_(5145+1_5146-1)dup")
(def dna-duplication4k :cdna)
(def dna-duplication4 (mut/map->DNADuplication {:coord-start [(coord/->CDNACoordinate 4071 nil 1)
                                                              (coord/->CDNACoordinate 4072 nil -1)]
                                                :coord-end [(coord/->CDNACoordinate 5145 nil 1)
                                                            (coord/->CDNACoordinate 5146 nil -1)]
                                                :ref nil}))

(def dna-duplication5s "(?_-30)_(12+1_13-1)dup")
(def dna-duplication5k :cdna)
(def dna-duplication5 (mut/map->DNADuplication {:coord-start [(coord/->UnknownCoordinate)
                                                              (coord/->CDNACoordinate 30 :up nil)]
                                                :coord-end [(coord/->CDNACoordinate 12 nil 1)
                                                            (coord/->CDNACoordinate 13 nil -1)]
                                                :ref nil}))

(def dna-duplication6s "(?_-1)_(*1_?)dup")
(def dna-duplication6k :cdna)
(def dna-duplication6 (mut/map->DNADuplication {:coord-start [(coord/->UnknownCoordinate)
                                                              (coord/->CDNACoordinate 1 :up nil)]
                                                :coord-end [(coord/->CDNACoordinate 1 :down nil)
                                                            (coord/->UnknownCoordinate)]
                                                :ref nil}))

(deftest format-dna-duplication-test
  (testing "returns a string expression of a DNA duplication"
    (are [m s] (= (mut/format m nil) s)
      dna-duplication1 dna-duplication1s
      dna-duplication2 dna-duplication2s
      dna-duplication3 dna-duplication3s
      dna-duplication4 dna-duplication4s
      dna-duplication5 dna-duplication5s
      dna-duplication6 dna-duplication6s)))

(deftest parse-dna-duplication-test
  (testing "returns a correct DNADuplication"
    (are [s k m] (= (mut/parse-dna-duplication s k) m)
      dna-duplication1s dna-duplication1k dna-duplication1
      dna-duplication2s dna-duplication2k dna-duplication2
      dna-duplication3s dna-duplication3k dna-duplication3
      dna-duplication4s dna-duplication4k dna-duplication4
      dna-duplication5s dna-duplication5k dna-duplication5
      dna-duplication6s dna-duplication6k dna-duplication6)))

;;; DNA - insertion

(def dna-insertion1s "5756_5757insAGG")
(def dna-insertion1k :genome)
(def dna-insertion1 (mut/map->DNAInsertion {:coord-start (coord/->GenomicCoordinate 5756)
                                            :coord-end (coord/->GenomicCoordinate 5757)
                                            :alt "AGG"}))

(def dna-insertion2s "123_124insL37425.1:23_361")
(def dna-insertion2k :genome)
(def dna-insertion2 (mut/map->DNAInsertion {:coord-start (coord/->GenomicCoordinate 123)
                                            :coord-end (coord/->GenomicCoordinate 124)
                                            :alt {:transcript "L37425.1"
                                                  :coord-start (coord/->GenomicCoordinate 23)
                                                  :coord-end (coord/->GenomicCoordinate 361)}}))

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
(def dna-inversion1 (mut/map->DNAInversion {:coord-start (coord/->GenomicCoordinate 1077)
                                            :coord-end (coord/->GenomicCoordinate 1080)}))

(def dna-inversion2s "77_80inv")
(def dna-inversion2k :cdna)
(def dna-inversion2 (mut/map->DNAInversion {:coord-start (coord/->CDNACoordinate 77 nil nil)
                                            :coord-end (coord/->CDNACoordinate 80 nil nil)}))

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
(def dna-conversion1 (mut/map->DNAConversion {:coord-start (coord/->GenomicCoordinate 333)
                                              :coord-end (coord/->GenomicCoordinate 590)
                                              :alt {:transcript nil
                                                    :kind nil
                                                    :coord-start (coord/->GenomicCoordinate 1844)
                                                    :coord-end (coord/->GenomicCoordinate 2101)}}))

(def dna-conversion2s "415_1655conAC096506.5:g.409_1683")
(def dna-conversion2k :genome)
(def dna-conversion2 (mut/map->DNAConversion {:coord-start (coord/->GenomicCoordinate 415)
                                              :coord-end (coord/->GenomicCoordinate 1655)
                                              :alt {:transcript "AC096506.5"
                                                    :kind :genome
                                                    :coord-start (coord/->GenomicCoordinate 409)
                                                    :coord-end (coord/->GenomicCoordinate 1683)}}))

(def dna-conversion3s "15_355conNM_004006.1:20_360")
(def dna-conversion3k :cdna)
(def dna-conversion3 (mut/map->DNAConversion {:coord-start (coord/->CDNACoordinate 15 nil nil)
                                              :coord-end (coord/->CDNACoordinate 355 nil nil)
                                              :alt {:transcript "NM_004006.1"
                                                    :kind nil
                                                    :coord-start (coord/->CDNACoordinate 20 nil nil)
                                                    :coord-end (coord/->CDNACoordinate 360 nil nil)}}))

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
(def dna-indel1 (mut/map->DNAIndel {:coord-start (coord/->GenomicCoordinate 6775)
                                    :coord-end nil
                                    :alt "GA"}))

(def dna-indel2s "145_147delinsTGG")
(def dna-indel2k :cdna)
(def dna-indel2 (mut/map->DNAIndel {:coord-start (coord/->CDNACoordinate 145 nil nil)
                                    :coord-end (coord/->CDNACoordinate 147 nil nil)
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
(def dna-repeated-seqs (mut/map->DNARepeatedSeqs {:coord-start (coord/->GenomicCoordinate 123)
                                                  :coord-end (coord/->GenomicCoordinate 124)
                                                  :ncopy 14}))

(deftest format-dna-repeated-seqs-test
  (testing "returns a string expression of a DNA repeated sequences"
    (are [m s] (= (mut/format m nil) s)
      dna-repeated-seqs dna-repeated-seqss)))

(deftest parse-dna-repeated-seqs-test
  (testing "returns a correct DNARepeatedSeqs"
    (are [s k m] (= (mut/parse-dna-repeated-seqs s k) m)
      dna-repeated-seqss dna-repeated-seqsk dna-repeated-seqs)))

;;; Protein mutations

;;; Protein - substitution

(def protein-substitution1s "Arg54Ser")
(def protein-substitution1ss "R54S")
(def protein-substitution1 (mut/map->ProteinSubstitution {:ref "Arg"
                                                          :coord (coord/->ProteinCoordinate 54)
                                                          :alt "Ser"}))

(def protein-substitution2s "Cys123=")
(def protein-substitution2ss "C123=")
(def protein-substitution2 (mut/map->ProteinSubstitution {:ref "Cys"
                                                          :coord (coord/->ProteinCoordinate 123)
                                                          :alt "Cys"}))

(deftest format-protein-substitution-test
  (testing "returns a string expression of a protein substitution"
    (is (= (mut/format protein-substitution1 nil) protein-substitution1s))
    (is (= (mut/format protein-substitution1 {:amino-acid-format :short}) protein-substitution1ss))
    (is (= (mut/format protein-substitution2 nil) protein-substitution2s))
    (is (= (mut/format protein-substitution2 {:amino-acid-format :short}) protein-substitution2ss))))

(deftest parse-protein-substitution-test
  (testing "returns a correct ProteinSubstitution"
    (is (= (mut/parse-protein-substitution protein-substitution1s) protein-substitution1))
    (is (= (mut/parse-protein-substitution protein-substitution1ss) protein-substitution1))
    (is (= (mut/parse-protein-substitution protein-substitution2s) protein-substitution2))
    (is (= (mut/parse-protein-substitution protein-substitution2ss) protein-substitution2))))

;;; Protein - deletion

(def protein-deletion1s "Ala3del")
(def protein-deletion1ss "A3del")
(def protein-deletion1 (mut/map->ProteinDeletion {:ref-start "Ala"
                                                  :coord-start (coord/->ProteinCoordinate 3)
                                                  :ref-end nil
                                                  :coord-end nil}))

(def protein-deletion2s "Cys76_Glu79del")
(def protein-deletion2ss "C76_E79del")
(def protein-deletion2 (mut/map->ProteinDeletion {:ref-start "Cys"
                                                  :coord-start (coord/->ProteinCoordinate 76)
                                                  :ref-end "Glu"
                                                  :coord-end (coord/->ProteinCoordinate 79)}))
(deftest format-protein-deletion-test
  (testing "returns a string expression of a protein deletion"
    (is (= (mut/format protein-deletion1 nil) protein-deletion1s))
    (is (= (mut/format protein-deletion1 {:amino-acid-format :short}) protein-deletion1ss))
    (is (= (mut/format protein-deletion2 nil) protein-deletion2s))
    (is (= (mut/format protein-deletion2 {:amino-acid-format :short}) protein-deletion2ss))))

(deftest parse-protein-deletion-test
  (testing "returns a correct ProteinDeletion"
    (is (= (mut/parse-protein-deletion protein-deletion1s) protein-deletion1))
    (is (= (mut/parse-protein-deletion protein-deletion1ss) protein-deletion1))
    (is (= (mut/parse-protein-deletion protein-deletion2s) protein-deletion2))
    (is (= (mut/parse-protein-deletion protein-deletion2ss) protein-deletion2))))

;;; Protein - duplication

(def protein-duplication1s "Ala3dup")
(def protein-duplication1ss "A3dup")
(def protein-duplication1 (mut/map->ProteinDuplication {:ref-start "Ala"
                                                        :coord-start (coord/->ProteinCoordinate 3)
                                                        :ref-end nil
                                                        :coord-end nil}))

(def protein-duplication2s "Ala3_Ser5dup")
(def protein-duplication2ss "A3_S5dup")
(def protein-duplication2 (mut/map->ProteinDuplication {:ref-start "Ala"
                                                        :coord-start (coord/->ProteinCoordinate 3)
                                                        :ref-end "Ser"
                                                        :coord-end (coord/->ProteinCoordinate 5)}))

(deftest format-protein-duplication-test
  (testing "returns a string expression of a protein duplication"
    (is (= (mut/format protein-duplication1 nil) protein-duplication1s))
    (is (= (mut/format protein-duplication1 {:amino-acid-format :short}) protein-duplication1ss))
    (is (= (mut/format protein-duplication2 nil) protein-duplication2s))
    (is (= (mut/format protein-duplication2 {:amino-acid-format :short}) protein-duplication2ss))))

(deftest parse-protein-duplication-test
  (testing "returns a correct ProteinDuplication"
    (is (= (mut/parse-protein-duplication protein-duplication1s) protein-duplication1))
    (is (= (mut/parse-protein-duplication protein-duplication1ss) protein-duplication1))
    (is (= (mut/parse-protein-duplication protein-duplication2s) protein-duplication2))
    (is (= (mut/parse-protein-duplication protein-duplication2ss) protein-duplication2))))

;;; Protein - insertion

(def protein-insertions "Lys23_Leu24insArgSerGln")
(def protein-insertionss "K23_L24insRSQ")
(def protein-insertion (mut/map->ProteinInsertion {:ref-start "Lys"
                                                   :coord-start (coord/->ProteinCoordinate 23)
                                                   :ref-end "Leu"
                                                   :coord-end (coord/->ProteinCoordinate 24)
                                                   :alts ["Arg" "Ser" "Gln"]}))

(deftest format-protein-insertion-test
  (testing "returns a string expression of a protein insertion"
    (is (= (mut/format protein-insertion nil) protein-insertions))
    (is (= (mut/format protein-insertion {:amino-acid-format :short}) protein-insertionss))))

(deftest parse-protein-insertion-test
  (testing "returns a correct ProteinInsertion"
    (is (= (mut/parse-protein-insertion protein-insertions) protein-insertion))
    (is (= (mut/parse-protein-insertion protein-insertionss) protein-insertion))))

;;; Protein - indel

(def protein-indel1s "Cys28delinsTrpVal")
(def protein-indel1ss "C28delinsWV")
(def protein-indel1 (mut/map->ProteinIndel {:ref-start "Cys"
                                            :coord-start (coord/->ProteinCoordinate 28)
                                            :ref-end nil
                                            :coord-end nil
                                            :alts ["Trp" "Val"]}))

(def protein-indel2s "Cys28_Lys29delinsTrp")
(def protein-indel2ss "C28_K29delinsW")
(def protein-indel2 (mut/map->ProteinIndel {:ref-start "Cys"
                                            :coord-start (coord/->ProteinCoordinate 28)
                                            :ref-end "Lys"
                                            :coord-end (coord/->ProteinCoordinate 29)
                                            :alts ["Trp"]}))

(deftest format-protein-indel-test
  (testing "returns a string expression of a protein indel"
    (is (= (mut/format protein-indel1 nil) protein-indel1s))
    (is (= (mut/format protein-indel1 {:amino-acid-format :short}) protein-indel1ss))
    (is (= (mut/format protein-indel2 nil) protein-indel2s))
    (is (= (mut/format protein-indel2 {:amino-acid-format :short}) protein-indel2ss))))

(deftest parse-protein-indel-test
  (testing "returns a correct ProteinIndel"
    (is (= (mut/parse-protein-indel protein-indel1s) protein-indel1))
    (is (= (mut/parse-protein-indel protein-indel1ss) protein-indel1))
    (is (= (mut/parse-protein-indel protein-indel2s) protein-indel2))
    (is (= (mut/parse-protein-indel protein-indel2ss) protein-indel2))))

;;; Protein - repeated sequences

(def protein-repeated-seqs1s "Ala2[10]")
(def protein-repeated-seqs1ss "A2[10]")
(def protein-repeated-seqs1 (mut/map->ProteinRepeatedSeqs {:ref-start "Ala"
                                                           :coord-start (coord/->ProteinCoordinate 2)
                                                           :ref-end nil
                                                           :coord-end nil
                                                           :ncopy 10
                                                           :ncopy-other nil}))

(def protein-repeated-seqs2s "Ala2[10];[11]")
(def protein-repeated-seqs2ss "A2[10];[11]")
(def protein-repeated-seqs2 (mut/map->ProteinRepeatedSeqs {:ref-start "Ala"
                                                           :coord-start (coord/->ProteinCoordinate 2)
                                                           :ref-end nil
                                                           :coord-end nil
                                                           :ncopy 10
                                                           :ncopy-other 11}))

(def protein-repeated-seqs3s "Arg65_Ser67[12]")
(def protein-repeated-seqs3ss "R65_S67[12]")
(def protein-repeated-seqs3 (mut/map->ProteinRepeatedSeqs {:ref-start "Arg"
                                                           :coord-start (coord/->ProteinCoordinate 65)
                                                           :ref-end "Ser"
                                                           :coord-end (coord/->ProteinCoordinate 67)
                                                           :ncopy 12
                                                           :ncopy-other nil}))

(deftest format-protein-repeated-seqs-test
  (testing "returns a string expression of a protein repeated-seqs"
    (is (= (mut/format protein-repeated-seqs1 nil) protein-repeated-seqs1s))
    (is (= (mut/format protein-repeated-seqs1 {:amino-acid-format :short}) protein-repeated-seqs1ss))
    (is (= (mut/format protein-repeated-seqs2 nil) protein-repeated-seqs2s))
    (is (= (mut/format protein-repeated-seqs2 {:amino-acid-format :short}) protein-repeated-seqs2ss))
    (is (= (mut/format protein-repeated-seqs3 nil) protein-repeated-seqs3s))
    (is (= (mut/format protein-repeated-seqs3 {:amino-acid-format :short}) protein-repeated-seqs3ss))))

(deftest parse-protein-repeated-seqs-test
  (testing "returns a correct ProteinRepeatedSeqs"
    (is (= (mut/parse-protein-repeated-seqs protein-repeated-seqs1s) protein-repeated-seqs1))
    (is (= (mut/parse-protein-repeated-seqs protein-repeated-seqs1ss) protein-repeated-seqs1))
    (is (= (mut/parse-protein-repeated-seqs protein-repeated-seqs2s) protein-repeated-seqs2))
    (is (= (mut/parse-protein-repeated-seqs protein-repeated-seqs2ss) protein-repeated-seqs2))
    (is (= (mut/parse-protein-repeated-seqs protein-repeated-seqs3s) protein-repeated-seqs3))
    (is (= (mut/parse-protein-repeated-seqs protein-repeated-seqs3ss) protein-repeated-seqs3))))

;;; Protein - frame shift

(def protein-frame-shift1s "Arg97ProfsTer23")
(def protein-frame-shift1 (mut/map->ProteinFrameShift {:ref "Arg"
                                                       :coord (coord/->ProteinCoordinate 97)
                                                       :alt "Pro"
                                                       :new-site "Ter23"}))

(def protein-frame-shift2s "Arg97fs")
(def protein-frame-shift2ss "R97fs")
(def protein-frame-shift2 (mut/map->ProteinFrameShift {:ref "Arg"
                                                       :coord (coord/->ProteinCoordinate 97)
                                                       :alt nil
                                                       :new-site nil}))

(def protein-frame-shift3s "Ile327Argfs*?")
(def protein-frame-shift3 (mut/map->ProteinFrameShift {:ref "Ile"
                                                       :coord (coord/->ProteinCoordinate 327)
                                                       :alt "Arg"
                                                       :new-site "*?"}))

(def protein-frame-shift4s "Gln151Thrfs*9")
(def protein-frame-shift4 (mut/map->ProteinFrameShift {:ref "Gln"
                                                       :coord (coord/->ProteinCoordinate 151)
                                                       :alt "Thr"
                                                       :new-site "*9"}))
(deftest format-protein-frame-shift-test
  (testing "returns a string expression of a protein frame-shift"
    (is (= (mut/format protein-frame-shift1 nil) protein-frame-shift1s))
    (is (= (mut/format protein-frame-shift2 nil) protein-frame-shift2s))
    (is (= (mut/format protein-frame-shift2 {:amino-acid-format :short}) protein-frame-shift2ss))
    (is (= (mut/format protein-frame-shift3 nil) protein-frame-shift3s))
    (is (= (mut/format protein-frame-shift4 nil) protein-frame-shift4s))))

(deftest parse-protein-frame-shift-test
  (testing "returns a correct ProteinFrameShift"
    (is (= (mut/parse-protein-frame-shift protein-frame-shift1s) protein-frame-shift1))
    (is (= (mut/parse-protein-frame-shift protein-frame-shift2s) protein-frame-shift2))
    (is (= (mut/parse-protein-frame-shift protein-frame-shift2ss) protein-frame-shift2))
    (is (= (mut/parse-protein-frame-shift protein-frame-shift3s) protein-frame-shift3))
    (is (= (mut/parse-protein-frame-shift protein-frame-shift4s) protein-frame-shift4))))

;;; Protein - extension

(def protein-extension1s "Met1ext-5")
(def protein-extension1ss "M1ext-5")
(def protein-extension1 (mut/map->ProteinExtension {:ref "Met"
                                                    :coord (coord/->ProteinCoordinate 1)
                                                    :alt nil
                                                    :new-site "-5"}))

(def protein-extension2s "Met1Valext-12")
(def protein-extension2ss "M1Vext-12")
(def protein-extension2 (mut/map->ProteinExtension {:ref "Met"
                                                    :coord (coord/->ProteinCoordinate 1)
                                                    :alt "Val"
                                                    :new-site "-12"}))

(def protein-extension3s "Ter110GlnextTer17")
(def protein-extension3 (mut/map->ProteinExtension {:ref "Ter"
                                                    :coord (coord/->ProteinCoordinate 110)
                                                    :alt "Gln"
                                                    :new-site "Ter17"}))

(deftest format-protein-extension-test
  (testing "returns a string expression of a protein extension"
    (is (= (mut/format protein-extension1 nil) protein-extension1s))
    (is (= (mut/format protein-extension1 {:amino-acid-format :short}) protein-extension1ss))
    (is (= (mut/format protein-extension2 nil) protein-extension2s))
    (is (= (mut/format protein-extension2 {:amino-acid-format :short}) protein-extension2ss))
    (is (= (mut/format protein-extension3 nil) protein-extension3s))))

(deftest parse-protein-extension-test
  (testing "returns a correct ProteinExtension"
    (is (= (mut/parse-protein-extension protein-extension1s) protein-extension1))
    (is (= (mut/parse-protein-extension protein-extension1ss) protein-extension1))
    (is (= (mut/parse-protein-extension protein-extension2s) protein-extension2))
    (is (= (mut/parse-protein-extension protein-extension2ss) protein-extension2))
    (is (= (mut/parse-protein-extension protein-extension3s) protein-extension3))))
