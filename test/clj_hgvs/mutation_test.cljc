(ns clj-hgvs.mutation-test
  (:require #?(:clj [clojure.test :refer :all]
               :cljs [cljs.test :refer-macros [deftest is testing]])
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
