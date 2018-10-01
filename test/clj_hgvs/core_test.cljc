(ns clj-hgvs.core-test
  (:require #?(:clj [clojure.test :refer :all]
               :cljs [cljs.test :refer-macros [deftest is are testing]])
            [clj-hgvs.coordinate :as coord]
            [clj-hgvs.core :as hgvs]
            [clj-hgvs.mutation :as mut]))

(deftest transcript?-test
  (testing "returns true if string is transcript"
    (are [s] (true? (#'clj-hgvs.core/transcript? s))
      "NC_000023.10"
      "NC_000023"
      "LRG_199"
      "LRG_199t1"
      "LRG_199p1"
      "NG_012232.1"
      "NM_004006.2"
      "NR_002196.1"
      "NP_003997.1"))
  (testing "returns false if string is not transcript"
    (are [s] (false? (#'clj-hgvs.core/transcript? s))
      "LRG_199.1"
      "NT_000023.10")))

(def hgvs1s "NM_005228.3:c.2361G>A")
(def hgvs1m {:transcript "NM_005228.3", :kind :cdna,
             :mutation (mut/map->DNASubstitution {:coord (coord/cdna-coordinate 2361)
                                                  :ref "G"
                                                  :type ">"
                                                  :alt "A"})})
(def hgvs1pm {:transcript "NM_005228.3", :kind "cdna",
              :mutation {:mutation "dna-substitution"
                         :coord (coord/plain (coord/cdna-coordinate 2361))
                         :ref "G"
                         :type ">"
                         :alt "A"}})

(def hgvs2s "c.2361G>A")
(def hgvs2m {:transcript nil, :kind :cdna,
             :mutation (mut/map->DNASubstitution {:coord (coord/cdna-coordinate 2361)
                                                  :ref "G"
                                                  :type ">"
                                                  :alt "A"})})

(def hgvs3s "g.[2376A>C;3103del]")
(def hgvs3m {:transcript nil, :kind :genome,
             :mutation (mut/dna-alleles
                        [(mut/map->DNASubstitution {:coord (coord/genomic-coordinate 2376)
                                                    :ref "A"
                                                    :type ">"
                                                    :alt "C"}),
                         (mut/map->DNADeletion {:coord-start (coord/genomic-coordinate 3103)
                                                :coord-end nil
                                                :ref nil})]
                        nil)})

(def hgvs4s "NC_000022.11:g.28703511delA")
(def hgvs4m {:transcript "NC_000022.11", :kind :genome,
             :mutation (mut/map->DNADeletion {:coord-start (coord/genomic-coordinate 28703511)
                                              :coord-end nil
                                              :ref "A"})})

(def hgvs5s "NM_004380.2:c.86-1G>T")
(def hgvs5m {:transcript "NM_004380.2", :kind :cdna,
             :mutation (mut/map->DNASubstitution {:coord (coord/cdna-coordinate 86 -1 nil)
                                                  :ref "G"
                                                  :type ">"
                                                  :alt "T"})})

(def hgvs6s "NM_000000.1:r.76a>c")
(def hgvs6m {:transcript "NM_000000.1", :kind :rna,
             :mutation (mut/map->RNASubstitution {:coord (coord/rna-coordinate 76 nil nil)
                                                  :ref "a"
                                                  :alt "c"})})

(def hgvs7s "NM_004006.1:r.0")
(def hgvs7m {:transcript "NM_004006.1", :kind :rna, :mutation (mut/no-rna)})

(def hgvs8s "LRG_199t1:r.?")
(def hgvs8m {:transcript "LRG_199t1", :kind :rna,
             :mutation (mut/rna-unknown-mutation)})

(def hgvs9s "NP_005219.2:p.Leu858Arg")
(def hgvs9ss "NP_005219.2:p.L858R")
(def hgvs9m {:transcript "NP_005219.2", :kind :protein,
             :mutation (mut/map->ProteinSubstitution {:ref "Leu"
                                                      :coord (coord/protein-coordinate 858)
                                                      :alt "Arg"})})

(def hgvs10s "NP_001096.1:p.Arg258=")
(def hgvs10m {:transcript "NP_001096.1", :kind :protein,
              :mutation (mut/map->ProteinSubstitution {:ref "Arg"
                                                       :coord (coord/protein-coordinate 258)
                                                       :alt "Arg"})})

(def hgvs11s "NP_001005735.1:p.Leu344Trpfs")
(def hgvs11m {:transcript "NP_001005735.1", :kind :protein,
              :mutation (mut/map->ProteinFrameShift {:ref "Leu"
                                                     :coord (coord/protein-coordinate 344)
                                                     :alt "Trp"
                                                     :new-ter-site nil})})

(deftest hgvs-test
  (testing "allows mutation records"
    (is (= (hgvs/hgvs "NM_005228.3" :cdna
                      (mut/map->DNASubstitution {:coord (coord/cdna-coordinate 2361)
                                                 :ref "G"
                                                 :type ">"
                                                 :alt "A"}))
           hgvs1m))
    (is (= (hgvs/hgvs nil :genome
                      (mut/dna-alleles
                       [(mut/map->DNASubstitution {:coord (coord/genomic-coordinate 2376)
                                                   :ref "A"
                                                   :type ">"
                                                   :alt "C"})
                        (mut/map->DNADeletion {:coord-start (coord/genomic-coordinate 3103)
                                               :coord-end nil
                                               :ref nil})]
                       nil))
           hgvs3m)))
  (testing "allows mutation strings"
    (is (= (hgvs/hgvs "NM_005228.3" :cdna "2361G>A") hgvs1m))
    (is (= (hgvs/hgvs nil :genome "[2376A>C;3103del]") hgvs3m))))

(deftest parse-test
  (testing "returns HGVS map"
    (are [s m] (= (hgvs/parse s) m)
      hgvs1s hgvs1m
      hgvs2s hgvs2m
      hgvs3s hgvs3m
      hgvs4s hgvs4m
      hgvs5s hgvs5m
      hgvs6s hgvs6m
      hgvs7s hgvs7m
      hgvs8s hgvs8m
      hgvs9s hgvs9m
      hgvs9ss hgvs9m
      hgvs10s hgvs10m
      hgvs11s hgvs11m))
  (testing "throws Exception when an illegal HGVS is passed"
    (are [x] (thrown? #?(:clj Exception, :cljs js/Error) (hgvs/parse x))
      ":2361G>A"
      "NM_005228.3:2361G>A"
      "NM_005228.3:z.2361G>A"
      "NM_005228.3:c.G>A"
      "NM_005228.3:c.2361G"
      ""
      nil)))

(deftest format-test
  (testing "returns HGVS string"
    (are [m s] (= (hgvs/format m) s)
      hgvs1m hgvs1s
      hgvs2m hgvs2s
      hgvs3m hgvs3s
      hgvs5m hgvs5s
      hgvs6m hgvs6s
      hgvs7m hgvs7s
      hgvs8m hgvs8s
      hgvs9m hgvs9s
      hgvs10m hgvs10s
      hgvs11m hgvs11s)
    (is (= (hgvs/format hgvs4m {:show-bases? true}) hgvs4s))
    (is (= (hgvs/format hgvs9m {:amino-acid-format :short}) hgvs9ss))))

(deftest plain-test
  (testing "returns a plain map representing HGVS"
    (is (= (hgvs/plain hgvs1m) hgvs1pm))))

(deftest restore-test
  (testing "restores a plain map to HGVS"
    (is (= (hgvs/restore hgvs1pm) hgvs1m))))

(deftest normalize-test
  (testing "returns a normalized HGVS string"
    (are [s ns] (= (hgvs/normalize s) ns)
      "NG_012232.1:g.19_21delTGC" "NG_012232.1:g.19_21del"
      "g.6775delTinsGA" "g.6775delinsGA"
      "p.I327Rfs*?" "p.Ile327ArgfsTer?"
      "p.*110Glnext*17" "p.Ter110Glnext*17"))
  (testing "throws exception"
    (are [s] (thrown? #?(:clj Exception, :cljs js/Error) (hgvs/normalize s))
      ":2361G>A"
      "NM_005228.3:2361G>A"
      ""
      nil)))
