(ns clj-hgvs.core-test
  (:require #?(:clj [clojure.test :refer :all]
               :cljs [cljs.test :refer-macros [deftest is are testing]])
            [clj-hgvs.coordinate :as coord]
            [clj-hgvs.core :as hgvs]
            [clj-hgvs.mutation :as mut]))

(deftest split-mutations-test
  (testing "splits multiple mutations, returning a vector"
    (are [s e] (= (#'hgvs/split-mutations s) e)
      "123456A>G" ["123456A>G"]
      "[123456A>G;345678G>C]" ["123456A>G" "345678G>C"]
      "[123456A>G];[345678G>C]" ["123456A>G" "345678G>C"]
      "112GAT[14]" ["112GAT[14]"]
      "[112GAT[14];113ATC[15]]" ["112GAT[14]" "113ATC[15]"]
      "[112GAT[14]];[113ATC[15]]" ["112GAT[14]" "113ATC[15]"]))
  (testing "throws Exception when an illegal string is passed"
    (are [s] (thrown? #?(:clj Error, :cljs js/Error) (#'hgvs/split-mutations s))
      "[123456A>G;345678G>C"
      "[123456A>G;[345678G>C]")))

(def hgvs1s "NM_005228.3:c.2361G>A")
(def hgvs1m {:transcript "NM_005228.3", :kind :cdna,
             :mutations [(mut/map->DNASubstitution {:coord (coord/cdna-coordinate 2361)
                                                    :ref "G"
                                                    :type ">"
                                                    :alt "A"})]})
(def hgvs1pm {:transcript "NM_005228.3", :kind "cdna",
              :mutations [{:mutation "dna-substitution"
                           :coord (coord/plain (coord/cdna-coordinate 2361))
                           :ref "G"
                           :type ">"
                           :alt "A"}]})

(def hgvs2s "c.2361G>A")
(def hgvs2m {:transcript nil, :kind :cdna,
             :mutations [(mut/map->DNASubstitution {:coord (coord/cdna-coordinate 2361)
                                                    :ref "G"
                                                    :type ">"
                                                    :alt "A"})]})

(def hgvs3s "g.[2376A>C;3103del]")
(def hgvs3m {:transcript nil, :kind :genome,
             :mutations [(mut/map->DNASubstitution {:coord (coord/genomic-coordinate 2376)
                                                    :ref "A"
                                                    :type ">"
                                                    :alt "C"}),
                         (mut/map->DNADeletion {:coord-start (coord/genomic-coordinate 3103)
                                                :coord-end nil
                                                :ref nil})]})

(def hgvs4s "NC_000022.11:g.28703511delA")
(def hgvs4m {:transcript "NC_000022.11", :kind :genome,
             :mutations [(mut/map->DNADeletion {:coord-start (coord/genomic-coordinate 28703511)
                                                :coord-end nil
                                                :ref "A"})]})

(def hgvs5s "NM_004380.2:c.86-1G>T")
(def hgvs5m {:transcript "NM_004380.2", :kind :cdna,
             :mutations [(mut/map->DNASubstitution {:coord (coord/cdna-coordinate 86 -1 nil)
                                                    :ref "G"
                                                    :type ">"
                                                    :alt "T"})]})

(def hgvs6s "NM_000000.1:r.76a>c")
(def hgvs6m {:transcript "NM_000000.1", :kind :rna,
             :mutations [(mut/map->RNASubstitution {:coord (coord/rna-coordinate 76 nil nil)
                                                    :ref "a"
                                                    :alt "c"})]})

(def hgvs7s "NP_005219.2:p.Leu858Arg")
(def hgvs7ss "NP_005219.2:p.L858R")
(def hgvs7m {:transcript "NP_005219.2", :kind :protein,
             :mutations [(mut/map->ProteinSubstitution {:ref "Leu"
                                                        :coord (coord/protein-coordinate 858)
                                                        :alt "Arg"})]})

(def hgvs8s "NP_001096.1:p.Arg258=")
(def hgvs8m {:transcript "NP_001096.1", :kind :protein,
             :mutations [(mut/map->ProteinSubstitution {:ref "Arg"
                                                        :coord (coord/protein-coordinate 258)
                                                        :alt "Arg"})]})

(def hgvs9s "NP_001005735.1:p.Leu344Trpfs")
(def hgvs9m {:transcript "NP_001005735.1", :kind :protein,
             :mutations [(mut/map->ProteinFrameShift {:ref "Leu"
                                                      :coord (coord/protein-coordinate 344)
                                                      :alt "Trp"
                                                      :new-ter-site nil})]})

(deftest hgvs-test
  (testing "allows mutation maps"
    (is (= (hgvs/hgvs "NM_005228.3" :cdna
                      (mut/map->DNASubstitution {:coord (coord/cdna-coordinate 2361)
                                                 :ref "G"
                                                 :type ">"
                                                 :alt "A"}))
           hgvs1m))
    (is (= (hgvs/hgvs nil :genome
                      (mut/map->DNASubstitution {:coord (coord/genomic-coordinate 2376)
                                                 :ref "A"
                                                 :type ">"
                                                 :alt "C"})
                      (mut/map->DNADeletion {:coord-start (coord/genomic-coordinate 3103)
                                             :coord-end nil
                                             :ref nil}))
           hgvs3m)))
  (testing "allows mutation string"
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
      hgvs7ss hgvs7m
      hgvs8s hgvs8m
      hgvs9s hgvs9m))
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
      hgvs9m hgvs9s)
    (is (= (hgvs/format hgvs4m {:show-bases? true}) hgvs4s))
    (is (= (hgvs/format hgvs7m {:amino-acid-format :short}) hgvs7ss))))

(deftest plain-test
  (testing "returns a plain map representing HGVS"
    (is (= (hgvs/plain hgvs1m) hgvs1pm))))

(deftest restore-test
  (testing "restores a plain map to HGVS"
    (is (= (hgvs/restore hgvs1pm) hgvs1m))))
