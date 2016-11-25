(ns clj-hgvs.core-test
  (:require #?(:clj [clojure.test :refer :all]
               :cljs [cljs.test :refer-macros [deftest is testing]])
            [clj-hgvs.core :as hgvs]
            [clj-hgvs.mutation :as mut]))

(def hgvs1s "NM_005228.3:c.2361G>A")
(def hgvs1m {:transcript "NM_005228.3", :kind :cdna,
             :mutations [(mut/map->CDNAMutation {:numbering "2361",
                                                 :type :substitution,
                                                 :ref "G",
                                                 :alt "A"})]})

(def hgvs2s "c.2361G>A")
(def hgvs2m {:transcript nil, :kind :cdna,
             :mutations [(mut/map->CDNAMutation {:numbering "2361",
                                                 :type :substitution,
                                                 :ref "G",
                                                 :alt "A"})]})

(def hgvs3s "g.[2376A>C;3103del]")
(def hgvs3m {:transcript nil, :kind :genome,
             :mutations [(mut/map->GenomeMutation {:numbering "2376",
                                                   :type :substitution,
                                                   :ref "A",
                                                   :alt "C"}),
                         (mut/map->GenomeMutation {:numbering "3103",
                                                   :type :deletion,
                                                   :ref nil,
                                                   :alt nil})]})

(def hgvs4s "NC_000022.11:g.28703511delA")
(def hgvs4m {:transcript "NC_000022.11", :kind :genome,
             :mutations [(mut/map->GenomeMutation {:numbering "28703511",
                                                   :type :deletion,
                                                   :ref nil,
                                                   :alt "A"})]})

(def hgvs5s "NM_004380.2:c.86-1G>T")
(def hgvs5m {:transcript "NM_004380.2", :kind :cdna,
             :mutations [(mut/map->CDNAMutation {:numbering "86-1",
                                                 :type :substitution,
                                                 :ref "G",
                                                 :alt "T"})]})

(def hgvs6s "NP_005219.2:p.Leu858Arg")
(def hgvs6ss "NP_005219.2:p.L858R")
(def hgvs6m {:transcript "NP_005219.2", :kind :protein,
             :mutations [(mut/map->ProteinSubstitution {:coord {:amino-acid "Leu"
                                                                :position 858}
                                                        :alt "Arg"})]})

(def hgvs7s "NP_001096.1:p.Arg258=")
(def hgvs7m {:transcript "NP_001096.1", :kind :protein,
             :mutations [(mut/map->ProteinSubstitution {:coord {:amino-acid "Arg"
                                                                :position 258}
                                                        :alt "Arg"})]})

(def hgvs8s "NP_001005735.1:p.Leu344Trpfs")
(def hgvs8m {:transcript "NP_001005735.1", :kind :protein,
             :mutations [(mut/map->ProteinFrameShift {:coord {:amino-acid "Leu"
                                                              :position 344}
                                                      :alt "Trp"
                                                      :new-site nil})]})

(deftest hgvs-test
  (testing "allows mutation maps"
    (is (= (hgvs/hgvs "NM_005228.3" :cdna
                      (mut/map->CDNAMutation {:numbering "2361",
                                              :type :substitution,
                                              :ref "G",
                                              :alt "A"}))
           hgvs1m))
    (is (= (hgvs/hgvs nil :genome
                      (mut/map->GenomeMutation {:numbering "2376",
                                                :type :substitution,
                                                :ref "A",
                                                :alt "C"})
                      (mut/map->GenomeMutation {:numbering "3103",
                                                :type :deletion,
                                                :ref nil,
                                                :alt nil}))
           hgvs3m)))
  (testing "allows mutation string"
    (is (= (hgvs/hgvs "NM_005228.3" :cdna "2361G>A") hgvs1m))
    (is (= (hgvs/hgvs nil :genome "[2376A>C;3103del]") hgvs3m))))

(deftest parse-test
  (testing "returns HGVS map"
    (is (= (hgvs/parse hgvs1s) hgvs1m))
    (is (= (hgvs/parse hgvs2s) hgvs2m))
    (is (= (hgvs/parse hgvs3s) hgvs3m))
    (is (= (hgvs/parse hgvs4s) hgvs4m))
    (is (= (hgvs/parse hgvs5s) hgvs5m))
    (is (= (hgvs/parse hgvs6s) hgvs6m))
    (is (= (hgvs/parse hgvs6ss) hgvs6m))
    (is (= (hgvs/parse hgvs7s) hgvs7m))
    (is (= (hgvs/parse hgvs8s) hgvs8m)))
  (testing "throws Exception when an illegal HGVS is passed"
    (is (thrown? #?(:clj Exception, :cljs js/Error) (hgvs/parse ":2361G>A")))
    (is (thrown? #?(:clj Exception, :cljs js/Error) (hgvs/parse "NM_005228.3:2361G>A")))
    (is (thrown? #?(:clj Exception, :cljs js/Error) (hgvs/parse "NM_005228.3:z.2361G>A")))
    (is (thrown? #?(:clj Exception, :cljs js/Error) (hgvs/parse "NM_005228.3:c.G>A")))
    (is (thrown? #?(:clj Exception, :cljs js/Error) (hgvs/parse "NM_005228.3:c.2361G")))
    (is (thrown? #?(:clj Exception, :cljs js/Error) (hgvs/parse "")))
    (is (thrown? #?(:clj Exception, :cljs js/Error) (hgvs/parse nil)))))

(deftest format-test
  (testing "returns HGVS string"
    (is (= (hgvs/format hgvs1m) hgvs1s))
    (is (= (hgvs/format hgvs2m) hgvs2s))
    (is (= (hgvs/format hgvs3m) hgvs3s))
    (is (= (hgvs/format hgvs4m) hgvs4s))
    (is (= (hgvs/format hgvs5m) hgvs5s))
    (is (= (hgvs/format hgvs6m) hgvs6s))
    (is (= (hgvs/format hgvs6m :amino-acid-format :short) hgvs6ss))
    (is (= (hgvs/format hgvs7m) hgvs7s))
    (is (= (hgvs/format hgvs8m) hgvs8s))))
