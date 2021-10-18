(ns clj-hgvs.core-test
  (:require #?(:clj [clojure.pprint :as pp])
            [clojure.spec.alpha :as s]
            [clojure.string :as string]
            [clojure.test :refer [are deftest is testing]]
            [clj-hgvs.coordinate :as coord]
            [clj-hgvs.core :as hgvs]
            [clj-hgvs.mutation :as mut]
            [clj-hgvs.repairer :as repairer]
            clj-hgvs.test-common))

(deftest transcript-spec-test
  (are [s] (s/valid? ::hgvs/transcript s)
    "NC_000023.10"
    "NC_000023"
    "LRG_199"
    "LRG_199t1"
    "LRG_199p1"
    "NG_012232.1"
    "NM_004006.2"
    "NR_002196.1"
    "NP_003997.1"
    "J01749.1"
    "ENST00000000001.1"
    "ENSP00000000001.11"
    "MGP_CBAJ_G00000000002.1")
  (are [s] (not (s/valid? ::hgvs/transcript s))
    "LRG_199.1"
    "NT_000023.10"
    "ENST0001"
    "ENSF00000000001.1"
    "MGP_CBAJ_P0000000002.1"))

(def hgvs1s "NM_005228.3:c.2361G>A")
(def hgvs1m (hgvs/map->HGVS
             {:transcript "NM_005228.3", :kind :coding-dna,
              :mutation (mut/map->DNASubstitution {:coord (coord/coding-dna-coordinate 2361)
                                                   :ref "G"
                                                   :type ">"
                                                   :alt "A"})}))
(def hgvs1pm {:transcript "NM_005228.3", :kind "coding-dna",
              :mutation {:mutation "dna-substitution"
                         :coord (coord/plain (coord/coding-dna-coordinate 2361))
                         :ref "G"
                         :type ">"
                         :alt "A"}})

(def hgvs2s "c.2361G>A")
(def hgvs2m (hgvs/map->HGVS
             {:transcript nil, :kind :coding-dna,
              :mutation (mut/map->DNASubstitution {:coord (coord/coding-dna-coordinate 2361)
                                                   :ref "G"
                                                   :type ">"
                                                   :alt "A"})}))

(def hgvs3s "g.[2376A>C;3103del]")
(def hgvs3m (hgvs/map->HGVS
             {:transcript nil, :kind :genome,
              :mutation (mut/dna-alleles
                         [(mut/map->DNASubstitution {:coord (coord/genomic-coordinate 2376)
                                                     :ref "A"
                                                     :type ">"
                                                     :alt "C"}),
                          (mut/map->DNADeletion {:coord-start (coord/genomic-coordinate 3103)
                                                 :coord-end nil
                                                 :ref nil})]
                         nil)}))

(def hgvs4s "NC_000022.11:g.28703511delA")
(def hgvs4m (hgvs/map->HGVS
             {:transcript "NC_000022.11", :kind :genome,
              :mutation (mut/map->DNADeletion {:coord-start (coord/genomic-coordinate 28703511)
                                               :coord-end nil
                                               :ref "A"})}))

(def hgvs5s "NM_004380.2:c.86-1G>T")
(def hgvs5m (hgvs/map->HGVS
             {:transcript "NM_004380.2", :kind :coding-dna,
              :mutation (mut/map->DNASubstitution {:coord (coord/coding-dna-coordinate 86 -1 nil)
                                                   :ref "G"
                                                   :type ">"
                                                   :alt "T"})}))

(def hgvs6s "NM_000000.1:r.76a>c")
(def hgvs6m (hgvs/map->HGVS
             {:transcript "NM_000000.1", :kind :rna,
              :mutation (mut/map->RNASubstitution {:coord (coord/rna-coordinate 76 nil nil)
                                                   :ref "a"
                                                   :alt "c"})}))

(def hgvs7s "NM_004006.1:r.0")
(def hgvs7m (hgvs/map->HGVS
             {:transcript "NM_004006.1", :kind :rna, :mutation (mut/no-rna)}))

(def hgvs8s "LRG_199t1:r.?")
(def hgvs8m (hgvs/map->HGVS
             {:transcript "LRG_199t1", :kind :rna,
              :mutation (mut/rna-unknown-mutation)}))

(def hgvs9s "NP_005219.2:p.Leu858Arg")
(def hgvs9ss "NP_005219.2:p.L858R")
(def hgvs9m (hgvs/map->HGVS
             {:transcript "NP_005219.2", :kind :protein,
              :mutation (mut/map->ProteinSubstitution {:ref "Leu"
                                                       :coord (coord/protein-coordinate 858)
                                                       :alt "Arg"})}))

(def hgvs10s "NP_001096.1:p.Arg258=")
(def hgvs10m (hgvs/map->HGVS
              {:transcript "NP_001096.1", :kind :protein,
               :mutation (mut/map->ProteinSubstitution {:ref "Arg"
                                                        :coord (coord/protein-coordinate 258)
                                                        :alt "Arg"})}))

(def hgvs11s "NP_001005735.1:p.Leu344Trpfs")
(def hgvs11m (hgvs/map->HGVS
              {:transcript "NP_001005735.1", :kind :protein,
               :mutation (mut/map->ProteinFrameShift {:ref "Leu"
                                                      :coord (coord/protein-coordinate 344)
                                                      :alt "Trp"
                                                      :new-ter-site nil})}))

(def hgvs12s "NC_012920.1:m.16563_13del")
(def hgvs12m (hgvs/map->HGVS
              {:transcript "NC_012920.1"
               :kind :mitochondria
               :mutation (mut/dna-deletion (coord/mitochondrial-coordinate 16563)
                                           (coord/mitochondrial-coordinate 13))}))

(def hgvs13s "J01749.1:o.4344_197dup")
(def hgvs13m (hgvs/map->HGVS
              {:transcript "J01749.1"
               :kind :circular-dna
               :mutation (mut/dna-duplication (coord/circular-dna-coordinate 4344)
                                              (coord/circular-dna-coordinate 197))}))

(def hgvs14s "ENST00000331920.11:c.*1+68G>A")
(def hgvs14m (hgvs/map->HGVS
              {:transcript "ENST00000331920.11"
               :kind :coding-dna
               :mutation (mut/map->DNASubstitution {:coord (coord/coding-dna-coordinate 1 68 :downstream)
                                                    :ref "G"
                                                    :type ">"
                                                    :alt "A"})}))
(deftest hgvs-test
  (testing "allows mutation records"
    (is (= (hgvs/hgvs "NM_005228.3" :coding-dna
                      (mut/map->DNASubstitution {:coord (coord/coding-dna-coordinate 2361)
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
    (is (= (hgvs/hgvs "NM_005228.3" :coding-dna "2361G>A") hgvs1m))
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
      hgvs11s hgvs11m
      hgvs12s hgvs12m
      hgvs13s hgvs13m
      hgvs14s hgvs14m))
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
      hgvs11m hgvs11s
      hgvs12m hgvs12s
      hgvs13m hgvs13s)
    (is (= (hgvs/format hgvs4m {:show-bases? true}) hgvs4s))
    (is (= (hgvs/format hgvs9m {:amino-acid-format :short}) hgvs9ss))))

(deftest equiv-test
  (testing "one arg"
    (is (true? (hgvs/== (hgvs/parse "NM_005228:c.2361G>A"))))
    (is (true? (hgvs/== nil))))

  (testing "two args"
    (are [s1 s2] (true? (hgvs/== (hgvs/parse s1) (hgvs/parse s2)))
      "NM_005228:c.2361G>A"   "NM_005228:c.2361G>A"
      "NM_005228:c.2361G>A"   "NM_005228.3:c.2361G>A"
      "NM_005228.3:c.2361G>A" "NM_005228.4:c.2361G>A"
      "c.2361G>A"             "c.2361G>A"
      "p.L858R"               "p.Leu858Arg"
      "p.K53Afs*9"            "p.K53Afs")
    (is (true? (hgvs/== nil nil)))
    (are [s1 s2] (false? (hgvs/== (hgvs/parse s1) (hgvs/parse s2)))
      "NM_005228:c.2361G>A" "NM_001346898:c.2361G>A"
      "NM_005228:c.2361G>A" "NM_005228:c.2361G>C"
      "NM_005228:c.2361G>A" "c.2361G>A"
      "p.L858R"             "p.L858M")
    (is (false? (hgvs/== (hgvs/parse "NM_005228:c.2361G>A") nil)))
    (is (false? (hgvs/== nil (hgvs/parse "NM_005228:c.2361G>A")))))

  (testing "more args"
    (is (true? (hgvs/== (hgvs/parse "NM_005228:c.2361G>A")
                        (hgvs/parse "NM_005228.3:c.2361G>A")
                        (hgvs/parse "NM_005228.4:c.2361G>A"))))
    (is (false? (hgvs/== (hgvs/parse "NM_005228:c.2361G>A")
                         (hgvs/parse "NM_005228.3:c.2361G>A")
                         (hgvs/parse "NM_005228.4:c.2361G>C"))))))

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
      "p.I327Rfs*?" "p.Ile327Argfs"
      "p.*110Glnext*17" "p.Ter110Glnext*17"))
  (testing "throws exception"
    (are [s] (thrown? #?(:clj Exception, :cljs js/Error) (hgvs/normalize s))
      ":2361G>A"
      "NM_005228.3:2361G>A"
      ""
      nil)))

(deftest repair-hgvs-str-test
  (are [s e] (let [s* (hgvs/repair-hgvs-str s)]
               (and (= s* e) (some? (hgvs/parse s*))))
    ;; remove-unnecessary-chars
    "NM_00001.1:c.123_124 delinsCGT" "NM_00001.1:c.123_124delinsCGT"
    "p.Asn3450Valfs*54\""            "p.Asn3450Valfs*54"

    ;; dedupe-colons
    "NM_00001.1::c.123T>G" "NM_00001.1:c.123T>G"

    ;; remove-trailing-period
    "c.123T>G." "c.123T>G"

    ;; upper-case-ins
    "c.673-43_673-33del11inscccagagc" "c.673-43_673-33delinsCCCAGAGC"

    ;; remove-gene-symbol
    "NM_000492.3(CFTR):c.1585-1G>A" "NM_000492.3:c.1585-1G>A"

    ;; dedupe-kinds
    "c.c.123T>G" "c.123T>G"

    ;; deletion->del
    "p.Leu128deletion" "p.Leu128del"

    ;; indel->delins
    "c.123_124indelCTGA" "c.123_124delinsCTGA"

    ;; frameshift->fs
    "p.Gln410frameshift" "p.Gln410fs"

    ;; stop->ter
    "p.Gln13Stop" "p.Gln13Ter"
    "p.Tyr9stop"  "p.Tyr9Ter"

    ;; substitutions->inv
    "c.1234_1235CA>TG"	"c.1234_1235inv"

    ;; substitution->delins
    "c.123G>CCACGTG" "c.123delGinsCCACGTG"

    ;; substitutions->delins
    "c.123_124>T"    "c.123_124delinsT"
    "c.123_124GC>AA" "c.123_124delGCinsAA"

    ;; remove-assembly
    "c.123_124delCT[hg19]" "c.123_124delCT"

    ;; remove-affected-count
    "c.123_234del14" "c.123_234del"
    "c.123_125dup3" "c.123_125dup"
    "c.123_125inv3" "c.123_125inv"

    ;; remove-same-end
    "c.123_123delA"       "c.123delA"
    "c.123_123dupA"       "c.123dupA"
    "c.123_123delAinsTAC" "c.123delAinsTAC"
    "g.123_123[14]"       "g.123[14]"

    ;; remove-alternative
    "c.1902dup(1897_1898insA)"  "c.1902dup"
    "c.1902dupA(1897_1898insA)" "c.1902dupA"

    ;; remove-inv-bases
    "c.123_124invTG" "c.123_124inv"

    ;; remove-del-count-from-delins
    "c.123_124del2insCTGA" "c.123_124delinsCTGA"

    ;; replace-repeated-seqs-parens1
    "c.112GAT(14)" "c.112GAT[14]"

    ;; replace-repeated-seqs-parens2
    "c.112GAT(14_16)" "c.112GAT[(14_16)]"

    ;; remove-genomic-bases-from-protein
    "p.G307S:GGC>AGC" "p.G307S"

    ;; remove-same-amino-acid
    "p.Phe269Phe="   "p.Phe269="
    "p.Gly413Gly=fs" "p.Gly413=fs"

    ;; replace-protein-no-change
    "p.*189*" "p.*189="

    ;; remove-extra-bases-from-protein-del
    "p.Phe53delPhe" "p.Phe53del"

    ;; protein-substitution->delins
    "p.N771>KL" "p.N771delinsKL"

    ;; protein-substitutions->delins
    "p.E746_S752>V" "p.E746_S752delinsV"

    ;; fix-protein-repeated-seqs-pos
    "p.348_349SerPro[4]" "p.Ser348_Pro349[4]"
    "p.68_70AAP[1]"      "p.A68_P70[1]"

    ;; frameshift-x->ter
    "p.G72AfsX13" "p.G72Afs*13"
    "p.Q94Hfsx?"  "p.Q94Hfs*?"

    ;; lower-case-fs
    "p.S896KFS*7" "p.S896Kfs*7"

    ;; remove-fs-greater
    "p.R123fs*>51" "p.R123fs*51"

    ;; protein-ter-substitution->ext
    "p.*320L" "p.*320Lext*?"

    ;; frameshift->ext
    "p.*833fs?" "p.*833ext*?"
    "p.*833fs"  "p.*833ext*?"

    ;; ext-ter->downstream
    "p.Ter397ThrextTer?" "p.Ter397Thrext*?"

    ;; ext->ins
    "p.S733_*734insS" "p.*734Sext*?")

  (are [s] (let [s* (hgvs/repair-hgvs-str s)]
             (and (= s* s) (some? (hgvs/parse s*))))
    "NM_005228.3:c.2361G>A"
    "c.2361G>A"
    "g.[2376A>C;3103del]"
    "NC_000022.11:g.28703511delA"
    "NM_004380.2:c.86-1G>T"
    "NM_000000.1:r.76a>c"
    "r.673-43_673-33delinscccagagc"
    "NM_004006.1:r.0"
    "LRG_199t1:r.?"
    "NP_005219.2:p.Leu858Arg"
    "NP_005219.2:p.L858R"
    "NP_001096.1:p.Arg258="
    "NP_001005735.1:p.Leu344Trpfs"
    "NC_012920.1:m.16563_13del"
    "J01749.1:o.4344_197dup"
    "ENST00000331920.11:c.34G>A")

  (let [lower-case-ext (fn [s kind]
                         (if (= kind :protein)
                           (string/replace s #"EXT" "ext")
                           s))
        my-repairers (conj repairer/built-in-repairers lower-case-ext)]
    (is (= (hgvs/repair-hgvs-str "NP_000000.1:p.*833EXT*?.")
           "NP_000000.1:p.*833EXT*?"))
    (is (= (hgvs/repair-hgvs-str "NP_000000.1:p.*833EXT*?." my-repairers)
           "NP_000000.1:p.*833ext*?")))

  (is (thrown? #?(:clj Exception, :cljs js/Error) (hgvs/repair-hgvs-str nil))))

#?(:clj (deftest print-test
          (let [expect (str "#clj-hgvs/hgvs \"" hgvs1s "\"")]
            (is (= (pr-str hgvs1m) expect))
            (is (= (pp/write hgvs1m :stream nil) expect)))))
