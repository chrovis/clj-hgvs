(ns clj-hgvs.repairer-test
  (:require [clojure.test :refer [are deftest]]
            [clj-hgvs.repairer :as repairer]))

(deftest infer-kind-test
  (are [s e] (= (repairer/infer-kind s) e)
    ;; valid HGVS
    "NC_000022.11:g.28703511delA" :genome
    "NC_012920.1:m.16563_13del"   :mitochondria
    "NM_005228.3:c.2361G>A"       :coding-dna
    "NM_005228:c.2361G>A"         :coding-dna
    "c.2361G>A"                   :coding-dna
    "J01749.1:o.4344_197dup"      :circular-dna
    "LRG_199t1:r.?"               :rna
    "NP_005219.2:p.L858R"         :protein

    ;; invalid HGVS
    "NC_000022.11:g.28703511deletionA" :genome
    "NM_000492.3(CFTR):c.1585-1G>A"    :coding-dna
    "c.123T>G."                        :coding-dna
    "c.673-43_673-33del11inscccagagc"  :coding-dna
    "c.c.123T>G"                       :coding-dna
    "c.123_124delCT[hg19]"             :coding-dna
    "c.1902dup(1897_1898insA)"         :coding-dna
    "p.Gln13Stop"                      :protein)

  (are [s] (nil? (repairer/infer-kind s))
    "2361G>A"
    "cc.2361G>A"
    "NC_000022.11:x.28703511delA"
    "p."
    ""))
