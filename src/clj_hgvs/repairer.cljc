(ns clj-hgvs.repairer
  "HGVS repair functions."
  (:require [clojure.string :as string]
            [clj-hgvs.internal :as intl]))

(defn infer-kind
  "Infers a kind from the provided string, returning the kind keyword."
  [s]
  (when-let [[_ kind] (re-matches #"(?:[^:]+:)?([gmcnorp])\..+" s)]
    (intl/->kind-keyword kind)))

;;; Repair rules

;; NM_00001.1:c.123_124 delinsCGT -> NM_00001.1:c.123_124delinsCGT
;; p.Asn3450Valfs*54" -> p.Asn3450Valfs*54
(defn ^:no-doc remove-unnecessary-chars
  [s _]
  (string/replace s #"[\s\"]" ""))

;; NM_00001.1::c.123T>G -> NM_00001.1:c.123T>G
(defn ^:no-doc dedupe-colons
  [s _]
  (string/replace s #"::" ":"))

;; c.123T>G. -> c.123T>G
(defn ^:no-doc remove-trailing-period
  [s _]
  (string/replace s #"\.$" ""))

;; c.673-43_673-33del11inscccagagccca -> c.673-43_673-33del11insCCCAGAGCCCA
(defn ^:no-doc upper-case-ins
  [s kind]
  (if (#{:genome :mitochondria :coding-dna :non-coding-dna :circular-dna} kind)
    (string/replace s
                    #"ins([acgnt]+)$"
                    (fn [[_ match]]
                      (str "ins" (string/upper-case match))))
    s))

;; NM_000492.3(CFTR):c.1585-1G>A -> NM_000492.3:c.1585-1G>A
(defn ^:no-doc remove-gene-symbol
  [s _]
  (string/replace s #"\([\dA-Z]+\):" ":"))

;; c.c.123T>G -> c.123T>G
(defn ^:no-doc dedupe-kinds
  [s _]
  (string/replace s #"([gmcnorp]\.){2,}" "$1"))

;; p.Leu128deletion -> p.Leu128del
(defn ^:no-doc deletion->del
  [s _]
  (string/replace s #"deletion" "del"))

;; c.123_124indelCTGA -> c.123_124delinsCTGA
(defn ^:no-doc indel->delins
  [s _]
  (string/replace s #"indel" "delins"))

;; p.Gln410frameshift -> p.Gln410fs
(defn ^:no-doc frameshift->fs
  [s kind]
  (if (= kind :protein)
    (string/replace s #"(F|f)rameshift" "fs")
    s))

;; p.Gln13Stop -> p.Gln13Ter
;; p.Tyr9stop -> p.Tyr9Ter
(defn ^:no-doc stop->ter
  [s kind]
  (if (= kind :protein)
    (string/replace s #"(S|s)top" "Ter")
    s))

;; c.123G>CCACGTG -> c.123delGinsCCACGTG
(defn ^:no-doc substitution->delins
  [s kind]
  (if (#{:genome :mitochondria :coding-dna :non-coding-dna :circular-dna} kind)
    (string/replace s #"([-\d+*]+)([A-Z]+)>([A-Z]{2,})" "$1del$2ins$3")
    s))

;; c.123_124>T -> c.123_124delinsT
;; c.123_124GC>AA -> c.123_124delGCinsAA
(defn ^:no-doc substitutions->delins
  [s kind]
  (if (#{:genome :mitochondria :coding-dna :non-coding-dna :circular-dna} kind)
    (string/replace s #"([-\d+*]+_[-\d+*]+)([A-Z]+)?>([A-Z]+)" "$1del$2ins$3")
    s))

;; c.233_234TC>CT	-> c.233_234delinsCT
(defn ^:no-doc substitutions->delins*
  [s kind]
  (if (= :coding-dna kind)
    (string/replace s #"([A-Z]?[-\d+*]+_[A-Z]?[-\d+*]+)([A-Z]+)?>([A-Z]+)" "$1delins$3")))

;; c.123_124delCT[hg19] -> c.123_124delCT
(defn ^:no-doc remove-assembly
  [s _]
  (string/replace s #"\[hg\d+\]$" ""))

;; c.123_234del14 -> c.123_234del
;; c.123_125dup3 -> c.123_125dup
;; c.123_125inv3 -> c.123_125inv
(defn ^:no-doc remove-affected-count
  [s kind]
  (if (#{:genome :mitochondria :coding-dna :non-coding-dna :circular-dna :rna}
       kind)
    (string/replace s #"(del|dup|inv)\d+$" "$1")
    s))

;; c.2210_2211CA>TG	-> c.2210_2211inv
(defn ^:no-doc basis->inv
  [s kind]
  (let [base-pairs {\A \T \T \A \C \G \G \C}
        [_ region b a] (re-find #"(\d+_\d+)([A-Z]+)>([A-Z]+)" s)]
    (if (= a
           (->> b
                (map base-pairs)
                reverse
                (apply str)))
      (str "c." region "inv")
      s)))

;; c.123_123delAinsTAC -> c.123delAinsTAC
;; g.123_123[14] -> g.123[14]
(defn ^:no-doc remove-same-end
  [s kind]
  (if (#{:genome :mitochondria :coding-dna :non-coding-dna :circular-dna :rna}
       kind)
    (string/replace s
                    #"([-\d+*]+)_([-\d+*]+)(del|dup|\[([\d\(\)_]+)\]$)"
                    (fn [[s start end mut]]
                      (if (= start end)
                        (str start mut)
                        s)))
    s))

;; c.1902dupA(1897_1898insA) -> c.1902dupA
(defn ^:no-doc remove-alternative
  [s kind]
  (if (#{:genome :mitochondria :coding-dna :non-coding-dna :circular-dna} kind)
    (string/replace s #"(dup[\dA-Z]*)\(.+\)" "$1")
    s))

;; c.123_124invTG -> c.123_124inv
(defn ^:no-doc remove-inv-bases
  [s kind]
  (if (#{:genome :mitochondria :coding-dna :non-coding-dna :circular-dna :rna}
       kind)
    (string/replace s #"inv[ACGNT]+$" "inv")
    s))

;; c.123_124del2insCTGA -> c.123_124delinsCTGA
(defn ^:no-doc remove-del-count-from-delins
  [s kind]
  (if (#{:genome :mitochondria :coding-dna :non-coding-dna :circular-dna :rna}
       kind)
    (string/replace s #"del\d+ins" "delins")
    s))

;; c.112GAT(14) -> c.112GAT[14]
(defn ^:no-doc replace-repeated-seqs-parens1
  [s kind]
  (if (#{:genome :mitochondria :coding-dna :non-coding-dna :circular-dna :rna}
       kind)
    (string/replace s #"\((\d+)\)$" "[$1]")
    s))

;; c.112GAT(14_16) -> c.112GAT[(14_16)]
(defn ^:no-doc replace-repeated-seqs-parens2
  [s kind]
  (if (#{:genome :mitochondria :coding-dna :non-coding-dna :circular-dna :rna}
       kind)
    (string/replace s #"(\(\d+_\d+\))$" "[$1]")
    s))

;; p.G307S:GGC>AGC -> p.G307S
(defn ^:no-doc remove-genomic-bases-from-protein
  [s kind]
  (if (= kind :protein)
    (string/replace s #":[A-Z]{3}>[A-Z]{3}$" "")
    s))

;; p.Phe269Phe= -> p.Phe269=
;; p.Gly413Gly=fs -> p.Gly413=fs
(defn ^:no-doc remove-same-amino-acid
  [s kind]
  (if (= kind :protein)
    (string/replace s #"[A-Z]([a-z]{2})?=" "=")
    s))

;; p.*189* -> p.*189=
(defn ^:no-doc replace-protein-no-change
  [s kind]
  (if (= kind :protein)
    (string/replace s
                    #"([A-Z*](?:[a-z]{2})?)(\d+)([A-Z*](?:[a-z]{2})?)$"
                    (fn [[s ref coord alt]]
                      (if (= ref alt)
                        (str ref coord "=")
                        s)))
    s))

;; p.Phe53delPhe -> p.Phe53del
(defn ^:no-doc remove-extra-bases-from-protein-del
  [s kind]
  (if (= kind :protein)
    (string/replace s #"del([A-Z]([a-z]{2})?)+$" "del")
    s))

;; p.N771>KL -> p.N771delinsKL
(defn ^:no-doc protein-substitution->delins
  [s kind]
  (if (= kind :protein)
    (string/replace s #"([0-9A-Za-z]+)>([A-Za-z]{2,})" "$1delins$2")
    s))

;; p.E746_S752>V -> p.E746_S752delinsV
(defn ^:no-doc protein-substitutions->delins
  [s kind]
  (if (= kind :protein)
    (string/replace s #"([0-9A-Za-z]+_[0-9A-Za-z]+)>([A-Za-z]+)" "$1delins$2")
    s))

;; p.348_349SerPro[4] -> p.Ser348_Pro349[4]
;; p.68_70AAP[1] -> p.A68_P70[1]
(defn ^:no-doc fix-protein-repeated-seqs-pos
  [s kind]
  (if (= kind :protein)
    (string/replace s
                    #"(\d+)_(\d+)([A-Z](?:[a-z]{2})?)(?:[A-Z](?:[a-z]{2})?)*([A-Z](?:[a-z]{2})?)(\[[\d()_]+\])"
                    "$3$1_$4$2$5")
    s))

;; p.G72AfsX13 -> p.G72Afs*13
;; p.Q94Hfsx? -> p.Q94Hfs*?
(defn ^:no-doc frameshift-x->ter
  [s kind]
  (if (= kind :protein)
    (string/replace s #"fs[Xx](\?|\d+)$" "fs*$1")
    s))

;; p.S896KFS*7 -> p.S896Kfs*7
(defn ^:no-doc lower-case-fs
  [s kind]
  (if (= kind :protein)
    (string/replace s #"FS((\*|Ter)(\?|\d+))$" "fs$1")
    s))

;; p.R123fs*>51 -> p.R123fs*51
(defn ^:no-doc remove-fs-greater
  [s kind]
  (if (= kind :protein)
    (string/replace s #"fs(\*|Ter)>(\?|\d+)" "fs$1$2")
    s))

;; p.*320L -> p.*320Lext*?
(defn ^:no-doc protein-ter-substitution->ext
  [s kind]
  (if (= kind :protein)
    (string/replace s #"((\*|Ter)\d+[A-Z]([a-z]{2})?)$" "$1ext*?")
    s))

;; p.*833fs? -> p.*833ext*?
(defn ^:no-doc frameshift->ext
  [s kind]
  (if (= kind :protein)
    (string/replace s
                    #"((?:\*|Ter)\d+(?:[A-Z](?:[a-z]{2})?)?)fs(?:\*|Ter)?(\?|\d+)?"
                    (fn [[_ ter new-site]]
                      (str ter "ext*" (or new-site "?"))))
    s))

;; p.Ter397ThrextTer? -> p.Ter397Thrext*?
(defn ^:no-doc ext-ter->downstream
  [s kind]
  (if (= kind :protein)
    (string/replace s #"extTer(\?|\d+)$" "ext*$1")
    s))

(def built-in-repairers
  "The built-in repair functions.

  A repair fn must take a HGVS string and an inferred kind, and return a
  repaired HGVS string. e.g.

    (defn lower-case-ext
      [s kind]
      (if (= kind :protein)
        (clojure.string/replace s #\"EXT\" \"ext\")
        s))

  The built-in rules are based on frequent mistakes in popular public-domain
  databases such as dbSNP and ClinVar."
  [remove-unnecessary-chars
   dedupe-colons
   remove-trailing-period
   upper-case-ins
   remove-gene-symbol
   dedupe-kinds
   deletion->del
   indel->delins
   frameshift->fs
   stop->ter
   substitution->delins
   substitutions->delins
   remove-assembly
   remove-affected-count
   remove-same-end
   remove-alternative
   remove-inv-bases
   remove-del-count-from-delins
   replace-repeated-seqs-parens1
   replace-repeated-seqs-parens2
   remove-genomic-bases-from-protein
   remove-same-amino-acid
   replace-protein-no-change
   remove-extra-bases-from-protein-del
   protein-substitution->delins
   protein-substitutions->delins
   fix-protein-repeated-seqs-pos
   frameshift-x->ter
   lower-case-fs
   remove-fs-greater
   protein-ter-substitution->ext
   frameshift->ext
   ext-ter->downstream])
