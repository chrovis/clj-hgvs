(ns clj-hgvs.mutation
  "Data structures and functions to handle HGVS mutations."
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.string :as string]
            [clojure.walk :as walk]
            [clj-hgvs.coordinate :as coord]
            [clj-hgvs.internal :refer [parse-long ->kind-keyword ->kind-str]]))

(declare parse-dna parse-rna parse-protein common-mutations?)

(def short-amino-acids
  "A list of single-letter amino acids."
  ["A"
   "R"
   "N"
   "D"
   "C"
   "Q"
   "E"
   "G"
   "H"
   "I"
   "L"
   "K"
   "M"
   "F"
   "P"
   "S"
   "T"
   "W"
   "Y"
   "V"
   "*"
   "X"])

(def long-amino-acids
  "A list of three-letter amino acids."
  ["Ala"
   "Arg"
   "Asn"
   "Asp"
   "Cys"
   "Gln"
   "Glu"
   "Gly"
   "His"
   "Ile"
   "Leu"
   "Lys"
   "Met"
   "Phe"
   "Pro"
   "Ser"
   "Thr"
   "Trp"
   "Tyr"
   "Val"
   "Ter"
   "Xaa"])

(def ^:private short-long-amino-acid-map
  (zipmap short-amino-acids long-amino-acids))

(def ^:private long-short-amino-acid-map
  (zipmap long-amino-acids short-amino-acids))

(defn ->long-amino-acid
  "Converts a single-letter amino acid into a three-letter one. s must be String
  or Character. Returns nil if s is not present in amino acid list."
  [s]
  (let [s (cond-> s (char? s) str)]
    (if ((set long-amino-acids) s)
      s
      (get short-long-amino-acid-map s))))

(defn ->short-amino-acid
  "Converts a three-letter amino acid into a single-letter one. s must be String
  or Character. Returns nil if s is not present in amino acid list."
  [s]
  (let [s (cond-> s (char? s) str)]
    (if ((set short-amino-acids) s)
      s
      (get long-short-amino-acid-map s))))

;; "1" => ["1"]
;; "1_2" => ["1" "2"]
;; "(1_2)" => ["(1_2)"]
;; "(1_2)_(3_4)" => ["(1_2)" "(3_4)"]
(defn- split-coord-range
  [s]
  (vec (condp re-find s
         #"\)_\(" (re-seq #"\([\*\-\+\d\?_]+\)" s)
         #"\)_" (next (re-find #"(\([\*\-\+\d\?_]+\))_([\*\-\+\d\?]+)" s))
         #"_\(" (next (re-find #"([\*\-\+\d\?]+)_(\([\*\-\+\d\?_]+\))" s))
         #"\([\*\-\+\d\?_]+\)" [s]
         (re-seq #"[\*\-\+\d\?]+" s))))

;; "[123G>A];[345del]" => ["123G>A" "345del"]
;; "2376[G>C];[G>C]" => ["2376G>C" "2376G>C"]
(defn- split-alleles-mutations
  [s]
  (let [[_ common alleles] (re-matches #"([^\[\];]+)(\[.+\])" s)]
    (if common
      (if (re-matches #"\[\d+\];\[\d+\]" alleles)
        (->> (string/split alleles #";")
             (mapv #(str common %)))
        (->> (re-seq #"\[(.+?)\]" alleles)
             (map second)
             (mapv #(str common %))))
      (mapv second (re-seq #"\[(.+?)\]" s)))))

;; "14" => 14
;; "(600_800)" => [600 800]
(defn- parse-ncopy
  [s]
  (if-let [[_ s1 s2] (re-matches #"\((\d+)_(\d+)\)" s)]
    [(parse-long s1) (parse-long s2)]
    (parse-long s)))

;; 14 => "14"
;; [600 800] => "(600_800)"
(defn- format-ncopy
  [x]
  (cond
    (vector? x) (str "(" (first x) "_" (second x) ")")
    (integer? x) (str x)))

(defn- plain-coords
  [m]
  (walk/prewalk (fn [x]
                  (if (and (vector? x) (satisfies? coord/Coordinate (second x)))
                    (update x 1 coord/plain)
                    x))
                m))

(defn- restore-coords
  [m]
  (walk/prewalk (fn [x]
                  (if (and (vector? x) (:coordinate (second x)))
                    (update x 1 coord/restore)
                    x))
                m))

(defprotocol Mutation
  (format [this] [this opts]
    "Returns a string representing the given mutation. The second argument is an
  optional map to specify style. See document of clj-hgvs.core/format for
  details of the option.")
  (plain [this] "Returns a plain map representing the given mutation."))

(defmulti restore
  "Restores a plain map to a suitable mutation record."
  {:arglists '([m])}
  :mutation)

;; SeparatelyFormat protocol is used for formatting alleles mutation, such as
;; c.2376[G>C];[G>C], including common part.
(defprotocol SeparatelyFormat
  (format-common [this opts])
  (format-unique [this opts]))

;;; DNA mutations

;; See http://varnomen.hgvs.org/bg-material/standards#dna
(defn- dna-bases?
  [s]
  (and (string? s) (some? (re-matches #"[ACGTBDHKMNRSVWY]+" s))))

(defn- coord-parser
  [kind]
  (case kind
    :genome coord/parse-genomic-coordinate
    :mitochondria coord/parse-mitochondrial-coordinate
    :cdna coord/parse-cdna-coordinate
    :ncdna coord/parse-ncdna-coordinate))

;;; DNA - substitution
;;;
;;; e.g. g.45576A>C
;;;      c.88+1G>T
;;;      c.123G=
;;;      c.85C=/>T
;;;      c.85C=//>T

(defrecord DNASubstitution [coord ref type alt]
  Mutation
  (format [this] (format this nil))
  (format [this _]
    (str (format-common this _) (format-unique this _)))
  (plain [this]
    (into {:mutation "dna-substitution"} (plain-coords this)))
  SeparatelyFormat
  (format-common [this _] (coord/format coord))
  (format-unique [this _] (str ref type alt)))

(defn dna-substitution
  "Constructor of DNASubstitution. Throws an exception if any input is illegal."
  ([coord ref typ] (dna-substitution coord ref typ nil))
  ([coord ref typ alt]
   {:pre [(satisfies? coord/Coordinate coord)
          (dna-bases? ref)
          (#{">" "=" "=/>" "=//>"} typ)
          (or (nil? alt) (dna-bases? alt))]}
   (DNASubstitution. coord ref typ alt)))

(def ^:private dna-substitution-re
  #"([\d\-\+\*\?]+)([A-Z])([>=/]+)([A-Z])?")

(defn parse-dna-substitution
  [s kind]
  (let [[_ coord ref typ alt] (re-matches dna-substitution-re s)
        parse-coord (coord-parser kind)]
    (dna-substitution (parse-coord coord) ref typ alt)))

(defmethod restore "dna-substitution"
  [m]
  (let [{:keys [coord ref type alt]} (restore-coords m)]
    (dna-substitution coord ref type alt)))

;;; DNA - deletion
;;;
;;; e.g. g.7del
;;;      g.6_8del
;;;      g.6_8delTGC
;;;      c.120_123+48del
;;;      c.(4071+1_4072-1)_(5145+1_5146-1)del
;;;      c.(?_-30)_(12+1_13-1)del
;;;      c.(?_-1)_(*1_?)del

(defrecord DNADeletion [coord-start coord-end ref]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [show-bases?] :or {show-bases? false}}]
    (apply str (flatten [(if (vector? coord-start)
                           ["(" (coord/format (first coord-start))
                            "_" (coord/format (second coord-start)) ")"]
                           (coord/format coord-start))
                         (if coord-end "_")
                         (if (vector? coord-end)
                           ["(" (coord/format (first coord-end))
                            "_" (coord/format (second coord-end)) ")"]
                           (some-> coord-end coord/format))
                         "del"
                         (if show-bases? ref)])))
  (plain [this]
    (into {:mutation "dna-deletion"} (plain-coords this))))

(defn dna-deletion
  "Constructor of DNADeletion. Throws an exception if any input is illegal."
  ([coord-start coord-end] (dna-deletion coord-start coord-end nil))
  ([coord-start coord-end ref]
   {:pre [(if (vector? coord-start)
            (every? #(satisfies? coord/Coordinate %) coord-start)
            (satisfies? coord/Coordinate coord-start))
          (or (nil? coord-end)
              (if (vector? coord-start)
                (every? #(satisfies? coord/Coordinate %) coord-end)
                (satisfies? coord/Coordinate coord-end)))
          (or (nil? ref) (dna-bases? ref))]}
   (DNADeletion. coord-start coord-end ref)))

(def ^:private dna-deletion-re
  #"([\(\)\*\-\+\d\?_]+)del([A-Z]+)?")

(defn parse-dna-deletion
  [s kind]
  (let [[_ coord ref] (re-matches dna-deletion-re s)
        parse-coord (coord-parser kind)
        [coord-s coord-e] (split-coord-range coord)]
    (dna-deletion (parse-coord coord-s)
                  (some-> coord-e parse-coord)
                  ref)))

(defmethod restore "dna-deletion"
  [m]
  (let [{:keys [coord-start coord-end ref]} (restore-coords m)]
    (dna-deletion coord-start coord-end ref)))

;;; DNA - duplication
;;;
;;; e.g. g.7dup
;;;      g.6_8dup
;;;      g.6_8dupTGC
;;;      c.120_123+48dup
;;;      c.(4071+1_4072-1)_(5145+1_5146-1)dup
;;;      c.(?_-30)_(12+1_13-1)dup
;;;      c.(?_-1)_(*1_?)dup

(defrecord DNADuplication [coord-start coord-end ref]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [show-bases?] :or {show-bases? false}}]
    (apply str (flatten [(if (vector? coord-start)
                           ["(" (coord/format (first coord-start))
                            "_" (coord/format (second coord-start)) ")"]
                           (coord/format coord-start))
                         (if (and (some? coord-end) (not= coord-start coord-end))
                           ["_"
                            (if (vector? coord-end)
                              ["(" (coord/format (first coord-end))
                               "_" (coord/format (second coord-end)) ")"]
                              (coord/format coord-end))])
                         "dup"
                         (if show-bases? ref)])))
  (plain [this]
    (into {:mutation "dna-duplication"} (plain-coords this))))

(defn dna-duplication
  "Constructor of DNADuplication. Throws an exception if any input is illegal."
  ([coord-start coord-end] (dna-duplication coord-start coord-end nil))
  ([coord-start coord-end ref]
   {:pre [(if (vector? coord-start)
            (every? #(satisfies? coord/Coordinate %) coord-start)
            (satisfies? coord/Coordinate coord-start))
          (or (nil? coord-end)
              (if (vector? coord-start)
                (every? #(satisfies? coord/Coordinate %) coord-end)
                (satisfies? coord/Coordinate coord-end)))
          (or (nil? ref) (dna-bases? ref))]}
   (DNADuplication. coord-start coord-end ref)))

(def ^:private dna-duplication-re
  #"([\(\)\*\-\+\d\?_]+)dup([A-Z]+)?")

(defn parse-dna-duplication
  [s kind]
  (let [[_ coord ref] (re-matches dna-duplication-re s)
        parse-coord (coord-parser kind)
        [coord-s coord-e] (split-coord-range coord)]
    (dna-duplication (parse-coord coord-s)
                     (some-> coord-e parse-coord)
                     ref)))

(defmethod restore "dna-duplication"
  [m]
  (let [{:keys [coord-start coord-end ref]} (restore-coords m)]
    (dna-duplication coord-start coord-end ref)))

;;; DNA - insertion
;;;
;;; e.g. g.5756_5757insAGG
;;;      g.123_124insL37425.1:23_361
;;;      g.122_123ins123_234inv (TODO)
;;;      g.122_123ins213_234invinsAins123_211inv (TODO)
;;;      g.549_550insN
;;;      g.1134_1135ins(100)
;;;      g.?_?insNC_000023.10:(12345_23456)_(34567_45678)

(defrecord DNAInsertion [coord-start coord-end alt]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [ins-format] :or {ins-format :auto}}]
    (apply str (flatten [(coord/format coord-start)
                         "_"
                         (coord/format coord-end)
                         "ins"
                         (cond
                           (string? alt) (case ins-format
                                           :auto (if (and (every? #(= % \N) alt)
                                                          (>= (count alt) 10))
                                                   (str "(" (count alt) ")")
                                                   alt)
                                           :bases alt
                                           :count (str "(" (count alt) ")"))
                           (map? alt) [(:transcript alt)
                                       ":"
                                       (coord/format (:coord-start alt))
                                       "_"
                                       (coord/format (:coord-end alt))])])))
  (plain [this]
    (into {:mutation "dna-insertion"} (plain-coords this))))

(defn dna-insertion
  "Constructor of DNAInsertion. Throws an exception if any input is illegal."
  [coord-start coord-end alt]
  {:pre [(satisfies? coord/Coordinate coord-start)
         (satisfies? coord/Coordinate coord-end)
         (or (dna-bases? alt) (map? alt))]}
  (DNAInsertion. coord-start coord-end alt))

(defn- parse-dna-insertion-alt
  [s kind]
  (or (re-matches #"[A-Z]+" s)
      (some-> (re-matches #"\((\d+)\)" s)
              (second)
              (parse-long)
              (repeat "N")
              (#(apply str %)))
      (let [[_ transcript coord] (re-matches #"(?:([^:]+):)([\(\)\*\-\+\d\?_]+)" s)
            parse-coord (coord-parser kind)
            [coord-s coord-e] (split-coord-range coord)]
        {:transcript transcript
         :coord-start (parse-coord coord-s)
         :coord-end (parse-coord coord-e)})))

(def ^:private dna-insertion-re
  #"([\d\-\+\*\?]+)_([\d\-\+\*\?]+)ins(.+)")

(defn parse-dna-insertion
  [s kind]
  (let [[_ coord-s coord-e alt] (re-matches dna-insertion-re s)
        parse-coord (coord-parser kind)]
    (dna-insertion (parse-coord coord-s)
                   (parse-coord coord-e)
                   (parse-dna-insertion-alt alt kind))))

(defmethod restore "dna-insertion"
  [m]
  (let [{:keys [coord-start coord-end alt]} (restore-coords m)]
    (dna-insertion coord-start coord-end alt)))

;;; DNA - inversion
;;;
;;; e.g. g.1077_1080inv
;;;      c.77_80inv

(defrecord DNAInversion [coord-start coord-end]
  Mutation
  (format [this] (format this nil))
  (format [this _]
    (str (coord/format coord-start)
         "_"
         (coord/format coord-end)
         "inv"))
  (plain [this]
    (into {:mutation "dna-inversion"} (plain-coords this))))

(defn dna-inversion
  "Constructor of DNAInversion. Throws an exception if any input is illegal."
  [coord-start coord-end]
  {:pre [(satisfies? coord/Coordinate coord-start)
         (satisfies? coord/Coordinate coord-end)
         (neg? (compare coord-start coord-end))]}
  (DNAInversion. coord-start coord-end))

(def ^:private dna-inversion-re
  #"([\d\-\+\*\?]+)_([\d\-\+\*\?]+)inv")

(defn parse-dna-inversion
  [s kind]
  (let [[_ coord-s coord-e] (re-matches dna-inversion-re s)
        parse-coord (coord-parser kind)]
    (dna-inversion (parse-coord coord-s) (parse-coord coord-e))))

(defmethod restore "dna-inversion"
  [m]
  (let [{:keys [coord-start coord-end]} (restore-coords m)]
    (dna-inversion coord-start coord-end)))

;;; DNA - conversion
;;;
;;; e.g. g.333_590con1844_2101
;;;      g.415_1655conAC096506.5:g.409_1683
;;;      c.15_355conNM_004006.1:20_360

(defrecord DNAConversion [coord-start coord-end alt]
  Mutation
  (format [this] (format this nil))
  (format [this _]
    (apply str (flatten [(coord/format coord-start)
                         "_"
                         (coord/format coord-end)
                         "con"
                         (if (:transcript alt) [(:transcript alt) ":"])
                         (if (:kind alt) [(->kind-str (:kind alt)) "."])
                         (coord/format (:coord-start alt))
                         "_"
                         (coord/format (:coord-end alt))])))
  (plain [this]
    (into {:mutation "dna-conversion"} (plain-coords this))))

(defn dna-conversion
  "Constructor of DNAConversion. Throws an exception if any input is illegal."
  [coord-start coord-end alt]
  {:pre [(satisfies? coord/Coordinate coord-start)
         (satisfies? coord/Coordinate coord-end)
         (neg? (compare coord-start coord-end))
         (map? alt)]}
  (DNAConversion. coord-start coord-end alt))

(def ^:private dna-conversion-re
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))con(.+)")

(def ^:private dna-conversion-alt-re
  #"(?:([^:]+):)?(?:([gmcn])\.)?([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))")

(defn- parse-dna-conversion-alt
  [s kind default-coord-parser]
  (let [[_ transcript kind coord-s coord-e] (re-matches dna-conversion-alt-re s)
        parse-alt-coord (if kind
                          (coord-parser (->kind-keyword kind))
                          default-coord-parser)]
    {:transcript transcript
     :kind (some-> kind ->kind-keyword)
     :coord-start (parse-alt-coord coord-s)
     :coord-end (parse-alt-coord coord-e)}))

(defn parse-dna-conversion
  [s kind]
  (let [[_ coord-s coord-e alt] (re-matches dna-conversion-re s)
        parse-coord (coord-parser kind)]
    (dna-conversion (parse-coord coord-s)
                    (parse-coord coord-e)
                    (parse-dna-conversion-alt alt kind parse-coord))))

(defmethod restore "dna-conversion"
  [m]
  (let [{:keys [coord-start coord-end alt]} (restore-coords m)]
    (dna-conversion coord-start coord-end alt)))

;;; DNA - indel
;;;
;;; e.g. g.6775delinsGA
;;;      g.6775delTinsGA
;;;      c.145_147delinsTGG

(defrecord DNAIndel [coord-start coord-end ref alt]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [show-bases?] :or {show-bases? false}}]
    (apply str (flatten [(coord/format coord-start)
                         (if (neg? (compare coord-start coord-end))
                           ["_" (coord/format coord-end)])
                         "del"
                         (if show-bases? ref)
                         "ins"
                         alt])))
  (plain [this]
    (into {:mutation "dna-indel"} (plain-coords this))))

(defn dna-indel
  "Constructor of DNAIndel. Throws an exception if any input is illegal."
  [coord-start coord-end ref alt]
  {:pre [(satisfies? coord/Coordinate coord-start)
         (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))
         (or (nil? ref) (dna-bases? ref))
         (dna-bases? alt)]}
  (DNAIndel. coord-start coord-end ref alt))

(def ^:private dna-indel-re
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))?del([A-Z]+)?ins([A-Z]+)")

(defn parse-dna-indel
  [s kind]
  (let [[_ coord-s coord-e ref alt] (re-matches dna-indel-re s)
        parse-coord (coord-parser kind)]
    (dna-indel (parse-coord coord-s) (some-> coord-e parse-coord) ref alt)))

(defmethod restore "dna-indel"
  [m]
  (let [{:keys [coord-start coord-end ref alt]} (restore-coords m)]
    (dna-indel coord-start coord-end ref alt)))

;;; DNA - alleles
;;;
;;; e.g. g.[123G>A;345del]
;;;      g.[123G>A];[345del]
;;;      c.2376[G>C];[G>C]
;;;      c.2376G>C(;)3103del (TODO)
;;;      c.2376[G>C];[(G>C)] (TODO)
;;;      c.2376[G>C];[=] (TODO)
;;;      c.[2376G>C];[?] (TODO)
;;;      c.[296T>G;476C>T;1083A>C];[296T>G;1083A>C]
;;;      c.[296T>G;476C>T];[476C>T](;)1083A>C (TODO)
;;;      c.[296T>G];[476C>T](;)1083G>C(;)1406del (TODO)
;;;      c.[NM_000167.5:94A>G;NM_004006.2:76A>C] (TODO)
;;;      g.123_124[14];[18]

(defrecord DNAAlleles [mutations1 mutations2]
  Mutation
  (format [this] (format this nil))
  (format [this opts]
    (if (common-mutations? (concat mutations1 mutations2))
      (str (format-common (first mutations1) opts)
           "[" (string/join ";" (map #(format-unique % opts) mutations1)) "]"
           (if (seq mutations2)
             (str ";[" (string/join ";" (map #(format-unique % opts) mutations2)) "]")))
      (str "[" (string/join ";" (map #(format % opts) mutations1)) "]"
           (if (seq mutations2)
             (str ";[" (string/join ";" (map #(format % opts) mutations2)) "]")))))
  (plain [this]
    {:mutation "dna-alleles"
     :mutations1 (mapv plain mutations1)
     :mutations2 (mapv plain mutations2)}))

(defn dna-alleles
  [mutations1 mutations2]
  {:pre [(every? #(satisfies? Mutation %) mutations1)
         (every? #(satisfies? Mutation %) mutations2)]}
  (DNAAlleles. mutations1 mutations2))

(defn parse-dna-alleles
  [s kind]
  (let [[mut1 mut2] (split-alleles-mutations s)]
    (dna-alleles (mapv #(parse-dna % kind) (string/split mut1 #";"))
                 (if mut2
                   (mapv #(parse-dna % kind) (string/split mut2 #";"))))))

(defmethod restore "dna-alleles"
  [m]
  (dna-alleles (mapv restore (:mutations1 m))
               (mapv restore (:mutations2 m))))

;;; DNA - repeated sequences
;;;
;;; e.g. g.123_124[14]
;;;      g.123TG[14]
;;;      c.-128_-126[(600_800)]

(defrecord DNARepeatedSeqs [coord-start coord-end ref ncopy]
  Mutation
  (format [this] (format this nil))
  (format [this opts]
    (str (format-common this opts) "[" (format-unique this opts) "]"))
  (plain [this]
    (into {:mutation "dna-repeated-seqs"} (plain-coords this)))
  SeparatelyFormat
  (format-common [this {:keys [range-format] :or {range-format :auto}}]
    (let [should-show-end? (neg? (compare coord-start coord-end))]
      (str (coord/format coord-start)
           (case range-format
             :auto (or ref (if should-show-end? (str "_" (coord/format coord-end))))
             :bases ref
             :coord (if should-show-end? (str "_" (coord/format coord-end)))))))
  (format-unique [this _]
    (format-ncopy ncopy)))

(defn dna-repeated-seqs
  "Constructor of DNARepeatedSeqs. Throws an exception if any input is illegal."
  [coord-start coord-end ref ncopy]
  {:pre [(satisfies? coord/Coordinate coord-start)
         (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))
         (or (nil? ref) (dna-bases? ref))
         (or (integer? ncopy) (vector? ncopy))]}
  (DNARepeatedSeqs. coord-start coord-end ref ncopy))

(def ^:private dna-repeated-seqs-re
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))?([A-Z]+)?\[([\d\(\)_]+)\]")

(defn parse-dna-repeated-seqs
  [s kind]
  (let [[_ coord-s coord-e ref ncopy] (re-matches dna-repeated-seqs-re s)
        parse-coord (coord-parser kind)]
    (dna-repeated-seqs (parse-coord coord-s)
                       (some-> coord-e parse-coord)
                       ref
                       (parse-ncopy ncopy))))

(defmethod restore "dna-repeated-seqs"
  [m]
  (let [{:keys [coord-start coord-end ref ncopy]} (restore-coords m)]
    (dna-repeated-seqs coord-start coord-end ref ncopy)))

(defn parse-dna
  "Parses a DNA mutation string s, returning a record implementing Mutation
  protocol. kind must be selected from :genome, :mitochondria, :cdna, or
  :ncdna."
  [s kind]
  ((condp re-find s
     #"\[.+;.+\]$" parse-dna-alleles
     #"del[A-Z]*ins" parse-dna-indel
     #"del" parse-dna-deletion
     #"dup" parse-dna-duplication
     #"ins" parse-dna-insertion
     #"inv" parse-dna-inversion
     #"con" parse-dna-conversion
     #"\[[\d\(\)_]+\]" parse-dna-repeated-seqs
     parse-dna-substitution)
   s kind))

;;; RNA mutations

;; See http://varnomen.hgvs.org/bg-material/standards#rna
(defn- rna-bases?
  [s]
  (and (string? s) (some? (re-matches #"[acgubdhkmnrsvwy]+" s))))

;;; RNA - substitution
;;;
;;; e.g. r.76a>c
;;;      r.-14g>c
;;;      r.*46u>a

(defrecord RNASubstitution [coord ref alt]
  Mutation
  (format [this] (format this nil))
  (format [this _]
    (str (format-common this _) (format-unique this _))
    (apply str (coord/format coord) ref ">" alt))
  (plain [this]
    (into {:mutation "rna-substitution"} (plain-coords this)))
  SeparatelyFormat
  (format-common [this _] (coord/format coord))
  (format-unique [this _] (str ref ">" alt)))

(defn rna-substitution
  "Constructor of RNASubstitution. Throws an exception if any input is illegal."
  [coord ref alt]
  {:pre [(satisfies? coord/Coordinate coord)
         (rna-bases? ref)
         (rna-bases? alt)]}
  (RNASubstitution. coord ref alt))

(def ^:private rna-substitution-re
  #"([\d\-\+\*]+)([a-z]?)>([a-z]?)")

(defn parse-rna-substitution
  [s]
  (let [[_ coord ref alt] (re-matches rna-substitution-re s)]
    (rna-substitution (coord/parse-rna-coordinate coord) ref alt)))

(defmethod restore "rna-substitution"
  [m]
  (let [{:keys [coord ref alt]} (restore-coords m)]
    (rna-substitution coord ref alt)))

;;; RNA - deletion
;;;
;;; e.g. r.7del
;;;      r.7delu
;;;      r.6_8del
;;;      r.(4072_5145)del (TODO)

(defrecord RNADeletion [coord-start coord-end ref]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [show-bases?] :or {show-bases? false}}]
    (str (coord/format coord-start)
         (if coord-end
           (str "_" (coord/format coord-end)))
         "del"
         (if show-bases? ref)))
  (plain [this]
    (into {:mutation "rna-deletion"} (plain-coords this))))

(defn rna-deletion
  "Constructor of DNAdeletion. Throws an exception if any input is illegal."
  ([coord-start coord-end] (rna-deletion coord-start coord-end nil))
  ([coord-start coord-end ref]
   {:pre [(satisfies? coord/Coordinate coord-start)
          (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))
          (or (nil? ref) (rna-bases? ref))]}
   (RNADeletion. coord-start coord-end ref)))

(def ^:private rna-deletion-re
  #"([\d\-\+\*]+)(?:_([\d\-\+\*]+))?del([a-z]+)?")

(defn parse-rna-deletion
  [s]
  (let [[_ coord-s coord-e ref] (re-matches rna-deletion-re s)]
    (rna-deletion (coord/parse-rna-coordinate coord-s)
                  (some-> coord-e coord/parse-rna-coordinate)
                  ref)))

(defmethod restore "rna-deletion"
  [m]
  (let [{:keys [coord-start coord-end ref]} (restore-coords m)]
    (rna-deletion coord-start coord-end ref)))

;;; RNA - duplication
;;;
;;; e.g. r.7dup
;;;      r.7dupu
;;;      r.6_8dup

(defrecord RNADuplication [coord-start coord-end ref]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [show-bases?] :or {show-bases? false}}]
    (str (coord/format coord-start)
         (if coord-end
           (str "_" (coord/format coord-end)))
         "dup"
         (if show-bases? ref)))
  (plain [this]
    (into {:mutation "rna-duplication"} (plain-coords this))))

(defn rna-duplication
  "Constructor of RNADuplication. Throws an exception if any input is illegal."
  ([coord-start coord-end] (rna-duplication coord-start coord-end nil))
  ([coord-start coord-end ref]
   {:pre [(satisfies? coord/Coordinate coord-start)
          (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))
          (or (nil? ref) (rna-bases? ref))]}
   (RNADuplication. coord-start coord-end ref)))

(def ^:private rna-duplication-re
  #"([\d\-\+\*]+)(?:_([\d\-\+\*]+))?dup([a-z]+)?")

(defn parse-rna-duplication
  [s]
  (let [[_ coord-s coord-e ref] (re-matches rna-duplication-re s)]
    (rna-duplication (coord/parse-rna-coordinate coord-s)
                     (some-> coord-e coord/parse-rna-coordinate)
                     ref)))

(defmethod restore "rna-duplication"
  [m]
  (let [{:keys [coord-start coord-end ref]} (restore-coords m)]
    (rna-duplication coord-start coord-end ref)))

;;; RNA - insertion
;;;
;;; e.g. r.756_757insacu
;;;      r.431_432ins(5)
;;;      r.123_124insL37425.1:23_361

(defrecord RNAInsertion [coord-start coord-end alt]
  Mutation
  (format [this] (format this nil))
  (format [this _]
    (str (coord/format coord-start)
         "_"
         (coord/format coord-end)
         "ins"
         (cond
           (map? alt) (str (:genbank alt) ":" (:coord-start alt) "_" (:coord-end alt))
           (re-matches #"n{2,}" alt) (str "(" (count alt) ")")
           :else alt)))
  (plain [this]
    (into {:mutation "rna-insertion"} (plain-coords this))))

(defn rna-insertion
  "Constructor of RNAInsertion. Throws an exception if any input is illegal."
  [coord-start coord-end alt]
  {:pre [(satisfies? coord/Coordinate coord-start)
         (satisfies? coord/Coordinate coord-end)
         (or (rna-bases? alt) (map? alt) (re-matches #"n{2,}" alt))]}
  (RNAInsertion. coord-start coord-end alt))

(defn- parse-rna-alt-n
  [s]
  (if-let [n (some-> (re-find #"\((\d)\)" s)
                     second
                     parse-long)]
    (apply str (repeat n "n"))))

(def ^:private genbank-re
  #"([A-Z]+[\d\.]+):(\d+)_(\d+)")

(defn- parse-rna-alt-genbank
  [s]
  (let [[_ genbank coord-s coord-e] (re-matches genbank-re s)]
    {:genbank genbank
     :coord-start (parse-long coord-s)
     :coord-end (parse-long coord-e)}))

(def ^:private rna-insertion-re
  #"([\d\-\+\*]+)_([\d\-\+\*]+)ins(.+)?")

(defn parse-rna-insertion
  [s]
  (let [[_ coord-s coord-e alt] (re-matches rna-insertion-re s)]
    (rna-insertion (coord/parse-rna-coordinate coord-s)
                   (some-> coord-e coord/parse-rna-coordinate)
                   (cond
                     (re-find #"[a-z]+" alt) alt
                     (re-find #"\(\d\)" alt) (parse-rna-alt-n alt)
                     :else (parse-rna-alt-genbank alt)))))

(defmethod restore "rna-insertion"
  [m]
  (let [{:keys [coord-start coord-end alt]} (restore-coords m)]
    (rna-insertion coord-start coord-end alt)))

;;; RNA - inversion
;;;
;;; e.g. r.177_180inv

(defrecord RNAInversion [coord-start coord-end]
  Mutation
  (format [this] (format this nil))
  (format [this _]
    (str (coord/format coord-start)
         "_"
         (coord/format coord-end)
         "inv"))
  (plain [this]
    (into {:mutation "rna-inversion"} (plain-coords this))))

(defn rna-inversion
  "Constructor of RNAInversion. Throws an exception if any input is illegal."
  [coord-start coord-end]
  {:pre [(satisfies? coord/Coordinate coord-start)
         (satisfies? coord/Coordinate coord-end)
         (neg? (compare coord-start coord-end))]}
  (RNAInversion. coord-start coord-end))

(def ^:private rna-inversion-re
  #"([\d\-\+\*]+)_([\d\-\+\*]+)inv")

(defn parse-rna-inversion
  [s]
  (let [[_ coord-s coord-e] (re-matches rna-inversion-re s)]
    (rna-inversion (coord/parse-rna-coordinate coord-s)
                   (coord/parse-rna-coordinate coord-e))))

(defmethod restore "rna-inversion"
  [m]
  (let [{:keys [coord-start coord-end]} (restore-coords m)]
    (rna-inversion coord-start coord-end)))

;;; RNA - conversion
;;;
;;; e.g. r.123_345con888_1110
;;;      r.415_1655conAC096506.5:409_1649

(defrecord RNAConversion [coord-start coord-end alt]
  Mutation
  (format [this] (format this nil))
  (format [this _]
    (str (coord/format coord-start)
         "_"
         (coord/format coord-end)
         "con"
         (some-> (:transcript alt) (str ":"))
         (coord/format (:coord-start alt))
         "_"
         (coord/format (:coord-end alt))))
  (plain [this]
    (into {:mutation "rna-conversion"} (plain-coords this))))

(defn rna-conversion
  "Constructor of RNAConversion. Throws an exception if any input is illegal."
  [coord-start coord-end alt]
  {:pre [(satisfies? coord/Coordinate coord-start)
         (satisfies? coord/Coordinate coord-end)
         (neg? (compare coord-start coord-end))
         (map? alt)]}
  (RNAConversion. coord-start coord-end alt))

(def ^:private rna-conversion-re
  #"([\d\-\+\*]+)_([\d\-\+\*]+)con(.+)")

(def ^:private rna-conversion-alt-re
  #"(?:([^:]+):)?([\d\-\+\*\?]+)_([\d\-\+\*\?]+)")

(defn- parse-rna-conversion-alt
  [s]
  (let [[_ transcript coord-s coord-e] (re-matches rna-conversion-alt-re s)]
    {:transcript transcript
     :coord-start (coord/parse-rna-coordinate coord-s)
     :coord-end (coord/parse-rna-coordinate coord-e)}))

(defn parse-rna-conversion
  [s]
  (let [[_ coord-s coord-e alt] (re-matches rna-conversion-re s)]
    (rna-conversion (coord/parse-rna-coordinate coord-s)
                    (coord/parse-rna-coordinate coord-e)
                    (parse-rna-conversion-alt alt))))

(defmethod restore "rna-conversion"
  [m]
  (let [{:keys [coord-start coord-end alt]} (restore-coords m)]
    (rna-conversion coord-start coord-end alt)))

;;; RNA - indel
;;;
;;; e.g. r.775delinsga
;;;      r.775deluinsga
;;;      r.775_777delinsc

(defrecord RNAIndel [coord-start coord-end ref alt]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [show-bases?] :or {show-bases? false}}]
    (str (coord/format coord-start)
         (if (neg? (compare coord-start coord-end))
           (str "_" (coord/format coord-end)))
         "del"
         (if show-bases? ref)
         "ins"
         alt))
  (plain [this]
    (into {:mutation "rna-indel"} (plain-coords this))))

(defn rna-indel
  "Constructor of RNAIndel. Throws an exception if any input is illegal."
  [coord-start coord-end ref alt]
  {:pre [(satisfies? coord/Coordinate coord-start)
         (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))
         (or (nil? ref) (rna-bases? ref))
         (rna-bases? alt)]}
  (RNAIndel. coord-start coord-end ref alt))

(def ^:private rna-indel-re
  #"([\d\-\+\*]+)(?:_([\d\-\+\*]+))?del([a-z]+)?ins([a-z]+)")

(defn parse-rna-indel
  [s]
  (let [[_ coord-s coord-e ref alt] (re-matches rna-indel-re s)]
    (rna-indel (coord/parse-rna-coordinate coord-s)
               (some-> coord-e coord/parse-rna-coordinate)
               ref
               alt)))

(defmethod restore "rna-indel"
  [m]
  (let [{:keys [coord-start coord-end ref alt]} (restore-coords m)]
    (rna-indel coord-start coord-end ref alt)))

;;; RNA - alleles
;;;
;;; e.g. r.[76a>u;103del]
;;;      r.[76a>u];[103del]
;;;      r.76[a>u];[a>u]
;;;      r.76a>u(;)103del (TODO)
;;;      r.76[a>u];[(a>u)] (TODO)
;;;      r.76[a>u];[=] (TODO)
;;;      r.[76a>u];[?] (TODO)
;;;      r.-124_-123[14];[18]

(defrecord RNAAlleles [mutations1 mutations2]
  Mutation
  (format [this] (format this nil))
  (format [this opts]
    (if (common-mutations? (concat mutations1 mutations2))
      (str (format-common (first mutations1) opts)
           "[" (string/join ";" (map #(format-unique % opts) mutations1)) "]"
           (if (seq mutations2)
             (str ";[" (string/join ";" (map #(format-unique % opts) mutations2)) "]")))
      (str "[" (string/join ";" (map #(format % opts) mutations1)) "]"
           (if (seq mutations2)
             (str ";[" (string/join ";" (map #(format % opts) mutations2)) "]")))))
  (plain [this]
    {:mutation "rna-alleles"
     :mutations1 (mapv plain mutations1)
     :mutations2 (mapv plain mutations2)}))

(defn rna-alleles
  [mutations1 mutations2]
  {:pre [(every? #(satisfies? Mutation %) mutations1)
         (every? #(satisfies? Mutation %) mutations2)]}
  (RNAAlleles. mutations1 mutations2))

(defn parse-rna-alleles
  [s]
  (let [[muts1 muts2] (split-alleles-mutations s)]
    (rna-alleles (mapv parse-rna (string/split muts1 #";"))
                 (if muts2
                   (mapv parse-rna (string/split muts2 #";"))))))

(defmethod restore "rna-alleles"
  [m]
  (rna-alleles (mapv restore (:mutations1 m))
               (mapv restore (:mutations2 m))))

;;; RNA - repeated sequences
;;;
;;; e.g. r.-124_-123[14]
;;;      r.-124ug[14]
;;;      r.-128_-126[(600_800)]

(defrecord RNARepeatedSeqs [coord-start coord-end ref ncopy]
  Mutation
  (format [this] (format this nil))
  (format [this opts]
    (str (format-common this opts) "[" (format-unique this opts) "]"))
  (plain [this]
    (into {:mutation "rna-repeated-seqs"} (plain-coords this)))
  SeparatelyFormat
  (format-common [this {:keys [range-format] :or {range-format :auto}}]
    (let [should-show-end? (neg? (compare coord-start coord-end))]
      (str (coord/format coord-start)
           (case range-format
             :auto (or ref (if should-show-end? (str "_" (coord/format coord-end))))
             :bases ref
             :coord (if should-show-end? (str "_" (coord/format coord-end)))))))
  (format-unique [this _]
    (format-ncopy ncopy)))

(defn rna-repeated-seqs
  "Constructor of RNARepeatedSeqs. Throws an exception if any input is illegal."
  [coord-start coord-end ref ncopy]
  {:pre [(satisfies? coord/Coordinate coord-start)
         (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))
         (or (nil? ref) (rna-bases? ref))
         (or (integer? ncopy) (vector? ncopy))]}
  (RNARepeatedSeqs. coord-start coord-end ref ncopy))

(def ^:private rna-repeated-seqs-re
  #"([\d\-\+\*]+)(?:_([\d\-\+\*]+))?([a-z]+)?\[([\d\(\)_]+)\]")

(defn parse-rna-repeated-seqs
  [s]
  (let [[_ coord-s coord-e ref ncopy] (re-matches rna-repeated-seqs-re s)]
    (rna-repeated-seqs (coord/parse-rna-coordinate coord-s)
                       (some-> coord-e coord/parse-rna-coordinate)
                       ref
                       (parse-ncopy ncopy))))

(defmethod restore "rna-repeated-seqs"
  [m]
  (let [{:keys [coord-start coord-end ref ncopy]} (restore-coords m)]
    (rna-repeated-seqs coord-start coord-end ref ncopy)))

(defn parse-rna
  "Parses a RNA mutation string s, returning a record implementing Mutation
  protocol."
  [s]
  ((condp re-find s
     #"\[.+;.+\]$" parse-rna-alleles
     #"delins" parse-rna-indel
     #"del" parse-rna-deletion
     #"dup" parse-rna-duplication
     #"ins" parse-rna-insertion
     #"inv" parse-rna-inversion
     #"con" parse-rna-conversion
     #"\[[\d\(\)_]+\]" parse-rna-repeated-seqs
     parse-rna-substitution)
   s))

;;; Protein mutations

(defn- amino-acid?
  [s]
  (or (some? ((set short-amino-acids) s))
      (some? ((set long-amino-acids) s))))

(defn- should-show-end?
  [ref-start coord-start ref-end coord-end]
  (and (some? ref-end)
       (neg? (compare coord-start coord-end))))

;;; Protein - substitution
;;;
;;; e.g. Arg54Ser
;;;      Trp26Ter
;;;      Trp26*
;;;      Cys123=

(defrecord ProteinSubstitution [ref coord alt]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    {:pre [(#{:long :short} amino-acid-format)]}
    (str (cond-> ref
           (= amino-acid-format :short) ->short-amino-acid)
         (coord/format coord)
         (if (= ref alt)
           "="
           (cond-> alt
             (= amino-acid-format :short) ->short-amino-acid))))
  (plain [this]
    (into {:mutation "protein-substitution"} (plain-coords this))))

(defn protein-substitution
  "Constructor of ProteinSubstitution. Throws an exception if any input is illegal."
  [ref coord alt]
  {:pre [(amino-acid? ref)
         (satisfies? coord/Coordinate coord)
         (amino-acid? alt)]}
  (ProteinSubstitution. ref coord alt))

(def ^:private protein-substitution-re
  #"([A-Z](?:[a-z]{2})?)(\d+)([A-Z\*=](?:[a-z]{2})?)")

(defn parse-protein-substitution
  [s]
  (let [[_ ref coord' alt] (re-matches protein-substitution-re s)]
    (protein-substitution (->long-amino-acid ref)
                          (coord/parse-protein-coordinate coord')
                          (case alt
                            "=" (->long-amino-acid ref)
                            "*" "Ter"
                            (->long-amino-acid alt)))))

(defmethod restore "protein-substitution"
  [m]
  (let [{:keys [ref coord alt]} (restore-coords m)]
    (protein-substitution ref coord alt)))

;;; Protein - deletion
;;;
;;; e.g. Ala3del
;;;      Cys76_Glu79del

(defrecord ProteinDeletion [ref-start coord-start ref-end coord-end]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (apply str (flatten [(cond-> ref-start
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-start)
                         (if (should-show-end? ref-start coord-start ref-end coord-end)
                           ["_"
                            (cond-> ref-end
                              (= amino-acid-format :short) ->short-amino-acid)
                            (coord/format coord-end)])
                         "del"])))
  (plain [this]
    (into {:mutation "protein-deletion"} (plain-coords this))))

(defn protein-deletion
  "Constructor of ProteinDeletion. Throws an exception if any input is illegal."
  ([ref-start coord-start] (protein-deletion ref-start coord-start nil nil))
  ([ref-start coord-start ref-end coord-end]
   {:pre [(amino-acid? ref-start)
          (satisfies? coord/Coordinate coord-start)
          (or (nil? ref-end) (amino-acid? ref-end))
          (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))]}
   (ProteinDeletion. ref-start coord-start ref-end coord-end)))

(def ^:private protein-deletion-re
  #"([A-Z](?:[a-z]{2})?)(\d+)(?:_([A-Z](?:[a-z]{2})?)(\d+))?del")

(defn parse-protein-deletion
  [s]
  (let [[_ ref-s coord-s ref-e coord-e] (re-matches protein-deletion-re s)]
    (protein-deletion (->long-amino-acid ref-s)
                      (coord/parse-protein-coordinate coord-s)
                      (->long-amino-acid ref-e)
                      (some-> coord-e coord/parse-protein-coordinate))))

(defmethod restore "protein-deletion"
  [m]
  (let [{:keys [ref-start coord-start ref-end coord-end]} (restore-coords m)]
    (protein-deletion ref-start coord-start ref-end coord-end)))

;;; Protein - duplication
;;;
;;; e.g. Ala3dup
;;;      Ala3_Ser5dup

(defrecord ProteinDuplication [ref-start coord-start ref-end coord-end]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (apply str (flatten [(cond-> ref-start
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-start)
                         (if (should-show-end? ref-start coord-start ref-end coord-end)
                           ["_"
                            (cond-> ref-end
                              (= amino-acid-format :short) ->short-amino-acid)
                            (coord/format coord-end)])
                         "dup"])))
  (plain [this]
    (into {:mutation "protein-duplication"} (plain-coords this))))

(defn protein-duplication
  "Constructor of ProteinDuplication. Throws an exception if any input is illegal."
  ([ref-start coord-start] (protein-duplication ref-start coord-start nil nil))
  ([ref-start coord-start ref-end coord-end]
   {:pre [(amino-acid? ref-start)
          (satisfies? coord/Coordinate coord-start)
          (or (nil? ref-end) (amino-acid? ref-end))
          (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))]}
   (ProteinDuplication. ref-start coord-start ref-end coord-end)))

(def ^:private protein-duplication-re
  #"([A-Z](?:[a-z]{2})?)(\d+)(?:_([A-Z](?:[a-z]{2})?)(\d+))?dup")

(defn parse-protein-duplication
  [s]
  (let [[_ ref-s coord-s ref-e coord-e] (re-matches protein-duplication-re s)]
    (protein-duplication (->long-amino-acid ref-s)
                         (coord/parse-protein-coordinate coord-s)
                         (->long-amino-acid ref-e)
                         (some-> coord-e coord/parse-protein-coordinate))))

(defmethod restore "protein-duplication"
  [m]
  (let [{:keys [ref-start coord-start ref-end coord-end]} (restore-coords m)]
    (protein-duplication ref-start coord-start ref-end coord-end)))

;;; Protein - insertion
;;;
;;; e.g. Lys23_Leu24insArgSerGln
;;;      Arg78_Gly79ins23

(defrecord ProteinInsertion [ref-start coord-start ref-end coord-end alts]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (apply str (flatten [(cond-> ref-start
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-start)
                         "_"
                         (cond-> ref-end
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-end)
                         "ins"
                         (if (every? #(= % "Xaa") alts)
                           (count alts)
                           (cond->> alts
                             (= amino-acid-format :short) (map ->short-amino-acid)))])))
  (plain [this]
    (into {:mutation "protein-insertion"} (plain-coords this))))

(defn protein-insertion
  "Constructor of ProteinInsertion. Throws an exception if any input is illegal."
  [ref-start coord-start ref-end coord-end alts]
  {:pre [(amino-acid? ref-start)
         (satisfies? coord/Coordinate coord-start)
         (amino-acid? ref-end)
         (satisfies? coord/Coordinate coord-end)
         (every? amino-acid? alts)]}
  (ProteinInsertion. ref-start coord-start ref-end coord-end alts))

(defn- parse-protein-insertion-alts
  [s]
  (condp re-matches s
    #"([A-Z]([a-z]{2})?)+" (mapv ->long-amino-acid (re-seq #"[A-Z](?:[a-z]{2})?" s))
    #"\d+" (vec (repeat (parse-long s) "Xaa"))))

(def ^:private protein-insertion-re
  #"([A-Z](?:[a-z]{2})?)(\d+)_([A-Z](?:[a-z]{2})?)(\d+)ins([\da-zA-Z]+)")

(defn parse-protein-insertion
  [s]
  (let [[_ ref-s coord-s ref-e coord-e alts] (re-matches protein-insertion-re s)]
    (protein-insertion (->long-amino-acid ref-s)
                       (coord/parse-protein-coordinate coord-s)
                       (->long-amino-acid ref-e)
                       (some-> coord-e coord/parse-protein-coordinate)
                       (parse-protein-insertion-alts alts))))

(defmethod restore "protein-insertion"
  [m]
  (let [{:keys [ref-start coord-start ref-end coord-end alts]} (restore-coords m)]
    (protein-insertion ref-start coord-start ref-end coord-end alts)))

;;; Protein - indel
;;;
;;; e.g. Cys28delinsTrpVal
;;;      Cys28_Lys29delinsTrp

(defrecord ProteinIndel [ref-start coord-start ref-end coord-end alts]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (apply str (flatten [(cond-> ref-start
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-start)
                         (if (should-show-end? ref-start coord-start ref-end coord-end)
                           ["_"
                            (cond-> ref-end
                              (= amino-acid-format :short) ->short-amino-acid)
                            (coord/format coord-end)])
                         "delins"
                         (cond->> alts
                           (= amino-acid-format :short) (map ->short-amino-acid))])))
  (plain [this]
    (into {:mutation "protein-indel"} (plain-coords this))))

(defn protein-indel
  "Constructor of ProteinIndel. Throws an exception if any input is illegal."
  [ref-start coord-start ref-end coord-end alts]
  {:pre [(amino-acid? ref-start)
         (satisfies? coord/Coordinate coord-start)
         (or (nil? ref-end) (amino-acid? ref-end))
         (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))
         (every? amino-acid? alts)]}
  (ProteinIndel. ref-start coord-start ref-end coord-end alts))

(def ^:private protein-indel-re
  #"([A-Z](?:[a-z]{2})?)(\d+)(?:_([A-Z](?:[a-z]{2})?)(\d+))?delins([A-Z][a-zA-Z]*)?")

(defn parse-protein-indel
  [s]
  (let [[_ ref-s coord-s ref-e coord-e alts] (re-matches protein-indel-re s)]
    (protein-indel (->long-amino-acid ref-s)
                   (coord/parse-protein-coordinate coord-s)
                   (->long-amino-acid ref-e)
                   (some-> coord-e coord/parse-protein-coordinate)
                   (mapv ->long-amino-acid (some->> alts (re-seq #"[A-Z](?:[a-z]{2})?"))))))

(defmethod restore "protein-indel"
  [m]
  (let [{:keys [ref-start coord-start ref-end coord-end alts]} (restore-coords m)]
    (protein-indel ref-start coord-start ref-end coord-end alts)))

;;; Protein - alleles
;;;
;;; e.g. p.[Ser73Arg;Asn603del]
;;;      p.[(Ser73Arg;Asn603del)] (TODO)
;;;      p.[Ser73Arg;(Asn603del)] (TODO)
;;;      p.[Ser73Arg];[Asn603del]
;;;      p.[(Ser73Arg)];[(Asn603del)] (TODO)
;;;      p.(Ser73Arg)(;)(Asn603del) (TODO)
;;;      p.[Ser73Arg];[Ser73=]
;;;      p.[Ser73Arg];[(?)] (TODO)
;;;      p.[Asn26His,Ala25_Gly29del] (TODO)
;;;      p.[Arg83=/Arg83Ser] (TODO)
;;;      p.[Arg83=//Arg83Ser] (TODO)
;;;      p.Ala2[10];[11]

(defrecord ProteinAlleles [mutations1 mutations2]
  Mutation
  (format [this] (format this nil))
  (format [this opts]
    (if (common-mutations? (concat mutations1 mutations2))
      (str (format-common (first mutations1) opts)
           "[" (string/join ";" (map #(format-unique % opts) mutations1)) "]"
           (if (seq mutations2)
             (str ";[" (string/join ";" (map #(format-unique % opts) mutations2)) "]")))
      (str "[" (string/join ";" (map #(format % opts) mutations1)) "]"
           (if (seq mutations2)
             (str ";[" (string/join ";" (map #(format % opts) mutations2)) "]")))))
  (plain [this]
    {:mutation "protein-alleles"
     :mutations1 (mapv plain mutations1)
     :mutations2 (mapv plain mutations2)}))

(defn protein-alleles
  [mutations1 mutations2]
  {:pre [(every? #(satisfies? Mutation %) mutations1)
         (every? #(satisfies? Mutation %) mutations2)]}
  (ProteinAlleles. mutations1 mutations2))

(defn parse-protein-alleles
  [s]
  (let [[mut1 mut2] (split-alleles-mutations s)]
    (protein-alleles (mapv parse-protein (string/split mut1 #";"))
                     (if mut2
                       (mapv parse-protein (string/split mut2 #";"))))))

(defmethod restore "protein-alleles"
  [m]
  (protein-alleles (mapv restore (:mutations1 m))
                   (mapv restore (:mutations2 m))))

;;; Protein - repeated sequences
;;;
;;; e.g. Ala2[10]
;;;      Arg65_Ser67[12]
;;;      (Gln18)[(70_80)]

(defrecord ProteinRepeatedSeqs [ref-start coord-start ref-end coord-end ncopy]
  Mutation
  (format [this] (format this nil))
  (format [this opts]
    (str (format-common this opts) "[" (format-unique this opts) "]"))
  (plain [this]
    (into {:mutation "protein-repeated-seqs"} (plain-coords this)))
  SeparatelyFormat
  (format-common [this {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (str (cond-> ref-start
           (= amino-acid-format :short) ->short-amino-acid)
         (coord/format coord-start)
         (if (should-show-end? ref-start coord-start ref-end coord-end)
           (str "_"
                (cond-> ref-end
                  (= amino-acid-format :short) ->short-amino-acid)
                (coord/format coord-end)))))
  (format-unique [this _]
    (format-ncopy ncopy)))

(defn protein-repeated-seqs
  "Constructor of ProteinRepeatedSeqs. Throws an exception if any input is illegal."
  [ref-start coord-start ref-end coord-end ncopy]
  {:pre [(amino-acid? ref-start)
         (satisfies? coord/Coordinate coord-start)
         (or (nil? ref-end) (amino-acid? ref-end))
         (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))
         (or (integer? ncopy) (vector? ncopy))]}
  (ProteinRepeatedSeqs. ref-start coord-start ref-end coord-end ncopy))

(def ^:private protein-repeated-seqs-re
  #"([A-Z](?:[a-z]{2})?)(\d+)(?:_([A-Z](?:[a-z]{2})?)(\d+))?\[([\d\(\)_]+)\]")

(defn parse-protein-repeated-seqs
  [s]
  (let [[_ ref-s coord-s ref-e coord-e ncopy] (re-matches protein-repeated-seqs-re s)]
    (protein-repeated-seqs (->long-amino-acid ref-s)
                           (coord/parse-protein-coordinate coord-s)
                           (->long-amino-acid ref-e)
                           (some-> coord-e coord/parse-protein-coordinate)
                           (parse-ncopy ncopy))))

(defmethod restore "protein-repeated-seqs"
  [m]
  (let [{:keys [ref-start coord-start ref-end coord-end ncopy]} (restore-coords m)]
    (protein-repeated-seqs ref-start coord-start ref-end coord-end ncopy)))

;;; Protein - frame shift
;;;
;;; e.g. Arg97ProfsTer23
;;;      Arg97fs
;;;      Pro661=fs
;;;      Ile327Argfs*?
;;;      Gln151Thrfs*9

(defrecord ProteinFrameShift [ref coord alt new-ter-site]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [amino-acid-format ter-format]
                 :or {amino-acid-format :long}}]
    (str (cond-> ref
           (= amino-acid-format :short) ->short-amino-acid)
         (coord/format coord)
         (if (= ref alt)
           "="
           (cond-> alt
             (= amino-acid-format :short) ->short-amino-acid))
         "fs"
         (if (some? new-ter-site)
           (str (cond
                  (= ter-format :short) (->short-amino-acid "Ter")
                  (= ter-format :long) "Ter"
                  (= amino-acid-format :short) (->short-amino-acid "Ter")
                  :else "Ter")
                (coord/format new-ter-site)))))
  (plain [this]
    (into {:mutation "protein-frame-shift"} (plain-coords this))))

(defn protein-frame-shift
  "Constructor of ProteinFrameShift. Throws an exception if any input is illegal."
  [ref coord alt new-ter-site]
  {:pre [(amino-acid? ref)
         (satisfies? coord/Coordinate coord)
         (or (nil? alt) (amino-acid? alt))
         (or (nil? new-ter-site) (satisfies? coord/Coordinate new-ter-site))]}
  (ProteinFrameShift. ref coord alt new-ter-site))

(def ^:private protein-frame-shift-re
  #"([A-Z](?:[a-z]{2})?)(\d+)([A-Z=](?:[a-z]{2})?)?fs(?:Ter|\*)?(\?|\d+)?")

(defn parse-protein-frame-shift
  [s]
  (let [[_ ref coord' alt new-ter-site] (re-matches protein-frame-shift-re s)]
    (protein-frame-shift (->long-amino-acid ref)
                         (coord/parse-protein-coordinate coord')
                         (->long-amino-acid (if (= alt "=") ref alt))
                         (some-> new-ter-site coord/parse-protein-coordinate))))

(defmethod restore "protein-frame-shift"
  [m]
  (let [{:keys [ref coord alt new-ter-site]} (restore-coords m)]
    (protein-frame-shift ref coord alt new-ter-site)))

;;; Protein - extension
;;;
;;; e.g. Met1ext-5
;;;      Met1Valext-12
;;;      Ter110Glnext*17
;;;      *110Glnext*17
;;;      Ter327Argext*?

(defrecord ProteinExtension [ref coord alt region new-site]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [amino-acid-format ter-format]
                 :or {amino-acid-format :long}}]
    (str (case ref
           "Met" (cond-> ref (= amino-acid-format :short) ->short-amino-acid)
           "Ter" (cond
                   (= ter-format :short) (->short-amino-acid ref)
                   (= ter-format :long) ref
                   (= amino-acid-format :short) (->short-amino-acid ref)
                   :else ref))
         (coord/format coord)
         (cond-> alt
           (= amino-acid-format :short) ->short-amino-acid)
         "ext"
         (coord/->region-str region)
         (coord/format new-site)))
  (plain [this]
    (into {:mutation "protein-extension"} (plain-coords this))))

(defn protein-extension
  "Constructor of ProteinExtension. Throws an exception if any input is illegal."
  [ref coord alt region new-site]
  {:pre [(amino-acid? ref)
         (satisfies? coord/Coordinate coord)
         (or (nil? alt) (amino-acid? alt))
         (#{:upstream :downstream} region)
         (satisfies? coord/Coordinate coord)]}
  (ProteinExtension. ref coord alt region new-site))

(def ^:private protein-extension-re
  #"([A-Z\*](?:[a-z]{2})?)(\d+)([A-Z](?:[a-z]{2})?)?ext(\-|\*)(\?|\d+)")

(defn parse-protein-extension
  [s]
  (let [[_ ref coord' alt region new-site] (re-matches protein-extension-re s)]
    (protein-extension (->long-amino-acid ref)
                       (coord/parse-protein-coordinate coord')
                       (->long-amino-acid alt)
                       (coord/->region-keyword region)
                       (coord/parse-protein-coordinate new-site))))

(defmethod restore "protein-extension"
  [m]
  (let [{:keys [ref coord alt region new-site]} (restore-coords m)]
    (protein-extension ref coord alt region new-site)))

(defn parse-protein
  "Parses a protein mutation string s, returning a record implementing Mutation
  protocol."
  [s]
  ((condp re-find s
     #"\[.+;.+\]$" parse-protein-alleles
     #"delins" parse-protein-indel
     #"del" parse-protein-deletion
     #"dup" parse-protein-duplication
     #"ins" parse-protein-insertion
     #"fs" parse-protein-frame-shift
     #"ext" parse-protein-extension
     #"\[[\d\(\)_]+\]" parse-protein-repeated-seqs
     parse-protein-substitution)
   s))

(defn- common-mutations?
  [mutations]
  (and (apply = (map type mutations))
       (condp instance? (first mutations)
         DNASubstitution (apply = (map :coord mutations))
         DNARepeatedSeqs (apply = (map #(dissoc % :ncopy) mutations))
         RNASubstitution (apply = (map :coord mutations))
         RNARepeatedSeqs (apply = (map #(dissoc % :ncopy) mutations))
         ProteinRepeatedSeqs (apply = (map #(dissoc % :ncopy) mutations))
         false)))
