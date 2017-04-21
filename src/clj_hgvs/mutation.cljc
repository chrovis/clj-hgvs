(ns clj-hgvs.mutation
  "Data structures and functions to handle HGVS mutations."
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.string :as string]
            [clj-hgvs.coordinate :as coord]
            [clj-hgvs.internal :refer [parse-long ->kind-keyword ->kind-str]]))

(def short-amino-acids
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
   "*"])

(def long-amino-acids
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
   "Ter"])

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

(defprotocol Mutation
  (format [this] [this opts]))

;;; DNA mutations

(defn- dna-bases?
  [s]
  (and (string? s) (some? (re-matches #"[ACGT]+" s))))

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
    (apply str (flatten [(coord/format coord)
                         ref
                         type
                         alt]))))

(defn dna-substitution
  ([coord ref typ] (dna-substitution coord ref typ nil))
  ([coord ref typ alt]
   {:pre [(satisfies? coord/Coordinate coord)
          (dna-bases? ref)
          (#{">" "=" "=/>" "=//>"} typ)
          (or (nil? alt) (dna-bases? alt))]}
   (DNASubstitution. coord ref typ alt)))

(def ^:private dna-substitution-re
  #"([\d\-\+]+)([A-Z])([>=/]+)([A-Z])?")

(defn parse-dna-substitution
  [s kind]
  (let [[_ coord ref typ alt] (re-matches dna-substitution-re s)
        parse-coord (coord-parser kind)]
    (dna-substitution (parse-coord coord) ref typ alt)))

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
                         (if show-bases? ref)]))))

(defn dna-deletion
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
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))?del([A-Z]+)?")

(defn- parse-dna-deletion*
  [s kind]
  (let [[_ coord-s coord-e ref] (re-matches dna-deletion-re s)
        parse-coord (coord-parser kind)]
    (dna-deletion (parse-coord coord-s)
                  (some-> coord-e parse-coord)
                  ref)))

(def ^:private dna-deletion-range-re
  #"\(([\d\-\+\?\*]+)_([\d\-\+\?\*]+)\)_\(([\d\-\+\?\*]+)_([\d\-\+\?\*]+)\)del([A-Z]+)?")

(defn- parse-dna-deletion-range
  [s kind]
  (let [[_ coord-s1 coord-s2 coord-e1 coord-e2 ref] (re-matches dna-deletion-range-re s)
        parse-coord (coord-parser kind)]
    (dna-deletion (mapv parse-coord [coord-s1 coord-s2])
                  (mapv parse-coord [coord-e1 coord-e2])
                  ref)))

(defn parse-dna-deletion
  [s kind]
  (if (re-matches #"\(.+\)_\(.+\)del[A-Z]*" s)
    (parse-dna-deletion-range s kind)
    (parse-dna-deletion* s kind)))

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
                         (if show-bases? ref)]))))

(defn dna-duplication
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
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))?dup([A-Z]+)?")

(defn- parse-dna-duplication*
  [s kind]
  (let [[_ coord-s coord-e ref] (re-matches dna-duplication-re s)
        parse-coord (coord-parser kind)]
    (dna-duplication (parse-coord coord-s)
                     (some-> coord-e parse-coord)
                     ref)))

(def ^:private dna-duplication-range-re
  #"\(([\d\-\+\?\*]+)_([\d\-\+\?\*]+)\)_\(([\d\-\+\?\*]+)_([\d\-\+\?\*]+)\)dup([A-Z]+)?")

(defn- parse-dna-duplication-range
  [s kind]
  (let [[_ coord-s1 coord-s2 coord-e1 coord-e2 ref] (re-matches dna-duplication-range-re s)
        parse-coord (coord-parser kind)]
    (dna-duplication (mapv parse-coord [coord-s1 coord-s2])
                     (mapv parse-coord [coord-e1 coord-e2])
                     ref)))

(defn parse-dna-duplication
  [s kind]
  (if (re-matches #"\(.+\)_\(.+\)dup[A-Z]*" s)
    (parse-dna-duplication-range s kind)
    (parse-dna-duplication* s kind)))

;;; DNA - insertion
;;;
;;; e.g. g.5756_5757insAGG
;;;      g.123_124insL37425.1:23_361
;;;      g.122_123ins123_234inv (TODO)
;;;      g.122_123ins213_234invinsAins123_211inv (TODO)

(defrecord DNAInsertion [coord-start coord-end alt]
  Mutation
  (format [this] (format this nil))
  (format [this _]
    (apply str (flatten [(coord/format coord-start)
                         "_"
                         (coord/format coord-end)
                         "ins"
                         (cond
                           (string? alt) alt
                           (map? alt) [(:transcript alt)
                                       ":"
                                       (coord/format (:coord-start alt))
                                       "_"
                                       (coord/format (:coord-end alt))])]))))

(defn dna-insertion
  [coord-start coord-end alt]
  {:pre [(satisfies? coord/Coordinate coord-start)
         (satisfies? coord/Coordinate coord-end)
         (or (dna-bases? alt) (map? alt))]}
  (DNAInsertion. coord-start coord-end alt))

(def ^:private dna-insertion-re
  #"([\d\-\+\*\?]+)_([\d\-\+\*\?]+)ins(.+)")

(def ^:private dna-insertion-alt-re
  #"(?:([^:]+):)([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))")

(defn parse-dna-insertion
  [s kind]
  (let [[_ coord-s coord-e alt] (re-matches dna-insertion-re s)
        parse-coord (coord-parser kind)]
    (dna-insertion (parse-coord coord-s)
                   (parse-coord coord-e)
                   (or (re-matches #"[A-Z]+" alt)
                       (let [[_ transcript coord-s coord-e] (re-matches dna-insertion-alt-re alt)]
                         {:transcript transcript
                          :coord-start (parse-coord coord-s)
                          :coord-end (parse-coord coord-e)})))))

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
         "inv")))

(defn dna-inversion
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
                         (coord/format (:coord-end alt))]))))

(defn dna-conversion
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
                         alt]))))

(defn dna-indel
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

;;; DNA - repeated sequences
;;;
;;; e.g. g.123_124[14]
;;;      g.123TG[14]
;;;      g.123_124[14];[18]

(defrecord DNARepeatedSeqs [coord-start coord-end ref ncopy ncopy-other]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [range-format] :or {range-format :auto}}]
    (let [should-show-end? (neg? (compare coord-start coord-end))]
      (str (coord/format coord-start)
           (case range-format
             :auto (or ref (if should-show-end? (str "_" (coord/format coord-end))))
             :bases ref
             :coord (if should-show-end? (str "_" (coord/format coord-end))))
           "[" ncopy "]"
           (if ncopy-other (str ";[" ncopy-other "]"))))))

(defn dna-repeated-seqs
  ([coord-start coord-end ref ncopy]
   (dna-repeated-seqs coord-start coord-end ref ncopy nil))
  ([coord-start coord-end ref ncopy ncopy-other]
   {:pre [(satisfies? coord/Coordinate coord-start)
          (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))
          (or (nil? ref) (dna-bases? ref))
          (integer? ncopy)
          (or (nil? ncopy-other) (integer? ncopy-other))]}
   (DNARepeatedSeqs. coord-start coord-end ref ncopy ncopy-other)))

(def ^:private dna-repeated-seqs-re
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))?([A-Z]+)?\[(\d+)\](?:;\[(\d+)\])?")

(defn parse-dna-repeated-seqs
  [s kind]
  (let [[_ coord-s coord-e ref ncopy1 ncopy2] (re-matches dna-repeated-seqs-re s)
        parse-coord (coord-parser kind)]
    (dna-repeated-seqs (parse-coord coord-s)
                       (some-> coord-e parse-coord)
                       ref
                       (parse-long ncopy1)
                       (some-> ncopy2 parse-long))))

(defn parse-dna
  [s kind]
  ((condp re-find s
     #"delins" parse-dna-indel
     #"del" parse-dna-deletion
     #"dup" parse-dna-duplication
     #"ins" parse-dna-insertion
     #"inv" parse-dna-inversion
     #"con" parse-dna-conversion
     #"\[\d+\]" parse-dna-repeated-seqs
     parse-dna-substitution)
   s kind))

;;; RNA mutations

(defn- rna-bases?
  [s]
  (and (string? s) (some? (re-matches #"[acgu]+" s))))

;;; RNA - substitution
;;;
;;; e.g. r.76a>c
;;;      r.-14g>c
;;;      r.*46u>a

(defrecord RNASubstitution [coord ref alt]
  Mutation
  (format [this] (format this nil))
  (format [this _]
    (apply str (coord/format coord) ref ">" alt)))

(defn rna-substitution
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
         (if show-bases? ref))))

(defn rna-deletion
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
         (if show-bases? ref))))

(defn rna-duplication
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
           :else alt))))

(defn rna-insertion
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
         "inv")))

(defn rna-inversion
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
         (coord/format (:coord-end alt)))))

(defn rna-conversion
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
         alt)))

(defn rna-indel
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

;;; RNA - repeated sequences
;;;
;;; e.g. r.-124_-123[14]
;;;      r.-124ug[14]
;;;      r.-124_-123[14];[18]

(defrecord RNARepeatedSeqs [coord-start coord-end ref ncopy ncopy-other]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [range-format] :or {range-format :auto}}]
    (let [should-show-end? (neg? (compare coord-start coord-end))]
      (str (coord/format coord-start)
           (case range-format
             :auto (or ref (if should-show-end? (str "_" (coord/format coord-end))))
             :bases ref
             :coord (if should-show-end? (str "_" (coord/format coord-end))))
           "[" ncopy "]"
           (if ncopy-other (str ";[" ncopy-other "]"))))))

(defn rna-repeated-seqs
  ([coord-start coord-end ref ncopy]
   (rna-repeated-seqs coord-start coord-end ref ncopy nil))
  ([coord-start coord-end ref ncopy ncopy-other]
   {:pre [(satisfies? coord/Coordinate coord-start)
          (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))
          (or (nil? ref) (rna-bases? ref))
          (integer? ncopy)
          (or (nil? ncopy-other) (integer? ncopy-other))]}
   (RNARepeatedSeqs. coord-start coord-end ref ncopy ncopy-other)))

(def ^:private rna-repeated-seqs-re
  #"([\d\-\+\*]+)(?:_([\d\-\+\*]+))?([a-z]+)?\[(\d+)\](?:;\[(\d+)\])?")

(defn parse-rna-repeated-seqs
  [s]
  (let [[_ coord-s coord-e ref ncopy1 ncopy2] (re-matches rna-repeated-seqs-re s)]
    (rna-repeated-seqs (coord/parse-rna-coordinate coord-s)
                       (some-> coord-e coord/parse-rna-coordinate)
                       ref
                       (parse-long ncopy1)
                       (some-> ncopy2 parse-long))))

(defn parse-rna
  [s]
  ((condp re-find s
     #"delins" parse-rna-indel
     #"del" parse-rna-deletion
     #"dup" parse-rna-duplication
     #"ins" parse-rna-insertion
     #"inv" parse-rna-inversion
     #"con" parse-rna-conversion
     #"\[\d+\]" parse-rna-repeated-seqs
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
             (= amino-acid-format :short) ->short-amino-acid)))))

(defn protein-substitution
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
                         "del"]))))

(defn protein-deletion
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
                         "dup"]))))

(defn protein-duplication
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

;;; Protein - insertion
;;;
;;; e.g. Lys23_Leu24insArgSerGln

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
                         (cond->> alts
                           (= amino-acid-format :short) (map ->short-amino-acid))]))))

(defn protein-insertion
  [ref-start coord-start ref-end coord-end alts]
  {:pre [(amino-acid? ref-start)
         (satisfies? coord/Coordinate coord-start)
         (amino-acid? ref-end)
         (satisfies? coord/Coordinate coord-end)
         (every? amino-acid? alts)]}
  (ProteinInsertion. ref-start coord-start ref-end coord-end alts))

(def ^:private protein-insertion-re
  #"([A-Z](?:[a-z]{2})?)(\d+)_([A-Z](?:[a-z]{2})?)(\d+)ins([A-Z][a-zA-Z]*)?")

(defn parse-protein-insertion
  [s]
  (let [[_ ref-s coord-s ref-e coord-e alts] (re-matches protein-insertion-re s)]
    (protein-insertion (->long-amino-acid ref-s)
                       (coord/parse-protein-coordinate coord-s)
                       (->long-amino-acid ref-e)
                       (some-> coord-e coord/parse-protein-coordinate)
                       (mapv ->long-amino-acid (some->> alts (re-seq #"[A-Z](?:[a-z]{2})?"))))))

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
                           (= amino-acid-format :short) (map ->short-amino-acid))]))))

(defn protein-indel
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

;;; Protein - repeated sequences
;;;
;;; e.g. Ala2[10]
;;;      Ala2[10];[11]
;;;      Arg65_Ser67[12]

(defrecord ProteinRepeatedSeqs [ref-start coord-start ref-end coord-end ncopy
                                ncopy-other]
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
                         "[" ncopy "]"
                         (if ncopy-other [";[" ncopy-other "]"])]))))

(defn protein-repeated-seqs
  ([ref-start coord-start ref-end coord-end ncopy]
   (protein-repeated-seqs ref-start coord-start ref-end coord-end ncopy nil))
  ([ref-start coord-start ref-end coord-end ncopy ncopy-other]
   {:pre [(amino-acid? ref-start)
          (satisfies? coord/Coordinate coord-start)
          (or (nil? ref-end) (amino-acid? ref-end))
          (or (nil? coord-end) (satisfies? coord/Coordinate coord-end))
          (integer? ncopy)
          (or (nil? ncopy-other) (integer? ncopy-other))]}
   (ProteinRepeatedSeqs. ref-start coord-start ref-end coord-end ncopy ncopy-other)))

(def ^:private protein-repeated-seqs-re
  #"([A-Z](?:[a-z]{2})?)(\d+)(?:_([A-Z](?:[a-z]{2})?)(\d+))?\[(\d+)\](?:;\[(\d+)\])?")

(defn parse-protein-repeated-seqs
  [s]
  (let [[_ ref-s coord-s ref-e coord-e ncopy1 ncopy2] (re-matches protein-repeated-seqs-re s)]
    (protein-repeated-seqs (->long-amino-acid ref-s)
                           (coord/parse-protein-coordinate coord-s)
                           (->long-amino-acid ref-e)
                           (some-> coord-e coord/parse-protein-coordinate)
                           (parse-long ncopy1)
                           (if ncopy2 (parse-long ncopy2)))))

;;; Protein - frame shift
;;;
;;; e.g. Arg97ProfsTer23
;;;      Arg97fs
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
         (cond-> alt
           (= amino-acid-format :short) ->short-amino-acid)
         "fs"
         (if (some? new-ter-site)
           (str (cond
                  (= ter-format :short) (->short-amino-acid "Ter")
                  (= ter-format :long) "Ter"
                  (= amino-acid-format :short) (->short-amino-acid "Ter")
                  :else "Ter")
                (coord/format new-ter-site))))))

(defn protein-frame-shift
  [ref coord alt new-ter-site]
  {:pre [(amino-acid? ref)
         (satisfies? coord/Coordinate coord)
         (or (nil? alt) (amino-acid? alt))
         (or (nil? new-ter-site) (satisfies? coord/Coordinate new-ter-site))]}
  (ProteinFrameShift. ref coord alt new-ter-site))

(def ^:private protein-frame-shift-re
  #"([A-Z](?:[a-z]{2})?)(\d+)([A-Z](?:[a-z]{2})?)?fs(?:Ter|\*)?(\?|\d+)?")

(defn parse-protein-frame-shift
  [s]
  (let [[_ ref coord' alt new-ter-site] (re-matches protein-frame-shift-re s)]
    (protein-frame-shift (->long-amino-acid ref)
                         (coord/parse-protein-coordinate coord')
                         (->long-amino-acid alt)
                         (some-> new-ter-site coord/parse-protein-coordinate))))

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
         (coord/format new-site))))

(defn protein-extension
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

(defn parse-protein
  [s]
  ((condp re-find s
     #"delins" parse-protein-indel
     #"del" parse-protein-deletion
     #"dup" parse-protein-duplication
     #"ins" parse-protein-insertion
     #"fs" parse-protein-frame-shift
     #"ext" parse-protein-extension
     #"\[\d+\]" parse-protein-repeated-seqs
     parse-protein-substitution)
   s))
