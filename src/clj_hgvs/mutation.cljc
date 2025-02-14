(ns clj-hgvs.mutation
  "Data structures and functions to handle HGVS mutations."
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.spec.alpha :as s]
            [clojure.string :as string]
            [clojure.walk :as walk]
            [clj-hgvs.coordinate :as coord]
            [clj-hgvs.internal :as intl]))

(declare parse-dna parse-rna parse-protein parse common-mutations?)

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

(def ^:private short-amino-acid-set (set short-amino-acids))

(def ^:private long-amino-acid-set (set long-amino-acids))

(def ^:private short-long-amino-acid-map
  (zipmap short-amino-acids long-amino-acids))

(def ^:private long-short-amino-acid-map
  (zipmap long-amino-acids short-amino-acids))

(defn ->long-amino-acid
  "Converts a single-letter amino acid into a three-letter one. s must be String
  or Character. Returns nil if s is not present in amino acid list."
  [s]
  (let [s (cond-> s (char? s) str)]
    (or (long-amino-acid-set s) (short-long-amino-acid-map s))))

(defn ->short-amino-acid
  "Converts a three-letter amino acid into a single-letter one. s must be String
  or Character. Returns nil if s is not present in amino acid list."
  [s]
  (let [s (cond-> s (char? s) str)]
    (or (short-amino-acid-set s) (long-short-amino-acid-map s))))

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
    [(intl/parse-long s1) (intl/parse-long s2)]
    (intl/parse-long s)))

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

(s/def ::mutation #(satisfies? Mutation %))

(s/def :clj-hgvs.mutation.plain-mutation/mutation string?)
(s/def ::plain-mutation
  (s/keys :req-un [:clj-hgvs.mutation.plain-mutation/mutation]))

(defmulti restore
  "Restores a plain map to a suitable mutation record."
  {:arglists '([m])}
  :mutation)

;; SeparatelyFormat protocol is used for formatting alleles mutation, such as
;; c.2376[G>C];[G>C], including common part.
(defprotocol SeparatelyFormat
  (format-common [this opts])
  (format-unique [this opts]))

(defprotocol Equivalence
  (equiv* [this o]))

(defn equiv
  [mutation1 mutation2]
  (if (satisfies? Equivalence mutation1)
    (equiv* mutation1 mutation2)
    (= mutation1 mutation2)))

;;; Uncertain mutation
;;;
;;; e.g. r.(?)
;;;      r.(306g>u)
;;;      p.(Arg2371Ser)

(defrecord UncertainMutation [mutation]
  Mutation
  (format [this] (format this nil))
  (format [this opts]
    (str "(" (format mutation opts) ")"))
  (plain [this]
    {:mutation "uncertain-mutation"
     :content-mutation (plain mutation)})

  Equivalence
  (equiv* [this o]
    (if (instance? UncertainMutation o)
      (equiv mutation (:mutation o))
      false)))

(s/def ::uncertain-mutation
  (s/and ::mutation (s/keys :req-un [::mutation])))

(defn uncertain-mutation
  "Constructor of UncertainMutation. Throws an exception if any input is
  illegal."
  [mutation]
  {:pre [(not (instance? UncertainMutation mutation))]
   :post [(intl/valid? ::uncertain-mutation %)]}
  (UncertainMutation. mutation))

(def ^:private uncertain-mutation-re
  #"\((\S+)\)")

(defn parse-uncertain-mutation
  [s kind]
  (let [[_ mut] (re-matches uncertain-mutation-re s)]
    (uncertain-mutation (parse mut kind))))

(defmethod restore "uncertain-mutation"
  [m]
  (uncertain-mutation (restore (:content-mutation m))))

;;; DNA mutations

;; See https://hgvs-nomenclature.org/stable/background/standards/#dna
(s/def ::dna-bases
  (s/and string? #(re-matches #"[ACGTBDHKMNRSVWY]+" %)))

(defn- coord-parser
  [kind]
  (case kind
    :genome coord/parse-genomic-coordinate
    :mitochondria coord/parse-mitochondrial-coordinate
    :coding-dna coord/parse-coding-dna-coordinate
    :non-coding-dna coord/parse-non-coding-dna-coordinate
    :circular-dna coord/parse-circular-dna-coordinate))

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
  (format-unique [this _]
    (str ref type (when-not (= type "=") alt))))

(s/def :clj-hgvs.mutation.dna-substitution/coord ::coord/coordinate)
(s/def :clj-hgvs.mutation.dna-substitution/ref ::dna-bases)
(s/def :clj-hgvs.mutation.dna-substitution/type #{">" "=" "=/>" "=//>"})
(s/def :clj-hgvs.mutation.dna-substitution/alt (s/nilable ::dna-bases))
(s/def ::dna-substitution
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.dna-substitution/coord
                          :clj-hgvs.mutation.dna-substitution/ref
                          :clj-hgvs.mutation.dna-substitution/type
                          :clj-hgvs.mutation.dna-substitution/alt])))

(defn dna-substitution
  "Constructor of DNASubstitution. Throws an exception if any input is illegal."
  ([coord ref typ] (dna-substitution coord ref typ nil))
  ([coord ref typ alt]
   {:post [(intl/valid? ::dna-substitution %)]}
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
    (apply str (flatten [(coord/format coord-start)
                         (if coord-end "_")
                         (some-> coord-end coord/format)
                         "del"
                         (if show-bases? ref)])))
  (plain [this]
    (into {:mutation "dna-deletion"} (plain-coords this)))

  Equivalence
  (equiv* [this o]
    (if (instance? DNADeletion o)
      (and (= coord-start (:coord-start o))
           (= coord-end (:coord-end o))
           (if (and ref (:ref o))
             (= ref (:ref o))
             true))
      false)))

(s/def :clj-hgvs.mutation.dna-deletion/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.dna-deletion/coord-end (s/nilable ::coord/coordinate))
(s/def :clj-hgvs.mutation.dna-deletion/ref (s/nilable ::dna-bases))
(s/def ::dna-deletion
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.dna-deletion/coord-start
                          :clj-hgvs.mutation.dna-deletion/coord-end
                          :clj-hgvs.mutation.dna-deletion/ref])))

(defn dna-deletion
  "Constructor of DNADeletion. Throws an exception if any input is illegal."
  ([coord-start coord-end] (dna-deletion coord-start coord-end nil))
  ([coord-start coord-end ref]
   {:pre [(or (nil? coord-end)
              (not (coord/comparable-coordinates? coord-start coord-end))
              (neg? (compare coord-start coord-end)))]
    :post [(intl/valid? ::dna-deletion %)]}
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
    (apply str (flatten [(coord/format coord-start)
                         (if (and (some? coord-end) (not= coord-start coord-end))
                           ["_"
                            (coord/format coord-end)])
                         "dup"
                         (if show-bases? ref)])))
  (plain [this]
    (into {:mutation "dna-duplication"} (plain-coords this)))

  Equivalence
  (equiv* [this o]
    (if (instance? DNADuplication o)
      (and (= coord-start (:coord-start o))
           (= coord-end (:coord-end o))
           (if (and ref (:ref o))
             (= ref (:ref o))
             true))
      false)))

(s/def :clj-hgvs.mutation.dna-duplication/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.dna-duplication/coord-end (s/nilable ::coord/coordinate))
(s/def :clj-hgvs.mutation.dna-duplication/ref (s/nilable ::dna-bases))
(s/def ::dna-duplication
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.dna-duplication/coord-start
                          :clj-hgvs.mutation.dna-duplication/coord-end
                          :clj-hgvs.mutation.dna-duplication/ref])))

(defn dna-duplication
  "Constructor of DNADuplication. Throws an exception if any input is illegal."
  ([coord-start coord-end] (dna-duplication coord-start coord-end nil))
  ([coord-start coord-end ref]
   {:pre [(or (nil? coord-end)
              (not (coord/comparable-coordinates? coord-start coord-end))
              (neg? (compare coord-start coord-end)))]
    :post [(intl/valid? ::dna-duplication %)]}
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
;;;      g.1134_1135insN[100]
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
                                                   (str "N[" (count alt) "]")
                                                   alt)
                                           :bases alt
                                           :count (str "N[" (count alt) "]"))
                           (map? alt) [(:transcript alt)
                                       ":"
                                       (coord/format (:coord-start alt))
                                       "_"
                                       (coord/format (:coord-end alt))])])))
  (plain [this]
    (into {:mutation "dna-insertion"} (plain-coords this))))

(s/def :clj-hgvs.mutation.dna-insertion/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.dna-insertion/coord-end ::coord/coordinate)
(s/def :clj-hgvs.mutation.dna-insertion/alt (s/or :dna-bases ::dna-bases
                                                  :ref-seq map?))
(s/def ::dna-insertion
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.dna-insertion/coord-start
                          :clj-hgvs.mutation.dna-insertion/coord-end
                          :clj-hgvs.mutation.dna-insertion/alt])))

(defn dna-insertion
  "Constructor of DNAInsertion. Throws an exception if any input is illegal."
  [coord-start coord-end alt]
  {:pre [(or (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))]
   :post [(intl/valid? ::dna-insertion %)]}
  (DNAInsertion. coord-start coord-end alt))

(defn- parse-dna-insertion-alt
  [s kind]
  (or (re-matches #"[A-Z]+" s)
      (some-> (re-matches #"N\[(\d+)\]" s)
              (second)
              (intl/parse-long)
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

(s/def :clj-hgvs.mutation.dna-inversion/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.dna-inversion/coord-end ::coord/coordinate)
(s/def ::dna-inversion
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.dna-inversion/coord-start
                          :clj-hgvs.mutation.dna-inversion/coord-end])))

(defn dna-inversion
  "Constructor of DNAInversion. Throws an exception if any input is illegal."
  [coord-start coord-end]
  {:pre [(or (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))]
   :post [(intl/valid? ::dna-inversion %)]}
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
                         (if (:kind alt) [(intl/->kind-str (:kind alt)) "."])
                         (coord/format (:coord-start alt))
                         "_"
                         (coord/format (:coord-end alt))])))
  (plain [this]
    (into {:mutation "dna-conversion"} (plain-coords this))))

(s/def :clj-hgvs.mutation.dna-conversion/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.dna-conversion/coord-end ::coord/coordinate)
(s/def ::dna-conversion
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.dna-conversion/coord-start
                          :clj-hgvs.mutation.dna-conversion/coord-end])))

(defn dna-conversion
  "Constructor of DNAConversion. Throws an exception if any input is illegal."
  [coord-start coord-end alt]
  {:pre [(or (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))
         (map? alt)]
   :post [(intl/valid? ::dna-conversion %)]}
  (DNAConversion. coord-start coord-end alt))

(def ^:private dna-conversion-re
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))con(.+)")

(def ^:private dna-conversion-alt-re
  #"(?:([^:]+):)?(?:([gmcn])\.)?([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))")

(defn- parse-dna-conversion-alt
  [s kind default-coord-parser]
  (let [[_ transcript kind coord-s coord-e] (re-matches dna-conversion-alt-re s)
        parse-alt-coord (if kind
                          (coord-parser (intl/->kind-keyword kind))
                          default-coord-parser)]
    {:transcript transcript
     :kind (some-> kind intl/->kind-keyword)
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
;;;      c.145_147delinsN[10]

(defrecord DNAIndel [coord-start coord-end ref alt]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [show-bases? ins-format] :or {show-bases? false ins-format :auto}}]
    (apply str (flatten [(coord/format coord-start)
                         (when (and coord-end
                                    (or (not (coord/comparable-coordinates? coord-start coord-end))
                                        (neg? (compare coord-start coord-end))))
                           ["_" (coord/format coord-end)])
                         "del"
                         (when show-bases? ref)
                         "ins"
                         (when (string? alt)
                           (case ins-format
                             :auto (if (and (every? #(= % \N) alt)
                                            (>= (count alt) 10))
                                     (str "N[" (count alt) "]")
                                     alt)
                             :bases alt
                             :count (str "N[" (count alt) "]")))])))
  (plain [this]
    (into {:mutation "dna-indel"} (plain-coords this)))

  Equivalence
  (equiv* [this o]
    (if (instance? DNAIndel o)
      (and (= coord-start (:coord-start o))
           (= coord-end (:coord-end o))
           (if (and ref (:ref o))
             (= ref (:ref o))
             true)
           (= alt (:alt o)))
      false)))

(s/def :clj-hgvs.mutation.dna-indel/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.dna-indel/coord-end (s/nilable ::coord/coordinate))
(s/def :clj-hgvs.mutation.dna-indel/ref (s/nilable ::dna-bases))
(s/def :clj-hgvs.mutation.dna-indel/alt ::dna-bases)
(s/def ::dna-indel
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.dna-indel/coord-start
                          :clj-hgvs.mutation.dna-indel/coord-end
                          :clj-hgvs.mutation.dna-indel/ref
                          :clj-hgvs.mutation.dna-indel/alt])))

(defn- parse-dna-indel-alt
  [s]
  (or (re-matches #"[A-Z]+" s)
      (some-> (re-matches #"N\[(\d+)\]" s)
              (second)
              (intl/parse-long)
              (repeat "N")
              (#(apply str %)))))

(defn dna-indel
  "Constructor of DNAIndel. Throws an exception if any input is illegal."
  [coord-start coord-end ref alt]
  {:pre [(or (nil? coord-end)
             (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))]
   :post [(intl/valid? ::dna-indel %)]}
  (DNAIndel. coord-start coord-end ref alt))

(def ^:private dna-indel-re
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))?del([A-Z]+)?ins(N\[\d+\]|[A-Z]+)")

(defn parse-dna-indel
  [s kind]
  (let [[_ coord-s coord-e ref alt] (re-matches dna-indel-re s)
        parse-coord (coord-parser kind)]
    (dna-indel (parse-coord coord-s)
               (some-> coord-e parse-coord)
               ref
               (parse-dna-indel-alt alt))))

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

(s/def :clj-hgvs.mutation.dna-alleles/mutations1 (s/coll-of ::mutation))
(s/def :clj-hgvs.mutation.dna-alleles/mutations2 (s/nilable (s/coll-of ::mutation)))
(s/def ::dna-alleles
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.dna-alleles/mutations1
                          :clj-hgvs.mutation.dna-alleles/mutations2])))

(defn dna-alleles
  [mutations1 mutations2]
  {:post [(intl/valid? ::dna-alleles %)]}
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
    (let [should-show-end? (and (some? coord-end)
                                (or (not (coord/comparable-coordinates? coord-start coord-end))
                                    (neg? (compare coord-start coord-end))))]
      (str (coord/format coord-start)
           (case range-format
             :auto (or ref (if should-show-end? (str "_" (coord/format coord-end))))
             :bases (or ref (throw (#?(:clj NullPointerException.
                                       :cljs js/Error.) "ref missing")))
             :coord (if (and (nil? coord-end) (> (count ref) 1))
                      (throw (#?(:clj NullPointerException.
                                 :cljs js/Error.) "coord-end missing"))
                      (if should-show-end? (str "_" (coord/format coord-end))))))))
  (format-unique [this _]
    (format-ncopy ncopy))

  Equivalence
  (equiv* [this o]
    (if (instance? DNARepeatedSeqs o)
      (and (= coord-start (:coord-start o))
           (cond
             (and ref (:ref o)) (= ref (:ref o))
             (and coord-end (:coord-end o)) (= coord-end (:coord-end o))
             (and coord-end (:ref o)) (= coord-end
                                         (coord/plus (:coord-start o)
                                                     (dec (count (:ref o)))))
             (and ref (:coord-end o)) (= (coord/plus coord-start
                                                     (dec (count ref)))
                                         (:coord-end o))
             :else true)
           (= ncopy (:ncopy o)))
      false)))

(s/def :clj-hgvs.mutation.dna-repeated-seqs/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.dna-repeated-seqs/coord-end (s/nilable ::coord/coordinate))
(s/def :clj-hgvs.mutation.dna-repeated-seqs/ref (s/nilable ::dna-bases))
(s/def :clj-hgvs.mutation.dna-repeated-seqs/ncopy (s/or :integer integer?
                                                        :vector  vector?))
(s/def ::dna-repeated-seqs
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.dna-repeated-seqs/coord-start
                          :clj-hgvs.mutation.dna-repeated-seqs/coord-end
                          :clj-hgvs.mutation.dna-repeated-seqs/ref
                          :clj-hgvs.mutation.dna-repeated-seqs/ncopy])))

(defn dna-repeated-seqs
  "Constructor of DNARepeatedSeqs. Throws an exception if any input is illegal."
  [coord-start coord-end ref ncopy]
  {:pre [(or (nil? coord-end)
             (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))]
   :post [(intl/valid? ::dna-repeated-seqs %)]}
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
  protocol. kind must be selected from :genome, :mitochondria, :coding-dna, or
  :non-coding-dna."
  [s kind]
  ((condp re-find s
     #"^\((\S+)\)$" parse-uncertain-mutation
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

(s/def ::dna-mutation
  (s/or :uncertain     ::uncertain-mutation
        :substitution  ::dna-substitution
        :deletion      ::dna-deletion
        :duplication   ::dna-duplication
        :insertion     ::dna-insertion
        :inversion     ::dna-inversion
        :conversion    ::dna-conversion
        :indel         ::dna-indel
        :alleles       ::dna-alleles
        :repeated-seqs ::dna-repeated-seqs))

;;; RNA mutations

;; See https://hgvs-nomenclature.org/stable/background/standards/#rna
(s/def ::rna-bases
  (s/and string? #(re-matches #"[acgubdhkmnrsvwy]+" %)))

;;; RNA - substitution
;;;
;;; e.g. r.76a>c
;;;      r.-14g>c
;;;      r.*46u>a

(defrecord RNASubstitution [coord ref alt]
  Mutation
  (format [this] (format this nil))
  (format [this _]
    (str (format-common this _) (format-unique this _)))
  (plain [this]
    (into {:mutation "rna-substitution"} (plain-coords this)))
  SeparatelyFormat
  (format-common [this _] (coord/format coord))
  (format-unique [this _] (str ref ">" alt)))

(s/def :clj-hgvs.mutation.rna-substitution/coord ::coord/coordinate)
(s/def :clj-hgvs.mutation.rna-substitution/ref ::rna-bases)
(s/def :clj-hgvs.mutation.rna-substitution/alt ::rna-bases)
(s/def ::rna-substitution
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.rna-substitution/coord
                          :clj-hgvs.mutation.rna-substitution/ref
                          :clj-hgvs.mutation.rna-substitution/alt])))

(defn rna-substitution
  "Constructor of RNASubstitution. Throws an exception if any input is illegal."
  [coord ref alt]
  {:post [(intl/valid? ::rna-substitution %)]}
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
    (into {:mutation "rna-deletion"} (plain-coords this)))

  Equivalence
  (equiv* [this o]
    (if (instance? RNADeletion o)
      (and (= coord-start (:coord-start o))
           (= coord-end (:coord-end o))
           (if (and ref (:ref o))
             (= ref (:ref o))
             true))
      false)))

(s/def :clj-hgvs.mutation.rna-deletion/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.rna-deletion/coord-end (s/nilable ::coord/coordinate))
(s/def :clj-hgvs.mutation.rna-deletion/ref (s/nilable ::rna-bases))
(s/def ::rna-deletion
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.rna-deletion/coord-start
                          :clj-hgvs.mutation.rna-deletion/coord-end
                          :clj-hgvs.mutation.rna-deletion/ref])))

(defn rna-deletion
  "Constructor of DNAdeletion. Throws an exception if any input is illegal."
  ([coord-start coord-end] (rna-deletion coord-start coord-end nil))
  ([coord-start coord-end ref]
   {:pre [(or (nil? coord-end)
              (not (coord/comparable-coordinates? coord-start coord-end))
              (neg? (compare coord-start coord-end)))]
    :post [(intl/valid? ::rna-deletion %)]}
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
    (into {:mutation "rna-duplication"} (plain-coords this)))

  Equivalence
  (equiv* [this o]
    (if (instance? RNADuplication o)
      (and (= coord-start (:coord-start o))
           (= coord-end (:coord-end o))
           (if (and ref (:ref o))
             (= ref (:ref o))
             true))
      false)))

(s/def :clj-hgvs.mutation.rna-duplication/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.rna-duplication/coord-end (s/nilable ::coord/coordinate))
(s/def :clj-hgvs.mutation.rna-duplication/ref (s/nilable ::rna-bases))
(s/def ::rna-duplication
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.rna-duplication/coord-start
                          :clj-hgvs.mutation.rna-duplication/coord-end
                          :clj-hgvs.mutation.rna-duplication/ref])))

(defn rna-duplication
  "Constructor of RNADuplication. Throws an exception if any input is illegal."
  ([coord-start coord-end] (rna-duplication coord-start coord-end nil))
  ([coord-start coord-end ref]
   {:pre [(or (nil? coord-end)
              (not (coord/comparable-coordinates? coord-start coord-end))
              (neg? (compare coord-start coord-end)))]
    :post [(intl/valid? ::rna-duplication %)]}
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
;;;      r.431_432insn[5]
;;;      r.123_124insL37425.1:23_361

(defrecord RNAInsertion [coord-start coord-end alt]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [ins-format] :or {ins-format :auto}}]
    (str (coord/format coord-start)
         "_"
         (coord/format coord-end)
         "ins"
         (cond
           (map? alt) (str (:genbank alt) ":" (:coord-start alt) "_" (:coord-end alt))
           (re-matches #"n{2,}" alt) (case ins-format
                                       :auto (if (>= (count alt) 10)
                                               (str "n[" (count alt) "]")
                                               alt)
                                       :bases alt
                                       :count (str "n[" (count alt) "]"))
           :else alt)))
  (plain [this]
    (into {:mutation "rna-insertion"} (plain-coords this))))

(s/def :clj-hgvs.mutation.rna-insertion/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.rna-insertion/coord-end ::coord/coordinate)
(s/def :clj-hgvs.mutation.rna-insertion/alt (s/or :rna-bases ::rna-bases
                                                  :ref-seq map?))
(s/def ::rna-insertion
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.rna-insertion/coord-start
                          :clj-hgvs.mutation.rna-insertion/coord-end
                          :clj-hgvs.mutation.rna-insertion/alt])))

(defn rna-insertion
  "Constructor of RNAInsertion. Throws an exception if any input is illegal."
  [coord-start coord-end alt]
  {:pre [(or (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))]
   :post [(intl/valid? ::rna-insertion %)]}
  (RNAInsertion. coord-start coord-end alt))

(defn- parse-rna-alt-n
  [s]
  (when-let [n (some-> (re-find #"n\[(\d+)\]" s)
                       second
                       intl/parse-long)]
    (apply str (repeat n "n"))))

(def ^:private genbank-re
  #"([A-Z]+[\d\.]+):(\d+)_(\d+)")

(defn- parse-rna-alt-genbank
  [s]
  (let [[_ genbank coord-s coord-e] (re-matches genbank-re s)]
    {:genbank genbank
     :coord-start (intl/parse-long coord-s)
     :coord-end (intl/parse-long coord-e)}))

(def ^:private rna-insertion-re
  #"([\d\-\+\*]+)_([\d\-\+\*]+)ins(.+)?")

(defn parse-rna-insertion
  [s]
  (let [[_ coord-s coord-e alt] (re-matches rna-insertion-re s)]
    (rna-insertion (coord/parse-rna-coordinate coord-s)
                   (some-> coord-e coord/parse-rna-coordinate)
                   (cond
                     (re-find #"n\[\d+\]" alt) (parse-rna-alt-n alt)
                     (re-find #"[a-z]+" alt) alt
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

(s/def :clj-hgvs.mutation.rna-inversion/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.rna-inversion/coord-end ::coord/coordinate)
(s/def ::rna-inversion
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.rna-inversion/coord-start
                          :clj-hgvs.mutation.rna-inversion/coord-end])))

(defn rna-inversion
  "Constructor of RNAInversion. Throws an exception if any input is illegal."
  [coord-start coord-end]
  {:pre [(or (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))]
   :post [(intl/valid? ::rna-inversion %)]}
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

(s/def :clj-hgvs.mutation.rna-conversion/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.rna-conversion/coord-end ::coord/coordinate)
(s/def :clj-hgvs.mutation.rna-conversion/alt map?)
(s/def ::rna-conversion
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.rna-conversion/coord-start
                          :clj-hgvs.mutation.rna-conversion/coord-end
                          :clj-hgvs.mutation.rna-conversion/alt])))

(defn rna-conversion
  "Constructor of RNAConversion. Throws an exception if any input is illegal."
  [coord-start coord-end alt]
  {:pre [(or (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))]
   :post [(intl/valid? ::rna-conversion %)]}
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
;;;      r.775_777delinsn[10]

(defrecord RNAIndel [coord-start coord-end ref alt]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [show-bases? ins-format] :or {show-bases? false ins-format :auto}}]
    (str (coord/format coord-start)
         (when (and coord-end
                    (or (not (coord/comparable-coordinates? coord-start coord-end))
                        (neg? (compare coord-start coord-end))))
           (str "_" (coord/format coord-end)))
         "del"
         (when show-bases? ref)
         "ins"
         (if (re-matches #"n{2,}" alt)
           (case ins-format
             :auto (if (>= (count alt) 10)
                     (str "n[" (count alt) "]")
                     alt)
             :bases alt
             :count (str "n[" (count alt) "]"))
           alt)))
  (plain [this]
    (into {:mutation "rna-indel"} (plain-coords this)))

  Equivalence
  (equiv* [this o]
    (if (instance? RNAIndel o)
      (and (= coord-start (:coord-start o))
           (= coord-end (:coord-end o))
           (if (and ref (:ref o))
             (= ref (:ref o))
             true)
           (= alt (:alt o)))
      false)))

(s/def :clj-hgvs.mutation.rna-indel/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.rna-indel/coord-end (s/nilable ::coord/coordinate))
(s/def :clj-hgvs.mutation.rna-indel/ref (s/nilable ::rna-bases))
(s/def :clj-hgvs.mutation.rna-indel/alt ::rna-bases)
(s/def ::rna-indel
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.rna-indel/coord-start
                          :clj-hgvs.mutation.rna-indel/coord-end
                          :clj-hgvs.mutation.rna-indel/ref
                          :clj-hgvs.mutation.rna-indel/alt])))

(defn rna-indel
  "Constructor of RNAIndel. Throws an exception if any input is illegal."
  [coord-start coord-end ref alt]
  {:pre [(or (nil? coord-end)
             (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))]
   :post [(intl/valid? ::rna-indel %)]}
  (RNAIndel. coord-start coord-end ref alt))

(def ^:private rna-indel-re
  #"([\d\-\+\*]+)(?:_([\d\-\+\*]+))?del([a-z]+)?ins([a-z\d\[\]]+)")

(defn parse-rna-indel
  [s]
  (let [[_ coord-s coord-e ref alt] (re-matches rna-indel-re s)]
    (rna-indel (coord/parse-rna-coordinate coord-s)
               (some-> coord-e coord/parse-rna-coordinate)
               ref
               (if (re-find #"n\[\d+\]" alt)
                 (parse-rna-alt-n alt)
                 alt))))

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

(s/def :clj-hgvs.mutation.rna-alleles/mutations1 (s/coll-of ::rna-mutation))
(s/def :clj-hgvs.mutation.rna-alleles/mutations2 (s/nilable (s/coll-of ::rna-mutation)))
(s/def ::rna-alleles
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.rna-alleles/mutations1
                          :clj-hgvs.mutation.rna-alleles/mutations2])))

(defn rna-alleles
  [mutations1 mutations2]
  {:post [(intl/valid? ::rna-alleles %)]}
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
    (let [should-show-end? (and (some? coord-end)
                                (or (not (coord/comparable-coordinates? coord-start coord-end))
                                    (neg? (compare coord-start coord-end))))]
      (str (coord/format coord-start)
           (case range-format
             :auto (or ref (if should-show-end? (str "_" (coord/format coord-end))))
             :bases (or ref (throw (#?(:clj NullPointerException.
                                       :cljs js/Error.) "ref missing")))
             :coord (if (and (nil? coord-end) (> (count ref) 1))
                      (throw (#?(:clj NullPointerException.
                                 :cljs js/Error.) "coord-end missing"))
                      (if should-show-end? (str "_" (coord/format coord-end))))))))
  (format-unique [this _]
    (format-ncopy ncopy))

  Equivalence
  (equiv* [this o]
    (if (instance? RNARepeatedSeqs o)
      (and (= coord-start (:coord-start o))
           (cond
             (and ref (:ref o)) (= ref (:ref o))
             (and coord-end (:coord-end o)) (= coord-end (:coord-end o))
             (and coord-end (:ref o)) (= coord-end
                                         (coord/plus (:coord-start o)
                                                     (dec (count (:ref o)))))
             (and ref (:coord-end o)) (= (coord/plus coord-start
                                                     (dec (count ref)))
                                         (:coord-end o))
             :else true)
           (= ncopy (:ncopy o)))
      false)))

(s/def :clj-hgvs.mutation.rna-repeated-seqs/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.rna-repeated-seqs/coord-end (s/nilable ::coord/coordinate))
(s/def :clj-hgvs.mutation.rna-repeated-seqs/ref (s/nilable ::rna-bases))
(s/def :clj-hgvs.mutation.rna-repeated-seqs/ncopy (s/or :integer integer?
                                                        :vector  vector?))
(s/def ::rna-repeated-seqs
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.rna-repeated-seqs/coord-start
                          :clj-hgvs.mutation.rna-repeated-seqs/coord-end
                          :clj-hgvs.mutation.rna-repeated-seqs/ref
                          :clj-hgvs.mutation.rna-repeated-seqs/ncopy])))

(defn rna-repeated-seqs
  "Constructor of RNARepeatedSeqs. Throws an exception if any input is illegal."
  [coord-start coord-end ref ncopy]
  {:pre [(or (nil? coord-end)
             (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))]
   :post [(intl/valid? ::rna-repeated-seqs %)]}
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

;;; RNA - no RNA detected
;;;
;;; e.g. r.0

(defrecord NoRNA []
  Mutation
  (format [this] (format this nil))
  (format [this _] "0")
  (plain [this] {:mutation "no-rna"}))

(s/def ::no-rna ::mutation)

(defn no-rna
  []
  (NoRNA.))

(defmethod restore "no-rna"
  [m]
  (no-rna))

;;; RNA - unknown mutation
;;;
;;; e.g. r.?

(defrecord RNAUnknownMutation []
  Mutation
  (format [this] (format this nil))
  (format [this _] "?")
  (plain [this] {:mutation "rna-unknown"}))

(s/def ::rna-unknown ::mutation)

(defn rna-unknown-mutation
  []
  (RNAUnknownMutation.))

(defmethod restore "rna-unknown"
  [m]
  (rna-unknown-mutation))

;;; RNA - no effect
;;;
;;; e.g. r.=

(defrecord RNANoEffect []
  Mutation
  (format [this] (format this {}))
  (format [_ _] "=")
  (plain [_] {:mutation "rna-no-effect"}))

(s/def ::rna-no-effect ::mutation)

(defn rna-no-effect
  []
  (RNANoEffect.))

(defmethod restore "rna-no-effect"
  [_]
  (rna-no-effect))

;;; RNA - splice affected
;;;
;;; e.g. r.spl

(defrecord RNASpliceAffected []
  Mutation
  (format [this] (format this {}))
  (format [_ _] "spl")
  (plain [_] {:mutation "rna-splice-affected"}))

(s/def ::rna-splice-affected ::mutation)

(defn rna-splice-affected
  []
  (RNASpliceAffected.))

(defmethod restore "rna-splice-affected"
  [_]
  (rna-splice-affected))

(defn parse-rna
  "Parses a RNA mutation string s, returning a record implementing Mutation
  protocol."
  [s]
  (case s
    "0" (no-rna)
    "?" (rna-unknown-mutation)
    "=" (rna-no-effect)
    "spl" (rna-splice-affected)
    ((condp re-find s
       #"^\((\S+)\)$" #(parse-uncertain-mutation % :rna)
       #"\[.+;.+\]$" parse-rna-alleles
       #"delins" parse-rna-indel
       #"del" parse-rna-deletion
       #"dup" parse-rna-duplication
       #"ins" parse-rna-insertion
       #"inv" parse-rna-inversion
       #"con" parse-rna-conversion
       #"\[[\d\(\)_]+\]" parse-rna-repeated-seqs
       parse-rna-substitution)
     s)))

(s/def ::rna-mutation
  (s/or :uncertain       ::uncertain-mutation
        :substitution    ::rna-substitution
        :deletion        ::rna-deletion
        :duplication     ::rna-duplication
        :substitution    ::rna-substitution
        :insertion       ::rna-insertion
        :inversion       ::rna-inversion
        :conversion      ::rna-conversion
        :indel           ::rna-indel
        :alleles         ::rna-alleles
        :repeated-seqs   ::rna-repeated-seqs
        :no-rna          ::no-rna
        :unknown         ::rna-unknown
        :no-effect       ::rna-no-effect
        :splice-affected ::rna-splice-affected))

;;; Protein mutations

(s/def ::short-amino-acid (set short-amino-acids))
(s/def ::long-amino-acid  (set long-amino-acids))

(s/def ::amino-acid (s/or :short ::short-amino-acid
                          :long  ::long-amino-acid))

(defn- should-show-end?
  [ref-start coord-start ref-end coord-end]
  (and (some? ref-end)
       (or (not (coord/comparable-coordinates? coord-start coord-end))
           (neg? (compare coord-start coord-end)))))

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
    (str (cond-> ref
           (= amino-acid-format :short) ->short-amino-acid)
         (coord/format coord)
         (if (= ref alt)
           "="
           (cond-> alt
             (= amino-acid-format :short) ->short-amino-acid))))
  (plain [this]
    (into {:mutation "protein-substitution"} (plain-coords this))))

(s/def :clj-hgvs.mutation.protein-substitution/ref ::amino-acid)
(s/def :clj-hgvs.mutation.protein-substitution/coord ::coord/coordinate)
(s/def :clj-hgvs.mutation.protein-substitution/alt ::amino-acid)
(s/def ::protein-substitution
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.protein-substitution/ref
                          :clj-hgvs.mutation.protein-substitution/coord
                          :clj-hgvs.mutation.protein-substitution/alt])))

(defn protein-substitution
  "Constructor of ProteinSubstitution. Throws an exception if any input is illegal."
  [ref coord alt]
  {:post [(intl/valid? ::protein-substitution %)]}
  (ProteinSubstitution. ref coord alt))

(def ^:private protein-substitution-re
  #"([A-Z*](?:[a-z]{2})?)(\d+)([A-Z*=](?:[a-z]{2})?)")

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

(s/def :clj-hgvs.mutation.protein-deletion/ref-start ::amino-acid)
(s/def :clj-hgvs.mutation.protein-deletion/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.protein-deletion/ref-end (s/nilable ::amino-acid))
(s/def :clj-hgvs.mutation.protein-deletion/coord-end (s/nilable ::coord/coordinate))
(s/def ::protein-deletion
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.protein-deletion/ref-start
                          :clj-hgvs.mutation.protein-deletion/coord-start
                          :clj-hgvs.mutation.protein-deletion/ref-end
                          :clj-hgvs.mutation.protein-deletion/coord-end])))

(defn protein-deletion
  "Constructor of ProteinDeletion. Throws an exception if any input is illegal."
  ([ref-start coord-start] (protein-deletion ref-start coord-start nil nil))
  ([ref-start coord-start ref-end coord-end]
   {:pre [(or (nil? coord-start)
              (not (coord/comparable-coordinates? coord-start coord-end))
              (neg? (compare coord-start coord-end)))]
    :post [(intl/valid? ::protein-deletion %)]}
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

(s/def :clj-hgvs.mutation.protein-duplication/ref-start ::amino-acid)
(s/def :clj-hgvs.mutation.protein-duplication/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.protein-duplication/ref-end (s/nilable ::amino-acid))
(s/def :clj-hgvs.mutation.protein-duplication/coord-end (s/nilable ::coord/coordinate))
(s/def ::protein-duplication
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.protein-duplication/ref-start
                          :clj-hgvs.mutation.protein-duplication/coord-start
                          :clj-hgvs.mutation.protein-duplication/ref-end
                          :clj-hgvs.mutation.protein-duplication/coord-end])))

(defn protein-duplication
  "Constructor of ProteinDuplication. Throws an exception if any input is illegal."
  ([ref-start coord-start] (protein-duplication ref-start coord-start nil nil))
  ([ref-start coord-start ref-end coord-end]
   {:pre [(or (nil? coord-end)
              (not (coord/comparable-coordinates? coord-start coord-end))
              (neg? (compare coord-start coord-end)))]
    :post [(intl/valid? ::protein-duplication %)]}
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
;;;      Arg78_Gly79insX[23]
;;;      Thr89_Ter90insSerProIle (T89_*90insSPI)

(defrecord ProteinInsertion [ref-start coord-start ref-end coord-end alts]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [amino-acid-format ins-format] :or {amino-acid-format :long ins-format :auto}}]
    (apply str (flatten [(cond-> ref-start
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-start)
                         "_"
                         (cond-> ref-end
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-end)
                         "ins"
                         (if (every? #(= % "Xaa") alts)
                           (let [alts (cond->> alts
                                        (= amino-acid-format :short) (map ->short-amino-acid))]
                             (case ins-format
                               :auto (if (>= (count alts) 10)
                                       (str "X[" (count alts) "]")
                                       alts)
                               :amino-acids alts
                               :count (str "X[" (count alts) "]")))
                           (cond->> alts
                             (= amino-acid-format :short) (map ->short-amino-acid)))])))
  (plain [this]
    (into {:mutation "protein-insertion"} (plain-coords this))))

(s/def :clj-hgvs.mutation.protein-insertion/ref-start ::amino-acid)
(s/def :clj-hgvs.mutation.protein-insertion/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.protein-insertion/ref-end ::amino-acid)
(s/def :clj-hgvs.mutation.protein-insertion/coord-end ::coord/coordinate)
(s/def :clj-hgvs.mutation.protein-insertion/alts (s/coll-of ::amino-acid))
(s/def ::protein-insertion
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.protein-insertion/ref-start
                          :clj-hgvs.mutation.protein-insertion/coord-start
                          :clj-hgvs.mutation.protein-insertion/ref-end
                          :clj-hgvs.mutation.protein-insertion/coord-end
                          :clj-hgvs.mutation.protein-insertion/alts])))

(defn protein-insertion
  "Constructor of ProteinInsertion. Throws an exception if any input is illegal."
  [ref-start coord-start ref-end coord-end alts]
  {:pre [(or (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))]
   :post [(intl/valid? ::protein-insertion %)]}
  (ProteinInsertion. ref-start coord-start ref-end coord-end alts))

(defn- parse-protein-insertion-alts
  [s]
  (condp re-matches s
    #"([A-Z*]([a-z]{2})?)+" (mapv ->long-amino-acid (re-seq #"[A-Z*](?:[a-z]{2})?" s))
    #"X\[\d+\]" (-> (re-find #"X\[(\d+)\]" s)
                    second
                    intl/parse-long
                    (repeat "Xaa")
                    vec)))

(def ^:private protein-insertion-re
  #"([A-Z](?:[a-z]{2})?)(\d+)_([A-Z*](?:[a-z]{2})?)(\d+)ins([\da-zA-Z*\[\]]+)")

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
;;;      Cys28_Lys29delinsX[10]

(defrecord ProteinIndel [ref-start coord-start ref-end coord-end alts]
  Mutation
  (format [this] (format this nil))
  (format [this {:keys [amino-acid-format ins-format] :or {amino-acid-format :long ins-format :auto}}]
    (apply str (flatten [(cond-> ref-start
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-start)
                         (when (should-show-end? ref-start coord-start ref-end coord-end)
                           ["_"
                            (cond-> ref-end
                              (= amino-acid-format :short) ->short-amino-acid)
                            (coord/format coord-end)])
                         "delins"
                         (if (every? #(= % "Xaa") alts)
                           (let [alts (cond->> alts
                                        (= amino-acid-format :short) (map ->short-amino-acid))]
                             (case ins-format
                               :auto (if (>= (count alts) 10)
                                       (str "X[" (count alts) "]")
                                       alts)
                               :amino-acids alts
                               :count (str "X[" (count alts) "]")))
                           (cond->> alts
                             (= amino-acid-format :short) (map ->short-amino-acid)))])))
  (plain [this]
    (into {:mutation "protein-indel"} (plain-coords this))))

(s/def :clj-hgvs.mutation.protein-indel/ref-start ::amino-acid)
(s/def :clj-hgvs.mutation.protein-indel/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.protein-indel/ref-end (s/nilable ::amino-acid))
(s/def :clj-hgvs.mutation.protein-indel/coord-end (s/nilable ::coord/coordinate))
(s/def :clj-hgvs.mutation.protein-indel/alts (s/coll-of ::amino-acid))
(s/def ::protein-indel
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.protein-indel/ref-start
                          :clj-hgvs.mutation.protein-indel/coord-start
                          :clj-hgvs.mutation.protein-indel/ref-end
                          :clj-hgvs.mutation.protein-indel/coord-end
                          :clj-hgvs.mutation.protein-indel/alts])))

(defn protein-indel
  "Constructor of ProteinIndel. Throws an exception if any input is illegal."
  [ref-start coord-start ref-end coord-end alts]
  {:pre [(or (nil? coord-end)
             (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))]
   :post [(intl/valid? ::protein-indel %)]}
  (ProteinIndel. ref-start coord-start ref-end coord-end alts))

(def ^:private protein-indel-re
  #"([A-Z](?:[a-z]{2})?)(\d+)(?:_([A-Z](?:[a-z]{2})?)(\d+))?delins([A-Z*][a-zA-Z*\[\]\d]*)?")

(defn parse-protein-indel
  [s]
  (let [[_ ref-s coord-s ref-e coord-e alts] (re-matches protein-indel-re s)]
    (protein-indel (->long-amino-acid ref-s)
                   (coord/parse-protein-coordinate coord-s)
                   (->long-amino-acid ref-e)
                   (some-> coord-e coord/parse-protein-coordinate)
                   (mapv ->long-amino-acid (some->> alts parse-protein-insertion-alts)))))

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

(s/def :clj-hgvs.mutation.protein-alleles/mutations1 (s/coll-of ::protein-mutation))
(s/def :clj-hgvs.mutation.protein-alleles/mutations2 (s/nilable (s/coll-of ::protein-mutation)))
(s/def ::protein-alleles
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.protein-alleles/mutations1
                          :clj-hgvs.mutation.protein-alleles/mutations2])))

(defn protein-alleles
  [mutations1 mutations2]
  {:post [(intl/valid? ::protein-alleles %)]}
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

(s/def :clj-hgvs.mutation.protein-repeated-seqs/ref-start ::amino-acid)
(s/def :clj-hgvs.mutation.protein-repeated-seqs/coord-start ::coord/coordinate)
(s/def :clj-hgvs.mutation.protein-repeated-seqs/ref-end (s/nilable ::amino-acid))
(s/def :clj-hgvs.mutation.protein-repeated-seqs/coord-end (s/nilable ::coord/coordinate))
(s/def :clj-hgvs.mutation.protein-repeated-seqs/ncopy (s/or :integer integer?
                                                            :vector  vector?))
(s/def ::protein-repeated-seqs
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.protein-repeated-seqs/ref-start
                          :clj-hgvs.mutation.protein-repeated-seqs/coord-start
                          :clj-hgvs.mutation.protein-repeated-seqs/ref-end
                          :clj-hgvs.mutation.protein-repeated-seqs/coord-end
                          :clj-hgvs.mutation.protein-repeated-seqs/ncopy])))

(defn protein-repeated-seqs
  "Constructor of ProteinRepeatedSeqs. Throws an exception if any input is illegal."
  [ref-start coord-start ref-end coord-end ncopy]
  {:pre [(or (nil? coord-end)
             (not (coord/comparable-coordinates? coord-start coord-end))
             (neg? (compare coord-start coord-end)))]
   :post [(intl/valid? ::protein-repeated-seqs %)]}
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
  (format [this {:keys [amino-acid-format show-ter-site? ter-format]
                 :or {amino-acid-format :long, show-ter-site? false, ter-format :long}}]
    (str (cond-> ref
           (= amino-acid-format :short) ->short-amino-acid)
         (coord/format coord)
         (if (= ref alt)
           "="
           (cond-> alt
             (= amino-acid-format :short) ->short-amino-acid))
         "fs"
         (if (and show-ter-site? (some? new-ter-site))
           (str (cond-> "Ter"
                  (= ter-format :short) ->short-amino-acid)
                (coord/format new-ter-site)))))
  (plain [this]
    (into {:mutation "protein-frame-shift"} (plain-coords this)))

  Equivalence
  (equiv* [this o]
    (if (instance? ProteinFrameShift o)
      (and (= ref (:ref o))
           (= coord (:coord o))
           (if (and alt (:alt o))
             (= alt (:alt o))
             true)
           (if (and (instance? clj_hgvs.coordinate.ProteinCoordinate new-ter-site)
                    (instance? clj_hgvs.coordinate.ProteinCoordinate (:new-ter-site o)))
             (= new-ter-site (:new-ter-site o))
             true))
      false)))

(s/def :clj-hgvs.mutation.protein-frame-shift/ref ::amino-acid)
(s/def :clj-hgvs.mutation.protein-frame-shift/coord ::coord/coordinate)
(s/def :clj-hgvs.mutation.protein-frame-shift/alt (s/nilable ::amino-acid))
(s/def :clj-hgvs.mutation.protein-frame-shift/new-ter-site (s/nilable ::coord/coordinate))
(s/def ::protein-frame-shift
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.protein-frame-shift/ref
                          :clj-hgvs.mutation.protein-frame-shift/coord
                          :clj-hgvs.mutation.protein-frame-shift/alt
                          :clj-hgvs.mutation.protein-frame-shift/new-ter-site])))

(defn protein-frame-shift
  "Constructor of ProteinFrameShift. Throws an exception if any input is illegal."
  [ref coord alt new-ter-site]
  {:post [(intl/valid? ::protein-frame-shift %)]}
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
                 :or {amino-acid-format :long, ter-format :long}}]
    (str (case ref
           "Met" (cond-> ref (= amino-acid-format :short) ->short-amino-acid)
           "Ter" (cond-> ref (= ter-format :short) ->short-amino-acid))
         (coord/format coord)
         (cond-> alt
           (= amino-acid-format :short) ->short-amino-acid)
         "ext"
         (coord/->region-str region)
         (coord/format new-site)))
  (plain [this]
    (into {:mutation "protein-extension"}
          (update (plain-coords this) :region name)))

  Equivalence
  (equiv* [this o]
    (if (instance? ProteinExtension o)
      (and (= ref (:ref o))
           (= coord (:coord o))
           (if (and alt (:alt o))
             (= alt (:alt o))
             true)
           (= region (:region o))
           (if (and (instance? clj_hgvs.coordinate.ProteinCoordinate new-site)
                    (instance? clj_hgvs.coordinate.ProteinCoordinate (:new-site o)))
             (= new-site (:new-site o))
             true))
      false)))

(s/def :clj-hgvs.mutation.protein-extension/ref ::amino-acid)
(s/def :clj-hgvs.mutation.protein-extension/coord ::coord/coordinate)
(s/def :clj-hgvs.mutation.protein-extension/alt (s/nilable ::amino-acid))
(s/def :clj-hgvs.mutation.protein-extension/region #{:upstream :downstream})
(s/def :clj-hgvs.mutation.protein-extension/new-site ::coord/coordinate)
(s/def ::protein-extension
  (s/and ::mutation
         (s/keys :req-un [:clj-hgvs.mutation.protein-extension/ref
                          :clj-hgvs.mutation.protein-extension/coord
                          :clj-hgvs.mutation.protein-extension/alt
                          :clj-hgvs.mutation.protein-extension/region
                          :clj-hgvs.mutation.protein-extension/new-site])))

(defn protein-extension
  "Constructor of ProteinExtension. Throws an exception if any input is illegal."
  [ref coord alt region new-site]
  {:post [(intl/valid? ::protein-extension %)]}
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
    (protein-extension ref coord alt (keyword region) new-site)))

;;; Protein - no protein detected
;;;
;;; e.g. p.0

(defrecord NoProtein []
  Mutation
  (format [this] (format this nil))
  (format [this _] "0")
  (plain [this] {:mutation "no-protein"}))

(s/def ::no-protein ::mutation)

(defn no-protein
  []
  (NoProtein.))

(defmethod restore "no-protein"
  [m]
  (no-protein))

;;; Protein - unknown mutation
;;;
;;; e.g. p.?
;;;      p.Met1?

(defrecord ProteinUnknownMutation []
  Mutation
  (format [this] (format this nil))
  (format [this _] "?")
  (plain [this] {:mutation "protein-unknown"}))

(s/def ::protein-unknown ::mutation)

(defn protein-unknown-mutation
  []
  (ProteinUnknownMutation.))

(defmethod restore "protein-unknown"
  [m]
  (protein-unknown-mutation))

;;; Protein - no effect
;;;
;;; e.g. p.=

(defrecord ProteinNoEffect []
  Mutation
  (format [this] (format this {}))
  (format [_ _] "=")
  (plain [_] {:mutation "protein-no-effect"}))

(s/def ::protein-no-effect ::mutation)

(defn protein-no-effect
  []
  (ProteinNoEffect.))

(defmethod restore "protein-no-effect"
  [_]
  (protein-no-effect))

(defn parse-protein
  "Parses a protein mutation string s, returning a record implementing Mutation
  protocol."
  [s]
  (case s
    "0" (no-protein)
    "=" (protein-no-effect)
    (condp re-find s
      #"^(M(et)?1)?\?$" (protein-unknown-mutation)
      #"^\((\S+)\)$" (parse-uncertain-mutation s :protein)
      #"\[.+;.+\]$" (parse-protein-alleles s)
      #"delins" (parse-protein-indel s)
      #"del" (parse-protein-deletion s)
      #"dup" (parse-protein-duplication s)
      #"ins" (parse-protein-insertion s)
      #"fs" (parse-protein-frame-shift s)
      #"ext" (parse-protein-extension s)
      #"\[[\d\(\)_]+\]" (parse-protein-repeated-seqs s)
      (parse-protein-substitution s))))

(s/def ::protein-mutation
  (s/or :uncertain     ::uncertain-mutation
        :substitution  ::protein-substitution
        :deletion      ::protein-deletion
        :duplication   ::protein-duplication
        :insertion     ::protein-insertion
        :indel         ::protein-indel
        :alleles       ::protein-alleles
        :repeated-seqs ::protein-repeated-seqs
        :frame-shift   ::protein-frame-shift
        :extension     ::protein-extension
        :no-protein    ::no-protein
        :unknown       ::protein-unknown
        :no-effect     ::protein-no-effect))

(defn parse
  [s kind]
  (let [parse* (case kind
                 (:genome :mitochondria :coding-dna :non-coding-dna :circular-dna) #(parse-dna % kind)
                 :rna parse-rna
                 :protein parse-protein)]
    (parse* s)))

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
