(ns clj-hgvs.mutation
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.string :as string]
            [clj-hgvs.coordinate :as coord]
            [clj-hgvs.internal :refer [parse-long ->kind-keyword ->kind-str]]))

(defn- ->mutation-type-keyword
  [s]
  (case s
    ">" :substitution
    "del" :deletion
    "ins" :insertion
    "delins" :indel
    "inv" :inversion
    "con" :conversion
    "fs" :frame-shift
    "ext" :extension
    "dup" :duplication))

(defn- ->mutation-str
  [k]
  (case k
    :substitution ">"
    :deletion "del"
    :insertion "ins"
    :indel "delins"
    :inversion "inv"
    :conversion "con"
    :frame-shift "fs"
    :extension "ext"
    :duplication "dup"))

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
  [s]
  (if ((set long-amino-acids) s)
    s
    (get short-long-amino-acid-map s)))

(defn ->short-amino-acid
  [s]
  (if ((set short-amino-acids) s)
    s
    (get long-short-amino-acid-map s)))

(defprotocol Mutation
  (format [this opts]))

;;; DNA mutations

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
  (format [_ _]
    (apply str (flatten [(coord/format coord)
                         ref
                         type
                         alt]))))

(def ^:private dna-substitution-re
  #"^([\d\-\+]+)([A-Z]+)([>=/]+)([A-Z]+)?$")

(defn parse-dna-substitution
  [s kind]
  (let [[_ coord ref type alt] (re-find dna-substitution-re s)
        parse-coord (coord-parser kind)]
    (map->DNASubstitution {:coord (parse-coord coord)
                           :ref ref
                           :type type
                           :alt alt})))

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
  (format [_ {:keys [show-bases?] :or {show-bases? false}}]
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

(def ^:private dna-deletion-re
  #"^([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))?del([A-Z]+)?")

(defn- parse-dna-deletion*
  [s kind]
  (let [[_ coord-s coord-e ref] (re-find dna-deletion-re s)
        parse-coord (coord-parser kind)]
    (map->DNADeletion {:coord-start (parse-coord coord-s)
                       :coord-end (some-> coord-e parse-coord)
                       :ref ref})))

(def ^:private dna-deletion-range-re
  #"^\(([\d\-\+\?\*]+)_([\d\-\+\?\*]+)\)_\(([\d\-\+\?\*]+)_([\d\-\+\?\*]+)\)del([A-Z]+)?")

(defn- parse-dna-deletion-range
  [s kind]
  (let [[_ coord-s1 coord-s2 coord-e1 coord-e2 ref] (re-find dna-deletion-range-re s)
        parse-coord (coord-parser kind)]
    (map->DNADeletion {:coord-start (mapv parse-coord [coord-s1 coord-s2])
                       :coord-end (mapv parse-coord [coord-e1 coord-e2])
                       :ref ref})))

(defn parse-dna-deletion
  [s kind]
  (if (re-find #"^\(.+\)_\(.+\)del[A-Z]*$" s)
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
  (format [_ {:keys [show-bases?] :or {show-bases? false}}]
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

(def ^:private dna-duplication-re
  #"^([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))?dup([A-Z]+)?")

(defn- parse-dna-duplication*
  [s kind]
  (let [[_ coord-s coord-e ref] (re-find dna-duplication-re s)
        parse-coord (coord-parser kind)]
    (map->DNADuplication {:coord-start (parse-coord coord-s)
                          :coord-end (some-> coord-e parse-coord)
                          :ref ref})))

(def ^:private dna-duplication-range-re
  #"^\(([\d\-\+\?\*]+)_([\d\-\+\?\*]+)\)_\(([\d\-\+\?\*]+)_([\d\-\+\?\*]+)\)dup([A-Z]+)?")

(defn- parse-dna-duplication-range
  [s kind]
  (let [[_ coord-s1 coord-s2 coord-e1 coord-e2 ref] (re-find dna-duplication-range-re s)
        parse-coord (coord-parser kind)]
    (map->DNADuplication {:coord-start (mapv parse-coord [coord-s1 coord-s2])
                          :coord-end (mapv parse-coord [coord-e1 coord-e2])
                          :ref ref})))

(defn parse-dna-duplication
  [s kind]
  (if (re-find #"^\(.+\)_\(.+\)dup[A-Z]*$" s)
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
  (format [_ _]
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

(def ^:private dna-insertion-re
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))ins(.+)")

(def ^:private dna-insertion-alt-re
  #"(?:([^:]+):)([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))")

(defn parse-dna-insertion
  [s kind]
  (let [[_ coord-s coord-e alt] (re-matches dna-insertion-re s)
        parse-coord (coord-parser kind)]
    (map->DNAInsertion
     {:coord-start (parse-coord coord-s)
      :coord-end (parse-coord coord-e)
      :alt (or (re-matches #"[A-Z]+" alt)
               (let [[_ transcript coord-s coord-e] (re-matches dna-insertion-alt-re alt)]
                 {:transcript transcript
                  :coord-start (parse-coord coord-s)
                  :coord-end (parse-coord coord-e)}))})))

;;; DNA - inversion
;;;
;;; e.g. g.1077_1080inv
;;;      c.77_80inv

(defrecord DNAInversion [coord-start coord-end]
  Mutation
  (format [_ _]
    (str (coord/format coord-start)
         "_"
         (coord/format coord-end)
         "inv")))

(def ^:private dna-inversion-re
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))inv")

(defn parse-dna-inversion
  [s kind]
  (let [[_ coord-s coord-e] (re-matches dna-inversion-re s)
        parse-coord (coord-parser kind)]
    (map->DNAInversion {:coord-start (parse-coord coord-s)
                        :coord-end (parse-coord coord-e)})))

;;; DNA - conversion
;;;
;;; e.g. g.333_590con1844_2101
;;;      g.415_1655conAC096506.5:g.409_1683
;;;      c.15_355conNM_004006.1:20_360

(defrecord DNAConversion [coord-start coord-end alt]
  Mutation
  (format [_ _]
    (apply str (flatten [(coord/format coord-start)
                         "_"
                         (coord/format coord-end)
                         "con"
                         (if (:transcript alt) [(:transcript alt) ":"])
                         (if (:kind alt) [(->kind-str (:kind alt)) "."])
                         (coord/format (:coord-start alt))
                         "_"
                         (coord/format (:coord-end alt))]))))

(def ^:private dna-conversion-re
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))con(.+)")

(def ^:private dna-conversion-alt-re
  #"(?:([^:]+):)?(?:([gmcn])\.)?([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))")

(defn parse-dna-conversion
  [s kind]
  (let [[_ coord-s coord-e alt] (re-matches dna-conversion-re s)
        parse-coord (coord-parser kind)]
    (map->DNAConversion
     {:coord-start (parse-coord coord-s)
      :coord-end (parse-coord coord-e)
      :alt (let [[_ transcript kind coord-s coord-e] (re-matches dna-conversion-alt-re alt)
                 parse-alt-coord (if kind
                                   (coord-parser (->kind-keyword kind))
                                   parse-coord)]
             {:transcript transcript
              :kind (some-> kind ->kind-keyword)
              :coord-start (parse-alt-coord coord-s)
              :coord-end (parse-alt-coord coord-e)})})))

;;; DNA - indel
;;;
;;; e.g. g.6775delinsGA
;;;      c.145_147delinsTGG

(defrecord DNAIndel [coord-start coord-end alt]
  Mutation
  (format [_ _]
    (apply str (flatten [(coord/format coord-start)
                         (if coord-end ["_" (coord/format coord-end)])
                         "delins"
                         alt]))))

(def ^:private dna-indel-re
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))?delins([A-Z]+)")

(defn parse-dna-indel
  [s kind]
  (let [[_ coord-s coord-e alt] (re-matches dna-indel-re s)
        parse-coord (coord-parser kind)]
    (map->DNAIndel {:coord-start (parse-coord coord-s)
                    :coord-end (some-> coord-e parse-coord)
                    :alt alt})))

;;; DNA - repeated sequences
;;;
;;; e.g. g.123_124[14]

(defrecord DNARepeatedSeqs [coord-start coord-end ncopy]
  Mutation
  (format [_ _]
    (apply str (flatten [(coord/format coord-start)
                         "_" (coord/format coord-end)
                         "[" ncopy "]"]))))

(def ^:private dna-repeated-seqs-re
  #"([\d\-\+\*\?]+)(?:_([\d\-\+\*\?]+))\[(\d+)\]")

(defn parse-dna-repeated-seqs
  [s kind]
  (let [[_ coord-s coord-e ncopy] (re-matches dna-repeated-seqs-re s)
        parse-coord (coord-parser kind)]
    (map->DNARepeatedSeqs {:coord-start (parse-coord coord-s)
                           :coord-end (parse-coord coord-e)
                           :ncopy (parse-long ncopy)})))

(defn parse-dna
  [s kind]
  ((cond
     (string/includes? s "delins") parse-dna-indel
     (string/includes? s "del") parse-dna-deletion
     (string/includes? s "dup") parse-dna-duplication
     (string/includes? s "ins") parse-dna-insertion
     (string/includes? s "inv") parse-dna-inversion
     (string/includes? s "con") parse-dna-conversion
     (re-find #"\[\d+\]" s) parse-dna-repeated-seqs
     :else parse-dna-substitution)
   s kind))

;;; RNA mutations

;;; RNA - substitution
;;;
;;; e.g. r.76a>c
;;;      r.-14g>c
;;;      r.*46u>a

(defrecord RNASubstitution [coord ref alt]
  Mutation
  (format [_ _]
    (apply str (coord/format coord) ref ">" alt)))

(def ^:private rna-substitution-re
  #"([\d\-\+\*]+)([a-z]?)>([a-z]?)")

(defn parse-rna-substitution
  [s]
  (let [[_ coord ref alt] (re-matches rna-substitution-re s)]
    (map->RNASubstitution {:coord (coord/parse-rna-coordinate coord)
                           :ref ref
                           :alt alt})))

;;; RNA - deletion
;;;
;;; e.g. r.7del
;;;      r.7delu
;;;      r.6_8del
;;;      r.(4072_5145)del (TODO)

(defrecord RNADeletion [coord-start coord-end ref]
  Mutation
  (format [_ {:keys [show-bases?] :or {show-bases? false}}]
    (str (coord/format coord-start)
         (if coord-end
           (str "_" (coord/format coord-end)))
         "del"
         (if show-bases? ref))))

(def ^:private rna-deletion-re
  #"([\d\-\+\*]+)(?:_([\d\-\+\*]+))?del([a-z]+)?")

(defn parse-rna-deletion
  [s]
  (let [[_ coord-s coord-e ref] (re-matches rna-deletion-re s)]
    (map->RNADeletion {:coord-start (coord/parse-rna-coordinate coord-s)
                       :coord-end (some-> coord-e coord/parse-rna-coordinate)
                       :ref ref})))

;;; RNA - duplication
;;;
;;; e.g. r.7dup
;;;      r.7dupu
;;;      r.6_8dup

(defrecord RNADuplication [coord-start coord-end ref]
  Mutation
  (format [_ {:keys [show-bases?] :or {show-bases? false}}]
    (str (coord/format coord-start)
         (if coord-end
           (str "_" (coord/format coord-end)))
         "dup"
         (if show-bases? ref))))

(def ^:private rna-duplication-re
  #"([\d\-\+\*]+)(?:_([\d\-\+\*]+))?dup([a-z]+)?")

(defn parse-rna-duplication
  [s]
  (let [[_ coord-s coord-e ref] (re-matches rna-duplication-re s)]
    (map->RNADuplication {:coord-start (coord/parse-rna-coordinate coord-s)
                          :coord-end (some-> coord-e coord/parse-rna-coordinate)
                          :ref ref})))

;;; RNA - insertion
;;;
;;; e.g. r.756_757insacu
;;;      r.431_432ins(5)
;;;      r.123_124insL37425.1:23_361

(defrecord RNAInsertion [coord-start coord-end alt]
  Mutation
  (format [_ _]
    (str (coord/format coord-start)
         "_"
         (coord/format coord-end)
         "ins"
         (cond
           (map? alt) (str (:genbank alt) ":" (:coord-start alt) "_" (:coord-end alt))
           (re-matches #"n{2,}" alt) (str "(" (count alt) ")")
           :else alt))))

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
    (map->RNAInsertion {:coord-start (coord/parse-rna-coordinate coord-s)
                        :coord-end (some-> coord-e coord/parse-rna-coordinate)
                        :alt (cond
                               (re-find #"[a-z]+" alt) alt
                               (re-find #"\(\d\)" alt) (parse-rna-alt-n alt)
                               :else (parse-rna-alt-genbank alt))})))

;;; RNA - inversion
;;;
;;; e.g. r.177_180inv

(defrecord RNAInversion [coord-start coord-end]
  Mutation
  (format [_ _]
    (str (coord/format coord-start)
         "_"
         (coord/format coord-end)
         "inv")))

(def ^:private rna-inversion-re
  #"([\d\-\+\*]+)_([\d\-\+\*]+)inv")

(defn parse-rna-inversion
  [s]
  (let [[_ coord-s coord-e] (re-matches rna-inversion-re s)]
    (map->RNAInversion {:coord-start (coord/parse-rna-coordinate coord-s)
                        :coord-end (coord/parse-rna-coordinate coord-e)})))

;;; RNA - conversion
;;;
;;; e.g. r.123_345con888_1110
;;;      r.415_1655conAC096506.5:409_1649

(defrecord RNAConversion [coord-start coord-end alt]
  Mutation
  (format [_ _]
    (str (coord/format coord-start)
         "_"
         (coord/format coord-end)
         "con"
         (some-> (:transcript alt) (str ":"))
         (coord/format (:coord-start alt))
         "_"
         (coord/format (:coord-end alt)))))

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
    (map->RNAConversion {:coord-start (coord/parse-rna-coordinate coord-s)
                         :coord-end (coord/parse-rna-coordinate coord-e)
                         :alt (parse-rna-conversion-alt alt)})))

;;; RNA - indel
;;;
;;; e.g. r.775delinsga
;;;      r.775_777delinsc

(defrecord RNAIndel [coord-start coord-end alt]
  Mutation
  (format [_ _]
    (str (coord/format coord-start)
         (some->> coord-end coord/format (str "_"))
         "delins"
         alt)))

(def ^:private rna-indel-re
  #"([\d\-\+\*]+)(?:_([\d\-\+\*]+))?delins([a-z]+)")

(defn parse-rna-indel
  [s]
  (let [[_ coord-s coord-e alt] (re-matches rna-indel-re s)]
    (map->RNAIndel {:coord-start (coord/parse-rna-coordinate coord-s)
                    :coord-end (some-> coord-e coord/parse-rna-coordinate)
                    :alt alt})))

;;; RNA - repeated sequences
;;;
;;; e.g. r.-124_-123[14]
;;;      r.-124ug[14]
;;;      r.-124_-123[14];[18]

(defrecord RNARepeatedSeqs [coord-start coord-end ref ncopy ncopy-other]
  Mutation
  (format [_ _]
    (str (coord/format coord-start)
         (some->> coord-end coord/format (str "_"))
         ref
         "[" ncopy "]"
         (if ncopy-other (str ";[" ncopy-other "]")))))

(def ^:private rna-repeated-seqs-re
  #"([\d\-\+\*]+)(?:_([\d\-\+\*]+))?([a-z]+)?\[(\d+)\](?:;\[(\d+)\])?")

(defn parse-rna-repeated-seqs
  [s]
  (let [[_ coord-s coord-e ref ncopy1 ncopy2] (re-matches rna-repeated-seqs-re s)]
    (map->RNARepeatedSeqs {:coord-start (coord/parse-rna-coordinate coord-s)
                           :coord-end (some-> coord-e coord/parse-rna-coordinate)
                           :ref ref
                           :ncopy (parse-long ncopy1)
                           :ncopy-other (some-> ncopy2 parse-long)})))

(defn parse-rna
  [s]
  ((cond
     (string/includes? s "delins") parse-rna-indel
     (string/includes? s "del") parse-rna-deletion
     (string/includes? s "dup") parse-rna-duplication
     (string/includes? s "ins") parse-rna-insertion
     (string/includes? s "inv") parse-rna-inversion
     (string/includes? s "con") parse-rna-conversion
     (re-find #"\[\d+\]" s) parse-rna-repeated-seqs
     :else parse-rna-substitution)
   s))

;;; Protein mutations

;;; Protein - substitution
;;;
;;; e.g. Arg54Ser
;;;      Trp26Ter
;;;      Trp26*
;;;      Cys123=

(defrecord ProteinSubstitution [coord ref alt]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    {:pre [(#{:long :short} amino-acid-format)]}
    (str (cond-> ref
           (= amino-acid-format :short) ->short-amino-acid)
         (coord/format coord)
         (if (= ref alt)
           "="
           (cond-> alt
             (= amino-acid-format :short) ->short-amino-acid)))))

(def ^:private protein-substitution-re
  #"^([A-Z](?:[a-z]{2})?)(\d+)([A-Z\*=](?:[a-z]{2})?)$")

(defn parse-protein-substitution
  [s]
  (let [[_ ref coord' alt] (re-find protein-substitution-re s)]
    (map->ProteinSubstitution {:coord (coord/parse-protein-coordinate coord')
                               :ref (->long-amino-acid ref)
                               :alt (case alt
                                      "=" (->long-amino-acid ref)
                                      "*" "Ter"
                                      (->long-amino-acid alt))})))

;;; Protein - deletion
;;;
;;; e.g. Ala3del
;;;      Cys76_Glu79del

(defrecord ProteinDeletion [ref-start coord-start ref-end coord-end]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (apply str (flatten [(cond-> ref-start
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-start)
                         (if ref-end
                           ["_"
                            (cond-> ref-end
                              (= amino-acid-format :short) ->short-amino-acid)
                            (coord/format coord-end)])
                         "del"]))))

(def ^:private protein-deletion-re
  #"^([A-Z](?:[a-z]{2})?)(\d+)(?:_([A-Z](?:[a-z]{2})?)(\d+))?del$")

(defn parse-protein-deletion
  [s]
  (let [[_ ref-s coord-s ref-e coord-e] (re-find protein-deletion-re s)]
    (map->ProteinDeletion {:ref-start (->long-amino-acid ref-s)
                           :coord-start (coord/parse-protein-coordinate coord-s)
                           :ref-end (->long-amino-acid ref-e)
                           :coord-end (some-> coord-e coord/parse-protein-coordinate)})))

;;; Protein - duplication
;;;
;;; e.g. Ala3dup
;;;      Ala3_Ser5dup

(defrecord ProteinDuplication [ref-start coord-start ref-end coord-end]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (apply str (flatten [(cond-> ref-start
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-start)
                         (if ref-end
                           ["_"
                            (cond-> ref-end
                              (= amino-acid-format :short) ->short-amino-acid)
                            (coord/format coord-end)])
                         "dup"]))))

(def ^:private protein-duplication-re
  #"^([A-Z](?:[a-z]{2})?)(\d+)(?:_([A-Z](?:[a-z]{2})?)(\d+))?dup$")

(defn parse-protein-duplication
  [s]
  (let [[_ ref-s coord-s ref-e coord-e] (re-find protein-duplication-re s)]
    (map->ProteinDuplication {:ref-start (->long-amino-acid ref-s)
                              :coord-start (coord/parse-protein-coordinate coord-s)
                              :ref-end (->long-amino-acid ref-e)
                              :coord-end (some-> coord-e coord/parse-protein-coordinate)})))

;;; Protein - insertion
;;;
;;; e.g. Lys23_Leu24insArgSerGln

(defrecord ProteinInsertion [ref-start coord-start ref-end coord-end alts]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (apply str (flatten [(cond-> ref-start
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-start)
                         (if ref-end
                           ["_"
                            (cond-> ref-end
                              (= amino-acid-format :short) ->short-amino-acid)
                            (coord/format coord-end)])
                         "ins"
                         (cond->> alts
                           (= amino-acid-format :short) (map ->short-amino-acid))]))))

(def ^:private protein-insertion-re
  #"^([A-Z](?:[a-z]{2})?)(\d+)(?:_([A-Z](?:[a-z]{2})?)(\d+))?ins([A-Z][a-zA-Z]*)?$")\

(defn parse-protein-insertion
  [s]
  (let [[_ ref-s coord-s ref-e coord-e alts] (re-find protein-insertion-re s)]
    (map->ProteinInsertion {:ref-start (->long-amino-acid ref-s)
                            :coord-start (coord/parse-protein-coordinate coord-s)
                            :ref-end (->long-amino-acid ref-e)
                            :coord-end (some-> coord-e coord/parse-protein-coordinate)
                            :alts (mapv ->long-amino-acid (some->> alts (re-seq #"[A-Z](?:[a-z]{2})?")))})))

;;; Protein - indel
;;;
;;; e.g. Cys28delinsTrpVal
;;;      Cys28_Lys29delinsTrp

(defrecord ProteinIndel [ref-start coord-start ref-end coord-end alts]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (apply str (flatten [(cond-> ref-start
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-start)
                         (if ref-end
                           ["_"
                            (cond-> ref-end
                              (= amino-acid-format :short) ->short-amino-acid)
                            (coord/format coord-end)])
                         "delins"
                         (cond->> alts
                           (= amino-acid-format :short) (map ->short-amino-acid))]))))

(def ^:private protein-indel-re
  #"^([A-Z](?:[a-z]{2})?)(\d+)(?:_([A-Z](?:[a-z]{2})?)(\d+))?delins([A-Z][a-zA-Z]*)?$")

(defn parse-protein-indel
  [s]
  (let [[_ ref-s coord-s ref-e coord-e alts] (re-find protein-indel-re s)]
    (map->ProteinIndel {:ref-start (->long-amino-acid ref-s)
                        :coord-start (coord/parse-protein-coordinate coord-s)
                        :ref-end (->long-amino-acid ref-e)
                        :coord-end (some-> coord-e coord/parse-protein-coordinate)
                        :alts (mapv ->long-amino-acid (some->> alts (re-seq #"[A-Z](?:[a-z]{2})?")))})))

;;; Protein - repeated sequences
;;;
;;; e.g. Ala2[10]
;;;      Ala2[10];[11]
;;;      Arg65_Ser67[12]

(defrecord ProteinRepeatedSeqs [ref-start coord-start ref-end coord-end ncopy
                                ncopy-other]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (apply str (flatten [(cond-> ref-start
                           (= amino-acid-format :short) ->short-amino-acid)
                         (coord/format coord-start)
                         (if ref-end
                           ["_"
                            (cond-> ref-end
                              (= amino-acid-format :short) ->short-amino-acid)
                            (coord/format coord-end)])
                         "[" ncopy "]"
                         (if ncopy-other [";[" ncopy-other "]"])]))))

(def ^:private protein-repeated-seqs-re
  #"^([A-Z](?:[a-z]{2})?)(\d+)(?:_([A-Z](?:[a-z]{2})?)(\d+))?\[(\d+)\](?:;\[(\d+)\])?")

(defn parse-protein-repeated-seqs
  [s]
  (let [[_ ref-s coord-s ref-e coord-e ncopy1 ncopy2] (re-find protein-repeated-seqs-re s)]
    (map->ProteinRepeatedSeqs {:ref-start (->long-amino-acid ref-s)
                               :coord-start (coord/parse-protein-coordinate coord-s)
                               :ref-end (->long-amino-acid ref-e)
                               :coord-end (some-> coord-e coord/parse-protein-coordinate)
                               :ncopy (parse-long ncopy1)
                               :ncopy-other (if ncopy2 (parse-long ncopy2))})))

;;; Protein - frame shift
;;;
;;; e.g. Arg97ProfsTer23
;;;      Arg97fs
;;;      Ile327Argfs*?
;;;      Gln151Thrfs*9

(defrecord ProteinFrameShift [ref coord alt new-site]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (str (cond-> ref
           (= amino-acid-format :short) ->short-amino-acid)
         (coord/format coord)
         (cond-> alt
           (= amino-acid-format :short) ->short-amino-acid)
         "fs"
         new-site)))

(def ^:private protein-frame-shift-re
  #"^([A-Z](?:[a-z]{2})?)(\d+)([A-Z](?:[a-z]{2})?)?fs(.+)?$")

(defn parse-protein-frame-shift
  [s]
  (let [[_ ref coord' alt new-site] (re-find protein-frame-shift-re s)]
    (map->ProteinFrameShift {:ref (->long-amino-acid ref)
                             :coord (coord/parse-protein-coordinate coord')
                             :alt (->long-amino-acid alt)
                             :new-site new-site})))

;;; Protein - extension
;;;
;;; e.g. Met1ext-5
;;;      Met1Valext-12
;;;      Ter110GlnextTer17

(defrecord ProteinExtension [ref coord alt new-site]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (str (cond-> ref
           (= amino-acid-format :short) ->short-amino-acid)
         (coord/format coord)
         (cond-> alt
           (= amino-acid-format :short) ->short-amino-acid)
         "ext"
         new-site)))

(def ^:private protein-extension-re
  #"^([A-Z](?:[a-z]{2})?)(\d+)([A-Z](?:[a-z]{2})?)?ext(.+)$")

(defn parse-protein-extension
  [s]
  (let [[_ ref coord' alt new-site] (re-find protein-extension-re s)]
    (map->ProteinExtension {:ref (->long-amino-acid ref)
                             :coord (coord/parse-protein-coordinate coord')
                            :alt (->long-amino-acid alt)
                            :new-site new-site})))

(defn parse-protein
  [s]
  ((cond
     (string/includes? s "delins") parse-protein-indel
     (string/includes? s "del") parse-protein-deletion
     (string/includes? s "dup") parse-protein-duplication
     (string/includes? s "ins") parse-protein-insertion
     (string/includes? s "fs") parse-protein-frame-shift
     (string/includes? s "ext") parse-protein-extension
     (re-find #"\[\d+\]" s) parse-protein-repeated-seqs
     :else parse-protein-substitution)
   s))
