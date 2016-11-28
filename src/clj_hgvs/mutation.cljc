(ns clj-hgvs.mutation
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.string :as string]
            [clj-hgvs.coordinate :as coord]
            [clj-hgvs.internal :refer [parse-long]]))

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
;;;
;;; dup: g.123_345dup
;;; ins: g.123_124insAGC
;;; inv: g.123_345inv
;;; con: g.123_345con888_1110 (TODO)
;;; delins: g.123_127delinsAG
;;; repeat: g.123_125[36], g.123GGC[36] (TODO)

;;; DNA - substitution
;;;
;;; e.g. 45576A>C
;;;      88+1G>T
;;;      76_77delinsTT
;;;      123G=
;;;      85C=/>T
;;;      85C=//>T

(defrecord DNASubstitution [start end ref alt]
  Mutation
  (format [_ _]

    ))

(def ^:private dna-substitution-re
  #"")

(defn parse-dna-substitution
  [s]

)

(defrecord GenomeMutation [numbering type ref alt]
  Mutation
  (format [_ _]
    (string/join [numbering ref (->mutation-str type) alt])))

(def ^:private genome-mutation-re
  #"^([\d_\-\+\*\?]+)([A-Z]+)?(>|del|dup|ins|delins|inv|con|fs|ext)([A-Z]+)?$")

(defn parse-genome
  [s]
  (let [[_ numbering ref type alt] (re-find genome-mutation-re s)]
    (map->GenomeMutation {:numbering numbering
                          :type (->mutation-type-keyword type)
                          :ref ref
                          :alt alt})))

;;; Mitochondria

(defrecord MitochondriaMutation [numbering type ref alt]
  Mutation
  (format [_ _]
    (string/join [numbering ref (->mutation-str type) alt])))

(def ^:private mitochondria-mutation-re
  #"^([\d_\-\+\*\?]+)([A-Z]+)?(>|del|dup|ins|delins|inv|con|fs|ext)([A-Z]+)?$")

(defn parse-mitochondria
  [s]
  (let [[_ numbering ref type alt] (re-find genome-mutation-re s)]
    (map->MitochondriaMutation {:numbering numbering
                                :type (->mutation-type-keyword type)
                                :ref ref
                                :alt alt})))

;;; Coding DNA

(defrecord CDNAMutation [numbering type ref alt]
  Mutation
  (format [_ _]
    (string/join [numbering ref (->mutation-str type) alt])))

(def ^:private coding-dna-mutation-re
  #"^([\d_\-\+\*\?]+)([A-Z]+)?(>|del|dup|ins|delins|inv|con|fs|ext)([A-Z]+)?$")

(defn parse-cdna
  [s]
  (let [[_ numbering ref type alt] (re-find genome-mutation-re s)]
    (map->CDNAMutation {:numbering numbering
                        :type (->mutation-type-keyword type)
                        :ref ref
                        :alt alt})))

;;; Non-coding DNA

(defrecord NCDNAMutation [numbering type ref alt]
  Mutation
  (format [_ _]
    (string/join [numbering ref (->mutation-str type) alt])))

(def ^:private non-coding-dna-mutation-re
  #"^([\d_\-\+\*\?]+)([A-Z]+)?(>|del|dup|ins|delins|inv|con|fs|ext)([A-Z]+)?$")

(defn parse-ncdna
  [s]
  (let [[_ numbering ref type alt] (re-find genome-mutation-re s)]
    (map->NCDNAMutation {:numbering numbering
                         :type (->mutation-type-keyword type)
                         :ref ref
                         :alt alt})))

;;; RNA

(defrecord RNAMutation [numbering type ref alt]
  Mutation
  (format [_ _]
    (string/join [numbering ref (->mutation-str type) alt])))

(def ^:private rna-mutation-re
  #"^([\d_\-\+\*\?]+)([A-Z]+)?(>|del|dup|ins|delins|inv|con|fs|ext)([A-Z]+)?$")

(defn parse-rna
  [s]
  (let [[_ numbering ref type alt] (re-find genome-mutation-re s)]
    (map->RNAMutation {:numbering numbering
                       :type (->mutation-type-keyword type)
                       :ref ref
                       :alt alt})))

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
  #"^([A-Z](?:[a-z]{2})?)(\d+)(?:_([A-Z](?:[a-z]{2})?)(\d+))?ins([A-Z][a-zA-Z]*)?$")

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
