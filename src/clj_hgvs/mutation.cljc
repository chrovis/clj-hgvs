(ns clj-hgvs.mutation
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.string :as string]))

(defn- parse-long
  [s]
  #?(:clj (Long/parseLong s)
     :cljs (js/parseInt s)))

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

;;; Genome
;;;
;;; dup: g.123_345dup
;;; ins: g.123_124insAGC
;;; inv: g.123_345inv
;;; con: g.123_345con888_1110 (TODO)
;;; delins: g.123_127delinsAG
;;; repeat: g.123_125[36], g.123GGC[36] (TODO)

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

(defn format-protein-coord
  "{:amino-acid \"Ala\", :position 3} => Ala3"
  [{:keys [amino-acid position]} amino-acid-format]
  {:pre [(#{:long :short} amino-acid-format)]}
  (str (cond-> amino-acid
         (= amino-acid-format :short) ->short-amino-acid)
       position))

(defn parse-protein-coord
  "Ala3 => {:amino-acid \"Ala\", :position 3}"
  [s]
  (if s
    (if-let [[_ amino-acid pos] (re-find #"^([A-Z](?:[a-z]{2})?)(\d+)$" s)]
      {:amino-acid (->long-amino-acid amino-acid)
       :position (parse-long pos)})))

;;; Protein - substitution
;;;
;;; e.g. Arg54Ser
;;;      Trp26Ter
;;;      Trp26*
;;;      Cys123=

(defrecord ProteinSubstitution [coord alt]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    {:pre [(#{:long :short} amino-acid-format)]}
    (str (format-protein-coord coord amino-acid-format)
         (if (= (:amino-acid coord) alt)
           "="
           (cond-> alt
             (= amino-acid-format :short) ->short-amino-acid)))))

(def ^:private protein-substitution-re
  #"^([A-Z](?:[a-z]{2})?\d+)([A-Z\*=](?:[a-z]{2})?)$")

(defn parse-protein-substitution
  [s]
  (let [[_ coord alt] (re-find protein-substitution-re s)
        coord' (parse-protein-coord coord)]
    (map->ProteinSubstitution {:coord coord'
                               :alt (case alt
                                      "=" (->long-amino-acid (:amino-acid coord'))
                                      "*" "Ter"
                                      (->long-amino-acid alt))})))

;;; Protein - deletion
;;;
;;; e.g. Ala3del
;;;      Cys76_Glu79del

(defrecord ProteinDeletion [start end]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (apply str
           (format-protein-coord start amino-acid-format)
           (if end "_")
           (format-protein-coord end amino-acid-format)
           "del")))

(def ^:private protein-deletion-re
  #"^([A-Z](?:[a-z]{2})?\d+)(?:_([A-Z](?:[a-z]{2})?\d+))?del$")

(defn parse-protein-deletion
  [s]
  (let [[_ start end] (re-find protein-deletion-re s)]
    (map->ProteinDeletion {:start (parse-protein-coord start)
                           :end (parse-protein-coord end)})))

;;; Protein - duplication
;;;
;;; e.g. Ala3dup
;;;      Ala3_Ser5dup

(defrecord ProteinDuplication [start end]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (str (format-protein-coord start amino-acid-format)
         (if end "_")
         (format-protein-coord end amino-acid-format)
         "dup")))

(def ^:private protein-duplication-re
  #"^([A-Z](?:[a-z]{2})?\d+)(?:_([A-Z](?:[a-z]{2})?\d+))?dup$")

(defn parse-protein-duplication
  [s]
  (let [[_ start end] (re-find protein-duplication-re s)]
    (map->ProteinDuplication {:start (parse-protein-coord start)
                              :end (parse-protein-coord end)})))

;;; Protein - insertion
;;;
;;; e.g. Lys23_Leu24insArgSerGln

(defrecord ProteinInsertion [start end alts]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (apply str
           (format-protein-coord start amino-acid-format)
           (if end "_")
           (format-protein-coord end amino-acid-format)
           "ins"
           (cond->> alts
             (= amino-acid-format :short) (map ->short-amino-acid)))))

(def ^:private protein-insertion-re
  #"^([A-Z](?:[a-z]{2})?\d+)(?:_([A-Z](?:[a-z]{2})?\d+))?ins([A-Z][a-zA-Z]*)?$")

(defn parse-protein-insertion
  [s]
  (let [[_ start end alts] (re-find protein-insertion-re s)]
    (map->ProteinInsertion {:start (parse-protein-coord start)
                            :end (parse-protein-coord end)
                            :alts (mapv ->long-amino-acid (some->> alts (re-seq #"[A-Z](?:[a-z]{2})?")))})))

;;; Protein - indel
;;;
;;; e.g. Cys28delinsTrpVal
;;;      Cys28_Lys29delinsTrp

(defrecord ProteinIndel [start end alts]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (apply str
           (format-protein-coord start amino-acid-format)
           (if end "_")
           (format-protein-coord end amino-acid-format)
           "delins"
           (cond->> alts
             (= amino-acid-format :short) (map ->short-amino-acid)))))

(def ^:private protein-indel-re
  #"^([A-Z](?:[a-z]{2})?\d+)(?:_([A-Z](?:[a-z]{2})?\d+))?delins([A-Z][a-zA-Z]*)?$")

(defn parse-protein-indel
  [s]
  (let [[_ start end alts] (re-find protein-indel-re s)]
    (map->ProteinIndel {:start (parse-protein-coord start)
                        :end (parse-protein-coord end)
                        :alts (mapv ->long-amino-acid (some->> alts (re-seq #"[A-Z](?:[a-z]{2})?")))})))

;;; Protein - repeated sequences
;;;
;;; e.g. Ala2[10]
;;;      Ala2[10];[11]
;;;      Arg65_Ser67[12]

(defrecord ProteinRepeatedSeqs [start end ncopy ncopy-other]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (str (format-protein-coord start amino-acid-format)
         (if end "_")
         (format-protein-coord end amino-acid-format)
         "[" ncopy "]"
         (if ncopy-other (str ";[" ncopy-other "]")))))

(def ^:private protein-repeated-seqs-re
  #"^([A-Z](?:[a-z]{2})?\d+)(?:_([A-Z](?:[a-z]{2})?\d+))?\[(\d+)\](?:;\[(\d+)\])?")

(defn parse-protein-repeated-seqs
  [s]
  (let [[_ start end ncopy1 ncopy2] (re-find protein-repeated-seqs-re s)]
    (map->ProteinRepeatedSeqs {:start (parse-protein-coord start)
                               :end (parse-protein-coord end)
                               :ncopy (parse-long ncopy1)
                               :ncopy-other (if ncopy2 (parse-long ncopy2))})))

;;; Protein - frame shift
;;;
;;; e.g. Arg97ProfsTer23
;;;      Arg97fs
;;;      Ile327Argfs*?
;;;      Gln151Thrfs*9

(defrecord ProteinFrameShift [coord alt new-site]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (str (format-protein-coord coord amino-acid-format)
         (cond-> alt
           (= amino-acid-format :short) ->short-amino-acid)
         "fs"
         new-site)))

(def ^:private protein-frame-shift-re
  #"^([A-Z](?:[a-z]{2})?\d+)([A-Z](?:[a-z]{2})?)?fs(.+)?$")

(defn parse-protein-frame-shift
  [s]
  (let [[_ coord alt new-site] (re-find protein-frame-shift-re s)]
    (map->ProteinFrameShift {:coord (parse-protein-coord coord)
                             :alt (->long-amino-acid alt)
                             :new-site new-site})))

;;; Protein - extension
;;;
;;; e.g. Met1ext-5
;;;      Met1Valext-12
;;;      Ter110GlnextTer17

(defrecord ProteinExtension [coord alt new-site]
  Mutation
  (format [_ {:keys [amino-acid-format] :or {amino-acid-format :long}}]
    (str (format-protein-coord coord amino-acid-format)
         (cond-> alt
           (= amino-acid-format :short) ->short-amino-acid)
         "ext"
         new-site)))

(def ^:private protein-extension-re
  #"^([A-Z](?:[a-z]{2})?\d+)([A-Z](?:[a-z]{2})?)?ext(.+)$")

(defn parse-protein-extension
  [s]
  (let [[_ coord alt new-site] (re-find protein-extension-re s)]
    (map->ProteinExtension {:coord (parse-protein-coord coord)
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
