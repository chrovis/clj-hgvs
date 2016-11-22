(ns clj-hgvs.mutation
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.string :as string]))

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
    "=" :unchanged
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
    :unchanged "="
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
   "V"])

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
   "Val"])

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

(defmulti parse (fn [_ kind] kind))

;; Genome
;;
;; dup: g.123_345dup
;; ins: g.123_124insAGC
;; inv: g.123_345inv
;; con: g.123_345con888_1110 (TODO)
;; delins: g.123_127delinsAG
;; repeat: g.123_125[36], g.123GGC[36] (TODO)

(defrecord GenomeMutation [numbering type ref alt]
  Mutation
  (format [_ _]
    (string/join [numbering ref (->mutation-str type) alt])))

(def ^:private genome-mutation-re
  #"^([\d_\-\+\*\?]+)([A-Z]+)?(>|del|dup|ins|delins|inv|con|fs|ext)([A-Z]+)?$")

(defmethod parse :genome
  [s _]
  (let [[_ numbering ref type alt] (re-find genome-mutation-re s)]
    (map->GenomeMutation {:numbering numbering
                          :type (->mutation-type-keyword type)
                          :ref ref
                          :alt alt})))

;; Mitochondria

(defrecord MitochondriaMutation [numbering type ref alt]
  Mutation
  (format [_ _]
    (string/join [numbering ref (->mutation-str type) alt])))

(def ^:private mitochondria-mutation-re
  #"^([\d_\-\+\*\?]+)([A-Z]+)?(>|del|dup|ins|delins|inv|con|fs|ext)([A-Z]+)?$")

(defmethod parse :mitochondria
  [s kind]
  (let [[_ numbering ref type alt] (re-find genome-mutation-re s)]
    (map->MitochondriaMutation {:numbering numbering
                                :type (->mutation-type-keyword type)
                                :ref ref
                                :alt alt})))

;; Coding DNA

(defrecord CodingDNAMutation [numbering type ref alt]
  Mutation
  (format [_ _]
    (string/join [numbering ref (->mutation-str type) alt])))

(def ^:private coding-dna-mutation-re
  #"^([\d_\-\+\*\?]+)([A-Z]+)?(>|del|dup|ins|delins|inv|con|fs|ext)([A-Z]+)?$")

(defmethod parse :coding-dna
  [s _]
  (let [[_ numbering ref type alt] (re-find genome-mutation-re s)]
    (map->CodingDNAMutation {:numbering numbering
                             :type (->mutation-type-keyword type)
                             :ref ref
                             :alt alt})))

;; Non-coding DNA

(defrecord NonCodingDNAMutation [numbering type ref alt]
  Mutation
  (format [_ _]
    (string/join [numbering ref (->mutation-str type) alt])))

(def ^:private non-coding-dna-mutation-re
  #"^([\d_\-\+\*\?]+)([A-Z]+)?(>|del|dup|ins|delins|inv|con|fs|ext)([A-Z]+)?$")

(defmethod parse :non-coding-dna
  [s _]
  (let [[_ numbering ref type alt] (re-find genome-mutation-re s)]
    (map->NonCodingDNAMutation {:numbering numbering
                             :type (->mutation-type-keyword type)
                             :ref ref
                             :alt alt})))

;; RNA

(defrecord RNAMutation [numbering type ref alt]
  Mutation
  (format [_ _]
    (string/join [numbering ref (->mutation-str type) alt])))

(def ^:private rna-mutation-re
  #"^([\d_\-\+\*\?]+)([A-Z]+)?(>|del|dup|ins|delins|inv|con|fs|ext)([A-Z]+)?$")

(defmethod parse :rna
  [s _]
  (let [[_ numbering ref type alt] (re-find genome-mutation-re s)]
    (map->RNAMutation {:numbering numbering
                       :type (->mutation-type-keyword type)
                       :ref ref
                       :alt alt})))

;; Protein
;;
;; del: p.Cys76_Glu79del (TODO)
;; ins: p.Lys23_Leu24insArgSerGln (TODO)
;; fs: Arg123LysfsTer34
;; ext: p.Met1ext-5, p.Ter110GlnextTer17

(defrecord ProteinMutation [numbering type ref alt rest]
  Mutation
  (format [_ {:keys [amino-acid-format]
              :or {amino-acid-format :long}}]
    (string/join [(cond-> ref
                    (= amino-acid-format :short) ->short-amino-acid)
                  numbering
                  (if (#{:unchanged :deletion :insertion :indel} type)
                    (->mutation-str type))
                  (cond-> alt
                    (= amino-acid-format :short) ->short-amino-acid)
                  (if (#{:frame-shift :extension} type)
                    (->mutation-str type))
                  rest])))

(def ^:private protein-mutation-re
  #"^([A-Z](?:[a-z]{2})?)(\d+)(=|del|dup|ins|delins)?([A-Z](?:[a-z]{2})?)?(?:(fs|ext)(.+)?)?$")

(defmethod parse :protein
  [s _]
  (let [[_ ref numbering type1 alt type2 rest*] (re-find protein-mutation-re s)]
    (map->ProteinMutation {:numbering numbering
                           :type (->mutation-type-keyword (or type1 type2 ">"))
                           :ref (->long-amino-acid ref)
                           :alt (->long-amino-acid alt)
                           :rest rest*})))
