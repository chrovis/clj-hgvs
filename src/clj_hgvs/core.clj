(ns clj-hgvs.core
  (:require [clojure.string :as string]))

(defn- ->mutation-type-keyword
  [s]
  (case s
    ">" :substitution
    "del" :deletion
    "ins" :insertion
    "inv" :inversion
    "con" :conversion
    "fs" :frame-shift
    "ext" :extension))

(defn- ->kind-keyword
  [s]
  (case s
    "g" :genome
    "m" :mitochondria
    "c" :coding-dna
    "n" :non-coding-dna
    "r" :rna
    "p" :protein))

(defn- split-mutations
  [s]
  (map #(string/replace % #"[\[\]]" "") (string/split s #";")))

(def ^:private mutation-re #"^([\d_\-\+\*\?]+)([a-zA-Z]+)?(>|del|dup|ins|inv|con|fs|ext)([a-zA-Z]+)$")

(defn- parse-mutation*
  [s]
  (let [[_ numbering ref type alt] (re-find mutation-re s)]
    {:numbering numbering
     :type (->mutation-type-keyword type)
     :ref ref
     :alt alt}))

(def ^:private protein-mutation-re #"^([a-zA-Z]+)(\d+)([a-zA-Z]+)$")

(defn- parse-protein-mutation
  [s]
  (let [[_ ref numbering alt] (re-find protein-mutation-re s)]
    {:numbering numbering
     :type :substitution
     :ref ref
     :alt alt}))

(defn- parse-mutation
  [s kind]
  ((cond
    (#{:genome :mitochondria :coding-dna :non-coding-dna :rna} kind) parse-mutation*
    (= kind :protein) parse-protein-mutation)
   s))

(defn hgvs
  [transcript kind mutation & mutations]
  {:transcript transcript
   :kind kind
   :mutations (cond
                (map? mutation) (cons mutation mutations)
                (string? mutation) (map #(parse-mutation % kind) (split-mutations mutation)))})

(def ^:private hgvs-re #"^(?:([^:]+):)?([gmcnrp])\.(.+)$")

(defn parse
  [s]
  (let [[_ transcript kind mutations] (re-find hgvs-re s)
        kind-k (->kind-keyword kind)]
    (hgvs transcript kind-k mutations)))
