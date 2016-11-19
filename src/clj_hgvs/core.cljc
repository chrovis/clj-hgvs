(ns clj-hgvs.core
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

(defn- ->kind-keyword
  [s]
  (case s
    "g" :genome
    "m" :mitochondria
    "c" :coding-dna
    "n" :non-coding-dna
    "r" :rna
    "p" :protein))

(defn- ->kind-str
  [k]
  (case k
    :genome "g"
    :mitochondria "m"
    :coding-dna "c"
    :non-coding-dna "n"
    :rna "r"
    :protein"p"))

(defn- split-mutations
  [s]
  (map #(string/replace % #"[\[\]]" "") (string/split s #";")))

(def ^:private mutation-re
  #"^([\d_\-\+\*\?]+)([a-zA-Z]+)?(>|del|dup|ins|inv|con|fs|ext)([a-zA-Z]+)?$")

(defn- parse-mutation*
  [s]
  (let [[_ numbering ref type alt] (re-find mutation-re s)]
    {:numbering numbering
     :type (->mutation-type-keyword type)
     :ref ref
     :alt alt}))

;; del: p.Cys76_Glu79del (TODO)
;; fs: Arg123LysfsTer34
;; ext: p.Met1ext-5, p.Ter110GlnextTer17

(def ^:private protein-mutation-re
  #"^([A-Z](?:[a-z]{2})?)(\d+)(=|del|dup|ins|delins)?([A-Z](?:[a-z]{2})?)?(?:(fs|ext)(.+)?)?$")

(defn- parse-protein-mutation
  [s]
  (let [[_ ref numbering type1 alt type2 rest*] (re-find protein-mutation-re s)]
    {:numbering numbering
     :type (->mutation-type-keyword (or type1 type2 ">"))
     :ref ref
     :alt alt
     :rest rest*}))

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
                (string? mutation) (mapv #(parse-mutation % kind) (split-mutations mutation)))})

(def ^:private hgvs-re #"^(?:([^:]+):)?([gmcnrp])\.(.+)$")

(defn parse
  [s]
  (let [[_ transcript kind mutations] (re-find hgvs-re s)
        kind-k (->kind-keyword kind)]
    (hgvs transcript kind-k mutations)))

(defn- format-transcript
  [transcript]
  (if transcript
    (str transcript ":")))

(defn- format-kind
  [kind]
  (str (->kind-str kind) "."))

(defn- format-mutation*
  [{:keys [numbering type ref alt]}]
  (string/join [numbering ref (->mutation-str type) alt]))

(defn- format-protein-mutation
  [{:keys [numbering type ref alt rest]}]
  (string/join [ref
                numbering
                (if (#{:unchanged :deletion :insertion :indel} type)
                  (->mutation-str type))
                alt
                (if (#{:frame-shift :extension} type)
                  (->mutation-str type))
                rest]))

(defn- format-mutation
  [mutation kind]
  ((cond
    (#{:genome :mitochondria :coding-dna :non-coding-dna :rna} kind) format-mutation*
    (= kind :protein) format-protein-mutation)
   mutation))

(defn- format-mutations
  [mutations kind]
  (let [multi? (> (count mutations) 1)]
    (apply str (flatten [(if multi? "[")
                         (string/join ";" (map #(format-mutation % kind) mutations))
                         (if multi? "]")]))))

(defn format
  [hgvs]
  (apply str [(format-transcript (:transcript hgvs))
              (format-kind (:kind hgvs))
              (format-mutations (:mutations hgvs) (:kind hgvs))]))
