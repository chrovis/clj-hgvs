(ns clj-hgvs.core
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.string :as string]
            [clj-hgvs.mutation :as mut]))

(defn- ->kind-keyword
  [s]
  (case s
    "g" :genome
    "m" :mitochondria
    "c" :cdna
    "n" :ncdna
    "r" :rna
    "p" :protein))

(defn- ->kind-str
  [k]
  (case k
    :genome "g"
    :mitochondria "m"
    :cdna "c"
    :ncdna "n"
    :rna "r"
    :protein"p"))

(defn- mutation-parser
  [kind]
  (case kind
    :genome mut/parse-genome
    :mitochondria mut/parse-mitochondria
    :cdna mut/parse-cdna
    :ncdna mut/parse-ncdna
    :rna mut/parse-rna
    :protein mut/parse-protein))

(defn- split-mutations
  [s]
  (map #(string/replace % #"[\[\]]" "") (string/split s #";")))

(defn hgvs
  [transcript kind mutation & mutations]
  {:transcript transcript
   :kind kind
   :mutations (cond
                (map? mutation) (cons mutation mutations)
                (string? mutation) (mapv (mutation-parser kind) (split-mutations mutation)))})

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

(defn format-mutations
  [mutations & {:keys [amino-acid-format] :or {amino-acid-format :long}}]
  (let [multi? (> (count mutations) 1)]
    (apply str (flatten [(if multi? "[")
                         (->> mutations
                              (map #(mut/format % {:amino-acid-format amino-acid-format}))
                              (string/join ";"))
                         (if multi? "]")]))))

(defn format
  [hgvs & opts]
  (apply str [(format-transcript (:transcript hgvs))
              (format-kind (:kind hgvs))
              (apply format-mutations (:mutations hgvs) opts)]))
