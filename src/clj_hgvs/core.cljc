(ns clj-hgvs.core
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.string :as string]
            [clj-hgvs.internal :refer [->kind-keyword ->kind-str]]
            [clj-hgvs.mutation :as mut]))

(defn- mutation-parser
  [kind]
  (cond
    (#{:genome :mitochondria :cdna :ncdna} kind) #(mut/parse-dna % kind)
    (= kind :rna) mut/parse-rna
    (= kind :protein) mut/parse-protein))

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
  [mutations & opts]
  (let [multi? (> (count mutations) 1)]
    (apply str (flatten [(if multi? "[")
                         (->> mutations
                              (map #(mut/format % opts))
                              (string/join ";"))
                         (if multi? "]")]))))

(defn format
  [hgvs & opts]
  (apply str [(format-transcript (:transcript hgvs))
              (format-kind (:kind hgvs))
              (apply format-mutations (:mutations hgvs) opts)]))
