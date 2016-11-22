(ns clj-hgvs.core
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.string :as string]
            [clj-hgvs.mutation :as mut]))

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

(defn hgvs
  [transcript kind mutation & mutations]
  {:transcript transcript
   :kind kind
   :mutations (cond
                (map? mutation) (cons mutation mutations)
                (string? mutation) (mapv #(mut/parse % kind) (split-mutations mutation)))})

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
                         (string/join ";" (map #(mut/format % opts) mutations))
                         (if multi? "]")]))))

(defn format
  [hgvs & opts]
  (apply str [(format-transcript (:transcript hgvs))
              (format-kind (:kind hgvs))
              (apply format-mutations (:mutations hgvs) opts)]))
