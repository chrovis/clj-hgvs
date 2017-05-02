(ns clj-hgvs.core
  "Functions to handle HGVS. See http://varnomen.hgvs.org/ for the detail HGVS
  nomenclature."
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

;; [123456A>G;345678G>C]
;; [123456A>G];[345678G>C]
;; 112GAT[14]
(defn- split-mutations
  [s]
  {:pre [(= (count (re-seq #"\[" s)) (count (re-seq #"\]" s)))]}
  (if (re-find #";" s)
    (condp re-find s
      #"\];\[" (mapv #(subs % 1 (dec (count %))) (string/split s #";"))
      #";" (string/split (subs s 1 (dec (count s))) #";"))
    [s]))

(defn hgvs
  "Constructor of HGVS map."
  [transcript kind mutation & mutations]
  {:transcript transcript
   :kind kind
   :mutations (cond
                (map? mutation) (cons mutation mutations)
                (string? mutation) (mapv (mutation-parser kind) (split-mutations mutation)))})

(def ^:private hgvs-re #"^(?:([^:]+):)?([gmcnrp])\.(.+)$")

(defn parse
  "Parses a HGVS string s, returning a map representing the HGVS."
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
  "Returns a string representing mutations. The second argument is an optional
  map to specify style. See document of clj-hgvs.core/format for details of the
  option."
  ([mutations] (format-mutations mutations nil))
  ([mutations opts]
   (let [multi? (> (count mutations) 1)]
     (apply str (flatten [(if multi? "[")
                          (->> mutations
                               (map #(mut/format % opts))
                               (string/join ";"))
                          (if multi? "]")])))))

(defn format
  "Returns a HGVS string representing the given HGVS map. The second argument is
  an optional map to specify style:
    :show-bases? - displays additional bases, e.g. g.6_8delTGC, default false.
    :range-format - range style, default :auto. <:auto|:bases|:coord>
    :amino-acid-format - amino acid style of protein HGVS, default :long. <:long|:short>
    :ter-format - ter codon style of protein frame shift. <:long|:short>"
  ([hgvs] (format hgvs nil))
  ([hgvs opts]
   (apply str [(format-transcript (:transcript hgvs))
               (format-kind (:kind hgvs))
               (format-mutations (:mutations hgvs) opts)])))

(defn plain
  "Returns a plain map representing the given HGVS. This function is useful for
  sending data through another codec. Use restore to retrieve original HGVS
  data."
  [hgvs]
  (-> hgvs
      (update :kind name)
      (update :mutations #(map mut/plain %))))

(defn restore
  "Restores a plain map to a suitable HGVS data structure. This function is
  useful for sending data through another codec."
  [m]
  (-> m
      (update :kind keyword)
      (update :mutations #(map mut/restore %))))
