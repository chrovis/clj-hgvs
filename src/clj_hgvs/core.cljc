(ns clj-hgvs.core
  "Functions to handle HGVS. See http://varnomen.hgvs.org/ for the detail HGVS
  nomenclature."
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clj-hgvs.internal :refer [->kind-keyword ->kind-str]]
            [clj-hgvs.mutation :as mut]))

(defn- transcript?
  [s]
  (or (some? (re-matches #"N(C|G|M|R|P)_\d+(\.\d+)?" s))
      (some? (re-matches #"LRG_\d+((t|p)\d+)?" s))))

(defn- kind?
  [k]
  (some? (#{:genome :mitochondria :cdna :ncdna :rna :protein} k)))

(defn- mutation-parser
  [kind]
  (cond
    (#{:genome :mitochondria :cdna :ncdna} kind) #(mut/parse-dna % kind)
    (= kind :rna) mut/parse-rna
    (= kind :protein) mut/parse-protein))

(defn hgvs
  "Constructor of HGVS map. Throws an exception if any input is illegal."
  [transcript kind mutation]
  {:pre [(or (nil? transcript) (transcript? transcript))
         (kind? kind)
         (or (satisfies? mut/Mutation mutation) (string? mutation))]}
  {:transcript transcript
   :kind kind
   :mutation (if (string? mutation)
               ((mutation-parser kind) mutation)
               mutation)})

(def ^:private hgvs-re #"^(?:([^:]+):)?([gmcnrp])\.(.+)$")

(defn parse
  "Parses a HGVS string s, returning a map representing the HGVS."
  [s]
  (let [[_ transcript kind mutation] (re-find hgvs-re s)
        kind-k (->kind-keyword kind)]
    (hgvs transcript kind-k mutation)))

(defn- format-transcript
  [transcript]
  (if transcript
    (str transcript ":")))

(defn- format-kind
  [kind]
  (str (->kind-str kind) "."))

(defn format
  "Returns a HGVS string representing the given HGVS map. The second argument is
  an optional map to specify style:
    :show-bases? - displays additional bases, e.g. g.6_8delTGC, default false.
    :ins-format - bases style of insertion, default :auto. <:auto|:bases|:count>
    :range-format - range style, default :auto. <:auto|:bases|:coord>
    :amino-acid-format - amino acid style of protein HGVS, default :long. <:long|:short>
    :ter-format - ter codon style of protein frame shift. <:long|:short>"
  ([hgvs] (format hgvs nil))
  ([hgvs opts]
   (apply str [(format-transcript (:transcript hgvs))
               (format-kind (:kind hgvs))
               (mut/format (:mutation hgvs) opts)])))

(defn plain
  "Returns a plain map representing the given HGVS. This function is useful for
  sending data through another codec. Use restore to retrieve original HGVS
  data."
  [hgvs]
  (-> hgvs
      (update :kind name)
      (update :mutation mut/plain)))

(defn restore
  "Restores a plain map to a suitable HGVS data structure. This function is
  useful for sending data through another codec."
  [m]
  (-> m
      (update :kind keyword)
      (update :mutation mut/restore)))
