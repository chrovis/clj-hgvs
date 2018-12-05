(ns clj-hgvs.core
  "Functions to handle HGVS. See http://varnomen.hgvs.org/ for the detail HGVS
  nomenclature."
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.spec.alpha :as s]
            [clj-hgvs.internal :as intl]
            [clj-hgvs.mutation :as mut]))

(defmacro with-validation-disabled
  "Disables validation within a scope."
  [& body]
  `(binding [intl/*validation-enabled* nil]
     ~@body))

(s/def ::transcript
  (s/and string? (s/or :N*_ #(re-matches #"N(C|G|M|R|P)_\d+(\.\d+)?" %)
                       :LRG_ #(re-matches #"LRG_\d+((t|p)\d+)?" %))))

(s/def ::kind #{:genome :mitochondria :cdna :ncdna :rna :protein})

(s/def :clj-hgvs.hgvs/transcript (s/nilable ::transcript))
(s/def ::hgvs (s/keys :req-un [:clj-hgvs.hgvs/transcript
                               ::kind
                               ::mut/mutation]))

(defn hgvs
  "Constructor of HGVS map. Throws an exception if any input is illegal."
  [transcript kind mutation]
  {:post [(intl/valid? ::hgvs %)]}
  {:transcript transcript
   :kind kind
   :mutation (if (string? mutation)
               (mut/parse mutation kind)
               mutation)})

(s/fdef hgvs
  :args (s/cat :transcript (s/nilable ::transcript)
               :kind       ::kind
               :mutation   (s/or :s string?
                                 :m ::mut/mutation))
  :ret  ::hgvs)

(def ^:private hgvs-re #"^(?:([^:]+):)?([gmcnrp])\.(.+)$")

(defn parse
  "Parses a HGVS string s, returning a map representing the HGVS."
  [s]
  (let [[_ transcript kind mutation] (re-find hgvs-re s)
        kind-k (intl/->kind-keyword kind)]
    (hgvs transcript kind-k mutation)))

(s/fdef parse
  :args (s/cat :s string?)
  :ret  ::hgvs)

(defn- format-transcript
  [transcript]
  (if transcript
    (str transcript ":")))

(defn- format-kind
  [kind]
  (str (intl/->kind-str kind) "."))

(defn format
  "Returns a HGVS string representing the given HGVS map. The second argument is
  an optional map to specify style:
    :show-bases? - displays additional bases, e.g. g.6_8delTGC, default false.
    :ins-format - bases style of insertion, default :auto. <:auto|:bases|:count>
    :range-format - range style, default :auto. <:auto|:bases|:coord>
    :amino-acid-format - amino acid style of protein HGVS, default :long.
                         <:long|:short>
    :show-ter-site? - displays a new ter codon site of protein frame shift,
                      default false.
    :ter-format - ter codon style of protein frame shift and extension, default
                  :long. <:long|:short>"
  ([hgvs] (format hgvs {}))
  ([hgvs opts]
   (apply str [(format-transcript (:transcript hgvs))
               (format-kind (:kind hgvs))
               (mut/format (:mutation hgvs) opts)])))

(s/def :clj-hgvs.spec.format-options/show-bases? #(or (true? %) (false? %)))
(s/def :clj-hgvs.spec.format-options/ins-format #{:auto :bases :count})
(s/def :clj-hgvs.spec.format-options/range-format #{:auto :bases :coord})
(s/def :clj-hgvs.spec.format-options/amino-acid-format #{:long :short})
(s/def :clj-hgvs.spec.format-options/show-ter-site? #(or (true? %) (false? %)))
(s/def :clj-hgvs.spec.format-options/ter-format #{:long :short})
(s/def ::format-options
  (s/keys :opt-un [:clj-hgvs.spec.format-options/show-bases?
                   :clj-hgvs.spec.format-options/ins-format
                   :clj-hgvs.spec.format-options/range-format
                   :clj-hgvs.spec.format-options/amino-acid-format
                   :clj-hgvs.spec.format-options/show-ter-site?
                   :clj-hgvs.spec.format-options/ter-format]))

(s/fdef format
  :args (s/cat :hgvs ::hgvs
               :opts (s/? (s/nilable ::format-options)))
  :ret  string?)

(defn plain
  "Returns a plain map representing the given HGVS. This function is useful for
  sending data through another codec. Use restore to retrieve original HGVS
  data."
  [hgvs]
  (-> hgvs
      (update :kind name)
      (update :mutation mut/plain)))

(s/def :clj-hgvs.plain-hgvs/transcript (s/nilable ::transcript))
(s/def :clj-hgvs.plain-hgvs/kind #{"genome" "mitochondria" "cdna" "ncdna" "rna" "protein"})
(s/def :clj-hgvs.plain-hgvs/mutation ::mut/plain-mutation)
(s/def ::plain-hgvs (s/keys :req-un [:clj-hgvs.plain-hgvs/transcript
                                     :clj-hgvs.plain-hgvs/kind
                                     :clj-hgvs.plain-hgvs/mutation]))

(s/fdef plain
  :args (s/cat :hgvs ::hgvs)
  :ret  ::plain-hgvs)

(defn restore
  "Restores a plain map to a suitable HGVS data structure. This function is
  useful for sending data through another codec."
  [m]
  (-> m
      (update :kind keyword)
      (update :mutation mut/restore)))

(s/fdef restore
  :args (s/cat :m ::plain-hgvs)
  :ret  ::hgvs)

(defn normalize
  "Reformats the HGVS string s, returning the normalized HGVS string. Default
  values are used for style options. Throws an exception when an illegal HGVS
  string is supplied."
  [s]
  (format (parse s)))

(s/fdef clj-hgvs.core/normalize
  :args (s/cat :s string?)
  :ret  string?)
