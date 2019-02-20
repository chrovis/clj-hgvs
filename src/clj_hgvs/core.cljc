(ns clj-hgvs.core
  "Main functions for handling HGVS. See http://varnomen.hgvs.org/ for the
  detail HGVS nomenclature."
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.spec.alpha :as s]
            [clj-hgvs.internal :as intl]
            [clj-hgvs.mutation :as mut]))

(defmacro with-validation-disabled
  "Disables validation within a scope.

  HGVS data are checked upon interpretation phases by default, but in the scope,
  all validations will be skipped. This macro may improve the performance on
  handling HGVS which validity is already known."
  [& body]
  `(binding [intl/*validation-enabled* nil]
     ~@body))

(defrecord HGVS [transcript kind mutation])

(s/def ::transcript
  (s/and string? (s/or :N*_ #(re-matches #"N(C|G|M|R|P)_\d+(\.\d+)?" %)
                       :LRG_ #(re-matches #"LRG_\d+((t|p)\d+)?" %))))

(s/def ::kind #{:genome :mitochondria :cdna :ncdna :rna :protein})

(s/def :clj-hgvs.hgvs/transcript (s/nilable ::transcript))
(s/def ::hgvs (s/keys :req-un [:clj-hgvs.hgvs/transcript
                               ::kind
                               ::mut/mutation]))

(defn hgvs
  "Creates HGVS data represented as a record.

  transcript is nilable for conventional reasons, but you should supply it if
  possible.

  kind must be one of :genome, :mitochondria, :cdna, :ncdna, :rna, and :protein.

  mutation must be a clj-hgvs.mutation record or string. The string mutation
  will be parsed by clj-hgvs.mutation/parse."
  [transcript kind mutation]
  {:post [(intl/valid? ::hgvs %)]}
  (map->HGVS {:transcript transcript
              :kind kind
              :mutation (if (string? mutation)
                          (mut/parse mutation kind)
                          mutation)}))

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
  "Returns a HGVS string representing the given HGVS map.

  The second argument is an optional map to specify style:

    {:show-bases?        Display additional bases, e.g. g.6_8delTGC, default
                         false.

     :ins-format         Bases style of insertion, default :auto.
                         <:auto|:bases|:count>

     :range-format       Range style, default :auto. <:auto|:bases|:coord>

     :amino-acid-format  Amino acid style of protein HGVS, default :long.
                         <:long|:short>

     :show-ter-site?     Display a new ter codon site of protein frame shift,
                         default false.

     :ter-format         Ter codon style of protein frame shift and extension,
                         default :long. <:long|:short>}"
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
  sending data through another codec. Use clj-hgvs.core/restore to retrieve
  original HGVS data."
  [hgvs]
  {:transcript (:transcript hgvs)
   :kind (name (:kind hgvs))
   :mutation (mut/plain (:mutation hgvs))})

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
      (update :mutation mut/restore)
      map->HGVS))

(s/fdef restore
  :args (s/cat :m ::plain-hgvs)
  :ret  ::hgvs)

(defn normalize
  "Reformats the HGVS string s, returning the normalized HGVS string.

  The default style will be used for reformatting. See clj-hgvs.core/format
  document for further details of the style options."
  [s]
  (format (parse s)))

(s/fdef clj-hgvs.core/normalize
  :args (s/cat :s string?)
  :ret  string?)

#?(:clj (def data-readers
          "Tagged literal support if loader does not find \"data_readers.clj\"."
          {'clj-hgvs/hgvs parse}))
