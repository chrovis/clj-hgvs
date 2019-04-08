(ns clj-hgvs.coordinate
  "Data structures and functions to handle HGVS coordinates."
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clojure.spec.alpha :as s]
            [clj-hgvs.internal :as intl]))

(declare parser parse-coordinate coding-dna-coordinate rna-coordinate)

(defprotocol Coordinate
  (format [this] "Returns a string representing the given coordinate.")
  (plain [this] "Returns a plain map representing the given coordinate."))

(defmulti restore
  "Restores a plain map to a suitable coordinate record."
  {:arglists '([m])}
  :coordinate)

(defprotocol ICodingDNACoordinate
  (in-exon? [this] "Returns true if the coordinate is located in exon, else false."))

(s/def ::coordinate #(satisfies? Coordinate %))

(s/def ::position (s/and integer? pos?))
(s/def ::offset (s/nilable integer?))
(s/def ::region (s/nilable #{:upstream :downstream}))

(defn ->region-keyword
  [s]
  (if s
    (case s
      "-" :upstream
      "*" :downstream
      "" nil)))

(defn ->region-str
  [k]
  (case k
    :upstream "-"
    :downstream "*"
    nil ""))

;;; unknown coordinate

(defrecord UnknownCoordinate []
  Coordinate
  (format [_] "?")
  (plain [_] {:coordinate "unknown"}))

(s/def ::unknown-coordinate ::coordinate)

(defn unknown-coordinate
  "Returns UnknownCoordinate instance."
  []
  (UnknownCoordinate.))

(defmethod restore "unknown"
  [m]
  (unknown-coordinate))

;;; uncertain coordinate
;;;
;;; e.g. (123456_234567)
;;;      (?_1)
;;;      4072-?

(defrecord UncertainCoordinate [start end]
  Coordinate
  (format [this]
    (let [[known :as knowns] (remove #(instance? UnknownCoordinate %) [start end])]
      (cond
        (= (count knowns) 2) (str "(" (format start) "_" (format end) ")")

        (and (satisfies? ICodingDNACoordinate known) (not (zero? (:offset known))))
        (str (->region-str (:region known))
             (:position known)
             (if (pos? (:offset known)) "+?" "-?"))

        :else (str "(" (format start) "_" (format end) ")"))))
  (plain [this]
    {:coordinate "uncertain"
     :start (plain start)
     :end (plain end)}))

(s/def :clj-hgvs.coordinate.uncertain-coordinate/start ::coordinate)
(s/def :clj-hgvs.coordinate.uncertain-coordinate/end ::coordinate)
(s/def ::uncertain-coordinate
  (s/and ::coordinate
         (s/keys :req-un [:clj-hgvs.coordinate.uncertain-coordinate/start
                          :clj-hgvs.coordinate.uncertain-coordinate/end])))

(defn uncertain-coordinate
  "Returns UncertainCoordinate instance having start and end. Throws an
  exception if any input is illegal."
  [start end]
  {:pre [(not (instance? UncertainCoordinate start))
         (not (instance? UncertainCoordinate end))
         (or (some #(instance? UnknownCoordinate %) [start end])
             (= (type start) (type end)))]
   :post [(intl/valid? ::uncertain-coordinate %)]}
  (UncertainCoordinate. start end))

;; (123456_234567)
(def ^:private uncertain-coordinate-re
  #"\(([\*\+\-\d\?]+)_([\*\+\-\d\?]+)\)")

;; 4072-?
(def ^:private uncertain-coordinate-offset-re
  #"(\-|\*)?(\d+)(\-|\+)\?")

(defn parse-uncertain-coordinate
  [s t]
  (if-let [[_ start end] (re-matches uncertain-coordinate-re s)]
    (let [parse (parser t)]
      (uncertain-coordinate (parse start) (parse end)))
    (if-let [[_ region position offset-direction]
             (re-matches uncertain-coordinate-offset-re s)]
      (let [coordinate (case t
                         :coding-dna coding-dna-coordinate
                         :rna rna-coordinate)]
        (if (= offset-direction "+")
          (uncertain-coordinate (coordinate (intl/parse-long position)
                                            1
                                            (->region-keyword region))
                                (unknown-coordinate))
          (uncertain-coordinate (unknown-coordinate)
                                (coordinate (intl/parse-long position)
                                            -1
                                            (->region-keyword region)))))
      (throw (#?(:clj Exception., :cljs js/Error.)
              (str "Unable to parse string: " s))))))

(defmethod restore "uncertain"
  [m]
  (uncertain-coordinate (restore (:start m)) (restore (:end m))))

(defn- uncertain-coordinate-string?
  [s]
  (or (some? (re-matches uncertain-coordinate-re s))
      (some? (re-matches uncertain-coordinate-offset-re s))))

;;; genomic coordinate

(defrecord GenomicCoordinate [position]
  #?(:clj Comparable, :cljs IComparable)
  (#?(:clj compareTo, :cljs -compare) [this o]
    (compare position (:position o)))
  Coordinate
  (format [_] (str position))
  (plain [this] (into {:coordinate "genomic"} this)))

(s/def ::genomic-coordinate
  (s/and ::coordinate (s/keys :req-un [::position])))

(defn genomic-coordinate
  "Returns GenomicCoordinate instance having position. Throws an exception if
  position is illegal."
  [position]
  {:post [(intl/valid? ::genomic-coordinate %)]}
  (GenomicCoordinate. position))

(defn parse-genomic-coordinate
  "Parses a coordinate string used in genomic mutations, returning a
  GenomicCoordinate or UnknownCoordinate."
  [s]
  (cond
    (uncertain-coordinate-string? s) (parse-uncertain-coordinate s :genomic)
    (= s "?") (unknown-coordinate)
    :else (genomic-coordinate (intl/parse-long s))))

(defmethod restore "genomic"
  [m]
  (genomic-coordinate (:position m)))

;;; mitochondrial coordinate

(defrecord MitochondrialCoordinate [position]
  #?(:clj Comparable, :cljs IComparable)
  (#?(:clj compareTo, :cljs -compare) [this o]
    (compare position (:position o)))
  Coordinate
  (format [_] (str position))
  (plain [this] (into {:coordinate "mitochondrial"} this)))

(s/def ::mitochondrial-coordinate
  (s/and ::coordinate (s/keys :req-un [::position])))

(defn mitochondrial-coordinate
  "Returns MitochondrialCoordinate instance having position. Throws an exception
  if position is illegal."
  [position]
  {:post [(intl/valid? ::mitochondrial-coordinate %)]}
  (MitochondrialCoordinate. position))

(defn parse-mitochondrial-coordinate
  "Parses a coordinate string used in mitochondrial mutations, returning a
  MitochondrialCoordinate or UnknownCoordinate."
  [s]
  (cond
    (uncertain-coordinate-string? s) (parse-uncertain-coordinate s :mitochondrial)
    (= s "?") (unknown-coordinate)
    :else (mitochondrial-coordinate (intl/parse-long s))))

(defmethod restore "mitochondrial"
  [m]
  (mitochondrial-coordinate (:position m)))

;;; coding DNA coordinate
;;;
;;; e.g. 3
;;;      -3
;;;      *3
;;;      87+3
;;;      88-1
;;;      -85+3
;;;      *37+3

(defrecord CodingDNACoordinate [position offset region]
  #?(:clj Comparable, :cljs IComparable)
  (#?(:clj compareTo, :cljs -compare) [this o]
    (let [{o-position :position, o-region :region} o
          offset (or offset 0)
          o-offset (or (:offset o) 0)]
      (if (= region o-region)
        (if (= position o-position)
          (compare offset o-offset)
          (if (= region :upstream)
            (compare o-position position)
            (compare position o-position)))
        (compare (get {:upstream -1, :downstream 1} region 0)
                 (get {:upstream -1, :downstream 1} o-region 0)))))
  Coordinate
  (format [_]
    (str (->region-str region)
         position
         (if-not (or (nil? offset) (zero? offset))
           (str (if (pos? offset) "+") offset))))
  (plain [this]
    (into {:coordinate "coding-dna"} (update this :region #(some-> % name))))
  ICodingDNACoordinate
  (in-exon? [this]
    (or (nil? offset) (zero? offset))))

(s/def ::coding-dna-coordinate
  (s/and ::coordinate (s/keys :req-un [::position ::offset ::region])))

(defn coding-dna-coordinate
  "Returns CodingDNACoordinate instance having position, offset, and region.
  Throws an exception if any input is illegal."
  ([position] (coding-dna-coordinate position 0 nil))
  ([position offset region]
   {:post [(intl/valid? ::coding-dna-coordinate %)]}
   (CodingDNACoordinate. position offset region)))

(def ^:private coding-dna-coordinate-re
  #"^(\-|\*)?(\d+)([\-\+]\d+)?$")

(defn parse-coding-dna-coordinate
  "Parses a coordinate string used in coding DNA mutations, returning a
  CodingDNACoordinate or UnknownCoordinate."
  [s]
  (cond
    (uncertain-coordinate-string? s) (parse-uncertain-coordinate s :coding-dna)
    (= s "?") (unknown-coordinate)
    :else (let [[_ region position offset] (re-find coding-dna-coordinate-re s)]
            (coding-dna-coordinate (intl/parse-long position)
                                   (if (some? offset)
                                     (intl/parse-long offset)
                                     0)
                                   (->region-keyword region)))))

(defmethod restore "coding-dna"
  [m]
  (coding-dna-coordinate (:position m) (:offset m) (keyword (:region m))))

;;; non-coding DNA coordinate

(defrecord NonCodingDNACoordinate [position]
  #?(:clj Comparable, :cljs IComparable)
  (#?(:clj compareTo, :cljs -compare) [this o]
    (compare position (:position o)))
  Coordinate
  (format [_] (str position))
  (plain [this] (into {:coordinate "non-coding-dna"} this)))

(s/def ::non-coding-dna-coordinate
  (s/and ::coordinate (s/keys :req-un [::position])))

(defn non-coding-dna-coordinate
  "Returns NonCodingDNACoordinate instance having position. Throws an exception
  if position is illegal."
  [position]
  {:post [(intl/valid? ::non-coding-dna-coordinate %)]}
  (NonCodingDNACoordinate. position))

(defn parse-non-coding-dna-coordinate
  "Parses a coordinate string used in non-coding DNA mutations, returning a
  NonCodingDNACoordinate or UnknownCoordinate."
  [s]
  (cond
    (uncertain-coordinate-string? s) (parse-uncertain-coordinate s :non-coding-dna)
    (= s "?") (unknown-coordinate)
    :else (non-coding-dna-coordinate (intl/parse-long s))))

(defmethod restore "non-coding-dna"
  [m]
  (non-coding-dna-coordinate (:position m)))

;;; RNA coordinate

(defrecord RNACoordinate [position offset region]
  #?(:clj Comparable, :cljs IComparable)
  (#?(:clj compareTo, :cljs -compare) [this o]
    (let [{o-position :position, o-region :region} o
          offset (or offset 0)
          o-offset (or (:offset o) 0)]
      (if (= region o-region)
        (if (= position o-position)
          (compare offset o-offset)
          (if (= region :upstream)
            (compare o-position position)
            (compare position o-position)))
        (compare (get {:upstream -1, :downstream 1} region 0)
                 (get {:upstream -1, :downstream 1} o-region 0)))))
  Coordinate
  (format [_]
    (str (->region-str region)
         position
         (if-not (or (nil? offset) (zero? offset))
           (str (if (pos? offset) "+") offset))))
  (plain [this]
    (into {:coordinate "rna"} (update this :region #(some-> % name))))
  ICodingDNACoordinate
  (in-exon? [this]
    (or (nil? offset) (zero? offset))))

(s/def ::rna-coordinate
  (s/and ::coordinate (s/keys :req-un [::position ::offset ::region])))

(defn rna-coordinate
  "Returns RNACoordinate instance having position, offset, and region. Throws an
  exception if any input is illegal."
  [position offset region]
  {:post [(intl/valid? ::rna-coordinate %)]}
  (RNACoordinate. position offset region))

(def ^:private rna-coordinate-re
  #"^(\-|\*)?(\d+)([\-\+]\d+)?$")

(defn parse-rna-coordinate
  "Parses a coordinate string used in RNA mutations, returning a RNACoordinate
  or UnknownCoordinate."
  [s]
  (if (= s "?")
    (unknown-coordinate)
    (let [[_ region position offset] (re-find rna-coordinate-re s)]
      (if (or (some? region) (some? offset))
        (rna-coordinate (intl/parse-long position)
                        (if (some? offset)
                          (intl/parse-long offset)
                          0)
                        (->region-keyword region))
        (rna-coordinate (intl/parse-long position) nil nil)))))

(defmethod restore "rna"
  [m]
  (rna-coordinate (:position m) (:offset m) (keyword (:region m))))

;;; protein coordinate

(defrecord ProteinCoordinate [position]
  #?(:clj Comparable, :cljs IComparable)
  (#?(:clj compareTo, :cljs -compare) [this o]
    (compare position (:position o)))
  Coordinate
  (format [_] (str position))
  (plain [this] (into {:coordinate "protein"} this)))

(s/def ::protein-coordinate
  (s/and ::coordinate (s/keys :req-un [::position])))

(defn protein-coordinate
  "Returns ProteinCoordinate instance having position. Throws an exception if
  position is illegal."
  [position]
  {:post [(intl/valid? ::protein-coordinate %)]}
  (ProteinCoordinate. position))

(defn parse-protein-coordinate
  "Parses a coordinate string used in protein mutations, returning a
  ProteinCoordinate or UnknownCoordinate."
  [s]
  (if (= s "?")
    (unknown-coordinate)
    (protein-coordinate (intl/parse-long s))))

(defmethod restore "protein"
  [m]
  (protein-coordinate (:position m)))

;;; general parser

(defn- parser
  [t]
  (case t
    :genomic parse-genomic-coordinate
    :mitochondrial parse-mitochondrial-coordinate
    :coding-dna parse-coding-dna-coordinate
    :non-coding-dna parse-non-coding-dna-coordinate
    :rna parse-rna-coordinate
    :protein parse-protein-coordinate))

(defn parse-coordinate
  [s t]
  ((parser t) s))

(defn comparable-coordinates?
  "Returns true if the two coordinates are comparable, else false."
  [coord1 coord2]
  (not-any? (some-fn #(instance? UnknownCoordinate %)
                     #(instance? UncertainCoordinate %))
            [coord1 coord2]))
