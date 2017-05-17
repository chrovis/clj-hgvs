(ns clj-hgvs.coordinate
  "Data structures and functions to handle HGVS coordinates."
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clj-hgvs.internal :refer [parse-long]]))

(declare parser parse-coordinate cdna-coordinate rna-coordinate)

(defprotocol Coordinate
  (format [this] "Returns a string representing the given coordinate.")
  (plain [this] "Returns a plain map representing the given coordinate."))

(defmulti restore
  "Restores a plain map to a suitable coordinate record."
  {:arglists '([m])}
  :coordinate)

(defprotocol ICDNACoordinate
  (in-exon? [this] "Returns true if the coordinate is located in exon, else false."))

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

(defrecord UncertainCoordinate [start end]
  Coordinate
  (format [this] (str "(" (format start) "_" (format end) ")"))
  (plain [this]
    {:coordinate "uncertain"
     :start (plain start)
     :end (plain end)}))

(defn uncertain-coordinate
  "Returns UncertainCoordinate instance having start and end. Throws an
  exception if any input is illegal."
  [start end]
  {:pre [(satisfies? Coordinate start)
         (satisfies? Coordinate end)
         (not (instance? UncertainCoordinate start))
         (not (instance? UncertainCoordinate end))
         (or (some #(instance? UnknownCoordinate %) [start end])
             (= (type start) (type end)))]}
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
                         :cdna cdna-coordinate
                         :rna rna-coordinate)]
        (if (= offset-direction "+")
          (uncertain-coordinate (coordinate (parse-long position)
                                            1
                                            (->region-keyword region))
                                (unknown-coordinate))
          (uncertain-coordinate (unknown-coordinate)
                                (coordinate (parse-long position)
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

(defn genomic-coordinate
  "Returns GenomicCoordinate instance having position. Throws an exception if
  position is illegal."
  [position]
  {:pre [(integer? position) (pos? position)]}
  (GenomicCoordinate. position))

(defn parse-genomic-coordinate
  "Parses a coordinate string used in genomic mutations, returning a
  GenomicCoordinate or UnknownCoordinate."
  [s]
  (cond
    (uncertain-coordinate-string? s) (parse-uncertain-coordinate s :genomic)
    (= s "?") (unknown-coordinate)
    :else (genomic-coordinate (parse-long s))))

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

(defn mitochondrial-coordinate
  "Returns MitochondrialCoordinate instance having position. Throws an exception
  if position is illegal."
  [position]
  {:pre [(integer? position) (pos? position)]}
  (MitochondrialCoordinate. position))

(defn parse-mitochondrial-coordinate
  "Parses a coordinate string used in mitochondrial mutations, returning a
  MitochondrialCoordinate or UnknownCoordinate."
  [s]
  (cond
    (uncertain-coordinate-string? s) (parse-uncertain-coordinate s :mitochondrial)
    (= s "?") (unknown-coordinate)
    :else (mitochondrial-coordinate (parse-long s))))

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

(defrecord CDNACoordinate [position offset region]
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
    (into {:coordinate "cdna"} (update this :region #(some-> % name))))
  ICDNACoordinate
  (in-exon? [this]
    (or (nil? offset) (zero? offset))))

(defn cdna-coordinate
  "Returns CDNACoordinate instance having position, offset, and region. Throws
  an exception if any input is illegal."
  ([position] (cdna-coordinate position 0 nil))
  ([position offset region]
   {:pre [(integer? position) (pos? position)
          (integer? offset)
          (or (nil? region) (#{:upstream :downstream} region))]}
   (CDNACoordinate. position offset region)))

(def ^:private cdna-coordinate-re
  #"^(\-|\*)?(\d+)([\-\+]\d+)?$")

(defn parse-cdna-coordinate
  "Parses a coordinate string used in coding DNA mutations, returning a
  CDNACoordinate or UnknownCoordinate."
  [s]
  (cond
    (uncertain-coordinate-string? s) (parse-uncertain-coordinate s :cdna)
    (= s "?") (unknown-coordinate)
    :else (let [[_ region position offset] (re-find cdna-coordinate-re s)]
            (cdna-coordinate (parse-long position)
                             (if (some? offset)
                               (parse-long offset)
                               0)
                             (->region-keyword region)))))

(defmethod restore "cdna"
  [m]
  (cdna-coordinate (:position m) (:offset m) (keyword (:region m))))

;;; non-coding DNA coordinate

(defrecord NCDNACoordinate [position]
  #?(:clj Comparable, :cljs IComparable)
  (#?(:clj compareTo, :cljs -compare) [this o]
    (compare position (:position o)))
  Coordinate
  (format [_] (str position))
  (plain [this] (into {:coordinate "ncdna"} this)))

(defn ncdna-coordinate
  "Returns NCDNACoordinate instance having position. Throws an exception if
  position is illegal."
  [position]
  {:pre [(integer? position) (pos? position)]}
  (NCDNACoordinate. position))

(defn parse-ncdna-coordinate
  "Parses a coordinate string used in non-coding DNA mutations, returning a
  NCDNACoordinate or UnknownCoordinate."
  [s]
  (cond
    (uncertain-coordinate-string? s) (parse-uncertain-coordinate s :ncdna)
    (= s "?") (unknown-coordinate)
    :else (ncdna-coordinate (parse-long s))))

(defmethod restore "ncdna"
  [m]
  (ncdna-coordinate (:position m)))

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
  ICDNACoordinate
  (in-exon? [this]
    (or (nil? offset) (zero? offset))))

(defn rna-coordinate
  "Returns RNACoordinate instance having position, offset, and region. Throws an
  exception if any input is illegal."
  [position offset region]
  {:pre [(integer? position) (pos? position)
         (or (nil? offset) (integer? offset))
         (or (nil? region) (#{:upstream :downstream} region))]}
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
        (rna-coordinate (parse-long position)
                        (if (some? offset)
                          (parse-long offset)
                          0)
                        (->region-keyword region))
        (rna-coordinate (parse-long position) nil nil)))))

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

(defn protein-coordinate
  "Returns ProteinCoordinate instance having position. Throws an exception if
  position is illegal."
  [position]
  {:pre [(integer? position) (pos? position)]}
  (ProteinCoordinate. position))

(defn parse-protein-coordinate
  "Parses a coordinate string used in protein mutations, returning a
  ProteinCoordinate or UnknownCoordinate."
  [s]
  (if (= s "?")
    (unknown-coordinate)
    (protein-coordinate (parse-long s))))

(defmethod restore "protein"
  [m]
  (protein-coordinate (:position m)))

;;; general parser

(defn- parser
  [t]
  (case t
    :genomic parse-genomic-coordinate
    :mitochondrial parse-mitochondrial-coordinate
    :cdna parse-cdna-coordinate
    :ncdna parse-ncdna-coordinate
    :rna parse-rna-coordinate
    :protein parse-protein-coordinate))

(defn parse-coordinate
  [s t]
  ((parser t) s))
