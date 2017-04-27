(ns clj-hgvs.coordinate
  "Data structures and functions to handle HGVS coordinates."
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clj-hgvs.internal :refer [parse-long]]))

(defprotocol Coordinate
  (format [this] "Returns a string representing the given coordinate."))

(defprotocol ICDNACoordinate
  (in-exon? [this] "Returns true if the coordinate is located in exon, else false."))

(defn ->region-keyword
  [s]
  (case s
    "-" :upstream
    "*" :downstream
    "" nil
    nil nil))

(defn ->region-str
  [k]
  (case k
    :upstream "-"
    :downstream "*"
    nil ""))

;;; unknown coordinate

(defrecord UnknownCoordinate []
  Coordinate
  (format [_] "?"))

(defn unknown-coordinate
  "Returns UnknownCoordinate instance."
  []
  (UnknownCoordinate.))

;;; genomic coordinate

(defrecord GenomicCoordinate [position]
  #?(:clj Comparable, :cljs IComparable)
  (#?(:clj compareTo, :cljs -compare) [this o]
    (compare position (:position o)))
  Coordinate
  (format [_]
    (str position)))

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
  (if (= s "?")
    (unknown-coordinate)
    (genomic-coordinate (parse-long s))))

;;; mitochondrial coordinate

(defrecord MitochondrialCoordinate [position]
  #?(:clj Comparable, :cljs IComparable)
  (#?(:clj compareTo, :cljs -compare) [this o]
    (compare position (:position o)))
  Coordinate
  (format [_]
    (str position)))

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
  (if (= s "?")
    (unknown-coordinate)
    (mitochondrial-coordinate (parse-long s))))

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
  (if (= s "?")
    (unknown-coordinate)
    (let [[_ region position offset] (re-find cdna-coordinate-re s)]
      (cdna-coordinate (parse-long position)
                       (if (some? offset)
                         (parse-long offset)
                         0)
                       (->region-keyword region)))))

;;; non-coding DNA coordinate

(defrecord NCDNACoordinate [position]
  #?(:clj Comparable, :cljs IComparable)
  (#?(:clj compareTo, :cljs -compare) [this o]
    (compare position (:position o)))
  Coordinate
  (format [_]
    (str position)))

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
  (if (= s "?")
    (unknown-coordinate)
    (ncdna-coordinate (parse-long s))))

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

;;; protein coordinate

(defrecord ProteinCoordinate [position]
  #?(:clj Comparable, :cljs IComparable)
  (#?(:clj compareTo, :cljs -compare) [this o]
    (compare position (:position o)))
  Coordinate
  (format [_]
    (str position)))

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
