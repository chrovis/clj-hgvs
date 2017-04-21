(ns clj-hgvs.coordinate
  "Data structures and functions to handle HGVS coordinates."
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clj-hgvs.internal :refer [parse-long]]))

(defprotocol Coordinate
  (format [this]))

(defprotocol ICDNACoordinate
  (in-exon? [this]))

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
  [position]
  {:pre [(integer? position) (pos? position)]}
  (GenomicCoordinate. position))

(defn parse-genomic-coordinate
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
  [position]
  {:pre [(integer? position) (pos? position)]}
  (MitochondrialCoordinate. position))

(defn parse-mitochondrial-coordinate
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
  ([position] (cdna-coordinate position 0 nil))
  ([position offset region]
   {:pre [(integer? position) (pos? position)
          (integer? offset)
          (or (nil? region) (#{:upstream :downstream} region))]}
   (CDNACoordinate. position offset region)))

(def ^:private cdna-coordinate-re
  #"^(\-|\*)?(\d+)([\-\+]\d+)?$")

(defn parse-cdna-coordinate
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
  [position]
  {:pre [(integer? position) (pos? position)]}
  (NCDNACoordinate. position))

(defn parse-ncdna-coordinate
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
  [position offset region]
  {:pre [(integer? position) (pos? position)
         (or (nil? offset) (integer? offset))
         (or (nil? region) (#{:upstream :downstream} region))]}
  (RNACoordinate. position offset region))

(def ^:private rna-coordinate-re
  #"^(\-|\*)?(\d+)([\-\+]\d+)?$")

(defn parse-rna-coordinate
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
  [position]
  {:pre [(integer? position) (pos? position)]}
  (ProteinCoordinate. position))

(defn parse-protein-coordinate
  [s]
  (if (= s "?")
    (unknown-coordinate)
    (protein-coordinate (parse-long s))))
