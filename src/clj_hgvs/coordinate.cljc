(ns clj-hgvs.coordinate
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clj-hgvs.internal :refer [parse-long]]))

(defprotocol Coordinate
  (format [this]))

(defrecord UnknownCoordinate []
  Coordinate
  (format [_] "?"))

(defn unknown-coordinate
  []
  (UnknownCoordinate.))

;;; genomic coordinate

(defrecord GenomicCoordinate [position]
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

(defrecord CDNACoordinate [position stream intron-offset]
  Coordinate
  (format [_]
    (str (case stream
           :up "-"
           :down "*"
           nil)
         position
         (if-not (or (nil? intron-offset) (zero? intron-offset))
           (str (if (pos? intron-offset) "+") intron-offset)))))

(defn cdna-coordinate
  [position stream intron-offset]
  {:pre [(integer? position) (pos? position)
         (or (nil? stream) (#{:up :down} stream))
         (or (nil? intron-offset) (integer? intron-offset))]}
  (CDNACoordinate. position stream intron-offset))

(def ^:private cdna-coordinate-re
  #"^(\-|\*)?(\d+)([\-\+]\d+)?$")

(defn parse-cdna-coordinate
  [s]
  (if (= s "?")
    (unknown-coordinate)
    (let [[_ stream position intron-offset] (re-find cdna-coordinate-re s)]
      (cdna-coordinate (parse-long position)
                       (case stream
                         "-" :up
                         "*" :down
                         nil)
                       (some-> intron-offset parse-long)))))

;;; non-coding DNA coordinate

(defrecord NCDNACoordinate [position]
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

(defrecord RNACoordinate [position stream intron-offset]
  Coordinate
  (format [_]
    (str (case stream
           :up "-"
           :down "*"
           nil)
         position
         (if (and intron-offset (pos? intron-offset))
           "+")
         intron-offset)))

(defn rna-coordinate
  [position stream intron-offset]
  {:pre [(integer? position) (pos? position)
         (or (nil? stream) (#{:up :down} stream))
         (or (nil? intron-offset) (integer? intron-offset))]}
  (RNACoordinate. position stream intron-offset))

(def ^:private rna-coordinate-re
  #"^(\-|\*)?(\d+)([\-\+]\d+)?$")

(defn parse-rna-coordinate
  [s]
  (let [[_ stream position intron-offset] (re-find rna-coordinate-re s)]
    (rna-coordinate (parse-long position)
                    (case stream
                      "-" :up
                      "*" :down
                      nil)
                    (some-> intron-offset parse-long))))

;;; protein coordinate

(defrecord ProteinCoordinate [position]
  Coordinate
  (format [_]
    (str position)))

;; TODO
(defn protein-coordinate
  [position]
  {:pre [(integer? position) (pos? position)]}
  (ProteinCoordinate. position))

(defn parse-protein-coordinate
  [s]
  (if (= s "?")
    (unknown-coordinate)
    (protein-coordinate (parse-long s))))
