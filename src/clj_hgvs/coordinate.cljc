(ns clj-hgvs.coordinate
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clj-hgvs.internal :refer [parse-long]]))

(defprotocol Coordinate
  (format [this]))

;;; genomic coordinate

(defrecord GenomicCoordinate [position]
  Coordinate
  (format [_]
    (str position)))

(defn parse-genomic-coordinate
  [s]
  (->GenomicCoordinate (parse-long s)))

;;; mitochondrial coordinate

(defrecord MitochondrialCoordinate [position]
  Coordinate
  (format [_]
    (str position)))

(defn parse-mitochondrial-coordinate
  [s]
  (->MitochondrialCoordinate (parse-long s)))

;;; TODO: coding DNA coordinate

(defrecord CDNACoordinate [])

;;; non-coding DNA coordinate

(defrecord NCDNACoordinate [position]
  Coordinate
  (format [_]
    (str position)))

(defn parse-ncdna-coordinate
  [s]
  (->NCDNACoordinate (parse-long s)))

;;; TODO: RNA coordinate

(defrecord RNACoordinate [])

;;; protein coordinate

(defrecord ProteinCoordinate [position]
  Coordinate
  (format [_]
    (str position)))

(defn parse-protein-coordinate
  [s]
  (->ProteinCoordinate (parse-long s)))
