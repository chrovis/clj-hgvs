(ns clj-hgvs.coordinate
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clj-hgvs.internal :refer [parse-long]]))

(defprotocol Coordinate
  (format [this]))

(defprotocol ICDNACoordinate
  (in-exon? [this]))

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

(defrecord CDNACoordinate [position offset region]
  Coordinate
  (format [_]
    (str (case region
           :upstream "-"
           :downstream "*"
           nil)
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
                       (case region
                         "-" :upstream
                         "*" :downstream
                         nil)))))

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

(defrecord RNACoordinate [position offset region]
  Coordinate
  (format [_]
    (str (case region
           :upstream "-"
           :downstream "*"
           nil)
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
  (let [[_ region position offset] (re-find rna-coordinate-re s)]
    (if (or (some? region) (some? offset))
      (rna-coordinate (parse-long position)
                      (if (some? offset)
                        (parse-long offset)
                        0)
                      (case region
                        "-" :upstream
                        "*" :downstream
                        nil))
      (rna-coordinate (parse-long position) nil nil))))

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
