(ns clj-hgvs.coordinate
  #?(:clj (:refer-clojure :exclude [format]))
  (:require [clj-hgvs.internal :refer [parse-long]]))

(defprotocol Coordinate
  (format [this]))

(defrecord UnknownCoordinate []
  Coordinate
  (format [_] "?"))

;;; genomic coordinate

(defrecord GenomicCoordinate [position]
  Coordinate
  (format [_]
    (str position)))

(defn parse-genomic-coordinate
  [s]
  (if (= s "?")
    (UnknownCoordinate.)
    (->GenomicCoordinate (parse-long s))))

;;; mitochondrial coordinate

(defrecord MitochondrialCoordinate [position]
  Coordinate
  (format [_]
    (str position)))

(defn parse-mitochondrial-coordinate
  [s]
  (if (= s "?")
    (UnknownCoordinate.)
   (->MitochondrialCoordinate (parse-long s))))

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
         (if (and intron-offset (pos? intron-offset))
           "+")
         intron-offset)))

(def ^:private cdna-coordinate-re
  #"^(\-|\*)?(\d+)([\-\+]\d+)?$")

(defn parse-cdna-coordinate
  [s]
  (if (= s "?")
    (UnknownCoordinate.)
    (let [[_ stream position intron-offset] (re-find cdna-coordinate-re s)]
      (map->CDNACoordinate {:position (parse-long position)
                            :stream (case stream
                                      "-" :up
                                      "*" :down
                                      nil)
                            :intron-offset (some-> intron-offset parse-long)}))))

;;; non-coding DNA coordinate

(defrecord NCDNACoordinate [position]
  Coordinate
  (format [_]
    (str position)))

(defn parse-ncdna-coordinate
  [s]
  (if (= s "?")
    (UnknownCoordinate.)
    (->NCDNACoordinate (parse-long s))))

;;; TODO: RNA coordinate

(defrecord RNACoordinate [])

(defn parse-rna-coordinate
  [s]
)

;;; protein coordinate

(defrecord ProteinCoordinate [position]
  Coordinate
  (format [_]
    (str position)))

(defn parse-protein-coordinate
  [s]
  (if (= s "?")
    (UnknownCoordinate.)
    (->ProteinCoordinate (parse-long s))))
