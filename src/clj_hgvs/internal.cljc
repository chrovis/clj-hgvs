(ns clj-hgvs.internal)

(defn parse-long
  [s]
  #?(:clj (Long/parseLong s)
     :cljs (js/parseInt s)))
