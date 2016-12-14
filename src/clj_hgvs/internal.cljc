(ns clj-hgvs.internal)

(defn parse-long
  [s]
  #?(:clj (Long/parseLong s)
     :cljs (js/parseInt s)))

(defn ->kind-keyword
  [s]
  (case s
    "g" :genome
    "m" :mitochondria
    "c" :cdna
    "n" :ncdna
    "r" :rna
    "p" :protein))

(defn ->kind-str
  [k]
  (case k
    :genome "g"
    :mitochondria "m"
    :cdna "c"
    :ncdna "n"
    :rna "r"
    :protein"p"))
