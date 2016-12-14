# clj-hgvs

HGVS library for Clojure(Script).

## Installation

```clojure
[clj-hgvs "0.1.0-SNAPSHOT"]
```

## Usage

```clojure
(require '[clj-hgvs.core :as hgvs])

;; `parse` parses HGVS text, returning HGVS map.
(def hgvs1 (hgvs/parse "NM_005228.3:c.2573T>G"))

hgvs1
=> {:transcript "NM_005228.3"
    :kind :cdna
    :mutations [#clj_hgvs.mutation.DNASubstitution
                {:coord-start #clj_hgvs.coordinate.CDNACoordinate {:position 2573
                                                                   :stream nil
                                                                   :intron-offset nil}
                 :coord-end nil
                 :ref "T"
                 :type ">"
                 :alt "G"}]}

;; `format` returns HGVS text.
(hgvs/format hgvs1)
=> "NM_005228.3:c.2573T>G"
```
