# clj-hgvs

HGVS library for Clojure(Script).

## Installation

```clojure
[clj-hgvs "0.1.0-SNAPSHOT"]
```

## Usage

```clojure
(require '[clj-hgvs.core :as hgvs])

;; `hgvs` is a constructor of HGVS map.
(def hgvs1 (hgvs/hgvs "NM_005228.3" :coding-dna {:numbering "2573"
                                                 :type :substitution
                                                 :ref "T"
                                                 :alt "G"}))

hgvs1
=> {:transcript "NM_005228.3"
    :kind :coding-dna
    :mutations ({:numbering "2573", :type :substitution, :ref "T", :alt "G"})}

;; `format` returns HGVS text.
(hgvs/format hgvs1)
=> "NM_005228.3:c.2573T>G"

;; `parse` parses HGVS text, returning HGVS map.
(hgvs/parse "NM_005228.3:c.2573T>G")
=> {:transcript "NM_005228.3"
    :kind :coding-dna
    :mutations ({:numbering "2573", :type :substitution, :ref "T", :alt "G"})}
```
