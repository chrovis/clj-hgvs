# clj-hgvs

Clojure(Script) library for handling [HGVS](http://varnomen.hgvs.org/).

## Features

clj-hgvs provides:

* Data structure of HGVS
* Parser of HGVS text
* Formatter to HGVS text

## Installation

With Leiningen/Boot:

```clojure
[clj-hgvs "0.1.0-SNAPSHOT"]
```

## Usage

```clojure
(require '[clj-hgvs.core :as hgvs])

;; `parse` parses HGVS text, returning HGVS map.
(def hgvs1 (hgvs/parse "NM_005228.3:c.2573T>G"))

hgvs1
;;=> {:transcript "NM_005228.3"
;;    :kind :cdna
;;    :mutations [#clj_hgvs.mutation.DNASubstitution
;;                {:coord #clj_hgvs.coordinate.CDNACoordinate {:position 2573
;;                                                             :offset 0
;;                                                             :region nil}
;;                 :ref "T"
;;                 :type ">"
;;                 :alt "G"}]}

;; `format` returns HGVS text.
(hgvs/format hgvs1)
;;=> "NM_005228.3:c.2573T>G"
```

## License

Copyright 2017 [Xcoo, Inc.](https://xcoo.jp/)

Licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).
