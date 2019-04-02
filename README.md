# clj-hgvs

[![Clojars Project](https://img.shields.io/clojars/v/clj-hgvs.svg)](https://clojars.org/clj-hgvs)
[![Build Status](https://travis-ci.org/chrovis/clj-hgvs.svg?branch=master)](https://travis-ci.org/chrovis/clj-hgvs)
[![codecov](https://codecov.io/gh/chrovis/clj-hgvs/branch/master/graph/badge.svg)](https://codecov.io/gh/chrovis/clj-hgvs)

Clojure(Script) library for handling [HGVS](http://varnomen.hgvs.org/).

## Features

clj-hgvs provides:

* Data structure of HGVS
* Parser of HGVS text
* Formatter to HGVS text

## Installation

Clojure CLI/deps.edn:

```clojure
clj-hgvs {:mvn/version "0.3.1"}
```

Leiningen/Boot:

```clojure
[clj-hgvs "0.3.1"]
```

To use clj-hgvs with Clojure 1.8, you must include a dependency on
[clojure-future-spec](https://github.com/tonsky/clojure-future-spec).

## Usage

### Documentation

[API Reference](https://chrovis.github.io/clj-hgvs/)

### Basics

```clojure
(require '[clj-hgvs.core :as hgvs])

;; `parse` parses HGVS text, returning HGVS map.
(def hgvs1 (hgvs/parse "NM_005228.3:c.2573T>G"))

hgvs1
;;=> {:transcript "NM_005228.3"
;;    :kind :cdna
;;    :mutation #clj_hgvs.mutation.DNASubstitution
;;              {:coord #clj_hgvs.coordinate.CDNACoordinate {:position 2573
;;                                                           :offset 0
;;                                                           :region nil}
;;               :ref "T"
;;               :type ">"
;;               :alt "G"}}

;; `format` returns HGVS text.
(hgvs/format hgvs1)
;;=> "NM_005228.3:c.2573T>G"
```

## License

Copyright 2017-2019 [Xcoo, Inc.](https://xcoo.jp/)

Licensed under the [Apache License, Version 2.0](LICENSE).
