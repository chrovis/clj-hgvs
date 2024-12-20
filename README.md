# clj-hgvs

[![Clojars Project](https://img.shields.io/clojars/v/clj-hgvs.svg)](https://clojars.org/clj-hgvs)
[![build](https://github.com/chrovis/clj-hgvs/actions/workflows/build.yml/badge.svg)](https://github.com/chrovis/clj-hgvs/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/chrovis/clj-hgvs/branch/master/graph/badge.svg)](https://codecov.io/gh/chrovis/clj-hgvs)

Clojure(Script) library for handling [HGVS](https://hgvs-nomenclature.org/).

## Features

clj-hgvs provides:

* Data structure for HGVS
* HGVS text parser
* HGVS text formatter

## Installation

Clojure CLI/deps.edn:

```clojure
clj-hgvs {:mvn/version "0.5.1"}
```

Leiningen/Boot:

```clojure
[clj-hgvs "0.5.1"]
```

To use clj-hgvs with Clojure 1.8, you must include a dependency on
[clojure-future-spec](https://github.com/tonsky/clojure-future-spec).

## Breaking Changes in 0.5.0

- Fix uncertain bases and amino acids format of ins and delins because HGVS nomenclature was updated.
  - DNA: `ins(10)` -> `insN[10]`
  - RNA: `ins(10)` -> `insn[10]`
  - Protein: `ins10` -> `insX[10]`

## Breaking Changes in 0.4.0

- HGVS data structure changes from map to record (`clj-hgvs.core/HGVS`).
- (n)cdna is renamed to (non-)coding-dna to avoid misunderstanding.

See [CHANGELOG](CHANGELOG.md) for more information.

## Usage

### Documentation

[API Reference](https://chrovis.github.io/clj-hgvs/)

### Basics

```clojure
(require '[clj-hgvs.core :as hgvs])

;; `parse` parses a HGVS text, returning a HGVS record.
(def hgvs1 (hgvs/parse "NM_005228.3:c.2573T>G"))

hgvs1
;;=> #clj_hgvs.core.HGVS
;;   {:transcript "NM_005228.3"
;;    :kind :coding-dna
;;    :mutation #clj_hgvs.mutation.DNASubstitution
;;              {:coord #clj_hgvs.coordinate.CodingDNACoordinate
;;                      {:position 2573
;;                       :offset 0
;;                       :region nil}
;;               :ref "T"
;;               :type ">"
;;               :alt "G"}}

;; `format` returns a HGVS text.
(hgvs/format hgvs1)
;;=> "NM_005228.3:c.2573T>G"
```

### Tagged Literal

`#clj-hgvs/hgvs` tagged literal is useful for easy and readable definition of a
HGVS data.

```clojure
#clj-hgvs/hgvs "NM_005228.3:c.2573T>G"
```

### Formatter Styles

`clj-hgvs.core/format` has various options for specifying HGVS styles.

```clojure
(hgvs/format #clj-hgvs/hgvs "NM_005228.3:c.2307_2308insGCCAGCGTG"
             {:ins-format :count})
;;=> "NM_005228.3:c.2307_2308ins(9)"

(hgvs/format #clj-hgvs/hgvs "p.Leu858Arg"
             {:amino-acid-format :short})
;;=> "p.L858R"
```

See [API reference](https://chrovis.github.io/clj-hgvs/clj-hgvs.core.html#var-format)
for all formatter options.

### Equivalence

`clj-hgvs.core/==` tests the fundamental equivalence of the given HGVS.

```clojure
(hgvs/== #clj-hgvs/hgvs "NM_005228:c.2361G>A"
         #clj-hgvs/hgvs "NM_005228.4:c.2361G>A")
;;=> true

(hgvs/== #clj-hgvs/hgvs "p.K53Afs*9"
         #clj-hgvs/hgvs "p.Lys53Alafs")
;;=> true

(hgvs/== #clj-hgvs/hgvs "p.L858R"
         #clj-hgvs/hgvs "p.L858M")
;;=> false
```

### Plain Map Representation

`plain` converts HGVS data to a plain map, and `restore` converts the map back
to the HGVS data. These functions are useful for sending HGVS data through
another codec.

```clojure
(hgvs/plain #clj-hgvs/hgvs "NM_005228.3:c.2573T>G")
;;=> {:transcript "NM_005228.3"
;;    :kind "coding-dna"
;;    :mutation {:mutation "dna-substitution"
;;               :coord {:coordinate "coding-dna"
;;                       :position 2573
;;                       :offset 0
;;                       :region nil}
;;               :ref "T"
;;               :type ">"
;;               :alt "G"}}

(hgvs/restore *1)
;;=> #clj-hgvs/hgvs "NM_005228.3:c.2573T>G"
```

### HGVS Repair

`repair-hgvs-str` attempts to repair an invalid HGVS text.

```clojure
(hgvs/repair-hgvs-str "c.123_124GC>AA")
;;=> "c.123_124delGCinsAA"
```

The repair rules are based on frequent mistakes in popular public-domain
databases such as dbSNP and ClinVar.

You may supply custom repair rules to the second argument:

```clojure
(require '[clojure.string :as string]
         '[clj-hgvs.repairer :as repairer])

(defn lower-case-ext
  [s kind]
  (if (= kind :protein)
    (string/replace s #"EXT" "ext")
    s))

(def my-repairers (conj repairer/built-in-repairers
                        lower-case-ext))

(hgvs/repair-hgvs-str "p.*833EXT*?" my-repairers)
;;=> "p.*833ext*?"
```

## License

Copyright 2017-2024 [Xcoo, Inc.](https://xcoo.jp/)

Licensed under the [Apache License, Version 2.0](LICENSE).
