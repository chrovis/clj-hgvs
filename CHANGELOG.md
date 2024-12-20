# Changelog

## [0.5.1] - 2024-12-05

### Fixed

- Fix protein-insertion regex. [#13](https://github.com/chrovis/clj-hgvs/pull/13)

## [0.5.0] - 2024-10-16

### BREAKING

- Fix uncertain bases and amino acids format of ins and delins because HGVS nomenclature was updated.
  - DNA: `ins(10)` -> `insN[10]`
  - RNA: `ins(10)` -> `insn[10]`
  - Protein: `ins10` -> `insX[10]`

### Fixed

- Fix uncertain insertion. [#10](https://github.com/chrovis/clj-hgvs/pull/10)

## [0.4.7] - 2023-04-25

### Fixed

- Fix format of uncertain protein mutation. [#8](https://github.com/chrovis/clj-hgvs/pull/8)

## [0.4.6] - 2022-07-27

### Changed

- Update dependencies. [#6](https://github.com/chrovis/clj-hgvs/pull/6)

## [0.4.5] - 2021-11-08

### Added

- Add a repairer for protein insertion to extension. [#4](https://github.com/chrovis/clj-hgvs/pull/4)
- Parse HGVS with Ensembl stable ID. [#5](https://github.com/chrovis/clj-hgvs/pull/5)

## [0.4.4] - 2020-09-28

### Added

- Repair repeated sequences which have the same end.
- Repair illegal substitutions to an inversion. [#3](https://github.com/chrovis/clj-hgvs/pull/3)

## [0.4.3] - 2020-03-23

### Added

- Add "fix-protein-repeated-seqs-pos" repairer.
- Make coordinates calculable.
- Add equiv fn.

## [0.4.2] - 2019-12-27

### Added

- Add HGVS repair function.

## [0.4.1] - 2019-12-09

### Changed

- Make order of start and end strict.

### Fixed

- Fix ter substitution.
- Fix ter insertion.

## [0.4.0] - 2019-08-01

### BREAKING

- HGVS data structure changes from map to record (`clj-hgvs.core/HGVS`).
- (n)cdna is renamed to (non-)coding-dna to avoid misunderstanding. See [#2](https://github.com/chrovis/clj-hgvs/issues/2) for more information.
    - kind: `:cdna`, `:ncdna` → `:coding-dna`, `:non-coding-dna`
    - coordinate (`clj-hgvs.coordinate`):
        - `CDNACoordinate`, `NCDNACoordinate` → `CodingDNACoordinate`, `NonCodingDNACoordinate`
        - `cdna-coordinate`, `ncdna-coordinate` → `coding-dna-coordinate`, `non-coding-dna-coordinate`
        - `parse-cdna-coordinate`, `parse-ncdna-coordinate` → `parse-coding-dna-coordinate`, `parse-non-coding-dna-coordinate`

### Added

- Add a tagged literal.
- Add custom print method.
- Support circular DNA. [#1](https://github.com/chrovis/clj-hgvs/issues/1)

### Changed

- Change HGVS map to record.
- Rename (n)cdna to (non-)coding-dna. [#2](https://github.com/chrovis/clj-hgvs/issues/2)

## [0.3.1] - 2018-12-12

### Added

- Add a macro to disable validation in a scope.
- Support parsing protein unknown alt "p.Met1?".

## [0.3.0] - 2018-11-15

### BREAKING

- clj-hgvs 0.3.0 uses clojure.spec for HGVS validation. To use clj-hgvs with
  clojure 1.8, you must include a dependency on
  [clojure-future-spec](https://github.com/tonsky/clojure-future-spec).
- A new ter codon site of protein frame shift is not displayed by default.
  Set `:show-ter-site?` option to `true` to display that.

### Added

- Add show-ter-site? option for frame shift format.

### Changed

- Use spec for validation.
- Refactor AA long-short conversion.

### Removed

- Drop clojure 1.7 support.

### Fixed

- Fix plain/restore fns of protein extension.

## [0.2.4] - 2018-10-02

## [0.2.3] - 2017-12-01

## [0.2.2] - 2017-08-29

## [0.2.1] - 2017-06-05

## [0.2.0] - 2017-05-18

## [0.1.3] - 2017-05-12

## [0.1.2] - 2017-05-08

## [0.1.1] - 2017-05-01

## 0.1.0 - 2017-04-17

[Unreleased]: https://github.com/chrovis/clj-hgvs/compare/0.5.1...HEAD
[0.5.1]: https://github.com/chrovis/clj-hgvs/compare/0.5.0...0.5.1
[0.5.0]: https://github.com/chrovis/clj-hgvs/compare/0.4.7...0.5.0
[0.4.7]: https://github.com/chrovis/clj-hgvs/compare/0.4.6...0.4.7
[0.4.6]: https://github.com/chrovis/clj-hgvs/compare/0.4.5...0.4.6
[0.4.5]: https://github.com/chrovis/clj-hgvs/compare/0.4.4...0.4.5
[0.4.4]: https://github.com/chrovis/clj-hgvs/compare/0.4.3...0.4.4
[0.4.3]: https://github.com/chrovis/clj-hgvs/compare/0.4.2...0.4.3
[0.4.2]: https://github.com/chrovis/clj-hgvs/compare/0.4.1...0.4.2
[0.4.1]: https://github.com/chrovis/clj-hgvs/compare/0.4.0...0.4.1
[0.4.0]: https://github.com/chrovis/clj-hgvs/compare/0.3.1...0.4.0
[0.3.1]: https://github.com/chrovis/clj-hgvs/compare/0.3.0...0.3.1
[0.3.0]: https://github.com/chrovis/clj-hgvs/compare/0.2.4...0.3.0
[0.2.4]: https://github.com/chrovis/clj-hgvs/compare/0.2.3...0.2.4
[0.2.3]: https://github.com/chrovis/clj-hgvs/compare/0.2.2...0.2.3
[0.2.2]: https://github.com/chrovis/clj-hgvs/compare/0.2.1...0.2.2
[0.2.1]: https://github.com/chrovis/clj-hgvs/compare/0.2.0...0.2.1
[0.2.0]: https://github.com/chrovis/clj-hgvs/compare/0.1.3...0.2.0
[0.1.3]: https://github.com/chrovis/clj-hgvs/compare/0.1.2...0.1.3
[0.1.2]: https://github.com/chrovis/clj-hgvs/compare/0.1.1...0.1.2
[0.1.1]: https://github.com/chrovis/clj-hgvs/compare/0.1.0...0.1.1
