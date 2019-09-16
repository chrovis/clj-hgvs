# Changelog

## [Unreleased]

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

[Unreleased]: https://github.com/chrovis/clj-hgvs/compare/0.4.0...HEAD
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
