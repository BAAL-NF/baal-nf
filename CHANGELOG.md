# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
### Changed

## [[0.6.1]](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf/compare/0.6.0...0.6.1) - 15.06.2021
### Added
- New config option to export RDS files for debugging purposes.
### Changed

## [[0.6.0]](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf/compare/0.5.3...0.6.0) - 5.02.2021
### Added
- Support for google life sciences API through gls profile
### Changed
- Fix bug preventing users from overriding default output directories
- Move to baal-nf-env 1.2.0, adding support for GLS
- Refactor container configuration for easier container updates in future

## [[0.5.3]](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf/compare/0.5.2...0.5.3) - 29.01.2021
### Changed
- Add PID isolation to eddie configuration profile to fix issues with child processes getting killed.

## [[0.5.2]](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf/compare/0.5.1...0.5.2) - 26.11.2020
### Changed
- Memory optimisations for BaalChIP
- Move to new version of baal-nf-env to fix issue with compressed bedfiles

## [[0.5.1]](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf/compare/0.5.0...0.5.1) - 26.11.2020
### Changed
- Support a mixture of compressed and uncompressed bedfiles

## [[0.5.0]](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf/compare/0.4.0...0.5.0) - 23.11.2020

### Added
- Support for downloading fastq-screen genomes automatically

### Changed
- Updated `Readme.md` with details of configuration profiles
- Cleaned up unused process labels

## [[0.4.0]](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf/compare/0.3.0...0.4.0) - 20.11.2020

### Changed
- Renamed main script to `main.nf` to conform with nf-core guidelines
- Linted using groovyLint
- Add enrichment analysis using GAT and then ensembl annotations
- Bump container versions to 1.1.0 for both containers

## [[0.3.0]](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf/compare/0.2.1...0.3.0) - 28.8.2020

### Added
- Support for automatic staging of the iGenomes hg19 genome
- Pipeline reporting is now part of the reports by default
- Support for staging directly from URLs and s3 buckets

### Changed
- New input format with explicit paired-end input
- Required configuration files are now (mostly) staged through nextflow channels, reducing the need for custom mount points
- Added test profile which stages a test set from github

## [[0.2.0]](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf/compare/0.1.1...0.2.1) - 28.8.2020

### Changed
- Fix bug causing incorrect rendering of Baal ChIP reports
- Move to absolute paths for all input files

## [[0.1.1]](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf/compare/0.1.0...0.1.1) - 26.8.2020

### Changed

- Add the ability for overlap_beds.py to deal with empty bed files gracefully

## [0.1.0] - 18.8.2020

### Added

- First full pipeline release.

