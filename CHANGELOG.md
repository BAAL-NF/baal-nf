# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

