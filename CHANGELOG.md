# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.1.2]
### Changed
- Updated description in manifest

## [v0.1.1]
### Fixed
* Default basecaller_basemod_threads value
* Undefined `colors` variable

## [v0.1.0]
### Added
* Workflow will now output pass and fail CRAM
    * Reads are separated into pass and fail based on their mean qscore as calculated by dorado
    * The threshold can be changed with `--qscore_filter`
### Changed
* Improved `--help` documentation
### Fixed
* Workflow will exit with "No files match pattern" if no suitable files are found to basecall
    * Ensure to set `--dorado_ext` to `fast5` or `pod5` as appropriate

## [v0.0.1]
- Initial release of wf-basecalling supporting the Dorado basecaller

## [v0.0.0]
- Initialised wf-basecalling from wf-template #30ff92d
