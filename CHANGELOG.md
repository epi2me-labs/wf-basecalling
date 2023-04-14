# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]
### Changed
- Updated to Oxford Nanopore Technologies PLC. Public License
  
## [v0.4.1]
### Fixed
- Input ref channel depleted after first alignment

## [v0.4.0]
### Changed
- Reference is no longer required for basecalling
    - CRAM files with no alignments will be generated if `--ref` is not provided
    - FASTQ may be output *instead* of CRAM by providing `--fastq_only`
- PG line for converting Dorado SAM output to uBAM is no longer written to output header
- Work directory is automatically cleaned up on successful completion to remove large intermediate files
    - Override this by including `cleanup = false` in a custom Nextflow configuration file
- Number of threads for merging is now configurable for advanced users

## [v0.3.0]
### Changed
- Updated Dorado to v0.2.1
- `--basecaller_cfg` and `--remora_cfg` are now validated against a list of models installed in the Dorado container
### Fixed
- Workflow no longer prints a confusing error when Dorado fails

## [v0.2.0]
### Added
* `--basecaller_args` may be used to provide custom arguments to the basecalling process
### Changed
- Updated Dorado to v0.1.1
    - Latest models are now v4.0.0
- Workflow prints a more helpful error when Dorado fails due to unknown model name

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
