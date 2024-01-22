# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.1.4]
### Changed
- Updated Dorado to v0.5.2 (see https://github.com/nanoporetech/dorado/releases/tag/v0.5.2)

## [v1.1.3]
### Changed
- Bumped memory directives for intense tasks to reduce likelihood of job failure
- Default to parallel GPU usage when using awsbatch profile
### Removed
- Runtime driver check in Dorado process, as this is no longer available in the Dorado image

## [v1.1.2]
### Changed
- Updated dorado version to v0.5.1 (see https://github.com/nanoporetech/dorado/releases/tag/v0.5.1)

## [v1.1.1]
### Added
- Reintroduced RNA002 models

## [v1.1.0]
### [Added]
- `--duplex` basecalling converts FAST5 to POD5 automatically
    - Converted POD5 files are deleted by default, use `--output_pod5` to output converted POD5 files to the workflow output directory.

### Changed
- Updated Dorado to v0.3.4 (see https://github.com/nanoporetech/dorado/releases/tag/v0.3.4)

## [v1.0.1]
### Fixed
- Workflow crashes with fast5 input
- Workflow fails early when trying to use FAST5 input with Dorado duplex

## [v1.0.0]
### Added
- RNA004 models
- R941 v3.3 5mCG 5hmCG models
- Duplex calling with option `--duplex`
    - Note that duplex calling is not optimised for streaming basecalling with `--watch_path` and may lead to lower duplex yield.
    - Duplex basecalling is currently not compatible with modified basecalling.

### Changed
- Updated Dorado to v0.3.2 (see https://github.com/nanoporetech/dorado/releases/tag/v0.3.2)
- Pascal architecture GPUs are now supported
- Bumped minimum required Nextflow version to 23.04.2
- Users no longer need to provide `--basecaller_cfg custom` and/or `--remora_cfg custom` to override models with `--basecaller_model_path` and/or `--remora_model_path` respectively.

### Fixed
- `bamstats` process very slow when `output_bam` has been selected

## [v0.7.2]
### Added
- v4.2 5mC and 6mA modification models

### Changed
- Updated Dorado to v0.3.1
- GPU tasks are limited to run in serial by default to avoid memory errors
    - Users in cluster and cloud environments where GPU devices are scheduled must use `-profile discrete_gpus` to parallelise GPU work
    - A warning will be printed if the workflow detects it is running non-local execution but the discrete_gpus profile is not enabled
    - Additional guidance on GPU support is provided in our Quickstart
- Bumped minimum required Nextflow version to 22.10.8

## [v0.7.1]
### Fixed
- Command not found on `cram_cache` step
- Typo in report that refers to the workflow as "wf-basecalling-report"

## [v0.7.0]
### Changed
- Updated Dorado to v0.3.0
- BAM may be output **instead** of CRAM by providing `--output_bam`
- `--help` message will list basecalling and modbasecalling models available for use with the workflow

### Added
- v4.2.0 models, which must be used for sequencing runs performed at new 5 kHz sampling rate
- v4.1.0 models replace v4.0.0 models and must be used for sequencing runs performed at 4 kHz sampling rate

### Removed
- v4.0.0 models

### Fixed
- Custom models were previously rejected by the workflow as `basecaller_cfg` and `remora_cfg` are validated against a list of basecalling models installed in the Dorado container.
    - Users should now provide `--basecaller_cfg custom` and/or `--remora_cfg custom` to override models with `--basecaller_model_path` and/or `--remora_model_path` respectively.
    - Providing `--basecaller_cfg custom` or `--remora_cfg custom` without the corresponding `--basecaller_model_path` or `--remora_model_path` will result in an error.

## [v0.6.0]
### Added
- Ability to watch the input path and process files as they become available in real time.

## [v0.5.2]
### Added
- Configuration for running demo data in AWS

## [v0.5.1]
### Fixed
- Missing models from list of valid models
- "dna_r9.4.1_e8_hac@v3.4_5mCG@v0" is now correctly referred to as "dna_r9.4.1_e8_hac@v3.3_5mCG@v0", to match the simplex model version
- "dna_r9.4.1_e8_sup@v3.4_5mCG@v0" is now correctly referred to as "dna_r9.4.1_e8_sup@v3.3_5mCG@v0", to match the simplex model version

## [v0.5.0]
### Changed
- Updated Dorado to v0.2.4
- Updated to Oxford Nanopore Technologies PLC. Public License

### Fixed
- Dorado image correctly ships with CUDA runtime library

## [v0.4.1]
### Fixed
- Input ref channel depleted after first alignment

## [v0.4.0]
### Changed
- Reference is no longer required for basecalling
    - CRAM files with no alignments will be generated if `--ref` is not provided
    - FASTQ may be output **instead** of CRAM by providing `--fastq_only`
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
- `--basecaller_args` may be used to provide custom arguments to the basecalling process

### Changed
- Updated Dorado to v0.1.1
    - Latest models are now v4.0.0
- Workflow prints a more helpful error when Dorado fails due to unknown model name

## [v0.1.2]
### Changed
- Updated description in manifest

## [v0.1.1]
### Fixed
- Default basecaller_basemod_threads value
- Undefined `colors` variable

## [v0.1.0]
### Added
- Workflow will now output pass and fail CRAM
    - Reads are separated into pass and fail based on their mean qscore as calculated by dorado
    - The threshold can be changed with `--qscore_filter`

### Changed
- Improved `--help` documentation

### Fixed
- Workflow will exit with "No files match pattern" if no suitable files are found to basecall
    - Ensure to set `--dorado_ext` to `fast5` or `pod5` as appropriate

## [v0.0.1]
* Initial release of wf-basecalling supporting the Dorado basecaller

## [v0.0.0]
* Initialised wf-basecalling from wf-template #30ff92d

