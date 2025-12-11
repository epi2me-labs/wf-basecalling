# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.5.6]
### Changed
- Updated to wf-template v5.6.3 to maintain compliance with our latest wf-template standard, changing:
  - Pipeline overview now appears before pipeline parameters in README.
  - ezCharts plotting library has been updated to 0.15.1, there are no user facing changes to plots.
  - Fastcat FASTQ pre-processing program has been updated to 0.24.2, it is more robust to malformed FASTQ input.
  - CHANGELOG to be compliant with our formatting rules.
- Updated to Dorado [v1.3.0](https://github.com/nanoporetech/dorado/releases/tag/v1.3.0).
  - Files output by `dorado demux` are now organized into subdirectories following the structure `demuxed/sample_id/run_id/bam_pass/barcode_id/file.bam`.

## [v1.5.5]
This patch release of `wf-basecalling` updates internal workflow naming, and does not affect any workflow outputs.
### Changed
- Removed workflow suffix from workflow title.

## [v1.5.4]
This patch update to wf-basecalling fixes the "Not a valid path value type" error observed by some users when demultiplexing, and updates Dorado to v1.1.1 to benefit from improved HAC basecalling speed. It is recommended that users update to this version to keep up to date with Dorado's model and basecalling performance improvements.

### Changed
- Memory requirements reduced to account for cloud instances expressing capacity in GiB.
- Updated Dorado to [v1.1.1](https://github.com/nanoporetech/dorado/releases/tag/v1.1.1):
    - Update to improved 5mC_5hmC and 5mCG_5hmCG modified base models for DNA v5.2.0 HAC and SUP.
    - Improve speed of HAC basecalling models on a range of Nvidia GPUs.

### Fixed
- `store_dir` parameter format incorrectly declared in the schema. This does not affect this workflow as it does not use the storeDir directive and has been changed to maintain compliance with our latest testing standard.
-  `Not a valid path value type: org.codehaus.groovy.runtime.NullObject (null)` when demuxing and outputting FASTQ without performing alignment. FASTQ is now correctly demuxed into the expected folder output. 

## [v1.5.3]
This release of wf-basecalling is made to support the development of other EPI2ME workflows that handle their own Poly(A) configuration. Users do not need to adopt this release.
### Changed
- Poly(A) configuration for dorado signal ingress now expects a channel
- Updated to wf-template v5.6.2 to maintain compliance with our latest wf-template standard: this does not impact the workflow.

## [v1.5.2]
This version of wf-basecalling updates Dorado to v1.0.2.
Since Dorado v1.0.0, the 4 kHz models have been deprecated. Users with 4 kHz data must use wf-basecalling v1.5.1.
Additionally, Dorado no longer supports FAST5 input. wf-basecalling will automatically convert FAST5 to POD5 but this is time consuming, users are strongly encouraged to provide POD5 files where possible.
### Changed
- Updated to Dorado [v1.0.2](https://github.com/nanoporetech/dorado/releases/tag/v1.0.2).
- All FAST5 inputs are converted to POD5 for basecalling as FAST5 files are no longer supported by Dorado, previously this was only done for duplex basecalling.
### Removed
- Support for 4 kHz models as they are not supported by Dorado >v1.0.0.

## [v1.5.1]
The updates in this release do not affect wf-basecalling but are required for EPI2ME workflows that make use of wf-basecalling, to maintain compliance with our latest wf-template standard.
Users do not need to update to this release.
### Changed
- Updated to wf-template v5.6.1, changing:
    - pre-commit configuration to resolve an internal dependency problem with flake8. This has no effect on the workflow.
    - Log banner art to say "EPI2ME" instead of "EPI2ME Labs" to match current branding. This has no effect on the workflow outputs.

## [v1.5.0]
This version of wf-basecalling updates Dorado to v0.9.5 which should improve the speed of basecalling on some GPU architectures.
Dorado v0.9.5 increases the minimum NVIDIA Driver requirement to 525.105.
### Changed
- Updated Dorado to [v0.9.5](https://github.com/nanoporetech/dorado/releases/tag/v0.9.5)
- Alignment uses the high quality long read preset (-x lr:hq) to reduce mapping time.
- Basecaller model options are now reverse version sorted in the workflow schema to ensure newer models appear at the top of drop-downs and listings of available models.

### Fixed
- "Input data problem" error in downstream workflows necessitating use of the override_basecaller_cfg option. Relevant metadata from the input XAM header is now retained after alignment to ensure that the basecaller configuration is automatically detected.

## [v1.4.7]
This maintenance release updates the models used for outputting results.
This release is to support our other workflows.
Users do not need to update to this release.
### Changed
- Removed Pydantic and autogenerated models from results_schema.yml and created a new model namespace using dataclasses.

## [v1.4.6]
This maintenance release updates the version of our code that plots our post-workflow reports.
This release is to support our other workflows.
Users do not need to update to this release.
### Changed
- Updated common Docker image to update ezcharts to v0.12.0. This improves formatting of plots in the report.
### Fixed
- Typo in schema.

## [v1.4.5]
### Changed
- Reconciled workflow with wf-template v5.3.4.

## [v1.4.4]
### Changed
- Updated Dorado to [v0.9.0](https://github.com/nanoporetech/dorado/releases/tag/v0.9.0)

## [v1.4.3]
### Changed
- Reconciled workflow with wf-template v5.3.3.

## [v1.4.2]
### Added
- q-score filter added to signal ingress.
### Changed
- Updated Dorado to [v0.8.3](https://github.com/nanoporetech/dorado/releases/tag/v0.8.3)
- Reconciled workflow with wf-template v5.3.1.

## [v1.4.1]
### Changed
- Reconciled workflow with wf-template v5.3.0.
- Updated Dorado to [v0.8.1](https://github.com/nanoporetech/dorado/releases/tag/v0.8.1)

## [v1.4.0]
### Added
- IGV configuration file with `--ref --igv` options and either `--output_fmt bam` or `--output_fmt cram`.
- Support for gzipped reference genomes.
- `output_fmt` selects the output format for basecalled and aligned files.
### Changed
- Updated Dorado to [v0.8.0](https://github.com/nanoporetech/dorado/releases/tag/v0.8.0)
- Reconciled workflow with wf-template v5.2.6.
- Do not emit the reference FASTA file.
- Collapse redundant RG and PG header lines when emitting BAM or CRAM.
### Fixed
- Workflow starting with `--duplex --barcode_kit`, despite duplex not supporting barcoding.
- Workflow crashing with `--ref {{ reference }} --barcode_kit`.
- Aligned reads will no longer be trimmed when demuxing to preserve mapping information.
- Workflow emits confusing warning about Bonito filtering when using Dorado.
### Removed
- `fastq_only` and `output_bam` options replaced by `output_fmt`.
    - `--output_fmt fastq` can be used to output unaligned FASTQ instead of unaligned CRAM.
    - `--output_fmt bam` can be used to output unaligned or aligned BAM instead of CRAM.

## [v1.3.0]
### Added
- Modified base calling with `--duplex`.
- APK 5.0.0 model.
### Changed
- Updated Dorado to v0.7.2 (see https://github.com/nanoporetech/dorado/releases/tag/v0.7.2)

## [v1.2.2]
### Changed
- Bug fix for downstream workflows and `--poly_a_config` which does not affect normal workflow use.

## [v1.2.1]
### Added
- Output channel for demuxed BAM files for downstream use.

## [v1.2.0]
### Added
- Support for `dorado demux` to demultiplex barcoded runs. Specify your `--barcode_kit` to activate demultiplexing.
- Support for poly(a) tail length estimation with `--poly_a_config`. You can configure by providng a TOML file to `--poly_a_config` which is described in detail [here](https://github.com/nanoporetech/dorado?tab=readme-ov-file#polya-tail-estimation)
### Changed
- Updated Dorado to v0.7.1 (see https://github.com/nanoporetech/dorado/releases/tag/v0.7.1)

## [v1.1.9]
### Fixed
- Report crashing when no data are present in the input pod5.
- Reconciled workflow with wf-template v5.1.3.
- Updated Dorado to v0.7.0 (see https://github.com/nanoporetech/dorado/releases/tag/v0.7.0)
- Added new DNA and RNA 5.0.0 models.

## [v1.1.8]
### Changed
- Updated Dorado to v0.6.0 (see https://github.com/nanoporetech/dorado/releases/tag/v0.6.0)

## [v1.1.7]
### Fixed
- Workflow accepting incompatible `--fastq_only` and `--duplex` options
- Dynamically updated report in `--watch_path` mode.

## [v1.1.6]
### Fixed
- qscore_filter inadvertently disabled in v1.1.5

## [v1.1.5]
### Changed
- Minor update to default resource requests on dorado task.
### Added
- Experimental feature switch.

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
### Added
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
Initial release of wf-basecalling supporting the Dorado basecaller

## [v0.0.0]
Initialised wf-basecalling from wf-template #30ff92d

