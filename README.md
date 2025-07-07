# Basecalling workflow

Helper workflow for basecalling nanopore sequencing data.



## Introduction

In brief this workflow can be used to perform: 

+ Basecalling of a directory of pod5 or fast5 signal data
+ Basecalling in Duplex mode
+ Modified basecalling
+ Basecalling in real time
+ Output basecalled sequences in various formats: FASTQ, CRAM or Unaligned BAM
+ If a reference is provided a sorted and indexed BAM or CRAM will be output
for basecalling a directory of `pod5` or `fast5` signal data with `dorado`
and aligning it with `minimap2` to produce a sorted, indexed CRAM.




## Compute requirements

Recommended requirements:

+ CPUs = 64
+ Memory = 256GB

Minimum requirements:

+ CPUs = 8
+ Memory = 64GB

Approximate run time: Variable depending on coverage, genome size, model of choice and GPU model.

ARM processor support: False




## Install and run


These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://docs.docker.com/get-started/)
or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
to provide isolation of the required software.
Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below.

It is not required to clone or download the git repository
in order to run the workflow.
More information on running EPI2ME workflows can
be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This will pull the repository in to the assets folder of
Nextflow and provide a list of all parameters
available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-basecalling --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull epi2me-labs/wf-basecalling
```

A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-basecalling/wf-basecalling-demo.tar.gz
tar -xzvf wf-basecalling-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using:
```
nextflow run epi2me-labs/wf-basecalling \
	--basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' \
	--dorado_ext 'pod5' \
	--input 'wf-basecalling-demo/input' \
	--ref 'wf-basecalling-demo/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
	--remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v2' \
	-profile standard
```

For further information about running a workflow on
the command line see https://labs.epi2me.io/wfquickstart/




## Related protocols

<!---Hyperlinks to any related protocols that are directly related to this workflow, check the community for any such protocols.--->

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Input example

This workflow accepts a folder containing FAST5 or POD5 files as input.
The folder may contain other folders of FAST5 or POD5 files and all files will be processed by the workflow.





## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| input | string | Directory containing FAST5 (or POD5) signal for basecalling. | This directory will be searched recursively. All FAST5 or POD5 files (depending on which extension you select in the Basecalling Options) in this directory or any subdirectory (no matter how deep) will be basecalled. |  |
| ref | string | Optional reference FASTA file to align basecalled reads to. | Without a reference, basecalls are output to unaligned CRAM. When using a reference, take care to retain this FASTA file as the output CRAM file cannot be read without the reference it was aligned to. |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all files. |  | output |
| sample_name | string | Sample name to prefix file names of workflow outputs. |  | SAMPLE |
| output_fmt | string | Desired file format of files created by basecalling and alignment. | FASTQ can only be output when a reference has not been provided. Aligned output will always be written to CRAM unless BAM is selected. | cram |
| igv | boolean | Visualize outputs in the EPI2ME IGV visualizer. | Enabling this option will visualize the output alignment files in the EPI2ME desktop app IGV visualizer. | False |


### Basecalling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| basecaller_cfg | string | Name of the model to use for converting signal. | Required for basecalling. The model list only shows models that are compatible with this workflow. |  |
| duplex | boolean | Run the basecaller in duplex mode. | By default the workflow conducts simplex basecalling. If you used a chemistry and flowcell combination that supports duplex reads, you may switch this option on. This option is incompatible with the watch_path option due to the way the input files must be traversed in order to find duplex pairs. | False |
| remora_cfg | string | Name of the model to use for calling modified bases. | Required for calling modified bases while basecalling. The model list only shows models that are compatible with this workflow. |  |
| dorado_ext | string | File extension for Dorado inputs. | Set this to fast5 if you have not converted your fast5 to pod5. It is recommended to [convert existing fast5 files to pod5 for use with Dorado](https://github.com/nanoporetech/pod5-file-format/blob/master/python/README.md#pod5-convert-from-fast5). | pod5 |
| poly_a_config | string | Provide this TOML file to turn on and configure dorado poly(A) calling. | This TOML file allows you to turn on and configure poly(A) tail calling options in dorado. This feature is described [here](https://github.com/nanoporetech/dorado?tab=readme-ov-file#polya-tail-estimation). |  |
| barcode_kit | string | Name of the kit to use for barcoding. Demultiplex the data. | Providing a kit here will instruct the workflow to demultiplex your 'pass' data to BAM files, which can be found in your output directory under the folder 'demuxed' in a structure reminiscent of MinKNOW. |  |


### Advanced basecalling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| output_pod5 | boolean | Save the converted POD5 when running in duplex with FAST5 inputs. | Dorado duplex only supports POD5 input. The workflow will automatically convert FAST5 input to POD5 when duplex calling. By default, converted POD5 are deleted to save disk space. Enabling this option will make the workflow output converted POD5 files to a subfolder within the output directory. | False |
| qscore_filter | number | Mean qscore by which to filter reads. Inclusive such that reads with score >= qscore_filter are kept. | The mean qscore of reads is calculated by dorado and rounded to an integer by dorado and stored as a tag in dorado's SAM output. The pipeline separates reads into pass and fail categories based on this SAM tag. | 10 |
| cuda_device | string | GPU device to use for basecalling [cuda:all]. | For local execution this can be used to pin GPU tasks to one (or more) specific GPU devices. Use cuda:all to use all available GPU devices, or cuda:idx[,idx,...] where idx is an index number(s) of GPU device(s) to use. | cuda:all |
| basecaller_model_path | string | Override the named basecalling model with a custom basecalling model. | For typical use, users should set --basecaller_cfg which will use a named model from inside the container. Experimental or custom basecallers will not be available in the container and can be loaded from the host with --basecaller_model_path. |  |
| remora_model_path | string | Override the named remora model with a custom remora model. | For typical use, users should set --remora_cfg which will use a named model from inside the container. Experimental or custom models will not be available in the container and can be loaded from the host with --remora_model_path. |  |
| basecaller_basemod_threads | number | Number of threads to use for base modification calling. | You must set this to > 0 when using a modbase aware model. Modbase calling does not require much additional CPU and should be set carefully when using GPU servers with a small number of CPUs per GPU. | 2 |
| basecaller_args | string | Additional command line arguments to pass to the basecaller process. |  |  |
| demux_args | string | Additional command line arguments to pass to the basecaller barcoding process. |  |  |


### Multiprocessing Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| ubam_map_threads | integer | Set max number of threads to use for aligning reads from uBAM (limited by config executor cpus) |  | 8 |
| ubam_sort_threads | integer | Set max number of threads to use for sorting and indexing aligned reads from uBAM (limited by config executor cpus) |  | 3 |
| ubam_bam2fq_threads | integer | Set max number of threads to use for uncompressing uBAM and generating FASTQ for alignment (limited by config executor cpus) |  | 1 |
| merge_threads | integer | Set max number of threads to use for merging BAM files (limited by config executor cpus) |  | 4 |
| stats_threads | integer | Set max number of threads to use for getting stats from output files. (limited by config executor cpus) |  | 4 |


### Real Time Analysis Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| watch_path | boolean | Enable to continuously watch the input directory for new input files. Reads will be analysed as they appear. | This option enables the use of Nextflow's directory watching feature to constantly monitor input directories for new files. As soon as files are written by an external process Nextflow will begin analysing these files. The workflow will accumulate data over time to produce an updating report. Real time analysis of duplex data may lead to lower duplex rates than what would have been obtained by running basecalling after sequencing. | False |
| read_limit | integer | Stop processing data when a particular number of reads have been analysed. | By default the workflow will run indefinitely when using the real time watch path option. This will set the upper bound on the number of reads that will be analysed before the workflow is automatically stopped and no more data is analysed. |  |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | wf-basecalling-report.html | Report summarising the work done by the basecalling workflow | per-sample |
| Simplex alignment file of passed reads | {{ alias }}.pass.simplex.{{ format }} | BAM or CRAM file of simplex reads for the sample that pass QC filtering. | per-sample |
| Duplex alignment file of passed reads | {{ alias }}.pass.duplex.{{ format }} | BAM or CRAM file of duplex reads for the sample that pass QC filtering. Created if duplex basecalling is requested. | per-sample |
| Simplex alignment file index of passed reads | {{ alias }}.pass.simplex.{{ format }}.{{ index_format }} | The index of the resulting BAM or CRAM file with the simplex reads that pass QC filtering. | per-sample |
| Duplex alignment file index of passed reads | {{ alias }}.pass.duplex.{{ format }}.{{ index_format }} | The index of the resulting BAM or CRAM file with the duplex reads that pass QC filtering. Created if duplex basecalling is requested. | per-sample |
| Simplex alignment file of failed reads | {{ alias }}.fail.simplex.{{ format }} | BAM or CRAM file of simplex reads for the sample that fail QC filtering. | per-sample |
| Duplex alignment file of failed reads | {{ alias }}.fail.duplex.{{ format }} | BAM or CRAM file of duplex reads for the sample that fail QC filtering. Created if duplex basecalling is requested. | per-sample |
| Simplex alignment file index of failed reads | {{ alias }}.fail.simplex.{{ format }}.{{ index_format }} | The index of the resulting BAM or CRAM file with the simplex reads that fail QC filtering. | per-sample |
| Duplex alignment file index of failed reads | {{ alias }}.fail.duplex.{{ format }}.{{ index_format }} | The index of the resulting BAM or CRAM file with the duplex reads that fail QC filtering. Created if duplex basecalling is requested. | per-sample |
| Index of the reference FASTA file | {{ ref }}.fai | Index of the reference FASTA file. | aggregated |
| JSON configuration file for IGV browser | igv.json | JSON configuration file to be loaded in IGV for visualising alignments against the reference genome. | aggregated |




## Pipeline overview

### 1. Prerequisites

The workflow uses [Dorado](https://github.com/nanoporetech/dorado) for basecalling which includes the use of [Remora](https://github.com/nanoporetech/remora) for modified basecalling.
Basecalling with `Dorado` requires an NVIDIA GPU with [Pascal architecture or newer](https://www.nvidia.com/en-gb/technologies/) and at least 8 GB of vRAM.

#### Windows

Windows should not be considered as a supported operating systems for wf-basecalling as we do not directly support configuration of accelerated computing through WSL2 and Docker.
Although we do not offer support, it is possible to set up Docker to use GPUs for most versions of Windows 11 and some versions of Windows 10 and we direct users to the [CUDA on WSL User Guide](https://docs.nvidia.com/cuda/wsl-user-guide/index.html).
Users should take note of the support constraints section to ensure their environment is suitable before following the guidance. **Do not install an NVIDIA driver into your WSL2 environment**.
Users are encouraged to download Dorado for Windows from the [Dorado GitHub repository](https://github.com/nanoporetech/dorado#installation).

#### MacOS

MacOS should not be considered as a supported operating systems for wf-basecalling as we do not support accelerated computing through Docker on MacOS.
On MacOS, GPU support through Docker remains in technical infancy. In addition, the containers we provide will not be able to leverage the M1 and M2 architecture and will not run as performantly as if Dorado had been run natively.
Users are encouraged to download Dorado for MacOS directly from the [Dorado GitHub repository](https://github.com/nanoporetech/dorado#installation).

#### Linux

When using Docker for accelerated computing on Linux, you will need the `nvidia-container-toolkit` installed.
If you observe the error "could not select device driver with capabilities gpu", you should follow the instructions to install `nvidia-container-toolkit` [here](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html#setting-up-nvidia-container-toolkit). You will need to follow the steps to:

- Setup the package repository and the GPG key (ignore the box about experimental releases)
- Update package listings
- Install nvidia-container-toolkit
- Configure the Docker daemon to recognize the NVIDIA Container Runtime
- Restart the Docker daemon to complete the installation after setting the default runtime

By default, workflows are configured to run GPU tasks in serial. That is, only one basecalling task will be run at a time. This is to prevent the GPU from running out of memory on local execution.
When running workflows on a cluster, or in a cloud where GPU resources are isolated from one another, users should specify `-profile discrete_gpus` as part of the command invocation. This will allow for parallel execution of GPU tasks.
You should ask your system administrator if you need to configure any additional options to leverage GPUs on your cluster. For example, you may need to provide a special string to the workflow's `--cuda_device` option to ensure tasks use the GPU assigned to them by the job scheduler.

### 2. Choosing a model

To select the relevant model see the `dorado` repository for a [table of available models](https://github.com/nanoporetech/dorado#available-basecalling-models) to select for `--basecaller_cfg` and `--remora_cfg`.

### 3. Aligning to a reference

The workflow can optionally perform the alignment of the basecalled data using [minimap2](https://github.com/lh3/minimap2) to a reference of choice, provided with the `--ref` option.
Additionally, the workflow will generate an IGV configuration file. This file allows the user to view the filtered aligned BAM in the EPI2ME Desktop Application in the Viewer tab.

### 4. Duplex calling

wf-basecalling supports [duplex calling](https://github.com/nanoporetech/dorado#duplex), which is enabled with the `--duplex` option. If you used a chemistry and flowcell combination that supported duplex reads, you should switch this option on. The resulting BAM/CRAM will quality filtered and then automatically split in separate BAM/CRAM files for the simplex and duplex reads.
Since `dorado duplex` requires the inputs to be in `pod5` format, the workflow will perform the conversion automatically using [pod5 convert fast5](https://github.com/nanoporetech/pod5-file-format/blob/master/python/pod5/README.md#pod5-convert-fast5). These files are normally deleted upon completion of the analysis, but can optionally be saved by the user by providing the `--output_pod5` option.

### 5. Real-time analysis

wf-basecalling can perform the basecalling as the pod5 files are generated. To enable this, provide the `--watch_path` option. The workflow will process the newly generated files as soon as they become available.

### 6. Barcode classification and demultiplexing

wf-basecalling can perform data demultiplexing by providing the appropriate barcoding kit with the `--barcode_kit` option.
This will generate a new `{{ out_dir }}/demuxed` directory, with one subfolder for each barcode and one additional `unclassified` folder for reads that cannot be demultiplexed. This option is not available for `dorado duplex`.
Please note that the demultiplexed reads will always be in BAM format, even when the user sets `--output_bam false`.



## Troubleshooting

* Duplex mode with wf-basecalling is reliant on internal optimisations to organise input files for better duplex rates, which is not possible when using streaming basecalling; therefore duplex combined with the `--watch_path` option could lead to lower duplex rates than what would be achieved running the algorithm after sequencing is completed.
* Renaming, moving or deleting the reference genome or the output directory from the location provided at runtime will cause IGV to not load anymore.



## FAQ's

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-basecalling/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).




## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



