# Basecalling workflow

This repository contains a [nextflow](https://www.nextflow.io/) workflow
for basecalling a directory of `pod5` or `fast5` signal data with `dorado`
and aligning it with `minimap2` to produce a sorted, indexed CRAM.





## Introduction

This workflow introduces users to [`Dorado`](https://github.com/nanoporetech/dorado),
which is now our standard basecaller. `dorado` is still under active development and
will be kept updated as new releases are made. We strongly encourage users to check
the CHANGELOG for breaking changes.




## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/latest/user-guide/) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit our website](https://labs.epi2me.io/wfindex).

### Workflow options

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-basecalling --help
```

to see the options and running examples for the workflow.

### Example command

```
nextflow run epi2me-labs/wf-basecalling \
    -profile singularity
    --input /path/to/my/fast5 \
    --dorado_ext fast5 \
    --ref /path/to/my/ref.fa \
    --out_dir /path/to/my/outputs \
    --basecaller_cfg "dna_r10.4.1_e8.2_400bps_hac@v4.0.0" \
    --basecaller_basemod_threads 2 \
    --remora_cfg "dna_r10.4.1_e8.2_400bps_hac@v4.0.0_5mCG@v2"
```

### Choosing a model

The `dorado` repository has [a table of available models](https://github.com/nanoporetech/dorado#available-basecalling-models) to choose for `--basecaller_cfg` and `--remora_cfg`.

### Duplex calling

wf-basecalling supports [duplex calling](https://github.com/nanoporetech/dorado#duplex), which is enabled with the `--duplex` option. If you used a chemistry and flowcell combination that supported duplex reads, you should switch this option on. 

There are some caveats to duplex calling, namely:
* Duplex mode cannot be used when calling modified bases. You must either run simplex basecalling with modified bases; or duplex calling without modified bases.
* Duplex mode with wf-basecalling is reliant on internal optimisations to organise input files for better duplex rates, which is not possible when using streaming basecalling; therefore duplex combined with the `--watch_path` option could lead to lower duplex rates than what would be achieved running the algorithm after sequencing is completed.

### Updating the workflow

It is recommended to keep this workflow updated to take advantage of the latest basecalling models with:

```
nextflow pull epi2me-labs/wf-basecalling
```

### Workflow outputs

The primary outputs of the workflow include:

* two sorted, indexed CRAMs of basecalls, aligned to the provided reference, with reads separated by their quality score
    * `<sample_name>.pass.cram` contains reads with `qscore >= threshold`
    * `<sample_name>.fail.cram` contains reads with `qscore < threshold`

If `duplex` calling is enabled, the workflow will generate four CRAMs instead:
* `<sample_name>.pass.duplex.cram` contains duplex reads + simplex reads not beloging to a pair with `qscore >= threshold`
* `<sample_name>.fail.duplex.cram` contains duplex reads + simplex reads not beloging to a pair with `qscore < threshold`
* `<sample_name>.pass.simplex.cram` contains simplex reads beloging to a pair with `qscore >= threshold`
* `<sample_name>.fail.simplex.cram` contains simplex reads beloging to a pair with `qscore < threshold`

Take care to retain the input reference as CRAM files cannot be read without the corresponding reference!


### Support for basecalling on GPU

This section will be kept up to date with latest advice for running our workflows on the GPU.

#### Prerequisites

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




## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [singularity](https://sylabs.io/singularity/)
* [dorado](https://github.com/nanoporetech/dorado/)
