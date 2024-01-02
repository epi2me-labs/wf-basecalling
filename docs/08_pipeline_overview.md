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

### 4. Duplex calling

wf-basecalling supports [duplex calling](https://github.com/nanoporetech/dorado#duplex), which is enabled with the `--duplex` option. If you used a chemistry and flowcell combination that supported duplex reads, you should switch this option on. The resulting BAM/CRAM will quality filtered and then automatically split in separate BAM/CRAM files for the simplex and duplex reads.
Since `dorado duplex` requires the inputs to be in `pod5` format, the workflow will perform the conversion automatically using [pod5 convert fast5](https://github.com/nanoporetech/pod5-file-format/blob/master/python/pod5/README.md#pod5-convert-fast5). These files are normally deleted upon completion of the analysis, but can optionally be saved by the user by providing the `--output_pod5` option.

### 5. Real-time analysis

wf-basecalling can perform the basecalling as the pod5 files are generated. To enable this, provide the `--watch_path` option. The workflow will process the newly generated files as soon as they become available.
