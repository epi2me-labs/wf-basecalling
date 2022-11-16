## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://sylabs.io/singularity/) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

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
    --basecaller_cfg "dna_r10.4.1_e8.2_400bps_hac@v3.5.2" \
    --basecaller_basemod_threads 2 \
    --remora_cfg "dna_r10.4.1_e8.2_400bps_hac@v3.5.2_5mCG@v2"
```

### Choosing a model

The `dorado` repository has [a table of available models](https://github.com/nanoporetech/dorado#available-basecalling-models) to choose for `--basecaller_cfg` and `--remora_cfg`.

### Updating the workflow

It is strongly recommended to keep this experimental workflow updated with:

```
nextflow pull epi2me-labs/wf-basecalling
```

Users are reminded that `dorado` is released for evaluation purposes only.
Users should consult the CHANGELOG to keep up to date with breaking changes.

### Workflow outputs

The primary outputs of the workflow include:

* a sorted, indexed CRAM of basecalls, aligned to the provided reference

