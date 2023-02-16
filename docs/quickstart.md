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
Please note that wf-basecalling does not currently support duplex stereo basecalling.

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

Take care to retain the input reference as CRAM files cannot be read without the corresponding reference!

