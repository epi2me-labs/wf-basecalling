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
| fastq_only | boolean | Output unaligned FASTQ instead of unaligned CRAM. | FASTQ can only be output when a reference has not been provided. Aligned output will always be written to CRAM even if fastq_only is set. | False |
| output_bam | boolean | Output unaligned BAM instead of unaligned CRAM. | Some downstream applications do not yet support CRAM and will require a BAM file. Enabling this option will output BAM instead of CRAM. You should only use this option if you know that it is needed. Output files will be larger than the corresponding CRAM files that would have been written if this option was not enabled. | False |


### Basecalling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| basecaller_cfg | string | Name of the model to use for converting signal. | Required for basecalling. The model list only shows models that are compatible with this workflow. |  |
| duplex | boolean | Run the basecaller in duplex mode. | By default, the workflow conducts simplex basecalling. If you used a chemistry and flowcell combination that supported duplex reads, you should switch this option on. Currently, duplex basecalling is not compatible with modified basecalling. Additionally, duplex basecalling within this workflow is reliant on internal optimisations to organise input files for better duplex rates, which is not possible when using streaming basecalling; therefore duplex combined with the watch_path option could lead to lower duplex rates. | False |
| remora_cfg | string | Name of the model to use for calling modified bases. | Required for calling modified bases while basecalling. The model list only shows models that are compatible with this workflow. |  |
| dorado_ext | string | File extension for Dorado inputs. | Set this to fast5 if you have not converted your fast5 to pod5. It is recommended to [convert existing fast5 files to pod5 for use with Dorado](https://github.com/nanoporetech/pod5-file-format/blob/master/python/README.md#pod5-convert-from-fast5). | pod5 |


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


### Generic options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| disable_ping | boolean | Disable workflow ping. |  | False |


