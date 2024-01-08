import ArgumentParser
import java.lang.Math
import java.util.regex.Pattern
import java.util.regex.Matcher

include {
    merge_calls as merge_pass_calls;
    merge_calls as merge_fail_calls;
    merge_calls_to_fastq as merge_pass_calls_to_fastq;
    merge_calls_to_fastq as merge_fail_calls_to_fastq;
} from './merge.nf'  // sigh. module aliasing seems easier than flow.


Map parse_arguments(Map arguments) {
    ArgumentParser parser = new ArgumentParser(
        args:[
            "input_path",
            "input_ref",
            "basecaller_model_name",
            "remora_model_name",
        ],
        kwargs:[
            "basecaller_model_path": null,
            "remora_model_path": null,
            "watch_path": false,
            "dorado_ext": "pod5",
            "output_bam": false,
            "fastq_only": false
        ],
        name: "signal_ingress")
    return parser.parse_args(arguments)
}

process dorado {
    label "wf_dorado"
    label "wf_basecalling"
    label "gpu"
    accelerator 1 // further configuration should be overloaded using withLabel:gpu
    cpus 8
    memory 48.GB // 48 should be plenty, if needed can push for more based on instance type: p3=61 GB RAM per GPU, p4=144 GB RAM per GPU
    input:
        tuple val(chunk_idx), path('*')
        tuple val(basecaller_cfg), path("dorado_model"), val(basecaller_model_override)
        tuple val(remora_cfg), path("remora_model"), val(remora_model_override)
    output:
        path("${chunk_idx}.ubam"), emit: ubams
        path("converted/*.pod5"), emit: converted_pod5s, optional: true
    script:
    def remora_model = remora_model_override ? "remora_model" : "\${DRD_MODELS_PATH}/${remora_cfg}"
    def remora_args = (params.basecaller_basemod_threads > 0 && (params.remora_cfg || remora_model_override)) && !params.duplex ? "--modified-bases-models ${remora_model}" : ''
    def model_arg = basecaller_model_override ? "dorado_model" : "\${DRD_MODELS_PATH}/${basecaller_cfg}"
    def basecaller_args = params.basecaller_args ?: ''
    def caller = params.duplex ? "duplex" : "basecaller"
    // CW-2569: delete pod5 is not required them to be emitted
    def signal_path = (params.duplex && params.dorado_ext == 'fast5') ? "converted/" : "."
    def delete_pod5s = !params.output_pod5 ? "rm -r converted/" : "echo 'No cleanup'"
    // we can't set maxForks dynamically, but we can detect it might be wrong!
    if (task.executor != "local" && task.maxForks == 1) {
        log.warn "Non-local workflow execution detected but GPU tasks are currently configured to run in serial, perhaps you should be using '-profile discrete_gpus' to parallelise GPU tasks for better performance?"
    }
    // If no pairs list, run vanilla duplex
    """
    # CW-2569: convert the pod5s contextually
    if [[ "${params.duplex}" == "true" && "${params.dorado_ext}" == "fast5" ]]; then
        pod5 convert fast5 ./*.fast5 --output ${signal_path} --threads ${task.cpus} --one-to-one ./
    fi

    # Run dorado on the new pod5s
    dorado ${caller} \
        ${model_arg} \
        ${signal_path} \
        ${remora_args} \
        ${basecaller_args} \
        --device ${params.cuda_device} | samtools view --no-PG -b -o ${chunk_idx}.ubam -

    # CW-2569: delete the pod5s, if emit not required.
    if [[ "${params.duplex}" == "true" && "${params.dorado_ext}" == "fast5" ]]; then
        ${delete_pod5s}
    fi
    """
}


process align_and_qsFilter {
    label "wf_basecalling"
    cpus {params.ubam_map_threads + params.ubam_sort_threads + params.ubam_bam2fq_threads}
    memory 32.GB
    input:
        path mmi_reference
        path reference
        path reads
    output:
        // NOTE merge does not need an index if merging with region/BED (https://github.com/samtools/samtools/blob/969d44990df7fa9c7bda3a7140a2c1d1bd8c62a0/bam_sort.c#L1256-L1272)
        // so we can save a few cycles and just output the CRAM
        path("${reads.baseName}.pass.cram"), emit: pass
        path("${reads.baseName}.fail.cram"), emit: fail
    script:
    """
    samtools bam2fq -@ ${params.ubam_bam2fq_threads} -T 1 ${reads} \
        | minimap2 -y -t ${params.ubam_map_threads} -ax map-ont ${mmi_reference} - \
        | samtools sort -@ ${params.ubam_sort_threads} \
        | samtools view -e '[qs] >= ${params.qscore_filter}' --output ${reads.baseName}.pass.cram --unoutput ${reads.baseName}.fail.cram -O CRAM --reference ${reference} -
    """
}


process qsFilter {
    label "wf_basecalling"
    cpus 2
    memory 16.GB
    input:
        path reads
    output:
        path("${reads.baseName}.pass.cram"), emit: pass
        path("${reads.baseName}.fail.cram"), emit: fail
    script:
    """
    samtools view -e '[qs] >= ${params.qscore_filter}' ${reads} --output ${reads.baseName}.pass.cram --unoutput ${reads.baseName}.fail.cram -O CRAM
    """
}


process make_mmi {
    label "wf_basecalling"
    cpus 4
    memory 16.GB
    input:
        path(ref)
    output:
        path("ref.mmi")
    script:
    """
    minimap2 -t ${task.cpus} -x map-ont -d ref.mmi ${ref}
    """
}



// Compute summaries from the raw unmapped bam files
process dorado_summary {
    label "wf_basecalling"
    cpus 1
    input:
        path xam // chunks are always CRAM
    output:
        path("${xam.baseName}.summary.tsv.gz"), emit: summary
    script:
    """
    dorado summary ${xam} | gzip -c > ${xam.baseName}.summary.tsv.gz
    """
}


// Create list of simplex reads
process combine_dorado_summaries {
    label "wf_basecalling"
    cpus 1
    input:
        path tsvs // individual summaries
    output:
        path("${params.sample_name}.summary.tsv.gz"), emit: summary
    shell:
    '''
    n=0
    for fname in !{tsvs}; do
        if [ $n -eq 0 ]; then
            n=1
            zcat $fname
        else
            zcat $fname | awk 'NR>1 {print}'
        fi
    done | gzip -c > !{params.sample_name}.summary.tsv.gz
    '''
}


workflow wf_dorado {
    take:
        arguments
    main:
        Map margs = parse_arguments(arguments)

        // determine output extentions
        def align_ext = "cram"
        def index_ext = "crai"
        if (margs.output_bam) {
            align_ext = "bam"
            index_ext = "bai"
        }
        output_exts = Channel.of([align_ext, index_ext]).collect()

        if (margs.input_ref) {
            if (margs.fastq_only) {
                log.warn "Ignoring request to output FASTQ as you have provided a reference for alignment."
            }
            // create value channel of ref by calling first
            ref = Channel.fromPath(margs.input_ref, checkIfExists: true).first()
        }
        else {
            ref = file("${projectDir}/data/OPTIONAL_FILE")
        }

        // Munge models
        // I didn't want to use the same trick from wf-humvar as I thought the models here are much larger
        // ...they aren't, but nevermind this is less hilarious than the humvar way
        basecaller_model = file("${projectDir}/data/OPTIONAL_FILE")
        def basecaller_model_override = false
        if (margs.basecaller_model_path) {
            basecaller_model = file(margs.basecaller_model_path, type: "dir", checkIfExists: true)
            basecaller_model_override = true
            log.warn "Overriding basecaller model with ${margs.basecaller_model_path}"
        }
        remora_model = file("${projectDir}/data/OPTIONAL_FILE")
        def remora_model_override = false
        if (margs.remora_model_path) {
            remora_model = file(margs.remora_model_path, type: "dir", checkIfExists: true)
            remora_model_override = true
            log.warn "Overriding remora model with ${margs.remora_model_path}"
        }

        Integer chunk_idx = 0
        String stop_filename = "STOP.${workflow.sessionId}.${margs.dorado_ext}" // use the sessionId so resume works
        existing_pod5_chunks = Channel
            .fromPath(margs.input_path + "**.${margs.dorado_ext}", checkIfExists: true)

        // Define naming patterns
        // The default naming is: PAN00000_pass__abc012d4_wxy789z0_0.pod5
        def pattern_long = /[A-Z]{3}[0-9]{5}_(pass|fail)__[0-9a-zA-Z]{8}_[A-Za-z0-9]{8}_[0-9]*/
        // The alternative naming is: PAN00000_abc012d4_0.pod5
        def pattern_short = /[A-Z]{3}[0-9]{5}_[0-9a-zA-Z]{8}_[0-9]*/

        // Check if it is a watched path
        if (margs.watch_path) {
            // watch input path for more pod5s
            if (params.read_limit) {
                log.warn "Watching ${margs.input_path} for new ${margs.dorado_ext} files, until ${params.read_limit} reads have been observed."
                log.warn "To stop this workflow early, create: ${margs.input_path}/${stop_filename}"
            }
            else {
                log.warn "Watching ${margs.input_path} for new ${margs.dorado_ext} files indefinitely."
                log.warn "To stop this workflow create: ${margs.input_path}/${stop_filename}"
            }
            // Warn if duplex requested during streaming
            if (params.duplex){
                log.warn "Streamed duplex calling cannot be optimized, and might lead to lower duplex rates."
            }
            // Define watch path
            watch_pod5_chunks = Channel
                .watchPath("${margs.input_path}/**.${margs.dorado_ext}")
                .until{ it.name == stop_filename }
            // Create list of pod5s to process
            ready_pod5_chunks = existing_pod5_chunks
                .concat(watch_pod5_chunks)
                .buffer(size:params.basecaller_chunk_size, remainder:true)
                .map { tuple(chunk_idx++, it) }
        } else {
            // Check if a pod5 matches the long pattern, the short pattern or
            // no pattern. After that, branch them in separate subchannel to
            // be processed separately to extract the pieces of information needed
            // to define the best grouping.
            existing_pod5_chunks
                .branch{
                    pattern_long: it.baseName ==~ pattern_long
                    pattern_short: it.baseName ==~ pattern_short
                    no_pattern: true
                }.set{branched_patterns}
            // Check if there are files with unknown naming patterns and log that
            if (params.duplex) {
                branched_patterns.no_pattern.first().subscribe { log.warn "Detected input files with an unknown naming pattern. These will be basecalled, but duplex rates may be impacted." }
            }
            // Split the name keeping the flow cell, run, pass/fail and pod5 index.
            branched_patterns.pattern_long
                // If the pattern is long, extract flow cell, run, pass/fail and index
                .map{filename -> 
                    fields = filename.baseName.split("_")
                    [fields[-1] as int, fields[0], fields[3], fields[1], filename]
                }
                .mix(
                    branched_patterns.pattern_short
                        // If the pattern is short, extract flow cell, run and index, and set all files as pass
                        .map{ filename ->
                            fields = filename.baseName.split("_")
                            [fields[-1] as int, fields[0], fields[1], 'pass', filename]
                        }
                )
                // Computed chunk index as floored pod5 index / chunk size value, and then concatenate them. 
                // If not pass is provided, treat all as pass. The chunk number goes first to perform chunk
                // clustering appropriately.
                .map{ pod5_index, cell_id, run_id, pass, pod5 ->
                    [Math.floor(pod5_index/params.basecaller_chunk_size), cell_id, run_id, pass, pod5]
                }
                // Group them by chunk, flowcell, run and pass/fail
                .groupTuple(
                    by: [0, 1, 2, 3]
                )
                // Emit only files for chunking
                .map{pod5_index, cell_id, run_id, pass, pod5s -> pod5s}
                // Add back pod5s not following naming pattern
                .mix(
                    branched_patterns.no_pattern
                        .buffer(size:params.basecaller_chunk_size, remainder: true)
                )
                // Replace the chunk number with a new progressive numbering
                .map { pod5s -> [chunk_idx++, pod5s] }
                // Set the channel
                .set{ready_pod5_chunks}
        }


        called_bams = dorado(
            ready_pod5_chunks,
            tuple(margs.basecaller_model_name, basecaller_model, basecaller_model_override),
            tuple(margs.remora_model_name, remora_model, remora_model_override)
        )

        // Compute summary
        if (params.dorado_ext == 'pod5' && params.duplex){
            dorado_summary(called_bams.ubams) | collect | combine_dorado_summaries    
            summary = combine_dorado_summaries.out.summary
        } else {
            summary = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE")
        }

        // Run filtering or mapping
        if (margs.input_ref) {
            // make mmi for faster alignment
            mmi_ref = make_mmi(ref)

            // align, qscore_filter and sort
            crams = align_and_qsFilter(mmi_ref, ref, called_bams.ubams)
        }
        else {
            // skip alignment and just collate pass and fail
            crams = qsFilter(called_bams.ubams)
        }

        // merge passes and fails
        // we've aliased the merge_calls process to save writing some unpleasant looking flow
        // FASTQ output can only be used when there is no input_ref
        if (margs.fastq_only && !margs.input_ref) {
            pass = merge_pass_calls_to_fastq(crams.pass.collect(), "pass")
            fail = merge_fail_calls_to_fastq(crams.fail.collect(), "fail")
        }
        else {
            pass = merge_pass_calls(ref, crams.pass.collect(), "pass", output_exts)
            fail = merge_fail_calls(ref, crams.fail.collect(), "fail", output_exts)
        }

    emit:
        chunked_pass_crams = crams.pass
        pass = pass
        fail = fail
        output_exts = output_exts
        summary = summary
        converted_pod5s = called_bams.converted_pod5s
}
