#!/usr/bin/env nextflow
import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList

nextflow.enable.dsl = 2

include { wf_dorado } from './lib/signal/ingress'
nextflow.preview.recursion=true

process getVersions {
    label "wf_basecalling"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    dorado --version 2>&1 | head -n1 | sed 's/^/dorado,/' >> versions.txt
    minimap2 --version | head -n 1 | sed 's/^/minimap2,/' >> versions.txt
    """
}


process getParams {
    label "wf_basecalling"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process cram_cache {
    label "wf_common"
    input:
        path reference
    output:
        path("ref_cache/"), emit: cram_cache
    script:
    """
    seq_cache_populate.pl -root ref_cache/ ${reference}
    """
}


process bamstats {
    label "wf_common"
    cpus params.stats_threads
    input:
        path "input.cram" // chunks are always CRAM
        path ref_cache

    output:
        path "bamstats.tsv", emit: stats
        path "stats.${task.index}.json", emit: json
    script:
        String ref_path = ref_cache.name.startsWith('OPTIONAL_FILE') ? '' : "export REF_PATH=${ref_cache}/%2s/%2s/%s"
    """
    ${ref_path}
    bamstats --threads=${task.cpus} -u input.cram > bamstats.tsv
    fastcat_histogram.py \
            --sample_id "${params.sample_name}" \
            bamstats.tsv "stats.${task.index}.json"
    """
}


// Scan step for accumulating fastcat stats
//
// Nextflow scan does a silly thing where it feeds back the growing list of
// historical outputs. We only ever need the most recent output (the "state").
process progressive_stats {
    label "wf_common"
    maxForks 1
    cpus 1
    input: 
        path fastcat_stats
    output:
        path("all_stats.${task.index}")
    script:
        def new_input = fastcat_stats instanceof BlankSeparatedList ? fastcat_stats.first() : fastcat_stats
        def state = fastcat_stats instanceof BlankSeparatedList ? fastcat_stats.last() : "NOSTATE"
        def output = "all_stats.${task.index}"
    """
    touch "${state}"
    add_jsons.py "${new_input}" "${state}" "${output}"
    """
}


// Split simplex reads belonging to a pair
process split_xam {
    label "wf_common"
    cpus 2
    input:
        tuple path(xam), path(xam_index)
        tuple val(align_ext), val(index_ext)
        path ref
    output:
        tuple path("${xam.baseName}.duplex.${align_ext}"), path("${xam.baseName}.duplex.${align_ext}.${index_ext}"), emit: xam_dx
        tuple path("${xam.baseName}.simplex.${align_ext}"), path("${xam.baseName}.simplex.${align_ext}.${index_ext}"), emit: xam_sx
    script:
        String reference = ref.name.startsWith('OPTIONAL_FILE') ? '' : "--reference ${ref}"
    """
    samtools view ${reference} \
        -@ ${task.cpus} \
        -O ${align_ext} \
        --tag dx:-1 \
        --unoutput ${xam.baseName}.duplex.${align_ext} \
        -o ${xam.baseName}.simplex.${align_ext} \
        ${xam}
    samtools index ${xam.baseName}.simplex.${align_ext}
    samtools index ${xam.baseName}.duplex.${align_ext}
    """
}


// Compute pairing statistics progressively, if duplex enabled
process pair_stats {
    label "wf_common"
    cpus 1
    input:
        path cram // chunks are always CRAM
        path ref_cache
    output:
        path("pairs.${task.index}.csv"), emit: csv
    script:
        String ref_path = ref_cache.name.startsWith('OPTIONAL_FILE') ? '' : "export REF_PATH=${ref_cache}/%2s/%2s/%s"
    """
    ${ref_path}
    duplex_stats.py ${cram} pairs.${task.index}.csv
    """
}


process progressive_pairings {
    label "wf_common"
    maxForks 1
    cpus 1
    input: 
        path pairings
    output:
        path("pairing_stats.${task.index}")
    script:
        def new_input = pairings instanceof BlankSeparatedList ? pairings.first() : pairings
        def state = pairings instanceof BlankSeparatedList ? pairings.last() : "NOSTATE"
        def output = "pairing_stats.${task.index}"
    """
    cat ${new_input} > pairing_stats.${task.index}
    if [ -e ${state} ]; then
        awk 'NR>1 {print \$0}' ${state} >> pairing_stats.${task.index}
    fi
    """
}


// Make reports
process makeReport {
    label "wf_common"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path per_read_stats
        path pairings
        path "versions/*"
        path "params.json"
    output:
        path "wf-basecalling-*.html"
    script:
        String report_name = "wf-basecalling-report.html"
        def report_pairings = params.duplex ? "--pairings ${pairings}" : ""
    """
    report.py $report_name \
        --sample_name $params.sample_name \
        --versions versions \
        --stats $per_read_stats \
        --params params.json \
        $report_pairings
    """
}


// watch path stop condition, if params.read_limit is met will inject a stop file in to input folder.
process stopCondition { 
    label "wf_common"
    cpus 1 
    publishDir params.input, mode: 'copy', pattern: "*"
    input:
        path json
        val (stop_filename)
    output:
        path "${stop_filename}", optional: true, emit: stop
    script:
        int threshold = params.read_limit
    """    
    #!/usr/bin/env python
    import json
    from pathlib import Path
    with open("$json") as json_file:
        state = json.load(json_file)
        total = 0 
        for k,v in state.items():
            total += v["total_reads"]
        if total >= $threshold:
            p = Path("$stop_filename")
            p.touch(exist_ok=False)
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_stream {
    // publish inputs to output directory
    label "wf_basecalling"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files."
    """
}

// Output the last report, once each of them finish
process output_last {
    // publish inputs to output directory
    label "wf_basecalling"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files."
    """
}

// CW-2569: Emit pod5s if requested, in a new directory
process output_pod5s {
    // publish inputs to output directory
    label "wf_basecalling"
    publishDir "${params.out_dir}/pod5s/", mode: 'copy', pattern: "*"
    input:
        path pod5s
    output:
        path pod5s
    """
    echo "Writing output files."
    """
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Map colors = NfcoreTemplate.logColours(params.monochrome_logs)

    Pinguscript.ping_start(nextflow, workflow, params)

    // Basecall
    // Ensure basecaller config is set
    if (!params.basecaller_cfg && !params.basecaller_model_path) {
        throw new Exception(colors.red + "You must provide a basecaller profile with --basecaller_cfg <profile>" + colors.reset)
    }
    if (params.basecaller_cfg && params.basecaller_model_path) {
        log.warn("--basecaller_cfg and --basecaller_model_path both provided. Custom remora model path (${params.basecaller_cfg}) will override enum choice (${params.basecaller_model_path}).")
    }
    if (params.remora_cfg && params.remora_model_path) {
        log.warn("--remora_cfg and --remora_model_path both provided. Custom remora model path (${params.remora_model_path}) will override enum choice (${params.remora_cfg}).")
    }
    if (params.duplex && params.dorado_ext != 'pod5') {
        log.warn("Duplex currently requires POD5 files and is not compatible with FAST5. The workflow will convert the FAST5 inputs to POD5 format automatically.")
    }

    // Ensure modbase threads are set if calling them
    if (params.remora_cfg || params.remora_model_path) {
        if (params.duplex) {
            throw new Exception(colors.red + "--duplex cannot call modified bases.\nUnset --remora_cfg/--remora_model_path to run duplex basecalling, or unset --duplex to run simplex modified basecalling." + colors.reset)
        }
        if (params.basecaller_basemod_threads == 0) {
            throw new Exception(colors.red + "--remora_cfg modbase aware config requires setting --basecaller_basemod_threads > 0" + colors.reset)
        }
    }
    // ring ring it's for you
    basecaller_out = wf_dorado([
        "input_path": params.input,
        "input_ref": params.ref,
        "basecaller_model_name": params.basecaller_cfg,
        "remora_model_name": params.remora_cfg,
        "basecaller_model_path": params.basecaller_model_path,
        "remora_model_path": params.remora_model_path,
        "watch_path": params.watch_path,
        "output_bam": params.output_bam,
        "dorado_ext": params.dorado_ext,
        "fastq_only": params.fastq_only,
    ])
    software_versions = getVersions()
    workflow_params = getParams()

    if (params.ref) {
        ref = Channel.fromPath(params.ref, checkIfExists: true).first()
    } else {
        ref = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE")
    }

    // create cram ref cache if there is a ref (basecaller always emit cram)
    if (params.ref) {
        ref_cache = cram_cache(ref)
    }
    else {
        ref_cache = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE")
    }

    // stream stats for report
    stat = bamstats(basecaller_out.chunked_pass_crams, ref_cache)
    stats = progressive_stats.scan(stat.json)
    // stream pair stats for report
    pairings = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE")
    if (params.duplex){
        // Separate the simplex reads belonging to a pair from the
        // duplex and simplex reads.
        // Save the simplex reads in a duplex in a separate xam file.
        split_xam(
            basecaller_out.pass.concat(
                basecaller_out.fail
            ),
            basecaller_out.output_exts,
            ref
            )

        // Create emission channel
        emit_xam = split_xam.out.xam_dx.flatten()
            .concat(split_xam.out.xam_sx.flatten())
            .concat(basecaller_out.summary)

        // Then, compute the stats on the duplex
        pairs = pair_stats(basecaller_out.chunked_pass_crams, ref_cache)
        pairings = progressive_pairings.scan(pairs.csv)
    } else {
        emit_xam = basecaller_out.pass.flatten()
            .concat(basecaller_out.fail.flatten())
    }
    // Make the report
    report = makeReport(stats, pairings, software_versions, workflow_params) | last | collect | output_last

    // dump out artifacts thanks for calling
    output_stream(emit_xam.concat(pairings.last(), software_versions, workflow_params, ref, ref_cache))

    // dump pod5s if requested
    if (params.duplex && params.dorado_ext == 'fast5' && params.output_pod5){
        output_pod5s(basecaller_out.converted_pod5s)
    }

    //  Stop file to input folder when read_limit stop condition is met.
    String stop_filename = "STOP.${workflow.sessionId}.${params.dorado_ext}"
    if (params.watch_path && params.read_limit){
        stopCondition(stats, stop_filename).first().subscribe {
            log.info "Creating STOP file: '$stop_filename'"
        }
    }

}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
