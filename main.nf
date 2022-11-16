#!/usr/bin/env nextflow
import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { wf_dorado } from './basecalling'

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


process makeReport {
    label "wf_basecalling"
    input:
        val metadata
        path "seqs.txt"
        path "versions/*"
        path "params.json"
    output:
        path "wf-template-*.html"
    script:
        report_name = "wf-template-" + params.report_name + '.html'
        def metadata = new JsonBuilder(metadata).toPrettyString()
    """
    echo '${metadata}' > metadata.json
    report.py $report_name \
        --versions versions \
        seqs.txt \
        --params params.json \
        --metadata metadata.json
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
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


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }
    if (workflow.profile == "conda") {
        throw new Exception(colors.red + "Sorry, wf-basecalling is not compatible with --profile conda, please use --profile standard (Docker) or --profile singularity." + colors.reset)
    }

    // Basecall
    // Ensure basecaller config is set
    if (!params.basecaller_cfg) {
        throw new Exception(colors.red + "You must provide a basecaller profile with --basecaller_cfg <profile>" + colors.reset)
    }
    // Ensure modbase threads are set if calling them
    if (params.remora_cfg && params.basecaller_basemod_threads == 0) {
        throw new Exception(colors.red + "--remora_cfg modbase aware config requires setting --basecaller_basemod_threads > 0" + colors.reset)
    }

    // ring ring it's for you
    basecaller_out = wf_dorado(params.input, file(params.ref))

    software_versions = getVersions()
    workflow_params = getParams()

    bam = basecaller_out.cram
    idx = basecaller_out.crai

    // dump out artifacts thanks for calling
    output(bam.concat(idx, software_versions, workflow_params))
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
