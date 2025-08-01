
// this process is shared by both the uCRAM and CRAM arms of the basecalling workflow
// for uCRAM the staged ref is the OPTIONAL_FILE, so we withhold the ref arg
process merge_calls {
    label "wf_basecalling"
    cpus params.merge_threads
    memory 16.GB
    input:
        path(ref)
        path(crams, stageAs: "filtered_*.cram")
        val(filetag)
        tuple val(align_ext), val(index_ext) // either [bam, bai] or [cram, crai]
    output:
        tuple path("${out_file}"), path("${out_index}")
    script:
    out_file = "${params.sample_name}.${filetag}.${align_ext}"  // nodef used in output 
    out_index = "${out_file}.${index_ext}"  // nodef used in output 
    def ref_arg = ref.name != "OPTIONAL_FILE" ? "--reference \"${ref}\"" : ""
    """
    samtools merge -c "${out_file}##idx##${out_index}" ${crams} --no-PG -O ${align_ext} --write-index ${ref_arg} --threads ${task.cpus}
    """
}

process merge_calls_to_fastq {
    label "wf_basecalling"
    cpus { params.merge_threads + params.ubam_bam2fq_threads }
    memory 16.GB
    input:
        path(crams, stageAs: "filtered_*.cram")
        val(filetag)
    output:
        path(fq_file)
    script:
    fq_file = "${params.sample_name}.${filetag}.fq.gz"  // nodef used in output 
    """
    samtools merge -c ${crams} --no-PG -O CRAM -@ ${params.merge_threads} -o - | samtools bam2fq -T 1 -@ ${params.ubam_bam2fq_threads} -0 "${fq_file}" -
    """
}
