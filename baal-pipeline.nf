#!/usr/bin/env nextflow

params.staging_root = "${workflow.launchDir}/test/staging/runs"
params.experiments = "${params.staging_root}/runs.csv"
params.report_dir = "${workflow.launchDir}/test/reports/"

Channel
    .fromPath(params.experiments)
    .splitCsv(header:true)
    .map({row -> tuple(row.cell_line, row.transcription_factor, row.path )})
    .flatMap({cell_line, transcription_factor, path -> file("${params.staging_root}/${path}/experiments/*", type: "dir")})
    .flatMap({path -> file("${path}/*/*.fastq.gz")})
    .map { file ->
        def sampleID = file.getParent().getName()
        return tuple(sampleID, file)
    }
    .into{ fastqc_input; trim_galore_input }


process fastQC {
    publishDir("${params.report_dir}/${sampleID}", mode: "move")

    input:
    set sampleID, file("${sampleID}.fastq.gz") from fastqc_input

    output:
    file "${sampleID}_fastqc*" into fastqc_out

    script:
    """
    fastqc ${sampleID}.fastq.gz
    """
}



process trimGalore {
    publishDir("${params.report_dir}/${sampleID}", mode: "move", pattern: "*report*")

    input:
    set sampleID, file("${sampleID}.fastq.gz") from trim_galore_input

    output:
    file '*trimmed*' into trim_galore_out
    file '*report*'

    script:
    """
    trim_galore ${sampleID}.fastq.gz
    """
}
