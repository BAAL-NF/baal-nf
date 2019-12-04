#!/usr/bin/env nextflow

params.staging_root = "${workflow.launchDir}/test/staging/runs"
params.experiments = "${params.staging_root}/runs.csv"
params.report_dir = "${workflow.launchDir}/test/reports/"
params.genome = "hg19"
params.picard_cmd = "picard"

Channel
    .fromPath(params.experiments)
    .splitCsv(header:true)
    .map({row -> tuple(row.cell_line, row.transcription_factor, row.path )})
    .flatMap({cell_line, transcription_factor, path -> file("${params.staging_root}/${path}/experiments/*", type: "dir")})
    .flatMap({folder -> file("${folder}/*/*.fastq.gz")})
    .map { file ->
        def sampleID = file.getParent().getName()
        def filename = file.name
        return tuple(sampleID, filename, file)
    }
    .into{ fastqc_input; for_trim_galore }

process fastQC {
    label 'fastq'
    publishDir("${params.report_dir}/${sampleID}", mode: "move")

    input:
    tuple sampleID, filename, file("${filename}") from fastqc_input

    output:
    file "*_fastqc*"

    script:
    """
    fastqc ${filename}
    """
}

// group multiple fastq files for the same run, since they will need
// to be processed together by trimGalore and subsequent tasks
trim_galore_input = for_trim_galore.groupTuple()

process trimGalore {
    label 'fastq'
    publishDir("${params.report_dir}/${sampleID}", mode: "move", pattern: "*report*")

    input:
    set sampleID, filename, file(files) from trim_galore_input

    output:
    set sampleID, file("${sampleID}_tg*") into trim_galore_out
    file '*report*'

    script:
    switch (files) {
         case nextflow.processor.TaskPath:
         return """trim_galore --basename ${sampleID}_tg ${files}"""

         case nextflow.util.BlankSeparatedList:
         return """trim_galore --basename ${sampleID}_tg --paired ${files}"""

         default:
         println("Error getting files for sample ${sampleID}, exiting")
         return "exit 1"
    }
}

trim_galore_out.into {
    fastq_screen_input
    bowtie_input
}

process fastqScreen {
    label 'fastq'
    publishDir("${params.report_dir}/${sampleID}", mode: "move", pattern: "*screen*")

    input:
    set sampleID, file(trimmed) from fastq_screen_input

    output:
    file trimmed into fastq_screen_out
    file '*screen*'

    script:
    switch (trimmed) {
        case nextflow.processor.TaskPath:
        return """fastq_screen ${trimmed}"""

        case nextflow.util.BlankSeparatedList:
        return """fastq_screen --paired ${trimmed}"""

        default:
        println("Error getting files for sample ${sampleID}, exiting")
        return "exit 1"
    }

}

process align {
    label 'fastq'
    memory '8GB'

    input:
    set sampleID, file(trimmed) from bowtie_input

    output:
    set sampleID, file("${sampleID}.sam") into aligned_sequences

    script:
    switch (trimmed) {
        case nextflow.processor.TaskPath:
        return """bowtie2 -x ${params.genome} -U ${trimmed} -S ${sampleID}.sam"""

        case nextflow.util.BlankSeparatedList:
        first = trimmed[0]
        second = trimmed[1]
        return """bowtie2 -x ${params.genome} -1 ${first} -2 ${second} -S ${sampleID}.sam"""

        default:
        println("Error getting files for sample ${sampleID}, exiting")
        return "exit 1"
    }
}

process sortAndCompress {
    label 'fastq'
    input:
    set sampleID, file(samfile) from aligned_sequences

    output:
    set sampleID, file("${samfile.baseName}.bam") into bamfiles

    """
    samtools sort ${samfile} -o ${samfile}.sorted
    samtools view -h -S -b ${samfile}.sorted > ${samfile.baseName}.bam
    """
}

process markDuplicates {
    label 'fastq'
    publishDir("${params.report_dir}/${sampleID}", mode: "move", pattern: "**.metrics")

    input:
    set sampleID, file(bamfile) from bamfiles

    output:
    set sampleID, file("${bamfile.baseName}_dedup.bam") into marked_bamfiles
    file "${bamfile.baseName}.metrics"

    """
    ${params.picard_cmd} MarkDuplicates I="${bamfile}" O="${bamfile.baseName}_dedup.bam" M="${bamfile.baseName}.metrics"
    """
}

process index {
    label 'fastq'
    input:
    set sampleID, file(bamfile) from marked_bamfiles

    output:
    set sampleID, file(bamfile), file("${bamfile}.bai") into indexed_bamfiles

    """
    samtools index ${bamfile}
    """
}

