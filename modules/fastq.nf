params.genome = "hg19"
params.report_dir = "${workflow.launchDir}/test/reports/"
params.picard_cmd = "picard"

process trimGalore {
    label 'fastq'
    publishDir("${params.report_dir}/${run}", mode: "move", pattern: "*report*")

    input:
    tuple run, file(fastq_files)

    output:
    tuple run, file("${run}_tg*"), emit: trimmed_fastq
    file '*report*'

    script:
    switch (fastq_files) {
         case nextflow.processor.TaskPath:
         return """trim_galore --basename ${run}_tg ${fastq_files}"""

         case nextflow.util.BlankSeparatedList:
         return """trim_galore --basename ${run}_tg --paired ${fastq_files}"""

         default:
         println("Error getting files for sample ${run}, exiting")
         return "exit 1"
    }
}

process align {
    label 'fastq'
    label 'bigmem'

    input:
    tuple run, file(trimmed)

    output:
    tuple run, file("${run}.sam")

    script:
    switch (trimmed) {
        case nextflow.processor.TaskPath:
        return """bowtie2 -x ${params.genome} -U ${trimmed} -S ${run}.sam"""

        case nextflow.util.BlankSeparatedList:
        first = trimmed[0]
        second = trimmed[1]
        return """bowtie2 -x ${params.genome} -1 ${first} -2 ${second} -S ${run}.sam"""

        default:
        println("Error getting files for sample ${run}, exiting")
        return "exit 1"
    }
}

process sortAndCompress {
    label 'fastq'
    input:
    tuple run, file(samfile)

    output:
    tuple run, file("${samfile.baseName}.bam")

    """
    samtools sort ${samfile} -o ${samfile}.sorted
    samtools view -h -S -b ${samfile}.sorted > ${samfile.baseName}.bam
    """
}

process markDuplicates {
    label 'fastq'
    publishDir("${params.report_dir}/${run}", mode: "move", pattern: "**.metrics")

    input:
    tuple run, file(bamfile)

    output:
    tuple run, file("${bamfile.baseName}_dedup.bam"), emit: bam_files
    file "${bamfile.baseName}.metrics"

    """
    ${params.picard_cmd} MarkDuplicates I="${bamfile}" O="${bamfile.baseName}_dedup.bam" M="${bamfile.baseName}.metrics"
    """
}

process index {
    label 'fastq'
    input:
    tuple run, file(bamfile)

    output:
    tuple run, file(bamfile), file("${bamfile}.bai")

    """
    samtools index ${bamfile}
    """
}

workflow create_bam {
    get:
    fastq_files

    main:
    fastq_files | align | sortAndCompress | markDuplicates
    markDuplicates.out.bam_files | index

    emit:
    index.out
}
