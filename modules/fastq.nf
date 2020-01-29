params.genome = "hg19"
params.report_dir = "${workflow.launchDir}/test/reports/"
params.picard_cmd = "picard"

process trimGalore {
    label 'fastq'

    input:
    tuple run, file(fastq_files)

    output:
    tuple run, file("${run}_tg*"), emit: trimmed_fastq
    tuple run, file('*report*'), emit: report

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

process createBam {
    label 'fastq'
    label 'bigmem'

    input:
    tuple run, file(trimmed)

    output:
    tuple run, file("${run}_dedup.bam"), emit: bamfile
    tuple run, file("**.metrics"), emit: report

    script:
    result = ""
    switch (trimmed) {
        case nextflow.processor.TaskPath:
        result += """bowtie2 -x ${params.genome} -U ${trimmed} -S ${run}.sam\n"""
        break

        case nextflow.util.BlankSeparatedList:
        first = trimmed[0]
        second = trimmed[1]
        result += """bowtie2 -x ${params.genome} -1 ${first} -2 ${second} -S ${run}.sam\n"""
        break

        default:
        println("Error getting files for sample ${run}, exiting")
        return "exit 1"
    }
    result += """samtools sort ${run}.sam -o ${run}.sam.sorted
    samtools view -h -S -b ${run}.sam.sorted > ${run}.bam
    ${params.picard_cmd} MarkDuplicates I="${run}.bam" O="${run}_dedup.bam" M="${run}.metrics"
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
    fastq_files | createBam
    createBam.out.bamfile | index

    emit:
    bamfile = index.out
    report = createBam.out.report
}
