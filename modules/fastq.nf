params.genome = "genome"
params.report_dir = "${workflow.launchDir}/reports/"
params.picard_cmd = "picard"

process trimGalore {
    label 'fastq'

    input:
    tuple run, file("${run}*.fastq.gz")

    output:
    tuple run, file("${run}_tg*"), emit: trimmed_fastq
    tuple run, file('*report*'), emit: report

    script:
    """
    FILES=(${run}*.fastq.gz)
    case \${#FILES[@]} in
        1)
        trim_galore --basename ${run}_tg \${FILES[0]}
        ;;

        2)
        trim_galore --basename ${run}_tg --paired \${FILES[0]} \${FILES[1]}
        ;;

        *)
        echo "Error getting files for sample ${run}, exiting"
        exit 1
    esac
    """
}

process createBam {
    label 'fastq'
    label 'bigmem'
    publishDir("${params.report_dir}/logs/bowtie2/", mode: "copy", pattern: "${run}.log")    
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple run, file(trimmed)
    file index_files

    output:
    tuple run, file("${run}_dedup.bam"), emit: bamfile
    tuple run, file("**.metrics"), emit: report
    file "${run}.log"

    script:
    // I hate the fact that I have to fall back on bash to check the length of this list.
    """
    export BOWTIE2_INDEXES=${index_files}
    FILES=(${trimmed})
    case \${#FILES[@]} in
        1)
        bowtie2 -x ${params.genome} -U \${FILES[0]} -S ${run}.sam 2> >(tee ${run}.log >&2)\n
        ;;

        2)
        bowtie2 -x ${params.genome} -1 \${FILES[0]} -2 \${FILES[1]} -S ${run}.sam 2> >(tee ${run}.log >&2)\n
        ;;

        *)
        echo Error getting files for sample ${run}, exiting
        exit 1
    esac
    samtools sort ${run}.sam -o ${run}.sam.sorted
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
    take:
    fastq_files

    main:
    index_ch = file(params.bowtie2_index, type: "dir")
    createBam(fastq_files, index_ch)
    createBam.out.bamfile | index

    emit:
    bamfile = index.out
    report = createBam.out.report
}
