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
    """
    FILES=(${fastq_files})
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
    
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple run, file(trimmed)

    output:
    tuple run, file("${run}_dedup.bam"), emit: bamfile
    tuple run, file("**.metrics"), emit: report

    script:
    // I hate the fact that I have to fall back on bash to check the length of this list.
    """
    FILES=(${trimmed})
    case \${#FILES[@]} in
        1)
        bowtie2 -x ${params.genome} -U \${FILES[0]} -S ${run}.sam\n
        ;;

        2)
        bowtie2 -x ${params.genome} -1 \${FILES[0]} -2 \${FILES[1]} -S ${run}.sam\n
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
    get:
    fastq_files

    main:
    fastq_files | createBam
    createBam.out.bamfile | index

    emit:
    bamfile = index.out
    report = createBam.out.report
}
