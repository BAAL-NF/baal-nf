// process fastq files and prepare them for use in BaalChIP
process trimGalore {
    label 'parallel'

    input:
    tuple val(run), path("${run}*.fastq.gz")

    output:
    tuple val(run), path("${run}_tg*"), emit: trimmed_fastq
    tuple val(run), path('*report*'), emit: report

    script:
    extra_args = ''
    if (task.cpus > 1) {
        extra_args += "-j ${task.cpus}"
    }

    """
    FILES=(${run}*.fastq.gz)
    case \${#FILES[@]} in
        1)
        trim_galore ${extra_args} --basename ${run}_tg \${FILES[0]}
        ;;

        2)
        trim_galore ${extra_args} --basename ${run}_tg --paired \${FILES[0]} \${FILES[1]}
        ;;

        *)
        echo "Error getting files for sample ${run}, exiting"
        exit 1
    esac
    """
}

process createBam {
    label 'bigmem'
    label 'parallel'

    publishDir("${params.report_dir}/logs/bowtie2/", mode: 'copy', pattern: "${run}.log")
    publishDir("${params.report_dir}/logs/dedup/", mode: 'copy', pattern: "${run}.dedup.log", enable: params.dedup_umi)
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(run), path(trimmed)
    file index_files

    output:
    tuple val(run), path("${run}_dedup.bam"), emit: bamfile
    tuple val(run), path('metrics/*'), emit: report
    file "${run}.log"

    script:
    // Check that we have one or two fastq files provided. Set a flag for whether we consequently have paired-end reads
    paired = (trimmed instanceof List)
    if (paired && (trimmed.size() > 2)) throw new Exception("Too many fastq files provided for mapping: ${trimmed.getClass()}")

    // If we have multiple cores, utilise bowtie2s threaded mode
    extra_args = ''
    if (task.cpus > 1) {
        extra_args += "-p ${task.cpus}"
    }
    
    script = "export BOWTIE2_INDEXES=${index_files}\n"

    if (paired) {
        script += "bowtie2 -x ${params.genome} ${extra_args} -1 ${trimmed[0]} -2 ${trimmed[1]} -S ${run}.sam 2> >(tee ${run}.log >&2)\n"
    } else {
        script += "bowtie2 -x ${params.genome} ${extra_args} -U ${trimmed} -S ${run}.sam 2> >(tee ${run}.log >&2)\n"
    }
    
    script += """\
              samtools sort ${run}.sam -o ${run}.sam.sorted
              samtools view -h -S -b ${run}.sam.sorted > ${run}.bam
              mkdir metrics
              """.stripIndent()

    if (params.dedup_umi) {
        // UMI tools needs the bam file to be indexed. We still need to index the deduplicated files afterwards.
        script += "samtools index ${run}.bam\n"
        script += "umi_tools dedup ${params.umi_tools_options} ${paired ? '--paired' : ''} --stdin=${run}.bam --output-stats=\"metrics/${run}\"  --log=${run}.dedup.log > ${run}_dedup.bam\n"
    } else {
        script += "picard MarkDuplicates I=\"${run}.bam\" O=\"${run}_dedup.bam\" M=\"metrics/${run}.metrics\"\n"
    }

    script
}

process index {
    input:
    tuple val(run), path(bamfile)

    output:
    tuple val(run), path(bamfile), path("${bamfile}.bai")

    """
    samtools index ${bamfile}
    """
}

workflow create_bam {
    take:
    fastq_files

    main:
    index_ch = file(params.bowtie2_index, type: 'dir', checkIfExists: true)
    createBam(fastq_files, index_ch)
    createBam.out.bamfile | index

    emit:
    bamfile = index.out
    report = createBam.out.report
}
