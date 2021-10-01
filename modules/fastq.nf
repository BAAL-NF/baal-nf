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
    if (params.dedup_umi) publishDir("${params.report_dir}/logs/dedup/", mode: 'copy', pattern: "${run}.dedup.log")
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(run), path(trimmed)
    file index_files

    output:
    tuple val(run), path("${run}_dedup.bam"), emit: bamfile
    tuple val(run), path('**.metrics'), emit: report
    file "${run}.log"

    script:
    extra_args = ''
    if (task.cpus > 1) {
        extra_args += "-p ${task.cpus}"
    }
    // I hate the fact that I have to fall back on bash to check the length of this list.
    result = """
    export BOWTIE2_INDEXES=${index_files}
    FILES=(${trimmed})
    case \${#FILES[@]} in
        1)
        bowtie2 -x ${params.genome} ${extra_args} -U \${FILES[0]} -S ${run}.sam 2> >(tee ${run}.log >&2)\n
        ;;

        2)
        bowtie2 -x ${params.genome} ${extra_args} -1 \${FILES[0]} -2 \${FILES[1]} -S ${run}.sam 2> >(tee ${run}.log >&2)\n
        ;;

        *)
        echo Error getting files for sample ${run}, exiting
        exit 1
    esac
    samtools sort ${run}.sam -o ${run}.sam.sorted
    samtools view -h -S -b ${run}.sam.sorted > ${run}.bam
    """

    if (params.dedup_umi) {
        result += "umi_tools dedup ${params.umi_tools_options} --stdin=${run.bam} --log=${run}.dedup.log > ${run}_dedup.bam"
    } else {
        result += "picard MarkDuplicates I=\"${run}.bam\" O=\"${run}_dedup.bam\" M=\"${run}.metrics\""
    }
    result
}

process mergeBeds {
     publishDir("${params.report_dir}/bed_files/", mode: 'copy')

     input:
     tuple(val(group_name), path(bedfiles))

     output:
     tuple(val(group_name), path("${group_name}.bed"))

     script:
     // Use zcat -f since it works with both gzipped and plaintext files
     """
     zcat -f ${bedfiles} | sort -k 1,1 -k2,2n | mergeBed > ${group_name}.bed
     """
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
