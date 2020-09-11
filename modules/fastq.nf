process trimGalore {
    label 'fastq'
    label 'parallel'

    input:
    tuple val(run), file("${run}*.fastq.gz")

    output:
    tuple val(run), file("${run}_tg*"), emit: trimmed_fastq
    tuple val(run), file('*report*'), emit: report

    script:
    extra_args = ""
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
    label 'fastq'
    label 'bigmem'
    label 'parallel'

    publishDir("${params.report_dir}/logs/bowtie2/", mode: "copy", pattern: "${run}.log")    
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(run), file(trimmed)
    file index_files

    output:
    tuple val(run), file("${run}_dedup.bam"), emit: bamfile
    tuple val(run), file("**.metrics"), emit: report
    file "${run}.log"

    script:
    extra_args = ""
    if (task.cpus > 1) {
        extra_args += "-p ${task.cpus}"
    }
    // I hate the fact that I have to fall back on bash to check the length of this list.
    """
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
    picard MarkDuplicates I="${run}.bam" O="${run}_dedup.bam" M="${run}.metrics"
    """

}

process mergeBeds {
     publishDir("${params.report_dir}/bed_files/", mode: "copy")
     label "fastq"
 
     input:
     tuple val(runs), val(group_name), val(antigens), val(experiments), file(bedfiles), file(snp_files), file(bamfiles), file(index_files)

     output:
     tuple val(runs), val(group_name), val(antigens), val(experiments), file("${group_name}.bed"), file(snp_files), file(bamfiles), file(index_files)
 
     script:
     cat_cmd = ("${bedfiles}".endsWith(".gz")) ? "zcat" : "cat"
     """
     ${cat_cmd} ${bedfiles} | sort -k 1,1 -k2,2n | mergeBed > ${group_name}.bed
     """
}

process index {
    label 'fastq'
    input:
    tuple val(run), file(bamfile)

    output:
    tuple val(run), file(bamfile), file("${bamfile}.bai")

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
