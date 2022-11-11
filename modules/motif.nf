process bamToBed {
    input:
    tuple val(antigen), path(bam_file)

    output:
    tuple val(antigen), path(bam_file), path("${bam_file.baseName}.bed")

    script:
    """
    bamToBed -i ${bam_file} | sort -k1,1 -k2n  > ${bam_file.baseName}.bed
    """
}

process profileMotifs {
    label 'parallel'
    label 'nopeak'

    publishDir("${params.report_dir}/motifs/${antigen}/profiles/", mode: "copy")

    input:
    tuple val(antigen), path(bam_file), path(bed_file)
    path genome

    output:
    tuple val(antigen), val("${bam_file}"),  path("profile_${bed_file}.csv")

    script:
    """
    noPeak PROFILE -t ${task.cpus} --reads ${bed_file} --genome ${genome} -k ${params.motif_kmer_length}
    """
}

process getFragmentSize {
    publishDir("${params.report_dir}/motifs/${antigen}/spp/", mode: "copy")
    input:
    tuple val(antigen), path(bam_file)

    output:
    tuple val(antigen), path(bam_file), path("${bam_file.baseName}.txt")

    script:
    """
    run_spp.R -c=${bam_file} -out=${bam_file.baseName}.txt
    """
}

process parseFragmentSize {
    input:
    tuple val(antigen), path(bam_file), path(spp_output)
    path parse_script

    output:
    tuple val(antigen), val("${bam_file}"), path("${bam_file.baseName}.fragment_size.txt")

    script:
    """
    python ${parse_script} ${spp_output} ${bam_file.baseName}.fragment_size.txt
    """
}

process getMotifs {
    label 'nopeak'

    publishDir("${params.report_dir}/motifs/${antigen}/motifs/", mode: "copy")

    input:
    tuple  val(antigen), val(bam_file), path(profile), path(fragment_size)

    output:
    tuple val(bam_file), path("${bam_file}.motifs.txt"), path("${bam_file}.kmers.txt")

    script:
    """
    FRAGMENT_SIZE=`cat ${fragment_size}`
    noPeak LOGO --strict --signal ${profile} --fraglen \$FRAGMENT_SIZE --export-kmers ${bam_file}.kmers.txt > ${bam_file}.motifs.txt
    """

}

workflow no_peak {
    take:
    input

    main:

    parse_script = file("${projectDir}/py/get_fragment_sizes.py")
    getFragmentSize(input) | set { spp_output }
    parseFragmentSize(spp_output, parse_script)
    
    genome = file(params.nopeak_index, type: 'dir', checkIfExists: true)
    bamToBed(input) | set { pileup }
    profileMotifs(pileup, genome)

    profileMotifs.out.join(parseFragmentSize.out, by:[0,1] ) | getMotifs
}
