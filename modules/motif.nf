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
    label 'bigmem'
    label 'nopeak'

    publishDir("${params.report_dir}/motifs/${antigen}/profiles/${kmer}mer/", mode: "copy")

    input:
    tuple val(antigen), path(bam_file), path(bed_file)
    path genome
    each kmer

    output:
    tuple val(antigen), val("${bam_file}"), val(kmer), path("profile_${bed_file}.csv")

    script:
    """
    noPeak PROFILE -t ${task.cpus} --reads ${bed_file} --genome ${genome} -k ${kmer}
    """
}

process getFragmentSize {
    publishDir("${params.report_dir}/motifs/${antigen}/spp/", mode: "copy")
    input:
    tuple val(antigen), path(bam_file)

    output:
    tuple val(antigen), path(bam_file), path("${bam_file.baseName}.txt"), env(quality)

    script:
    """
    run_spp.R -c=${bam_file} -out=${bam_file.baseName}.txt

    # Keep only top fragment size in this file
    sed -r 's/,[^\t]+//g' ${bam_file.baseName}.txt > ${bam_file.baseName}.filt.txt

    # Export quality score as variable in the output channel
    quality=\$( awk '{print \$11}' ${bam_file.baseName}.filt.txt )
    """
}

process parseFragmentSize {
    input:
    tuple val(antigen), path(bam_file), path(spp_output), val(quality)
    path parse_script

    output:
    tuple val(antigen), val("${bam_file}"), path("${bam_file.baseName}.fragment_size.txt"), val(quality)

    script:
    """
    python ${parse_script} ${spp_output} ${bam_file.baseName}.fragment_size.txt
    """
}

process getMotifs {
    label 'nopeak'

    publishDir("${params.report_dir}/motifs/${antigen}/motifs/${kmer}mer/", mode: "copy")

    input:
    tuple  val(antigen), val(bam_file), val(kmer), path(profile), path(fragment_size), val(quality)

    output:
    tuple val(antigen), val(bam_file), path("${bam_file}.motifs.txt"), path("${bam_file}.kmers.txt")

    when:
    quality != "-2" && quality != "-1"

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

    // Create value channel for sensitivity analysis
    kmer_lengths = Channel.of(6,8)

    parse_script = file("${projectDir}/py/get_fragment_sizes.py")
    getFragmentSize(input) | set { spp_output }
    parseFragmentSize(spp_output, parse_script)
    
    genome = file(params.nopeak_index, type: 'dir', checkIfExists: true)
    bamToBed(input) | set { pileup }
    profileMotifs(pileup, genome, kmer_lengths)

    profileMotifs.out.join(parseFragmentSize.out, by:[0,1] ) | getMotifs

    emit:
    motifs = getMotifs.out
}
