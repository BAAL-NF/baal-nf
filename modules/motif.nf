process bamToBed {
    input:
    file bam_file

    output:
    tuple val("${bam_file}"), file("${bam_file.baseName}.bed")

    script:
    """
    bamToBed -i ${bam_file} | sort -k1,1 -k2n  > ${bam_file.baseName}.bed
    """
}

process mergeBeds {
     publishDir("${params.report_dir}/bed_files/", mode: 'copy')

     input:
     tuple val(group_name), path(bedfiles)

     output:
     tuple val(group_name), path("${group_name}.bed")

     script:
     // Use zcat -f since it works with both gzipped and plaintext files
     """
     zcat -f ${bedfiles} | sort -k 1,1 -k2,2n | mergeBed > ${group_name}.bed
     """
}

process profileMotifs {
    container 'oalmelid/nopeak:latest'
    publishDir("${params.report_dir}/motifs/profiles/", mode: "copy")

    input:
    tuple val(bam_file), path(bed_file)
    path genome

    output:
    tuple val(bam_file), path("profile_${bed_file}.csv")

    script:
    // FIXME: add a nopeak script to container for easy running
    // FIXME: copy utility scripts to container
    """
    java -jar /usr/local/lib/noPeak.jar PROFILE --reads ${bed_file} --genome ${genome} -k ${params.motif_kmer_length}
    """
}

process getFragmentSize {
    publishDir("${params.report_dir}/spp/", mode: "copy")
    input:
    path(bam_file)

    output:
    tuple val("${bam_file}"), path("${bam_file.baseName}.txt")

    script:
    """
    run_spp.R -c=${bam_file} -out=${bam_file.baseName}.txt
    """
}

process getMotifs {
    container 'oalmelid/nopeak:latest'
    publishDir("${params.report_dir}/motifs/motifs/", mode: "copy")

    input:
    tuple val(bam_file), path(profile), path(fragment_size)

    output:
    tuple val(bam_file), path("${bam_file}.motifs.txt"), path("${bam_file}.kmers.txt")

    script:
    """
    FRAG_SIZE=`cat ${fragment_size} | cut -f3 | cut -d, -f1`
    java -jar /usr/local/lib/noPeak.jar LOGO --strict --signal ${profile} --fraglen \$FRAG_SIZE --export-kmers ${bam_file}.kmers.txt > ${bam_file}.motifs.txt
    """

}

workflow no_peak {
    take:
    input
    main:
    
    input | map { antigen, bam_file -> bam_file } | set { bam_files }
    //| groupTuple() | mergeBeds
    getFragmentSize(bam_files) | set { fragment_sizes }
    
    genome = file(params.nopeak_index, type: 'dir', checkIfExists: true)
    bamToBed(bam_files) | set { pileup }
    profileMotifs(pileup, genome)

    profileMotifs.out.join(fragment_sizes) | getMotifs
}