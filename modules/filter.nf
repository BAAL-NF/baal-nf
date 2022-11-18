// Filter bQTLs by concordance with motif score

process pullMotifs {
    label 'nopeak_utils'

    input:
    val antigen
    path pullMotifsScript

    output:
    tuple val(antigen), path("*.jaspar"), optional: true

    script:
    """
    python ${pullMotifsScript} ${antigen}
    """
}

process getGenomepy {
    label 'nopeak_utils'
    storeDir params.genomepy_cache

    output:
    path "${params.assembly}"

    script:
    """
    genomepy install --annotation ${params.assembly} --genomes_dir .
    """
}

process filterByMotifs {
    label 'nopeak_utils'

    input:
    tuple val(antigen), path(motifs), path(snps)

    output:
    path "*.csv"

    script:
    """
    echo "run filter"
    """
}


workflow filter_snps {
    take:
    asbs
    motifs
    groups

    main:
    
    motifs
        .groupTuple()
        .map { antigen, bam_files, motifs, kmers -> antigen }
        .set { tf }
    jaspar = pullMotifs(tf, file("${projectDir}/py/pull_jaspar_motifs.py"))

    genomepy_idx = getGenomepy()

    groups
        .map { group_name, runs, antigens, snp_files, bam_files, index_files -> [group_name, *antigens.unique()] }
        .join(asbs)
        .groupTuple(by: 1)
        .set { group_metadata }

}
