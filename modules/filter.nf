// Filter bQTLs by concordance with motif score

process pullMotifs {
    label 'nopeak_utils'
    storeDir params.jaspar_cache

    input:
    tuple val(antigen), path(snps)
    path pullMotifsScript

    output:
    tuple val(antigen), path(snps), path("${antigen}/*.jaspar"), optional: true

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
    publishDir("${params.report_dir}/asb_filt/${antigen}/", mode: 'copy')

    input:
    tuple val(antigen), path(snps), path(motifs)
    path index
    path filterScript

    output:
    path "*.csv"

    script:
    """
    python ${filterScript} --assembly ${params.assembly} --genomepy_dir ${index} --tf_asb ${antigen} --tf_motifs ${antigen}
    """
}


workflow filter_snps {
    take:
    asbs
    groups

    main:

    groups
        .map { group_name, runs, antigens, snp_files, bam_files, index_files -> [group_name, *antigens.unique()] }
        .join(asbs)
        .groupTuple(by: 1)
        .map { groups, antigen, asbs -> [antigen, asbs] }
        .set { tf_asbs } 
    asbs_motifs = pullMotifs(tf_asbs, file("${projectDir}/py/pull_jaspar_motifs.py"))
    genomepy_idx = getGenomepy()

    filt = filterByMotifs(asbs_motifs, genomepy_idx, file("${projectDir}/py/filter_snps_by_motifs.py"))
}
