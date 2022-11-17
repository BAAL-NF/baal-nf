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

workflow filter_snps {
    take:
    asbs
    motifs

    main:
    
    motifs
        .groupTuple()
        .map { antigen, bam_files, motifs, kmers -> antigen }
        .set { tf }
    jaspar = pullMotifs(tf, file("${projectDir}/py/pull_jaspar_motifs.py"))

    jaspar.view()

}
