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

process filterByJASPARMotifs {
    label 'nopeak_utils'
    label 'moremem'
    publishDir("${params.report_dir}/asb_filt/JASPAR/${antigen}/", mode: 'copy')

    input:
    tuple val(antigen), path(snps), path(motifs)
    path index
    path filterScript

    output:
    path "JASPAR_Motif_metadata_${antigen}.csv", emit: motif_dat
    tuple val(antigen), path("${antigen}_ASBs_JASPAR_Motifs_AR_score_diff.csv"), emit: filt_asb

    script:
    """
    python ${filterScript} --assembly ${params.assembly} --genomepy_dir ${index} --tf_asb ${antigen} --tf_motifs ${antigen}
    """
}

process filterByNoPeakMotifs {
    label 'nopeak_utils'
    label 'moremem'
    publishDir("${params.report_dir}/asb_filt/NoPeak/${antigen}/", mode: 'copy')

    input:
    tuple val(antigen), path(motifs), path(asbs)
    path index
    path filterScript

    output:
    path "NoPeak_Motif_metadata_${antigen}.csv", emit: motif_dat, optional: true
    tuple val(antigen), path("${antigen}_ASBs_NoPeak_Motifs_AR_score_diff.csv"), emit: filt_asb, optional: true

    script:
    """
    python ${filterScript} --assembly ${params.assembly} --k ${params.k} --genomepy_dir ${index} --tf_asb ${antigen} --tf_motifs ${antigen} 
    """
}

process compileMotifInformation {
    label 'nopeak_utils'
    publishDir("${params.report_dir}/asb_filt/compiled/", mode: 'copy')
    
    input:
    tuple val(antigen), path(asbs), val(jaspar), val(nopeak)
    path compileScript
    
    output:
    path "${antigen}.withPeaks.withMotifs.csv"
    
    script:
    """
    # because missing values are set to null when joining these channels, we need to pass these in as a value channel and then symlink the files in order for the process to run properly
	if [[ -e ${jaspar} && -e ${nopeak} ]]; then
		ln -s ${jaspar} .
		ln -s ${nopeak} .
		python ${compileScript} --tf ${antigen} --mode allMotifs
		echo "Running script with all motifs from JASPAR and NoPeak for ${antigen}"
	elif [[ -e ${jaspar} ]]; then 
		echo "Running script with all motifs from JASPAR, as there are  no NoPeak motifs for ${antigen}"
		ln -s ${jaspar} .
		python ${compileScript} --tf ${antigen} --mode noNoPeak
	elif [[ -e ${nopeak} ]]; then 
		echo "Running script with all motifs from NoPeak, as there are no JASPAR motifs for ${antigen}"
		ln -s ${nopeak} .
		python ${compileScript} --tf ${antigen} --mode noJaspar
	else
		echo "Running script with no motifs from NoPeak or JASPAR for ${antigen}"
		python ${compileScript} --tf ${antigen} --mode noMotifs
	fi
    """
}


workflow filter_snps {
    take:
    asbs
    groups
    motifs

    main:

    groups
        .map { group_name, runs, antigens, snp_files, bam_files, index_files -> [group_name, *antigens.unique()] }
        .join(asbs)
        .groupTuple(by: 1)
        .map { groups, antigen, asbs -> [antigen, asbs] }
        .set { tf_asbs } 
    asbs_motifs = pullMotifs(tf_asbs, file("${projectDir}/py/pull_jaspar_motifs.py"))
    genomepy_idx = getGenomepy()

    // Filter by JASPAR motifs
    jaspar_filt = filterByJASPARMotifs(asbs_motifs, genomepy_idx, file("${projectDir}/py/filter_snps_by_JASPAR_motifs.py"))
    
    // Filter by NoPeak motifs
    motifs
        .groupTuple()
        .map { antigen, bam_files, motifs, kmers -> [antigen, motifs] }
        .join(tf_asbs)
        .set { nopeak_motifs }

    nopeak_filt = filterByNoPeakMotifs(nopeak_motifs, genomepy_idx, file("${projectDir}/py/filter_snps_by_NoPeak_motifs.py"))

    // join all annotated asb channels
    tf_asbs
        .join(jaspar_filt.filt_asb, remainder: true)
        .join(nopeak_filt.filt_asb, remainder: true)
        .set { all_asbs }
     
    compileMotifInformation(all_asbs, file("${projectDir}/py/compile_asb_information.py"))

}
