// Filter bQTLs by concordance with motif score

process pullMotifs {
    label 'nopeak_utils'
    storeDir params.jaspar_cache

    input:
    val antigen
    path pullMotifsScript

    output:
    tuple val(antigen), path("${antigen}/*.jaspar"), optional: true

    script:
    """
    python ${pullMotifsScript} ${antigen}
    """

    stub:
    """
    mkdir ${antigen}
    touch ${antigen}/motif1.jaspar
    touch ${antigen}/motif2.jaspar
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

    stub:
    """
    mkdir ${params.assembly}
    touch ${params.assembly}/README.txt
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
    path "JASPAR_Motif_metadata_${antigen}.csv", emit: motif_dat, optional: true
    tuple val(antigen), path("${antigen}_ASBs_JASPAR_Motifs_AR_score_diff.csv"), emit: filt_asb, optional: true

    script:
    """
    NEW_CACHE=\$PWD/cache
    mkdir -p \$NEW_CACHE
    export XDG_CACHE_HOME=\$NEW_CACHE
    echo "Using \$XDG_CACHE_HOME for cache"

    python ${filterScript} --assembly ${params.assembly} --genomepy_dir ${index} --tf_asb ${antigen} --tf_motifs ${antigen}
    """

    stub:
    """
    touch JASPAR_Motif_metadata_${antigen}.csv
    touch ${antigen}_ASBs_JASPAR_Motifs_AR_score_diff.csv
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
    NEW_CACHE=\$PWD/cache
    mkdir -p \$NEW_CACHE
    export XDG_CACHE_HOME=\$NEW_CACHE
    echo "Using \$XDG_CACHE_HOME for cache"

    python ${filterScript} --assembly ${params.assembly} --k ${params.k} --genomepy_dir ${index} --tf_asb ${antigen} --tf_motifs ${antigen} 
    """

    stub:
    """
    touch NoPeak_Motif_metadata_${antigen}.csv
    touch ${antigen}_ASBs_NoPeak_Motifs_AR_score_diff.csv
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

    stub:
    """
    touch ${antigen}.withPeaks.withMotifs.csv
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
        .multiMap { 
            groups, antigen, asbs -> 
            tf_asbs: [antigen, asbs]
            tfs: antigen 
        }
        .set { out_ch }
    jaspar_motifs = pullMotifs(out_ch.tfs, file("${projectDir}/py/pull_jaspar_motifs.py"))
    asbs_motifs = out_ch.tf_asbs.join(jaspar_motifs)

    genomepy_idx = getGenomepy()

    // Filter by JASPAR motifs
    jaspar_filt = filterByJASPARMotifs(asbs_motifs, genomepy_idx, file("${projectDir}/py/filter_snps_by_JASPAR_motifs.py"))
  
    // Filter by NoPeak motifs
    motifs
        .groupTuple()
        .map { antigen, bam_files, motifs, kmers -> [antigen, motifs] }
        .join(out_ch.tf_asbs)
        .set { nopeak_motifs }

    nopeak_filt = filterByNoPeakMotifs(nopeak_motifs, genomepy_idx, file("${projectDir}/py/filter_snps_by_NoPeak_motifs.py"))

    // join all annotated asb channels
    out_ch.tf_asbs
        .join(jaspar_filt.filt_asb, remainder: true)
        .join(nopeak_filt.filt_asb, remainder: true)
        .set { all_asbs }
     
    compileMotifInformation(all_asbs, file("${projectDir}/py/compile_asb_information.py"))

}
