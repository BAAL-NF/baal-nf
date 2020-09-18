// Analysis performed after BaalChIP has been run
process overlapPeaks {
    label 'python'
    publishDir("${params.report_dir}/asb", mode: 'copy')

    input:
    tuple val(group_name), file(asb_file), file(bed_file)
    path script

    output:
    tuple val(group_name), file("${group_name}.withPeaks.csv")

    script:
    """
    python ${script} ${asb_file} ${bed_file} ${group_name}.withPeaks.csv
    """
}

process makeGatBedFiles {
    label 'python'

    input:
    tuple val(group), file(asb_file), file(snp_file)
    file bedfile_script
    
    output:
    tuple val(group), file("foreground.bed"), file("background.bed")

    script:
    """
    python ${bedfile_script} ${asb_file} ${snp_file}
    """
}

process runGat {
    label 'python'
    publishDir("${params.report_dir}/enrichment/", mode: "copy")
    
    input:
    tuple val(group), file(foreground), file(background)
    file annotations

    output:
    file "${group}.tsv"

    script:
    """
    gat-run.py --segments ${foreground} --annotations=${annotations} --workspace=${background}\
               --ignore-segment-tracks --num-samples=${params.gat_num_samples} --log=gat.log > ${group}.tsv
    """
}

workflow process_results {
    take:
    baal_results
    snp_files

    main:
    overlapPeaks(baal_results, file("${projectDir}/py/overlap_beds.py"))

    if (params.run_gat){  
        makeGatBedFiles(overlapPeaks.out.join(snp_files), 
                        file("${projectDir}/py/make_gat_bedfiles.py"))
        runGat(makeGatBedFiles.out, file("${params.annotation_file}"))
    }
}