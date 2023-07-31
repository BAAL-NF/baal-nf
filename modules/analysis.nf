// Analysis performed after BaalChIP has been run
process overlapPeaks {
    publishDir(params.baal_output_dir, mode: 'copy')

    input:
    tuple val(group_name), path(asb_file), path(bed_file)
    path script

    output:
    tuple val(group_name), path("${group_name}.withPeaks.csv")

    script:
    """
    python ${script} ${asb_file} ${bed_file} ${group_name}.withPeaks.csv
    """

    stub:
    """
    touch ${group_name}.withPeaks.csv
    """
}

process makeGatBedFiles {
    // the script will exit with status EX_DATAERR if
    // the output from BaalChIP contains no significant SNPs
    errorStrategy { task.exitStatus == 65 ? 'ignore' : 'terminate' }

    input:
    tuple val(group), path(asb_file), path(snp_file)
    file bedfile_script
    
    output:
    tuple val(group), path("foreground.bed"), path("background.bed")

    script:
    """
    python ${bedfile_script} ${asb_file} ${snp_file}
    """

    stub:
    """
    touch foreground.bed
    touch background.bed
    """
}

process runGat {
    publishDir("${params.gat_output_dir}/", mode: "copy")

    input:
    tuple val(group), path(foreground), path(background)
    file annotations

    output:
    tuple val(group), file("${group}.tsv")

    script:
    """
    gat-run.py --segments ${foreground} --annotations=${annotations} --workspace=${background}\
               --ignore-segment-tracks --num-samples=${params.gat_num_samples} --log=gat.log > ${group}.tsv
    """

    stub:
    """
    touch ${group}.tsv
    """
}

workflow process_results {
    take:
    baal_results
    snp_files

    main:
    overlapPeaks(baal_results, file("${projectDir}/py/overlap_beds.py"))
   
    if (params.run_gat) {
        makeGatBedFiles(overlapPeaks.out.join(snp_files),
                        file("${projectDir}/py/make_gat_bedfiles.py"))
        gat = runGat(makeGatBedFiles.out, file("${params.annotation_file}", checkIfExists: true))
    } else {
        gat = Channel.value()
    }

    emit:
    overlap_peaks = overlapPeaks.out
    gat = gat
}

workflow create_report {
    take:
    asb_report
    multiqc_results
    overlap_peaks_results
    gat_results

    main:
    multiqc_flat = multiqc_results.map { key, values -> [key] + values.collect({report -> "${params.multiqc_report_dir}/${key}/${report.name}"}) }
    asb_output = asb_report.map{ key, report -> [key, "${params.baal_report_dir}/${report.name}"]}
    overlap_peaks_results = overlap_peaks_results.map{ key, report -> [key, "${params.baal_output_dir}/${report.name}"] }

    header = ["key", "baal_report", "asb", "multiqc_data", "multiqc_report", "sample_file"]
    combined_results = asb_output.join(overlap_peaks_results).join(multiqc_flat)
    
    if (params.run_gat) {
        gat_results = gat_results.map({key, report -> [key, "${params.gat_output_dir}/${report.name}"]})
        combined_results = combined_results.join(gat_results)
        header += ["gat_results"]
    }

    header_str = header.join(",")
    combined_results.collectFile({ line -> ["pipeline_report.csv", line.join(",")] }, 
                                    storeDir: params.pipeline_report_dir,
                                    newLine: true,
                                    seed: header_str)
}
