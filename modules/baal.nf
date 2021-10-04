// Run BaalChIP ASB detection

// Create the baalChIP sample file
// This is objectively a horrible way to write to a file, and feels silly.
// I haven't found a better way to do it because nextflow doesn't give access to the working directory
// when you run an exec command
process createSampleFile {
    publishDir("${params.report_dir}/samples/", mode:'copy')

    input:
    tuple(val(group_name), path(bed_file), val(runs), val(antigens),
          path(snp_paths), path(bamfiles), path(index_files))

    output:
    tuple val(group_name), path("${group_name}.tsv")

    script:
        output = "cat << EOF > ${group_name}.tsv\n"
        output += "group_name\ttarget\treplicate_number\tbam_name\tbed_name\tsampleID\n"
        0.upto(runs.size() - 1) {
            i ->
            replicate = i + 1
            output +=
                "${group_name}\t${antigens[i]}\t${replicate}\t${bamfiles[i].name}\t${bed_file.name}\t${runs[i]}\n"
        }
        output += 'EOF\n'
        output
}

process baalProcessBams {
    label 'baal_chip'

    input:
    tuple val(group_name), path(bed_file), path(snp_file), path(bamfiles), path(index_files), path(sample_file)

    output:
    tuple val(group_name), path('process_bams.rds'), path(snp_file), path(bed_file)

    script:
    script = """
    #!/usr/bin/env Rscript
    library(BaalChIP)

    data(blacklist_hg19)
    data(pickrell2011cov1_hg19)
    data(UniqueMappability50bp_hg19)

    samplesheet <- '${sample_file}'

    hets <- c('${group_name}' = '${snp_file}')

    res <- new('BaalChIP', samplesheet=samplesheet, hets=hets)
    res <- alleleCounts(res, min_base_quality=10, min_mapq=15, all_hets=TRUE)
    res <- QCfilter(res,
                    RegionsToFilter=list('blacklist'=blacklist_hg19,
                                         'highcoverage'=pickrell2011cov1_hg19),
                    RegionsToKeep=list('UniqueMappability'=UniqueMappability50bp_hg19))
    res <- mergePerGroup(res)
    """
    
    if( ! params.dedup_umi ) script += "res <- filter1allele(res)\n"

    script += "saveRDS(res, file='process_bams.rds')\n"
    script    
}

process baalGetASB {
    publishDir(params.baal_report_dir, mode: 'copy', pattern: "${group_name}.html")
    // Optionally export the object files
    if (params.save_baal_objects) publishDir(params.baal_object_dir, mode: 'copy', pattern: "*.rds")

    label 'baal_chip'
    label 'parallel'

    input:
    tuple val(group_name), path('process_bams.rds'), path(snp_file), path(bed_file)
    path report_md, stageAs: 'baal_report.Rmd'

    output:
    tuple val(group_name), path("${group_name}.csv"), path(bed_file), emit: asb
    tuple val(group_name), path("${group_name}.html"), emit: report
    path("${group_name}.rds")

    script:
    """
    #!/usr/bin/env Rscript
    library(BaalChIP)
    library(knitr)
    library(rmarkdown)
    # Read in hets from file
    res <- readRDS("process_bams.rds")
    res <- getASB(res, Iter=5000, conf_level=c(${params.confidence_levels.join(',')}), cores=${task.cpus}, clusterType = "PSOCK")
    saveRDS(res, "${group_name}.rds")
    report <- BaalChIP.report(res)

    # Write ASB results to CSV file
    write.csv(report[["${group_name}"]], "${group_name}.csv")

    # generate final report
    knit("baal_report.Rmd")
    render("baal_report.md", output_format="all", output_file="${group_name}")
    """
}

workflow run_baal {
    take:
    baal_groups

    main:
    samples = createSampleFile(baal_groups)
    baal_groups
        .map { 
            group_name, bed_file, runs, antigens, snp_files, bamfiles, index_files -> 
            [ group_name, bed_file, snp_files, bamfiles, index_files ] }
        .join(samples) | baalProcessBams
    baalGetASB(baalProcessBams.out, file("${projectDir}/doc/baal_report.Rmd"))

    emit:
    asb = baalGetASB.out.asb
    report = baalGetASB.out.report
}
