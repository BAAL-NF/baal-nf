// Run BaalChIP ASB detection

// Create the baalChIP sample file
// This is objectively a horrible way to write to a file, and feels silly.
// I haven't found a better way to do it because nextflow doesn't give access to the working directory
// when you run an exec command
process createSampleFile {
    publishDir("${params.report_dir}/samples/", mode:'copy')

    input:
    tuple(val(group_name), file(bed_file), val(runs), val(antigens), val(experiments),
          file(snp_files), file(bamfiles), file(index_files))

    output:
    tuple val(group_name), file(bed_file), file(snp_files), file(bamfiles), file(index_files), file("${group_name}.tsv")

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
    tuple val(group_name), file(bed_file), file(snp_file), file(bamfiles), file(index_files), file(sample_file)

    output:
    tuple val(group_name), file('process_bams.rds'), file(snp_file), file(bed_file)

    script:
    """
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
    res <- filter1allele(res)
    saveRDS(res, file='process_bams.rds')
    """
}

process baalGetASB {
    publishDir("${params.report_dir}/baal_reports", mode: 'copy', pattern: "${group_name}.html")

    label 'baal_chip'
    label 'parallel'
    label 'bigmem'

    input:
    tuple val(group_name), file('process_bams.rds'), file(snp_file), file(bed_file)
    path report_md, stageAs: 'baal_report.Rmd'

    output:
    tuple val(group_name), file('*.csv'), file(bed_file), emit: asb
    file("${group_name}.html")

    script:
    """
    #!/usr/bin/env Rscript
    library(BaalChIP)
    library(knitr)
    library(rmarkdown)
    # Read in hets from file
    res <- readRDS("process_bams.rds")
    res <- getASB(res, Iter=5000, conf_level=0.95, cores=${task.cpus})
    saveRDS(res, "final.rds")
    report <- BaalChIP.report(res)

    # Write ASB results to CSV file
    for (group in names(report)) {
            write.csv(report[[group]], paste(group,".csv", sep=""))
    }

    # generate final report
    knit("baal_report.Rmd")
    render("baal_report.md", output_format="all", output_file="${group_name}")
    """
}

workflow run_baal {
    take:
    baal_groups

    main:
    baal_groups | createSampleFile | baalProcessBams
    baalGetASB(baalProcessBams.out, file("${projectDir}/doc/baal_report.Rmd"))

    emit:
    baalGetASB.out.asb
}
