params.report_dir = "${workflow.launchDir}/test/reports/"
params.mpiflags = ""

// This is objectively a horrible way to write to a file, and feels silly.
// I haven't found a better way to do it because nextflow doesn't give access to the working directory
// when you run an exec command
process createSampleFile {
    publishDir("${params.report_dir}/samples/", mode:"copy")

    input:
    tuple runs, group_name, antigens, experiments, bed_file, snp_files, bamfiles, index_files

    output:
    tuple group_name, bed_file, snp_files, bamfiles, index_files, file("${group_name}.tsv")

    script:
        output = "cat << EOF > ${group_name}.tsv\n"
        output += "group_name\ttarget\treplicate_number\tbam_name\tbed_name\tsampleID\n"
        0.upto(runs.size()-1, {
            replicate = it + 1
            output += "${group_name}\t${antigens[it]}\t${replicate}\t${bamfiles[it].name}\t${bed_file.name}\t${runs[it]}\n"
        })
        output += "EOF\n"
        output
}

process baalProcessBams {
    label "baal_chip"
    label "retry_mem"

    errorStrategy { (task.attempt < 3) ? "retry" : "ignore"}

    input:
    tuple group_name, file(bed_file), file(snp_file), file(bamfiles), file(index_files), file(sample_file)

    output:
    tuple group_name, file("process_bams.rds"), file(snp_file), file(bed_file)

    script:
    """
    #!/usr/bin/env Rscript
    library(BaalChIP)

    data(blacklist_hg19)
    data(pickrell2011cov1_hg19)
    data(UniqueMappability50bp_hg19)

    samplesheet <- "${sample_file}"

    hets <- c("${group_name}" = "${snp_file}")

    res <- new("BaalChIP", samplesheet=samplesheet, hets=hets)
    res <- alleleCounts(res, min_base_quality=10, min_mapq=15, all_hets=TRUE)
    res <- QCfilter(res,
                    RegionsToFilter=list("blacklist"=blacklist_hg19,
                                         "highcoverage"=pickrell2011cov1_hg19),
                    RegionsToKeep=list("UniqueMappability"=UniqueMappability50bp_hg19))
    res <- mergePerGroup(res)
    res <- filter1allele(res)
    saveRDS(res, file="process_bams.rds")
    """
}

process baalGetASB {
    publishDir("${params.report_dir}/baal_reports", mode: "copy", pattern: "${group_name}.html")
    errorStrategy { (task.attempt < 5) ? "retry" : "ignore"}

    label "baal_chip"
    label "parallel"
    label "bigmem"

    input:
    tuple group_name, file("process_bams.rds"), file(snp_file), file(bed_file)

    output:
    tuple group_name, file("*.csv"), file(bed_file), emit: asb
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

    # set the knitr working directory to the current working directory.
    opts_knit$set(root.dir = getwd())
    # generate final report
    knit("${workflow.projectDir}/doc/baal_report.Rmd")
    render("baal_report.md", output_format="all", output_file="${group_name}")
    """
}

process overlapPeaks {
    label "python"
    publishDir("${params.report_dir}/asb", mode: "copy")

    input:
    tuple group_name, file(asb_file), file(bed_file)

    output:
    file("${group_name}.withPeaks.csv")

    script:
    """
    python ${workflow.projectDir}/py/overlap_beds.py ${asb_file} ${bed_file} ${group_name}.withPeaks.csv 
    """
}

workflow run_baal {
    take:
    baal_groups

    main:
    baal_groups | createSampleFile | baalProcessBams | baalGetASB 
    baalGetASB.out.asb | overlapPeaks

    emit:
    overlapPeaks.out
}
