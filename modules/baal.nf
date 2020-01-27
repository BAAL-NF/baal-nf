params.report_dir = "${workflow.launchDir}/test/reports/"
params.mpiflags = ""

// This is objectively a horrible way to write to a file, and feels silly.
// I haven't found a better way to do it because nextflow doesn't give access to the working directory
// when you run an exec command
process createSampleFile {
    publishDir("${params.report_dir}/samples/")

    input:
    tuple runs, group_name, antigens, experiments, bedfiles, snp_files, bamfiles, index_files

    output:
    tuple group_name, bedfiles, snp_files, bamfiles, index_files, file("${group_name}.tsv")

    script:
        output = "cat << EOF > ${group_name}.tsv\n"
        output += "group_name\ttarget\treplicate_number\tbam_name\tbed_name\tsampleID\n"
        0.upto(runs.size()-1, {
            replicate = it + 1
            output += "${group_name}\t${antigens[it]}\t${replicate}\t${bamfiles[it].name}\t${bedfiles[it].name}\t${runs[it]}\n"
        })
        output += "EOF\n"
        output
}

process baalProcessBams {
    label "baal_chip"

    input:
    tuple group_name, file(bedfiles), file(snp_file), file(bamfiles), file(index_files), file(sample_file)

    output:
    tuple file("process_bams.rds"), file(snp_file)

    script:
    """
    #!/usr/bin/env Rscript
    library(BaalChIP)

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
    publishDir("${params.report_dir}/asb", mode: "copy")
    errorStrategy "retry"
    maxRetries 3

    label "baal_chip"
    label "mpi"
    label "bigmem"

    input:
    tuple file("process_bams.rds"), file(snp_file)

    output:
    file "*.csv"

    script:
    """
    #!/usr/bin/env mpirun -np 1 ${params.mpiflags} Rscript
    library(BaalChIP)
    # Read in hets from file
    res <- readRDS("process_bams.rds")
    res <- getASB(res, Iter=5000, conf_level=0.95)
    saveRDS(res, "final.rds")
    report <- BaalChIP.report(res)

    for (group in names(report)) {
            write.csv(report[[group]], paste(group,".csv", sep=""))
    }
    """
}

workflow run_baal {
    get:
    baal_groups

    main:
    baal_groups | createSampleFile

    sample_files = createSampleFile.out.map({
        group_name, bedfiles, snp_files, bamfiles, index_files, sample_file ->
        [group_name, bedfiles.unique(), snp_files.unique(), bamfiles, index_files, sample_file]
    })

    sample_files | baalProcessBams | baalGetASB

    emit:
    baalProcessBams.out
}
