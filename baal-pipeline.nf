#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.staging_root = "${workflow.launchDir}/test/staging/runs"
params.experiments = "${params.staging_root}/runs.csv"
params.report_dir = "${workflow.launchDir}/test/reports/"
params.genome = "hg19"
params.picard_cmd = "picard"
params.mpiflags = ""
params.fastqc_conf = ""
params.fastq_screen_conf = ""

// key: "${row.cell_line}_${row.transcription_factor}"

workflow filter_fastq_before {
    include filter_fastq as pre_filter_fastq from "./modules/qc.nf" params(report_dir: params.report_dir, fastqc_conf: params.fastqc_conf)
    get:
    fastq_list

    main:
    result = pre_filter_fastq(fastq_list, "before")

    emit:
    result
}

workflow filter_fastq_after {
    include filter_fastq as post_filter_fastq from "./modules/qc.nf" params(report_dir: params.report_dir, fastqc_conf: params.fastqc_conf)
    include fastq_screen from "./modules/qc.nf" params(fastq_screen_conf: params.fastq_screen_conf)

    get:
    fastq_list

    main:
    post_filter_fastq(fastq_list, "after")
    fastq_list | fastq_screen
    result = post_filter_fastq.out.join(fastq_screen.out)

    emit:
    result
}

workflow {
    include "./modules/fastq.nf" params(genome: params.genome, report_dir: params.report_dir, picard_cmd: params.picard_cmd)
    include "./modules/baal.nf" params(report_dir: params.report_dir)

    Channel
        .fromPath(params.experiments)
        .splitCsv(header:true)
        .map({row -> tuple(row.run,
                           "${row.cell_line}_${row.transcription_factor}",
                           row.transcription_factor,
                           row.experiment,
                           file("${params.staging_root}/${row.fastq_folder}*.fastq.gz"),
                           file("${params.staging_root}/${row.bed_file}"),
                           file("${params.staging_root}/${row.snp_list}"))})
        .fork {
            run, group, transcription_factor, experiment, fastq_files, bed_file, snp_file ->
            fastq: [run, fastq_files]
            metadata: [run, group, transcription_factor, experiment, bed_file, snp_file]
        }
        .set { srr_ch }

    filter_fastq_before(srr_ch.fastq) | trimGalore
    trimGalore.out.trimmed_fastq | filter_fastq_after
    // Once filtering is done, we should be able to count the  number of fastq
    // files that will actually go into our analysis
    srr_ch.metadata.join(filter_fastq_after.out)
        .groupTuple(by: 1)
        .map({ runs, tag, transcription_factor, experiments, bed_files, snp_files, fastq_files ->
            num_samples = runs.size()
            [ runs, groupKey(tag, num_samples), transcription_factor, experiments, bed_files, snp_files, fastq_files ]
        })
        .transpose()
        .fork {
          run, key, transcription_factor, experiment, bed_file, snp_file, fastq_files ->
          fastq: [ run, fastq_files ]
          metadata: [ run, key, transcription_factor, experiment, bed_file, snp_file ]
        }
        .set { filtered_data }

    filtered_data.fastq | create_bam
    bam_files = filtered_data.metadata
                .join(create_bam.out)
                .groupTuple(by: 1)
    bam_files | run_baal
}



// baal_sample_file_ch
//     .map({
//         runs, group_name, cell_lines, antigens, experiments, bedfiles, snp_files, bamfiles, index_files, sample_file ->
//         [runs, group_name, cell_lines, antigens, experiments, bedfiles.unique(), snp_files[0], bamfiles, index_files, sample_file]
//     })
//     .set { baal_process_bams_ch }

// process baalProcessBams {
//     label "baal_chip"

//     input:
//     set runs, group_name, cell_lines, antigens, experiments, file(bedfiles), file(snp_file), file(bamfiles), file(index_files), file(sample_file) from baal_process_bams_ch

//     output:
//     set file("process_bams.rds"), file(snp_file) into get_asb_ch

//     script:
//     """
//     #!/usr/bin/env Rscript
//     library(BaalChIP)

//     samplesheet <- "${sample_file}"

//     hets <- c("${group_name}" = "${snp_file}")

//     res <- new("BaalChIP", samplesheet=samplesheet, hets=hets)
//     res <- alleleCounts(res, min_base_quality=10, min_mapq=15, all_hets=TRUE)
//     res <- QCfilter(res,
//                     RegionsToFilter=list("blacklist"=blacklist_hg19,
//                                          "highcoverage"=pickrell2011cov1_hg19),
//                     RegionsToKeep=list("UniqueMappability"=UniqueMappability50bp_hg19))
//     res <- mergePerGroup(res)
//     res <- filter1allele(res)
//     saveRDS(res, file="process_bams.rds")
//     """
// }

// process baalGetASB {
//     publishDir("${params.report_dir}/asb")

//     label "baal_chip"
//     label "mpi"

//     input:
//     set file("process_bams.rds"), file(snp_file) from get_asb_ch

//     output:
//     file "*.csv"

//     script:
//     """
//     #!/usr/bin/env mpirun -np 1 ${params.mpiflags} Rscript
//     library(BaalChIP)
//     # Read in hets from file
//     res <- readRDS("process_bams.rds")
//     res <- getASB(res, Iter=5000, conf_level=0.95)
//     saveRDS(res, "final.rds")
//     report <- BaalChIP.report(res)

//     for (group in names(report)) {
//             write.csv(report[[group]], paste(group,".csv", sep=""))
//     }
//     """
// }
