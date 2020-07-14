#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.staging_root = "${workflow.projectDir}/test/staging/runs"
params.experiments = "${params.staging_root}/runs.csv"
params.report_dir = "${workflow.projectDir}/test/reports/"
params.genome = "hg19"
params.picard_cmd = "picard"
params.mpiflags = ""
params.fastqc_conf_pre = "${workflow.projectDir}/data/before_limits.txt"
params.fastqc_conf_post = "${workflow.projectDir}/data/after_limits.txt"
params.fastq_screen_conf = ""
params.run_baal = true

// Import a CSV file with all sample identifiers
workflow import_samples {
    main:
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
        .set{ srr_ch }


    emit:
    fastq = srr_ch.fastq
    metadata = srr_ch.metadata
}

// Count fastq files that have made it through filtering, create a group key, and separate out the fastq files
// into their own channel for the sake of further processing.
workflow count_fastq {
    get:
    metadata
    fastq_files

    main:
    metadata
        .join(fastq_files)
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

    emit:
    fastq = filtered_data.fastq
    metadata =  filtered_data.metadata
}

// Filtering stages. These are done before and after TrimGalore.
workflow filter_fastq_before {
    include filter_fastq as pre_filter_fastq from "./modules/qc.nf" params(report_dir: params.report_dir, fastqc_conf: params.fastqc_conf_pre)
    get:
    fastq_list

    main:
    pre_filter_fastq(fastq_list)

    emit:
    result = pre_filter_fastq.out.fastq_list
    report = pre_filter_fastq.out.report
}

// After TrimGalore we also run FastQ-screen to check for contamination
workflow filter_fastq_after {
    include filter_fastq as post_filter_fastq from "./modules/qc.nf" params(report_dir: params.report_dir, fastqc_conf: params.fastqc_conf_post)
    include fastq_screen from "./modules/qc.nf" params(fastq_screen_conf: params.fastq_screen_conf)

    get:
    fastq_list

    main:
    post_filter_fastq(fastq_list)
    fastq_list | fastq_screen

    result = post_filter_fastq.out.fastq_list.join(fastq_screen.out.result)
    report = post_filter_fastq.out.report
                .mix(fastq_screen.out.report)

    emit:
    result = result
    report = report
}

process reportFastQC {
    publishDir("${params.report_dir}/preFastQC/${key}/", mode: "copy")

    input:
    tuple key, file(reports)

    output:
    file reports

    script:
    "exit 0"
}

process mergeBeds {
     publishDir(params.out_dir, mode: "move")
 
     input:
     tuple runs, group_name, antigens, experiments, bedfiles, snp_files, bamfiles, index_files

     output:
     tuple runs, group_name, antigens, experiments, file("${group_name}.bed"), snp_files, bamfiles, index_files
 
     script:
     """
     cat ${bedfiles} | sort -k 1,1 -k2,2n | mergeBed > ${group_name}.bed
     """
}

workflow {
    include "./modules/fastq.nf" params(genome: params.genome, report_dir: params.report_dir, picard_cmd: params.picard_cmd)
    include "./modules/baal.nf" params(report_dir: params.report_dir, mpiflags: params.mpiflags)
    include multi_qc from "./modules/qc.nf"  params(report_dir: params.report_dir)

    // Load CSV file
    import_samples()

    // Pre-trimming filtering step
    filter_fastq_before(import_samples.out.fastq)

    import_samples.out.metadata
        .join(filter_fastq_before.out.report)
        .groupTuple(by:1)
        .map{run, group, transcription_factor, experiment, bed_file, snp_file, reports -> [group, reports.flatten()]} | reportFastQC

    // Adapter trimming
    trimGalore(filter_fastq_before.out.result)

    // Rerun fastQC and fastqScreen
    filter_fastq_after(trimGalore.out.trimmed_fastq)

    // Once filtering is done, we should be able to count the  number of fastq
    // files that will actually go into our analysis
    count_fastq(import_samples.out.metadata, filter_fastq_after.out.result)

    // Create BAM files for each SRR
    count_fastq.out.fastq | create_bam

    // Generate multiQC reports per TF/Cell line group.
    reports = filter_fastq_after.out.report
                .mix(trimGalore.out.report)
                .mix(create_bam.out.report)
                .groupTuple()
                .map { key, files -> [key, files.flatten() ] }

    multi_qc(import_samples.out.metadata, reports)

    // Regroup the bam files with their associated metadata and run baal chip
    bam_files = count_fastq.out.metadata
                .join(create_bam.out.bamfile)
                .groupTuple(by: 1)
		.mergeBeds()

    if (params.run_baal) {
        bam_files | run_baal
    }
}
