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

workflow {
    include "./modules/fastq.nf" params(genome: params.genome, report_dir: params.report_dir, picard_cmd: params.picard_cmd)
    include "./modules/baal.nf" params(report_dir: params.report_dir, mpiflags: params.mpiflags)
    include multi_qc from "./modules/qc.nf"  params(report_dir: params.report_dir)

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

    filter_fastq_before(srr_ch.fastq)

    trimGalore(filter_fastq_before.out.result)
    filter_fastq_after(trimGalore.out.trimmed_fastq)

    // Once filtering is done, we should be able to count the  number of fastq
    // files that will actually go into our analysis
    srr_ch.metadata.join(filter_fastq_after.out.result)
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


    // Generate multiQC reports per TF/Cell line group.
    reports = filter_fastq_after.out.report
                .mix(trimGalore.out.report)
                .mix(create_bam.out.report)
                .groupTuple()
                .map { key, files -> [key, files.flatten() ] }

    multi_qc(srr_ch.metadata, reports)


    bam_files = filtered_data.metadata
                .join(create_bam.out.bamfile)
                .groupTuple(by: 1)
    bam_files | run_baal
}
