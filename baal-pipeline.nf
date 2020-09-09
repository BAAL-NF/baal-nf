#!/usr/bin/env nextflow
nextflow.preview.dsl=2

if (params.sample_file.isEmpty()) {
    println("Error: no sample file provided")
    exit 1
}

// Import a CSV file with all sample identifiers
workflow import_samples {
    main:
    Channel
        .fromPath(params.sample_file, checkIfExists: true)
        .splitCsv(header:true)
        .map({row -> tuple(
            row.run,
            "${row.cell_line}_${row.transcription_factor}",
            row.transcription_factor,
            row.experiment,
            (([row.fastq_1, row.fastq_2] - "").collect { path -> file(path, checkIfExists: true) }),
            file("${row.bed_file}", checkIfExists: true),
            file("${row.snp_list}", checkIfExists: true))}
    ).multiMap {
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
    take:
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
        .multiMap {
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
    include { filter_fastq as fastqc_before_trimming } from "./modules/qc.nf" addParams(fastqc_conf: params.fastqc_conf_pre, report_dir: params.report_dir)
    take:
    fastq_list

    main:
    fastqc_before_trimming(fastq_list)
    result = fastqc_before_trimming.out.fastq_list.join(fastq_list)

    emit:
    result = result
    report = fastqc_before_trimming.out.report
}

// After TrimGalore we also run FastQ-screen to check for contamination
workflow filter_fastq_after {
    include { filter_fastq as fastqc_after_trimming; fastq_screen } from "./modules/qc.nf" addParams(fastqc_conf: params.fastqc_conf_post, report_dir: params.report_dir)

    take:
    fastq_list

    main:
    fastqc_after_trimming(fastq_list)
    fastq_screen(fastq_list)

    result = fastqc_after_trimming.out.fastq_list
                                      .join(fastq_screen.out.result)
                                      .join(fastq_list)
    report = fastqc_after_trimming.out.report
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
     publishDir("${params.report_dir}/bed_files/", mode: "copy")
     label "fastq"
 
     input:
     tuple runs, group_name, antigens, experiments, file(bedfiles), snp_files, bamfiles, index_files

     output:
     tuple runs, group_name, antigens, experiments, file("${group_name}.bed"), snp_files, bamfiles, index_files
 
     script:
     """
     cat ${bedfiles} | sort -k 1,1 -k2,2n | mergeBed > ${group_name}.bed
     """
}

workflow {
    include { trimGalore; create_bam } from "./modules/fastq.nf" params(bowtie2_index: params.bowtie2_index, genome: params.genome, report_dir: params.report_dir)
    include { run_baal } from "./modules/baal.nf" params(report_dir: params.report_dir)
    include { multi_qc } from "./modules/qc.nf"  params(report_dir: params.report_dir)

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
        .map {
                   runs, group_name, antigens, experiments, bedfiles, snp_files, bamfiles, index_files -> 
           [runs, group_name, antigens, experiments, bedfiles.unique(), snp_files.unique(), bamfiles, index_files]
        }

    if (params.run_baal) {
        bam_files | mergeBeds | run_baal
    }
}
