#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

if (params.sample_file.isEmpty()) {
    println('Error: no sample file provided')
    exit 1
}

if (workflow.profile.split(',').contains('conda') && (! params.baal_chip_env)) {
    println('Error: conda configuration profile selected, but baal_chip_env not specified.')
    println('Please specify the path to an anaconda environment containing BaalChIP.')
    println('See https://git.ecdf.ed.ac.uk/oalmelid/baal-nf/blob/master/README.md#configuration-options for details.')
    exit 1
}

// Import a CSV file with all sample identifiers
workflow import_samples {
    main:
    Channel
        .fromPath(params.sample_file, checkIfExists: true)
        .splitCsv(header:true)
        .map {
            row -> tuple(
                row.run,
                "${row.cell_line}_${row.transcription_factor}",
                row.transcription_factor,
                row.background_run,
                ([row.fastq_1, row.fastq_2] - '').collect {
                    path -> file(path, checkIfExists: true) 
                },
                ([row.background_1, row.background_2] - '').collect {
                    path -> file(path, checkIfExists: true)
                },
                file("${row.bed_file}", checkIfExists: true),
                file("${row.snp_list}", checkIfExists: true))
        } .multiMap {
            run, group, transcription_factor, background_run, fastq_files, background_files, bed_file, snp_file ->
            fastq: [run, group, fastq_files, "true"]
            background: [background_run, group, background_files, "false"]
            metadata: [run, group, transcription_factor, background_run, bed_file, snp_file]
        }
        .set { srr_ch }

    emit:
    fastq = srr_ch.fastq
    background = srr_ch.background.unique()
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
        .map {
            runs, tag, transcription_factor, bed_files, snp_files, fastq_files ->
            num_samples = runs.size()
            [ runs, groupKey(tag, num_samples), transcription_factor, bed_files, snp_files, fastq_files ]
        }
        .transpose()
        .multiMap {
            run, key, transcription_factor, bed_file, snp_file, fastq_files ->
            fastq: [ run, fastq_files ]
            metadata: [ run, key, transcription_factor, bed_file, snp_file ]
        }
        .set { filtered_data }

    emit:
    fastq = filtered_data.fastq
    metadata =  filtered_data.metadata
}

// Filtering stages. These are done before and after TrimGalore.
workflow filter_fastq_before {
    include { filter_fastq as fastqc_before_trimming } from './modules/qc.nf' \
        addParams(fastqc_conf: params.fastqc_conf_pre)
    take:
    fastq_list

    main:
    fastqc_before_trimming(fastq_list)
    result = fastqc_before_trimming.out.fastq_list.join(fastq_list)

    emit:
    result
    report = fastqc_before_trimming.out.report
}

// After TrimGalore we also run FastQ-screen to check for contamination
workflow filter_fastq_after {
    include { filter_fastq as fastqc_after_trimming; fastq_screen } from './modules/qc.nf' \
        addParams(fastqc_conf: params.fastqc_conf_post)

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
    result
    report
}

// Dummy process to export FastQC files from the pre-trimming run
process reportFastQC {
    publishDir("${params.report_dir}/preFastQC/${key}/", mode: 'copy')

    input:
    tuple val(key), path(reports)

    output:
    file reports

    script:
    'exit 0'
}

workflow {
    include { trimGalore; create_bam; mergeBeds } from './modules/fastq.nf'
    include { run_baal } from './modules/baal.nf'
    include { multi_qc } from './modules/qc.nf'
    include { process_results; create_report } from './modules/analysis.nf'

    // Load CSV file
    import_samples()

    // Pre-trimming filtering step - run on both samples and gDNA files
    filter_fastq_before(import_samples.out.fastq.mix(import_samples.out.background))

    import_samples.out.metadata
        .join(filter_fastq_before.out.report)
        .groupTuple(by: 1)
        .map {
            run, group, transcription_factor, bed_file, snp_file, reports ->
            [group, reports.flatten()]
        } | reportFastQC

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

    // Maybe use the count_fastq metadata in stead
    multi_qc(count_fastq.out.metadata, reports)

    // Regroup the bam files with their associated metadata and run baal chip
    count_fastq.out.metadata
        .join(create_bam.out.bamfile)
        .groupTuple(by: 1)
        .multiMap {
            runs, group_name, antigens, bedfiles, snp_files, bamfiles, index_files -> 
            
            snp_files : [group_name, snp_files.unique()]
            bed_files : [group_name, bedfiles.unique()]
            baal_files : [group_name, runs, antigens, snp_files.unique(), bamfiles, index_files] }
        .set { group_ch }

    mergeBeds(group_ch.bed_files)
    run_baal(mergeBeds.out.join(group_ch.baal_files))

    process_results(run_baal.out.asb, group_ch.snp_files)

    create_report(run_baal.out.report, multi_qc.out, process_results.out.overlap_peaks, process_results.out.gat)
}
