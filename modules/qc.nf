params.report_dir = "${workflow.launchDir}/test/reports/"
params.fastqc_conf = ""

process fastQC {
    label 'fastq'
    publishDir("${params.report_dir}/${run}/${report_subdir}", pattern: "*fastqc*",  mode: "copy")

    input:
    tuple run, file(fastq_files)
    val report_subdir

    output:
    tuple file("*_fastqc.zip"), file("*_fastqc.html"), run, file(fastq_files)

    script:
    config_args=""

    if (params.fastqc_conf) {
        config_args += "-l ${params.fastqc_conf}"
    }

    """
    fastqc ${config_args} ${fastq_files}
    """
}

process getFastqcResult {
    input:
    tuple file(report_zip), file(html_report), identifier, fastq_files

    output:
    tuple stdout, identifier, fastq_files

    script:
    script = ""
    for (report_file in report_zip) {
        script += """unzip ${report_file} ${report_file.baseName}/summary.txt >/dev/null
                     cat ${report_file.baseName}/summary.txt | awk '/^FAIL/{print \"FAIL\"}' | uniq
                  """
    }
    script
}

workflow filter_fastq {
    get:
    fastq_list
    report_subdir

    main:
    result = fastQC(fastq_list, report_subdir) | getFastqcResult | filter{ it[0].isEmpty() } | map{ it -> it[1..-1] }

    emit:
    result
}


process fastqScreen {
    label 'fastq'
    publishDir("${params.report_dir}/${run}", mode: "copy", pattern: "*screen*")

    input:
    tuple run, file(trimmed)

    output:
    tuple run, file("*screen.txt"), emit: screening_result
    file '*screen*'

    script:
    optargs = ""

    if (params.fastq_screen_conf){
        optargs += "--conf ${params.fastq_screen_conf}"
    }

    switch (trimmed) {
        case nextflow.processor.TaskPath:
        return """fastq_screen ${optargs} ${trimmed}"""

        case nextflow.util.BlankSeparatedList:
        return """fastq_screen ${optargs} --paired ${trimmed}"""

        default:
        println("Error getting files for sample ${run}, exiting")
        return "exit 1"
    }
}

workflow fastq_screen {
    get:
    fastq_ch

    main:
    result = fastqScreen(fastq_ch)

    emit:
    result.screening_result
}
