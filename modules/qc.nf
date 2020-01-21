params.report_dir = "${workflow.launchDir}/test/reports/"
params.fastqc_conf = ""
params.fastq_screen_conf = ""
params.max_acceptable_unmapped = 90

process fastQC {
    label 'fastq'
    publishDir("${params.report_dir}/${run}/${report_subdir}", pattern: "*fastqc*",  mode: "copy")

    input:
    tuple run, file(fastq_files)
    val report_subdir

    output:
    tuple file("*_fastqc.zip"), file("*_fastqc.html"), run, file(fastq_files)

    script:
    config_args=[]

    if (params.fastqc_conf) {
        config_args += ["-l", "${params.fastqc_conf}"]
    }

    config_cmd = config_args.join(" ")

    """
    fastqc ${config_cmd} ${fastq_files}
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
        script += """
                  unzip ${report_file} ${report_file.baseName}/summary.txt >/dev/null
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
    options = []

    if (params.fastq_screen_conf){
        options += ["--conf", "${params.fastq_screen_conf}"]
    }
    optargs = options.join(" ")

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

process getFastqScreenResult {
    label "python"

    input:
    tuple run, file(screening_result)

    output:
    tuple run, stdout

    script:
    """
    #!/usr/bin/env python
    import pandas as pd

    # Files are provided as a space separated list by nextflow, so split on spaces
    files = "${screening_result}".split(" ")

    result = True
    for file in files:
        screening = pd.read_csv(file, sep='\t', skiprows=1, skipfooter=2, engine='python')
        result = result and (screening[screening.Genome != "Human"]["%Unmapped"].min() > ${params.max_acceptable_unmapped})

    print("pass" if result else "fail")
    """

}

workflow fastq_screen {
    get:
    fastq_ch

    main:
    screen = fastqScreen(fastq_ch)
    result = getFastqScreenResult(screen.screening_result) | filter { it[1].strip() == "pass" } | map { it -> it [0] }

    emit:
    result
}
