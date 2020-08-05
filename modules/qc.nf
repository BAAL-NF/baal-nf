params.report_dir = "${workflow.launchDir}/test/reports/"
params.fastqc_conf = ""
params.fastq_screen_conf = ""
params.max_acceptable_unmapped = 90

process fastQC {
    label 'fastq'

    input:
    tuple run, file(fastq_files)

    output:
    tuple file("*_fastqc.zip"), file("*_fastqc.html"), run, file(fastq_files)

    script:
    config_args=[]

    if (params.fastqc_conf) {
        config_args += ["-l", "${params.fastqc_conf}"]
    }

    config_cmd = config_args.join(" ")
    """
    # FASTQC run ${run}
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

    main:
    result = fastQC(fastq_list) | getFastqcResult | filter{ it[0].isEmpty() } | map{ it -> it[1..-1] }

    report = fastQC.out.map {
        zip, html, run, fastq_files -> [run, [zip, html].flatten()]
    }

    emit:
      fastq_list = result
      report = report
}


process fastqScreen {
    label 'fastq'

    input:
    tuple run, file(trimmed)

    output:
    tuple run, file("*screen.txt"), emit: screening_result
    tuple run, file('*screen*'), emit: report

    script:
    options = []

    if (params.fastq_screen_conf){
        options += ["--conf", "${params.fastq_screen_conf}"]
    }
    optargs = options.join(" ")

    """
    FILES=(${trimmed})
    case \${#FILES[@]} in

        1)
        fastq_screen ${optargs} \${FILES[0]}
        ;;

        2)
        fastq_screen ${optargs} --paired \${FILES[0]} \${FILES[1]}
        ;;

        *)
        echo "Error getting files for sample ${run}, exiting"
        exit 1
    esac
    """
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
    result = getFastqScreenResult(screen.screening_result) | filter { it[1] =~ /pass/ } | map { it -> it [0] }

    emit:
    result = result
    report = screen.report
}


process multiQC {
    publishDir("${params.report_dir}/multiQC/${key}", mode: "move")

    label "fastq"
    input:
    tuple key, file(results)

    output:
    file("multiqc_*")

    script:
    """
    multiqc ${results}
    """
}

workflow multi_qc {
    get:
    metadata
    report

    main:
    metadata
        .join(report)
        .groupTuple(by: 1)
        .map {
            run, group, transcription_factor, experiment, bed_file, snp_file, reports -> [group, reports.flatten()]
        } | multiQC

    emit:
    multiQC.out
}
