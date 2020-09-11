params.max_acceptable_unmapped = 90

process fastQC {
    label 'fastq'

    input:
    tuple val(run), file("${run}*.fastq.gz")
    file fastqc_conf

    output:
    tuple val(run), file("*_fastqc.zip"), file("*_fastqc.html")

    script:
    """
    # FASTQC run ${run}
    fastqc -l ${fastqc_conf} ${run}*.fastq.gz
    """
}

process getFastqcResult {
    input:
    tuple val(run), file(report_zip), file(html_report)

    output:
    tuple stdout, val(run)

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
    take:
    fastq_list

    main:
    ch_fastqc_conf = file(params.fastqc_conf)
    filtered = fastQC(fastq_list, ch_fastqc_conf) | getFastqcResult | filter{ it[0].isEmpty() } | map{ it -> it[1..-1] }

    report = fastQC.out.map {
        run, zip, html -> [run, [zip, html].flatten()]
    }

    emit:
      fastq_list = filtered
      report = report
}


process fastqScreen {
    label 'fastq'
    label 'parallel'

    input:
    tuple val(run), path("${run}*.fastq.gz")
    file fastq_screen_conf
    output:
    tuple val(run), file("*screen.txt"), emit: screening_result
    tuple val(run), file('*screen*'), emit: report

    script:
    options = ["--aligner", "bowtie2"]
    options += ["--conf", "${params.fastq_screen_conf}"]
    
    if (task.cpus > 1) {
        options += ["--threads", "${task.cpus}"]
    }
    
    optargs = options.join(" ")

    """
    FILES=(${run}*.fastq.gz)
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
    tuple val(run), file(screening_result)

    output:
    tuple val(run), stdout

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
    take:
    fastq_ch

    main:
    ch_fastq_screen_conf = file(params.fastq_screen_conf)
    screen = fastqScreen(fastq_ch, ch_fastq_screen_conf)
    result = getFastqScreenResult(screen.screening_result) | filter { it[1] =~ /pass/ } | map { it -> it [0] }

    emit:
    result = result
    report = screen.report
}


process multiQC {
    publishDir("${params.report_dir}/multiQC/${key}", mode: "move")

    label "fastq"
    input:
    tuple val(key), file(results)

    output:
    file("multiqc_*")

    script:
    """
    multiqc ${results}
    """
}

workflow multi_qc {
    take:
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
