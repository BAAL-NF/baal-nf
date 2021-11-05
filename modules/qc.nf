// Initial QC prior to processing
params.max_acceptable_unmapped = 90

process fastQC {
    input:
    tuple val(run), val(group), val(not_background), path("${run}*.fastq.gz")
    file fastqc_conf

    output:
    tuple val(run), val(group), val(not_background), path('*_fastqc.zip'), path('*_fastqc.html')

    script:
    """
    # FASTQC run ${run}
    fastqc -l ${fastqc_conf} ${run}*.fastq.gz
    """
}

process getFastqcResult {
    input:
    tuple val(run), val(group), val(not_background), path(report_zip), path(html_report)

    output:
    tuple stdout, val(run)

    script:
    script = ''
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
    ch_fastqc_conf = file(params.fastqc_conf, checkIfExists: true)
    fastq_list = fastQC(fastq_list, ch_fastqc_conf) |
                 getFastqcResult |
                 filter { result -> result[0].isEmpty() } |
                 map { result -> result[1..-1] }

    report = fastQC.out.map {
        run, group, not_background, zip, html -> [run, group, not_background, [zip, html].flatten()]
    }

    emit:
      fastq_list
      report
}

process fetchFastqScreenFiles {
    storeDir params.fastq_screen_cache

    output:
    path "FastQ_Screen_Genomes"

    shell:
    '''
    fastq_screen --get_genomes
    # Strip out absolute path to accomodate staging into working directory
    sed -i "s#${PWD}/##g" FastQ_Screen_Genomes/fastq_screen.conf
    '''
}

process fastqScreen {
    label 'parallel'

    input:
    tuple val(run), val(group), val(not_background), path("${run}*.fastq.gz")
    path fastq_screen_folder
    output:
    tuple val(run), path('*screen.txt'), emit: screening_result
    tuple val(run), path('*screen*'), emit: report

    script:
    options = ['--aligner', 'bowtie2']
    options += ['--conf', "${fastq_screen_folder}/fastq_screen.conf"]

    if (task.cpus > 1) {
        options += ['--threads', "${task.cpus}"]
    }

    optargs = options.join(' ')

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
    input:
    tuple val(run), path(screening_result)

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
        if (screening[screening.Genome != "Human"]["%Unmapped"].min() < ${params.max_acceptable_unmapped}):
            result = False
            break

    print("pass" if result else "fail")
    """
}

workflow fastq_screen {
    take:
    fastq_ch

    main:
    // Fetch fastq-screen files if cache directory is empty
    // This will only work if the directory contains a folder called FastQ_Screen_Genomes
    // That folder will in turn need to contain the file fastq_screen.conf
    screen = fastqScreen(fastq_ch, fetchFastqScreenFiles())
    result = getFastqScreenResult(screen.screening_result) |
             filter { result -> result[1] =~ /pass/ } |
             map { result -> result [0] }

    emit:
    result
    report = screen.report
}

process multiQC {
    publishDir("${params.multiqc_report_dir}/${key}", mode: 'copy')

    input:
    tuple val(key), path(results)

    output:
    tuple val(key), path('multiqc_*')

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
            run, group, transcription_factor, bed_file, snp_file, reports -> [group, reports.flatten()]
        } | multiQC

    emit:
    multiQC.out
}
