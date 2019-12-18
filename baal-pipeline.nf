#!/usr/bin/env nextflow

params.staging_root = "${workflow.launchDir}/test/staging/runs"
params.experiments = "${params.staging_root}/runs.csv"
params.report_dir = "${workflow.launchDir}/test/reports/"
params.genome = "hg19"
params.picard_cmd = "picard"
params.fastq_screen_conf = ""
params.mpiflags = ""

Channel
    .fromPath(params.experiments)
    .splitCsv(header:true)
    .map({row -> tuple("${row.cell_line}_${row.transcription_factor}",
                        row.cell_line,
                        row.transcription_factor,
                        row.experiment,
                        row.run,
                        file("${params.staging_root}/${row.fastq_folder}*.fastq.gz"),
                        file("${params.staging_root}/${row.bed_file}"),
                        file("${params.staging_root}/${row.snp_list}"))})
    .groupTuple()
    .map({tag, cell_lines, transcription_factors, experiments, runs, fastq_files, bed_files, snp_files ->
        num_samples = runs.size()
        [ groupKey(tag, num_samples), cell_lines, transcription_factors, experiments, runs, fastq_files, bed_files, snp_files ]
    })
    .transpose()
    .transpose()
    .fork { key, cell_line, transcription_factor, experiment, run, fastq_file, bed_file, snp_file ->
        filename = fastq_file.name
        srr_ch: tuple(run, filename, fastq_file)
        baal_ch: tuple(run, key, cell_line, transcription_factor, experiment, bed_file, snp_file)
    }
    .set { samples }

samples.srr_ch.into{ fastqc_input; for_trim_galore }

process fastQC {
    label 'fastq'
    publishDir("${params.report_dir}/${sampleID}", mode: "move")

    input:
    // I don't like having to pass a filename in here to coerce nextflow
    // into linking the file. I'll have to look for a better solution to this.
    tuple sampleID, filename, file("${filename}") from fastqc_input

    output:
    file "*_fastqc*"

    script:
    """
    fastqc ${filename}
    """
}

// group multiple fastq files for the same run, since they will need
// to be processed together by trimGalore and subsequent tasks
trim_galore_input = for_trim_galore.groupTuple()

process trimGalore {
    label 'fastq'
    publishDir("${params.report_dir}/${sampleID}", mode: "move", pattern: "*report*")

    input:
    set sampleID, filename, file(files) from trim_galore_input

    output:
    set sampleID, file("${sampleID}_tg*") into trim_galore_out
    file '*report*'

    script:
    switch (files) {
         case nextflow.processor.TaskPath:
         return """trim_galore --basename ${sampleID}_tg ${files}"""

         case nextflow.util.BlankSeparatedList:
         return """trim_galore --basename ${sampleID}_tg --paired ${files}"""

         default:
         println("Error getting files for sample ${sampleID}, exiting")
         return "exit 1"
    }
}

trim_galore_out.into {
    fastq_screen_input
    bowtie_input
}

// process fastqScreen {
//     label 'fastq'
//     publishDir("${params.report_dir}/${sampleID}", mode: "move", pattern: "*screen*")

//     input:
//     set sampleID, file(trimmed) from fastq_screen_input

//     output:
//     file trimmed into fastq_screen_out
//     file '*screen*'

//     script:
//     optargs = ''

//     if (params.fastq_screen_conf != null && !params.fastq_screen_conf.isEmpty()){
//         optargs += "--conf ${params.fastq_screen_conf}"
//     }

//     switch (trimmed) {
//         case nextflow.processor.TaskPath:
//         return """fastq_screen ${optargs} ${trimmed}"""

//         case nextflow.util.BlankSeparatedList:
//         return """fastq_screen ${optargs} --paired ${trimmed}"""

//         default:
//         println("Error getting files for sample ${sampleID}, exiting")
//         return "exit 1"
//     }

// }

process align {
    label 'fastq'
    label 'bigmem'

    input:
    set sampleID, file(trimmed) from bowtie_input

    output:
    set sampleID, file("${sampleID}.sam") into aligned_sequences

    script:
    switch (trimmed) {
        case nextflow.processor.TaskPath:
        return """bowtie2 -x ${params.genome} -U ${trimmed} -S ${sampleID}.sam"""

        case nextflow.util.BlankSeparatedList:
        first = trimmed[0]
        second = trimmed[1]
        return """bowtie2 -x ${params.genome} -1 ${first} -2 ${second} -S ${sampleID}.sam"""

        default:
        println("Error getting files for sample ${sampleID}, exiting")
        return "exit 1"
    }
}

process sortAndCompress {
    label 'fastq'
    input:
    set sampleID, file(samfile) from aligned_sequences

    output:
    set sampleID, file("${samfile.baseName}.bam") into bamfiles

    """
    samtools sort ${samfile} -o ${samfile}.sorted
    samtools view -h -S -b ${samfile}.sorted > ${samfile.baseName}.bam
    """
}

process markDuplicates {
    label 'fastq'
    publishDir("${params.report_dir}/${sampleID}", mode: "move", pattern: "**.metrics")

    input:
    set sampleID, file(bamfile) from bamfiles

    output:
    set sampleID, file("${bamfile.baseName}_dedup.bam") into marked_bamfiles
    file "${bamfile.baseName}.metrics"

    """
    ${params.picard_cmd} MarkDuplicates I="${bamfile}" O="${bamfile.baseName}_dedup.bam" M="${bamfile.baseName}.metrics"
    """
}

process index {
    label 'fastq'
    input:
    set sampleID, file(bamfile) from marked_bamfiles

    output:
    set sampleID, file(bamfile), file("${bamfile}.bai") into indexed_bamfiles

    """
    samtools index ${bamfile}
    """
}

samples.baal_ch
    .join(indexed_bamfiles)
    .groupTuple(by: 1)
    .set { baal_file_ch }

process createSampleFile {
    publishDir("test/reports/samples/")

    input:
    set runs, group_name, cell_lines, antigens, experiments, bedfiles, snp_files, bamfiles, index_files from baal_file_ch

    output:
    set runs, group_name, cell_lines, antigens, experiments, bedfiles, snp_files, bamfiles, index_files, file("${group_name}.tsv") into baal_sample_file_ch

    script:
        output = "cat << EOF > ${group_name}.tsv\n"
        output += "group_name\ttarget\treplicate_number\tbam_name\tbed_name\tsampleID\n"
        0.upto(runs.size()-1, {
            replicate = it + 1
            output += "${group_name}\t${antigens[it]}\t${replicate}\t${bamfiles[it].name}\t${bedfiles[it].name}\t${runs[it]}\n"
        })
        output += "EOF\n"
        output
}

baal_sample_file_ch
    .map({
        runs, group_name, cell_lines, antigens, experiments, bedfiles, snp_files, bamfiles, index_files, sample_file ->
        [runs, group_name, cell_lines, antigens, experiments, bedfiles.unique(), snp_files[0], bamfiles, index_files, sample_file]
    })
    .set { baal_process_bams_ch }

process baalProcessBams {
    label "baal_chip"

    input:
    set runs, group_name, cell_lines, antigens, experiments, file(bedfiles), file(snp_file), file(bamfiles), file(index_files), file(sample_file) from baal_process_bams_ch

    output:
    set file("process_bams.rds"), file(snp_file) into get_asb_ch

    script:
    """
    #!/usr/bin/env Rscript
    library(BaalChIP)

    samplesheet <- "${sample_file}"

    hets <- c("${group_name}" = "${snp_file}")

    res <- new("BaalChIP", samplesheet=samplesheet, hets=hets)
    res <- alleleCounts(res, min_base_quality=10, min_mapq=15, all_hets=TRUE)
    res <- QCfilter(res,
                    RegionsToFilter=list("blacklist"=blacklist_hg19,
                                         "highcoverage"=pickrell2011cov1_hg19),
                    RegionsToKeep=list("UniqueMappability"=UniqueMappability50bp_hg19))
    res <- mergePerGroup(res)
    res <- filter1allele(res)
    saveRDS(res, file="process_bams.rds")
    """
}

process baalGetASB {
    publishDir("${params.report_dir}/asb")

    label "baal_chip"
    label "mpi"

    input:
    set file("process_bams.rds"), file(snp_file) from get_asb_ch

    output:
    file "*.csv"

    script:
    """
    #!/usr/bin/env mpirun -np 1 ${params.mpiflags} Rscript
    library(BaalChIP)
    # Read in hets from file
    res <- readRDS("process_bams.rds")
    res <- getASB(res, Iter=5000, conf_level=0.95)
    saveRDS(res, "final.rds")
    report <- BaalChIP.report(res)

    for (group in names(report)) {
            write.csv(report[[group]], paste(group,".csv", sep=""))
    }
    """
}
