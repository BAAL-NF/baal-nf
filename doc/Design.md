# Pipeline Design

# Background and Motivation

# Design Decisions
- Coarse-grained QC
- How to handle replicates

# How to add features to the pipeline

# Future Work

- Migrate BaalChIP to MCMC for whole monte carlo run
- Fix BaalChIP unit tests
- Motif calling and filtering
- Refactor pipeline to reduce number of parameters being passed around

# Stages

The pipeline consists of a number of stages. See the flowchart at the end of this document for a brief overview.

## QC

### Initial FastQC

Processes:
- `qc::fastQC`
- `qc::getFastqcResult`

This stage is a coarse filter, intended to throw out input that for whatever reason is unparseable or of so low quality that it's not fit for further processing. Most FastQC filters are disabled. The ones that are used are

- `quality_sequence`
- `tile`
- `sequence_length`

Full configuration can be found in [data/before_limits.txt](/data/before_limits.txt).

`getFastqcResult` is a simple bash script that will print FAIL to stdout if one or more filters have failed.
The `filter_fastq` workflow uses this to discard any fastq files that do not meet the threshold.

### TrimGalore

Process: `fastq::trimGalore`.

This runs `TrimGalore` on the fastq files, to perform best-guess adapter trimming.
No optional arguments are passed, but `TrimGalore` is called appropriately for single-stranded and double-stranded fastq files.

Trim galore produces a report, which is collated alongside other reports and provided to multiQC later.

### Fetching FastqScreen files

Process: `qc::fetchFastqScreenFiles`

This process stores the reference genomes required by `FastqScreen` in a configurable folder as designated by `fastq_screen_cache`.
This step is only run once per user.
If the folder already exists, this process will not run.

### FastQ-screen

Process: `qc::fastqScreen`

`FastqScreen` looks for sample contamination by randomly sampling reads and seeing if they match reference genomes from a wide variety of species.
It produces a report which is collated by MultiQC.

### Filtering based on FastqScreen results

Process: `qc::getFastqScreenResult`

This is a simple python script that looks through all non-human reference genomes in the FastQ-screen report.
If, for any genome, less than a configurable percentage (`max_acceptable_unmapped`) of the input reads are not mapped (i.e. if a large amount of the input reads are successfully mapped to this genome), the input file is discarded as likely contaminated. 

### FastQC

`FastQC` and filtering is done exactly as for [the initial fastqc filtering](#initial-fastqc), but with a different configuration file, namely [data/after_limits.txt](/data/after_limits.txt)

These are the same processes, the workflow is imported twice into [`main.nf`](/main.nf) as `filter_fastq_before` and `filter_fastq_after`.

## Mapping

### Bowtie2

Process: `fastq::createBam`
Maps using the `hg19` reference genome from igenomes by default.
Note that if you change the reference genome you will also need to modify the use of BaalChIP to use different blacklists when doing variant calling.

### picard

Process: `fastq::markDuplicates`
Marks duplicate reads.

## BaalChIP

### Create Config File

Process:`baal::createSampleFile`

Creates a configuration file for use by BaalChIP (see e.g [the example files included with BaalChiP](https://github.com/InesdeSantiago/BaalChIP/blob/master/inst/test/exampleChIP.tsv)).


### Filter Reads and Count Instance of Each Allele

Process: `baal::baalProcessBams`

The first of two stages using the BaalChIP package. The following is done in this step

- Discards all reads below a certain read quality threshold
- Filters out reads overlapping blacklisted areas
- Filters out duplicate reads
- Constructs a table of heterozygous SNPs, with counts of how many reads have been observed for each allele
- Filters out any SNPs for which only one of the two alleles have been observed

### Get Allele-specific Binding Sites

Process: `baal::getAsb`

Runs a markov-chain monte carlo simulation based on the results from the previous step, along with input allelic ratios.
This outputs the final allele-specific binding results.

Also produces a BaalChIP report, in the form of a rendered PDF.

## Post-processing

### multiQC

Collate all QC reports except for the initial run of FastQC into a single report. Inputs from
- FastQC
- Picard
- Bowtie2
- FastQ-Screen

### GAT

Genome Analysis Toolkit. This runs enrichment analysis, comparing called ASB sites to all heterozygous SNPs in a given cell line, and looking for enrichments in e.g. enhancers, promoters, promoter-flanking regions, etc.

## Motif Calling




# Output Files and Structure

# Debugging tips and tricks



# Flowchart

```mermaid
flowchart TD
    p0((Channel.fromPath))
    p1([splitCsv])
    p2([map])
    p3([multiMap])
    p4(( ))
    p5[filter_fastq_before:fastqc_before_trimming:fastQC]
    p6[filter_fastq_before:fastqc_before_trimming:getFastqcResult]
    p7([filter])
    p8([map])
    p9([map])
    p10([join])
    p11([join])
    p12([groupTuple])
    p13([map])
    p14[reportFastQC]
    p15(( ))
    p16[trimGalore]
    p17(( ))
    p18[filter_fastq_after:fastqc_after_trimming:fastQC]
    p19[filter_fastq_after:fastqc_after_trimming:getFastqcResult]
    p20([filter])
    p21([map])
    p22([map])
    p23[filter_fastq_after:fastq_screen:fetchFastqScreenFiles]
    p24[filter_fastq_after:fastq_screen:fastqScreen]
    p25[filter_fastq_after:fastq_screen:getFastqScreenResult]
    p26([filter])
    p27([map])
    p28([join])
    p29([join])
    p30([mix])
    p31([join])
    p32([groupTuple])
    p33([map])
    p34([transpose])
    p35([multiMap])
    p36(( ))
    p37[create_bam:createBam]
    p38(( ))
    p39[create_bam:index]
    p40([mix])
    p41([mix])
    p42([groupTuple])
    p43([map])
    p44([join])
    p45([groupTuple])
    p46([map])
    p47[multi_qc:multiQC]
    p48([join])
    p49([groupTuple])
    p50([multiMap])
    p51([map])
    p52[no_peak:getFragmentSize]
    p53(( ))
    p54[no_peak:parseFragmentSize]
    p55[no_peak:bamToBed]
    p56(( ))
    p57[no_peak:profileMotifs]
    p58([join])
    p59[no_peak:getMotifs]
    p60(( ))
    p61[mergeBeds]
    p62([join])
    p63[run_baal:createSampleFile]
    p64([map])
    p65([join])
    p66[run_baal:baalProcessBams]
    p67(( ))
    p68[run_baal:baalGetASB]
    p69(( ))
    p70(( ))
    p71[process_results:overlapPeaks]
    p72([join])
    p73(( ))
    p74[process_results:makeGatBedFiles]
    p75(( ))
    p76[process_results:runGat]
    p77([map])
    p78([map])
    p79([map])
    p80([join])
    p81([join])
    p82([map])
    p83([join])
    p84([collectFile])
    p85(( ))
    p0 --> p1
    p1 --> p2
    p2 --> p3
    p3 -->|fastq_list| p5
    p3 -->|metadata| p11
    p4 -->|fastqc_conf| p5
    p5 --> p6
    p6 --> p7
    p7 --> p8
    p8 -->|fastq_list| p10
    p5 --> p9
    p9 -->|report| p11
    p3 -->|fastq_list| p10
    p10 -->|result| p16
    p11 --> p12
    p12 --> p13
    p13 --> p14
    p14 -->|reports| p15
    p16 -->|fastq_ch| p18
    p16 --> p40
    p17 -->|fastqc_conf| p18
    p18 --> p19
    p19 --> p20
    p20 --> p21
    p21 -->|fastq_list| p28
    p18 --> p22
    p22 -->|report| p30
    p23 --> p24
    p16 -->|fastq_ch| p24
    p24 --> p25
    p24 -->|report| p30
    p25 --> p26
    p26 --> p27
    p27 -->|result| p28
    p28 --> p29
    p16 -->|fastq_ch| p29
    p29 -->|fastq_files| p31
    p30 -->|report| p40
    p3 -->|metadata| p31
    p31 --> p32
    p32 --> p33
    p33 --> p34
    p34 --> p35
    p35 -->|fastq_files| p37
    p35 -->|metadata| p44
    p36 -->|index_files| p37
    p37 --> p39
    p37 -->|report| p41
    p37 --> p38
    p39 --> p48
    p40 --> p41
    p41 --> p42
    p42 --> p43
    p43 -->|report| p44
    p44 --> p45
    p45 --> p46
    p46 --> p47
    p47 -->|multiqc_results| p77
    p35 -->|metadata| p48
    p48 -->|merged_data| p49
    p49 --> p50
    p50 -->|snp_files| p72
    p50 --> p62
    p50 --> p61
    p48 -->|merged_data| p51
    p51 -->|input| p52
    p52 --> p54
    p53 -->|parse_script| p54
    p54 --> p58
    p51 -->|input| p55
    p55 --> p57
    p56 -->|genome| p57
    p57 --> p58
    p58 --> p59
    p59 --> p60
    p61 --> p62
    p62 -->|baal_groups| p63
    p63 --> p65
    p62 -->|baal_groups| p64
    p64 --> p65
    p65 --> p66
    p66 --> p68
    p67 -->|report_md| p68
    p68 -->|baal_results| p71
    p68 -->|asb_report| p78
    p68 --> p69
    p70 -->|script| p71
    p71 -->|overlap_peaks_results| p72
    p72 --> p74
    p73 -->|bedfile_script| p74
    p74 --> p76
    p75 -->|annotations| p76
    p76 -->|gat_results| p82
    p77 -->|multiqc_flat| p81
    p78 -->|asb_output| p80
    p71 -->|overlap_peaks_results| p79
    p79 -->|overlap_peaks_results| p80
    p80 --> p81
    p81 -->|combined_results| p83
    p82 -->|gat_results| p83
    p83 -->|combined_results| p84
    p84 --> p85
```