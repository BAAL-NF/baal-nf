This project contains a pipeline for calculating allele-specific binding affinity from a pre-downloaded set of fastq files using BaalChIP.

# Requirements
- If using anaconda, you will need to have bioconda in your channel configuration
```
conda config --add channels defaults
conda config --add channels bioconda
```

# Design

Baal-NF is a pipeline for computing allele-specific binding variants (ASBVs) in cancer and non-cancer cell lines.
It features quality checking and reporting using fastqc, fastq-screen and multiqc, and uses the [BaalChIP](https://github.com/InesdeSantiago/BaalChIP) R package for the detection of ASBVs.
The pipeline is run using nextflow, and currently supports execution using native environments, anaconda, docker and singularity.

For more details on the design and each process in the pipeline, see [the design docs](doc/Design.md)

# Usage and Configuration

The hg19 index for bowtie2 is automatically downloaded from the AWS iGenomes s3 bucket. Should you wish to download it manually, set the bowtie2_index parameter in your pipeline configuration or on the command line to point to the folder containing the required index files, i.e.

```groovy
params.bowtie2_index = "/path/to/bowtie2/index"
params.genome = "index_file_basename" // e.g. hg19
```
**Warning**: Note that BaalChIP uses blacklists based on hg19, and currently no other reference genomes are supported. You must use hg19 for all input files to this pipeline.

In a similar manner, reference genomes for `fastq-screen` are downloaded automatically unless their location are provided in `params.fastq_screen_cache`. The location provided in this variable must contain a folder titled `FastQ_Screen_Genomes`, which again must contain a copy of `fastq_screen.conf`. This is the folder structure which is automatically generated by `fastq-screen` when run the `--get_genomes` command line option.

## Configuration profiles

The pipeline currently supports docker, singularity and anaconda for dependency management. We recommend using docker or singularity for reproducibility, though anaconda may be used if your computing infrastructure does not support containers. 

If using anaconda, you will need to manually create an anaconda environment and install our modified version of BaalChIP manually from https://git.ecdf.ed.ac.uk/oalmelid/BaalChIP. The path to this environment should then be provided as the parameter `baal_chip_env` (see below).

The following configuration profiles are available for general use

| Name | Description |
| :----- | :----- |
| docker | Use the docker container engine |
| singularity | Use the singularity container engine |
| conda | Use anaconda environments |
| gls | Execute using the [Google Life Sciences API](https://www.nextflow.io/docs/latest/google.html#cloud-life-sciences) |
| test | Run the pipeline test set, can be used in combination with any of the above environment profiles |
| eddie | Platform-specific profile for the University of Edinburgh eddie cluster |

## Input file format

baal-nf requires a listing of input files in CSV format, with the following named columns

| Column name | Example | Description |
| :----- | :----- | :----- |
| cell_line | GM12878 | Cell line |
| transcription_factor | ESR1 | Transcription Factor |
| run | SRR123456 | Unique identifier of sequencing run |
| fastq_1 | `/scratch/mydata/experiment/GM12878_ESR1_1.fastq.gz` | FASTQ file 1, use this field if using single-end sequencing |
| fastq_2 | `/scratch/mydata/experiment/GM12878_ESR1_2.fastq.gz` | FASTQ file 2, leave blank if using single-end sequencing |
| bed_file | `/scratch/mydata/bedfiles/GM12878_ESR1.bed` | Path to bed file containing peak calls for the sequencing run |
| snp_list | `/scratch/het_snps/GM12878_hetSNP.txt` | TSV file containing het SNPs and RAF in the [format expected by BaalChIP](https://github.com/InesdeSantiago/BaalChIP/blob/master/inst/test/GM12891_hetSNP.txt) |

All files may be either local paths, HTTP URLs or S3 uris. Nextflow will automatically stage any files not present in the local filesystem.

Other columns may be included in the input file for bookkeeping purposes. `baal-nf` will ignore these fields in its output.

## Configuration options

Some configuration can be set using nextflow's usual custom parameters, either on the command line using double-dashed command line options, or in a nextflow configuration file in the `params` scope. The full list of configuration options is as follows

| Parameter | Default | Description | Required? | 
| :----- | :----- | :----- | :-- |
| sample_file | No default| Path to input file as specified [in the previous section](#input-file-format) | Yes 
| fastq_screen_cache | `"${launchDir}/assets/fastq_screen"` | Path to folder containing fastq screen reference genomes, as generated by running the command fastq_screen --get_genomes. The folder must contain a folder named `FastQ_Screen_Genomes`, which again must contain a copy of `fastq_screen.conf`. The pipeline will automatically download these genomes if not specified. | No
| bowtie2_index | `"s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/"` | Location of bowtie2 index files if using local cache | No
| confidence_levels | `[0.95, 0.99]` | Credible intervals to be output by BaalChIP. The first entry will be used by BaalChIP to assess whether a given site is an ASB site or not | Yes
| dedup_umi | `false` | Use UMI_tools to deduplicate mapped files based on UMIs. This also disables BaalChIP's 1-allele filter | No
| umi_tools_options | `""` | Additional command-line parameters for umi_tools dedup | No
| genome | `"genome"`| Name of the reference genome used for mapping. This should correspond to the file name for your local copy of hg19, if changed. | No 
| fastqc_conf_pre | `"${workflow.projectDir}/data/before_limits.txt"` | `fastqc` configuration used for pre-screening| No
| fastqc_conf_post | `"${workflow.projectDir}/data/after_limits.txt"`| `fastqc` configuration used after adapter trimming | No
| report_dir | `"${workflow.launchDir}/reports/"`| Directory to place all reports in, defaults to a subfolder named `reports` in the launch directory. | No
| run_gat | `true` | Whether to run GAT enrichment analysis against the ENSEMBL genome annotations | No
| baal_chip_env | null | Only required when using the conda profile. Path to conda environment containing BaalChIP | Yes, if using the conda profile
| save_baal_objects | false | Export BaalChIP objects as RDS files | No
| motif_kmer_length | 8 | kmer length used by NoPeak when doing motif calling. 8-12 is recommended. | No

# Output Files and Structure

All output files are in the folder specied by `params.report_dir`, by default this is the subdirectory `reports` within the launch directory for the pipeline.
The report directory contains the following subfolders.

## `asb`

This folder contains the output of BaalChIP.
See [the BaalChIP documentation](https://git.ecdf.ed.ac.uk/oalmelid/BaalChIP/-/blob/dev/vignettes/BaalChIP.Rmd) for details on this.

## `baal_objects`

Only available when `params.save_baal_objects` is set to `true`.
RDS files produced by BaalChIP when counting reads overlapping each allele of a SNP.
These files are mainly useful for debugging.

## `baal_reports`

BaalChIP run summary as a pdf.
This contains information about corrections for RAF and CNV, as well as details on filtering done by BaalChIP.

## `bed_files`

The result of taking the union of all peak call bed files in the source directory.

## `enrichment`

Output from GAT, these calculate the relative enrichment of ASB sites in various regions of the genome, compared to a background of all heterozygous SNPs in the given cell line.

## `logs`

Logs from bowtie2 and umi_tools (`dedup`).

## `motifs`

Motif calls and kmer profiles as produced by NoPeak.

## `multiQC`

Collated QC reports from multiQC.
This report combines the reports from fastqc, fastq-screen, bowtie2 and picard into one report.

## `pipeline_info`

Nextflow reports, including pipeline performance statistics, timelines, a trace report, a diagram of the DAG and a CSV file summarising the completed sets of cell lines and transcription factors.

## `preFastQC`

FastQC reports for the initial QC run.

## `samples`

Sample files provided to BaalChIP. These contain details on which sequencing runs were included in the ASB calculation.