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
The pipeline is run using nextflow, and currently supports native and anaconda environments only.

![Pipeline flow](img/baal_pipeline.png)

# Usage and Configuration

The hg19 index for bowtie2 is automatically downloaded from the AWS iGenomes s3 bucket. Should you wish to download it manually, set the bowtie2_index parameter in your pipeline configuration or on the command line to point to the folder containing the required index files, i.e.

```groovy
params.bowtie2_index = "/path/to/bowtie2/index"
params.genome = "index_file_basename" // e.g. hg19
```
**Warning**: Note that BaalChIP uses blacklists based on hg19, and currently no other reference genomes are supported. You must use hg19 for all input files to this pipeline.

You will also need to have the BaalChIP package installed in your native environment, or in an anaconda environment in order to run. We use a modified version of BaalChIP currently available at (https://git.ecdf.ed.ac.uk/oalmelid/BaalChIP).
If your cluster supports docker or singularity, use the docker and singularity profiles. At present, baal-nf has no way to handle reference genomes and configuration for fastq-screen, so you will need to provide a path to a fastq-screen configuration file, and also make sure your container engine mounts the required indexes, using e.g.

```groovy
singularity.runOptions = "--bind /path/to/fastq-screen/indexes"
```
or
```groovy
docker.runOptions = "--mount type=bind,source=/path/to/fastq-screen/indexes,target=/path/to/fastq-screen/indexes,readonly"
```

## Input file format

baal-nf requires a listing of input files in CSV format, with the following named columns

| Column name | Example | Description |
| :----- | :----- | :----- |
| cell_line | GM12878 | Cell line |
| transcription_factor | ESR1 | Transcription Factor |
| experiment | SRX123456 | Unique identifier of the sequencing experiment|
| run | SRR123456 | Unique identifier of sequencing run |
| fastq_folder | `/scratch/mydata/folder/with/fastq/files` | Folder containing fastq files for this sequencing run |
| bed_file | `/scratch/mydata/bedfiles/GM12878_ESR1.bed` | Path to bed file containing peak calls for the sequencing run |
| snp_list | `/scratch/het_snps/GM12878_hetSNP.txt` | TSV file containing het SNPs and RAF in the [format expected by BaalChIP](https://github.com/InesdeSantiago/BaalChIP/blob/master/inst/test/GM12891_hetSNP.txt) |

All files may be either local paths, HTTP URLs or S3 uris. Nextflow will automatically stage any files not present in the local filesystem.

## Configuration options

Some configuration can be set using nextflow's usual custom parameters, either on the command line using double-dashed command line options, or in a nextflow configuration file in the `params` scope. The full list of configuration options is as follows

| Parameter | Default | Description | Required? | 
| :----- | :----- | :----- | :-- |
| sample_file | No default| Path to input file as specified [in the previous section](##input-file-format) | Yes 
| fastq_screen_conf | No default | Fastq-screen configuration file | Yes
| bowtie2_index | `"s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/"` | Location of bowtie2 index files if using local cache | No
| genome | `"genome"`| Name of the reference genome used for mapping. This should correspond to the file name for your local copy of hg19, if changed. | No 
| fastqc_conf_pre | `"${workflow.projectDir}/data/before_limits.txt"` | `fastqc` configuration used for pre-screening| No
| fastqc_conf_post | `"${workflow.projectDir}/data/after_limits.txt"`| `fastqc` configuration used after adapter trimming | No
| report_dir | `"${workflow.launchDir}/reports/"`| Directory to place all reports in, defaults to a subfolder named `reports` in the launch directory. | No
| run_gat | `true` | Whether to run GAT enrichment analysis against the ENSEMBL genome annotations | No

# ToDo

- Automatic export to SQL or SQLite database
- Fix singularity/docker mount points for fastq screen so pipeline can run without the user having to configure mount points
- Configurable genomes
