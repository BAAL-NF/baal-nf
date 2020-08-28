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

Due to their prohibitive size, you must download the hg19 index for bowtie2 and provide it to the pipeline by setting the BOWTIE2_INDEXES environment variable, i.e.

```groovy
BOWTIE2_INDEXES = "/path/to/bowtie2/index"
```

You will also need to have the BaalChIP package installed in your native environment, or in an anaconda environment in order to run. We use a modified version of BaalChIP currently available at (https://git.ecdf.ed.ac.uk/oalmelid/BaalChIP).
If your cluster supports docker or singularity, use the docker and singularity profiles. At present, baal-nf doesn't correctly handle source files, so you will need to add both the directory where nextflow has stored baal-nf and the directory with any configuration files, reference genomes and similar to the singularity or docker mount options manually, using e.g.

```groovy
singularity.runOptions = "--bind /path/to/baal-nf --bind /path/to/working/dir"
```

We intend to fix this in a future release.

# ToDo

- Automatic fetching of fastq files from the ENA archive
- Automatic export to SQL or SQLite database
- Fix singularity/docker mount points so pipeline can run
- Add a small testing set to the project
- Configurable genomes
- Automatic download of genomes
