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

You will need to have the BaalChIP package installed in your native environment, or in an anaconda environment in order to run. We use a modified version of BaalChIP currently available at (https://git.ecdf.ed.ac.uk/oalmelid/BaalChIP)

# ToDo

- Automatic fetching of fastq files from the ENA archive
- Automatic export to SQL or SQLite database
