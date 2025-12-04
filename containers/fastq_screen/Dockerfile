FROM mambaorg/micromamba:1.5.8

COPY --chown=$MAMBA_USER:$MAMBA_USER env/environment.yml /env/env.yaml

RUN micromamba install -n base -y -f /env/env.yaml \
 && micromamba clean --all --yes

# Run apt as root to add low-level utilities
USER root
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      bzip2 tar gzip xz-utils ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# go back to the original non-root user
USER $MAMBA_USER

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV CONDA_ENV_NAME=base
ENV PATH="/opt/conda/bin:${PATH}"

# Patch for FastQScreen --get_genomes which has a broken link
RUN find /opt/conda -type f -name "fastq_screen" -exec \
    sed -i "s|www.bioinformatics.babraham.ac.uk/projects/fastq_screen/genome_locations.txt|https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/genome_locations.txt|g" {} \;
