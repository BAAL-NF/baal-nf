# syntax=docker/dockerfile:experimental
FROM r-base:4.0.2

COPY env /env
COPY scripts /scripts

RUN apt-get update &&\
    cat /env/apt_packages | xargs apt-get install -y &&\
    apt-get clean &&\
    apt-get autoremove -y &&\
    rm -rf /tmp/* &&\
    rm -rf /var/lib/apt/lists/*

# Install R prerequisites
RUN Rscript --verbose scripts/install_prerequisites.R

# Hack to make sure we install BAAL_CHIP every time, since docker is 
# blind to upstream repo changes.
ARG BUILD_FORCE=1
# Install BaalChIP
RUN --mount=type=secret,id=gitlab,dst=/etc/gitlab.sh \
    . /etc/gitlab.sh &&\
    RMPI_TYPE="OPENMPI" Rscript --verbose /scripts/install_baal_chip.R &&\
    rm -rf /tmp/*

