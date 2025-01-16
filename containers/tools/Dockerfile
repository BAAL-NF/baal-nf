FROM nfcore/base:1.12.1

RUN apt-get update && apt-get -y install unzip

COPY env/* /env/
RUN conda env update -n root --file /env/environment.yml && conda clean -a
