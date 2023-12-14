# Dockerfile to build DL2 Converter for Eventdisplay
FROM python:3.10-slim
ENV LC_ALL=C

LABEL maintainer.name="Eventdisplay Team"
LABEL maintainer.email="gernot.maier@desy.de"
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y --no-install-recommends wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV PATH /conda/bin:$PATH
ENV PYTHONPATH=$PYTHONPATH:"/"

COPY Converters/DL2/environment.yml environment.yml

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /conda && \
    rm ~/miniconda.sh
RUN conda env update -n base --file environment.yml && \
    conda clean --all

WORKDIR /workdir/

SHELL ["/bin/bash", "-c"]

RUN source /root/.bashrc && \
    conda init bash
