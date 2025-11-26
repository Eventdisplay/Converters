# Dockerfile to build DL3-IRFs Converter for Eventdisplay
FROM rootproject/root:6.30.02-ubuntu22.04
ENV LC_ALL=C

LABEL maintainer.name="Eventdisplay Team"
LABEL maintainer.email="gernot.maier@desy.de"

# force ubuntu to use bash for /bin/sh
RUN yes no | dpkg-reconfigure dash

## Basic packages
RUN apt-get update && \
    apt-get install --no-install-recommends -y git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /data/
ENV DATA=/data/

## Install and compile cfitsio from GitHub at tag cfitsio-4.6.3
RUN git clone https://github.com/HEASARC/cfitsio.git && \
    cd cfitsio && \
    git checkout cfitsio-4.6.3 && \
    ./configure --prefix=${DATA}/ && \
    make && \
    make install && \
    make clean && \
    cd .. && rm -rf cfitsio

ENV LD_LIBRARY_PATH="${DATA}/lib:${LD_LIBRARY_PATH}"

# Install and compile DL3 converter
RUN git clone https://github.com/Eventdisplay/Converters.git && \
    cd Converters && \
    cd DL3-IRFs && \
    make convertSensitivityFilesToFITS
