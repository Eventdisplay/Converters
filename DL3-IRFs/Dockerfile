From ubuntu:20.04

# force ubunto to use bash for /bin/sh
# (otherwise trouble with many evndisp scripts)
RUN yes no | dpkg-reconfigure dash

## Basic packages
RUN apt-get update && \
  apt-get install -y \
  bash \
  bzip2 \
  gcc \
  g++ \
  git \
  libblas-dev \
  libbz2-dev \
  libfreetype6 \
  libxpm4 \
  libtiff-dev \
  libxtst6 \
  make \
  unzip \
  vim \
  wget \
  zstd
RUN apt-get clean
WORKDIR /data/

# Install ROOT
RUN wget https://root.cern/download/root_v6.22.00.Linux-ubuntu20-x86_64-gcc9.3.tar.gz && \
    tar -xzf root_v6.22.00.Linux-ubuntu20-x86_64-gcc9.3.tar.gz && \
    rm -f root_v6.22.00.Linux-ubuntu20-x86_64-gcc9.3.tar.gz
ENV ROOTSYS /data/root
ENV LD_LIBRARY_PATH "${ROOTSYS}/lib/:${LD_LIBRARY_PATH}"
ENV PATH ${ROOTSYS}/bin:${PATH}

# Install cfitsio
RUn wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.49.tar.gz && \
    tar -xzf cfitsio-3.49.tar.gz && \
    cd cfitsio-3.49 && \
    ./configure --prefix=${WORKDIR}/ && \
    make && \
    make install 

ENV LD_LIBRARY_PATH "${WORKDIR}/lib:${LD_LIBRARY_PATH}"
