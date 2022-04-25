From rootproject/root:latest

# force ubunto to use bash for /bin/sh
RUN yes no | dpkg-reconfigure dash

## Basic packages
RUN apt-get update && apt-get install -y \
   git
RUN apt-get clean
WORKDIR /data/
ENV DATA /data/

# Install and compile cfitsio
RUN wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.49.tar.gz && \
    tar -xzf cfitsio-3.49.tar.gz && \
    cd cfitsio-3.49 && \
    ./configure --prefix=${DATA}/ && \
    make && \
    make install && \
    make clean && \
    cd .. && rm -f cfitsio-3.49.tar.gz

ENV LD_LIBRARY_PATH "${DATA}/lib:${LD_LIBRARY_PATH}"

# Install and compile DL3 converter
RUN git clone https://github.com/Eventdisplay/Converters.git && \
    cd Converters && \
    cd DL3-IRFs && \
    make convertSensitivityFilesToFITS
