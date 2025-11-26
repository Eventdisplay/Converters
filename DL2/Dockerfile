# Dockerfile to build DL2 Converter for Eventdisplay
FROM mambaorg/micromamba:latest
ENV LC_ALL=C

LABEL maintainer.name="Eventdisplay Team"
LABEL maintainer.email="gernot.maier@desy.de"
ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /workdir/

COPY Converters/DL2/environment.yml /tmp/environment.yml

# Create environment using micromamba
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean -a -y

ENV PATH="/opt/conda/bin:$PATH"
