# Base image
FROM ubuntu:14.04

# No interaction during install
ARG DEBIAN_FRONTEND=noninteractive

# Install essentials
RUN apt-get -qqy update && \
    apt-get -qqy install \
    build-essential \
    libopenmpi-dev \
    openmpi-bin \
    libgsl0-dev \
    git \
    curl

# install python and jupyter
RUN apt-get update && apt-get install -qqy \
    python \
    python-dev

RUN curl -fSsL -O https://bootstrap.pypa.io/get-pip.py && \
    python get-pip.py && \
    rm get-pip.py

RUN pip --no-cache-dir install \
    setuptools \
    numpy \
    jupyter

# Load fastpm into docker
WORKDIR /fastpm
COPY . .

# Build fastpm
RUN cd /fastpm && \
	echo "CC = mpicc\nOPENMP = -fopenmp\nCPPFLAGS = -DFASTPM_FFT_PRECISION=32\nOPTIMIZE = -O3 -g" > Makefile.local && \
	make

# Expose an address for jupyter notebooks
EXPOSE 8888

# Folder which can be used for connecting data directories
WORKDIR /workspace

ENTRYPOINT ["/bin/bash"]
