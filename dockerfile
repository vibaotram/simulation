FROM continuumio/miniconda3:latest

# System packages
RUN apt-get update && apt-get install -yq wget autoconf automake bzip2 cmake make gcc perl \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev g++ \
    build-essential libtool git pkg-config openjdk-17-jdk \
    r-base r-base-dev libxml2-dev r-cran-devtools rsync tabix

# SLIM 4.3
RUN wget https://github.com/MesserLab/SLiM/releases/download/v4.3/SLiM.zip && \
    unzip SLiM.zip && \
    mkdir build && \
    cmake ../SLiM && \
    make slim && \
    make install slim

RUN apt-get install -yq parallel

# samtools 1.21
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 -O samtools-1.21.tar.bz2
RUN tar xvjf samtools-1.21.tar.bz2
RUN cd samtools-1.21 && \
    make && \
    make prefix=/usr/local/bin install && \
    ln -s /usr/local/bin/bin/samtools /usr/bin/samtools

# bcftools 1.21
RUN wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 -O bcftools-1.21.tar.bz2
RUN tar -xjvf bcftools-1.21.tar.bz2
RUN cd bcftools-1.21 && \
    make && \
    make prefix=/usr/local/bin install && \
    ln -s /usr/local/bin/bin/bcftools /usr/bin/bcftools
