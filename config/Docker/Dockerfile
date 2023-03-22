FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND noninteractive
ENV TZ America/Chicago

# Install required packages
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    tar \
    gzip \
    bzip2 \
    make \
    automake \
    gcc \
    g++ \
    git \
    r-base \
    python-is-python3 \
    python3-pip \
    perl \
    bcftools \
    vcftools \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install python packages
RUN pip3 install --upgrade pip && \
    pip3 install numpy pandas scipy cython pysam jupyter

# Install tidyverse package in R
RUN R -e 'install.packages("tidyverse")'

# Install bbmap
RUN wget https://downloads.sourceforge.net/project/bbmap/BBMap_38.90.tar.gz && \
    tar xvzf BBMap_38.90.tar.gz && \
    rm -f BBMap_38.90.tar.gz && \
    mv bbmap /usr/local/bin/

# Install vsearch
RUN wget https://github.com/torognes/vsearch/releases/download/v2.17.0/vsearch-2.17.0-linux-x86_64.tar.gz && \
    tar xvzf vsearch-2.17.0-linux-x86_64.tar.gz && \
    rm -f vsearch-2.17.0-linux-x86_64.tar.gz && \
    mv vsearch-2.17.0-linux-x86_64/bin/vsearch /usr/local/bin/

# Install fastp
RUN wget https://github.com/OpenGene/fastp/archive/v0.20.1.tar.gz && \
    tar xvzf v0.20.1.tar.gz && \
    rm -f v0.20.1.tar.gz && \
    cd fastp-0.20.1 && \
    make && \
    mv fastp /usr/local/bin/ && \
    cd .. && \
    rm -rf fastp-0.20.1

# Install fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    rm -f fastqc_v0.11.9.zip && \
    chmod +x FastQC/fastqc && \
    ln -s /FastQC/fastqc /usr/local/bin/

# Install multiQC
RUN pip3 install multiqc

# Install BWA
RUN wget https://github.com/lh3/bwa/archive/v0.7.17.tar.gz && \
    tar xvzf v0.7.17.tar.gz && \
    rm -f v0.7.17.tar.gz && \
    cd bwa-0.7.17 && \
    make && \
    mv bwa /usr/local/bin/ && \
    cd .. && \
    rm -rf bwa-0.7.17

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar xvjf samtools-1.13.tar.bz2 && \
    rm -f samtools-1.13.tar.bz2 && \
    cd samtools-1.13 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.13

# Install ANGSD
# RUN git clone https://github.com/ANGSD/angsd.git && \
#     cd angsd && \
#     make HTSSRC=/usr/include/ && \
#     make install && \
#     cd ..

# Install popscle
# RUN git clone --recursive https://github.com/cooplab/popscle.git && \
#     cd popscle && \
#     make

# Install EEMS
# RUN wget https://github.com/evolbioinf/eems/archive/refs/tags/v2.1.7.tar.gz && \
#     tar -xvzf v2.1.7.tar.gz && \
#     cd eems-2.1.7 && \
#     make

# Add tools to PATH
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/local/lib/R/site-library/ 
#:/usr/local/bin/EEMS

# Set the working directory
WORKDIR /usr/src/app

# Set the entrypoint command
CMD ["/bin/bash"]
