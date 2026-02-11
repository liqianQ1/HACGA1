FROM ubuntu:22.04

LABEL maintainer="HACGA1 Docker"
ENV LANG=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

# --------------------------------------------------
# Install system dependencies
# --------------------------------------------------
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    unzip \
    git \
    default-mysql-client \
    perl \
    default-jdk \
    bzip2 \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    ca-certificates \
    gnupg \
    dirmngr \
    sed \
    vim \
    && rm -rf /var/lib/apt/lists/*
    
# --------------------------------------------------
# Install Miniconda
# --------------------------------------------------
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda clean -ya

ENV PATH="/opt/conda/bin:${PATH}"


# --------------------------------------------------
# Conda environment
# --------------------------------------------------
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

RUN conda install -n base -c conda-forge mamba -y && \
    mamba env create -f /tmp/environment.yml

ENV PATH="/opt/conda/envs/ann/bin:${PATH}"


# --------------------------------------------------
# PASA global config
# --------------------------------------------------

ENV PASA_HOME=/opt/conda/envs/ann/opt/pasa-2.5.2

COPY annotation_gene /pipeline/annotation_gene

# --------------------------------------------------
# GeneMark-ES
# --------------------------------------------------
COPY gmes_linux_64_4 /software/gmes_linux_64_4
ENV GMES=/software/gmes_linux_64_4
ENV PATH=/software/gmes_linux_64_4:$PATH

# --------------------------------------------------
# Annotation pipeline
# --------------------------------------------------
COPY environment.yml /tmp/environment.yml

WORKDIR /pipeline

# --------------------------------------------------
# Entrypoint
# --------------------------------------------------
COPY docker-entrypoint.sh /docker-entrypoint.sh
RUN chmod +x /docker-entrypoint.sh

ENTRYPOINT ["/docker-entrypoint.sh"]
CMD []

