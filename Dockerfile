FROM ubuntu:22.04

LABEL maintainer="HACGA1 Docker"
ENV LANG=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

# --------------------------------------------------
# System dependencies
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
    gosu \
    && rm -rf /var/lib/apt/lists/*

# --------------------------------------------------
# Install Miniforge
# --------------------------------------------------
RUN wget -q https://mirrors.tuna.tsinghua.edu.cn/github-release/conda-forge/miniforge/LatestRelease/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

ENV PATH="/opt/conda/bin:${PATH}"

# --------------------------------------------------
# Create base non-root user (fallback user)
# --------------------------------------------------
RUN groupadd -g 1000 hacga && \
    useradd -m -u 1000 -g 1000 -s /bin/bash hacga

# --------------------------------------------------
# Install Conda Environment
# --------------------------------------------------
COPY environment.yml /tmp/environment.yml

RUN mamba env create -f /tmp/environment.yml -y && \
    conda clean -afy

RUN chmod -R 777 /opt/conda/envs/ann/opt

RUN chmod -R 777 /opt/conda/envs/ann/config

ENV PATH="/opt/conda/envs/ann/bin:${PATH}"

# --------------------------------------------------
# Software installation
# --------------------------------------------------
COPY gmes_linux_64 /software/gmes_linux_64_4
ENV GMES=/software/gmes_linux_64_4
ENV PATH=/software/gmes_linux_64_4:$PATH

COPY annotation_gene /pipeline/annotation_gene
WORKDIR /pipeline

RUN chmod -R 777 /pipeline /software

ENV PASA_HOME=/opt/conda/envs/ann/opt/pasa-2.5.2

# --------------------------------------------------
# Entrypoint
# --------------------------------------------------
COPY docker-entrypoint.sh /docker-entrypoint.sh
RUN chmod +x /docker-entrypoint.sh

# --------------------------------------------------
# Runtime as root (needed for dynamic UID mapping)
# --------------------------------------------------
USER root

ENTRYPOINT ["/docker-entrypoint.sh"]
CMD []
