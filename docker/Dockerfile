#------------------------------------------------------------------------------------#
#
#   Created by Alem Gusinac, last modified at 15-10-2026
#
#   It builds on top of a slim ubuntu (70MB) base
#
#------------------------------------------------------------------------------------#

FROM ubuntu:22.04

# set environment without graphical interface
ENV DEBIAN_FRONTEND=noninteractive

#------------------------------------------------------------------------------------#
# 1. Installation of essential linux packages
#------------------------------------------------------------------------------------#

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    gnupg \
    gfortran \
    make \
    cmake \
    git \
    wget curl \
    pandoc \
    rustc \
    cargo \
    pkg-config \
    autoconf \
    automake \
    libnode-dev \
    libblas-dev \
    liblapack-dev \
    libbz2-dev \
    libuv1-dev \
    libx11-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    libfreetype6-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff-dev \
    libicu-dev \
    libfontconfig1-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

#------------------------------------------------------------------------------------#
# 2. R-base:latest setup
#------------------------------------------------------------------------------------#
# Setting up R-base:latest, includes public key setup and addition of r-base version to apt depository
RUN install -d -m 0755 /etc/apt/keyrings && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
      | gpg --dearmor -o /etc/apt/keyrings/cran-archive-keyring.gpg && \
    echo "deb [signed-by=/etc/apt/keyrings/cran-archive-keyring.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
      > /etc/apt/sources.list.d/cran-r.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends r-base && \
    rm -rf /var/lib/apt/lists/*

# Install ggplot2 3.5.2 explicitly before other packages
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')"

# Copy requirements
COPY ./docker/install2.r .
COPY ./docker/installBioc.r .

# Required package for install2.r
RUN R -e "install.packages('docopt', dependencies=TRUE)"

# Essential R dependencies
RUN Rscript install2.r --error --skipinstalled \
    pak \
    BiocManager \
    b64 \
    cli \
    DT \
    downloadthis \
    rmarkdown \
    optparse \
    && rm -rf /tmp/downloaded_packages

RUN Rscript installBioc.r --error --skipinstalled \
    rhdf5 \
    && rm -rf /tmp/downloaded_packages

# Install via github issue: https://github.com/agusinac/OmicFlow/issues/17
RUN R -e "pak::pak('agusinac/OmicFlow')"

# Make autoFlow.R directly callable from /usr/local/bin
RUN ln -s "$(Rscript -e 'cat(system.file("exec", "autoFlow.R", package = "OmicFlow"))')" /usr/local/bin/autoflow \
    && chmod +x /usr/local/bin/autoflow

#------------------------------------------------------------------------------------#
# 3. non-root user
#------------------------------------------------------------------------------------#

ARG USERNAME=docker
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME