#------------------------------------------------------------------------------------#
#
#   Created by Alem Gusinac, last modified at 24-10-2025
#
#   It builds on top of a slim ubuntu (70MB) base
#
#------------------------------------------------------------------------------------#

FROM ubuntu:20.04

# set environment without graphical interface
ENV DEBIAN_FRONTEND=noninteractive

#------------------------------------------------------------------------------------#
# 1. Installation of essential linux packages
#------------------------------------------------------------------------------------#

# Install dependencies
RUN apt-get update && apt-get install -y \
    make \
    cmake \
    git \
    wget curl \
    pandoc \
    rustc \
    cargo \
    pgp \
    pkg-config \
    autoconf automake \
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
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" >> /etc/apt/sources.list
# apt update, upgrade and clean after above set-up is required!
# Otherwise it gives an error of broken pipes.
RUN apt-get update && apt-get upgrade -y && apt-get clean
RUN apt-get update && apt-get install -y r-base

# Install ggplot2 3.5.2 explicitly before other packages
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('ggplot2', version='3.5.2', repos='https://cloud.r-project.org')"

# Copy requirements
COPY install2.r .
COPY installBioc.r .

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