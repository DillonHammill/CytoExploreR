FROM ubuntu:16.04

# liprotobuf version 2.6.0
RUN apt-get update
RUN apt-get install libprotobuf-dev -y

# install dependencies
RUN apt-get update && apt-get install -y\
    autoconf \
    automake \
    libtool \
    libxml2 \
    libhdf5-dev

# R 
FROM rocker/r-devel

# Install X11 dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libx11-6 \
    libxss1 \
    libxt6 \
    libxext6 \
    libsm6 \
    libice6 \
    xdg-utils \
  && rm -rf /var/lib/apt/lists/*

# RStudio
FROM rocker/rstudio

# Tidyverse & devtools (for speed)
FROM rocker/tidyverse

# install packages
RUN R -e "BiocManager::install(version = 'devel')" 
RUN R -e "BiocManager::install('flowCore')" 
RUN R -e "BiocManager::install('flowWorkspace')" 
RUN R -e "BiocManager::install('openCyto')" 
RUN R -e "devtools::install_github('DillonHammill/CytoExploreRData')" 
RUN R -e "devtools::install_github('DillonHammill/CytoExploreR')"