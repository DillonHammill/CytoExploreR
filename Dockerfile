FROM ubuntu:16.04

# R and RStudio
FROM rocker/verse:devel

RUN apt-get update
RUN apt-get install libprotobuf-dev -y

# install dependencies
RUN apt-get update && apt-get install -y\
    autoconf \
    automake \
    libtool \
    libxml2 \
    libhdf5-dev

# install packages
RUN R -e "install.packages('BiocManager')" 
RUN R -e "BiocManager::install(version = 'devel')" 
RUN R -e "BiocManager::install('flowCore')" 
RUN R -e "BiocManager::install('flowWorkspace')" 
RUN R -e "BiocManager::install('openCyto')" 
RUN R -e "devtools::install_github('DillonHammill/CytoExploreRData')" 
RUN R -e "devtools::install_github('DillonHammill/CytoExploreR')"