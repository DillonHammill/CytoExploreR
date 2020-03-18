FROM ubuntu:16.04

FROM rocker/verse:devel

RUN apt-get update
RUN apt-get install libprotobuf-dev -y

RUN apt-get update && apt-get install -y\
    autoconf \
    automake \
    libtool \
    libxml2 \
    libhdf5-dev

# Install dependencies 
RUN R -e "library(devtools)"
RUN R -e "install.packages('BiocManager')"
RUN R -e "library(BiocManager)"
RUN R -e "install(version = 'devel')"
RUN R -e "install('flowCore')"
RUN R -e "install('flowWorkspace')"
RUN R -e "install('openCyto')"

# Install CytoExploreRData and CytoExploreR
RUN R -e "install_github('DillonHammill/CytoExploreRData')"
RUN R -e "install_github('DillonHammill/CytoExploreR')"