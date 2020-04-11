# UBUNTU OS IMAGE
FROM ubuntu:16.04

# RSTUDIO 
FROM rocker/rstudio

# INSTALL DEPENDENCIES
RUN apt-get update && apt-get install --no-install-recommends -y\
    libprotobuf-dev \
    autoconf \
    automake \
    libtool \
    libxml2-dev \
    libhdf5-dev \
    zlib1g-dev \
    libpng-dev \
    tcl \
    wget \
    unzip \
    make \
    g++

# DOWNLOAD FFTW
RUN wget http://www.fftw.org/fftw-3.3.8.tar.gz && \
    tar -zxvf fftw-3.3.8.tar.gz && \
    rm -rf fftw-3.3.8.tar.gz && \
    cd fftw-3.3.8 && \
    ./configure && \
    make && \
    make install

# FIT-SNE - (path to fast_tsne -> FIt-SNE-master/fast_tsne.R)
RUN cd /. && \
    wget --no-check-certificate -O FIt-SNE.zip https://github.com/KlugerLab/FIt-SNE/archive/master.zip && \
    unzip FIt-SNE.zip && \
    rm -rf FIt-SNE.zip && \
    cd FIt-SNE-master && \
    g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm -Wno-address-of-packed-member && \
    cd /.

# INSTALL DEPENDENCIES
RUN R -e "install.packages('remotes')" && \
    R -e "install.packages('latticeExtra')" && \
    R -e "install.packages('XML')" && \
    R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install('flowCore')" && \
    R -e "BiocManager::install('CytoML')" && \
    R -e "BiocManager::install('flowWorkspace')" && \
    R -e "BiocManager::install('openCyto')" && \
    R -e "remotes::install_github('RGLab/RProtoBufLib')" && \
    R -e "remotes::install_github('RGLab/cytolib')" && \
    R -e "remotes::install_github('RGLab/flowCore')" && \
    R -e "remotes::install_github('RGLab/flowWorkspace')" && \
    R -e "remotes::install_github('RGLab/flowStats')" && \
    R -e "remotes::install_github('RGLab/openCyto')" && \
    R -e "remotes::install_github('DillonHammill/CytoExploreRData')" && \
    R -e "remotes::install_github('DillonHammill/CytoExploreR')"

