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
    R -e "remotes::install_github('RGLab/RProtoBufLib@ce730aaf6613694aaff06df19dbee8448d379a73')" && \
    R -e "remotes::install_github('RGLab/cytolib@3d4045b15e86bdffb776b26e060fc5046180c4c5')" && \
    R -e "remotes::install_github('RGLab/flowCore@cd4446e9cf3c4e567829dc42773e323a18243b92')" && \
    R -e "remotes::install_github('RGLab/flowWorkspace@7a72eed6d90d67098be3676378c241a411267557')" && \
    R -e "remotes::install_github('RGLab/flowStats@f4e2acdd8a6bfab06d372529a6f79fd8ec6b3085')" && \
    R -e "remotes::install_github('RGLab/openCyto@5f7cf813369cf619f218231394d651042954e643')" && \
    R -e "remotes::install_github('DillonHammill/CytoExploreRData')" && \
    R -e "remotes::install_github('DillonHammill/CytoExploreR')"

