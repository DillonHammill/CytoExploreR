# UBUNTU OS IMAGE
FROM ubuntu:16.04

# REQUIRE LIBPROTOBUF 2.6.0
RUN apt-get update
RUN apt-get install libprotobuf-dev -y

# INSTALL DEPENDENCIES
RUN apt-get update && apt-get install -y\
    autoconf \
    automake \
    libtool \
    libxml2 \
    libhdf5-dev \
    wget \
    unzip \
    make \
    g++

# FIT-SNE ZIP
RUN wget -O FIt-SNE.zip https://github.com/KlugerLab/FIt-SNE/archive/master.zip

# UNZIP & REMOVE FIt-SNE ZIP
RUN unzip FIt-SNE.zip
RUN rm -rf FIt-SNE.zip

# DOWNLOAD FFTW
RUN wget http://www.fftw.org/fftw-3.3.8.tar.gz

# EXTRACT & REMOVE TARBALL
RUN tar -zxvf fftw-3.3.8.tar.gz
RUN rm -rf fftw-3.3.8.tar.gz

# CHANGE DIRECTORY
RUN cd fftw-3.3.8

# BUILD FFTW LIBRARY
RUN ./configure
RUN make
RUN make install

# FIT-SNE - (path to fast_tsne -> FIt-SNE-master/fast_tsne.R)
RUN cd /.
RUN cd FIt-SNE-master
RUN g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm -Wno-address-of-packed-member
RUN cd /.

# R 
FROM rocker/verse:devel

# INSTALL X11 DEPENDENCIES
#RUN apt-get update && apt-get install -y --no-install-recommends \
#    libx11-6 \
#    libxss1 \apt 
#    libxt6 \
#    libxext6 \
#    libsm6 \
#    libice6 \
#    xdg-utils \
#  && rm -rf /var/lib/apt/lists/*

# INSTALL REQUIRED PACKAGES
RUN R -e "BiocManager::install(version = 'devel')" 
RUN R -e "BiocManager::install('flowCore')" 
RUN R -e "BiocManager::install('flowWorkspace')" 
RUN R -e "BiocManager::install('openCyto')" 
RUN R -e "devtools::install_github('DillonHammill/CytoExploreRData')" 
RUN R -e "devtools::install_github('DillonHammill/CytoExploreR')"