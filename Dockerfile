# BIOCONDUCTOR ROCKER
FROM bioconductor/bioconductor_docker:latest

# DEPENDENCIES
RUN apt-get update && apt-get install --no-install-recommends -y \
    libopenblas-dev \
    autoconf \
    libltdl-dev \
    zlib1g-dev \
    wget \
    unzip \
    make \
    g++

# PYTHON
# RUN apt-get install -y python3 python3-dev python3-pip python3-venv

# FFTW
RUN wget http://www.fftw.org/fftw-3.3.8.tar.gz && \
    tar -zxvf fftw-3.3.8.tar.gz && \
    rm -rf fftw-3.3.8.tar.gz && \
    cd fftw-3.3.8 && \
    ./configure && \
    make && \
    make install

# FIT-SNE
RUN cd /. && \
    wget --no-check-certificate -O FIt-SNE.zip https://github.com/KlugerLab/FIt-SNE/archive/master.zip && \
    unzip FIt-SNE.zip && \
    rm -rf FIt-SNE.zip && \
    cd FIt-SNE-master && \
    g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm -Wno-address-of-packed-member && \
    cd /. && \
    mv /FIt-SNE-master home/rstudio/.FIt-SNE

# R DEPENDENCIES
RUN R -e "options('timeout' = 999999)" && \
    R -e "install.packages('BiocManager')" &&\
    R -e "BiocManager::install(c('cytolib', 'flowCore', 'flowWorkspace', 'openCyto', 'ggcyto', 'CytoML', 'flowWorkspaceData'))" &&\
    R -e "devtools::install_github('RGLab/cytoqc')" &&\
    R -e "devtools::install_github('DillonHammill/openCyto', force = TRUE)" && \
    R -e "devtools::install_github('DillonHammill/CytoExploreRData')" && \
    R -e "devtools::install_github('DillonHammill/DataEditR')" && \
    R -e "devtools::install_github('DillonHammill/HeatmapR')" && \
    R -e "devtools::install_github('DillonHammill/CytoExploreR', ref = 'refine')" && \
    R -e "install.packages(c('dplyr', 'tidyr', 'stringr', 'magrittr', 'forcats'))" && \
    R -e "install.packages(c('overlapping', 'FNN'))"

# RETICULATE
RUN R -e "install.packages('reticulate')"
    R -e "reticulate::py_discover_config()" && \
    R -e "reticulate::virtualenv_create('cytoexplorer')" && \
    R -e "reticulate::use_virtualenv('cytoexplorer')" && \
    R -e "reticulate::py_install(c('numpy', 'pacmap', 'openTSNE'), 'cytoexplorer', pip = TRUE)"
