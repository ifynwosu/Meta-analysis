FROM bioconductor/bioconductor_docker:RELEASE_3_17

#install required packages
RUN R -e 'BiocManager::install(c("BiocManager","tidyverse", "GEOquery", "janitor, "tools"), force = TRUE)'

WORKDIR /1_prepare_data

