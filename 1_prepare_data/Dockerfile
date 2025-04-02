FROM bioconductor/bioconductor_docker:RELEASE_3_20

#install required packages
RUN R -e 'BiocManager::install(c("tidyverse", "GEOquery", "janitor", "tools"), force = TRUE)'

WORKDIR /1_prepare_data

