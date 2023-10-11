FROM bioconductor/bioconductor_docker:RELEASE_3_17

#install required packages
RUN R -e 'BiocManager::install(c("tidyverse", "tools", "janitor", "DExMA"), force = TRUE)'

WORKDIR /4_PR_metaanalysis
