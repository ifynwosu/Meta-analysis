FROM bioconductor/bioconductor_docker:RELEASE_3_17

#install required packages
RUN R -e 'BiocManager::install(c("BiocManager", "GEOquery",  "SCAN.UPC", "RankProd"), force = TRUE)'
RUN R -e 'install.packages(c("devtools", "tidyverse", "metafor", "rjson", "tools"), force = TRUE)'

WORKDIR /1_run_metaanalysis


