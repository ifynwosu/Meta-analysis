FROM python:3.13

#install required packages
RUN pip install scikit-learn==1.6.1
RUN pip install pandas==2.2.3

FROM bioconductor/bioconductor_docker:RELEASE_3_20

#install required packages
RUN R -e 'BiocManager::install(c("tidyverse", "tools", "janitor"), force = TRUE)'

WORKDIR /3_cross_validation
