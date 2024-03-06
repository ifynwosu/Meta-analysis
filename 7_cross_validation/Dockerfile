FROM python:3.12

#install required packages
RUN pip install scikit-learn==1.3.1
RUN pip install pandas==2.1.1

FROM bioconductor/bioconductor_docker:RELEASE_3_17

#install required packages
RUN R -e 'BiocManager::install(c("tidyverse", "tools", "janitor"), force = TRUE)'

WORKDIR /7_cross_validation
