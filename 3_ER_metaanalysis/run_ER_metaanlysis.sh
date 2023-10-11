#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/metaanalysis_er_status_03 .

#######################################################
# Run docker command
#######################################################

dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/3_ER_metaanalysis \
    -v $(pwd)/../../Meta_Analysis/1_prepare_data:/prepare_data \
    -v $(pwd)/../../Meta_Analysis/Data:/Data \
    inwosu/metaanalysis_er_status_03"

time $dockerCommand Rscript scripts/run_ER_metaanalysis.R

# $dockerCommand bash
