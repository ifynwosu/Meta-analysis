#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/metaanalysis_cross_validation_03 .

#######################################################
# Run docker command
#######################################################

dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/3_cross_validation \
    -v $(pwd)/../Data:/Data \
    inwosu/metaanalysis_cross_validation_03"

time $dockerCommand Rscript scripts/prepare_cross_val_data.R
time $dockerCommand python3 scripts/cross_val.py
time $dockerCommand Rscript scripts/plot_graphs.R

# $dockerCommand bash
