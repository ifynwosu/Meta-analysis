#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/metaanalysis_race_02 .

#######################################################
# Run docker command
#######################################################

dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/2_race_metaanalysis \
    -v $(pwd)/../1_prepare_data:/prepare_data \
    -v $(pwd)/../Data:/Data \
    inwosu/metaanalysis_race_02"

time $dockerCommand Rscript scripts/run_race_metaanalysis.R

# $dockerCommand bash
