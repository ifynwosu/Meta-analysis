#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/metaanalysis_tri_neg_status_06 .

#######################################################
# Run docker command
#######################################################

dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/6_triple_negative \
    -v $(pwd)/../1_prepare_data:/prepare_data \
    -v $(pwd)/../Data:/Data \
    inwosu/metaanalysis_trip_neg_status_06"

time $dockerCommand Rscript scripts/run_tri_neg_metaanalysis.R

# $dockerCommand bash
