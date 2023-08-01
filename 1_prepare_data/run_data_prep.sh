#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/metaanalysis_01 .

#######################################################
# Run docker command
#######################################################

dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/1_prepare_data \
    -v $(pwd)/../../curated_data/Data:/Data \
    inwosu/metaanalysis_01"

time $dockerCommand Rscript scripts/data_prep.R

# $dockerCommand bash