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
    -v $(pwd):/1_run_metaanalysis \
    -v $(pwd)/../../curated_data/Data:/Data \
    inwosu/metaanalysis_01"

time $dockerCommand Rscript scripts/run_metaanalysis.R

# $dockerCommand bash