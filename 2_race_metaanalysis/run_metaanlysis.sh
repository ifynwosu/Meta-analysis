#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/metaanalysis_02 .

#######################################################
# Run docker command
#######################################################

dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/2_race_metaanalysis \
    -v $(pwd)/../../curated_data/Data:/Data \
    inwosu/metaanalysis_02"

time $dockerCommand Rscript scripts/run_race_metaanalysis.R

# $dockerCommand bash