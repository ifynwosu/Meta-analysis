#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/metaanalysis_cross_validation_07 .

#######################################################
# Run docker command
#######################################################

dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/7_cross_validation \
    -v $(pwd)/../Data:/Data \
    inwosu/metaanalysis_cross_validation_07"

time $dockerCommand python3 scripts/run_all_cross_val.py

# $dockerCommand bash
