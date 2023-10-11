#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/metaanalysis_cross_validation_06 .

#######################################################
# Run docker command
#######################################################

dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/6_cross_validation \
    -v $(pwd)/../../Meta_Analysis/Data:/Data \
    inwosu/metaanalysis_cross_validation_06"

time $dockerCommand python3 scripts/run_all_cross_val.py

# $dockerCommand bash
