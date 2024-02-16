#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/metaanalysis_pr_status_04 .

#######################################################
# Run docker command
#######################################################

dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/4_PR_metaanalysis \
    -v $(pwd)/../1_prepare_data:/prepare_data \
    -v $(pwd)/../Data:/Data \
    inwosu/metaanalysis_pr_status_04"

time $dockerCommand Rscript scripts/run_PR_metaanalysis.R

# $dockerCommand bash
