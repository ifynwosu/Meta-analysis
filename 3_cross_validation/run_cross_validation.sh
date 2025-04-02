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

# prepare_cross_val_data runtime    29m1.144s

# cross_val runtimes
# ER_status  1473m58.590s
# PR_status  1321m41.232s
# HER2_status 1072m48.164s
 