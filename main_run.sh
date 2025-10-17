#!/bin/bash

set -eu
DATE=$(date +%F)
NXF_CONFIG=./nextflow.config
# Options: hyperqueue, docker
NXF_PROFILE='hyperqueue'
# Options:  rnaseq_count, prep_genome, or sra_download
NXF_ENTRY='sra_download'
# The output prefix on filenames for reports/logs
REPORT=${1:-"rnaseq_star_counts"}

# Set Debug > 0 to increase verbosity in nextflow logs
export NXF_DEBUG=2

# Ensure nextflow home variable is on a shared file system
# https://www.nextflow.io/docs/latest/reference/env-vars.html#nextflow-settings
export NXF_HOME="$HOME/nextflow"
mkdir -p $NXF_HOME

# Nextflow run to execute the workflow
PREFIX="${REPORT}_${DATE}"
nextflow -c ${NXF_CONFIG} \
    -log reports/${PREFIX}_nextflow.log \
    run main.nf \
    -e.NXF_HOME=$NXF_HOME \
    -cache TRUE \
    -entry ${NXF_ENTRY} \
    -profile ${NXF_PROFILE} \
    -with-report reports/${PREFIX}.html \
    -with-dag dag/${PREFIX}_dag.html \
    -with-trace reports/${PREFIX}_trace.txt 
