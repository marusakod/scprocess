#!/bin/bash

echo "running snakemake with the following inputs:"
echo "  snakefile:          $1"
echo "  config file:        $2"
echo "  extra arguments:    $3"
echo

snakefile=$1
configfile=$2
extraargs=$3

# get directory of this file, change to it
SCPROCESS_DIR="$(dirname "$(readlink -f "$0")")"
cd $SCPROCESS_DIR

echo "starting scprocess setup"

SM_ARGS=" "

if [ "$extraargs" != 0 ]
  then
   SM_ARGS+="$extraargs"
fi

# run snakemake

echo snakemake ${SM_ARGS} \
  --snakefile $snakefile \
  --configfile $configfile \
  --rerun-incomplete \
  --verbose \
  --use-conda \

 snakemake ${SM_ARGS} \
  --snakefile $snakefile \
  --configfile $configfile \
  --rerun-incomplete \
  --verbose \
  --use-conda \

