#!/bin/bash

echo "running snakemake with the following inputs:"
echo "  snakefile:          $1"
echo "  config file:        $2"
echo "  project directory:  $3"
echo "  rule to run:        $4"
echo "  extra arguments:    $5"
echo

snakefile=$1
configfile=$2
proj_dir=$3
rule=$4
extraargs=$5

# get directory of this file, change to it
SCPROCESS_DIR="$(dirname "$(readlink -f "$0")")"
cd $SCPROCESS_DIR

echo "starting snakemake rule" $rule

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
  $rule

 snakemake ${SM_ARGS} \
  --snakefile $snakefile \
  --configfile $configfile \
  --rerun-incomplete \
  --verbose \
  --use-conda \
  $rule
