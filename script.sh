#!/bin/bash
#
# Script for generating simulation data.
# Author: Nathan Malamud
#

# SCRIPT CONFIG - - - //

INPUT_FILE="cariaco.txt"

NUM_SETS=2
RUNS_PER_SIM=30

# Never do too many threads -- this causes issues in data output.
NUM_THREADS=3

OUTPUT_DIR=data

# - - - - - - - - - - //

# WARNING: RUNNING THIS SCRIPT WILL OVERWRITE THE EXISTING DATA DIRECTORY
rm -rf $OUTPUT_DIR 

COUNTER=0

for SET in $(seq 1 $NUM_SETS); do

    # Run 'NUM_THREADS' simulations in parallel
    for i in $(seq 1 $NUM_THREADS); do

        ((COUNTER++))
        printf -v LABEL "%02d" $COUNTER

        # Execute simulation
        python3 main.py \
            --runs=$RUNS_PER_SIM \
            --out="$OUTPUT_DIR/batch_$LABEL" \
            --debug \
            < $INPUT_FILE \
            & \
    done

    wait # until current set is complete before moving on
done
