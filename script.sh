#!/bin/bash
#
# Example of how to use this code to generate simulation data.
# Author: Nathan Malamud
#

# SCRIPT CONFIG - - - //

MIN_SIZE=5 
MAX_SIZE=15

SETS_PER_SIZE=10
RUNS_PER_SIM=50

NUM_THREADS=12

OUTPUT_DIR=data

# - - - - - - - - - - //

# WARNING: RUNNING THIS SCRIPT WILL OVERWRITE THE EXISTING DATA DIRECTORY
rm -rf $OUTPUT_DIR 

COUNTER=0

for SIZE in $(seq $MIN_SIZE $MAX_SIZE); do

    for SET in $(seq 1 $SETS_PER_SIZE); do

        # Run 'NUM_THREADS' simulations in parallel
        for i in $(seq 1 $NUM_THREADS); do

                ((COUNTER++))
                printf -v LABEL "%02d" $COUNTER

                python3 main.py -k=$SIZE --runs=$RUNS_PER_SIM --out="$OUTPUT_DIR/batch_$LABEL" &
        done

        wait # until current set is complete before moving on
    done

done
