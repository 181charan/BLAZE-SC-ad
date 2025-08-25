#!/bin/bash

# Get the current working directory
WORK_DIR=$(pwd)
echo "Info: Working directory: $WORK_DIR"

#Ensure that perf_runs directory exists. If not throw an error and exit.
if [ ! -d "perf_runs" ]; then
    echo "Error: perf_runs directory does not exist. Please run build.sh and compile BLAZE first, and run perf_tests.sh before attempting to generate graphs."
    exit 1
fi

python graphs.py $WORK_DIR/perf_runs