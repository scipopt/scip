#!/bin/bash

# This script will work only if there is one .eval file in results/
# check/results should only contain results from one run

cd check/
if [ -z "$TESTSET" ]; then
    echo "Missing testset information. Aborting..."
    exit 1
fi
./evalcheck_cluster.sh -R results/check.$TESTSET.*.eval
