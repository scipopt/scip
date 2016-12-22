#!/bin/bash

# in the jenkins build instructions, the TESTSET variable should be set
cd check/
if [ -z "$TESTSET" ]; then
    echo "Missing testset information. Aborting..."
    exit 1
fi
./evalcheck_cluster.sh -R results/check.$TESTSET.*.eval
