#!/bin/bash

# This script will work only if there is one .eval file in results/
# check/results should only contain results from one run

./evalcheck_cluster.sh results/*.eval
rbcli up results/*.{err,out,set}
