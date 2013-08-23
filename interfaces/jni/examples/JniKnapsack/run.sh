#!/bin/bash

CURPATH=`pwd` 

# set class path
export CLASSPATH="$CURPATH/classes:../../lib/scip.jar"

# set library path
export LD_LIBRARY_PATH="$CURPATH/../../lib"

# run JniKnapsack example
java JniKnapsack
