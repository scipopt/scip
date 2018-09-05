#!/usr/bin/env bash
#
# This scripts generates/updates the src/depend.* files for SCIP
#

make depend

pushd applications/Scheduler > /dev/null
make depend
popd > /dev/null
