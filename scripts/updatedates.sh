#!/usr/bin/env bash
#
# This bash script updates all copyrights in files that are under version
# control (git) and have a ZIB copyright.
#
# You just have to run this script. There is nothing to adjust.
# The correct year is detected through the 'date' function
#
# Note that not all files (usually scripts) contain a copyright. A copyright is only
# needed for those files which are part of a SCIP distribution (see makedist.sh)

set -euo pipefail

NEWYEAR=`date +"%Y"`

for f in `git ls-files` ; do

# TODO list files which may miss a copyright or have a copyright in wrong format

# skip files that are not text or do not have copyright string
grep -Iq "Copyright (C) .* Konrad-Zuse-Zentrum" $f || continue

#echo $f
# fixup copyright date
sed -i "s/Copyright (C) \([0-9]*\)-[0-9]* Konrad-Zuse-Zentrum/Copyright (C) \1-$NEWYEAR Konrad-Zuse-Zentrum/g" $f

done
