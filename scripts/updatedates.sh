#!/usr/bin/env bash
#
# This bash script updates all copyrights in files that are under version
# control (git) and have a ZIB copyright.
#
# You just have to run this script. There is nothing to adjust.
# The correct year is detected through the 'date' function
#
# Note that not all files (usually scripts) contain a copyright. A copyright is only
# needed for those files which are part of a SCIP distribution (see scripts/makedist.sh)
#
# USAGE: ./scripts/updatedates.sh

set -euo pipefail

NEWYEAR=`date +"%Y"`

echo "Updating copyright of all files under version control and list findings of possibly incorrect copyright string."

for f in `git ls-files` ; do

    # skip binary files
    grep -Iq . $f || continue

    # skip this file
    [[ $f =~ "updatedates.sh" ]] && continue

    # skip symbolic links
    if `test -L $f`; then
        echo "Skipping symbolic link $f"
        continue
    fi

    # process files with ZIB copyright string that do not include current year
    if grep -o 'Copyright (c) [0-9]*-[0-9]* Zuse Institute Berlin (ZIB)' $f | grep -vq $NEWYEAR ; then
        echo "Updating $f"
        sed -i "s/Copyright (c) \([0-9]*\)-[0-9]* Zuse Institute Berlin (ZIB)/Copyright (c) \1-$NEWYEAR Zuse Institute Berlin (ZIB)/g" $f
    fi

    # print matches for lines that have "Copyright" and "Zuse" but are not a valid ZIB copyright
    grep -iH "Copyright.*Zuse" $f | grep -v "Copyright (c) [0-9]*-$NEWYEAR Zuse Institute Berlin (ZIB)" || true

done

sed -i "s/\([0-9]*\)-[0-9]* by Zuse Institute Berlin (ZIB)/\1-$NEWYEAR by Zuse Institute Berlin (ZIB)/" doc/scipfooter.html
