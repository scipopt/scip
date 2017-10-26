#!/bin/sh

# moves all $USER files :

# for SoPlex
mv -n /usr/local/tmp/$USER.* $TARGETPATH/

# for SCIP
mv -n /usr/local/tmp/$USER-tmpdir/$USER.* $TARGETPATH/
