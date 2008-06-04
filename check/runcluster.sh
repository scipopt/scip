#!/bin/bash
uname -a
uname -a >&2
echo @01 $FILENAME ===========
echo @01 $FILENAME =========== >&2
echo -----------------------------
date
date >&2
echo -----------------------------
date +"@03 %s"
$SCIPPATH/../$BINNAME < $TMPFILE
date +"@04 %s"
echo -----------------------------
date
echo -----------------------------
date >&2
echo
echo =ready=
rm -f $TMPFILE
chmod g+r $OUTFILE
chmod g+r $ERRFILE
chmod g+r $SETFILE
