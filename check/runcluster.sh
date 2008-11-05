#!/bin/bash
OUTFILE="/scratch/scip_output.out"
ERRFILE=$BASENAME.err
TMPFILE=$BASENAME.tmp

uname -a                            > $OUTFILE
uname -a                            >> $ERRFILE
echo @01 $FILENAME ===========      >> $OUTFILE 
echo @01 $FILENAME ===========      >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
date                                >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date +"@03 %s"                      >> $OUTFILE
$SCIPPATH/../$BINNAME < $TMPFILE    >> $OUTFILE 2>>$ERRFILE
date +"@04 %s"                      >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $ERRFILE
echo                                >> $OUTFILE
echo =ready=                        >> $OUTFILE
rm -f $TMPFILE
chmod g+r $OUTFILE
chmod g+r $ERRFILE
chmod g+r $BASENAME.set

mv $OUTFILE $BASENAME.out
