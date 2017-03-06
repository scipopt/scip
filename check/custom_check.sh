#!/usr/bin/env bash

export LANG=C

FILES=$@
RESFILES=""
REDRESFILES=""

for FILE in $FILES
do
    DIR=`dirname $FILE`

    # strip all possible file extensions
    BASE=`basename $FILE .out`
    BASE=`basename $BASE .err`
    BASE=`basename $BASE .set`
    BASE=`basename $BASE .res`
    BASE=`basename $BASE .tex`
    BASE=`basename $BASE .pav`

    # define the paths to all necessary source files
    OUTFILE=$DIR/$BASE.out 
    ERRFILE=$DIR/$BASE.err
    TEXFILE=$DIR/$BASE.tex
    PAVFILE=$DIR/$BASE.pav
    
    # define the path to the target file
    RESFILE=$DIR/$BASE.1_normal.res
    REDRESFILE=$DIR/$BASE.2_woBranchTime.res

    awk -f check.awk -v "TEXFILE=$TEXFILE" -v "PAVFILE=$PAVFILE" -v "ERRFILE=$ERRFILE" $OUTFILE | tee $RESFILE
    RESFILES="$RESFILES $RESFILE"

    awk -f custom_check.awk -v "TEXFILE=$TEXFILE" -v "PAVFILE=$PAVFILE" -v "ERRFILE=$ERRFILE" $OUTFILE | tee $REDRESFILE
    REDRESFILES="$REDRESFILES $REDRESFILE"
done

echo $RESFILES
echo $REDRESFILES

./allcmpres.sh $RESFILES | tee $DIR/comparisonStandard
./allcmpres.sh $REDRESFILES | tee $DIR/comparisonWOBranchTime
