#!/usr/bin/env bash

export LANG=C

# if called from another directory, change to the directory of this script
# this is needed to find the awk scripts used
cd $(dirname $0)

# first parameter is the directory into which the .res files shoule be written
TARGET_DIR=$1
# all other parameters are handled as problem files
FILES=${@:2}

RESFILES=""
REDRESFILES=""

VERSION=""

for FILE in $FILES
do
    echo "File '$FILE'"
    NEXT_VERSION=$(echo "$FILE" | rev | cut -d. -f10-11 | rev) # separated by '.', the version is on the 10th and 11th place from the end
    if [[ $VERSION = "" ]]
    then
        VERSION=$NEXT_VERSION
    elif [[ $NEXT_VERSION != $VERSION ]]
    then
        echo "ERROR: File '$FILE' has not the same version as the previous files. Expected: $VERSION, Was: $NEXT_VERSION"
        exit 1
    fi
    echo
done

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

./allcmpres.sh $RESFILES | tee $TARGET_DIR/$VERSION.normal.cmpres
./allcmpres.sh $REDRESFILES | tee $TARGET_DIR/$VERSION.woBranchTime.cmpres

# clean up afterwards, as those res files are not needed anymore
rm $RESFILES
rm $REDRESFILES
