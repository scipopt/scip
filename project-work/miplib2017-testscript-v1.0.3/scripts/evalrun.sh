#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

AWKARGS=""
FILES=""

# construct paths
MIPLIBPATH=`pwd`

# Additional parameter
LINTOL=1e-4   # absolut tolerance for checking linear constraints and objective value
INTTOL=1e-4   # absolut tolerance for checking integrality constraints

. scripts/eval_functions.sh

for i in $@
do
    echo $i
    if test ! -e $i
    then
        AWKARGS="$AWKARGS $i"
    else
        FILES="$FILES $i"
    fi
done

for i in $FILES
do
    NAME=$(basename $i .out)
    DIR=$(dirname $i)
    OUTFILE=$DIR/$NAME.out
    METAFILE=$DIR/$NAME.meta
    RESFILE=$DIR/$NAME.res

    if [[ -e ${METAFILE} ]]
    then
        echo "Create results for '${NAME}'"
        TSTNAME=$(readMetaData ${METAFILE} "TestName")
        SOLVER=$(readMetaData ${METAFILE} "Solver")

        SOLUFILE=$MIPLIBPATH/testsets/$TSTNAME.solu

         # check if a solution  file/link exists
        if test ! -e $SOLUFILE
        then
            SOLUFILE=$MIPLIBPATH/testsets/miplib2017.solu

            if test ! -e $SOLUFILE
            then
                echo "Warning: solution file/link <$TSTNAME.solu> or <miplib2017.solu> does not exist in <testsets> folder; therefore, no consistency check"
                SOLUFILE=""
            fi
        fi

        awk -f scripts/parse.awk -f scripts/parse_$SOLVER.awk -v "LINTOL=$LINTOL" $SOLUFILE $OUTFILE | tee $RESFILE
    else
        echo "Cannot parse the out files without the meta file: '${METAFILE}'."
    fi
done
