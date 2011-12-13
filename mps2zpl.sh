#!/bin/bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#@file    mps2zpl.sh
#@brief   converts mps file to zpl file
#@author  Thorsten Koch
#@author  Kati Wolter
#

FILENAME=$1

# script requires one argument, the (gezipped) mps file name
ARGS=1                   
if [ $# -ne "$ARGS" ]
then
    echo "ERROR: wrong argument list, call <mps2zpl.sh filename.mps(.gz)>"
    exit;
fi

# check if file exists
if test ! -e $FILENAME
then
    echo "ERROR: file <$FILENAME> does not exist"
    exit;
fi

# check if file is (gzipped) mps file
if [[ "$FILENAME" != *.mps ]] && [[ "$FILENAME" != *.mps.gz ]]
then
    echo "ERROR: file <$FILENAME> is not of type *.mps or *.mps.gz"
    exit;
fi

# unzip gzipped mps file and get name of mps file
if [[ "$FILENAME" == *.mps.gz ]]
then
    gunzip $FILENAME
    MPSNAME=$(echo "$FILENAME" | sed 's/\.gz//')
else
    MPSNAME=$FILENAME
fi

# convert mps file to zpl file
ZPLNAME=$(echo "$MPSNAME"   | sed 's/\.mps/\.zpl/')
awk -f mps2zpl.awk $MPSNAME > $ZPLNAME

## call zimpl to test zpl file (creates lp and tbl file in current directory)
#HARDTIMELIMIT=1800
#HARDMEMLIMIT=`expr 5000 \* 1024`
#bash -c " ulimit -s 100000; ulimit -t $HARDTIMELIMIT s; ulimit -v $HARDMEMLIMIT k; ulimit -f 200000; zimpl $ZPLNAME"

# zip mps file again
if [[ "$FILENAME" == *.mps.gz ]]
then
    gzip $MPSNAME
    #gzip $ZPLNAME
fi


