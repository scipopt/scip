#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      *
#*                                                                           *
#*  Licensed under the Apache License, Version 2.0 (the "License");          *
#*  you may not use this file except in compliance with the License.         *
#*  You may obtain a copy of the License at                                  *
#*                                                                           *
#*      http://www.apache.org/licenses/LICENSE-2.0                           *
#*                                                                           *
#*  Unless required by applicable law or agreed to in writing, software      *
#*  distributed under the License is distributed on an "AS IS" BASIS,        *
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
#*  See the License for the specific language governing permissions and      *
#*  limitations under the License.                                           *
#*                                                                           *
#*  You should have received a copy of the Apache-2.0 license                *
#*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

SCIPDIR=..
OUTFILE=`pwd`/check_scip_defs.out

cd $SCIPDIR
rm -f $OUTFILE

DEFS=`grep -h '#ifdef' src/scip/* | awk ' // { printf("%s\n", $2); }' | sort -u | awk '// { printf("%s ", $1); }'`

# IMPLINTSARECONT, BETTERWEIGHTFORDEMANDNODES, MAKECUTINTEGRAL, and SEPARATEROWS are defined by default
# ALLDIFFERENT and WITHEQKNAPSACK are excluded because they need non-default plugins
EXCLUDES=" __GNUC__ __cplusplus CPX_VERSION_VERSION ___DEBUG _MSC_VER NDEBUG ___NDEBUG NO_NEXTAFTER NO_RAND_R SCIP_NO_SIGACTION NO_STRERROR_R SCIP_NO_STRTOK_R SCIP_ROUNDING_FE SCIP_ROUNDING_FP SCIP_ROUNDING_MS SCIP_DEBUG_SOLUTION SORTTPL_BACKWARDS SORTTPL_EXPANDNAME SORTTPL_FIELD1TYPE SORTTPL_FIELD2TYPE SORTTPL_FIELD3TYPE SORTTPL_FIELD4TYPE SORTTPL_FIELD5TYPE SORTTPL_FIELD6TYPE SORTTPL_INDCOMP SORTTPL_NAME SORTTPL_PTRCOMP __sun _WIN32 SCIP_WITH_GMP SCIP_WITH_READLINE SCIP_WITH_ZIMPL SCIP_WITH_ZLIB QUADCONSUPGD_PRIORITY LINCONSUPGD_PRIORITY ALLDIFFERENT WITHEQKNAPSACK IMPLINTSARECONT BETTERWEIGHTFORDEMANDNODES MAKECUTINTEGRAL SEPARATEROWS SCIP_MORE_DEBUG "

echo "gathering defines . . ." |& tee -a $OUTFILE

for i in $DEFS
do
    if [[ ! $EXCLUDES =~ " $i " ]]; then
        INCLUDES="$INCLUDES$i "
    else
        echo "- excluding $i" |& tee -a $OUTFILE
    fi
done

for i in $INCLUDES
do
    USRDEFS="$USRDEFS-D$i "
    echo "- including $i" |& tee -a $OUTFILE
done
echo |& tee -a $OUTFILE

LPSS=(none spx cpx grb) # xprs msk clp qso)
OPTS=(dbg opt)

for i in ${LPSS[@]}
do
    for k in ${OPTS[@]}
    do
        echo "make LPS=$i OPT=$k USRCFLAGS=\"$USRDEFS\" clean"
        echo
        make LPS=$i OPT=$k USRCFLAGS="$USRDEFS" clean
        echo

        echo "make LPS=$i OPT=$k USRCFLAGS=\"$USRDEFS\""
        echo
        make LPS=$i OPT=$k USRCFLAGS="$USRDEFS"  || exit 1
        echo


        if [[ ! $i == "grb" ]]
        then
            echo "make LPS=$i OPT=$k USRCFLAGS=\"$USRDEFS\" CC=g++ CFLAGS=\"\" ZIMPL=false clean"
            echo
            make LPS=$i OPT=$k USRCFLAGS="$USRDEFS" CC=g++ CFLAGS="" ZIMPL=false clean
            echo

            echo "make LPS=$i OPT=$k USRCFLAGS=\"$USRDEFS\" CC=g++ CFLAGS=\"\"" ZIMPL=false
            echo
            make LPS=$i OPT=$k USRCFLAGS="$USRDEFS" CC=g++ CFLAGS="" ZIMPL=false || exit 1
            echo
        fi


        echo "make LPS=$i OPT=$k COMP=clang USRCFLAGS=\"$USRDEFS\" clean"
        echo
        make LPS=$i OPT=$k COMP=clang USRCFLAGS="$USRDEFS" clean
        echo

        echo "make LPS=$i OPT=$k COMP=clang USRCFLAGS=\"$USRDEFS\""
        echo
        make LPS=$i OPT=$k COMP=clang USRCFLAGS="$USRDEFS"  || exit 1
        echo


        if [[ ! $i == "grb" ]]
        then
            echo "make LPS=$i OPT=$k COMP=clang USRCFLAGS=\"$USRDEFS\" CC=\"clang++ -x c++\" CFLAGS=\"\" ZIMPL=false clean"
            echo
            make LPS=$i OPT=$k COMP=clang USRCFLAGS="$USRDEFS" CC="clang++ -x c++" CFLAGS="" ZIMPL=false clean
            echo

            echo "make LPS=$i OPT=$k COMP=clang USRCFLAGS=\"$USRDEFS\" CC=\"clang++ -x c++\" CFLAGS=\"\"" ZIMPL=false
            echo
            make LPS=$i OPT=$k COMP=clang USRCFLAGS="$USRDEFS" CC="clang++ -x c++" CFLAGS="" ZIMPL=false || exit 1
            echo
        fi

#if test "$?" = 0
#then
#    echo
#
#    echo "make LPS=$i OPT=$k USRCFLAGS=\"$USRDEFS\" test"
#    echo
#    make LPS=$i OPT=$k USRCFLAGS="$USRDEFS" test
#    echo
#fi
    done
done |& tee -a $OUTFILE
