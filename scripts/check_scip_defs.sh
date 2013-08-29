#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

SCIPDIR=..
OUTFILE=`pwd`/check_scip_defs.out

cd $SCIPDIR
rm -f $OUTFILE

DEFS=`grep -h 'ifdef' src/scip/* | awk ' // { printf("%s\n", $2); }' | sort -u | awk '// { printf("%s ", $1); }'`
# echo "DEFS = \"$DEFS\""
# echo

# IMPLINTSARECONT, BETTERWEIGHTFORDEMANDNODES, MAKECUTINTEGRAL, and SEPARATEROWS are defined by default
# ALLDIFFERENT is excluded because it needs cons_alldifferent
EXCLUDES=" __cplusplus CPX_VERSION_VERSION ___DEBUG _MSC_VER NDEBUG ___NDEBUG NO_NEXTAFTER NO_RAND_R NO_SIGACTION NO_STRERROR_R NO_STRTOK_R ROUNDING_FE ROUNDING_FP ROUNDING_MS SCIP_DEBUG_SOLUTION SORTTPL_BACKWARDS SORTTPL_EXPANDNAME SORTTPL_FIELD1TYPE SORTTPL_FIELD2TYPE SORTTPL_FIELD3TYPE SORTTPL_FIELD4TYPE SORTTPL_FIELD5TYPE SORTTPL_FIELD6TYPE SORTTPL_INDCOMP SORTTPL_NAME SORTTPL_PTRCOMP __sun _WIN32 WITH_GMP WITH_READLINE WITH_ZIMPL WITH_ZLIB QUADCONSUPGD_PRIORITY LINCONSUPGD_PRIORITY ALLDIFFERENT IMPLINTSARECONT BETTERWEIGHTFORDEMANDNODES MAKECUTINTEGRAL SEPARATEROWS SCIP_MORE_DEBUG "
# echo "EXCLUDES = \"$EXCLUDES\""
# echo

echo
echo "gathering defines . . ."
echo
for i in $DEFS
do
    if [[ ! $EXCLUDES =~ " $i " ]]; then
	INCLUDES="$INCLUDES$i "
    else
	echo "- excluding $i"
    fi
done
echo

for i in $INCLUDES
do
    USRDEFS="$USRDEFS-D$i "
    echo "- including $i"
done
echo

LPSS=(grb cpx spx none) # spx132 xprs msk clp qso)
OPTS=(dbg opt prf opt-gccold)

# first ensure that links are available so we don't fail because of that when making with -j later
for i in ${LPSS[@]}
do
    for k in ${OPTS[@]}
    do
	echo "make LPS=$i OPT=$k USRCFLAGS=\"$USRDEFS\" links"
	echo
	make LPS=$i OPT=$k USRCFLAGS="$USRDEFS" links
	echo
    done
done

for i in ${LPSS[@]}
do
    for k in ${OPTS[@]}
    do
	echo "make LPS=$i OPT=$k USRCFLAGS=\"$USRDEFS\" clean"
	echo
	make LPS=$i OPT=$k USRCFLAGS="$USRDEFS" clean
	echo

	echo "make LPS=$i OPT=$k USRCFLAGS=\"$USRDEFS\"" -j
	echo
	make LPS=$i OPT=$k USRCFLAGS="$USRDEFS" -j || exit 1
	echo

	echo "make LPS=$i OPT=$k USRCFLAGS=\"$USRDEFS\" CC=g++ CFLAGS=\"\" ZIMPL=false clean"
	echo
	make LPS=$i OPT=$k USRCFLAGS="$USRDEFS" CC=g++ CFLAGS="" ZIMPL=false clean
	echo

	echo "make LPS=$i OPT=$k USRCFLAGS=\"$USRDEFS\" CC=g++ CFLAGS=\"\"" ZIMPL=false -j
	echo
	make LPS=$i OPT=$k USRCFLAGS="$USRDEFS" CC=g++ CFLAGS="" ZIMPL=false -j || exit 1
	echo

#	if test "$?" = 0
#	then
#	    echo
#
#	    echo "make LPS=$i OPT=$k USRCFLAGS=\"$USRDEFS\" test"
#	    echo
#	    make LPS=$i OPT=$k USRCFLAGS="$USRDEFS" test
#	    echo
#	fi
    done
done |& tee -a $OUTFILE
