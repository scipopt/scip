#!/usr/bin/env bash
#
# run with bash -e makeallclean.sh to stop on errors
#

APPLICATIONS=(Coloring CycleClustering MinIISC PolySCIP Ringpacking Scheduler STP)
LPSOLVERS=(spx2 cpx none)
OPTS=(opt dbg)
SHARED=(true false)

# determine architecture
ARCH=`uname -m | \
    sed \
    -e 's/sun../sparc/' \
    -e 's/i.86/x86/' \
    -e 's/i86pc/x86/' \
    -e 's/[0-9]86/x86/' \
    -e 's/amd64/x86_64/' \
    -e 's/IP../mips/' \
    -e 's/9000..../hppa/' \
    -e 's/Power\ Macintosh/ppc/' \
    -e 's/00........../pwr4/'`
OSTYPE=`uname -s | tr '[:upper:]' '[:lower:]' | \
    sed \
    -e 's/cygwin.*/cygwin/' \
    -e 's/irix../irix/' \
    -e 's/windows.*/windows/' \
    -e 's/mingw.*/mingw/'`

for APPLICATION in ${APPLICATIONS[@]}
do
    echo
    echo
    echo ===== $APPLICATION =====
    echo
    cd $APPLICATION
    echo
    for OPT in ${OPTS[@]}
    do
	for LPS in ${LPSOLVERS[@]}
	do
	    for SHAREDVAL in ${SHARED[@]}
	    do
		echo make OPT=$OPT LPS=$LPS SHARED=$SHAREDVAL clean
		if (! make OPT=$OPT LPS=$LPS SHARED=$SHAREDVAL clean )
		then
		    exit $STATUS
		fi
		echo
	    done
	done
    done
    cd -
done
