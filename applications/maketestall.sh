#!/usr/bin/env bash
#
# run all tests of applications
#
# run with bash -e maketestall.sh to stop on errors
#

APPLICATIONS=(Coloring Scheduler STP)
LPSOLVERS=(spx2)
OPTS=(dbg)

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
	    LPILIB=../../lib/liblpi$LPS.$OSTYPE.$ARCH.gnu.$OPT.a
	    if test -e $LPILIB
            then
		SCIPLIB=../../lib/libscip.$OSTYPE.$ARCH.gnu.$OPT.a
		if test -e $SCIPLIB
		then
		    echo make OPT=$OPT LPS=$LPS
		    if (! make OPT=$OPT LPS=$LPS )
		    then
			exit $STATUS
		    fi
		    echo
		    echo make OPT=$OPT LPS=$LPS test
		    if (! make OPT=$OPT LPS=$LPS test )
	            then
			exit $STATUS
		    fi
		else
		    echo $SCIPLIB" does not exist - skipping combination ("$OPT", "$LPS")"
		fi
            else
		echo $LPILIB" does not exist - skipping combination ("$OPT", "$LPS")"
            fi
	    echo
	done
    done
    cd -
done
