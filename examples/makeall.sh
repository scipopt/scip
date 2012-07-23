#!/usr/bin/env bash
#
# run with bash -e makeall.sh to stop on errors
#

EXAMPLES=(Coloring Binpacking Eventhdlr LOP MIPSolver Queens SamplePricer SamplePricer_C Scheduler Testslack TSP VRP CallableLibrary)
LPSOLVERS=(clp cpx none spx)
OPTS=(opt dbg)

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

for EXAMPLE in ${EXAMPLES[@]}
do
    echo
    echo
    echo ===== $EXAMPLE =====
    echo
    cd $EXAMPLE
    echo
    for OPT in ${OPTS[@]}
    do
#   currently do not run depend to detect errors
#        echo make OPT=$OPT LPS=none ZIMPL=false depend
#	if (! make OPT=$OPT LPS=none ZIMPL=false depend)
#	then
#	    exit
#        fi
	for LPS in ${LPSOLVERS[@]}
	do
	    LPILIB=../../lib/liblpi$LPS.$OSTYPE.$ARCH.gnu.$OPT.a
	    if test -e $LPILIB
            then
		SCIPLIB=../../lib/libscip.$OSTYPE.$ARCH.gnu.$OPT.a
		if test -e $SCIPLIB
		then
		    echo make OPT=$OPT LPS=$LPS clean
		    if (! make OPT=$OPT LPS=$LPS clean )
		    then
			exit $STATUS
		    fi
		    echo
		    echo make OPT=$OPT LPS=$LPS
		    if (! make OPT=$OPT LPS=$LPS )
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
