#!/usr/bin/env bash
#
# run with bash -e makeall.sh to stop on errors
#

APPLICATIONS=(Coloring CycleClustering MinIISC PolySCIP Ringpacking Scheduler STP)
LPSOLVERS=(clp cpx none spx spx2)
OPTS=(opt dbg)
LIBTYPE=(static shared)

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
#   currently do not run depend to detect errors
#        echo make OPT=$OPT LPS=none ZIMPL=false depend
#	if (! make OPT=$OPT LPS=none ZIMPL=false depend)
#	then
#	    exit
#        fi
	for LPS in ${LPSOLVERS[@]}
	do
	    for TYPE in ${LIBTYPE[@]}
	    do
		if test "$TYPE" = "shared"
		then
		    SHAREDVAL="true"
		    LIBEXT="so"
		else
		    SHAREDVAL="false"
		    LIBEXT="a"
		fi

		SCIPLIB=../../lib/$TYPE/libscip.$OSTYPE.$ARCH.gnu.$OPT.$LPS.$LIBEXT
		if test -e $SCIPLIB
		then
		    echo make OPT=$OPT LPS=$LPS SHARED=$SHAREDVAL clean
		    if (! make OPT=$OPT LPS=$LPS SHARED=$SHAREDVAL clean )
		    then
			exit $STATUS
		    fi
		    echo
		    echo make OPT=$OPT LPS=$LPS SHARED=$SHAREDVAL
		    if (! make OPT=$OPT LPS=$LPS SHARED=$SHAREDVAL )
		    then
			exit $STATUS
		    fi
		else
		    echo $SCIPLIB" does not exist - skipping combination ("$OPT", "$LPS", "$TYPE")"
		fi
		echo
	    done
	done
    done
    cd -
done
