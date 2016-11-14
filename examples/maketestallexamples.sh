#!/usr/bin/env bash
#
# Run all tests of examples, arguments are passed to the make test command.
# Parameter "-q" turns on the quiet mode, i.e., does not output the logging of the programs.
#
#

EXAMPLES=(Binpacking CallableLibrary Eventhdlr GMI LOP MIPSolver Queens TSP VRP)
LPSOLVERS=(spx2)
OPTS=(dbg)

echo "Running all tests on examples."

# parse command line
MAKEARGS=""
QUIET=0
for i in $@
do
    if test "$i" = "-q"
    then
	echo "Quiet mode."
	QUIET=1
    else
	MAKEARGS="$MAKEARGS $i"
    fi
done


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

# prepare log file
echo "" > exampletestsummary.log

# pretest
for OPT in ${OPTS[@]}
do
    for LPS in ${LPSOLVERS[@]}
    do
	LPILIB=../lib/liblpi$LPS.$OSTYPE.$ARCH.gnu.$OPT.a
	if test ! -e $LPILIB
	then
	    echo "Error: "$LPILIB" does not exist, please compile SCIP with OPT="$OPT" and LPS="$LPS"." >> ../applicationtestsummary.log
	    echo "Error: "$LPILIB" does not exist, please compile SCIP with OPT="$OPT" and LPS="$LPS"."
	    exit 1
	fi
	SCIPLIB=../lib/libscip.$OSTYPE.$ARCH.gnu.$OPT.a
	if test ! -e $SCIPLIB
	then
	    echo "Error: "$SCIPLIB" does not exist, please compile SCIP with OPT="$OPT" and LPS="$LPS"." >> ../applicationtestsummary.log
	    echo "Error: "$SCIPLIB" does not exist, please compile SCIP with OPT="$OPT" and LPS="$LPS"."
	    exit 1
	fi
    done
done

# run tests
for EXAMPLE in ${EXAMPLES[@]}
do
    echo
    echo
    echo ===== $EXAMPLE =====
    echo ===== $EXAMPLE ===== >> exampletestsummary.log
    echo
    cd $EXAMPLE
    for OPT in ${OPTS[@]}
    do
	for LPS in ${LPSOLVERS[@]}
	do
	    echo make OPT=$OPT LPS=$LPS $MAKEARGS
	    if (! make OPT=$OPT LPS=$LPS $MAKEARGS )
	    then
		echo "Making "$EXAMPLE" failed." >> ../exampletestsummary.log
		exit $STATUS
	    else
		echo "Making "$EXAMPLE" successful." >> ../exampletestsummary.log
	    fi
	    echo
	    if test $QUIET = 1
	    then
		echo make OPT=$OPT LPS=$LPS $MAKEARGS test
		if ( ! make OPT=$OPT LPS=$LPS $MAKEARGS test > /dev/null )
		then
		    echo "Testing "$EXAMPLE" failed."
		    echo "Testing "$EXAMPLE" failed." >> ../exampletestsummary.log
		    exit $STATUS
		fi
	    else
		echo make OPT=$OPT LPS=$LPS $MAKEARGS test
		if ( ! make OPT=$OPT LPS=$LPS $MAKEARGS test )
		then
		    echo "Testing "$EXAMPLE" failed."
		    echo "Testing "$EXAMPLE" failed." >> ../exampletestsummary.log
		    exit $STATUS
		fi
	    fi
	    echo "Testing "$EXAMPLE" successful."
	    echo "Testing "$EXAMPLE" successful." >> ../exampletestsummary.log

	    # find most recently changed result file and display it
	    if test -d check/results
	    then
		RESFILE=`find check/results/*.res -type f -printf '%T@ %p\n' | sort -n | tail -1 | cut -f2- -d" "`
		if test -e $RESFILE
		then
		    cat $RESFILE >> ../exampletestsummary.log
		fi
	    fi
	    echo
	    echo >> ../exampletestsummary.log
	done
    done
    cd - > /dev/null
done

echo
echo
echo ===== Summary =====

cat exampletestsummary.log
rm -f exampletestsummary.log
