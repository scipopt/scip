#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

SOLVER=$1
BINNAME=$2
NAME=$3
TIMELIMIT=$4
MEMLIMIT=$5   # The memory limit (in MB)
SOLFILE=$6
THREADS=$7
MIPGAP=$8

# import some useful functions that are reused by other scripts
. $(dirname "${BASH_SOURCE[0]}")/run_functions.sh

TMPFILE=check.$SOLVER.tmp

echo > $TMPFILE
echo > $SOLFILE

# disable log file
echo set logfile '*'                   >> $TMPFILE

# set threads to given value
if test $THREADS != 0
then
    echo set threads $THREADS          >> $TMPFILE
fi

# set mipgap to given value
echo set mip tolerances mipgap $MIPGAP >> $TMPFILE

# set timing to wall-clock time and pass time limit
echo set clocktype 2                   >> $TMPFILE
echo set timelimit $TIMELIMIT          >> $TMPFILE

# set the available memory
echo set mip limits treemem $MEMLIMIT  >> $TMPFILE

# use deterministic mode (warning if not possible)
echo set parallel 1                    >> $TMPFILE

# set display options
# nothing for CPLEX

# read, optimize, display statistics, write solution, and exit
echo read $NAME                        >> $TMPFILE
echo change sense 0                    >> $TMPFILE  # to identify obj sense
echo                                   >> $TMPFILE
echo optimize                          >> $TMPFILE
echo write $SOLFILE.sol                >> $TMPFILE
echo set logfile $SOLFILE              >> $TMPFILE  # create file to mark successful run
echo display solution objective        >> $TMPFILE
echo set logfile '*'                   >> $TMPFILE
echo quit                              >> $TMPFILE

rm -f $SOLFILE.sol
$BINNAME < $TMPFILE

if test -f $SOLFILE
then
    # translate CPLEX solution format into format for solution checker.
    #  The SOLFILE format is a very simple format where in each line
    #  we have a <variable, value> pair, separated by spaces.
    #  A variable name of =obj= is used to store the objective value
    #  of the solution, as computed by the solver. A variable name of
    #  =infeas= can be used to indicate that an instance is infeasible.
    if test -f $SOLFILE.sol
    then
        $(dirname "${BASH_SOURCE[0]}")/cpx2solchecker.py $SOLFILE.sol > $SOLFILE
    else
	    echo "=infeas=" > $SOLFILE
    fi
fi
rm -f $SOLFILE.sol

# Replace remaining HTML entities
unescapeHTMLentities $SOLFILE

# remove tmp file
rm $TMPFILE

# remove CPLEX log
rm cplex.log
