#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

SOLVER=$1
BINNAME=$2
NAME=$3
TIMELIMIT=$4
MEMLIMIT=$5   # The memory limit (in MB) # this is used for each SCIP environment
SOLFILE=$6
THREADS=$7
MIPGAP=$8
TIMEOUT=`expr $4 + 30`

FIBERSETTINGS=check.$SOLVER.$THREADS.tmp
SCIPSETTINGS=check.$SOLVER.$THREADS.set

if test $THREADS != 0
then
    MEMLIMIT=`echo "$MEMLIMIT/$THREADS" | bc -l` # If you want to use MEMLIMIT for while solver threads except LC
fi

# add one thread for the load coordinator

echo > $FIBERSETTINGS
echo > $SCIPSETTINGS
echo > $SOLFILE

# set threads to given value
# nothing to be done for SCIP

# set mipgap to given value
echo limits/absgap = $MIPGAP           >> $SCIPSETTINGS

# set the available memory, but this is used for each SCIP environment
echo limits/memory = $MEMLIMIT        >> $SCIPSETTINGS

# set timing to wall-clock time and pass time limit
echo TimeLimit = $TIMELIMIT           >> $FIBERSETTINGS

# use deterministic mode (warning if not possible)
echo "WARNING: fscip can not run in deterministic mode"
# nothing to be done for SCIP

# run fiber SCIP : TIMEOUT should be more than TIMELIMIT
echo timeout $TIMELIMIT $BINNAME $FIBERSETTINGS $NAME -sr $SCIPSETTINGS -s $SCIPSETTINGS -fsol $SOLFILE.sol -sth $THREADS SCIP_MEMLIMIT = $MEMLIMIT
if test $THREADS != 0
then
    timeout $TIMEOUT $BINNAME $FIBERSETTINGS $NAME -sr $SCIPSETTINGS -s $SCIPSETTINGS -fsol $SOLFILE.sol -sth $THREADS
else
    timeout $TIMEOUT $BINNAME $FIBERSETTINGS $NAME -sr $SCIPSETTINGS -s $SCIPSETTINGS -fsol $SOLFILE.sol
fi

if test -f $SOLFILE.sol
then
   # translate SCIP solution format into format for solution checker. The
   # SOLFILE format is a very simple format where in each line we have a
   # <variable, value> pair, separated by spaces.  A variable name of
   # =obj= is used to store the objective value of the solution, as
   # computed by the solver. A variable name of =infeas= can be used to
   # indicate that an instance is infeasible.
    sed ' /solution status:/d;
          /\[ Final Solution \]/d;
         s/objective value:/=obj=/g;
         s/infinity/1e+20/g;
         s/No Solution//g' $SOLFILE.sol > $SOLFILE
fi

# remove all temporary files
rm -f $FIBERSETTINGS $SCIPSETTINGS $SOLFILE.sol
