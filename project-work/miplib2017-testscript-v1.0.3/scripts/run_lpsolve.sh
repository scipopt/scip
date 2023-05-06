#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

SOLVER=$1
BINNAME=$2
NAME=$3
TIMELIMIT=$4
MEMLIMIT=$5   # The memory limit (in MB) # currently unused
SOLFILE=$6
THREADS=$7
MIPGAP=$8

echo > $SOLFILE

# set threads to given value
# nothing to be done for lpsolve

# set mipgap to given value
PARAMETERS=$PARAMETERS" -ga $MIPGAP"

# set timing to wall-clock time and pass time limit
PARAMETERS=$PARAMETERS" -timeout $TIMELIMIT"

# use deterministic mode (warning if not possible)
# nothing to be done for GLPK

# check for lp format
MPS=`echo $NAME | grep "\.mps"`
if test $MPS
then
    PARAMETERS=$PARAMETERS" -fmps"
fi

# set display options
PARAMETERS=$PARAMETERS" -v4"

# set solution output
PARAMETERS=$PARAMETERS" -S2"

# run lpsolve
$BINNAME $PARAMETERS $NAME | tee $SOLFILE

if test -e $SOLFILE
then
    # translate lpsolve solution format into format for solution checker. The
    # SOLFILE format is a very simple format where in each line we have a
    # <variable, value> pair, separated by spaces.  A variable name of
    # =obj= is used to store the objective value of the solution, as
    # computed by the solver. A variable name of =infeas= can be used to
    # indicate that an instance is infeasible.
    awk '
    BEGIN {
       start = 0;
    }
    //{
       if( start )
          printf("%s\n", $0);
    }
    /^Actual values of the variables:/ {
       start = 1;
    }
    /^Optimal solution/ {
       printf("=obj= %s\n", $3);
    }' $SOLFILE > $SOLFILE.tmp

    mv $SOLFILE.tmp $SOLFILE
fi
