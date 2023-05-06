#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# make SOLVER=nuopt TEST=benchmark TIME=3600 test

SOLVER=$1
BINNAME=$2
NAME=$3
TIMELIMIT=$4
MEMLIMIT=$5   # The memory limit (in MB) # currently unused
SOLFILE=$6
THREADS=$7
MIPGAP=$8

PRMFILE=nuopt.prm
RESULT=solver.sol


# begin of parameter file
echo begin                             >  $PRMFILE

# set threads to given value
# nothing to be done for NUOPT

# set mipgap to given value
echo branch:gaptol=$MIPGAP             >> $PRMFILE

# set timing to wall-clock time and pass time limit
echo branch:maxtim=$TIMELIMIT          >> $PRMFILE

# set memory limit
echo branch:maxmem=$MEMLIMIT           >> $PRMFILE

# use deterministic mode (warning if not possible)
# nothing to be done for NUOPT

# set display options
#echo branch:out=debug                  >> $PRMFILE
#echo outputh:name=$RESULT              >> $PRMFILE

# end of parameter file
echo end                               >>  $PRMFILE

# read, optimize and exit
gunzip -dc $NAME | $BINNAME


# check the status
INFEASIBLE=`grep "infeasible" $RESULT | wc -l`
if [ $INFEASIBLE = 0 ]
then
# translate NUOPT solution format into format for solution checker. The
# SOLFILE format is a very simple format where in each line we have a
# <variable, value> pair, separated by spaces.  A variable name of
# =obj= is used to store the objective value of the solution, as
# computed by the solver. A variable name of =infeas= can be used to
# indicate that an instance is infeasible.
    OBJ=`grep "VALUE_OF_OBJECTIVE" $RESULT | sed 's/VALUE_OF_OBJECTIVE//g'`
    echo "=obj=                             " $OBJ > $SOLFILE
    grep "V#" $RESULT | awk '{printf "%s %50s\n", $3, $4}' >> $SOLFILE
else
    echo "=infeas=" > $SOLFILE
fi
