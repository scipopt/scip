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


echo $SOLFILE

# set threads to given value
# set mipgap to given value
# set timing to wall-clock time and pass time limit
# read, optimize, display statistics, write solution, and exit
echo $BINNAME  -solfile $SOLFILE -threads $THREADS -time $TIMELIMIT $NAME

$BINNAME -solfile $SOLFILE -threads $THREADS -time $TIMELIMIT $NAME

# remove MIPCL log file
rm *_mipcl.log
