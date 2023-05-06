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

TMPFILE=check.$SOLVER.tmp
MPSFILE=check.$SOLVER.mps

cp $NAME $MPSFILE.gz
gunzip $MPSFILE.gz

echo > $SOLFILE

$BINNAME $MPSFILE

rm $MPSFILE
