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

echo > $TMPFILE
echo > $SOLFILE

# set mipgap to given value
echo set gap_limit $MIPGAP             >> $TMPFILE

# set timing to wall-clock time and
# changing clock type is not possible

# pass time limit
echo set time_limit $TIMELIMIT         >> $TMPFILE

# use deterministic mode (warning if not possible)
# nothing to be done for SCIP

# read, optimize, display statistics, write solution, and exit
echo load $NAME                        >> $TMPFILE
echo solve                             >> $TMPFILE
echo display stats                     >> $TMPFILE
echo display solution                  >> $TMPFILE
echo quit                              >> $TMPFILE

$BINNAME < $TMPFILE | tee $SOLFILE

# translate SYMPHONY solution format into format for solution checker. The
# SOLFILE format is a very simple format where in each line we have a
# <variable, value> pair, separated by spaces.  A variable name of
# =obj= is used to store the objective value of the solution, as
# computed by the solver. A variable name of =infeas= can be used to
# indicate that an instance is infeasible.
awk '
BEGIN {
  start = 0;
}
/^\+\+\+\+\+\+/{
   next;
}
/^SYMPHONY:/ {
   next;
}

//{
  if( start )
     printf("%s\n", $0);
}
/^Nonzero column names and values in the solution/ {
   start = 1;
}
/^Solution Cost:/ {
   printf("=obj= %s\n", $3);
}' $SOLFILE > $TMPFILE

mv $TMPFILE $SOLFILE
