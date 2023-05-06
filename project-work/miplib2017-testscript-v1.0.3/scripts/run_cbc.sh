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

echo > $SOLFILE

if test $THREADS != 0
then
    $BINNAME -import $NAME -sec $TIMELIMIT -threads $THREADS -ratio $MIPGAP -timeMode elapsed -solve -solution $SOLFILE
else
    $BINNAME -import $NAME -sec $TIMELIMIT -ratio $MIPGAP -solve -solution $SOLFILE
fi

# translate CBC solution format into format for solution checker.
#  The SOLFILE format is a very simple format where in each line
#  we have a <variable, value> pair, separated by spaces.
#  A variable name of =obj= is used to store the objective value
#  of the solution, as computed by the solver. A variable name of
#  =infeas= can be used to indicate that an instance is infeasible.
awk '
   /^Stopped/ {
      if( NF > 7 )
	 exit;

      printf ("=obj= %s \n", $7);
      next;
   }
   /^Optimal/ {
      if( $5 == "-" )
         printf ("=obj= %s \n", $8);
      else
         printf ("=obj= %s \n", $5);
      next;
   }
   /^Infeasible/ {
      printf ("=infeas= \n");
      exit;
   }
   /^Integer/ {
      if( $2 == "infeasible")
	 printf ("=infeas= \n");
      exit;
   }
   //{
      printf ("%s %s \n", $2, $3);
   }' $SOLFILE | tee $TMPFILE

cp $SOLFILE $SOLFILE.sol
mv $TMPFILE $SOLFILE
