#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$BINNAME.$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
THREADS=$8
FEASTOL=$9
DISPFREQ=${10}
CONTINUE=${11}

# check if all variables defined (by checking the last one)
if test -z $CONTINUE
then
    echo Skipping test since not all variables are defined
    echo "TSTNAME       = $TSTNAME"
    echo "BINNAME       = $BINNAME"
    echo "SETNAMES      = $SETNAME"
    echo "BINID         = $BINID"
    echo "TIMELIMIT     = $TIMELIMIT"
    echo "NODELIMIT     = $NODELIMIT"
    echo "MEMLIMIT      = $MEMLIMIT"
    echo "THREADS       = $THREADS"
    echo "FEASTOL       = $FEASTOL"
    echo "DISPFREQ      = $DISPFREQ"
    echo "CONTINUE      = $CONTINUE"
    exit 1;
fi

# get current SCIP path
SCIPPATH=`pwd`

if test ! -e $SCIPPATH/results
then
    mkdir $SCIPPATH/results
fi

OUTFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.out
ERRFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.err
RESFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.res
TEXFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.tex
TMPFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.tmp
SETFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.prm

SETTINGS=$SCIPPATH/../settings/$SETNAME.set

# check if the settings file exists
if test $SETNAME != "default"
then
    if test ! -e $SETTINGS
    then
        echo skipping test due to non-existence of settings file $SETTINGS
        exit
    fi
fi

if test "$CONTINUE" = "true"
then
    MVORCP=cp
else
    MVORCP=mv
fi

DATEINT=`date +"%s"`
if test -e $OUTFILE
then
    $MVORCP $OUTFILE $OUTFILE.old-$DATEINT
fi
if test -e $ERRFILE
then
    $MVORCP $ERRFILE $ERRFILE.old-$DATEINT
fi

if test "$CONTINUE" = "true"
then
    LASTPROB=`awk -f getlastprob.awk $OUTFILE`
    echo Continuing benchmark. Last solved instance: $LASTPROB
    echo "" >> $OUTFILE
    echo "----- Continuing from here. Last solved: $LASTPROB -----" >> $OUTFILE
    echo "" >> $OUTFILE
else
    LASTPROB=""
fi



uname -a >>$OUTFILE
uname -a >>$ERRFILE
date >>$OUTFILE
date >>$ERRFILE

# we add 100% to the hard time limit and additional 600 seconds in case of small time limits
# NOTE: the jobs should have a hard running time of more than 5 minutes; if not so, these
#       jobs get automatically assigned in the "express" queue; this queue has only 4 CPUs
#       available
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 600\` + $TIMELIMIT`

# we add 10% to the hard memory limit and additional 100mb to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``
HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024`

EVALFILE=$SCIPPATH/results/check.$TSTNAME"_"$THREADS.$BINID.$QUEUE.$SETNAME.eval
echo > $EVALFILE

for i in `cat testset/$TSTNAME.test` DONE
do
  if test "$i" = "DONE"
  then
      break
  fi

  if test "$LASTPROB" = ""
  then
      LASTPROB=""
      # check if problem instance exists
      if test -f $SCIPPATH/$i
      then
          if test -e $SETFILE
          then
	      rm -f $SETFILE
          fi
          echo @01 $i ===========
          echo @01 $i ===========                 >> $ERRFILE

          echo > $TMPFILE
          echo ""                                 > $TMPFILE
          if test $SETNAME != "default"
          then
              echo "non-default settings not yet supported"
          fi
          echo maxtime = -$TIMELIMIT              >> $TMPFILE
          # use wallclock time
          echo cputime = 0                        >> $TMPFILE
          echo miprelstop = $MIPGAP               >> $TMPFILE
          if test $FEASTOL != "default"
          then
              echo miptol = $FEASTOL              >> $TMPFILE
          fi
          echo maxnode = $NODELIMIT               >> $TMPFILE
          echo threads = $THREADS                 >> $TMPFILE

          # the following should enforce the memory limit, but still it is only
          # a soft limit and a temporary file is also written
          echo treememorylimit = $MEMLIMIT        >> $TMPFILE
          echo treememorysavingtarget = 0.0       >> $TMPFILE

          echo format \"timelimit %g\" \$maxtime  >> $TMPFILE
          echo format \"mipgap %g\" \$miprelstop  >> $TMPFILE
          echo format \"feastol %g\" \$miptol     >> $TMPFILE
          echo format \"nodelimit %g\" \$maxnode  >> $TMPFILE
          echo format \"memlimit %g\" \$treememorylimit >> $TMPFILE
          echo format \"percentmemtofile %g\" \$treememorysavingtarget >> $TMPFILE

          echo readprob $SCIPPATH/$i              >> $TMPFILE
          echo time mipoptimize                   >> $TMPFILE
	  echo echo simplexiter:                  >> $TMPFILE
	  echo simplexiter                        >> $TMPFILE
          echo quit                               >> $TMPFILE
          echo -----------------------------
          date
          date >>$ERRFILE
          echo -----------------------------
          date +"@03 %s"
          bash -c "ulimit -t $HARDTIMELIMIT s; ulimit -v $HARDMEMLIMIT k; ulimit -f 1000000; $BINNAME < $TMPFILE" 2>>$ERRFILE
          date +"@04 %s"
          echo -----------------------------
          date
          date >>$ERRFILE
          echo -----------------------------
          echo =ready=
      else
          echo @02 FILE NOT FOUND: $i ===========
          echo @02 FILE NOT FOUND: $i =========== >>$ERRFILE
      fi
  else
      echo skipping $i
      if test "$LASTPROB" = "$i"
      then
          LASTPROB=""
      fi
  fi
done | tee -a $OUTFILE

rm -f $TMPFILE

date >>$OUTFILE
date >>$ERRFILE

./evalcheck_xpress.sh $OUTFILE
