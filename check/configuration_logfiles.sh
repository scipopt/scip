#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

### configures the right test output files such as the .eval, the .tmp and the .set
### files to run a test on.
### the invoking script should pass "init" as argument to signalize that
### files need to be reset

### environment variables declared in this script
### OUTFILE - the name of the (sticked together) output file
### ERRFILE - the name of the (sticked together) error file
### EVALFILE - evaluation file to glue single output and error files together
### OBJECTIVEVAL - the optimal or best-know objective value for this instance
### SHORTPROBNAME - the basename of $INSTANCE without file extension
### FILENAME - the basename of the local files (.out, .tmp, and .err)
### SKIPINSTANCE - should the instance be skipped because it was already evaluated in a previous setting?
### BASENAME - $SCIPPATH/$OUTPUTDIR/$FILENAME cf. FILENAME argument
### TMPFILE  - the batch file name to pass for solver instructions
### SETFILE  - the name of the settings file to save solver settings to

### environment variables passed as arguments to this script
INIT=$1      # should log files be initialized (this overwrite or copy/move some existing log files)
COUNT=$2     # the instance count as part of the filename
INSTANCE=$3  # the name of the instance
BINID=$4     # the ID of the binary to use
PERMUTE=$5   # the number of permutations to use - 0 for no permutation
SEEDS=$6     # the number of random seeds - 0 only default seeds
SETNAME=$7   # the name of the setting
TSTNAME=$8   # the name of the testset
CONTINUE=$9  # should test continue an existing run
QUEUE=${10}  # the queue name
p=${11}      # the index of the current permutation
s=${12}      # shift of the global random seed
THREADS=${13} # the number of threads

# common naming scheme for eval files
EVALFILE=$SCIPPATH/$OUTPUTDIR/check.$TSTNAME.$BINID.$QUEUE.$SETNAME

# if number of threads is larger than 1, add postfix
if test $THREADS -gt 1
then
    EVALFILE=$EVALFILE"-t"$THREADS
fi

# if seed is positive, add postfix
SEED=`expr $s + $GLBSEEDSHIFT`
if test $SEED -gt 0
then
    EVALFILE=$EVALFILE"-s"$SEED
fi

# if permutation is positive, add postfix
if test $p -gt 0
then
    EVALFILE=$EVALFILE"-p"$p
fi

OUTFILE=$EVALFILE.out
ERRFILE=$EVALFILE.err

# add .eval extension to evalfile
EVALFILE=$EVALFILE.eval

# create meta file
if test -e $EVALFILE
then
    fname=$SCIPPATH/$OUTPUTDIR/`basename $EVALFILE .eval`.meta
    if ! test -e $fname
    then
        echo @Permutation $p > $fname
        echo @Seed $SEED >> $fname
        echo @Settings $SETNAME >> $fname
        echo @TstName $TSTNAME >> $fname
        echo @BinName $BINNAME >> $fname
        echo @NodeLimit $NODELIMIT >> $fname
        echo @MemLimit $MEMLIMIT >> $fname
        echo @Threads $THREADS >> $fname
        echo @FeasTol $FEASTOL >> $fname
        echo @Queue $QUEUE >> $fname
        echo @Exclusive $EXCLUSIVE >> $fname
    fi
fi

if test "$INSTANCE" = "DONE"
then
    return
fi

# reset files if flag is set to 'init'
if test $INIT = "true"
then
    #reset the eval file
    echo > $EVALFILE

    #mv existing out and error files
    if test "$CONTINUE" = "true"
    then
        MVORCP=cp
    else
        MVORCP=mv
    fi
    DATEINT=`date +"%s"`
    for FILE in OUTFILE ERRFILE
    do
        if test -e $FILE
        then
            $MVORCP $FILE $FILE.old-$DATEINT
        fi
    done
fi


# filter all parseable file format extensions
SHORTPROBNAME=`basename $INSTANCE .gz`
for EXTENSION in .mps .lp .opb .gms .pip .zpl .cip .fzn .osil .wbo .cnf .difflist
do
    SHORTPROBNAME=`basename $SHORTPROBNAME $EXTENSION`
done

# get objective value from solution file
# we do this here to have it available for all solvers, even though it is not really related to logfiles
if test -e "$SOLUFILE"
then
    # get the objective value from the solution file: grep for the instance name and only use entries with an optimal or best known value;
    # if there are multiple entries for this instance in the solution file, sort them by objective value and take the objective value
    # written in the last line, i.e., the largest value;
    # as a double-check, we do the same again, but reverse the sorting to get the smallest value
    OBJECTIVEVAL=`grep " $SHORTPROBNAME " $SOLUFILE | grep -e =opt= -e =best= | sort -k 3 -g | tail -n 1 | awk '{print $3}'`
    CHECKOBJECTIVEVAL=`grep " $SHORTPROBNAME " $SOLUFILE | grep -e =opt= -e =best= | sort -k 3 -g -r | tail -n 1 | awk '{print $3}'`

    # if largest and smalles reference value given in the solution file differ by more than 1e-04, stop because of this inconsistency
    if awk -v n1="$OBJECTIVEVAL" -v n2="$CHECKOBJECTIVEVAL" 'BEGIN { exit (n1 <= n2 + 0.0001 && n2 <= n1 + 0.0001) }' /dev/null;
    then
	echo "Exiting test because objective value in solu file is inconsistent: $OBJECTIVEVAL vs. $CHECKOBJECTIVEVAL"
        exit
    fi
else
    OBJECTIVEVAL=""
fi
#echo "Reference value $OBJECTIVEVAL $SOLUFILE"

NEWSHORTPROBNAME=`echo $SHORTPROBNAME | cut -c1-25`
SHORTPROBNAME=$NEWSHORTPROBNAME

#define file name for temporary log file
FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTPROBNAME.$BINID.$QUEUE.$SETNAME

# if number of threads is larger than 1, add postfix
if test $THREADS -gt 1
then
    FILENAME=$FILENAME"-t"$THREADS
fi

# if seed is positive, add postfix
if test $SEED -gt 0
then
    FILENAME=$FILENAME"-s"$SEED
fi

# if permutation is positive, add postfix
if test $p -gt 0
then
    FILENAME=$FILENAME"-p"$p
fi

SKIPINSTANCE="false"
# in case we want to continue we check if the job was already performed
if test "$CONTINUE" = "true" && test -e $OUTPUTDIR/$FILENAME.out
then
    echo skipping file $INSTANCE due to existing output file $OUTPUTDIR/$FILENAME.out
    SKIPINSTANCE="true"
fi

# configure global names TMPFILE (batch file) and SETFILE to save settings to
BASENAME=$SCIPPATH/$OUTPUTDIR/$FILENAME
TMPFILE=$BASENAME.tmp
SETFILE=$BASENAME.set

# even if we decide to skip this instance, we write the basename to the eval file
echo $BASENAME >> $EVALFILE
