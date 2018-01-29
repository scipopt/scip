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
### SHORTPROBNAME - the basename of $INSTANCE without file extension
### FILENAME - the basename of the local files (.out, .tmp, and .err)
### EVALFILE - evaluation file to glue single output and error files together
### SKIPINSTANCE - should the instance be skipped because it was already evaluated in a previous setting?
### BASENAME - $SCIPPATH/$OUTPUTDIR/$FILENAME cf. FILENAME argument
### TMPFILE  - the batch file name to pass for solver instructions
### SETFILE  - the name of the settings file to save solver settings to

### environment variables passed as arguements to this script
INIT=$1      # should log files be initialized (this overwrite or copy/move some existing log files)
COUNT=$2     # the instance count as part of the filename
INSTANCE=$3  # the name of the instance
BINID=$4     # the ID of the binary to use
PERMUTE=$5   # the number of permutations to use - 0 for no permutation
SEEDS=$6     # the number of random seeds - 0 only default seeds
SETNAME=$7   # the name of the setting
TSTNAME=$8   # the name of the testset
CONTINUE=$9  # should test continue an existing run
QUEUE=${10}    # the queue name
p=${11}      # the index of the current permutation
s=${12}      # shift of the global random seed

OUTFILE=$OUTPUTDIR/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.out
ERRFILE=$OUTPUTDIR/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.err

# if number of permutations is positive, add postfix
if test $p -gt 0
then
    # if number of seeds is positive, add postfix
    if test $s -gt 0
    then
        EVALFILE=$SCIPPATH/$OUTPUTDIR/check.$TSTNAME.$BINID.$QUEUE.$SETNAME"-s"$s"-p"$p.eval
    else
        EVALFILE=$SCIPPATH/$OUTPUTDIR/check.$TSTNAME.$BINID.$QUEUE.$SETNAME"-p"$p.eval
    fi
else
    # if number of seeds is positive, add postfix
    if test $s -gt 0
    then
        EVALFILE=$SCIPPATH/$OUTPUTDIR/check.$TSTNAME.$BINID.$QUEUE.$SETNAME"-s"$s.eval
    else
        EVALFILE=$SCIPPATH/$OUTPUTDIR/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.eval
    fi
fi


if test -e $EVALFILE
then
    fname=$SCIPPATH/$OUTPUTDIR/`basename $EVALFILE .eval`.meta
    if ! test -e $fname
    then
        echo @Permutation $p > $fname
        echo @Seed $s >> $fname
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
NEWSHORTPROBNAME=`echo $SHORTPROBNAME | cut -c1-25`
SHORTPROBNAME=$NEWSHORTPROBNAME

# if number of permutations is positive, add postfix
if test $p -gt 0
then
    # if number of seeds is positive, add postfix
    if test $s -gt 0
    then
        FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTPROBNAME.$BINID.$QUEUE.$SETNAME-"s"$s-"p"$p
    else
        FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTPROBNAME.$BINID.$QUEUE.$SETNAME-"p"$p
    fi
else
    # if number of seeds is positive, add postfix
    if test $s -gt 0
    then
        FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTPROBNAME.$BINID.$QUEUE.$SETNAME-"s"$s
    else
        FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTPROBNAME.$BINID.$QUEUE.$SETNAME
    fi
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
# even if we decide skip this instance, we write the basename to the eval file
echo $BASENAME >> $EVALFILE
