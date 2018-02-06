#!/bin/bash
#
# This script runs a regression test for SCIP
#
# As input it optionally expects the parameter LPS - default is LPS=spx
#

SCRIPTNAME="makeregression.sh"

SCIPDIR=`pwd`

LPS=$1

# check which LP solver to use
if [ -z ${LPS} ];
then
    LPS=spx
fi

# email variables
ADMINEMAIL=miltenberger@zib.de
LPIPDEVELOPERSEMAIL=lpip-developers@zib.de
EMAILTO=$LPIPDEVELOPERSEMAIL
EMAILFROM="Git <git@zib.de>"
HOSTNAME=`hostname`

# get current time stamp of the script  file
SCRIPTTIMESTAMP=`stat -c %Y $SCRIPTNAME`

GITHASHFILE="check/results/githashorder"
GITHASH=0

# test run variables
MEM=14000
TIME=3600
LOCK=false
CONTINUE=false
OPTS=(dbg opt)
TESTS=(short MMM bugs SAP-MMP)

# find last git hash for which a complete test was performed
if [ -f $GITHASHFILE ];
then
    LASTGITHASH=`tail -1 $GITHASHFILE`
else
    LASTGITHASH=0
fi

# send mail to admin to indicate that makeregression.sh has (re)started
SUBJECT="[$HOSTNAME] (Re)Start $SCRIPTNAME"
echo "Time stamp $SCRIPTTIMESTAMP" | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL


NEWSCRIPTTIMESTAMP=$SCRIPTTIMESTAMP

# continue testing if makeregression.sh did not change
while [ $NEWSCRIPTTIMESTAMP -eq $SCRIPTTIMESTAMP ]
do
    # backup current download lists of SCIP and SoPlex
#     YEAR=`date +%Y`
#     MONTH=`date +%m`
#     DAY=`date +%d`
#     cp /www/Abt-Optimization/scip/counter/users.dat ~/download-counter/scip-$YEAR-$MONTH-$DAY-users.dat
#     cp /www/Abt-Optimization/soplex/counter/users.dat ~/download-counter/soplex-$YEAR-$MONTH-$DAY-users.dat
#     echo "created backups of download statistics"

    for OPT in ${OPTS[@]}
    do
        for TEST in ${TESTS[@]}
        do
            # pull and compile SoPlex
            cd ../soplex
            git pull
            make
            make OPT=dbg

            # pull new version from git repository
            cd ../scip
            git pull

            # to clean the local git repository perform a git status (to avoid a dirty hash)
            git status

            # get current git hash
            GITHASH=`git describe --always --dirty  | sed -re 's/^.+-g//'`

            # compile SCIP in debug or opt mode
            make OPT=$OPT LPSOPT=$OPT VERSION=$GITHASH LPS=$LPS ZIMPL=true clean
            make OPT=$OPT LPSOPT=$OPT VERSION=$GITHASH LPS=$LPS ZIMPL=true

            # run test
            make OPT=$OPT VERSION=$GITHASH LPS=$LPS TIME=$TIME LOCK=$LOCK CONTINUE=$CONTINUE TEST=$TEST MEM=$MEM test

            # remove tex and pav file
            rm -f check/results/check.$TEST.*$GITHASH*tex
            rm -f check/results/check.$TEST.*$GITHASH*pav

            BASEFILE="check/results/check.$TEST.*$GITHASH.*.$OPT.$LPS"

            # check if fail occurs
            NFAILS=`grep -c fail $BASEFILE.*res`

            # construct string which shows the destination of the out, err, and res files
            ERRORFILE=`ls $BASEFILE.*.err`
            OUTFILE=`ls $BASEFILE.*.out`
            RESFILE=`ls $BASEFILE.*.res`
            DESTINATION="$SCIPDIR/$OUTFILE \n$SCIPDIR/$ERRORFILE \n$SCIPDIR/$RESFILE"

            # check read fails
            if [ $NFAILS -gt 0 ];
            then
                SUBJECT="[FAIL] [$HOSTNAME] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH] $TEST"
                ERRORINSTANCES=`grep fail $BASEFILE.*.res`
                echo -e "$ERRORINSTANCES \n$DESTINATION" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
            fi

            # in any case send a mail to admin
            SUBJECT="[$HOSTNAME] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH] $TEST"
            RESULTS=`tail -8 $BASEFILE.*.res`
            echo -e "$RESULTS \n$DESTINATION" | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL
        done
    done

    # last test was performed for the current hash
    echo $GITHASH >> $GITHASHFILE

    LASTGITHASH=$GITHASH

    # wait until new version is available
    while [ "$GITHASH" == "$LASTGITHASH" ]
    do
        echo "wait for new version"
        sleep 60

        # pull new version form git repository
        git pull

        # to clean the local git repository perform a git status (to avoid a dirty hash)
        git status

        NEWSCRIPTTIMESTAMP=`stat -c %Y $SCRIPTNAME`

        # get current git hash
        GITHASH=`git describe --always --dirty  | sed -re 's/^.+-g//'`
    done
done

# send email to admin to indicate  that the script stopped
SUBJECT="[$HOSTNAME] $SCRIPTNAME [LPS=$LPS] has stopped - restart it"
echo "Time stamp $NEWSCRIPTTIMESTAMP" | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL
