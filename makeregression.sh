#!/bin/bash
# 
# This scripts a regression for SCIP 
#
# As input it expect optional the path to SCIP. If no path is given the
# current directory is assumed to be the SCIP path. This optional argument
# is need to be able to restart the script via a cron job
#


SCRIPTNAME="makeregression.sh"
SCIPDIR=$1

# check if the SCIP directory is given
if [ -z ${SCIPDIR} ];
then
    SCIPDIR=`pwd`
fi

# move into the SCIP directory; this is necessary due the cron job
cd $SCIPDIR

# email variables
ADMINEMAIL=heinz@zib.de
LPIPDEVELOPERSEMAIL=lpip-developers@zib.de
EMAILFROM="Git <git@zib.de>"
HOSTNAME=`hostname`

# get current time stamp of the script  file
SCRIPTTIMESTAMP=`stat -c %Y $SCRIPTNAME`

CRONTABFILE="crontabfile"

GITHASHFILE="check/results/githashorder"
GITHASH=0

# test run variables
MEM=6544
TIME=3600
LOCK=false
CONTINUE=false
OPTS=(dbg opt)
TESTS=(short MMM bugs)

# first delete cron jobs if one exists
crontab -r

# check if the script exist; if not we have the stop 
if [ ! -f $SCRIPTNAME ];
then
    SUBJECT="[$HOSTNAME] killed $SCRIPTNAME"
    echo "kill script due to not exists of file $SCRIPTNAME in directory $SCIPDIR" 
    echo "killed" | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL
    exit;
fi

# if file named "kill" exist, stop the regression test
if [ -f "kill" ];
then
    SUBJECT="[$HOSTNAME] killed $SCRIPTNAME"
    echo "kill script due to exists of file \"kill\" (rm kill)" 
    echo "killed" | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL
    exit;
fi

# find last git hash for which a complete test was performed
if [ -f $GITHASHFILE ];
then
    LASTGITHASH=`tail -1 $GITHASHFILE`
else
    LASTGITHASH=0
fi

# send mail to admin to indicate that makeregression.sh (re)start 
SUBJECT="[$HOSTNAME] (Re)Start $SCRIPTNAME"
echo "Time stamp $SCRIPTTIMESTAMP" | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL

# pull new version form git repository
git pull

# to clean the local git repository perform a git st (to avoid a dirty hash)
git st

# get maybe new time stamp for makeregression.sh 
NEWSCRIPTTIMESTAMP=`stat -c %Y $SCRIPTNAME`

# get current git hash
GITHASH=`git describe --always --dirty  | sed -re 's/^.+-g//'`

# continue testing if makeregression.sh did not change
while [ $NEWSCRIPTTIMESTAMP -eq $SCRIPTTIMESTAMP ]
do
    # backup current download lists of SCIP and SoPlex
    YEAR=`date +%Y`
    MONTH=`date +%m`
    DAY=`date +%d`
    cp /www/Abt-Optimization/scip/counter/users.dat ~/download-counter/scip-$YEAR-$MONTH-$DAY-users.dat
    cp /www/Abt-Optimization/soplex/counter/users.dat ~/download-counter/soplex-$YEAR-$MONTH-$DAY-users.dat
    echo "created backups of download statistics"

    # if file named "kill" exist, stop the regression test
    if [ -f "kill" ];
    then
	SUBJECT="[$HOSTNAME] killed $SCRIPTNAME"
	echo "kill script due to exists of file \"kill\" (rm kill)" 
	echo "killed" | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL
	exit;
    fi

    echo "get current SoPlex version"
    cd ../soplex
    git pull
    make
    cd ../scip

    echo "run regression test"

    for OPT in ${OPTS[@]}
    do
        # compile SCIP in debug or opt mode
        make OPT=$OPT VERSION=$GITHASH ZIMPL=true clean
        make OPT=$OPT VERSION=$GITHASH ZIMPL=true
	
        for TEST in ${TESTS[@]}
        do
            # run test
            make OPT=$OPT VERSION=$GITHASH TIME=$TIME LOCK=$LOCK CONTINUE=$CONTINUE TEST=$TEST MEM=$MEM test

            # remove tex and pav file
            rm -f check/results/check.$TEST.*$GITHASH*tex
            rm -f check/results/check.$TEST.*$GITHASH*pav

            # check if fail occurs
            NFAILS=`grep -c fail check/results/check.$TEST.*$GITHASH.*.$OPT.*res`
            NABORTS=`grep -c abort check/results/check.$TEST.*$GITHASH.*.$OPT.*res`
            NREADERRORS=`grep -c readerror check/results/check.$TEST.*$GITHASH.*.$OPT.*res`

	    # only send fail mail to the group if the occurs on MMM or short
	    if [ "$TEST" == "bugs" ]
	    then
		RECEIVER=$ADMINEMAIL
	    else
		EMAILTO=$LPIPDEVELOPERSEMAIL
	    fi

	    # construct string which shows the destination of the out, err, and res files
	    ERRORFILE=`ls check/results/check.$TEST.*$GITHASH.*.$OPT.*.err`
	    OUTFILE=`ls check/results/check.$TEST.*$GITHASH.*.$OPT.*.out`
	    RESFILE=`ls check/results/check.$TEST.*$GITHASH.*.$OPT.*.res`
	    DESTINATION="$SCIPDIR/$OUTFILE \n$SCIPDIR/$ERRORFILE \n$SCIPDIR/$RESFILE"

            # check read fails
            if [ $NFAILS -gt 0 ];
            then
                SUBJECT="[FAIL] [$HOSTNAME] [OPT=$OPT] [GITHASH: $GITHASH] $TEST"
		ERRORINSTANCES=`grep fail check/results/check.$TEST.*$GITHASH.*.$OPT.*.res`
                echo -e "$ERRORINSTANCES \n$DESTINATION" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
            fi

            # check read errors
            if [ $NREADERRORS -gt 0 ];
            then
                SUBJECT="[READERROR] [$HOSTNAME] [OPT=$OPT] [GITHASH: $GITHASH] $TEST"
		ERRORINSTANCES=`grep readerror check/results/check.$TEST.*$GITHASH.*.$OPT.*.res`
                echo -e "$ERRORINSTANCES \n$DESTINATION" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
            fi

            # check aborts
            if [ $NABORTS -gt 0 ];
            then
                SUBJECT="[ABORT] [$HOSTNAME] [OPT=$OPT] [GITHASH: $GITHASH] $TEST"
		ERRORINSTANCES=`grep abort check/results/check.$TEST.*$GITHASH.*.$OPT.*.res`
                echo -e "$ERRORINSTANCES \n$DESTINATION" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
            fi

            # check performance in opt mode
            if [ "$OPT" == "opt" ];
            then
                NOK=`grep -c ok check/results/check.$TEST.*$GITHASH.*.$OPT.*res`
                NSOLVED=`grep -c solved check/results/check.$TEST.*$GITHASH.*.$OPT.*res`
                NTIMEOUTS=`grep -c timeouts check/results/check.$TEST.*$GITHASH.*.$OPT.*res`

		if [ -f "check/results/check.$TEST.*$LASTGITHASH.*.$OPT.*res" ];
		then
                    NLASTTIMEOUTS=`grep -c timeouts check/results/check.$TEST.*$LASTGITHASH.*.$OPT.*res`

                    # check time outs
                    if [ $NTIMEOUTS -gt NLASTTIMEOUTS ];
                    then
			SUBJECT="[TIMEOUT] [$HOSTNAME] [OPT=$OPT] [GITHASH: $GITHASH] $TEST"
			grep timeouts check/results/check.$TEST.*$GITHASH.*.$OPT.*.res | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
                    fi
		fi
            fi

            # in any case send a mail to admin
            SUBJECT="[$HOSTNAME] [OPT=$OPT] [GITHASH: $GITHASH] $TEST"
	    RESULTS=`tail -8 check/results/check.$TEST.*$GITHASH.*.$OPT.*.res`
            echo -e "$RESULTS \n$DESTINATION" | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL
        done
    done
    
    # all test are performed for the current hash
    echo $GITHASH >> $GITHASHFILE
    
    echo $GITHASH
    LASTGITHASH=$GITHASH
    
    # wait until new version is available
    while [ "$GITHASH" == "$LASTGITHASH" ]
    do
        echo "wait for new version"
	sleep 60

        # pull new version form git repository
        git pull
	
        # to clean the local git repository perform a git st (to avoid a dirty hash)
	git st

	NEWSCRIPTTIMESTAMP=`stat -c %Y $SCRIPTNAME`
        
        # get current git hash
        GITHASH=`git describe --always --dirty  | sed -re 's/^.+-g//'`
    done
done

# create crontab file for restarting the script
echo "* * * * * $SCIPDIR/$SCRIPTNAME $SCIPDIR"  > $CRONTABFILE

# set cron job to to restart script
crontab $CRONTABFILE

# remove cron tab file 
rm -f $CRONTABFILE

# send email to admin to indicate  that the script stopped
SUBJECT="[$HOSTNAME] Stop $SCRIPTNAME"
echo "Time stamp $NEWSCRIPTTIMESTAMP" | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL
