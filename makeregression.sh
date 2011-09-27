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

# email variables
ADMINEMAIL=heinz@zib.de
EMAILTO=lpip-developers@zib.de
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
CONTINUE=true
OPTS=(dbg opt)
TESTS=(short miplib2010)

# move into the SCIP directory; this is necessary due the cron job 
cd $SCIPDIR

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
echo $SCRIPTTIMESTAMP | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL

# pull new version form git repository
git pull

# get maybe new time stamp for makeregression.sh 
NEWSCRIPTTIMESTAMP=`stat -c %Y $SCRIPTNAME`

# get current git hash
GITHASH=`git describe --always --dirty`

# continue testing if makeregression.sh did not change
while [ $NEWSCRIPTTIMESTAMP -eq $SCRIPTTIMESTAMP ]
do
    # if file named "kill" exist, stop the regression test
    if [ -f "kill" ];
    then
	SUBJECT="[$HOSTNAME] killed $SCRIPTNAME"
	echo "kill script due to exists of file \"kill\" (rm kill)" 
	echo "killed" | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL
	exit;
    fi

    echo "run regression test"

    for OPT in ${OPTS[@]}
    do
        # compile SCIP in debug or opt mode
        make OPT=$OPT VERSION=$GITHASH ZIMPL=false clean
        make OPT=$OPT VERSION=$GITHASH ZIMPL=false -j2
	
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

            # check read fails
            if [ $NFAILS -gt 0 ];
            then
                SUBJECT="[FAIL] [$HOSTNAME] [OPT=$OPT] [GITHASH: $GITHASH] $TEST"
                grep fail check/results/check.$TEST.*$GITHASH.*.$OPT.*.res | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
            fi

            # check read errors
            if [ $NREADERRORS -gt 0 ];
            then
                SUBJECT="[READERROR] [$HOSTNAME] [OPT=$OPT] [GITHASH: $GITHASH] $TEST"
                grep readerror check/results/check.$TEST.*$GITHASH.*.$OPT.*.res | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
            fi

            # check aborts
            if [ $NABORTS -gt 0 ];
            then
                SUBJECT="[ABORT] [$HOSTNAME] [OPT=$OPT] [GITHASH: $GITHASH] $TEST"
                grep abort check/results/check.$TEST.*$GITHASH.*.$OPT.*.res | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
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
            tail -8 check/results/check.$TEST.*$GITHASH.*.$OPT.*.res | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL
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
	
	NEWSCRIPTTIMESTAMP=`stat -c %Y $SCRIPTNAME`
        
        # get current git hash
        GITHASH=`git describe --always --dirty`
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
echo $NEWSCRIPTTIMESTAMP | mailx -s "$SUBJECT" -r "$EMAILFROM" $ADMINEMAIL
