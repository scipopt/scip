#!/bin/bash
LASTGITHASH=0
GITHASH=0
TIME=3600
LOCK=false
CONTINUE=true
OPTS=(opt dbg)
TESTS=(short miplib2010)

# pull new version form git repository
git pull

# get current git hash
GITHASH=`git describe --always --dirty`

# check if a new version exits
while [ true ]
do
    # wait until new version is available
    while [ "$GITHASH" = "$LASTGITHASH" ]
    do
        echo "wait for new version"
        sleep 60

        # pull new version form git repository
        git pull

        # get current git hash
        GITHASH=`git describe --always --dirty`
    done

    echo $GITHASH
    LASTGITHASH=$GITHASH

    echo "run regression test"

    for OPT in ${OPTS[@]}
    do
        # compile SCIP in debug and opt mode
        make OPT=$OPT VERSION=$GITHASH ZIMPL=false -j2

        for TEST in ${TESTS[@]}
        do
            # run test
            make OPT=$OPT VERSION=$GITHASH TIME=$TIME LOCK=$LOCK CONTINUE=$CONTINUE TEST=$TEST test

            # remove tex and pav file
            rm -f check/results/check.$TEST.*$GITHASH*tex
            rm -f check/results/check.$TEST.*$GITHASH*pav

            # check if fail occurs
            NFAILS=`grep -c fail check/results/check.$TEST.*$GITHASH.*.$OPT.*res`
            NABORTS=`grep -c abort check/results/check.$TEST.*$GITHASH.*.$OPT.*res`
            NREADERRORS=`grep -c readerror check/results/check.$TEST.*$GITHASH.*.$OPT.*res`

            HOSTNAME=`hostname`
            EMAILTO=lpip-developers@zib.de
            EMAILFROM="Git <git@zib.de>"

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

            if [ -f "check/results/check.$TEST.*$LASTGITHASH.*.$OPT.*res" ]
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

        # in any case send a mail to heinz@zib.de
        EMAILTO=heinz@zib.de
        SUBJECT="[$HOSTNAME] [OPT=$OPT] [GITHASH: $GITHASH] $TEST"
        tail check/results/check.$TEST.*$GITHASH.*.$OPT.*.res | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
        done
    done
done
