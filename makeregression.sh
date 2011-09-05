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
	make OPT=$OPT VERSION=$GITHASH ZIMPL=false


	for TEST in ${TESTS[@]}
	do
	    # run test 
	    make OPT=$OPT VERSION=$GITHASH TIME=$TIME LOCK=$LOCK CONTINUE=$CONTINUE TEST=$TEST test
    
            # remove tex and pav file
	    rm -f check/results/check.$TEST.*$GITHASH*tex
	    rm -f check/results/check.$TEST.*$GITHASH*pav

            # check if fail occurs
	    NFAILS=`grep -c fail check/results/check.$TEST.*$GITHASH*$OPT*res`
echo $NFAILS
	    if [ $NFAILS -gt 0 ]; 
	    then
		HOSTNAME=`hostname`
		SUBJECT="[Bug] [$HOSTNAME] [$GITHASH] $TEST $NFAILS fails"
		EMAILTO=heinz@zib.de
		EMAILFROM=`whoami`
		#metasend -z -b -t $EMAILTO -F $EMAILFROM -s "$SUBJECT" -f "bla" -e quoted-printable -m text/plain
		echo "hallo" | mailx -s "$SUBJECT" $EMAILTO
	    else
		HOSTNAME=`hostname`
		SUBJECT="[OK] [$HOSTNAME] [$GITHASH] $TEST "
		EMAILTO=heinz@zib.de
		EMAILFROM=`whoami`
		#metasend -z -b -t $EMAILTO -F $EMAILFROM -s "$SUBJECT" -f "bla" -e quoted-printable -m text/plain
		echo "hallo" | mailx -s "$SUBJECT" $EMAILTO
	    fi
	done
    done
done
