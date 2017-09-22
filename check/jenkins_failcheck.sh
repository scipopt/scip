#!/bin/bash
#
# This script uploads and checks for fails in a SCIP run.
# Sends an email if errors are detected. Is not meant to be use directly,
# but to be called by jenkins_check_results.sh.
# Note: TESTSET, GITHASH, etc are read from the environment, see
# jenkins_check_results.sh

sleep 5

# we use a name that is unique per test sent to the cluster (a jenkins job
# can have several tests sent to the cluster, that is why the jenkins job
# name (i.e, the directory name) is not enough)
DATABASE="/nfs/OPTI/adm_timo/databases/${PWD##*/}_${TESTSET}_$SETTING.txt"
TMPDATABASE="$DATABASE.tmp"
STILLFAILING="${DATABASE}_SF.tmp"
RBDB="${DATABASE}_rbdb.txt"
OUTPUT="${DATABASE}_output.tmp"
touch ${STILLFAILING}
touch ${RBDB}

# the first time, the file might not exists so we create it
# Even more, we have to write something to it, since otherwise
# the awk scripts below won't work (NR and FNR will not be different)
echo "Preparing database."
if ! [[ -s $DATABASE ]]; then  # check that file exists and has size larger that 0
  echo "Instance Fail_reason Branch Testset Setting Opt_mode LPS" > $DATABASE
fi

EMAILFROM="adm_timo <timo-admin@zib.de>"
EMAILTO="adm_timo <timo-admin@zib.de>"

# SCIP check files are check.TESTSET.VERSION.otherstuff.SETTING.{out,err,res,meta}
BASEFILE="check/results/check.$TESTSET.*.$SETTING"

# evaluate the run and upload it to rubberband
echo "Evaluating the run and uploading it to rubberband."
cd check/
PERF_MAIL=""
if [ ${PERFORMANCE} == 'performance' ]; then
  ./evalcheck_cluster.sh -R results/check.$TESTSET.*.$SETTING[.0-9]*eval >> ${OUTPUT}
  NEWRBID=`echo $OUTPUT | grep "rubberband.zib" |sed -e 's|https://rubberband.zib.de/result/||'`
  OLDRBID=`tail $RBDB -n 1`
  PERF_MAIL=`echo "The results of the weekly performance runs are ready. Take a look at https://rubberband.zib.de/result/${NEWRBID}?compare=${OLDRBID}"`
  echo $NEWRBID >> $RBDB
  rm ${OUTPUT}
else
  ./evalcheck_cluster.sh -R results/check.$TESTSET.*.$SETTING[.0-9]*eval
fi
cd ..

# construct string which shows the destination of the out, err, and res files
SCIPDIR=`pwd`
ERRORFILE=`ls $BASEFILE[.0-9]*err`
OUTFILE=`ls $BASEFILE[.0-9]*out`
RESFILE=`ls $BASEFILE[.0-9]*res`
DESTINATION="$SCIPDIR/$OUTFILE \n$SCIPDIR/$ERRORFILE \n$SCIPDIR/$RESFILE"

# check for fixed instances
echo "Checking for fixed instances."
RESOLVEDINSTANCES=`awk '
NR != FNR {
  if( $3 == "'$GITBRANCH'" && $4 == "'$TESTSET'" && $5 == "'$SETTING'" && $6 == "'$OPT'" && $7 == "'$LPS'")
  {
     if( $0 in bugs == 0 )
        print "Previously failing instance " $1 " with error " $2 " does not fail anymore"
     else
        print $0 >> "'$TMPDATABASE'"
  }
  else
     print $0 >> "'$TMPDATABASE'"
  next;
}
## fail instances for this configuration
/fail/ {
  failmsg=$13; for(i=14;i<=NF;i++){failmsg=failmsg"_"$i;}
  errorstring=$1 " " failmsg " '$GITBRANCH' '$TESTSET' '$SETTING' '$OPT' '$LPS'";
  bugs[errorstring]
}' $RESFILE $DATABASE`
mv $TMPDATABASE $DATABASE

# check if fail occurs
NFAILS=`grep -c fail $RESFILE`

# if there are fails check for new fails and send email with information if needed
if [ $NFAILS -gt 0 ]; then
  echo "Detected ${NFAILS} fails."
  ERRORINSTANCES=`awk '
  ## read all known bugs
  NR == FNR {known_bugs[$0]; next}
  /fail/ {
     ## get the fail error and build string with the format of "database"
     failmsg=$13; for(i=14;i<=NF;i++){failmsg=failmsg"_"$i;}
     errorstring=$1 " " failmsg " '$GITBRANCH' '$TESTSET' '$SETTING' '$OPT' '$LPS'";
     ## if error is not in "database", add it and print it in ERRORINSTANCES to send email
     if( errorstring in known_bugs == 0 )
     {
        print errorstring >> "'$DATABASE'";
        print $0;
     }
     else # these are instances that failed before
     {
        print $1 " " failmsg >> "'$STILLFAILING'"; # only report the name of the instance and the fail message
     }
  }' $DATABASE $RESFILE`
  STILLFAILINGDB=`cat ${STILLFAILING}`

  # check if there are errors (string non empty)
  if [ -n "$ERRORINSTANCES" ]; then
     echo "Found new errors, sending emails."
     SUBJECT="FAIL [BRANCH: $GITBRANCH] [TESTSET: $TESTSET] [SETTING=$SETTING] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH]"
     echo -e "There are newly failed instances:\n${ERRORINSTANCES}\n\nThe follwing instances are still failing:\n${STILLFAILINGDB}\n\nThe files can be found here:\n$DESTINATION\n\nPlease note that they might be deleted soon" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
  else
     echo "No new errors, sending no emails."
  fi
else
  echo "No fails detected."
fi

# send email if there are fixed instances
if [ -n "$RESOLVEDINSTANCES" ]; then
   SUBJECT="FIX [BRANCH: $GITBRANCH] [TESTSET: $TESTSET] [SETTING=$SETTING] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH]"
   echo -e "Congratulations!\n\nThe following errors have been fixed:\n${RESOLVEDINSTANCES}\n\nThe following instances are still failing:\n${STILLFAILINGDB}\n\n" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
fi
rm ${STILLFAILING}

if [ ${PERFORMANCE} == 'performance' ]; then
   SUBJECT="WEEKLYPERF [BRANCH: $GITBRANCH] [TESTSET: $TESTSET] [SETTING=$SETTING] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH]"
   echo -e "${PERF_MAIL}" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
fi
