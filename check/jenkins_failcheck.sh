#!/bin/bash
#
# This script uploads and checks for fails in a SCIP run.
# Sends an email if errors are detected.
# Note: TESTSET, GITHAS, etc are read from the environment, see
# jenkins_check_results.sh

sleep 5

DATABASE="/nfs/OPTI/bzfserra/jenkins/known_bugs.txt"

EMAILFROM="adm_timo <timo-admin@zib.de>"
EMAILTO="adm_timo <timo-admin@zib.de>"

# SCIP check files are check.TESTSET.VERSION.otherstuff.SETTING.{out,err,res}
BASEFILE="check/results/check.$TESTSET.*.$SETTING"

# evaluate the run and upload it to rubberband
cd check/
./evalcheck_cluster.sh -R results/check.$TESTSET.*.$SETTING.eval
cd ..

# check if fail occurs
NFAILS=`grep -c fail $BASEFILE.res`

# construct string which shows the destination of the out, err, and res files
SCIPDIR=`pwd`
ERRORFILE=`ls $BASEFILE.err`
OUTFILE=`ls $BASEFILE.out`
RESFILE=`ls $BASEFILE.res`
DESTINATION="$SCIPDIR/$OUTFILE \n$SCIPDIR/$ERRORFILE \n$SCIPDIR/$RESFILE"

# if there are fails check for new fails and send email with information if needed
if [ $NFAILS -gt 0 ]; then
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
  }' $DATABASE $BASEFILE.res`

  # check if there are errors (string non empty)
  if [ -n "$ERRORINSTANCES" ]; then
     SUBJECT="FAIL [BRANCH: $GITBRANCH] [TESTSET: $TESTSET] [SETTING=$SETTING] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH]"
     echo -e "$ERRORINSTANCES \n\nThe files can be found here:\n$DESTINATION\n\nPlease note that the files might be deleted soon" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
  fi
fi
