#!/bin/bash
#
# This script uploads and checks for fails in a SCIP run.
# Sends an email if errors are detected.
# Note: TESTSET, GITHAS, etc are read from the environment, see
# jenkins_check_results.sh

sleep 5

EMAILFROM="adm_timo <timo-admin@zib.de>"
EMAILTO="adm_timo <timo-admin@zib.de>"

BASEFILE="check/results/check.$TESTSET"

# evaluate the run and upload it to rubberband
cd check/
./evalcheck_cluster.sh -R results/check.$TESTSET.*.eval
cd ..

# check if fail occurs
NFAILS=`grep -c fail $BASEFILE.*res`

# construct string which shows the destination of the out, err, and res files
SCIPDIR=`pwd`
ERRORFILE=`ls $BASEFILE.*.err`
OUTFILE=`ls $BASEFILE.*.out`
RESFILE=`ls $BASEFILE.*.res`
DESTINATION="$SCIPDIR/$OUTFILE \n$SCIPDIR/$ERRORFILE \n$SCIPDIR/$RESFILE"

# if there are fails send email with information
if [ $NFAILS -gt 0 ];
then
  SUBJECT="FAIL [BRANCH: $GITBRANCH] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH]"
  ERRORINSTANCES=`grep fail $BASEFILE.*.res`
  echo -e "$ERRORINSTANCES \n\nThe files can be found here:\n$DESTINATION\n\nPlease note that the files might be deleted soon" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
fi
