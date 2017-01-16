#!/bin/bash
#
# This script checks for fails in a SCIP run and sends emails if errors are detected.
# Otherwise, it uploads to ruberband.
#
# TODO: maybe clean it up a little bit
# TODO: check whether there is more information that we would like to have.
# TODO: we can also mention where the outfiles are, but beware that jenkins deletes every day.

sleep 5

EMAILFROM="adm_timo <timo-admin@zib.de>"
EMAILTO="adm_timo <timo-admin@zib.de>"

BASEFILE="check/results/check.$TESTSET"

cd check/
./evalcheck_cluster.sh -r results/check.$TESTSET.*.eval
cd ..

# check if fail occurs
NFAILS=`grep -c fail $BASEFILE.*res`

# construct string which shows the destination of the out, err, and res files
ERRORFILE=`ls $BASEFILE.*.err`
OUTFILE=`ls $BASEFILE.*.out`
RESFILE=`ls $BASEFILE.*.res`
#DESTINATION="$SCIPDIR/$OUTFILE \n$SCIPDIR/$ERRORFILE \n$SCIPDIR/$RESFILE"

# check read fails
if [ $NFAILS -gt 0 ];
then
  SUBJECT="[FAIL] [] [OPT=] [LPS=] [GITHASH: ] "
  ERRORINSTANCES=`grep fail $BASEFILE.*.res`
  echo -e "$ERRORINSTANCES \n this is a test email" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
else
  #TODO: use ruberband
  echo -e "uploading to ruberband... not yet"
fi
