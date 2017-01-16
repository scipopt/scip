#!/bin/bash
#
# This script checks for fails in a SCIP run and sends emails if errors are detected.
# Otherwise, it uploads to ruberband.

sleep 5

EMAILFROM="adm_timo <timo-admin@zib.de>"
EMAILTO="adm_timo <timo-admin@zib.de>"

BASEFILE="check/results/check.$TESTSET"

# evaluate the run
cd check/
./evalcheck_cluster.sh results/check.$TESTSET.*.eval
cd ..

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
  SUBJECT="[FAIL] [GITBRANCH: $GITBRANCH] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH]"
  ERRORINSTANCES=`grep fail $BASEFILE.*.res`
  echo -e "$ERRORINSTANCES \nThe files can be found here:\n$DESTINATION\n Please note that the files might be deleted soon" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
else
  echo -e "Uploading to ruberband"
  cd check/
  #./evalcheck_cluster.sh -R results/check.$TESTSET.*.eval
  cd ..
fi
