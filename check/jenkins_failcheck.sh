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
DATABASE="/nfs/OPTI/adm_timo/databases/${PWD##*/}_${TESTSET}_${SETTING}_${LPS}.txt"
TMPDATABASE="$DATABASE.tmp"
STILLFAILING="${DATABASE}_SF.tmp"
RBDB="/nfs/OPTI/adm_timo/databases/rbdb/${PWD##*/}_${TESTSET}_${SETTING}_${LPS}_rb.txt"
OUTPUT="${DATABASE}_output.tmp"
touch ${STILLFAILING}

# the first time, the file might not exists so we create it
# Even more, we have to write something to it, since otherwise
# the awk scripts below won't work (NR and FNR will not be different)
echo "Preparing database."
if ! [[ -s $DATABASE ]]; then  # check that file exists and has size larger that 0
  echo "Instance Fail_reason Branch Testset Setting Opt_mode LPS" > $DATABASE
fi

EMAILFROM="adm_timo <timo-admin@zib.de>"
EMAILTO="adm_timo <timo-admin@zib.de>"

########################
# FIND evalfile prefix #
########################

# SCIP check files are check.TESTSET.SCIPVERSION.otherstuff.SETTING.{out,err,res,meta} (SCIPVERSION is of the form scip-VERSION)
BASEFILE="check/results/check.${TESTSET}.${SCIPVERSION}.*.${SETTING}"
EVALFILE=`ls ${BASEFILE}*eval`
# if no evalfile was found --> check if this is fscip output
if [ "${EVALFILE}" == "" ]; then
    echo "Ignore previous ls error; looking again for eval file"
    BASEFILE="check/results/check.${TESTSET}.fscip.*.${SETTING}"
    EVALFILE=`ls ${BASEFILE}*eval`
fi

# if still no evalfile was found --> send an email informing that something is wrong and exit
if [ "${EVALFILE}" == "" ]; then
    echo "Couldn't find eval file, sending email"
    SUBJECT="ERROR [BRANCH: $GITBRANCH] [TESTSET: $TESTSET] [SETTING=$SETTING] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH]"
    echo -e "Aborting because the .eval file cannot be found.\nTried " | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
    exit 1
fi

# if more than one evalfile was found --> something is wrong, send an email
if [ `wc -w <<< ${EVALFILE}` -gt 1 ]; then
    echo "More than one eval file found; sending email"
    SUBJECT="ERROR [BRANCH: $GITBRANCH] [TESTSET: $TESTSET] [SETTING=$SETTING] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH]"
    echo -e "Aborting because there were more than one .eval files found:\n${EVALFILE}\n\nAfter fixing this run\ncd `pwd`\nPERFORMANCE=$PERFORMANCE SCIPVERSION=$SCIPVERSION SETTING=$SETTING LPS=$LPS GITHASH=$GITHASH OPT=$OPT TESTSET=$TESTSET GITBRANCH=$GITBRANCH ./check/jenkins_failcheck.sh\n" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
    exit 1
fi

############################################
# Process evalfile and upload to ruberband #
############################################

# evaluate the run and upload it to rubberband
echo "Evaluating the run and uploading it to rubberband."
cd check/
PERF_MAIL=""
if [ "${PERFORMANCE}" == "performance" ]; then
  ./evalcheck_cluster.sh -R ../${EVALFILE} > ${OUTPUT}
  cat ${OUTPUT}
  NEWRBID=`cat $OUTPUT | grep "rubberband.zib" |sed -e 's|https://rubberband.zib.de/result/||'`
  OLDRBID=`tail $RBDB -n 1`
  PERF_MAIL=`echo "The results of the weekly performance runs are ready. Take a look at https://rubberband.zib.de/result/${NEWRBID}?compare=${OLDRBID}"`
  echo $NEWRBID >> $RBDB
  rm ${OUTPUT}
else
  ./evalcheck_cluster.sh -T ../${EVALFILE}
fi
cd ..

# Store paths of err out res and set file
SCIPDIR=`pwd`
ERRFILE=SCIPDIR/`ls $BASEFILE*err`
OUTFILE=SCIPDIR/`ls $BASEFILE*out`
RESFILE=SCIPDIR/`ls $BASEFILE*res`
SETFILE=SCIPDIR/`ls $BASEFILE*set`

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


###################
# Check for fails #
###################

# if there are fails; process them and send email when there are new ones
NFAILS=`grep -c fail $RESFILE`
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

  # check if there are new fails!
  if [ -n "$ERRORINSTANCES" ]; then
      ###################
      ## Process fails ##
      ###################

      # get SCIP's header
      SCIP_HEADER=`awk '
      BEGIN{printLines=0;}

      /^SCIP version/ {printLines=1;}
      printLines > 0 && /^$/ {printLines+=1;}
      printLines > 0 {print $0}
      {
          if ( printLines == 3 ){
              exit 0;
          }
      }' $OUTFILE`

      # Get assertions and instance where they were generated
      ERRORS_INFO=`echo "${ERRORINSTANCES}" | awk '
      BEGIN{searchAssert=0; idx=1; delete failed[0]; delete human[0]}

      NR==FNR && /fail \(abort\)/ {failed[length(failed)+1]=$1; next;}

      # searching for an assert, but we see different instance name -> could not find assert for failed[idx],
      # continue with next failed (search for instance output, not assert!
      # failed[idx]"\\." --> name of the failed instance followed by a "."; this because the scripts keep the last characters of the name!
      searchAssert == 1 && /^@01/ && !($0 ~ failed[idx]"\\.") { human[failed[idx]]=1; searchAssert=0; idx+=1;}
      searchAssert == 1 && /Assertion.*failed\.$/ {print instance "\n" $0 "\n"; searchAssert=0; idx+=1;}

      /^@01/ && $0 ~ failed[idx] {searchAssert=1; instance=$0}
      END{ if( length(human) > 0 ){ print "The following fails need human inspection:"; for(key in human){ print key} }}
      ' - $ERRFILE`

      ###############
      # ERROR EMAIL #
      ###############
      echo "Found new errors, sending emails."
      SUBJECT="FAIL [BRANCH: $GITBRANCH] [TESTSET: $TESTSET] [SETTING=$SETTING] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH]"
      echo -e "There are newly failed instances.
      The instances run with the following SCIP version and setting file:
      ${SCIP_HEADER}
      ${SETFILE}

      Here is a list of the instances and the assertion that fails, if any.
      ${ERRORS_INFO}

      Here is the complete list of new fails:
      ${ERRORINSTANCES}

      The follwing instances are still failing:
      ${STILLFAILINGDB}

      Finally, the err, out and res file can be found here:
      $ERRFILE
      $OUTFILE
      $RESFILE

      Please note that they might be deleted soon" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
  else
      echo "No new errors, sending no emails."
  fi
else
  echo "No fails detected."
fi

# send email if there are fixed instances
if [ -n "$RESOLVEDINSTANCES" ]; then
   SUBJECT="FIX [BRANCH: $GITBRANCH] [TESTSET: $TESTSET] [SETTING=$SETTING] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH]"
   echo -e "Congratulations, see bottom for fixed instances!
   The following instances are still failing:
   ${STILLFAILINGDB}
   The err, out and res file can be found here:
   $ERRFILE
   $OUTFILE
   $RESFILE
   The following errors have been fixed:
   ${RESOLVEDINSTANCES}" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
fi
rm ${STILLFAILING}

if [ "${PERFORMANCE}" == "performance" ]; then
   SUBJECT="WEEKLYPERF [BRANCH: $GITBRANCH] [TESTSET: $TESTSET] [SETTING=$SETTING] [OPT=$OPT] [LPS=$LPS] [GITHASH: $GITHASH]"
   echo -e "${PERF_MAIL}" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
fi
