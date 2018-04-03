#!/bin/bash -x
#
# This script uploads and checks for fails in a SCIP run.
# Sends an email if errors are detected. Is not meant to be use directly,
# but to be called by jenkins_check_results.sh.
# Note: TESTSET, GITHASH, etc are read from the environment, see
# jenkins_check_results.sh

echo "This is jenkins_failcheck.sh running."

sleep 5

#################################
# AWK scripts for later         #
#################################
read -d '' awkscript_findasserts << 'EOF'
# call: awk "$awkscript_findasserts" $ERRORINSTANCES $ERRFILE (or pipe a string with the error instances into awk if they are not in a file)

# init variables
BEGIN {
    searchAssert=0;
}

# read failed instances into array as keys
# the names in the errorinstances have to be the full names of the instances
# set failed[instname] to
#     1 if instance found in err file once, and no assert found
#     2 if instance found in err file once, also found assert
#     3 if instance found more than once in errfile
NR==FNR && /fail.*abort/ {
    failed[$1]=0; next;
}

# find instances in errorfile
NR!=FNR && /^@01/ {

    # if we were looking for an assertion in currinstname, we are now at the beginning of the error output
    # of another instance. Therefore, we didn't find an assertion and this instance needs human inspection.
    if (searchAssert == 1) { searchAssert=0; }

    # get instancename (copied from check.awk)
    n = split($2, a, "/"); m = split(a[n], b, "."); currinstname = b[1];
    if( b[m] == "gz" || b[m] == "z" || b[m] == "GZ" || b[m] == "Z" ) { m--; }
    for( i = 2; i < m; ++i ) { currinstname = currinstname "." b[i]; }

    instancestr = $2;

    # adjust searchAssert
    if (currinstname in failed) {
        searchAssert=1;
        if ( failed[currinstname] == 0 ) {
            failed[currinstname]=1
        } else {
            failed[currinstname]=3
        }
    }
}

# find assertions in errorfile
NR!=FNR && searchAssert == 1 && /Assertion.*failed.$/ {
    print "";
    print instancestr
    for(i=2;i<=NF;i++){printf "%s ", $i}; print "";
    if ( failed[currinstname] < 3 ) {
        failed[currinstname]=2
    } else {
        failed[currinstname]=3
    }
    searchAssert=0;
}

# print results
END {
    print "";
    count = 0;
    for( key in failed ) {
        if( failed[key] == 1 ) {
            if( count == 0 ) {
                print "The following fails need human inspection, because there was no Assertion found:";
            }
            print key;
            count = count + 1;
        }
    }
    count = 0;
    for( key in failed ) {
        if( failed[key] == 3 ) {
            if( count == 0 ) {
                print "The following instances do not appear in the .err file, or some other error appeared:";
            }
            print key;
            count = count + 1;
        }
    }
}
EOF

read -d '' awkscript_checkfixedinstances << 'EOF'
# call: awk $AWKARGS "$awkscript_checkfixedinstances" $RESFILE $DATABASE
# prints fixed instances
# afterwards $TMPDATABASE contains all still failing bugs, not new ones!

# read fail instances for this configuration from resfile
NR == FNR && /fail/ {
    failmsg=$13; for(i=14;i<=NF;i++){ failmsg=failmsg"_"$i; }
    errorstring=$1 " " failmsg " " GITBRANCH " " TESTSET " " SETTING " " OPT " " LPS " " PERM;
    bugs[errorstring]
    next;
}

# read from database
NR != FNR {
    if( $3 == GITBRANCH && $4 == TESTSET && $5 == SETTING && $6 == OPT && $7 == LPS ) {
        if (!( $0 in bugs )) {
            # if current line from database matches our settings and
            # it is not in the set of failed instances from this run it was fixed
            print "Previously failing instance " $1 " with error " $2 " does not fail anymore"
            next;
        }
    }
    # write record into the database for next time
    print $0 >> TMPDATABASE
}
EOF

read -d '' awkscript_readknownbugs << 'EOF'
# call: awk $AWKARGS "$awkscript_readknownbugs" $DATABASE $RESFILE
# append new fails to DATABASE, also print them
# write all still failing bugs into STILLFAILING

# read known bugs from database
NR == FNR {known_bugs[$0]; next}

# find fails in resfile
/fail/ {
    # get the fail error and build string in database format
    failmsg=$13; for(i=14;i<=NF;i++){failmsg=failmsg"_"$i;}
    errorstring=$1 " " failmsg " " GITBRANCH " " TESTSET " " SETTING " " OPT " " LPS " " PERM;

    if (!( errorstring in known_bugs )) {
        # if error is not known, add it and print it to ERRORINSTANCES for sending mail later
        print errorstring >> DATABASE;
        print $0;
    } else {
        # if error is known, then instance failed before with same settings
        # only report the name and the fail message of the instance
        print $1 " " failmsg >> STILLFAILING;
    }
}
EOF

read -d '' awkscript_scipheader << 'EOF'
# call: awk "$awkscript_scipheader" $OUTFILE
# prints current scipheader from OURFILE

BEGIN{printLines=0;}

/^SCIP version/ {printLines=1;}
printLines > 0 && /^$/ {printLines+=1;}

# UG does not have blank lines after the header.
printLines > 0 && /^Default LC presolving/ {exit 0;}

# the third blank line marks the end of the header
printLines > 0 {print $0}
{
    if ( printLines == 3 ){
        exit 0;
    }
}
EOF
#################################
# End of AWK Scripts            #
#################################

# The RBDB database has the form: timestamp_of_testrun rubberbandid p=PERM s=SEED
RBDB="/nfs/OPTI/adm_timo/databases/rbdb/${PWD##*/}_${TESTSET}_${SETTING}_${LPS}_rbdb.txt"
touch $RBDB
OLDTIMESTAMP=`tail -n 1 ${RBDB}|cut -d ' ' -f 1`
NEWTIMESTAMP=`date '+%F-%H-%M'`

PERM=0
while [ $PERM -le $PERMUTE ]; do
  # get ending given by permutation
  if [ "${PERM}" == "0" ]; then
      PERM_ENDING="."
  else
      PERM_ENDING="-p${PERM}."
  fi

  # we use a name that is unique per test sent to the cluster (a jenkins job
  # can have several tests sent to the cluster, that is why the jenkins job
  # name (i.e, the directory name) is not enough)
  DATABASE="/nfs/OPTI/adm_timo/databases/${PWD##*/}_${TESTSET}_${SETTING}_${LPS}${PERM_ENDING}txt"
  TMPDATABASE="$DATABASE.tmp"
  STILLFAILING="${DATABASE}_SF.tmp"
  OUTPUT="${DATABASE}_output.tmp"
  touch ${STILLFAILING}

  SUBJECTINFO="[BRANCH: $GITBRANCH] [TESTSET: $TESTSET] [SETTING: $SETTING] [OPT: $OPT] [LPS: $LPS] [GITHASH: $GITHASH] [PERM: $PERM]"

  AWKARGS="-v GITBRANCH=$GITBRANCH -v TESTSET=$TESTSET -v SETTING=$SETTING -v OPT=$OPT -v LPS=$LPS -v DATABASE=$DATABASE -v TMPDATABASE=$TMPDATABASE -v STILLFAILING=$STILLFAILING -v PERM=$PERM"
  echo $AWKARGS

  # the first time, the file might not exists so we create it
  # Even more, we have to write something to it, since otherwise
  # the awk scripts below won't work (NR and FNR will not be different)
  echo "Preparing database."
  if ! [[ -s $DATABASE ]]; then  # check that file exists and has size larger that 0
    echo "Instance Fail_reason Branch Testset Setting Opt_mode LPS" > $DATABASE
  fi

  EMAILFROM="adm_timo <timo-admin@zib.de>"
  EMAILTO="adm_timo <timo-admin@zib.de>"

  #################
  # FIND evalfile #
  #################

  # SCIP check files are check.TESTSET.SCIPVERSION.otherstuff.SETTING.{out,err,res,meta} (SCIPVERSION is of the form scip-VERSION)
  BASEFILE="check/results/check.${TESTSET}.${SCIPVERSION}.*.${SETTING}${PERM_ENDING}"
  EVALFILE=`ls ${BASEFILE}*eval`
  # if no evalfile was found --> check if this is fscip output
  if [ "${EVALFILE}" == "" ]; then
      echo "Ignore previous ls error; looking again for eval file"
      BASEFILE="check/results/check.${TESTSET}.fscip.*.${SETTING}" # we do not use permutations with fiber scip
      EVALFILE=`ls ${BASEFILE}*eval`
  fi

  # if still no evalfile was found --> send an email informing that something is wrong and exit
  if [ "${EVALFILE}" == "" ]; then
      echo "Couldn't find eval file, sending email"
      SUBJECT="ERROR ${SUBJECTINFO}"
      echo -e "Aborting because the .eval file cannot be found.\n\nTried:\n${BASEFILE}.eval\ncheck/results/check.${TESTSET}.${SCIPVERSION}.*.${SETTING}${PERM_ENDING}*eval\n\nDirectory: `pwd`.\n" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
      exit 1
  fi

  # if more than one evalfile was found --> something is wrong, send an email
  if [ `wc -w <<< ${EVALFILE}` -gt 1 ]; then
      echo "More than one eval file found; sending email"
      SUBJECT="ERROR ${SUBJECTINFO}"
      echo -e "Aborting because there were more than one .eval files found:\n${EVALFILE}\n\nAfter fixing this run\ncd `pwd`\nPERFORMANCE=$PERFORMANCE SCIPVERSION=$SCIPVERSION SETTING=$SETTING LPS=$LPS GITHASH=$GITHASH OPT=$OPT TESTSET=$TESTSET GITBRANCH=$GITBRANCH PERM=$PERM ./check/jenkins_failcheck.sh\n" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
      exit 1
  fi

  # at this point we have exactly one evalfile
  BASENAME=${EVALFILE%.*} # remove extension
  WORKINGDIR=`pwd -P`

  # Store paths of err out res and set file
  ERRFILE="${WORKINGDIR}/${BASENAME}.err"
  OUTFILE="${WORKINGDIR}/${BASENAME}.out"
  RESFILE="${WORKINGDIR}/${BASENAME}.res"
  SETFILE="${WORKINGDIR}/${BASENAME}.set"
  EVALFILE="${WORKINGDIR}/${BASENAME}.eval"

  ############################################
  # Process evalfile and upload to ruberband #
  ############################################

  # evaluate the run and upload it to rubberband
  echo "Evaluating the run and uploading it to rubberband."
  cd check/
  PERF_MAIL=""
  if [ "${PERFORMANCE}" == "performance" ]; then
    ./evalcheck_cluster.sh -R ${EVALFILE} > ${OUTPUT}
    NEWRBID=`cat $OUTPUT | grep "rubberband.zib" |sed -e 's|https://rubberband.zib.de/result/||'`
    echo "${NEWTIMESTAMP} ${NEWRBID} p=${PERM}" >> $RBDB
  else
    ./evalcheck_cluster.sh -r "-v useshortnames=0" ${EVALFILE} > ${OUTPUT}
  fi
  cat ${OUTPUT}
  rm ${OUTPUT}
  cd ..

  # check for fixed instances
  echo "Checking for fixed instances."
  RESOLVEDINSTANCES=`awk $AWKARGS "$awkscript_checkfixedinstances" $RESFILE $DATABASE`
  echo "Temporary database: $TMPDATABASE\n"
  mv $TMPDATABASE $DATABASE


  ###################
  # Check for fails #
  ###################

  # if there are fails; process them and send email when there are new ones
  NFAILS=`grep -c fail $RESFILE`
  if [ $NFAILS -gt 0 ]; then
    echo "Detected ${NFAILS} fails."
    ## read all known bugs
    ERRORINSTANCES=`awk $AWKARGS "$awkscript_readknownbugs" $DATABASE $RESFILE`
    STILLFAILINGDB=`cat ${STILLFAILING}`

    # check if there are new fails!
    if [ -n "$ERRORINSTANCES" ]; then
        ###################
        ## Process fails ##
        ###################

        # get SCIP's header
        SCIP_HEADER=`awk "$awkscript_scipheader" $OUTFILE`

        if [ "${PERFORMANCE}" != "performance" ]; then
            # Get assertions and instance where they were generated
            ERRORS_INFO=`echo "${ERRORINSTANCES}" | awk "$awkscript_findasserts" - ${ERRFILE}`
        fi

        ###############
        # ERROR EMAIL #
        ###############
        echo "Found new errors, sending emails."
        SUBJECT="FAIL ${SUBJECTINFO}"
        echo -e "There are newly failed instances.
The instances run with the following SCIP version and setting file:

\`\`\`
BRANCH: $GITBRANCH

SCIP HEADER:
${SCIP_HEADER}

SETTINGS FILE:
${SETFILE}
\`\`\`

Here is a list of the instances and the assertion that fails (fails with _fail (abort)_), if any:
${ERRORS_INFO}

Here is the complete list of new fails:
${ERRORINSTANCES}

The following instances are still failing:
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
     #########################
     # RESOLVED ERRORS EMAIL #
     #########################
     SUBJECT="FIX ${SUBJECTINFO}"
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

  PERM=$((PERM + 1))
done

# construct the rubberband link
if [ "${PERFORMANCE}" == "performance" ]; then
  # collect all ids with timestamps OLDTIMESTAMP NEWTIMESTAMP in RBIDS
  RBDB_STRS=`grep -e "\(${OLDTIMESTAMP}\|${NEWTIMESTAMP}\)" ${RBDB}|cut -d ' ' -f 2`

  i=0
  while read -r line; do
      arr=($line)
      RBIDS[$i]=${arr[-1]}
      ((i++))
  done <<< "${RBDB_STRS}"

  IDSTR=$(printf ",%s" "${RBIDS[@]}")
  IDSTR=${IDSTR:1}

  URLSTR=`echo ${IDSTR} | sed 's/,/?compare=/'`

  PERF_MAIL=`echo "The results of the weekly performance runs are ready. Take a look at https://rubberband.zib.de/result/${URLSTR}
"`

  PERM=0
  while [ $PERM -le $PERMUTE ]; do
    LASTWEEK=`grep -e ${OLDTIMESTAMP} ${RBDB}|grep p=$PERM|cut -d ' ' -f 2`
    THISWEEK=`grep -e ${NEWTIMESTAMP} ${RBDB}|grep p=$PERM|cut -d ' ' -f 2`
    if [ "${LASTWEEK}" != "" ]; then
      if [ "${THISWEEK}" != "" ]; then
        PERF_MAIL="${PERF_MAIL}
Compare permutation ${PERM}: https://rubberband.zib.de/result/${LASTWEEK}?compare=${THISWEEK}"
      fi
    fi
    PERM=$((PERM + 1))
  done

  SUBJECT="WEEKLYPERF ${SUBJECTINFO}"
  echo -e "$PERF_MAIL" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
fi
