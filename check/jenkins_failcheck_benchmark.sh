#!/bin/bash -x
#
# This script uploads and checks for fails in a SCIP run.
# Sends an email if errors are detected. Is not meant to be use directly,
# but to be called by jenkins_check_results_benchmark.sh.
# Note: TESTSET etc are read from the environment, see
# jenkins_check_results_benchmark.sh
# Careful: WE ARE IN check/ !

echo "This is jenkins_failcheck_benchmark.sh running."

sleep 5

# evaluate the run and upload it to rubberband
echo "Evaluating the runs and uploading them to rubberband."

# SCIP check files are check.clusterbench.SCIPVERSION.otherstuff.SETTING.{out,err,res,meta} (SCIPVERSION is of the form scip-VERSION)
EVALFILES=`ls ${OUTPUTDIR}*.eval`
OUTPUT="clusterbenchmark_output.tmp"

# declare the array
declare -A RBIDS

# evaluate with evalcheck_cluster.sh, that also uploads to rubberband, and collect rubberbandids
i=0
for EVALFILE in ${EVALFILES};
do
  ./evalcheck_cluster.sh -R ${EVALFILE} > ${OUTPUT}
  cat ${OUTPUT}
  RBID=`cat $OUTPUT | grep "rubberband.zib" |sed -e 's|https://rubberband.zib.de/result/||'`

  # save the rubberbandid only if RBID is not empty
  if [ "${RBID}" != "" ]; then
    RBIDS[$i]=${RBID}
    ((i++))
  fi
done

rm ${OUTPUT}

# last benchmarkrun is in database
RBDB="/nfs/OPTI/adm_timo/databases/rbdb/clusterbenchmark_rb.txt"
touch $RBDB
LASTRUN=`tail $RBDB -n 1`

# construct rubberband link and mailtext
IDSTR=$(printf ",%s" "${RBIDS[@]:1}")
IDSTR=${IDSTR:1}
IDSTR="${RBIDS[0]}?compare=${IDSTR}${LASTRUN}"

# save comma seperated string of ids in database
SVESTR=$(printf ",%s" "${RBIDS[@]}")
echo $SVESTR >> $RBDB

MAILTEXT="The results of the clusterbenchmark are ready. Take a look at https://rubberband.zib.de/result/${IDSTR}"

# send email with cluster results
EMAILFROM="adm_timo <timo-admin@zib.de>"
EMAILTO="adm_timo <timo-admin@zib.de>"

SUBJECT="CLUSTERBENCHMARK"
echo -e "${MAILTEXT}" | mailx -s "$SUBJECT" -r "$EMAILFROM" $EMAILTO
