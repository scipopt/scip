#!/bin/bash -x
#
# This script is not meant to be used directly, but to be called.
# Usage: see clusterbench.sh
#
# This script glues the benchmark testrun files and sends an email if errors are detected.
#
# Careful: WE ARE IN check/ !


# reminder from clusterbench.sh: OUTPUTDIR="${CB_OUTPUTDIR}/${Q}_${n}"

CB_BASENAME="clusterbench_${CB_ID}"

# prepare database
# last benchmarkrun is in database
RBDB="/nfs/OPTI/adm_timo/databases/rbdb/clusterbenchmark_rb.txt"
touch ${RBDB}

cd ${CB_OUTPUTDIR}

for Q in $(echo ${QUEUE} | sed -e 's/,/ /g'); do
  QBASENAME="${Q}_${CB_BASENAME}"
  QDIRS="${Q}/*"
  # glue the .err and .out files
  #   the .err file is just the concatenation of individual .err files
  cat ${QDIRS}/check*.err > ${QBASENAME}.err
  #   the .out file is the concatenation of individual .meta and .out files
  #   where the .meta file comes directly before its respective .out file
  ls ${QDIRS}/check*.{meta,out} | sort | xargs cat > ${QBASENAME}.out

  # copy one settings file for rubberband upload, as all should be created equal by clusterbench.sh
  ls ${QDIRS}/check*.set | head -n 1 | xargs cat > ${QBASENAME}.set

  OUTPUT="jenkins_collect_benchmark_${CB_ID}.tmp"

  # upload to rubberband, and collect rubberbandid
  rbcli up ${QBASENAME}.* > ${OUTPUT}
  RBID=`cat $OUTPUT | grep "rubberband.zib" |sed -e 's|https://rubberband.zib.de/result/||'`
  # only save id if run by adm_timo
  if [ "${USER}" == "adm_timo" ]; then
    echo "${CB_ID} ${RBID} queue=${Q}" >> ${RBDB}
  else
    echo "Queue=${Q}: https://rubberband.zib.de/result/${RBID}"
  fi

  # clean up
  rm ${OUTPUT}
done

# if this is run by adm_timo then we want to save the rubberband id in the database
# and send an email
# otherwise just print the rubberband url in slurmlog
if [ "${USER}" == "adm_timo" ]; then
  # get the last clusterbench id
  OLDCB_ID=$(grep -P "queue=$Q($| )" ${RBDB} |tail -n 2|head -n 1|cut -d ' ' -f 2)

  MAILTEXT="The results of the clusterbenchmark are ready."

  # add a comparison for all queues
  for Q in $(echo ${QUEUE} | sed -e 's/,/ /g'); do
    LASTWEEK=$(grep -e ${OLDCB_ID} ${RBDB}|grep -P "queue=$Q($| )" |cut -d ' ' -f 2)
    THISWEEK=$(grep -e ${CB_ID} ${RBDB}|grep -P "queue=$Q($| )" |cut -d ' ' -f 2)
    if [ "${LASTWEEK}" != "" ]; then
      if [ "${THISWEEK}" != "" ]; then
        # construct rubberband link and mailtext
        MAILTEXT="${MAILTEXT}
Compare queue ${Q}: https://rubberband.zib.de/result/${LASTWEEK}?compare=${THISWEEK}"
      fi
    fi
  done

  # send email with cluster results
  EMAILFROM="adm_timo <timo-admin@zib.de>"
  EMAILTO="<timo-admin@zib.de>"

  SUBJECT="CLUSTERBENCHMARK"
  echo -e "${MAILTEXT}" | mailx -s "${SUBJECT}" -r "${EMAILFROM}" ${EMAILTO}
fi

