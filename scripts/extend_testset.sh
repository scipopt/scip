#!/bin/bash  -e

# usage:
#    TESTSETFILE=awesome.test CHECKFILESBASE=check.scip bash extend_testset.sh

# CHECKFILESBASE : basename of testrunfiles
# TESTSETFILE    : testsetfile

echo "Scanning ${CHECKFILESBASE}.* for instances to add to ${TESTSETFILE}."

NEWTESTSETFILE=${TESTSETFILE}_merged_$(date "+%Y%m%d%H%M%S")
cp ${TESTSETFILE} ${NEWTESTSETFILE}

GOODINSTANCES=$(awk '{if (NR > 3) { if (substr($1,1,1) == "-") { exit; }; if (($13 == "ok") || ($13 == "gaplimit") || ($13 == "better") || ($13 == "solved")) { print $1; }}}' ${CHECKFILESBASE}.res)
PATHS=$(grep @01 ${CHECKFILESBASE}.err)
for i in ${GOODINSTANCES}; do
  INSTANCEPATH=$(echo "${PATHS}" | grep "/$i\." | cut -d " " -f 2 | grep -o "/check/.*" | cut -c 8- )
  TEST=$(grep "${INSTANCEPATH}" ${TESTSETFILE})
  if [ "${TEST}" = "" ]; then
    echo "Adding ${INSTANCEPATH} to testsetfile."
    echo "${INSTANCEPATH}" >> ${NEWTESTSETFILE}
  fi
done
cat ${NEWTESTSETFILE} | sort > ${NEWTESTSETFILE}.tmp
mv ${NEWTESTSETFILE}.tmp ${NEWTESTSETFILE}
echo "Done, you can find you new testsetfile here '${NEWTESTSETFILE}'."
