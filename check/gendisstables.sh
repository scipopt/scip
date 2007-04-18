#!/bin/sh

TEXINCFILE="disstables/detailedtables.tex"

SETTINGS="mem/default branch/default brscore/default nodesel/default childsel/default prop/default sepa/default sepa/sepaonly sepa/sepalocal cutsel/default heur/default heur/heurrounding heur/heurdiving heur/heurobjdiving heur/heurimprovement presol/default conf/default"
TESTSETS="miplib coral milp enlight alu fctp acc fc arcset mik cls"

rm -f gendisstables.out

# count valid test sets
NTESTSETS=0
for testset in $TESTSETS
do
  for atamtuerkres in results/check.diss_atamtuerk_${testset}*.*.default.res
  do
    break
  done
  for normalres in results/check.diss_${testset}*.*.default.res
  do
    break
  done
  
  if [ -f $atamtuerkres ]
  then
      NTESTSETS=$(($NTESTSETS+1))
  elif [ -f $normalres ]
  then
      NTESTSETS=$(($NTESTSETS+1))
  fi
done

echo processing $NTESTSETS test sets

echo "% tables generated automatically" > $TEXINCFILE
echo "" >> $TEXINCFILE

#echo "% setting names" >> $TEXINCFILE
#for setgroupname in $SETTINGS
#do
#  setname=`dirname $setgroupname`
#  echo "\\newcommand{\\setting${setname}}{${setname}\\xspace}" >> $TEXINCFILE
#done
#echo "" >> $TEXINCFILE

for setgroupname in $SETTINGS
do
  setname=`dirname $setgroupname`
  groupname=`basename $setgroupname`

  TEXSUMMARYFILE="disstables/Table_${setname}_${groupname}_summary.tex"
  rm -f $TEXSUMMARYFILE
  SUMMARYHEADER=$NTESTSETS

  echo "" >> $TEXINCFILE
  echo "\\header${setname}${groupname}" >> $TEXINCFILE
  echo "" >> $TEXINCFILE

  for testset in $TESTSETS
  do
    for atamtuerkres in results/check.diss_atamtuerk_${testset}*.*.default.res
    do
      break
    done
    for normalres in results/check.diss_${testset}*.*.default.res
    do
      break
    done

    if [ -f $atamtuerkres ]
    then
	ATAMTUERK="atamtuerk_"
    elif [ -f $normalres ]
    then
	ATAMTUERK=""
    fi
    testsetfile=${ATAMTUERK}${testset}

    echo "SET: $setname  GROUP: $groupname  TESTSET: $testset  TESTSETFILE: $testsetfile"

    for defaultres in results/check.diss_${testsetfile}*.*.default.res
    do
      break
    done
    for setres in results/check.diss_${testsetfile}*.*.${setname}_*.res
    do
      break
    done

    if [ -f $defaultres ]
    then
	if [ -f $setres ]
	then
	    disscmpres.sh onlygroup=$groupname texincfile=$TEXINCFILE texfile="disstables/Table_${setname}_${groupname}_${testset}.tex" texsummaryfile="$TEXSUMMARYFILE" texsummaryheader=$SUMMARYHEADER textestset="$testset" results/check.diss_${testsetfile}*.*.default.res results/check.diss_${testsetfile}*.*.${setname}_*.res >> gendisstables.out
	fi
    fi
    echo "" >> $TEXINCFILE
    echo "\\clearpage${setname}${groupname}${testset}" >> $TEXINCFILE
    echo "" >> $TEXINCFILE
    SUMMARYHEADER=0
  done
done
