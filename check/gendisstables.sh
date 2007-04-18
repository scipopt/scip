#!/bin/sh

TEXINCFILE="disstables/detailedtables.tex"

# settings: (primary class)/subgroup
SETTINGS="mem/default branch/default brscore/default nodesel/default childsel/default prop/default sepa/default sepa/sepaonly sepa/sepalocal cutsel/default heur/default heur/heurrounding heur/heurdiving heur/heurobjdiving heur/heurimprovement presol/default conf/default"

# test sets: set/(weight in total score), weight=0 means weight=#instances
TESTSETS="miplib/0 coral/0 milp/0 enlight/3 alu/3 fctp/3 acc/3 fc/3 arcset/3 mik/3 cls/3"

rm -f gendisstables.out

# count valid test sets
NTESTSETS=0
for testsetweight in $TESTSETS
do
  testset=`dirname $testsetweight`

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

  TEXSUMMARYBASE="disstables/Table_${setname}_${groupname}_summary"
  TEXSUMMARYFILE="${TEXSUMMARYBASE}.tex"
  rm -f $TEXSUMMARYFILE
  SUMMARYHEADER=$NTESTSETS

  echo "" >> $TEXINCFILE
  echo "\\header${setname}${groupname}" >> $TEXINCFILE
  echo "" >> $TEXINCFILE

  for testsetweight in $TESTSETS
  do
    testset=`dirname $testsetweight`
    weight=`basename $testsetweight`

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
	    disscmpres.sh onlygroup=$groupname texincfile=$TEXINCFILE texfile="disstables/Table_${setname}_${groupname}_${testset}.tex" texsummaryfile="$TEXSUMMARYFILE" texsummaryheader=$SUMMARYHEADER texsummaryweight=$weight textestset="$testset" results/check.diss_${testsetfile}*.*.default.res results/check.diss_${testsetfile}*.*.${setname}_*.res >> gendisstables.out
	fi
    fi
    echo "" >> $TEXINCFILE
    echo "\\clearpage${setname}${groupname}${testset}" >> $TEXINCFILE
    echo "" >> $TEXINCFILE
    SUMMARYHEADER=0
  done

  # generate "TOTAL" row in summary file
  gawk -v texsummaryfiletime="${TEXSUMMARYBASE}_time.tex" -v texsummaryfilenodes="${TEXSUMMARYBASE}_nodes.tex" -v texcolorlimit=5 '
function abs(x)
{
   return x < 0 ? -x : x;
}
function floor(x)
{
   return (x == int(x) ? x : (x < 0 ? int(x-1) : int(x)));
}
function texcompstr(val, x,s,t)
{
   x = floor(100*(val-1.0)+0.5);
   s = "";
   t = "";
   if( x < 0 )
   {
      if( x <= -texcolorlimit )
      {
         s = "\\textcolor{red}{\\raisebox{0.25ex}{\\tiny $-$}";
         t = "}";
      }
      else
         s = "\\raisebox{0.25ex}{\\tiny $-$}";
   }
   else if( x > 0 )
   {
      if( x >= +texcolorlimit )
      {
         s = "\\textcolor{blue}{\\raisebox{0.25ex}{\\tiny $+$}";
         t = "}";
      }
      else
         s = "\\raisebox{0.25ex}{\\tiny $+$}";
   }

   return sprintf("{\\bf%s%d%s}", s, abs(x), t);
}
BEGIN { nsolver = 0; }
/^% =geom=/ {
  if( !counted[$3] ) {
    nsolver++;
    counted[$3] = 1;
    solver[nsolver] = $3;
  }
  nvals[$3]++;
  timegeom[$3,nvals[$3]] = $4;
  nodesgeom[$3,nvals[$3]] = $5;
  weight[$3,nvals[$3]] = $6;
}
END {
  printf("\\addlinespace[0.4\\defaultaddspace]\n") >> texsummaryfiletime;
  printf("\\addlinespace[0.4\\defaultaddspace]\n") >> texsummaryfilenodes;
  printf("& \\textbf{total}") >> texsummaryfiletime;
  printf("& \\textbf{total}") >> texsummaryfilenodes;
  for( i = 1; i <= nsolver; i++ )
  {
    s = solver[i];
    tottimelogsum = 0.0;
    totnodeslogsum = 0.0;
    cnt = 0;
    for( j = 1; j <= nvals[s]; j++ )
    {
      tottimelogsum += log(timegeom[s,j]) * weight[s,j];
      totnodeslogsum += log(nodesgeom[s,j]) * weight[s,j];
      cnt += weight[s,j];
    }
    tottimegeom = exp(tottimelogsum/cnt);
    totnodesgeom = exp(totnodeslogsum/cnt);
    printf("& %s", texcompstr(tottimegeom)) >> texsummaryfiletime;
    printf("& %s", texcompstr(totnodesgeom)) >> texsummaryfilenodes;
  }
  printf("\\\\\n") >> texsummaryfiletime;
  printf("\\\\\n") >> texsummaryfilenodes;
  printf("\\addlinespace[0.4\\defaultaddspace]\n") >> texsummaryfiletime;
  printf("\\addlinespace[0.4\\defaultaddspace]\n") >> texsummaryfilenodes;
}' $TEXSUMMARYFILE
done
