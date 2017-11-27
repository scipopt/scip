#!/usr/bin/env bash

set -eu

declare -A OPTVAL
OPTVAL['magic']=988540
OPTVAL['alkyl']=-1.765011

MAINDIR=$(cd `dirname $0`; pwd)
GAMSDIR="${MAINDIR}/gams"

rm -rf test
mkdir -p test
cd test

cat > subsys.txt <<EOF
SCIPTEST 111 0 0001020304 1 0 2 MIP NLP CNS DNLP RMINLP MINLP QCP MIQCP RMIQCP
gmsgenus.run
gmsgenux.out
${MAINDIR}/$1 scp 1
EOF

for m in ${!OPTVAL[@]}
do
   echo "Running model $m"
   echo "* Trace Record Definition" > $m.trc
   echo "* GamsSolve" >> $m.trc
   echo "* ModelStatus SolverStatus ObjectiveValue ObjectiveValueEstimate" >> $m.trc

   "${GAMSDIR}/gamslib" -q $m
   "${GAMSDIR}/gams" $m subsys=subsys.txt solvelink=5 optcr=0 trace=$m.trc traceopt=3 lo=4
done

echo
echo "Checking outcome"
fail=0
for m in ${!OPTVAL[@]}
do
   result=(`tail -1 $m.trc`)
   echo "Model $m (optval ${OPTVAL[$m]}): model status ${result[0]}, solve status ${result[1]}, primal bound ${result[2]}, dual bound ${result[3]}"

   if [[ ${result[0]} != 1 ]] ; then
      echo "Unexpected model status ${result[0]}, should be 1"
      (( fail++ )) || true
   fi

   if [[ ${result[1]} != 1 ]] ; then
      echo "Unexpected solver status ${result[1]}, should be 1"
      (( fail++ )) || true
   fi

   diff=`echo "scale=20; (${result[2]} - (${OPTVAL[$m]}))^2 < 10^-12" | bc`
   if [[ ${diff} != 1 ]] ; then
      echo "Unexpected primal bound ${result[2]}, should be ${OPTVAL[$m]}"
      (( fail++ )) || true
   fi

   diff=`echo "scale=20; (${result[3]} - (${OPTVAL[$m]}))^2 < 10^-12" | bc`
   if [[ ${diff} != 1 ]] ; then
      echo "Unexpected dual bound ${result[3]}, should be ${OPTVAL[$m]}"
      (( fail++ )) || true
   fi

done

echo "Number of fails: $fail"
exit $fail
