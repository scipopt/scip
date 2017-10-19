#!/bin/bash

# benchmarks the performance of all nodes of a given slurm queue: on each node run exclusively a predefined test set
# with a specified SCIP binary; to distinguish the runs, create artificial empty settings files such that comparison is
# easily possible with standard workflow

if test -z $USER
then
    echo $USER
    exit
fi

if [ ! "$(readlink -f `pwd`)" = "$(dirname $(readlink -f ${BASH_SOURCE[0]}))" ]; then
  echo 'cd to scip/check in order to run script'
  exit
fi;

if [ $# -ge 3 ] || [ $# -eq 0 ]; then
  echo 'usage: ./clusterbench.sh QUEUES EXECUTABLE'
  echo '       note: queues can be separated by comma'
  exit
fi;

QUEUE=$1
export EXECUTABLE=`pwd`/$2
DATETIME=`date '+%Y%m%d-%H%M%S'`

echo "queue      :" $QUEUE
echo "executable :" $EXECUTABLE
echo "date-time  :" $DATETIME
echo

# multiple queues are split by a comma (hopefully, no one ever attempts to define slurm partitions with a comma ;))
split_queue=(`echo $QUEUE | tr "," "\n"`)

declare -a opts
declare -a nodes

for Q in ${split_queue[@]}; do
  unset opts
  unset nodes

  echo $Q

  # evaluate sinfo for this partition and collect the nodes
  SINFO=$(`echo sinfo -p ${Q}` )
  LSINFO=(`echo ${SINFO[@]} | wc -w`)

  #echo SINFO is ${SINFO[@]}
  #echo LSINFO is $LSINFO
  SINFOSTR=(`echo ${SINFO[@]} | tr " " "\n"`)
  for i in $(seq 1 $LSINFO); do

    # no clue why this is needed, n is the current element
    n=(`echo ${SINFOSTR[i]} | sed 's/\*//'`)

    # modulo 6 because the table has six columns -> every first entry of a row gets checked
    if [[ $(($i % 6)) -eq 0 ]]; then
      if [ "$n" = "$Q" ] ; then

        # opts are listed as comma-separated list, split on comma and append to opts
        # we must not split simply on a comma, because this will spoil everything if
        # a comma is contained in brackets: optc-07-[01,05-16]
        # the sed-command replaces occurences of commas outside brackets ], ... [
        # by a white space
        moreopts=(`echo ${SINFOSTR[$i+5]} | sed -e 's/\(]\)\(,\)\([^[]*\)/\1 \3/g'`)
        echo "Nodes  :" ${moreopts[@]}
        opts=( ${opts[@]} ${moreopts[@]} );
      fi;
    fi;
  done

  # Expand compressed slurm array
  for nodelist in ${opts[@]}; do
    subnodelist=(`echo $nodelist | tr "[]" "\n"`)
    hostprefix=${subnodelist[0]}

    # my favorite one: single host ids without brackets; add the host id directly and continue
    if [ "${hostprefix}" = "${subnodelist[*]}" ];  then
      nodes=( ${nodes[*]} $hostprefix )
      continue
    fi

    # get the list content as array
    hostnumberlist=(`echo ${subnodelist[1]} | tr "," "\n"`)

    # loop through list of host numbers. these could be single numbers OR a range of numbers
    for hostnumbers in ${hostnumberlist[@]};
    do
        # is hostnumbers a range?
        if `echo ${hostnumbers} | grep "-" 1>/dev/null 2>&1`; then
          e=(`echo $hostnumbers | tr "-" "\n"`)
          s=${e[0]}
          l=${e[1]}

          # loop over the range and append every number to nodes
          for i in $(seq $s $l); do

            #make sure that single-digit numbers have a leading zero
            printf -v inputNo "%02d" $i
            nodes=( ${nodes[*]} $hostprefix${inputNo} );
          done
        else
          #hostnumbers is a single host, append it directly
          nodes=( ${nodes[*]} $hostprefix$hostnumbers )
        fi;
    done
    hostprefix=""
  done
  echo "#Nodes :" ${#nodes[@]}
  echo "start testing ..."
  echo

  # we need to switch to the check directorys of th global SCIP installation
  cd ..

  # check if the test set file exists
  TEST=clusterbench
  if [ ! -e "./check/testset/$TEST.test" ]; then
      echo Cannot find $TEST.test in `pwd`/check/testset/
      exit 1
  fi;

  echo "make testcluster"
  echo "     EXECUTABLE=$EXECUTABLE"
  echo "     QUEUE=$Q"
  echo "     EXCLUSIVE=true"
  echo "     TEST=$TEST"
  echo "     TIME=600"
  echo "     OUTPUTDIR=results/clusterbench"
  echo

  # execute the testing script on all nodes
  for n in ${nodes[*]}
  do
    # create an empty settings file named by the slurm queue, node, date, and time
    SETTINGS=$Q-$n-$DATETIME
    mkdir -p settings
    mkdir -p check/results/clusterbench
    touch settings/$SETTINGS.set

    # run full test set on each node
    echo "     CLUSTERNODES=$n"
    echo "     SETTINGS=$SETTINGS"
    make testcluster EXECUTABLE=$EXECUTABLE QUEUE=$Q CLUSTERNODES=$n EXCLUSIVE=true TEST=$TEST SETTINGS=$SETTINGS TIME=601 OUTPUTDIR=results/clusterbench

    # artificial settings file may only be removed after executable was started
  done
done

