#!/bin/bash

# moves all $USER.* files and directories to directory $TARGETPATH on all machines

if test -z $USER
then
    echo $USER
    exit
fi

if [ $# -gt 5  ] || [ $# -eq 0 ]; then
  echo 'usage: ./collectnodes.sh ACCOUNT QUEUE TARGETPATH USESSH'
  echo '       note: queues can be separated by comma'
  exit
fi;

ACCOUNT=$1
QUEUE=$2
TARGETPATH=$3
USESSH=$4

echo "account:" $ACCOUNT
echo "queue  :" $QUEUE
echo "target :" $TARGETPATH

if test ! -e $TARGETPATH;
then
    echo "target does not exist"
    exit
fi

declare -a opts


# multiple queues are split by a comma (hopefully, no one ever attempts to define slurm partitions with a comma ;))
split_queue=(`echo $QUEUE | tr "," "\n"`)

for Q in ${split_queue[@]}; do

  # evaluate sinfo for this partition and collect the nodes
  SINFO=$(`echo sinfo -p ${Q}` )
  LSINFO=(`echo ${SINFO[@]} | wc -w`)

  echo SINFO is ${SINFO[@]}
  echo LSINFO is $LSINFO
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
        echo moreopts: ${moreopts[@]}
        opts=( ${opts[@]} ${moreopts[@]} );
      fi;
    fi;
  done
done

echo opts is ${opts[@]}
declare -a nodes

# Expand compressed slurm array
for nodelist in ${opts[@]}; do
  subnodelist=(`echo $nodelist | tr "[]" "\n"`)
  echo subnodelist is ${subnodelist[@]}

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
          nodes=( ${nodes[*]} $hostprefix${inputNo} ); done
      else

      #hostnumbers is a single host, append it directly
      nodes=( ${nodes[*]} $hostprefix$hostnumbers )
      fi;
  done
  hostprefix=""

done

echo "nodes  :" ${nodes[*]}
echo Number of nodes: ${#nodes[@]}

export TARGETPATH=$TARGETPATH

for i in ${nodes[@]}
do
  if [ "$USESSH" != "true" ];
  then
    /usr/local/bin/ZIBopti-wol $i 1>/dev/null
    sbatch --job-name=$i-grab -p $QUEUE -w $i -A $ACCOUNT --output=/dev/null collectnode.sh
  else
    ssh -oStrictHostKeyChecking=no -i $HOME/.ssh/id_rsa $i "mv -nv /usr/local/tmp/${USER}-tmpdir/${USER}.* ${TARGETPATH}/"
  fi
done
