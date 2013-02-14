#!/bin/bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# Bash-shell version of schulz script from GAMS performance tools (http://www.gamsworld.org/performance/schulz.htm)
#
# Supervises processes and sends given signals via kill when elapsed time exceeds given thresholds
#
# Usage: ./schulz.sh <process> <timelimits> <signals> <sleep>
#  <process>     is a regular expression which determines the commands that should be supervised (default: ^gms)
#  <timelimits>  is a colon separated list of timelimits (default: 120:180:240)
#  <signals>     is a colon separated list of signals (default: 2:1:9)
#  <sleep>       is the number of seconds to wait before checking the process list again (default: 60)
#
#  Whenever a process which command matches the first argument exceeds one of the given timelimits,
#  the corresponding signal is send to that process.
#
#  That is, for the default values, for each process which command that starts with ^gms the following happens:
#  If it is running for more than 240 seconds, signal 9 (SIGKILL) is send.
#  If it is running for more than 280 seconds, but less-or-equal 240 seconds, signal 1 (SIGHUP) is send.
#  If it is running for more than 120 seconds, but less-or-equal 180 seconds, signal 2 (SIGINT) is send.
#  If it is running for less-or-equal 120 seconds, nothing happens.
#
#  The list of timelimits and signals need to have the same length.

watch=${1:-^gms}
ress=${2:-120\:180\:240}
sigs=${3:-2\:1\:9}
sleepsec=${4:-60}

# split ress and sigs into arrays
IFS=":"
j=1
for r in $ress
do
  resl[$j]=$r
  ((j++))
done

j=1
for s in $sigs
do
  sigl[$j]=$s
  ((j++))
done

# check if number of timelimits equals number of signals
if test ${#resl[@]} -ne ${#sigl[@]}
then
  echo "List of times and list of signals need to have equal length."
  exit 1
fi

echo "--- " `date` ": Start watching $watch   sleep = $sleepsec seconds"
for (( j=1; j <= ${#resl[@]}; j++ ))
do
  echo "  Threshold ${resl[$j]}s -> Signal ${sigl[$j]}"
done

# given a date in format [[[dd-]hh:]mm:]ss,
# - calculates corresponding number of seconds and stores them in $secs
# - checks if one of the timelimits is exceeded and stores corresponding signal in $signal
#   if no timelimit is exceeded, $signal is set to 0
function getsignal() {
  secs=${1: -2:2}
  mins=${1: -5:2}
  hours=${1: -8:2}
  days=${1: -11:2}
  
  #remove 0's from beginning
  days=${days##0}
  hours=${hours##0}
  mins=${mins##0}
  secs=${secs##0}
  
  #echo $days "d" $hours "h" $mins "m" $secs "s"
  
  (( secs = secs + 60 * mins + 3600 * hours + 86400 * days ))
  
  signal=0
  for (( j=${#resl[@]}; j ; j-- ))
  do
    if test $secs -gt ${resl[$j]}
    then
      signal=${sigl[$j]}
      break
    fi
  done
}

# run forever
for ((;;))
do

sleep $sleepsec || exit

# get list of processes to watch in the format "pid,elapsedtime,command#pid,elapsedtime,command#..."
ap=`ps -Ao comm,pid,etime | grep -E "$watch" | awk '{ printf("%s,%s,%s#",$2,$3,tolower($1)) }'`

echo "--- " `date` ":"

# process list of processes to watch
IFS="#"
for i in $ap
do
  # get pid, elapsed time of process, and command
  pid=`echo $i | awk -F"," '{print $1}'`
  processtime=`echo $i | awk -F"," '{print $2}'`
  command=`echo $i | awk -F"," '{print $3}'`
  
  # check if we need to send a signal
  getsignal $processtime
  if test $signal -ne 0
  then
    echo "Signal $signal to $command with pid $pid running for $secs seconds"
    kill -$signal $pid > /dev/null
  else
    echo "$pid:$command:$secs"
  fi
done

done
