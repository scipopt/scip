#!/usr/bin/env bash
LIMIT=$1
QUEUE=$2

while true
do
  ALLQUEUE=`qstat | grep -c ""`
  QUEUE=`qstat -q $QUEUE | grep -c $USER`
  RUNNING=`qstat -u $USER | grep -c " R "`

  # display current user load and total load 
  echo jobs in progress: $RUNNING / $QUEUE "("$ALLQUEUE")"
 
  if test $ALLQUEUE -le 1990
      then
      if test $QUEUE -le 400
	  then
	  break
      fi

      if test $ALLQUEUE -le $LIMIT
	  then
	  break
      fi
  fi

  sleep 30
done
