#!/bin/bash
LIMIT=$1

while true
do
  ALLQUEUE=`qstat | grep -c ""`
  QUEUE=`qstat | grep -c $USER`
  RUNNING=`qstat -u $USER | grep -c " R "`

  # display current user load and total load 
  echo jobs in progress: $RUNNING / $QUEUE "("$ALLQUEUE")"
  
  if test $ALLQUEUE -le $LIMIT
      then
      break
  fi

  sleep 30
done
