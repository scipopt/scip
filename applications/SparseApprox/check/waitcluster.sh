#!/usr/bin/env bash
LIMIT=$1
QUEUE=$2
LIMITQUEUE=$3

while true
do
  ALLQUEUED=`qstat | grep -c ""`
  QUEUED=`qstat -u $USER | grep -c " $QUEUE "`
  RUNNING=`qstat -u $USER | grep -c " R "`

  # display current user load and total load 
  echo jobs in progress: $RUNNING / $QUEUED "("$ALLQUEUED")"
 
  if test $ALLQUEUED -le 1990
      then

      if test $QUEUED -le $LIMITQUEUE
	  then
	  break
      fi

      if test $ALLQUEUED -le $LIMIT
	  then
	  break
      fi
  fi

  sleep 30
done
