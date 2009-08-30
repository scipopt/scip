#!/bin/bash
while true
do
  QUEUE=`qstat | grep -c $USER`
  RUNNING=`qstat -u $USER | grep -c " R "`
  echo jobs in progress: $RUNNING / $QUEUE

  if test $QUEUE -le 20
  then
     break
  fi

  sleep 30
done
