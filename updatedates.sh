#!/usr/bin/env bash
#
# This bash script updates all copyrights in the SCIP files and posted
# those files which contain a copyright which do not have the right format
#
# You just have to run this script. There is nothing to adjust. 
# The correct year is detected through the 'date' function 
#
# Note that not all files (usually scripts) contain a copyright. A copyright is only 
# needed for those files which are part of a SCIP distribution (see makedist.sh)
#
# This bash script also changes the copyrights of the SCIP examples
#

NEWYEAR=`date +"%Y"`
LASTYEAR=`expr $NEWYEAR - 1`

DIRECTORIES=(check doc src src/* examples examples/* examples/*/src examples/*/doc)
EXTENSIONS=(sh awk h c hpp cpp)
EXTRAFILE=(Makefile)

echo ""
echo "This script reports *all* files which have not a correct COPYRIGHT."
echo "Only files which are included in the distribution need a COPYRIGHT (see makedist.sh)"
echo ""

for DIRECTORY in ${DIRECTORIES[@]}
do
  for EXTENSION in ${EXTENSIONS[@]}
  do
    for FILE in $DIRECTORY/*.$EXTENSION $DIRECTORY/$EXTRAFILE $EXTRAFILE
    do
      if test -f $FILE
      then
	  # check if the file has already the new copyright 
	  COUNT1=`grep -c 2002-$NEWYEAR $FILE`
	  COUNT2=`grep -c 1996-$NEWYEAR $FILE`

	  if test "$COUNT1" != "$COUNT2"
	  then
	      continue
	  fi
	  
	  # check if the file has a correct old copyright 
	  COUNT1=`grep -c 2002-$LASTYEAR $FILE`
	  COUNT2=`grep -c 1996-$LASTYEAR $FILE`

	  if test "$COUNT1" == "$COUNT2"
	  then
	      # post those files which have a wrong old copyright
	      echo "COPYRIGHT ERROR --------------------> $FILE"
	      grep "2002-2" $FILE
	      grep "1996-2" $FILE	  
	  else
	      echo $FILE
  
	      mv $FILE $FILE.olddate
	      sed 's!2002-'$LASTYEAR'!2002-'$NEWYEAR'!g 
                   s!1996-'$LASTYEAR'!1996-'$NEWYEAR'!g' $FILE.olddate > $FILE
	  fi
      fi
    done
  done
done

