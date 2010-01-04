#!/bin/sh

OLDYEAR=2008
NEWYEAR=2010

DIRECTORIES=(check doc src src/scip src/objscip src/blockmemshell src/tclique examples examples/*/src examples/*/doc)
EXTENSIONS=(sh awk h c hpp cpp)

for DIRECTORY in ${DIRECTORIES[@]}
do
  for EXTENSION in ${EXTENSIONS[@]}
  do
    for FILE in $DIRECTORY/*.$EXTENSION 
    do
      if test -f $FILE
      then
	  # check if the file has a correct old date 

	  COUNT1=`grep -c 2002-$NEWYEAR $FILE`
	  COUNT2=`grep -c 1996-$NEWYEAR $FILE`

	  if test "$COUNT1" != "$COUNT2"
	  then
	      continue
	  fi
	  
	  COUNT1=`grep -c 2002-$OLDYEAR $FILE`
	  COUNT2=`grep -c 1996-$OLDYEAR $FILE`

	  if test "$COUNT1" == "$COUNT2"
	  then
	      echo "DATE ERROR --------------------> $FILE"
	      grep "2002-2" $FILE
	      grep "1996-2" $FILE	  
	  else
	      echo $FILE
  
	      mv $FILE $FILE.olddate
	      sed 's!2002-'$OLDYEAR'!2002-'$NEWYEAR'!g 
s!1996-'$OLDYEAR'!1996-'$NEWYEAR'!g' $FILE.olddate > $FILE
	  fi
      fi
    done
  done
done

