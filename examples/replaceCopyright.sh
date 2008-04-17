#!/bin/sh
mkdir tmp
for i in */src 
  do
  for j in {$i/*.h,$i/*.c,$i/*.cpp,$i/*.hpp}; do
      file=`basename $j`
      if [ -f $i/$file ] 
	  then
	  echo $i/$file
	  sed 's/              2002-2007 Konrad-Zuse-Zentrum/Copyright (C) 2002-2008 Konrad-Zuse-Zentrum/g' $i/$file > tmp/$file
	  mv tmp/$file $i/$file
      fi
  done
done
for i in *
  do
  for j in {INSTALL,Makefile}; do
      if [ -f $i/$j ] 
	  then
	  echo $i/$j
	  sed 's/              2002-2007 Konrad-Zuse-Zentrum/Copyright (C) 2002-2008 Konrad-Zuse-Zentrum/g' $i/$j > tmp/$j
	  mv tmp/$j $i/$j
      fi
  done
done
rmdir tmp
