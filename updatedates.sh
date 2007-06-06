#!/bin/sh
#for i in Makefile INSTALL check/*.sh check/*.awk doc/*.c src/*.c src/*.cpp src/*.h src/blockmemshell/*.c src/blockmemshell/*.h src/objscip/*.cpp src/objscip/*.h src/scip/*.c src/scip/*.h src/tclique/*.c src/tclique/*.h
for i in check/*.sh check/*.awk doc/*.c
do
echo $i
mv $i $i.olddate
sed 's!2002-2006!2002-2007!g
s!1996-2006!1996-2007!g' $i.olddate > $i
done
