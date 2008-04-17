#!/bin/sh
for i in Makefile INSTALL check/*.sh check/*.awk doc/*.c src/*.c src/*.cpp src/*.h src/blockmemshell/*.c src/blockmemshell/*.h src/objscip/*.cpp src/objscip/*.h src/scip/*.c src/scip/*.cpp src/scip/*.h src/tclique/*.c src/tclique/*.h
do
if [ -f $i ]
then
echo $i
mv $i $i.olddate
sed 's!2002-2007!2002-2008!g
s!1996-2007!1996-2008!g' $i.olddate > $i
fi
done
