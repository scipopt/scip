#! /bin/bash

i=0

mkdir -p newfiles
rm newfiles/*

while read line
do
    i=$((i + 1))
    echo $i $line
    doxygroup=`echo $line | awk '{print $1}'`
    headername=`echo $line | awk '{print $2}'`
    echo $doxygroup
    echo $headername
    echo   "/${doxygroup}/,/@}/p"
    header=newfiles/scip_${headername}.h
    echo >> $header
    sed -n "/${doxygroup}/,/@}/p" ../../src/scip/scip.h >> $header
    echo >> $header
done < groups_names.list

for i in newfiles/*.h
do
    header=$i
    module=newfiles/$(basename $i .h).c

    ./parse_functions.py ../../src/scip/scip.c $header > tmp.list
    sed -n -f tmp.list ../../src/scip/scip.c > $module
done

