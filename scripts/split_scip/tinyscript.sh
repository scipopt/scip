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
    export headername
    echo   "/${doxygroup}/,/@}/p"
    export header=newfiles/scip_${headername}.h
    export doxygengroup=PUBLICCOREAPI
    export doxygenbrief="public methods for $headername"
    export guard="__SCIP_SCIP_`echo $headername | tr [:lower:] [:upper:]`_H__"

    envsubst < header_head_template.h > $header
    sed -n "/${doxygroup}/,/@}/p" ../../src/scip/scip.h >> $header
    cat header_foot.h >> $header
done < groups_names.list

for i in newfiles/*.h
do
    header=$i
    module=newfiles/$(basename $i .h).c

    ./parse_functions.py ../../src/scip/scip.c $header > tmp.list
    sed -n -f tmp.list ../../src/scip/scip.c > $module
done

