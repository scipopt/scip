#! /bin/bash

i=0

mkdir -p newfiles/scip
rm newfiles/scip/*

function getFileLen {
    filename=$1
    len=`wc -l $filename | awk '{print $1}'`
    echo $len
}

while read line
do
    i=$((i + 1))
    echo $i $line
    doxygroup=`echo $line | awk '{print $1}'`
    headername=`echo $line | awk '{print $2}'`
    echo $doxygroup
    export headername
    echo   "/${doxygroup}/,/@}/p"
    export header=newfiles/scip/scip_${headername}.h

    if [ ! -f $header ]
    then
        export doxygengroup=PUBLICCOREAPI
        export doxygenbrief="public methods for $headername"
        export guard="__SCIP_SCIP_`echo $headername | tr [:lower:] [:upper:]`_H__"

        envsubst < header_head_template.h > $header
    else
        #
        # delete the footer from the previous round with an inline sed
        #
        len=`getFileLen $header`
        footlen=`getFileLen header_foot.h`
        startline=`expr $len - $footlen + 2`

        sed -i "${startline},$$d" $header
    fi

    sed -n "/${doxygroup}/,/@}/p" ../../src/scip/scip.h >> $header
    cat header_foot.h >> $header
done < groups_names.list

#
# use the write gaps functionality to write the gap file
#
./parse_functions.py ../../src/scip/scip.{c,h} --write_gaps

for i in newfiles/scip/*.h
do
    header=$i
    module=newfiles/scip/$(basename $i .h).c

    echo $module

    ./parse_functions.py ../../src/scip/scip.c $header > tmp.list

    modulebasename=$(basename $module)
    export modulebasename
    envsubst < preamble.c > $module

    cat all_includes.c >> $module

    if [ $header = newfiles/scip/scip_numerics.h ]
    then
        echo >> $module
        sed -n '/In debug mode, the following methods are implemented as function calls to ensure/,/undef SCIPgetHugeValue/p' ../../src/scip/scip.c >> $module
        echo >> $module
    fi

    if [ $header = newfiles/scip/scip_expr.h ]
    then
        cat >>$module <<EOL
/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) (val >= infty1 ? infty2 : val)

EOL
    fi

    if [ $header = newfiles/scip/scip_nlp.h ]
    then
        sed -n '/method to call, when the priority of an NLPI was changed/,/^\}/p' ../../src/scip/scip.c >> $module
    fi

    sed -n -f tmp.list ../../src/scip/scip.c >> $module

    if [ $header = newfiles/scip/scip_expr.h ]
    then
        echo >> $module
        echo '#undef infty2infty' >> $module
    fi
done

rm newfiles/scip/scip_bandit.c

./organizeincludes.sh

#
# replace doxygen brief descriptions
#
./descriptions.sh
