#!/usr/bin/env bash
#
# This script removes all whitespace in otherwise empty lines.
#
# There is nothing to adjust. Afterward please check the changes via "git
# diff" to make sure that the replacement worked

DIRECTORY="src/scip"

echo ""
echo "This script removes all whitespace in otherwise empty lines"
echo ""

DIRECTORIES=(src/scip src/blockmemshell src/dijkstra src/lpi src/objscip src/tclique src/xml)

# collect all files
for DIR in ${DIRECTORIES[@]}
do
    FILES="$FILES $DIR/*.h $DIR/*.c"
done

CNT=0;
for FILE in ${FILES[@]}
do
    if test -f $FILE
    then
	echo $FILE
	COUNT=`grep -c -h "^[ \t]\([ \t]*\)$" $FILE`
	CNT=`expr $CNT + $COUNT`
	mv $FILE $FILE.oldwhitespace
	sed 's/^\([ \t][ \t]*\)$//g' $FILE.oldwhitespace > $FILE
	rm $FILE.oldwhitespace
    fi
done
echo ""
echo "Replaced "$CNT" lines."
