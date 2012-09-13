#!/usr/bin/env bash
#
# This bash script updates all "extern" of public methods with "EXTERN".
# The latter is a macro define in def.h and is needed to easily create
# libraries for Windows. In case of any other system this macro "EXTERN" is
# replaces by "extern".
#
# This scripts should be run before every release. You just have to run
# this script as it is from the main directory of SCIP.
#
# > scripts/updateextern.sh
#
# There is nothing to adjust. Afterward please check the changes via "git
# diff" to make sure that the replacement worked

DIRECTORY="src/scip"
EXTRAFILES=(scip)
PLUGINTYPES=(branch cons dialog disp event heur message nodesel prop presol pricer reader relax sepa)

echo ""
echo "This script changes all "extern" to "EXTERN" in the header of public callable methods"
echo ""

FILES="src/scip/scip*.h src/scip/pub_*.h src/nlpi/pub_expr.h src/nlpi/nlpi.h src/nlpi/exprinterpret.h  src/blockmemshell/memory.h"

# collect all header files related to plugins
for PLUGINTYPE in ${PLUGINTYPES[@]}
do
    for FILE in $DIRECTORY/$PLUGINTYPE"_"*.h
    do
	# ignore xyz plugin templates
	if test "$FILE" = "$DIRECTORY/$PLUGINTYPE"_"xyz.h"
	then
	    continue
	fi

	if test -f $FILE
	then
	    FILES="$FILES $FILE"
	fi
    done
done

# collect all public header files
for FILE in $DIRECTORY/pub_*.h
do
    if test -t $FILE
    then
	FILES="$FILES $FILE"
    fi
done

# collect the header of the extra files
for EXTRAFILE in ${EXTRAFILES[@]}
do
    if test -f $DIRECTORY/$EXTRAFILE
    then
	FILES="$FILES $DIRECTORY/$EXTRAFILE"
    fi
done

# replace in all collected header files the lines which start and end with
# "extern" by "EXTERN"; note that we also consider trailing white spaces
for FILE in ${FILES[@]}
do
    if test -f $FILE
    then
	COUNT=`grep -c "^extern\([ ]*\)$" $FILE`

	if test $COUNT -gt 0
	then
	    echo "--> replaced $COUNT extern in $FILE"

	    mv $FILE $FILE.oldextern
	    sed 's/^extern\([ ]*\)$/EXTERN/g' $FILE.oldextern > $FILE
	fi
    fi
done
