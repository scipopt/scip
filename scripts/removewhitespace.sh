#!/usr/bin/env bash
#
# This script removes all trailing whitespace.
#
# There is nothing to adjust. Afterward please check the changes via "git
# diff" to make sure that the replacement worked

echo ""
echo "This script removes all whitespace in otherwise empty lines in src/"
echo ""

# collect all files
for DIR in src/*
do
    test -d $DIR || continue

    # external codes for which we do not want to modify whitespace
    case $DIR in src/amplmp | src/cppad | src/dejavu | src/nauty | src/tinycthread ) continue ;; esac

    FILES="$FILES $DIR/*.h $DIR/*.c"
done

for FILE in ${FILES[@]}
do
    if test -f $FILE
    then
        echo $FILE
        sed -i -e 's/\([ \t][ \t]*\)$//g' $FILE
    fi
done
