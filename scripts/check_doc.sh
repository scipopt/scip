#! /bin/bash

# filters doxygen documentation files for occurrences of unresolved references
#
# usage: check_doc.sh <DIRECTORYNAME>
#
# DIRECTORYNAME should be the root directory of doxygen-generated documentation, e.g., 'doc/html'

# get file list after excluding source files
FILELIST=`ls $1/*.* -a | grep "_source" -v`

# lists all files in directory which do not contain source in their name,
# but '@ref' in the file with their number of matches
GREPRESULT=`grep "@ref" -l $FILELIST | xargs grep "@ref" -cH`

# gives them out
for RESULT in ${GREPRESULT[@]}
do
    # check length of RESULT variable: If grep command finds no match, GREPRESULT is '0'
    if [ "${#RESULT}" -gt "2" ]
    then
        echo === Warning: unresolved references @ref in: $RESULT
    fi
done

# lists all files in directory which do not contain source in their name,
# but '@ref' in the file and classifies the unresolved references by their matching
echo === Classes of unresolved references:
echo `grep -hP "@ref [A-Za-z][^ .<>]*" -o $FILELIST | sort | uniq`
