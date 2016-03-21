# regex for all canonical files in given directory
CHECKDIR="$1*.*"

# lists all files in directory which does not contain source in their name, but @ref in the file with their number of matches
GREPRESULT=`ls $CHECKDIR -a | grep "source" -v | xargs grep "@ref" -l | xargs grep "@ref" -c`

# gives them out
for RESULT in ${GREPRESULT[@]}
do
    if [ "${#RESULT}" -gt "2" ]
    then
        echo WARNING: unresolved references @ref in: $RESULT
    fi
done

# lists all files in directory which does not contain source in their name, but @ref in the file and classifies the unresolved references by their matching
echo classes of unresolved references:`ls $CHECKDIR -a | grep "source" -v | xargs grep -hP "@ref [A-Z][^ .<>]*" -o | sort | uniq | tee fails.txt`