#! /bin/bash

#
# train search completion approximations from SCIP log files with tree size estimation reports
#
# This script processes an input directory with SCIP log files to train
# two search completion approximations; a regression forest and a linear regression, the so-called monotone regression.
# The training is performed within an R script, and the output can be input into SCIP
# to change the appearance of the search completion display column.
#
# The only input specification is the name of a directory with SCIP log files that contain periodic tree size estimation reports.
#
# Please see the README.md for installation instructions, further information on creating the log file output, and
# help in interpreting the results
#
# Note that the training only happens on instances that have been reported solved
#
# Ensure that different settings are not mixed within one directory
#
echo "Usage: ${0} <logfile-directory> [outputdir]"

# read log file directory
if [ -z ${1} ]
then
   echo "Exiting because missing logfile directory (mandatory)"
   exit 1
else
   logfiledirectory=${1}
fi

# read optional output directory
if [ -z ${2} ]
then
   outputdir="output"
else
   outputdir=${2}
fi
echo
echo "Storing all relevant information in directory ${outputdir}/"
echo
mkdir -p ${outputdir}

# search for solved problems
 grep -l "problem is solved" ${logfiledirectory}/*.out > ${outputdir}/solvedproblems.list

 echo
 echo "Number of solved problems: $(cat ${outputdir}/solvedproblems.list | wc -l), stored in ${outputdir}/solvedproblems.list"
 echo

# parse reports of all log files of all solved problems
 for i in `cat ${outputdir}/solvedproblems.list`
 do
    awk -f $(dirname ${0})/parse_logfile.awk $i
 done | sed -e "s/progress/tree-weight/g" > ${outputdir}/table1.csv

# parse actual total time and number of nodes from the log files
 for i in `cat ${outputdir}/solvedproblems.list`
 do
    actual=$(grep "Tree Data\|Estimation Tree" $i |\
                tail -n 1 | \
                grep -oP "(\d+) nodes" | \
                sed 's/ nodes//g')
    time=$(grep -i "Total Time" $i | \
             grep -oP '[^ ]+$' \
            )
    echo "$i,$actual,$time"
done > ${outputdir}/actual.csv

echo
echo "Running ./train.R"
echo
# call R training script
./train.R ${outputdir}

echo
echo "Training completed. Find more information and learned methods in directory ${outputdir}/"
ls -l ${outputdir}