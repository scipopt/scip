#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

NEEDED_PARAMS=(HARDMEMLIMITMB MEMLIMITMB MIPLIBPATH JOBSET SOLVER BINNAME TIMELIMIT THREADS MIPGAP LINTOL INTTOL OUTFILEPATTERN SOLFILEPATTERN )

# check that each of the needed parameters is set
MISSING="false"
for PARAM in ${NEEDED_PARAMS[@]}
do
    if [[ -z ${!PARAM} ]]
    then
        echo "There is no value for variable '${PARAM}'."
        MISSING="true"
    fi
done
if [[ ${MISSING} == "true" ]]
then
    echo "At least one of the needed enviroment variables is missing."
    exit -1
fi

# import some useful functions that make this script less cluttered
. ${MIPLIBPATH}/scripts/run_functions.sh

function printDateTime() {
    echo "-----------------------------"
    date
    echo "-----------------------------"
}

# Copies and then deletes the working directory and all of its content on exit
function cleanUp() {
    for RELATIVEOUTFILE in ${OUT_FILES}
    do
        echo "Appending temporary out data from ${WORKINGDIR}/${RELATIVEOUTFILE} to ${MIPLIBPATH}/${RELATIVEOUTFILE}"
        cat ${WORKINGDIR}/${RELATIVEOUTFILE} >> ${MIPLIBPATH}/${RELATIVEOUTFILE}
    done
    for REL_SOL_FILE in ${SOL_FILES}
    do
        echo "Appending temporary solution data from ${WORKINGDIR}/${RELATIVESOLFILE} to ${MIPLIBPATH}/${RELATIVESOLFILE}"
        cat ${WORKINGDIR}/${RELATIVESOLFILE} >> ${MIPLIBPATH}/${RELATIVESOLFILE}
    done
    if [[ -e ${WORKINGDIR} ]]
    then
        rm -rf ${WORKINGDIR}
    fi
}

# copy the output and solution and remove the temporary files on exit
trap cleanUp EXIT

# the temp directoryand files are needed in case this script is executed on a remote cluster
WORKINGDIR=$(mktemp -d)
echo "Created temporary directory ${WORKINGDIR}"
OUT_FILES=""
SOL_FILES=""

# convert hard memory limit to kilo bytes
HARDMEMLIMITKB=$(expr ${HARDMEMLIMITMB} \* 1024)

for INSTANCE in ${JOBSET}
do
    # get the paths to the relative out and sol file
    INSTANCENAME=$(stripPathAndInstanceExtension ${INSTANCE})
    RELATIVEOUTFILE=$(replaceInPattern ${OUTFILEPATTERN} ${INSTANCENAME})
    RELATIVESOLFILE=$(replaceInPattern ${SOLFILEPATTERN} ${INSTANCENAME})
    
    TMPOUTFILE="${WORKINGDIR}/${RELATIVEOUTFILE}"
    TMPSOLFILE="${WORKINGDIR}/${RELATIVESOLFILE}"
    
    echo "Writing out data temporarily to ${TMPOUTFILE}"
    echo "Writing solution data temporarily to ${TMPSOLFILE}"
    
    OUT_FILES="${OUT_FILES} ${RELATIVEOUTFILE}"
    SOL_FILES="${SOL_FILES} ${RELATIVESOLFILE}"
    
    # create empty files in the given path while creating all necessary parent directories via install
    install -Dm 640 /dev/null ${TMPOUTFILE}
    install -Dm 640 /dev/null ${TMPSOLFILE}

    # post system information and current time into the output file
    uname -a >> ${TMPOUTFILE}
    date >> ${TMPOUTFILE}
    cat /proc/cpuinfo | grep 'model name' | uniq | sed 's/model name/CPU/g' >> ${TMPOUTFILE}

    # post MIPLIB script version
    echo "MIPLIB script version ${VERSION}" >> ${TMPOUTFILE}

    # print the hard memory limit (in kb) into the output file
    echo "hard mem limit: ${HARDMEMLIMITKB} k" >> ${TMPOUTFILE}

    if [[ -f ${MIPLIBPATH}/${INSTANCE} ]]
    then
        echo "@01 ${INSTANCE} ==========="
        printDateTime
        TIMESTART=$(date +"%s")
        echo "@03 ${TIMESTART}"
        bash -c " ulimit -v ${HARDMEMLIMITKB} k; cd ${WORKINGDIR}; ${MIPLIBPATH}/scripts/run_${SOLVER}.sh ${SOLVER} ${MIPLIBPATH}/bin/${BINNAME} ${MIPLIBPATH}/${INSTANCE} ${TIMELIMIT} ${MEMLIMITMB} ${TMPSOLFILE} ${THREADS} ${MIPGAP}"
        echo
        TIMEEND=`date +"%s"`
        echo "@04 ${TIMEEND}"
        echo "@05 ${TIMELIMIT}"

        # if we have a (non-empty) solution and a checker: check the solution
        if [[ -e ${TMPSOLFILE} ]]
        then
            if [[ -e ${MIPLIBPATH}/checker/bin/solchecker ]]
            then
                echo ""
                bash -c " ${MIPLIBPATH}/checker/bin/solchecker ${MIPLIBPATH}/${INSTANCE} ${TMPSOLFILE} ${LINTOL} ${INTTOL}"
                echo ""
            else
                echo "solution checker not found!"
            fi
        else
            echo "Solution file is empty."
        fi
        printDateTime
        echo ""
        echo "=ready="
    else
        echo "@02 FILE NOT FOUND: ${INSTANCE} ==========="
    fi 2>&1 | tee -a ${TMPOUTFILE}
done
