#!/usr/bin/env bash

#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# Script parameter
VERSION=${1}      # (string) MIPLIB script version
BINNAME=${2}      # (string) name of the binary in the 'bin' dir
TSTNAME=${3}      # (string) name of the testset to execute
TIMELIMIT=${4}    # (int > 0) in seconds
MEMLIMITMB=${5}   # (int > 0) in MB
THREADS=${6}      # (int > 0) number of threads to be used by the solver
PERMUTE=${7}      # (int >= 0) permutation to run

# Additional parameter
LINTOL=1e-4 # absolut tolerance for checking linear constraints and objective value
INTTOL=1e-4 # absolut tolerance for checking integrality constraints
MIPGAP=0.0  # Note that the MIP gap (gap between primal and dual solution) is not
            # uniqely defined through all solvers. For example, there is a difference
            # between SCIP and CPLEX. All solver, however, have the some behaviour in
            # case of a MIP gap of 0.0.

# check if all variables defined (by checking the last one)
if [[ -z "${PERMUTE}" ]]
then
    echo "Skipping test since not all variables are defined"
    echo "VERSION      = ${VERSION}"
    echo "BINNAME      = ${BINNAME}"
    echo "TSTNAME      = ${TSTNAME}"
    echo "TIMELIMIT    = ${TIMELIMIT}"
    echo "MEMLIMITMB   = ${MEMLIMITMB}"
    echo "THREADS      = ${THREADS}"
    echo "PERMUTE      = ${PERMUTE}"
    exit -1
fi

# import some useful functions that make this script less cluttered
. $(dirname "${BASH_SOURCE[0]}")/run_functions.sh

# construct paths
MIPLIBPATH=$(pwd)
BINPATH=${MIPLIBPATH}/bin
CHECKERPATH=${MIPLIBPATH}/checker
RESULTSDIR=results
RESULTSPATH=${MIPLIBPATH}/${RESULTSDIR}
SOLUTIONDIR=${RESULTSDIR}/solutions
SOLUTIONPATH=${MIPLIBPATH}/${SOLUTIONDIR}
SCRIPTPATH=${MIPLIBPATH}/scripts
TSTPATH=${MIPLIBPATH}/testsets

# check if the solver link (binary) exists
if [[ ! -e "${BINPATH}/${BINNAME}" ]]
then
    echo "ERROR: solver link '${BINNAME}' does not exist in 'bin' folder; see bin/README"
    exit -1
fi

# check if the test set file/link exists
if [[ ! -e "${TSTPATH}/${TSTNAME}.test" ]]
then
    echo "ERROR: test set file/link '${TSTNAME}.test' does not exist in 'testset' folder"
    exit -1
fi

SOLUFILE=${TSTPATH}/${TSTNAME}.solu

# check if a solution  file/link exists
if [[ ! -e ${SOLUFILE} ]]
then
    echo "Warning: solution file/link '${TSTNAME}.solu' does not exist in 'testset' folder; therefore, no consistency check"
    SOLUFILE=""
fi

# check if the result folder exist. if not create the result folder
if [[ ! -e ${RESULTSPATH} ]]
then
    mkdir ${RESULTSPATH}
fi

# check if the solution folder exist. if not create the solution folder
if [[ ! -e $SOLUTIONPATH ]]
then
    mkdir ${SOLUTIONPATH}
fi

if [[ ${PERMUTE} -gt 0 ]]
then
    # construct name of output and results file
    BASENAME="${RESULTSDIR}/${TSTNAME}.${BINNAME}.${THREADS}threads.${TIMELIMIT}s.p${PERMUTE}"
    SOLFILEPATTERN="${SOLUTIONDIR}/${TSTNAME}.{INSTANCENAME}.${BINNAME}.${THREADS}threads.${TIMELIMIT}s.p${PERMUTE}.sol"
else
    # construct name of output and results file
    BASENAME="${RESULTSDIR}/${TSTNAME}.${BINNAME}.${THREADS}threads.${TIMELIMIT}s"
    SOLFILEPATTERN="${SOLUTIONDIR}/${TSTNAME}.{INSTANCENAME}.${BINNAME}.${THREADS}threads.${TIMELIMIT}s.sol"
fi

METAFILE=${BASENAME}.meta
OUTFILE=${BASENAME}.out
RESFILE=${BASENAME}.res

# we add a small memory overhead for the cluster and file system, such that a small violation
# does not result in a killed process
HARDMEMLIMITMB=$(( ${MEMLIMITMB} + ${MEMLIMITMB}/10 + 100 ))

writeMetaFile ${METAFILE} ${PERMUTE} ${TSTNAME} ${BINNAME} ${TIMELIMIT} ${MEMLIMITMB} ${THREADS}
echo > ${MIPLIBPATH}/${OUTFILE}
echo > ${MIPLIBPATH}/${RESFILE}

# grep solver name
SOLVER=$(stripversion ${BINNAME})

# loop over all instance names which are listed in the test set file name
for i in $(cat ${TSTPATH}/${TSTNAME}.test)
do
    # check if the current instance exists
    if [[ -f ${i} ]]
    then
        #permute problem
        if [[ ${PERMUTE} -gt 0 ]]
        then
            INSTANCE=$(permuteInstance ${i} ${PERMUTE} ${BINPATH})
        else
            INSTANCE=${i}
        fi

        FILENAME=$(basename ${INSTANCE} ".gz")

        # export the variables used by runjob.sh
        export HARDMEMLIMITMB
        export MEMLIMITMB
        export MIPLIBPATH
        export JOBSET=${INSTANCE}
        export SOLVER
        export BINNAME
        export TIMELIMIT
        export THREADS
        export MIPGAP
        export LINTOL
        export INTTOL
        export OUTFILEPATTERN=${OUTFILE}
        export SOLFILEPATTERN
        export VERSION

        ${SCRIPTPATH}/runjob.sh
    else
        echo "Instance ${i} not found, skipping it."
    fi
done

awk -f ${SCRIPTPATH}/parse.awk -f  ${SCRIPTPATH}/parse_${SOLVER}.awk -v "LINTOL=${LINTOL}" ${SOLUFILE} ${OUTFILE} | tee ${RESFILE}
