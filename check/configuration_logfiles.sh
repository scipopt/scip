#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      *
#*                                                                           *
#*  Licensed under the Apache License, Version 2.0 (the "License");          *
#*  you may not use this file except in compliance with the License.         *
#*  You may obtain a copy of the License at                                  *
#*                                                                           *
#*      http://www.apache.org/licenses/LICENSE-2.0                           *
#*                                                                           *
#*  Unless required by applicable law or agreed to in writing, software      *
#*  distributed under the License is distributed on an "AS IS" BASIS,        *
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
#*  See the License for the specific language governing permissions and      *
#*  limitations under the License.                                           *
#*                                                                           *
#*  You should have received a copy of the Apache-2.0 license                *
#*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

### configures the right test output files such as the .eval, the .tmp and the .set
### files to run a test on.
### the invoking script should pass "init" as argument to signalize that
### files need to be reset

### environment variables declared in this script
### OUTFILE - the name of the (sticked together) output file
### ERRFILE - the name of the (sticked together) error file
### EVALFILE - evaluation file to glue single output and error files together
### CHECKSETFILE  - the name of the basic settings file starting with 'check.'
### OBJECTIVEVAL - the optimal or best-know objective value for this instance
### SHORTPROBNAME - the basename of "${INSTANCE}" without file extension
### FILENAME - the basename of the local files (.out, .tmp, and .err)
### SKIPINSTANCE - should the instance be skipped because it was already evaluated in a previous setting?
### BASENAME - "${SCIPPATH}/${OUTPUTDIR}/${FILENAME}" cf. FILENAME argument
### TMPFILE  - the batch file name to pass for solver instructions
### SETFILE  - the name of the settings file to save solver settings to

### environment variables passed as arguments to this script
INIT="${1}"          # should log files be initialized (this overwrite or copy/move some existing log files)
COUNT="${2}"         # the instance count as part of the filename
INSTANCE="${3}"      # the name of the instance
BINID="${4}"         # the ID of the binary to use
PERMUTE="${5}"       # the number of permutations to use - 0 for no permutation
SEEDS="${6}"         # the number of random seeds - 0 only default seeds
SETNAME="${7}"       # the name of the setting
TSTNAME="${8}"       # the name of the testset
CONTINUE="${9}"      # should test continue an existing run
QUEUE="${10}"        # the queue name
p="${11}"            # shift of the global permutation seed
s="${12}"            # shift of the global random seed
THREADS="${13}"      # the number of threads
GLBSEEDSHIFT="${14}" # the global seed shift
STARTPERM="${15}"    # the starting permutation

# common naming scheme for eval files
CHECKBASENAME="${SCIPPATH}/${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${QUEUE}.${SETNAME}"

# if number of threads is larger than 1, add postfix
if test "${THREADS}" -gt 1
then
    CHECKBASENAME="${CHECKBASENAME}-t${THREADS}"
fi

# if seed is positive, add postfix
SEED=$((s + GLBSEEDSHIFT))
if test "${SEED}" -gt 0
then
    CHECKBASENAME="${CHECKBASENAME}-s${SEED}"
fi

# if permutation is positive, add postfix
PERM=$((p + STARTPERM))
if test "${PERM}" -gt 0
then
    CHECKBASENAME="${CHECKBASENAME}-p${PERM}"
fi

OUTFILE="${CHECKBASENAME}.out"
ERRFILE="${CHECKBASENAME}.err"
CHECKSETFILE="${CHECKBASENAME}.set"
EVALFILE="${CHECKBASENAME}.eval"
METAFILE="${CHECKBASENAME}.meta"

# create meta file
if ! test -e "${METAFILE}"
then
    echo "@Permutation ${PERM}"             >  "${METAFILE}"
    echo "@Seed ${SEED}"                    >> "${METAFILE}"
    echo "@Settings ${SETNAME}"             >> "${METAFILE}"
    echo "@TstName ${TSTNAME}"              >> "${METAFILE}"
    echo "@BinName ${BINNAME}"              >> "${METAFILE}"
    echo "@NodeLimit ${NODELIMIT}"          >> "${METAFILE}"
    echo "@MemLimit ${MEMLIMIT}"            >> "${METAFILE}"
    echo "@Threads ${THREADS}"              >> "${METAFILE}"
    echo "@FeasTol ${FEASTOL}"              >> "${METAFILE}"
    echo "@Queue ${QUEUE}"                  >> "${METAFILE}"
    echo "@Exclusive ${EXCLUSIVE}"          >> "${METAFILE}"
    if [ "${CLUSTERBENCHMARK}" == "yes" ]; then
        echo "@QueueNode ${CB_QUEUENODE}"   >> "${METAFILE}"
        echo "@ClusterBenchmarkID ${CB_ID}" >> "${METAFILE}"
    fi
fi

if test "${INSTANCE}" = "DONE"
then
    return
fi

# reset files if flag is set to 'init'
if test "${INIT}" = "true"
then
    #reset the eval file
    echo > "${EVALFILE}"

    #mv existing out and error files
    if test "${CONTINUE}" = "true"
    then
        MVORCP=cp
    else
        MVORCP=mv
    fi
    DATEINT=$(date +"%s")
    for FILE in OUTFILE ERRFILE CHECKSETFILE
    do
        if test -e "${FILE}"
        then
            "${MVORCP}" "${FILE}" "${FILE}.old-${DATEINT}"
        fi
    done
fi


# filter all parseable file format extensions
SHORTPROBNAME=$(basename "${INSTANCE}" .gz)
for EXTENSION in .mps .lp .opb .gms .pip .zpl .cip .fzn .osil .wbo .cnf .difflist .cbf .dat-s
do
    SHORTPROBNAME=$(basename "${SHORTPROBNAME}" "${EXTENSION}")
done

# get objective value from solution file
# we do this here to have it available for all solvers, even though it is not really related to logfiles
if test -e "${SOLUFILE}"
then
    # get the objective value from the solution file: grep for the instance name and only use entries with an optimal or best known value;
    # if there are multiple entries for this instance in the solution file, sort them by objective value and take the objective value
    # written in the last line, i.e., the largest value;
    # as a double-check, we do the same again, but reverse the sorting to get the smallest value
    OBJECTIVEVAL=$(grep " ${SHORTPROBNAME} " "${SOLUFILE}" | grep -e =opt= -e =best= | sort -k 3 -g | tail -n 1 | awk '{print $3}')
    CHECKOBJECTIVEVAL=$(grep " ${SHORTPROBNAME} " "${SOLUFILE}" | grep -e =opt= -e =best= | sort -k 3 -g -r | tail -n 1 | awk '{print $3}')

    # if largest and smalles reference value given in the solution file differ by more than 1e-04, stop because of this inconsistency
    if awk -v n1="${OBJECTIVEVAL}" -v n2="${CHECKOBJECTIVEVAL}" 'BEGIN { exit (n1 <= n2 + 0.0001 && n2 <= n1 + 0.0001) }' /dev/null;
    then
        echo "Exiting test because objective value for instance ${SHORTPROBNAME} in solu file ${SOLUFILE} is inconsistent: ${OBJECTIVEVAL} vs. ${CHECKOBJECTIVEVAL}"
        exit 1
    fi
else
    OBJECTIVEVAL=""
fi

NEWSHORTPROBNAME=$(echo "${SHORTPROBNAME}" | cut -c1-25)
SHORTPROBNAME="${NEWSHORTPROBNAME}"

#define file name for temporary log file
FILENAME="${USER}.${TSTNAME}.${COUNT}"_"${SHORTPROBNAME}.${BINID}.${QUEUE}.${SETNAME}"

# if number of threads is larger than 1, add postfix
if test "${THREADS}" -gt 1
then
    FILENAME="${FILENAME}-t${THREADS}"
fi

# if seed is positive, add postfix
if test "${SEED}" -gt 0
then
    FILENAME="${FILENAME}-s${SEED}"
fi

# if permutation is positive, add postfix
if test "${PERM}" -gt 0
then
    FILENAME="${FILENAME}-p${PERM}"
fi

SKIPINSTANCE="false"
# in case we want to continue we check if the job was already performed
if test "${CONTINUE}" = "true" && test -e "${OUTPUTDIR}/${FILENAME}".out
then
    echo "skipping file ${INSTANCE} due to existing output file ${OUTPUTDIR}/${FILENAME}.out"
    SKIPINSTANCE="true"
fi

# configure global names TMPFILE (batch file) and SETFILE to save settings to
BASENAME="${SCIPPATH}/${OUTPUTDIR}/${FILENAME}"
TMPFILE="${BASENAME}.tmp"
SETFILE="${BASENAME}.set"

# even if we decide to skip this instance, we write the basename to the eval file
echo "${OUTPUTDIR}/${FILENAME}" >> "${EVALFILE}"
