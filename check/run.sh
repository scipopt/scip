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

# The script executing EXECNAME on one instance and producing the logfiles.
# Can be executed either locally or on a cluster node.
# Is to be invoked inside a 'check(_cluster)*.sh' script.

# absolut tolerance for checking linear constraints and objective value
LINTOL=1e-04
# absolut tolerance for checking integrality constraints
INTTOL=1e-04

TMPPATH="${CLIENTTMPDIR}/${USER}-tmpdir"

# check if tmp-path exists
if test ! -d "${TMPPATH}"
then
    mkdir "${TMPPATH}"
    echo "Creating directory ${TMPPATH} for temporary outfile."
fi

OUTFILE="${TMPPATH}/${BASENAME}.out"
ERRFILE="${TMPPATH}/${BASENAME}.err"
SOLFILE="${TMPPATH}/${BASENAME}.sol"
DATFILE="${TMPPATH}/${BASENAME}.dat"
TMPFILE="${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.tmp"

uname -a                            > "${OUTFILE}"
uname -a                            > "${ERRFILE}"

# function to copy back the results and delete temporary files
function cleanup {
    mv "${OUTFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.out"
    mv "${ERRFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.err"
    if [ -d "${ERRFILE}.rr" ]
    then
        mv -f "${ERRFILE}.rr" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.rr"
    fi
    # move a possible data file
    if [ -f "${DATFILE}" ] ;
    then
        mv "${DATFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.dat"
    fi
    rm -f "${TMPFILE}"
    rm -f "${SOLFILE}"
}

# ensure TMPFILE is deleted and results are copied when exiting (normally or due to abort/interrupt)
trap cleanup EXIT

# only wait for optimi to be mounted in run.sh if you are on a recognized computer at zib
OPTHOST=$(uname -n | sed 's/.zib.de//g' | sed 's/portal//g' | tr -cd '[:alpha:]')

# file on optimi to check for
case "$OPTHOST" in
  opt | optc ) OPTIMIFILE=/nfs/optimi/QUOTAS ;;
  htccmp     ) OPTIMIFILE=/data/optimi/optimi/QUOTAS ;;
esac

# check if the scripts runs a *.zib.de host and we have recognized the host (OPTIMIFILE being set)
if $(hostname -f | grep -q zib.de) && [ -n "$OPTIMIFILE" ];
then
    # access optimi once to force a mount
    ls ${OPTIMIFILE} >/dev/null 2>&1

    # check if optimi is mounted
    MOUNTED=0

    # count number of fails and abort after 10 min to avoid an endless loop
    FAILED=0

    while [ "${MOUNTED}" -ne 1 ]
    do
        # stop if the system does not mount optimi for ~10 minutes
        if [ "${FAILED}" -eq 600 ]
        then
            exit 1
        fi

        if [ -f ${OPTIMIFILE} ] ;
        then
            MOUNTED=1
        else
            ((FAILED++))
            echo "optimi is not mounted yet, waiting 1 second" >> "${ERRFILE}"
            sleep 1
        fi
    done
fi

echo                                        >> "${OUTFILE}"
if test $(uname) == Linux                   # -b does not work with top on macOS
then
    top -b -n 1 | head -n 15                >> "${OUTFILE}"
fi
echo                                        >> "${OUTFILE}"
echo "hard time limit: ${HARDTIMELIMIT}"    >> "${OUTFILE}"
echo "hard mem limit: ${HARDMEMLIMIT}"      >> "${OUTFILE}"
echo                                        >> "${OUTFILE}"
echo "SLURM jobID: ${SLURM}_JOB_ID"         >> "${OUTFILE}"
echo                                        >> "${OUTFILE}"
echo "@01 ${FILENAME} ==========="          >> "${OUTFILE}"
echo "@01 ${FILENAME} ==========="          >> "${ERRFILE}"
echo "-----------------------------"        >> "${OUTFILE}"
date                                        >> "${OUTFILE}"
date                                        >> "${ERRFILE}"
echo "-----------------------------"        >> "${OUTFILE}"
date +"@03 %s"                              >> "${OUTFILE}"
echo "@05 ${TIMELIMIT}"                     >> "${OUTFILE}"

#if we use a debugger command, we need to replace the errfile place holder by the actual err-file for logging
#and if we run on the cluster we want to use srun with CPU binding which is defined by the check_cluster script
EXECNAME="${EXECNAME/ERRFILE_PLACEHOLDER/${ERRFILE}}"
EXECNAME="${SRUN}${EXECNAME/RRTRACEFOLDER_PLACEHOLDER/${ERRFILE}}"
echo "${EXECNAME}"            >> "${ERRFILE}"
if test -e "${TMPFILE}"
then
    eval "${EXECNAME}"            < "${TMPFILE}" 2>>"${ERRFILE}"  | tee -a "${OUTFILE}"
else
    eval "${EXECNAME}"                           2>>"${ERRFILE}"  | tee -a "${OUTFILE}"
fi
retcode=${PIPESTATUS[0]}
if test "${retcode}" != 0
then
    echo "${EXECNAME} returned with error code ${retcode}." >>"${ERRFILE}"
fi

if test -e "${SOLFILE}"
then
    # translate SCIP solution format into format for solution checker. The
    # SOLFILE format is a very simple format where in each line we have a
    # <variable, value> pair, separated by spaces.  A variable name of
    # =obj= is used to store the objective value of the solution, as
    # computed by the solver. A variable name of =infeas= can be used to
    # indicate that an instance is infeasible.
    sed ' /solution status:/d;
    s/objective value:/=obj=/g;
    s/infinity/1e+20/g;
    s/no solution available//g' "${SOLFILE}" > "${TMPFILE}"
    mv "${TMPFILE}" "${SOLFILE}"

    # check if the link to the solution checker exists
    if test -f "${CHECKERPATH}/bin/solchecker"
    then
        echo
        "${SHELL}" -c " ${CHECKERPATH}/bin/solchecker ${FILENAME} ${SOLFILE} ${LINTOL} ${INTTOL}" 2>> "${ERRFILE}" | tee -a "${OUTFILE}"
        echo
    fi
fi

date +"@04 %s"                        >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${ERRFILE}"
echo                                  >> "${OUTFILE}"
echo "=ready="                        >> "${OUTFILE}"
