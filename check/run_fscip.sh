#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      *
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

# The script executing fscip on one instance and producing the logfiles.
# Can be executed either locally or on a cluster node.
# Is to be invoked inside a 'check(_cluster)*.sh' script.

# absolut tolerance for checking linear constraints and objective value
LINTOL=1e-04
# absolut tolerance for checking integrality constraints
INTTOL=1e-04

TMPDIR="${CLIENTTMPDIR}/${USER}-tmpdir"

# check if tmp-path exists
if test ! -d "${TMPDIR}"
then
    echo "Creating directory ${TMPDIR} for temporary outfile"
    mkdir -p "${TMPDIR}"
fi

OUTFILE="${TMPDIR}/${BASENAME}.out"
ERRFILE="${TMPDIR}/${BASENAME}.err"
SOLFILE="${TMPDIR}/${BASENAME}.sol"
SETFILE="${TMPDIR}/${BASENAME}.prm"
SCIPSETFILE="${TMPDIR}/${BASENAME}.set"
TMPFILE="${TMPDIR}/${BASENAME}.tmp"

SETTINGS="${SOLVERPATH}/../settings/${SETNAME}.prm"
if test ! -e "${SETTINGS}"
then
    SETTINGS="${SOLVERPATH}/../../ug/settings/${SETNAME}.prm"
fi
SCIPSETTINGS="${SOLVERPATH}/../settings/${SETNAME}.set"
if test ! -e "${SCIPSETTINGS}"
then
    SCIPSETTINGS="${SOLVERPATH}/../../ug/settings/${SETNAME}.set"
fi

uname -a                            > "${OUTFILE}"
uname -a                            > "${ERRFILE}"

# function to copy back the results and delete temporary files
function cleanup {
    echo "=ready="                      >> "${OUTFILE}"
    echo                                >> "${OUTFILE}"
    echo "start moving files"           >> "${OUTFILE}"
    date                                >> "${OUTFILE}"
    echo                                >> "${OUTFILE}"

    mv "${ERRFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.err"
    if [ -d "${ERRFILE}.rr" ]
    then
        mv -f "${ERRFILE}.rr" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.rr"
    fi

    rm -f "${SOLFILE}"
    rm -f "${SETFILE}"
    rm -f "${SCIPSETFILE}"

    echo                                >> "${OUTFILE}"
    echo "--- FINISH ---"               >> "${OUTFILE}"
    date                                >> "${OUTFILE}"
    echo                                >> "${OUTFILE}"

    mv "${OUTFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.out"
}

# ensure TMPFILE is deleted and results are copied when exiting (normally or due to abort/interrupt)
trap cleanup EXIT

# check if the scripts runs a *.zib.de Linux host
if $(hostname -f | grep -q zib.de) && [ $(uname) == Linux ]
then
    # access optimi once to force a mount
    OPTIMIFILE=/data/optimi/optimi/kombadon/IP/miplib2003/10teams.mps.gz
    ls $OPTIMIFILE >/dev/null 2>&1

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

        if [ -f $OPTIMIFILE ]
        then
            MOUNTED=1
        else
            ((FAILED++))
            echo "optimi is not mounted yet ($OPTIMIFILE is no file), waiting 1 second"
            sleep 1
        fi
    done
fi

echo "start checking mount"                        >> "${OUTFILE}"
date                                               >> "${OUTFILE}"
echo                                               >> "${OUTFILE}"

echo "${SETTINGS} is used as UG parameters"        >> "${OUTFILE}"
echo "${SCIPSETTINGS} is used as SCIP parameters"  >> "${OUTFILE}"
# gnerate UG parameter file
cp "${SETTINGS}" "${SETFILE}"
echo "TimeLimit = ${TIMELIMIT}"                    >> "${SETFILE}"
# gnerate SCIP parameter file
cp "${SCIPSETTINGS}" "${SCIPSETFILE}"

echo "start printing some stats"                   >> "${OUTFILE}"
date                                               >> "${OUTFILE}"
echo                                               >> "${OUTFILE}"

echo                                               >> "${OUTFILE}"
if test $(uname) == Linux ; then   # -b does not work with top on macOS
    top -b -n 1 | head -n 15                       >> "${OUTFILE}"
fi
echo                                               >> "${OUTFILE}"
echo "hard time limit: ${HARDTIMELIMIT}          " >> "${OUTFILE}"
echo "hard mem limit: ${HARDMEMLIMIT}"             >> "${OUTFILE}"
echo                                               >> "${OUTFILE}"
echo "SLURM jobID: ${SLURM}_JOB_ID"                >> "${OUTFILE}"
echo                                               >> "${OUTFILE}"
echo "@01 ${FILENAME} ==========="                 >> "${OUTFILE}"
echo "@01 ${FILENAME} ==========="                 >> "${ERRFILE}"
echo "-----------------------------"               >> "${OUTFILE}"
date                                               >> "${OUTFILE}"
date                                               >> "${ERRFILE}"
echo "-----------------------------"               >> "${OUTFILE}"
date +"@03 %s"                                     >> "${OUTFILE}"
echo "@05 ${TIMELIMIT}"                            >> "${OUTFILE}"

echo                                               >> "${OUTFILE}"
echo "execute binary"                              >> "${OUTFILE}"
date                                               >> "${OUTFILE}"
echo                                               >> "${OUTFILE}"

#if we use a debugger command, we need to replace the errfile place holder by the actual err-file for logging
#and if we run on the cluster we want to use srun with CPU binding which is defined by the check_cluster script
EXECNAME="${EXECNAME/ERRFILE_PLACEHOLDER/${ERRFILE}}"
EXECNAME="${SRUN}${EXECNAME/RRTRACEFOLDER_PLACEHOLDER/${ERRFILE}}"
echo "running the following command from $(pwd):" "${EXECNAME}"  "${SETFILE}" "${FILENAME}" -sth "${THREADS}" -fsol "${SOLFILE}" -sr "${SCIPSETFILE}" -s "${SCIPSETFILE}" >> "${OUTFILE}"
eval "${EXECNAME}"  "${SETFILE}" "${FILENAME}" -sth "${THREADS}" -fsol "${SOLFILE}" -sr "${SCIPSETFILE}" -s "${SCIPSETFILE}" 2>> "${ERRFILE}"  | tee -a "${OUTFILE}"
retcode="${PIPESTATUS[0]}"
if test "${retcode}" != 0
then
    echo "${EXECNAME} returned with error code ${retcode}." >> "${ERRFILE}"
fi

echo                                >> "${OUTFILE}"
echo "call solution checker"        >> "${OUTFILE}"
date                                >> "${OUTFILE}"
echo                                >> "${OUTFILE}"

if test -e "${SOLFILE}"
then
    # translate FiberSCIP solution format into SCIP solution format and workaround infinite values for solution checker
    sed ' /\[ Final Solution \]/d;
    s/No Solution/no solution available/g;
    s/infinity/1e+20/g' "${SOLFILE}" > "${TMPFILE}"
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
