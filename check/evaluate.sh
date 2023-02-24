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

# Depending on the solver used for the testrun calls corresponding 'check_*.awk'
# on the testrun files and writes the output in a .res file.
#
# Is to be invoked by 'evalcheck(_cluster).sh'

export LANG=C
export LC_NUMERIC=C

FILE="${1}"

BASENAME=$(basename "${FILE}" .out)
DIR=$(dirname ${FILE})
OUTFILE="${DIR}/${BASENAME}.out"
ERRFILE="${DIR}/${BASENAME}.err"
RESFILE="${DIR}/${BASENAME}.res"
TEXFILE="${DIR}/${BASENAME}.tex"
PAVFILE="${DIR}/${BASENAME}.pav"

# detect test set
TSTNAME=$(echo "${BASENAME}" | sed 's/check.\([a-zA-Z0-9_-]*\).*/\1/g')

# detect test used solver
SOLVER=$(echo "${BASENAME}" | sed 's/check.\([a-zA-Z0-9_-]*\).\([a-zA-Z0-9_]*\).*/\2/g')

echo "Testset ${TSTNAME}"
echo "Solver  ${SOLVER}"

#search for test file and check if we use a ttest or a test file
if [ -f "instancedata/testsets/${TSTNAME}.ttest" ]
then
    TESTFILE="instancedata/testsets/${TSTNAME}.ttest"
elif [ -f "instancedata/testsets/${TSTNAME}.test" ]
then
    TESTFILE="instancedata/testsets/${TSTNAME}.test"
elif [ -f "testset/${TSTNAME}.ttest" ]
then
    TESTFILE="testset/${TSTNAME}.ttest"
elif [ -f "testset/${TSTNAME}.test" ]
then
    TESTFILE="testset/${TSTNAME}.test"
fi

# call method to obtain solution file
# defines the following environment variable: SOLUFILE
. ./configuration_solufile.sh "${TSTNAME}"

echo "${OUTFILE}"

# the variable AWKARGS needs to be without quotation marks in all awk calls!
if test "${SOLVER}" = "gurobi_cl"
then
    awk -f check_gurobi.awk -v "TEXFILE=${TEXFILE}" -v "PAVFILE=${PAVFILE}" ${AWKARGS} "${TESTFILE}" "${SOLUFILE}" "${OUTFILE}" | tee "${RESFILE}"
elif test  "${SOLVER}" = "cplex"
then
    awk -f check_cplex.awk -v "TEXFILE=${TEXFILE}" ${AWKARGS} "${SOLUFILE}" "${OUTFILE}" | tee "${RESFILE}"
elif test  "${SOLVER}" = "xpress"
then
    awk -f check_xpress.awk -v "TEXFILE=${TEXFILE}" ${AWKARGS} "${SOLUFILE}" "${OUTFILE}" | tee "${RESFILE}"
elif test  "${SOLVER}" = "mosek"
then
    awk -f check_mosek.awk -v "TEXFILE=${TEXFILE}" ${AWKARGS} "${SOLUFILE}" "${OUTFILE}" | tee "${RESFILE}"
elif test  "${SOLVER}" = "glpk"
then
    awk -f check_glpk.awk -v "TEXFILE=${TEXFILE}" ${AWKARGS} "${SOLUFILE}" "${OUTFILE}" | tee "${RESFILE}"
elif test  "${SOLVER}" = "symphony"
then
    awk -f check_symphony.awk -v "TEXFILE=${TEXFILE}" ${AWKARGS} "${SOLUFILE}" "${OUTFILE}" | tee "${RESFILE}"
elif test  "${SOLVER}" = "cbc"
then
    awk -f check_cbc.awk -v "TEXFILE=${TEXFILE}" ${AWKARGS} "${TESTFILE}" "${SOLUFILE}" "${OUTFILE}" | tee "${RESFILE}"
    # we should not check for SOLVER = scip here, because check.awk needs also to be called for examples with other names
else
    awk -f check.awk -v "TEXFILE=${TEXFILE}" -v "PAVFILE=${PAVFILE}" -v "ERRFILE=${ERRFILE}" ${AWKARGS} "${TESTFILE}" "${SOLUFILE}" "${OUTFILE}" | tee "${RESFILE}"
fi
