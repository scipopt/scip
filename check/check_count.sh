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

# Count feasible solutions with scip.
#
# To be invoked by Makefile 'make testcount'.

TSTNAME=${1}
BINNAME=${2}
SETNAME=${3}
BINID=${4}
TIMELIMIT=${5}
NODELIMIT=${6}
MEMLIMIT=${7}
FEASTOL=${8}
DISPFREQ=${9}
CONTINUE=${10}
LOCK=${11}
VERSION=${12}
LPS=${13}
OUTPUTDIR=${14}

SETDIR=../settings

if test ! -e ${OUTPUTDIR}
then
    mkdir ${OUTPUTDIR}
fi
if test ! -e locks
then
    mkdir locks
fi

LOCKFILE=locks/count.${TSTNAME}.${SETNAME}.${VERSION}.${LPS}.lock
RUNFILE=locks/count.${TSTNAME}.${SETNAME}.${VERSION}.${LPS}.run.${BINID}
DONEFILE=locks/count.${TSTNAME}.${SETNAME}.${VERSION}.${LPS}.done

OUTFILE=${OUTPUTDIR}/checkcount.${TSTNAME}.${BINID}.${SETNAME}.out
ERRFILE=${OUTPUTDIR}/checkcount.${TSTNAME}.${BINID}.${SETNAME}.err
RESFILE=${OUTPUTDIR}/checkcount.${TSTNAME}.${BINID}.${SETNAME}.res
TEXFILE=${OUTPUTDIR}/checkcount.${TSTNAME}.${BINID}.${SETNAME}.tex
TMPFILE=${OUTPUTDIR}/checkcount.${TSTNAME}.${BINID}.${SETNAME}.tmp
SETFILE=${OUTPUTDIR}/checkcount.${TSTNAME}.${BINID}.${SETNAME}.set

SETTINGS=${SETDIR}/${SETNAME}.set

if test "${LOCK}" = "true"
then
    if test -e ${DONEFILE}
    then
        echo skipping test due to existing done file ${DONEFILE}
        exit 1
    fi
    if test -e ${LOCKFILE}
    then
        if test -e ${RUNFILE}
        then
            echo continuing aborted run with run file ${RUNFILE}
        else
            echo skipping test due to existing lock file ${LOCKFILE}
            exit 1
        fi
    fi
    date > ${LOCKFILE}
    date > ${RUNFILE}
fi

if test ! -e ${OUTFILE}
then
    CONTINUE=false
fi

if test "${CONTINUE}" = "true"
then
    MVORCP=cp
else
    MVORCP=mv
fi

DATEINT=$(date +"%s")
if test -e ${OUTFILE}
then
    ${MVORCP} ${OUTFILE} ${OUTFILE}.old-${DATEINT}
fi
if test -e ${ERRFILE}
then
    ${MVORCP} ${ERRFILE} ${ERRFILE}.old-${DATEINT}
fi

if test "${CONTINUE}" = "true"
then
    LASTPROB=$(awk -f getlastprob.awk ${OUTFILE})
    echo Continuing benchmark. Last solved instance: ${LASTPROB}
    echo "" >> ${OUTFILE}
    echo "----- Continuing from here. Last solved: ${LASTPROB} -----" >> ${OUTFILE}
    echo "" >> ${OUTFILE}
else
    LASTPROB=""
fi

uname -a >>${OUTFILE}
uname -a >>${ERRFILE}
date >>${OUTFILE}
date >>${ERRFILE}

# we add 10% to the hard time limit and additional 10 seconds in case of small time limits
HARDTIMELIMIT=$(((TIMELIMIT + 10) + (TIMELIMIT / 10)))

# we add 10% to the hard memory limit and additional 100mb to the hard memory limit
HARDMEMLIMIT=$(((MEMLIMIT + 100) + (MEMLIMIT / 10)))
HARDMEMLIMIT=$((HARDMEMLIMIT * 1024))

echo "hard time limit: ${HARDTIMELIMIT} s" >>${OUTFILE}
echo "hard mem limit: ${HARDMEMLIMIT} k" >>${OUTFILE}

for i in $(cat testset/${TSTNAME}.test) DONE
do
    if test "${i}" = "DONE"
    then
        date > ${DONEFILE}
        break
    fi

    if test "${LASTPROB}" = ""
    then
        LASTPROB=""
        if test -f ${i}
        then
            echo @01 ${i} ===========
            echo @01 ${i} ===========                >> ${ERRFILE}
            echo > ${TMPFILE}
            echo set emphasis benchmark            >> ${TMPFILE} # avoid switching to dfs etc. - better abort with memory error; this has to be first
            if test ${SETTINGS} != "default"
            then
                echo set load ${SETTINGS}            >> ${TMPFILE}
            else
                echo set emphasis count            >> ${TMPFILE}
            fi
            if test ${FEASTOL} != "default"
            then
                echo set numerics feastol ${FEASTOL} >> ${TMPFILE}
            fi
            echo set limits time ${TIMELIMIT}        >> ${TMPFILE}
            echo set limits nodes ${NODELIMIT}       >> ${TMPFILE}
            echo set limits memory ${MEMLIMIT}       >> ${TMPFILE}
            echo set timing clocktype 1            >> ${TMPFILE}
            echo set display freq ${DISPFREQ}        >> ${TMPFILE}
            if test "${LPS}" == "none"
            then
                echo set lp solvefreq -1           >> ${TMPFILE} # avoid solving LPs in case of LPS=none
            fi
            echo set save ${SETFILE}                 >> ${TMPFILE}
            echo read ${i}                           >> ${TMPFILE}
            echo count                             >> ${TMPFILE}
            echo display statistics                >> ${TMPFILE}
            echo quit                              >> ${TMPFILE}

#            if test "${LPS}" == "cpx"
#            then
#                waitcplex.sh # ??????????????????
#            fi

            echo -----------------------------
            date
            date >>${ERRFILE}
            echo -----------------------------
            date +"@03 %s"
            bash -c "ulimit -t ${HARDTIMELIMIT}; ulimit -v ${HARDMEMLIMIT}; ulimit -f 200000; ../${BINNAME} < ${TMPFILE}" 2>>${ERRFILE}
            date +"@04 %s"
            echo -----------------------------
            date
            date >>${ERRFILE}
            echo -----------------------------
            echo
            echo =ready=
        else
            echo @02 FILE NOT FOUND: ${i} ===========
            echo @02 FILE NOT FOUND: ${i} =========== >>${ERRFILE}
        fi
    else
        echo skipping ${i}
        if test "${LASTPROB}" = "${i}"
        then
            LASTPROB=""
        fi
    fi
done | tee -a ${OUTFILE}

rm -f ${TMPFILE}

date >>${OUTFILE}
date >>${ERRFILE}

if test -e ${DONEFILE}
then
    ./evalcheck_count.sh ${OUTFILE}

    if test "${LOCK}" = "true"
    then
        rm -f ${RUNFILE}
    fi
fi
