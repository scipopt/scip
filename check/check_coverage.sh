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

# Start a local scip coverage report.
#
# To be invoked by Makefile 'make coverage'.

BINNAME="${1}"
BINID="${2}"
MEMLIMIT="${3}"
VERSION="${4}"
LPS="${5}"
OBJDIR="${6}"

# settings relative to results folder
SETDIR=../coverage/settings

if ! test -e check/results
then
    mkdir -p check/results
fi

if ! test -e obj/O.linux.x86_64.gnu.gcov/lib/scip/src
then
    ln -s ../../../../src obj/O.linux.x86_64.gnu.gcov/lib/scip
fi
if ! test -e obj/O.linux.x86_64.gnu.gcov/lib/tclique/src
then
    ln -s ../../../../src obj/O.linux.x86_64.gnu.gcov/lib/tclique
fi
if ! test -e obj/O.linux.x86_64.gnu.gcov/lib/dijkstra/src
then
    ln -s ../../../../src obj/O.linux.x86_64.gnu.gcov/lib/dijkstra
fi
if ! test -e obj/O.linux.x86_64.gnu.gcov/lib/blockmemshell/src
then
    ln -s ../../../../src obj/O.linux.x86_64.gnu.gcov/lib/blockmemshell
fi
if test ! -e obj/O.linux.x86_64.gnu.gcov/lib/xml/src
then
    ln -s ../../../../src obj/O.linux.x86_64.gnu.gcov/lib/xml
fi

if ! test -e coverage/gcov
then
    mkdir -p coverage/gcov
fi

lcov -d obj/O.linux.x86_64.gnu.gcov -z

TESTSETS=(coverage)
SETTINGS=(default allaggr bivariate convertinttobin fullstrongbfs heurlprows leastinf nlpdiving1-relprop nlpdiving3 nlpdivingsolvesubmip presolaggr presoloff randomhybrid sepaaggr undercover2 allfull
    cgmip feaspump20 heuraggr intobj mostinf nlpdiving2 nlpdiving4 oddcycle presolfast pscost rensnlp undercover1 nologicor setppclifting proporbit oddcycleliftheur oddcyclelift indicatoralterlp
    indicatorlogicor indicatorsepa cgmipstrong cgmipviol dynamic alldisp uct_breadth_dualval cloud)

TESTSETS=(shortshort)
SETTINGS=(default)

cd check;

for TSTNAME in ${TESTSETS[@]}
do
    for SETNAME in ${SETTINGS[@]}
    do

    TIMELIMIT=180
    NODELIMIT=2100000000
    DISPFREQ=10000

    OUTFILE="results/check.${TSTNAME}.${BINID}.${SETNAME}.out"
    ERRFILE="results/check.${TSTNAME}.${BINID}.${SETNAME}.err"
    RESFILE="results/check.${TSTNAME}.${BINID}.${SETNAME}.res"
    TEXFILE="results/check.${TSTNAME}.${BINID}.${SETNAME}.tex"
    TMPFILE="results/check.${TSTNAME}.${BINID}.${SETNAME}.tmp"
    SETFILE="results/check.${TSTNAME}.${BINID}.${SETNAME}.set"

    SETTINGS="${SETDIR}/${SETNAME}.set"

    DATEINT=$(date +"%s")
    if test -e "${OUTFILE}"
    then
        mv "${OUTFILE}" "${OUTFILE}.old-${DATEINT}"
    fi
    if test -e "${ERRFILE}"
    then
        mv "${ERRFILE}" "${ERRFILE}.old-${DATEINT}"
    fi

    uname -a >> "${OUTFILE}"
    uname -a >> "${ERRFILE}"
    date     >> "${OUTFILE}"
    date     >> "${ERRFILE}"

    # we add 10% to the hard time limit and additional 10 seconds in case of small time limits
    HARDTIMELIMIT=$(((TIMELIMIT + 10) + (TIMELIMIT / 10)))

    # we add 10% to the hard memory limit and additional 100mb to the hard memory limit
    HARDMEMLIMIT=$(((MEMLIMIT + 1000) + (MEMLIMIT / 10)))
    HARDMEMLIMIT=$((HARDMEMLIMIT * 1024))

    echo "hard time limit: "${HARDTIMELIMIT}" s"          >> "${OUTFILE}"
    echo "hard mem limit: ${HARDMEMLIMIT} k"              >> "${OUTFILE}"

    for i in $(cat "testset/${TSTNAME}.test")
    do
        if test -f "${i}"
        then
            echo "@01 ${i} ==========="
            echo "@01 ${i} ==========="                   >> "${ERRFILE}"
            echo                                          > "${TMPFILE}"
	    echo "set emphasis benchmark"                 >> "${TMPFILE}" # avoid switching to dfs etc. - better abort with memory error; this has to be first
            if test "${SETNAME}" != "default"
            then
                echo "set load ${SETTINGS}"               >>  "${TMPFILE}"
            fi
            echo "set limits time ${TIMELIMIT}"           >> "${TMPFILE}"
            echo "set limits nodes ${NODELIMIT}"          >> "${TMPFILE}"
            echo "set limits memory ${MEMLIMIT}"          >> "${TMPFILE}"
            echo "set lp advanced threads 1"              >> "${TMPFILE}"
            echo "set timing clocktype 1"                 >> "${TMPFILE}"
            echo "set display verblevel 4"                >> "${TMPFILE}"
            echo "set display freq ${DISPFREQ}"           >> "${TMPFILE}"
            echo "set visual vbcfilename vbctest.vbc"     >> "${TMPFILE}"
            echo "set save ${SETFILE}"                    >> "${TMPFILE}"
            echo "read ${i}"                              >> "${TMPFILE}"
            echo "write problem cipreadparsetest.cip"     >> "${TMPFILE}"
            echo "write problem gamswritetest.gms"        >> "${TMPFILE}"
            echo "presolve"                               >> "${TMPFILE}"
            echo "write transproblem gamswritetest.gms"   >> "${TMPFILE}"
            echo "write transproblem ppmtest.ppm"         >> "${TMPFILE}"
            echo "write transproblem pbmtest.pbm"         >> "${TMPFILE}"
            echo "write transproblem ccgtest.ccg"         >> "${TMPFILE}"
            echo "optimize"                               >> "${TMPFILE}"
            echo "display statistics"                     >> "${TMPFILE}"
            echo "write solution soltest.sol"             >> "${TMPFILE}"
            echo "checksol"                               >> "${TMPFILE}"
            echo "quit"                                   >> "${TMPFILE}"
            echo "-----------------------------"
            date
            date                                          >> "${ERRFILE}"
            echo "-----------------------------"
            date +"@03 %s"
            bash -c " ulimit -t ${HARDTIMELIMIT} s; ulimit -v ${HARDMEMLIMIT} k; ulimit -f 200000; ../${BINNAME} < ${TMPFILE}" 2>> "${ERRFILE}"
            date +"@04 %s"
            echo "-----------------------------"
            date
            date                                          >> "${ERRFILE}"
            echo "-----------------------------"
            echo
            echo "=ready="

            cp soltest.sol soltest.fix

            echo "@01 ${i} ==========="
            echo "@01 ${i} ==========="                   >> "${ERRFILE}"
            echo > "${TMPFILE}"
	    echo "set emphasis benchmark"                 >> "${TMPFILE}" # avoid switching to dfs etc. - better abort with memory error; this has to be first
            if test "${SETNAME}" != "default"
            then
                echo "set load ${SETTINGS}"               >> "${TMPFILE}"
            fi
            echo "set limits time ${TIMELIMIT}"           >> "${TMPFILE}"
            echo "set limits nodes ${NODELIMIT}"          >> "${TMPFILE}"
            echo "set limits memory ${MEMLIMIT}"          >> "${TMPFILE}"
            echo "set lp advanced threads 1"              >> "${TMPFILE}"
            echo "set timing clocktype 1"                 >> "${TMPFILE}"
            echo "set display verblevel 4"                >> "${TMPFILE}"
            echo "set display freq ${DISPFREQ}"           >> "${TMPFILE}"
            echo "set save ${SETFILE}"                    >> "${TMPFILE}"
            echo "read cipreadparsetest.cip"              >> "${TMPFILE}"
            echo "read soltest.sol"                       >> "${TMPFILE}"
            echo "read soltest.fix"                       >> "${TMPFILE}"
            echo "optimize"                               >> "${TMPFILE}"
            echo "display statistics"                     >> "${TMPFILE}"
            echo "checksol"                               >> "${TMPFILE}"
            echo "quit"                                   >> "${TMPFILE}"
            echo "-----------------------------"
            date
            date                                          >> "${ERRFILE}"
            echo "-----------------------------"
            date +"@03 %s"
            bash -c " ulimit -t ${HARDTIMELIMIT} s; ulimit -v ${HARDMEMLIMIT} k; ulimit -f 200000; ../${BINNAME} < "${TMPFILE}"" 2>>"${ERRFILE}"
            date +"@04 %s"
            echo "-----------------------------"
            date
            date                                          >> "${ERRFILE}"
            echo "-----------------------------"
            echo
            echo "=ready="

            echo "@01 ${i} ==========="
            echo "@01 ${i} ==========="                   >> "${ERRFILE}"
            echo                                           > "${TMPFILE}"
	    echo "set emphasis benchmark"                 >> "${TMPFILE}" # avoid switching to dfs etc. - better abort with memory error; this has to be first
            if test "${SETNAME}" != "default"
            then
                echo "set load ${SETTINGS}"               >>  "${TMPFILE}"
            fi
            echo "set limits time ${TIMELIMIT}"           >> "${TMPFILE}"
            echo "set limits nodes ${NODELIMIT}"          >> "${TMPFILE}"
            echo "set limits memory ${MEMLIMIT}"          >> "${TMPFILE}"
            echo "set lp advanced threads 1"              >> "${TMPFILE}"
            echo "set timing clocktype 1"                 >> "${TMPFILE}"
            echo "set display verblevel 4"                >> "${TMPFILE}"
            echo "set display freq ${DISPFREQ}"           >> "${TMPFILE}"
            echo "set save ${SETFILE}"                    >> "${TMPFILE}"
            echo "read ${i}"                              >> "${TMPFILE}"
            echo "presolve"                               >> "${TMPFILE}"
            echo "read soltest.sol"                       >> "${TMPFILE}"
            echo "read soltest.fix"                       >> "${TMPFILE}"
            if test "${i}" = "instances/MIP/bell5.mps"
            then
                echo "read bell5.sol"                     >> "${TMPFILE}"
            fi
            echo "presolve"                               >> "${TMPFILE}"
            echo "write transproblem lpwritetest.lp"      >> "${TMPFILE}"
            echo "write transproblem mpswritetest.mps"    >> "${TMPFILE}"
            echo "write transproblem opbtest.opb"         >> "${TMPFILE}"
            echo "write transproblem piptest.pip"         >> "${TMPFILE}"
            echo "write transproblem rlptest.rlp"         >> "${TMPFILE}"
            echo "read mpswritetest.mps"                  >> "${TMPFILE}"
            echo "read lpwritetest.lp"                    >> "${TMPFILE}"
            echo "read rlptest.rlp"                       >> "${TMPFILE}"
            echo "read opbtest.opb"                       >> "${TMPFILE}"
            echo "read piptest.pip"                       >> "${TMPFILE}"
            echo "display statistics"                     >> "${TMPFILE}"
            echo "checksol"                               >> "${TMPFILE}"
            echo "quit"                                   >> "${TMPFILE}"
            echo "-----------------------------"
            date
            date                                          >> "${ERRFILE}"
            echo "-----------------------------"
            date +"@03 %s"
            bash -c " ulimit -t ${HARDTIMELIMIT} s; ulimit -v ${HARDMEMLIMIT} k; ulimit -f 200000; ../"${BINNAME}" < "${TMPFILE}"" 2>> "${ERRFILE}"
            date +"@04 %s"
            echo "-----------------------------"
            date
            date                                          >> "${ERRFILE}"
            echo "-----------------------------"
            echo
            echo "=ready="

            rm -f cipreadparsetest.cip
            rm -f gamswritetest.gms
            rm -f ppmtest.ppm
            rm -f pbmtest.pbm
            rm -f ccgtest.ccg
            rm -f soltest.sol
            rm -f soltest.fix
            rm -f vbctest.vbc
            rm -f piptest.pip
            rm -f opbtest.opb
            rm -f rlptest.rlp
            rm -f lpwritetest.lp
            rm -f mpswritetest.mps
        else
            echo "@02 FILE NOT FOUND: ${i} ==========="
            echo "@02 FILE NOT FOUND: ${i} ===========" >> "${ERRFILE}"
        fi
    done | tee -a "${OUTFILE}"

    rm -f "${TMPFILE}"

    date >> "${OUTFILE}"
    date >> "${ERRFILE}"

    ./evalcheck.sh "${OUTFILE}"
    done
done

cd ..;

OUTFILE="results/check.commands.${BINID}.out"
ERRFILE="results/check.commands.${BINID}.err"


bash -c " ulimit -t ${HARDTIMELIMIT} s; ulimit -v ${HARDMEMLIMIT} k; ulimit -f 200000; "${BINNAME}" < coverage/commands.bat''" 2>> "${ERRFILE}" | tee -a "${OUTFILE}"

#make LPS=spx OPT=gcov TIME=180 testexamples

#cd ../gcg-coverage;
#make OPT=gcov LPS=spx TEST=cpmp50s TIME=600 SETTINGS=presolve
#make OPT=gcov LPS=spx TEST=cpmp50s TIME=600 SETTINGS=presolve test
#cd ../scip-coverage

lcov -d obj/O.linux.x86_64.gnu.gcov -c >coverage/gcov/scip.capture
genhtml -o coverage/html coverage/gcov/scip.capture
