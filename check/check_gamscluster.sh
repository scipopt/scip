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

# Start a testrun with gams on the cluster.
#
# To be invoked by Makefile 'make testgamscluster'.

TSTNAME="${1}"
GAMSBIN="${2}"
SOLVER=${3^^}
SETNAME="${4}"
TIMELIMIT="${6}"
NODELIMIT="${7}"
MEMLIMIT="${8}"
GAPLIMIT=${9:-0}
THREADS="${10}"
CONTINUE=${11,,}
CONVERTSCIP="${12}"
QUEUETYPE="${13}"
QUEUE="${14}"
PPN="${15}"
CLIENTTMPDIR="${16}"
NOWAITCLUSTER="${17}"
EXCLUSIVE="${18}"
# set this to 1 if you want the scripts to (try to) pass a best known primal bound (from .solu file) to the GAMS solver
SETCUTOFF="${19}"

# set this to 1 if you want the scripts to (try to) pass a best known solution (from .gdx file) to the GAMS solver
PASSSTARTSOL=0

# set this to 1 to keep solutions in .gdx files
KEEPSOLS=0

# set this to 1 to run the solver through Examiner2
EXAMINER=0

# check all variables defined
if [ -z "${EXCLUSIVE}" ]
then
    echo "Skipping test since not all variables are defined."
    exit 1;
fi

# check if queuetype has been defined
if test "${QUEUETYPE}" = ""
then
    echo "Skipping test since the queuetype has not been defined."
    exit
fi

# if we run locally, use hostname as "queue"-name
if test "${QUEUETYPE}" = "local"
then
    QUEUE="${HOSTNAME}"
fi

# check if gams system is available
if ! which "${GAMSBIN}" > /dev/null 2>&1
then
    echo "No GAMS system available: ${GAMSBIN} does not work. Abort."
    exit 1
fi
BINNAME=$(echo "${GAMSBIN}" | sed -e 's/[^A-Za-z0-9_.]/_/g')

# check if gams solver has been specified
if test -z "${SOLVER}"
then
    echo "GAMSSOLVER not specified. Abort."
    exit 1
fi

# check if testset file is available
if test ! -r "testset/${TSTNAME}.test"
then
    echo "Testset file testset/${TSTNAME}.test not available. Abort."
    exit 1
fi

if test ! -e results
then
    mkdir results
fi
if test ! -e locks
then
    mkdir locks
fi

EVALFILE="results/check.${TSTNAME}.${BINNAME}.${SOLVER}.${QUEUE}.${SETNAME}.eval"
SETFILE="results/check.${TSTNAME}.${BINNAME}.${SOLVER}.${QUEUE}.${SETNAME}.set"
SCHFILE="results/check.${TSTNAME}.${BINNAME}.${SOLVER}.${QUEUE}.${SETNAME}.sch"
OPTFILE="$(pwd)/results/check.${TSTNAME}.${BINNAME}.${SOLVER}.${QUEUE}.${SETNAME}.opt"
GMSDIR="$(pwd)/results/check.${TSTNAME}.${BINNAME}.${SOLVER}.${QUEUE}.${SETNAME}.gms"
SOLDIR="$(pwd)/results/check.${TSTNAME}.${BINNAME}.${SOLVER}.${QUEUE}.${SETNAME}.sol"

# additional environment variables needed by finishgamscluster.sh at the end (or when trap is setup)
export GMSDIR="${GMSDIR}"
export EVALFILE="${EVALFILE}"

echo > "${EVALFILE}"

# we add 50% to the time limit and additional 10 seconds in case of small time limits
HARDTIMELIMIT=$(((TIMELIMIT + 10) + (TIMELIMIT / 2)))

# echo "hard time limit: ${HARDTIMELIMIT} s" >> ${OUTFILE}

# set pf4=0 to get no default upper bounds on integer variables
# set domlim to infinity to not stop on function evaluation errors
## do not use scratch files (solvelink=5) if possible, GAMS resets to solvelink=2 if not supported by solver
# set logoption=3 to get default output to stdout
# listing file: append mode, print step summary, disable solution printing, disable rows and columns output, disable page control
GAMSOPTS="pf4=0 domlim=9999999" # solvelink=5
GAMSOPTS="${GAMSOPTS} logoption=3 stepsum=1 solprint=0 limcol=0 limrow=0 pc=2 pw=255" #appendout=1
GAMSOPTS="${GAMSOPTS} reslim=${TIMELIMIT}"
GAMSOPTS="${GAMSOPTS} nodlim=${NODELIMIT}"
GAMSOPTS="${GAMSOPTS} workspace=${MEMLIMIT}"
GAMSOPTS="${GAMSOPTS} optcr=${GAPLIMIT}"
GAMSOPTS="${GAMSOPTS} threads=${THREADS}"

# set SBB option to overwrite NLP solver status files in each node instead of appending
if test "${SOLVER}" = SBB
then
    GAMSOPTS="${GAMSOPTS} integer4=2"
fi

# create directory for solutions, if KEEPSOLS is true
if test "${KEEPSOLS}" = 1
then
    mkdir -p "${SOLDIR}"
    GAMSOPTS="${GAMSOPTS} gdxcompress=1"
fi

# setup solver option file
# copy optionfile <solvername>.opt to ${OPTFILE}
if test "${SETNAME}" != "default"
then
    SETDIR=$(cd ../settings ; pwd)
    if test -f "${SETDIR}/${SETNAME}.gamsset"
    then
        SETTINGS="${SETDIR}/${SETNAME}.gamsset"
        cp "${SETTINGS}" "${OPTFILE}"
    else
        echo "${m} settings file ${SETDIR}/${SETNAME}.gamsset not found"
        exit 1
    fi
fi

# add information on solver and limits for eval script
echo "SOLVER,${SOLVER}" >> "${EVALFILE}"
echo "TIMELIMIT,${TIMELIMIT}" >> "${EVALFILE}"
echo "NODELIMIT,${NODELIMIT}" >> "${EVALFILE}"
echo "GAPLIMIT,${GAPLIMIT}" >> "${EVALFILE}"
echo "SETTINGS,${SETNAME}" >> "${EVALFILE}"

# save gams command line options in ${SETFILE}
echo "${GAMSOPTS}" > "${SETFILE}"

# check for SCIP for converting test instances into GAMS files
# set CONVERTSCIP to empty if no SCIP found or feature is disabled
if test -z "${CONVERTSCIP}"
then
    if test -x ../bin/scip
    then
        CONVERTSCIP=../bin/scip
    else
        CONVERTSCIP=
    fi
elif test "${CONVERTSCIP}" = no
then
    CONVERTSCIP=
else
    if test ! -x "${CONVERTSCIP}"
    then
        echo "${CONVERTSCIP} for converting into gams format is not an executable file. Abort."
        exit 1
    fi
fi

# get name of solver executable
solverexe=$(which "${GAMSBIN}")
solverexe=$(dirname "${solverexe}")
solverexe=$(grep -A 2 ^${SOLVER} ${solverexe}/"gmscmpun.txt" | tail -1)
if test -z "${solverexe}"
then
    echo "${SOLVER} does not seem to be a GAMS solver (does not appear in gmscmpun.txt). Abort."
    exit 1
fi

# if run locally, run schulz to make sure solvers do not run forever
# also setup what happens on exit
if test "${QUEUETYPE}" = local ; then
    # signal 2 (sigint, ^C) when 5 seconds above timelimit
    # signal 1 (sighup) when at hard timelimit
    # signal 9 (sigkill) when at hard timelimit + 60 seconds
    # set sleepseconds small enough so that we check at least once between ${TIMELIMIT} and ${HARDTIMELIMIT}, which have difference at least ${TIMELIMIT}/10,
    # and at least once between ${HARDTIMELIMIT} and ${HARDTIMELIMIT} + 60
    (( sleepsec = TIMELIMIT > 600 ? 60 : (TIMELIMIT > 30 ? TIMELIMIT / 10 - 2 : 1) ))
    ./schulz.sh "^${solverexe}" "$((TIMELIMIT + 5)):${HARDTIMELIMIT}:$((HARDTIMELIMIT + 60))" "2:1:9" "${sleepsec}" > "${SCHFILE}" 2>&1 &
    schulzpid=$!

    # kill schulz on exit and call finishgamscluster script
    trap "echo 'Finishing up.'; kill ${schulzpid}; ./finishgamscluster.sh" EXIT
fi

# if cutoff should be passed, check for solu file
if test "${SETCUTOFF}" = 1 || test "${SETCUTOFF}" = true ; then
    # call method to obtain solution file
    # defines the following environment variable: SOLUFILE
    . ./configuration_solufile.sh "${TSTNAME}"

    if test "${SOLUFILE}" = ""
    then
        echo "Warning: SETCUTOFF=1 set, but no .solu file available"
        SETCUTOFF=0
    fi
fi

#define clusterqueue, which might not be the QUEUE, cause this might be an alias for a bunch of QUEUEs
CLUSTERQUEUE="${QUEUE}"

NICE=""
ACCOUNT="optimi_integer"

CONSTRAINT=""
if test "${CLUSTERQUEUE}" = "Gold6338"
then
    CONSTRAINT="Gold6338"
    CLUSTERQUEUE="big"
elif test "${CLUSTERQUEUE}" = "Gold6342"
then
    CONSTRAINT="Gold6342"
    CLUSTERQUEUE="big"
elif test "${CLUSTERQUEUE}" = "M640v2"
then
    CONSTRAINT="Gold5222"
    CLUSTERQUEUE="opt_int"
elif test "${CLUSTERQUEUE}" = "M640"
then
    CONSTRAINT="Gold5122"
    CLUSTERQUEUE="opt_int"
fi

# check if the slurm blades should be used exclusively
if test "${EXCLUSIVE}" = "true"
then
    EXCLUSIVE=" --exclusive"
else
    EXCLUSIVE=""
fi

# counter to define file names for a test set uniquely
COUNT=1
# dependencies of finish-job, if run through slurm
FINISHDEPEND=afterany

for i in $(cat "testset/${TSTNAME}.test") DONE
do
    if test "${i}" = "DONE"
    then
        break
    fi
    if test -f "${i}"
    then
        # the cluster queue has an upper bound of 2000 jobs; if this limit is
        # reached the submitted jobs are dumped; to avoid that we check the total
        # load of the cluster and wait until it is save (total load not more than
        # 1900 jobs) to submit the next job.
        if test "${NOWAITCLUSTER}" != "1"
        then
            ./waitcluster.sh 1600 "${QUEUE}" 200
        fi

        SHORTFILENAME=$(basename "${i}" .gz)
        SHORTFILENAME=$(basename "${SHORTFILENAME}" .mps)
        SHORTFILENAME=$(basename "${SHORTFILENAME}" .lp)
        SHORTFILENAME=$(basename "${SHORTFILENAME}" .opb)
        SHORTFILENAME=$(basename "${SHORTFILENAME}" .gms)
        SHORTFILENAME=$(basename "${SHORTFILENAME}" .pip)
        SHORTFILENAME=$(basename "${SHORTFILENAME}" .zpl)

        FILENAME="${USER}.${TSTNAME}.${COUNT}"_"${SHORTFILENAME}.${BINNAME}.${SOLVER}.${QUEUE}.${SETNAME}"
        BASENAME="results/${FILENAME}"

        TMPFILE="${BASENAME}.tmp"

        echo "${BASENAME}" >> "${EVALFILE}"

        COUNT=$((COUNT + 1))

    # in case we want to continue we check if the job was already performed
    if test "${CONTINUE}" != "false"
    then
        if test -e "results/${FILENAME}.out"
        then
            echo "skipping file ${i} due to existing output file ${FILENAME}.out"
            continue
        fi
    fi

    GMSFILE=$(basename "${i}")
    if test "${i:0:1}" == "/" ; then
        INPUTDIR=$(dirname "${i}")
    else
        INPUTDIR=$(pwd)/$(dirname "${i}")
    fi
    case ${GMSFILE} in
        *.gms | *.gms.gz )
            ;;
        *.gms.z )
            echo "Temporarily decompressing ${i}."
            GMSFILE=${GMSFILE/%.gz/}
            GMSFILE=${GMSFILE/%.z/}
            INPUTDIR="${GMSDIR}"
            mkdir -p "${INPUTDIR}"
            gunzip -c "${i}" > "${INPUTDIR}/${GMSFILE}"
            ;;
        *.gms.bz2 )
            echo "Temporarily decompressing ${i}."
            GMSFILE=${GMSFILE/.bz2/}
            INPUTDIR="${GMSDIR}"
            mkdir -p "${INPUTDIR}"
            bunzip2 -c "${i}" > "${INPUTDIR}/${GMSFILE}"
            ;;
        * )
            if test -n "${CONVERTSCIP}"
            then
                echo "Attempt converting ${i} to GAMS format using ${CONVERTSCIP}."
                INPUTDIR="${GMSDIR}"
                mkdir -p "${INPUTDIR}"
                # construct gams file name as basename of model, compressor extensions stripped, everything starting at last '.' removed, and .gms appened
                GMSFILE="$(echo "${GMSFILE}" | sed -e 's/\.bz2$//' -e 's/\.gz$//' -e 's/\.z$//' -e 's/\.[^.]\{1,\}$//').gms"
                "${CONVERTSCIP}" -s "" -c "set reading gmsreader freeints TRUE r ${i} w genprob ${INPUTDIR}/${GMSFILE} q"
            else
                echo "Warning: Instance may not be in GAMS format, but conversion is disabled. Expect trouble!"
            fi
            ;;
    esac

    # get modeltype of instance and convert into capitals
    MODTYPE=$(zgrep -i '^[^*]' "${INPUTDIR}/${GMSFILE}" | tr '[a-z]' '[A-Z]' | sed -n -e '/SOLVE.* USING/s/\(.* USING [ %]*\)\([^ %]*\)\(.*\)/\2/p')
    #echo "Modeltype: ${MODTYPE:-UNKNOWN!}"

    if test -z "${MODTYPE}"
    then
        echo "Could not recognize model type. Using SOLVER."
        MODTYPE=SOLVER
    fi

    if test "${KEEPSOLS}" = 1
    then
        GDXFILE="gdx=${SOLDIR}/${GMSFILE/%.gms*/.gdx}"
    fi

    if test "${SETCUTOFF}" = 1 || test "${SETCUTOFF}" = true
    then
        export CUTOFF=$(grep ${GMSFILE/%.gms*/} "${SOLUFILE}" | grep -v =feas= | grep -v =inf= | tail -n 1 | awk '{print $3}')
    fi

    # additional environment variables needed by rungamscluster.sh
    export BASENAME="${FILENAME}"
    export FILENAME="${i}"
    export GAMSBIN
    export GAMSOPTS
    export GMSFILE
    export INPUTDIR
    export MODTYPE
    export SOLVER
    export GDXFILE
    export CLIENTTMPDIR
    export PASSSTARTSOL
    export EXAMINER
    export SETTINGS

    case "${QUEUETYPE}" in
        srun )
            # slurm prefers to have a memory limit set
            # we add 10% to the hard memory limit and additional 100MB to the memory limit
            HARDMEMLIMIT=$(((MEMLIMIT + 100) + (MEMLIMIT / 10)))

        sbatchret=$(sbatch --job-name="${SHORTFILENAME}" --mem="${HARDMEMLIMIT}" --time="0:${HARDTIMELIMIT}" -p "${CLUSTERQUEUE}" --constraint="${CONSTRAINT}" -A "${ACCOUNT}" ${EXCLUSIVE} ${NICE} --output=/dev/null rungamscluster.sh)
        echo "${sbatchret}"
        FINISHDEPEND="${FINISHDEPEND}":$(echo "${sbatchret}" | cut -d " " -f 4)
        ;;
    qsub )
        # -V to copy all environment variables
        qsub -l walltime="${HARDTIMELIMIT}" -l nodes=1:ppn="${PPN}" -N "${SHORTFILENAME}" -V -q "${QUEUE}" -o /dev/null -e /dev/null rungamscluster.sh
        ;;
    local )
        echo "Running ${GMSFILE} locally."
        ./rungamscluster.sh
        ;;
    esac
else
    echo "input file "${i}" not found!"
    fi
done

# add job that calls finishgamscluster.sh for slurm runs
# for local runs, finishgamscluster.sh is called via the trap above, so it's also called if we interrupt with a Ctrl+C or similar
#TODO call finishgamscluster also in qsub runs
case "${QUEUETYPE}" in
    srun )
        sbatch --job-name=GAMSFINISH --mem=1000 -p "${CLUSTERQUEUE}" -A "${ACCOUNT}" --output=/dev/null -d "${FINISHDEPEND}" finishgamscluster.sh
        echo
        squeue -p "${CLUSTERQUEUE}"
        ;;
esac
