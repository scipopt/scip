#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      *
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
#
# Call via Makefile 'make testcluster' for a clusterrun with scip,
#                   'make testclustercbc'                    cbc
#                   'make testclustercpx'                    cplex
#                   'make testclustergurobi'                 gurobi
#                   'make testclusterxpress'                 xpress
#
# The queue is passed via "${QUEUE}" (possibly defined in a local makefile in scip/make/local).
#
# For each run, we can specify the number of nodes reserved for a run via "${PPN}". If tests runs
# with valid time measurements should be executed, this number should be chosen in such a way
# that a job is run on a single computer, i.e., in general, "${PPN}" should equal the number of cores
# of each computer. Of course, the value depends on the specific computer/queue.
#
# To get the result files call "./evalcheck_cluster.sh
# "${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${SETNAME}.eval" in directory check/
# This leads to result files
#  - "${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${SETNAME}.out"
#  - "${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${SETNAME}.res"
#  - "${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${SETNAME}.err"
#
# To get verbose output from Slurm, have SRUN_FLAGS="-v -v" set in your environment.

TSTNAME="${1}"
BINNAME="${2}"
SETNAMES="${3}"
BINID="${4}"
OUTPUTDIR="${5}"
TIMELIMIT="${6}"
NODELIMIT="${7}"
MEMLIMIT="${8}"
THREADS="${9}"
FEASTOL="${10}"
LPS="${11}"
DISPFREQ="${12}"
CONTINUE="${13}"
QUEUETYPE="${14}"
QUEUE="${15}"
PPN="${16}"
CLIENTTMPDIR="${17}"
NOWAITCLUSTER="${18}"
EXCLUSIVE="${19}"
PERMUTE="${20}"
SEEDS="${21}"
GLBSEEDSHIFT="${22}"
STARTPERM="${23}"
DEBUGTOOL="${24}"
REOPT="${25}"
OPTCOMMAND="${26}"
SETCUTOFF="${27}"
VISUALIZE="${28}"
CLUSTERNODES="${29}"
EXCLUDENODES="${30}"
SLURMACCOUNT="${31}"
PYTHON="${32}"
EMPHBENCHMARK="${33}"
CLOCKTYPE="${34}"
WITHCERTIFICATE="${35}"
KEEPSOL="${36}"

# check if all variables defined (by checking the last one)
if test -z "${KEEPSOL}"
then
    echo Skipping test since not all variables are defined
    echo "TSTNAME       = ${TSTNAME}"
    echo "BINNAME       = ${BINNAME}"
    echo "SETNAMES      = ${SETNAMES}"
    echo "BINID         = ${BINID}"
    echo "OUTPUTDIR     = ${OUTPUTDIR}"
    echo "TIMELIMIT     = ${TIMELIMIT}"
    echo "NODELIMIT     = ${NODELIMIT}"
    echo "MEMLIMIT      = ${MEMLIMIT}"
    echo "THREADS       = ${THREADS}"
    echo "FEASTOL       = ${FEASTOL}"
    echo "LPS           = ${LPS}"
    echo "DISPFREQ      = ${DISPFREQ}"
    echo "CONTINUE      = ${CONTINUE}"
    echo "QUEUETYPE     = ${QUEUETYPE}"
    echo "QUEUE         = ${QUEUE}"
    echo "PPN           = ${PPN}"
    echo "CLIENTTMPDIR  = ${CLIENTTMPDIR}"
    echo "NOWAITCLUSTER = ${NOWAITCLUSTER}"
    echo "EXCLUSIVE     = ${EXCLUSIVE}"
    echo "PERMUTE       = ${PERMUTE}"
    echo "SEEDS         = ${SEEDS}"
    echo "GLBSEEDSHIFT  = ${GLBSEEDSHIFT}"
    echo "STARTPERM     = ${STARTPERM}"
    echo "DEBUGTOOL     = ${DEBUGTOOL}"
    echo "REOPT         = ${REOPT}"
    echo "OPTCOMMAND    = ${OPTCOMMAND}"
    echo "SETCUTOFF     = ${SETCUTOFF}"
    echo "VISUALIZE     = ${VISUALIZE}"
    echo "CLUSTERNODES  = ${CLUSTERNODES}"
    echo "EXCLUDENODES  = ${EXCLUDENODES}"
    echo "SLURMACCOUNT  = ${SLURMACCOUNT}"
    echo "PYTHON        = ${PYTHON}"
    echo "EMPHBENCHMARK = ${EMPHBENCHMARK}"
    echo "CLOCKTYPE     = ${CLOCKTYPE}"
    echo "WITHCERTIFICATE = ${WITHCERTIFICATE}"
    echo "KEEPSOL       = ${KEEPSOL}"
    exit 1;
fi

# configure cluster-related environment variables
# defines the following environment variables: NICE, ACCOUNT, CLUSTERQUEUE, CONSTRAINT
. ./configuration_cluster.sh "${QUEUE}" "${PPN}" "${EXCLUSIVE}" "${QUEUETYPE}"

# the srun queue requires a format duration HH:MM:SS (and optionally days),
# whereas the qsub requires the memory limit in kB
if test "${QUEUETYPE}" != "qsub"
then
    TIMEFORMAT="format"
    MEMFORMAT="MB"
else
    TIMEFORMAT="sec"
    MEMFORMAT="B"
fi
# call routines for creating the result directory, checking for existence
# of passed settings, etc
# defines the following environment variables: SCIPPATH, SETTINGSLIST, SOLUFILE, HARDMEMLIMIT, DEBUGTOOLCMD, INSTANCELIST,
#                                              TIMELIMLIST, HARDTIMELIMLIST
. ./configuration_set.sh "${BINNAME}" "${TSTNAME}" "${SETNAMES}" "${TIMELIMIT}" "${TIMEFORMAT}" "${MEMLIMIT}" "${MEMFORMAT}" "${DEBUGTOOL}" "${SETCUTOFF}"

# ensure THREADS_SAFE is always available (used by AUTO mode and srun steps)
THREADS_SAFE=$(( THREADS > 0 ? THREADS : 1 ))

# query node topology to compute resource requests per instance
SINFO_CMD="sinfo --noheader -N -p ${CLUSTERQUEUE} -O"
if test -n "${CONSTRAINT}"
then
    SINFO_FILTER="grep ${CONSTRAINT}"
else
    SINFO_FILTER="cat"
fi
NODE_SOCKETS=$(${SINFO_CMD} sockets,features | ${SINFO_FILTER} | awk '{print $1}' | sort -n | head -1 | tr -d ' ')
NODE_CORES=$(${SINFO_CMD} cores,features | ${SINFO_FILTER} | awk '{print $1}' | sort -n | head -1 | tr -d ' ')

if test "${NODE_SOCKETS:-0}" -gt 0 && test "${NODE_CORES:-0}" -gt 0
then
    NODE_CPUS=$(( NODE_CORES * NODE_SOCKETS ))
else
    NODE_CPUS=1
    echo "Warning: could not determine node topology for partition ${CLUSTERQUEUE}; assuming ${NODE_CPUS} core"
fi

# compute AUTO instances-per-node
if test "${AUTO_PPN_PENDING}" -eq 1
then
    NODE_MEM_MB=$(${SINFO_CMD} memory,features | ${SINFO_FILTER} | awk '{print $1}' | sort -n | head -1 | tr -d ' ')
    if test "${NODE_MEM_MB:-0}" -le 0
    then
        echo "Warning: could not determine node memory; falling back to AUTO_PPN=1"
    else
        PPN_CPU=$(( NODE_CPUS / THREADS_SAFE ))
        PPN_MEM=$(( NODE_MEM_MB / HARDMEMLIMIT ))
        AUTO_PPN=$(( PPN_CPU < PPN_MEM ? PPN_CPU : PPN_MEM ))
        if test "${AUTO_PPN}" -lt 1
        then
            AUTO_PPN=1
        fi
    fi
    # apply user-specified upper bound
    if test "${PPN}" -gt 0 && test "${AUTO_PPN}" -gt "${PPN}"
    then
        AUTO_PPN="${PPN}"
    fi
fi

# compute NODE_FLAGS once (CLUSTERNODES/EXCLUDENODES don't change per instance)
if test "${CLUSTERNODES}" = "all" && test "${EXCLUDENODES}" = "none"
then
    NODE_FLAGS=""
elif test "${CLUSTERNODES}" != "all" && test "${EXCLUDENODES}" = "none"
then
    NODE_FLAGS="-w ${CLUSTERNODES}"
elif test "${CLUSTERNODES}" = "all" && test "${EXCLUDENODES}" != "none"
then
    NODE_FLAGS="-x ${EXCLUDENODES}"
else
    NODE_FLAGS="-w ${CLUSTERNODES} -x ${EXCLUDENODES}"
fi

# at the first time, some files need to be initialized. set to "" after the innermost loop
# finished the first time
INIT="true"

# counter to define file names for a test set uniquely
COUNT=0

# auto mode: compute srun prefix for within-job steps
if test "${AUTO_PPN_PENDING}" -eq 1
then
    AUTO_SRUN="srun --exact -n 1 -c ${THREADS_SAFE} --mem=${HARDMEMLIMIT} --cpu_bind=verbose,cores --propagate=STACK ${SRUN_FLAGS}"
fi

# auto mode batch accumulator
AUTO_BATCH_COUNT=0          # number of instances accumulated so far
AUTO_BATCH_TIME=""          # max HARDTIMELIMIT for the group
AUTO_BATCH_JOBNAME=""       # JOBNAME of first instance in group (for sbatch --job-name)
AUTO_BATCH_ARGS=()          # positional arguments for run_group.sh (14 per instance)
AUTO_BATCH_NAMES=""         # space-separated instance names for echo

# flush the current auto mode batch as a single exclusive sbatch job
flush_auto_batch() {
    if test "${AUTO_BATCH_COUNT}" -eq 0
    then
        return
    fi
    GROUP_MEM=$(( AUTO_BATCH_COUNT * HARDMEMLIMIT ))
    echo sbatch --job-name="${AUTO_BATCH_JOBNAME}" --constraint="${CONSTRAINT}" --mem="${GROUP_MEM}" -n "${AUTO_BATCH_COUNT}" -c "${THREADS_SAFE}" --distribution=cyclic:cyclic -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${AUTO_BATCH_TIME}" --cpu-freq=medium-medium:Performance --exclusive ${NODE_FLAGS} --output=/dev/null run_group.sh
    echo "instances:${AUTO_BATCH_NAMES}"
    sbatch --job-name="${AUTO_BATCH_JOBNAME}" --constraint="${CONSTRAINT}" --mem="${GROUP_MEM}" -n "${AUTO_BATCH_COUNT}" -c "${THREADS_SAFE}" --distribution=cyclic:cyclic -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${AUTO_BATCH_TIME}" --cpu-freq=medium-medium:Performance --exclusive ${NODE_FLAGS} --output=/dev/null run_group.sh "${AUTO_SRUN}" "${AUTO_BATCH_COUNT}" "${AUTO_BATCH_ARGS[@]}"
    # reset accumulator
    AUTO_BATCH_COUNT=0
    AUTO_BATCH_JOBNAME=""
    AUTO_BATCH_TIME=""
    AUTO_BATCH_ARGS=()
    AUTO_BATCH_NAMES=""
}

# loop over permutations
# loop over testset
for idx in "${!INSTANCELIST[@]}"
do
    # retrieve instance and timelimits from arrays set in the configuration_set.sh script
    INSTANCE=${INSTANCELIST[${idx}]}
    TIMELIMIT=${TIMELIMLIST[${idx}]}
    HARDTIMELIMIT=${HARDTIMELIMLIST[${idx}]}
    # increase the index for the instance tried to solve, even if the filename does not exist
    COUNT=$((COUNT + 1))

    # we need the DONE keyword for the check.sh script to automatically run evalcheck, here it is not needed
    if test "${INSTANCE}" = "DONE"
    then
        continue
    fi

    # run different random seeds
    for ((s = 0; ${s} <= ${SEEDS}; s++))
    do
        # permute transformed problem
        for ((p = 0; ${p} <= ${PERMUTE}; p++))
        do
            # the cluster queue has an upper bound of 2000 jobs; if this limit is
            # reached the submitted jobs are dumped; to avoid that we check the total
            # load of the cluster and wait until it is safe (total load not more than
            # 1600 jobs) to submit the next job.
            if test "${NOWAITCLUSTER}" -eq "0" && test "${QUEUETYPE}" = "qsub"
            then
                ./waitcluster.sh 1600 "${QUEUE}" 200
            elif test "${NOWAITCLUSTER}" -eq "0"
            then
                echo "waitcluster does not work on slurm cluster"
            fi

            # loop over settings
            for SETNAME in "${SETTINGSLIST[@]}"
            do
                # infer the names of all involved files from the arguments
                # defines the following environment variables: OUTFILE, ERRFILE, EVALFILE, OBJECTIVEVAL, SHORTPROBNAME,
                #                                              FILENAME, SKIPINSTANCE, BASENAME, TMPFILE, SETFILE
                . ./configuration_logfiles.sh "${INIT}" "${COUNT}" "${INSTANCE}" "${BINID}" "${PERMUTE}" "${SEEDS}" "${SETNAME}" "${TSTNAME}" "${CONTINUE}" "${QUEUE}" \
                                              "${p}" "${s}" "${THREADS}" "${GLBSEEDSHIFT}" "${STARTPERM}"

                # skip instance if log file is present and we want to continue a previously launched test run
                if test "${SKIPINSTANCE}" = "true"
                then
                    continue
                fi
                # check if binary exists. The second condition checks whether there is a binary of that name directly available
                # independent of whether it is a symlink, file in the working directory, or application in the path
                if test -e "${SCIPPATH}/../${BINNAME}"
                then
                    EXECNAME="${DEBUGTOOLCMD}${SCIPPATH}/../${BINNAME}"
                elif type "${BINNAME}" >/dev/null 2>&1
                then
                    EXECNAME="${DEBUGTOOLCMD}${BINNAME}"
                fi

                # use specified python version if the binary ends with ".py"
                EXT="${BINNAME##*.}"
                if test "${EXT}" = "py"
                then
                    EXECNAME="${PYTHON} ${BINNAME}"
                fi

                # find out the solver that should be used
                SOLVER=$(stripversion "${BINNAME}")

                # xpress executable is called optimizer
                if test "${SOLVER}" = "optimizer"
                then
                    SOLVER="xpress"
                fi

                CONFFILE="configuration_tmpfile_setup_${SOLVER}.sh"

                # call tmp file configuration for the solver
                # this may modify the EXECNAME environment variable
                . ./"${CONFFILE}" "${INSTANCE}" "${SCIPPATH}" "${TMPFILE}" "${SETNAME}" "${SETFILE}" "${THREADS}" "${SETCUTOFF}" \
                            "${FEASTOL}" "${TIMELIMIT}" "${MEMLIMIT}" "${NODELIMIT}" "${LPS}" "${DISPFREQ}" "${REOPT}" "${OPTCOMMAND}" \
                            "${CLIENTTMPDIR}" "${FILENAME}" "${VISUALIZE}" "${SOLUFILE}" "${EMPHBENCHMARK}" "${CLOCKTYPE}" \
                            "${WITHCERTIFICATE}" "${KEEPSOL}"

                JOBNAME="$(capitalize ${SOLVER})${SHORTPROBNAME}"
                # additional environment variables needed by run.sh
                export TARGETFREQ
                export SOLVERPATH="${SCIPPATH}"
                # this looks wrong but is totally correct
                export BASENAME="${FILENAME}"
                export FILENAME="${INSTANCE}"
                export CLIENTTMPDIR
                export OUTPUTDIR
                export HARDTIMELIMIT
                export HARDMEMLIMIT
                export CHECKERPATH="${SCIPPATH}/solchecker"
                export SETFILE
                export TIMELIMIT
                export EXECNAME
                export VIPRCHECKNAME=viprchk
                export VIPRCOMPNAME=viprcomp
                export VIPRCOMPRESSNAME=viprttn

                if test "${SLURMACCOUNT}" == "default"
                then
                    SLURMACCOUNT="${ACCOUNT}"
                fi

                WRITESETTINGS="false"
                if test "${INIT}" = "true"
                then
                    if test "${SOLVER}" = "scip"
                    then
                        WRITESETTINGS="true"
                        echo -e "#!/usr/bin/env bash \n ${EXECNAME} -s ${SETTINGS} -c 'set save ${CHECKSETFILE} quit'" > write-settings.sh
                    fi
                fi

                # check queue type
                if test  "${QUEUETYPE}" = "srun"
                then
                    if test "${WRITESETTINGS}" = "true"
                    then
                        sbatch --job-name=write-settings -n 1 -c "${THREADS_SAFE}" --mem="${HARDMEMLIMIT}" -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=medium-medium:Performance ${EXCLUSIVE} --output=/dev/null write-settings.sh
                    fi

                    if test "${AUTO_PPN_PENDING}" -eq 1
                    then
                        # auto mode: accumulate instance arguments into current batch
                        AUTO_BATCH_ARGS+=("${SOLVERPATH}" "${BASENAME}" "${FILENAME}" "${CLIENTTMPDIR}" "${OUTPUTDIR}" \
                            "${HARDTIMELIMIT}" "${HARDMEMLIMIT}" "${CHECKERPATH}" "${SETFILE}" "${TIMELIMIT}" \
                            "${EXECNAME}" "${VIPRCHECKNAME}" "${VIPRCOMPNAME}" "${VIPRCOMPRESSNAME}")
                        AUTO_BATCH_COUNT=$(( AUTO_BATCH_COUNT + 1 ))
                        AUTO_BATCH_NAMES="${AUTO_BATCH_NAMES} ${SHORTPROBNAME}"
                        if test -z "${AUTO_BATCH_JOBNAME}"
                        then
                            AUTO_BATCH_JOBNAME="${JOBNAME}"
                        fi
                        if [[ -z "${AUTO_BATCH_TIME}" || "${HARDTIMELIMIT}" > "${AUTO_BATCH_TIME}" ]]
                        then
                            AUTO_BATCH_TIME="${HARDTIMELIMIT}"
                        fi
                        # flush when batch is full
                        if test "${AUTO_BATCH_COUNT}" -ge "${AUTO_PPN}"
                        then
                            flush_auto_batch
                        fi
                    else
                        export SRUN="srun --exact -n 1 -c ${THREADS_SAFE} --mem=${HARDMEMLIMIT} --propagate=STACK --cpu_bind=verbose,cores ${SRUN_FLAGS}"
                        if test "${CLUSTERNODES}" = "all" && test "${EXCLUDENODES}" = "none"
                        then
                            echo sbatch --job-name="${JOBNAME}" --constraint="${CONSTRAINT}" -n 1 -c "${THREADS_SAFE}" --mem="${HARDMEMLIMIT}" -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=medium-medium:Performance ${EXCLUSIVE} --output=/dev/null run.sh
                            sbatch --job-name="${JOBNAME}" --constraint="${CONSTRAINT}" -n 1 -c "${THREADS_SAFE}" --mem="${HARDMEMLIMIT}" -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=medium-medium:Performance ${EXCLUSIVE} --output=/dev/null run.sh
                        elif test "${CLUSTERNODES}" != "all" && test "${EXCLUDENODES}" = "none"
                        then
                            echo sbatch --job-name="${JOBNAME}" --constraint="${CONSTRAINT}" -n 1 -c "${THREADS_SAFE}" --mem="${HARDMEMLIMIT}" -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=medium-medium:Performance ${EXCLUSIVE} -w "${CLUSTERNODES}" --output=/dev/null run.sh
                            sbatch --job-name="${JOBNAME}" --constraint="${CONSTRAINT}" -n 1 -c "${THREADS_SAFE}" --mem="${HARDMEMLIMIT}" -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=medium-medium:Performance ${EXCLUSIVE} -w "${CLUSTERNODES}" --output=/dev/null run.sh
                        elif test "${CLUSTERNODES}" = "all" && test "${EXCLUDENODES}" != "none"
                        then
                            echo sbatch --job-name="${JOBNAME}" --constraint="${CONSTRAINT}" -n 1 -c "${THREADS_SAFE}" --mem="${HARDMEMLIMIT}" -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=medium-medium:Performance ${EXCLUSIVE} -x "${EXCLUDENODES}" --output=/dev/null run.sh
                            sbatch --job-name="${JOBNAME}" --constraint="${CONSTRAINT}" -n 1 -c "${THREADS_SAFE}" --mem="${HARDMEMLIMIT}" -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=medium-medium:Performance ${EXCLUSIVE} -x "${EXCLUDENODES}" --output=/dev/null run.sh
                        else
                            echo sbatch --job-name="${JOBNAME}" --constraint="${CONSTRAINT}" -n 1 -c "${THREADS_SAFE}" --mem="${HARDMEMLIMIT}" -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=medium-medium:Performance ${EXCLUSIVE} -w "${CLUSTERNODES}" -x "${EXCLUDENODES}" --output=/dev/null run.sh
                            sbatch --job-name="${JOBNAME}" --constraint="${CONSTRAINT}" -n 1 -c "${THREADS_SAFE}" --mem="${HARDMEMLIMIT}" -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=medium-medium:Performance ${EXCLUSIVE} -w "${CLUSTERNODES}" -x "${EXCLUDENODES}" --output=/dev/null run.sh
                        fi
                    fi
                else
                    if test "${WRITESETTINGS}" = "true"
                    then
                        qsub -l walltime="${HARDTIMELIMIT}" -l nodes=1:ppn=$PPN -N write-settings \
                            -V -q "${CLUSTERQUEUE}" -o /dev/null -e /dev/null write-settings.sh
                        rm write-settings.sh
                    fi

                    # -V to copy all environment variables
                    qsub -l walltime="${HARDTIMELIMIT}" -l nodes=1:ppn="${PPN}" -N "${JOBNAME}" \
                        -V -q "${CLUSTERQUEUE}" -o /dev/null -e /dev/null run.sh

                fi

                if test "$WRITESETTINGS" = "true"
                then
                    rm write-settings.sh
                fi

            done # end for SETNAME
        done # end for PERMUTE
    done #end for SEEDS

    # after the first termination of the set loop, no file needs to be initialized anymore
    INIT="false"
done # end for TSTNAME

# flush any remaining instances in auto mode
if test "${AUTO_PPN_PENDING}" -eq 1 && test "${AUTO_BATCH_COUNT}" -gt 0
then
    flush_auto_batch
fi
