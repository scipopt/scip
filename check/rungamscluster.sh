#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# The script executing gams on one instance and producing the logfiles.
# Can be executed either locally or on a cluster node.
# Is to be invoked inside a 'check(_cluster)*.sh' script.

# check if tmp-path exists
if test ! -d "${CLIENTTMPDIR}"
then
    echo "Skipping test since the path for the tmp-dir does not exist."
    exit
fi

OUTFILE="${CLIENTTMPDIR}/${BASENAME}.out"
ERRFILE="${CLIENTTMPDIR}/${BASENAME}.err"
LSTFILE="${CLIENTTMPDIR}/${BASENAME}.lst"
TRCFILE="${CLIENTTMPDIR}/${BASENAME}.trc"
EXMFILE="${CLIENTTMPDIR}/${BASENAME}.exm"
WORKDIR="${CLIENTTMPDIR}/${BASENAME}.scr"
OPTDIR="${CLIENTTMPDIR}/${BASENAME}.opt"

# setup scratch directory
mkdir -p "${WORKDIR}"
GAMSOPTS="${GAMSOPTS} curdir=${WORKDIR}"

# ensure scratch directory is deleted and results are copied when exiting (normally or due to abort/interrupt)
trap "
rm -rf "${WORKDIR}" "${OPTDIR}";
test -e "${OUTFILE}" && mv "${OUTFILE}" "results/${BASENAME}.out"
test -e "${LSTFILE}" && mv "${LSTFILE}" "results/${BASENAME}.lst"
test -e "${ERRFILE}" && mv "${ERRFILE}" "results/${BASENAME}.err"
test -e "${TRCFILE}" && mv "${TRCFILE}" "results/${BASENAME}.trc"
test -e "${EXMFILE}" && mv "${EXMFILE}" "results/${BASENAME}.exm"
" EXIT


# create directory "${OPTDIR}" and put optionfile <solvername>.opt there
if test -n "${SETTINGS}"
then
    if test -d "${OPTDIR}"
    then
        rm -f "${OPTDIR}"/*
    else
        mkdir -p "${OPTDIR}"
    fi
    # replace all "${var}" by their value w.r.t. the current environment
    awk '{while(match($0,"[$]{[^}]*}")) {var=substr($0,RSTART+2,RLENGTH -3);gsub("[$]{"var"}",ENVIRON[var])}}1' "${SETTINGS}" > ${OPTDIR}/${SOLVER,,}.opt
    GAMSOPTS="${GAMSOPTS} optdir="${OPTDIR}" optfile=1"
fi

# initialize trace file
echo "* Trace Record Definition" > "${TRCFILE}"
echo "* GamsSolve" >> "${TRCFILE}"
echo "* InputFileName,ModelType,SolverName,OptionFile,Direction,NumberOfEquations,NumberOfVariables,NumberOfDiscreteVariables,NumberOfNonZeros,NumberOfNonlinearNonZeros," >> "${TRCFILE}"
echo "* ModelStatus,SolverStatus,ObjectiveValue,ObjectiveValueEstimate,SolverTime,ETSolver,NumberOfIterations,NumberOfNodes" >> "${TRCFILE}"

# setup examiner option file
if test "${EXAMINER}" = 1
then
    mkdir -p "${OPTDIR}"
    echo "subsolver ${SOLVER,,}" > "${OPTDIR}/examiner2.opt"
    if test -n "${SETTINGS}"
    then
        echo "subsolveropt 1" >> "${OPTDIR}/examiner2.opt"
    else
        GAMSOPTS="${GAMSOPTS} optdir=${OPTDIR} optfile=1"
    fi
    echo "scaled yes"           >> "${OPTDIR}/examiner2.opt"
    echo "unscaled yes"         >> "${OPTDIR}/examiner2.opt"
    echo "examinesolupoint yes" >> "${OPTDIR}/examiner2.opt"
    echo "examinesolvpoint yes" >> "${OPTDIR}/examiner2.opt"
    echo "examinegamspoint no"  >> "${OPTDIR}/examiner2.opt"
    echo "examineinitpoint no"  >> "${OPTDIR}/examiner2.opt"
    echo "trace ${EXMFILE}"     >> "${OPTDIR}/examiner2.opt"
    echo "tracestyle 1"         >> "${OPTDIR}/examiner2.opt"

    #Examiner currently ignores Trace Record Definition, so no sense in initializing it
    #cp -f "${TRCFILE}" "${EXMFILE}"
fi

# setup gams file that sets cutoff
# this only works for models that include %gams.u1% and where the model name is m (e.g., MINLPLib instances)
if test -n "${CUTOFF}"
then
    echo "m.cutoff = ${CUTOFF};" > "${WORKDIR}/include.u1"
    GAMSOPTS="${GAMSOPTS} u1=${WORKDIR}/include.u1"
fi

# add commands to .u1 file to read start solutions from available gdx files, if any
if test "${PASSSTARTSOL}" = 1
then
    for sol in ${INPUTDIR}/${GMSFILE/%.gms*/}*.gdx
    do
        if test -e "${sol}"
        then
            # create .u1 file if not existing yet and add u1 command to GAMS options
            if test ! -e "${WORKDIR}/include.u1"
            then
                touch "${WORKDIR}/include.u1"
                GAMSOPTS="${GAMSOPTS} u1=${WORKDIR}/include.u1"
            fi
            echo "execute_loadpoint '${sol}';" >> "${WORKDIR}/include.u1"
        fi
    done
fi

if test "${EXAMINER}" = 1 ; then
    ACTSOLVER=EXAMINER2
else
    ACTSOLVER="${SOLVER}"
fi

uname -a                                 > "${OUTFILE}"
uname -a                                 > "${ERRFILE}"
echo "@01 ${FILENAME} ==========="      >> "${OUTFILE}"
echo "@01 ${FILENAME} ==========="      >> "${ERRFILE}"
echo "-----------------------------"    >> "${OUTFILE}"
date                                    >> "${OUTFILE}"
date                                    >> "${ERRFILE}"
echo "-----------------------------"    >> "${OUTFILE}"
date +"@03 %s"                          >> "${OUTFILE}"

# run GAMS and check return code
"${GAMSBIN}" "${GMSFILE}" "${GAMSOPTS}" output="${LSTFILE}" inputdir="${INPUTDIR}" "${MODTYPE}"="${ACTSOLVER}" "${GDXFILE}" traceopt=3 trace="${TRCFILE}" >> "${OUTFILE}" 2>> "${ERRFILE}"
gamsrc=$?
if test "${gamsrc}" != 0
then
    echo "GAMS returned with error code ${gamsrc}. There was some problem." >>"${ERRFILE}"
    #TODO write 13/13 into trace file
fi

date +"@04 %s"                        >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${ERRFILE}"
echo                                  >> "${OUTFILE}"
echo "=ready="                        >> "${OUTFILE}"
