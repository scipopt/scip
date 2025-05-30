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

export LANG=C

REMOVE=0
AWKARGS=""
FILES=""

for i in $@
do
    if test ! -e "${i}"
    then
        if test "${i}" = "-r"
        then
            REMOVE=1
        else
            AWKARGS="${AWKARGS} ${i}"
        fi
    else
        FILES="${FILES} ${i}"
    fi
done

for FILE in ${FILES}
do

    DIR=$(dirname "${FILE}")
    EVALFILE=$(basename "${FILE}" .eval)
    EVALFILE=$(basename "${EVALFILE}" .out)

    OUTFILE="${DIR}/${EVALFILE}.out"
    ERRFILE="${DIR}/${EVALFILE}.err"
    TRCFILE="${DIR}/${EVALFILE}.trc"
    LSTFILE="${DIR}/${EVALFILE}.lst"
    RESFILE="${DIR}/${EVALFILE}.res"
    TEXFILE="${DIR}/${EVALFILE}.tex"
    PAVFILE="${DIR}/${EVALFILE}.pav"
    EXMFILE="${DIR}/${EVALFILE}.exm"

    # check if the eval file exists; if this is the case construct the overall solution files
    if test -e "${DIR}/${EVALFILE}.eval"
    then
        # in case an output file exists, copy it away to save the results
        DATEINT=$(date +"%s")
        if test -e "${OUTFILE}"
        then
            cp "${OUTFILE}" "${OUTFILE}.old-${DATEINT}"
        fi
        if test -e "${ERRFILE}"
        then
            cp "${ERRFILE}" "${ERRFILE}.old-${DATEINT}"
        fi
        if test -e "${TRCFILE}"
        then
            cp "${TRCFILE}" "${TRCFILE}.old-${DATEINT}"
        fi
        if test -e "${EXMFILE}"
        then
            cp "${EXMFILE}" "${EXMFILE}.old-${DATEINT}"
            rm "${EXMFILE}"
        fi
        if test -e "${LSTFILE}"
        then
            cp "${LSTFILE}" "${LSTFILE}.old-${DATEINT}"
        fi

        echo > "${OUTFILE}"
        echo > "${ERRFILE}"
        # initialize gams trace file
        echo "* Trace Record Definition" > "${TRCFILE}"
        echo "* GamsSolve" >> "${TRCFILE}"
        echo "* InputFileName,ModelType,SolverName,OptionFile,Direction,NumberOfEquations,NumberOfVariables,NumberOfDiscreteVariables,NumberOfNonZeros,NumberOfNonlinearNonZeros," >> "${TRCFILE}"
        echo "* ModelStatus,SolverStatus,ObjectiveValue,ObjectiveValueEstimate,SolverTime,ETSolver,NumberOfIterations,NumberOfNodes" >> "${TRCFILE}"
        echo "*" >> "${TRCFILE}"
        echo > "${LSTFILE}"

        echo "create overall output, error, and trace file for ${EVALFILE}"

        for i in $(cat "${DIR}/${EVALFILE}.eval") DONE
        do
            if test "${i}" = "DONE"
            then
                break
            fi

            case "${i}" in SOLVER,* )
                SOLVER=${i:7}
                ;;
            esac

            # pass auxiliary lines about solver and limits form eval file to trace file
            case "${i}" in *,* )
                echo "* ${i}" >> "${TRCFILE}"
                continue
                ;;
            esac

            FILE="${i}.out"
            if test -e "${FILE}"
            then
                cat "${FILE}" >> "${OUTFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            else
                echo "Missing ${i}"
            fi

            FILE="${i}.err"
            if test -e "${FILE}"
            then
                cat "${FILE}" >> "${ERRFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            fi

            FILE="${i}.trc"
            if test -e "${FILE}"
            then
                grep -v "^*" "${FILE}" | sed -e "s/EXAMINER2/${SOLVER}/" >> "${TRCFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            else
                echo "Missing ${i}"
            fi

            FILE="${i}.exm"
            if test -e "${FILE}"
            then
                test -e "${EXMFILE}" || grep "^*" "${FILE}" > "${EXMFILE}"
                grep -v "^*" "${FILE}" >> "${EXMFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            fi

            FILE="${i}.lst"
            if test -e "${FILE}"
            then
                cat "${FILE}" >> "${LSTFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            fi
        done

        if test "${REMOVE}" = "1"
        then
            rm -f "${DIR}/${EVALFILE}.eval"
        fi
    fi

    # check if the out file exists
    if test -e "${DIR}/${EVALFILE}.out"
    then
        echo "create results for ${EVALFILE}"

        # detect test set
        TSTNAME=$(echo "${EVALFILE}" | sed 's/check.\([a-zA-Z0-9_-]*\).*/\1/g')
        echo "Testset ${TSTNAME}"

        if test -f "testset/${TSTNAME}.test"
        then
            TESTFILE="testset/${TSTNAME}.test"
        else
            TESTFILE=""
        fi

        # call method to obtain solution file
        # defines the following environment variable: SOLUFILE
        . ./configuration_solufile.sh "${TSTNAME}"

        # the variable AWKARGS needs to be without quotation marks here
        awk -f check_gams.awk -v "TEXFILE=${TEXFILE}" -v "PAVFILE=${PAVFILE}" ${AWKARGS} "${SOLUFILE}" "${TRCFILE}" | tee "${RESFILE}"
    fi
done
