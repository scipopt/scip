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

# Evaluates a testrun and concatenates the individual logfiles, possibly uploads to rubberband.
#
# If the environment variable RBCLI_TAG is set, then it will be passed as tags to rbcli.
#
# Usage: from folder 'check' call
# ./evalcheck_cluster.sh results/check.*.eval


export LANG=C
export LC_NUMERIC=C

REMOVE=0
UPLOAD=0
APPEND=0
EXPIRE=0
AWKARGS=""
FILES=""

for i in $@
do
    if test ! -e "${i}"
    then
        if test "${i}" = "-r"
        then
            REMOVE=1
        elif test "${i}" = "-U"
        then
            UPLOAD=1
        elif test "${i}" = "-E"
        then
            UPLOAD=1
            EXPIRE=1
        elif test "${i}" = "-R"
        then
            REMOVE=1
            UPLOAD=1
        elif test "${i}" = "-T"
        then
            REMOVE=1
            UPLOAD=1
            EXPIRE=1
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
    SETFILE="${DIR}/${EVALFILE}.set"
    METAFILE="${DIR}/${EVALFILE}.meta"

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

        echo > "${OUTFILE}"
        echo > "${ERRFILE}"
        echo ""
        echo "create overall output and error file for ${EVALFILE}"

        # check first if all out and err files exist for this eval-file.
        NMISSING=0
        for i in $(cat "${DIR}/${EVALFILE}.eval") DONE
        do
            if test "${i}" = "DONE"
            then
                break
            fi

            for extension in out err
            do
                FILE="${i}.${extension}"
                if ! test -e "${FILE}"
                then
                    echo "Missing ${FILE}"
                    ((NMISSING++))
                fi
            done
        done

        if [ "${NMISSING}" -gt 0 -a "${REMOVE}" -eq 1 ]
        then
            echo "Exiting because ${NMISSING} out/err file$([ ${NMISSING} -gt 1 ] && echo "s are" || echo " is" ) missing, please rerun without the REMOVE flag"
            exit
        fi

        for i in $(cat "${DIR}/${EVALFILE}.eval") DONE
        do
            if test "${i}" = "DONE"
            then
                break
            fi

            FILE="${i}.out"
            if test -e "${FILE}"
            then
                cat "${FILE}" >> "${OUTFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            else
                echo "@01 ${FILE} ==MISSING=="  >> "${OUTFILE}"
                echo                            >> "${OUTFILE}"
            fi

            FILE="${i}.err"
            if test -e "${FILE}"
            then
                cat "${FILE}" >> "${ERRFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            else
                echo "@01 ${FILE} ==MISSING=="  >> "${ERRFILE}"
                echo                            >> "${ERRFILE}"
            fi

            FILE="${i}.set"
            if test -e "${FILE}"
            then
                if ! test -e $SETFILE
                then
                    cp "${FILE}" "${SETFILE}"
                fi
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            fi

            FILE="${i}".tmp
            if test -e "${FILE}"
            then
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
    if test -e "${OUTFILE}"
    then
        echo "create results for ${EVALFILE}"

        # run check.awk (or the solver specialization) to evaluate the outfile
        . ./evaluate.sh "${OUTFILE}"

        # upload results to rubberband.zib.de
        if test "${UPLOAD}" = "1"
        then
            if test "${EXPIRE}" = "1"
            then
                RB_EXP_DATE=$(date '+%Y-%b-%d' -d "+6 weeks")
                echo "rbcli -e ${RB_EXP_DATE} up ${OUTFILE} ${ERRFILE} ${SETFILE} ${METAFILE}"
                rbcli -e "${RB_EXP_DATE}" up "${OUTFILE}" "${ERRFILE}" "${SETFILE}" "${METAFILE}"
            else
                if test -z "${RBCLI_TAG}"
                then
                    echo "rbcli up ${OUTFILE} ${ERRFILE} ${SETFILE} ${METAFILE}"
                    rbcli up "${OUTFILE}" "${ERRFILE}" "${SETFILE}" "${METAFILE}"
                else
                    echo "rbcli --tags ${RBCLI_TAG} up ${OUTFILE} ${ERRFILE} ${SETFILE} ${METAFILE}"
                    rbcli --tags "${RBCLI_TAG}" up "${OUTFILE}" "${ERRFILE}" "${SETFILE}" "${METAFILE}"
                fi
            fi
        fi
    fi
done
