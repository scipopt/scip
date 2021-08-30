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

# Evaluates one or more testrun by each checking the concatenated outfile.
# Is to be invoked inside a 'check*.sh' script.
#
# Usage: from folder 'check' call
# ./evalcheck.sh results/check.*.eval

AWKARGS=""
FILES=""

for i in $@
do
    if test ! -e "${i}"
    then
        AWKARGS="${AWKARGS} ${i}"
    else
        FILES="${FILES} ${i}"
    fi
done

for FILE in ${FILES}
do
    # run check.awk (or the solver specialization) to evaluate the outfile
    . ./evaluate.sh "${FILE}"
done
