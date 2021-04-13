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

# Cleans up after gams testrun, calls 'evalcheck_gamscluster.sh'.
# To be invoked by 'check_gamscluster.sh'.

# Input environment variables
# GMSDIR needs to be defined, the corresponding directory will be deleted.

# New environment variables defined by this script: None

if test -z "${GMSDIR}"
then
    echo "Error: finishgamscluster.sh called with empty GMSDIR variable."
    exit 0
fi

if test -d "${GMSDIR}"
then
    rm "${GMSDIR}/*"
    rmdir "${GMSDIR}"
fi

if test -z "${EVALFILE}"
then
    echo "Error: finishgamscluster.sh called with empty EVALFILE variable."
    exit 0
fi

./evalcheck_gamscluster.sh "${EVALFILE}"
