#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

if test -z "$GMSDIR"
then
  echo "Error: finishgamscluster.sh called with empty GMSDIR variable."
  exit 0
fi

if test -d "$GMSDIR"
then
  rm $GMSDIR/*
  rmdir $GMSDIR
fi

if test -z "$EVALFILE"
then
  echo "Error: finishgamscluster.sh called with empty EVALFILE variable."
  exit 0
fi

./evalcheck_gamscluster.sh $EVALFILE
