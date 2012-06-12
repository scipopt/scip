#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# the list was selected from
# grep -h 'ifdef' src/scip/* | awk ' // { printf("-D%s \\\n", $2); }' | sort -u
SCIPDEFS="-DALLDIFFERENT -DBETTERWEIGHTFORDEMANDNODES -DBRANCHLP -DBRANCHTOLINEARITY -DCHECKCONSARRAYS -DCHECKIMPLINBILINEAR -DCOUNTNETWORKVARIABLETYPES -DDEBUG_LEXDUAL -DDECONVEXIFY -DFIXSIMPLEVALUE -DGMLGATEPRINTING -DIMPLINTSARECONT -DLPI_CLP_DEBUG_WRITE_FILES -DMAKECUTINTEGRAL -DMAKEINTCUTINTEGRAL -DMCF_DEBUG -DOUTPUTGRAPH -DPRINT_MATRIX -DSCIP_CONFGRAPH -DSCIP_DEBUG -DSCIP_DEBUG_LP_INTERFACE -DSCIP_ENABLE_IISCHECK -DSCIP_MORE_DEBUG -DSCIP_OUPUT -DSCIP_OUTPUT -DSCIP_STATISTIC -DSCIP_WRITEPROB -DSEPARATEROWS -DSHOW_SCI -DSTRONGBRANCH_RESTOREBASIS -DTIEBREAKING -DTIMEEVENTEXEC -DTYPESTATISTICS -DUSECAPACITYFORTIEBREAKING -DUSEFLOWFORTIEBREAKING -DVARUSES -DWITH_BOUNDFLIPPING -DWITHDECOMPOSE -DWITH_LPSCHECK -DWITH_PRINTORIGCONSTYPES -DWITH_STATISTICS -DWITHSUBROOTS -DXPRS_SOLUTIONFILE -DZEROHALF__PRINT_STATISTICS"

SCIPDIR=..

cd $SCIPDIR
make OPT=dbg USRCFLAGS="$SCIPDEFS" clean
make OPT=dbg USRCFLAGS="$SCIPDEFS"
