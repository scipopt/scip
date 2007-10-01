/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: probdata_LOP.h,v 1.1 2007/10/01 13:41:53 bzfpfets Exp $"

#ifndef __HCP_PROBDATA_LOP__
#define __HCP_PROBDATA_LOP__

#include <scip/scip.h>


extern
SCIP_RETCODE LOPcreateProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of file to read */
   );

extern
SCIP_RETCODE LOPgenerateModel(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
SCIP_RETCODE LOPevalSolution(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
int LOPgetNElements(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
