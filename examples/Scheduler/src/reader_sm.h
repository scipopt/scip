/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_sm.h
 * @brief  scheduling problem file reader for RCPSP format
 * @author Michael Bastubbe
 * @author Stefan Heinz
 *
 * This reader is capabale of parsing resource-constrained project scheduling problem (RCPSP) instances. The <a
 * href="http://129.187.106.231/psplib/datasm.html">PSPlib</a> provides several instances set.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_SM_H__
#define __SCIP_READER_SM_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the sm file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderSm(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a cumulative scheduling problem */
extern
SCIP_RETCODE SCIPcreateSchedulingProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           problemname,        /**< problem name */
   const char**          jobnames,           /**< job names, or NULL */
   const char**          resourcenames,      /**< resource names, or NULL */
   int**                 demands,            /**< demand matrix resource job demand */
   SCIP_DIGRAPH*         precedencegraph,    /**< direct graph to store the precedence conditions */
   int*                  durations,          /**< array to store the processing for each job */
   int*                  capacities,         /**< array to store the different capacities */
   int                   njobs,              /**< number of jobs to be parsed */
   int                   nresources          /**< number of capacities to be parsed */
   );

#ifdef __cplusplus
}
#endif

#endif
