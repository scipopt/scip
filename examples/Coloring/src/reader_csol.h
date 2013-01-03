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

/**@file   reader_csol.h
 * @brief  file reader and writer for vertex coloring solutions
 * @author Gerald Gamrath
 *
 * This file implements the reader and writer for coloring solution files.
 *
 * These files have the following structure:@n The first line contains the name of the problem, the
 * number of colors used in the solution, and - optional - the name of the algorithm that computed
 * this solution.  The second line lists the colors of the nodes, separated by spaces. It is sorted
 * increasingly by the node indices. The numbers for the colors start with 0.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_CSOL_H__
#define __SCIP_READER_CSOL_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the csol file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderCsol(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
