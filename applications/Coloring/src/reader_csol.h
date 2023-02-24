/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
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
SCIP_RETCODE SCIPincludeReaderCsol(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
