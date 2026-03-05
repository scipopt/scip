/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   reader_bd.h
 * @ingroup FILEREADERS
 * @brief  Benders file reader - this is used to read multiple files that form the master and subproblems of a
 *         Benders' decomposition
 * @author Stephen J. Maher
 *
 * The Benders file reader is designed to accept multiple files that correspond to the master and subproblems of a
 * problem decomposed using Benders' decomposition. The instance files will be read in and then the default Benders'
 * decomposition algorithm will be executed.
 *
 * The input file format has the keywords `MASTER` and `SUBPROBLEMS`. The line immediately after the `MASTER` keyword is
 * the instance file for the master problem. The `SUBPROBLEMS` keyword must be followed (on the same line) by the
 * number of subproblems, i.e. `SUBPROBLEMS 10`. The instance files for the subproblems follow the `SUBPROBLEMS`
 * keyword.
 *
 * An example of a benders file looks as follows.
 *
 * @verbinclude multi-zone.bd
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_BD_H__
#define __SCIP_READER_BD_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the benders file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeReaderBenders(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
