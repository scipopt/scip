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

/**@file   reader_cyc.h
 * @brief  file reader for cycle clustering instances
 * @author Leon Eifler
 *
 * This file implements the reader for the cycle clustering problem. The data is read from a matrix, entries separated
 * by whitespace. The first line in the file has to be of the form "# p nstates ncluster",
 * where nstates is the size of the matrix and ncluster is the number of clusters that should be used.
 * The file has to have the ending ".cyc" to be recognized by the reader.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_CYC_H__
#define __SCIP_READER_CYC_H__

#include "scip/scip.h"
#include "tclique/tclique.h"
#include "scip/cons_setppc.h"
#include "scip/type_cons.h"
#include "scip/scip.h"
#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the cyc file reader into SCIP */
SCIP_RETCODE SCIPincludeReaderCyc(
   SCIP*             scip                 /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif



#endif
