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

/**@file   vardata_binpacking.h
 * @brief  Variable data containing the ids of constraints in which the variable appears
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file implements the handling of the variable data which is attached to each file. See SCIP_VarData and \ref BINPACKING_PRICER.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_VARDATA_BINPACKING__
#define __SCIP_VARDATA_BINPACKING__

#include "scip/scip.h"

/** create variable data */
SCIP_RETCODE SCIPvardataCreateBinpacking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   int*                  consids,            /**< array of constraints ids */
   int                   nconss              /**< number of constraints */
   );

/** get number of constraints */
int SCIPvardataGetNConsids(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

/** returns sorted constraint id array */
int* SCIPvardataGetConsids(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

/** creates variable */
SCIP_RETCODE SCIPcreateVarBinpacking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   );

/** prints vardata to file stream */
void SCIPvardataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata,            /**< variable data */
   FILE*                 file                /**< the text file to store the information into */
   );

#endif
