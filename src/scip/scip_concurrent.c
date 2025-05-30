/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   scip_concurrent.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for concurrent solving mode
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/concsolver.h"
#include "scip/debug.h"
#include "scip/pub_message.h"
#include "scip/scip_concurrent.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/syncstore.h"

/** creates a concurrent solver type and includes it in SCIP.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPincludeConcsolverType(
   SCIP*                               scip,                       /**< SCIP data structure */
   const char*                         name,                       /**< name of concurrent_solver */
   SCIP_Real                           prefpriodefault,            /**< the default preferred priority of this concurrent solver type */
   SCIP_DECL_CONCSOLVERCREATEINST      ((*concsolvercreateinst)),  /**< data copy method of concurrent solver */
   SCIP_DECL_CONCSOLVERDESTROYINST     ((*concsolverdestroyinst)), /**< data copy method of concurrent solver */
   SCIP_DECL_CONCSOLVERINITSEEDS       ((*concsolverinitseeds)),   /**< initialize random seeds of concurrent solver */
   SCIP_DECL_CONCSOLVEREXEC            ((*concsolverexec)),        /**< execution method of concurrent solver */
   SCIP_DECL_CONCSOLVERCOPYSOLVINGDATA ((*concsolvercopysolvdata)),/**< method to copy solving data */
   SCIP_DECL_CONCSOLVERSTOP            ((*concsolverstop)),        /**< terminate solving in concurrent solver */
   SCIP_DECL_CONCSOLVERSYNCWRITE       ((*concsolversyncwrite)),   /**< synchronization method of concurrent solver */
   SCIP_DECL_CONCSOLVERSYNCREAD        ((*concsolversyncread)),    /**< synchronization method of concurrent solver */
   SCIP_DECL_CONCSOLVERTYPEFREEDATA    ((*concsolvertypefreedata)),/**< method to free data of concurrent solver type */
   SCIP_CONCSOLVERTYPEDATA*            data                        /**< the concurent solver type's data */
   )
{
   SCIP_CONCSOLVERTYPE* concsolvertype;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeConcsolverType", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether concurrent solver type is already present */
   if( SCIPfindConcsolverType(scip, name) != NULL )
   {
      SCIPerrorMessage("concurrent solver type <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPconcsolverTypeCreate(&concsolvertype, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, prefpriodefault, concsolvercreateinst, concsolverdestroyinst,
         concsolverinitseeds, concsolverexec, concsolvercopysolvdata,
         concsolverstop, concsolversyncwrite, concsolversyncread,
         concsolvertypefreedata, data) );

   SCIP_CALL( SCIPsetIncludeConcsolverType(scip->set, concsolvertype) );

   return SCIP_OKAY;
}

/** returns the concurrent solver type with the given name, or NULL if not existing */
SCIP_CONCSOLVERTYPE* SCIPfindConcsolverType(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of concurrent_solver */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindConcsolverType(scip->set, name);
}

/** returns the array of included concurrent solver types */
SCIP_CONCSOLVERTYPE** SCIPgetConcsolverTypes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->concsolvertypes;
}

/** returns the number of included concurrent solver types */
int SCIPgetNConcsolverTypes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nconcsolvertypes;
}

/** Constructs the parallel interface to execute processes concurrently.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPconstructSyncstore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPconstructSyncstore", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPsyncstoreCreate(&scip->syncstore) );

   return SCIP_OKAY;
}

/** releases the current parallel interface
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *       - \ref SCIP_STAGE_FREE
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPfreeSyncstore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPfreeSyncstore", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsyncstoreRelease(&(scip->syncstore)) );

   return SCIP_OKAY;
}

/** Gets the parallel interface to execute processes concurrently.
 *
 *  @return the \ref SCIP_SYNCSTORE parallel interface pointer to submit jobs for concurrent processing.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_SYNCSTORE* SCIPgetSyncstore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSyncstore", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->syncstore;
}
