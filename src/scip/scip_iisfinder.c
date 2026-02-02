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

/**@file   scip_iisfinder.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for IIS plugins
 * @author Mark Turner
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/debug.h"
#include "scip/iisfinder.h"
#include "scip/pub_message.h"
#include "scip/scip_iisfinder.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"

/** creates an IIS finder and includes it in SCIP
 *
 *  @note this method has all IIS finder callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludeIISfinderBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeIISfinder(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of IIS finder */
   const char*           desc,               /**< description of IIS finder */
   int                   priority,           /**< priority of the IIS finder */
   SCIP_Bool             enable,             /**< whether the IIS finder should be enabled */
   SCIP_DECL_IISFINDERCOPY ((*iisfindercopy)), /**< copy method of IIS finder or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_IISFINDERFREE ((*iisfinderfree)), /**< destructor of IIS finder */
   SCIP_DECL_IISFINDEREXEC ((*iisfinderexec)), /**< IIS finder execution method */
   SCIP_IISFINDERDATA*   iisfinderdata       /**< IIS finder data */
   )
{
   SCIP_IISFINDER* iisfinder;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeIISfinder", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether the IIS is already present */
   if( SCIPfindIISfinder(scip, name) != NULL )
   {
      SCIPerrorMessage("IIS <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPiisfinderCreate(&iisfinder, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         enable, iisfindercopy, iisfinderfree, iisfinderexec, iisfinderdata) );
   SCIP_CALL( SCIPsetIncludeIISfinder(scip->set, iisfinder) );

   return SCIP_OKAY;
}

/** Creates an IIS finder and includes it in SCIP with its most fundamental callbacks.
 *
 *  All non-fundamental (or optional) callbacks as, e.g., copy and free callbacks, will be set to NULL. Optional
 *  callbacks can be set via specific setter functions, see SCIPsetIISfinderCopy(), and SCIPsetIISfinderFree()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeIISfinder() instead
 */
SCIP_RETCODE SCIPincludeIISfinderBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IISFINDER**      iisfinder,          /**< reference to an IIS finder, or NULL */
   const char*           name,               /**< name of IIS finder */
   const char*           desc,               /**< description of IIS finder */
   int                   priority,           /**< priority of the IIS finder in standard mode */
   SCIP_Bool             enable,             /**< whether the IIS finder should be enabled */
   SCIP_DECL_IISFINDEREXEC((*iisfinderexec)), /**< IIS finder execution method */
   SCIP_IISFINDERDATA*   iisfinderdata       /**< IIS finder data */
   )
{
    SCIP_IISFINDER* iisfinderptr;

    SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeIISfinderBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

    /* check whether the IIS is already present */
    if( SCIPfindIISfinder(scip, name) != NULL )
    {
        SCIPerrorMessage("IIS finder <%s> already included.\n", name);
        return SCIP_INVALIDDATA;
    }

    SCIP_CALL( SCIPiisfinderCreate(&iisfinderptr, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
          enable, NULL, NULL, iisfinderexec, iisfinderdata) );
    SCIP_CALL( SCIPsetIncludeIISfinder(scip->set, iisfinderptr) );

    if( iisfinder != NULL )
        *iisfinder = iisfinderptr;

    return SCIP_OKAY;
}

/** sets copy method of IIS finder */
SCIP_RETCODE SCIPsetIISfinderCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_DECL_IISFINDERCOPY ((*iisfindercopy)) /**< copy method of IIS finder or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
    SCIP_CALL( SCIPcheckStage(scip, "SCIPsetIISfinderCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

    assert(iisfinder != NULL);

    SCIPiisfinderSetCopy(iisfinder, iisfindercopy);

    return SCIP_OKAY;
}

/** sets destructor method of IIS finder */
SCIP_RETCODE SCIPsetIISfinderFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_DECL_IISFINDERFREE ((*iisfinderfree)) /**< destructor of IIS finder */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetIISfinderFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(iisfinder != NULL);

   SCIPiisfinderSetFree(iisfinder, iisfinderfree);

   return SCIP_OKAY;
}

/** the execution method that iterates over the IIS finder plugins */
SCIP_RETCODE SCIPgenerateIIS(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPiisGenerate(scip->set) );

   return SCIP_OKAY;
}

/** returns the IIS finder of the given name, or NULL if not existing */
SCIP_IISFINDER* SCIPfindIISfinder(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of the IIS finder */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindIISfinder(scip->set, name);
}

/** returns the array of currently available IIS finders */
SCIP_IISFINDER** SCIPgetIISfinders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortIISfinders(scip->set);

   return scip->set->iisfinders;
}

/** returns the number of currently available IIS finders */
int SCIPgetNIISfinders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->niisfinders;
}

/** sets the priority of an IIS finder */
SCIP_RETCODE SCIPsetIISfinderPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   int                   priority            /**< new priority of the IIS finder */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert( iisfinder != NULL );

   SCIPiisfinderSetPriority(iisfinder, scip->set, priority);

   return SCIP_OKAY;
}

/** Gets the IIS storage.
 *
 *  @return the \ref SCIP_IIS iis storage.
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
SCIP_IIS* SCIPgetIIS(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetIIS", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->iis;
}
