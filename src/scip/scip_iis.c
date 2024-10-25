/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   scip_iis.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for IIS plugins
 * @author Mark Turner
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/debug.h"
#include "scip/iis.h"
#include "scip/pub_message.h"
#include "scip/scip_iis.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"

/** creates an IIS and includes it in SCIP
 *
 *  @note this method has all IIS callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludeIISBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeIIS(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of IIS */
   const char*           desc,               /**< description of IIS */
   int                   priority,           /**< priority of the IIS */
   SCIP_DECL_IISCOPY     ((*iiscopy)),        /**< copy method of IIS or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_IISFREE     ((*iisfree)),       /**< destructor of IIS */
   SCIP_DECL_IISGENERATE ((*iisgenerate)),   /**< IIS generation method */
   SCIP_IISDATA*         iisdata             /**< IIS data */
   )
{
   SCIP_IIS* iis;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeIIS", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether the IIS is already present */
   if( SCIPfindIIS(scip, name) != NULL )
   {
      SCIPerrorMessage("IIS <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPiisCreate(&iis, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         iiscopy, iisfree, iisgenerate, iisdata) );
   SCIP_CALL( SCIPsetIncludeIIS(scip->set, iis) );

   return SCIP_OKAY;
}

/** Creates an IIS and includes it in SCIP with its most fundamental callbacks.
 *
 *  All non-fundamental (or optional) callbacks as, e.g., copy and free callbacks, will be set to NULL. Optional
 *  callbacks can be set via specific setter functions, see SCIPsetIISCopy(), and SCIPsetIISFree()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeIIS() instead
 */
SCIP_RETCODE SCIPincludeIISBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IIS**            iis,                /**< reference to an IIS, or NULL */
   const char*           name,               /**< name of IIS */
   const char*           desc,               /**< description of IIS */
   int                   priority,           /**< priority of the IIS in standard mode */
   SCIP_DECL_IISGENERATE((*iisgenerate)),    /**< IIS generation method */
   SCIP_IISDATA*         iisdata             /**< IIS data */
   )
{
    SCIP_IIS* iisptr;

    SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeIISBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

    /* check whether the IIS is already present */
    if( SCIPfindIIS(scip, name) != NULL )
    {
        SCIPerrorMessage("IIS <%s> already included.\n", name);
        return SCIP_INVALIDDATA;
    }

    SCIP_CALL( SCIPiisCreate(&iisptr, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
          NULL, NULL, iisgenerate, iisdata) );
    SCIP_CALL( SCIPsetIncludeIIS(scip->set, iisptr) );

    if( iis != NULL )
        *iis = iisptr;

    return SCIP_OKAY;
}

/** sets copy method of IIS */
SCIP_RETCODE SCIPsetIISCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_DECL_IISCOPY     ((*iiscopy))        /**< copy method of IIS or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
    SCIP_CALL( SCIPcheckStage(scip, "SCIPsetIISCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

    assert(iis != NULL);

    SCIPiisSetCopy(iis, iiscopy);

    return SCIP_OKAY;
}

/** sets destructor method of IIS */
SCIP_RETCODE SCIPsetIISFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_DECL_IISFREE     ((*iisfree))        /**< destructor of IIS */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetIISFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(iis != NULL);

   SCIPiisSetFree(iis, iisfree);

   return SCIP_OKAY;
}

/** the execution method that iterates over the IIS plugins */
SCIP_RETCODE SCIPgenerateIIS(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPiisGenerate(scip->set) );
   
   return SCIP_OKAY;
}

/** returns the IIS of the given name, or NULL if not existing */
SCIP_IIS* SCIPfindIIS(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of the IIS */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindIIS(scip->set, name);
}

/** returns the array of currently available IISs */
SCIP_IIS** SCIPgetIIS(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortIIS(scip->set);

   return scip->set->iis;
}

/** returns the number of currently available IISs */
int SCIPgetNIIS(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->niis;
}

/** sets the priority of an IIS */
SCIP_RETCODE SCIPsetIISPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IIS*             iis,                /**< IIS */
   int                   priority            /**< new priority of the IIS */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPiisSetPriority(iis, scip->set, priority);

   return SCIP_OKAY;
}
