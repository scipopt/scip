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

/**@file   iis.c
 * @ingroup OTHER_CFILES
 * @brief  methods for IIS
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/set.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/scip.h"
#include "scip/iis.h"

#include "scip/struct_iis.h"


/** method to call, when the priority of an IIS was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdIISPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetIISPriority() to mark the IIS unsorted */
   SCIP_CALL( SCIPsetIISPriority(scip, (SCIP_IIS*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** internal method for creating a IIS */
static
SCIP_RETCODE doIISCreate(
   SCIP_IIS**            iis,                /**< pointer to store IIS */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of iis */
   const char*           desc,               /**< description of iis */
   int                   priority,           /**< priority of the iis */
   SCIP_DECL_IISCOPY     ((*iiscopy)),       /**< copy method of IIS or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_IISFREE     ((*iisfree)),       /**< destructor of IIS */
   SCIP_DECL_IISGENERATE ((*iisgenerate)),   /**< IIS generation method */
   SCIP_IISDATA*         iisdata             /**< IIS data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(iis != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(iisgenerate != NULL);

   SCIP_ALLOC( BMSallocMemory(iis) );
   BMSclearMemory(*iis);

   SCIP_ALLOC( BMSduplicateMemoryArray(&(*iis)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*iis)->desc, desc, strlen(desc)+1) );
   (*iis)->priority = priority;
   (*iis)->iiscopy = iiscopy;
   (*iis)->iisfree = iisfree;
   (*iis)->iisgenerate = iisgenerate;
   (*iis)->iisdata = iisdata;

   /* create clocks */
   SCIP_CALL( SCIPclockCreate(&(*iis)->iistime, SCIP_CLOCKTYPE_DEFAULT) );

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "iis/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of iis generation rule <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*iis)->priority, FALSE, priority, INT_MIN/4, INT_MAX/2,
         paramChgdIISPriority, (SCIP_PARAMDATA*)(*iis)) ); /*lint !e740*/

   return SCIP_OKAY;
}


/** creates an IIS */
SCIP_RETCODE SCIPiisCreate(
   SCIP_IIS**            iis,                /**< pointer to store IIS */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of IIS */
   const char*           desc,               /**< description of IIS */
   int                   priority,           /**< priority of the IIS in standard mode */
   SCIP_DECL_IISCOPY     ((*iiscopy)),       /**< copy method of IIS or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_IISFREE     ((*iisfree)),       /**< destructor of IIS */
   SCIP_DECL_IISGENERATE ((*iisgenerate)),   /**< IIS generation method */
   SCIP_IISDATA*         iisdata             /**< IIS data */
   )
{
   assert(iis != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(iisgenerate != NULL);

   SCIP_CALL_FINALLY( doIISCreate(iis, set, messagehdlr, blkmem, name, desc, priority,
         iiscopy, iisfree, iisgenerate, iisdata), (void) SCIPiisFree(iis, set) );

   return SCIP_OKAY;
}

/** gets name of IIS */
const char* SCIPiisGetName(
   SCIP_IIS*             iis                 /**< IIS */
   )
{
   assert(iis != NULL);

   return iis->name;
}

/** calls IIS generation method */
SCIP_RETCODE SCIPiisGenerate(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;
   SCIP_RESULT result = SCIP_DIDNOTFIND;

   /* sort the iis's by priority */
   SCIPsetSortIIS(set);

   /* try all IIS generators until one succeeds */
   for( i = 0; i < set->niis && result == SCIP_DIDNOTFIND; ++i )
   {
      SCIP_IIS* iis;

      iis = set->iis[i];

      assert(iis != NULL);

      /* start timing */
      SCIPclockStart(iis->iistime, set);

      SCIP_CALL( iis->iisgenerate(set->scip, iis, &result) );

      /* stop timing */
      SCIPclockStop(iis->iistime, set);

      assert(result == SCIP_SUCCESS || result == SCIP_DIDNOTFIND);
   }

   return SCIP_OKAY;
}

/** gets description of IIS */
const char* SCIPiisGetDesc(
   SCIP_IIS*             iis                 /**< IIS */
   )
{
   assert(iis != NULL);

   return iis->desc;
}

/** copies the given IIS to a new scip */
SCIP_RETCODE SCIPiisCopyInclude(
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(iis != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( iis->iiscopy != NULL )
   {
      SCIPsetDebugMsg(set, "including IIS %s in subscip %p\n", SCIPiisGetName(iis), (void*)set->scip);
      SCIP_CALL( iis->iiscopy(set->scip, iis) );
   }
   return SCIP_OKAY;
}

/** frees memory of IIS */
SCIP_RETCODE SCIPiisFree(
   SCIP_IIS**            iis,                /**< pointer to IIS data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(iis != NULL);

   if( *iis == NULL )
      return SCIP_OKAY;

   assert(set != NULL);

   /* call destructor of IIS */
   if( (*iis)->iisfree != NULL )
   {
      SCIP_CALL( (*iis)->iisfree(set->scip, *iis) );
   }

   /* free clocks */
   SCIPclockFree(&(*iis)->iistime);

   BMSfreeMemoryArrayNull(&(*iis)->name);
   BMSfreeMemoryArrayNull(&(*iis)->desc);
   BMSfreeMemory(iis);

   return SCIP_OKAY;
}

/** gets user data of IIS */
SCIP_IISDATA* SCIPiisGetData(
   SCIP_IIS*             iis                 /**< IIS */
   )
{
   assert(iis != NULL);

   return iis->iisdata;
}

/** sets user data of IIS; user has to free old data in advance! */
void SCIPiisSetData(
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_IISDATA*         iisdata             /**< new IIS user data */
   )
{
   assert(iis != NULL);

   iis->iisdata = iisdata;
}

/** gets priority of IIS */
int SCIPiisGetPriority(
   SCIP_IIS*             iis                 /**< IIS */
   )
{
   assert(iis != NULL);

   return iis->priority;
}

/** enables or disables all clocks of @p iis, depending on the value of the flag */
void SCIPiisEnableOrDisableClocks(
   SCIP_IIS*             iis,                /**< the IIS for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the IIS be enabled? */
   )
{
   assert(iis != NULL);

   SCIPclockEnableOrDisable(iis->iistime, enable);
}


/* new callback/method setter methods */

/** sets copy method of IIS */
void SCIPiisSetCopy(
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_DECL_IISCOPY     ((*iiscopy))        /**< copy method of IIS or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   assert(iis != NULL);

   iis->iiscopy = iiscopy;
}

/** sets destructor method of IIS */
void SCIPiisSetFree(
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_DECL_IISFREE     ((*iisfree))        /**< destructor of IIS */
   )
{
   assert(iis != NULL);

   iis->iisfree = iisfree;
}

/** sets priority of IIS */
void SCIPiisSetPriority(
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the IIS */
   )
{
   assert(iis != NULL);
   assert(set != NULL);

   iis->priority = priority;
   set->iissorted = FALSE;
}

/** gets time in seconds used in this IIS */
SCIP_Real SCIPiisGetTime(
   SCIP_IIS*             iis                 /**< IIS */
   )
{
   assert(iis != NULL);

   return SCIPclockGetTime(iis->iistime);
}

/** get number of times the IIS was called */
SCIP_Longint SCIPiisGetNCalls(
   SCIP_IIS*             iis                 /**< IIS */
   )
{
  assert(iis != NULL);

  return iis->ncalls;
}

/** compares two IIS's w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPiisComp)
{  /*lint --e{715}*/
   return ((SCIP_IIS*)elem2)->priority - ((SCIP_IIS*)elem1)->priority;
}
