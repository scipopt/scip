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

/**@file   iisfinder.c
 * @ingroup OTHER_CFILES
 * @brief  methods for IIS finders
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/set.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/iisfinder.h"
#include "scip/iisfinder_greedy.h"

#include "scip/struct_iisfinder.h"


/** method to call, when the priority of an IIS finder was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdIISfinderPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetIISPriority() to mark the IIS unsorted */
   SCIP_CALL( SCIPsetIISfinderPriority(scip, (SCIP_IISFINDER*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** internal method for creating the subscip that will hold the IIS */
static
SCIP_RETCODE createSubscipIIS(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_IIS*             iis,                /**< pointer to store IIS */
   SCIP_Real             timelim,            /**< timelimit */
   SCIP_Longint          nodelim             /**< nodelimit */
   )
{
   SCIP_VAR** vars;
   SCIP_Bool success;
   int nvars;
   int i;

   /* Create the subscip used for storing the IIS */
   if( iis->subscip != NULL )
   {
      SCIPdebugMsg(set->scip, "An IIS for this problem already exists. Removing it before starting search procedure again.\n");

      /* free sub-SCIP */
      SCIP_CALL( SCIPiisReset(&iis) );
   }

   assert( iis->subscip == NULL );
   assert( iis->varsmap == NULL );
   assert( iis->conssmap == NULL );

   /* create a new SCIP instance */
   SCIP_CALL( SCIPcreate(&(iis->subscip)) );

   /* Create hash maps */
   SCIP_CALL( SCIPhashmapCreate(&(iis->varsmap), SCIPblkmem(set->scip), SCIPgetNOrigVars(set->scip)) );
   SCIP_CALL( SCIPhashmapCreate(&(iis->conssmap), SCIPblkmem(set->scip), SCIPgetNOrigConss(set->scip)) );

   /* create problem in sub-SCIP */
   SCIP_CALL( SCIPcopyOrig(set->scip, iis->subscip, iis->varsmap, iis->conssmap, "iis", TRUE, FALSE, FALSE, &success) );

   if( ! success )
      return SCIP_ERROR;

   /* Remove the objective */
   vars = SCIPgetOrigVars(iis->subscip);
   nvars = SCIPgetNOrigVars(iis->subscip);
   for( i = 0; i < nvars; i++ )
      SCIP_CALL( SCIPchgVarObj(iis->subscip, vars[i], 0.0 ) );

   /* copy parameter settings */
   /** @todo: Do we really want to copy the parameter settings?
   SCIP_CALL( SCIPcopyParamSettings(set->scip, iis->subscip) );
#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(iis->subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(iis->subscip, "display/freq", 100000000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(iis->subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(iis->subscip, "timing/statistictiming", FALSE) );
#endif
   SCIP_CALL( SCIPsetIntParam(iis->subscip, "limits/bestsol", 1) );
   SCIP_CALL( SCIPsetRealParam(iis->subscip, "limits/time", timelim - SCIPclockGetTime(iis->iistime)) );
   SCIP_CALL( SCIPsetLongintParam(iis->subscip, "limits/nodes", nodelim) );

   return SCIP_OKAY;
}

/** checks the problem for trivial infeasibility reasons, e.g. contradicting bounds */
static
SCIP_RETCODE checkTrivialInfeas(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_Bool*            trivial             /**< pointer to store whether the problem is trivially infeasible */
   )
{
   SCIP_CONS** origconss;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   int nconss;
   int nvars;
   int i;

   assert( trivial != NULL );
   *trivial = FALSE;
   nvars = SCIPgetNOrigVars(scip);
   vars = SCIPgetOrigVars(scip);
   for( i = 0; i < nvars; i++ )
   {
      if( SCIPisGT(scip, SCIPvarGetLbOriginal(vars[i]), SCIPvarGetUbOriginal(vars[i])) )
      {
         *trivial = TRUE;
         break;
      }
   }

   if( *trivial )
   {
      nconss = SCIPgetNOrigConss(scip);
      origconss = SCIPgetOrigConss(scip);
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &conss, origconss, nconss) );
      for( i = 0; i < nconss; i++ )
      {
         SCIP_CALL( SCIPdelCons(scip, conss[i]) );
      }
      SCIPfreeBlockMemoryArray(scip, &conss, nconss);
   }
   return SCIP_OKAY;
}

/** internal method for creating an IIS finder */
static
SCIP_RETCODE doIISfinderCreate(
   SCIP_IISFINDER**      iisfinder,          /**< pointer to store IIS finder */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of IIS finder */
   const char*           desc,               /**< description of IIS finder */
   int                   priority,           /**< priority of the IIS finder */
   SCIP_DECL_IISFINDERCOPY ((*iisfindercopy)), /**< copy method of IIS finder or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_IISFINDERFREE ((*iisfinderfree)), /**< destructor of IIS finder */
   SCIP_DECL_IISFINDEREXEC ((*iisfinderexec)), /**< IIS finder execution method */
   SCIP_IISFINDERDATA*   iisfinderdata       /**< IIS finder data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(iisfinder != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(iisfinderexec != NULL);

   SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, iisfinder) );

   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*iisfinder)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*iisfinder)->desc, desc, strlen(desc)+1) );
   (*iisfinder)->priority = priority;
   (*iisfinder)->iisfindercopy = iisfindercopy;
   (*iisfinder)->iisfinderfree = iisfinderfree;
   (*iisfinder)->iisfinderexec = iisfinderexec;
   (*iisfinder)->iisfinderdata = iisfinderdata;

   /* create clocks */
   SCIP_CALL( SCIPclockCreate(&(*iisfinder)->iisfindertime, SCIP_CLOCKTYPE_DEFAULT) );

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "iis/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of iis generation rule <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*iisfinder)->priority, FALSE, priority, INT_MIN/4, INT_MAX/2,
         paramChgdIISfinderPriority, (SCIP_PARAMDATA*)(*iisfinder)) ); /*lint !e740*/

   return SCIP_OKAY;
}


/** creates an IIS finder */
SCIP_RETCODE SCIPiisfinderCreate(
   SCIP_IISFINDER**      iisfinder,          /**< pointer to store IIS finder */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of IIS finder */
   const char*           desc,               /**< description of IIS finder */
   int                   priority,           /**< priority of the IIS finder in standard mode */
   SCIP_DECL_IISFINDERCOPY ((*iisfindercopy)), /**< copy method of IIS finder or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_IISFINDERFREE ((*iisfinderfree)), /**< destructor of IIS finder */
   SCIP_DECL_IISFINDEREXEC ((*iisfinderexec)), /**< IIS finder execution method */
   SCIP_IISFINDERDATA*   iisfinderdata       /**< IIS finder data */
   )
{
   assert(iisfinder != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(iisfinderexec != NULL);

   SCIP_CALL_FINALLY( doIISfinderCreate(iisfinder, set, messagehdlr, blkmem, name, desc, priority,
         iisfindercopy, iisfinderfree, iisfinderexec, iisfinderdata), (void) SCIPiisfinderFree(iisfinder, set, blkmem) );

   return SCIP_OKAY;
}

/** gets name of IIS finder */
const char* SCIPiisfinderGetName(
   SCIP_IISFINDER*       iisfinder           /**< IIS finder */
   )
{
   assert(iisfinder != NULL);

   return iisfinder->name;
}

/** calls IIS finder generation method */
SCIP_RETCODE SCIPiisGenerate(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_CONS** origconss;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_IIS* iis;
   int i;
   int j;
   int nconss;
   int nvars;
   SCIP_RESULT result = SCIP_DIDNOTFIND;
   SCIP_Bool silent;
   SCIP_Bool removebounds;
   SCIP_Bool minimal;
   SCIP_Bool stopafterone;
   SCIP_Bool removeunusedvars;
   SCIP_Bool trivial;
   SCIP_Bool islinear;
   SCIP_Real timelim;
   SCIP_Longint nodelim;

   /* start timing */
   SCIP_CALL( SCIPgetRealParam(set->scip, "iis/time", &timelim) );
   SCIP_CALL( SCIPgetLongintParam(set->scip, "iis/nodes", &nodelim) );
   SCIP_CALL( SCIPgetBoolParam(set->scip, "iis/silent", &silent) );

   /* sort the iis finders by priority */
   SCIPsetSortIISfinders(set);

   /* Get the IIS data. */
   iis = SCIPgetIIS(set->scip);

   /* Create the subscip used for storing the IIS */
   SCIP_CALL( createSubscipIIS(set, iis, timelim, nodelim) );

   SCIPclockStart(iis->iistime, set);

   /* If the model is not yet shown to be infeasible then check for infeasibility */
   if( SCIPgetStage(set->scip) == SCIP_STAGE_PROBLEM && set->iisfinder_checkinitfeas )
   {
      SCIP_CALL( SCIPsolve(iis->subscip) );
      if( SCIPgetStage(iis->subscip) == SCIP_STAGE_SOLVED )
      {
         switch( SCIPgetStatus(iis->subscip) )
         {
            case SCIP_STATUS_TIMELIMIT:
            case SCIP_STATUS_NODELIMIT:
            case SCIP_STATUS_TOTALNODELIMIT:
            case SCIP_STATUS_STALLNODELIMIT:
            case SCIP_STATUS_MEMLIMIT:
            case SCIP_STATUS_GAPLIMIT:
            case SCIP_STATUS_SOLLIMIT:
            case SCIP_STATUS_RESTARTLIMIT:
            case SCIP_STATUS_DUALLIMIT:
               SCIPinfoMessage(iis->subscip, NULL, "Some limit reached. Failed to prove infeasibility of initial problem.\n");
               SCIPclockStop(iis->iistime, set);
               return SCIP_OKAY;

            case SCIP_STATUS_INFORUNBD:
               SCIPinfoMessage(iis->subscip, NULL, "Initial problem is infeasible or Unbounded. Not performing IIS generation.\n");
               SCIPclockStop(iis->iistime, set);
               return SCIP_OKAY;

            case SCIP_STATUS_USERINTERRUPT:
            case SCIP_STATUS_TERMINATE:
               SCIPdebugMsg(iis->subscip, "User interrupt. Stopping.\n");
               SCIPclockStop(iis->iistime, set);
               return SCIP_OKAY;

            case SCIP_STATUS_BESTSOLLIMIT:
            case SCIP_STATUS_OPTIMAL:
            case SCIP_STATUS_UNBOUNDED:
            case SCIP_STATUS_PRIMALLIMIT:
               SCIPinfoMessage(iis->subscip, NULL, "Initial problem is feasible. No IIS exists.\n");
               SCIPclockStop(iis->iistime, set);
               return SCIP_OKAY;

            case SCIP_STATUS_INFEASIBLE:
               break;

            case SCIP_STATUS_UNKNOWN:
            default:
               SCIPinfoMessage(iis->subscip, NULL, "Unexpected return status %d. Failing to show (in)feasibility of initial problem. Exiting ...\n", SCIPgetStatus(iis->subscip));
               SCIPclockStop(iis->iistime, set);
               return SCIP_OKAY;
         }
      }
      else
      {
         SCIPinfoMessage(iis->subscip, NULL, "Initial solve to verify infeasibility failed.\n");
         SCIPclockStop(iis->iistime, set);
         return SCIP_OKAY;
      }
      iis->nnodes += SCIPgetNTotalNodes(iis->subscip);
      SCIP_CALL( SCIPfreeTransform(iis->subscip) );
      iis->valid = TRUE;
   }
   else if( SCIPgetStage(set->scip) == SCIP_STAGE_SOLVED )
   {
      if( SCIPgetStatus(set->scip) != SCIP_STATUS_INFEASIBLE )
      {
         SCIPinfoMessage(iis->subscip, NULL, "Initial problem does not have an infeasible status. Not performing an IIS algorithm.\n");
         SCIPclockStop(iis->iistime, set);
         return SCIP_OKAY;
      }
      iis->valid = TRUE;
   }

   /* Check for trivial infeasibility reasons */
   SCIP_CALL( checkTrivialInfeas(iis->subscip, &trivial) );

   /* Try all IIS generators */
   SCIP_CALL( SCIPgetBoolParam(set->scip, "iis/stopafterone", &stopafterone) );
   SCIP_CALL( SCIPgetBoolParam(set->scip, "iis/removebounds", &removebounds) );
   if( !trivial )
   {
      for( i = 0; i < set->niisfinders; ++i )
      {
         SCIP_IISFINDER* iisfinder;
         iisfinder = set->iisfinders[i];
         assert(iis != NULL);

         /* Recreate the subscip if one of the IIS finder algorithms has produced an invalid IS */
         if( ! iis->valid )
            SCIP_CALL( createSubscipIIS(set, iis, timelim - SCIPclockGetTime(iis->iistime), nodelim) );

         /* start timing */
         SCIPclockStart(iisfinder->iisfindertime, set);

         SCIPdebugMsg(set->scip, "----- STARTING IIS FINDER %s -----\n", iisfinder->name);
         SCIP_CALL( iisfinder->iisfinderexec(iis, iisfinder, timelim, nodelim, removebounds, silent, &result) );

         /* stop timing */
         SCIPclockStop(iisfinder->iisfindertime, set);

         assert( result == SCIP_SUCCESS || result == SCIP_DIDNOTFIND || result == SCIP_DIDNOTRUN );

         if( timelim - SCIPclockGetTime(iis->iistime) <= 0 || (nodelim != -1 && iis->nnodes > nodelim) )
         {
            SCIPdebugMsg(set->scip, "Time or node limit hit. Stopping Search.\n");
            break;
         }

         if( (stopafterone && (result == SCIP_SUCCESS)) || (iis->irreducible == TRUE) )
            break;
      }
   }

   /* Ensure the problem is irreducible if requested */
   SCIP_CALL( SCIPgetBoolParam(set->scip, "iis/minimal", &minimal) );
   if( !iis->irreducible && minimal && !(timelim - SCIPclockGetTime(iis->iistime) <= 0 || (nodelim != -1 && iis->nnodes > nodelim)) && !trivial )
   {
      SCIP_RANDNUMGEN* randnumgen;

      SCIPdebugMsg(set->scip, "----- STARTING GREEDY DELETION ALGORITHM WITH BATCHSIZE=1. ATTEMPT TO ENSURE IRREDUCIBILITY -----\n");

      if( !(iis->valid) )
         SCIP_CALL( createSubscipIIS(set, iis, timelim, nodelim) );

      SCIP_CALL( SCIPcreateRandom(set->scip, &randnumgen, 0x5EED, TRUE) );
      SCIP_CALL( SCIPexecIISfinderGreedy(iis, timelim, nodelim, removebounds, silent, randnumgen, 1e+20, FALSE, TRUE, TRUE, TRUE, -1L, 1, 1.0, &result) );
      SCIPfreeRandom(set->scip, &randnumgen);
      assert( result == SCIP_SUCCESS || result == SCIP_DIDNOTFIND || result == SCIP_DIDNOTRUN );
   }

   /* Remove redundant constraints that potentially are left over from indicator constraints,
    * that is constraints containing a variables with no bounds and that only features in a single constraint */
   nconss = SCIPgetNOrigConss(iis->subscip);
   origconss = SCIPgetOrigConss(iis->subscip);
   SCIP_CALL( SCIPduplicateBlockMemoryArray(iis->subscip, &conss, origconss, nconss) );
   for( i = 0; i < nconss; ++i )
   {
      islinear = strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])), "linear") == 0;
      if( islinear && SCIPconsGetNUses(conss[i]) <= 1)
      {
         nvars = SCIPgetNVarsLinear(iis->subscip, conss[i]);
         vars = SCIPgetVarsLinear(iis->subscip, conss[i]);
         for( j = 0; j < nvars; ++j )
         {
            if( SCIPvarGetNUses(vars[j]) <= 2 && SCIPisInfinity(iis->subscip, -SCIPvarGetLbOriginal(vars[j])) && SCIPisInfinity(iis->subscip, SCIPvarGetUbOriginal(vars[j])) )
            {
               SCIP_CALL( SCIPdelCons(iis->subscip, conss[i]) );
               break;
            }
         }
      }
   }
   SCIPfreeBlockMemoryArray(iis->subscip, &conss, nconss);

   if( !silent )
      SCIPiisfinderInfoMessage(iis, FALSE);

   SCIP_CALL( SCIPgetBoolParam(set->scip, "iis/removeunusedvars", &removeunusedvars) );
   if( removeunusedvars )
   {
      SCIP_Bool deleted;

      nvars = SCIPgetNOrigVars(iis->subscip);
      vars = SCIPgetOrigVars(iis->subscip);
      // TODO: Ask someone if this is safe. It is passing basic tests.
      for( i = nvars - 1; i >= 0; i-- )
      {
         if( SCIPvarGetNUses(vars[i]) <= 1 && SCIPvarGetLbOriginal(vars[i]) <= SCIPvarGetUbOriginal(vars[i]) )
         {
            SCIP_CALL( SCIPdelVar(iis->subscip, vars[i], &deleted) );
            assert( deleted );
         }
      }
   }

   /* stop timing */
   SCIPclockStop(iis->iistime, set);

   return SCIP_OKAY;
}

/** gets description of the IIS finder */
const char* SCIPiisfinderGetDesc(
   SCIP_IISFINDER*       iisfinder           /**< IIS finder */
   )
{
   assert(iisfinder != NULL);

   return iisfinder->desc;
}

/** copies the given IIS finder to a new scip */
SCIP_RETCODE SCIPiisfinderCopyInclude(
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(iisfinder != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( iisfinder->iisfindercopy != NULL )
   {
      SCIPsetDebugMsg(set, "including IIS finder %s in subscip %p\n", SCIPiisfinderGetName(iisfinder), (void*)set->scip);
      SCIP_CALL( iisfinder->iisfindercopy(set->scip, iisfinder) );
   }
   return SCIP_OKAY;
}

/** frees memory of IIS finder */
SCIP_RETCODE SCIPiisfinderFree(
   SCIP_IISFINDER**      iisfinder,          /**< IIS finder */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(iisfinder != NULL);

   if( *iisfinder == NULL )
      return SCIP_OKAY;

   assert(set != NULL);

   /* call destructor of IIS */
   if( (*iisfinder)->iisfinderfree != NULL )
   {
      SCIP_CALL( (*iisfinder)->iisfinderfree(set->scip, *iisfinder) );
   }

   /* free clocks */
   SCIPclockFree(&(*iisfinder)->iisfindertime);

   BMSfreeBlockMemoryArrayNull(blkmem, &(*iisfinder)->name, strlen((*iisfinder)->name)+1);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*iisfinder)->desc, strlen((*iisfinder)->desc)+1);
   BMSfreeBlockMemory(blkmem, iisfinder);

   return SCIP_OKAY;
}

/** gets user data of IIS finder */
SCIP_IISFINDERDATA* SCIPiisfinderGetData(
   SCIP_IISFINDER*       iisfinder           /**< IIS finder */
   )
{
   assert(iisfinder != NULL);

   return iisfinder->iisfinderdata;
}

/** sets user data of IIS finder; user has to free old data in advance! */
void SCIPiisfinderSetData(
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_IISFINDERDATA*   iisfinderdata       /**< new IIS finder user data */
   )
{
   assert(iisfinder != NULL);

   iisfinder->iisfinderdata = iisfinderdata;
}

/** gets priority of IIS finder */
int SCIPiisfinderGetPriority(
   SCIP_IISFINDER*       iisfinder           /**< IIS finder */
   )
{
   assert(iisfinder != NULL);

   return iisfinder->priority;
}

/** enables or disables all clocks of @p iisfinder, depending on the value of the flag */
void SCIPiisfinderEnableOrDisableClocks(
   SCIP_IISFINDER*       iisfinder,          /**< the IIS finder for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the IIS be enabled? */
   )
{
   assert(iisfinder != NULL);

   SCIPclockEnableOrDisable(iisfinder->iisfindertime, enable);
}

/** sets copy method of IIS finder */
void SCIPiisfinderSetCopy(
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_DECL_IISFINDERCOPY ((*iisfindercopy)) /**< copy method of IIS finder or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   assert(iisfinder != NULL);

   iisfinder->iisfindercopy = iisfindercopy;
}

/** sets destructor method of IIS finder */
void SCIPiisfinderSetFree(
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_DECL_IISFINDERFREE ((*iisfinderfree)) /**< destructor of IIS finder */
   )
{
   assert(iisfinder != NULL);

   iisfinder->iisfinderfree = iisfinderfree;
}

/** sets priority of IIS finder */
void SCIPiisfinderSetPriority(
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the IIS finder */
   )
{
   assert(iisfinder != NULL);
   assert(set != NULL);

   iisfinder->priority = priority;
   set->iisfinderssorted = FALSE;
}

/** gets time in seconds used in this IIS finder */
SCIP_Real SCIPiisfinderGetTime(
   SCIP_IISFINDER*       iisfinder           /**< IIS finder */
   )
{
   assert(iisfinder != NULL);

   return SCIPclockGetTime(iisfinder->iisfindertime);
}

/** prints output line during IIS calculations */
void SCIPiisfinderInfoMessage(
   SCIP_IIS*             iis,                /**< pointer to the IIS */
   SCIP_Bool             printheaders        /**< whether the headers should be printed instead of the info */
   )
{
   SCIP_VAR** vars;
   SCIP* scip;
   int i;
   int nvars;
   int nbounds = 0;
   const char* valid = iis->valid ? "yes" : "no";
   const char* irreducible = iis->irreducible ? "yes" : "no";

   scip = SCIPiisGetSubscip(iis);
   assert(scip != NULL);

   if( printheaders )
   {
      SCIPinfoMessage(scip, NULL, "  conss | bounds |  valid | irreducible |  nodes | time  \n");
      return;
   }
   if( iis->niismessagecalls % 15 == 0 )
      SCIPinfoMessage(scip, NULL, "  conss | bounds |  valid | irreducible |  nodes | time  \n");
   iis->niismessagecalls += 1;

   nvars = SCIPgetNOrigVars(scip);
   vars = SCIPgetOrigVars(scip);
   for( i = 0; i < nvars; i++ )
   {
      if( !SCIPisInfinity(scip, -SCIPvarGetLbOriginal(vars[i])) )
         ++nbounds;
      if( !SCIPisInfinity(scip, SCIPvarGetUbOriginal(vars[i])) )
         ++nbounds;
   }
   SCIPinfoMessage(scip, NULL, "%7d |%7d |%7s | %11s |%7lld | %7f\n", SCIPgetNOrigConss(scip), nbounds, valid, irreducible, SCIPiisGetNNodes(iis), SCIPiisGetTime(iis));
}

/** creates and captures a new IIS */
SCIP_RETCODE SCIPiisCreate(
   SCIP_IIS**            iis,                /**< pointer to return the created IIS */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(iis != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, iis) );

   (*iis)->subscip = NULL;
   (*iis)->subscip = NULL;
   (*iis)->varsmap = NULL;
   (*iis)->conssmap = NULL;
   SCIP_CALL( SCIPclockCreate(&(*iis)->iistime, SCIP_CLOCKTYPE_DEFAULT) );
   (*iis)->niismessagecalls = 0;
   (*iis)->nnodes = 0;
   (*iis)->valid = FALSE;
   (*iis)->irreducible = FALSE;

   return SCIP_OKAY;
}

/** releases an IIS */
SCIP_RETCODE SCIPiisFree(
   SCIP_IIS**            iis,                /**< pointer to the IIS */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(iis != NULL);
   if( *iis == NULL )
      return SCIP_OKAY;

   if( (*iis)->subscip != NULL )
   {
      SCIP_CALL( SCIPfree(&((*iis)->subscip)) );
      (*iis)->subscip = NULL;
   }

   if( (*iis)->varsmap != NULL )
   {
      SCIPhashmapFree(&((*iis)->varsmap));
      (*iis)->varsmap = NULL;
   }

   if( (*iis)->conssmap != NULL )
   {
      SCIPhashmapFree(&((*iis)->conssmap));
      (*iis)->conssmap = NULL;
   }

   SCIPclockFree(&(*iis)->iistime);

   BMSfreeBlockMemory(blkmem, iis);
   *iis = NULL;

   return SCIP_OKAY;
}

/** reset an IIS (in case one exists from a previous solve) */
SCIP_RETCODE SCIPiisReset(
   SCIP_IIS**           iis                  /**< pointer to the IIS */
   )
{
   assert(iis != NULL);
   if( *iis == NULL )
      return SCIP_OKAY;

   if( (*iis)->subscip != NULL )
   {
      SCIP_CALL( SCIPfree(&((*iis)->subscip)) );
      (*iis)->subscip = NULL;
   }

   if( (*iis)->varsmap != NULL )
   {
      SCIPhashmapFree(&((*iis)->varsmap));
      (*iis)->varsmap = NULL;
   }

   if( (*iis)->conssmap != NULL )
   {
      SCIPhashmapFree(&((*iis)->conssmap));
      (*iis)->conssmap = NULL;
   }

   SCIPclockReset((*iis)->iistime);
   (*iis)->niismessagecalls = 0;
   (*iis)->nnodes = 0;
   (*iis)->valid = FALSE;
   (*iis)->irreducible = FALSE;

   return SCIP_OKAY;
}

/** gets time in seconds used in the IIS */
SCIP_Real SCIPiisGetTime(
   SCIP_IIS*             iis                 /**< IIS */
   )
{
   assert( iis != NULL );

   return SCIPclockGetTime(iis->iistime);
}

/** Gets whether the IIS subscip is valid. */
SCIP_Bool SCIPiisGetValid(
   SCIP_IIS*             iis                 /**< IIS data structure */
   )
{
   assert( iis != NULL );

   return iis->valid;
}

/** Gets whether the IIS subscip is irreducible. */
SCIP_Bool SCIPiisGetIrreducible(
   SCIP_IIS*             iis                 /**< IIS data structure */
   )
{
   assert( iis != NULL );

   return iis->irreducible;
}

/** Gets the number of nodes in the IIS solve. */
SCIP_Longint SCIPiisGetNNodes(
   SCIP_IIS*             iis                 /**< IIS data structure */
   )
{
   assert( iis != NULL );

   return iis->nnodes;
}

/** Sets the flag that states whether the IIS subscip is valid. */
void SCIPiisSetValid(
   SCIP_IIS*             iis,                /**< IIS data structure */
   SCIP_Bool             valid               /**< The new validity status of the IIS */
   )
{
   iis->valid = valid;
}

/** Sets the flag that states whether the IIS subscip is irreducible. */
void SCIPiisSetIrreducible(
   SCIP_IIS*             iis,                /**< IIS data structure */
   SCIP_Bool             irreducible         /**< The new irreducible status of the IIS */
)
{
   iis->irreducible = irreducible;
}

/** Increments the number of nodes in the IIS solve. */
void SCIPiisAddNNodes(
   SCIP_IIS*             iis,                /**< IIS data structure */
   SCIP_Longint          nnodes              /**< The number of nodes to add to the IIS */
   )
{
   iis->nnodes += nnodes;
}

/** get the subscip of an IIS */
SCIP* SCIPiisGetSubscip(
   SCIP_IIS*            iis                  /**< pointer to the IIS */
   )
{
   return iis->subscip;
}

/** compares two IIS finders w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPiisfinderComp)
{  /*lint --e{715}*/
   return ((SCIP_IISFINDER*)elem2)->priority - ((SCIP_IISFINDER*)elem1)->priority;
}
