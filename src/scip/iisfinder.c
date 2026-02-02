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

/**@file   iisfinder.c
 * @ingroup OTHER_CFILES
 * @brief  methods for IIS finders
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/set.h"
#include "scip/clock.h"
#include "scip/misc.h"
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
   SCIP_Longint          nodelim,            /**< nodelimit */
   SCIP_Bool*            success             /**< whether the created subscip is complete */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(set != NULL);
   assert(iis != NULL);
   assert(success != NULL);

   *success = FALSE;

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
   SCIP_CALL( SCIPcopyOrig(set->scip, iis->subscip, iis->varsmap, iis->conssmap, "iis", TRUE, FALSE, TRUE, success) );

   if( !(*success) )
      return SCIP_OKAY;

   /* Remove the objective */
   vars = SCIPgetOrigVars(iis->subscip);
   nvars = SCIPgetNOrigVars(iis->subscip);
   for( i = 0; i < nvars; i++ )
      SCIP_CALL( SCIPchgVarObj(iis->subscip, vars[i], 0.0 ) );

   /* copy parameter settings */
   /** @todo: Do we really want to copy the parameter settings? */
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
   SCIP_CALL( SCIPsetSubscipsOff(iis->subscip, TRUE) );
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
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Real maxactivity;
   SCIP_Real minactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool success;
   int violatingcons;
   int nconss;
   int nvars;
   int i;
   int j;

   assert( trivial != NULL );
   *trivial = FALSE;

   /* Check for contradicting bounds */
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

   /* Check for linear max (min) activities that do not respect their lhs (rhs) */
   violatingcons = -1;
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);
   for( i = 0; i < nconss; ++i )
   {
      /**@todo generalize activity evaluation */
      /* Skip the constraint if it is not linear */
      if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])), "linear") != 0 )
         continue;

      /* Get variable information */
      nvars = SCIPgetNVarsLinear(scip, conss[i]);
      vars = SCIPgetVarsLinear(scip, conss[i]);
      coefs = SCIPgetValsLinear(scip, conss[i]);

      /* Check the left-hand side */
      lhs = SCIPconsGetLhs(scip, conss[i], &success);
      assert( success );

      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* Compute the maximum activity */
         maxactivity = 0.0;
         for( j = 0; j < nvars; ++j )
            maxactivity += coefs[j] * (coefs[j] >= 0.0 ? SCIPvarGetUbOriginal(vars[j]) : SCIPvarGetLbOriginal(vars[j]));

         /* Is the violation (maxactivity < lhs) true? */
         if( SCIPisSumLT(scip, maxactivity, lhs) )
         {
            *trivial = TRUE;
            violatingcons = i;
            break;
         }
      }

      /* Check the right-hand side */
      rhs = SCIPconsGetRhs(scip, conss[i], &success);
      assert( success );

      if( !SCIPisInfinity(scip, rhs) )
      {
         /* Compute the minimum activity */
         minactivity = 0.0;
         for( j = 0; j < nvars; ++j )
            minactivity += coefs[j] * (coefs[j] >= 0.0 ? SCIPvarGetLbOriginal(vars[j]) : SCIPvarGetUbOriginal(vars[j]));

         /* Is the violation (rhs < minactivity) true? */
         if( SCIPisSumLT(scip, rhs, minactivity) )
         {
            *trivial = TRUE;
            violatingcons = i;
            break;
         }
      }
   }

   /* Delete all constraints not relevant to the infeasibility */
   if( *trivial )
   {
      for( i = nconss - 1; i >= 0; i-- )
      {
         if( i == violatingcons )
            continue;
         SCIP_CALL( SCIPdelCons(scip, conss[i]) );
      }
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
   SCIP_Bool             enable,             /**< whether the IIS finder should be enabled */
   SCIP_DECL_IISFINDERCOPY ((*iisfindercopy)), /**< copy method of IIS finder or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_IISFINDERFREE ((*iisfinderfree)), /**< destructor of IIS finder */
   SCIP_DECL_IISFINDEREXEC ((*iisfinderexec)), /**< IIS finder execution method */
   SCIP_IISFINDERDATA*   iisfinderdata       /**< IIS finder data */
   )
{
   char priorityparamname[SCIP_MAXSTRLEN];
   char priorityparamdesc[SCIP_MAXSTRLEN];
   char enableparamname[SCIP_MAXSTRLEN];
   char enableparamdesc[SCIP_MAXSTRLEN];

   assert(iisfinder != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(iisfinderexec != NULL);

   SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, iisfinder) );

   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*iisfinder)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*iisfinder)->desc, desc, strlen(desc)+1) );
   (*iisfinder)->priority = priority;
   (*iisfinder)->enable = enable;
   (*iisfinder)->iisfindercopy = iisfindercopy;
   (*iisfinder)->iisfinderfree = iisfinderfree;
   (*iisfinder)->iisfinderexec = iisfinderexec;
   (*iisfinder)->iisfinderdata = iisfinderdata;

   /* create clocks */
   SCIP_CALL( SCIPclockCreate(&(*iisfinder)->iisfindertime, SCIP_CLOCKTYPE_DEFAULT) );

   /* add parameters */
   (void) SCIPsnprintf(priorityparamname, SCIP_MAXSTRLEN, "iis/%s/priority", name);
   (void) SCIPsnprintf(priorityparamdesc, SCIP_MAXSTRLEN, "priority of iis generation rule <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, priorityparamname, priorityparamdesc,
         &(*iisfinder)->priority, FALSE, priority, INT_MIN/4, INT_MAX/2,
         paramChgdIISfinderPriority, (SCIP_PARAMDATA*)(*iisfinder)) ); /*lint !e740*/
   (void) SCIPsnprintf(enableparamname, SCIP_MAXSTRLEN, "iis/%s/enable", name);
   (void) SCIPsnprintf(enableparamdesc, SCIP_MAXSTRLEN, "whether the iis finder <%s> should be enabled", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, enableparamname, enableparamdesc,
         &(*iisfinder)->enable, FALSE, enable, NULL, NULL) );

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
   SCIP_Bool             enable,             /**< whether the IIS finder should be enabled */
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

   SCIP_CALL_FINALLY( doIISfinderCreate(iisfinder, set, messagehdlr, blkmem, name, desc, priority, enable,
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
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_IIS* iis;
   SCIP_RESULT result = SCIP_DIDNOTFIND;
   SCIP_RETCODE retcode;
   SCIP_Real timelim;
   SCIP_Longint nodelim;
   SCIP_Bool silent;
   SCIP_Bool makeirreducible;
   SCIP_Bool stopafterone;
   SCIP_Bool removeunusedvars;
   SCIP_Bool trivial;
   SCIP_Bool islinear;
   SCIP_Bool success;
   int nconss;
   int nvars;
   int nbounds;
   int i;
   int j;

   /* exact mode is not supported */
   if( set->exact_enable )
   {
      SCIPinfoMessage(set->scip, NULL, "IIS generation does not yet support exact mode.\n");
      return SCIP_OKAY;
   }

   /* sort the iis finders by priority */
   SCIPsetSortIISfinders(set);

   /* Get the IIS data. */
   iis = SCIPgetIIS(set->scip);
   SCIP_CALL( SCIPiisReset(&iis) );
   SCIP_CALL( SCIPgetRealParam(set->scip, "iis/time", &timelim) );
   SCIP_CALL( SCIPgetLongintParam(set->scip, "iis/nodes", &nodelim) );

   /* Create the subscip used for storing the IIS */
   SCIP_CALL( createSubscipIIS(set, iis, timelim, nodelim, &success) );

   if( !success )
   {
      SCIPinfoMessage(iis->subscip, NULL, "Error copying  original problem instance. IIS generation suspended.\n");
      return SCIP_OKAY;
   }

   SCIPclockStart(iis->iistime, set);

   /* If the model is not yet shown to be infeasible then check for infeasibility */
   if( SCIPgetStage(set->scip) == SCIP_STAGE_PROBLEM )
   {
      retcode = SCIPsolve(iis->subscip);

      if( retcode != SCIP_OKAY )
      {
         SCIPinfoMessage(iis->subscip, NULL, "Error proving infeasibility of initial problem. IIS generation suspended.\n");
         SCIPclockStop(iis->iistime, set);
         return SCIP_OKAY;
      }

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
      iis->infeasible = TRUE;
   }
   else if( SCIPgetStage(set->scip) == SCIP_STAGE_SOLVED )
   {
      if( SCIPgetStatus(set->scip) != SCIP_STATUS_INFEASIBLE )
      {
         SCIPinfoMessage(iis->subscip, NULL, "Initial problem does not have an infeasible status. Not performing an IIS algorithm.\n");
         SCIPclockStop(iis->iistime, set);
         return SCIP_OKAY;
      }
      iis->infeasible = TRUE;
   }
   else
   {
      SCIPinfoMessage(iis->subscip, NULL, "Initial problem is neither in problem stage or infeasible.\n");
      return SCIP_OKAY;
   }

   /* Check for trivial infeasibility reasons */
   SCIP_CALL( checkTrivialInfeas(iis->subscip, &trivial) );
   if( trivial )
      SCIPiisSetSubscipIrreducible(iis, TRUE);

   /* Try all IIS generators */
   SCIP_CALL( SCIPgetBoolParam(set->scip, "iis/stopafterone", &stopafterone) );
   if( !trivial )
   {
      for( i = 0; i < set->niisfinders; ++i )
      {
         SCIP_IISFINDER* iisfinder;
         iisfinder = set->iisfinders[i];
         assert( iis->infeasible );

         /* skip disabled IIS finders */
         if( !iisfinder->enable )
            continue;

         /* start timing */
         SCIPclockStart(iisfinder->iisfindertime, set);

         SCIPdebugMsg(iis->subscip, "----- STARTING IIS FINDER %s -----\n", iisfinder->name);
         SCIP_CALL( iisfinder->iisfinderexec(iis, iisfinder, &result) );
         assert( result == SCIP_SUCCESS || result == SCIP_DIDNOTFIND || result == SCIP_DIDNOTRUN );

         /* stop timing */
         SCIPclockStop(iisfinder->iisfindertime, set);

         /* recreate the initial subscip if the IIS finder has produced an invalid infeasible subsystem */
         if( !iis->infeasible )
         {
            SCIP_CALL( createSubscipIIS(set, iis, timelim, nodelim, &success) );
            assert(success);
         }

         if( timelim - SCIPclockGetTime(iis->iistime) <= 0 || (nodelim != -1 && iis->nnodes > nodelim) )
         {
            SCIPdebugMsg(iis->subscip, "Time or node limit hit. Stopping Search.\n");
            break;
         }

         if( (stopafterone && (result == SCIP_SUCCESS)) || (iis->irreducible == TRUE) )
            break;
      }
   }

   /* Ensure the problem is irreducible if requested */
   SCIP_CALL( SCIPgetBoolParam(set->scip, "iis/irreducible", &makeirreducible) );
   if( !iis->irreducible && makeirreducible && !(timelim - SCIPclockGetTime(iis->iistime) <= 0 || (nodelim != -1 && iis->nnodes > nodelim)) && !trivial )
   {
      SCIPdebugMsg(iis->subscip, "----- STARTING GREEDY SINGLETON DELETION ALGORITHM. ATTEMPT TO ENSURE IRREDUCIBILITY -----\n");

      assert( iis->infeasible );
      SCIP_CALL( SCIPiisGreedyMakeIrreducible(iis) );
      assert( iis->infeasible );
   }

   /* Remove redundant constraints that potentially are left over from indicator constraints,
    * that is constraints containing a variables with no bounds and that only features in a single constraint */
   nconss = SCIPgetNOrigConss(iis->subscip);
   conss = SCIPgetOrigConss(iis->subscip);
   for( i = nconss -1; i >= 0; --i )
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

   SCIP_CALL( SCIPgetBoolParam(set->scip, "iis/removeunusedvars", &removeunusedvars) );
   if( removeunusedvars )
   {
      SCIP_Bool deleted;

      nvars = SCIPgetNOrigVars(iis->subscip);
      vars = SCIPgetOrigVars(iis->subscip);
      for( i = nvars - 1; i >= 0; i-- )
      {
         if( SCIPvarGetNUses(vars[i]) <= 1 && SCIPvarIsDeletable(vars[i]) && SCIPvarGetLbOriginal(vars[i]) <= SCIPvarGetUbOriginal(vars[i]) )
         {
            SCIP_CALL( SCIPdelVar(iis->subscip, vars[i], &deleted) );
            assert( deleted );
         }
      }
   }

   SCIP_CALL( SCIPgetBoolParam(set->scip, "iis/silent", &silent) );
   if( !silent )
      SCIPiisfinderInfoMessage(iis, FALSE);

   /* stop timing */
   SCIPclockStop(iis->iistime, set);

   if( !silent )
   {
      nbounds = 0;
      SCIPinfoMessage(set->scip, NULL, "\n");
      if( iis->infeasible && iis->irreducible )
         SCIPinfoMessage(set->scip, NULL, "IIS Status            : irreducible infeasible subsystem (IIS) found\n");
      else if( iis->infeasible )
         SCIPinfoMessage(set->scip, NULL, "IIS Status            : infeasible subsystem (IS) found\n");
      else
         SCIPinfoMessage(set->scip, NULL, "IIS Status            : failed to find an infeasible subsystem (IS)\n");
      if( iis->irreducible )
         SCIPinfoMessage(set->scip, NULL, "IIS irreducible       : yes\n");
      else
         SCIPinfoMessage(set->scip, NULL, "IIS irreducible       : no\n");
      assert( iis->subscip != NULL );
      SCIPinfoMessage(set->scip, NULL, "Generation Time (sec) : %.2f\n", SCIPiisGetTime(iis));
      SCIPinfoMessage(set->scip, NULL, "Generation Nodes      : %lld\n", SCIPiisGetNNodes(iis));
      SCIPinfoMessage(set->scip, NULL, "Num. Cons. in IIS     : %d\n", SCIPgetNOrigConss(iis->subscip));
      nvars = SCIPgetNOrigVars(iis->subscip);
      vars = SCIPgetOrigVars(iis->subscip);
      for( i = 0; i < nvars; i++ )
      {
         if( !SCIPisInfinity(iis->subscip, -SCIPvarGetLbOriginal(vars[i])) )
            ++nbounds;
         if( !SCIPisInfinity(iis->subscip, SCIPvarGetUbOriginal(vars[i])) )
            ++nbounds;
      }
      SCIPinfoMessage(set->scip, NULL, "Num. Vars. in IIS     : %d\n", nvars);
      SCIPinfoMessage(set->scip, NULL, "Num. Bounds in IIS    : %d\n", nbounds);
   }

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

/** returns whether IIS finder is enabled */
SCIP_Bool SCIPiisfinderGetEnable(
   SCIP_IISFINDER*       iisfinder           /**< IIS finder */
   )
{
   assert(iisfinder != NULL);

   return iisfinder->enable;
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

/** enables/disables IIS finder */
void SCIPiisfinderSetEnable(
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_Bool             enable              /**< whether the IIS finder should be enabled */
   )
{
   assert(iisfinder != NULL);

   iisfinder->enable = enable;
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
   const char* infeasible = iis->infeasible ? "yes" : "no";

   scip = SCIPiisGetSubscip(iis);
   assert(scip != NULL);

   if( printheaders || (iis->niismessagecalls % 15 == 0) )
   {
      SCIPinfoMessage(scip, NULL, "time(s)| node  | cons  | vars  | bounds| infeasible\n");
      if( printheaders )
         return;
   }
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
   SCIPinfoMessage(scip, NULL, "%7.1f|%7lld|%7d|%7d|%7d| %10s\n", SCIPiisGetTime(iis), SCIPiisGetNNodes(iis), SCIPgetNOrigConss(scip), nvars, nbounds, infeasible);
}

/** creates and captures a new IIS */
SCIP_RETCODE SCIPiisCreate(
   SCIP_IIS**            iis,                /**< pointer to return the created IIS */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(iis != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, iis) );

   (*iis)->subscip = NULL;
   (*iis)->subscip = NULL;
   (*iis)->varsmap = NULL;
   (*iis)->conssmap = NULL;
   SCIP_CALL( SCIPrandomCreate(&((*iis)->randnumgen), blkmem, SCIPsetInitializeRandomSeed(set, 0x5EED)) );
   SCIP_CALL( SCIPclockCreate(&(*iis)->iistime, SCIP_CLOCKTYPE_DEFAULT) );
   (*iis)->niismessagecalls = 0;
   (*iis)->nnodes = 0;
   (*iis)->infeasible = FALSE;
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

   if( (*iis)->randnumgen != NULL )
   {
      SCIPrandomFree(&((*iis)->randnumgen), blkmem);
      (*iis)->randnumgen = NULL;
   }

   SCIPclockFree(&(*iis)->iistime);

   BMSfreeBlockMemory(blkmem, iis);
   *iis = NULL;

   return SCIP_OKAY;
}

/** reset an IIS (in case one exists from a previous solve) */
SCIP_RETCODE SCIPiisReset(
   SCIP_IIS**            iis                 /**< pointer to the IIS */
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
   (*iis)->infeasible = FALSE;
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

/** Gets whether the IIS subscip is currently infeasible. */
SCIP_Bool SCIPiisIsSubscipInfeasible(
   SCIP_IIS*             iis                 /**< IIS data structure */
   )
{
   assert( iis != NULL );

   return iis->infeasible;
}

/** Gets whether the IIS subscip is irreducible. */
SCIP_Bool SCIPiisIsSubscipIrreducible(
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

/** Sets the flag that states whether the IIS subscip is currently infeasible. */
void SCIPiisSetSubscipInfeasible(
   SCIP_IIS*             iis,                /**< IIS data structure */
   SCIP_Bool             infeasible          /**< The new infeasibility status of the IIS subscip */
   )
{
   assert( iis != NULL );
   iis->infeasible = infeasible;
}

/** Sets the flag that states whether the IIS subscip is irreducible. */
void SCIPiisSetSubscipIrreducible(
   SCIP_IIS*             iis,                /**< IIS data structure */
   SCIP_Bool             irreducible         /**< The new irreducible status of the IIS */
   )
{
   assert( iis != NULL );
   iis->irreducible = irreducible;
}

/** Increments the number of nodes in the IIS solve. */
void SCIPiisAddNNodes(
   SCIP_IIS*             iis,                /**< IIS data structure */
   SCIP_Longint          nnodes              /**< The number of nodes to add to the IIS */
   )
{
   assert( iis != NULL );
   iis->nnodes += nnodes;
}

/** get the randnumgen of the IIS */
SCIP_RANDNUMGEN* SCIPiisGetRandnumgen(
   SCIP_IIS*             iis                 /**< pointer to the IIS */
   )
{
   assert( iis != NULL );
   return iis->randnumgen;
}

/** get the subscip of an IIS */
SCIP* SCIPiisGetSubscip(
   SCIP_IIS*             iis                 /**< pointer to the IIS */
   )
{
   assert( iis != NULL );
   return iis->subscip;
}

/** compares two IIS finders w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPiisfinderComp)
{  /*lint --e{715}*/
   return ((SCIP_IISFINDER*)elem2)->priority - ((SCIP_IISFINDER*)elem1)->priority;
}
