/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   conflict_resolution.c
 * @ingroup OTHER_CFILES
 * @brief  methods and datastructures for resolution-based conflict analysis
 * @author Gioni Mexi
 *
 * @todo
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG

#include "lpi/lpi.h"
#include "scip/conflict_resolution.h"
#include "scip/conflict_graphanalysis.h"
#include "scip/conflict_dualproofanalysis.h"
#include "scip/clock.h"
#include "scip/conflict.h"
#include "scip/cons.h"
#include "scip/cons_linear.h"
#include "scip/cuts.h"
#include "scip/history.h"
#include "scip/lp.h"
#include "scip/presolve.h"
#include "scip/prob.h"
#include "scip/prop.h"
#include "scip/pub_conflict.h"
#include "scip/pub_cons.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_paramset.h"
#include "scip/pub_prop.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/scip_conflict.h"
#include "scip/scip_cons.h"
#include "scip/scip_mem.h"
#include "scip/scip_sol.h"
#include "scip/scip_var.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/struct_conflict.h"
#include "scip/struct_lp.h"
#include "scip/struct_prob.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_tree.h"
#include "scip/struct_var.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/visual.h"
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif

/** resolve the conflict constraint with the reason */
static
SCIP_RETCODE resolveWithReason(
   SCIP_RESOLUTIONSET*   conflictresolutionset,    /**< conflict resolution set */
   SCIP_RESOLUTIONSET*   reasonresolutionset,      /**< reason resolution set */
   BMS_BLKMEM*           blkmem,                   /**< block memory */
   int                   residx                   /**< index of variable to resolve */
   )
{
   /* assert that the coefficients of residx have different signs in the two resolution sets */
   /* Do resolution ( maybe use SCIP_AGGRROW ? ) */

   return SCIP_OKAY;
}

/** weaken variables in the reason */
static
SCIP_RETCODE weakenVarReason(
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_SET*             set,
   SCIP_VAR*             var,
   int                   pos,
   SCIP_Bool*            valid
   )
{
   assert(resolutionset != NULL);
   assert(var != NULL);
   assert(pos >= 0 && pos < resolutionset->nnz);
   assert(valid != NULL);

   *valid = TRUE;
   /* weaken with upper bound */
   if( resolutionset->vals[pos] > 0.0 )
   {
      resolutionset->lhs -= resolutionset->vals[pos] * SCIPvarGetUbGlobal(var);
   }
   /* weaken with lower bound */
   else
   {
      assert(resolutionset->vals[pos] < 0.0);
      resolutionset->lhs -= resolutionset->vals[pos] * SCIPvarGetLbGlobal(var);
   }

   --resolutionset->nnz;

   resolutionset->vals[pos] = resolutionset->vals[resolutionset->nnz];
   resolutionset->inds[pos] = resolutionset->inds[resolutionset->nnz];
   resolutionset->vals[resolutionset->nnz] = 0.0;
   resolutionset->inds[resolutionset->nnz] = 0;

   if( SCIPsetIsInfinity(set, -resolutionset->lhs) )
      *valid = FALSE;
}

/** calculates the slack of a given set of bounds and coefficients */
static
SCIP_RETCODE getSlack(
   SCIP_RESOLUTIONSET**  resolutionset,      /**< resolution set */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to resolve */
   )
{
   /* @todo */
   return SCIP_OKAY;
}

/** calculates the minimal activity of a given set of bounds and coefficients */
static
SCIP_Real getMinActivity(
   )
{
   /* @todo */
   SCIP_Real QUAD(minact);
   return QUAD_TO_DBL(minact);
}

/** calculates the maximal activity of a given set of bounds and coefficients */
static
SCIP_Real getMaxActivity(
   )
{
   /* @todo */
   SCIP_Real QUAD(maxact);
   return QUAD_TO_DBL(maxact);
}

/** creates a resolution constraint and tries to add it to the storage */
static
SCIP_RETCODE createAndAddResolutionCons(
   )
   {
      /* @todo */
      return SCIP_OKAY;
   }

/* create resolution constraints out of resolution sets */
SCIP_RETCODE SCIPconflictFlushResolutionSets(
   )
   {
      /* @todo */
      /* Calls static method createAndAddResolutionCons() */
      return SCIP_OKAY;
   }

/** creates a resolutionset */
static
SCIP_RETCODE resolutionsetCreate(
   SCIP_RESOLUTIONSET**  resolutionset,      /**< resolution set */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(resolutionset != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, resolutionset) );
   (*resolutionset)->vals = NULL;
   (*resolutionset)->inds = NULL;
   (*resolutionset)->lhs = 0.0;
   (*resolutionset)->origlhs = 0.0;
   (*resolutionset)->origrhs = 0.0;
   (*resolutionset)->nnz = 0;
   (*resolutionset)->size = 0;
   (*resolutionset)->validdepth = 0;
   (*resolutionset)->row = NULL;
   (*resolutionset)->propside = SCIP_PROPSIDE_UNKNOWN;
   (*resolutionset)->conflicttype = SCIP_CONFTYPE_UNKNOWN;

   return SCIP_OKAY;
}

/** frees a resolutionset */
void SCIPresolutionsetFree(
   SCIP_RESOLUTIONSET**  resolutionset,      /**< resolution set */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(resolutionset != NULL);
   assert(*resolutionset != NULL);
   assert(blkmem != NULL);

   BMSfreeBlockMemoryArrayNull(blkmem, &(*resolutionset)->vals, (*resolutionset)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*resolutionset)->inds, (*resolutionset)->size);
   BMSfreeBlockMemory(blkmem, resolutionset);
   (*resolutionset) = NULL;
}


/** adds given data as row to the resolution set */
static
SCIP_RETCODE resolutionsetAddSparseData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Real*            vals,               /**< variable coefficients */
   int*                  inds,               /**< variable array */
   int                   nnz,                /**< size of variable and coefficient array */
   SCIP_Real             lhs,                /**< left-hand side of resolution set */
   SCIP_Real             origrhs,            /**< right-hand side of the row */
   SCIP_Real             origlhs,            /**< left-hand side of the row */
   SCIP_ROW*             row                 /**< pointer to row */
   )
{

   assert(resolutionset != NULL);
   assert(blkmem != NULL);

   if( resolutionset->size == 0 )
   {
      assert(resolutionset->vals == NULL);
      assert(resolutionset->inds == NULL);

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &resolutionset->vals, vals, nnz) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &resolutionset->inds, inds, nnz) );

      resolutionset->size = nnz;
   }
   else
   {
      int i;

      assert(resolutionset->vals != NULL);
      assert(resolutionset->inds != NULL);

      if( resolutionset->size < nnz )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &resolutionset->vals, resolutionset->size, nnz) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &resolutionset->inds, resolutionset->size, nnz) );
         resolutionset->size = nnz;
      }

      for( i = 0; i < nnz; i++ )
      {
         resolutionset->vals[i] = vals[i];
         resolutionset->inds[i] = inds[i];
      }
   }
   resolutionset->lhs = lhs;
   resolutionset->origrhs = origrhs;
   resolutionset->origlhs = origlhs;
   resolutionset->nnz = nnz;

   resolutionset->row = row;

   return SCIP_OKAY;
}
/** add a row to the resolutionset */
static
SCIP_RETCODE reasonResolutionsetFromRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change to resolve */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row                 /**< row to add */
   )
{
   SCIP_COL** cols;
   SCIP_VAR* var;
   SCIP_Real* vals;
   int* inds;
   SCIP_Real lhs;
   SCIP_Real origlhs;
   SCIP_Real origrhs;
   int nnz;
   int varidx;
   int i;
   SCIP_Bool changesign;

   assert(resolutionset != NULL);
   assert(set != NULL);

   nnz = SCIProwGetNNonz(row);
   assert(nnz > 0);

   origlhs = SCIProwGetLhs(row);
   origrhs = SCIProwGetRhs(row);

   SCIP_CALL( SCIPsetAllocBufferArray(set, &inds, nnz) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, nnz) );

   var = SCIPbdchginfoGetVar(bdchginfo);
   varidx = SCIPvarGetProbindex(var);


   vals = SCIProwGetVals(row);
   cols = SCIProwGetCols(row);
   for( i = 0; i < nnz; i++ )
   {
      var = SCIPcolGetVar(cols[i]);
      inds[i] = SCIPvarGetProbindex(var);
      if ( (inds[i] == varidx) )
      {
         if ( (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER && vals[i] < 0) || (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER && vals[i] > 0) )
         {
            changesign = TRUE;
         }
         else if ( (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER && vals[i] > 0) || (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER && vals[i] < 0) )
         {
            changesign = FALSE;
         }
         else
         {
            assert(FALSE);
         }
      }
   }
   if ( changesign )
   {
      lhs = -origrhs;
      for( i = 0; i < nnz; i++ )
      {
         vals[i] = -vals[i];
      }
   }


   assert(inds != NULL);
   assert(vals != NULL);


   SCIP_CALL( resolutionsetAddSparseData(scip, resolutionset, blkmem, vals, inds, nnz, lhs, origrhs, origlhs, row) );

   SCIPsetFreeBufferArray(set, &inds);
   SCIPsetFreeBufferArray(set, &vals);

   return SCIP_OKAY;
}

/** add a row to the resolutionset */
static
SCIP_RETCODE conflictResolutionsetFromRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change to resolve */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row                 /**< row to add */
   )
{
   SCIP_COL** cols;
   SCIP_VAR* var;
   SCIP_Real* vals;
   int* inds;
   SCIP_Real lhs;
   SCIP_Real origlhs;
   SCIP_Real origrhs;
   int nnz;
   int varidx;
   int i;
   SCIP_Bool changesign;

   assert(resolutionset != NULL);
   assert(set != NULL);

   nnz = SCIProwGetNNonz(row);
   assert(nnz > 0);

   origlhs = SCIProwGetLhs(row);
   origrhs = SCIProwGetRhs(row);

   SCIP_CALL( SCIPsetAllocBufferArray(set, &inds, nnz) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, nnz) );

   var = SCIPbdchginfoGetVar(bdchginfo);
   varidx = SCIPvarGetProbindex(var);


   vals = SCIProwGetVals(row);
   cols = SCIProwGetCols(row);
   for( i = 0; i < nnz; i++ )
   {
      var = SCIPcolGetVar(cols[i]);
      inds[i] = SCIPvarGetProbindex(var);
      if ( (inds[i] == varidx) )
      {
         if ( (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER && vals[i] > 0) || (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER && vals[i] < 0) )
         {
            changesign = TRUE;
         }
         else if ( (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER && vals[i] < 0) || (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER && vals[i] > 0) )
         {
            changesign = FALSE;
         }
         else
         {
            assert(FALSE);
         }
      }
   }
   if ( changesign )
   {
      lhs = -origrhs;
      for( i = 0; i < nnz; i++ )
      {
         vals[i] = -vals[i];
      }
   }


   assert(inds != NULL);
   assert(vals != NULL);


   SCIP_CALL( resolutionsetAddSparseData(scip, resolutionset, blkmem, vals, inds, nnz, lhs, origrhs, origlhs, row) );

   SCIPsetFreeBufferArray(set, &inds);
   SCIPsetFreeBufferArray(set, &vals);

   return SCIP_OKAY;
}

/** resets the data structure of a resolutionset */
static
void resolutionSetClear(
   SCIP_RESOLUTIONSET*        resolutionset            /**< resolution set */
   )
{
   assert(resolutionset != NULL);

   resolutionset->nnz = 0;
   resolutionset->lhs = 0.0;
   resolutionset->origlhs = 0.0;
   resolutionset->origrhs = 0.0;
   resolutionset->validdepth = 0;
   resolutionset->conflicttype = SCIP_CONFTYPE_UNKNOWN;
   resolutionset->propside = SCIP_PROPSIDE_UNKNOWN;
}

/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound() and
 *  SCIPconflictAddRelaxedBound(), and on success, calls the conflict handlers to create a conflict constraint;
 *  afterwards the conflict queue and the conflict set is cleared
 */

SCIP_RETCODE conflictAnalyzeResolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONS*            cons,               /**< constraint that detected the conflict */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool             mustresolve,        /**< should the conflict set only be used, if a resolution was applied? */
   int*                  nconss,             /**< pointer to store the number of generated conflict constraints */
   int*                  nconfvars           /**< pointer to store the number of variables in generated conflict constraints */
   )
{
   SCIP_RESOLUTIONSET *conflictresolutionset;
   SCIP_RESOLUTIONSET *reasonresolutionset;
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BDCHGINFO** firstuips;
   SCIP_CONFTYPE conftype;
   SCIP_ROW* conflictrow;
   SCIP_ROW* reasonrow;
   SCIP_CONS* reasoncon;

   int nfirstuips;
   int focusdepth;
   int currentdepth;
   int maxvaliddepth;
   int resolvedepth;
   int nresolutions;
   int lastconsnresolutions;
   int lastconsresoldepth;

   assert(cons != NULL);
   assert(conflict != NULL);
   assert(conflict->conflictset != NULL);
   assert(conflict->conflictset->nbdchginfos >= 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(0 <= validdepth && validdepth <= SCIPtreeGetCurrentDepth(tree));
   assert(nconss != NULL);
   assert(nconfvars != NULL);


   focusdepth = SCIPtreeGetFocusDepth(tree);
   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);
   assert(focusdepth <= currentdepth);

   resolvedepth = ((set->conf_fuiplevels >= 0 && set->conf_fuiplevels <= currentdepth)
      ? currentdepth - set->conf_fuiplevels + 1 : 0);
   assert(0 <= resolvedepth && resolvedepth <= currentdepth + 1);

   /* if we must resolve at least one bound change, find the first UIP at least in the last depth level */
   if( mustresolve )
      resolvedepth = MIN(resolvedepth, currentdepth);

   SCIP_CALL( resolutionsetCreate(&conflictresolutionset, blkmem) );

   SCIPsetDebugMsg(set, "analyzing conflict with %d+%d conflict candidates and starting conflict constraint of size %d in depth %d (resolvedepth=%d)\n",
      SCIPpqueueNElems(conflict->forcedbdchgqueue), SCIPpqueueNElems(conflict->bdchgqueue),
      conflictresolutionset->size, currentdepth, resolvedepth);

   *nconss = 0;
   *nconfvars = conflictresolutionset->size;

   /* check, whether local conflicts are allowed; however, don't generate conflict constraints that are only valid in the
    * probing path and not in the problem tree (i.e. that exceed the focusdepth)
    */
   maxvaliddepth = (set->conf_allowlocal ? MIN(currentdepth-1, focusdepth) : 0);
   if( validdepth > maxvaliddepth )
      return SCIP_OKAY;

   /* allocate temporary memory for storing first UIPs (in each depth level, at most two bound changes can be flagged
    * as UIP, namely a binary and a non-binary bound change)
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &firstuips, 2*(currentdepth+1)) ); /*lint !e647*/

   /* process all bound changes in the conflict candidate queue */
   nresolutions = 0;
   lastconsnresolutions = (mustresolve ? 0 : -1);
   lastconsresoldepth = (mustresolve ? currentdepth : INT_MAX);
   bdchginfo = conflictFirstCand(conflict);
   nfirstuips = 0;

   SCIP_CALL( resolutionsetCreate(&conflictresolutionset, blkmem) );


   /* get the corresponding conflict row */
   conflictrow = SCIPconsGetRow(scip, cons);


   if (conflictrow == NULL)
      return SCIP_OKAY;

   /* get the resolution set of the conflict row */
   SCIP_CALL( conflictResolutionsetFromRow(scip, conflictresolutionset, bdchginfo, set, blkmem, conflictrow) );
   if ( SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER )
   {
      reasoncon = SCIPbdchginfoGetInferCons(bdchginfo);
      SCIP_CALL( resolutionsetCreate(&reasonresolutionset, blkmem) );

      /* get the corresponding conflict row */
      reasonrow = SCIPconsGetRow(scip, reasoncon);

      if (reasonrow == NULL)
         return SCIP_OKAY;

      /* get the resolution set of the conflict row */
      SCIP_CALL( reasonResolutionsetFromRow(scip, reasonresolutionset, bdchginfo, set, blkmem, reasonrow) );

   }
   else
   {
      return SCIP_OKAY;
   }
   /* @todo assert negative slack of conflict constraint */



   /* check if the initial reason on debugging solution */
   SCIP_CALL( SCIPdebugCheckConflictFrontier(blkmem, set, tree->path[validdepth], \
         NULL, conflict->conflictset->bdchginfos, conflict->conflictset->relaxedbds, conflict->conflictset->nbdchginfos, \
         conflict->bdchgqueue, conflict->forcedbdchgqueue) ); /*lint !e506 !e774*/

   /* main loop */
   while( bdchginfo != NULL && validdepth <= maxvaliddepth )
   {
      break;
      /* go through conflict constraint to get variable at the latest depth level
       * if the reason for the latest is a propagation of a constraint we extract the reason
       * so that we can apply generalized resolution */

      /*  assert slack for reason >= 0 */
      /* go through variables and weaken + coef. tightening till the resolved constraint has  negative slack */

   }
   return SCIP_OKAY;
}


/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound(), and on success, calls the
 *  conflict handlers to create a conflict constraint; updates statistics for propagation conflict analysis
 */
SCIP_RETCODE SCIPconflictAnalyzeResolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONS*            cons,               /**< constraint that detected the conflict */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   int nconss;
   int nconfvars;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(prob != NULL);

   if( success != NULL )
      *success = FALSE;

   /* @todo add option in SCIPconflictApplicable for resolution conflict analysis*/
   /* check if resolution conflict analysis is applicable */
   if( !SCIPconflictApplicable(set) )
      return SCIP_OKAY;

   /* @todo check, if the conflict constraint will get too large with high probability */

   SCIPsetDebugMsg(set, "resolution based conflict analysis after infeasible propagation in depth %d\n", SCIPtreeGetCurrentDepth(tree));

   /* start timing */
   SCIPclockStart(conflict->resanalyzetime, set);

   conflict->nrescalls++;

   /* analyze the conflict set, and create a conflict constraint on success */
   SCIP_CALL( conflictAnalyzeResolution(scip, conflict, cons, blkmem, set, stat, prob, tree, FALSE, validdepth, TRUE, &nconss, &nconfvars) );
   conflict->nressuccess += (nconss > 0 ? 1 : 0);
   conflict->nresconfconss += nconss;
   conflict->nresconfvariables += nconfvars;
   if( success != NULL )
      *success = (nconss > 0);

   /* stop timing */
   SCIPclockStop(conflict->resanalyzetime, set);

   return SCIP_OKAY;
}
