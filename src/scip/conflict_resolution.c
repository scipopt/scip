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
// #define SCIP_DEBUG

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
#include "scip/scip_prob.h"
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
#include "scip/scip_numerics.h"

#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif


/** perform activity based coefficient tigthening on a row defined with a left hand side; returns TRUE if the row
 *  is redundant due to acitvity bounds
 */

static
SCIP_Bool tightenCoefLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROB*            prob,
   SCIP_Bool             localbounds,        /**< do we use local bounds? */
   SCIP_Real*            rowcoefs,           /**< array of the non-zero coefficients in the row */
   SCIP_Real*            rowlhs,             /**< the left hand side of the row */
   int*                  rowinds,            /**< array of indices of variables with a non-zero coefficient in the row */
   int*                  rownnz,             /**< the number of non-zeros in the row */
   int*                  nchgcoefs           /**< number of changed coefficients */
   )
{

   int i;
   int nintegralvars;
   SCIP_VAR** vars;
   SCIP_Real* absvals;
   SCIP_Real QUAD(minacttmp);
   SCIP_Real minact;
   SCIP_Real maxabsval = 0.0;
   SCIP_Bool redundant = FALSE;

   assert(nchgcoefs != NULL);

   QUAD_ASSIGN(minacttmp, 0.0);

   vars = SCIPprobGetVars(prob);
   nintegralvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &absvals, *rownnz) );

   assert(nchgcoefs != NULL);
   *nchgcoefs = 0;

   for( i = 0; i < *rownnz; ++i )
   {
      assert(rowinds[i] >= 0);
      assert(vars[rowinds[i]] != NULL);

      if( rowcoefs[i] > 0.0 )
      {
         SCIP_Real lb = localbounds ? SCIPvarGetLbLocal(vars[rowinds[i]]) : SCIPvarGetLbGlobal(vars[rowinds[i]]);

         if( SCIPisInfinity(scip, -lb) )
            goto TERMINATE;

         if( rowinds[i] < nintegralvars )
         {
            maxabsval = MAX(maxabsval, rowcoefs[i]);
            absvals[i] = rowcoefs[i];
         }
         else
            absvals[i] = 0.0;

         SCIPquadprecSumQD(minacttmp, minacttmp, lb * rowcoefs[i]);
      }
      else
      {
         SCIP_Real ub = localbounds ? SCIPvarGetUbLocal(vars[rowinds[i]]) : SCIPvarGetUbGlobal(vars[rowinds[i]]);

         if( SCIPisInfinity(scip, ub) )
            goto TERMINATE;

         if( rowinds[i] < nintegralvars )
         {
            maxabsval = MAX(maxabsval, -rowcoefs[i]);
            absvals[i] = -rowcoefs[i];
         }
         else
            absvals[i] = 0.0;

         SCIPquadprecSumQD(minacttmp, minacttmp, ub * rowcoefs[i]);
      }
   }

   minact = QUAD_TO_DBL(minacttmp);

   /* row is redundant if minact is infinity */
   if (SCIPisInfinity(scip, minact) )
   {
      redundant = TRUE;
      goto TERMINATE;
   }
   /* no coefficients can be tightened */
   if (SCIPisInfinity(scip, -minact) )
   {
      goto TERMINATE;
   }

   /* row is redundant in activity bounds */
   if( SCIPisFeasGE(scip, minact, *rowlhs) )
   {
      redundant = TRUE;
      goto TERMINATE;
   }

   /* terminate, because coefficient tightening cannot be performed; also excludes the case in which no integral variable is present */
   // for lhs terminate if amin + maxabsval < rowlhs
   if( SCIPisLT(scip, minact + maxabsval, *rowlhs) )
      goto TERMINATE;

   SCIPsortDownRealRealInt(absvals, rowcoefs, rowinds, *rownnz);
   SCIPfreeBufferArray(scip, &absvals);

   /* loop over the integral variables and try to tighten the coefficients; see cons_linear for more details */
   for( i = 0; i < *rownnz; ++i )
   {
      /* due to sorting, we can exit if we reached a continuous variable: all further integral variables have 0 coefficents anyway */
      if( rowinds[i] >= nintegralvars )
         break;

      assert(SCIPvarIsIntegral(vars[rowinds[i]]));

      if( rowcoefs[i] < 0.0 && SCIPisGE(scip, minact - rowcoefs[i], *rowlhs) )
      {
         SCIP_Real coef = minact - (*rowlhs);
         SCIP_Real ub = localbounds ? SCIPvarGetUbLocal(vars[rowinds[i]]) : SCIPvarGetUbGlobal(vars[rowinds[i]]);

         if( coef > rowcoefs[i] )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumDD(delta, coef, -rowcoefs[i]);
            SCIPquadprecProdQD(delta, delta, ub);

            SCIPquadprecSumQD(tmp, delta, *rowlhs);
            SCIPdebugPrintf("tightened coefficient from %g to %g; lhs changed from %g to %g; the bounds are [%g,%g]\n",
               rowcoefs[i], coef, (*rowlhs), QUAD_TO_DBL(tmp),
               localbounds ? SCIPvarGetUbLocal(vars[rowinds[i]]) : SCIPvarGetUbGlobal(vars[rowinds[i]]), ub);

            *rowlhs = QUAD_TO_DBL(tmp);

            assert(!SCIPisPositive(scip, coef));

            ++(*nchgcoefs);

            if( SCIPisNegative(scip, coef) )
            {
               SCIPquadprecSumQQ(minacttmp, minacttmp, delta);
               minact = QUAD_TO_DBL(minacttmp);
               rowcoefs[i] = coef;
            }
            else
            {
               --(*rownnz);
               rowinds[i] = rowinds[*rownnz];
               rowcoefs[i] = rowcoefs[*rownnz];
               continue;
            }
         }
      }
      else if( rowcoefs[i] > 0.0 && SCIPisGE(scip, minact + rowcoefs[i], *rowlhs) )
      {
         SCIP_Real coef = (*rowlhs) - minact;
         SCIP_Real lb = localbounds ? SCIPvarGetLbLocal(vars[rowinds[i]]) : SCIPvarGetLbGlobal(vars[rowinds[i]]);

         if( coef < rowcoefs[i] )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumDD(delta, coef, -rowcoefs[i]);
            SCIPquadprecProdQD(delta, delta, lb);

            SCIPquadprecSumQD(tmp, delta, *rowlhs);
            SCIPdebugPrintf("tightened coefficient from %g to %g; lhs changed from %g to %g; the bounds are [%g,%g]\n",
               rowcoefs[i], coef, (*rowlhs), QUAD_TO_DBL(tmp), lb,
               localbounds ? SCIPvarGetUbLocal(vars[rowinds[i]]) : SCIPvarGetUbGlobal(vars[rowinds[i]]));

            *rowlhs = QUAD_TO_DBL(tmp);

            assert(!SCIPisNegative(scip, coef));

            ++(*nchgcoefs);

            if( SCIPisPositive(scip, coef) )
            {
               SCIPquadprecSumQQ(minacttmp, minacttmp, delta);
               minact = QUAD_TO_DBL(minacttmp);
               rowcoefs[i] = coef;
            }
            else
            {
               --(*rownnz);
               rowinds[i] = rowinds[*rownnz];
               rowcoefs[i] = rowcoefs[*rownnz];
               continue;
            }
         }
      }
      else /* due to sorting we can stop completely if the precondition was not fulfilled for this variable */
         break;
   }

  TERMINATE:
   SCIPfreeBufferArrayNull(scip, &absvals);

   return redundant;
}

/** returns next conflict analysis candidate from the candidate queue without removing it */
static
SCIP_BDCHGINFO* conflictFirstCand(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   SCIP_BDCHGINFO* bdchginfo;

   assert(conflict != NULL);

   if( SCIPpqueueNElems(conflict->resforcedbdchgqueue) > 0 )
   {
      /* get next potential candidate */
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueFirst(conflict->resforcedbdchgqueue));

      /* check if this candidate is valid */
      if( bdchginfoIsInvalid(conflict, bdchginfo) )
      {
         SCIPdebugMessage("bound change info [%d:<%s> %s %g] is invaild -> pop it from the force queue\n", SCIPbdchginfoGetDepth(bdchginfo),
            SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo));

         /* pop the invalid bound change info from the queue */
         (void)(SCIPpqueueRemove(conflict->resforcedbdchgqueue));

         /* call method recursively to get next conflict analysis candidate */
         bdchginfo = conflictFirstCand(conflict);
      }
   }
   else
   {
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueFirst(conflict->resbdchgqueue));

      /* check if this candidate is valid */
      if( bdchginfo != NULL && bdchginfoIsInvalid(conflict, bdchginfo) )
      {
         SCIPdebugMessage("bound change info [%d:<%s> %s %g] is invaild -> pop it from the queue\n", SCIPbdchginfoGetDepth(bdchginfo),
            SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo));

         /* pop the invalid bound change info from the queue */
         (void)(SCIPpqueueRemove(conflict->resbdchgqueue));
         /* call method recursively to get next conflict analysis candidate */
         bdchginfo = conflictFirstCand(conflict);
      }
   }
   assert(bdchginfo == NULL || !SCIPbdchginfoIsRedundant(bdchginfo));

   return bdchginfo;
}

/** clean up the queue of bound changes. To be called after each resolution step
 * in case signs of variables are changed
 */
static
void conflictCleanUpbdchgqueue(
   SCIP_CONFLICT*        conflict,            /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RESOLUTIONSET*   resolutionset       /**< resolution set */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_VAR* var;

   int nelems;
   int deleted;
   int i;

   assert(conflict != NULL);

   /* @todo this is inefficient. Think of something better */
   if( SCIPpqueueNElems(conflict->resforcedbdchgqueue) > 0 )
   {
      for( i = 0; i < SCIPpqueueNElems(conflict->resforcedbdchgqueue); ++i )
      {
         bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->resforcedbdchgqueue)[i]);
         SCIPpqueueDelPos(conflict->resforcedbdchgqueue, i);
      }
   }
   nelems = SCIPpqueueNElems(conflict->resbdchgqueue);
   deleted = 0;
   /* @todo this is inefficient. Think of something better */
   if( nelems > 0 )
   {
      for( i = 0; i < nelems - deleted; ++i )
      {
         int j;
         SCIP_Bool idxinrow;
         SCIP_Real val;
         int idxvar;
         bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->resbdchgqueue)[i]);
         var = bdchginfo->var;
         idxvar = SCIPvarGetProbindex(var);
         idxinrow = FALSE;
         for( j = 0; j < resolutionset->nnz; j++ )
         {
            if (resolutionset->inds[j] == idxvar)
            {
               idxinrow = TRUE;
               val = resolutionset->vals[j];
               break;
            }
         }

         if ( !idxinrow || (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER && val < 0) || (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER && val > 0) )
         {
            SCIPsetDebugMsg(set, " -> Remove bound change <%s> %s %.15g [status:%d, type:%d, depth:%d, pos:%d, reason:<%s>, info:%d] \n",
               SCIPvarGetName(var),
               SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
               SCIPbdchginfoGetNewbound(bdchginfo),
               SCIPvarGetStatus(var), SCIPvarGetType(var),
               SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
               SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "branch"
               : (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
                  ? SCIPconsGetName(SCIPbdchginfoGetInferCons(bdchginfo))
                  : (SCIPbdchginfoGetInferProp(bdchginfo) != NULL ? SCIPpropGetName(SCIPbdchginfoGetInferProp(bdchginfo))
                     : "none")),
               SCIPbdchginfoGetChgtype(bdchginfo) != SCIP_BOUNDCHGTYPE_BRANCHING ? SCIPbdchginfoGetInferInfo(bdchginfo) : -1);

               SCIPpqueueDelPos(conflict->resbdchgqueue, i);
               deleted++;
               i--;
         }
      }
   }

}

/** removes and returns next conflict analysis candidate from the candidate queue */
static
SCIP_BDCHGINFO* conflictRemoveCand(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_VAR* var;

   assert(conflict != NULL);

   if( SCIPpqueueNElems(conflict->resforcedbdchgqueue) > 0 )
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueRemove(conflict->resforcedbdchgqueue));
   else
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueRemove(conflict->resbdchgqueue));

   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   /* if we have a candidate this one should be valid for the current conflict analysis */
   assert(!bdchginfoIsInvalid(conflict, bdchginfo));

   /* mark the bound change to be no longer in the conflict (it will be either added again to the conflict set or
    * replaced by resolving, which might add a weaker change on the same bound to the queue)
    */
   var = SCIPbdchginfoGetVar(bdchginfo);
   if( SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER )
   {
      var->conflictlbcount = 0;
      var->conflictrelaxedlb = SCIP_REAL_MIN;
   }
   else
   {
      assert(SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER);
      var->conflictubcount = 0;
      var->conflictrelaxedub = SCIP_REAL_MAX;
   }

   return bdchginfo;
}


/** return TRUE if generalized resolution conflict analysis is applicable */
SCIP_Bool SCIPconflictResolutionApplicable(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   /* check, if generalized resolution conflict analysis is enabled */
   if( !set->conf_enable || !set->conf_usegeneralres )
      return FALSE;

   return TRUE;
}


/** gets number of conflict constraints detected in resolution conflict analysis */
SCIP_Longint SCIPconflictGetNResConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nresconfconss;
}

/** gets number of calls to resolution conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNResSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nressuccess;
}

/** gets number of calls to resolution conflict analysis */
SCIP_Longint SCIPconflictGetNResCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nrescalls;
}

/** creates a resolution set */
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
   (*resolutionset)->conflicttype = SCIP_CONFTYPE_UNKNOWN;

   return SCIP_OKAY;
}

/** creates and clears the resolution set */
SCIP_RETCODE SCIPconflictInitResolutionset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(conflict != NULL);
   assert(blkmem != NULL);

   SCIP_CALL( resolutionsetCreate(&conflict->resolutionset, blkmem) );
   SCIP_CALL( resolutionsetCreate(&conflict->reasonset, blkmem) );

   return SCIP_OKAY;
}

/** frees a resolution set */
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

/** resets the data structure of a resolution set */
static
void resolutionSetClear(
   SCIP_RESOLUTIONSET*   resolutionset       /**< resolution set */
   )
{
   assert(resolutionset != NULL);

   resolutionset->nnz = 0;
   resolutionset->lhs = 0.0;
   resolutionset->origlhs = 0.0;
   resolutionset->origrhs = 0.0;
   resolutionset->validdepth = 0;
   resolutionset->conflicttype = SCIP_CONFTYPE_UNKNOWN;
}

/** weaken variables in the reason */
static
SCIP_RETCODE weakenVarReason(
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_SET*             set,
   SCIP_VAR*             var,
   int                   pos
   )
{
   assert(resolutionset != NULL);
   assert(var != NULL);
   assert(pos >= 0 && pos < resolutionset->nnz);

   /* weaken with upper bound */
   if( SCIPsetIsGT(set, resolutionset->vals[pos], 0.0) )
   {
      resolutionset->lhs -= resolutionset->vals[pos] * SCIPvarGetUbGlobal(var);
   }
   /* weaken with lower bound */
   /* @todo check this */
   else
   {
      assert( SCIPsetIsLT(set, resolutionset->vals[pos], 0.0) );
      resolutionset->lhs -= resolutionset->vals[pos] * SCIPvarGetLbGlobal(var);
   }

   --resolutionset->nnz;

   resolutionset->vals[pos] = resolutionset->vals[resolutionset->nnz];
   resolutionset->inds[pos] = resolutionset->inds[resolutionset->nnz];

   return SCIP_OKAY;
}

/* Removes a variable with zero coefficient in the resolutionset */
static
void resolutionsetRemoveZeroVar(
   SCIP_RESOLUTIONSET*   resolutionset,
   SCIP_SET*             set,
   int                   pos
   )
{
   assert(resolutionset != NULL);
   assert(pos >= 0 && pos < resolutionset->nnz);
   assert(SCIPsetIsZero(set, resolutionset->vals[pos]));

   --resolutionset->nnz;
   resolutionset->vals[pos] = resolutionset->vals[resolutionset->nnz];
   resolutionset->inds[pos] = resolutionset->inds[resolutionset->nnz];
}

/** return the values of variable coefficients in the resolutionset */
static
SCIP_Real* resolutionsetGetVals(
   SCIP_RESOLUTIONSET*   resolutionset       /**< resolution set */
   )
{
   assert(resolutionset != NULL);

   return resolutionset->vals;
}

/** return the left-hand side of the resolutionset */
static
SCIP_Real resolutionsetGetLhs(
   SCIP_RESOLUTIONSET*   resolutionset       /**< resolution set */
   )
{
   assert(resolutionset != NULL);

   return resolutionset->lhs;
}

/** returns the number of non zeros in the resolutionset */
static
int resolutionsetGetNNzs(
   SCIP_RESOLUTIONSET*   resolutionset       /**< resolution set */
   )
{
   assert(resolutionset != NULL);

   return resolutionset->nnz;
}

/** calculates the slack of a given set of bounds and coefficients */
static
SCIP_Real getSlack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_RESOLUTIONSET*   resolutionset       /**< resolution set */
   )
{
   SCIP_VAR** vars;
   SCIP_Real QUAD(slack);
   int i;

   vars = SCIPprobGetVars(prob);
   assert(vars != NULL);

   QUAD_ASSIGN(slack, 0.0);

   for( i = 0; i < resolutionsetGetNNzs(resolutionset); i++ )
   {
      SCIP_Real coef;
      SCIP_Real QUAD(delta);
      int v;
      v = resolutionset->inds[i];

      assert(SCIPvarGetProbindex(vars[v]) == v);

      coef = resolutionset->vals[i];

      if( coef > 0.0 )
      {
         SCIP_Real bnd;
         int nubchgs;
         SCIP_BDCHGIDX *bdchgidx;

         nubchgs = SCIPvarGetNBdchgInfosUb(vars[v]);
         bdchgidx = SCIPvarGetLastBdchgIndex(vars[v]);

         /* @todo bdchgidx should be the one from the currently resolving variable! */
         if ( nubchgs > 0 )
         {
            bnd = SCIPgetVarUbAtIndex(scip, vars[v], bdchgidx, TRUE);
         }
         /* @todo if no bound change for the given variable should we use the local or global upper bound? */
         else bnd = SCIPvarGetUbLocal(vars[v]);

         SCIPquadprecProdDD(delta, coef, bnd);
      }
      else
      {
         SCIP_Real bnd;
         int nlbchgs;
         SCIP_BDCHGIDX *bdchgidx;

         nlbchgs = SCIPvarGetNBdchgInfosLb(vars[v]);
         bdchgidx = SCIPvarGetLastBdchgIndex(vars[v]);
         /* @todo bdchgidx should be the one from the currently resolving variable! */
         if ( nlbchgs > 0 )
         {
            bnd = SCIPgetVarLbAtIndex(scip, vars[v], bdchgidx, TRUE);
         }
         /* @todo if no bound change for the given variable should we use the local or global upper bound? */
         else bnd = SCIPvarGetLbLocal(vars[v]);

         SCIPquadprecProdDD(delta, coef, bnd);
      }
      SCIPquadprecSumQQ(slack, slack, delta);
   }
   SCIPquadprecSumQD(slack, slack, -resolutionset->lhs);
   return QUAD_TO_DBL(slack);
}

/** for every variable in the row, except the inferred variable, add bound changes */
static
SCIP_RETCODE addConflictBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROB*            prob,               /**< problem */
   SCIP_SET*             set,                /**< set */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< conflicting resolution set */
   SCIP_VAR*             infervar,           /**< variable that was deduced, or NULL */
   SCIP_BDCHGIDX*        inferbdchgidx       /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{

   /* scan through the row and add bound changes that make the constraint infeasible */
   /* @todo refine this by considering the row activity and stop adding bounds when infeasibility is detected */
   SCIP_VAR** vars;
   int i;

   vars = SCIPprobGetVars(prob);
   assert(vars != NULL);

   for( i = 0; i < resolutionsetGetNNzs(resolutionset); i++ )
   {
      SCIP_Real coef;
      int v;
      v = resolutionset->inds[i];

      assert(SCIPvarGetProbindex(vars[v]) == v);

      coef = resolutionset->vals[i];

      if( coef > 0.0 )
      {
         SCIP_Real bnd;
         if ( SCIPvarGetNBdchgInfosUb(vars[v]) > 0 )
         {
            /* @todo should the bdchgindex be inferbdchgidx or SCIPvarGetLastBdchgIndex(vars[v]) */
            bnd = SCIPgetVarUbAtIndex(scip, vars[v], inferbdchgidx, TRUE);
            if ( SCIPsetIsLT(set, bnd, SCIPvarGetUbGlobal(vars[v])) )
            {
               SCIP_CALL( SCIPaddConflictUb(scip, vars[v], inferbdchgidx) );
            }
         }
      }
      else
      {
         SCIP_Real bnd;
         if ( SCIPvarGetNBdchgInfosLb(vars[v]) > 0 )
         {
            /* @todo should the bdchgindex be inferbdchgidx or SCIPvarGetLastBdchgIndex(vars[v]) */
            bnd = SCIPgetVarLbAtIndex(scip, vars[v], inferbdchgidx, TRUE);
            if ( SCIPsetIsGT(set, bnd, SCIPvarGetLbGlobal(vars[v])) )
            {
               SCIP_CALL( SCIPaddConflictLb(scip, vars[v], inferbdchgidx) );
            }
         }
      }
   }
   return SCIP_OKAY;
}

/** creates a resolution constraint and tries to add it to the storage */
static
SCIP_RETCODE createAndAddResolutionCons(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set to add to the tree */
   int                   insertdepth,        /**< depth level at which the conflict set should be added */
   SCIP_Bool*            success             /**< pointer to store whether the addition was successful */
   )
{

   SCIP_VAR** vars;
   SCIP_VAR** consvars;
   SCIP_CONS* cons;
   SCIP_CONS* upgdcons;

   char consname[SCIP_MAXSTRLEN];

   SCIP_Real* vals;
   SCIP_Real lhs;
   int i;

   vars = SCIPprobGetVars(transprob);
   assert(vars != NULL);

   /* @todo use conflict->resolutionsets for multiple conflicts per round */

   vals = resolutionsetGetVals(resolutionset);

   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, resolutionsetGetNNzs(resolutionset)) );

   lhs = resolutionsetGetLhs(resolutionset);

   for( i = 0; i < resolutionsetGetNNzs(resolutionset); ++i )
   {
      consvars[i] = vars[resolutionset->inds[i]];
   }



   /* create a constraint out of the conflict set */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "confres%" SCIP_LONGINT_FORMAT, conflict->nresconfconss);
   SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, resolutionsetGetNNzs(resolutionset), consvars, vals, lhs, SCIPsetInfinity(set),
         FALSE, set->conf_separate, FALSE, FALSE, TRUE, (SCIPnodeGetDepth(tree->path[resolutionset->validdepth]) > 0 ), FALSE, set->conf_dynamic, set->conf_removable, FALSE) );

   /* try to automatically convert a linear constraint into a more specific and more specialized constraint */
   SCIP_CALL( SCIPupgradeConsLinear(scip, cons, &upgdcons) );
   if( upgdcons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      cons = upgdcons;
   }

   /* add conflict to SCIP */
   SCIP_CALL( SCIPaddConflict(scip, tree->path[insertdepth], cons, tree->path[resolutionset->validdepth], SCIP_CONFTYPE_RESOLUTION, FALSE) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &consvars);


   return SCIP_OKAY;
}

/** create resolution constraints out of resolution sets */
SCIP_RETCODE SCIPconflictFlushResolutionSets(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,         /**< clique table data structure */
   SCIP_RESOLUTIONSET*   resolutionset       /**< resolution set to add to the tree */
   )
{

   int focusdepth;
   int maxsize;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(tree != NULL);

   /* @todo use conflict->resolutionsets for multiple conflicts per round */

   focusdepth = SCIPtreeGetFocusDepth(tree);
   assert(focusdepth <= SCIPtreeGetCurrentDepth(tree));
   assert(SCIPtreeGetCurrentDepth(tree) == tree->pathlen-1);

   /* calculate the maximal size of each accepted conflict set */
   maxsize = conflictCalcMaxsize(set, transprob);

   SCIPsetDebugMsg(set, "flushing %d resolution sets at focus depth %d (maxsize: %d)\n",
      1, focusdepth, maxsize);

   /* @todo insert the conflict sets at the corresponding nodes */
   /* @todo add a check for the number of resolution sets */

   assert(resolutionset != NULL);
   assert(0 <= resolutionset->validdepth);
   assert( getSlack(scip, set, transprob, resolutionset) < 0 );

   /* if the resolution set is empty, the node and its sub tree in the conflict set's valid depth can be
    * cut off completely
    * @todo should check this! May be wrong
    */
   if( resolutionsetGetNNzs(resolutionset) == 0 )
   {
      SCIPsetDebugMsg(set, " -> empty resolution set in depth %d cuts off sub tree at depth %d\n",
         focusdepth, resolutionset->validdepth);

      SCIP_CALL( SCIPnodeCutoff(tree->path[resolutionset->validdepth], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
      return SCIP_OKAY;
   }
   /* @todo long conflicts can be used for repropagation */
   else
   {
      SCIP_Bool success;

      /* @todo use the right insert depth and not valid depth */
      SCIP_CALL( createAndAddResolutionCons(conflict, blkmem, scip, set, stat, transprob, origprob, \
                     tree, reopt, lp, branchcand, eventqueue, cliquetable, resolutionset, resolutionset->validdepth, &success) );
      SCIPsetDebugMsg(set, " -> resolution set added (cdpt:%d, fdpt:%d, insert:%d, valid:%d, conf: -, reprop: - , len:%d):\n",
         SCIPtreeGetCurrentDepth(tree), SCIPtreeGetFocusDepth(tree),
         resolutionset->validdepth, resolutionset->validdepth, resolutionsetGetNNzs(resolutionset));

   }

   return SCIP_OKAY;
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
   SCIP_ROW*             row,                /**< pointer to row */
   SCIP_Bool             reverse             /**< reverse coefficients */

   )
{
   int i;

   assert(resolutionset != NULL);
   assert(blkmem != NULL);

   if( resolutionset->size == 0 )
   {
      assert(resolutionset->vals == NULL);
      assert(resolutionset->inds == NULL);

      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &resolutionset->vals, nnz) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &resolutionset->inds, nnz) );
      resolutionset->size = nnz;
   }
   else
   {
      assert(resolutionset->vals != NULL);
      assert(resolutionset->inds != NULL);

      if( resolutionset->size < nnz )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &resolutionset->vals, resolutionset->size, nnz) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &resolutionset->inds, resolutionset->size, nnz) );
         resolutionset->size = nnz;
      }
   }

   if ( reverse )
   {
      for( i = 0; i < nnz; i++ )
      {
         resolutionset->vals[i] = -vals[i];
         resolutionset->inds[i] = inds[i];
      }
   }
   else
   {
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

   return SCIP_OKAY;
}

/** resolve the conflict constraint with the reason */
static
SCIP_RETCODE resolveWithReason(
   SCIP*                 scip,                     /**< SCIP */
   SCIP_SET*             set,                      /**< global SCIP settings */
   SCIP_RESOLUTIONSET*   conflictresolutionset,    /**< conflict resolution set */
   SCIP_RESOLUTIONSET*   reasonresolutionset,      /**< reason resolution set */
   BMS_BLKMEM*           blkmem,                   /**< block memory */
   int                   residx,                   /**< index of variable to resolve */
   SCIP_Bool             wasresolved,              /**< resolution has been applied */
   SCIP_Bool*            success                   /**< apply resolution */
   )
{
   int i;
   SCIP_Real coefconf;
   SCIP_Real coefreas;

   int newsize;
   int cidx;
   int previousnnz;

   SCIP_Bool idxinconflict;
   SCIP_Bool idxinreason;

   idxinconflict = FALSE;
   idxinreason = FALSE;
   *success = FALSE;

   for( i = 0; i < resolutionsetGetNNzs(conflictresolutionset); i++ )
   {
      if (conflictresolutionset->inds[i] == residx)
      {
         idxinconflict = TRUE;
         coefconf = conflictresolutionset->vals[i];
         // assert(wasresolved || !SCIPsetIsZero(set, coefconf));
         break;
      }
   }

   for( i = 0; i < resolutionsetGetNNzs(reasonresolutionset); i++ )
   {
      if (reasonresolutionset->inds[i] == residx)
      {
         idxinreason = TRUE;
         coefreas = reasonresolutionset->vals[i];
         break;
      }
   }

   assert(wasresolved || !SCIPsetIsZero(set, coefconf));
   assert(!SCIPsetIsZero(set, coefreas));

   if ( wasresolved )
   {
      if ( SCIPsetIsZero(set, coefconf) )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }
      assert((idxinconflict && idxinreason));
      /** @todo this should be replace by an assert. We have to update the current bdchginfo and
       * nextbdchginfo before adding the assert. At the moment we just stop */
      if (SCIPsetIsGE(set, coefconf * coefreas, 0.0))
      {
         *success = FALSE;
         return SCIP_OKAY;
      }
   }

   assert(coefconf * coefreas < 0);
   assert((idxinconflict && idxinreason));

   /* @todo rethink this */
   /* in case of integer coefficients we can use scm for smaller coefficients */
   if ( SCIPsetIsEQ(set, coefconf, floor(coefconf)) && SCIPsetIsEQ(set, coefreas, floor(coefreas)))
   {
      SCIP_Longint scm;
      scm = SCIPcalcSmaComMul(fabs(coefconf), fabs(coefreas));
      coefconf = scm/fabs(coefconf);
      coefreas = scm/fabs(coefreas);
   }
   else
   {
      SCIP_Real tmp;
      tmp = coefconf;
      coefconf = coefreas;
      coefreas = tmp;
   }

   SCIP_UNUSED(idxinconflict);
   SCIP_UNUSED(idxinreason);

   newsize = resolutionsetGetNNzs(conflictresolutionset) + resolutionsetGetNNzs(reasonresolutionset);
   if ( conflictresolutionset->size < newsize )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictresolutionset->inds, conflictresolutionset->size, newsize ) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictresolutionset->vals, conflictresolutionset->size, newsize ) );
      conflictresolutionset->size = newsize;
   }

   cidx = 0;
   previousnnz = resolutionsetGetNNzs(conflictresolutionset);

   /* multiply conflict by coefreas */
   for( i = 0; i < resolutionsetGetNNzs(conflictresolutionset); i++ )
   {
      conflictresolutionset->vals[i] = fabs(coefconf) * conflictresolutionset->vals[i];
   }
   /* multiply reason by coefconf */
   for( i = 0; i < resolutionsetGetNNzs(reasonresolutionset); i++ )
   {
      reasonresolutionset->vals[i] = fabs(coefreas) * reasonresolutionset->vals[i];
   }


   i = 0;
   /* add conflict and reason resolution sets */
   while ( i < resolutionsetGetNNzs(reasonresolutionset) )
   {
      if (cidx >= previousnnz)
      {
         conflictresolutionset->inds[conflictresolutionset->nnz] = reasonresolutionset->inds[i];
         conflictresolutionset->vals[conflictresolutionset->nnz] = reasonresolutionset->vals[i];
         conflictresolutionset->nnz++;
         i++;
      }
      else if (reasonresolutionset->inds[i] == conflictresolutionset->inds[cidx])
      {
         conflictresolutionset->vals[cidx] = conflictresolutionset->vals[cidx] + reasonresolutionset->vals[i];
         cidx++;
         i++;
      }
      else if (reasonresolutionset->inds[i] > conflictresolutionset->inds[cidx])
      {
         conflictresolutionset->vals[cidx] = conflictresolutionset->vals[cidx];
         cidx++;
      }
      else if (reasonresolutionset->inds[i] < conflictresolutionset->inds[cidx])
      {
         conflictresolutionset->inds[conflictresolutionset->nnz] = reasonresolutionset->inds[i];
         conflictresolutionset->vals[conflictresolutionset->nnz] = reasonresolutionset->vals[i];
         conflictresolutionset->nnz++;
         i++;
      }
   }
   conflictresolutionset->lhs = fabs(coefconf) * conflictresolutionset->lhs + fabs(coefreas) * reasonresolutionset->lhs;

   for( i = 0; i < resolutionsetGetNNzs(conflictresolutionset); i++ )
   {
      if (SCIPsetIsZero(set, conflictresolutionset->vals[i] ))
         resolutionsetRemoveZeroVar(conflictresolutionset, set, i);
   }

   SCIPsetDebugMsg(set, "Nonzeros in resolved constraint: %d \n", resolutionsetGetNNzs(conflictresolutionset));

   /* sort for linear time resolution */
   SCIPsortIntReal(conflictresolutionset->inds, conflictresolutionset->vals, resolutionsetGetNNzs(conflictresolutionset));
   *success = TRUE;

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
   SCIP_Bool isincon;

   assert(resolutionset != NULL);
   assert(set != NULL);

   nnz = SCIProwGetNNonz(row);
   assert(nnz > 0);

   origlhs = SCIProwGetLhs(row) - SCIProwGetConstant(row);
   origrhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);

   var = SCIPbdchginfoGetVar(bdchginfo);
   varidx = SCIPvarGetProbindex(var);

   vals = SCIProwGetVals(row);
   cols = SCIProwGetCols(row);

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &inds, nnz) );

   isincon = FALSE;
   for( i = 0; i < nnz; i++ )
   {
      var = SCIPcolGetVar(cols[i]);
      inds[i] = SCIPvarGetProbindex(var);
      if ( inds[i] == varidx )
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
         isincon = TRUE;
      }
   }
   assert(isincon);
   SCIP_UNUSED(isincon);

   if ( changesign )
      lhs = -origrhs;
   else
      lhs = origlhs;

   assert(inds != NULL);
   assert(vals != NULL);

   SCIP_CALL( resolutionsetAddSparseData(scip, resolutionset, blkmem, vals, inds, nnz, lhs, origrhs, origlhs, row, changesign) );

   BMSfreeBlockMemoryArray(blkmem, &inds, nnz);

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
   SCIP_Bool isincon;
   assert(resolutionset != NULL);
   assert(set != NULL);

   nnz = SCIProwGetNNonz(row);
   assert(nnz > 0);

   origlhs = SCIProwGetLhs(row) - SCIProwGetConstant(row);
   origrhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);

   var = SCIPbdchginfoGetVar(bdchginfo);
   varidx = SCIPvarGetProbindex(var);

   vals = SCIProwGetVals(row);
   cols = SCIProwGetCols(row);

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &inds, nnz) );

   isincon = FALSE;
   for( i = 0; i < nnz; i++ )
   {
      var = SCIPcolGetVar(cols[i]);
      inds[i] = SCIPvarGetProbindex(var);
      if ( inds[i] == varidx )
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
         isincon = TRUE;
      }
   }
   assert(isincon);
   SCIP_UNUSED(isincon);

   if ( changesign )
      lhs = -origrhs;
   else
      lhs = origlhs;

   assert(inds != NULL);
   assert(vals != NULL);

   SCIP_CALL( resolutionsetAddSparseData(scip, resolutionset, blkmem, vals, inds, nnz, lhs, origrhs, origlhs, row, changesign) );

   BMSfreeBlockMemoryArray(blkmem, &inds, nnz);
   return SCIP_OKAY;
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
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
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
   SCIP_BDCHGINFO* nextbdchginfo;
   SCIP_BDCHGIDX* bdchgidx;
   SCIP_VAR** vars;
   SCIP_ROW* conflictrow;
   SCIP_ROW* reasonrow;
   SCIP_CONS* reasoncon;

   int bdchgdepth;
   int focusdepth;
   int currentdepth;
   int maxvaliddepth;
   int i;
   SCIP_Bool addconstraint;
   int nchgcoefs;
   int nressteps;

   SCIP_VAR* vartoresolve;
   int residx;

   SCIP_Real conflictslack;
   SCIP_Real reasonslack;
   SCIP_Bool successresolution;

   assert(cons != NULL);
   assert(conflict != NULL);
   assert(conflict->conflictset != NULL);
   assert(conflict->conflictset->nbdchginfos >= 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(0 <= validdepth && validdepth <= SCIPtreeGetCurrentDepth(tree));
   assert(nconss != NULL);
   assert(nconfvars != NULL);

   vars = SCIPprobGetVars(transprob);
   assert(vars != NULL);

   resolutionSetClear(conflict->resolutionset);
   resolutionSetClear(conflict->reasonset);
   conflictresolutionset = conflict->resolutionset;
   reasonresolutionset = conflict->reasonset;

   focusdepth = SCIPtreeGetFocusDepth(tree);
   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);
   assert(focusdepth <= currentdepth);

   *nconss = 0;
   *nconfvars = 0;
   /* check, whether local conflicts are allowed; however, don't generate conflict constraints that are only valid in the
    * probing path and not in the problem tree (i.e. that exceed the focusdepth)
    */
   maxvaliddepth = (set->conf_allowlocal ? MIN(currentdepth-1, focusdepth) : 0);
   if( validdepth > maxvaliddepth )
      return SCIP_OKAY;

   /* get the corresponding conflict row */
   conflictrow = SCIPconsGetRow(scip, cons);

   /* if no row exists conflict analysis is not applicable */
   if (conflictrow == NULL)
   {
      SCIPsetDebugMsg(set, "Conflict analysis not applicable since no row is available for the conflict constraint \n");
      return SCIP_OKAY;
   }

   SCIPsetDebugMsg(set, "Conflict constraint: %s \n", SCIPconsGetName(cons));

   /* last bound change that led to infeasibility */
   bdchginfo = conflictFirstCand(conflict);

   if (bdchginfo == NULL)
   {
      SCIPsetDebugMsg(set, "Conflict analysis not applicable since bdchginfo is empty \n");
      return SCIP_OKAY;
   }

   /* remove the last bound change */
   bdchginfo = conflictRemoveCand(conflict);
   bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
   bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);

   vartoresolve = bdchginfo->var;
   residx = SCIPvarGetProbindex(vartoresolve);

SCIPsetDebugMsg(set, " -> First bound change to resolve <%s> %s %.15g [status:%d, type:%d, depth:%d, pos:%d, reason:<%s>, info:%d] \n",
   SCIPvarGetName(vartoresolve),
   SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
   SCIPbdchginfoGetNewbound(bdchginfo),
   SCIPvarGetStatus(vartoresolve), SCIPvarGetType(vartoresolve),
   SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
   SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "branch"
   : (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
      ? SCIPconsGetName(SCIPbdchginfoGetInferCons(bdchginfo))
      : (SCIPbdchginfoGetInferProp(bdchginfo) != NULL ? SCIPpropGetName(SCIPbdchginfoGetInferProp(bdchginfo))
         : "none")),
   SCIPbdchginfoGetChgtype(bdchginfo) != SCIP_BOUNDCHGTYPE_BRANCHING ? SCIPbdchginfoGetInferInfo(bdchginfo) : -1);


   /* if the bound change was not infered by a (reason) constraint then conflict analysis is not applicable */
   if ( SCIPbdchginfoGetChgtype(bdchginfo) != SCIP_BOUNDCHGTYPE_CONSINFER )
   {
      SCIPsetDebugMsg(set, "Conflict analysis not applicable since no reason row exists \n");
      return SCIP_OKAY;
   }

   conflictresolutionset->validdepth = validdepth;

   /* get the resolution set of the conflict row */
   SCIP_CALL( conflictResolutionsetFromRow(scip, conflictresolutionset, bdchginfo, set, blkmem, conflictrow) );


   conflictslack = getSlack(scip, set, transprob, conflictresolutionset);

   /* apply coefficient tightening for the conflict constraint */
   /* @todo add parameter for this */
   tightenCoefLhs(scip, transprob, FALSE, conflictresolutionset->vals, &conflictresolutionset->lhs, conflictresolutionset->inds, &conflictresolutionset->nnz, &nchgcoefs);

   /* sort for linear time resolution */
   SCIPsortIntReal(conflictresolutionset->inds, conflictresolutionset->vals, resolutionsetGetNNzs(conflictresolutionset));

   /* slack should be negative. Are there cases where this is not the case? */
   /* @todo assert(conflictslack < 0); */
   if ( conflictslack >= 0 )
      return SCIP_OKAY;

   nressteps = 0;
   addconstraint = FALSE;
   /* main loop */
   while( bdchginfo != NULL && validdepth <= maxvaliddepth )
   {

      if ( SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER )
      {
         reasoncon = SCIPbdchginfoGetInferCons(bdchginfo);

         /* get the corresponding reason row */
         reasonrow = SCIPconsGetRow(scip, reasoncon);

         /* if no row exists conflict analysis is not applicable */
         if (reasonrow == NULL)
            goto TERMINATE;
         /* get the resolution set of the conflict row */
         SCIP_CALL( reasonResolutionsetFromRow(scip, reasonresolutionset, bdchginfo, set, blkmem, reasonrow) );

         /* @todo Treat negated variables differently to avoid this */
         if (SCIPisInfinity(scip, -reasonresolutionset->lhs) || SCIPisInfinity(scip, reasonresolutionset->lhs))
            goto TERMINATE;
         reasonslack = getSlack(scip, set, transprob, reasonresolutionset);
         SCIPsetDebugMsg(set, "conflict resolution set: nnzs: %d, slack: %f \n", resolutionsetGetNNzs(conflictresolutionset), conflictslack);
         SCIPsetDebugMsg(set, "reason resolution set: nnzs: %d, slack: %f \n", resolutionsetGetNNzs(reasonresolutionset), reasonslack);

         /* @todo assert(reasonslack >= 0); */
         i = 0;
         while ( i < resolutionsetGetNNzs(reasonresolutionset) )
         {
            double currentslack;
            SCIP_VAR* vartoweaken;
            SCIP_Bool success;

            currentslack = reasonslack;
            success = FALSE;
            vartoweaken = vars[reasonresolutionset->inds[i]];

            if ( reasonresolutionset->inds[i] != residx )
            {
               if( reasonresolutionset->vals[i] > 0.0 )
               {
                  int nubchgs;
                  nubchgs = SCIPvarGetNBdchgInfosUb(vars[reasonresolutionset->inds[i]]);
                  if ( nubchgs == 0 )
                  {
                     weakenVarReason(reasonresolutionset, set, vartoweaken, i);
                     success = TRUE;
                  }
               }
               else
               {
                  int nlbchgs;
                  nlbchgs = SCIPvarGetNBdchgInfosLb(vars[reasonresolutionset->inds[i]]);

                  assert( reasonresolutionset->vals[i] < 0.0);
                  if ( nlbchgs == 0 )
                  {
                     weakenVarReason(reasonresolutionset, set, vartoweaken, i);
                     success = TRUE;
                  }
               }
            }
            if (success)
            {
               tightenCoefLhs(scip, transprob, FALSE, reasonresolutionset->vals, &reasonresolutionset->lhs, reasonresolutionset->inds, &reasonresolutionset->nnz, &nchgcoefs);
               reasonslack = getSlack(scip, set, transprob, reasonresolutionset);
               assert(SCIPsetIsLE(set, reasonslack, currentslack) );
               if ( SCIPsetIsLT(set, reasonslack, currentslack))
                  SCIPsetDebugMsg(set, "Reduced slack of reason from %f to %f\n", currentslack, reasonslack);
            }
            else
            i++;
         }
         /* sort for linear time resolution*/
         SCIPsortIntReal(reasonresolutionset->inds, reasonresolutionset->vals, resolutionsetGetNNzs(reasonresolutionset));

         /* terminate if the sum of slacks is positive */
         if ( conflictslack + reasonslack >= 0 )
            goto TERMINATE;

#ifdef SCIP_DEBUG
      {
         int v;

         SCIPsetDebugMsgPrint(set, "Conflict row: ");
         for( i = 0; i < resolutionsetGetNNzs(conflictresolutionset); i++ )
         {
            v = conflictresolutionset->inds[i];
            assert(SCIPvarGetProbindex(vars[v]) == v);
            SCIPsetDebugMsgPrint(set, "%f<%s> ", conflictresolutionset->vals[i], SCIPvarGetName(vars[v]));
         }
         SCIPsetDebugMsgPrint(set, ">= %f \n", conflictresolutionset->lhs);

         SCIPsetDebugMsgPrint(set, "Reason row: ");
         for( i = 0; i < resolutionsetGetNNzs(reasonresolutionset); i++ )
         {
            v = reasonresolutionset->inds[i];
            assert(SCIPvarGetProbindex(vars[v]) == v);
            SCIPsetDebugMsgPrint(set, "%f<%s> ", reasonresolutionset->vals[i], SCIPvarGetName(vars[v]));
         }
         SCIPsetDebugMsgPrint(set, ">= %f \n", reasonresolutionset->lhs);

         SCIPsetDebugMsg(set, " -> Resolve bound change <%s> %s %.15g [status:%d, type:%d, depth:%d, pos:%d, reason:<%s>, info:%d] \n",
            SCIPvarGetName(vartoresolve),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo),
            SCIPvarGetStatus(vartoresolve), SCIPvarGetType(vartoresolve),
            SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
            SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "branch"
            : (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
               ? SCIPconsGetName(SCIPbdchginfoGetInferCons(bdchginfo))
               : (SCIPbdchginfoGetInferProp(bdchginfo) != NULL ? SCIPpropGetName(SCIPbdchginfoGetInferProp(bdchginfo))
                  : "none")),
            SCIPbdchginfoGetChgtype(bdchginfo) != SCIP_BOUNDCHGTYPE_BRANCHING ? SCIPbdchginfoGetInferInfo(bdchginfo) : -1);
      }
#endif
         /* resolution step */
         resolveWithReason(scip, set, conflictresolutionset, reasonresolutionset, blkmem, residx, addconstraint, &successresolution);
         nressteps++;

         if (!successresolution)
            goto TERMINATE;

         addconstraint = TRUE;

#ifdef SCIP_DEBUG
      {
         int v;
         SCIPsetDebugMsgPrint(set, "Resolved row: ");
         for( i = 0; i < resolutionsetGetNNzs(conflictresolutionset); i++ )
         {
            v = conflictresolutionset->inds[i];
            assert(SCIPvarGetProbindex(vars[v]) == v);
            SCIPsetDebugMsgPrint(set, "%f<%s> ", conflictresolutionset->vals[i], SCIPvarGetName(vars[v]));
         }
         SCIPsetDebugMsgPrint(set, ">= %f \n", conflictresolutionset->lhs);
      }
#endif
      }
      else
      {
         goto TERMINATE;
      }

      conflictslack = getSlack(scip, set, transprob, conflictresolutionset);
      tightenCoefLhs(scip, transprob, FALSE, conflictresolutionset->vals, &conflictresolutionset->lhs, conflictresolutionset->inds, &conflictresolutionset->nnz, &nchgcoefs);
      if (nchgcoefs > 0)
      {
         SCIPsetDebugMsg(set, "Tightened %d coefficients in the resolved constraint \n", nchgcoefs);
      }
      /* sort for linear time resolution */
      SCIPsortIntReal(conflictresolutionset->inds, conflictresolutionset->vals, resolutionsetGetNNzs(conflictresolutionset));

      SCIPsetDebugMsg(set, "Slack of resolved row: %f \n", conflictslack);

      /* terminate after at most nressteps resolution iterations */
      if (set->conf_maxnumressteps > 0 && nressteps >= set->conf_maxnumressteps)
         goto TERMINATE;

      conflictCleanUpbdchgqueue(conflict, set, conflictresolutionset);
      addConflictBounds(scip, transprob, set, conflictresolutionset, vartoresolve, bdchgidx);

      bdchginfo = conflictFirstCand(conflict);
      if (bdchginfo == NULL)
         goto TERMINATE;
      bdchginfo = conflictRemoveCand(conflict);
      bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
      vartoresolve = bdchginfo->var;
      residx = SCIPvarGetProbindex(vartoresolve);

      nextbdchginfo = conflictFirstCand(conflict);

      /* if this is a UIP we stop resolving with the reason */
      if( nextbdchginfo == NULL || SCIPbdchginfoGetDepth(nextbdchginfo) != bdchgdepth )
         goto TERMINATE;

      SCIPsetDebugMsg(set, " -> Next bound change <%s> %s %.15g [status:%d, type:%d, depth:%d, pos:%d, reason:<%s>, info:%d] \n",
      SCIPvarGetName(vartoresolve),
      SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
      SCIPbdchginfoGetNewbound(bdchginfo),
      SCIPvarGetStatus(vartoresolve), SCIPvarGetType(vartoresolve),
      SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
      SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "branch"
      : (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
         ? SCIPconsGetName(SCIPbdchginfoGetInferCons(bdchginfo))
         : (SCIPbdchginfoGetInferProp(bdchginfo) != NULL ? SCIPpropGetName(SCIPbdchginfoGetInferProp(bdchginfo))
            : "none")),
      SCIPbdchginfoGetChgtype(bdchginfo) != SCIP_BOUNDCHGTYPE_BRANCHING ? SCIPbdchginfoGetInferInfo(bdchginfo) : -1);

      SCIPsetDebugMsg(set, " -> after next bound change <%s> %s %.15g [depth:%d, pos:%d, reason:<%s>] \n",
         " ",
         SCIPbdchginfoGetBoundtype(nextbdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
         SCIPbdchginfoGetNewbound(nextbdchginfo),
         SCIPbdchginfoGetDepth(nextbdchginfo), SCIPbdchginfoGetPos(bdchginfo),
         SCIPbdchginfoGetChgtype(nextbdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "branch"
         : (SCIPbdchginfoGetChgtype(nextbdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
            ? SCIPconsGetName(SCIPbdchginfoGetInferCons(nextbdchginfo))
            : (SCIPbdchginfoGetInferProp(nextbdchginfo) != NULL ? SCIPpropGetName(SCIPbdchginfoGetInferProp(nextbdchginfo))
               : "none")));


   }

  TERMINATE:
   if ( addconstraint && SCIPsetIsLT(set, conflictslack, 0.0) )
   {
      SCIPconflictFlushResolutionSets(conflict, blkmem, scip, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, cliquetable, conflictresolutionset);
      (*nconss)++;
      (*nconfvars) = resolutionsetGetNNzs(conflictresolutionset);
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
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{

   int nconss;
   int nconfvars;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(origprob != NULL);
   assert(transprob != NULL);

   if( success != NULL )
      *success = FALSE;

   /* check if generalized resolution conflict analysis is applicable */
   if( !SCIPconflictResolutionApplicable(set) )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "resolution based conflict analysis after infeasible propagation in depth %d\n", SCIPtreeGetCurrentDepth(tree));

   /* start timing */
   SCIPclockStart(conflict->resanalyzetime, set);

   conflict->nrescalls++;
   conflict->bdchgonlyresqueue = TRUE;
   /* analyze the conflict set, and create a conflict constraint on success */
   SCIP_CALL( conflictAnalyzeResolution(scip, conflict, cons, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, cliquetable, FALSE, validdepth, TRUE, &nconss, &nconfvars) );

   conflict->nressuccess += (nconss > 0 ? 1 : 0);
   conflict->nresconfconss += nconss;
   conflict->nresconfvariables += nconfvars;
   if( success != NULL )
      *success = (nconss > 0);

   SCIPpqueueClear(conflict->resbdchgqueue);
   SCIPpqueueClear(conflict->resforcedbdchgqueue);
   conflict->bdchgonlyresqueue = FALSE;

   /* stop timing */
   SCIPclockStop(conflict->resanalyzetime, set);
   SCIPsetDebugMsg(set, "resolution based conflict analysis added %d constraints \n \n", nconss);

   return SCIP_OKAY;
}
