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
 * @brief  methods and datastructures for generalized resolution-based conflict
 * analysis
 * @author Gioni Mexi
 *
 * @todo Description of the algorithm
 * @todo implement resolution for propagators that are linearisable
 * @todo use weakening with division and/orcMIR to reduce the slack of the
 * reason constraint
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

// #define SCIP_STATISTIC
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

#define BOUNDSWITCH                0.51 /**< threshold for bound switching - see cuts.c */
#define POSTPROCESS               FALSE /**< apply postprocessing to the cut - see cuts.c */
#define USEVBDS                   FALSE /**< use variable bounds - see cuts.c */
#define ALLOWLOCAL                FALSE /**< allow to generate local cuts - see cuts. */
#define MINFRAC                   0.05  /**< minimal fractionality of floor(rhs) - see cuts.c */
#define MAXFRAC                   0.999 /**< maximal fractionality of floor(rhs) - see cuts.c */

/** creates a copy of the given resolution set, allocating an additional amount of memory */
static
SCIP_RETCODE resolutionsetCopy(
   SCIP_RESOLUTIONSET**  targetresolutionset,/**< pointer to store the resolution set */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_RESOLUTIONSET*   sourceresolutionset /**< source resolution set */
   )
{
   int targetsize;

   assert(targetresolutionset != NULL);
   assert(sourceresolutionset != NULL);

   targetsize = sourceresolutionset->nnz;
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, targetresolutionset) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*targetresolutionset)->inds, targetsize) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*targetresolutionset)->vals, targetsize) );

   /* copy all data from source to target */
   BMScopyMemoryArray((*targetresolutionset)->inds, sourceresolutionset->inds, targetsize);
   BMScopyMemoryArray((*targetresolutionset)->vals, sourceresolutionset->vals, targetsize);

   (*targetresolutionset)->nnz = targetsize;
   (*targetresolutionset)->size = targetsize;
   (*targetresolutionset)->lhs = sourceresolutionset->lhs;
   (*targetresolutionset)->origlhs = sourceresolutionset->origlhs;
   (*targetresolutionset)->origrhs = sourceresolutionset->origrhs;
   (*targetresolutionset)->slack = sourceresolutionset->slack;
   (*targetresolutionset)->coefquotient = sourceresolutionset->coefquotient;
   (*targetresolutionset)->nnz = sourceresolutionset->nnz;
   (*targetresolutionset)->validdepth = sourceresolutionset->validdepth;
   (*targetresolutionset)->conflictdepth = sourceresolutionset->conflictdepth;
   (*targetresolutionset)->repropdepth = sourceresolutionset->repropdepth;
   (*targetresolutionset)->conflicttype = sourceresolutionset->conflicttype;

   return SCIP_OKAY;
}

/** replaces a resolution set by another; allocate an additional amount of memory if needed */
static
SCIP_RETCODE resolutionsetReplace(
   SCIP_RESOLUTIONSET*   targetresolutionset,/**< pointer to store the resolution set */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_RESOLUTIONSET*   sourceresolutionset /**< source resolution set */
   )
{
   int sourcesize;
   int targetsize;

   assert(targetresolutionset != NULL);
   assert(sourceresolutionset != NULL);

   sourcesize = sourceresolutionset->size;
   targetsize = targetresolutionset->size;

   /* allocate additional memory for the inds and vals arrays if needed */
   if( targetsize < sourcesize )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &targetresolutionset->vals, targetsize, sourcesize) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &targetresolutionset->inds, targetsize, sourcesize) );
      }
   targetresolutionset->size = MAX(sourcesize, targetsize);
   /* copy all data from source to target */
   BMScopyMemoryArray(targetresolutionset->inds, sourceresolutionset->inds, sourcesize);
   BMScopyMemoryArray(targetresolutionset->vals, sourceresolutionset->vals, sourcesize);

   targetresolutionset->nnz = sourceresolutionset->nnz;
   targetresolutionset->lhs = sourceresolutionset->lhs;
   targetresolutionset->origlhs = sourceresolutionset->origlhs;
   targetresolutionset->origrhs = sourceresolutionset->origrhs;
   targetresolutionset->slack = sourceresolutionset->slack;
   targetresolutionset->coefquotient = sourceresolutionset->coefquotient;
   targetresolutionset->nnz = sourceresolutionset->nnz;
   targetresolutionset->validdepth = sourceresolutionset->validdepth;
   targetresolutionset->conflictdepth = sourceresolutionset->conflictdepth;
   targetresolutionset->repropdepth = sourceresolutionset->repropdepth;
   targetresolutionset->conflicttype = sourceresolutionset->conflicttype;

   return SCIP_OKAY;
}

/** resizes resolutionsets array to be able to store at least num entries */
static
SCIP_RETCODE conflictEnsureResolutionsetsMem(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);

   if( num > conflict->resolutionsetssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conflict->resolutionsets, newsize) );
      conflict->resolutionsetssize = newsize;
   }
   assert(num <= conflict->resolutionsetssize);

   return SCIP_OKAY;
}

/** add a resolutionset to the list of all resolutionsets */
static
SCIP_RETCODE conflictInsertResolutionset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RESOLUTIONSET**  resolutionset       /**< resolutionset to add */
   )
{
   assert(conflict != NULL);
   assert(resolutionset != NULL);

   /* insert resolution into the resolutionsets array */
   SCIP_CALL( conflictEnsureResolutionsetsMem(conflict, set, conflict->nresolutionsets + 1) );

   SCIPsetDebugMsg(set, "inserting resolution set (valid depth: %d, conf depth: %d, reprop depth: %d):\n",
                   (*resolutionset)->validdepth, (*resolutionset)->conflictdepth, (*resolutionset)->repropdepth);

   conflict->resolutionsets[conflict->nresolutionsets] = *resolutionset;
   ++conflict->nresolutionsets;

   *resolutionset = NULL; /* ownership of pointer is now in the resolutionsets array */

   return SCIP_OKAY;
}

/** perform activity based coefficient tightening on a row defined with a left hand side; returns if the row
 *  is redundant due to activity bounds
 */
static
SCIP_RETCODE tightenCoefLhs(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Bool             localbounds,        /**< do we use local bounds? */
   SCIP_Real*            rowcoefs,           /**< array of the non-zero coefficients in the row */
   SCIP_Real*            rowlhs,             /**< the left hand side of the row */
   int*                  rowinds,            /**< array of indices of variables with a non-zero coefficient in the row */
   int*                  rownnz,             /**< the number of non-zeros in the row */
   int*                  nchgcoefs,          /**< number of changed coefficients */
   SCIP_Bool*            redundant           /**< pointer to store whether the row is redundant */
   )
{
   /* @todo slack update during coefficient tightening */
   int i;
   int nintegralvars;
   SCIP_VAR** vars;
   SCIP_Real* absvals;
   SCIP_Real QUAD(minacttmp);
   SCIP_Real minact;
   SCIP_Real maxabsval = 0.0;

   assert(nchgcoefs != NULL);

   QUAD_ASSIGN(minacttmp, 0.0);

   vars = SCIPprobGetVars(prob);
   nintegralvars = SCIPgetNVars(set->scip) - SCIPgetNContVars(set->scip);
   SCIP_CALL_ABORT( SCIPallocBufferArray(set->scip, &absvals, *rownnz) );

   assert(nchgcoefs != NULL);
   *nchgcoefs = 0;

   if (redundant != NULL)
      *redundant = FALSE;

   for( i = 0; i < *rownnz; ++i )
   {
      assert(rowinds[i] >= 0);
      assert(vars[rowinds[i]] != NULL);

      if( rowcoefs[i] > 0.0 )
      {
         SCIP_Real lb = localbounds ? SCIPvarGetLbLocal(vars[rowinds[i]]) : SCIPvarGetLbGlobal(vars[rowinds[i]]);

         if( SCIPisInfinity(set->scip, -lb) )
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

         if( SCIPisInfinity(set->scip, ub) )
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
   if (SCIPisInfinity(set->scip, minact) )
   {
      if (redundant != NULL)
         *redundant = TRUE;
      goto TERMINATE;
   }
   /* no coefficients can be tightened */
   if (SCIPisInfinity(set->scip, -minact) )
   {
      goto TERMINATE;
   }

   /* row is redundant in activity bounds */
   if( SCIPisFeasGE(set->scip, minact, *rowlhs) )
   {
      if (redundant != NULL)
         *redundant = TRUE;
      goto TERMINATE;
   }

   /* terminate, because coefficient tightening cannot be performed; also excludes the case in which no integral variable is present */
   /* for lhs terminate if minact + maxabsval < rowlhs */
   if( SCIPisLT(set->scip, minact + maxabsval, *rowlhs) )
      goto TERMINATE;

   SCIPsortDownRealRealInt(absvals, rowcoefs, rowinds, *rownnz);
   SCIPfreeBufferArray(set->scip, &absvals);

   /* loop over the integral variables and try to tighten the coefficients; see cons_linear for more details */
   for( i = 0; i < *rownnz; ++i )
   {
      /* due to sorting, we can exit if we reached a continuous variable: all further integral variables have 0 coefficents anyway */
      if( rowinds[i] >= nintegralvars )
         break;

      assert(SCIPvarIsIntegral(vars[rowinds[i]]));

      if( rowcoefs[i] < 0.0 && SCIPisGE(set->scip, minact - rowcoefs[i], *rowlhs) )
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

            assert(!SCIPisPositive(set->scip, coef));

            ++(*nchgcoefs);

            if( SCIPisNegative(set->scip, coef) )
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
      else if( rowcoefs[i] > 0.0 && SCIPisGE(set->scip, minact + rowcoefs[i], *rowlhs) )
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

            assert(!SCIPisNegative(set->scip, coef));

            ++(*nchgcoefs);

            if( SCIPisPositive(set->scip, coef) )
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
   SCIPfreeBufferArrayNull(set->scip, &absvals);

   return SCIP_OKAY;
}

/* returns whether a bound change is resolvable or not */
static
SCIP_Bool bdchginfoIsResolvable(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to check */
   )
{
   SCIP_CONSHDLR *conshdlr;
   SCIP_BOUNDCHGTYPE bdchgtype;
   const char* conshdlrname;

   bdchgtype = SCIPbdchginfoGetChgtype(bdchginfo);
   if (bdchgtype == SCIP_BOUNDCHGTYPE_BRANCHING)
      return FALSE;
   else if (bdchgtype == SCIP_BOUNDCHGTYPE_PROPINFER)
   {
      /* todo some propagators can be resolved */
      return FALSE;
   }
   assert(bdchgtype == SCIP_BOUNDCHGTYPE_CONSINFER);
   conshdlr = SCIPconsGetHdlr(SCIPbdchginfoGetInferCons(bdchginfo));
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return TRUE;
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      return TRUE;
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      return TRUE;
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      return TRUE;
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      return TRUE;
   }
   return FALSE;
}

/** returns whether there exists a resolvable bound change or not */
static
SCIP_Bool existsResolvablebdchginfo(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   int i;

   /* loop through bound change and check if there exists a resolvable bound change */
   for( i = 0; i < SCIPpqueueNElems(conflict->resbdchgqueue); ++i )
   {
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->resbdchgqueue)[i]);
      if (bdchginfoIsResolvable(bdchginfo))
         return TRUE;
   }
   return FALSE;
}

/** returns next conflict analysis candidate from the candidate queue without removing it */
static
SCIP_BDCHGINFO* conflictFirstCand(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   SCIP_BDCHGINFO* bdchginfo;

   assert(conflict != NULL);

   bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueFirst(conflict->resbdchgqueue));

   /* check if this candidate is valid */
   if( bdchginfo != NULL && bdchginfoIsInvalid(conflict, bdchginfo) )
   {
      SCIPdebugMessage("bound change info [%d:<%s> %s %g] is invalid -> pop it from the queue\n",
         SCIPbdchginfoGetDepth(bdchginfo),
         SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
         SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
         SCIPbdchginfoGetNewbound(bdchginfo));

      /* pop the invalid bound change info from the queue */
      (void)(SCIPpqueueRemove(conflict->resbdchgqueue));
      /* call method recursively to get next conflict analysis candidate */
      bdchginfo = conflictFirstCand(conflict);
   }
   assert(bdchginfo == NULL || !SCIPbdchginfoIsRedundant(bdchginfo));

   return bdchginfo;
}

/** clean up the queue of bound changes. To be called after each resolution step
 * in case signs of variables are changed which means that some bdchgs may not be relevant any more
 */
static
void conflictCleanUpbdchgqueue(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
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

   nelems = SCIPpqueueNElems(conflict->resbdchgqueue);
   deleted = 0;

   /* @todo this is inefficient */
   if( nelems > 0 )
   {
      for( i = 0; i < nelems; ++i )
      {
         int j;
         SCIP_Bool idxinrow;
         SCIP_Real val;
         int idxvar;

         bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->resbdchgqueue)[i - deleted]);
         var = bdchginfo->var;
         idxvar = SCIPvarGetProbindex(var);
         idxinrow = FALSE;
         val = 0.0;
         for( j = 0; j < resolutionset->nnz; j++ )
         {
            if (resolutionset->inds[j] == idxvar)
            {
               idxinrow = TRUE;
               val = resolutionset->vals[j];
               assert(!SCIPsetIsZero(set, val));
               break;
            }
         }

         if ( !idxinrow || (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER && val < 0) ||
            (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER && val > 0) )
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

               /* mark the bound change to be no longer in the conflict (it will be either added again to the resolution set or
               * replaced by resolving, which might add a weaker change on the same bound to the queue)
               */
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

               SCIPpqueueDelPos(conflict->resbdchgqueue, i - deleted);
               deleted++;
               if (i - deleted >= nelems) break;
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
   (*resolutionset)->slack = 0.0;
   (*resolutionset)->coefquotient = 0.0;
   (*resolutionset)->nnz = 0;
   (*resolutionset)->size = 0;
   (*resolutionset)->validdepth = 0;
   (*resolutionset)->conflictdepth = 0;
   (*resolutionset)->repropdepth = 0;
   (*resolutionset)->conflicttype = SCIP_CONFTYPE_UNKNOWN;
   (*resolutionset)->usescutoffbound = FALSE;

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
   SCIP_CALL( resolutionsetCreate(&conflict->prevresolutionset, blkmem) );

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
   resolutionset->slack = 0.0;
   resolutionset->coefquotient = 0.0;
   resolutionset->validdepth = 0;
   resolutionset->conflictdepth = 0;
   resolutionset->repropdepth = 0;
   resolutionset->conflicttype = SCIP_CONFTYPE_UNKNOWN;
   resolutionset->usescutoffbound = FALSE;
}

/** weaken variables in the reason */
static
void weakenVarReason(
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to weaken */
   int                   pos                 /**< position of coefficient in resolutionset */
   )
{
   assert(resolutionset != NULL);
   assert(var != NULL);
   assert(pos >= 0 && pos < resolutionset->nnz);

   SCIPdebugMessage("weaken variable <%s> in reason \n", SCIPvarGetName(var));
   /* weaken with global upper bound */
   if( SCIPsetIsGT(set, resolutionset->vals[pos], 0.0) )
   {
      assert( SCIPsetIsEQ(set, SCIPvarGetUbGlobal(var), SCIPvarGetUbLocal(var)) );
      resolutionset->lhs -= resolutionset->vals[pos] * SCIPvarGetUbGlobal(var);
   }
   /* weaken with global lower bound */
   else
   {
      assert( SCIPsetIsLT(set, resolutionset->vals[pos], 0.0) );
      assert( SCIPsetIsEQ(set, SCIPvarGetLbGlobal(var), SCIPvarGetLbLocal(var)) );
      resolutionset->lhs -= resolutionset->vals[pos] * SCIPvarGetLbGlobal(var);
   }

   --resolutionset->nnz;
   resolutionset->vals[pos] = resolutionset->vals[resolutionset->nnz];
   resolutionset->inds[pos] = resolutionset->inds[resolutionset->nnz];
}

/* Removes a variable with zero coefficient in the resolutionset */
static
void resolutionsetRemoveZeroVar(
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   pos                 /**< position of coefficient in resolutionset */
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

/* returns if the variable index is in the conflict resolution set */
static
SCIP_Bool varIdxInResolutionset(
   SCIP_RESOLUTIONSET*   resolutionset,      /**< conflict resolution set */
   int                   varidx              /**< variable index to check */
   )
{
   int i;

   assert(resolutionset != NULL);
   assert(varidx >= 0);

   for( i = 0; i < resolutionsetGetNNzs(resolutionset); ++i )
   {
      if( resolutionset->inds[i] == varidx )
         return TRUE;
   }
   return FALSE;
}

/* returns if the variable index is in the indices array */
static
SCIP_Bool varIdxInArray(
   int*                  inds,               /**< array of variable indices */
   int                   ninds,              /**< number of variable indices in array */
   int                   varidx              /**< variable index to check */
   )
{
   int i;

   assert(inds != NULL);
   assert(varidx >= 0);

   for( i = 0; i < ninds; ++i )
   {
      if( inds[i] == varidx )
         return TRUE;
   }
   return FALSE;
}

/** returns the quotient of the largest and smallest value in an array */
static
SCIP_Real getQuotLargestSmallestCoef(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            vals,               /**< array of values */
   int                   nnz                 /**< number of nonzeros */
   )
   {
      int i;
      SCIP_Real minval;
      SCIP_Real maxval;

      assert( vals != NULL);

      if ( nnz == 0 )
         return 0.0;

      minval = SCIPsetInfinity(set);
      maxval = -SCIPsetInfinity(set);

      for ( i = 0; i < nnz; i++)
      {
         minval = MIN(minval, vals[i]);
         maxval = MAX(maxval, vals[i]);
      }
      return REALABS(maxval / minval);
   }

/** calculates the slack (maxact - rhs) for a resolutionset given a set of bounds and coefficients */
static
SCIP_Real getSlack(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_BDCHGIDX *       currbdchgidx,       /**< index of current bound change */
   SCIP_Real*            fixbounds,          /**< array of fixed bounds */
   int*                  fixinds             /**< array of indices of fixed variables */
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
      SCIP_Real bound;
      SCIP_Real QUAD(delta);
      int v;
      v = resolutionset->inds[i];

      assert(SCIPvarGetProbindex(vars[v]) == v);

      coef = resolutionset->vals[i];
      bound = 0.0;
      /* get the latest bound change before currbdchgidx */
      if( coef > 0.0 )
      {
         if ( fixinds != NULL && fixinds[v] == 1 ) /* if the variable is fixed */
         {
            bound = fixbounds[v];
         }
         else
         {
            bound = SCIPgetVarUbAtIndex(set->scip, vars[v], currbdchgidx, TRUE);
         }
         SCIPquadprecProdDD(delta, coef, bound);
      }
      else
      {
         if (fixinds != NULL && fixinds[v] == -1) /* if the variable is fixed */
         {
            bound = fixbounds[v];
         }
         else
         {
            bound = SCIPgetVarLbAtIndex(set->scip, vars[v], currbdchgidx, TRUE);
         }
         SCIPquadprecProdDD(delta, coef, bound);
      }
      SCIPquadprecSumQQ(slack, slack, delta);
   }
   SCIPquadprecSumQD(slack, slack, -resolutionset->lhs);
   return QUAD_TO_DBL(slack);
}

/** return the coefficient of a variable in the resolution set */
static
SCIP_Bool getCoefInResolutionSet(
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   int                   varidx,             /**< index of variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   )
{
   int i;
   SCIP_Bool found;

   assert(resolutionset != NULL);
   assert(resolutionset->nnz > 0);

   found = FALSE;
   *coef = 0.0;
   for( i = 0; i < resolutionset->nnz; i++ )
   {
      if( resolutionset->inds[i] == varidx )
      {
         found = TRUE;
         *coef = resolutionset->vals[i];
         break;
      }
   }
   return found;
}

/** reference solution besed on the conflict resolution set to use in cMIR */
static
SCIP_RETCODE computeReferenceSolutionConflict(
   SCIP_RESOLUTIONSET*   resolutionset,      /**< conflict resolution set */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SOL*             sol,                /**< solution to use as reference */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change information */
   SCIP_Bool*            success             /**< pointer to store whether the reference solution was successfully computed */
   )
{
   SCIP_VAR** vars;
   SCIP_BDCHGIDX* bdchgidx;
   int nvars;
   int i;

   assert(resolutionset != NULL);
   assert(set != NULL);
   assert(prob != NULL);
   assert(sol != NULL);

   *success = TRUE;
   vars = SCIPprobGetVars(prob);
   assert(vars != NULL);
   nvars = SCIPprobGetNVars(prob);


   bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
   /* initialize with average solution */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[i], SCIPvarGetAvgSol(vars[i])) );
   }

   /* set all variables that are part of the resolution set to their active local bounds */
   for( i = 0; i < resolutionsetGetNNzs(resolutionset); i++ )
   {
      int v;
      SCIP_Real val;
      SCIP_Real lb;
      SCIP_Real ub;

      v = resolutionset->inds[i];
      lb = SCIPvarGetLbGlobal(vars[v]);
      ub = SCIPvarGetUbGlobal(vars[v]);
      /* take the negation of the value since the aggregation of the resolution set is a <= constraint */
      val = -resolutionset->vals[i];

      /* stop if both bounds are infinite */
      if( SCIPsetIsInfinity(set, -lb) && SCIPsetIsInfinity(set, ub) )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }
      else if( !SCIPsetIsInfinity(set, -lb) || !SCIPsetIsInfinity(set, ub) )
      {

         SCIP_Real meanbound;

         /* take the negation of the value since the aggregation is a <= constraint */
         meanbound = ( ub + lb ) / 2.0;

         if( val > 0.0 )
         {
            SCIP_Real locallb;

            locallb = SCIPvarGetLbAtIndex(vars[v], bdchgidx, TRUE);
            if( SCIPsetIsGE(set, locallb, meanbound) )
               SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[v], SCIPvarGetUbGlobal(vars[v])) );
            else
               SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[v], SCIPvarGetLbGlobal(vars[v])) );
         }
         else
         {
            SCIP_Real localub;

            localub = SCIPvarGetUbAtIndex(vars[v], bdchgidx, TRUE);
            if( SCIPsetIsGE(set, localub, meanbound) )
               SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[v], SCIPvarGetUbGlobal(vars[v])) );
            else
               SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[v], SCIPvarGetLbGlobal(vars[v])) );
         }
      }
      else if( !SCIPsetIsInfinity(set, -lb) )
      {
         SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[v], lb) );
      }
      else
      {
         assert(!SCIPsetIsInfinity(set, ub));
         SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[v], ub) );
      }
   }

   return SCIP_OKAY;
}

/** reference solution besed on the reason resolution set to use in cMIR */
static
SCIP_RETCODE computeReferenceSolutionReason(
   SCIP_RESOLUTIONSET*   resolutionset,      /**< conflict resolution set */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SOL*             sol,                /**< solution to use as reference */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change information */
   SCIP_Bool*            success             /**< pointer to store whether the reference solution was successfully computed */
   )
{
   SCIP_VAR** vars;
   SCIP_BDCHGIDX* bdchgidx;
   int nvars;
   int i;

   assert(resolutionset != NULL);
   assert(set != NULL);
   assert(prob != NULL);
   assert(sol != NULL);

   *success = TRUE;
   vars = SCIPprobGetVars(prob);
   assert(vars != NULL);
   nvars = SCIPprobGetNVars(prob);
   bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);

   /* initialize with average solution */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[i], SCIPvarGetAvgSol(vars[i])) );
   }

   /* set all variables that are part of the resolution set to their active local bounds */
   for( i = 0; i < resolutionsetGetNNzs(resolutionset); i++ )
   {
      int v;
      SCIP_Real val;
      SCIP_Real locallb;
      SCIP_Real localub;

      v = resolutionset->inds[i];
      locallb = SCIPvarGetLbAtIndex(vars[v], bdchgidx, TRUE);
      localub = SCIPvarGetUbAtIndex(vars[v], bdchgidx, TRUE);
      /* take the negation of the value since the aggregation of the resolution set is a <= constraint */
      val = resolutionset->vals[i];

      /* assert that not both bounds are infinite since if they were no propagation would have happened */
      assert( !(SCIPsetIsInfinity(set, -locallb) && SCIPsetIsInfinity(set, localub)) );

      if( vars[v] == SCIPbdchginfoGetVar(bdchginfo) )
      {
         SCIP_Real solval;

         if( val > 0.0 )
         {
            /* it the case of >= constraints we get tight propagation at the point (coef * ub - slack) / coef */
            solval = MIN( SCIPvarGetUbGlobal(vars[v]), ( val * SCIPvarGetUbAtIndex(vars[v], bdchgidx, TRUE) - resolutionset->slack ) / val );
            /* this can only happen if the reason was a negated clique in the knapsack constraint handler */
            if( SCIPsetIsLT(set, solval, SCIPvarGetLbGlobal(vars[v])) )
            {
               assert(strcmp(SCIPconshdlrGetName(bdchginfo->inferencedata.reason.cons->conshdlr), "knapsack") == 0);
               solval = SCIPvarGetUbAtIndex(vars[v], bdchgidx, TRUE);
            }
         }
         else
         {
            assert(val < 0.0);
            solval = MAX( SCIPvarGetLbGlobal(vars[v]), ( val * SCIPvarGetLbAtIndex(vars[v], bdchgidx, TRUE) - resolutionset->slack ) / val );
            if( SCIPsetIsGT(set, solval, SCIPvarGetUbGlobal(vars[v])) )
            {
               assert(strcmp(SCIPconshdlrGetName(bdchginfo->inferencedata.reason.cons->conshdlr), "knapsack") == 0);
               solval = SCIPvarGetLbAtIndex(vars[v], bdchgidx, TRUE);
            }
         }

         SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[v], solval) );
         SCIPsetDebugMsg(set, " reference point value for resolving variable <%s> in reason: %f \n", SCIPvarGetName(vars[v]), solval);
      }

      else if( val > 0.0 )
      {
         assert(!SCIPsetIsInfinity(set, -locallb));
         SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[v], locallb) );
      }
      else
      {
         assert(val < 0.0);
         assert(!SCIPsetIsInfinity(set, localub));
         SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[v], localub) );
      }
   }
   return SCIP_OKAY;
}

/** calculates efficacy of a given aggregation row w.r.t. a given reference point */
static
SCIP_Real aggrRowGetEfficacy(
   SCIP_AGGRROW*         aggrrow,            /**< aggregation row */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_SOL*             sol                 /**< reference point */
   )
{
   SCIP_VAR** vars;
   SCIP_Real activity;
   SCIP_Real norm;
   int* inds;
   int nnz;
   int i;

   assert(set != NULL);

   vars = SCIPprobGetVars(transprob);
   assert(vars != NULL);

   activity = 0.0;
   nnz = SCIPaggrRowGetNNz(aggrrow);
   inds = SCIPaggrRowGetInds(aggrrow);

   for( i = 0; i < nnz; i++ )
   {
      SCIP_Real val;
      int v = inds[i];

      assert(SCIPvarGetProbindex(vars[v]) == v);

      val = SCIPaggrRowGetProbvarValue(aggrrow, v);
      activity += val * SCIPsolGetVal(sol, set, stat, vars[v]);
   }

   norm = SCIPaggrRowCalcEfficacyNorm(set->scip, aggrrow);
   return (activity - SCIPaggrRowGetRhs(aggrrow)) / MAX(1e-6, norm);
}

/** calculates a c-MIR cut from the coefficients of the resolution set
 */
static
SCIP_RETCODE computecMIRfromResolutionSet(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Real*            cutcoefs,           /**< the coefficients of the MIR cut */
   int*                  cutinds,            /**< the variable indices of the MIR cut */
   SCIP_Real*            cutrhs,             /**< the RHS of the MIR cut */
   int*                  cutnnz,             /**< the number of non-zeros in the cut */
   SCIP_Bool             isconflict,         /**< distinguish between reason and conflict constraint */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change information */
   SCIP_Bool*            success             /**< was the MIR cut successfully computed? */
   )
{
   SCIP_AGGRROW* aggrrow;
   SCIP_SOL* refsol;

   SCIP_Real* rowvals;
   int* rowinds;

   SCIP_Real resolutionefficacy;
   SCIP_Real cutefficacy;
   int nnz;
   int i;

   SCIP_Bool islocal;
   SCIP_Bool cutsuccess;

   *success = FALSE;

   /* creating the aggregation row. There will be only a single row in this aggregation, since it is only used to
    * compute the MIR coefficients
    */
   SCIP_CALL( SCIPaggrRowCreate(set->scip, &aggrrow) );

   /* All values must be negated since the aggregation row requires a RHS, and resolution sets are computed with a LHS */
   SCIP_CALL( SCIPallocBufferArray(set->scip, &rowvals, resolutionset->nnz) );
   SCIP_CALL( SCIPallocBufferArray(set->scip, &rowinds, resolutionset->nnz) );

   assert(!SCIPisInfinity(set->scip, resolutionset->lhs) && !SCIPisInfinity(set->scip, -resolutionset->lhs));

   nnz = 0;
   for( i = 0; i < resolutionset->nnz; i++ )
   {
      if ( !SCIPsetIsZero(set, resolutionset->vals[i]) )
      {
         rowinds[i] = resolutionset->inds[i];
         rowvals[i] = -resolutionset->vals[i];
         nnz++;
      }
   }

   assert(resolutionsetGetNNzs(resolutionset) == nnz);
   if ( nnz > 0 )
   {
      /* create the aggregation row */
      SCIP_CALL( SCIPaggrRowAddCustomCons(set->scip, aggrrow, rowinds, rowvals, nnz, -resolutionset->lhs, 1.0, 1, FALSE) );

#ifdef SCIP_DEBUG
{
      SCIP_VAR** vars;
      vars = SCIPprobGetVars(prob);

      for( i = 0; i < SCIPaggrRowGetNNz(aggrrow); i++ )
      {
         SCIP_Real aggrrow_val;
         aggrrow_val = SCIPaggrRowGetProbvarValue(aggrrow, SCIPvarGetProbindex(vars[resolutionset->inds[i]]));
         assert(SCIPsetIsEQ(set, -resolutionset->vals[i], aggrrow_val));
      }
}
#endif

      /* create reference solution */
      SCIP_CALL( SCIPcreateSol(set->scip, &refsol, NULL) );

      if ( isconflict )
      {
         /* compute the reference point */
         SCIP_CALL( computeReferenceSolutionConflict(resolutionset, set, prob, stat, tree, refsol, bdchginfo, success ) );
      }
      else
      {
         /* compute the cut efficacy */
         SCIP_CALL( computeReferenceSolutionReason(resolutionset, set, prob, stat, tree, refsol, bdchginfo, success ) );
      }

      resolutionefficacy = aggrRowGetEfficacy(aggrrow, prob, set, stat, refsol);
      SCIPdebugMessage("efficacy of resolution set: %f \n", resolutionefficacy);

      cutefficacy = 0.0;

      /* start timing flowcover */
      SCIPclockStart(conflict->resflowcovertime, set);

      /* apply flow cover */
      SCIP_CALL( SCIPcalcFlowCover(set->scip, refsol, POSTPROCESS, BOUNDSWITCH, ALLOWLOCAL, aggrrow, \
            cutcoefs, cutrhs, cutinds, cutnnz, &cutefficacy, NULL, &islocal, &cutsuccess) );

      SCIPclockStop(conflict->resflowcovertime, set);
      *success = cutsuccess;
      conflict->nresflowcovercalls += 1;
      if( cutsuccess )
         conflict->nresflowcover += 1;

      /* start timing MIR */
      SCIPclockStart(conflict->resmirtime, set);

      /* apply MIR */
      if( set->conf_applycmir || !isconflict )
      {
         SCIP_CALL( SCIPcutGenerationHeuristicCMIR(set->scip, refsol, POSTPROCESS, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, TRUE,  \
            INT_MAX, NULL, NULL, MINFRAC, MAXFRAC, aggrrow, cutcoefs, cutrhs, cutinds, cutnnz, &cutefficacy, NULL, \
            &islocal, &cutsuccess) );
      }
      else if( set->conf_applysimplemir )
      {
         SCIP_CALL( SCIPcalcMIR(set->scip, refsol, POSTPROCESS, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, \
            FALSE, NULL, NULL, MINFRAC, MAXFRAC, 1.0, aggrrow, cutcoefs, cutrhs, cutinds, cutnnz, &cutefficacy, NULL, \
            &islocal, &cutsuccess) );
      }

      SCIPclockStop(conflict->resmirtime, set);

      conflict->nresmircalls += 1;
      if( cutsuccess )
         conflict->nresmir += 1;
      *success = (*success || cutsuccess);

      /* try to tighten the coefficients of the cut */
      if( (*success) && !islocal )
      {
         SCIP_Bool redundant;
         int nchgcoefs;

         redundant = SCIPcutsTightenCoefficients(set->scip, FALSE, cutcoefs, cutrhs, cutinds, cutnnz, &nchgcoefs);

         (*success) = !redundant;
      }

      SCIP_CALL( SCIPfreeSol(set->scip, &refsol) );
   }
   /* freeing the local memory */
   SCIPfreeBufferArray(set->scip, &rowinds);
   SCIPfreeBufferArray(set->scip, &rowvals);
   SCIPaggrRowFree(set->scip, &aggrrow);

   return SCIP_OKAY;
}


/** for every variable in the row, except the inferred variable, add bound changes */
static
SCIP_RETCODE addConflictBounds(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_BDCHGIDX*        inferbdchgidx       /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{

   /* scan through the row and add bound changes that make the constraint infeasible */
   /* todo stop adding bounds when infeasibility is detected */
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
            bnd = SCIPgetVarUbAtIndex(set->scip, vars[v], inferbdchgidx, FALSE);
            if ( SCIPsetIsLT(set, bnd, SCIPvarGetUbGlobal(vars[v])) )
            {
               SCIP_CALL( SCIPaddConflictUb(set->scip, vars[v], inferbdchgidx) );
            }
         }
      }
      else
      {
         SCIP_Real bnd;
         if ( SCIPvarGetNBdchgInfosLb(vars[v]) > 0 )
         {
            bnd = SCIPgetVarLbAtIndex(set->scip, vars[v], inferbdchgidx, FALSE);
            if ( SCIPsetIsGT(set, bnd, SCIPvarGetLbGlobal(vars[v])) )
            {
               SCIP_CALL( SCIPaddConflictLb(set->scip, vars[v], inferbdchgidx) );
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
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
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

   vals = resolutionsetGetVals(resolutionset);

   SCIP_CALL( SCIPallocBufferArray(set->scip, &consvars, resolutionsetGetNNzs(resolutionset)) );

   lhs = resolutionsetGetLhs(resolutionset);

   for( i = 0; i < resolutionsetGetNNzs(resolutionset); ++i )
   {
      consvars[i] = vars[resolutionset->inds[i]];
   }

   /* create a constraint out of the conflict set */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "confres_%" SCIP_LONGINT_FORMAT, conflict->nresconfconss);
   SCIP_CALL( SCIPcreateConsLinear(set->scip, &cons, consname, resolutionsetGetNNzs(resolutionset), consvars, vals,
              lhs, SCIPsetInfinity(set), FALSE, set->conf_separesolution, FALSE, FALSE, TRUE, (SCIPnodeGetDepth(tree->path[resolutionset->validdepth]) > 0 ),
              FALSE, set->conf_dynamic, set->conf_removable, FALSE) );

   /* try to automatically convert a linear constraint into a more specific and more specialized constraint */
   SCIP_CALL( SCIPupgradeConsLinear(set->scip, cons, &upgdcons) );
   if( upgdcons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(set->scip, &cons) );
      cons = upgdcons;
   }
   /* chck if the constraint is valid for the dubug solution */
   SCIP_CALL( SCIPdebugCheckConss(set->scip, &cons, 1) );

   /* add conflict to SCIP */
   SCIP_CALL( SCIPaddConflict(set->scip, tree->path[insertdepth], cons, tree->path[resolutionset->validdepth], SCIP_CONFTYPE_RESOLUTION, conflict->resolutionset->usescutoffbound) );
   *success = TRUE;
   /* free temporary memory */
   SCIPfreeBufferArray(set->scip, &consvars);

   return SCIP_OKAY;
}

/** create resolution constraints out of resolution sets */
SCIP_RETCODE SCIPconflictFlushResolutionSets(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set to add to the tree */
   SCIP_Bool*            success             /**< true if the conflict is added to the problem */

   )
{

   int focusdepth;
   int maxsize;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(tree != NULL);

   focusdepth = SCIPtreeGetFocusDepth(tree);
   assert(focusdepth <= SCIPtreeGetCurrentDepth(tree));
   assert(SCIPtreeGetCurrentDepth(tree) == tree->pathlen-1);

   /* todo Maybe remove the maxsize for individual constraints? */
   /* calculate the maximal size of each accepted conflict set */
   maxsize = transprob->nvars / 2;

   SCIPsetDebugMsg(set, "flushing %d resolution sets at focus depth %d (vd: %d, cd: %d, rd: %d, maxsize: %d)\n",
      1, focusdepth, resolutionset->validdepth, resolutionset->conflictdepth, resolutionset->repropdepth, maxsize);

   assert(resolutionset != NULL);
   assert(resolutionset->validdepth == 0);

   *success = FALSE;
   /* do not add long conflicts */
   if( resolutionsetGetNNzs(resolutionset) > maxsize )
   {
      SCIPsetDebugMsg(set, " -> resolution set is too long: %d > %d nnzs\n", resolutionsetGetNNzs(resolutionset), maxsize);
      return SCIP_OKAY;
   }
   /* if the resolution set is empty and the lhs negative, the node and its sub tree in the conflict set's valid depth
    *  can be cut off completely
    */
   else if( resolutionsetGetNNzs(resolutionset) == 0 && SCIPsetIsLT(set, resolutionset->lhs, 0.0))
   {
      SCIPsetDebugMsg(set, " -> empty resolution set with lhs %f in depth %d cuts off sub tree at depth %d\n",
         resolutionset->lhs, focusdepth, resolutionset->validdepth);

      SCIP_CALL( SCIPnodeCutoff(tree->path[resolutionset->validdepth], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
      return SCIP_OKAY;
   }
   /* @todo do not create if relaxation only variable exist */
   else
   {
      /* @todo use the right insert depth and not valid depth */
      SCIP_CALL( createAndAddResolutionCons(conflict, blkmem, set, stat, transprob, origprob, \
                     tree, reopt, lp, cliquetable, resolutionset, resolutionset->validdepth, success) );
      SCIPsetDebugMsg(set, " -> resolution set added (cdpt:%d, fdpt:%d, insert:%d, valid:%d, conf: %d, reprop: %d , len:%d):\n",
                     SCIPtreeGetCurrentDepth(tree), SCIPtreeGetFocusDepth(tree),
                     resolutionset->validdepth, resolutionset->validdepth, resolutionset->conflictdepth,
                     resolutionset->repropdepth, resolutionsetGetNNzs(resolutionset));

   }

   return SCIP_OKAY;
}

/** adds given data as row to the resolution set */
static
SCIP_RETCODE resolutionsetAddSparseData(
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Real*            vals,               /**< variable coefficients */
   int*                  inds,               /**< variable array */
   int                   nnz,                /**< size of variable and coefficient array */
   SCIP_Real             lhs,                /**< left-hand side of resolution set */
   SCIP_Real             origrhs,            /**< right-hand side of the row */
   SCIP_Real             origlhs,            /**< left-hand side of the row */
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

/** compute scale for the reason constraint */
static
SCIP_Real computeScaleReason(
   SCIP_SET*             set,                   /**< global SCIP settings */
   SCIP_RESOLUTIONSET*   conflictresolutionset, /**< conflict resolution set */
   SCIP_RESOLUTIONSET*   reasonresolutionset,   /**< reason resolution set */
   int                   residx,                /**< index of variable to resolve */
   SCIP_Bool             wasresolved            /**< resolution has been applied */
   )
{
   SCIP_Real coefconf;
   SCIP_Real coefreas;
   SCIP_Real scale;

   SCIP_Bool idxinconflict;
   SCIP_Bool idxinreason;

   idxinconflict = FALSE;
   idxinreason = FALSE;

   coefconf = 0.0;
   coefreas = 0.0;

   /* find in the conflict resolution set the coefficient of the variable we are resolving */
   idxinconflict = getCoefInResolutionSet(conflictresolutionset, residx, &coefconf);

   /* find in the reason resolution set the coefficient of the variable we are resolving */
   idxinreason = getCoefInResolutionSet(reasonresolutionset, residx, &coefreas);

   assert(wasresolved || !SCIPsetIsZero(set, coefconf));
   assert(!SCIPsetIsZero(set, coefreas));

   if ( wasresolved )
   {
      if ( SCIPsetIsZero(set, coefconf) )
      {
         SCIPsetDebugMsg(set, "Coefficient of resolvent in conflict is zero");
         return SCIPsetInfinity(set);
      }
      assert((idxinconflict && idxinreason));
      if (SCIPsetIsGE(set, coefconf * coefreas, 0.0))
      {
         SCIPsetDebugMsg(set, "Coefficient of resolvent has the same sign in both conflict and reason");
         return SCIPsetInfinity(set);
      }
   }

   assert(coefconf * coefreas < 0);
   assert((idxinconflict && idxinreason));

   scale = REALABS( coefconf / coefreas );

   SCIP_UNUSED(idxinconflict);
   SCIP_UNUSED(idxinreason);

   return scale;

}

/** reduce a resolution set:
 * - Iteratively weaken variables from the resolution set that do not affect the
 *   slack. Then apply coefficient tightening to reduce the slack.
 * - We weaken a variable if:
 *    * it is free
 *    * it has a positive coefficient and its local upper bound is equal to the
 *      global upper bound
 *    * it has a negative coefficient and its local lower bound is equal to the
 *      global lower bound
 * - For the reason resolution set (isreason = TRUE) we weaken as long as the linear combination of
 *   slack_conflict and slack_reason is positive. i. e. slack_conflict + scale * slack_reason >= 0
*/
static
SCIP_RETCODE ReduceResolutionSet(
   SCIP_CONFLICT *       conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_BDCHGIDX*        currbdchgidx,       /**< index of current bound change */
   int                   residx,             /**< index of variable to resolve */
   int*                  nvarsweakened,      /**< number of weakened variables */
   SCIP_Bool             isreason,           /**< distinguish between reason and conflict constraint */
   SCIP_Real*            fixbounds,          /**< array of fixed bounds */
   int*                  fixinds,            /**< array of indices of fixed variables */
   SCIP_Bool             wasresolved         /**< whether resolution has been applied */
   )
{
   SCIP_VAR** vars;
   SCIP_Real previousslack;
   SCIP_Real scale;
   SCIP_RESOLUTIONSET* resolutionset;
   int i;
   SCIP_Bool applytightening;
   int nchgcoefs;

   assert(conflict != NULL);


   if( isreason )
      resolutionset = conflict->reasonset;
   else
      resolutionset = conflict->resolutionset;

   vars = SCIPprobGetVars(prob);
   i = 0;
   *nvarsweakened = 0;
   applytightening = FALSE;
   previousslack = resolutionset->slack;

   while ( i < resolutionsetGetNNzs(resolutionset) )
   {
      SCIP_Bool varwasweakened;
      SCIP_VAR* vartoweaken;

      varwasweakened = FALSE;
      vartoweaken = vars[resolutionset->inds[i]];

      if ( resolutionset->inds[i] != residx )
      {
         if( resolutionset->vals[i] > 0.0 )
         {
            int nubchgs;
            nubchgs = SCIPvarGetNBdchgInfosUb(vars[resolutionset->inds[i]]);
            if ( nubchgs == 0 )
            {
               weakenVarReason(resolutionset, set, vartoweaken, i);
               varwasweakened = TRUE;
               applytightening = TRUE;
               ++(*nvarsweakened);
            }
         }
         else
         {
            int nlbchgs;
            nlbchgs = SCIPvarGetNBdchgInfosLb(vars[resolutionset->inds[i]]);

            assert( resolutionset->vals[i] < 0.0);
            if ( nlbchgs == 0 )
            {
               weakenVarReason(resolutionset, set, vartoweaken, i);
               varwasweakened = TRUE;
               applytightening = TRUE;
               ++(*nvarsweakened);
            }
         }
      }
      if (varwasweakened)
      {
         if ( isreason && !set->conf_weakenreasonall && (set->conf_batchcoeftight > 0) &&
            (*nvarsweakened % set->conf_batchcoeftight == 0) )
         {
            SCIP_CALL( tightenCoefLhs(set, prob, FALSE, resolutionset->vals, &resolutionset->lhs,
                           resolutionset->inds, &resolutionset->nnz, &nchgcoefs, NULL) );
            applytightening = FALSE;
            if (nchgcoefs > 0)
            {
               previousslack = resolutionset->slack;
               /* todo the update of the slack should be included in tightenCoefLhs */
               resolutionset->slack = getSlack(set, prob, resolutionset, currbdchgidx, fixbounds, fixinds);
               scale = computeScaleReason(set, conflict->resolutionset, conflict->reasonset, residx, wasresolved);
               if (SCIPsetIsLT(set, conflict->resolutionset->slack + scale * resolutionset->slack, 0.0))
               {
                  /** if the linear combination of the slack of the conflict and the slack of the reason is negative
                   * we can stop weakening since the slack is subadditive and the resolvent is guaranteed to have
                   * negative slack
                   */
                  break;
               }
            }
         }
      }
      else
         i++;
   }
   if (applytightening)
   {
      SCIP_CALL( tightenCoefLhs(set, prob, FALSE, resolutionset->vals, &resolutionset->lhs, resolutionset->inds,
                     &resolutionset->nnz, &nchgcoefs, NULL) );
      /* recompute slack */
      resolutionset->slack = getSlack(set, prob, resolutionset, currbdchgidx, fixbounds, fixinds);
      assert(SCIPsetIsRelLE(set, resolutionset->slack, previousslack));

   }

   SCIP_UNUSED(previousslack);

#ifdef SCIP_DEBUG
   {
      if ( isreason )
         SCIPsetDebugMsg(set, " Slack of reason after weakening %d variables is: %f \n", *nvarsweakened, resolutionset->slack);
      else
         SCIPsetDebugMsg(set, " Slack of resolved constraint after weakening %d variables is: %f \n", *nvarsweakened, resolutionset->slack);
   }
#endif
   /* sort for linear time resolution*/
   SCIPsortIntReal(resolutionset->inds, resolutionset->vals, resolutionsetGetNNzs(resolutionset));

   return SCIP_OKAY;
}

/** Apply the resolution step: conflict + scale * reason */
static
SCIP_RETCODE resolveWithReason(
   SCIP_SET*             set,                   /**< global SCIP settings */
   SCIP_RESOLUTIONSET*   conflictresolutionset, /**< conflict resolution set */
   SCIP_RESOLUTIONSET*   reasonresolutionset,   /**< reason resolution set */
   BMS_BLKMEM*           blkmem,                /**< block memory */
   SCIP_PROB*            transprob,             /**< transformed problem */
   int                   residx,                /**< index of variable to resolve */
   SCIP_Bool             wasresolved,           /**< resolution has been applied */
   SCIP_Bool*            success                /**< apply resolution */
   )
{
   int i;
   SCIP_Real scale;
   SCIP_Real largestcoef;
   SCIP_Real smallestcoef;

   int newsize;
   int cidx;
   int previousnnz;
   int newnnz;

   *success = FALSE;

   scale = computeScaleReason(set, conflictresolutionset, reasonresolutionset, residx, wasresolved);
   /* stop if the scale becomes too large */
   if ( SCIPsetIsGE(set, scale,  set->conf_generalresminmaxquot) )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   newsize = resolutionsetGetNNzs(conflictresolutionset) + resolutionsetGetNNzs(reasonresolutionset);
   if ( conflictresolutionset->size < newsize )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictresolutionset->inds, conflictresolutionset->size, newsize ) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictresolutionset->vals, conflictresolutionset->size, newsize ) );
      conflictresolutionset->size = newsize;
   }

   cidx = 0;
   previousnnz = resolutionsetGetNNzs(conflictresolutionset);

   i = 0;
   /* add conflict and reason resolution sets */
   while ( i < resolutionsetGetNNzs(reasonresolutionset) )
   {
      if (cidx >= previousnnz)
      {
         conflictresolutionset->inds[conflictresolutionset->nnz] = reasonresolutionset->inds[i];
         conflictresolutionset->vals[conflictresolutionset->nnz] = scale * reasonresolutionset->vals[i];
         conflictresolutionset->nnz++;
         i++;
      }
      else if (reasonresolutionset->inds[i] == conflictresolutionset->inds[cidx])
      {
         /* @todo quadprecision when adding coefficients? */
         conflictresolutionset->vals[cidx] = conflictresolutionset->vals[cidx] + scale * reasonresolutionset->vals[i];
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
         conflictresolutionset->vals[conflictresolutionset->nnz] = scale * reasonresolutionset->vals[i];
         conflictresolutionset->nnz++;
         i++;
      }
   }
   conflictresolutionset->lhs = conflictresolutionset->lhs + scale * reasonresolutionset->lhs;

   newnnz = resolutionsetGetNNzs(conflictresolutionset);

   largestcoef = -SCIPsetInfinity(set);
   smallestcoef = SCIPsetInfinity(set);
   /* remove coefficients that are almost zero (10^-9 tolerance), loop backwards */
   for( i = newnnz - 1; i >= 0 ; i-- )
   {
      if (SCIPsetIsZero(set, conflictresolutionset->vals[i] ))
      {
         resolutionsetRemoveZeroVar(conflictresolutionset, set, i);
      }
      else
      {
         smallestcoef = MIN(smallestcoef, conflictresolutionset->vals[i]);
         largestcoef = MAX(largestcoef, conflictresolutionset->vals[i]);
      }
   }
   SCIPsetDebugMsg(set, "Nonzeros in resolved constraint: %d \n", resolutionsetGetNNzs(conflictresolutionset));
   SCIPdebug(resolutionsetPrintRow(conflictresolutionset, set, transprob, 3));

   /* check if the quotient of coefficients in the resolvent exceeds the max allowed quotient */
   conflictresolutionset->coefquotient = (conflictresolutionset->nnz > 0) ? fabs(largestcoef / smallestcoef) : 0.0;
   if ( SCIPsetIsGT(set, conflictresolutionset->coefquotient, set->conf_generalresminmaxquot) )
   {
      *success = FALSE;
      SCIPsetDebugMsg(set, "Quotient %f exceeds max allowed quotient", (conflictresolutionset->nnz > 0) ? fabs(largestcoef / smallestcoef) : 0.0);
      return SCIP_OKAY;
   }

   /* sort for linear time resolution */
   SCIPsortIntReal(conflictresolutionset->inds, conflictresolutionset->vals, resolutionsetGetNNzs(conflictresolutionset));
   *success = TRUE;

   return SCIP_OKAY;
}

/** add a row to the resolutionset */
static
SCIP_RETCODE reasonResolutionsetFromRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row,                /**< row to add */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to resolve */
   )
{
   SCIP_COL** cols;
   SCIP_VAR* vartoresolve;
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

   vartoresolve = SCIPbdchginfoGetVar(bdchginfo);
   varidx = SCIPvarGetProbindex(vartoresolve);

   vals = SCIProwGetVals(row);
   cols = SCIProwGetCols(row);

   /* @todo buffer mem (or avoid it by using the cols) */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &inds, nnz) );

   isincon = FALSE;
   changesign = FALSE;
   for( i = 0; i < nnz; i++ )
   {
      SCIP_VAR* var;

      var = SCIPcolGetVar(cols[i]);
      inds[i] = SCIPvarGetProbindex(var);
      if ( inds[i] == varidx )
      {
         assert(var == vartoresolve);
         if ( (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER && vals[i] < 0)
           || (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER && vals[i] > 0) )
         {
            if (strcmp(SCIPconshdlrGetName(bdchginfo->inferencedata.reason.cons->conshdlr), "knapsack") != 0)
               assert(!SCIPsetIsInfinity(set, origrhs));
            changesign = TRUE;
         }
         else if ( (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER && vals[i] > 0)
                || (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER && vals[i] < 0) )
         {
            if (strcmp(SCIPconshdlrGetName(bdchginfo->inferencedata.reason.cons->conshdlr), "knapsack") != 0)
               assert(!SCIPsetIsInfinity(set, -origlhs));
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

   SCIP_CALL( resolutionsetAddSparseData(resolutionset, blkmem, vals, inds, nnz, lhs, origrhs, origlhs, changesign) );

   BMSfreeBlockMemoryArray(blkmem, &inds, nnz);

   return SCIP_OKAY;
}

/** add a row to the resolution set */
static
SCIP_RETCODE conflictResolutionsetFromRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row,                /**< row to add */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to resolve */
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

      /* todo buffer mem (or avoid it by using the cols) */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &inds, nnz) );

   isincon = FALSE;
   changesign = FALSE;
   for( i = 0; i < nnz; i++ )
   {
      var = SCIPcolGetVar(cols[i]);
      inds[i] = SCIPvarGetProbindex(var);
      if ( inds[i] == varidx )
      {
         if ( (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER && vals[i] > 0) ||
              (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER && vals[i] < 0) )
         {
            changesign = TRUE;
         }
         else if ( (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER && vals[i] < 0) ||
                   (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER && vals[i] > 0) )

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

   SCIP_CALL( resolutionsetAddSparseData(resolutionset, blkmem, vals, inds, nnz, lhs, origrhs, origlhs, changesign) );

   BMSfreeBlockMemoryArray(blkmem, &inds, nnz);
   return SCIP_OKAY;
}

#ifdef SCIP_DEBUG
/** prints a resolution set in debug mode */
void resolutionsetPrintRow(
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set to print */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem */
   int                   type                /**< row type (0: initial, 1:conflict, 2:reason, 3:resolved, 4:cmir) */
   )
{
      SCIP_VAR** vars;
      int v;
      int i;

      vars = SCIPprobGetVars(transprob);
      assert(vars != NULL);

      switch( type )
      {
      case 0:
         SCIPsetDebugMsgPrint(set, "Initial row: ");
         break;
      case 1:
         SCIPsetDebugMsgPrint(set, "Conflict row: ");
         break;
      case 2:
         SCIPsetDebugMsgPrint(set, "Reason row: ");
         break;
      case 3:
         SCIPsetDebugMsgPrint(set, "Resolved row: ");
         break;
      case 4:
         SCIPsetDebugMsgPrint(set, "c-MIR row: ");
         break;
      default:
         SCIPsetDebugMsgPrint(set, "Resolution set row: ");
         break;
      }
      for( i = 0; i < resolutionsetGetNNzs(resolutionset); i++ )
      {
         v = resolutionset->inds[i];
         assert(SCIPvarGetProbindex(vars[v]) == v);
         SCIPsetDebugMsgPrint(set, "%f<%s> ", resolutionset->vals[i], SCIPvarGetName(vars[v]));
      }
      SCIPsetDebugMsgPrint(set, ">= %f \n", resolutionset->lhs);
}

/** print a single bound change in debug mode
*/
static
void printSingleBoundChange(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to print */
)
{
      SCIP_VAR* var;
      var = SCIPbdchginfoGetVar(bdchginfo);
      SCIPsetDebugMsg(set, " -> Bound change <%s> %s %.15g [status:%d, type:%d, depth:%d, pos:%d, reason:<%s>, info:%d] \n",
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
}
/** prints all bound changes in the queue in debug mode
 */
static
void printAllBoundChanges(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_VAR* var;
   int i;

   assert(conflict != NULL);

   SCIPsetDebugMsg(set, " -> Bound changes in queue: \n");
   for( i = 0; i < SCIPpqueueNElems(conflict->resbdchgqueue); ++i )
   {
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->resbdchgqueue)[i]);
      var = bdchginfo->var;
      SCIPsetDebugMsg(set, " -> Bound change %d: <%s> %s %.15g [status:%d, type:%d, depth:%d, pos:%d, reason:<%s>, info:%d] \n",
      i, SCIPvarGetName(var),
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
   }
   SCIPsetDebugMsg(set, " -> End of bound changes in queue. \n");
}

/* print the type of the non resolvable reason in debug mode */
static
void printNonResolvableReasonType(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to check */
   )
{
   SCIP_BOUNDCHGTYPE bdchgtype;

   bdchgtype = SCIPbdchginfoGetChgtype(bdchginfo);
   if (bdchgtype == SCIP_BOUNDCHGTYPE_BRANCHING)
   {
      SCIPsetDebugMsg(set, " -> Not resolvable bound change: branching \n");
   }
   else if (bdchgtype == SCIP_BOUNDCHGTYPE_PROPINFER)
   {
      SCIP_PROP* reasonprop;
      reasonprop = SCIPbdchginfoGetInferProp(bdchginfo);

      /* todo check why the propagator can be none */
      SCIPsetDebugMsg(set, " -> Not resolvable bound change: propagation %s \n",
      reasonprop != NULL ? SCIPpropGetName(reasonprop) : "none");
   }
   else
   {
      assert(bdchgtype == SCIP_BOUNDCHGTYPE_CONSINFER);
      SCIPsetDebugMsg(set, " -> Not resolvable bound change: constraint %s \n", SCIPconsGetName(SCIPbdchginfoGetInferCons(bdchginfo)));
   }
}

#endif

/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound() and
 *  SCIPconflictAddRelaxedBound(), and on success, creates and possibly adds a linear constraint
 *  that explains the infeasibility; afterwards the resolution set(s) are cleared
 */
SCIP_RETCODE conflictAnalyzeResolution(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_ROW*             initialconflictrow, /**< row of constraint that detected the conflict */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool             infeasibleLP,       /**< does the conflict originate from an infeasible LP? */
   SCIP_Bool             pseudoobj,          /**< does the conflict originate from a violated pseudo objective bound? */
   int*                  nconss,             /**< pointer to store the number of generated conflict constraints */
   int*                  nconfvars           /**< pointer to store the number of variables in generated conflict constraints */
   )
{
   SCIP_RESOLUTIONSET *conflictresolutionset;
   SCIP_RESOLUTIONSET *reasonresolutionset;
   SCIP_RESOLUTIONSET *prevconflictresolutionset;
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BDCHGINFO* nextbdchginfo;
   SCIP_BDCHGIDX* bdchgidx;
   SCIP_ROW* reasonrow;
   SCIP_CONS* reasoncon;

   int bdchgdepth;
   int focusdepth;
   int currentdepth;
   int maxvaliddepth;
   int maxsize;
   int nchgcoefs;
   int nressteps;
   int nresstepslast;
   int nfuips;
   SCIP_Real* cutcoefs;
   SCIP_Real* fixbounds;
   int* cutinds;
   int* fixinds;
   SCIP_Real conflictslack;
   SCIP_Real reasonslack;
   SCIP_Bool successresolution;
   SCIP_Bool usescutoffbound;
   int i;

   SCIP_VAR* vartoresolve;
   int residx;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(0 <= validdepth && validdepth <= SCIPtreeGetCurrentDepth(tree));
   assert(nconss != NULL);
   assert(nconfvars != NULL);

   usescutoffbound = conflict->resolutionset->usescutoffbound;

   resolutionSetClear(conflict->resolutionset);
   resolutionSetClear(conflict->reasonset);
   resolutionSetClear(conflict->prevresolutionset);
   conflictresolutionset = conflict->resolutionset;
   conflictresolutionset->usescutoffbound = usescutoffbound;
   reasonresolutionset = conflict->reasonset;
   prevconflictresolutionset = conflict->prevresolutionset;

   focusdepth = SCIPtreeGetFocusDepth(tree);
   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);
   assert(focusdepth <= currentdepth);

   *nconss = 0;
   *nconfvars = 0;
   /** check, whether local conflicts are allowed; however, don't generate
    * conflict constraints that are only valid in the probing path and not in
    * the problem tree (i.e. that exceed the focusdepth)
    */
   maxvaliddepth = (set->conf_resallowlocal ? MIN(currentdepth-1, focusdepth) : 0);
   if( validdepth > maxvaliddepth )
      return SCIP_OKAY;

   /* if no row exists conflict analysis is not applicable */
   if (initialconflictrow == NULL)
   {
      SCIPsetDebugMsg(set, "Conflict analysis not applicable since no row is available \n");
      return SCIP_OKAY;
   }

   /* calculate the maximal size of each accepted conflict set */
   maxsize = transprob->nvars / 2;
   if( SCIProwGetNNonz(initialconflictrow) > maxsize )
   {
      SCIPsetDebugMsg(set, "Number of nonzeros in conflict is larger than maxsize %d > %d\n",
                      SCIProwGetNNonz(initialconflictrow), maxsize);
      return SCIP_OKAY;
   }

   SCIPsetDebugMsg(set, "Initial conflict Row: %s \n", SCIProwGetName(initialconflictrow));

   /* last bound change that led to infeasibility */
   bdchginfo = conflictFirstCand(conflict);

   if ( bdchginfo == NULL )
   {
      SCIPsetDebugMsg(set, "Conflict analysis not applicable since resbdchginfo is empty \n");
      return SCIP_OKAY;
   }

#ifdef SCIP_DEBUG
{
   printAllBoundChanges(conflict, set);
}
#endif

   /* if no bound change was infered by a resolvable constraint then we terminate */
   if ( !existsResolvablebdchginfo(conflict) )
   {
      SCIPsetDebugMsg(set, "Conflict analysis not applicable since no resolvable bounds exist \n");
      return SCIP_OKAY;
   }

   /* remove the last bound change */
   bdchginfo = conflictRemoveCand(conflict);
   bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
   bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);

   vartoresolve = bdchginfo->var;
   residx = SCIPvarGetProbindex(vartoresolve);

   /* check if the variable we are resolving is active */
   assert(SCIPvarIsActive(vartoresolve));
#ifdef SCIP_DEBUG
   {
      SCIPsetDebugMsgPrint(set, " First bound change to resolve \n");
      printSingleBoundChange(set, bdchginfo);
   }
#endif

   conflictresolutionset->validdepth = validdepth;

   /* get the resolution set of the conflict row */
   SCIP_CALL( conflictResolutionsetFromRow(set, blkmem, initialconflictrow, conflictresolutionset, bdchginfo) );
   SCIP_CALL( conflictResolutionsetFromRow(set, blkmem, initialconflictrow, prevconflictresolutionset, bdchginfo) );

   SCIPdebug(resolutionsetPrintRow(conflictresolutionset, set, transprob, 0));

   /* weaken or just apply coefficient tightening to the conflict constraint */
   if ( set->conf_weakenconflict )
   {
      int nvarsweakened;
      SCIP_CALL( ReduceResolutionSet(conflict, set, transprob, bdchgidx, residx, &nvarsweakened, FALSE, fixbounds, fixinds, FALSE) );
      SCIPsetDebugMsg(set, "Weakened %d variables, new conflict slack %f \n", nvarsweakened, conflictresolutionset->slack);
      SCIPdebug(resolutionsetPrintRow(conflictresolutionset, set, transprob, 0));
   }
   else
   {
      /* always apply coefficient tightening. It can only reduce the slack even further */
      SCIP_CALL( tightenCoefLhs(set, transprob, FALSE, conflictresolutionset->vals, &conflictresolutionset->lhs,
                     conflictresolutionset->inds, &conflictresolutionset->nnz, &nchgcoefs, NULL) );
      conflictslack = getSlack(set, transprob, conflictresolutionset, bdchgidx, NULL, NULL);
      /* todo assert that the slack is not larger than before? */
      conflictresolutionset->slack = conflictslack;
   }

   /* the slack should be negative */
   if ( SCIPsetIsGE(set, conflictresolutionset->slack, 0.0) )
   {
      /** The only cases where this may not be true is if the conflict is found:
       *  - by a negated clique in the knapsack constraint handler
       *  - by propagating a ranged row
       */
      SCIP_CONSHDLR* conshdlr;
      SCIPsetDebugMsg(set, "Slack of conflict constraint is not negative \n");

      assert(!infeasibleLP && !pseudoobj);

      conshdlr = SCIProwGetOriginConshdlr(initialconflictrow);
      SCIPsetDebugMsg(set, "%s",SCIPconshdlrGetName(conshdlr));
      /* relaxed assertion */
      assert(strcmp(SCIPconshdlrGetName(conshdlr), "knapsack") == 0 || strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0);
      return SCIP_OKAY;
   }

   conflictresolutionset->coefquotient = getQuotLargestSmallestCoef(set, conflictresolutionset->vals, conflictresolutionset->nnz);

   /* initialize indices and bounds for the unresolvable bound changes */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &fixinds, SCIPprobGetNVars(transprob)) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &fixbounds, SCIPprobGetNVars(transprob)) );
   /** set value in fixed indices to 0 to indicate that they are not set
    * if a variable is set at an upper bound, then the value is 1
    * if a variable is set at a lower bound, then the value is -1
    */
   for( i = 0; i < SCIPprobGetNVars(transprob); ++i )
      fixinds[i] = 0;

   nressteps = 0;
   nresstepslast = 0;
   nfuips = 0;

   assert(conflict->nresolutionsets == 0);

   if( set->conf_applycmirreason || set->conf_applycmir || set->conf_applysimplemir)
   {
      SCIP_CALL( SCIPsetAllocBufferArray(set, &cutcoefs, SCIPprobGetNVars(transprob)) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &cutinds, SCIPprobGetNVars(transprob)) );
   }

   SCIPstatisticPrintf("Start Statistics \n");
   SCIPstatisticPrintf("ConflictSTAT: %d %d %f %f\n", nressteps, resolutionsetGetNNzs(conflictresolutionset),
                                                      conflictresolutionset->coefquotient, conflictresolutionset->slack);
   resolutionsetReplace(prevconflictresolutionset, blkmem, conflictresolutionset);

   /** main loop: All-FUIP RESOLUTION
    * --------------------------------
    * - (we already have the initial conflict row and the first bound change to
    *   resolve)
    * - apply coefficient tightening to the conflict row
    * - Ignore & Continue: if we can't explain the bound change, i.e. the reason
    *   is a branching or non-linear then we ignore it and continue with the
    *   next bound change. We have to ignore all other bound changes from
    *   for this variable
    * - if the bound change is resolvable:
    *   * get the reason row for the bound change
    *   * apply coefficient tightening to the reason (maybe also cMIR?)
    *   * take the linear combination of the conflict row and the reason row
    *   * apply coefficient tightening to the resolved row (maybe also cMIR?)
    *       - if there is no other bound change in the queue from the same depth level
    *         then we are at a UIP -> keep this constraint and continue
    */
   while( TRUE )
   {
      /** check if the bound change is resolvable. If it is not, we can ignore the bound change and continue
       * with the next one
       */
      if ( !bdchginfoIsResolvable(bdchginfo) && set->conf_fixandcontinue)
      {
#ifdef SCIP_DEBUG
         {
            printNonResolvableReasonType(set, bdchginfo);
         }
#endif
            if ( SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER )
            {
               reasoncon = SCIPbdchginfoGetInferCons(bdchginfo);
               /* we resolve only with  globally valid constraints */
               if(!SCIPconsIsGlobal(reasoncon))
               {
                  goto TERMINATE;
               }
            }

         if( existsResolvablebdchginfo(conflict) )
         {
            SCIP_BOUNDTYPE boundtype;
            SCIP_BOUNDCHGTYPE bdchgtype;

            /* if a bound for the variable has already been ignored then abort */
            if( fixinds[SCIPvarGetProbindex(SCIPbdchginfoGetVar(bdchginfo))] != 0 )
               goto TERMINATE;

            boundtype = SCIPbdchginfoGetBoundtype(bdchginfo);
            bdchgtype = SCIPbdchginfoGetChgtype(bdchginfo);
            /* ignore the bound change and continue */
            fixinds[SCIPvarGetProbindex(SCIPbdchginfoGetVar(bdchginfo))] = boundtype == SCIP_BOUNDTYPE_UPPER ? 1 : -1;
            fixbounds[SCIPvarGetProbindex(SCIPbdchginfoGetVar(bdchginfo))] = SCIPbdchginfoGetNewbound(bdchginfo);
            SCIPsetDebugMsgPrint(set, "ignoring-fixing variable %s to %f \n", SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
                  SCIPbdchginfoGetNewbound(bdchginfo));


            /* extract latest bound change from queue */
            bdchginfo = conflictRemoveCand(conflict);

            if( bdchginfo == NULL )
               goto TERMINATE;
            /* todo if we still have not resolved and reached a branching decision we can continue for the 2-FUIP */
            if( nressteps == 0 && bdchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
            {
               bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);
            }

            bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
            vartoresolve = SCIPbdchginfoGetVar(bdchginfo);

            SCIPdebug(resolutionsetPrintRow(conflictresolutionset, set, transprob, 1));

            assert(SCIPsetIsLT(set, getSlack(set, transprob, conflictresolutionset, bdchgidx, fixbounds, fixinds), 0.0));
            /* check if the variable we are resolving is active */
            assert(SCIPvarIsActive(vartoresolve));

            /* todo now its easy to continue for All-FUIP */

#ifdef SCIP_DEBUG
         {
            SCIPsetDebugMsgPrint(set, " First bound change in queue (next one to resolve) \n");
            printSingleBoundChange(set, bdchginfo);
         }
#endif
         }
         /* if no bound change was infered by a resolvable constraint then we terminate */
         else
            goto TERMINATE;
      }
      else if( !bdchginfoIsResolvable(bdchginfo) )
      {

         goto TERMINATE;
      }
      else
      {
         residx = SCIPvarGetProbindex(vartoresolve);
         reasoncon = SCIPbdchginfoGetInferCons(bdchginfo);

         if(!SCIPconsIsGlobal(reasoncon))
         {
            goto TERMINATE;
         }
         /* get the corresponding reason row */
         reasonrow = SCIPconsCreateRow(set->scip, reasoncon);

         assert(reasonrow != NULL);

         /* get the resolution set of the conflict row */
         SCIP_CALL( reasonResolutionsetFromRow(set, blkmem, reasonrow, reasonresolutionset, bdchginfo) );

         /* this happens if the reason is a negated clique found in the knapsack constraint handler */
         /* todo handle this case (should be hard)*/
         if (strcmp(SCIPconshdlrGetName(reasoncon->conshdlr), "knapsack") != 0)
         {
            assert(!SCIPsetIsInfinity(set, -reasonresolutionset->lhs) || !SCIPsetIsInfinity(set, reasonresolutionset->lhs));
         }
         else if(SCIPsetIsInfinity(set, -reasonresolutionset->lhs) || SCIPsetIsInfinity(set, reasonresolutionset->lhs))
         {
            SCIPdebug(resolutionsetPrintRow(reasonresolutionset, set, transprob, 2));
            goto TERMINATE;
         }
         reasonslack = getSlack(set, transprob, reasonresolutionset, bdchgidx, fixbounds, fixinds);
         reasonresolutionset->slack = reasonslack;

         /* Apply cmir after each iteration to strengthen the reason constraint */
         if( SCIPsetIsGT(set, reasonslack, 0.0) && set->conf_applycmirreason )
         {
            int cutnnz;
            SCIP_Real cutrhs;
            SCIP_Bool success;

            cutnnz = 0;
            /* todo rethink about the reference point and the scaling in SCIPcalcMIR */
            SCIP_CALL( computecMIRfromResolutionSet(conflict, set, reasonresolutionset, transprob, stat, tree, cutcoefs,
                                                    cutinds, &cutrhs, &cutnnz, FALSE, bdchginfo, &success) ); /*lint !e644*/

            if( success && varIdxInArray(cutinds, cutnnz, residx))
            {
               SCIP_RESOLUTIONSET* cutresolutionset;
               SCIP_CALL( resolutionsetCreate(&cutresolutionset, blkmem) );
               SCIP_CALL( resolutionsetAddSparseData(cutresolutionset, blkmem, cutcoefs, cutinds, cutnnz, -cutrhs,
                                          reasonresolutionset->origrhs, reasonresolutionset->origlhs, TRUE) );

               /* replace the current resolution set by the cut if the slack of the cut is smaller */
               cutresolutionset->slack = getSlack(set, transprob, cutresolutionset, bdchgidx, fixbounds, fixinds);
               if ( SCIPisLT(set->scip, cutresolutionset->slack, reasonslack) )
               {

                  SCIPdebug(resolutionsetPrintRow(reasonresolutionset, set, transprob, 2));
                  SCIPsetDebugMsg(set, "replacing reason resolution set by cMIR cut resolution set: new slack %f < old slack %f \n",
                                         cutresolutionset->slack, reasonslack);
                  SCIP_CALL( resolutionsetReplace(reasonresolutionset, blkmem, cutresolutionset) );
                  reasonresolutionset->slack = reasonslack;
               }
               SCIPsortIntReal(reasonresolutionset->inds, reasonresolutionset->vals, resolutionsetGetNNzs(reasonresolutionset));
               SCIPresolutionsetFree(&cutresolutionset, blkmem);
            }
         }

         /* sort for linear time resolution */
         SCIPsortIntReal(conflictresolutionset->inds, conflictresolutionset->vals, resolutionsetGetNNzs(conflictresolutionset));
         SCIPsortIntReal(reasonresolutionset->inds, reasonresolutionset->vals, resolutionsetGetNNzs(reasonresolutionset));

         SCIPsetDebugMsg(set, "conflict resolution set: nnzs: %d, slack: %f \n", resolutionsetGetNNzs(conflictresolutionset), conflictslack);
         SCIPsetDebugMsg(set, "reason resolution set: nnzs: %d, slack: %f \n", resolutionsetGetNNzs(reasonresolutionset), reasonslack);
         SCIPdebug(resolutionsetPrintRow(reasonresolutionset, set, transprob, 2));


         if ( set->conf_weakenreason && reasonresolutionset->slack > 0.0 )
         {
            SCIP_Real scale;
            scale = computeScaleReason(set, conflictresolutionset, reasonresolutionset, residx, (conflict->nresolutionsets > 0));
            if ( conflictresolutionset->slack + scale * reasonresolutionset->slack > 0.0 )
            {
               int nvarsweakened;
               SCIP_CALL( ReduceResolutionSet(conflict, set, transprob, bdchgidx, residx, &nvarsweakened, TRUE, fixbounds, fixinds, (nressteps > 0)) );
               SCIPsetDebugMsg(set, "Weakened %d variables, new reason slack %f \n", nvarsweakened, reasonresolutionset->slack);
               SCIPdebug(resolutionsetPrintRow(reasonresolutionset, set, transprob, 2));
            }
         }

         /** Unfortunately we cannot guarrante that the slack becomes zero after reducing the reason
          *  TIll now there are two major problems:
          *    - Knapsack constraints that use negated cliques in the propagation
          *    - Ranged row propagation (gcd argument)
          */
#ifdef SCIP_DEBUG
      {
         printSingleBoundChange(set, bdchginfo);
      }
#endif
         SCIPsetDebugMsg(set, " Resolving conflict and reason on variable <%s>\n", SCIPvarGetName(vartoresolve));

         /* resolution step */
         SCIP_CALL( resolveWithReason(set, conflictresolutionset, reasonresolutionset, blkmem, transprob,
                                      residx, (nressteps > 0), &successresolution) );
         SCIPstatisticPrintf("ConflictSTAT: %d %d %f %f\n", nressteps, resolutionsetGetNNzs(conflictresolutionset),
                             conflictresolutionset->coefquotient, conflictresolutionset->slack);

         /* if the resolution failed then we can still return the latest resolved constraint and terminate*/
         if (!successresolution)
         {
            if ( nressteps >= 1 && nresstepslast != nressteps )
            {
               /* add the previous conflict in the list of resolution sets */
               SCIP_RESOLUTIONSET* tmpconflictresolutionset;
               SCIP_CALL( resolutionsetCopy(&tmpconflictresolutionset, blkmem, prevconflictresolutionset) );
               SCIP_CALL( conflictInsertResolutionset(conflict, set, &tmpconflictresolutionset) );
            }
            goto TERMINATE;
         }
         conflictslack = getSlack(set, transprob, conflictresolutionset, bdchgidx, fixbounds, fixinds);

         conflictresolutionset->slack = conflictslack;

         SCIPsetDebugMsg(set, "Slack of resolved row: %f \n", conflictslack);

         if (SCIPsetIsGE(set, conflictslack, 0.0))
         {
            if ( nressteps >= 1 && nresstepslast != nressteps )
            {
               /* add the previous conflict in the list of resolution sets */
               SCIP_RESOLUTIONSET* tmpconflictresolutionset;
               SCIP_CALL( resolutionsetCopy(&tmpconflictresolutionset, blkmem, prevconflictresolutionset) );
               SCIP_CALL( conflictInsertResolutionset(conflict, set, &tmpconflictresolutionset) );
            }
            goto TERMINATE;
         }

         nressteps++;


         SCIP_CALL( tightenCoefLhs(set, transprob, FALSE, conflictresolutionset->vals, &conflictresolutionset->lhs, conflictresolutionset->inds, &conflictresolutionset->nnz, &nchgcoefs, NULL) );
         if (nchgcoefs > 0)
         {
            SCIP_Real previousslack;
            previousslack = conflictresolutionset->slack;
            conflictresolutionset->slack = getSlack(set, transprob, conflictresolutionset, bdchgidx, fixbounds, fixinds);
            SCIPsetDebugMsg(set, "Tightened %d coefficients in the resolved constraint, old slack %f, new slack %f \n", nchgcoefs, previousslack,conflictresolutionset->slack);
            assert(SCIPsetIsRelLE(set, conflictresolutionset->slack, previousslack));
         }
         /* sort for linear time resolution */
         SCIPsortIntReal(conflictresolutionset->inds, conflictresolutionset->vals, resolutionsetGetNNzs(conflictresolutionset));
         if (SCIPsetIsLT(set, conflictresolutionset->slack, 0.0))
         {
            /* Apply cmir after each iteration to strengthen the conflict constraint */
            if( set->conf_applycmir || set->conf_applysimplemir)
            {
               int cutnnz;
               SCIP_Real cutrhs;
               SCIP_Bool success;

               assert(cutinds != NULL);
               assert(cutcoefs != NULL);

               cutnnz = 0;
               SCIP_CALL( computecMIRfromResolutionSet(conflict, set, conflictresolutionset, transprob, stat, tree,
                                      cutcoefs, cutinds, &cutrhs, &cutnnz, TRUE, bdchginfo, &success) ); /*lint !e644*/

               if( success )
               {
                  SCIP_RESOLUTIONSET* cutresolutionset;
                  SCIP_CALL( resolutionsetCreate(&cutresolutionset, blkmem) );
                  SCIP_CALL( resolutionsetAddSparseData(cutresolutionset, blkmem, cutcoefs, cutinds, cutnnz, -cutrhs,
                                             conflictresolutionset->origrhs, conflictresolutionset->origlhs, TRUE) );
                  /* replace the current resolution set by the cut if the slack of the cut is negative */
                  cutresolutionset->slack = getSlack(set, transprob, cutresolutionset, bdchgidx, fixbounds, fixinds);
                  if ( SCIPisLT(set->scip, cutresolutionset->slack, 0.0) )
                  {
                     SCIP_CALL( resolutionsetReplace(conflictresolutionset, blkmem, cutresolutionset) );
                  }
                  SCIPsortIntReal(conflictresolutionset->inds, conflictresolutionset->vals, resolutionsetGetNNzs(conflictresolutionset));
                  SCIPresolutionsetFree(&cutresolutionset, blkmem);
               }
            }
            SCIP_CALL( resolutionsetReplace(prevconflictresolutionset, blkmem, conflictresolutionset) );
         }

         /* terminate after at most nressteps resolution iterations */
         if( set->conf_maxnumressteps > 0 && nressteps >= set->conf_maxnumressteps )
         {
            if( nresstepslast != nressteps )
            {
               /* add the previous conflict in the list of resolution sets */
               SCIP_RESOLUTIONSET* tmpconflictresolutionset;
               SCIP_CALL( resolutionsetCopy(&tmpconflictresolutionset, blkmem, prevconflictresolutionset) );
               SCIP_CALL( conflictInsertResolutionset(conflict, set, &tmpconflictresolutionset) );
            }
            goto TERMINATE;
         }

         /** clean up the queue of bound changes, e.g.
          * if variables get canceled during resolution
          */
         conflictCleanUpbdchgqueue(conflict, set, conflictresolutionset);

         SCIP_CALL( addConflictBounds(set, transprob, conflictresolutionset, bdchgidx) );

         /* get the next bound change */
         bdchginfo = conflictFirstCand(conflict);

         /* if no bound change exists the we should stop */
         /* in case we have applies resolution steps we keep the last conflict constraint */
         if( bdchginfo == NULL )
         {
            if ( nressteps >= 1 && nresstepslast != nressteps )
            {
               /* add the previous conflict in the list of resolution sets */
               SCIP_RESOLUTIONSET* tmpconflictresolutionset;
               SCIP_CALL( resolutionsetCopy(&tmpconflictresolutionset, blkmem, prevconflictresolutionset) );
               SCIP_CALL( conflictInsertResolutionset(conflict, set, &tmpconflictresolutionset) );
            }
            goto TERMINATE;
         }
         bdchginfo = conflictRemoveCand(conflict);
         bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
         vartoresolve = SCIPbdchginfoGetVar(bdchginfo);
         bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);

         residx = SCIPvarGetProbindex(vartoresolve);

         assert(varIdxInResolutionset(conflictresolutionset, residx));
         /* get the bound change before bdchginfo */
         nextbdchginfo = conflictFirstCand(conflict);

         /* check if the variable we are resolving is active */
         assert(SCIPvarIsActive(vartoresolve));

         /* when at an UIP add the previous conflict in the list of resolution sets */
         if( nextbdchginfo == NULL || SCIPbdchginfoGetDepth(nextbdchginfo) != bdchgdepth  )
         {
            if (nresstepslast != nressteps)
            {
               SCIPsetDebugMsgPrint(set, " reached UIP in depth %d \n", bdchgdepth);
               /* add the previous conflict in the list of resolution sets */
               prevconflictresolutionset->conflictdepth = bdchgdepth;
               prevconflictresolutionset->repropdepth = (nextbdchginfo == NULL) ? 0 : SCIPbdchginfoGetDepth(nextbdchginfo);
               SCIP_RESOLUTIONSET* tmpconflictresolutionset;
               SCIP_CALL( resolutionsetCopy(&tmpconflictresolutionset, blkmem, prevconflictresolutionset) );
               SCIP_CALL( conflictInsertResolutionset(conflict, set, &tmpconflictresolutionset) );
               nresstepslast = nressteps;
               nfuips ++;
            }
         }

#ifdef SCIP_DEBUG
      {
         printAllBoundChanges(conflict, set);
         SCIPsetDebugMsgPrint(set, " First bound change in queue (next one to resolve) \n");
         printSingleBoundChange(set, bdchginfo);

         if( nextbdchginfo != NULL )
         {
            SCIPsetDebugMsgPrint(set, " Next bound change in queue \n");
            printSingleBoundChange(set, nextbdchginfo);
         }
      }
#endif
      }

   }

  TERMINATE:
   SCIPstatisticPrintf("End Statistics \n");
   SCIPsetDebugMsg(set, "Total number of resolution sets found %d\n", conflict->nresolutionsets);
   if ( conflict->nresolutionsets > 0 )
   {
      int nconstoadd;

      /* add the conflict constraints to the current node
      (till now we generate only global constraints, i.e. add constraints to the root node)*/
      nconstoadd = (set->conf_resolutioncons > 0) ? MIN(set->conf_resolutioncons, conflict->nresolutionsets) : conflict->nresolutionsets;

      for( i = 0; i < nconstoadd; i++ )
      {
         SCIP_RESOLUTIONSET* resolutionset;
         resolutionset = conflict->resolutionsets[i];
         assert(SCIPsetIsLT(set, resolutionset->slack, 0.0));

         if ( SCIPsetIsLT(set, conflict->resolutionsets[i]->coefquotient, set->conf_generalresminmaxquot) )
         {
            SCIP_Bool success;
            SCIP_CALL( SCIPconflictFlushResolutionSets(conflict, blkmem, set, stat, transprob, origprob, tree, reopt,
                                            lp, cliquetable, resolutionset, &success) );
            if( success )
            {
               (*nconss)++;
               (*nconfvars) = resolutionsetGetNNzs(resolutionset);
            }
         }
      }
   }

   if( set->conf_applycmirreason || set->conf_applycmir || set->conf_applysimplemir)
   {
      SCIPsetFreeBufferArray(set, &cutinds);
      SCIPsetFreeBufferArray(set, &cutcoefs);
   }


   SCIPsetFreeBufferArray(set, &fixinds);
   SCIPsetFreeBufferArray(set, &fixbounds);

   /* free all resolution sets */
   for( i = 0; i < conflict->nresolutionsets; i++ )
      SCIPresolutionsetFree(&conflict->resolutionsets[i], blkmem);
   conflict->nresolutionsets = 0;

   return SCIP_OKAY;
}

/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound(), and on success,
 * creates a linear constraint that explains the infeasibility
 */
SCIP_RETCODE SCIPconflictAnalyzeResolution(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_ROW*             initialconflictrow, /**< row of constraint that detected the conflict */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool             infeasibleLP,       /**< does the conflict originate from an infeasible LP? */
   SCIP_Bool             pseudoobj,          /**< does the conflict originate from a violated pseudo objective bound? */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{

   SCIP_VAR** vars;
   int nconss;
   int nvars;
   int nconfvars;
   int i;

   /* arrays to store variable information related to conflict analysis */
   int tmp_count;
   SCIP_Real* tmp_conflictlb;
   SCIP_Real* tmp_conflictub;
   SCIP_Real* tmp_conflictrelaxedlb;
   SCIP_Real* tmp_conflictrelaxedub;
   int* tmp_conflictlbcount;
   int* tmp_conflictubcount;

   /* check if generalized resolution conflict analysis is applicable */
   if( !SCIPconflictResolutionApplicable(set) )
      return SCIP_OKAY;

   vars = SCIPprobGetVars(transprob);
   nvars = SCIPprobGetNVars(transprob);

   SCIPallocBufferArray(set->scip, &tmp_conflictlb, nvars); /*lint !e522*/
   SCIPallocBufferArray(set->scip, &tmp_conflictub, nvars); /*lint !e522*/
   SCIPallocBufferArray(set->scip, &tmp_conflictrelaxedlb, nvars); /*lint !e522*/
   SCIPallocBufferArray(set->scip, &tmp_conflictrelaxedub, nvars); /*lint !e522*/
   SCIPallocBufferArray(set->scip, &tmp_conflictlbcount, nvars); /*lint !e522*/
   SCIPallocBufferArray(set->scip, &tmp_conflictubcount, nvars); /*lint !e522*/

   assert(conflict != NULL);
   assert(set != NULL);
   assert(origprob != NULL);
   assert(transprob != NULL);

   if( success != NULL )
      *success = FALSE;

   SCIPsetDebugMsg(set, "resolution based conflict analysis after infeasible propagation in depth %d\n", SCIPtreeGetCurrentDepth(tree));

   tmp_count = conflict->count;

   for (i = 0; i < nvars; i++)
   {
      tmp_conflictlb[i] = vars[i]->conflictlb;
      tmp_conflictub[i] = vars[i]->conflictub;
      tmp_conflictrelaxedlb[i] = vars[i]->conflictrelaxedlb;
      tmp_conflictrelaxedub[i] = vars[i]->conflictrelaxedub;
      tmp_conflictlbcount[i] = vars[i]->conflictlbcount;
      tmp_conflictubcount[i] = vars[i]->conflictubcount;
   }
   /* start timing */
   SCIPclockStart(conflict->resanalyzetime, set);

   conflict->nrescalls++;

   /* setting this to true adds bound changes only to the resolution bdchg queue */
   conflict->bdchgonlyresqueue = TRUE;

   /* analyze the conflict set, and create a conflict constraint on success */
   SCIP_CALL( conflictAnalyzeResolution(conflict, blkmem, set, stat, transprob, origprob, tree, reopt, lp, \
          cliquetable, initialconflictrow, validdepth, infeasibleLP, pseudoobj, &nconss, &nconfvars) );

   conflict->nressuccess += (nconss > 0 ? 1 : 0);
   conflict->nresconfconss += nconss;
   conflict->nresconfvariables += nconfvars;
   if( success != NULL )
      *success = (nconss > 0);

   /* Set variable information related to conflict analysis to the values before using generalized resolution */
   for (i = 0; i < nvars; i++)
   {
      vars[i]->conflictlb = tmp_conflictlb[i];
      vars[i]->conflictub = tmp_conflictub[i];
      vars[i]->conflictrelaxedlb = tmp_conflictrelaxedlb[i];
      vars[i]->conflictrelaxedub = tmp_conflictrelaxedub[i];
      vars[i]->conflictlbcount = tmp_conflictlbcount[i];
      vars[i]->conflictubcount = tmp_conflictubcount[i];
   }
   conflict->count = tmp_count;

   /* free data */
   SCIPfreeBufferArray(set->scip, &tmp_conflictlb);
   SCIPfreeBufferArray(set->scip, &tmp_conflictub);
   SCIPfreeBufferArray(set->scip, &tmp_conflictrelaxedlb);
   SCIPfreeBufferArray(set->scip, &tmp_conflictrelaxedub);
   SCIPfreeBufferArray(set->scip, &tmp_conflictlbcount);
   SCIPfreeBufferArray(set->scip, &tmp_conflictubcount);

   /* free all resolutionsets */
   for( i = 0; i < conflict->nresolutionsets; i++ )
      SCIPresolutionsetFree(&conflict->resolutionsets[i], blkmem);

   conflict->nresolutionsets = 0;
   conflict->resolutionminslack = 0.0;
   conflict->bdchgonlyresqueue = FALSE;

   /* clear the bound change queues */
   SCIPpqueueClear(conflict->resbdchgqueue);
   SCIPpqueueClear(conflict->resforcedbdchgqueue);

   /* stop timing */
   SCIPclockStop(conflict->resanalyzetime, set);
   SCIPsetDebugMsg(set, "resolution based conflict analysis added %d constraints \n \n", nconss);

   return SCIP_OKAY;
}
