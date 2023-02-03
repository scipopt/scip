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

/* todo think about these parameters when using MIR */
#define BOUNDSWITCH                0.51 /**< threshold for bound switching - see cuts.c */
#define POSTPROCESS               FALSE /**< apply postprocessing to the cut - see cuts.c */
#define USEVBDS                   FALSE /**< use variable bounds - see cuts.c */
#define ALLOWLOCAL                FALSE /**< allow to generate local cuts - see cuts. */
#define MINFRAC                   0.05  /**< minimal fractionality of floor(rhs) - see cuts.c */
#define MAXFRAC                   0.999 /**< maximal fractionality of floor(rhs) - see cuts.c */

#define EPS                       1e-06

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
   (*targetresolutionset)->usescutoffbound = sourceresolutionset->usescutoffbound;
   (*targetresolutionset)->isbinary = sourceresolutionset->isbinary;

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
   targetresolutionset->usescutoffbound = sourceresolutionset->usescutoffbound;
   targetresolutionset->isbinary = sourceresolutionset->isbinary;

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

   /*@todo compute activity in a seperate method */
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

   /* propagating constraint cannot be redundant */
   assert(!SCIPisGE(set->scip, minact, *rowlhs));

   /* terminate, because coefficient tightening cannot be performed; also
   excludes the case in which no integral variable is present */
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
         SCIP_Real newcoef = minact - (*rowlhs);
         SCIP_Real ub = localbounds ? SCIPvarGetUbLocal(vars[rowinds[i]]) : SCIPvarGetUbGlobal(vars[rowinds[i]]);

         assert(SCIPisGE(set->scip, newcoef + EPS, rowcoefs[i]));
         assert(!SCIPisPositive(set->scip, newcoef));

         if( newcoef > rowcoefs[i] )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(newlhs);

            SCIPquadprecSumDD(delta, newcoef, -rowcoefs[i]);
            SCIPquadprecProdQD(delta, delta, ub);

            SCIPquadprecSumQD(newlhs, delta, *rowlhs);
            SCIPdebugPrintf("tightened coefficient from %g to %g; lhs changed from %g to %g; the bounds are [%g,%g]\n",
               rowcoefs[i], newcoef, (*rowlhs), QUAD_TO_DBL(newlhs),
               localbounds ? SCIPvarGetUbLocal(vars[rowinds[i]]) : SCIPvarGetUbGlobal(vars[rowinds[i]]), ub);

            *rowlhs = QUAD_TO_DBL(newlhs);

            ++(*nchgcoefs);

            /* todo check if newcoef = 0 is possible */
            /* This should not be possible, because then the constraint would be redundant */
            if( SCIPisNegative(set->scip, newcoef) )
            {
               SCIPquadprecSumQQ(minacttmp, minacttmp, delta);
               minact = QUAD_TO_DBL(minacttmp);
               rowcoefs[i] = newcoef;
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
         SCIP_Real newcoef = (*rowlhs) - minact;
         SCIP_Real lb = localbounds ? SCIPvarGetLbLocal(vars[rowinds[i]]) : SCIPvarGetLbGlobal(vars[rowinds[i]]);

         assert(SCIPisLE(set->scip, newcoef, rowcoefs[i] + EPS));
         assert(!SCIPisNegative(set->scip, newcoef));

         if( newcoef < rowcoefs[i] )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(newlhs);

            SCIPquadprecSumDD(delta, newcoef, -rowcoefs[i]);
            SCIPquadprecProdQD(delta, delta, lb);

            SCIPquadprecSumQD(newlhs, delta, *rowlhs);
            SCIPdebugPrintf("tightened coefficient from %g to %g; lhs changed from %g to %g; the bounds are [%g,%g]\n",
               rowcoefs[i], newcoef, (*rowlhs), QUAD_TO_DBL(newlhs), lb,
               localbounds ? SCIPvarGetUbLocal(vars[rowinds[i]]) : SCIPvarGetUbGlobal(vars[rowinds[i]]));

            *rowlhs = QUAD_TO_DBL(newlhs);

            ++(*nchgcoefs);

            /* todo check if newcoef = 0 is possible */
            /* This should not be possible, because then the constraint would be redundant */
            if( SCIPisPositive(set->scip, newcoef) )
            {
               SCIPquadprecSumQQ(minacttmp, minacttmp, delta);
               minact = QUAD_TO_DBL(minacttmp);
               rowcoefs[i] = newcoef;
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

   /* todo If a reason constraint propagates on a variable x and another
      variable y is free then tightening the coefficient of y is not possible */
  TERMINATE:
   SCIPfreeBufferArrayNull(set->scip, &absvals);

   return SCIP_OKAY;
}

/* check if the resolution set has a relaxation only variable */
static
SCIP_Bool hasRelaxationOnlyVar(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_RESOLUTIONSET*   resolutionset       /**< resolution set */
   )
{
   SCIP_VAR** vars;
   int i;

   assert(set != NULL);
   assert(prob != NULL);
   assert(resolutionset != NULL);

   vars = SCIPprobGetVars(prob);
   /* loop over resolution set */
   for( i = 0; i < resolutionset->nnz; ++i )
   {
      SCIP_VAR* var;

      var = vars[resolutionset->inds[i]];
      assert(var != NULL);

      if( SCIPvarIsRelaxationOnly(var) )
         return TRUE;
   }
   return FALSE;
}

/* check if a resolution set has only binary variables */
static
SCIP_Bool isBinaryResolutionSet(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_RESOLUTIONSET*   resolutionset       /**< resolution set */
   )
{
   SCIP_VAR** vars;
   int i;

   assert(set != NULL);
   assert(prob != NULL);
   assert(resolutionset != NULL);

   vars = SCIPprobGetVars(prob);
   /* loop over resolution set */
   for( i = 0; i < resolutionset->nnz; ++i )
   {
      SCIP_VAR* var;

      var = vars[resolutionset->inds[i]];
      assert(var != NULL);

      if( !SCIPvarIsBinary(var) )
         return FALSE;
   }
   return TRUE;
}

/** perform Normalized Chvatal-Gomory with a divisor d on a constraint: \sum_i \alpha_i x_i \geq b, \] where $\alpha_i , b \in R$
 *   where $x_k$ is a binary variable with positive coefficient $d = \alpha_k > 0$
 *   The CG Cut for $C$ with divisor $d$ is given by
 *   \[\sum_{i \in a_i \ge 0} ceil(\frac{\alpha_i}{d}} x_i + \sum_{j \in a_j < 0} floor(\frac{\alpha_j}{d}) x_j
 *     \geq ceil(\frac{b}{d} - \sum_{j \in a_j < 0}\frac{\alpha_j}{d}) +   \sum_{j \in a_j < 0} floor(\frac{\alpha_j}{d}) ,\]
 */
static
SCIP_RETCODE ChvatalGomoryLhs(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_Real             divisor             /**< the divisor of the row */
   )
{
   /* todo check instance 10teams where a local bound does not seem to have a bdchginfo */
   SCIP_Real negcoefsum;
   int negcoeffloorsum;

   assert(set != NULL);
   assert(prob != NULL);
   assert(resolutionset != NULL);
   assert(resolutionset->inds != NULL);
   assert(resolutionset->vals != NULL);
   assert(resolutionset->nnz > 0);

   assert(SCIPsetIsGT(set, divisor, 0.0));

   SCIPsetDebugMsg(set, "Normalized Chvatal-Gomory on constraint with LHS %f and divisor %f\n" , resolutionset->lhs, divisor);

   /* todo extend Chvatal-Gomory for constraints with general integer variables */
   assert(isBinaryResolutionSet(set, prob, resolutionset));

   negcoefsum = 0.0;
   negcoeffloorsum = 0;
   /* loop over all non zeros and divide each row element by the divisor */
   for( int i = 0; i < resolutionset->nnz; ++i )
   {
      SCIP_Real newcoef;

      newcoef = resolutionset->vals[i] / divisor;
      /* ceil the coefficient if it is positive and floor it if it is negative */
      if ( SCIPsetIsGE(set, newcoef, 0.0) )
      {
         resolutionset->vals[i] = SCIPsetCeil(set, newcoef);
      }
      else
      {
         negcoefsum += newcoef;
         negcoeffloorsum += SCIPsetFloor(set, newcoef);
         resolutionset->vals[i] = SCIPsetFloor(set, newcoef);
      }
   }
   /** new lhs: $ceil(\frac{b}{d} - \sum_{j \in a_j < 0}\frac{\alpha_j}{d}) +
    *  \sum_{j \in a_j < 0} floor(\frac{\alpha_j}{d}) $
    */
   resolutionset->lhs = SCIPsetCeil(set, resolutionset->lhs / divisor - negcoefsum) + negcoeffloorsum;

   /* remove variables with zero coefficient. Loop backwards */
   for( int i = resolutionset->nnz - 1; i >= 0; --i )
   {
      if( SCIPsetIsZero(set, resolutionset->vals[i]) )
      {
         --resolutionset->nnz;
         resolutionset->inds[i] = resolutionset->inds[resolutionset->nnz];
         resolutionset->vals[i] = resolutionset->vals[resolutionset->nnz];
      }
   }
   return SCIP_OKAY;
}

/** perform MIR on a constraint: \sum_i \alpha_i x_i \geq b, \] where $\alpha_i , b \in R$
 *   and $x_k$ is a binary variable with positive coefficient $d = \alpha_k > 0$
 *   The MIR Cut of $C$ with divisor $d$ is the constraint
 *   \[\sum_{i \in I_1} ceil(\frac{\alpha_i}{d}} x_i + \sum_{j \in I_2}\left(floor(\frac{\alpha_j}{d})+\frac{f(\alpha_j/d)}{f(b/d)}\right) x_j
 *     \geq ceil(\frac{b}{d}),\]
 *   where \[I_1 = \{i\, : \, f(a_i/d)\geq f(b/d) \text{ or } f(a_i/d) \in \Z\}, \] \[I_2 = \{j\, : \, f(a_i/d) < f(b/d) \text{ and } f(a_i/d) \notin \Z\},\]
 *   and $f(\cdot) = \cdot - \floor{\cdot}$.
 */
static
SCIP_RETCODE MirLhs(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_RESOLUTIONSET*   resolutionset,      /**< resolution set */
   SCIP_Real             divisor             /**< the divisor of the row */
   )
{
   SCIP_Real fraclhs;

   assert(set != NULL);
   assert(prob != NULL);
   assert(resolutionset != NULL);
   assert(resolutionset->inds != NULL);
   assert(resolutionset->vals != NULL);
   assert(resolutionset->nnz > 0);

   assert(SCIPsetIsGT(set, divisor, 0.0));

   SCIPsetDebugMsg(set, "MIR on constraint with LHS %f and divisor %f\n" , resolutionset->lhs, divisor);

   /* divide the lhs by the divisor and subtract the sum of fractionalities */
   resolutionset->lhs = resolutionset->lhs / divisor;
   fraclhs = resolutionset->lhs - SCIPsetFloor(set, resolutionset->lhs);
   resolutionset->lhs = SCIPsetCeil(set, resolutionset->lhs);
   SCIPsetDebugMsg(set, "New lhs %f with fractionality %f\n" , resolutionset->lhs, fraclhs);

   /* todo extend MIR for continuous and general integer variables */
   assert(isBinaryResolutionSet(set, prob, resolutionset));

   /* loop over all non zeros and divide each row element by the divisor */
   for( int i = 0; i < resolutionset->nnz; ++i )
   {

      SCIP_Real newcoef;
      SCIP_Real frac;

      newcoef = resolutionset->vals[i] / divisor;
      frac = newcoef - SCIPsetFloor(set, newcoef);
      /* ceil the coefficient if it is positive and floor plus quotient of
      fractionalities  if it is negative */
      if ( SCIPsetIsGE(set, frac, fraclhs) )
      {
         resolutionset->vals[i] = SCIPsetCeil(set, newcoef);
      }
      else
      {
         resolutionset->vals[i] = SCIPsetFloor(set, newcoef) + frac / fraclhs;
      }
   }

   /* remove variables with zero coefficient. Loop backwards */
   for( int i = resolutionset->nnz - 1; i >= 0; --i )
   {
      if( SCIPsetIsZero(set, resolutionset->vals[i]) )
      {
         --resolutionset->nnz;
         resolutionset->inds[i] = resolutionset->inds[resolutionset->nnz];
         resolutionset->vals[i] = resolutionset->vals[resolutionset->nnz];
      }
   }
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
      /* todo also other propagators can be resolved */
      if (SCIPbdchginfoGetInferProp(bdchginfo) == NULL)
         return FALSE;
      else if( strcmp(SCIPpropGetName(SCIPbdchginfoGetInferProp(bdchginfo)), "pseudoobj") == 0 )
         return TRUE;
      else
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
   else if( strcmp(conshdlrname, "orbisack") == 0 )
   {
      return TRUE;
   }
   else if( strcmp(conshdlrname, "orbitope") == 0 )
   {
      return TRUE;
   }
   else if( strcmp(conshdlrname, "and") == 0 )
   {
      return TRUE;
   }
   if( strcmp(conshdlrname, "xor") == 0 )
   {
      return TRUE;
   }
   if( strcmp(conshdlrname, "or") == 0 )
   {
      return TRUE;
   }
   if( strcmp(conshdlrname, "bounddisjunction") == 0 )
   {
      return TRUE;
   }
   return FALSE;
}

/** returns if the we can extract the reason bound changed by reverse propagation
 *  todo at the moment not all constraints/propagators are supported
 */
static
SCIP_Bool reasonIsLinearizable(
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
      /* todo also other propagators can be resolved */
      if (SCIPbdchginfoGetInferProp(bdchginfo) == NULL)
         return FALSE;
      else if( strcmp(SCIPpropGetName(SCIPbdchginfoGetInferProp(bdchginfo)), "pseudoobj") == 0 )
         return TRUE;
      else
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
      /* todo the problem here are the relaxed bounds */
      return FALSE;
   }
   else if( strcmp(conshdlrname, "orbisack") == 0 )
   {
      return TRUE;
   }
   else if( strcmp(conshdlrname, "orbitope") == 0 )
   {
      return TRUE;
   }
   else if( strcmp(conshdlrname, "and") == 0 )
   {
      return TRUE;
   }
   if( strcmp(conshdlrname, "xor") == 0 )
   {
      return TRUE;
   }
   if( strcmp(conshdlrname, "or") == 0 )
   {
      return TRUE;
   }
   if( strcmp(conshdlrname, "bounddisjunction") == 0 )
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
void cleanBdchgQueue(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RESOLUTIONSET*   resolutionset       /**< resolution set */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_VAR* var;

   int i;/*lint !e850*/

   assert(conflict != NULL);

   for( i = SCIPpqueueNElems(conflict->resbdchgqueue) - 1; i >= 0; --i )/*lint !e850*/
   {
      int j;
      SCIP_Bool idxinrow;
      SCIP_Real val;
      int idxvar;

      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->resbdchgqueue)[i]);

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

            /* todo fix this hack */
            if ( i == SCIPpqueueNElems(conflict->resbdchgqueue))
               SCIPpqueueDelPos(conflict->resbdchgqueue, i);
            else if ( i != SCIPpqueueNElems(conflict->resbdchgqueue) )
            {
               SCIPpqueueDelPos(conflict->resbdchgqueue, i);
               i = SCIPpqueueNElems(conflict->resbdchgqueue);/*lint !e850*/
            }
      }
   }

}/*lint !e850*/

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

/** gets number of calls that resolution conflict analysis stopped for an unknown reason*/
SCIP_Longint SCIPconflictGetNResUnkTerm(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nrescalls - conflict->nressuccess - conflict->ncorrectaborts;
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
   (*resolutionset)->isbinary = FALSE;

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
   SCIP_CALL( resolutionsetCreate(&conflict->resolvedresolutionset, blkmem) );

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
   resolutionset->isbinary = FALSE;
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
      resolutionset->lhs -= resolutionset->vals[pos] * SCIPvarGetUbGlobal(var);
   }
   /* weaken with global lower bound */
   else
   {
      assert( SCIPsetIsLT(set, resolutionset->vals[pos], 0.0) );
      resolutionset->lhs -= resolutionset->vals[pos] * SCIPvarGetLbGlobal(var);
   }

   --resolutionset->nnz;
   resolutionset->vals[pos] = resolutionset->vals[resolutionset->nnz];
   resolutionset->inds[pos] = resolutionset->inds[resolutionset->nnz];
   SCIPsortIntReal(resolutionset->inds, resolutionset->vals, resolutionset->nnz);
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

/* returns the index of a variable in the conflict resolution set */
static
int getVarIdxInResolutionset(
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
         return i;
   }
   return -1;
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
         SCIP_Real bnd;
         SCIP_Real tightval;

         if( val > 0.0 )
         {
            bnd = SCIPvarGetUbGlobal(vars[v]);
            tightval = (val * SCIPvarGetUbAtIndex(vars[v], bdchgidx, TRUE) - resolutionset->slack ) / val;
            /* it the case of >= constraints we get tight propagation at the point (coef * ub - slack) / coef */
            solval = MIN( bnd, tightval );
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
            bnd = SCIPvarGetLbGlobal(vars[v]);
            tightval = (val * SCIPvarGetLbAtIndex(vars[v], bdchgidx, TRUE) - resolutionset->slack ) / val;
            solval = MAX( bnd, tightval );
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
      SCIP_VAR** vars;
      vars = SCIPprobGetVars(prob);

      for( i = 0; i < SCIPaggrRowGetNNz(aggrrow); i++ )
      {
         assert(SCIPsetIsEQ(set, -resolutionset->vals[i], SCIPaggrRowGetProbvarValue(aggrrow, SCIPvarGetProbindex(vars[resolutionset->inds[i]]))));
      }
      SCIP_UNUSED(vars);
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
SCIP_RETCODE updateBdchgQueue(
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
               SCIP_CALL( SCIPaddConflictUb(set->scip, vars[v], inferbdchgidx, FALSE) );
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
               SCIP_CALL( SCIPaddConflictLb(set->scip, vars[v], inferbdchgidx, FALSE) );
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
}/*lint !e715*/

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

   /* calculate the maximal size of each accepted conflict set */
   maxsize = (int) (set->conf_maxvarsfracres * transprob->nvars);

   SCIPsetDebugMsg(set, "flushing %d resolution sets at focus depth %d (vd: %d, cd: %d, rd: %d, maxsize: %d)\n",
      1, focusdepth, resolutionset->validdepth, resolutionset->conflictdepth, resolutionset->repropdepth, maxsize);

   assert(resolutionset != NULL);
   assert(resolutionset->validdepth == 0);

   *success = FALSE;
   /* do not add long conflicts */
   if( resolutionsetGetNNzs(resolutionset) > maxsize )
   {
      conflict->ncorrectaborts++;
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
      conflict->ncorrectaborts++;
      SCIP_CALL( SCIPnodeCutoff(tree->path[resolutionset->validdepth], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
      return SCIP_OKAY;
   }
   /* if the conflict has a relaxation only variable it is not generated at the moment */
   else if( hasRelaxationOnlyVar(set, transprob, resolutionset) )
   {
      conflict->ncorrectaborts++;
      SCIPsetDebugMsg(set, " -> resolution set has relaxation only variable \n");
      return SCIP_OKAY;
   }
   else
   {
      /* @todo use the right insert depth and not valid depth */
      SCIP_CALL( createAndAddResolutionCons(conflict, blkmem, set, stat, transprob, origprob, \
                     tree, reopt, lp, cliquetable, resolutionset, resolutionset->validdepth, success) );
      conflict->nappliedglbresconss++;
      SCIPsetDebugMsg(set, " -> resolution set added (cdpt:%d, fdpt:%d, insert:%d, valid:%d, conf: %d, reprop: %d , len:%d):\n",
                     SCIPtreeGetCurrentDepth(tree), SCIPtreeGetFocusDepth(tree),
                     resolutionset->validdepth, resolutionset->validdepth, resolutionset->conflictdepth,
                     resolutionset->repropdepth, resolutionsetGetNNzs(resolutionset));
   }

   return SCIP_OKAY;
}/*lint !e715*/

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
   int                   residx                 /**< index of variable to resolve */
   )
{
   SCIP_Real coefconf;
   SCIP_Real coefreas;
   SCIP_Real scale;

   SCIP_Bool idxinconflict;
   SCIP_Bool idxinreason;

   coefconf = 0.0;
   coefreas = 0.0;

   /* find in the conflict resolution set the coefficient of the variable we are resolving */
   idxinconflict = getCoefInResolutionSet(conflictresolutionset, residx, &coefconf);

   /* find in the reason resolution set the coefficient of the variable we are resolving */
   idxinreason = getCoefInResolutionSet(reasonresolutionset, residx, &coefreas);

   assert((idxinconflict && idxinreason));
   assert(!SCIPsetIsZero(set, coefreas) && !SCIPsetIsZero(set, coefconf));

   assert(coefconf * coefreas < 0);

   scale = REALABS( coefconf / coefreas );

   SCIP_UNUSED(idxinconflict);
   SCIP_UNUSED(idxinreason);

   return scale;

}

/** compute the resolved resolution set conflict + scale * reason */
static
SCIP_RETCODE rescaleAndResolve(
   SCIP_SET*             set,                   /**< global SCIP settings */
   SCIP_RESOLUTIONSET*   conflictresolutionset, /**< conflict resolution set */
   SCIP_RESOLUTIONSET*   reasonresolutionset,   /**< reason resolution set */
   SCIP_RESOLUTIONSET*   resolvedresolutionset, /**< reason resolution set */
   BMS_BLKMEM*           blkmem,                /**< block memory */
   int                   residx,                /**< index of variable to resolve */
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

   scale = computeScaleReason(set, conflictresolutionset, reasonresolutionset, residx);

   /* stop if the scale becomes too large */
   if ( SCIPsetIsGE(set, scale,  set->conf_generalresminmaxquot) )
   {
      return SCIP_OKAY;
   }

   SCIP_CALL( resolutionsetReplace(resolvedresolutionset, blkmem, conflictresolutionset) );

   newsize = resolutionsetGetNNzs(resolvedresolutionset) + resolutionsetGetNNzs(reasonresolutionset);
   if ( resolvedresolutionset->size < newsize )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &resolvedresolutionset->inds, resolvedresolutionset->size, newsize ) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &resolvedresolutionset->vals, resolvedresolutionset->size, newsize ) );
      resolvedresolutionset->size = newsize;
   }

   cidx = 0;
   previousnnz = resolutionsetGetNNzs(resolvedresolutionset);

   i = 0;
   /* add conflict and reason resolution sets */
   while ( i < resolutionsetGetNNzs(reasonresolutionset) )
   {
      if (cidx >= previousnnz)
      {
         resolvedresolutionset->inds[resolvedresolutionset->nnz] = reasonresolutionset->inds[i];
         resolvedresolutionset->vals[resolvedresolutionset->nnz] = scale * reasonresolutionset->vals[i];
         resolvedresolutionset->nnz++;
         i++;
      }
      else if (reasonresolutionset->inds[i] == resolvedresolutionset->inds[cidx])
      {
         /* @todo quadprecision when adding coefficients? */
         resolvedresolutionset->vals[cidx] = resolvedresolutionset->vals[cidx] + scale * reasonresolutionset->vals[i];
         cidx++;
         i++;
      }
      else if (reasonresolutionset->inds[i] > resolvedresolutionset->inds[cidx])
      {
         resolvedresolutionset->vals[cidx] = resolvedresolutionset->vals[cidx];
         cidx++;
      }
      else if (reasonresolutionset->inds[i] < resolvedresolutionset->inds[cidx])
      {
         resolvedresolutionset->inds[resolvedresolutionset->nnz] = reasonresolutionset->inds[i];
         resolvedresolutionset->vals[resolvedresolutionset->nnz] = scale * reasonresolutionset->vals[i];
         resolvedresolutionset->nnz++;
         i++;
      }
   }
   resolvedresolutionset->lhs = resolvedresolutionset->lhs + scale * reasonresolutionset->lhs;

   newnnz = resolutionsetGetNNzs(resolvedresolutionset);

   largestcoef = -SCIPsetInfinity(set);
   smallestcoef = SCIPsetInfinity(set);
   /* remove coefficients that are almost zero (10^-9 tolerance), loop backwards */
   for( i = newnnz - 1; i >= 0 ; i-- )
   {
      if (SCIPsetIsZero(set, resolvedresolutionset->vals[i] ))
      {
         resolutionsetRemoveZeroVar(resolvedresolutionset, set, i);
      }
      else
      {
         smallestcoef = MIN(smallestcoef, resolvedresolutionset->vals[i]);
         largestcoef = MAX(largestcoef, resolvedresolutionset->vals[i]);
      }
   }
   SCIPsetDebugMsg(set, "Nonzeros in resolved constraint: %d \n", resolutionsetGetNNzs(resolvedresolutionset));

   /* check if the quotient of coefficients in the resolvent exceeds the max allowed quotient */
   resolvedresolutionset->coefquotient = (resolvedresolutionset->nnz > 0) ? fabs(largestcoef / smallestcoef) : 0.0;
   if ( SCIPsetIsGT(set, resolvedresolutionset->coefquotient, set->conf_generalresminmaxquot) )
   {
      SCIPsetDebugMsg(set, "Quotient %f exceeds max allowed quotient", (resolvedresolutionset->nnz > 0) ? fabs(largestcoef / smallestcoef) : 0.0);
      return SCIP_OKAY;
   }

   /* sort for linear time resolution */
   SCIPsortIntReal(resolvedresolutionset->inds, resolvedresolutionset->vals, resolutionsetGetNNzs(resolvedresolutionset));
   *success = TRUE;

   return SCIP_OKAY;
}


/** Apply Division based reduction:
 * - Iteratively weaken variables from the resolution set that do not affect the
 *   slack and are not divisible by the .
 * - Then apply MIR or CG to reduce the slack.
 * - We weaken a variable if:
 *   * it is not divisible by the coefficient in the conflict resolution set AND
 *    * it is free (@todo favor this case), or
 *    * it has a positive coefficient and its local upper bound is equal to the
 *      global upper bound, or
 *    * it has a negative coefficient and its local lower bound is equal to the
 *      global lower bound
 * - We weaken as long as the resolved slack is positive.
 * Remark: For binary variables, the slack of the resolvent should become negative
 *         at the latest when no more variables can be weakened.
*/
static
SCIP_RETCODE DivisionBasedReduction(
   SCIP_CONFLICT *       conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_BDCHGIDX*        currbdchgidx,       /**< index of current bound change */
   int                   residx,             /**< index of variable to resolve */
   int*                  nvarsweakened,      /**< number of weakened variables */
   SCIP_Real*            fixbounds,          /**< array of fixed bounds */
   int*                  fixinds,            /**< array of indices of fixed variables */
   SCIP_Bool*            successresolution   /**< pointer to store whether the resolution was successful */
   )
{
   SCIP_VAR** vars;
   SCIP_Real previousslack;
   SCIP_RESOLUTIONSET* reasonset;
   int i;
   SCIP_Bool applydivision;
   int idxinreason;
   SCIP_Real coefinreason;

   assert(conflict != NULL);

   reasonset = conflict->reasonset;

   vars = SCIPprobGetVars(prob);
   i = 0;
   *nvarsweakened = 0;
   applydivision = TRUE;
   previousslack = reasonset->slack;

   idxinreason = getVarIdxInResolutionset(reasonset, residx);
   assert(idxinreason >= 0);

   coefinreason = fabs(reasonset->vals[idxinreason]);
   SCIP_CALL( rescaleAndResolve(set, conflict->resolutionset, reasonset, conflict->resolvedresolutionset, blkmem,
                        residx, successresolution) );

   /* todo extend Chvatal-Gomory and MIR for constraints with general integer variables */
   /* MIR can also be used in the presence of continuous variables */
   if (!isBinaryResolutionSet(set, prob, reasonset))
   {
      SCIPsetDebugMsg(set, "Normalized Chvatal-Gomory and MIR are only implemented for binary constraints\n");
      return SCIP_OKAY;
   }

   while ( SCIPsetIsGE(set, getSlack(set, prob, conflict->resolvedresolutionset, currbdchgidx, fixbounds, fixinds), 0.0) && applydivision )
   {
      applydivision = FALSE;

      /* loop over all variables and weaken one by one */
      for( i = 0; i < resolutionsetGetNNzs(reasonset); ++i )
      {
         SCIP_Bool varwasweakened;
         SCIP_VAR* vartoweaken;

         varwasweakened = FALSE;
         vartoweaken = vars[reasonset->inds[i]];

         if ( reasonset->inds[i] != residx )
         {
            if( reasonset->vals[i] > 0.0 && !SCIPsetIsZero(set, fmod(reasonset->vals[i], coefinreason)) )
            {
               SCIP_Real ub;
               ub = SCIPgetVarUbAtIndex(set->scip, vartoweaken, currbdchgidx, TRUE);

               if( SCIPsetIsEQ(set, ub, SCIPvarGetUbGlobal(vartoweaken)) )
               {
                  weakenVarReason(reasonset, set, vartoweaken, i);
                  varwasweakened = TRUE;
                  applydivision = TRUE;
                  ++(*nvarsweakened);
               }
            }
            else if( reasonset->vals[i] < 0.0 && !SCIPsetIsZero(set, fmod(reasonset->vals[i], coefinreason)) )
            {
               SCIP_Real lb;
               lb = SCIPgetVarLbAtIndex(set->scip, vartoweaken, currbdchgidx, TRUE);
               if( SCIPsetIsEQ(set, lb, SCIPvarGetLbGlobal(vartoweaken)) )
               {
                  weakenVarReason(reasonset, set, vartoweaken, i);
                  varwasweakened = TRUE;
                  applydivision = TRUE;
                  ++(*nvarsweakened);
               }
            }
         }
         if (varwasweakened)
         {
            if ( !set->conf_weakenreasonall && (set->conf_batchcoeftight > 0) &&
               (*nvarsweakened % set->conf_batchcoeftight == 0) )
            {
               SCIP_RESOLUTIONSET *reducedreason;
               resolutionsetCopy(&reducedreason, blkmem, reasonset);

               /* apply the chosen reduction technique */
               if (set->conf_reductiontechnique == 'd')
                  SCIP_CALL( ChvatalGomoryLhs(set, prob, reducedreason, coefinreason) );
               else
                  SCIP_CALL( MirLhs(set, prob, reducedreason, coefinreason) );

               SCIPsortIntReal(reducedreason->inds, reducedreason->vals, resolutionsetGetNNzs(reducedreason));
               /* todo the update of the slack should be done in the division/mir algorithm */
               reducedreason->slack = getSlack(set, prob, reducedreason, currbdchgidx, fixbounds, fixinds);
               SCIPsetDebugMsg(set, "slack before tightening: %g \n",  previousslack);
               SCIPsetDebugMsg(set, "slack after tightening: %g \n", reducedreason->slack);
               SCIP_CALL( rescaleAndResolve(set, conflict->resolutionset, reducedreason, conflict->resolvedresolutionset, blkmem, residx, successresolution) );
               SCIPresolutionsetFree(&reducedreason, blkmem);
               break;
            }
         }
      }
   }
   /* weaken all and apply strengthening once */
   if ( set->conf_weakenreasonall && *nvarsweakened > 0 )
   {
      SCIP_RESOLUTIONSET *reducedreason;
      resolutionsetCopy(&reducedreason, blkmem, reasonset);

      /* apply the chosen reduction technique */
      if (set->conf_reductiontechnique == 'd')
         SCIP_CALL( ChvatalGomoryLhs(set, prob, reducedreason, coefinreason) );
      else
         SCIP_CALL( MirLhs(set, prob, reducedreason, coefinreason) );

      SCIPsortIntReal(reducedreason->inds, reducedreason->vals, resolutionsetGetNNzs(reducedreason));

      /* todo the update of the slack should be included in the division algorithms */
      reducedreason->slack = getSlack(set, prob, reducedreason, currbdchgidx, fixbounds, fixinds);
      SCIPsetDebugMsg(set, "slack before tightening: %g \n",  previousslack);
      SCIPsetDebugMsg(set, "slack after tightening: %g \n", reducedreason->slack);
      SCIP_CALL( rescaleAndResolve(set, conflict->resolutionset, reducedreason, conflict->resolvedresolutionset, blkmem, residx, successresolution) );
      SCIPresolutionsetFree(&reducedreason, blkmem);
   }

   SCIP_UNUSED(previousslack);

   return SCIP_OKAY;
}

/** Apply tightening based reduction:
 * - Iteratively weaken variables from the resolution set that do not affect the
 *   slack. Then apply coefficient tightening to reduce the slack.
 * - We weaken a variable if:
 *    * it is free (@todo favor this case)
 *    * it has a positive coefficient and its local upper bound is equal to the
 *      global upper bound
 *    * it has a negative coefficient and its local lower bound is equal to the
 *      global lower bound
 * - For the reason resolution set (isreason = TRUE) we weaken as long as the linear combination of
 *   slack_conflict and slack_reason is positive. i. e. slack_conflict + scale * slack_reason >= 0
*/
static
SCIP_RETCODE tighteningBasedReduction(
   SCIP_CONFLICT *       conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_BDCHGIDX*        currbdchgidx,       /**< index of current bound change */
   int                   residx,             /**< index of variable to resolve */
   int*                  nvarsweakened,      /**< number of weakened variables */
   SCIP_Real*            fixbounds,          /**< array of fixed bounds */
   int*                  fixinds,            /**< array of indices of fixed variables */
   SCIP_Bool*            successresolution   /**< pointer to store whether the resolution was successful */
   )
{
   SCIP_VAR** vars;
   SCIP_Real previousslack;
   SCIP_RESOLUTIONSET* resolutionset;
   int i;
   SCIP_Bool applytightening;
   int nchgcoefs;

   assert(conflict != NULL);

   resolutionset = conflict->reasonset;

   vars = SCIPprobGetVars(prob);
   i = 0;
   *nvarsweakened = 0;
   applytightening = TRUE;
   previousslack = resolutionset->slack;

   SCIP_CALL( rescaleAndResolve(set, conflict->resolutionset, resolutionset, conflict->resolvedresolutionset, blkmem,
                        residx, successresolution) );

   while ( SCIPsetIsGE(set, getSlack(set, prob, conflict->resolvedresolutionset, currbdchgidx, fixbounds, fixinds), 0.0) &&         applytightening )
   {
      applytightening = FALSE;
      for( i = 0; i < resolutionsetGetNNzs(resolutionset); ++i )
      {
         SCIP_Bool varwasweakened;
         SCIP_VAR* vartoweaken;

         varwasweakened = FALSE;
         vartoweaken = vars[resolutionset->inds[i]];

         if ( resolutionset->inds[i] != residx )
         {
            if( resolutionset->vals[i] > 0.0 )
            {
               SCIP_Real ub;

               ub = SCIPgetVarUbAtIndex(set->scip, vartoweaken, currbdchgidx, TRUE);

               if( SCIPsetIsEQ(set, ub, SCIPvarGetUbGlobal(vartoweaken)) )
               {
                  weakenVarReason(resolutionset, set, vartoweaken, i);
                  varwasweakened = TRUE;
                  applytightening = TRUE;
                  ++(*nvarsweakened);
               }
            }
            else
            {
               SCIP_Real lb;

               lb = SCIPgetVarLbAtIndex(set->scip, vartoweaken, currbdchgidx, TRUE);

               if( SCIPsetIsEQ(set, lb, SCIPvarGetLbGlobal(vartoweaken)) )
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
            if ( !set->conf_weakenreasonall && (set->conf_batchcoeftight > 0) &&
               (*nvarsweakened % set->conf_batchcoeftight == 0) )
            {
               SCIP_CALL( tightenCoefLhs(set, prob, FALSE, resolutionset->vals, &resolutionset->lhs,
                              resolutionset->inds, &resolutionset->nnz, &nchgcoefs, NULL) );
               previousslack = resolutionset->slack;
               SCIPsortIntReal(resolutionset->inds, resolutionset->vals, resolutionsetGetNNzs(resolutionset));

               /* todo the update of the slack should be included in tightenCoefLhs */
               resolutionset->slack = getSlack(set, prob, resolutionset, currbdchgidx, fixbounds, fixinds);
               SCIPsetDebugMsg(set, "slack after tightening: %g \n", resolutionset->slack);
               assert(SCIPsetIsLE(set, resolutionset->slack, previousslack + EPS));
               SCIP_CALL( rescaleAndResolve(set, conflict->resolutionset, resolutionset, conflict->resolvedresolutionset, blkmem,
                                    residx, successresolution) );
               break;
            }
         }
      }
   }
   if ( set->conf_weakenreasonall && *nvarsweakened > 0 )
   {
      SCIP_CALL( tightenCoefLhs(set, prob, FALSE, resolutionset->vals, &resolutionset->lhs,
                     resolutionset->inds, &resolutionset->nnz, &nchgcoefs, NULL) );
      previousslack = resolutionset->slack;
      SCIPsortIntReal(resolutionset->inds, resolutionset->vals, resolutionsetGetNNzs(resolutionset));

      /* todo the update of the slack should be included in tightenCoefLhs */
      resolutionset->slack = getSlack(set, prob, resolutionset, currbdchgidx, fixbounds, fixinds);
      SCIPsetDebugMsg(set, "slack after tightening: %g \n", resolutionset->slack);
      assert(SCIPsetIsLE(set, resolutionset->slack, previousslack + EPS));
      SCIP_CALL( rescaleAndResolve(set, conflict->resolutionset, resolutionset, conflict->resolvedresolutionset, blkmem,
                           residx, successresolution) );
   }

   SCIP_UNUSED(previousslack);

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
      SCIPsetDebugMsgPrint(set, " -> Bound change <%s> %s %.15g [status:%d, type:%d, depth:%d, pos:%d, reason:<%s>, info:%d] \n",
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

/** tries to resolve given bound change */
static
SCIP_RETCODE reasonBoundChanges(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change to resolve */
   SCIP_Real             relaxedbd,          /**< the relaxed bound */
   int                   validdepth,         /**< minimal depth level at which the conflict is valid */
   SCIP_Bool*            resolved            /**< pointer to store whether the bound change was resolved */
   )
{
   SCIP_VAR* actvar;
   SCIP_CONS* infercons;
   SCIP_PROP* inferprop;
   SCIP_RESULT result;

#ifndef NDEBUG

   /* store the current size of the conflict queues */
   assert(conflict != NULL);
#else
   assert(conflict != NULL);
#endif

   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   *resolved = FALSE;

   actvar = SCIPbdchginfoGetVar(bdchginfo);
   assert(actvar != NULL);
   assert(SCIPvarIsActive(actvar));

   /* check, if the bound change can and should be resolved:
    *  - the reason must be either a global constraint
    *  - @todo later a propagator
    */
   switch( SCIPbdchginfoGetChgtype(bdchginfo) )
   {
   case SCIP_BOUNDCHGTYPE_CONSINFER:
      infercons = SCIPbdchginfoGetInferCons(bdchginfo);
      assert(infercons != NULL);

      if( SCIPconsIsGlobal(infercons) || SCIPconsGetValidDepth(infercons) <= validdepth )
      {
         SCIP_VAR* infervar;
         int inferinfo;
         SCIP_BOUNDTYPE inferboundtype;
         SCIP_BDCHGIDX* bdchgidx;

         /* resolve bound change by asking the constraint that infered the bound to put all bounds that were
          * the reasons for the conflicting bound change on the priority queue
          */
         infervar = SCIPbdchginfoGetInferVar(bdchginfo);
         inferinfo = SCIPbdchginfoGetInferInfo(bdchginfo);
         inferboundtype = SCIPbdchginfoGetInferBoundtype(bdchginfo);
         bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
         assert(infervar != NULL);

         SCIPsetDebugMsg(set, "getting reason for <%s> %s %g(%g) [status:%d, type:%d, depth:%d, pos:%d]: <%s> %s %g [cons:<%s>(%s), info:%d]\n",
            SCIPvarGetName(actvar),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo), relaxedbd,
            SCIPvarGetStatus(actvar), SCIPvarGetType(actvar),
            SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
            SCIPvarGetName(infervar),
            inferboundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPgetVarBdAtIndex(set->scip, infervar, inferboundtype, bdchgidx, TRUE),
            SCIPconsGetName(infercons),
            SCIPconsIsGlobal(infercons) ? "global" : "local",
            inferinfo);

         /* in case the inference variables is not an active variables, we need to transform the relaxed bound */
         if( actvar != infervar )
         {
            SCIP_VAR* var;
            SCIP_Real scalar;
            SCIP_Real constant;

            assert(SCIPvarGetStatus(infervar) == SCIP_VARSTATUS_AGGREGATED
               || SCIPvarGetStatus(infervar) == SCIP_VARSTATUS_NEGATED
               || (SCIPvarGetStatus(infervar) == SCIP_VARSTATUS_MULTAGGR && SCIPvarGetMultaggrNVars(infervar) == 1));

            scalar = 1.0;
            constant = 0.0;

            var = infervar;

            /* transform given varibale to active varibale */
            SCIP_CALL( SCIPvarGetProbvarSum(&var, set, &scalar, &constant) );
            assert(var == actvar);

            relaxedbd *= scalar;
            relaxedbd += constant;
         }

         SCIP_CALL( SCIPconsResolvePropagation(infercons, set, infervar, inferinfo, inferboundtype, bdchgidx, relaxedbd, TRUE, &result) );
         *resolved = (result == SCIP_SUCCESS);
      }
      break;

   case SCIP_BOUNDCHGTYPE_PROPINFER:
      inferprop = SCIPbdchginfoGetInferProp(bdchginfo);
      if( inferprop != NULL )
      {
         SCIP_VAR* infervar;
         int inferinfo;
         SCIP_BOUNDTYPE inferboundtype;
         SCIP_BDCHGIDX* bdchgidx;

         /* resolve bound change by asking the propagator that infered the bound to put all bounds that were
          * the reasons for the conflicting bound change on the priority queue
          */
         infervar = SCIPbdchginfoGetInferVar(bdchginfo);
         inferinfo = SCIPbdchginfoGetInferInfo(bdchginfo);
         inferboundtype = SCIPbdchginfoGetInferBoundtype(bdchginfo);
         bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
         assert(infervar != NULL);

         SCIPsetDebugMsg(set, "getting reason for <%s> %s %g(%g) [status:%d, depth:%d, pos:%d]: <%s> %s %g [prop:<%s>, info:%d]\n",
            SCIPvarGetName(actvar),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo), relaxedbd,
            SCIPvarGetStatus(actvar), SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
            SCIPvarGetName(infervar),
            inferboundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPgetVarBdAtIndex(set->scip, infervar, inferboundtype, bdchgidx, TRUE),
            SCIPpropGetName(inferprop), inferinfo);

         SCIP_CALL( SCIPpropResolvePropagation(inferprop, set, infervar, inferinfo, inferboundtype, bdchgidx, relaxedbd, TRUE, &result) );
         *resolved = (result == SCIP_SUCCESS);
      }
      break;

   case SCIP_BOUNDCHGTYPE_BRANCHING:
      assert(!(*resolved));
      break;

   default:
      SCIPerrorMessage("invalid bound change type <%d>\n", SCIPbdchginfoGetChgtype(bdchginfo));
      return SCIP_INVALIDDATA;
   }

   SCIPsetDebugMsg(set, "resolving status: %u\n", *resolved);

#ifndef NDEBUG

   /* in case the bound change was not resolved, the separate conflict queue should have zero elements */
   assert((*resolved) || (SCIPpqueueNElems(conflict->separatebdchgqueue) == 0));
#endif

   return SCIP_OKAY;
}

/* get a conflict resolution set from bound changes */
static
SCIP_RETCODE getClauseConflictSet(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       currbdchginfo,      /**< bound change to resolve */
   SCIP_Bool*            success             /**< pointer to store whether we could find an  initial conflict */
)
{

   *success = FALSE;
   if (!SCIPvarIsBinary(SCIPbdchginfoGetVar(currbdchginfo)))
   {
      conflict->ncorrectaborts++;
      SCIPsetDebugMsg(set, "Could not obtain an initial conflict set \n");
      return SCIP_OKAY;
   }
   else
   {
      SCIP_Bool isbinary;
      SCIP_Real lhs;

      /* todo stop if current bound change is on a non-binary */
      isbinary = TRUE;
      lhs = 1.0;

      /** given the set of bound changes that renders infeasibility create a non-good cut
       *  as initial conflict. E.g. if x = 1, y = 1, and z = 0 leads to infeasibility,
       *  then the initial conflict constraint is (1 - x) + (1 - y) + z >= 1
       */
      for(int i = 0; i < SCIPpqueueNElems(conflict->resbdchgqueue); i++)
      {
         SCIP_BDCHGINFO* bdchginfo;
         bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->resbdchgqueue)[i]);

         if( !SCIPvarIsBinary(SCIPbdchginfoGetVar(bdchginfo)) )
         {
            isbinary = FALSE;
            break;
         }
         if (SCIPbdchginfoGetNewbound(bdchginfo) >= 0.5)
            lhs--;
      }

      if( isbinary )
      {
         conflict->resolutionset->nnz = SCIPpqueueNElems(conflict->resbdchgqueue) + 1;
         conflict->resolutionset->lhs = lhs;
         conflict->resolutionset->origlhs = lhs;
         conflict->resolutionset->origrhs = SCIPsetInfinity(set);

         if( conflict->resolutionset->size == 0 )
         {
            assert(conflict->resolutionset->vals == NULL);
            assert(conflict->resolutionset->inds == NULL);

            SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflict->resolutionset->vals, conflict->resolutionset->nnz) );
            SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflict->resolutionset->inds, conflict->resolutionset->nnz) );
            conflict->resolutionset->size = conflict->resolutionset->nnz;
         }

         else if( conflict->resolutionset->size < conflict->resolutionset->nnz )
         {
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflict->resolutionset->vals, conflict->resolutionset->size, conflict->resolutionset->nnz) );
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflict->resolutionset->inds, conflict->resolutionset->size, conflict->resolutionset->nnz) );
            conflict->resolutionset->size = conflict->resolutionset->nnz;
         }

         /* for the current bound change */
         if (SCIPbdchginfoGetNewbound(currbdchginfo) > 0.5)
         {
            conflict->resolutionset->vals[0] = -1.0;
            conflict->resolutionset->lhs -= 1.0;
         }
         else
            conflict->resolutionset->vals[0] = 1.0;
         conflict->resolutionset->inds[0] = SCIPvarGetProbindex(SCIPbdchginfoGetVar(currbdchginfo));

         for( int i = 0; i < conflict->resolutionset->nnz - 1; i++ )
         {
            SCIP_BDCHGINFO* bdchginfo;
            bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->resbdchgqueue)[i]);
            conflict->resolutionset->vals[i+1] = SCIPbdchginfoGetNewbound(bdchginfo) > 0.5 ? -1.0 : 1.0;
            conflict->resolutionset->inds[i+1] = SCIPvarGetProbindex(SCIPbdchginfoGetVar(bdchginfo));
         }
         *success = TRUE;
      }
      else
      {
         conflict->ncorrectaborts++;
         SCIPsetDebugMsg(set, "Could not obtain an initial conflict set \n");
         return SCIP_OKAY;
      }
   }
   return SCIP_OKAY;
}

/* get a resolution set from bound changes */
static
SCIP_RETCODE getClauseReasonSet(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       currbdchginfo,      /**< bound change to resolve */
   SCIP_Real             relaxedbd,          /**< the relaxed bound */
   int                   validdepth,         /**< minimal depth level at which the conflict is valid */
   SCIP_Bool*            success             /**< pointer to store whether we could find a reason*/
)
{

   /* if the current bound change is on a non-binary variable then we cannot find a linear reason */
   if( !SCIPvarIsBinary(SCIPbdchginfoGetVar(currbdchginfo)) || !reasonIsLinearizable(currbdchginfo))
   {
      conflict->ncorrectaborts++;
      *success = FALSE;
      return SCIP_OKAY;
   }
   /* make sure that the separate bound change queue is empty */
   if (SCIPpqueueNElems(conflict->separatebdchgqueue) != 0)
   {
      SCIPpqueueClear(conflict->separatebdchgqueue);
   }
   assert(SCIPpqueueNElems(conflict->separatebdchgqueue) == 0);

   SCIP_CALL( reasonBoundChanges(conflict, set, currbdchginfo, relaxedbd, validdepth, success) );
   if ( SCIPpqueueNElems(conflict->separatebdchgqueue) == 0 )
   {
      SCIPsetDebugMsg(set, "Could not obtain a reason set \n");
      *success = FALSE;
      return SCIP_OKAY;
   }
   else
   {
      SCIP_Bool isbinary;
      SCIP_Real lhs;

      *success = FALSE;
      isbinary = TRUE;
      lhs = 1.0;

      /** given the set of bound changes that leads to propagation of the current
       *  bound change, create a clause reason resolution set
       *  E.g. if x = 1, y = 1, leads to z = 0 then the reason
       *  constraint is (1 - x) + (1 - y) + (1 - z) >= 1
       */
      for(int i = 0; i < SCIPpqueueNElems(conflict->separatebdchgqueue); i++)
      {
         SCIP_BDCHGINFO* bdchginfo;
         bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->separatebdchgqueue)[i]);

         if( !SCIPvarIsBinary(SCIPbdchginfoGetVar(bdchginfo)) )
         {
            isbinary = FALSE;
            break;
         }
         if (SCIPbdchginfoGetNewbound(bdchginfo) == 1.0)
            lhs--;
      }

      if( isbinary )
      {
         conflict->reasonset->nnz = SCIPpqueueNElems(conflict->separatebdchgqueue) + 1;
         conflict->reasonset->lhs = lhs;
         conflict->reasonset->origlhs = lhs;
         conflict->reasonset->origrhs = SCIPsetInfinity(set);

         if( conflict->reasonset->size == 0 )
         {
            assert(conflict->reasonset->vals == NULL);
            assert(conflict->reasonset->inds == NULL);

            /* todo the next line is a temporay fix for the vector size */
            SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflict->reasonset->vals, conflict->reasonset->nnz) );
            SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflict->reasonset->inds, conflict->reasonset->nnz) );
            conflict->reasonset->size = conflict->reasonset->nnz;
         }

         else if( conflict->reasonset->size < conflict->reasonset->nnz )
         {
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflict->reasonset->vals, conflict->reasonset->size, conflict->reasonset->nnz) );
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflict->reasonset->inds, conflict->reasonset->size, conflict->reasonset->nnz) );
            conflict->reasonset->size = conflict->reasonset->nnz;
         }
         /* add the variable we are resolving and update lhs */
         conflict->reasonset->vals[0] = SCIPbdchginfoGetNewbound(currbdchginfo) > 0.5 ? 1.0 : -1.0;
         conflict->reasonset->inds[0] = SCIPvarGetProbindex(SCIPbdchginfoGetVar(currbdchginfo));
         conflict->reasonset->lhs += SCIPbdchginfoGetNewbound(currbdchginfo) > 0.5 ? 0.0 : -1.0;

         for( int i = 0; i < conflict->reasonset->nnz - 1; i++ )
         {
            SCIP_BDCHGINFO* bdchginfo;
            bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->separatebdchgqueue)[i]);
            conflict->reasonset->vals[i+1] = SCIPbdchginfoGetNewbound(bdchginfo) > 0.5 ? -1.0 : 1.0;
            conflict->reasonset->inds[i+1] = SCIPvarGetProbindex(SCIPbdchginfoGetVar(bdchginfo));
         }
         *success = TRUE;
      }
      else
      {
         conflict->ncorrectaborts++;
         SCIPsetDebugMsg(set, "Could not obtain a reason set \n");
         *success = FALSE;
         return SCIP_OKAY;
      }
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE getReasonRow(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       currbdchginfo,      /**< bound change to resolve */
   int                   residx,             /**< index of the bound change to resolve */
   int                   validdepth,         /**< minimal depth level at which the conflict is valid */
   SCIP_Real*            fixbounds,          /**< array of fixed bounds */
   int*                  fixinds,            /**< array of indices of fixed variables */
   SCIP_Bool*            success             /**< pointer to store whether we could get a linear reason */
)
{
   assert(success !=  NULL);

   *success = FALSE;
   if (bdchginfoIsResolvable(currbdchginfo) && SCIPbdchginfoGetChgtype(currbdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER)
   {
         SCIP_CONS* reasoncon;
         SCIP_ROW* reasonrow;

         reasoncon = SCIPbdchginfoGetInferCons(currbdchginfo);

         if(!SCIPconsIsGlobal(reasoncon))
         {
            SCIPsetDebugMsg(set, "Reason constraint is not global \n");
            return SCIP_OKAY;
         }

         /* get the corresponding reason row */
         reasonrow = SCIPconsCreateRow(set->scip, reasoncon);

         /* in case of orbitope-, orbisack-, and-constaints we construct a linearized clause as reason */
         if( reasonrow == NULL )
         {
               SCIP_CALL( getClauseReasonSet(conflict, blkmem,  set, currbdchginfo, SCIPbdchginfoGetRelaxedBound(currbdchginfo), validdepth, success) );
               if (*success)
               {
                  assert(SCIPsetIsZero(set, getSlack(set, prob, conflict->reasonset, SCIPbdchginfoGetIdx(currbdchginfo), fixbounds, fixinds)));
                  conflict->reasonset->slack = 0.0;
                  return SCIP_OKAY;
               }
               else
                  return SCIP_OKAY;
         }

         /* get the resolution set of the reason row */
         *success = TRUE;
         SCIP_CALL( reasonResolutionsetFromRow(set, blkmem, reasonrow, conflict->reasonset, currbdchginfo) );
         /* it may happen that some specialized propagation took place and the real reason is not the constraint
            e.g. negated cliques in cons_knapsack or ranged row propagation in cons_linear. */

         /* this happens if the reason is a negated clique found in the knapsack constraint handler */
         if (strcmp(SCIPconshdlrGetName(reasoncon->conshdlr), "knapsack") != 0)
         {
            assert(!SCIPsetIsInfinity(set, -conflict->reasonset->lhs) || !SCIPsetIsInfinity(set, conflict->reasonset->lhs));
         }
         else if(SCIPsetIsInfinity(set, -conflict->reasonset->lhs) || SCIPsetIsInfinity(set, conflict->reasonset->lhs))
         {
            /* to be able to continue we construct a linearized clause as reason */
            SCIP_CALL( getClauseReasonSet(conflict, blkmem, set, currbdchginfo, SCIPbdchginfoGetRelaxedBound(currbdchginfo), validdepth, success) );
            if (*success)
            {
               assert(SCIPsetIsZero(set, getSlack(set, prob, conflict->reasonset, SCIPbdchginfoGetIdx(currbdchginfo), fixbounds, fixinds)));
               conflict->reasonset->slack = 0.0;
            }
            return SCIP_OKAY;

         }
         conflict->reasonset->slack = getSlack(set, prob, conflict->reasonset, SCIPbdchginfoGetIdx(currbdchginfo), fixbounds, fixinds);

         /* If the slack is greater than 0, we check that the reason actually
         propagated the variable we resolve. It propagates a variable x_i if
         (slack - a_i * (oldbound - newbound) is smaller than 0 */
         if (SCIPsetIsGT(set, conflict->reasonset->slack, 0.0))
         {
            SCIP_VAR* var;
            SCIP_BDCHGIDX* currbdchgidx;
            SCIP_Real coef;
            SCIP_Real boundusedinslack;

            currbdchgidx = SCIPbdchginfoGetIdx(currbdchginfo);
            var = SCIPbdchginfoGetVar(currbdchginfo);
            assert(var != NULL);
            assert(SCIPvarGetProbindex(var) == residx);

            coef = conflict->reasonset->vals[getVarIdxInResolutionset(conflict->reasonset, residx)];
            boundusedinslack = coef > 0 ? SCIPgetVarUbAtIndex(set->scip, var, currbdchgidx, TRUE) : SCIPgetVarLbAtIndex(set->scip, var, currbdchgidx, TRUE);

            if (!SCIPsetIsLT(set, conflict->reasonset->slack - coef * ( boundusedinslack - SCIPbdchginfoGetOldbound(currbdchginfo) ) , 0.0))
            {

               assert( (strcmp(SCIPconshdlrGetName(reasoncon->conshdlr), "knapsack") == 0) || (strcmp(SCIPconshdlrGetName(reasoncon->conshdlr), "linear") == 0) );
               SCIP_CALL( getClauseReasonSet(conflict, blkmem, set, currbdchginfo, SCIPbdchginfoGetRelaxedBound(currbdchginfo), validdepth, success) );
               if (*success)
               {
                  assert(SCIPsetIsZero(set, getSlack(set, prob, conflict->reasonset, currbdchgidx, fixbounds, fixinds)));
                  conflict->reasonset->slack = 0.0;
               }
               return SCIP_OKAY;
            }
         }
   }
   else if (bdchginfoIsResolvable(currbdchginfo) && SCIPbdchginfoGetChgtype(currbdchginfo) == SCIP_BOUNDCHGTYPE_PROPINFER)
   {
      SCIP_CALL( getClauseReasonSet(conflict, blkmem,  set, currbdchginfo, SCIPbdchginfoGetRelaxedBound(currbdchginfo), validdepth, success) );
      if (*success)
      {
         assert(SCIPsetIsZero(set, getSlack(set, prob, conflict->reasonset, SCIPbdchginfoGetIdx(currbdchginfo), fixbounds, fixinds)));
         conflict->reasonset->slack = 0.0;
      }
         return SCIP_OKAY;

   }
   else
   {
      SCIPsetDebugMsg(set, "Could not obtain a reason row \n");
      *success = FALSE;
   }
   return SCIP_OKAY;
}
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
   SCIP_RESOLUTIONSET *resolvedresolutionset;
   SCIP_RESOLUTIONSET *prevconflictresolutionset;
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BDCHGINFO* nextbdchginfo;
   SCIP_BDCHGIDX* bdchgidx;
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
   resolutionSetClear(conflict->resolvedresolutionset);
   conflictresolutionset = conflict->resolutionset;
   conflictresolutionset->usescutoffbound = usescutoffbound;
   reasonresolutionset = conflict->reasonset;
   prevconflictresolutionset = conflict->prevresolutionset;
   resolvedresolutionset = conflict->resolvedresolutionset;

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
   {
      conflict->ncorrectaborts++;
      return SCIP_OKAY;
   }

   /* calculate the maximal size of each accepted conflict set */
   maxsize = (int) (set->conf_maxvarsfracres * transprob->nvars);

   /* last bound change that led to infeasibility */
   bdchginfo = conflictFirstCand(conflict);

   /* if no bound change exists or none of them was infered by a resolvable
   constraint then we terminate */
   /* todo we should always be able to resolve a bound change. For this we
   need to implement callbacks in all constraint handlers and propagators  */
   if ( bdchginfo == NULL || !existsResolvablebdchginfo(conflict) )
   {
#ifdef SCIP_DEBUG
      /* if at least one bound change is in the queue, print them all */
      if(bdchginfo != NULL)
         printAllBoundChanges(conflict, set);
#endif
      conflict->ncorrectaborts++;
      SCIPsetDebugMsg(set, "Conflict analysis not applicable since no resolvable bounds exist \n");
      return SCIP_OKAY;
   }

   /* remove the last bound change */
   bdchginfo = conflictRemoveCand(conflict);
   bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
   bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);

   vartoresolve = bdchginfo->var;
   residx = SCIPvarGetProbindex(vartoresolve);

   /* get the next bound change */
   nextbdchginfo = conflictFirstCand(conflict);

   /* check if the variable we are resolving is active */
   assert(SCIPvarIsActive(vartoresolve));

   /* the valid depth is always 0 (global case) */
   conflictresolutionset->validdepth = validdepth;

   if (initialconflictrow != NULL)
   {
      if( SCIProwGetNNonz(initialconflictrow) > maxsize )
      {
         SCIPsetDebugMsg(set, "Number of nonzeros in conflict is larger than maxsize %d > %d\n",
                        SCIProwGetNNonz(initialconflictrow), maxsize);
         conflict->ncorrectaborts++;
         return SCIP_OKAY;
      }
      SCIPsetDebugMsg(set, "Initial conflict Row: %s \n", SCIProwGetName(initialconflictrow));
      /* get the resolution set of the conflict row */
      SCIP_CALL( conflictResolutionsetFromRow(set, blkmem, initialconflictrow, conflictresolutionset, bdchginfo) );

      /* compute the slack */
      conflictresolutionset->slack = getSlack(set, transprob, conflictresolutionset, bdchgidx, NULL, NULL);

   }

   /** if no row exists create the resolution set (if possible) from the bound changes that lead to infeasibility
    * Moreover the slack should be negative. If not (for a good reason) create again the resolution set
    * (if possible) from the bound changes.
    * The only cases where this may not be true is if the conflict is found:
    *  - by a negated clique in the knapsack constraint handler
    *  - by propagating a ranged row
    */
   if ( initialconflictrow == NULL || SCIPsetIsGE(set, conflictresolutionset->slack, 0.0) )
   {

      SCIP_Bool success;
      SCIP_CONSHDLR* conshdlr;
      SCIPsetDebugMsg(set, "Slack of conflict constraint is not negative \n");

      assert(!infeasibleLP && !pseudoobj);
      if (initialconflictrow != NULL)
      {
         conshdlr = SCIProwGetOriginConshdlr(initialconflictrow);
         SCIPsetDebugMsg(set, "%s",SCIPconshdlrGetName(conshdlr));
         /* relaxed assertion */
         assert(strcmp(SCIPconshdlrGetName(conshdlr), "knapsack") == 0 || strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0);
      }

      getClauseConflictSet(conflict, blkmem, set, bdchginfo, &success);
      if (!success)
      {
         SCIPsetDebugMsg(set, "Initial conflict could not be retrieved \n");
         return SCIP_OKAY;
      }
      assert( getSlack(set, transprob, conflictresolutionset, bdchgidx, NULL, NULL) == -1.0);
      conflictresolutionset->slack = -1.0;
   }

   SCIP_CALL( resolutionsetReplace(prevconflictresolutionset, blkmem, conflictresolutionset) );
   SCIP_CALL( resolutionsetReplace(resolvedresolutionset, blkmem, conflictresolutionset) );

   SCIPdebug(resolutionsetPrintRow(conflictresolutionset, set, transprob, 0));

   /* Apply coefficient tightening to the conflict constraint should never hurt */
   SCIP_CALL( tightenCoefLhs(set, transprob, FALSE, conflictresolutionset->vals, &conflictresolutionset->lhs,
                  conflictresolutionset->inds, &conflictresolutionset->nnz, &nchgcoefs, NULL) );

   if (nchgcoefs > 0)
      conflictslack = getSlack(set, transprob, conflictresolutionset, bdchgidx, NULL, NULL);
   else
      conflictslack = conflictresolutionset->slack;

   assert(SCIPsetIsGE(set, conflictresolutionset->slack + EPS, conflictslack));
   conflictresolutionset->slack = conflictslack;

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


   /** main loop: All-FUIP RESOLUTION
    * --------------------------------
    * - (we already have the initial conflict row and the first bound change to
    *   resolve)
    * - apply coefficient tightening to the conflict row
    * - Ignore & Continue: if we can't explain the bound change, i.e. the reason
    *   is a branching or non-linear then we ignore it and continue with the
    *   next bound change. We have to ignore all other bound changes for
    *   this variable ( @todo resolve with no-good)
    * - if the bound change is resolvable:
    *   * get the reason row for the bound change
    *   * apply coefficient tightening to the reason (maybe also cMIR?)
    *   * take the linear combination of the conflict row and the reason row
    *   * apply coefficient tightening to the resolved row (maybe also cMIR?)
    *       - if there is no other bound change in the queue from the same depth level
    *         then we are at a UIP -> keep this constraint and continue
    */
   while( TRUE )  /*lint !e716*/
   {
#ifdef SCIP_DEBUG
      {
         SCIPsetDebugMsgPrint(set, " Number of applied resolution steps %d \n", nressteps);
         printAllBoundChanges(conflict, set);
         SCIPsetDebugMsgPrint(set, " Current bound change already removed from the queue: \n");
         printSingleBoundChange(set, bdchginfo);
         if( nextbdchginfo != NULL )
         {
            SCIPsetDebugMsgPrint(set, " First element in the bound change queue: \n");
            printSingleBoundChange(set, nextbdchginfo);
         }
         else
         {
            SCIPsetDebugMsgPrint(set, " The bound change queue is empty\n");
         }
      }
#endif
      /** check if the bound change is resolvable. If it is not, we can ignore the bound change and continue
       * with the next one
       */
      if ( !bdchginfoIsResolvable(bdchginfo) && set->conf_fixandcontinue)
      {
#ifdef SCIP_DEBUG
            printNonResolvableReasonType(set, bdchginfo);
#endif
            if ( SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER )
            {
               reasoncon = SCIPbdchginfoGetInferCons(bdchginfo);
               /* we resolve only with globally valid constraints */
               if(!SCIPconsIsGlobal(reasoncon))
               {
                  conflict->ncorrectaborts++;
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

            /* get next bound change from queue */
            nextbdchginfo = conflictFirstCand(conflict);

            /* if we still have not resolved and the current bound change is a
            branching decision we can continue */
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
         }
         /* if no bound change was infered by a resolvable constraint then we terminate */
         else
         {
            conflict->ncorrectaborts++;
            goto TERMINATE;
         }
      }
      else if( !bdchginfoIsResolvable(bdchginfo) && existsResolvablebdchginfo(conflict) )
      {
         /* todo this should be written as a function since it is already used above */
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

         /* get next bound change from queue */
         nextbdchginfo = conflictFirstCand(conflict);

         /* if no resolution has been applied yet, and the bound change is a
         branching decision, we can ignore it and continue with the next
         bound change */
         if( nressteps == 0 && bdchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
         {
            bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);
         }
         else
         {
            conflict->ncorrectaborts++;
            goto TERMINATE;
         }
         bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
         vartoresolve = SCIPbdchginfoGetVar(bdchginfo);

         SCIPdebug(resolutionsetPrintRow(conflictresolutionset, set, transprob, 1));

         assert(SCIPsetIsLT(set, getSlack(set, transprob, conflictresolutionset, bdchgidx, fixbounds, fixinds), 0.0));
         /* check if the variable we are resolving is active */
         assert(SCIPvarIsActive(vartoresolve));
      }
      else
      {
         SCIP_Bool obtainedreason;
         residx = SCIPvarGetProbindex(vartoresolve);

         /* get reason row of the latest bdchginfo */
         SCIP_CALL( getReasonRow(conflict, blkmem, transprob, set, bdchginfo, residx, validdepth, fixbounds, fixinds, &obtainedreason) );
         if( !obtainedreason )
         {
            conflict->ncorrectaborts++;
            SCIPsetDebugMsgPrint(set, "Could not obtain reason row for bound change \n");
            goto TERMINATE;
         }
         reasonslack = reasonresolutionset->slack;

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

         /** todo add the following: if the slack of the reason is negative, we
          * can restart conflict analysis with the reason as conflict constraint
          * (this may happen in depending on the order of propagation)
          */

         if ( !set->conf_weakenreason || set->conf_reductiontechnique == 'o')
         {
            SCIPsetDebugMsg(set, " Applying resolution to remove variable <%s>\n", SCIPvarGetName(vartoresolve));

            SCIP_CALL( rescaleAndResolve(set, conflictresolutionset, reasonresolutionset, resolvedresolutionset, blkmem,
                                 residx, &successresolution) );
            if (!successresolution)
            {
               conflict->ncorrectaborts++;
               goto TERMINATE;
            }
         }
         else
         {
            int nvarsweakened;
            SCIPsetDebugMsg(set, " Applying tightening based reduction with resolving variable <%s>\n", SCIPvarGetName(vartoresolve));

            if ( set->conf_reductiontechnique == 'c' )
               SCIP_CALL( tighteningBasedReduction(conflict, set, transprob, blkmem, bdchgidx, residx, &nvarsweakened, fixbounds, fixinds, &successresolution ) );
             else if ( set->conf_reductiontechnique == 'd' || set->conf_reductiontechnique == 'm' )
               SCIP_CALL(DivisionBasedReduction(conflict, set, transprob, blkmem, bdchgidx, residx, &nvarsweakened, fixbounds, fixinds, &successresolution ) );
         }

         SCIP_CALL( resolutionsetReplace(conflictresolutionset, blkmem, resolvedresolutionset) );

         SCIPstatisticPrintf("ConflictSTAT: %d %d %f %f\n", nressteps, resolutionsetGetNNzs(conflictresolutionset),
                             conflictresolutionset->coefquotient, conflictresolutionset->slack);


         conflictslack = getSlack(set, transprob, conflictresolutionset, bdchgidx, fixbounds, fixinds);
         SCIPdebug(resolutionsetPrintRow(conflictresolutionset, set, transprob, 3));
         conflictresolutionset->slack = conflictslack;

         SCIPsetDebugMsg(set, "Slack of resolved row: %f \n", conflictslack);

         /** Unfortunately we cannot guarrante that the slack becomes zero after reducing the reason (even if we have only binary variables)
          *  TIll now there are two major problems:
          *    - Knapsack constraints that use negated cliques in the propagation
          *    - Ranged row propagation (gcd argument)
          */

         /* check that we fail for a valid reason */
         if (SCIPsetIsGE(set, conflictslack, 0.0))
         {
            /* todo either remove the member isbinary in the resolution sets or update it */
            if ( set->conf_reductiontechnique == 'o' || !isBinaryResolutionSet(set, transprob, conflictresolutionset))
               conflict->ncorrectaborts++;
            goto TERMINATE;
         }

         nressteps++;

         /* apply coefficient tightening to the resolved constraint should never hurt */
         SCIP_CALL( tightenCoefLhs(set, transprob, FALSE, conflictresolutionset->vals, &conflictresolutionset->lhs, conflictresolutionset->inds, &conflictresolutionset->nnz, &nchgcoefs, NULL) );
         if (nchgcoefs > 0)
         {
            SCIP_Real previousslack;
            previousslack = conflictresolutionset->slack;
            conflictresolutionset->slack = getSlack(set, transprob, conflictresolutionset, bdchgidx, fixbounds, fixinds);
            SCIPsetDebugMsg(set, "Tightened %d coefficients in the resolved constraint, old slack %f, new slack %f \n", nchgcoefs, previousslack,conflictresolutionset->slack);
            assert(SCIPsetIsLE(set, conflictresolutionset->slack, previousslack + EPS));
         }
         /* sort for linear time resolution */
         /* todo this can be removed if we do not apply coefficient tightening */
         SCIPsortIntReal(conflictresolutionset->inds, conflictresolutionset->vals, resolutionsetGetNNzs(conflictresolutionset));
         if (SCIPsetIsLT(set, conflictresolutionset->slack, 0.0))
         {
            /* Apply cmir after each iteration to strengthen the conflict constraint */
            if( set->conf_applycmir || set->conf_applysimplemir )
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
         /* By default conf_maxnumressteps is -1 -> we do not stop early */
         if( set->conf_maxnumressteps > 0 && nressteps >= set->conf_maxnumressteps )
         {
            goto TERMINATE;
         }


         /** clean up the queue of bound changes, e.g.
          * if variables get canceled during resolution
          */
         cleanBdchgQueue(conflict, set, conflictresolutionset);

         SCIP_CALL( updateBdchgQueue(set, transprob, conflictresolutionset, bdchgidx) );

         /* get the next bound change */
         bdchginfo = conflictFirstCand(conflict);

         /* if no bound change exists the we should stop */
         /* in case we have applied resolution steps we keep the last conflict constraint */
         if( bdchginfo == NULL )
         {
            /* this can happen only if we already have resolved and some
            and it means that we have already reached a FUIP */
            SCIPsetDebugMsgPrint(set, " reached UIP in depth %d \n", bdchgdepth);
            /* add the previous conflict in the list of resolution sets */
            prevconflictresolutionset->conflictdepth = bdchgdepth;
            prevconflictresolutionset->repropdepth = (nextbdchginfo == NULL) ? 0 : SCIPbdchginfoGetDepth(nextbdchginfo);
            SCIP_RESOLUTIONSET* tmpconflictresolutionset;
            SCIP_CALL( resolutionsetCopy(&tmpconflictresolutionset, blkmem, prevconflictresolutionset) );
            SCIP_CALL( conflictInsertResolutionset(conflict, set, &tmpconflictresolutionset) );
            goto TERMINATE;
         }

         bdchginfo = conflictRemoveCand(conflict);
         bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
         vartoresolve = SCIPbdchginfoGetVar(bdchginfo);
         bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);

         residx = SCIPvarGetProbindex(vartoresolve);

         assert(getVarIdxInResolutionset(conflictresolutionset, residx) >= 0);
         /* get the bound change before bdchginfo */
         nextbdchginfo = conflictFirstCand(conflict);

         /* check if the variable we are resolving is active */
         assert(SCIPvarIsActive(vartoresolve));

         /* when at an UIP add the previous conflict in the list of resolution sets */
         if( nextbdchginfo == NULL || SCIPbdchginfoGetDepth(nextbdchginfo) != bdchgdepth  )
         {
            assert( nresstepslast != nressteps );
            SCIPsetDebugMsgPrint(set, " reached UIP in depth %d \n", bdchgdepth);
            /* add the previous conflict in the list of resolution sets */
            prevconflictresolutionset->conflictdepth = bdchgdepth;
            prevconflictresolutionset->repropdepth = (nextbdchginfo == NULL) ? 0 : SCIPbdchginfoGetDepth(nextbdchginfo);
            SCIP_RESOLUTIONSET* tmpconflictresolutionset;
            SCIP_CALL( resolutionsetCopy(&tmpconflictresolutionset, blkmem, prevconflictresolutionset) );
            SCIP_CALL( conflictInsertResolutionset(conflict, set, &tmpconflictresolutionset) );
            nresstepslast = nressteps;
            nfuips ++;
            /* stop after conf_resfuiplevels UIPs */
            if (set->conf_resfuiplevels > 0 && nfuips >= set->conf_resfuiplevels)
               goto TERMINATE;
         }
      }

   }

  TERMINATE:

   /* if resolution fails at some point we can still return add the latest valid
   conflict in the list of resolution set */
   if ( nressteps >= 1 && nresstepslast != nressteps )
   {
      SCIP_RESOLUTIONSET* tmpconflictresolutionset;
      SCIP_CALL( resolutionsetCopy(&tmpconflictresolutionset, blkmem, prevconflictresolutionset) );
      SCIP_CALL( conflictInsertResolutionset(conflict, set, &tmpconflictresolutionset) );
   }

   SCIPstatisticPrintf("End Statistics \n");
   SCIPsetDebugMsg(set, "Total number of resolution sets found %d\n", conflict->nresolutionsets);
   SCIPsetDebugMsg(set, "Total number of FUIPS found %d\n", nfuips);

   if ( conflict->nresolutionsets > 0 )
   {
      // int nconstoadd;

      /** add the conflict constraints to the current node
       * (all operations are done only on global constraints -> add conflict constraints to the root node)
       */
      // nconstoadd = (set->conf_resolutioncons > 0) ? MIN(set->conf_resolutioncons, conflict->nresolutionsets) : conflict->nresolutionsets;

      for( i = 0; i < nfuips; i++ )
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

   SCIPsetDebugMsg(set, "Starting resolution based conflict analysis after infeasible propagation in depth %d\n",
                   SCIPtreeGetCurrentDepth(tree));

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
   SCIPpqueueClear(conflict->separatebdchgqueue);

   /* stop timing */
   SCIPclockStop(conflict->resanalyzetime, set);
   SCIPsetDebugMsg(set, "resolution based conflict analysis added %d constraints \n \n", nconss);

   return SCIP_OKAY;
}
