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

/**@file   conflict_dualproofanalysis.c
 * @ingroup OTHER_CFILES
 * @brief  internal methods for dual proof conflict analysis
 * @author Timo Berthold
 * @author Jakob Witzig
 *
 * In dual proof analysis, an infeasible LP relaxation is analysed.
 * Using the dual solution, a valid constraint is derived that is violated
 * by all values in the domain. This constraint is added to the problem
 * and can then be used for domain propagation. More details can be found in [1]
 *
 * [1] J. Witzig, T. Berthold, en S. Heinz, ‘Computational aspects of infeasibility analysis in mixed integer programming’,
 * Math. Prog. Comp., mrt. 2021, doi: 10.1007/s12532-021-00202-0.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "lpi/lpi.h"
#include "scip/clock.h"
#include "scip/conflict_general.h"
#include "scip/conflict_dualproofanalysis.h"
#include "scip/conflictstore.h"
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

#define BOUNDSWITCH                0.51 /**< threshold for bound switching - see cuts.c */
#define POSTPROCESS               FALSE /**< apply postprocessing to the cut - see cuts.c */
#define USEVBDS                   FALSE /**< use variable bounds - see cuts.c */
#define ALLOWLOCAL                FALSE /**< allow to generate local cuts - see cuts. */
#define MINFRAC                   0.05  /**< minimal fractionality of floor(rhs) - see cuts.c */
#define MAXFRAC                   0.999 /**< maximal fractionality of floor(rhs) - see cuts.c */


/*
 * Proof Sets
 */

/** return the char associated with the type of the variable */
static
char varGetChar(
   SCIP_VAR*             var                 /**< variable */
   )
{
   SCIP_VARTYPE vartype = SCIPvarGetType(var);

   return (!SCIPvarIsIntegral(var) ? 'C' :
          (vartype == SCIP_VARTYPE_BINARY ? 'B' :
          (vartype == SCIP_VARTYPE_INTEGER ? 'I' : 'M')));
}

/** resets the data structure of a proofset */
static
void proofsetClear(
   SCIP_PROOFSET*        proofset            /**< proof set */
   )
{
   assert(proofset != NULL);

   proofset->nnz = 0;
   proofset->rhs = 0.0;
   proofset->validdepth = 0;
   proofset->conflicttype = SCIP_CONFTYPE_UNKNOWN;
}

/** creates a proofset */
static
SCIP_RETCODE proofsetCreate(
   SCIP_PROOFSET**       proofset,           /**< proof set */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(proofset != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, proofset) );
   (*proofset)->vals = NULL;
   (*proofset)->inds = NULL;
   (*proofset)->rhs = 0.0;
   (*proofset)->nnz = 0;
   (*proofset)->size = 0;
   (*proofset)->validdepth = 0;
   (*proofset)->conflicttype = SCIP_CONFTYPE_UNKNOWN;

   return SCIP_OKAY;
}

/** frees a proofset */
void SCIPproofsetFree(
   SCIP_PROOFSET**       proofset,           /**< proof set */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(proofset != NULL);
   assert(*proofset != NULL);
   assert(blkmem != NULL);

   BMSfreeBlockMemoryArrayNull(blkmem, &(*proofset)->vals, (*proofset)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*proofset)->inds, (*proofset)->size);
   BMSfreeBlockMemory(blkmem, proofset);
   (*proofset) = NULL;
}

#ifdef SCIP_DEBUG
/** print a proof set */
void proofsetPrint(
   SCIP_PROOFSET*        proofset,
   SCIP_SET*             set,
   SCIP_PROB*            transprob
   )
{
   SCIP_VAR** vars;
   int i;

   assert(proofset != NULL);

   vars = SCIPprobGetVars(transprob);
   assert(vars != NULL);

   printf("proofset: ");
   for( i = 0; i < proofset->nnz; i++ )
      printf("%+.15g <%s> ", proofset->vals[i], SCIPvarGetName(vars[proofset->inds[i]]));
   printf(" <= %.15g\n", proofset->rhs);
}
#endif

/** return the indices of variables in the proofset */
static
int* proofsetGetInds(
   SCIP_PROOFSET*        proofset            /**< proof set */
   )
{
   assert(proofset != NULL);

   return proofset->inds;
}

/** return coefficient of variable in the proofset with given probindex */
static
SCIP_Real* proofsetGetVals(
   SCIP_PROOFSET*        proofset            /**< proof set */
   )
{
   assert(proofset != NULL);

   return proofset->vals;
}

/** return the right-hand side if a proofset */
static
SCIP_Real proofsetGetRhs(
   SCIP_PROOFSET*        proofset            /**< proof set */
   )
{
   assert(proofset != NULL);

   return proofset->rhs;
}

/** returns the number of variables in the proofset */
int SCIPproofsetGetNVars(
   SCIP_PROOFSET*        proofset            /**< proof set */
   )
{
   assert(proofset != NULL);

   return proofset->nnz;
}

/** returns the number of variables in the proofset */
static
SCIP_CONFTYPE proofsetGetConftype(
   SCIP_PROOFSET*        proofset            /**< proof set */
   )
{
   assert(proofset != NULL);

   return proofset->conflicttype;
}

/** adds given data as aggregation row to the proofset */
static
SCIP_RETCODE proofsetAddSparseData(
   SCIP_PROOFSET*        proofset,           /**< proof set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Real*            vals,               /**< variable coefficients */
   int*                  inds,               /**< variable array */
   int                   nnz,                /**< size of variable and coefficient array */
   SCIP_Real             rhs                 /**< right-hand side of the aggregation row */
   )
{
   assert(proofset != NULL);
   assert(blkmem != NULL);

   if( proofset->size == 0 )
   {
      assert(proofset->vals == NULL);
      assert(proofset->inds == NULL);

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &proofset->vals, vals, nnz) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &proofset->inds, inds, nnz) );

      proofset->size = nnz;
   }
   else
   {
      int i;

      assert(proofset->vals != NULL);
      assert(proofset->inds != NULL);

      if( proofset->size < nnz )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &proofset->vals, proofset->size, nnz) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &proofset->inds, proofset->size, nnz) );
         proofset->size = nnz;
      }

      for( i = 0; i < nnz; i++ )
      {
         proofset->vals[i] = vals[i];
         proofset->inds[i] = inds[i];
      }
   }

   proofset->rhs = rhs;
   proofset->nnz = nnz;

   return SCIP_OKAY;
}

/** adds an aggregation row to the proofset */
static
SCIP_RETCODE proofsetAddAggrrow(
   SCIP_PROOFSET*        proofset,           /**< proof set */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_AGGRROW*         aggrrow             /**< aggregation row to add */
   )
{
   SCIP_Real* vals;
   int* inds;
   int nnz;
   int i;

   assert(proofset != NULL);
   assert(set != NULL);

   inds = SCIPaggrRowGetInds(aggrrow);
   assert(inds != NULL);

   nnz = SCIPaggrRowGetNNz(aggrrow);
   assert(nnz > 0);

   SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, nnz) );

   for( i = 0; i < nnz; i++ )
   {
      vals[i] = SCIPaggrRowGetProbvarValue(aggrrow, inds[i]);
   }

   SCIP_CALL( proofsetAddSparseData(proofset, blkmem, vals, inds, nnz, SCIPaggrRowGetRhs(aggrrow)) );

   SCIPsetFreeBufferArray(set, &vals);

   return SCIP_OKAY;
}

/** Removes a given variable @p var from position @p pos from the proofset and updates the right-hand side according
 *  to sign of the coefficient, i.e., rhs -= coef * bound, where bound = lb if coef >= 0 and bound = ub, otherwise.
 *
 *  @note: The list of non-zero indices and coefficients will be updated by swapping the last non-zero index to @p pos.
 */
static
void proofsetCancelVarWithBound(
   SCIP_PROOFSET*        proofset,
   SCIP_SET*             set,
   SCIP_VAR*             var,
   int                   pos,
   SCIP_Bool*            valid
   )
{
   assert(proofset != NULL);
   assert(var != NULL);
   assert(pos >= 0 && pos < proofset->nnz);
   assert(valid != NULL);

   *valid = TRUE;

   /* cancel with lower bound */
   if( proofset->vals[pos] > 0.0 )
   {
      proofset->rhs -= proofset->vals[pos] * SCIPvarGetLbGlobal(var);
   }
   /* cancel with upper bound */
   else
   {
      assert(proofset->vals[pos] < 0.0);
      proofset->rhs -= proofset->vals[pos] * SCIPvarGetUbGlobal(var);
   }

   --proofset->nnz;

   proofset->vals[pos] = proofset->vals[proofset->nnz];
   proofset->inds[pos] = proofset->inds[proofset->nnz];
   proofset->vals[proofset->nnz] = 0.0;
   proofset->inds[proofset->nnz] = 0;

   if( SCIPsetIsInfinity(set, proofset->rhs) )
      *valid = FALSE;
}

/** creates and clears the proofset */
SCIP_RETCODE SCIPconflictInitProofset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(conflict != NULL);
   assert(blkmem != NULL);

   SCIP_CALL( proofsetCreate(&conflict->proofset, blkmem) );

   return SCIP_OKAY;
}

/** resizes proofsets array to be able to store at least num entries */
static
SCIP_RETCODE conflictEnsureProofsetsMem(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);

   if( num > conflict->proofsetssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conflict->proofsets, newsize) );
      conflict->proofsetssize = newsize;
   }
   assert(num <= conflict->proofsetssize);

   return SCIP_OKAY;
}

/** add a proofset to the list of all proofsets */
static
SCIP_RETCODE conflictInsertProofset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROOFSET*        proofset            /**< proof set to add */
   )
{
   assert(conflict != NULL);
   assert(proofset != NULL);

   /* insert proofset into the sorted proofsets array */
   SCIP_CALL( conflictEnsureProofsetsMem(conflict, set, conflict->nproofsets + 1) );

   conflict->proofsets[conflict->nproofsets] = proofset;
   ++conflict->nproofsets;

   return SCIP_OKAY;
}

/** tighten the bound of a singleton variable in a constraint
 *
 *  if the bound is contradicting with a global bound we cannot tighten the bound directly.
 *  in this case we need to create and add a constraint of size one such that propagating this constraint will
 *  enforce the infeasibility.
 */
static
SCIP_RETCODE tightenSingleVar(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_TREE*            tree,               /**< tree data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             val,                /**< coefficient of the variable */
   SCIP_Real             rhs,                /**< rhs of the constraint */
   SCIP_CONFTYPE         prooftype,          /**< type of the proof */
   int                   validdepth          /**< depth where the bound change is valid */
   )
{
   SCIP_Real newbound;
   SCIP_Bool applyglobal;
   SCIP_BOUNDTYPE boundtype;

   assert(tree != NULL);
   assert(validdepth >= 0);

   applyglobal = (validdepth <= SCIPtreeGetEffectiveRootDepth(tree));

   /* if variable and coefficient are integral the rhs can be rounded down */
   if( SCIPvarIsIntegral(var) && SCIPsetIsIntegral(set, val) )
      newbound = SCIPsetFeasFloor(set, rhs)/val;
   else
      newbound = rhs/val;

   boundtype = (val > 0.0 ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
   SCIPvarAdjustBd(var, set, boundtype, &newbound);

   /* skip numerical unstable bound changes */
   if( applyglobal
      && ((boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsLE(set, newbound, SCIPvarGetLbGlobal(var)))
       || (boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsGE(set, newbound, SCIPvarGetUbGlobal(var)))) )
   {
      return SCIP_OKAY;
   }

   /* the new bound contradicts a global bound, we can cutoff the root node immediately */
   if( applyglobal
      && ((boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsGT(set, newbound, SCIPvarGetUbGlobal(var)))
       || (boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsLT(set, newbound, SCIPvarGetLbGlobal(var)))) )
   {
      SCIPsetDebugMsg(set, "detected global infeasibility at var <%s>: locdom=[%g,%g] glbdom=[%g,%g] new %s bound=%g\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var),
            SCIPvarGetUbLocal(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var),
            (boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper"), newbound);
      SCIP_CALL( SCIPnodeCutoff(tree->path[0], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
   }
   else
   {
      if( lp->strongbranching || !applyglobal )
      {
         SCIP_CONS* cons;
         SCIP_Real conslhs;
         SCIP_Real consrhs;
         char name[SCIP_MAXSTRLEN];

         SCIPsetDebugMsg(set, "add constraint <%s>[%c] %s %g to node #%lld in depth %d\n",
               SCIPvarGetName(var), varGetChar(var), boundtype == SCIP_BOUNDTYPE_UPPER ? "<=" : ">=", newbound,
               SCIPnodeGetNumber(tree->path[validdepth]), validdepth);

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "pc_fix_%s", SCIPvarGetName(var));

         if( boundtype == SCIP_BOUNDTYPE_UPPER )
         {
            conslhs = -SCIPsetInfinity(set);
            consrhs = newbound;
         }
         else
         {
            conslhs = newbound;
            consrhs = SCIPsetInfinity(set);
         }

         SCIP_CALL( SCIPcreateConsLinear(set->scip, &cons, name, 0, NULL, NULL, conslhs, consrhs,
               FALSE, FALSE, FALSE, FALSE, TRUE, !applyglobal, FALSE, TRUE, TRUE, FALSE) );

         SCIP_CALL( SCIPaddCoefLinear(set->scip, cons, var, 1.0) );

         if( applyglobal )
         {
            SCIP_CALL( SCIPprobAddCons(transprob, set, stat, cons) );
         }
         else
         {
            SCIP_CALL( SCIPnodeAddCons(tree->path[validdepth], blkmem, set, stat, tree, cons) );
         }

         SCIP_CALL( SCIPconsRelease(&cons, blkmem, set) );
      }
      else
      {
         assert(applyglobal);

         SCIPsetDebugMsg(set, "change global %s bound of <%s>[%c]: %g -> %g\n",
               (boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper"),
               SCIPvarGetName(var), varGetChar(var),
               (boundtype == SCIP_BOUNDTYPE_LOWER ? SCIPvarGetLbGlobal(var) : SCIPvarGetUbGlobal(var)),
               newbound);

         SCIP_CALL( SCIPnodeAddBoundchg(tree->path[0], blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, \
               eventqueue, cliquetable, var, newbound, boundtype, FALSE) );

         /* mark the node in the validdepth to be propagated again */
         SCIPnodePropagateAgain(tree->path[0], set, stat, tree);
      }
   }

   if( applyglobal )
      ++conflict->nglbchgbds;
   else
      ++conflict->nlocchgbds;

   if( prooftype == SCIP_CONFTYPE_INFEASLP || prooftype == SCIP_CONFTYPE_ALTINFPROOF )
   {
      ++conflict->dualproofsinfnnonzeros; /* we count a global bound reduction as size 1 */
      ++conflict->ndualproofsinfsuccess;
      ++conflict->ninflpsuccess;

      if( applyglobal )
         ++conflict->ndualproofsinfglobal;
      else
         ++conflict->ndualproofsinflocal;
   }
   else
   {
      ++conflict->dualproofsbndnnonzeros; /* we count a global bound reduction as size 1 */
      ++conflict->ndualproofsbndsuccess;
      ++conflict->nboundlpsuccess;

      if( applyglobal )
         ++conflict->ndualproofsbndglobal;
      else
         ++conflict->ndualproofsbndlocal;
   }

   return SCIP_OKAY;
}

/** calculates the minimal activity of a given set of bounds and coefficients */
static
SCIP_Real getMinActivity(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_Real*            coefs,              /**< coefficients in sparse representation */
   int*                  inds,               /**< non-zero indices */
   int                   nnz,                /**< number of non-zero indices */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables (or NULL for global bounds) */
   SCIP_Real*            curvarubs           /**< current upper bounds of active problem variables (or NULL for global bounds) */
   )
{
   SCIP_VAR** vars;
   SCIP_Real QUAD(minact);
   int i;

   assert(coefs != NULL);
   assert(inds != NULL);

   vars = SCIPprobGetVars(transprob);
   assert(vars != NULL);

   QUAD_ASSIGN(minact, 0.0);

   for( i = 0; i < nnz; i++ )
   {
      SCIP_Real val;
      SCIP_Real QUAD(delta);
      int v = inds[i];

      assert(SCIPvarGetProbindex(vars[v]) == v);

      val = coefs[i];

      if( val > 0.0 )
      {
         SCIP_Real bnd;

         assert(curvarlbs == NULL || !SCIPsetIsInfinity(set, -curvarlbs[v]));

         bnd = (curvarlbs == NULL ? SCIPvarGetLbGlobal(vars[v]) : curvarlbs[v]);

         if( SCIPsetIsInfinity(set, -bnd) )
            return -SCIPsetInfinity(set);

         SCIPquadprecProdDD(delta, val, bnd);
      }
      else
      {
         SCIP_Real bnd;

         assert(curvarubs == NULL || !SCIPsetIsInfinity(set, curvarubs[v]));

         bnd = (curvarubs == NULL ? SCIPvarGetUbGlobal(vars[v]) : curvarubs[v]);

         if( SCIPsetIsInfinity(set, bnd) )
            return -SCIPsetInfinity(set);

         SCIPquadprecProdDD(delta, val, bnd);
      }

      /* update minimal activity */
      SCIPquadprecSumQQ(minact, minact, delta);
   }

   /* check whether the minmal activity is infinite */
   if( SCIPsetIsInfinity(set, QUAD_TO_DBL(minact)) )
      return SCIPsetInfinity(set);
   if( SCIPsetIsInfinity(set, -QUAD_TO_DBL(minact)) )
      return -SCIPsetInfinity(set);

   return QUAD_TO_DBL(minact);
}

/** calculates the minimal activity of a given set of bounds and coefficients */
static
SCIP_Real getMaxActivity(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_Real*            coefs,              /**< coefficients in sparse representation */
   int*                  inds,               /**< non-zero indices */
   int                   nnz,                /**< number of non-zero indices */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables (or NULL for global bounds) */
   SCIP_Real*            curvarubs           /**< current upper bounds of active problem variables (or NULL for global bounds) */
   )
{
   SCIP_VAR** vars;
   SCIP_Real QUAD(maxact);
   int i;

   assert(coefs != NULL);
   assert(inds != NULL);

   vars = SCIPprobGetVars(transprob);
   assert(vars != NULL);

   QUAD_ASSIGN(maxact, 0.0);

   for( i = 0; i < nnz; i++ )
   {
      SCIP_Real val;
      SCIP_Real QUAD(delta);
      int v = inds[i];

      assert(SCIPvarGetProbindex(vars[v]) == v);

      val = coefs[i];

      if( val < 0.0 )
      {
         SCIP_Real bnd;

         assert(curvarlbs == NULL || !SCIPsetIsInfinity(set, -curvarlbs[v]));

         bnd = (curvarlbs == NULL ? SCIPvarGetLbGlobal(vars[v]) : curvarlbs[v]);

         if( SCIPsetIsInfinity(set, -bnd) )
            return SCIPsetInfinity(set);

         SCIPquadprecProdDD(delta, val, bnd);
      }
      else
      {
         SCIP_Real bnd;

         assert(curvarubs == NULL || !SCIPsetIsInfinity(set, curvarubs[v]));

         bnd = (curvarubs == NULL ? SCIPvarGetUbGlobal(vars[v]) : curvarubs[v]);

         if( SCIPsetIsInfinity(set, bnd) )
            return SCIPsetInfinity(set);

         SCIPquadprecProdDD(delta, val, bnd);
      }

      /* update maximal activity */
      SCIPquadprecSumQQ(maxact, maxact, delta);
   }

   /* check whether the maximal activity got infinite */
   if( SCIPsetIsInfinity(set, QUAD_TO_DBL(maxact)) )
      return SCIPsetInfinity(set);
   if( SCIPsetIsInfinity(set, -QUAD_TO_DBL(maxact)) )
      return -SCIPsetInfinity(set);

   return QUAD_TO_DBL(maxact);
}

/** propagate a long proof */
static
SCIP_RETCODE propagateLongProof(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_TREE*            tree,               /**< tree data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real*            coefs,              /**< coefficients in sparse representation */
   int*                  inds,               /**< non-zero indices */
   int                   nnz,                /**< number of non-zero indices */
   SCIP_Real             rhs,                /**< right-hand side */
   SCIP_CONFTYPE         conflicttype,       /**< type of the conflict */
   int                   validdepth          /**< depth where the proof is valid */
   )
{
   SCIP_VAR** vars;
   SCIP_Real minact;
   int i;

   assert(coefs != NULL);
   assert(inds != NULL);
   assert(nnz >= 0);

   vars = SCIPprobGetVars(transprob);
   minact = getMinActivity(set, transprob, coefs, inds, nnz, NULL, NULL);

   /* we cannot find global tightenings */
   if( SCIPsetIsInfinity(set, -minact) )
      return SCIP_OKAY;

   for( i = 0; i < nnz; i++ )
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_Real resminact;
      SCIP_Real lb;
      SCIP_Real ub;
      int pos;

      pos = inds[i];
      val = coefs[i];
      var = vars[pos];
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      assert(!SCIPsetIsZero(set, val));

      resminact = minact;

      /* we got a potential new upper bound */
      if( val > 0.0 )
      {
         SCIP_Real newub;

         resminact -= (val * lb);
         newub = (rhs - resminact)/val;

         if( SCIPsetIsInfinity(set, newub) )
            continue;

         /* we cannot tighten the upper bound */
         if( SCIPsetIsGE(set, newub, ub) )
            continue;

         SCIP_CALL( tightenSingleVar(conflict, set, stat, tree, blkmem, origprob, transprob, reopt, lp, branchcand, \
               eventqueue, cliquetable, var, val, rhs-resminact, conflicttype, validdepth) );
      }
      /* we got a potential new lower bound */
      else
      {
         SCIP_Real newlb;

         resminact -= (val * ub);
         newlb = (rhs - resminact)/val;

         if( SCIPsetIsInfinity(set, -newlb) )
            continue;

         /* we cannot tighten the lower bound */
         if( SCIPsetIsLE(set, newlb, lb) )
            continue;

         SCIP_CALL( tightenSingleVar(conflict, set, stat, tree, blkmem, origprob, transprob, reopt, lp, branchcand, \
               eventqueue, cliquetable, var, val, rhs-resminact, conflicttype, validdepth) );
      }

      /* the minimal activity should stay unchanged because we tightened the bound that doesn't contribute to the
       * minimal activity
       */
      assert(SCIPsetIsEQ(set, minact, getMinActivity(set, transprob, coefs, inds, nnz, NULL, NULL)));
   }

   return SCIP_OKAY;
}

/** creates a proof constraint and tries to add it to the storage */
static
SCIP_RETCODE createAndAddProofcons(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict pool data */
   SCIP_PROOFSET*        proofset,           /**< proof set */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< tree data */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_CONS* cons;
   SCIP_CONS* upgdcons;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int* inds;
   SCIP_Real rhs;
   SCIP_Real fillin;
   SCIP_Real globalminactivity;
   SCIP_Bool applyglobal;
   SCIP_Bool toolong;
   SCIP_Bool contonly;
   SCIP_Bool hasrelaxvar;
   SCIP_CONFTYPE conflicttype;
   char name[SCIP_MAXSTRLEN];
   int nnz;
   int i;

   assert(conflict != NULL);
   assert(conflictstore != NULL);
   assert(proofset != NULL);
   assert(proofset->validdepth == 0 || proofset->validdepth < SCIPtreeGetFocusDepth(tree));

   nnz = SCIPproofsetGetNVars(proofset);

   if( nnz == 0 )
      return SCIP_OKAY;

   vars = SCIPprobGetVars(transprob);

   rhs = proofsetGetRhs(proofset);
   assert(!SCIPsetIsInfinity(set, rhs));

   coefs = proofsetGetVals(proofset);
   assert(coefs != NULL);

   inds = proofsetGetInds(proofset);
   assert(inds != NULL);

   conflicttype = proofsetGetConftype(proofset);

   applyglobal = (proofset->validdepth <= SCIPtreeGetEffectiveRootDepth(tree));

   if( applyglobal )
   {
      SCIP_Real globalmaxactivity = getMaxActivity(set, transprob, coefs, inds, nnz, NULL, NULL);

      /* check whether the alternative proof is redundant */
      if( SCIPsetIsLE(set, globalmaxactivity, rhs) )
         return SCIP_OKAY;

      /* check whether the constraint proves global infeasibility */
      globalminactivity = getMinActivity(set, transprob, coefs, inds, nnz, NULL, NULL);
      if( SCIPsetIsGT(set, globalminactivity, rhs) )
      {
         SCIPsetDebugMsg(set, "detect global infeasibility: minactivity=%g, rhs=%g\n", globalminactivity, rhs);

         SCIP_CALL( SCIPnodeCutoff(tree->path[proofset->validdepth], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );

         goto UPDATESTATISTICS;
      }
   }

   if( set->conf_minmaxvars >= nnz  )
      toolong = FALSE;
   else
   {
      SCIP_Real maxnnz;

      if( transprob->startnconss < 100 )
         maxnnz = 0.85 * transprob->nvars;
      else
         maxnnz = (SCIP_Real)transprob->nvars;

      fillin = nnz;
      if( conflicttype == SCIP_CONFTYPE_INFEASLP || conflicttype == SCIP_CONFTYPE_ALTINFPROOF )
      {
         fillin += SCIPconflictstoreGetNDualInfProofs(conflictstore) * SCIPconflictstoreGetAvgNnzDualInfProofs(conflictstore);
         fillin /= (SCIPconflictstoreGetNDualInfProofs(conflictstore) + 1.0);
         toolong = (fillin > MIN(2.0 * stat->avgnnz, maxnnz));
      }
      else
      {
         assert(conflicttype == SCIP_CONFTYPE_BNDEXCEEDING || conflicttype == SCIP_CONFTYPE_ALTBNDPROOF);

         fillin += SCIPconflictstoreGetNDualBndProofs(conflictstore) * SCIPconflictstoreGetAvgNnzDualBndProofs(conflictstore);
         fillin /= (SCIPconflictstoreGetNDualBndProofs(conflictstore) + 1.0);
         toolong = (fillin > MIN(1.5 * stat->avgnnz, maxnnz));
      }

      toolong = (toolong && (nnz > set->conf_maxvarsfac * transprob->nvars));
   }

   /* don't store global dual proofs that are too long / have too many non-zeros */
   if( toolong )
   {
      if( applyglobal )
      {
         SCIP_CALL( propagateLongProof(conflict, set, stat, reopt, tree, blkmem, origprob, transprob, lp, branchcand,
               eventqueue, cliquetable, coefs, inds, nnz, rhs, conflicttype, proofset->validdepth) );
      }
      return SCIP_OKAY;
   }

   /* check if conflict contains variables that are invalid after a restart to label it appropriately */
   hasrelaxvar = FALSE;
   contonly = TRUE;
   for( i = 0; i < nnz && (!hasrelaxvar || contonly); ++i )
   {
      hasrelaxvar |= SCIPvarIsRelaxationOnly(vars[inds[i]]);

      if( SCIPvarIsIntegral(vars[inds[i]]) )
         contonly = FALSE;
   }

   if( !applyglobal && contonly )
      return SCIP_OKAY;

   if( conflicttype == SCIP_CONFTYPE_INFEASLP || conflicttype == SCIP_CONFTYPE_ALTINFPROOF )
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "dualproof_inf_%" SCIP_LONGINT_FORMAT, conflict->ndualproofsinfsuccess);
   else if( conflicttype == SCIP_CONFTYPE_BNDEXCEEDING || conflicttype == SCIP_CONFTYPE_ALTBNDPROOF )
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "dualproof_bnd_%" SCIP_LONGINT_FORMAT, conflict->ndualproofsbndsuccess);
   else
      return SCIP_INVALIDCALL;

   SCIP_CALL( SCIPcreateConsLinear(set->scip, &cons, name, 0, NULL, NULL, -SCIPsetInfinity(set), rhs,
         FALSE, FALSE, FALSE, FALSE, TRUE, !applyglobal,
         FALSE, TRUE, TRUE, FALSE) );

   for( i = 0; i < nnz; i++ )
   {
      int v = inds[i];
      SCIP_CALL( SCIPaddCoefLinear(set->scip, cons, vars[v], coefs[i]) );
   }

   /* do not upgrade linear constraints of size 1 */
   if( nnz > 1 )
   {
      upgdcons = NULL;
      /* try to automatically convert a linear constraint into a more specific and more specialized constraint */
      SCIP_CALL( SCIPupgradeConsLinear(set->scip, cons, &upgdcons) );
      if( upgdcons != NULL )
      {
         SCIP_CALL( SCIPreleaseCons(set->scip, &cons) );
         cons = upgdcons;

         if( conflicttype == SCIP_CONFTYPE_INFEASLP )
            conflicttype = SCIP_CONFTYPE_ALTINFPROOF;
         else if( conflicttype == SCIP_CONFTYPE_BNDEXCEEDING )
            conflicttype = SCIP_CONFTYPE_ALTBNDPROOF;
      }
   }

   /* mark constraint to be a conflict */
   SCIPconsMarkConflict(cons);

   /* add constraint to storage */
   if( conflicttype == SCIP_CONFTYPE_INFEASLP || conflicttype == SCIP_CONFTYPE_ALTINFPROOF )
   {
      /* add dual proof to storage */
      SCIP_CALL( SCIPconflictstoreAddDualraycons(conflictstore, cons, blkmem, set, stat, transprob, reopt, hasrelaxvar) );
   }
   else
   {
      SCIP_Real scale = 1.0;
      SCIP_Bool updateside = FALSE;

      /* In some cases the constraint could not be updated to a more special type. However, it is possible that
       * constraint got scaled. Therefore, we need to be very careful when updating the lhs/rhs after the incumbent
       * solution has improved.
       */
      if( conflicttype == SCIP_CONFTYPE_BNDEXCEEDING )
      {
         SCIP_Real side;

#ifndef NDEBUG
         SCIP_CONSHDLR* conshdlr = SCIPconsGetHdlr(cons);

         assert(conshdlr != NULL);
         assert(strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0);
#endif
         side = SCIPgetLhsLinear(set->scip, cons);

         if( !SCIPsetIsInfinity(set, -side) )
         {
            if( SCIPsetIsZero(set, side) )
            {
               scale = -1.0;
            }
            else
            {
               scale = proofsetGetRhs(proofset) / side;
               assert(SCIPsetIsNegative(set, scale));
            }
         }
         else
         {
            side = SCIPgetRhsLinear(set->scip, cons);
            assert(!SCIPsetIsInfinity(set, side));

            if( SCIPsetIsZero(set, side) )
            {
               scale = 1.0;
            }
            else
            {
               scale = proofsetGetRhs(proofset) / side;
               assert(SCIPsetIsPositive(set, scale));
            }
         }
         updateside = TRUE;
      }

      /* add dual proof to storage */
      SCIP_CALL( SCIPconflictstoreAddDualsolcons(conflictstore, cons, blkmem, set, stat, transprob, reopt, scale, updateside, hasrelaxvar) );
   }

   if( applyglobal ) /*lint !e774*/
   {
      /* add the constraint to the global problem */
      SCIP_CALL( SCIPprobAddCons(transprob, set, stat, cons) );
   }
   else
   {
      SCIP_CALL( SCIPnodeAddCons(tree->path[proofset->validdepth], blkmem, set, stat, tree, cons) );
   }

   SCIPsetDebugMsg(set, "added proof-constraint to node %p (#%lld) in depth %d (nproofconss %d)\n",
         (void*)tree->path[proofset->validdepth], SCIPnodeGetNumber(tree->path[proofset->validdepth]),
         proofset->validdepth,
         (conflicttype == SCIP_CONFTYPE_INFEASLP || conflicttype == SCIP_CONFTYPE_ALTINFPROOF)
            ? SCIPconflictstoreGetNDualInfProofs(conflictstore) : SCIPconflictstoreGetNDualBndProofs(conflictstore));

   /* release the constraint */
   SCIP_CALL( SCIPreleaseCons(set->scip, &cons) );

  UPDATESTATISTICS:
   /* update statistics */
   if( conflicttype == SCIP_CONFTYPE_INFEASLP || conflicttype == SCIP_CONFTYPE_ALTINFPROOF )
   {
      conflict->dualproofsinfnnonzeros += nnz;
      if( applyglobal ) /*lint !e774*/
         ++conflict->ndualproofsinfglobal;
      else
         ++conflict->ndualproofsinflocal;
      ++conflict->ndualproofsinfsuccess;
   }
   else
   {
      assert(conflicttype == SCIP_CONFTYPE_BNDEXCEEDING || conflicttype == SCIP_CONFTYPE_ALTBNDPROOF);
      conflict->dualproofsbndnnonzeros += nnz;
      if( applyglobal ) /*lint !e774*/
         ++conflict->ndualproofsbndglobal;
      else
         ++conflict->ndualproofsbndlocal;

      ++conflict->ndualproofsbndsuccess;
   }
   return SCIP_OKAY;
}

/** create proof constraints out of proof sets */
SCIP_RETCODE SCIPconflictFlushProofset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(conflict != NULL);

   if( proofsetGetConftype(conflict->proofset) != SCIP_CONFTYPE_UNKNOWN )
   {
      /* only one variable has a coefficient different to zero, we add this bound change instead of a constraint */
      if( SCIPproofsetGetNVars(conflict->proofset) == 1 )
      {
         SCIP_VAR** vars;
         SCIP_Real* coefs;
         int* inds;
         SCIP_Real rhs;

         vars = SCIPprobGetVars(transprob);

         coefs = proofsetGetVals(conflict->proofset);
         inds = proofsetGetInds(conflict->proofset);
         rhs = proofsetGetRhs(conflict->proofset);

         SCIP_CALL( tightenSingleVar(conflict, set, stat, tree, blkmem, origprob, transprob, reopt, lp, \
               branchcand, eventqueue, cliquetable, vars[inds[0]], coefs[0], rhs, conflict->proofset->conflicttype,
               conflict->proofset->validdepth) );
      }
      else
      {
         SCIP_Bool skipinitialproof = FALSE;

         /* prefer an infeasibility proof
          *
          * todo: check whether this is really what we want
          */
         if( set->conf_prefinfproof && proofsetGetConftype(conflict->proofset) == SCIP_CONFTYPE_BNDEXCEEDING )
         {
            int i;

            for( i = 0; i < conflict->nproofsets; i++ )
            {
               if( proofsetGetConftype(conflict->proofsets[i]) == SCIP_CONFTYPE_INFEASLP )
               {
                  skipinitialproof = TRUE;
                  break;
               }
            }
         }

         if( !skipinitialproof )
         {
            /* create and add the original proof */
            SCIP_CALL( createAndAddProofcons(conflict, conflictstore, conflict->proofset, set, stat, origprob, transprob, \
                  tree, reopt, lp, branchcand, eventqueue, cliquetable, blkmem) );
         }
      }

      /* clear the proof set anyway */
      proofsetClear(conflict->proofset);
   }

   if( conflict->nproofsets > 0 )
   {
      int i;

      for( i = 0; i < conflict->nproofsets; i++ )
      {
         assert(conflict->proofsets[i] != NULL);
         assert(proofsetGetConftype(conflict->proofsets[i]) != SCIP_CONFTYPE_UNKNOWN);

         /* only one variable has a coefficient different to zero, we add this bound change instead of a constraint */
         if( SCIPproofsetGetNVars(conflict->proofsets[i]) == 1 )
         {
            SCIP_VAR** vars;
            SCIP_Real* coefs;
            int* inds;
            SCIP_Real rhs;

            vars = SCIPprobGetVars(transprob);

            coefs = proofsetGetVals(conflict->proofsets[i]);
            inds = proofsetGetInds(conflict->proofsets[i]);
            rhs = proofsetGetRhs(conflict->proofsets[i]);

            SCIP_CALL( tightenSingleVar(conflict, set, stat, tree, blkmem, origprob, transprob, reopt, lp,
                  branchcand, eventqueue, cliquetable, vars[inds[0]], coefs[0], rhs,
                  conflict->proofsets[i]->conflicttype, conflict->proofsets[i]->validdepth) );
         }
         else
         {
            /* create and add proof constraint */
            SCIP_CALL( createAndAddProofcons(conflict, conflictstore, conflict->proofsets[i], set, stat, origprob, \
                  transprob, tree, reopt, lp, branchcand, eventqueue, cliquetable, blkmem) );
         }
      }

      /* free all proofsets */
      for( i = 0; i < conflict->nproofsets; i++ )
         SCIPproofsetFree(&conflict->proofsets[i], blkmem);

      conflict->nproofsets = 0;
   }

   return SCIP_OKAY;
}



#ifdef SCIP_DEBUG
/** print violation for debugging */
static
void debugPrintViolationInfo(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             minact,             /**< min activity */
   SCIP_Real             rhs,                /**< right hand side */
   const char*           infostr             /**< additional info for this debug message, or NULL */
   )
{
   SCIPsetDebugMsg(set, "-> %sminact=%.15g rhs=%.15g violation=%.15g\n",infostr != NULL ? infostr : "" , minact, rhs, minact - rhs);
}
#else
#define debugPrintViolationInfo(...) /**/
#endif

/** apply coefficient tightening */
static
void tightenCoefficients(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROOFSET*        proofset,           /**< proof set */
   int*                  nchgcoefs,          /**< pointer to store number of changed coefficients */
   SCIP_Bool*            redundant           /**< pointer to store whether the proof set is redundant */
   )
{
#ifdef SCIP_DEBUG
   SCIP_Real absmax = 0.0;
   SCIP_Real absmin = SCIPsetInfinity(set);
   int i;

   for( i = 0; i < proofset->nnz; i++ )
   {
      absmax = MAX(absmax, REALABS(proofset->vals[i]));
      absmin = MIN(absmin, REALABS(proofset->vals[i]));
   }
#endif

   (*redundant) = SCIPcutsTightenCoefficients(set->scip, FALSE, proofset->vals, &proofset->rhs, proofset->inds, &proofset->nnz, nchgcoefs);

#ifdef SCIP_DEBUG
   {
      SCIP_Real newabsmax = 0.0;
      SCIP_Real newabsmin = SCIPsetInfinity(set);

      for( i = 0; i < proofset->nnz; i++ )
      {
         newabsmax = MAX(newabsmax, REALABS(proofset->vals[i]));
         newabsmin = MIN(newabsmin, REALABS(proofset->vals[i]));
      }

      SCIPsetDebugMsg(set, "coefficient tightening: [%.15g,%.15g] -> [%.15g,%.15g] (nnz: %d, nchg: %d rhs: %.15g)\n",
            absmin, absmax, newabsmin, newabsmax, proofsetGetNVars(proofset), *nchgcoefs, proofsetGetRhs(proofset));
      printf("coefficient tightening: [%.15g,%.15g] -> [%.15g,%.15g] (nnz: %d, nchg: %d rhs: %.15g)\n",
            absmin, absmax, newabsmin, newabsmax, proofsetGetNVars(proofset), *nchgcoefs, proofsetGetRhs(proofset));
   }
#endif
}

/** try to generate alternative proofs by applying subadditive functions */
static
SCIP_RETCODE separateAlternativeProofs(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< tree data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_AGGRROW*         proofrow,           /**< proof rows data */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   SCIP_CONFTYPE         conflicttype        /**< type of the conflict */
   )
{
   SCIP_VAR** vars;
   SCIP_SOL* refsol;
   SCIP_Real* cutcoefs;
   SCIP_Real cutefficacy;
   SCIP_Real cutrhs;
   SCIP_Real proofefficacy;
   SCIP_Real efficacynorm;
   SCIP_Bool islocal;
   SCIP_Bool cutsuccess;
   SCIP_Bool success;
   SCIP_Bool infdelta;
   int* cutinds;
   int* inds;
   int cutnnz;
   int nnz;
   int nvars;
   int i;

   vars = SCIPprobGetVars(transprob);
   nvars = SCIPprobGetNVars(transprob);

   inds = SCIPaggrRowGetInds(proofrow);
   nnz = SCIPaggrRowGetNNz(proofrow);

   proofefficacy = SCIPaggrRowGetMinActivity(set, transprob, proofrow, curvarlbs, curvarubs, &infdelta);

   if( infdelta )
      return SCIP_OKAY;

   proofefficacy -= SCIPaggrRowGetRhs(proofrow);

   efficacynorm = SCIPaggrRowCalcEfficacyNorm(set->scip, proofrow);
   proofefficacy /= MAX(1e-6, efficacynorm);

   /* create reference solution */
   SCIP_CALL( SCIPcreateSol(set->scip, &refsol, NULL) );

   /* initialize with average solution */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPsolSetVal(refsol, set, stat, tree, vars[i], SCIPvarGetAvgSol(vars[i])) );
   }

   /* set all variables that are part of the proof to its active local bound */
   for( i = 0; i < nnz; i++ )
   {
      SCIP_Real val = SCIPaggrRowGetProbvarValue(proofrow, inds[i]);

      if( val > 0.0 )
      {
         SCIP_CALL( SCIPsolSetVal(refsol, set, stat, tree, vars[inds[i]], curvarubs[inds[i]]) );
      }
      else
      {
         SCIP_CALL( SCIPsolSetVal(refsol, set, stat, tree, vars[inds[i]], curvarlbs[inds[i]]) );
      }
   }

   SCIP_CALL( SCIPsetAllocBufferArray(set, &cutcoefs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cutinds, nvars) );

   cutnnz = 0;
   cutefficacy = -SCIPsetInfinity(set);

   /* apply flow cover */
   SCIP_CALL( SCIPcalcFlowCover(set->scip, refsol, POSTPROCESS, BOUNDSWITCH, ALLOWLOCAL, proofrow, \
         cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy, NULL, &islocal, &cutsuccess) );
   success = cutsuccess;

   /* apply MIR */
   SCIP_CALL( SCIPcutGenerationHeuristicCMIR(set->scip, refsol, POSTPROCESS, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, INT_MAX, \
         NULL, NULL, MINFRAC, MAXFRAC, proofrow, cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy, NULL, \
         &islocal, &cutsuccess) );
   success = (success || cutsuccess);

   /* replace the current proof */
   if( success && !islocal && SCIPsetIsPositive(set, cutefficacy) && cutefficacy * nnz > proofefficacy * cutnnz )
   {
      SCIP_PROOFSET* alternativeproofset;
      SCIP_Bool redundant;
      int nchgcoefs;

      SCIP_CALL( proofsetCreate(&alternativeproofset, blkmem) );
      alternativeproofset->conflicttype = (conflicttype == SCIP_CONFTYPE_INFEASLP ? SCIP_CONFTYPE_ALTINFPROOF : SCIP_CONFTYPE_ALTBNDPROOF);

      SCIP_CALL( proofsetAddSparseData(alternativeproofset, blkmem, cutcoefs, cutinds, cutnnz, cutrhs) );

      /* apply coefficient tightening */
      tightenCoefficients(set, alternativeproofset, &nchgcoefs, &redundant);

      if( !redundant )
      {
         SCIP_CALL( conflictInsertProofset(conflict, set, alternativeproofset) );
      }
      else
      {
         SCIPproofsetFree(&alternativeproofset, blkmem);
      }
   }  /*lint !e438*/

   SCIPsetFreeBufferArray(set, &cutinds);
   SCIPsetFreeBufferArray(set, &cutcoefs);

   SCIP_CALL( SCIPfreeSol(set->scip, &refsol) );

   return SCIP_OKAY;
}

/** tighten a given infeasibility proof a^Tx <= b with minact > b w.r.t. local bounds
 *
 *  1) Apply cut generating functions
 *    - c-MIR
 *    - Flow-cover
 *    - TODO: implement other subadditive functions
 *  2) Remove continuous variables contributing with its global bound
 *    - TODO: implement a variant of non-zero-cancellation
 */
static
SCIP_RETCODE tightenDualproof(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< tree data */
   SCIP_AGGRROW*         proofrow,           /**< aggregated row representing the proof */
   int                   validdepth,         /**< depth where the proof is valid */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   SCIP_Bool             initialproof        /**< do we analyze the initial reason of infeasibility? */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int* inds;
   SCIP_PROOFSET* proofset;
   SCIP_Bool valid;
   SCIP_Bool redundant;
   int nnz;
   int nchgcoefs;
   int nbinvars;
   int ncontvars;
   int nintvars;
   int i;

   assert(conflict->proofset != NULL);
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);

   vars = SCIPprobGetVars(transprob);
   nbinvars = 0;
   nintvars = 0;
   ncontvars = 0;

   inds = SCIPaggrRowGetInds(proofrow);
   nnz = SCIPaggrRowGetNNz(proofrow);

   /* count number of binary, integer, and continuous variables */
   for( i = 0; i < nnz; i++ )
   {
      assert(SCIPvarGetProbindex(vars[inds[i]]) == inds[i]);

      if( SCIPvarIsBinary(vars[inds[i]]) )
         ++nbinvars;
      else if( SCIPvarIsIntegral(vars[inds[i]]) )
         ++nintvars;
      else
         ++ncontvars;
   }

   SCIPsetDebugMsg(set, "start dual proof tightening:\n");
   SCIPsetDebugMsg(set, "-> tighten dual proof: nvars=%d (bin=%d, int=%d, cont=%d)\n",
         nnz, nbinvars, nintvars, ncontvars);
   debugPrintViolationInfo(set, aggrRowGetMinActivity(set, transprob, proofrow, curvarlbs, curvarubs, NULL), SCIPaggrRowGetRhs(proofrow), NULL);

   /* try to find an alternative proof of local infeasibility that is stronger */
   if( set->conf_sepaaltproofs )
   {
      SCIP_CALL( separateAlternativeProofs(conflict, set, stat, transprob, tree, blkmem, proofrow, curvarlbs, curvarubs,
            conflict->conflictset->conflicttype) );
   }

   if( initialproof )
      proofset = conflict->proofset;
   else
   {
      SCIP_CALL( proofsetCreate(&proofset, blkmem) );
   }

   /* start with a proofset containing all variables with a non-zero coefficient in the dual proof */
   SCIP_CALL( proofsetAddAggrrow(proofset, set, blkmem, proofrow) );
   proofset->conflicttype = conflict->conflictset->conflicttype;
   proofset->validdepth = validdepth;

   /* get proof data */
   vals = proofsetGetVals(proofset);
   inds = proofsetGetInds(proofset);
   nnz = SCIPproofsetGetNVars(proofset);

#ifndef NDEBUG
   for( i = 0; i < nnz; i++ )
   {
      int idx = inds[i];
      if( vals[i] > 0.0 )
         assert(!SCIPsetIsInfinity(set, -curvarlbs[idx]));
      if( vals[i] < 0.0 )
         assert(!SCIPsetIsInfinity(set, curvarubs[idx]));
   }
#endif

   /* remove continuous variable contributing with their global bound
    *
    * todo: check whether we also want to do that for bound exceeding proofs, but then we cannot update the
    *       conflict anymore
    */
   if( proofset->conflicttype == SCIP_CONFTYPE_INFEASLP )
   {
      /* remove all continuous variables that have equal global and local bounds (ub or lb depend on the sign)
       * from the proof
       */

      for( i = 0; i < nnz && nnz > 1; )
      {
         SCIP_Real val;
         int idx = inds[i];

         assert(vars[idx] != NULL);

         val = vals[i];
         assert(!SCIPsetIsZero(set, val));

         /* skip integral variables */
         if( SCIPvarGetType(vars[idx]) != SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(vars[idx]) != SCIP_VARTYPE_IMPLINT )
         {
            i++;
            continue;
         }
         else
         {
            SCIP_Real glbbd;
            SCIP_Real locbd;

            /* get appropriate global and local bounds */
            glbbd = (val < 0.0 ? SCIPvarGetUbGlobal(vars[idx]) : SCIPvarGetLbGlobal(vars[idx]));
            locbd = (val < 0.0 ? curvarubs[idx] : curvarlbs[idx]);

            if( !SCIPsetIsEQ(set, glbbd, locbd) )
            {
               i++;
               continue;
            }

            SCIPsetDebugMsg(set, "-> remove continuous variable <%s>: glb=[%g,%g], loc=[%g,%g], val=%g\n",
                  SCIPvarGetName(vars[idx]), SCIPvarGetLbGlobal(vars[idx]), SCIPvarGetUbGlobal(vars[idx]),
                  curvarlbs[idx], curvarubs[idx], val);

            proofsetCancelVarWithBound(proofset, set, vars[idx], i, &valid);
            assert(valid); /* this should be always fulfilled at this place */

            --nnz;
         }
      }
   }

   /* apply coefficient tightening to initial proof */
   tightenCoefficients(set, proofset, &nchgcoefs, &redundant);

   /* it can happen that the constraints is almost globally redundant w.r.t to the maximal activity,
    * e.g., due to numerics. in this case, we want to discard the proof
    */
   if( redundant )
   {
#ifndef NDEBUG
      SCIP_Real eps = MIN(0.01, 10.0*set->num_feastol);
      assert(proofset->rhs - getMaxActivity(set, transprob, proofset->vals, proofset->inds, proofset->nnz, NULL, NULL) < eps);
#endif
      if( initialproof )
      {
         proofsetClear(proofset);
      }
      else
      {
         SCIPproofsetFree(&proofset, blkmem);
      }
   }
   else
   {
      if( !initialproof )
      {
         SCIP_CALL( conflictInsertProofset(conflict, set, proofset) );
      }

      if( nchgcoefs > 0 )
      {
         if( proofset->conflicttype == SCIP_CONFTYPE_INFEASLP )
            proofset->conflicttype = SCIP_CONFTYPE_ALTINFPROOF;
         else if( proofset->conflicttype == SCIP_CONFTYPE_BNDEXCEEDING )
            proofset->conflicttype = SCIP_CONFTYPE_ALTBNDPROOF;
      }
   }

   return SCIP_OKAY;
}

/** perform conflict analysis based on a dual unbounded ray
 *
 *  given an aggregation of rows lhs <= a^Tx such that lhs > maxactivity. if the constraint has size one we add a
 *  bound change instead of the constraint.
 */
SCIP_RETCODE SCIPconflictAnalyzeDualProof(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< tree data */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_AGGRROW*         proofrow,           /**< aggregated row representing the proof */
   int                   validdepth,         /**< valid depth of the dual proof */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   SCIP_Bool             initialproof,       /**< do we analyze the initial reason of infeasibility? */
   SCIP_Bool*            globalinfeasible,   /**< pointer to store whether global infeasibility could be proven */
   SCIP_Bool*            success             /**< pointer to store success result */
   )
{
   SCIP_Real rhs;
   SCIP_Real minact;
   SCIP_Bool infdelta;
   int nnz;

   assert(set != NULL);
   assert(transprob != NULL);
   assert(validdepth >= 0);
   assert(validdepth == 0 || validdepth < SCIPtreeGetFocusDepth(tree));

   /* get sparse data */
   nnz = SCIPaggrRowGetNNz(proofrow);
   rhs = SCIPaggrRowGetRhs(proofrow);

   *globalinfeasible = FALSE;
   *success = FALSE;

   /* get minimal activity w.r.t. local bounds */
   minact = SCIPaggrRowGetMinActivity(set, transprob, proofrow, curvarlbs, curvarubs, &infdelta);

   if( infdelta )
      return SCIP_OKAY;

   /* only run is the proof proves local infeasibility */
   if( SCIPsetIsFeasLE(set, minact, rhs) )
      return SCIP_OKAY;

   /* if the farkas-proof is empty, the node and its sub tree can be cut off completely */
   if( nnz == 0 )
   {
      SCIPsetDebugMsg(set, " -> empty farkas-proof in depth %d cuts off sub tree at depth %d\n", SCIPtreeGetFocusDepth(tree), validdepth);

      SCIP_CALL( SCIPnodeCutoff(tree->path[validdepth], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );

      *globalinfeasible = TRUE;
      *success = TRUE;

      ++conflict->ndualproofsinfsuccess;

      return SCIP_OKAY;
   }

   /* try to enforce the constraint based on a dual ray */
   SCIP_CALL( tightenDualproof(conflict, set, stat, blkmem, transprob, tree, proofrow, validdepth,
      curvarlbs, curvarubs, initialproof) );

   if( *globalinfeasible )
   {
      SCIPsetDebugMsg(set, "detect global: cutoff root node\n");
      SCIP_CALL( SCIPnodeCutoff(tree->path[0], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
      *success = TRUE;

      ++conflict->ndualproofsinfsuccess;
   }

   return SCIP_OKAY;
}
