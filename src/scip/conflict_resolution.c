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

   /* @todo check this */
   --resolutionset->nnz;
   resolutionset->vals[pos] = resolutionset->vals[resolutionset->nnz];
   resolutionset->inds[pos] = resolutionset->inds[resolutionset->nnz];
}

/** return the indices of variables in the resolutionset */
static
int* resolutionsetGetInds(
   SCIP_RESOLUTIONSET*   resolutionset       /**< resolution set */
   )
{
   assert(resolutionset != NULL);

   return resolutionset->inds;
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

/** return the left-hand side of row from which we derive the resolutionset */
static
SCIP_Real resolutionsetGetOrigLhs(
   SCIP_RESOLUTIONSET*   resolutionset  /**< resolution set */
   )
{
   assert(resolutionset != NULL);

   return resolutionset->origlhs;
}

/** return the right-hand side of row from which we derive the resolutionset */
static
SCIP_Real resolutionsetGetOrigRhs(
   SCIP_RESOLUTIONSET*   resolutionset       /**< resolution set */
   )
{
   assert(resolutionset != NULL);

   return resolutionset->origrhs;
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

   SCIPsetDebugMsg(set, "resolution set with LHS: %f \n", resolutionset->lhs);
   SCIPsetDebugMsg(set, "resolution set with orig LHS: %f, and orig RHS: %f \n", resolutionset->origlhs,resolutionset->origrhs);

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
         if ( nubchgs > 0 )
         {
            bnd = SCIPgetVarUbAtIndex(scip, vars[v], bdchgidx, 1);
            SCIPsetDebugMsg(set, "Variable <%s>, coef: %f, UB changes: %d, UB used for slack is %f \n", SCIPvarGetName(vars[v]), coef, nubchgs, bnd);
         }
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
         if ( nlbchgs > 0 )
         {
            bnd = SCIPgetVarLbAtIndex(scip, vars[v], bdchgidx, 1);
            SCIPsetDebugMsg(set, "Variable <%s>, coef: %f, LB changes: %d, LB used for slack is %f \n", SCIPvarGetName(vars[v]), coef, nlbchgs, bnd);
         }
         else bnd = SCIPvarGetLbLocal(vars[v]);

         SCIPquadprecProdDD(delta, coef, bnd);
      }
      SCIPquadprecSumQQ(slack, slack, delta);
   }
   SCIPquadprecSumQD(slack, slack, -resolutionset->lhs);
   return QUAD_TO_DBL(slack);
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
   vars = SCIPprobGetVars(transprob);
   assert(vars != NULL);

   /* @todo */
   SCIP_VAR** consvars;
   SCIP_Real* vals;
   SCIP_Real lhs;
   int i;

   vals = resolutionset->vals;

   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, resolutionsetGetNNzs(resolutionset)) );

   lhs = resolutionsetGetLhs(resolutionset);

   for( i = 0; i < resolutionsetGetNNzs(resolutionset); ++i )
   {
      consvars[i] = vars[resolutionset->inds[i]];
   }

   SCIP_CONS* cons;
   SCIP_CONS* upgdcons;

   char consname[SCIP_MAXSTRLEN];


   /* create a constraint out of the conflict set */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "confres%" SCIP_LONGINT_FORMAT, conflict->nresconfconss);
   SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, resolutionsetGetNNzs(resolutionset), consvars, vals, lhs, SCIPinfinity(scip),
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
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(tree != NULL);

      int focusdepth;
#ifndef NDEBUG
      int currentdepth;
#endif
      int cutoffdepth;
      int repropdepth;
      int maxconflictsets;
      int maxsize;
      int i;

      /* calculate the maximal number of conflict sets to accept, and the maximal size of each accepted conflict set */
      maxconflictsets = (set->conf_maxconss == -1 ? INT_MAX : set->conf_maxconss);
      maxsize = conflictCalcMaxsize(set, transprob);

      focusdepth = SCIPtreeGetFocusDepth(tree);
#ifndef NDEBUG
      currentdepth = SCIPtreeGetCurrentDepth(tree);
      assert(focusdepth <= currentdepth);
      assert(currentdepth == tree->pathlen-1);
#endif

      SCIPsetDebugMsg(set, "flushing %d resolution sets at focus depth %d (maxsize: %d)\n",
         1, focusdepth, maxsize);

   /* insert the conflict sets at the corresponding nodes */
   cutoffdepth = INT_MAX;
   repropdepth = INT_MAX;

   /* @todo add a check for the number of resolution sets */
   /* @todo loop over all resolution sets */
   SCIP_RESOLUTIONSET* resolutionset;
   resolutionset = conflict->resolutionset;

   assert(resolutionset != NULL);
   assert(0 <= resolutionset->validdepth);
   assert( getSlack(scip, set, transprob, resolutionset) < 0 );


   /* if the resolution set is empty, the node and its sub tree in the conflict set's valid depth can be
      * cut off completely
      */
   if( resolutionsetGetNNzs(resolutionset) == 0 )
   {
      SCIPsetDebugMsg(set, " -> empty resolution set in depth %d cuts off sub tree at depth %d\n",
         focusdepth, resolutionset->validdepth);

      SCIP_CALL( SCIPnodeCutoff(tree->path[resolutionset->validdepth], set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
      cutoffdepth = resolutionset->validdepth;
      return SCIP_OKAY;
   }
   /* if the conflict set is too long, use the conflict set only if it decreases the repropagation depth */
   if( resolutionsetGetNNzs(resolutionset) > maxsize )
   {
      SCIPsetDebugMsg(set, " -> conflict set is too long: %d > %d nonzeros\n", resolutionsetGetNNzs(resolutionset), maxsize);
      /* @todo keep constraint for repropagation */
   }
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

   /* @todo free the conflict */

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
   int                   residx                    /**< index of variable to resolve */
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
   /* @todo add asserts to make sure that residx exist in bot sets! */
   for( i = 0; i < conflictresolutionset->nnz; i++ )
   {
      if (conflictresolutionset->inds[i] == residx)
      {
         idxinconflict = TRUE;
         coefconf = conflictresolutionset->vals[i];
         break;
      }
   }

   for( i = 0; i < reasonresolutionset->nnz; i++ )
   {
      if (reasonresolutionset->inds[i] == residx)
      {
         idxinreason = TRUE;
         coefreas = reasonresolutionset->vals[i];
         break;
      }
   }
   SCIPsetDebugMsg(set, "Nonzeros in conflict constraint: %d \n", conflictresolutionset->nnz);
   SCIPsetDebugMsg(set, "Nonzeros in reason constraint: %d \n", reasonresolutionset->nnz);

   assert(coefconf * coefreas < 0);
   assert(idxinconflict && idxinreason);

   newsize = conflictresolutionset->nnz + reasonresolutionset->nnz;
   SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictresolutionset->inds, conflictresolutionset->size, newsize ) );
   SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictresolutionset->vals, conflictresolutionset->size, newsize ) );
   conflictresolutionset->size = newsize;

   cidx = 0;
   previousnnz = conflictresolutionset->nnz;

   /* multiply conflict by coefreas */
   for( i = 0; i < conflictresolutionset->nnz; i++ )
   {
      conflictresolutionset->vals[i] = fabs(coefreas) * conflictresolutionset->vals[i];
   }
   /* multiply reason by coefconf */
   for( i = 0; i < reasonresolutionset->nnz; i++ )
   {
      reasonresolutionset->vals[i] = fabs(coefconf) * reasonresolutionset->vals[i];
   }


   i = 0;
   /* add conflict and reason resolution sets */
   while ( i < reasonresolutionset->nnz )
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
   conflictresolutionset->lhs = fabs(coefreas) * conflictresolutionset->lhs + fabs(coefconf) * reasonresolutionset->lhs;

   // for( i = 0; i < conflictresolutionset->nnz; i++ )
   // {
   //    if (SCIPsetIsZero(set, conflictresolutionset->vals[i] ))
   //       resolutionsetRemoveZeroVar(conflictresolutionset, set, i);
   // }

   SCIPsetDebugMsg(set, "Nonzeros in resolved constraint: %d \n", conflictresolutionset->nnz);

   /* sort for linear time resolution */
   SCIPsortIntReal(conflictresolutionset->inds, conflictresolutionset->vals, conflictresolutionset->nnz);

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
   SCIP_Real* changedvals;
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

   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &changedvals, vals, nnz) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &inds, nnz) );

   isincon = FALSE;
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
         isincon = TRUE;
      }
   }
   assert(isincon);

   if ( changesign )
   {
      lhs = -origrhs;
      for( i = 0; i < nnz; i++ )
      {
         changedvals[i] = -vals[i];
      }
   }
   else lhs = origlhs;

   assert(inds != NULL);
   assert(vals != NULL);

   SCIP_CALL( resolutionsetAddSparseData(scip, resolutionset, blkmem, changedvals, inds, nnz, lhs, origrhs, origlhs, row) );

   BMSfreeBlockMemoryArray(blkmem, &changedvals, nnz);
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
   SCIP_Real* changedvals;
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

   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &changedvals, vals, nnz) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &inds, nnz) );

   isincon = FALSE;
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
         isincon = TRUE;
      }
   }
   assert(isincon);

   if ( changesign )
   {
      lhs = -origrhs;
      for( i = 0; i < nnz; i++ )
      {
         changedvals[i] = -vals[i];
      }
   }
   else lhs = origlhs;

   assert(inds != NULL);
   assert(vals != NULL);

   SCIP_CALL( resolutionsetAddSparseData(scip, resolutionset, blkmem, changedvals, inds, nnz, lhs, origrhs, origlhs, row) );

   BMSfreeBlockMemoryArray(blkmem, &changedvals, nnz);
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
   SCIP_ROW* conflictrow;
   SCIP_ROW* reasonrow;
   SCIP_CONS* reasoncon;

   int focusdepth;
   int currentdepth;
   int maxvaliddepth;

   SCIP_VAR* vartoresolve;
   int residx;

   SCIP_Real conflictslack;
   SCIP_Real reasonslack;

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

   /* if the bound change was not infered by a constraint (reason) then conflict analysis is not applicable */
   if ( SCIPbdchginfoGetChgtype(bdchginfo) != SCIP_BOUNDCHGTYPE_CONSINFER )
   {
      SCIPsetDebugMsg(set, "Conflict analysis not applicable since no reason row exists \n");
      return SCIP_OKAY;
   }

   SCIP_CALL( resolutionsetCreate(&conflictresolutionset, blkmem) );
   conflictresolutionset->validdepth = validdepth;

   /* get the resolution set of the conflict row */
   SCIP_CALL( conflictResolutionsetFromRow(scip, conflictresolutionset, bdchginfo, set, blkmem, conflictrow) );

   /* sort for linear time resolution*/
   SCIPsortIntReal(conflictresolutionset->inds, conflictresolutionset->vals, conflictresolutionset->nnz);

   conflictslack = getSlack(scip, set, transprob, conflictresolutionset);

   /* slack should be negative */
   /* @todo assert(conflictslack < 0); */
   if ( conflictslack >= 0 )
   {
      SCIPresolutionsetFree(&conflictresolutionset, blkmem);
      return SCIP_OKAY;
   }

   /* @todo apply coefficient tightening for conflict set */

   if ( SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER )
   {
      reasoncon = SCIPbdchginfoGetInferCons(bdchginfo);

      /* get the corresponding reason row */
      reasonrow = SCIPconsGetRow(scip, reasoncon);

      /* if no row exists conflict analysis is not applicable */
      if (reasonrow == NULL)
      {
         SCIPresolutionsetFree(&conflictresolutionset, blkmem);
         return SCIP_OKAY;
      }

      SCIP_CALL( resolutionsetCreate(&reasonresolutionset, blkmem) );

      /* get the resolution set of the conflict row */
      SCIP_CALL( reasonResolutionsetFromRow(scip, reasonresolutionset, bdchginfo, set, blkmem, reasonrow) );
      /* sort for linear time resolution*/
      SCIPsortIntReal(reasonresolutionset->inds, reasonresolutionset->vals, reasonresolutionset->nnz);
      reasonslack = getSlack(scip, set, transprob, reasonresolutionset);

      /* @todo assert(reasonslack >= 0); */
      /* @todo apply coefficient tightening for reason set */

   }
   else
   {
      SCIPresolutionsetFree(&conflictresolutionset, blkmem);
      return SCIP_OKAY;
   }
   /* @todo apply weakining and coef. tightening is  conflictslack + reasonslack >= 0 */
   /* for now apply resolution only if conflictslack + reasonslack < 0 */
   if ( conflictslack + reasonslack >= 0 )
   {
      SCIPresolutionsetFree(&conflictresolutionset, blkmem);
      SCIPresolutionsetFree(&reasonresolutionset, blkmem);
      return SCIP_OKAY;
   }

   SCIPsetDebugMsg(set, "Slack of conflict: %f, Slack of reason %f \n", conflictslack, reasonslack);
   vartoresolve = bdchginfo->var;
   residx = SCIPvarGetProbindex(vartoresolve);
   resolveWithReason(scip, set, conflictresolutionset, reasonresolutionset, blkmem, residx);

   conflictslack = getSlack(scip, set, transprob, conflictresolutionset);

   if ( SCIPsetIsLT(set, conflictslack, 0.0) )
   {
      /* @todo add constraint */
      conflict->resolutionset = conflictresolutionset;

      /* @todo call flush from the main solving loop! */
      SCIPconflictFlushResolutionSets(conflict, blkmem, scip, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, cliquetable);
      (*nconss)++;
      (*nconfvars) = resolutionsetGetNNzs(conflictresolutionset);
      SCIPresolutionsetFree(&conflictresolutionset, blkmem);
      SCIPresolutionsetFree(&reasonresolutionset, blkmem);
      SCIPresolutionsetFree(&conflict->resolutionset, blkmem);
      return SCIP_OKAY;
   }
   SCIPsetDebugMsg(set, "Slack of resolved row: %f \n", conflictslack);

   /* main loop */
   while( bdchginfo != NULL && validdepth <= maxvaliddepth )
   {
      break;
      /* @todo add bound changes in the conflict queue ( will not work well with graph conflict analysis */
      /* @todo go through variables and weaken + coef. tightening till the resolved constraint has negative slack */
   }

   SCIPresolutionsetFree(&conflictresolutionset, blkmem);
   SCIPresolutionsetFree(&reasonresolutionset, blkmem);
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

   /* @todo check, if the conflict constraint will get too large with high probability */

   SCIPsetDebugMsg(set, "resolution based conflict analysis after infeasible propagation in depth %d\n", SCIPtreeGetCurrentDepth(tree));

   /* start timing */
   SCIPclockStart(conflict->resanalyzetime, set);

   conflict->nrescalls++;

   /* analyze the conflict set, and create a conflict constraint on success */
   SCIP_CALL( conflictAnalyzeResolution(scip, conflict, cons, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, cliquetable, FALSE, validdepth, TRUE, &nconss, &nconfvars) );

   conflict->nressuccess += (nconss > 0 ? 1 : 0);
   conflict->nresconfconss += nconss;
   conflict->nresconfvariables += nconfvars;
   if( success != NULL )
      *success = (nconss > 0);

   /* stop timing */
   SCIPclockStop(conflict->resanalyzetime, set);

   return SCIP_OKAY;
}
