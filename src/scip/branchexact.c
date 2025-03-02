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

/**@file   branchexact.c
 * @ingroup OTHER_CFILES
 * @brief  methods for branching rules and branching candidate storage (exact SCIP version)
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "blockmemshell/memory.h"
#include "scip/branch.h"
#include "scip/branchexact.h"
#include "scip/certificate.h"
#include "scip/def.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/prob.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/scip.h"
#include "scip/struct_branch.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/visual.h"


/** ensures, that lpcands array can store at least num entries */
static
SCIP_RETCODE ensureLpcandsSize(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(branchcand->nlpcands <= branchcand->lpcandssize);

   if( num > branchcand->lpcandssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&branchcand->lpcands, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&branchcand->lpcandssol, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&branchcand->lpcandsfrac, newsize) );
      branchcand->lpcandssize = newsize;
   }
   assert(num <= branchcand->lpcandssize);

   return SCIP_OKAY;
}

/** calculates branching candidates for LP solution branching (fractional variables) */
static
SCIP_RETCODE branchcandCalcLPCandsExact(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_COL** cols;
   SCIP_VAR* var;
   SCIP_COL* col;
   SCIP_Rational* tmp;
   SCIP_Real primsol;
   SCIP_Real frac;
   SCIP_VARTYPE vartype;
   SCIP_IMPLINTTYPE impltype;
   int branchpriority;
   int ncols;
   int c;
   int insertpos;

   assert(branchcand != NULL);
   assert(stat != NULL);
   assert(branchcand->validlpcandslp <= stat->lpcount);
   assert(lp != NULL);
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);
   /* this should only be called when fp-methods did not find any candidates */
   assert(branchcand->nlpcands == 0);
   assert(branchcand->validlpcandslp == stat->lpcount);

   SCIPsetDebugMsg(set, "calculating LP branching candidates exactly: validlp=%" SCIP_LONGINT_FORMAT ", lpcount=%" SCIP_LONGINT_FORMAT "\n",
      branchcand->validlpcandslp, stat->lpcount);

   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      branchcand->lpmaxpriority = INT_MIN / 2;
      branchcand->nlpcands = 0;
      branchcand->npriolpcands = 0;
      branchcand->npriolpbins = 0;
      branchcand->nimpllpfracs = 0;
      branchcand->validlpcandslp = stat->lpcount;

      SCIPsetDebugMsg(set, " LP is unbounded -> no branching candidates\n");
      return SCIP_OKAY;
   }

   /* always recalculate */
   SCIPsetDebugMsg(set, " -> recalculating LP branching candidates exactly\n");

   cols = SCIPlpGetCols(lp);
   ncols = SCIPlpGetNCols(lp);

   /* construct the LP branching candidate set, moving the candidates with maximal priority to the front */
   SCIP_CALL( ensureLpcandsSize(branchcand, set, ncols) );

   branchcand->lpmaxpriority = INT_MIN / 2;
   branchcand->nlpcands = 0;
   branchcand->nimpllpfracs = 0;
   branchcand->npriolpcands = 0;
   branchcand->npriolpbins = 0;

   SCIP_CALL( RatCreateBuffer(set->buffer, &tmp) );

   for( c = 0; c < ncols; ++c )
   {
      col = cols[c];
      assert(col != NULL);
      assert(col->lppos == c);
      assert(col->lpipos >= 0);

      primsol = SCIPcolGetPrimsol(col);
      assert(primsol < SCIP_INVALID);
      assert(SCIPsetIsInfinity(set, -col->lb) || SCIPsetIsFeasGE(set, primsol, col->lb));
      assert(SCIPsetIsInfinity(set, col->ub) || SCIPsetIsFeasLE(set, primsol, col->ub));

      var = col->var;
      assert(var != NULL);
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(var) == col);

      /* LP branching candidates are fractional binary and integer variables; implicit variables are kept at the end
      * of the candidates array for some rounding heuristics
      */
      vartype = SCIPvarGetType(var);
      impltype = SCIPvarGetImplType(var);
      if( vartype == SCIP_VARTYPE_CONTINUOUS && impltype == SCIP_IMPLINTTYPE_NONE )
         continue;

      /* ignore fixed variables (due to numerics, it is possible, that the LP solution of a fixed integer variable
      * (with large fixed value) is fractional in terms of absolute feasibility measure)
      */
      /** @todo exip: this should be fine, right? */
      if( SCIPvarGetLbLocal(var) >= SCIPvarGetUbLocal(var) - 0.5 )
         continue;

      RatSetReal(tmp, primsol);

      /* check, if the LP solution value is fractional */
      if( RatIsIntegral(tmp) )
         continue;

      RatGetFrac(tmp, tmp);
      frac = RatApproxReal(tmp);

      /* insert candidate in candidate list */
      branchpriority = SCIPvarGetBranchPriority(var);
      insertpos = branchcand->nlpcands + branchcand->nimpllpfracs;
      assert(insertpos < branchcand->lpcandssize);

      if( impltype != SCIP_IMPLINTTYPE_NONE )
         branchpriority = INT_MIN;

      assert(impltype != SCIP_IMPLINTTYPE_NONE || branchpriority >= INT_MIN/2);
      /* ensure that implicit variables are stored at the end of the array */
      if( impltype == SCIP_IMPLINTTYPE_NONE && branchcand->nimpllpfracs > 0 )
      {
         assert(branchcand->lpcands[branchcand->nlpcands] != NULL
               && SCIPvarIsImpliedIntegral(branchcand->lpcands[branchcand->nlpcands]));

         branchcand->lpcands[insertpos] = branchcand->lpcands[branchcand->nlpcands];
         branchcand->lpcandssol[insertpos] = branchcand->lpcandssol[branchcand->nlpcands];
         branchcand->lpcandsfrac[insertpos] = branchcand->lpcandsfrac[branchcand->nlpcands];

         insertpos = branchcand->nlpcands;
      }

      if( branchpriority > branchcand->lpmaxpriority )
      {
         /* candidate has higher priority than the current maximum:
         * move it to the front and declare it to be the single best candidate
         */
         if( insertpos != 0 )
         {
            branchcand->lpcands[insertpos] = branchcand->lpcands[0];
            branchcand->lpcandssol[insertpos] = branchcand->lpcandssol[0];
            branchcand->lpcandsfrac[insertpos] = branchcand->lpcandsfrac[0];
            insertpos = 0;
         }
         branchcand->npriolpcands = 1;
         assert(!SCIPvarIsImpliedIntegral(var));
         branchcand->npriolpbins = (vartype == SCIP_VARTYPE_BINARY ? 1 : 0);
         branchcand->lpmaxpriority = branchpriority;
      }
      else if( branchpriority == branchcand->lpmaxpriority )
      {
         /* candidate has equal priority as the current maximum:
         * move away the first non-maximal priority candidate, move the current candidate to the correct
         * slot (binaries first) and increase the number of maximal priority candidates
         */
         if( insertpos != branchcand->npriolpcands )
         {
            branchcand->lpcands[insertpos] = branchcand->lpcands[branchcand->npriolpcands];
            branchcand->lpcandssol[insertpos] = branchcand->lpcandssol[branchcand->npriolpcands];
            branchcand->lpcandsfrac[insertpos] = branchcand->lpcandsfrac[branchcand->npriolpcands];
            insertpos = branchcand->npriolpcands;
         }
         branchcand->npriolpcands++;

         assert(!SCIPvarIsImpliedIntegral(var));
         if( vartype == SCIP_VARTYPE_BINARY )
         {
            if( insertpos != branchcand->npriolpbins )
            {
               branchcand->lpcands[insertpos] = branchcand->lpcands[branchcand->npriolpbins];
               branchcand->lpcandssol[insertpos] = branchcand->lpcandssol[branchcand->npriolpbins];
               branchcand->lpcandsfrac[insertpos] = branchcand->lpcandsfrac[branchcand->npriolpbins];
               insertpos = branchcand->npriolpbins;
            }
            branchcand->npriolpbins++;
         }
      }
      /* insert variable at the correct position of the candidates storage */
      branchcand->lpcands[insertpos] = var;
      branchcand->lpcandssol[insertpos] = primsol;
      branchcand->lpcandsfrac[insertpos] = frac;

      /* increase the counter depending on the variable type */
      if( impltype == SCIP_IMPLINTTYPE_NONE )
         branchcand->nlpcands++;
      else
         branchcand->nimpllpfracs++;

      SCIPsetDebugMsg(set, " -> candidate %d: var=<%s>, sol=%g, frac=%g, prio=%d (max: %d) -> pos %d\n",
         branchcand->nlpcands, SCIPvarGetName(var), primsol, frac, branchpriority, branchcand->lpmaxpriority,
         insertpos);
   }

#ifndef NDEBUG
   /* in debug mode we assert that the variables are positioned correctly (binaries and integers first,
   * implicit integers last)
   */
   for( c = 0; c < branchcand->nlpcands + branchcand->nimpllpfracs; ++c )
   {
      assert(SCIPvarIsImpliedIntegral(branchcand->lpcands[c]) == (c >= branchcand->nlpcands));
   }
#endif

   branchcand->validlpcandslp = stat->lpcount;

   assert(0 <= branchcand->npriolpcands && branchcand->npriolpcands <= branchcand->nlpcands);

   SCIPsetDebugMsg(set, " -> %d fractional variables (%d of maximal priority)\n", branchcand->nlpcands, branchcand->npriolpcands);

   RatFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** branches on a variable x; unlike the fp-version this will also branch x <= floor(x'), x >= ceil(x')
 * if x' is very close to being integral at one of its bounds;
 * in the fp version this case would be branched in the middle of the domain;
 * if x' is almost integral but not at a bound, this will branch (x <= x'-1, x == x', x >= x'+1);
 * not meant for branching on a continuous variables
 */
SCIP_RETCODE SCIPtreeBranchVarExact(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_NODE**           downchild,          /**< pointer to return the left child with variable rounded down, or NULL */
   SCIP_NODE**           eqchild,            /**< pointer to return the middle child with variable fixed, or NULL */
   SCIP_NODE**           upchild             /**< pointer to return the right child with variable rounded up, or NULL */
   )
{
   SCIP_NODE* node;
   SCIP_Real priority;
   SCIP_Real estimate;

   SCIP_Real downub;
   SCIP_Real uplb;
   SCIP_Real val;

   assert(tree != NULL);
   assert(set != NULL);
   assert(var != NULL);

   /* initialize children pointer */
   if( downchild != NULL )
      *downchild = NULL;
   if( eqchild != NULL )
      *eqchild = NULL;
   if( upchild != NULL )
      *upchild = NULL;

   var = SCIPvarGetProbvar(var);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIPerrorMessage("cannot branch on fixed or multi-aggregated variable <%s>\n", SCIPvarGetName(var));
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   /* ensure, that branching on continuous variables will only be performed when a branching point is given. */
   if( !SCIPvarIsIntegral(var) )
   {
      SCIPerrorMessage("Cannot branch exactly on continuous variable %s.\n", SCIPvarGetName(var));
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   assert(SCIPvarIsActive(var));
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPsetIsFeasIntegral(set, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsFeasIntegral(set, SCIPvarGetUbLocal(var)));
   assert(SCIPsetIsLT(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

   /* update the information for the focus node before creating children */
   SCIP_CALL( SCIPvisualUpdateChild(stat->visual, set, stat, tree->focusnode) );

   /* get value of variable in current LP or pseudo solution */
   val = SCIPvarGetSol(var, tree->focusnodehaslp);

   /* avoid branching on infinite values in pseudo solution */
   if( SCIPsetIsInfinity(set, -val) || SCIPsetIsInfinity(set, val) )
   {
      val = SCIPvarGetWorstBoundLocal(var);

      /* if both bounds are infinite, choose zero as branching point */
      if( SCIPsetIsInfinity(set, -val) || SCIPsetIsInfinity(set, val) )
      {
         assert(SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)));
         assert(SCIPsetIsInfinity(set, SCIPvarGetUbLocal(var)));
         val = 0.0;
      }
   }

   assert(SCIPsetIsFeasGE(set, val, SCIPvarGetLbLocal(var)));
   assert(SCIPsetIsFeasLE(set, val, SCIPvarGetUbLocal(var)));
   /* see comment in SCIPbranchVarVal */
   assert(SCIPvarIsIntegral(var));

   /* create child nodes with x <= floor(x'), and x >= ceil(x') */
   downub = floor(val);
   uplb = downub + 1.0;
   SCIPsetDebugMsg(set, "fractional branch on variable <%s> with value %g, root value %g, priority %d (current lower bound: %g)\n",
      SCIPvarGetName(var), val, SCIPvarGetRootSol(var), SCIPvarGetBranchPriority(var), SCIPnodeGetLowerbound(tree->focusnode));

   /* perform the branching;
    * set the node selection priority in a way, s.t. a node is preferred whose branching goes in the same direction
    * as the deviation from the variable's root solution
    */
   if( downub != SCIP_INVALID )    /*lint !e777*/
   {
      /* create child node x <= downub */
      priority = SCIPtreeCalcNodeselPriority(tree, set, stat, var, SCIP_BRANCHDIR_DOWNWARDS, downub);
      /* if LP solution is cutoff in child, compute a new estimate
       * otherwise we cannot expect a direct change in the best solution, so we keep the estimate of the parent node */
      if( SCIPsetIsGT(set, val, downub) )
         estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, downub);
      else
         estimate = SCIPnodeGetEstimate(tree->focusnode);
      SCIPsetDebugMsg(set, " -> creating child: <%s> <= %g (priority: %g, estimate: %g)\n",
         SCIPvarGetName(var), downub, priority, estimate);
      SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, priority, estimate) );

      /* update branching information in certificate, if certificate is active */
      SCIP_CALL( SCIPcertificateUpdateBranchingData(set, stat->certificate, stat, transprob, lp, tree, node, var, SCIP_BOUNDTYPE_UPPER, downub) );

      SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
            eventqueue, eventfilter, NULL, var, downub, SCIP_BOUNDTYPE_UPPER, FALSE) );
      /* output branching bound change to visualization file */
      SCIP_CALL( SCIPvisualUpdateChild(stat->visual, set, stat, node) );

      if( downchild != NULL )
         *downchild = node;
   }

   if( uplb != SCIP_INVALID )    /*lint !e777*/
   {
      /* create child node with x >= uplb */
      priority = SCIPtreeCalcNodeselPriority(tree, set, stat, var, SCIP_BRANCHDIR_UPWARDS, uplb);
      if( SCIPsetIsLT(set, val, uplb) )
         estimate = SCIPtreeCalcChildEstimate(tree, set, stat, var, uplb);
      else
         estimate = SCIPnodeGetEstimate(tree->focusnode);
      SCIPsetDebugMsg(set, " -> creating child: <%s> >= %g (priority: %g, estimate: %g)\n",
         SCIPvarGetName(var), uplb, priority, estimate);
      SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, priority, estimate) );

      /* update branching information in certificate, if certificate is active */
      SCIP_CALL( SCIPcertificateUpdateBranchingData(set, stat->certificate, stat, transprob, lp, tree, node, var, SCIP_BOUNDTYPE_LOWER, uplb) );

      SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
            eventqueue, eventfilter, NULL, var, uplb, SCIP_BOUNDTYPE_LOWER, FALSE) );
      /* output branching bound change to visualization file */
      SCIP_CALL( SCIPvisualUpdateChild(stat->visual, set, stat, node) );

      if( upchild != NULL )
         *upchild = node;
   }

   return SCIP_OKAY;
}

/** calls branching rules to branch on an LP solution; if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *  if the branch priority of an unfixed variable is larger than the maximal branch priority of the fractional
 *  variables, pseudo solution branching is applied on the unfixed variables with maximal branch priority
 */
SCIP_RETCODE SCIPbranchExecLPExact(
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   int i;
   int nalllpcands;  /* sum of binary, integer, and implicit branching candidates */

   assert(branchcand != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* calculate branching candidates */
   SCIP_CALL( branchcandCalcLPCandsExact(branchcand, set, stat, lp) );
   assert(0 <= branchcand->npriolpcands && branchcand->npriolpcands <= branchcand->nlpcands);
   assert((branchcand->npriolpcands == 0) == (branchcand->nlpcands == 0));

   SCIPsetDebugMsg(set, "branching on LP solution with %d (+%d) fractional (+implicit fractional) variables (%d of maximal priority)\n",
      branchcand->nlpcands, branchcand->nimpllpfracs, branchcand->npriolpcands);

   nalllpcands = branchcand->nlpcands + branchcand->nimpllpfracs;
   /* do nothing, if no fractional variables exist */
   if( nalllpcands == 0 )
      return SCIP_OKAY;

   /* if there is a non-fixed variable with higher priority than the maximal priority of the fractional candidates,
    * use pseudo solution branching instead
    */
   if( branchcand->pseudomaxpriority > branchcand->lpmaxpriority )
   {
      SCIP_CALL( SCIPbranchExecPseudo(blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
            eventfilter, cutoffbound, allowaddcons, result) );
      assert(*result != SCIP_DIDNOTRUN && *result != SCIP_DIDNOTFIND);
      return SCIP_OKAY;
   }

   /* it does not make sense to call the normal branching rules, due to assumed very small fractionalities,
    * SCIP is not designed to branch on such values. So we simply branch on the first possible variable
    */
   for( i = 0; i < branchcand->nlpcands && *result != SCIP_BRANCHED; i++ )
   {
      SCIP_VAR* branchvar;

      branchvar = branchcand->lpcands[i];

#ifndef NDEBUG
      {
         SCIP_Rational* tmp;
         SCIP_Real branchval;

         SCIP_CALL( RatCreateBuffer(set->buffer, &tmp) );
         branchval = branchcand->lpcandssol[i];
         RatSetReal(tmp, branchval);
         assert(!RatIsIntegral(tmp));
         RatFreeBuffer(set->buffer, &tmp);
      }
#endif

      SCIP_CALL( SCIPtreeBranchVarExact(tree, reopt, blkmem, set, stat, transprob, origprob, lp,
            branchcand, eventqueue, eventfilter, branchvar, NULL, NULL, NULL) );
      *result = SCIP_BRANCHED;
   }
   /* reset the validlpcandslp to recalculate the branchcands for normal branching */
   branchcand->validlpcandslp = -1;

   return SCIP_OKAY;
}
