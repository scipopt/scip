/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: solve.c,v 1.165 2005/02/03 17:50:45 bzfpfend Exp $"

/**@file   solve.c
 * @brief  main solving loop and node processing
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "def.h"
#include "set.h"
#include "stat.h"
#include "buffer.h"
#include "clock.h"
#include "vbc.h"
#include "interrupt.h"
#include "misc.h"
#include "event.h"
#include "lp.h"
#include "var.h"
#include "prob.h"
#include "sol.h"
#include "primal.h"
#include "tree.h"
#include "pricestore.h"
#include "sepastore.h"
#include "cutpool.h"
#include "solve.h"
#include "scip.h"
#include "branch.h"
#include "conflict.h"
#include "cons.h"
#include "disp.h"
#include "heur.h"
#include "nodesel.h"
#include "pricer.h"
#include "relax.h"
#include "sepa.h"
#include "prop.h"


#define MAXNLPERRORS  10                /**< maximal number of LP error loops in a single node */


/** returns whether the solving process will be / was stopped before proving optimality;
 *  if the solving process was stopped, stores the reason as status in stat
 */
Bool SCIPsolveIsStopped(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   )
{
   if( SCIPinterrupted() )
      stat->status = SCIP_STATUS_USERINTERRUPT;
   else if( set->limit_nodes >= 0 && stat->nnodes >= set->limit_nodes )
      stat->status = SCIP_STATUS_NODELIMIT;
   else if( SCIPclockGetTime(stat->solvingtime) >= set->limit_time )
      stat->status = SCIP_STATUS_TIMELIMIT;
   else if( SCIPgetMemUsed(set->scip) >= set->limit_memory*1024.0*1024.0 )
      stat->status = SCIP_STATUS_MEMLIMIT;
   else if( SCIPgetStage(set->scip) >= SCIP_STAGE_SOLVING && SCIPsetIsLT(set, SCIPgetGap(set->scip), set->limit_gap) )
      stat->status = SCIP_STATUS_GAPLIMIT;
   else if( set->limit_solutions >= 0 && SCIPgetStage(set->scip) >= SCIP_STAGE_PRESOLVED
      && SCIPgetNSolsFound(set->scip) >= set->limit_solutions )
      stat->status = SCIP_STATUS_SOLLIMIT;
   else if( set->limit_bestsol >= 0 && SCIPgetStage(set->scip) >= SCIP_STAGE_PRESOLVED
      && SCIPgetNBestSolsFound(set->scip) >= set->limit_bestsol )
      stat->status = SCIP_STATUS_BESTSOLLIMIT;

   return (stat->status != SCIP_STATUS_UNKNOWN);
}

/** applies domain propagation on current node */
RETCODE SCIPpropagateDomains(
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   int              depth,              /**< depth level to use for propagator frequency checks */
   int              maxproprounds,      /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   NODE* node;
   RESULT result;
   Bool propagain;
   int propround;
   int h;
   int p;

   assert(set != NULL);
   assert(tree != NULL);
   assert(cutoff != NULL);

   node = SCIPtreeGetCurrentNode(tree);
   assert(node != NULL);
   assert(SCIPnodeIsActive(node));
   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_FOCUSNODE
      || SCIPnodeGetType(node) == SCIP_NODETYPE_REFOCUSNODE
      || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);

   /* adjust maximal number of propagation rounds */
   if( maxproprounds == -1 )
      maxproprounds = INT_MAX;
   else if( maxproprounds == 0 )
      maxproprounds = (depth == 0 ? set->prop_maxroundsroot : set->prop_maxrounds);

   debugMessage("domain propagation of node %p in depth %d (using depth %d, maxrounds %d)\n", 
      node, SCIPnodeGetDepth(node), depth, maxproprounds);

   /* propagate as long new bound changes were found and the maximal number of propagation rounds is not exceeded */
   *cutoff = FALSE;
   propround = 0;
   do
   {
      propround++;
      propagain = FALSE;

      /* sort propagators */
      SCIPsetSortProps(set);

      /* call additional propagators with nonnegative priority */
      for( p = 0; p < set->nprops && !(*cutoff); ++p )
      {
         if( SCIPpropGetPriority(set->props[p]) < 0 )
            continue;

         CHECK_OKAY( SCIPpropExec(set->props[p], set, stat, depth, &result) );
         propagain = propagain || (result == SCIP_REDUCEDDOM);
         *cutoff = *cutoff || (result == SCIP_CUTOFF);
      }

      /* propagate constraints */
      for( h = 0; h < set->nconshdlrs && !(*cutoff); ++h )
      {
         CHECK_OKAY( SCIPconshdlrPropagate(set->conshdlrs[h], blkmem, set, stat, prob, depth, &result) );
         propagain = propagain || (result == SCIP_REDUCEDDOM);
         *cutoff = *cutoff || (result == SCIP_CUTOFF);
      }

      /* call additional propagators with negative priority */
      for( p = 0; p < set->nprops && !(*cutoff); ++p )
      {
         if( SCIPpropGetPriority(set->props[p]) >= 0 )
            continue;

         CHECK_OKAY( SCIPpropExec(set->props[p], set, stat, depth, &result) );
         propagain = propagain || (result == SCIP_REDUCEDDOM);
         *cutoff = *cutoff || (result == SCIP_CUTOFF);
      }
   }
   while( propagain && !(*cutoff) && propround < maxproprounds );

   /* mark the node to be completely propagated in the current repropagation subtree level */
   SCIPnodeMarkPropagated(node, tree);

   return SCIP_OKAY;
}

/** strengthens variable's bounds looking at reduced costs */
static
RETCODE redcostStrengthening(
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   NODE* node;
   COL** cols;
   VAR* var;
   int* cstat;
   Real lpobjval;
   Real redcost;
   Real oldlb;
   Real oldub;
   Real newbd;
   Bool strengthen;
   int ncols;
   int c;

   assert(set != NULL);
   assert(tree != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(primal != NULL);
   assert(SCIPsetIsLT(set, SCIPlpGetObjval(lp, set), primal->cutoffbound) );

   /* we cannot apply reduced cost fixing, if we want to solve exactly */
   /**@todo implement reduced cost fixing with interval arithmetics */
   if( set->misc_exactsolve )
      return SCIP_OKAY;

   /* reduced cost strengthening can only be applied, if we have an upper bound on the LP value */
   if( SCIPsetIsInfinity(set, primal->cutoffbound) )
      return SCIP_OKAY;

   /* get LP columns */
   cols = SCIPlpGetCols(lp);
   ncols = SCIPlpGetNCols(lp);
   if( ncols == 0 )
      return SCIP_OKAY;

   stat->nredcoststrcalls++;
   node = SCIPtreeGetFocusNode(tree);
   assert(node != NULL);
   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_FOCUSNODE);

   /* start redcost strengthening timer */
   SCIPclockStart(stat->redcoststrtime, set);

   /* get temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &cstat, ncols) );

   /* get basis status for columns and LP objective value */
   CHECK_OKAY( SCIPlpGetBase(lp, cstat, NULL) );
   lpobjval = SCIPlpGetObjval(lp, set);

   /* check reduced costs for non-basic columns */
   for( c = 0; c < ncols; ++c )
   {
      switch( cstat[c] )
      {
      case SCIP_BASESTAT_LOWER:
         redcost = SCIPcolGetRedcost(cols[c], stat, lp);
         assert(!SCIPsetIsFeasNegative(set, redcost));
         if( SCIPsetIsFeasPositive(set, redcost) )
         {
            var = SCIPcolGetVar(cols[c]);
            oldlb = SCIPvarGetLbLocal(var);
            oldub = SCIPvarGetUbLocal(var);
            assert(SCIPsetIsEQ(set, oldlb, SCIPcolGetLb(cols[c])));
            assert(SCIPsetIsEQ(set, oldub, SCIPcolGetUb(cols[c])));
            if( SCIPsetIsFeasLT(set, oldlb, oldub) )
            {
               /* calculate reduced cost based bound */
               newbd = (primal->cutoffbound - lpobjval) / redcost + oldlb;

               /* check, if new bound is good enough:
                *  - integer variables: take all possible strengthenings
                *  - continuous variables: strengthening must cut part of the variable's dynamic range, and
                *                          at least 20% of the current domain
                */
               if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
               {
                  SCIPvarAdjustUb(var, set, &newbd);
                  strengthen = (newbd < oldub - 0.5);
               }
               else
                  strengthen = (newbd < cols[c]->maxprimsol && newbd <= 0.2 * oldlb + 0.8 * oldub);

               if( strengthen )
               {
                  /* strengthen upper bound */
                  debugMessage("redcost strengthening upper bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
                     SCIPvarGetName(var), oldlb, oldub, oldlb, newbd, primal->cutoffbound, lpobjval, redcost);
                  CHECK_OKAY( SCIPnodeAddBoundchg(node, blkmem, set, stat, tree, lp, branchcand, eventqueue,
                        var, newbd, SCIP_BOUNDTYPE_UPPER, FALSE) );
                  stat->nredcoststrfound++;
               }
            }
         }
         break;

      case SCIP_BASESTAT_BASIC:
         break;

      case SCIP_BASESTAT_UPPER:
         redcost = SCIPcolGetRedcost(cols[c], stat, lp);
         assert(!SCIPsetIsFeasPositive(set, redcost));
         if( SCIPsetIsFeasNegative(set, redcost) )
         {
            var = SCIPcolGetVar(cols[c]);
            oldlb = SCIPvarGetLbLocal(var);
            oldub = SCIPvarGetUbLocal(var);
            assert(SCIPsetIsEQ(set, oldlb, SCIPcolGetLb(cols[c])));
            assert(SCIPsetIsEQ(set, oldub, SCIPcolGetUb(cols[c])));
            if( SCIPsetIsFeasLT(set, oldlb, oldub) )
            {
               /* calculate reduced cost based bound */
               newbd = (primal->cutoffbound - lpobjval) / redcost + oldub;

               /* check, if new bound is good enough:
                *  - integer variables: take all possible strengthenings
                *  - continuous variables: strengthening must cut part of the variable's dynamic range, and
                *                          at least 20% of the current domain
                */
               if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
               {
                  SCIPvarAdjustLb(var, set, &newbd);
                  strengthen = (newbd > oldlb + 0.5);
               }
               else
                  strengthen = (newbd > cols[c]->minprimsol && newbd >= 0.8 * oldlb + 0.2 * oldub);

               /* check, if new bound is good enough: at least 20% strengthening for continuous variables */
               if( strengthen )
               {
                  /* strengthen lower bound */
                  debugMessage("redcost strengthening lower bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
                     SCIPvarGetName(var), oldlb, oldub, newbd, oldub, primal->cutoffbound, lpobjval, redcost);
                  CHECK_OKAY( SCIPnodeAddBoundchg(node, blkmem, set, stat, tree, lp, branchcand, eventqueue,
                        var, newbd, SCIP_BOUNDTYPE_LOWER, FALSE) );
                  stat->nredcoststrfound++;
               }
            }
         }
         break;

      case SCIP_BASESTAT_ZERO:
         assert(SCIPsetIsFeasZero(set, SCIPcolGetRedcost(cols[c], stat, lp)));
         break;

      default:
         errorMessage("invalid basis state\n");
         return SCIP_INVALIDDATA;
      }
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &cstat);

   /* stop redcost strengthening activation timer */
   SCIPclockStop(stat->redcoststrtime, set);

   return SCIP_OKAY;
}

/** returns whether the given variable with the old LP solution value should lead to an update of the pseudo cost entry */
static
Bool isPseudocostUpdateValid(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   Real             oldlpsolval         /**< solution value of variable in old LP */
   )
{
   Real newlpsolval;

   assert(var != NULL);

   /* if the old LP solution value is unknown, the pseudo cost update cannot be performed */
   if( oldlpsolval >= SCIP_INVALID )
      return FALSE;

   /* the bound change on the given variable was responsible for the gain in the dual bound, if the variable's
    * old solution value is outside the current bounds, and the new solution value is equal to the bound
    * closest to the old solution value
    */

   /* find out, which of the current bounds is violated by the old LP solution value */
   if( SCIPsetIsLT(set, oldlpsolval, SCIPvarGetLbLocal(var)) )
   {
      newlpsolval = SCIPvarGetLPSol(var);
      return SCIPsetIsEQ(set, newlpsolval, SCIPvarGetLbLocal(var));
   }
   else if( SCIPsetIsGT(set, oldlpsolval, SCIPvarGetUbLocal(var)) )
   {
      newlpsolval = SCIPvarGetLPSol(var);
      return SCIPsetIsEQ(set, newlpsolval, SCIPvarGetUbLocal(var));
   }
   else
      return FALSE;
}

/** pseudo cost flag stored in the variables to mark them for the pseudo cost update */
enum PseudocostFlag
{
   PSEUDOCOST_NONE     = 0,             /**< variable's bounds were not changed */
   PSEUDOCOST_IGNORE   = 1,             /**< bound changes on variable should be ignored for pseudo cost updates */
   PSEUDOCOST_UPDATE   = 2              /**< pseudo cost value of variable should be updated */
};

/** updates the variable's pseudo cost values after the node's initial LP was solved */
static
RETCODE updatePseudocost(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< LP data */
   )
{
   NODE* focusnode;
   int actdepth;

   assert(lp != NULL);
   assert(tree != NULL);
   assert(tree->path != NULL);

   focusnode = SCIPtreeGetFocusNode(tree);
   assert(SCIPnodeIsActive(focusnode));
   assert(SCIPnodeGetType(focusnode) == SCIP_NODETYPE_FOCUSNODE);
   actdepth = SCIPnodeGetDepth(focusnode);
   assert(tree->path[actdepth] == focusnode);

   if( lp->solved && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && tree->focuslpfork != NULL )
   {
      BOUNDCHG** updates;
      NODE* node;
      VAR* var;
      Real weight;
      Real lpgain;
      int nupdates;
      int nvalidupdates;
      int d;
      int i;

      assert(SCIPnodeIsActive(tree->focuslpfork));
      assert(tree->path[tree->focuslpfork->depth] == tree->focuslpfork);

      /* get a buffer for the collected bound changes; start with a size twice as large as the number of nodes between
       * current node and LP fork
       */
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &updates, 2*(actdepth - tree->focuslpfork->depth)) );
      nupdates = 0;
      nvalidupdates = 0;

      /* search the nodes from LP fork down to current node for bound changes in between; move in this direction, 
       * because the bound changes closer to the LP fork are more likely to have a valid LP solution information
       * attached; collect the bound changes for pseudo cost value updates and mark the corresponding variables such
       * that they are not updated twice in case of more than one bound change on the same variable
       */
      for( d = tree->focuslpfork->depth+1; d <= actdepth; ++d )
      {
         node = tree->path[d];

         if( node->domchg != NULL )
         {
            BOUNDCHG* boundchgs;
            int nboundchgs;

            boundchgs = node->domchg->domchgbound.boundchgs;
            nboundchgs = node->domchg->domchgbound.nboundchgs;
            for( i = 0; i < nboundchgs; ++i )
            {
               if( boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
               {
                  var = boundchgs[i].var;
                  assert(var != NULL);
                  if( var->pseudocostflag == PSEUDOCOST_NONE )
                  {
                     /* remember the bound change and mark the variable */
                     CHECK_OKAY( SCIPsetReallocBufferArray(set, &updates, nupdates+1) );
                     updates[nupdates] = &boundchgs[i];
                     nupdates++;

                     /* check, if the bound change would lead to a valid pseudo cost update */
                     if( isPseudocostUpdateValid(var, set, boundchgs[i].data.branchingdata.lpsolval) )
                     {
                        var->pseudocostflag = PSEUDOCOST_UPDATE;
                        nvalidupdates++;
                     }
                     else
                        var->pseudocostflag = PSEUDOCOST_IGNORE;
                  }
               }
            }
         }
      }

      /* update the pseudo cost values and reset the variables' flags; assume, that the responsibility for the dual gain
       * is equally spread on all bound changes that lead to valid pseudo cost updates
       */
      weight = nvalidupdates > 0 ? 1.0 / (Real)nvalidupdates : 1.0;
      lpgain = SCIPlpGetObjval(lp, set) - tree->focuslpfork->lowerbound;
      lpgain = MAX(lpgain, 0.0);
      for( i = 0; i < nupdates; ++i )
      {
         assert(updates[i]->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING);
         var = updates[i]->var;
         assert(var != NULL);
         assert(var->pseudocostflag != PSEUDOCOST_NONE);
         if( var->pseudocostflag == PSEUDOCOST_UPDATE )
         {
            debugMessage("updating pseudocosts of <%s>: sol: %g -> %g, LP: %e -> %e => gain=%g, weight: %g\n",
               SCIPvarGetName(var), updates[i]->data.branchingdata.lpsolval, SCIPvarGetLPSol(var),
               tree->focuslpfork->lowerbound, SCIPlpGetObjval(lp, set), lpgain, weight);
            CHECK_OKAY( SCIPvarUpdatePseudocost(var, set, stat, 
                  SCIPvarGetLPSol(var) - updates[i]->data.branchingdata.lpsolval, lpgain, weight) );
         }
         var->pseudocostflag = PSEUDOCOST_NONE;
      }

      /* free the buffer for the collected bound changes */
      SCIPsetFreeBufferArray(set, &updates);
   }

   return SCIP_OKAY;
}

/** updates lower bound of node using lower bound of LP */
static
RETCODE nodeUpdateLowerboundLP(
   NODE*            node,               /**< node to set lower bound for */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< LP data */
   )
{
   Real lpobjval;
   
   assert(set != NULL);

   if( set->misc_exactsolve )
   {
      CHECK_OKAY( SCIPlpGetProvedLowerbound(lp, set, &lpobjval) );
   }
   else
      lpobjval = SCIPlpGetObjval(lp, set);

   SCIPnodeUpdateLowerbound(node, stat, lpobjval);

   return SCIP_OKAY;
}

/** constructs the LP of the root node */
static
RETCODE initRootLP(
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICESTORE*      pricestore,         /**< pricing storage */
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   VAR* var;
   int v;
   int h;

   assert(lp != NULL);
   assert(SCIPlpGetNCols(lp) == 0);
   assert(SCIPlpGetNRows(lp) == 0);
   assert(lp->nremoveablecols == 0);
   assert(lp->nremoveablerows == 0);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* inform pricing and separation storage, that LP is now filled with initial data */
   SCIPpricestoreStartInitialLP(pricestore);
   SCIPsepastoreStartInitialLP(sepastore);

   /* add all initial variables to LP */
   debugMessage("init root LP: initial columns\n");
   for( v = 0; v < prob->nvars; ++v )
   {
      var = prob->vars[v];
      assert(SCIPvarGetProbindex(var) >= 0);

      if( SCIPvarIsInitial(var) )
      {
         CHECK_OKAY( SCIPpricestoreAddVar(pricestore, blkmem, set, lp, var, 0.0, TRUE) );
      }
   }
   assert(lp->nremoveablecols == 0);
   CHECK_OKAY( SCIPpricestoreApplyVars(pricestore, blkmem, set, stat, prob, tree, lp) );

   /* add LP relaxations of all initial constraints to LP */
   debugMessage("init root LP: initial rows\n");
   for( h = 0; h < set->nconshdlrs; ++h )
   {
      CHECK_OKAY( SCIPconshdlrInitLP(set->conshdlrs[h], blkmem, set, stat, prob) );
   }
   CHECK_OKAY( SCIPsepastoreApplyCuts(sepastore, blkmem, set, stat, tree, lp, branchcand, eventqueue, cutoff) );

   /* inform pricing and separation storage, that initial LP setup is now finished */
   SCIPpricestoreEndInitialLP(pricestore);
   SCIPsepastoreEndInitialLP(sepastore);

   return SCIP_OKAY;
}

/** load and solve the initial LP of a node */
static
RETCODE solveNodeInitialLP(
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICESTORE*      pricestore,         /**< pricing storage */
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   Bool*            lperror             /**< pointer to store whether an unresolved error in LP solving occured */
   )
{
   NODE* focusnode;
   EVENT event;

   assert(stat != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(cutoff != NULL);
   assert(lperror != NULL);

   *cutoff = FALSE;
   *lperror = FALSE;

   focusnode = SCIPtreeGetFocusNode(tree);
   assert(focusnode != NULL);
   assert(SCIPnodeGetType(focusnode) == SCIP_NODETYPE_FOCUSNODE);

   /* load the LP into the solver and load the LP state */
   debugMessage("loading LP\n");
   CHECK_OKAY( SCIPtreeLoadLP(tree, blkmem, set, stat, lp) );
   
   /* init root node LP */
   if( SCIPnodeGetDepth(focusnode) == 0 )
   {
      CHECK_OKAY( initRootLP(blkmem, set, stat, prob, tree, lp, pricestore, sepastore, branchcand, eventqueue, cutoff) );
   }

   /* only the bounds and the constraint list may have been changed since the LP state node
    * -> dual simplex is applicable as first solver
    */
   debugMessage("node: solve initial LP\n");
   CHECK_OKAY( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, TRUE, FALSE, lperror) );
   assert(lp->flushed);
   assert(lp->solved || *lperror);

   if( !(*lperror) )
   {
      /* issue FIRSTLPSOLVED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_FIRSTLPSOLVED) );
      CHECK_OKAY( SCIPeventChgNode(&event, focusnode) );
      CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );
      
      /* update pseudo cost values */
      CHECK_OKAY( updatePseudocost(set, stat, tree, lp) );
   }

   return SCIP_OKAY;
}

/** solve the current LP of a node with a price-and-cut loop */
static
RETCODE priceAndCutLoop(
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICESTORE*      pricestore,         /**< pricing storage */
   SEPASTORE*       sepastore,          /**< separation storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   CONFLICT*        conflict,           /**< conflict analysis data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   CONSHDLR**       conshdlrs_sepa,     /**< constraint handlers sorted by separation priority */
   Bool             initialloop,        /**< is this the first price-and-cut loop call at the current node? */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   Bool*            unbounded,          /**< pointer to store whether an unbounded ray was found in the LP */
   Bool*            lperror             /**< pointer to store whether an unresolved error in LP solving occured */
   )
{
   NODE* focusnode;
   RESULT result;
   EVENT event;
   Real loclowerbound;
   Real glblowerbound;
   Bool separate;
   Bool mustprice;
   Bool mustsepar;
   Bool root;
   int maxseparounds;
   int nsepastallrounds;
   int maxnsepastallrounds;
   int actdepth;
   int npricedcolvars;
   int h;
   int s;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(pricestore != NULL);
   assert(sepastore != NULL);
   assert(cutpool != NULL);
   assert(primal != NULL);
   assert(set->nconshdlrs == 0 || conshdlrs_sepa != NULL);
   assert(cutoff != NULL);
   assert(unbounded != NULL);
   assert(lperror != NULL);

   focusnode = SCIPtreeGetFocusNode(tree);
   assert(focusnode != NULL);
   assert(SCIPnodeGetType(focusnode) == SCIP_NODETYPE_FOCUSNODE);
   actdepth = SCIPnodeGetDepth(focusnode);
   root = (actdepth == 0);

   /* check, if we want to separate at this node */
   loclowerbound = SCIPnodeGetLowerbound(focusnode);
   glblowerbound = SCIPtreeGetLowerbound(tree, set);
   separate = SCIPsetIsLE(set, loclowerbound - glblowerbound, 
      set->sepa_maxbounddist * (primal->cutoffbound - glblowerbound));

   /* get maximal number of separation rounds */
   maxseparounds = (root ? set->sepa_maxroundsroot : set->sepa_maxrounds);
   if( maxseparounds == -1 )
      maxseparounds = INT_MAX;
   if( initialloop && set->sepa_maxaddrounds >= 0 )
      maxseparounds = MIN(maxseparounds, stat->nseparounds + set->sepa_maxaddrounds);
   maxnsepastallrounds = set->sepa_maxstallrounds;
   if( maxnsepastallrounds == -1 )
      maxnsepastallrounds = INT_MAX;

   /* solve initial LP of price-and-cut loop */
   debugMessage("node: solve LP with price and cut\n");
   CHECK_OKAY( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, TRUE, FALSE, lperror) );
   assert(lp->flushed);
   assert(lp->solved);

   /* price-and-cut loop */
   npricedcolvars = prob->ncolvars;
   mustprice = TRUE;
   mustsepar = separate;
   *cutoff = FALSE;
   *unbounded = FALSE;
   nsepastallrounds = 0;
   while( !(*cutoff) && !(*lperror) && (mustprice || mustsepar) )
   {
      debugMessage("-------- node solving loop --------\n");
      assert(lp->flushed);
      assert(lp->solved);

      /* if the LP is unbounded, we don't need to price */
      mustprice = mustprice && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

      /* if all the variables are already in the LP, we don't need to price */
      mustprice = mustprice && !SCIPprobAllColsInLP(prob, set, lp);

      /* pricing (has to be done completely to get a valid lower bound) */
      while( !(*cutoff) && !(*lperror) && mustprice )
      {
         assert(lp->flushed);
         assert(lp->solved);
         assert(SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

         /* price problem variables */
         debugMessage("problem variable pricing\n");
         assert(SCIPpricestoreGetNVars(pricestore) == 0);
         assert(SCIPpricestoreGetNBoundResets(pricestore) == 0);
         CHECK_OKAY( SCIPpricestoreAddProbVars(pricestore, blkmem, set, stat, prob, tree, lp, branchcand, eventqueue) );
         npricedcolvars = prob->ncolvars;

         /* if no additional problem variables were found, call external pricers to create additional problem variables */
         if( SCIPpricestoreGetNVars(pricestore) == 0 )
         {
            Bool enoughvars;
            int p;

            debugMessage("external variable pricing\n");

            /* sort pricer algorithms by priority */
            SCIPsetSortPricers(set);
            
            /* call external pricer algorithms, that are active for the current problem */
            enoughvars = FALSE;
            for( p = 0; p < set->nactivepricers && !enoughvars; ++p )
            {
               CHECK_OKAY( SCIPpricerExec(set->pricers[p], set, prob, lp) );
               enoughvars = enoughvars || (SCIPpricestoreGetNVars(pricestore) >= SCIPsetGetPriceMaxvars(set, root)/2);
            }
         }

         /* apply the priced variables to the LP */
         CHECK_OKAY( SCIPpricestoreApplyVars(pricestore, blkmem, set, stat, prob, tree, lp) );
         assert(SCIPpricestoreGetNVars(pricestore) == 0);
         assert(!lp->flushed || lp->solved);
         mustprice = !lp->flushed || (prob->ncolvars != npricedcolvars);
         mustsepar = mustsepar || !lp->flushed;
         
         /* after adding columns, the LP should be primal feasible such that primal simplex is applicable;
          * if LP was infeasible, we have to use dual simplex
          */
         debugMessage("pricing: solve LP\n");
         CHECK_OKAY( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, TRUE, FALSE, lperror) );
         assert(lp->flushed);
         assert(lp->solved || *lperror);

         /* reset bounds temporarily set by pricer to their original values */
         debugMessage("pricing: reset bounds\n");
         CHECK_OKAY( SCIPpricestoreResetBounds(pricestore, blkmem, set, stat, lp, branchcand, eventqueue) );
         assert(SCIPpricestoreGetNVars(pricestore) == 0);
         assert(SCIPpricestoreGetNBoundResets(pricestore) == 0);
         assert(!lp->flushed || lp->solved);
         mustprice = mustprice || !lp->flushed || (prob->ncolvars != npricedcolvars);
         mustsepar = mustsepar || !lp->flushed;

         /* solve LP again after resetting bounds (with dual simplex) */
         debugMessage("pricing: solve LP after resetting bounds\n");
         CHECK_OKAY( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, FALSE, FALSE, lperror) );
         assert(lp->flushed);
         assert(lp->solved || *lperror);

         /* increase pricing round counter */
         stat->npricerounds++;

         /* display node information line for root node */
         if( root && mustprice && (VERBLEVEL)set->disp_verblevel >= SCIP_VERBLEVEL_FULL )
         {
            CHECK_OKAY( SCIPdispPrintLine(set, stat, NULL, TRUE) );
         }

         /* if the LP is unbounded, we can stop pricing */
         mustprice = mustprice && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY);
      }
      assert(lp->flushed);
      assert(lp->solved);

      /* update lower bound w.r.t. the the LP solution */
      if( !(*lperror) )
      {
         CHECK_OKAY( nodeUpdateLowerboundLP(focusnode, set, stat, lp) );
         debugMessage(" -> new lower bound: %g (LP status: %d, LP obj: %g)\n",
            SCIPnodeGetLowerbound(focusnode), SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set));
      }
      else
      {
         debugMessage(" -> error solving LP. keeping old bound: %g\n", SCIPnodeGetLowerbound(focusnode));
      }

      /* display node information line for root node */
      if( root && (VERBLEVEL)set->disp_verblevel >= SCIP_VERBLEVEL_HIGH )
      {
         CHECK_OKAY( SCIPdispPrintLine(set, stat, NULL, TRUE) );
      }

      /* if the LP is infeasible or exceeded the objective limit, we don't need to separate cuts */
      mustsepar = mustsepar
         && separate
         && stat->nseparounds < maxseparounds
         && nsepastallrounds < maxnsepastallrounds
         && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_INFEASIBLE)
         && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OBJLIMIT)
         && SCIPsetIsLT(set, SCIPnodeGetLowerbound(focusnode), primal->cutoffbound);

      /* separation and reduced cost strengthening
       * (needs not to be done completely, because we just want to increase the lower bound)
       */
      if( !(*cutoff) && !(*lperror) && mustsepar )
      {
         Longint olddomchgcount;
         Real oldlpobjval;
         Bool separateagain;
         Bool enoughcuts;

         assert(lp->flushed);
         assert(lp->solved);
         assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);

         olddomchgcount = stat->domchgcount;
         oldlpobjval = SCIPlpGetObjval(lp, set);

         mustsepar = FALSE;
         enoughcuts = (SCIPsetGetSepaMaxcuts(set, root) == 0);

         /* global cut pool separation */
         if( !enoughcuts )
         {
            if( (set->sepa_poolfreq == 0 && actdepth == 0)
               || (set->sepa_poolfreq > 0 && actdepth % set->sepa_poolfreq == 0) )
            {
               debugMessage("global cut pool separation\n");
               assert(SCIPsepastoreGetNCuts(sepastore) == 0);
               CHECK_OKAY( SCIPcutpoolSeparate(cutpool, blkmem, set, stat, lp, sepastore, root, &result) );
               *cutoff = *cutoff || (result == SCIP_CUTOFF);
               enoughcuts = enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
            }
         }
         assert(lp->flushed);
         assert(lp->solved);
         assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);

         /* constraint separation */
         debugMessage("constraint separation\n");

         /* separate constraints and LP */
         separateagain = TRUE;
         while( !(*cutoff) && !(*lperror) && !enoughcuts && separateagain
            && lp->solved
            && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY) )
         {
            separateagain = FALSE;

            /* sort separators by priority */
            SCIPsetSortSepas(set);

            /* call LP separators with nonnegative priority */
            for( s = 0; s < set->nsepas && !(*cutoff) && !separateagain && !enoughcuts && lp->flushed; ++s )
            {
               if( SCIPsepaGetPriority(set->sepas[s]) < 0 )
                  continue;

               CHECK_OKAY( SCIPsepaExec(set->sepas[s], set, stat, sepastore, actdepth, &result) );
               *cutoff = *cutoff || (result == SCIP_CUTOFF);
               separateagain = separateagain || (result == SCIP_CONSADDED);
               enoughcuts = enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
            }

            /* try separating constraints of the constraint handlers until the cut pool is at least half full;
             * we have to stop after a domain reduction was found, because they go directly into the LP and invalidate
             * the current solution
             */
            for( h = 0; h < set->nconshdlrs && !(*cutoff) && !enoughcuts && lp->flushed; ++h )
            {
               CHECK_OKAY( SCIPconshdlrSeparate(conshdlrs_sepa[h], blkmem, set, stat, prob, sepastore, actdepth,
                     &result) );
               *cutoff = *cutoff || (result == SCIP_CUTOFF);
               separateagain = separateagain || (result == SCIP_CONSADDED);
               enoughcuts = enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
            }

            /* call LP separators with negative priority */
            for( s = 0; s < set->nsepas && !(*cutoff) && !separateagain && !enoughcuts && lp->flushed; ++s )
            {
               if( SCIPsepaGetPriority(set->sepas[s]) >= 0 )
                  continue;

               CHECK_OKAY( SCIPsepaExec(set->sepas[s], set, stat, sepastore, actdepth, &result) );
               *cutoff = *cutoff || (result == SCIP_CUTOFF);
               separateagain = separateagain || (result == SCIP_CONSADDED);
               enoughcuts = enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
            }

            if( !(*cutoff) && !lp->flushed )
            {
               /* solve LP (with dual simplex) */
               debugMessage("separation: solve LP\n");
               CHECK_OKAY( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, TRUE, FALSE, lperror) );
               assert(lp->flushed);
               assert(lp->solved || lperror);
               separateagain = TRUE;
            }
         }
         assert(*cutoff || *lperror || (lp->flushed && lp->solved));
         assert(!lp->solved
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT);

         if( *cutoff || *lperror
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT )
         {
            /* the found cuts are of no use, because the node is infeasible anyway (or we have an error in the LP) */
            CHECK_OKAY( SCIPsepastoreClearCuts(sepastore, blkmem, set, lp) );
         }
         else
         {
            /* apply reduced cost bound strengthening */
            if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && !mustprice && prob->ncolvars == npricedcolvars )
            {
               if( (set->prop_redcostfreq == 0 && root)
                  || (set->prop_redcostfreq > 0 && actdepth % set->prop_redcostfreq == 0) )
               {
                  CHECK_OKAY( redcostStrengthening(blkmem, set, stat, primal, tree, lp, branchcand, eventqueue) );
               }
            }
        
            /* apply found cuts */
            CHECK_OKAY( SCIPsepastoreApplyCuts(sepastore, blkmem, set, stat, tree, lp, branchcand, eventqueue, cutoff) );

            if( !(*cutoff) )
            {
               /* if a new bound change (e.g. a cut with only one column) was found, propagate domains again */
               if( stat->domchgcount != olddomchgcount )
               {
                  CHECK_OKAY( SCIPpropagateDomains(blkmem, set, stat, prob, tree, 
                        SCIPtreeGetCurrentDepth(tree), 0, cutoff) );
               }
               
               mustprice = mustprice || !lp->flushed || (prob->ncolvars != npricedcolvars);
               mustsepar = !lp->flushed;
               
               if( !(*cutoff) )
               {
                  Real lpobjval;

                  /* solve LP (with dual simplex) */
                  debugMessage("separation: solve LP\n");
                  CHECK_OKAY( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, TRUE, FALSE, lperror) );
                  assert(lp->flushed);
                  assert(lp->solved || lperror);

                  lpobjval = SCIPlpGetObjval(lp, set);
                  if( SCIPsetIsFeasGT(set, lpobjval, oldlpobjval) )
                     nsepastallrounds = 0;
                  else
                     nsepastallrounds++;
               }
            }
         }
         assert(*cutoff || *lperror || (lp->flushed && lp->solved)); /* cutoff: LP may be unsolved due to bound changes */

         /* increase separation round counter */
         stat->nseparounds++;
      }
   }

   /* update lower bound w.r.t. the the LP solution */
   if( *cutoff )
   {
      SCIPnodeUpdateLowerbound(focusnode, stat, primal->cutoffbound);
   }
   else if( !(*lperror) )
   {
      assert(lp->flushed);
      assert(lp->solved);

      CHECK_OKAY( nodeUpdateLowerboundLP(focusnode, set, stat, lp) );

      /* issue LPSOLVED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_LPSOLVED) );
      CHECK_OKAY( SCIPeventChgNode(&event, focusnode) );
      CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );

      /* analyze an infeasible LP (not necessary in the root node) */
      if( !set->misc_exactsolve && !root
         && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT) )
      {
         CHECK_OKAY( SCIPconflictAnalyzeLP(conflict, blkmem, set, stat, prob, tree, lp, NULL) );
      }

      /* check for unboundness */
      if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
      {
         assert(root); /* this can only happen in the root node */
         *unbounded = TRUE;
      }
   }
   debugMessage(" -> final lower bound: %g (LP status: %d, LP obj: %g)\n",
      SCIPnodeGetLowerbound(focusnode), SCIPlpGetSolstat(lp), *cutoff ? SCIPsetInfinity(set) : SCIPlpGetObjval(lp, set));

   return SCIP_OKAY;
}

/** solves the current node's LP in a price-and-cut loop */
static
RETCODE solveNodeLP(
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICESTORE*      pricestore,         /**< pricing storage */
   SEPASTORE*       sepastore,          /**< separation storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   CONFLICT*        conflict,           /**< conflict analysis data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   CONSHDLR**       conshdlrs_sepa,     /**< constraint handlers for separating constraints, sorted by priority */
   Bool             initiallpsolved,    /**< was the initial LP already solved? */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            unbounded,          /**< pointer to store TRUE, if an unbounded ray was found in the LP */
   Bool*            lperror             /**< pointer to store TRUE, if an unresolved error in LP solving occured */
   )
{
   Longint nlpiterations;
   int nlps;

   assert(stat != NULL);
   assert(tree != NULL);
   assert(SCIPtreeHasFocusNodeLP(tree));
   assert(cutoff != NULL);
   assert(lperror != NULL);
   assert(cutoff != NULL);
   assert(unbounded != NULL);
   assert(lperror != NULL);
   assert(*cutoff == FALSE);
   assert(*unbounded == FALSE);
   assert(*lperror == FALSE);

   nlps = stat->nlps;
   nlpiterations = stat->nlpiterations;

   if( !initiallpsolved )
   {
      /* load and solve the initial LP of the node */
      CHECK_OKAY( solveNodeInitialLP(blkmem, set, stat, prob, tree, lp, pricestore, sepastore, 
            branchcand, eventfilter, eventqueue, cutoff, lperror) );
      assert(*cutoff || *lperror || (lp->flushed && lp->solved));
      debugMessage("price-and-cut-loop: initial LP status: %d, LP obj: %g\n", 
         SCIPlpGetSolstat(lp), *cutoff ? SCIPsetInfinity(set) : SCIPlpGetObjval(lp, set));

      /* update initial LP iteration counter */
      stat->ninitlps += stat->nlps - nlps;
      stat->ninitlpiterations += stat->nlpiterations - nlpiterations;
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);

   if( !(*cutoff) && !(*lperror) )
   {
      /* solve the LP with price-and-cut*/
      CHECK_OKAY( priceAndCutLoop(blkmem, set, stat, prob, primal, tree, lp, pricestore, sepastore, cutpool, 
            branchcand, conflict, eventfilter, eventqueue, conshdlrs_sepa, initiallpsolved, cutoff, unbounded, lperror) );
   }
   assert(*cutoff || *lperror || (lp->flushed && lp->solved));

   /* update node's LP iteration counter */
   stat->nnodelps += stat->nlps - nlps;
   stat->nnodelpiterations += stat->nlpiterations - nlpiterations;

#if 0 /**@todo check if this is valid and useful (faster (but not so well?) strong branching -> see 10teams.mps) */
   /*???????????????????????????????*/
   if( !(*cutoff) && !(*lperror) && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      /* clean up newly created part of LP to keep only necessary columns and rows */
      CHECK_OKAY( SCIPlpCleanupNew(lp, blkmem, set, stat, SCIPtreeGetFocusDepth(tree) == 0) );
      
      /* resolve LP after cleaning up */
      if( !lp->solved || !lp->flushed )
      {
         debugMessage("resolving LP after cleanup\n");
         CHECK_OKAY( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, FALSE, TRUE, lperror) );
      }
   }
#endif
   
   return SCIP_OKAY;
}

/** calls relaxators */
static
RETCODE solveNodeRelax(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   int              depth,              /**< depth of current node */
   Bool             beforelp,           /**< should the relaxators with non-negative or negative priority be called? */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            propagateagain,     /**< pointer to store TRUE, if domain propagation should be applied again */
   Bool*            solvelpagain        /**< pointer to store TRUE, if the node's LP has to be solved again */
   )
{
   RESULT result;
   int r;

   assert(set != NULL);
   assert(cutoff != NULL);
   assert(solvelpagain != NULL);
   assert(propagateagain != NULL);
   assert(!(*cutoff));

   for( r = 0; r < set->nrelaxs && !(*cutoff); ++r )
   {
      if( beforelp != (SCIPrelaxGetPriority(set->relaxs[r]) >= 0) )
         continue;

      CHECK_OKAY( SCIPrelaxExec(set->relaxs[r], set, stat, depth, &result) );

      switch( result )
      {  
      case SCIP_CUTOFF:
         *cutoff = TRUE;
         break;

      case SCIP_CONSADDED:
         *solvelpagain = TRUE;   /* the separation for new constraints should be called */
         *propagateagain = TRUE; /* the propagation for new constraints should be called */
         break;

      case SCIP_REDUCEDDOM:
         *solvelpagain = TRUE;
         *propagateagain = TRUE;
         break;

      case SCIP_SEPARATED:
         *solvelpagain = TRUE;
         break;

      case SCIP_SUCCESS:
      case SCIP_SUSPENDED:
      case SCIP_DIDNOTRUN:
         break;

      default:
         errorMessage("invalid result code <%d> of relaxator <%s>\n", result, SCIPrelaxGetName(set->relaxs[r]));
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** enforces constraints by branching, separation, or domain reduction */
static
RETCODE enforceConstraints(
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   SEPASTORE*       sepastore,          /**< separation storage */
   CONSHDLR**       conshdlrs_enfo,     /**< constraint handlers for enforcing constraints, sorted by priority */
   Bool*            branched,           /**< pointer to store whether a branching was created */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            infeasible,         /**< pointer to store TRUE, if the LP/pseudo solution is infeasible */
   Bool*            propagateagain,     /**< pointer to store TRUE, if domain propagation should be applied again */
   Bool*            solvelpagain        /**< pointer to store TRUE, if the node's LP has to be solved again */
   )
{
   RESULT result;
   Real pseudoobjval;
   Bool resolved;
   Bool objinfeasible;
   int h;

   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(SCIPtreeGetFocusNode(tree) != NULL);
   assert(prob != NULL);
   assert(conshdlrs_enfo != NULL);
   assert(branched != NULL);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(propagateagain != NULL);
   assert(solvelpagain != NULL);
   assert(!(*cutoff));
   assert(!(*propagateagain));
   assert(!(*solvelpagain));

   *branched = FALSE;
   /**@todo avoid checking the same pseudosolution twice */

   /* enforce constraints by branching, applying additional cutting planes (if LP is being processed),
    * introducing new constraints, or tighten the domains
    */
   debugMessage("enforcing constraints on %s solution\n", SCIPtreeHasFocusNodeLP(tree) ? "LP" : "pseudo");

   /* check, if the solution is infeasible anyway due to it's objective value */
   if( SCIPtreeHasFocusNodeLP(tree) )
      objinfeasible = FALSE;
   else
   {
      pseudoobjval = SCIPlpGetPseudoObjval(lp, set);
      objinfeasible = SCIPsetIsLT(set, pseudoobjval, SCIPnodeGetLowerbound(SCIPtreeGetFocusNode(tree)));
   }

   /* enforce constraints until a handler resolved an infeasibility with cutting off the node, branching, 
    * reducing a domain, or separating a cut
    * if a constraint handler introduced new constraints to enforce his constraints, the newly added constraints
    * have to be enforced themselves
    */
   resolved = FALSE;
   for( h = 0; h < set->nconshdlrs && !resolved; ++h )
   {
      assert(SCIPsepastoreGetNCuts(sepastore) == 0); /* otherwise, the LP should have been resolved first */

      if( SCIPtreeHasFocusNodeLP(tree) )
      {
         assert(lp->flushed);
         assert(lp->solved);
         assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
         CHECK_OKAY( SCIPconshdlrEnforceLPSol(conshdlrs_enfo[h], blkmem, set, stat, prob, tree, sepastore, &result) );
      }
      else
      {
         CHECK_OKAY( SCIPconshdlrEnforcePseudoSol(conshdlrs_enfo[h], blkmem, set, stat, prob, tree, objinfeasible, 
               &result) );
         if( SCIPsepastoreGetNCuts(sepastore) != 0 )
         {
            errorMessage("pseudo enforcing method of constraint handler <%s> separated cuts\n",
               SCIPconshdlrGetName(conshdlrs_enfo[h]));
            return SCIP_INVALIDRESULT;
         }
      }
      debugMessage("enforcing of <%s> returned result %d\n", SCIPconshdlrGetName(conshdlrs_enfo[h]), result);

      switch( result )
      {  
      case SCIP_CUTOFF:
         assert(tree->nchildren == 0);
         *cutoff = TRUE;
         *infeasible = TRUE;
         resolved = TRUE;
         break;

      case SCIP_CONSADDED:
         assert(tree->nchildren == 0);
         *infeasible = TRUE;
         *solvelpagain = TRUE;   /* the separation for new constraints should be called */
         *propagateagain = TRUE; /* the propagation for new constraints should be called */
         resolved = TRUE;
         break;

      case SCIP_REDUCEDDOM:
         assert(tree->nchildren == 0);
         *infeasible = TRUE;
         *solvelpagain = TRUE;
         *propagateagain = TRUE;
         resolved = TRUE;
         break;

      case SCIP_SEPARATED:
         assert(tree->nchildren == 0);
         assert(SCIPsepastoreGetNCuts(sepastore) > 0);
         *infeasible = TRUE;
         *solvelpagain = TRUE;
         resolved = TRUE;
         break;

      case SCIP_BRANCHED:
         assert(tree->nchildren >= 1);
         assert(!SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         *infeasible = TRUE;
         *branched = TRUE;
         resolved = TRUE;
         break;

      case SCIP_SOLVELP:
         assert(!SCIPtreeHasFocusNodeLP(tree));
         assert(tree->nchildren == 0);
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         *infeasible = TRUE;
         *solvelpagain = TRUE;
         resolved = TRUE;
         SCIPtreeSetFocusNodeLP(tree, TRUE); /* the node's LP must be solved */
         break;

      case SCIP_INFEASIBLE:
         assert(tree->nchildren == 0);
         assert(!SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         *infeasible = TRUE;
         break;

      case SCIP_FEASIBLE:
         assert(tree->nchildren == 0);
         assert(!SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         break;

      case SCIP_DIDNOTRUN:
         assert(tree->nchildren == 0);
         assert(!SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         assert(objinfeasible);
         *infeasible = TRUE;
         break;

      default:
         errorMessage("invalid result code <%d> from enforcing method of constraint handler <%s>\n",
            result, SCIPconshdlrGetName(conshdlrs_enfo[h]));
         return SCIP_INVALIDRESULT;
      }  /*lint !e788*/

      /* the enforcement method may add a primal solution, after which the LP status could be set to
       * objective limit reached
       */
      if( SCIPtreeHasFocusNodeLP(tree) && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT )
      {
         *cutoff = TRUE;
         *infeasible = TRUE;
         resolved = TRUE;
      }

      assert(!(*branched) || (resolved && !(*cutoff) && *infeasible && !(*propagateagain) && !(*solvelpagain)));
      assert(!(*cutoff) || (resolved && !(*branched) && *infeasible && !(*propagateagain) && !(*solvelpagain)));
      assert(*infeasible || (!resolved && !(*branched) && !(*cutoff) && !(*propagateagain) && !(*solvelpagain)));
      assert(!(*propagateagain) || (resolved && !(*branched) && !(*cutoff) && *infeasible));
      assert(!(*solvelpagain) || (resolved && !(*branched) && !(*cutoff) && *infeasible));
   }
   assert(!objinfeasible || *infeasible);
   assert(resolved == (*branched || *cutoff || *propagateagain || *solvelpagain));
   assert(*cutoff || *solvelpagain || SCIPsepastoreGetNCuts(sepastore) == 0);

   debugMessage(" -> enforcing result: branched=%d, cutoff=%d, infeasible=%d, propagateagain=%d, solvelpagain=%d, resolved=%d\n",
      *branched, *cutoff, *infeasible, *propagateagain, *solvelpagain, resolved);

   return SCIP_OKAY;
}

/** applies the cuts stored in the separation store, or clears the store if the node can be cut off */
static
RETCODE applyCuts(
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            cutoff,             /**< pointer to whether the node can be cut off */
   Bool*            propagateagain,     /**< pointer to store TRUE, if domain propagation should be applied again */
   Bool*            solvelpagain        /**< pointer to store TRUE, if the node's LP has to be solved again */
   )
{
   assert(stat != NULL);
   assert(cutoff != NULL);
   assert(propagateagain != NULL);
   assert(solvelpagain != NULL);

   if( *cutoff )
   {
      /* the found cuts are of no use, because the node is infeasible anyway (or we have an error in the LP) */
      CHECK_OKAY( SCIPsepastoreClearCuts(sepastore, blkmem, set, lp) );
   }
   else if( SCIPsepastoreGetNCuts(sepastore) > 0 )
   {
      Longint olddomchgcount;

      olddomchgcount = stat->domchgcount;
      CHECK_OKAY( SCIPsepastoreApplyCuts(sepastore, blkmem, set, stat, tree, lp, branchcand, eventqueue, cutoff) );
      *propagateagain = *propagateagain || (stat->domchgcount != olddomchgcount);
      *solvelpagain = TRUE;
   }

   return SCIP_OKAY;
}

/** updates the cutoff, propagateagain, and solverelaxagain status of the current solving loop */
static
void updateLoopStatus(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   int              depth,              /**< depth of current node */
   Bool*            cutoff,             /**< pointer to whether the node can be cut off */
   Bool*            propagateagain,     /**< pointer to store TRUE, if domain propagation should be applied again */
   Bool*            solverelaxagain     /**< pointer to store TRUE, if at least one relaxator should be called again */
   )
{
   NODE* focusnode;
   int r;

   assert(set != NULL);
   assert(stat != NULL);
   assert(cutoff != NULL);
   assert(propagateagain != NULL);
   assert(solverelaxagain != NULL);

   /* check, if the path was cutoff */
   *cutoff = *cutoff || (tree->cutoffdepth <= depth);
   
   /* check, if the focus node should be repropagated */
   focusnode = SCIPtreeGetFocusNode(tree);
   *propagateagain = *propagateagain || SCIPnodeIsPropagatedAgain(focusnode);

   /* check, if one of the external relaxations should be solved again */
   for( r = 0; r < set->nrelaxs && !(*solverelaxagain); ++r )
      *solverelaxagain = !SCIPrelaxIsSolved(set->relaxs[r], stat);
}

/** solves the focus node */
static
RETCODE solveNode(
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICESTORE*      pricestore,         /**< pricing storage */
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   CONFLICT*        conflict,           /**< conflict analysis data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   CONSHDLR**       conshdlrs_sepa,     /**< constraint handlers for separating constraints, sorted by priority */
   CONSHDLR**       conshdlrs_enfo,     /**< constraint handlers for enforcing constraints, sorted by priority */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   Bool*            unbounded,          /**< pointer to store whether the focus node is unbounded */
   Bool*            infeasible,         /**< pointer to store whether the focus node's solution is infeasible */
   Bool*            restart             /**< should solving process be started again with presolving? */
   )
{
   NODE* focusnode;
   Longint lastdomchgcount;
   int lastlpcount;
   int actdepth;
   int nlperrors;
   Bool focusnodehaslp;
   Bool initiallpsolved;
   Bool solverelaxagain;
   Bool solvelpagain;
   Bool propagateagain;
   Bool branched;

   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);
   assert(primal != NULL);
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);
   assert(cutoff != NULL);
   assert(unbounded != NULL);
   assert(infeasible != NULL);
   assert(restart != NULL);

   *cutoff = FALSE;
   *unbounded = FALSE;
   *infeasible = FALSE;
   *restart = FALSE;

   focusnode = SCIPtreeGetFocusNode(tree);
   assert(focusnode != NULL);
   assert(SCIPnodeGetType(focusnode) == SCIP_NODETYPE_FOCUSNODE);
   actdepth = SCIPnodeGetDepth(focusnode);

   debugMessage("Processing node %lld in depth %d, %d siblings\n", stat->nnodes, actdepth, tree->nsiblings);
   debugMessage("current pseudosolution: obj=%g", SCIPlpGetPseudoObjval(lp, set));
   debug(SCIPprobPrintPseudoSol(prob, set));
   
   /* check, if we want to solve the LP at the selected node:
    * - solve the LP, if the lp solve depth and frequency demand solving
    * - solve the root LP, if the LP solve frequency is set to 0
    * - solve the root LP, if there are continuous variables present
    * - don't solve the node if its cut off by the pseudo objective value anyway
    */
   focusnodehaslp = (set->lp_solvedepth == -1 || actdepth <= set->lp_solvedepth);
   focusnodehaslp = focusnodehaslp && (set->lp_solvefreq >= 1 && actdepth % set->lp_solvefreq == 0);
   focusnodehaslp = focusnodehaslp || (actdepth == 0 && set->lp_solvefreq == 0);
   focusnodehaslp = focusnodehaslp || (actdepth == 0 && prob->ncontvars > 0);
   focusnodehaslp = focusnodehaslp && SCIPsetIsLT(set, SCIPlpGetPseudoObjval(lp, set), primal->cutoffbound);
   SCIPtreeSetFocusNodeLP(tree, focusnodehaslp);

   /* external node solving loop:
    *  - propagate domains
    *  - solve LP
    *  - enforce constraints
    * if a constraint handler adds constraints to enforce its own constraints, both, propagation and LP solving
    * is applied again (if applicable on current node); however, if the new constraints don't have the enforce flag set,
    * it is possible, that the current infeasible solution is not cut off; in this case, we have to declare the solution
    * infeasible and perform a branching
    */
   lastdomchgcount = stat->domchgcount;
   lastlpcount = stat->lpcount;
   initiallpsolved = FALSE;
   nlperrors = 0;
   stat->npricerounds = 0;
   stat->nseparounds = 0;
   solverelaxagain = TRUE;
   solvelpagain = TRUE;
   propagateagain = TRUE;
   while( !(*cutoff) && (solverelaxagain || solvelpagain || propagateagain) && nlperrors < MAXNLPERRORS && !(*restart) )
   {
      Real pseudoobjval;
      Bool lperror;
      int r;

      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      lperror = FALSE;

      /* domain propagation */
      if( propagateagain && !(*cutoff) )
      {
         propagateagain = FALSE;
         CHECK_OKAY( SCIPpropagateDomains(blkmem, set, stat, prob, tree, SCIPtreeGetCurrentDepth(tree), 0, cutoff) );

         /* check, if the path was cutoff */
         *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);
      }
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      /* if the LP should be resolved, all relaxations should also be resolved */
      /**@todo if the LP modification methods of the relax interface is implemented, we can remove this and give
       *       total control to the relaxators
       */
      if( solvelpagain )
      {
         solverelaxagain = TRUE;
         for( r = 0; r < set->nrelaxs; ++r )
            SCIPrelaxMarkUnsolved(set->relaxs[r]);
      }

      /* solve external relaxations with non-negative priority */
      if( solverelaxagain && !(*cutoff) )
      {
         CHECK_OKAY( solveNodeRelax(set, stat, actdepth, TRUE, cutoff, &propagateagain, &solvelpagain) );

         /* check, if the path was cutoff */
         *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);

         /* apply found cuts */
         CHECK_OKAY( applyCuts(blkmem, set, stat, tree, lp, sepastore, branchcand, eventqueue, 
               cutoff, &propagateagain, &solvelpagain) );
      }
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      /* check, if we want to solve the LP at this node */
      if( solvelpagain && !(*cutoff) )
      {
         solvelpagain = FALSE;
         if( SCIPtreeHasFocusNodeLP(tree) )
         {
            /* solve the node's LP */
            CHECK_OKAY( solveNodeLP(blkmem, set, stat, prob, primal, tree, lp, pricestore, sepastore,
                  cutpool, branchcand, conflict, eventfilter, eventqueue, conshdlrs_sepa,
                  initiallpsolved, cutoff, unbounded, &lperror) );
            initiallpsolved = TRUE;
            debugMessage(" -> LP status: %d, LP obj: %g, iter: %lld, count: %d\n", 
               SCIPlpGetSolstat(lp), *cutoff ? SCIPsetInfinity(set) : SCIPlpGetObjval(lp, set),
               stat->nlpiterations, stat->lpcount);
            
            /* check, if the path was cutoff */
            *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);
            
            /* if an error occured during LP solving, switch to pseudo solution */
            if( lperror )
            {
               SCIPtreeSetFocusNodeLP(tree, FALSE);
               nlperrors++;
               infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) unresolved numerical troubles in LP %d -- using pseudo solution instead (loop %d)\n", 
                  stat->nnodes, stat->nlps, nlperrors);
            }
               
            /* if we solve exactly, the LP claims to be infeasible but the infeasibility could not be proved,
             * we have to forget about the LP and use the pseudo solution instead
             */
            if( !(*cutoff) && set->misc_exactsolve && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE
               && SCIPnodeGetLowerbound(focusnode) < primal->cutoffbound )
            {
               if( SCIPbranchcandGetNPseudoCands(branchcand) == 0 && prob->ncontvars > 0 )
               {
                  errorMessage("(node %lld) could not prove infeasibility of LP %d, all variables are fixed, %d continuous vars\n",
                     stat->nnodes, stat->nlps, prob->ncontvars);
                  errorMessage("(node %lld)  -> have to call PerPlex() (feature not yet implemented)\n", stat->nnodes);
                  /**@todo call PerPlex */
                  return SCIP_LPERROR;
               }
               else
               {
                  SCIPtreeSetFocusNodeLP(tree, FALSE);
                  infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                     "(node %lld) could not prove infeasibility of LP %d -- using pseudo solution (%d unfixed vars) instead\n", 
                     stat->nnodes, stat->nlps, SCIPbranchcandGetNPseudoCands(branchcand));
               }
            }
         }
      }
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);
      assert(*cutoff || !SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));

      /* solve external relaxations with negative priority */
      if( solverelaxagain && !(*cutoff) )
      {
         solverelaxagain = FALSE;
         CHECK_OKAY( solveNodeRelax(set, stat, actdepth, FALSE, cutoff, &propagateagain, &solvelpagain) );

         /* check, if the path was cutoff */
         *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);

         /* apply found cuts */
         CHECK_OKAY( applyCuts(blkmem, set, stat, tree, lp, sepastore, branchcand, eventqueue, 
               cutoff, &propagateagain, &solvelpagain) );
      }
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      /* update lower bound w.r.t. the pseudo solution */
      pseudoobjval = SCIPlpGetPseudoObjval(lp, set);
      SCIPnodeUpdateLowerbound(focusnode, stat, pseudoobjval);
      debugMessage(" -> new lower bound: %g [%g] (pseudoobj: %g [%g])\n",
         SCIPnodeGetLowerbound(focusnode), SCIPprobExternObjval(prob, set, SCIPnodeGetLowerbound(focusnode)),
         pseudoobjval, SCIPprobExternObjval(prob, set, pseudoobjval));

      /* check for infeasible node by bounding */
      if( !(*cutoff)
         && (SCIPnodeGetLowerbound(focusnode) >= primal->cutoffbound
            || (!set->misc_exactsolve && SCIPsetIsGE(set, SCIPnodeGetLowerbound(focusnode), primal->cutoffbound))) )
      {
         debugMessage("node is cut off by bounding (lower=%g, upper=%g)\n",
            SCIPnodeGetLowerbound(focusnode), primal->cutoffbound);
         SCIPnodeUpdateLowerbound(focusnode, stat, primal->cutoffbound);
         *cutoff = TRUE;

         /* call pseudo conflict analysis, if the node is cut off due to the pseudo objective value */
         if( pseudoobjval >= primal->cutoffbound )
         {
            CHECK_OKAY( SCIPconflictAnalyzePseudo(conflict, blkmem, set, stat, prob, tree, lp, NULL) );
         }
      }

      /* update the cutoff, propagateagain, and solverelaxagain status of current solving loop */
      updateLoopStatus(set, stat, tree, actdepth, cutoff, &propagateagain, &solverelaxagain);

      /* enforce constraints */
      if( !(*cutoff) && !solverelaxagain && !solvelpagain && !propagateagain )
      {
         /* if the solution changed since the last enforcement, we have to completely reenforce it; otherwise, we
          * only have to enforce the additional constraints added in the last enforcement, but keep the infeasible
          * flag TRUE in order to not declare the infeasible solution feasible due to disregarding the already
          * enforced constraints
          */
         if( lastdomchgcount != stat->domchgcount || lastlpcount != stat->lpcount )
         {
            lastdomchgcount = stat->domchgcount;
            lastlpcount = stat->lpcount;
            *infeasible = FALSE;
         }
         
         /* call constraint enforcement */
         CHECK_OKAY( enforceConstraints(blkmem, set, stat, prob, tree, lp, sepastore, conshdlrs_enfo,
               &branched, cutoff, infeasible, &propagateagain, &solvelpagain) );
         assert(branched == (tree->nchildren > 0));
         assert(!branched || (!(*cutoff) && *infeasible && !propagateagain && !solvelpagain));
         assert(!(*cutoff) || (!branched && *infeasible && !propagateagain && !solvelpagain));
         assert(*infeasible || (!branched && !(*cutoff) && !propagateagain && !solvelpagain));
         assert(!propagateagain || (!branched && !(*cutoff) && *infeasible));
         assert(!solvelpagain || (!branched && !(*cutoff) && *infeasible));

         /* apply found cuts */
         CHECK_OKAY( applyCuts(blkmem, set, stat, tree, lp, sepastore, branchcand, eventqueue, 
               cutoff, &propagateagain, &solvelpagain) );

         /* update the cutoff, propagateagain, and solverelaxagain status of current solving loop */
         updateLoopStatus(set, stat, tree, actdepth, cutoff, &propagateagain, &solverelaxagain);
      }
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      /* if the node is infeasible, but no constraint handler could resolve the infeasibility
       * -> branch on LP or the pseudo solution
       * -> e.g. select non-fixed binary or integer variable x with value x', create three
       *    sons: x <= x'-1, x = x', and x >= x'+1.
       *    In the left and right branch, the current solution is cut off. In the middle
       *    branch, the constraints can hopefully reduce domains of other variables to cut
       *    off the current solution.
       * In LP branching, we cannot allow adding constraints, because this does not necessary change the LP and can
       * therefore lead to an infinite loop.
       */
      if( *infeasible && !(*cutoff) && !(*unbounded) && !solverelaxagain && !solvelpagain && !propagateagain && !branched )
      {
         RESULT result;
         int nlpcands;

         if( SCIPtreeHasFocusNodeLP(tree) )
         {
            CHECK_OKAY( SCIPbranchcandGetLPCands(branchcand, set, stat, lp, NULL, NULL, NULL, &nlpcands, NULL) );
         }
         else
            nlpcands = 0;
         
         if( nlpcands > 0 )
         {
            /* branch on LP solution */
            debugMessage("infeasibility in depth %d was not resolved: branch on LP solution with %d fractionals\n",
               SCIPnodeGetDepth(focusnode), nlpcands);
            CHECK_OKAY( SCIPbranchExecLP(blkmem, set, stat, tree, lp, sepastore, branchcand, eventqueue, 
                  primal->cutoffbound, FALSE, &result) );
            assert(result != SCIP_DIDNOTRUN);
         }
         else
         {
            /* branch on pseudo solution */
            debugMessage("infeasibility in depth %d was not resolved: branch on pseudo solution with %d unfixed integers\n",
               SCIPnodeGetDepth(focusnode), SCIPbranchcandGetNPseudoCands(branchcand));
            CHECK_OKAY( SCIPbranchExecPseudo(blkmem, set, stat, tree, lp, branchcand, eventqueue, 
                  primal->cutoffbound, TRUE, &result) );
         }

         switch( result )
         {  
         case SCIP_CUTOFF:
            assert(tree->nchildren == 0);
            *cutoff = TRUE;
            break;
         case SCIP_CONSADDED:
            assert(tree->nchildren == 0);
            if( nlpcands > 0 )
            {
               errorMessage("LP branching rule added constraint, which was not allowed this time\n");
               return SCIP_INVALIDRESULT;
            }
            solvelpagain = TRUE;
            propagateagain = TRUE;
            break;
         case SCIP_REDUCEDDOM:
            assert(tree->nchildren == 0);
            solvelpagain = TRUE;
            propagateagain = TRUE;
            break;
         case SCIP_SEPARATED:
            assert(tree->nchildren == 0);
            assert(SCIPsepastoreGetNCuts(sepastore) > 0);
            solvelpagain = TRUE;
            break;
         case SCIP_BRANCHED:
            assert(tree->nchildren >= 1);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            branched = TRUE;
            break;
         case SCIP_DIDNOTRUN:
            /* all integer variables in the infeasible solution are fixed,
             * - if no continuous variables exist and all variables are known, the infeasible pseudo solution is completely
             *   fixed, and the node can be cut off
             * - if at least one continuous variable exist or we do not know all variables due to external pricers, we
             *   cannot resolve the infeasibility by branching -> solve LP (and maybe price in additional variables)
             */
            assert(tree->nchildren == 0);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            assert(SCIPbranchcandGetNPseudoCands(branchcand) == 0);

            if( prob->ncontvars == 0 && set->nactivepricers == 0 )
               *cutoff = TRUE;
            else
            {
               assert(!SCIPtreeHasFocusNodeLP(tree)); /* feasible LP solutions with all integers fixed must be feasible */
            
               /* solve the LP in the next loop */
               SCIPtreeSetFocusNodeLP(tree, TRUE);
               solvelpagain = TRUE;
            }
            break;
         default:
            errorMessage("invalid result code <%d> from SCIPbranchLP() or SCIPbranchPseudo()\n", result);
            return SCIP_INVALIDRESULT;
         }  /*lint !e788*/
         assert(*cutoff || solvelpagain || propagateagain || branched); /* something must have been done */
         assert(!(*cutoff) || (!solvelpagain && !propagateagain && !branched));
         assert(!solvelpagain || (!(*cutoff) && !branched));
         assert(!propagateagain || (!(*cutoff) && !branched));
         assert(!branched || (!solvelpagain && !propagateagain));
         assert(branched == (tree->nchildren > 0));

         /* apply found cuts */
         CHECK_OKAY( applyCuts(blkmem, set, stat, tree, lp, sepastore, branchcand, eventqueue, 
               cutoff, &propagateagain, &solvelpagain) );

         /* update the cutoff, propagateagain, and solverelaxagain status of current solving loop */
         updateLoopStatus(set, stat, tree, actdepth, cutoff, &propagateagain, &solverelaxagain);
      }

      /* check for immediate restart */
      *restart = *restart
         || (actdepth == 0 && set->presol_restartbdchgs > 0 && stat->nrootboundchgsrun >= set->presol_restartbdchgs);

      debugMessage("node solving iteration finished: cutoff=%d, propagateagain=%d, solverelaxagain=%d, solvelpagain=%d, nlperrors=%d, restart=%d\n",
         *cutoff, propagateagain, solverelaxagain, solvelpagain, nlperrors, *restart);
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);

   /* check for too many LP errors */
   if( nlperrors >= MAXNLPERRORS )
   {
      errorMessage("(node %lld) unresolved numerical troubles in LP %d -- aborting\n", stat->nnodes, stat->nlps);
      return SCIP_LPERROR;
   }

   /* check for final restart */
   *restart = *restart || (actdepth == 0 && set->presol_restartbdchgs >= 0 && stat->nrootboundchgsrun > 0);

   /* remember root LP solution */
   if( actdepth == 0 && !(*cutoff) && !(*restart) && !(*unbounded) )
   {
      SCIPprobStoreRootSol(prob, SCIPtreeHasFocusNodeLP(tree));
   }

   /* check for cutoff */
   if( *cutoff )
   {
      debugMessage("node is cut off\n");
      SCIPnodeUpdateLowerbound(focusnode, stat, primal->cutoffbound);
      *infeasible = TRUE;
      *restart = FALSE;
      *unbounded = FALSE;
   }

   return SCIP_OKAY;
}

/** calls primal heuristics */
static
RETCODE primalHeuristics(
   SET*             set,                /**< global SCIP settings */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   NODE*            nextnode,           /**< next node that will be processed, or NULL if no more nodes left */
   Bool*            foundsol            /**< pointer to store whether a solution has been found */
   )
{
   RESULT result;
   Bool plunging;
   int ndelayedheurs;
   int depth;
   int lpforkdepth;
   int h;

   assert(set != NULL);
   assert(tree != NULL);
   assert((nextnode == NULL) == (SCIPtreeGetNNodes(tree) == 0));
   assert(foundsol != NULL);

   *foundsol = FALSE;

   /* nothing to do, if no heuristics are available, or if the branch-and-bound process is finished */
   if( set->nheurs == 0 || nextnode == NULL )
      return SCIP_OKAY;

   /* sort heuristics by priority, but move the delayed heuristics to the front */
   SCIPsetSortHeurs(set);

   /* we are in plunging mode iff the next node is a sibling or a child, and no leaf */
   assert(SCIPnodeGetType(nextnode) == SCIP_NODETYPE_SIBLING
      || SCIPnodeGetType(nextnode) == SCIP_NODETYPE_CHILD
      || SCIPnodeGetType(nextnode) == SCIP_NODETYPE_LEAF);
   plunging = (SCIPnodeGetType(nextnode) != SCIP_NODETYPE_LEAF);
   depth = SCIPtreeGetFocusDepth(tree);
   lpforkdepth = tree->focuslpfork != NULL ? SCIPnodeGetDepth(tree->focuslpfork) : -1;

   /* call heuristics */
   ndelayedheurs = 0;
   for( h = 0; h < set->nheurs; ++h )
   {
      CHECK_OKAY( SCIPheurExec(set->heurs[h], set, primal, depth, lpforkdepth, SCIPtreeHasFocusNodeLP(tree), plunging,
            &ndelayedheurs, &result) );
      *foundsol = *foundsol || (result == SCIP_FOUNDSOL);
   }
   assert(0 <= ndelayedheurs && ndelayedheurs <= set->nheurs);

   return SCIP_OKAY;
}

/** if feasible, adds current solution to the solution storage */
static
RETCODE addCurrentSolution(
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   EVENTFILTER*     eventfilter         /**< event filter for global (not variable dependent) events */
   )
{
   SOL* sol;
   Bool foundsol;

   /* found a feasible solution */
   if( SCIPtreeHasFocusNodeLP(tree) )
   {
      /* start clock for LP solutions */
      SCIPclockStart(stat->lpsoltime, set);
      
      /* add solution to storage */
      CHECK_OKAY( SCIPsolCreateLPSol(&sol, blkmem, set, stat, primal, tree, lp, NULL) );
      if( set->misc_exactsolve )
      {
         /* if we want to solve exactly, we have to check the solution exactly again */
         CHECK_OKAY( SCIPprimalTrySolFree(primal, blkmem, set, stat, prob, tree, lp, eventfilter, &sol,
               TRUE, TRUE, TRUE, &foundsol) );
      }
      else
      {
         CHECK_OKAY( SCIPprimalAddSolFree(primal, blkmem, set, stat, prob, tree, lp, eventfilter, &sol, &foundsol) );
      }
      if( foundsol )
         stat->nlpsolsfound++;

      /* stop clock for LP solutions */
      SCIPclockStop(stat->lpsoltime, set);
   }
   else
   {
      /* start clock for pseudo solutions */
      SCIPclockStart(stat->pseudosoltime, set);

      /* add solution to storage */
      CHECK_OKAY( SCIPsolCreatePseudoSol(&sol, blkmem, set, stat, primal, tree, lp, NULL) );
      if( set->misc_exactsolve )
      {
         /* if we want to solve exactly, we have to check the solution exactly again */
         CHECK_OKAY( SCIPprimalTrySolFree(primal, blkmem, set, stat, prob, tree, lp, eventfilter, &sol,
               TRUE, TRUE, TRUE, &foundsol) );
      }
      else
      {
         CHECK_OKAY( SCIPprimalAddSolFree(primal, blkmem, set, stat, prob, tree, lp, eventfilter, &sol, &foundsol) );
      }

      /* stop clock for pseudo solutions */
      SCIPclockStop(stat->pseudosoltime, set);

      if( foundsol )
         stat->npssolsfound++;
   }

   return SCIP_OKAY;
}

/** main solving loop */
RETCODE SCIPsolveCIP(
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   MEM*             mem,                /**< block memory pools */
   PROB*            prob,               /**< transformed problem after presolve */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICESTORE*      pricestore,         /**< pricing storage */
   SEPASTORE*       sepastore,          /**< separation storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   CONFLICT*        conflict,           /**< conflict analysis data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            restart             /**< should solving process be started again with presolving? */
   )
{
   CONSHDLR** conshdlrs_sepa;
   CONSHDLR** conshdlrs_enfo;
   NODESEL* nodesel;
   NODE* focusnode;
   NODE* nextnode;
   EVENT event;
   int nnodes;
   int depth;
   Bool cutoff;
   Bool unbounded;
   Bool infeasible;
   Bool foundsol;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(pricestore != NULL);
   assert(sepastore != NULL);
   assert(branchcand != NULL);
   assert(cutpool != NULL);
   assert(primal != NULL);
   assert(eventfilter != NULL);
   assert(eventqueue != NULL);
   assert(restart != NULL);

   /* check for immediate restart (if problem solving marked to be restarted was aborted) */
   *restart = (SCIPtreeGetCurrentDepth(tree) == 0 && set->presol_restartbdchgs >= 0 && stat->nrootboundchgsrun > 0);

   /* switch status to UNKNOWN */
   stat->status = SCIP_STATUS_UNKNOWN;

   /* sort constraint handlers by priorities */
   ALLOC_OKAY( duplicateMemoryArray(&conshdlrs_sepa, set->conshdlrs, set->nconshdlrs) );
   ALLOC_OKAY( duplicateMemoryArray(&conshdlrs_enfo, set->conshdlrs, set->nconshdlrs) );
   SCIPbsortPtr((void**)conshdlrs_sepa, set->nconshdlrs, SCIPconshdlrCompSepa);
   SCIPbsortPtr((void**)conshdlrs_enfo, set->nconshdlrs, SCIPconshdlrCompEnfo);

   nextnode = NULL;

   while( !SCIPsolveIsStopped(set, stat) && !(*restart) )
   {
      assert(SCIPbufferGetNUsed(set->buffer) == 0);

      foundsol = FALSE;
      infeasible = FALSE;

      do
      {
         /* update the memory saving flag, switch algorithms respectively */
         SCIPstatUpdateMemsaveMode(stat, set, mem);

         /* get the current node selector */
         nodesel = SCIPsetGetNodesel(set, stat);
         
         /* inform tree about the current node selector */
         CHECK_OKAY( SCIPtreeSetNodesel(tree, set, stat, nodesel) );
      
         /* the next node was usually already selected in the previous solving loop before the primal heuristics were
          * called, because they need to know, if the next node will be a child/sibling (plunging) or not;
          * if the heuristics found a new best solution that cut off some of the nodes, the node selector must be called
          * again, because the selected next node may be invalid due to cut off
          */
         if( nextnode == NULL )
         {
            /* select next node to process */
            CHECK_OKAY( SCIPnodeselSelect(nodesel, set, &nextnode) );
         }
         focusnode = nextnode;
         nextnode = NULL;
         assert(SCIPbufferGetNUsed(set->buffer) == 0);
         
         /* start node activation timer */
         SCIPclockStart(stat->nodeactivationtime, set);
         
         /* focus selected node */
         CHECK_OKAY( SCIPnodeFocus(&focusnode, blkmem, set, stat, prob, primal, tree, lp, branchcand,
               eventfilter, eventqueue, &cutoff) );
         if( cutoff )
            stat->ndelayedcutoffs++;

         /* stop node activation timer */
         SCIPclockStop(stat->nodeactivationtime, set);

         assert(SCIPbufferGetNUsed(set->buffer) == 0);
      }
      while( cutoff ); /* select new node, if the current one was located in a cut off subtree */

      assert(SCIPtreeGetCurrentNode(tree) == focusnode);
      assert(SCIPtreeGetFocusNode(tree) == focusnode);

      /* if no more node was selected, we finished optimization */
      if( focusnode == NULL )
      {
         assert(SCIPtreeGetNNodes(tree) == 0);
         break;
      }

      /* update maxdepth and node count statistics */
      depth = SCIPnodeGetDepth(focusnode);
      stat->maxdepth = MAX(stat->maxdepth, depth);
      stat->maxtotaldepth = MAX(stat->maxtotaldepth, depth);
      stat->nnodes++;
      stat->ntotalnodes++;
      
      /* issue NODEFOCUSED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEFOCUSED) );
      CHECK_OKAY( SCIPeventChgNode(&event, focusnode) );
      CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );

      /* solve focus node */
      CHECK_OKAY( solveNode(blkmem, set, stat, prob, primal, tree, lp, pricestore, sepastore, branchcand, cutpool,
            conflict, eventfilter, eventqueue, conshdlrs_sepa, conshdlrs_enfo,
            &cutoff, &unbounded, &infeasible, restart) );
      assert(!cutoff || infeasible);
      assert(SCIPbufferGetNUsed(set->buffer) == 0);
      assert(SCIPtreeGetCurrentNode(tree) == focusnode);
      assert(SCIPtreeGetFocusNode(tree) == focusnode);

      /* check for restart */
      if( !(*restart) )
      {
         /* change color of node in VBC output */
         SCIPvbcSolvedNode(stat->vbc, stat, focusnode);

         /* check, if the current solution is feasible */
         if( !infeasible )
         {
            assert(!SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));
            assert(!cutoff);

            /* node solution is feasible: add it to the solution store */
            CHECK_OKAY( addCurrentSolution(blkmem, set, stat, prob, primal, tree, lp, eventfilter) );

            /* issue NODEFEASIBLE event */
            CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEFEASIBLE) );
            CHECK_OKAY( SCIPeventChgNode(&event, focusnode) );
            CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );
         }
         else if( !unbounded )
         {
            /* node solution is not feasible */
            if( tree->nchildren == 0 )
            {
               /* issue NODEINFEASIBLE event */
               CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEINFEASIBLE) );

               /* increase the cutoff counter of the branching variable */
               if( stat->lastbranchvar != NULL )
               {
                  CHECK_OKAY( SCIPvarIncNCutoffs(stat->lastbranchvar, stat, stat->lastbranchdir) );
               }
               /**@todo if last branching variable is unknown, retrieve it from the nodes' boundchg arrays */
            }
            else
            {
               /* issue NODEBRANCHED event */
               CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEBRANCHED) );
            }
            CHECK_OKAY( SCIPeventChgNode(&event, focusnode) );
            CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );
         }
         assert(SCIPbufferGetNUsed(set->buffer) == 0);

         /* if no branching was created, the node was not cut off, but it's lower bound is still smaller than
          * the cutoff bound, we have to branch on a non-fixed variable;
          * this can happen, if we want to solve exactly, the current solution was declared feasible by the
          * constraint enforcement, but in exact solution checking it was found out to be infeasible;
          * in this case, no branching would have been generated by the enforcement of constraints, but we
          * have to further investigate the current sub tree
          */
         if( !cutoff && !unbounded && tree->nchildren == 0 && SCIPnodeGetLowerbound(focusnode) < primal->cutoffbound )
         {
            RESULT result;

            assert(set->misc_exactsolve);

            do
            {
               result = SCIP_DIDNOTRUN;
               if( SCIPbranchcandGetNPseudoCands(branchcand) == 0 )
               {
                  if( prob->ncontvars > 0 )
                  {
                     /**@todo call PerPlex */
                     errorMessage("cannot branch on all-fixed LP -- have to call PerPlex instead\n");
                  }
               }
               else
               {
                  CHECK_OKAY( SCIPbranchExecPseudo(blkmem, set, stat, tree, lp, branchcand, eventqueue, 
                        primal->cutoffbound, FALSE, &result) );
                  assert(result != SCIP_DIDNOTRUN);
               }
            }
            while( result == SCIP_REDUCEDDOM );
         }
         assert(SCIPbufferGetNUsed(set->buffer) == 0);

         /* select node to process in next solving loop; the primal heuristics need to know whether a child/sibling
          * (plunging) will be selected as next node or not
          */
         CHECK_OKAY( SCIPnodeselSelect(nodesel, set, &nextnode) );
         assert(SCIPbufferGetNUsed(set->buffer) == 0);

         /* call primal heuristics */
         nnodes = SCIPtreeGetNNodes(tree);
         CHECK_OKAY( primalHeuristics(set, primal, tree, nextnode, &foundsol) );
         assert(SCIPbufferGetNUsed(set->buffer) == 0);

         /* if the heuristics found a new best solution that cut off some of the nodes, the node selector must be called
          * again, because the selected next node may be invalid due to cut off
          */
         assert(!tree->cutoffdelayed);
         if( nnodes != SCIPtreeGetNNodes(tree) )
            nextnode = NULL;
      }

      /* display node information line */
      CHECK_OKAY( SCIPdispPrintLine(set, stat, NULL, (SCIPnodeGetDepth(focusnode) == 0) && infeasible && !foundsol) );

      debugMessage("Processing of node in depth %d finished. %d siblings, %d children, %d leaves left\n", 
         SCIPnodeGetDepth(focusnode), tree->nsiblings, tree->nchildren, SCIPtreeGetNLeaves(tree));
      debugMessage("**********************************************************************\n");
   }
   assert(SCIPbufferGetNUsed(set->buffer) == 0);

   debugMessage("Problem solving finished (restart=%d)\n", *restart);

   /* if the current node is the only remaining node, and if its lower bound exceeds the upper bound, we have
    * to delete it manually in order to get to the SOLVED stage instead of thinking, that only the gap limit
    * was reached (this may happen, if the current node is the one defining the global lower bound and a
    * feasible solution with the same value was found at this node)
    */
   if( tree->focusnode != NULL && SCIPtreeGetNNodes(tree) == 0
      && SCIPsetIsGE(set, tree->focusnode->lowerbound, primal->cutoffbound) )
   {
      focusnode = NULL;
      CHECK_OKAY( SCIPnodeFocus(&focusnode, blkmem, set, stat, prob, primal, tree, lp, branchcand, 
            eventfilter, eventqueue, &cutoff) );
   }

   /* free sorted constraint handler arrays */
   freeMemoryArrayNull(&conshdlrs_sepa);
   freeMemoryArrayNull(&conshdlrs_enfo);

   /* check if we finised solving */
   if( SCIPtreeGetNNodes(tree) == 0 && SCIPtreeGetCurrentNode(tree) == NULL )
   {
      /* set the solution status */
      if( unbounded )
      {
         /* switch status to UNBOUNDED */
         stat->status = SCIP_STATUS_UNBOUNDED;
      }
      else if( primal->nsols == 0 )
      {
         /* switch status to INFEASIBLE */
         stat->status = SCIP_STATUS_INFEASIBLE;
      }
      else
      {
         /* switch status to INFEASIBLE */
         stat->status = SCIP_STATUS_OPTIMAL;
      }
   }

   return SCIP_OKAY;
}

