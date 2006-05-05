/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: solve.c,v 1.215 2006/05/05 13:55:25 bzfpfend Exp $"

/**@file   solve.c
 * @brief  main solving loop and node processing
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/buffer.h"
#include "scip/clock.h"
#include "scip/vbc.h"
#include "scip/interrupt.h"
#include "scip/misc.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/sol.h"
#include "scip/primal.h"
#include "scip/tree.h"
#include "scip/pricestore.h"
#include "scip/sepastore.h"
#include "scip/cutpool.h"
#include "scip/solve.h"
#include "scip/scip.h"
#include "scip/branch.h"
#include "scip/conflict.h"
#include "scip/cons.h"
#include "scip/disp.h"
#include "scip/heur.h"
#include "scip/nodesel.h"
#include "scip/pricer.h"
#include "scip/relax.h"
#include "scip/sepa.h"
#include "scip/prop.h"


#define MAXNLPERRORS  10                /**< maximal number of LP error loops in a single node */


/** returns whether the solving process will be / was stopped before proving optimality;
 *  if the solving process was stopped, stores the reason as status in stat
 */
SCIP_Bool SCIPsolveIsStopped(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(set != NULL);
   assert(stat != NULL);

   if( SCIPinterrupted() || stat->userinterrupt )
   {
      stat->status = SCIP_STATUS_USERINTERRUPT;
      stat->userinterrupt = FALSE;
   }
   else if( set->limit_nodes >= 0 && stat->nnodes >= set->limit_nodes )
      stat->status = SCIP_STATUS_NODELIMIT;
   else if( set->limit_stallnodes >= 0 && stat->nnodes >= stat->bestsolnode + set->limit_stallnodes )
      stat->status = SCIP_STATUS_STALLNODELIMIT;
   else if( SCIPclockGetTime(stat->solvingtime) >= set->limit_time )
      stat->status = SCIP_STATUS_TIMELIMIT;
   else if( SCIPgetMemUsed(set->scip) >= set->limit_memory*1024.0*1024.0 )
      stat->status = SCIP_STATUS_MEMLIMIT;
   else if( set->stage >= SCIP_STAGE_SOLVING && SCIPsetIsLT(set, SCIPgetGap(set->scip), set->limit_gap) )
      stat->status = SCIP_STATUS_GAPLIMIT;
   else if( set->limit_solutions >= 0 && set->stage >= SCIP_STAGE_PRESOLVED
      && SCIPgetNSolsFound(set->scip) >= set->limit_solutions )
      stat->status = SCIP_STATUS_SOLLIMIT;
   else if( set->limit_bestsol >= 0 && set->stage >= SCIP_STAGE_PRESOLVED
      && SCIPgetNBestSolsFound(set->scip) >= set->limit_bestsol )
      stat->status = SCIP_STATUS_BESTSOLLIMIT;

   return (stat->status != SCIP_STATUS_UNKNOWN);
}

/** applies one round of propagation */
static
SCIP_RETCODE propagationRound(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth level to use for propagator frequency checks */
   SCIP_Bool             fullpropagation,    /**< should all constraints be propagated (or only new ones)? */
   SCIP_Bool             onlydelayed,        /**< should only delayed propagators be called? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a propagator was delayed */
   SCIP_Bool*            propagain,          /**< pointer to store whether propagation should be applied again */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_RESULT result;
   int i;

   assert(set != NULL);
   assert(delayed != NULL);
   assert(propagain != NULL);
   assert(cutoff != NULL);

   *delayed = FALSE;
   *propagain = FALSE;

   /* sort propagators */
   SCIPsetSortProps(set);

   /* call additional propagators with nonnegative priority */
   for( i = 0; i < set->nprops && !(*cutoff); ++i )
   {
      if( SCIPpropGetPriority(set->props[i]) < 0 )
         continue;

      if( onlydelayed && !SCIPpropWasDelayed(set->props[i]) )
         continue;

      SCIP_CALL( SCIPpropExec(set->props[i], set, stat, depth, onlydelayed, &result) );
      *delayed = *delayed || (result == SCIP_DELAYED);
      *propagain = *propagain || (result == SCIP_REDUCEDDOM);
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      if( *cutoff )
      {
         SCIPdebugMessage(" -> propagator <%s> detected cutoff\n", SCIPpropGetName(set->props[i]));
      }

      /* if we work off the delayed propagators, we stop immediately if a reduction was found */
      if( onlydelayed && result == SCIP_REDUCEDDOM )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* propagate constraints */
   for( i = 0; i < set->nconshdlrs && !(*cutoff); ++i )
   {
      if( onlydelayed && !SCIPconshdlrWasPropagationDelayed(set->conshdlrs[i]) )
         continue;

      SCIP_CALL( SCIPconshdlrPropagate(set->conshdlrs[i], blkmem, set, stat, depth, fullpropagation, onlydelayed,
            &result) );
      *delayed = *delayed || (result == SCIP_DELAYED);
      *propagain = *propagain || (result == SCIP_REDUCEDDOM);
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      if( *cutoff )
      {
         SCIPdebugMessage(" -> constraint handler <%s> detected cutoff in propagation\n", 
            SCIPconshdlrGetName(set->conshdlrs[i]));
      }

      /* if we work off the delayed propagators, we stop immediately if a reduction was found */
      if( onlydelayed && result == SCIP_REDUCEDDOM )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* call additional propagators with negative priority */
   for( i = 0; i < set->nprops && !(*cutoff); ++i )
   {
      if( SCIPpropGetPriority(set->props[i]) >= 0 )
         continue;

      if( onlydelayed && !SCIPpropWasDelayed(set->props[i]) )
         continue;

      SCIP_CALL( SCIPpropExec(set->props[i], set, stat, depth, onlydelayed, &result) );
      *delayed = *delayed || (result == SCIP_DELAYED);
      *propagain = *propagain || (result == SCIP_REDUCEDDOM);
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      if( *cutoff )
      {
         SCIPdebugMessage(" -> propagator <%s> detected cutoff\n", SCIPpropGetName(set->props[i]));
      }

      /* if we work off the delayed propagators, we stop immediately if a reduction was found */
      if( onlydelayed && result == SCIP_REDUCEDDOM )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

/** applies domain propagation on current node */
static
SCIP_RETCODE propagateDomains(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   depth,              /**< depth level to use for propagator frequency checks */
   int                   maxproprounds,      /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   SCIP_Bool             fullpropagation,    /**< should all constraints be propagated (or only new ones)? */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_NODE* node;
   SCIP_Bool delayed;
   SCIP_Bool propagain;
   int propround;

   assert(set != NULL);
   assert(tree != NULL);
   assert(depth >= 0);
   assert(cutoff != NULL);

   node = SCIPtreeGetCurrentNode(tree);
   assert(node != NULL);
   assert(SCIPnodeIsActive(node));
   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_FOCUSNODE
      || SCIPnodeGetType(node) == SCIP_NODETYPE_REFOCUSNODE
      || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);

   /* adjust maximal number of propagation rounds */
   if( maxproprounds == 0 )
      maxproprounds = (depth == 0 ? set->prop_maxroundsroot : set->prop_maxrounds);
   if( maxproprounds == -1 )
      maxproprounds = INT_MAX;

   SCIPdebugMessage("domain propagation of node %p in depth %d (using depth %d, maxrounds %d)\n",
      node, SCIPnodeGetDepth(node), depth, maxproprounds);

   /* propagate as long new bound changes were found and the maximal number of propagation rounds is not exceeded */
   *cutoff = FALSE;
   propround = 0;
   propagain = TRUE;
   while( propagain && !(*cutoff) && propround < maxproprounds )
   {
      propround++;

      /* perform the propagation round by calling the propagators and constraint handlers */
      SCIP_CALL( propagationRound(blkmem, set, stat, depth, fullpropagation, FALSE, &delayed, &propagain, cutoff) );

      /* if the propagation will be terminated, call the delayed propagators */
      while( delayed && (!propagain || propround >= maxproprounds) && !(*cutoff) )
      {
         /* call the delayed propagators and constraint handlers */
         SCIP_CALL( propagationRound(blkmem, set, stat, depth, fullpropagation, TRUE, &delayed, &propagain, cutoff) );
      }

      /* if a reduction was found, we want to do another full propagation round (even if the propagator only claimed
       * to have done a domain reduction without applying a domain change)
       */
      fullpropagation = TRUE;
   }

   /* mark the node to be completely propagated in the current repropagation subtree level */
   SCIPnodeMarkPropagated(node, tree);

   return SCIP_OKAY;
}

/** applies domain propagation on current node and flushes the conflict storage afterwards */
SCIP_RETCODE SCIPpropagateDomains(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   int                   depth,              /**< depth level to use for propagator frequency checks */
   int                   maxproprounds,      /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   /* apply domain propagation */
   SCIP_CALL( propagateDomains(blkmem, set, stat, tree, depth, maxproprounds, TRUE, cutoff) );

   /* flush the conflict set storage */
   SCIP_CALL( SCIPconflictFlushConss(conflict, blkmem, set, stat, prob, tree) );

   return SCIP_OKAY;
}

/** returns whether the given variable with the old LP solution value should lead to an update of the pseudo cost entry */
static
SCIP_Bool isPseudocostUpdateValid(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             oldlpsolval         /**< solution value of variable in old LP */
   )
{
   SCIP_Real newlpsolval;

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
typedef enum PseudocostFlag PSEUDOCOSTFLAG;

/** updates the variable's pseudo cost values after the node's initial LP was solved */
static
SCIP_RETCODE updatePseudocost(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< LP data */
   )
{
   SCIP_NODE* focusnode;
   int actdepth;

   assert(lp != NULL);
   assert(tree != NULL);
   assert(tree->path != NULL);

   focusnode = SCIPtreeGetFocusNode(tree);
   assert(SCIPnodeIsActive(focusnode));
   assert(SCIPnodeGetType(focusnode) == SCIP_NODETYPE_FOCUSNODE);
   actdepth = SCIPnodeGetDepth(focusnode);
   assert(tree->path[actdepth] == focusnode);

   if( lp->solved && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && tree->focuslpstatefork != NULL )
   {
      SCIP_BOUNDCHG** updates;
      SCIP_NODE* node;
      SCIP_VAR* var;
      SCIP_Real weight;
      SCIP_Real lpgain;
      int nupdates;
      int nvalidupdates;
      int d;
      int i;

      assert(SCIPnodeIsActive(tree->focuslpstatefork));
      assert(tree->path[tree->focuslpstatefork->depth] == tree->focuslpstatefork);

      /* get a buffer for the collected bound changes; start with a size twice as large as the number of nodes between
       * current node and LP fork
       */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &updates, 2*(actdepth - tree->focuslpstatefork->depth)) );
      nupdates = 0;
      nvalidupdates = 0;

      /* search the nodes from LP fork down to current node for bound changes in between; move in this direction,
       * because the bound changes closer to the LP fork are more likely to have a valid LP solution information
       * attached; collect the bound changes for pseudo cost value updates and mark the corresponding variables such
       * that they are not updated twice in case of more than one bound change on the same variable
       */
      for( d = tree->focuslpstatefork->depth+1; d <= actdepth; ++d )
      {
         node = tree->path[d];

         if( node->domchg != NULL )
         {
            SCIP_BOUNDCHG* boundchgs;
            int nboundchgs;

            boundchgs = node->domchg->domchgbound.boundchgs;
            nboundchgs = node->domchg->domchgbound.nboundchgs;
            for( i = 0; i < nboundchgs; ++i )
            {
               /* we even collect redundant bound changes, since they were not redundant in the LP branching decision
                * and therefore should be regarded in the pseudocost updates
                */ 
               if( (SCIP_BOUNDCHGTYPE)boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
               {
                  var = boundchgs[i].var;
                  assert(var != NULL);
                  if( (PSEUDOCOSTFLAG)var->pseudocostflag == PSEUDOCOST_NONE )
                  {
                     /* remember the bound change and mark the variable */
                     SCIP_CALL( SCIPsetReallocBufferArray(set, &updates, nupdates+1) );
                     updates[nupdates] = &boundchgs[i];
                     nupdates++;

                     /* check, if the bound change would lead to a valid pseudo cost update */
                     if( isPseudocostUpdateValid(var, set, boundchgs[i].data.branchingdata.lpsolval) )
                     {
                        var->pseudocostflag = PSEUDOCOST_UPDATE; /*lint !e641*/
                        nvalidupdates++;
                     }
                     else
                        var->pseudocostflag = PSEUDOCOST_IGNORE; /*lint !e641*/
                  }
               }
            }
         }
      }

      /* update the pseudo cost values and reset the variables' flags; assume, that the responsibility for the dual gain
       * is equally spread on all bound changes that lead to valid pseudo cost updates
       */
      weight = nvalidupdates > 0 ? 1.0 / (SCIP_Real)nvalidupdates : 1.0;
      lpgain = (SCIPlpGetObjval(lp, set) - tree->focuslpstatefork->lowerbound) * weight;
      lpgain = MAX(lpgain, 0.0);
      for( i = 0; i < nupdates; ++i )
      {
         assert((SCIP_BOUNDCHGTYPE)updates[i]->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING);
         var = updates[i]->var;
         assert(var != NULL);
         assert((PSEUDOCOSTFLAG)var->pseudocostflag != PSEUDOCOST_NONE);
         if( (PSEUDOCOSTFLAG)var->pseudocostflag == PSEUDOCOST_UPDATE )
         {
            SCIPdebugMessage("updating pseudocosts of <%s>: sol: %g -> %g, LP: %e -> %e => gain=%g, weight: %g\n",
               SCIPvarGetName(var), updates[i]->data.branchingdata.lpsolval, SCIPvarGetLPSol(var),
               tree->focuslpstatefork->lowerbound, SCIPlpGetObjval(lp, set), lpgain, weight);
            SCIP_CALL( SCIPvarUpdatePseudocost(var, set, stat,
                  SCIPvarGetLPSol(var) - updates[i]->data.branchingdata.lpsolval, lpgain, weight) );
         }
         var->pseudocostflag = PSEUDOCOST_NONE; /*lint !e641*/
      }

      /* free the buffer for the collected bound changes */
      SCIPsetFreeBufferArray(set, &updates);
   }

   return SCIP_OKAY;
}

/** updates lower bound of node using lower bound of LP */
static
SCIP_RETCODE nodeUpdateLowerboundLP(
   SCIP_NODE*            node,               /**< node to set lower bound for */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< LP data */
   )
{
   SCIP_Real lpobjval;

   assert(set != NULL);

   if( set->misc_exactsolve )
   {
      SCIP_CALL( SCIPlpGetProvedLowerbound(lp, set, &lpobjval) );
   }
   else
      lpobjval = SCIPlpGetObjval(lp, set);

   SCIPnodeUpdateLowerbound(node, stat, lpobjval);

   return SCIP_OKAY;
}

/** constructs the initial LP of the current node */
static
SCIP_RETCODE initLP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool             root,               /**< is this the initial root LP? */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_VAR* var;
   int v;
   int h;

   assert(set != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* at the root node, we have to add the initial variables as columns */
   if( root )
   {
      assert(SCIPlpGetNCols(lp) == 0);
      assert(SCIPlpGetNRows(lp) == 0);
      assert(lp->nremoveablecols == 0);
      assert(lp->nremoveablerows == 0);

      /* inform pricing storage, that LP is now filled with initial data */
      SCIPpricestoreStartInitialLP(pricestore);

      /* add all initial variables to LP */
      SCIPdebugMessage("init LP: initial columns\n");
      for( v = 0; v < prob->nvars; ++v )
      {
         var = prob->vars[v];
         assert(SCIPvarGetProbindex(var) >= 0);

         if( SCIPvarIsInitial(var) )
         {
            SCIP_CALL( SCIPpricestoreAddVar(pricestore, blkmem, set, lp, var, 0.0, TRUE) );
         }
      }
      assert(lp->nremoveablecols == 0);
      SCIP_CALL( SCIPpricestoreApplyVars(pricestore, blkmem, set, stat, prob, tree, lp) );

      /* inform pricing storage, that initial LP setup is now finished */
      SCIPpricestoreEndInitialLP(pricestore);
   }

   /* inform separation storage, that LP is now filled with initial data */
   SCIPsepastoreStartInitialLP(sepastore);

   /* add LP relaxations of all initial constraints to LP */
   SCIPdebugMessage("init LP: initial rows\n");
   for( h = 0; h < set->nconshdlrs; ++h )
   {
      SCIP_CALL( SCIPconshdlrInitLP(set->conshdlrs[h], blkmem, set, stat) );
   }
   SCIP_CALL( SCIPsepastoreApplyCuts(sepastore, blkmem, set, stat, tree, lp, branchcand, eventqueue, cutoff) );

   /* inform separation storage, that initial LP setup is now finished */
   SCIPsepastoreEndInitialLP(sepastore);

   return SCIP_OKAY;
}

/** constructs the LP of the current node, but does not load the LP state and warmstart information  */
SCIP_RETCODE SCIPconstructCurrentLP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_Bool initroot;

   assert(tree != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   if( !SCIPtreeIsFocusNodeLPConstructed(tree) )
   {
      /* load the LP into the solver and load the LP state */
      SCIPdebugMessage("loading LP\n");
      SCIP_CALL( SCIPtreeLoadLP(tree, blkmem, set, lp, &initroot) );
      assert(initroot || SCIPnodeGetDepth(SCIPtreeGetFocusNode(tree)) > 0);
      assert(SCIPtreeIsFocusNodeLPConstructed(tree));
      
      /* setup initial LP relaxation of node */
      SCIP_CALL( initLP(blkmem, set, stat, prob, tree, lp, pricestore, sepastore, branchcand, eventqueue, initroot,
            cutoff) );
   }

   return SCIP_OKAY;
}

/** load and solve the initial LP of a node */
static
SCIP_RETCODE solveNodeInitialLP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved error in LP solving occured */
   )
{
   assert(stat != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(cutoff != NULL);
   assert(lperror != NULL);
   assert(SCIPtreeGetFocusNode(tree) != NULL);
   assert(SCIPnodeGetType(SCIPtreeGetFocusNode(tree)) == SCIP_NODETYPE_FOCUSNODE);

   *cutoff = FALSE;
   *lperror = FALSE;

   /* load the LP into the solver */
   SCIP_CALL( SCIPconstructCurrentLP(blkmem, set, stat, prob, tree, lp, pricestore, sepastore, branchcand, eventqueue,
         cutoff) );
   if( *cutoff )
      return SCIP_OKAY;

   /* load the LP state */
   SCIP_CALL( SCIPtreeLoadLPState(tree, blkmem, set, stat, lp) );

   /* solve initial LP */
   SCIPdebugMessage("node: solve initial LP\n");
   SCIP_CALL( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, TRUE, FALSE, lperror) );
   assert(lp->flushed);
   assert(lp->solved || *lperror);

   if( !(*lperror) )
   {
      SCIP_EVENT event;

      /* issue FIRSTLPSOLVED event */
      SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_FIRSTLPSOLVED) );
      SCIP_CALL( SCIPeventChgNode(&event, SCIPtreeGetFocusNode(tree)) );
      SCIP_CALL( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );

      /* update pseudo cost values */
      SCIP_CALL( updatePseudocost(set, stat, tree, lp) );
   }

   return SCIP_OKAY;
}

/** calls primal heuristics */
static
SCIP_RETCODE primalHeuristics(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            nextnode,           /**< next node that will be processed, or NULL if no more nodes left */
   SCIP_Bool             nodesolved,         /**< is the current node already solved? */
   SCIP_Bool             inlploop,           /**< are we currently in the LP solving loop? */
   SCIP_Bool*            foundsol            /**< pointer to store whether a solution has been found */
   )
{
   SCIP_RESULT result;
   SCIP_Bool plunging;
   SCIP_Longint oldnbestsolsfound;
   int ndelayedheurs;
   int depth;
   int lpstateforkdepth;
   int h;

   assert(set != NULL);
   assert(primal != NULL);
   assert(tree != NULL);
   assert(!nodesolved || (nextnode == NULL) == (SCIPtreeGetNNodes(tree) == 0));
   assert(!nodesolved || !inlploop);
   assert(foundsol != NULL);

   *foundsol = FALSE;

   /* nothing to do, if no heuristics are available, or if the branch-and-bound process is finished */
   if( set->nheurs == 0 || (nodesolved && nextnode == NULL) )
      return SCIP_OKAY;

   /* sort heuristics by priority, but move the delayed heuristics to the front */
   SCIPsetSortHeurs(set);

   /* we are in plunging mode iff the next node is a sibling or a child, and no leaf */
   assert(!nodesolved
      || SCIPnodeGetType(nextnode) == SCIP_NODETYPE_SIBLING
      || SCIPnodeGetType(nextnode) == SCIP_NODETYPE_CHILD
      || SCIPnodeGetType(nextnode) == SCIP_NODETYPE_LEAF);
   plunging = (!nodesolved || SCIPnodeGetType(nextnode) != SCIP_NODETYPE_LEAF);
   depth = SCIPtreeGetFocusDepth(tree);
   lpstateforkdepth = tree->focuslpstatefork != NULL ? SCIPnodeGetDepth(tree->focuslpstatefork) : -1;

   /* call heuristics */
   ndelayedheurs = 0;
   oldnbestsolsfound = primal->nbestsolsfound;
   for( h = 0; h < set->nheurs; ++h )
   {
      SCIP_CALL( SCIPheurExec(set->heurs[h], set, primal, depth, lpstateforkdepth, SCIPtreeHasFocusNodeLP(tree),
            plunging, nodesolved, inlploop, &ndelayedheurs, &result) );
   }
   assert(0 <= ndelayedheurs && ndelayedheurs <= set->nheurs);

   *foundsol = (primal->nbestsolsfound > oldnbestsolsfound);

   return SCIP_OKAY;
}

/** applies one round of LP separation */
static
SCIP_RETCODE separationRoundLP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   int                   actdepth,           /**< current depth in the tree */
   SCIP_Bool             onlydelayed,        /**< should only delayed separators be called? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a separator was delayed */
   SCIP_Bool*            enoughcuts,         /**< pointer to store whether enough cuts have been found this round */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_RESULT result;
   int i;
   SCIP_Bool consadded;
   SCIP_Bool root;

   assert(set != NULL);
   assert(lp != NULL);
   assert(set->conshdlrs_sepa != NULL);
   assert(delayed != NULL);
   assert(enoughcuts != NULL);
   assert(cutoff != NULL);

   *delayed = FALSE;
   *enoughcuts = FALSE;
   consadded = FALSE;
   root = (actdepth == 0);
   
   /* sort separators by priority */
   SCIPsetSortSepas(set);

   /* call LP separators with nonnegative priority */
   for( i = 0; i < set->nsepas && !(*cutoff) && !(*enoughcuts) && lp->flushed; ++i )
   {
      if( SCIPsepaGetPriority(set->sepas[i]) < 0 )
         continue;

      if( onlydelayed && !SCIPsepaWasLPDelayed(set->sepas[i]) )
         continue;

      SCIPdebugMessage(" -> executing separator <%s> with priority %d\n", 
         SCIPsepaGetName(set->sepas[i]), SCIPsepaGetPriority(set->sepas[i]));
      SCIP_CALL( SCIPsepaExecLP(set->sepas[i], set, stat, sepastore, actdepth, onlydelayed, &result) );
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      consadded = consadded || (result == SCIP_CONSADDED);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
      *delayed = *delayed || (result == SCIP_DELAYED);
      if( *cutoff )
      {
         SCIPdebugMessage(" -> separator <%s> detected cutoff\n", SCIPsepaGetName(set->sepas[i]));
      }

      /* if we work off the delayed separators, we stop immediately if a cut was found */
      if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED) )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* try separating constraints of the constraint handlers */
   for( i = 0; i < set->nconshdlrs && !(*cutoff) && !(*enoughcuts) && lp->flushed; ++i )
   {
      if( onlydelayed && !SCIPconshdlrWasLPSeparationDelayed(set->conshdlrs_sepa[i]) )
         continue;

      SCIPdebugMessage(" -> executing separation of constraint handler <%s> with priority %d\n", 
		       SCIPconshdlrGetName(set->conshdlrs_sepa[i]), SCIPconshdlrGetSepaPriority(set->conshdlrs_sepa[i]));
      SCIP_CALL( SCIPconshdlrSeparateLP(set->conshdlrs_sepa[i], blkmem, set, stat, sepastore, actdepth, onlydelayed,
            &result) );
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      consadded = consadded || (result == SCIP_CONSADDED);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
      *delayed = *delayed || (result == SCIP_DELAYED);
      if( *cutoff )
      {
         SCIPdebugMessage(" -> constraint handler <%s> detected cutoff in separation\n",
            SCIPconshdlrGetName(set->conshdlrs_sepa[i]));
      }

      /* if we work off the delayed separators, we stop immediately if a cut was found */
      if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED) )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* call LP separators with negative priority */
   for( i = 0; i < set->nsepas && !(*cutoff) && !(*enoughcuts) && lp->flushed; ++i )
   {
      if( SCIPsepaGetPriority(set->sepas[i]) >= 0 )
         continue;

      if( onlydelayed && !SCIPsepaWasLPDelayed(set->sepas[i]) )
         continue;

      SCIPdebugMessage(" -> executing separator <%s> with priority %d\n", 
		       SCIPsepaGetName(set->sepas[i]), SCIPsepaGetPriority(set->sepas[i]));
      SCIP_CALL( SCIPsepaExecLP(set->sepas[i], set, stat, sepastore, actdepth, onlydelayed, &result) );
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      consadded = consadded || (result == SCIP_CONSADDED);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
      *delayed = *delayed || (result == SCIP_DELAYED);
      if( *cutoff )
      {
         SCIPdebugMessage(" -> separator <%s> detected cutoff\n", SCIPsepaGetName(set->sepas[i]));
      }

      /* if we work off the delayed separators, we stop immediately if a cut was found */
      if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED) )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* process the constraints that were added during this separation round */
   while( consadded )
   {
      consadded = FALSE;

      for( i = 0; i < set->nconshdlrs && !(*cutoff) && !(*enoughcuts) && lp->flushed; ++i )
      {
	 if( onlydelayed && !SCIPconshdlrWasLPSeparationDelayed(set->conshdlrs_sepa[i]) )
	    continue;

	 SCIPdebugMessage(" -> executing separation of constraint handler <%s> with priority %d\n", 
			  SCIPconshdlrGetName(set->conshdlrs_sepa[i]), SCIPconshdlrGetSepaPriority(set->conshdlrs_sepa[i]));
	 SCIP_CALL( SCIPconshdlrSeparateLP(set->conshdlrs_sepa[i], blkmem, set, stat, sepastore, actdepth, onlydelayed,
					   &result) );
	 *cutoff = *cutoff || (result == SCIP_CUTOFF);
	 consadded = consadded || (result == SCIP_CONSADDED);
	 *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
	 *delayed = *delayed || (result == SCIP_DELAYED);
         if( *cutoff )
         {
            SCIPdebugMessage(" -> constraint handler <%s> detected cutoff in separation\n",
               SCIPconshdlrGetName(set->conshdlrs_sepa[i]));
         }

	 /* if we work off the delayed separators, we stop immediately if a cut was found */
	 if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED) )
	 {
	    *delayed = TRUE;
	    return SCIP_OKAY;
	 }
      }
   }

   return SCIP_OKAY;
}

/** applies one round of separation on the given primal solution */
static
SCIP_RETCODE separationRoundSol(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   int                   actdepth,           /**< current depth in the tree */
   SCIP_Bool             onlydelayed,        /**< should only delayed separators be called? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a separator was delayed */
   SCIP_Bool*            enoughcuts,         /**< pointer to store whether enough cuts have been found this round */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_RESULT result;
   int i;
   SCIP_Bool consadded;
   SCIP_Bool root;

   assert(set != NULL);
   assert(set->conshdlrs_sepa != NULL);
   assert(delayed != NULL);
   assert(enoughcuts != NULL);
   assert(cutoff != NULL);

   *delayed = FALSE;
   *enoughcuts = FALSE;
   consadded = FALSE;
   root = (actdepth == 0);

   /* sort separators by priority */
   SCIPsetSortSepas(set);

   /* call separators with nonnegative priority */
   for( i = 0; i < set->nsepas && !(*cutoff) && !(*enoughcuts); ++i )
   {
      if( SCIPsepaGetPriority(set->sepas[i]) < 0 )
         continue;

      if( onlydelayed && !SCIPsepaWasSolDelayed(set->sepas[i]) )
         continue;

      SCIP_CALL( SCIPsepaExecSol(set->sepas[i], set, stat, sepastore, sol, actdepth, onlydelayed, &result) );
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      consadded = consadded || (result == SCIP_CONSADDED);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
      *delayed = *delayed || (result == SCIP_DELAYED);
      if( *cutoff )
      {
         SCIPdebugMessage(" -> separator <%s> detected cutoff\n", SCIPsepaGetName(set->sepas[i]));
      }

      /* if we work off the delayed separators, we stop immediately if a cut was found */
      if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED) )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* try separating constraints of the constraint handlers */
   for( i = 0; i < set->nconshdlrs && !(*cutoff) && !(*enoughcuts); ++i )
   {
      if( onlydelayed && !SCIPconshdlrWasSolSeparationDelayed(set->conshdlrs_sepa[i]) )
         continue;

      SCIP_CALL( SCIPconshdlrSeparateSol(set->conshdlrs_sepa[i], blkmem, set, stat, sepastore, sol, actdepth, onlydelayed,
            &result) );
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      consadded = consadded || (result == SCIP_CONSADDED);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
      *delayed = *delayed || (result == SCIP_DELAYED);
      if( *cutoff )
      {
         SCIPdebugMessage(" -> constraint handler <%s> detected cutoff in separation\n",
            SCIPconshdlrGetName(set->conshdlrs_sepa[i]));
      }

      /* if we work off the delayed separators, we stop immediately if a cut was found */
      if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED) )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* call separators with negative priority */
   for( i = 0; i < set->nsepas && !(*cutoff) && !(*enoughcuts); ++i )
   {
      if( SCIPsepaGetPriority(set->sepas[i]) >= 0 )
         continue;

      if( onlydelayed && !SCIPsepaWasSolDelayed(set->sepas[i]) )
         continue;

      SCIP_CALL( SCIPsepaExecSol(set->sepas[i], set, stat, sepastore, sol, actdepth, onlydelayed, &result) );
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      consadded = consadded || (result == SCIP_CONSADDED);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
      *delayed = *delayed || (result == SCIP_DELAYED);
      if( *cutoff )
      {
         SCIPdebugMessage(" -> separator <%s> detected cutoff\n", SCIPsepaGetName(set->sepas[i]));
      }

      /* if we work off the delayed separators, we stop immediately if a cut was found */
      if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED) )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* process the constraints that were added during this separation round */
   while( consadded )
   {
      consadded = FALSE;

      for( i = 0; i < set->nconshdlrs && !(*cutoff) && !(*enoughcuts); ++i )
      {
	 if( onlydelayed && !SCIPconshdlrWasSolSeparationDelayed(set->conshdlrs_sepa[i]) )
	    continue;

	 SCIP_CALL( SCIPconshdlrSeparateSol(set->conshdlrs_sepa[i], blkmem, set, stat, sepastore, sol, actdepth, onlydelayed,
					    &result) );
	 *cutoff = *cutoff || (result == SCIP_CUTOFF);
	 consadded = consadded || (result == SCIP_CONSADDED);
	 *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
	 *delayed = *delayed || (result == SCIP_DELAYED);
         if( *cutoff )
         {
            SCIPdebugMessage(" -> constraint handler <%s> detected cutoff in separation\n",
               SCIPconshdlrGetName(set->conshdlrs_sepa[i]));
         }

	 /* if we work off the delayed separators, we stop immediately if a cut was found */
	 if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED) )
	 {
	    *delayed = TRUE;
	    return SCIP_OKAY;
	 }
      }
   }

   return SCIP_OKAY;
}

/** applies one round of separation on the given primal solution or on the LP solution */
SCIP_RETCODE SCIPseparationRound(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   int                   actdepth,           /**< current depth in the tree */
   SCIP_Bool             onlydelayed,        /**< should only delayed separators be called? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a separator was delayed */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_Bool enoughcuts;

   assert(delayed != NULL);
   assert(cutoff != NULL);

   *delayed = FALSE;
   *cutoff = FALSE;
   enoughcuts = FALSE;

   if( sol == NULL )
   {
      /* apply a separation round on the LP solution */
      SCIP_CALL( separationRoundLP(blkmem, set, stat, lp, sepastore, actdepth, onlydelayed, delayed, &enoughcuts, cutoff) );
   }
   else
   {
      /* apply a separation round on the given primal solution */
      SCIP_CALL( separationRoundSol(blkmem, set, stat, sepastore, sol, actdepth, onlydelayed, delayed, &enoughcuts, cutoff) );
   }

   return SCIP_OKAY;
}

/** solves the current LP completely with pricing in new variables */
SCIP_RETCODE SCIPpriceLoop(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool             pretendroot,        /**< should the pricers be called as if we are at the root node? */
   SCIP_Bool             displayinfo,        /**< should info lines be displayed after each pricing round? */
   int                   maxpricerounds,     /**< maximal number of pricing rounds (-1: no limit);
                                              *   a finite limit means that the LP might not be solved to optimality! */
   int*                  npricedcolvars,     /**< pointer to store number of column variables after problem vars were priced */
   SCIP_Bool*            mustsepa,           /**< pointer to store TRUE if a separation round should follow */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved error in LP solving occured */
   )
{
   int npricerounds;
   SCIP_Bool mustprice;

   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(npricedcolvars != NULL);
   assert(mustsepa != NULL);
   assert(lperror != NULL);

   *npricedcolvars = prob->ncolvars;
   *lperror = FALSE;

   /* if the LP is unbounded, we don't need to price */
   mustprice = (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

   /* if all the variables are already in the LP, we don't need to price */
   mustprice = mustprice && !SCIPprobAllColsInLP(prob, set, lp);

   /* check if infinite number of pricing rounds should be used */
   if( maxpricerounds == -1 )
      maxpricerounds = INT_MAX;

   /* pricing (has to be done completely to get a valid lower bound) */
   npricerounds = 0;
   while( !(*lperror) && mustprice && npricerounds < maxpricerounds )
   {
      SCIP_Bool enoughvars;
      int p;

      assert(lp->flushed);
      assert(lp->solved);
      assert(SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

      /* price problem variables */
      SCIPdebugMessage("problem variable pricing\n");
      assert(SCIPpricestoreGetNVars(pricestore) == 0);
      assert(SCIPpricestoreGetNBoundResets(pricestore) == 0);
      SCIP_CALL( SCIPpricestoreAddProbVars(pricestore, blkmem, set, stat, prob, tree, lp, branchcand, eventqueue) );
      *npricedcolvars = prob->ncolvars;

      /* call external pricers to create additional problem variables */
      SCIPdebugMessage("external variable pricing\n");
      
      /* sort pricer algorithms by priority */
      SCIPsetSortPricers(set);
      
      /* call external pricer algorithms, that are active for the current problem */
      enoughvars = (SCIPpricestoreGetNVars(pricestore) >= SCIPsetGetPriceMaxvars(set, pretendroot)/2);
      for( p = 0; p < set->nactivepricers && !enoughvars; ++p )
      {
         SCIP_CALL( SCIPpricerExec(set->pricers[p], set, prob, lp, pricestore) );
         enoughvars = enoughvars || (SCIPpricestoreGetNVars(pricestore) >= SCIPsetGetPriceMaxvars(set, pretendroot)/2);
      }
   
      /* apply the priced variables to the LP */
      SCIP_CALL( SCIPpricestoreApplyVars(pricestore, blkmem, set, stat, prob, tree, lp) );
      assert(SCIPpricestoreGetNVars(pricestore) == 0);
      assert(!lp->flushed || lp->solved);
      mustprice = !lp->flushed || (prob->ncolvars != *npricedcolvars);
      *mustsepa = *mustsepa || !lp->flushed;

      /* after adding columns, the LP should be primal feasible such that primal simplex is applicable;
       * if LP was infeasible, we have to use dual simplex
       */
      SCIPdebugMessage("pricing: solve LP\n");
      SCIP_CALL( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, TRUE, FALSE, lperror) );
      assert(lp->flushed);
      assert(lp->solved || *lperror);

      /* reset bounds temporarily set by pricer to their original values */
      SCIPdebugMessage("pricing: reset bounds\n");
      SCIP_CALL( SCIPpricestoreResetBounds(pricestore, blkmem, set, stat, lp, branchcand, eventqueue) );
      assert(SCIPpricestoreGetNVars(pricestore) == 0);
      assert(SCIPpricestoreGetNBoundResets(pricestore) == 0);
      assert(!lp->flushed || lp->solved || *lperror);
      mustprice = mustprice || !lp->flushed || (prob->ncolvars != *npricedcolvars);
      *mustsepa = *mustsepa || !lp->flushed;

      /* solve LP again after resetting bounds (with dual simplex) */
      SCIPdebugMessage("pricing: solve LP after resetting bounds\n");
      SCIP_CALL( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, FALSE, FALSE, lperror) );
      assert(lp->flushed);
      assert(lp->solved || *lperror);

      /* increase pricing round counter */
      stat->npricerounds++;
      npricerounds++;

      /* display node information line */
      if( displayinfo && mustprice && (SCIP_VERBLEVEL)set->disp_verblevel >= SCIP_VERBLEVEL_FULL )
      {
         SCIP_CALL( SCIPdispPrintLine(set, stat, NULL, TRUE) );
      }

      /* if the LP is unbounded, we can stop pricing */
      mustprice = mustprice && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY);
   }
   assert(lp->flushed);
   assert(lp->solved || *lperror);

   return SCIP_OKAY;
}

/** solve the current LP of a node with a price-and-cut loop */
static
SCIP_RETCODE priceAndCutLoop(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cut pool */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool             initialloop,        /**< is this the first price-and-cut loop call at the current node? */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            unbounded,          /**< pointer to store whether an unbounded ray was found in the LP */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved error in LP solving occured */
   )
{
   SCIP_NODE* focusnode;
   SCIP_RESULT result;
   SCIP_EVENT event;
   SCIP_Real loclowerbound;
   SCIP_Real glblowerbound;
   SCIP_Real stalllpobjval;
   SCIP_Bool separate;
   SCIP_Bool mustprice;
   SCIP_Bool mustsepa;
   SCIP_Bool delayedsepa;
   SCIP_Bool root;
   int maxseparounds;
   int nsepastallrounds;
   int maxnsepastallrounds;
   int stallnfracs;
   int actdepth;
   int npricedcolvars;

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
   SCIPdebugMessage("node: solve LP with price and cut\n");
   SCIP_CALL( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, TRUE, FALSE, lperror) );
   assert(lp->flushed);
   assert(lp->solved || *lperror);

   /* price-and-cut loop */
   npricedcolvars = prob->ncolvars;
   mustprice = TRUE;
   mustsepa = separate;
   delayedsepa = FALSE;
   *cutoff = FALSE;
   *unbounded = FALSE;
   nsepastallrounds = 0;
   stalllpobjval = SCIP_REAL_MAX;
   stallnfracs = INT_MAX;
   while( !(*cutoff) && !(*lperror) && (mustprice || mustsepa || delayedsepa) )
   {
      SCIP_Bool foundsol;

      SCIPdebugMessage("-------- node solving loop --------\n");
      assert(lp->flushed);
      assert(lp->solved);

      /* solve the LP with pricing in new variables */
      if( mustprice )
      {
         SCIP_CALL( SCIPpriceLoop(blkmem, set, stat, prob, tree, lp, pricestore, branchcand, eventqueue,
               root, root, -1, &npricedcolvars, &mustsepa, lperror) );
         mustprice = FALSE;
      }
      assert(lp->flushed);
      assert(lp->solved || *lperror);

      /* update lower bound w.r.t. the the LP solution */
      if( !(*lperror) )
      {
         SCIP_CALL( nodeUpdateLowerboundLP(focusnode, set, stat, lp) );
         SCIPdebugMessage(" -> new lower bound: %g (LP status: %d, LP obj: %g)\n",
            SCIPnodeGetLowerbound(focusnode), SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set));
      }
      else
      {
         SCIPdebugMessage(" -> error solving LP. keeping old bound: %g\n", SCIPnodeGetLowerbound(focusnode));
      }

      /* call primal heuristics that are applicable during node processing loop */
      if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         SCIP_CALL( primalHeuristics(set, primal, tree, NULL, FALSE, TRUE, &foundsol) );
      }
      else
         foundsol = FALSE;

      /* display node information line for root node */
      if( !foundsol && root && (SCIP_VERBLEVEL)set->disp_verblevel >= SCIP_VERBLEVEL_HIGH )
      {
         SCIP_CALL( SCIPdispPrintLine(set, stat, NULL, TRUE) );
      }

      /* check, if we exceeded the separation round limit */
      mustsepa = mustsepa
         && stat->nseparounds < maxseparounds
         && nsepastallrounds < maxnsepastallrounds;

      /* if separators were delayed, we want to apply a final separation round with the delayed separators */
      delayedsepa = delayedsepa && !mustsepa; /* if regular separation applies, we ignore delayed separators */
      mustsepa = mustsepa || delayedsepa;

      /* if the LP is infeasible or exceeded the objective limit, we don't need to separate cuts */
      if( !separate
         || (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE)
         || (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT)
         || SCIPsetIsGE(set, SCIPnodeGetLowerbound(focusnode), primal->cutoffbound) )
      {
         mustsepa = FALSE;
         delayedsepa = FALSE;
      }

      /* separation and reduced cost strengthening
       * (needs not to be done completely, because we just want to increase the lower bound)
       */
      if( !(*cutoff) && !(*lperror) && mustsepa )
      {
         SCIP_Longint olddomchgcount;
         SCIP_Bool enoughcuts;

         assert(lp->flushed);
         assert(lp->solved);
         assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);

         olddomchgcount = stat->domchgcount;

         mustsepa = FALSE;
         enoughcuts = (SCIPsetGetSepaMaxcuts(set, root) == 0);

         /* global cut pool separation */
         if( !enoughcuts && !delayedsepa )
         {
            if( (set->sepa_poolfreq == 0 && actdepth == 0)
               || (set->sepa_poolfreq > 0 && actdepth % set->sepa_poolfreq == 0) )
            {
               SCIPdebugMessage("global cut pool separation\n");
               assert(SCIPsepastoreGetNCuts(sepastore) == 0);
               SCIP_CALL( SCIPcutpoolSeparate(cutpool, blkmem, set, stat, lp, sepastore, root, &result) );
               *cutoff = *cutoff || (result == SCIP_CUTOFF);
               enoughcuts = enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetSepaMaxcuts(set, root));
               if( *cutoff )
               {
                  SCIPdebugMessage(" -> global cut pool detected cutoff\n");
               }
            }
         }
         assert(lp->flushed);
         assert(lp->solved);
         assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);

         /* constraint separation */
         SCIPdebugMessage("constraint separation\n");

         /* separate constraints and LP */
         if( !(*cutoff) && !(*lperror) && !enoughcuts && lp->solved
            && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY) )
         {
            /* apply a separation round */
            SCIP_CALL( separationRoundLP(blkmem, set, stat, lp, sepastore, actdepth, delayedsepa,
                  &delayedsepa, &enoughcuts, cutoff) );

            /* if bound changes were applied in the separation round, we have to resolve the LP */
            if( !(*cutoff) && !lp->flushed )
            {
               /* solve LP (with dual simplex) */
               SCIPdebugMessage("separation: solve LP\n");
               SCIP_CALL( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, TRUE, FALSE, lperror) );
               assert(lp->flushed);
               assert(lp->solved || *lperror);
               delayedsepa = FALSE;
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
            SCIP_CALL( SCIPsepastoreClearCuts(sepastore, blkmem, set, lp) );
         }
         else
         {
            /* apply found cuts */
            SCIP_CALL( SCIPsepastoreApplyCuts(sepastore, blkmem, set, stat, tree, lp, branchcand, eventqueue, cutoff) );

            if( !(*cutoff) )
            {
               /* if a new bound change (e.g. a cut with only one column) was found, propagate domains again */
               if( stat->domchgcount != olddomchgcount )
               {
                  /* propagate domains */
                  SCIP_CALL( propagateDomains(blkmem, set, stat, tree, SCIPtreeGetCurrentDepth(tree), 0, FALSE, cutoff) );

                  /* in the root node, remove redundant rows permanently from the LP */
                  if( root )
                  {
                     SCIP_CALL( SCIPlpFlush(lp, blkmem, set) );
                     SCIP_CALL( SCIPlpRemoveRedundantRows(lp, blkmem, set, stat) );
                  }
               }

               mustprice = mustprice || !lp->flushed || (prob->ncolvars != npricedcolvars);
               mustsepa = !lp->flushed;

               if( !(*cutoff) )
               {
                  SCIP_Real lpobjval;

                  /* solve LP (with dual simplex) */
                  SCIPdebugMessage("separation: solve LP\n");
                  SCIP_CALL( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, TRUE, FALSE, lperror) );
                  assert(lp->flushed);
                  assert(lp->solved || *lperror);

                  if( !(*lperror) && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
                  {
                     SCIP_Real objreldiff;
                     int nfracs;

		     if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
		     {
			SCIP_CALL( SCIPbranchcandGetLPCands(branchcand, set, stat, lp, NULL, NULL, NULL, &nfracs, NULL) );
		     }
		     else
			nfracs = INT_MAX;
                     lpobjval = SCIPlpGetObjval(lp, set);
                     objreldiff = SCIPrelDiff(lpobjval, stalllpobjval);
                     if( objreldiff > 1e-03 || nfracs <= (0.9 - 0.1 * nsepastallrounds) * stallnfracs )
                     {
                        nsepastallrounds = 0;
                        stalllpobjval = lpobjval;
                        stallnfracs = nfracs;
                     }
                     else
                        nsepastallrounds++;
                  }
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

      SCIP_CALL( nodeUpdateLowerboundLP(focusnode, set, stat, lp) );

      /* issue LPSOLVED event */
      SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_LPSOLVED) );
      SCIP_CALL( SCIPeventChgNode(&event, focusnode) );
      SCIP_CALL( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );

      /* analyze an infeasible LP (not necessary in the root node) */
      if( !set->misc_exactsolve && !root
         && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT) )
      {
         SCIP_CALL( SCIPconflictAnalyzeLP(conflict, blkmem, set, stat, prob, tree, lp, NULL) );
      }

      /* check for unboundness */
      if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
      {
         assert(root); /* this can only happen in the root node */
         *unbounded = TRUE;
      }
   }
   SCIPdebugMessage(" -> final lower bound: %g (LP status: %d, LP obj: %g)\n",
      SCIPnodeGetLowerbound(focusnode), SCIPlpGetSolstat(lp),
      *cutoff ? SCIPsetInfinity(set) : *lperror ? -SCIPsetInfinity(set) : SCIPlpGetObjval(lp, set));

   return SCIP_OKAY;
}

/** solves the current node's LP in a price-and-cut loop */
static
SCIP_RETCODE solveNodeLP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cut pool */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool             initiallpsolved,    /**< was the initial LP already solved? */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            unbounded,          /**< pointer to store TRUE, if an unbounded ray was found in the LP */
   SCIP_Bool*            lperror             /**< pointer to store TRUE, if an unresolved error in LP solving occured */
   )
{
   SCIP_Longint nlpiterations;
   int nlps;

   assert(stat != NULL);
   assert(tree != NULL);
   assert(SCIPtreeHasFocusNodeLP(tree));
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
      SCIP_CALL( solveNodeInitialLP(blkmem, set, stat, prob, tree, lp, pricestore, sepastore,
            branchcand, eventfilter, eventqueue, cutoff, lperror) );
      assert(*cutoff || *lperror || (lp->flushed && lp->solved));
      SCIPdebugMessage("price-and-cut-loop: initial LP status: %d, LP obj: %g\n",
         SCIPlpGetSolstat(lp),
         *cutoff ? SCIPsetInfinity(set) : *lperror ? -SCIPsetInfinity(set) : SCIPlpGetObjval(lp, set));

      /* update initial LP iteration counter */
      stat->ninitlps += stat->nlps - nlps;
      stat->ninitlpiterations += stat->nlpiterations - nlpiterations;
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);

   if( !(*cutoff) && !(*lperror) )
   {
      /* solve the LP with price-and-cut*/
      SCIP_CALL( priceAndCutLoop(blkmem, set, stat, prob, primal, tree, lp, pricestore, sepastore, cutpool,
            branchcand, conflict, eventfilter, eventqueue, initiallpsolved, cutoff, unbounded, lperror) );
   }
   assert(*cutoff || *lperror || (lp->flushed && lp->solved));

   /* update node's LP iteration counter */
   stat->nnodelps += stat->nlps - nlps;
   stat->nnodelpiterations += stat->nlpiterations - nlpiterations;

   return SCIP_OKAY;
}

/** calls relaxators */
static
SCIP_RETCODE solveNodeRelax(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             beforelp,           /**< should the relaxators with non-negative or negative priority be called? */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            propagateagain,     /**< pointer to store TRUE, if domain propagation should be applied again */
   SCIP_Bool*            solvelpagain        /**< pointer to store TRUE, if the node's LP has to be solved again */
   )
{
   SCIP_RESULT result;
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

      SCIP_CALL( SCIPrelaxExec(set->relaxs[r], set, stat, depth, &result) );

      switch( result )
      {
      case SCIP_CUTOFF:
         *cutoff = TRUE;
         SCIPdebugMessage(" -> relaxator <%s> detected cutoff\n", SCIPrelaxGetName(set->relaxs[r]));
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
         SCIPerrorMessage("invalid result code <%d> of relaxator <%s>\n", result, SCIPrelaxGetName(set->relaxs[r]));
         return SCIP_INVALIDRESULT;
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** enforces constraints by branching, separation, or domain reduction */
static
SCIP_RETCODE enforceConstraints(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_Bool*            branched,           /**< pointer to store whether a branching was created */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the LP/pseudo solution is infeasible */
   SCIP_Bool*            propagateagain,     /**< pointer to store TRUE, if domain propagation should be applied again */
   SCIP_Bool*            solvelpagain        /**< pointer to store TRUE, if the node's LP has to be solved again */
   )
{
   SCIP_RESULT result;
   SCIP_Real pseudoobjval;
   SCIP_Bool resolved;
   SCIP_Bool objinfeasible;
   int h;

   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(SCIPtreeGetFocusNode(tree) != NULL);
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
   SCIPdebugMessage("enforcing constraints on %s solution\n", SCIPtreeHasFocusNodeLP(tree) ? "LP" : "pseudo");

   /* check, if the solution is infeasible anyway due to it's objective value */
   if( SCIPtreeHasFocusNodeLP(tree) )
      objinfeasible = FALSE;
   else
   {
      pseudoobjval = SCIPlpGetPseudoObjval(lp, set);
      objinfeasible = SCIPsetIsLT(set, pseudoobjval, SCIPnodeGetLowerbound(SCIPtreeGetFocusNode(tree)));
   }

   /* during constraint enforcemenst, generated cuts should enter the LP in any case; otherwise, a constraint handler
    * would fail to enforce its constraints if it relies on the modification of the LP relaxation
    */
   SCIPsepastoreStartForceCuts(sepastore);

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
         SCIP_CALL( SCIPconshdlrEnforceLPSol(set->conshdlrs_enfo[h], blkmem, set, stat, tree, sepastore, &result) );
      }
      else
      {
         SCIP_CALL( SCIPconshdlrEnforcePseudoSol(set->conshdlrs_enfo[h], blkmem, set, stat, tree, objinfeasible,
               &result) );
         if( SCIPsepastoreGetNCuts(sepastore) != 0 )
         {
            SCIPerrorMessage("pseudo enforcing method of constraint handler <%s> separated cuts\n",
               SCIPconshdlrGetName(set->conshdlrs_enfo[h]));
            return SCIP_INVALIDRESULT;
         }
      }
      SCIPdebugMessage("enforcing of <%s> returned result %d\n", SCIPconshdlrGetName(set->conshdlrs_enfo[h]), result);

      switch( result )
      {
      case SCIP_CUTOFF:
         assert(tree->nchildren == 0);
         *cutoff = TRUE;
         *infeasible = TRUE;
         resolved = TRUE;
         SCIPdebugMessage(" -> constraint handler <%s> detected cutoff in enforcement\n", 
            SCIPconshdlrGetName(set->conshdlrs_enfo[h]));
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
         SCIPerrorMessage("invalid result code <%d> from enforcing method of constraint handler <%s>\n",
            result, SCIPconshdlrGetName(set->conshdlrs_enfo[h]));
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
         SCIPdebugMessage(" -> LP exceeded objective limit\n");
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

   /* deactivate the cut forcing of the constraint enforcement */
   SCIPsepastoreEndForceCuts(sepastore);

   SCIPdebugMessage(" -> enforcing result: branched=%d, cutoff=%d, infeasible=%d, propagateagain=%d, solvelpagain=%d, resolved=%d\n",
      *branched, *cutoff, *infeasible, *propagateagain, *solvelpagain, resolved);

   return SCIP_OKAY;
}

/** applies the cuts stored in the separation store, or clears the store if the node can be cut off */
static
SCIP_RETCODE applyCuts(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool*            cutoff,             /**< pointer to whether the node can be cut off */
   SCIP_Bool*            propagateagain,     /**< pointer to store TRUE, if domain propagation should be applied again */
   SCIP_Bool*            solvelpagain        /**< pointer to store TRUE, if the node's LP has to be solved again */
   )
{
   assert(stat != NULL);
   assert(cutoff != NULL);
   assert(propagateagain != NULL);
   assert(solvelpagain != NULL);

   if( *cutoff )
   {
      /* the found cuts are of no use, because the node is infeasible anyway (or we have an error in the LP) */
      SCIP_CALL( SCIPsepastoreClearCuts(sepastore, blkmem, set, lp) );
   }
   else if( SCIPsepastoreGetNCuts(sepastore) > 0 )
   {
      SCIP_Longint olddomchgcount;

      olddomchgcount = stat->domchgcount;
      SCIP_CALL( SCIPsepastoreApplyCuts(sepastore, blkmem, set, stat, tree, lp, branchcand, eventqueue, cutoff) );
      *propagateagain = *propagateagain || (stat->domchgcount != olddomchgcount);
      *solvelpagain = TRUE;
   }

   return SCIP_OKAY;
}

/** updates the cutoff, propagateagain, and solverelaxagain status of the current solving loop */
static
void updateLoopStatus(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   depth,              /**< depth of current node */
   SCIP_Bool*            cutoff,             /**< pointer to whether the node can be cut off */
   SCIP_Bool*            propagateagain,     /**< pointer to store TRUE, if domain propagation should be applied again */
   SCIP_Bool*            solverelaxagain     /**< pointer to store TRUE, if at least one relaxator should be called again */
   )
{
   SCIP_NODE* focusnode;
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
SCIP_RETCODE solveNode(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cut pool */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            unbounded,          /**< pointer to store whether the focus node is unbounded */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the focus node's solution is infeasible */
   SCIP_Bool*            restart             /**< should solving process be started again with presolving? */
   )
{
   SCIP_NODE* focusnode;
   SCIP_Longint lastdomchgcount;
   SCIP_Real restartfac;
   int lastlpcount;
   int actdepth;
   int nlperrors;
   SCIP_Bool focusnodehaslp;
   SCIP_Bool initiallpsolved;
   SCIP_Bool solverelaxagain;
   SCIP_Bool solvelpagain;
   SCIP_Bool propagateagain;
   SCIP_Bool fullpropagation;
   SCIP_Bool branched;
   SCIP_Bool forcedlpsolve;

   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);
   assert(primal != NULL);
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);
   assert(SCIPconflictGetNConflicts(conflict) == 0);
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

   SCIPdebugMessage("Processing node %"SCIP_LONGINT_FORMAT" in depth %d, %d siblings\n",
      stat->nnodes, actdepth, tree->nsiblings);
   SCIPdebugMessage("current pseudosolution: obj=%g\n", SCIPlpGetPseudoObjval(lp, set));
   /*debug(SCIPprobPrintPseudoSol(prob, set));*/

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
    *  - solve SCIP_LP
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
   fullpropagation = TRUE;
   forcedlpsolve = FALSE;
   while( !(*cutoff) && (solverelaxagain || solvelpagain || propagateagain) && nlperrors < MAXNLPERRORS && !(*restart) )
   {
      SCIP_Real pseudoobjval;
      SCIP_Bool lperror;
      SCIP_Bool solverelax;
      SCIP_Bool solvelp;
      SCIP_Bool propagate;
      int r;

      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      lperror = FALSE;
      solverelax = solverelaxagain;
      solverelaxagain = FALSE;
      solvelp = solvelpagain;
      solvelpagain = FALSE;
      propagate = propagateagain;
      propagateagain = FALSE;

      /* domain propagation */
      if( propagate && !(*cutoff) )
      {
         SCIP_Bool lpwasflushed;

         lpwasflushed = lp->flushed;

         SCIP_CALL( propagateDomains(blkmem, set, stat, tree, SCIPtreeGetCurrentDepth(tree), 0, fullpropagation, cutoff) );
         fullpropagation = FALSE;

         /* check, if the path was cutoff */
         *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);

         /* if the LP was flushed and is now no longer flushed, a bound change occurred, and the LP has to be resolved */
         solvelp = solvelp || (lpwasflushed && !lp->flushed);
      }
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      /* if the LP should be resolved, all relaxations should also be resolved */
      /**@todo if the LP modification methods of the relax interface is implemented, we can remove this and give
       *       total control to the relaxators
       */
      if( solvelp )
      {
         solverelax = TRUE;
         for( r = 0; r < set->nrelaxs; ++r )
            SCIPrelaxMarkUnsolved(set->relaxs[r]);
      }

      /* solve external relaxations with non-negative priority */
      if( solverelax && !(*cutoff) )
      {
         SCIP_CALL( solveNodeRelax(set, stat, actdepth, TRUE, cutoff, &propagateagain, &solvelpagain) );

         /* check, if the path was cutoff */
         *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);

         /* apply found cuts */
         SCIP_CALL( applyCuts(blkmem, set, stat, tree, lp, sepastore, branchcand, eventqueue,
               cutoff, &propagateagain, &solvelpagain) );
      }
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      /* check, if we want to solve the LP at this node */
      if( solvelp && !(*cutoff) )
      {
         if( SCIPtreeHasFocusNodeLP(tree) )
         {
            /* solve the node's LP */
            SCIP_CALL( solveNodeLP(blkmem, set, stat, prob, primal, tree, lp, pricestore, sepastore,
                  cutpool, branchcand, conflict, eventfilter, eventqueue, initiallpsolved, cutoff, unbounded, &lperror) );
            initiallpsolved = TRUE;
            SCIPdebugMessage(" -> LP status: %d, LP obj: %g, iter: %"SCIP_LONGINT_FORMAT", count: %d\n",
               SCIPlpGetSolstat(lp),
               *cutoff ? SCIPsetInfinity(set) : lperror ? -SCIPsetInfinity(set) : SCIPlpGetObjval(lp, set),
               stat->nlpiterations, stat->lpcount);

            /* check, if the path was cutoff */
            *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);

            /* if an error occured during LP solving, switch to pseudo solution */
            if( lperror )
            {
               if( forcedlpsolve )
               {
                  SCIPerrorMessage("(node %"SCIP_LONGINT_FORMAT") unresolved numerical troubles in LP %d cannot be dealt with\n",
                     stat->nnodes, stat->nlps);
                  return SCIP_LPERROR;
               }
               SCIPtreeSetFocusNodeLP(tree, FALSE);
               nlperrors++;
               SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %"SCIP_LONGINT_FORMAT") unresolved numerical troubles in LP %d -- using pseudo solution instead (loop %d)\n",
                  stat->nnodes, stat->nlps, nlperrors);
            }

            /* if we solve exactly, the LP claims to be infeasible but the infeasibility could not be proved,
             * we have to forget about the LP and use the pseudo solution instead
             */
            if( !(*cutoff) && !lperror && set->misc_exactsolve && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE
               && SCIPnodeGetLowerbound(focusnode) < primal->cutoffbound )
            {
               if( SCIPbranchcandGetNPseudoCands(branchcand) == 0 && prob->ncontvars > 0 )
               {
                  SCIPerrorMessage("(node %"SCIP_LONGINT_FORMAT") could not prove infeasibility of LP %d, all variables are fixed, %d continuous vars\n",
                     stat->nnodes, stat->nlps, prob->ncontvars);
                  SCIPerrorMessage("(node %"SCIP_LONGINT_FORMAT")  -> have to call PerPlex() (feature not yet implemented)\n", stat->nnodes);
                  /**@todo call PerPlex */
                  return SCIP_LPERROR;
               }
               else
               {
                  SCIPtreeSetFocusNodeLP(tree, FALSE);
                  SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                     "(node %"SCIP_LONGINT_FORMAT") could not prove infeasibility of LP %d -- using pseudo solution (%d unfixed vars) instead\n",
                     stat->nnodes, stat->nlps, SCIPbranchcandGetNPseudoCands(branchcand));
               }
            }
         }
      }
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);
      assert(*cutoff || !SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));

      /* solve external relaxations with negative priority */
      if( solverelax && !(*cutoff) )
      {
         SCIP_CALL( solveNodeRelax(set, stat, actdepth, FALSE, cutoff, &propagateagain, &solvelpagain) );

         /* check, if the path was cutoff */
         *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);

         /* apply found cuts */
         SCIP_CALL( applyCuts(blkmem, set, stat, tree, lp, sepastore, branchcand, eventqueue,
               cutoff, &propagateagain, &solvelpagain) );
      }
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      /* update lower bound w.r.t. the pseudo solution */
      pseudoobjval = SCIPlpGetPseudoObjval(lp, set);
      SCIPnodeUpdateLowerbound(focusnode, stat, pseudoobjval);
      SCIPdebugMessage(" -> new lower bound: %g [%g] (pseudoobj: %g [%g])\n",
         SCIPnodeGetLowerbound(focusnode), SCIPprobExternObjval(prob, set, SCIPnodeGetLowerbound(focusnode)),
         pseudoobjval, SCIPprobExternObjval(prob, set, pseudoobjval));

      /* check for infeasible node by bounding */
      if( !(*cutoff)
         && (SCIPnodeGetLowerbound(focusnode) >= primal->cutoffbound
            || (!set->misc_exactsolve && SCIPsetIsGE(set, SCIPnodeGetLowerbound(focusnode), primal->cutoffbound))) )
      {
         SCIPdebugMessage("node is cut off by bounding (lower=%g, upper=%g)\n",
            SCIPnodeGetLowerbound(focusnode), primal->cutoffbound);
         SCIPnodeUpdateLowerbound(focusnode, stat, primal->cutoffbound);
         *cutoff = TRUE;

         /* call pseudo conflict analysis, if the node is cut off due to the pseudo objective value */
         if( pseudoobjval >= primal->cutoffbound )
         {
            SCIP_CALL( SCIPconflictAnalyzePseudo(conflict, blkmem, set, stat, prob, tree, lp, NULL) );
         }
      }

      /* update the cutoff, propagateagain, and solverelaxagain status of current solving loop */
      updateLoopStatus(set, stat, tree, actdepth, cutoff, &propagateagain, &solverelaxagain);

      /* enforce constraints */
      branched = FALSE;
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
         SCIP_CALL( enforceConstraints(blkmem, set, stat, tree, lp, sepastore,
               &branched, cutoff, infeasible, &propagateagain, &solvelpagain) );
         assert(branched == (tree->nchildren > 0));
         assert(!branched || (!(*cutoff) && *infeasible && !propagateagain && !solvelpagain));
         assert(!(*cutoff) || (!branched && *infeasible && !propagateagain && !solvelpagain));
         assert(*infeasible || (!branched && !(*cutoff) && !propagateagain && !solvelpagain));
         assert(!propagateagain || (!branched && !(*cutoff) && *infeasible));
         assert(!solvelpagain || (!branched && !(*cutoff) && *infeasible));

         /* apply found cuts */
         SCIP_CALL( applyCuts(blkmem, set, stat, tree, lp, sepastore, branchcand, eventqueue,
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
      forcedlpsolve = FALSE;
      if( *infeasible && !(*cutoff) && !(*unbounded) && !solverelaxagain && !solvelpagain && !propagateagain && !branched )
      {
         SCIP_RESULT result;
         int nlpcands;

         if( SCIPtreeHasFocusNodeLP(tree) )
         {
            SCIP_CALL( SCIPbranchcandGetLPCands(branchcand, set, stat, lp, NULL, NULL, NULL, &nlpcands, NULL) );
         }
         else
            nlpcands = 0;

         if( nlpcands > 0 )
         {
            /* branch on LP solution */
            SCIPdebugMessage("infeasibility in depth %d was not resolved: branch on LP solution with %d fractionals\n",
               SCIPnodeGetDepth(focusnode), nlpcands);
            SCIP_CALL( SCIPbranchExecLP(blkmem, set, stat, tree, lp, sepastore, branchcand, eventqueue,
                  primal->cutoffbound, FALSE, &result) );
            assert(result != SCIP_DIDNOTRUN);
         }
         else
         {
            /* branch on pseudo solution */
            SCIPdebugMessage("infeasibility in depth %d was not resolved: branch on pseudo solution with %d unfixed integers\n",
               SCIPnodeGetDepth(focusnode), SCIPbranchcandGetNPseudoCands(branchcand));
            SCIP_CALL( SCIPbranchExecPseudo(blkmem, set, stat, tree, lp, branchcand, eventqueue,
                  primal->cutoffbound, TRUE, &result) );
         }

         switch( result )
         {
         case SCIP_CUTOFF:
            assert(tree->nchildren == 0);
            *cutoff = TRUE;
            SCIPdebugMessage(" -> branching rule detected cutoff\n");
            break;
         case SCIP_CONSADDED:
            assert(tree->nchildren == 0);
            if( nlpcands > 0 )
            {
               SCIPerrorMessage("LP branching rule added constraint, which was not allowed this time\n");
               return SCIP_INVALIDRESULT;
            }
            solverelaxagain = TRUE;
            solvelpagain = TRUE;
            propagateagain = TRUE;
            break;
         case SCIP_REDUCEDDOM:
            assert(tree->nchildren == 0);
            solverelaxagain = TRUE;
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
            {
               *cutoff = TRUE;
               SCIPdebugMessage(" -> cutoff because all variables are fixed in current node\n");
            }
            else
            {
               assert(!SCIPtreeHasFocusNodeLP(tree)); /* feasible LP solutions with all integers fixed must be feasible */

               /* solve the LP in the next loop */
               SCIPtreeSetFocusNodeLP(tree, TRUE);
               solvelpagain = TRUE;
               forcedlpsolve = TRUE; /* this LP must be solved without error - otherwise we have to abort */
            }
            break;
         default:
            SCIPerrorMessage("invalid result code <%d> from SCIPbranchLP() or SCIPbranchPseudo()\n", result);
            return SCIP_INVALIDRESULT;
         }  /*lint !e788*/
         assert(*cutoff || solvelpagain || propagateagain || branched); /* something must have been done */
         assert(!(*cutoff) || (!solvelpagain && !propagateagain && !branched));
         assert(!solvelpagain || (!(*cutoff) && !branched));
         assert(!propagateagain || (!(*cutoff) && !branched));
         assert(!branched || (!solvelpagain && !propagateagain));
         assert(branched == (tree->nchildren > 0));

         /* apply found cuts */
         SCIP_CALL( applyCuts(blkmem, set, stat, tree, lp, sepastore, branchcand, eventqueue,
               cutoff, &propagateagain, &solvelpagain) );

         /* update the cutoff, propagateagain, and solverelaxagain status of current solving loop */
         updateLoopStatus(set, stat, tree, actdepth, cutoff, &propagateagain, &solverelaxagain);
      }

      /* check for immediate restart */
      *restart = *restart
         || (actdepth == 0
            && (set->presol_maxrestarts == -1 || stat->nruns <= set->presol_maxrestarts)
            && set->presol_restartfac > 0.0
            && stat->nrootintfixingsrun > set->presol_restartfac * (prob->nvars - prob->ncontvars));

      SCIPdebugMessage("node solving iteration finished: cutoff=%d, propagateagain=%d, solverelaxagain=%d, solvelpagain=%d, nlperrors=%d, restart=%d\n",
         *cutoff, propagateagain, solverelaxagain, solvelpagain, nlperrors, *restart);
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);
   assert(*cutoff || SCIPconflictGetNConflicts(conflict) == 0);

   /* flush the conflict set storage */
   SCIP_CALL( SCIPconflictFlushConss(conflict, blkmem, set, stat, prob, tree) );

   /* check for too many LP errors */
   if( nlperrors >= MAXNLPERRORS )
   {
      SCIPerrorMessage("(node %"SCIP_LONGINT_FORMAT") unresolved numerical troubles in LP %d -- aborting\n", stat->nnodes, stat->nlps);
      return SCIP_LPERROR;
   }

   /* check for final restart */
   restartfac = set->presol_subrestartfac;
   if( actdepth == 0 )
      restartfac = MIN(restartfac, set->presol_restartfac);
   *restart = *restart
      || ((set->presol_maxrestarts == -1 || stat->nruns <= set->presol_maxrestarts)
         && stat->nrootintfixingsrun > restartfac * (prob->nvars - prob->ncontvars));

   /* remember root LP solution */
   if( actdepth == 0 && !(*cutoff) && !(*restart) && !(*unbounded) )
   {
      SCIPprobStoreRootSol(prob, set, stat, lp, SCIPtreeHasFocusNodeLP(tree));
   }

   /* check for cutoff */
   if( *cutoff )
   {
      SCIPdebugMessage("node is cut off\n");
      SCIPnodeUpdateLowerbound(focusnode, stat, primal->cutoffbound);
      *infeasible = TRUE;
      *restart = FALSE;
      *unbounded = FALSE;
   }

   return SCIP_OKAY;
}

/** if feasible, adds current solution to the solution storage */
static
SCIP_RETCODE addCurrentSolution(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_EVENTFILTER*     eventfilter         /**< event filter for global (not variable dependent) events */
   )
{
   SCIP_SOL* sol;
   SCIP_Bool foundsol;

   /* found a feasible solution */
   if( SCIPtreeHasFocusNodeLP(tree) )
   {
      /* start clock for LP solutions */
      SCIPclockStart(stat->lpsoltime, set);

      /* add solution to storage */
      SCIP_CALL( SCIPsolCreateLPSol(&sol, blkmem, set, stat, primal, tree, lp, NULL) );
      if( set->misc_exactsolve )
      {
         /* if we want to solve exactly, we have to check the solution exactly again */
         SCIP_CALL( SCIPprimalTrySolFree(primal, blkmem, set, stat, prob, tree, lp, eventfilter, &sol,
               TRUE, TRUE, TRUE, &foundsol) );
      }
      else
      {
         SCIP_CALL( SCIPprimalAddSolFree(primal, blkmem, set, stat, prob, tree, lp, eventfilter, &sol, &foundsol) );
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
      SCIP_CALL( SCIPsolCreatePseudoSol(&sol, blkmem, set, stat, primal, tree, lp, NULL) );
      if( set->misc_exactsolve )
      {
         /* if we want to solve exactly, we have to check the solution exactly again */
         SCIP_CALL( SCIPprimalTrySolFree(primal, blkmem, set, stat, prob, tree, lp, eventfilter, &sol,
               TRUE, TRUE, TRUE, &foundsol) );
      }
      else
      {
         SCIP_CALL( SCIPprimalAddSolFree(primal, blkmem, set, stat, prob, tree, lp, eventfilter, &sol, &foundsol) );
      }

      /* stop clock for pseudo solutions */
      SCIPclockStop(stat->pseudosoltime, set);

      if( foundsol )
         stat->npssolsfound++;
   }

   return SCIP_OKAY;
}

/** main solving loop */
SCIP_RETCODE SCIPsolveCIP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_MEM*             mem,                /**< block memory pools */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cut pool */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool*            restart             /**< should solving process be started again with presolving? */
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODE* focusnode;
   SCIP_NODE* nextnode;
   SCIP_EVENT event;
   SCIP_Real restartfac;
   int nnodes;
   int depth;
   SCIP_Bool cutoff;
   SCIP_Bool unbounded;
   SCIP_Bool infeasible;
   SCIP_Bool foundsol;

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
   restartfac = set->presol_subrestartfac;
   if( SCIPtreeGetCurrentDepth(tree) == 0 )
      restartfac = MIN(restartfac, set->presol_restartfac);
   *restart = (set->presol_maxrestarts == -1 || stat->nruns <= set->presol_maxrestarts)
      && stat->nrootintfixingsrun > restartfac * (prob->nvars - prob->ncontvars);

   /* switch status to UNKNOWN */
   stat->status = SCIP_STATUS_UNKNOWN;

   nextnode = NULL;
   unbounded = FALSE;

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
         SCIP_CALL( SCIPtreeSetNodesel(tree, set, stat, nodesel) );

         /* the next node was usually already selected in the previous solving loop before the primal heuristics were
          * called, because they need to know, if the next node will be a child/sibling (plunging) or not;
          * if the heuristics found a new best solution that cut off some of the nodes, the node selector must be called
          * again, because the selected next node may be invalid due to cut off
          */
         if( nextnode == NULL )
         {
            /* select next node to process */
            SCIP_CALL( SCIPnodeselSelect(nodesel, set, &nextnode) );
         }
         focusnode = nextnode;
         nextnode = NULL;
         assert(SCIPbufferGetNUsed(set->buffer) == 0);

         /* start node activation timer */
         SCIPclockStart(stat->nodeactivationtime, set);

         /* focus selected node */
         SCIP_CALL( SCIPnodeFocus(&focusnode, blkmem, set, stat, prob, primal, tree, lp, branchcand, conflict,
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
      SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEFOCUSED) );
      SCIP_CALL( SCIPeventChgNode(&event, focusnode) );
      SCIP_CALL( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );

      /* call primal heuristics that should be applied before the node was solved */
      SCIP_CALL( primalHeuristics(set, primal, tree, nextnode, FALSE, FALSE, &foundsol) );
      assert(SCIPbufferGetNUsed(set->buffer) == 0);

      /* solve focus node */
      SCIP_CALL( solveNode(blkmem, set, stat, prob, primal, tree, lp, pricestore, sepastore, branchcand, cutpool,
            conflict, eventfilter, eventqueue, &cutoff, &unbounded, &infeasible, restart) );
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
            SCIP_CALL( addCurrentSolution(blkmem, set, stat, prob, primal, tree, lp, eventfilter) );

            /* issue NODEFEASIBLE event */
            SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEFEASIBLE) );
            SCIP_CALL( SCIPeventChgNode(&event, focusnode) );
            SCIP_CALL( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );
         }
         else if( !unbounded )
         {
            /* node solution is not feasible */
            if( tree->nchildren == 0 )
            {
               /* issue NODEINFEASIBLE event */
               SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEINFEASIBLE) );

               /* increase the cutoff counter of the branching variable */
               if( stat->lastbranchvar != NULL )
               {
                  SCIP_CALL( SCIPvarIncNCutoffs(stat->lastbranchvar, stat, stat->lastbranchdir) );
               }
               /**@todo if last branching variable is unknown, retrieve it from the nodes' boundchg arrays */
            }
            else
            {
               /* issue NODEBRANCHED event */
               SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEBRANCHED) );
            }
            SCIP_CALL( SCIPeventChgNode(&event, focusnode) );
            SCIP_CALL( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );
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
            SCIP_RESULT result;

            assert(set->misc_exactsolve);

            do
            {
               result = SCIP_DIDNOTRUN;
               if( SCIPbranchcandGetNPseudoCands(branchcand) == 0 )
               {
                  if( prob->ncontvars > 0 )
                  {
                     /**@todo call PerPlex */
                     SCIPerrorMessage("cannot branch on all-fixed LP -- have to call PerPlex instead\n");
                  }
               }
               else
               {
                  SCIP_CALL( SCIPbranchExecPseudo(blkmem, set, stat, tree, lp, branchcand, eventqueue,
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
         SCIP_CALL( SCIPnodeselSelect(nodesel, set, &nextnode) );
         assert(SCIPbufferGetNUsed(set->buffer) == 0);

         /* call primal heuristics that should be applied after the relaxation was solved */
         nnodes = SCIPtreeGetNNodes(tree);
         SCIP_CALL( primalHeuristics(set, primal, tree, nextnode, TRUE, FALSE, &foundsol) );
         assert(SCIPbufferGetNUsed(set->buffer) == 0);

         /* if the heuristics found a new best solution that cut off some of the nodes, the node selector must be called
          * again, because the selected next node may be invalid due to cut off
          */
         assert(!tree->cutoffdelayed);
         if( nnodes != SCIPtreeGetNNodes(tree) )
            nextnode = NULL;
      }

      /* display node information line */
      SCIP_CALL( SCIPdispPrintLine(set, stat, NULL, (SCIPnodeGetDepth(focusnode) == 0) && infeasible && !foundsol) );

      SCIPdebugMessage("Processing of node %"SCIP_LONGINT_FORMAT" in depth %d finished. %d siblings, %d children, %d leaves left\n",
         stat->nnodes, SCIPnodeGetDepth(focusnode), tree->nsiblings, tree->nchildren, SCIPtreeGetNLeaves(tree));
      SCIPdebugMessage("**********************************************************************\n");
   }
   assert(SCIPbufferGetNUsed(set->buffer) == 0);

   SCIPdebugMessage("Problem solving finished (restart=%d)\n", *restart);

   /* if the current node is the only remaining node, and if its lower bound exceeds the upper bound, we have
    * to delete it manually in order to get to the SOLVED stage instead of thinking, that only the gap limit
    * was reached (this may happen, if the current node is the one defining the global lower bound and a
    * feasible solution with the same value was found at this node)
    */
   if( tree->focusnode != NULL && SCIPtreeGetNNodes(tree) == 0
      && SCIPsetIsGE(set, tree->focusnode->lowerbound, primal->cutoffbound) )
   {
      focusnode = NULL;
      SCIP_CALL( SCIPnodeFocus(&focusnode, blkmem, set, stat, prob, primal, tree, lp, branchcand, conflict,
            eventfilter, eventqueue, &cutoff) );
   }

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
