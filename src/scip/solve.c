/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: solve.c,v 1.132 2004/09/07 18:22:20 bzfpfend Exp $"

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
#include "sepa.h"


#define MAXNLPERRORS  10                /**< maximal number of LP error loops in a single node */


/** returns whether the solving process will be / was stopped before proving optimality */
Bool SCIPsolveIsStopped(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   )
{
   return (SCIPinterrupted()
      || (set->nodelimit >= 0 && stat->nnodes >= set->nodelimit)
      || SCIPclockGetTime(stat->solvingtime) >= set->timelimit
      || (SCIPgetMemUsed(set->scip) >= set->memlimit*1024.0*1024.0)
      || (SCIPstage(set->scip) >= SCIP_STAGE_SOLVING && SCIPsetIsLT(set, SCIPgetGap(set->scip), set->gaplimit))
      || (set->sollimit >= 0 && SCIPstage(set->scip) >= SCIP_STAGE_PRESOLVED &&
         SCIPgetNSolsFound(set->scip) >= set->sollimit)
      || (set->bestsollimit >= 0 && SCIPstage(set->scip) >= SCIP_STAGE_PRESOLVED &&
         SCIPgetNBestSolsFound(set->scip) >= set->bestsollimit));
}

/** outputs the reason for termination */
void SCIPsolvePrintStopReason(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   if( file == NULL )
      file = stdout;

   if( SCIPinterrupted() )
      fprintf(file, "user interrupt");
   else if( set->nodelimit >= 0 && stat->nnodes >= set->nodelimit )
      fprintf(file, "node limit reached");
   else if( SCIPclockGetTime(stat->solvingtime) >= set->timelimit )
      fprintf(file, "time limit reached");
   else if( SCIPgetMemUsed(set->scip) >= set->memlimit*1024.0*1024.0 )
      fprintf(file, "memory limit reached");
   else if( SCIPstage(set->scip) >= SCIP_STAGE_SOLVING && SCIPsetIsLT(set, SCIPgetGap(set->scip), set->gaplimit) )
      fprintf(file, "gap limit reached");
   else if( set->sollimit >= 0 && SCIPstage(set->scip) >= SCIP_STAGE_PRESOLVED
      && SCIPgetNSolsFound(set->scip) >= set->sollimit )
      fprintf(file, "solution limit reached");
   else if( set->bestsollimit >= 0 && SCIPstage(set->scip) >= SCIP_STAGE_PRESOLVED
      && SCIPgetNBestSolsFound(set->scip) >= set->bestsollimit )
      fprintf(file, "solution improvement limit reached");
}

/** applies domain propagation on current node */
RETCODE SCIPpropagateDomains(
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   NODE* node;
   RESULT result;
   Bool propagain;
   int maxproprounds;
   int propround;
   int h;

   assert(set != NULL);
   assert(tree != NULL);
   assert(cutoff != NULL);

   debugMessage("domain propagation\n");

   node = SCIPtreeGetCurrentNode(tree);
   assert(node != NULL);
   assert(SCIPnodeIsActive(node));

   maxproprounds = (SCIPnodeGetDepth(node) == 0 ? set->maxproproundsroot : set->maxproprounds);

   *cutoff = FALSE;
   propround = 0;
   do
   {
      propround++;
      propagain = FALSE;
      for( h = 0; h < set->nconshdlrs && !(*cutoff); ++h )
      {
         CHECK_OKAY( SCIPconshdlrPropagate(set->conshdlrs[h], memhdr, set, stat, prob, SCIPnodeGetDepth(node),
               &result) );
         propagain = propagain || (result == SCIP_REDUCEDDOM);
         *cutoff = *cutoff || (result == SCIP_CUTOFF);
      }
   }
   while( propagain && !(*cutoff) && propround < maxproprounds );

   return SCIP_OKAY;
}

/** strengthens variable's bounds looking at reduced costs */
static
RETCODE redcostStrengthening(
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
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
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(primal != NULL);
   assert(SCIPsetIsLT(set, SCIPlpGetObjval(lp, set), primal->cutoffbound) );

   /* we cannot apply reduced cost fixing, if we want to solve exactly */
   /**@todo implement reduced cost fixing with interval arithmetics */
   if( set->exactsolve )
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
   node = SCIPtreeGetActiveNode(tree);
   assert(node != NULL);
   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_ACTNODE);

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
                  CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, tree, lp, branchcand, eventqueue,
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
                  CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, tree, lp, branchcand, eventqueue,
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
   NODE* actnode;
   int actdepth;

   assert(lp != NULL);
   assert(tree != NULL);
   assert(tree->path != NULL);

   actnode = SCIPtreeGetActiveNode(tree);
   assert(SCIPnodeIsActive(actnode));
   assert(SCIPnodeGetType(actnode) == SCIP_NODETYPE_ACTNODE);
   actdepth = SCIPnodeGetDepth(actnode);
   assert(tree->path[actdepth] == actnode);

   if( lp->solved && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && tree->actlpfork != NULL )
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

      assert(SCIPnodeIsActive(tree->actlpfork));
      assert(tree->path[tree->actlpfork->depth] == tree->actlpfork);

      /* get a buffer for the collected bound changes; start with a size twice as large as the number of nodes between
       * current node and LP fork
       */
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &updates, 2*(actdepth - tree->actlpfork->depth)) );
      nupdates = 0;
      nvalidupdates = 0;

      /* search the nodes from LP fork down to current node for bound changes in between; move in this direction, 
       * because the bound changes closer to the LP fork are more likely to have a valid LP solution information
       * attached; collect the bound changes for pseudo cost value updates and mark the corresponding variables such
       * that they are not updated twice in case of more than one bound change on the same variable
       */
      for( d = tree->actlpfork->depth+1; d <= actdepth; ++d )
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
      lpgain = SCIPlpGetObjval(lp, set) - tree->actlpfork->lowerbound;
      lpgain = MAX(lpgain, 0.0);
      for( i = 0; i < nupdates; ++i )
      {
         assert(updates[i]->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING);
         var = updates[i]->var;
         assert(var->pseudocostflag != PSEUDOCOST_NONE);
         if( var->pseudocostflag == PSEUDOCOST_UPDATE )
         {
            debugMessage("updating pseudocosts of <%s>: sol: %g -> %g, LP: %e -> %e => gain=%g, weight: %g\n",
               SCIPvarGetName(var), updates[i]->data.branchingdata.lpsolval, SCIPvarGetLPSol(var),
               tree->actlpfork->lowerbound, SCIPlpGetObjval(lp, set), lpgain, weight);
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

   if( set->exactsolve )
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
   MEMHDR*          memhdr,             /**< block memory buffers */
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
         CHECK_OKAY( SCIPpricestoreAddVar(pricestore, memhdr, set, lp, var, 0.0, TRUE) );
      }
   }
   assert(lp->nremoveablecols == 0);
   CHECK_OKAY( SCIPpricestoreApplyVars(pricestore, memhdr, set, stat, prob, tree, lp) );

   /* add LP relaxations of all initial constraints to LP */
   debugMessage("init root LP: initial rows\n");
   for( h = 0; h < set->nconshdlrs; ++h )
   {
      CHECK_OKAY( SCIPconshdlrInitLP(set->conshdlrs[h], memhdr, set, stat, prob) );
   }
   CHECK_OKAY( SCIPsepastoreApplyCuts(sepastore, memhdr, set, stat, tree, lp, branchcand, eventqueue, cutoff) );

   /* inform pricing and separation storage, that initial LP setup is now finished */
   SCIPpricestoreEndInitialLP(pricestore);
   SCIPsepastoreEndInitialLP(sepastore);

   return SCIP_OKAY;
}

/** load and solve the initial LP of a node */
static
RETCODE solveNodeInitialLP(
   MEMHDR*          memhdr,             /**< block memory buffers */
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
   NODE* actnode;
   EVENT event;

   assert(stat != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(cutoff != NULL);
   assert(lperror != NULL);

   *cutoff = FALSE;
   *lperror = FALSE;

   actnode = SCIPtreeGetActiveNode(tree);
   assert(actnode != NULL);
   assert(SCIPnodeGetType(actnode) == SCIP_NODETYPE_ACTNODE);

   /* load the LP into the solver and load the LP state */
   debugMessage("loading LP\n");
   CHECK_OKAY( SCIPtreeLoadLP(tree, memhdr, set, stat, lp) );
   
   /* init root node LP */
   if( SCIPnodeGetDepth(actnode) == 0 )
   {
      CHECK_OKAY( initRootLP(memhdr, set, stat, prob, tree, lp, pricestore, sepastore, branchcand, eventqueue, cutoff) );
   }

   /* only the bounds and the constraint list may have been changed since the LP state node
    * -> dual simplex is applicable as first solver
    */
   debugMessage("node: solve initial LP\n");
   CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, -1, TRUE, lperror) );
   assert(lp->solved || *lperror);

   if( !(*lperror) )
   {
      /* issue FIRSTLPSOLVED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_FIRSTLPSOLVED) );
      CHECK_OKAY( SCIPeventChgNode(&event, actnode) );
      CHECK_OKAY( SCIPeventProcess(&event, memhdr, set, NULL, NULL, NULL, eventfilter) );
      
      /* update pseudo cost values */
      CHECK_OKAY( updatePseudocost(set, stat, tree, lp) );
   }

   return SCIP_OKAY;
}

/** solve the current LP of a node with a price-and-cut loop */
static
RETCODE priceAndCutLoop(
   MEMHDR*          memhdr,             /**< block memory buffers */
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
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   Bool*            lperror             /**< pointer to store whether an unresolved error in LP solving occured */
   )
{
   NODE* actnode;
   RESULT result;
   EVENT event;
   Bool mustprice;
   Bool mustsepar;
   Bool root;
   int actdepth;
   int ncolvars;
   int h;
   int s;

   assert(set != NULL);
   assert(memhdr != NULL);
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
   assert(lperror != NULL);

   actnode = SCIPtreeGetActiveNode(tree);
   assert(actnode != NULL);
   assert(SCIPnodeGetType(actnode) == SCIP_NODETYPE_ACTNODE);
   actdepth = SCIPnodeGetDepth(actnode);
   root = (actdepth == 0);

   /* solve initial LP of price-and-cut loop */
   debugMessage("node: solve LP with price and cut\n");
   CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, -1, TRUE, lperror) );
   assert(lp->solved);

   /* display node information line for root node */
   if( root && (VERBLEVEL)set->verblevel >= SCIP_VERBLEVEL_HIGH )
   {
      CHECK_OKAY( SCIPdispPrintLine(set, stat, TRUE) );
   }

   /* price-and-cut loop */
   ncolvars = prob->ncolvars;
   mustprice = TRUE;
   mustsepar = TRUE;
   *cutoff = FALSE;
   while( !(*cutoff) && !(*lperror) && (mustprice || mustsepar) )
   {
      debugMessage("-------- node solving loop --------\n");
      assert(lp->solved);

      /* if the LP is unbounded, we don't need to price */
      mustprice = mustprice && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDED);
      /* if all the variables are already in the LP, we don't need to price */
      mustprice = mustprice && (prob->nvars > SCIPlpGetNCols(lp) || set->nactivepricers > 0);

      /* pricing (has to be done completely to get a valid lower bound) */
      while( !(*cutoff) && !(*lperror) && mustprice )
      {
         assert(lp->solved);
         assert(lp->lpsolstat != SCIP_LPSOLSTAT_UNBOUNDED);

         /* price problem variables */
         debugMessage("problem variable pricing\n");
         assert(SCIPpricestoreGetNVars(pricestore) == 0);
         assert(SCIPpricestoreGetNBoundResets(pricestore) == 0);
         CHECK_OKAY( SCIPpricestoreAddProbVars(pricestore, memhdr, set, stat, prob, tree, lp, branchcand, eventqueue) );
         ncolvars = prob->ncolvars;

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
               enoughvars = enoughvars || (SCIPpricestoreGetNVars(pricestore) >= SCIPsetGetMaxpricevars(set, root)/2);
            }
         }

         /* apply the priced variables to the LP */
         CHECK_OKAY( SCIPpricestoreApplyVars(pricestore, memhdr, set, stat, prob, tree, lp) );
         assert(SCIPpricestoreGetNVars(pricestore) == 0);
         mustprice = !lp->solved || (prob->ncolvars != ncolvars);
         mustsepar = mustsepar || !lp->solved;
         
         /* after adding columns, the LP should be primal feasible such that primal simplex is applicable;
          * if LP was infeasible, we have to use dual simplex
          */
         debugMessage("pricing: solve LP\n");
         CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, -1, TRUE, lperror) );
         assert(lp->solved || *lperror);

         /* reset bounds temporarily set by pricer to their original values */
         debugMessage("pricing: reset bounds\n");
         CHECK_OKAY( SCIPpricestoreResetBounds(pricestore, memhdr, set, stat, lp, branchcand, eventqueue) );
         assert(SCIPpricestoreGetNVars(pricestore) == 0);
         assert(SCIPpricestoreGetNBoundResets(pricestore) == 0);
         mustprice = mustprice || !lp->solved || (prob->ncolvars != ncolvars);
         mustsepar = mustsepar || !lp->solved;

         /* solve LP again after resetting bounds (with dual simplex) */
         debugMessage("pricing: solve LP after resetting bounds\n");
         CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, -1, FALSE, lperror) );
         assert(lp->solved || *lperror);

         /* increase pricing round counter */
         stat->npricerounds++;

         /* display node information line for root node */
         if( root && mustprice && (VERBLEVEL)set->verblevel >= SCIP_VERBLEVEL_FULL )
         {
            CHECK_OKAY( SCIPdispPrintLine(set, stat, TRUE) );
         }

         /* if the LP is unbounded, we can stop pricing */
         mustprice = mustprice && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDED);
      }
      assert(lp->solved);

      /* update lower bound w.r.t. the the LP solution */
      if( !(*lperror) )
      {
         CHECK_OKAY( nodeUpdateLowerboundLP(actnode, set, stat, lp) );
         debugMessage(" -> new lower bound: %g (LP status: %d, LP obj: %g)\n",
            SCIPnodeGetLowerbound(actnode), SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set));
      }
      else
      {
         debugMessage(" -> error solving LP. keeping old bound: %g\n", SCIPnodeGetLowerbound(actnode));
      }

      /* if the LP is infeasible or exceeded the objective limit, we don't need to separate cuts */
      mustsepar = mustsepar
         && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_INFEASIBLE)
         && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OBJLIMIT)
         && SCIPsetIsLT(set, SCIPnodeGetLowerbound(actnode), primal->cutoffbound);

      /* separation and reduced cost strengthening
       * (needs not to be done completely, because we just want to increase the lower bound)
       */
      if( !(*cutoff) && !(*lperror) && mustsepar )
      {
         Bool separateagain;
         Bool enoughcuts;

         assert(lp->solved);
         assert(lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL);

         mustsepar = FALSE;
         enoughcuts = (SCIPsetGetMaxsepacuts(set, root) == 0);

         /* global cut pool separation */
         if( !enoughcuts )
         {
            debugMessage("global cut pool separation\n");
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            CHECK_OKAY( SCIPcutpoolSeparate(cutpool, memhdr, set, stat, lp, sepastore, root, &result) );
            *cutoff = *cutoff || (result == SCIP_CUTOFF);
            enoughcuts = enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetMaxsepacuts(set, root));
         }
         assert(lp->solved);
         assert(lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL);

         /* constraint separation */
         debugMessage("constraint separation\n");

         /* separate constraints and LP */
         separateagain = TRUE;
         while( !(*cutoff) && !(*lperror) && !enoughcuts && separateagain
            && lp->solved && lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
         {
            separateagain = FALSE;

            /* try separating constraints of the constraint handlers until the cut pool is at least half full;
             * we have to stop after a domain reduction was found, because they go directly into the LP and invalidate
             * the current solution
             */
            for( h = 0; h < set->nconshdlrs && !(*cutoff) && !enoughcuts && lp->solved; ++h )
            {
               CHECK_OKAY( SCIPconshdlrSeparate(conshdlrs_sepa[h], memhdr, set, stat, prob, sepastore, actdepth,
                     &result) );
               *cutoff = *cutoff || (result == SCIP_CUTOFF);
               separateagain = separateagain || (result == SCIP_CONSADDED);
               enoughcuts = enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetMaxsepacuts(set, root));
            }

            /* sort separators by priority */
            SCIPsetSortSepas(set);

            /* call LP separators */
            for( s = 0; s < set->nsepas && !(*cutoff) && !separateagain && !enoughcuts && lp->solved; ++s )
            {
               CHECK_OKAY( SCIPsepaExec(set->sepas[s], set, stat, sepastore, actdepth, &result) );
               *cutoff = *cutoff || (result == SCIP_CUTOFF);
               separateagain = separateagain || (result == SCIP_CONSADDED);
               enoughcuts = enoughcuts || (SCIPsepastoreGetNCutsStored(sepastore) >= 2*SCIPsetGetMaxsepacuts(set, root));
            }

            if( !(*cutoff) && !lp->solved )
            {
               /* solve LP (with dual simplex) */
               debugMessage("separation: solve LP\n");
               CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, -1, TRUE, lperror) );
               separateagain = TRUE;
            }
         }
         assert(*cutoff || *lperror || lp->solved);
         assert(!lp->solved
            || lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL
            || lp->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE
            || lp->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT);

         if( *cutoff || *lperror || SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OPTIMAL )
         {
            /* the found cuts are of no use, because the node is infeasible anyway (or we have an error in the LP) */
            CHECK_OKAY( SCIPsepastoreClearCuts(sepastore, memhdr, set, lp) );
         }
         else
         {
            Longint olddomchgcount;

            olddomchgcount = stat->domchgcount;

            /* apply reduced cost bound strengthening */
            if( lp->solved )
            {
               /* check, if we want to apply reduced cost fixing at current node */
               if( (set->redcostfreq == 0 && root) || (set->redcostfreq > 0 && actdepth % set->redcostfreq == 0) )
               {
                  CHECK_OKAY( redcostStrengthening(memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue) );
               }
            }
        
            /* apply found cuts */
            CHECK_OKAY( SCIPsepastoreApplyCuts(sepastore, memhdr, set, stat, tree, lp, branchcand, eventqueue, cutoff) );

            if( !(*cutoff) )
            {
               /* if at least one of the cut was a bound change, propagate domains again */
               if( stat->domchgcount != olddomchgcount )
               {
                  CHECK_OKAY( SCIPpropagateDomains(memhdr, set, stat, prob, tree, cutoff) );
               }
               
               mustprice = mustprice || !lp->solved || prob->ncolvars != ncolvars;
               mustsepar = !lp->solved;
               
               if( !(*cutoff) )
               {
                  /* solve LP (with dual simplex) */
                  debugMessage("separation: solve LP\n");
                  CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, -1, TRUE, lperror) );
               }
            }
         }
         assert(*cutoff || *lperror || lp->solved); /* after cutoff, the LP may be unsolved due to bound changes */

         /* increase separation round counter */
         stat->nseparounds++;

         /* display node information line for root node */
         if( root && mustsepar )
         {
            CHECK_OKAY( SCIPdispPrintLine(set, stat, TRUE) );
         }
      }
   }

   /* update lower bound w.r.t. the the LP solution */
   if( *cutoff )
   {
      SCIPnodeUpdateLowerbound(actnode, stat, primal->cutoffbound);
   }
   else if( !(*lperror) )
   {
      assert(lp->solved);

      CHECK_OKAY( nodeUpdateLowerboundLP(actnode, set, stat, lp) );

      /* issue LPSOLVED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_LPSOLVED) );
      CHECK_OKAY( SCIPeventChgNode(&event, actnode) );
      CHECK_OKAY( SCIPeventProcess(&event, memhdr, set, NULL, NULL, NULL, eventfilter) );

      /* analyze an infeasible LP (not necessary in the root node) */
      if( !set->exactsolve && !root
         && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT) )
      {
         CHECK_OKAY( SCIPconflictAnalyzeLP(conflict, memhdr, set, stat, prob, tree, lp, NULL) );
      }
   }
   debugMessage(" -> final lower bound: %g (LP status: %d, LP obj: %g)\n",
      SCIPnodeGetLowerbound(actnode), SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set));

   return SCIP_OKAY;
}

/** solves the current node's LP in a price-and-cut loop */
static
RETCODE solveNodeLP(
   MEMHDR*          memhdr,             /**< block memory buffers */
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
   Bool*            lperror             /**< pointer to store TRUE, if an unresolved error in LP solving occured */
   )
{
   Longint nlpiterations;
   int nlps;

   assert(stat != NULL);
   assert(tree != NULL);
   assert(SCIPtreeHasActiveNodeLP(tree));
   assert(cutoff != NULL);
   assert(lperror != NULL);
   assert(*cutoff == FALSE);
   assert(*lperror == FALSE);

   nlps = stat->nlps;
   nlpiterations = stat->nlpiterations;

   if( !initiallpsolved )
   {
      /* load and solve the initial LP of the node */
      CHECK_OKAY( solveNodeInitialLP(memhdr, set, stat, prob, tree, lp, pricestore, sepastore, 
            branchcand, eventfilter, eventqueue, cutoff, lperror) );
      assert(*cutoff || *lperror || lp->solved);
      debugMessage("price-and-cut-loop: initial LP status: %d, LP obj: %g\n", 
         SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set));
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);

   if( !(*cutoff) && !(*lperror) )
   {
      /* solve the LP with price-and-cut*/
      CHECK_OKAY( priceAndCutLoop(memhdr, set, stat, prob, primal, tree, lp, pricestore, sepastore, cutpool, 
            branchcand, conflict, eventfilter, eventqueue, conshdlrs_sepa, cutoff, lperror) );
   }
   assert(*cutoff || *lperror || lp->solved);

   /* update node's LP iteration counter */
   stat->nnodelps += stat->nlps - nlps;
   stat->nnodelpiterations += stat->nlpiterations - nlpiterations;

   return SCIP_OKAY;
}

/** enforces constraints by branching, separation, or domain reduction */
static
RETCODE enforceConstraints(
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   CONSHDLR**       conshdlrs_enfo,     /**< constraint handlers for enforcing constraints, sorted by priority */
   Bool*            cutoff,             /**< pointer to store whether node can be cut off */
   Bool*            infeasible,         /**< pointer to store whether LP/pseudo solution is infeasible */
   Bool*            solvelpagain,       /**< pointer to store whether the node's LP has to be solved again */
   Bool*            propagateagain      /**< pointer to store whether domain propagation should be applied again */
   )
{
   RESULT result;
   Real pseudoobjval;
   Longint olddomchgcount;
   Bool resolved;
   Bool objinfeasible;
   int h;

   assert(set != NULL);
   assert(tree != NULL);
   assert(SCIPtreeGetActiveNode(tree) != NULL);
   assert(prob != NULL);
   assert(primal != NULL);
   assert(conshdlrs_enfo != NULL);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(solvelpagain != NULL);
   assert(propagateagain != NULL);

   *cutoff = FALSE;
   *infeasible = FALSE;
   *solvelpagain = FALSE;
   *propagateagain = FALSE;

   /**@todo avoid checking the same pseudosolution twice */

   /* enforce constraints by branching, applying additional cutting planes (if LP is being processed),
    * introducing new constraints, or tighten the domains
    */
   debugMessage("enforcing constraints on %s solution\n", SCIPtreeHasActiveNodeLP(tree) ? "LP" : "pseudo");

   /* check, if the solution is infeasible anyway due to it's objective value */
   if( SCIPtreeHasActiveNodeLP(tree) )
      objinfeasible = FALSE;
   else
   {
      pseudoobjval = SCIPlpGetPseudoObjval(lp, set);
      objinfeasible = SCIPsetIsLT(set, pseudoobjval, SCIPnodeGetLowerbound(SCIPtreeGetActiveNode(tree)));
   }

   /* enforce constraints until a handler resolved an infeasibility with cutting off the node, branching, 
    * reducing a domain, or separating a cut
    * if a constraint handler introduced new constraints to enforce his constraints, the newly added constraints
    * have to be enforced themselves
    */
   resolved = FALSE;
   for( h = 0; h < set->nconshdlrs && !resolved; ++h )
   {
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      if( SCIPtreeHasActiveNodeLP(tree) )
      {
         assert(lp->solved);
         assert(lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL);
         CHECK_OKAY( SCIPconshdlrEnforceLPSol(conshdlrs_enfo[h], memhdr, set, stat, prob, tree, sepastore, &result) );
      }
      else
      {
         CHECK_OKAY( SCIPconshdlrEnforcePseudoSol(conshdlrs_enfo[h], memhdr, set, stat, prob, tree, objinfeasible, 
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
      case SCIP_DIDNOTRUN:
         assert(tree->nchildren == 0);
         assert(!SCIPtreeHasActiveNodeLP(tree) || lp->solved);
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         assert(objinfeasible);
         *infeasible = TRUE;
         break;

      case SCIP_FEASIBLE:
         assert(tree->nchildren == 0);
         assert(!SCIPtreeHasActiveNodeLP(tree) || lp->solved);
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         break;

      case SCIP_INFEASIBLE:
         assert(tree->nchildren == 0);
         assert(!SCIPtreeHasActiveNodeLP(tree) || lp->solved);
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         *infeasible = TRUE;
         break;

      case SCIP_CUTOFF:
         assert(tree->nchildren == 0);
         assert(!SCIPtreeHasActiveNodeLP(tree) || lp->solved);

         /* the found cuts are of no use, because the node is infeasible anyway */
         CHECK_OKAY( SCIPsepastoreClearCuts(sepastore, memhdr, set, lp) );

         *cutoff = TRUE;
         *infeasible = TRUE;
         resolved = TRUE;
         break;

      case SCIP_SEPARATED:
         assert(tree->nchildren == 0);
         assert(SCIPsepastoreGetNCuts(sepastore) > 0);

         /* apply found cuts */
         olddomchgcount = stat->domchgcount;
         CHECK_OKAY( SCIPsepastoreApplyCuts(sepastore, memhdr, set, stat, tree, lp, branchcand, eventqueue, cutoff) );

         *infeasible = TRUE;
         *solvelpagain = !(*cutoff);
         *propagateagain = !(*cutoff) && (stat->domchgcount != olddomchgcount);
         resolved = TRUE;
         break;

      case SCIP_REDUCEDDOM:
         assert(tree->nchildren == 0);
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         *infeasible = TRUE;
         *solvelpagain = TRUE;
         *propagateagain = TRUE;
         resolved = TRUE;
         break;

      case SCIP_CONSADDED:
         assert(tree->nchildren == 0);
         assert(!SCIPtreeHasActiveNodeLP(tree) || lp->solved);
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         *infeasible = TRUE;
         *solvelpagain = TRUE; /* the separation for new constraints should be called */
         *propagateagain = TRUE;
         resolved = TRUE;
         break;

      case SCIP_BRANCHED:
         assert(tree->nchildren >= 1);
         assert(!SCIPtreeHasActiveNodeLP(tree) || lp->solved);
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         *infeasible = TRUE;
         resolved = TRUE;
         break;

      case SCIP_SOLVELP:
         assert(!SCIPtreeHasActiveNodeLP(tree));
         assert(tree->nchildren == 0);
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         *infeasible = TRUE;
         *solvelpagain = TRUE;
         resolved = TRUE;
         SCIPtreeSetActiveNodeLP(tree, TRUE); /* the node's LP must be solved */
         break;

      default:
         errorMessage("invalid result code <%d> from enforcing method of constraint handler <%s>\n",
            result, SCIPconshdlrGetName(conshdlrs_enfo[h]));
         return SCIP_INVALIDRESULT;
      }  /*lint !e788*/

      /* the enforcement method may add a primal solution, after which the LP status could be set to
       * objective limit reached
       */
      if( SCIPtreeHasActiveNodeLP(tree) && lp->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT )
      {
         *cutoff = TRUE;
         *infeasible = TRUE;
         *solvelpagain = FALSE;
         *propagateagain = FALSE;
         resolved = TRUE;
      }

      assert(!(*cutoff) || (resolved && *infeasible && !(*solvelpagain) && !(*propagateagain)));
      assert(!(*solvelpagain) || (resolved && *infeasible));
      assert(!(*propagateagain) || (resolved && *infeasible));
   }
   assert(!objinfeasible || *infeasible);
   debugMessage(" -> enforcing result: infeasible=%d, solvelpagain=%d, propagateagain=%d, resolved=%d\n",
      *infeasible, *solvelpagain, *propagateagain, resolved);
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);

   if( *infeasible && !resolved )
   {
      /* the node is infeasible, but no constraint handler could resolve the infeasibility
       * -> branch on the pseudo solution
       * -> e.g. select non-fixed binary or integer variable x with value x', create three
       *    sons: x <= x'-1, x = x', and x >= x'+1.
       *    In the left and right branch, the current solution is cut off. In the middle
       *    branch, the constraints can hopefully reduce domains of other variables to cut
       *    off the current solution.
       */

      assert(!(*cutoff));
      assert(!(*solvelpagain));
   
      /* branch on pseudo solution */
      CHECK_OKAY( SCIPbranchPseudo(set->scip, &result) );

      switch( result )
      {  
      case SCIP_CUTOFF:
         assert(tree->nchildren == 0);
         *cutoff = TRUE;
         resolved = TRUE;
         break;
      case SCIP_BRANCHED:
         assert(tree->nchildren >= 1);
         resolved = TRUE;
         break;
      case SCIP_REDUCEDDOM:
         assert(tree->nchildren == 0);
         *solvelpagain = TRUE;
         *propagateagain = TRUE;
         resolved = TRUE;
         break;
      case SCIP_DIDNOTRUN:
         /* all integer variables in the infeasible pseudo solution are fixed,
          * - if no continuous variables exist and all variables are known, the infeasible pseudo solution is completely
          *   fixed, and the node is infeasible
          * - if at least one continuous variable exist or we do not know all variables due to external pricers, we cannot
          *   resolve the infeasibility by branching -> solve LP (and maybe price in additional variables)
          */
         assert(tree->nchildren == 0);
         assert(SCIPbranchcandGetNPseudoCands(branchcand) == 0);

         if( prob->ncontvars == 0 && set->nactivepricers == 0 )
         {
            *cutoff = TRUE;
            resolved = TRUE;
         }
         else
         {
            assert(!SCIPtreeHasActiveNodeLP(tree));
            
            /* solve the LP in the next loop */
            SCIPtreeSetActiveNodeLP(tree, TRUE);
            *solvelpagain = TRUE;
         }
         break;

      default:
         errorMessage("invalid result code <%d> from SCIPbranchPseudo()\n", result);
         abort();
      }  /*lint !e788*/
   }
   assert(*infeasible || !resolved);

   return SCIP_OKAY;
}

/** solves the active node */
static
RETCODE solveNode(
   MEMHDR*          memhdr,             /**< block memory buffers */
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
   Bool*            infeasible,         /**< pointer to store whether the active node's solution is infeasible */
   Bool*            restart             /**< should solving process be started again with presolving? */
   )
{
   NODE* actnode;
   Bool actnodehaslp;
   Bool initiallpsolved;
   Bool solvelpagain;
   Bool propagateagain;
   int actdepth;
   int nlperrors;

   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);
   assert(primal != NULL);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(restart != NULL);

   *cutoff = FALSE;
   *infeasible = FALSE;
   *restart = FALSE;

   actnode = SCIPtreeGetActiveNode(tree);
   assert(actnode != NULL);
   assert(SCIPnodeGetType(actnode) == SCIP_NODETYPE_ACTNODE);
   actdepth = SCIPnodeGetDepth(actnode);

   debugMessage("Processing node %lld in depth %d, %d siblings\n", stat->nnodes, actdepth, tree->nsiblings);
   debugMessage("current pseudosolution: obj=%g", SCIPlpGetPseudoObjval(lp, set));
   debug(SCIPprobPrintPseudoSol(prob, set));
   
   /* check, if we want to solve the LP at the selected node:
    * - solve the LP, if the lp solve depth and frequency demand solving
    * - solve the root LP, if the LP solve frequency is set to 0
    * - solve the root LP, if there are continuous variables present
    * - don't solve the node if its cut off by the pseudo objective value anyway
    */
   actnodehaslp = (set->lpsolvedepth == -1 || actdepth <= set->lpsolvedepth);
   actnodehaslp = actnodehaslp && (set->lpsolvefreq >= 1 && actdepth % set->lpsolvefreq == 0);
   actnodehaslp = actnodehaslp || (actdepth == 0 && set->lpsolvefreq == 0);
   actnodehaslp = actnodehaslp || (actdepth == 0 && prob->ncontvars > 0);
   actnodehaslp = actnodehaslp && SCIPsetIsLT(set, SCIPlpGetPseudoObjval(lp, set), primal->cutoffbound);
   SCIPtreeSetActiveNodeLP(tree, actnodehaslp);

   /* external node solving loop */
   initiallpsolved = FALSE;
   nlperrors = 0;
   stat->npricerounds = 0;
   stat->nseparounds = 0;
   solvelpagain = TRUE;
   propagateagain = TRUE;
   while( (solvelpagain || propagateagain) && nlperrors < MAXNLPERRORS && !(*restart) )
   {
      Real pseudoobjval;
      Bool lperror;

      *infeasible = FALSE;
      lperror = FALSE;

      /* domain propagation */
      if( propagateagain )
      {
         propagateagain = FALSE;
         CHECK_OKAY( SCIPpropagateDomains(memhdr, set, stat, prob, tree, cutoff) );
         if( *cutoff )
         {
            debugMessage("domain propagation determined that node is infeasible\n");
            SCIPnodeUpdateLowerbound(actnode, stat, primal->cutoffbound);
            *infeasible = TRUE;
            break;
         }
      }

      /* check, if we want to solve the LP at this node */
      if( solvelpagain && SCIPtreeHasActiveNodeLP(tree) )
      {
         solvelpagain = FALSE;

         /* solve the node's LP */
         CHECK_OKAY( solveNodeLP(memhdr, set, stat, prob, primal, tree, lp, pricestore, sepastore,
               cutpool, branchcand, conflict, eventfilter, eventqueue, conshdlrs_sepa,
               initiallpsolved, cutoff, &lperror) );
         initiallpsolved = TRUE;
         debugMessage(" -> LP status: %d, LP obj: %g\n", SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set));

         /* if an error occured during LP solving, switch to pseudo solution */
         if( lperror )
         {
            SCIPtreeSetActiveNodeLP(tree, FALSE);
            nlperrors++;
            infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
               "(node %lld) unresolved numerical troubles in LP %d -- using pseudo solution instead (loop %d)\n", 
               stat->nnodes, stat->nlps, nlperrors);
         }
               
         /* if we solve exactly, the LP claims to be infeasible but the infeasibility could not be proved,
          * we have to forget about the LP and use the pseudo solution instead
          */
         if( set->exactsolve && lp->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE
            && SCIPnodeGetLowerbound(actnode) < primal->cutoffbound )
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
               SCIPtreeSetActiveNodeLP(tree, FALSE);
               infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) could not prove infeasibility of LP %d -- using pseudo solution (%d unfixed vars) instead\n", 
                  stat->nnodes, stat->nlps, SCIPbranchcandGetNPseudoCands(branchcand));
            }
         }
      }
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      /* update lower bound w.r.t. the pseudo solution */
      pseudoobjval = SCIPlpGetPseudoObjval(lp, set);
      SCIPnodeUpdateLowerbound(actnode, stat, pseudoobjval);
      debugMessage(" -> new lower bound: %g [%g] (pseudoobj: %g [%g])\n",
         SCIPnodeGetLowerbound(actnode), SCIPprobExternObjval(prob, set, SCIPnodeGetLowerbound(actnode)),
         pseudoobjval, SCIPprobExternObjval(prob, set, pseudoobjval));

      /* call pseudo conflict analysis, if the node is cut off due to the pseudo objective value */
      if( pseudoobjval >= primal->cutoffbound )
      {
         CHECK_OKAY( SCIPconflictAnalyzePseudo(conflict, memhdr, set, stat, prob, tree, lp, NULL) );
      }
      
      /* check for infeasible node by bounding */
      if( *cutoff
         || SCIPnodeGetLowerbound(actnode) >= primal->cutoffbound
         || (!set->exactsolve && SCIPsetIsGE(set, SCIPnodeGetLowerbound(actnode), primal->cutoffbound)) )
      {
         debugMessage("node is infeasible (lower=%g, upper=%g)\n", SCIPnodeGetLowerbound(actnode), primal->cutoffbound);
         SCIPnodeUpdateLowerbound(actnode, stat, primal->cutoffbound);
         *infeasible = TRUE;
         break;
      }
      assert(!SCIPtreeHasActiveNodeLP(tree) || lp->solved);

      /* enforce constraints */
      CHECK_OKAY( enforceConstraints(memhdr, set, stat, prob, primal, tree, lp, sepastore, branchcand, eventqueue,
            conshdlrs_enfo, cutoff, infeasible, &solvelpagain, &propagateagain) );
      assert(!solvelpagain || !(*cutoff));

      /* check for immediate restart */
      if( actdepth == 0 && !(*cutoff) && set->restartbdchgs > 0 && stat->nrootboundchgsrun >= set->restartbdchgs )
         *restart = TRUE;
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);

   /* check for too many LP errors */
   if( nlperrors >= MAXNLPERRORS )
   {
      errorMessage("(node %lld) unresolved numerical troubles in LP %d -- aborting\n", stat->nnodes, stat->nlps);
      return SCIP_LPERROR;
   }

   /* check for final restart */
   if( actdepth == 0 && !(*cutoff) && set->restartbdchgs >= 0 && stat->nrootboundchgsrun > 0 )
      *restart = TRUE;

   /* remember root LP solution */
   if( actdepth == 0 && !(*restart) )
   {
      SCIPprobStoreRootSol(prob, SCIPtreeHasActiveNodeLP(tree));
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
   depth = SCIPtreeGetActiveDepth(tree);
   lpforkdepth = tree->actlpfork != NULL ? SCIPnodeGetDepth(tree->actlpfork) : -1;

   /* call heuristics */
   ndelayedheurs = 0;
   for( h = 0; h < set->nheurs; ++h )
   {
      CHECK_OKAY( SCIPheurExec(set->heurs[h], set, primal, depth, lpforkdepth, SCIPtreeHasActiveNodeLP(tree), plunging,
            &ndelayedheurs, &result) );
      *foundsol = *foundsol || (result == SCIP_FOUNDSOL);
   }
   assert(0 <= ndelayedheurs && ndelayedheurs <= set->nheurs);

   return SCIP_OKAY;
}

/** if feasible, adds current solution to the solution storage */
static
RETCODE addCurrentSolution(
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTFILTER*     eventfilter         /**< event filter for global (not variable dependent) events */
   )
{
   SOL* sol;
   Bool foundsol;

   /* found a feasible solution */
   if( SCIPtreeHasActiveNodeLP(tree) )
   {
      /* start clock for LP solutions */
      SCIPclockStart(stat->lpsoltime, set);
      
      /* add solution to storage */
      CHECK_OKAY( SCIPsolCreateLPSol(&sol, memhdr, set, stat, primal, tree, lp, NULL) );
      if( set->exactsolve )
      {
         /* if we want to solve exactly, we have to check the solution exactly again */
         CHECK_OKAY( SCIPprimalTrySolFree(primal, memhdr, set, stat, prob, tree, lp, eventfilter, &sol,
               TRUE, TRUE, &foundsol) );
      }
      else
      {
         CHECK_OKAY( SCIPprimalAddSolFree(primal, memhdr, set, stat, prob, tree, lp, eventfilter, &sol, &foundsol) );
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
      CHECK_OKAY( SCIPsolCreatePseudoSol(&sol, memhdr, set, stat, primal, tree, lp, NULL) );
      if( set->exactsolve )
      {
         /* if we want to solve exactly, we have to check the solution exactly again */
         CHECK_OKAY( SCIPprimalTrySolFree(primal, memhdr, set, stat, prob, tree, lp, eventfilter, &sol,
               TRUE, TRUE, &foundsol) );
      }
      else
      {
         CHECK_OKAY( SCIPprimalAddSolFree(primal, memhdr, set, stat, prob, tree, lp, eventfilter, &sol, &foundsol) );
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
   MEMHDR*          memhdr,             /**< block memory buffers */
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
   NODE* actnode;
   NODE* nextnode;
   EVENT event;
   int nnodes;
   int depth;
   Bool cutoff;
   Bool infeasible;
   Bool foundsol;

   assert(set != NULL);
   assert(memhdr != NULL);
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

   *restart = FALSE;

   /* sort constraint handlers by priorities */
   ALLOC_OKAY( duplicateMemoryArray(&conshdlrs_sepa, set->conshdlrs, set->nconshdlrs) );
   ALLOC_OKAY( duplicateMemoryArray(&conshdlrs_enfo, set->conshdlrs, set->nconshdlrs) );
   SCIPbsortPtr((void**)conshdlrs_sepa, set->nconshdlrs, SCIPconshdlrCompSepa);
   SCIPbsortPtr((void**)conshdlrs_enfo, set->nconshdlrs, SCIPconshdlrCompEnfo);

   nextnode = NULL;

   while( !SCIPsolveIsStopped(set, stat) && !(*restart) )
   {
      assert(SCIPbufferGetNUsed(set->buffer) == 0);

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
      actnode = nextnode;
      nextnode = NULL;
      assert(SCIPbufferGetNUsed(set->buffer) == 0);

      /* update plunging and backtracking statistics */
      if( actnode != NULL )
      {
         if( SCIPnodeGetType(actnode) == SCIP_NODETYPE_LEAF )
         {
            stat->plungedepth = 0;
            stat->nbacktracks++;
            debugMessage("selected leaf node, lowerbound=%g, plungedepth=%d\n", actnode->lowerbound, stat->plungedepth);
         }
         else if( SCIPnodeGetType(actnode) == SCIP_NODETYPE_CHILD )
         {
            stat->plungedepth++;
            debugMessage("selected child node, lowerbound=%g, plungedepth=%d\n", actnode->lowerbound, stat->plungedepth);
         }
         else
         {
            assert(SCIPnodeGetType(actnode) == SCIP_NODETYPE_SIBLING);
            debugMessage("selected sibling node, lowerbound=%g, plungedepth=%d\n", actnode->lowerbound, stat->plungedepth);
         }
      }

      /* start node activation timer */
      SCIPclockStart(stat->nodeactivationtime, set);

      /* delay events in node activation */
      CHECK_OKAY( SCIPeventqueueDelay(eventqueue) );

      /* activate selected node */
      CHECK_OKAY( SCIPnodeActivate(actnode, memhdr, set, stat, tree, lp, branchcand, eventqueue, primal->cutoffbound) );
      assert(SCIPtreeGetCurrentNode(tree) == actnode);
      assert(SCIPtreeGetActiveNode(tree) == actnode);

      /* process the delayed events */
      CHECK_OKAY( SCIPeventqueueProcess(eventqueue, memhdr, set, primal, lp, branchcand, eventfilter) );

      /* stop node activation timer */
      SCIPclockStop(stat->nodeactivationtime, set);
      assert(SCIPbufferGetNUsed(set->buffer) == 0);

      /* if no more node was selected, we finished optimization */
      if( actnode == NULL )
         break;

      /* if node was cut off due to SCIPnodeCutoff() call, continue with the next node */
      if( actnode->cutoff )
         continue;

      /* update maxdepth and node count statistics */
      assert(SCIPnodeGetType(actnode) == SCIP_NODETYPE_ACTNODE);
      depth = SCIPnodeGetDepth(actnode);
      stat->maxdepth = MAX(stat->maxdepth, depth);
      stat->maxtotaldepth = MAX(stat->maxtotaldepth, depth);
      stat->nnodes++;
      stat->ntotalnodes++;
      
      /* issue NODEACTIVATED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEACTIVATED) );
      CHECK_OKAY( SCIPeventChgNode(&event, actnode) );
      CHECK_OKAY( SCIPeventProcess(&event, memhdr, set, NULL, NULL, NULL, eventfilter) );

      /* solve active node */
      CHECK_OKAY( solveNode(memhdr, set, stat, prob, primal, tree, lp, pricestore, sepastore, branchcand, cutpool,
            conflict, eventfilter, eventqueue, conshdlrs_sepa, conshdlrs_enfo, &cutoff, &infeasible, restart) );
      assert(!cutoff || infeasible);
      assert(SCIPbufferGetNUsed(set->buffer) == 0);
      assert(SCIPtreeGetCurrentNode(tree) == actnode);
      assert(SCIPtreeGetActiveNode(tree) == actnode);

      /* check for restart */
      if( !(*restart) )
      {
         /* change color of node in VBC output */
         SCIPvbcSolvedNode(stat->vbc, stat, actnode);

         /* check, if the current solution is feasible */
         if( !infeasible )
         {
            assert(!SCIPtreeHasActiveNodeLP(tree) || lp->solved);
            assert(!cutoff);

            /* node solution is feasible: add it to the solution store */
            CHECK_OKAY( addCurrentSolution(memhdr, set, stat, prob, primal, tree, lp, branchcand, eventfilter) );

            /* issue NODEFEASIBLE event */
            CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEFEASIBLE) );
            CHECK_OKAY( SCIPeventChgNode(&event, actnode) );
            CHECK_OKAY( SCIPeventProcess(&event, memhdr, set, NULL, NULL, NULL, eventfilter) );
         }
         else
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
            CHECK_OKAY( SCIPeventChgNode(&event, actnode) );
            CHECK_OKAY( SCIPeventProcess(&event, memhdr, set, NULL, NULL, NULL, eventfilter) );
         }
         assert(SCIPbufferGetNUsed(set->buffer) == 0);

         /* if no branching was created, the node was not cut off, but it's lower bound is still smaller than
          * the cutoff bound, we have to branch on a non-fixed variable;
          * this can happen, if we want to solve exactly, the current solution was declared feasible by the
          * constraint enforcement, but in exact solution checking it was found out to be infeasible;
          * in this case, no branching would have been generated by the enforcement of constraints, but we
          * have to further investigate the current sub tree
          */
         if( !cutoff && tree->nchildren == 0 && SCIPnodeGetLowerbound(actnode) < primal->cutoffbound )
         {
            RESULT result;

            assert(set->exactsolve);

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
                  CHECK_OKAY( SCIPbranchExecPseudo(memhdr, set, stat, tree, lp, branchcand, eventqueue, &result) );
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
      CHECK_OKAY( SCIPdispPrintLine(set, stat, (SCIPnodeGetDepth(actnode) == 0) && infeasible && !foundsol) );

      debugMessage("Processing of node in depth %d finished. %d siblings, %d children, %d leaves left\n", 
         SCIPnodeGetDepth(actnode), tree->nsiblings, tree->nchildren, SCIPtreeGetNLeaves(tree));
      debugMessage("**********************************************************************\n");
   }
   assert(SCIPbufferGetNUsed(set->buffer) == 0);

   debugMessage("Problem solving finished\n");

   /* free sorted constraint handler arrays */
   freeMemoryArrayNull(&conshdlrs_sepa);
   freeMemoryArrayNull(&conshdlrs_enfo);

   return SCIP_OKAY;
}

