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
#pragma ident "@(#) $Id: solve.c,v 1.94 2004/03/22 17:46:42 bzfpfend Exp $"

/**@file   solve.c
 * @brief  main solving loop and node processing
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "def.h"
#include "set.h"
#include "stat.h"
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
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   )
{
   return (SCIPinterrupted()
      || (set->nodelimit >= 0 && stat->nnodes >= set->nodelimit)
      || SCIPclockGetTime(stat->solvingtime) >= set->timelimit
      || (SCIPgetMemUsed(set->scip) >= set->memlimit*1024.0*1024.0)
      || (SCIPstage(set->scip) >= SCIP_STAGE_SOLVING && SCIPsetIsLT(set, SCIPgetGap(set->scip), set->gaplimit))
      || (set->sollimit >= 0 && SCIPstage(set->scip) >= SCIP_STAGE_PRESOLVED &&
         SCIPgetNSolsFound(set->scip) >= set->sollimit));
}

/** outputs the reason for termination */
void SCIPsolvePrintStopReason(
   const SET*       set,                /**< global SCIP settings */
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
}

/** domain propagation */
static
RETCODE propagateDomains(
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   RESULT result;
   Bool propagain;
   int h;

   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(cutoff != NULL);

   debugMessage("domain propagation\n");
   *cutoff = FALSE;
   do
   {
      propagain = FALSE;
      for( h = 0; h < set->nconshdlrs && !(*cutoff); ++h )
      {
         CHECK_OKAY( SCIPconshdlrPropagate(set->conshdlrs[h], memhdr, set, stat, prob, SCIPnodeGetDepth(tree->actnode),
                        &result) );
         propagain = propagain || (result == SCIP_REDUCEDDOM);
         *cutoff = *cutoff || (result == SCIP_CUTOFF);
      }
   }
   while( propagain && !(*cutoff) );

   return SCIP_OKAY;
}

/** returns whether the given variable with the old LP solution value should lead to an update of the history entry */
static
Bool isHistoryUpdateValid(
   VAR*             var,                /**< problem variable */
   const SET*       set,                /**< global SCIP settings */
   Real             oldlpsolval         /**< solution value of variable in old LP */
   )
{
   Real newlpsolval;

   assert(var != NULL);

   /* if the old LP solution value is unknown, the history update cannot be performed */
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

/** history flag stored in the variables to mark them for the history update */
enum HistoryFlag
{
   HISTORY_NONE     = 0,                /**< variable's bounds were not changed */
   HISTORY_IGNORE   = 1,                /**< bound changes on variable should be ignored for history updates */
   HISTORY_UPDATE   = 2                 /**< history value of variable should be updated */
};

/** updates the variable's history values after the node's initial LP was solved */
static
RETCODE updateHistory(
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< LP data */
   )
{
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnode->active);
   assert(tree->path != NULL);
   assert(tree->path[tree->actnode->depth] == tree->actnode);
   assert(lp != NULL);

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

      assert(tree->actlpfork->active);
      assert(tree->path[tree->actlpfork->depth] == tree->actlpfork);

      /* get a buffer for the collected bound changes; start with a size twice as large as the number of nodes between
       * current node and LP fork
       */
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &updates, tree->actnode->depth - tree->actlpfork->depth) );
      nupdates = 0;
      nvalidupdates = 0;

      /* search the nodes from LP fork down to current node for bound changes in between; move in this direction, 
       * because the bound changes closer to the LP fork are more likely to have a valid LP solution information
       * attached; collect the bound changes for history value updates and mark the corresponding variables such
       * that they are not updated twice in case of more than one bound change on the same variable
       */
      for( d = tree->actlpfork->depth+1; d <= tree->actnode->depth; ++d )
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
                  if( var->historyflag == HISTORY_NONE )
                  {
                     /* remember the bound change and mark the variable */
                     CHECK_OKAY( SCIPsetReallocBufferArray(set, &updates, nupdates+1) );
                     updates[nupdates] = &boundchgs[i];
                     nupdates++;

                     /* check, if the bound change would lead to a valid history update */
                     if( isHistoryUpdateValid(var, set, boundchgs[i].data.branchingdata.lpsolval) )
                     {
                        var->historyflag = HISTORY_UPDATE;
                        nvalidupdates++;
                     }
                     else
                        var->historyflag = HISTORY_IGNORE;
                  }
               }
            }
         }
      }

      /* update the history values and reset the variables' flags; assume, that the responsibility for the dual gain
       * is equally spread on all bound changes that lead to valid history updates
       */
      weight = nvalidupdates > 0 ? 1.0 / (Real)nvalidupdates : 1.0;
      lpgain = tree->actnode->lowerbound - tree->actlpfork->lowerbound;
      for( i = 0; i < nupdates; ++i )
      {
         assert(updates[i]->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING);
         var = updates[i]->var;
         assert(var->historyflag != HISTORY_NONE);
         if( var->historyflag == HISTORY_UPDATE )
         {
            SCIPvarUpdateLPHistory(var, set, stat, 
               SCIPvarGetLPSol(var) - updates[i]->data.branchingdata.lpsolval, lpgain, weight);
         }
         var->historyflag = HISTORY_NONE;
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
   const SET*       set,                /**< global SCIP settings */
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

   SCIPnodeUpdateLowerbound(node, lpobjval);

   return SCIP_OKAY;
}

/** constructs the LP of the root node */
static
RETCODE initRootLP(
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< LP data */
   PRICESTORE*      pricestore,         /**< pricing storage */
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   VAR* var;
   COL* col;
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
      CHECK_OKAY( SCIPconshdlrInitLP(set->conshdlrs[h], memhdr, set, prob) );
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
   const SET*       set,                /**< global SCIP settings */
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
   EVENT event;

   assert(stat != NULL);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(lp != NULL);
   assert(cutoff != NULL);
   assert(lperror != NULL);

   *cutoff = FALSE;
   *lperror = FALSE;

   /* load the LP into the solver and load the LP state */
   debugMessage("loading LP\n");
   CHECK_OKAY( SCIPtreeLoadLP(tree, memhdr, set, stat, lp) );
   
   /* init root node LP */
   if( SCIPnodeGetDepth(tree->actnode) == 0 )
   {
      assert(stat->lpcount == 0);
      CHECK_OKAY( initRootLP(memhdr, set, stat, prob, tree, lp, pricestore, sepastore, branchcand, eventqueue, cutoff) );
   }

   /* only the bounds and the constraint list may have been changed since the LP state node
    * -> dual simplex is applicable as first solver
    */
   debugMessage("node: solve initial LP\n");
   CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, TRUE, lperror) );
   assert(lp->solved || *lperror);

   if( !(*lperror) )
   {
      /* issue FIRSTLPSOLVED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_FIRSTLPSOLVED) );
      CHECK_OKAY( SCIPeventChgNode(&event, tree->actnode) );
      CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, eventfilter) );
      
      /* update history values */
      CHECK_OKAY( updateHistory(set, stat, tree, lp) );
   }

   return SCIP_OKAY;
}

/** solve the current LP of a node with price and cut */
static
RETCODE solveNodeLP(
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICESTORE*      pricestore,         /**< pricing storage */
   SEPASTORE*       sepastore,          /**< separation storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   CONFLICT*        conflict,           /**< conflict analysis data */
   LPCONFLICT*      lpconflict,         /**< conflict analysis data for infeasible LP conflicts */
   PRIMAL*          primal,             /**< primal data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   CONSHDLR**       conshdlrs_sepa,     /**< constraint handlers sorted by separation priority */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   Bool*            lperror             /**< pointer to store whether an unresolved error in LP solving occured */
   )
{
   RESULT result;
   EVENT event;
   Bool mustprice;
   Bool mustsepar;
   Bool root;
   int h;
   int s;

   assert(set != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(lp != NULL);
   assert(pricestore != NULL);
   assert(sepastore != NULL);
   assert(cutpool != NULL);
   assert(primal != NULL);
   assert(set->nconshdlrs == 0 || conshdlrs_sepa != NULL);
   assert(cutoff != NULL);
   assert(lperror != NULL);

   root = (SCIPnodeGetDepth(tree->actnode) == 0);

   /* solve initial LP of price-and-cut loop */
   debugMessage("node: solve LP with price and cut\n");
   CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, TRUE, lperror) );
   assert(lp->solved);

   /* display node information line for root node */
   if( SCIPnodeGetDepth(tree->actnode) == 0 && (VERBLEVEL)set->verblevel >= SCIP_VERBLEVEL_HIGH )
   {
      CHECK_OKAY( SCIPdispPrintLine(set, stat, TRUE) );
   }

   /* price-and-cut loop */
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
         mustprice = !lp->solved;
         mustsepar = mustsepar || !lp->solved;
         
         /* after adding columns, the LP should be primal feasible such that primal simplex is applicable;
          * if LP was infeasible, we have to use dual simplex
          */
         debugMessage("pricing: solve LP\n");
         CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, TRUE, lperror) );
         assert(lp->solved || *lperror);

         /* reset bounds temporarily set by pricer to their original values */
         debugMessage("pricing: reset bounds\n");
         CHECK_OKAY( SCIPpricestoreResetBounds(pricestore, memhdr, set, stat, lp, branchcand, eventqueue) );
         assert(SCIPpricestoreGetNVars(pricestore) == 0);
         assert(SCIPpricestoreGetNBoundResets(pricestore) == 0);
         mustprice = mustprice || !lp->solved;
         mustsepar = mustsepar || !lp->solved;

         /* solve LP again after resetting bounds (with dual simplex) */
         debugMessage("pricing: solve LP after resetting bounds\n");
         CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, FALSE, lperror) );
         assert(lp->solved || *lperror);

         /* increase pricing round counter */
         stat->npricerounds++;

         /* display node information line for root node */
         if( SCIPnodeGetDepth(tree->actnode) == 0 && mustprice && (VERBLEVEL)set->verblevel >= SCIP_VERBLEVEL_FULL )
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
         CHECK_OKAY( nodeUpdateLowerboundLP(tree->actnode, set, lp) );
         debugMessage(" -> new lower bound: %g (LP status: %d, LP obj: %g)\n",
            tree->actnode->lowerbound, SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set));
      }
      else
      {
         debugMessage(" -> error solving LP. keeping old bound: %g\n", tree->actnode->lowerbound);
      }

      /* if the LP is infeasible or exceeded the objective limit, we don't need to separate cuts */
      mustsepar = mustsepar
         && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_INFEASIBLE)
         && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OBJLIMIT)
         && SCIPsetIsLT(set, tree->actnode->lowerbound, primal->cutoffbound);

      /* separation (has not to be done completely, because we just want to increase the lower bound) */
      if( !(*cutoff) && !(*lperror) && mustsepar )
      {
         Bool separateagain;
         Bool enoughcuts;

         assert(lp->solved);

         enoughcuts = FALSE;

         /* global cut pool separation */
         debugMessage("global cut pool separation\n");
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         CHECK_OKAY( SCIPcutpoolSeparate(cutpool, memhdr, set, stat, lp, sepastore, root, &result) );
         *cutoff = *cutoff || (result == SCIP_CUTOFF);
         enoughcuts = enoughcuts || (SCIPsepastoreGetNCuts(sepastore) >= SCIPsetGetMaxsepacuts(set, root)/2);

         /* constraint separation */
         debugMessage("constraint separation\n");

         /* reset the constraint handler's separation calls */
         for( h = 0; h < set->nconshdlrs; ++h )
            SCIPconshdlrResetSepa(set->conshdlrs[h]);
         
         /* separate constraints */
         separateagain = TRUE;
         while( !(*cutoff) && separateagain )
         {
            /* try separating constraints of the constraint handlers until the cut pool is at least half full */
            while( !(*cutoff) && separateagain && !enoughcuts )
            {
               separateagain = FALSE;
               for( h = 0; h < set->nconshdlrs && !(*cutoff) && !enoughcuts; ++h )
               {
                  CHECK_OKAY( SCIPconshdlrSeparate(conshdlrs_sepa[h], memhdr, set, stat, prob, sepastore, 
                                 SCIPnodeGetDepth(tree->actnode), &result) );
                  *cutoff = *cutoff || (result == SCIP_CUTOFF);
                  separateagain = separateagain || (result == SCIP_CONSADDED);
                  enoughcuts = enoughcuts || (SCIPsepastoreGetNCuts(sepastore) >= SCIPsetGetMaxsepacuts(set, root)/2);
               }
            }
            
            /* separate LP, if no cuts have been found by the constraint handlers */
            if( SCIPsepastoreGetNCuts(sepastore) == 0 )
            {
               /* sort separators by priority */
               SCIPsetSortSepas(set);

               /* call separators */
               for( s = 0; s < set->nsepas && !(*cutoff) && !separateagain && !enoughcuts; ++s )
               {
                  CHECK_OKAY( SCIPsepaExec(set->sepas[s], set, stat, sepastore, SCIPnodeGetDepth(tree->actnode), &result) );
                  *cutoff = *cutoff || (result == SCIP_CUTOFF);
                  separateagain = separateagain || (result == SCIP_CONSADDED);
                  enoughcuts = enoughcuts || (SCIPsepastoreGetNCuts(sepastore) >= SCIPsetGetMaxsepacuts(set, root)/2);
               }
            }
         }

         if( *cutoff )
         {
            /* the found cuts are of no use, because the node is infeasible anyway */
            CHECK_OKAY( SCIPsepastoreClearCuts(sepastore, memhdr, set, lp) );
         }
         else
         {
            Longint oldnboundchanges;

            oldnboundchanges = stat->nboundchanges;

            /* apply found cuts */
            CHECK_OKAY( SCIPsepastoreApplyCuts(sepastore, memhdr, set, stat, tree, lp, branchcand, eventqueue, cutoff) );

            if( !(*cutoff) )
            {
               /* if at least one of the cut was a bound change, propagate domains again */
               if( stat->nboundchanges != oldnboundchanges )
               {
                  CHECK_OKAY( propagateDomains(memhdr, set, stat, prob, tree, cutoff) );
               }
               
               mustprice = mustprice || !lp->solved;
               mustsepar = !lp->solved;
               
               if( !(*cutoff) )
               {
                  /* solve LP (with dual simplex) */
                  debugMessage("separation: solve LP\n");
                  CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, TRUE, lperror) );
               }
            }
         }
         assert(*cutoff || *lperror || lp->solved); /* after cutoff, the LP may be unsolved due to bound changes */

         /* increase separation round counter */
         stat->nseparounds++;

         /* display node information line for root node */
         if( SCIPnodeGetDepth(tree->actnode) == 0 )
         {
            CHECK_OKAY( SCIPdispPrintLine(set, stat, TRUE) );
         }
      }
   }

   /* update lower bound w.r.t. the the LP solution */
   if( *cutoff )
   {
      SCIPnodeUpdateLowerbound(tree->actnode, primal->cutoffbound);
   }
   else if( !(*lperror) )
   {
      assert(lp->solved);

      CHECK_OKAY( nodeUpdateLowerboundLP(tree->actnode, set, lp) );

      /* issue LPSOLVED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_LPSOLVED) );
      CHECK_OKAY( SCIPeventChgNode(&event, tree->actnode) );
      CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, eventfilter) );

      /* analyze an infeasible LP */
      if( !set->exactsolve
         && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT) )
      {
         CHECK_OKAY( SCIPlpconflictAnalyze(lpconflict, memhdr, set, stat, prob, tree, lp, conflict, NULL) );
      }
   }
   debugMessage(" -> final lower bound: %g (LP status: %d, LP obj: %g)\n",
      tree->actnode->lowerbound, SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set));

   return SCIP_OKAY;
}

/** strengthens variable's bounds looking at reduced costs */
static
RETCODE redcostStrengthening(
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   PRIMAL*          primal,             /**< primal data */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   COL** cols;
   VAR* var;
   int* cstat;
   Real lpobjval;
   Real redcost;
   Real oldlb;
   Real oldub;
   Real newbd;
   int ncols;
   int c;

   assert(tree != NULL);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(primal != NULL);
   assert(!SCIPsetIsInfinity(set, primal->cutoffbound));
   assert(SCIPsetIsLT(set, SCIPlpGetObjval(lp, set), primal->cutoffbound) );

   /* we cannot apply reduced cost fixing, if we want to solve exactly */
   /**@todo implement reduced cost fixing with interval arithmetics */
   if( set->exactsolve )
      return SCIP_OKAY;

   /* get LP columns */
   cols = SCIPlpGetCols(lp);
   ncols = SCIPlpGetNCols(lp);
   if( ncols == 0 )
      return SCIP_OKAY;

   stat->nredcoststrcalls++;

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
               SCIPvarAdjustUb(var, set, &newbd);

               /* check, if new bound is good enough */
               if( SCIPsetIsUbBetter(set, newbd, oldub) )
               {
                  /* strengthen upper bound */
                  debugMessage("redcost strengthening upper bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
                     SCIPvarGetName(var), oldlb, oldub, oldlb, newbd, primal->cutoffbound, lpobjval, redcost);
                  CHECK_OKAY( SCIPnodeAddBoundchg(tree->actnode, memhdr, set, stat, tree, lp, branchcand, eventqueue,
                                 var, newbd, SCIP_BOUNDTYPE_UPPER) );
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
               SCIPvarAdjustLb(var, set, &newbd);

               /* check, if new bound is good enough */
               if( SCIPsetIsLbBetter(set, newbd, oldlb) )
               {
                  /* strengthen lower bound */
                  debugMessage("redcost strengthening lower bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
                     SCIPvarGetName(var), oldlb, oldub, newbd, oldub, primal->cutoffbound, lpobjval, redcost);
                  CHECK_OKAY( SCIPnodeAddBoundchg(tree->actnode, memhdr, set, stat, tree, lp, branchcand, eventqueue,
                                 var, newbd, SCIP_BOUNDTYPE_LOWER) );
                  stat->nredcoststrfound++;
               }
            }
         }
         break;

      default:
         errorMessage("invalid basis state\n");
         return SCIP_INVALIDDATA;
      }
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &cstat);

   /**@todo check, if the resolve in next solveNodeLP() call is too expensive in contrast to resolving directly
    *       after reduced cost fixing
    */
#if 0 /* resolving is done in the next solveNodeLP() call */
   /* resolve the LP (the optimal solution should stay the same) */
   CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, FALSE) );
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(EPSZ(SCIPsetRelDiff(set, SCIPlpGetObjval(lp, set), lpobjval), 100.0*(set)->feastol));
   /*assert(SCIPsetIsFeasEQ(set, SCIPlpGetObjval(lp, set), lpobjval));*/
#endif

   /* stop redcost strengthening activation timer */
   SCIPclockStop(stat->redcoststrtime, set);

   return SCIP_OKAY;
}

/** solves the current node's LP in a price-and-cut loop */
static
RETCODE priceAndCutLoop(
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICESTORE*      pricestore,         /**< pricing storage */
   SEPASTORE*       sepastore,          /**< separation storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   CONFLICT*        conflict,           /**< conflict analysis data */
   LPCONFLICT*      lpconflict,         /**< conflict analysis data for infeasible LP conflicts */
   PRIMAL*          primal,             /**< primal data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
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
   assert(tree->actnodehaslp);
   assert(tree->actnode != NULL);
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
      assert(lp->solved || lperror);
      debugMessage("price-and-cut-loop: initial LP status: %d, LP obj: %g\n", 
         SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set));
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);

   if( !(*cutoff) && !(*lperror) )
   {
      /* solve-redcost-propagate loop */
      do
      {
         /* continue solving the LP with price and cut */
         CHECK_OKAY( solveNodeLP(memhdr, set, stat, prob, tree, lp, pricestore, sepastore, cutpool, 
                        conflict, lpconflict, primal, branchcand, eventfilter, eventqueue, conshdlrs_sepa,
                        cutoff, lperror) );
         assert(*cutoff || *lperror || lp->solved);
         debugMessage("price-and-cut-loop: LP status: %d, LP obj: %g\n", SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set));

         /* reduced cost bound strengthening */
         if( !(*cutoff) && !(*lperror)
            && ((set->redcostfreq == 0 && tree->actnode->depth == 0)
               || (set->redcostfreq >= 1 && (tree->actnode->depth % set->redcostfreq) == 0))
            && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL
            && !SCIPsetIsInfinity(set, primal->cutoffbound)
            && SCIPsetIsLT(set, SCIPlpGetObjval(lp, set), primal->cutoffbound) )
         {
            Longint oldnboundchanges;
                        
            oldnboundchanges = stat->nboundchanges;
                        
            /* apply reduced cost bound strengthening */
            CHECK_OKAY( redcostStrengthening(memhdr, set, stat, prob, tree, lp, branchcand, primal, eventqueue) );
                        
            /* apply further domain propagation, if addional bounds have been changed */
            if( stat->nboundchanges != oldnboundchanges )
            {
               CHECK_OKAY( propagateDomains(memhdr, set, stat, prob, tree, cutoff) );
               if( *cutoff )
               {
                  debugMessage("redcost domain propagation determined that node is infeasible\n");
                  SCIPnodeUpdateLowerbound(tree->actnode, primal->cutoffbound);
               }
            }
         }
      }
      while( !(*cutoff) && !(*lperror) && !lp->solved );
   }

   /* update node's LP iteration counter */
   stat->nnodelps += stat->nlps - nlps;
   stat->nnodelpiterations += stat->nlpiterations - nlpiterations;

   return SCIP_OKAY;
}

/** enforces constraints by branching, separation, or domain reduction */
static
RETCODE enforceConstraints(
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   PRIMAL*          primal,             /**< primal data */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   CONSHDLR**       conshdlrs_enfo,     /**< constraint handlers for enforcing constraints, sorted by priority */
   Bool*            cutoff,             /**< pointer to store whether node can be cut off */
   Bool*            infeasible,         /**< pointer to store whether LP/pseudo solution is infeasible */
   Bool*            solveagain          /**< pointer to store whether node has to be solved again */
   )
{
   RESULT result;
   Real pseudoobjval;
   Bool resolved;
   Bool enforceagain;
   Bool objinfeasible;
   int h;

   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(prob != NULL);
   assert(primal != NULL);
   assert(conshdlrs_enfo != NULL);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(solveagain != NULL);

   /* enforce constraints: branching, separating, reducing domains */
   assert(!(*cutoff));
   assert(!(*infeasible));
   assert(!(*solveagain));

   /**@todo avoid checking the same pseudosolution twice */

   /* enforce constraints by branching, applying additional cutting planes (if LP is being processed),
    * introducing new constraints, or tighten the domains
    */
   debugMessage("enforcing constraints on %s solution\n", tree->actnodehaslp ? "LP" : "pseudo");

   /* reset the constraint handler's enforcement calls */
   for( h = 0; h < set->nconshdlrs; ++h )
      SCIPconshdlrResetEnfo(set->conshdlrs[h]);

   /* check, if the solution is infeasible anyway due to it's objective value */
   if( tree->actnodehaslp )
      objinfeasible = FALSE;
   else
   {
      pseudoobjval = SCIPlpGetPseudoObjval(lp, set);
      objinfeasible = SCIPsetIsLT(set, pseudoobjval, tree->actnode->lowerbound);
   }

   /* enforce constraints until a handler resolved an infeasibility with cutting off the node, branching, 
    * reducing a domain, or separating a cut
    * if a constraint handler introduced new constraints to enforce his constraints, the newly added constraints
    * have to be enforced themselves
    */
   do
   {
      enforceagain = FALSE;
      resolved = FALSE;
      for( h = 0; h < set->nconshdlrs && !resolved; ++h )
      {
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);

         if( tree->actnodehaslp )
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
            assert(!tree->actnodehaslp || lp->solved);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            assert(objinfeasible);
            *infeasible = TRUE;
            break;

         case SCIP_FEASIBLE:
            assert(tree->nchildren == 0);
            assert(!tree->actnodehaslp || lp->solved);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            break;

         case SCIP_INFEASIBLE:
            assert(tree->nchildren == 0);
            assert(!tree->actnodehaslp || lp->solved);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            *infeasible = TRUE;
            break;

         case SCIP_CUTOFF:
            assert(tree->nchildren == 0);
            assert(!tree->actnodehaslp || lp->solved);

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
            CHECK_OKAY( SCIPsepastoreApplyCuts(sepastore, memhdr, set, stat, tree, lp, branchcand, eventqueue, cutoff) );

            *infeasible = TRUE;
            *solveagain = !(*cutoff);
            resolved = TRUE;
            break;

         case SCIP_REDUCEDDOM:
            assert(tree->nchildren == 0);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            *infeasible = TRUE;
            *solveagain = TRUE;
            resolved = TRUE;
            break;

         case SCIP_CONSADDED:
            assert(tree->nchildren == 0);
            assert(!tree->actnodehaslp || lp->solved);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            *infeasible = TRUE;
            resolved = TRUE;
            enforceagain = TRUE; /* the newly added constraints have to be enforced themselves */
            break;

         case SCIP_BRANCHED:
            assert(tree->nchildren >= 1);
            assert(!tree->actnodehaslp || lp->solved);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            *infeasible = TRUE;
            resolved = TRUE;
            break;

         case SCIP_SOLVELP:
            assert(!tree->actnodehaslp);
            assert(tree->nchildren == 0);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            assert(!enforceagain);
            *infeasible = TRUE;
            *solveagain = TRUE;
            resolved = TRUE;
            tree->actnodehaslp = TRUE; /* the node's LP must be solved */
            break;

         default:
            errorMessage("invalid result code <%d> from enforcing method of constraint handler <%s>\n",
               result, SCIPconshdlrGetName(conshdlrs_enfo[h]));
            return SCIP_INVALIDRESULT;
         }  /*lint !e788*/
         assert(!(*solveagain) || (resolved && *infeasible));
      }
      assert(!objinfeasible || *infeasible);
      debugMessage(" -> enforcing result: infeasible=%d, solveagain=%d, resolved=%d, enforceagain=%d\n",
         *infeasible, *solveagain, resolved, enforceagain);
   }
   while( enforceagain );
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
      assert(!(*solveagain));
   
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
         *solveagain = TRUE;
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

         if( prob->ncontvars == 0 && set->nactivepricers == 0 )
         {
            *cutoff = TRUE;
            resolved = TRUE;
         }
         else
         {
            assert(!tree->actnodehaslp);
            
            /* solve the LP in the next loop */
            tree->actnodehaslp = TRUE;
            *solveagain = TRUE;
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

/** solves the current node */
static
RETCODE solveNode(
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICESTORE*      pricestore,         /**< pricing storage */
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   CONFLICT*        conflict,           /**< conflict analysis data */
   LPCONFLICT*      lpconflict,         /**< conflict analysis data for infeasible LP conflicts */
   PRIMAL*          primal,             /**< primal data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   CONSHDLR**       conshdlrs_sepa,     /**< constraint handlers for separating constraints, sorted by priority */
   CONSHDLR**       conshdlrs_enfo,     /**< constraint handlers for enforcing constraints, sorted by priority */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   Bool*            infeasible          /**< pointer to store whether the current node's solution is infeasible */
   )
{
   Bool initiallpsolved;
   Bool solveagain;
   int nlperrors;

   assert(set != NULL);
   assert(prob != NULL);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(primal != NULL);
   assert(cutoff != NULL);
   assert(infeasible != NULL);

   debugMessage("current pseudosolution: obj=%g", SCIPlpGetPseudoObjval(lp, set));
   debug(SCIPprobPrintPseudoSol(prob, set));
   
   /* check, if we want to solve the LP at the selected node:
    * - solve the LP, if the lp solve depth and frequency demand solving
    * - solve the root LP, if the LP solve frequency is set to 0
    * - solve the root LP, if there are continuous variables present
    * - don't solve the node if its cut off by the pseudo objective value anyway
    */
   tree->actnodehaslp = (set->lpsolvedepth == -1 || SCIPnodeGetDepth(tree->actnode) <= set->lpsolvedepth);
   tree->actnodehaslp = tree->actnodehaslp
      && (set->lpsolvefreq >= 1 && SCIPnodeGetDepth(tree->actnode) % set->lpsolvefreq == 0);
   tree->actnodehaslp = tree->actnodehaslp || (SCIPnodeGetDepth(tree->actnode) == 0 && set->lpsolvefreq == 0);
   tree->actnodehaslp = tree->actnodehaslp || (SCIPnodeGetDepth(tree->actnode) == 0 && prob->ncontvars > 0);
   tree->actnodehaslp = tree->actnodehaslp && SCIPsetIsLT(set, SCIPlpGetPseudoObjval(lp, set), primal->cutoffbound);
         
   /* external node solving loop */
   initiallpsolved = FALSE;
   nlperrors = 0;
   do
   {
      Real pseudoobjval;
      Bool lperror;

      *infeasible = FALSE;
      solveagain = FALSE;
      lperror = FALSE;

      stat->npricerounds = 0;
      stat->nseparounds = 0;

      /* check, if we want to solve the LP at this node */
      if( tree->actnodehaslp )
      {
         /* solve the LP with price and cut */
         CHECK_OKAY( priceAndCutLoop(memhdr, set, stat, prob, tree, lp, pricestore, sepastore,
                        cutpool, conflict, lpconflict, primal, branchcand, eventfilter, eventqueue, conshdlrs_sepa,
                        initiallpsolved, cutoff, &lperror) );
         initiallpsolved = TRUE;
         debugMessage(" -> LP status: %d, LP obj: %g\n", SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set));

         /* if an error occured during LP solving, switch to pseudo solution */
         if( lperror )
         {
            tree->actnodehaslp = FALSE;
            nlperrors++;
            infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
               "(node %lld) unresolved numerical troubles in LP %d -- using pseudo solution instead (loop %d)\n", 
               stat->nnodes, stat->nlps, nlperrors);
         }
               
         /* if we solve exactly, the LP claims to be infeasible but the infeasibility could not be proved,
          * we have to forget about the LP and use the pseudo solution instead
          */
         if( set->exactsolve && lp->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE
            && tree->actnode->lowerbound < primal->cutoffbound )
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
               tree->actnodehaslp = FALSE;
               infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) could not prove infeasibility of LP %d -- using pseudo solution (%d unfixed vars) instead\n", 
                  stat->nnodes, stat->nlps, SCIPbranchcandGetNPseudoCands(branchcand));
            }
         }
      }
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      /* update lower bound w.r.t. the pseudo solution */
      pseudoobjval = SCIPlpGetPseudoObjval(lp, set);
      SCIPnodeUpdateLowerbound(tree->actnode, pseudoobjval);
      debugMessage(" -> new lower bound: %g (pseudoobj: %g)\n", tree->actnode->lowerbound, pseudoobjval);
            
      /* check for infeasible node by bounding */
      if( *cutoff
         || tree->actnode->lowerbound >= primal->cutoffbound
         || (!set->exactsolve && SCIPsetIsGE(set, tree->actnode->lowerbound, primal->cutoffbound)) )
      {
         debugMessage("node is infeasible (lower=%g, upper=%g)\n", tree->actnode->lowerbound, primal->cutoffbound);
         SCIPnodeUpdateLowerbound(tree->actnode, primal->cutoffbound);
         *infeasible = TRUE;
      }
      else
      {
         Longint oldnboundchanges;

         assert(!tree->actnodehaslp || lp->solved);

         oldnboundchanges = stat->nboundchanges;

         /* enforce constraints */
         CHECK_OKAY( enforceConstraints(memhdr, set, stat, prob, tree, lp, sepastore, branchcand, primal, eventqueue,
                        conshdlrs_enfo, cutoff, infeasible, &solveagain) );
         assert(!solveagain || !(*cutoff));

         /* apply further domain propagation, if addional bounds have been changed and the node is not yet finished */
         if( solveagain && stat->nboundchanges != oldnboundchanges )
         {
            CHECK_OKAY( propagateDomains(memhdr, set, stat, prob, tree, cutoff) );
            if( *cutoff )
            {
               debugMessage("enforced domain propagation determined that node is infeasible\n");
               SCIPnodeUpdateLowerbound(tree->actnode, primal->cutoffbound);
               *infeasible = TRUE;
               solveagain = FALSE;
            }
         }
      }
      assert(!solveagain || !(*cutoff));
   }
   while( solveagain && nlperrors < MAXNLPERRORS );
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);

   if( nlperrors >= MAXNLPERRORS )
   {
      errorMessage("(node %lld) unresolved numerical troubles in LP %d -- aborting\n", stat->nnodes, stat->nlps);
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** calls primal heuristics */
static
RETCODE primalHeuristics(
   SET*             set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   PRIMAL*          primal,             /**< primal data */
   Bool*            foundsol            /**< pointer to store whether a solution has been found */
   )
{
   RESULT result;
   int depth;
   int lpforkdepth;
   int h;

   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(foundsol != NULL);

   *foundsol = FALSE;

   if( SCIPtreeGetNNodes(tree) > 0 )
   {
      /* sort heuristics by priority */
      SCIPsetSortHeurs(set);

      /* call heuristics */
      depth = SCIPnodeGetDepth(tree->actnode);
      lpforkdepth = tree->actlpfork != NULL ? SCIPnodeGetDepth(tree->actlpfork) : -1;
      for( h = 0; h < set->nheurs; ++h )
      {
         CHECK_OKAY( SCIPheurExec(set->heurs[h], set, primal, depth, lpforkdepth, tree->actnodehaslp, &result) );
         *foundsol = *foundsol || (result == SCIP_FOUNDSOL);
      }
   }

   return SCIP_OKAY;
}

/** if feasible, adds current solution to the solution storage */
static
RETCODE addCurrentSolution(
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   PRIMAL*          primal,             /**< primal data */
   EVENTFILTER*     eventfilter         /**< event filter for global (not variable dependent) events */
   )
{
   SOL* sol;
   Bool foundsol;

   /* found a feasible solution */
   if( tree->actnodehaslp )
   {
      /* start clock for LP solutions */
      SCIPclockStart(stat->lpsoltime, set);
      
      /* add solution to storage */
      CHECK_OKAY( SCIPsolCreateLPSol(&sol, memhdr, set, stat, tree, lp, NULL) );
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
      CHECK_OKAY( SCIPsolCreatePseudoSol(&sol, memhdr, set, stat, tree, lp, NULL) );
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
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICESTORE*      pricestore,         /**< pricing storage */
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   CONFLICT*        conflict,           /**< conflict analysis data */
   LPCONFLICT*      lpconflict,         /**< conflict analysis data for infeasible LP conflicts */
   PRIMAL*          primal,             /**< primal data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   CONSHDLR** conshdlrs_sepa;
   CONSHDLR** conshdlrs_enfo;
   NODESEL* nodesel;
   NODE* actnode;
   EVENT event;
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

   /* sort constraint handlers by priorities */
   ALLOC_OKAY( duplicateMemoryArray(&conshdlrs_sepa, set->conshdlrs, set->nconshdlrs) );
   ALLOC_OKAY( duplicateMemoryArray(&conshdlrs_enfo, set->conshdlrs, set->nconshdlrs) );
   SCIPbsortPtr((void**)conshdlrs_sepa, set->nconshdlrs, SCIPconshdlrCompSepa);
   SCIPbsortPtr((void**)conshdlrs_enfo, set->nconshdlrs, SCIPconshdlrCompEnfo);

   while( !SCIPsolveIsStopped(set, stat) )
   {
      /* update the memory saving flag, switch algorithms respectively */
      SCIPstatUpdateMemsaveMode(stat, set, mem);

      /* get the current node selector */
      nodesel = SCIPsetGetNodesel(set, stat);

      /* inform tree about the current node selector */
      CHECK_OKAY( SCIPtreeSetNodesel(tree, set, stat, nodesel) );

      /* select next node to process */
      CHECK_OKAY( SCIPnodeselSelect(nodesel, set, &actnode) );

      if( actnode != NULL )
      {
         int depth;

         /* update statistics */
         depth = SCIPnodeGetDepth(actnode);
         stat->maxdepth = MAX(stat->maxdepth, depth);
         stat->nnodes++;
         if( SCIPnodeGetType(actnode) == SCIP_NODETYPE_LEAF )
         {
            stat->plungedepth = 0;
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
      assert(tree->actnode == actnode);

      /* process the delayed events */
      CHECK_OKAY( SCIPeventqueueProcess(eventqueue, memhdr, set, lp, branchcand, eventfilter) );

      /* issue NODEACTIVATED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEACTIVATED) );
      CHECK_OKAY( SCIPeventChgNode(&event, actnode) );
      CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, eventfilter) );

      /* stop node activation timer */
      SCIPclockStop(stat->nodeactivationtime, set);

      /* if no more node was selected, we finished optimization */
      if( actnode == NULL )
         break;

      debugMessage("Processing node %lld in depth %d, %d siblings\n", 
         stat->nnodes, SCIPnodeGetDepth(actnode), tree->nsiblings);

      /* domain propagation */
      CHECK_OKAY( propagateDomains(memhdr, set, stat, prob, tree, &cutoff) );

      if( cutoff )
      {
         debugMessage("domain propagation determined that node is infeasible\n");
         SCIPnodeUpdateLowerbound(tree->actnode, primal->cutoffbound);
         tree->actnodehaslp = FALSE;
         infeasible = TRUE;
      }
      else
      {
         CHECK_OKAY( solveNode(memhdr, set, stat, prob, tree, lp, pricestore, sepastore, branchcand, cutpool,
                        conflict, lpconflict, primal, eventfilter, eventqueue, conshdlrs_sepa, conshdlrs_enfo,
                        &cutoff, &infeasible) );
      }
      assert(!cutoff || infeasible);

      /* change color of node in VBC output */
      SCIPvbcSolvedNode(stat->vbc, stat, tree->actnode);

      /* check, if the current solution is feasible */
      if( !infeasible )
      {
         assert(!tree->actnodehaslp || lp->solved);
         assert(!cutoff);

         /* node solution is feasible: add it to the solution store */
         CHECK_OKAY( addCurrentSolution(memhdr, set, stat, prob, tree, lp, branchcand, primal, eventfilter) );

         /* issue NODEFEASIBLE event */
         CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEFEASIBLE) );
         CHECK_OKAY( SCIPeventChgNode(&event, actnode) );
         CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, eventfilter) );
      }
      else
      {
         /* node solution is not feasible */
         if( tree->nchildren == 0 )
         {
            /* issue NODEINFEASIBLE event */
            CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEINFEASIBLE) );
         }
         else
         {
            /* issue NODEBRANCHED event */
            CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEBRANCHED) );
         }
         CHECK_OKAY( SCIPeventChgNode(&event, actnode) );
         CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, eventfilter) );
      }

      /* if no branching was created, the node was not cut off, but it's lower bound is still smaller than
       * the cutoff bound, we have to branch on a non-fixed variable;
       * this can happen, if we want to solve exactly, the current solution was declared feasible by the
       * constraint enforcement, but in exact solution checking it was found out to be infeasible;
       * in this case, no branching would have been generated by the enforcement of constraints, but we
       * have to further investigate the current sub tree
       */
      if( !cutoff && tree->nchildren == 0 && tree->actnode->lowerbound < primal->cutoffbound )
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

      /* call primal heuristics */
      CHECK_OKAY( primalHeuristics(set, tree, primal, &foundsol) );

      /* display node information line */
      CHECK_OKAY( SCIPdispPrintLine(set, stat, (SCIPnodeGetDepth(actnode) == 0) && infeasible && !foundsol) );
      
      debugMessage("Processing of node in depth %d finished. %d siblings, %d children, %d leaves left\n", 
         SCIPnodeGetDepth(actnode), tree->nsiblings, tree->nchildren, SCIPtreeGetNLeaves(tree));
      debugMessage("**********************************************************************\n");
   }

   debugMessage("Problem solving finished\n");

   /* free sorted constraint handler arrays */
   freeMemoryArrayNull(&conshdlrs_sepa);
   freeMemoryArrayNull(&conshdlrs_enfo);

   return SCIP_OKAY;
}

