/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: solve.c,v 1.71 2003/11/27 17:48:48 bzfpfend Exp $"

/**@file   solve.c
 * @brief  main solving loop and node processing
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "var.h"
#include "cons.h"
#include "disp.h"
#include "clock.h"
#include "interrupt.h"
#include "solve.h"


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
         CHECK_OKAY( SCIPconshdlrPropagate(set->conshdlrs[h], memhdr, set, prob, SCIPnodeGetDepth(tree->actnode),
                        &result) );
         propagain = propagain || (result == SCIP_REDUCEDDOM);
         *cutoff = *cutoff || (result == SCIP_CUTOFF);
      }
   }
   while( propagain && !(*cutoff) );

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
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   VAR* var;
   COL* col;
   int v;
   int h;

   assert(lp != NULL);
   assert(lp->ncols == 0);
   assert(lp->nrows == 0);
   assert(lp->nremoveablecols == 0);
   assert(lp->nremoveablerows == 0);

   debugMessage("init root LP\n");

   /* add all unremoveable variables to LP */
   for( v = 0; v < prob->nvars; ++v )
   {
      var = prob->vars[v];
      assert(SCIPvarGetProbIndex(var) >= 0);
      
      if( !SCIPvarIsRemoveable(var) )
      {
         /* transform variable into column variable, if needed */
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
         {
            CHECK_OKAY( SCIPvarColumn(var, memhdr, set, stat) );
         }
         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
         
         col = SCIPvarGetCol(var);
         assert(col != NULL);
         assert(col->lppos == -1);
         assert(col->lpipos == -1);
         debugMessage("adding initial variable <%s>\n", SCIPvarGetName(var));
         CHECK_OKAY( SCIPlpAddCol(lp, set, col) );
      }
   }
   assert(lp->nremoveablecols == 0);

   debugMessage("adding initial rows\n");

   /* add LP relaxations of all initial constraints to LP */
   for( h = 0; h < set->nconshdlrs; ++h )
   {
      CHECK_OKAY( SCIPconshdlrInitLP(set->conshdlrs[h], memhdr, set, prob, sepastore) );
   }

   /* apply found cuts */
   CHECK_OKAY( SCIPsepastoreApplyCuts(sepastore, memhdr, set, stat, tree, lp, branchcand, eventqueue) );

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
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   EVENT event;

   assert(stat != NULL);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(lp != NULL);

   /* load the LP into the solver and load the LP state */
   debugMessage("loading LP\n");
   CHECK_OKAY( SCIPtreeLoadLP(tree, memhdr, set, lp) );
   
   /* init root node LP */
   if( SCIPnodeGetDepth(tree->actnode) == 0 )
   {
      assert(stat->lpcount == 0);
      CHECK_OKAY( initRootLP(memhdr, set, stat, prob, tree, lp, sepastore, branchcand, eventqueue) );
   }

   /* only the bounds and the constraint list may have been changed since the LP state node
    * -> dual simplex is applicable as first solver
    */
   debugMessage("node: solve initial LP\n");
   CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, TRUE) );
   assert(lp->solved);

   /* issue FIRSTLPSOLVED event */
   CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_FIRSTLPSOLVED) );
   CHECK_OKAY( SCIPeventChgNode(&event, tree->actnode) );
   CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );

   return SCIP_OKAY;
}

/** solve the actual LP of a node with price and cut */
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
   LPCONFLICT*      lpconflict,         /**< conflict analysis data for infeasible LP conflicts */
   PRIMAL*          primal,             /**< primal data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   CONSHDLR**       conshdlrs_sepa,     /**< constraint handlers sorted by separation priority */
   Bool*            cutoff              /**< pointer to store whether the node can be cut off */
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

   root = (SCIPnodeGetDepth(tree->actnode) == 0);

   debugMessage("node: solve LP with price and cut\n");
   CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, TRUE) );
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
   while( !(*cutoff) && (mustprice || mustsepar) )
   {
      debugMessage("-------- node solving loop --------\n");
      assert(lp->solved);

      /* if the LP is unbounded, we don't need to price */
      mustprice = mustprice && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDED);
      /* if all the variables are already in the LP, we don't need to price */
      mustprice = mustprice && (prob->nvars > lp->ncols || set->npricers > 0);

      /* pricing (has to be done completely to get a valid lower bound) */
      while( !(*cutoff) && mustprice )
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
            
            /* call external pricer algorithms */
            enoughvars = FALSE;
            for( p = 0; p < set->npricers && !enoughvars; ++p )
            {
               CHECK_OKAY( SCIPpricerExec(set->pricers[p], set, prob, lp) );
               enoughvars = enoughvars || (SCIPpricestoreGetNVars(pricestore) >= SCIPsetGetMaxpricevars(set, root)/2);
            }
         }

         /* apply the priced variables to the LP */
         CHECK_OKAY( SCIPpricestoreApplyVars(pricestore, memhdr, set, stat, prob, tree, lp, branchcand, eventqueue) );
         assert(SCIPpricestoreGetNVars(pricestore) == 0);
         mustprice = !lp->solved;
         mustsepar = mustsepar || !lp->solved;
         
         /* after adding columns, the LP should be primal feasible such that primal simplex is applicable;
          * if LP was infeasible, we have to use dual simplex
          */
         debugMessage("pricing: solve LP\n");
         CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, TRUE) );
         assert(lp->solved);

         /* reset bounds temporarily set by pricer to their original values */
         debugMessage("pricing: reset bounds\n");
         CHECK_OKAY( SCIPpricestoreResetBounds(pricestore, memhdr, set, stat, tree, lp, branchcand, eventqueue) );
         assert(SCIPpricestoreGetNVars(pricestore) == 0);
         assert(SCIPpricestoreGetNBoundResets(pricestore) == 0);
         mustprice = mustprice || !lp->solved;
         mustsepar = mustsepar || !lp->solved;

         /* solve LP again after resetting bounds (with dual simplex) */
         debugMessage("pricing: solve LP after resetting bounds\n");
         CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, FALSE) );
         assert(lp->solved);

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
      tree->actnode->lowerbound = SCIPlpGetObjval(lp);
      tree->actnode->lowerbound = MIN(tree->actnode->lowerbound, primal->upperbound);
      debugMessage(" -> new lower bound: %g\n", tree->actnode->lowerbound);

      /* if the LP is infeasible or exceeded the objective limit, we don't need to separate cuts */
      mustsepar = mustsepar
         && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_INFEASIBLE)
         && (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OBJLIMIT);

      /* separation (has not to be done completely, because we just want to increase the lower bound) */
      if( !(*cutoff) && mustsepar )
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
                  CHECK_OKAY( SCIPconshdlrSeparate(conshdlrs_sepa[h], memhdr, set, prob, sepastore, 
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
            CHECK_OKAY( SCIPsepastoreApplyCuts(sepastore, memhdr, set, stat, tree, lp, branchcand, eventqueue) );

            /* if at least one of the cut was a bound change, propagate domains again */
            if( stat->nboundchanges != oldnboundchanges )
            {
               CHECK_OKAY( propagateDomains(memhdr, set, prob, tree, cutoff) );
            }

            mustprice = mustprice || !lp->solved;
            mustsepar = !lp->solved;

            if( !(*cutoff) )
            {
               /* solve LP (with dual simplex) */
               debugMessage("separation: solve LP\n");
               CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, TRUE) );
            }
         }
         assert(*cutoff || lp->solved); /* after cutoff, the LP may be unsolved due to bound changes */

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
      tree->actnode->lowerbound = primal->upperbound;
   }
   else
   {
      assert(lp->solved);

      tree->actnode->lowerbound = SCIPlpGetObjval(lp);
      tree->actnode->lowerbound = MIN(tree->actnode->lowerbound, primal->upperbound);

      /* issue LPSOLVED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_LPSOLVED) );
      CHECK_OKAY( SCIPeventChgNode(&event, tree->actnode) );
      CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );
   }
   debugMessage(" -> new lower bound: %g\n", tree->actnode->lowerbound);

   /* analyze an infeasible LP */
   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE
      || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT )
   {
      CHECK_OKAY( SCIPlpconflictAnalyze(lpconflict, set, prob, lp, NULL) );
   }

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
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   PRIMAL*          primal,             /**< primal data */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   VAR* var;
   int* cstat;
   Real lpobjval;
   Real redcost;
   Real oldlb;
   Real oldub;
   Real newbd;
   int c;

   assert(tree != NULL);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(primal != NULL);
   assert(!SCIPsetIsInfinity(set, primal->upperbound));
   assert(SCIPsetIsLT(set, SCIPlpGetObjval(lp), primal->upperbound) );

   stat->nredcoststrcalls++;

   /* start redcost strengthening timer */
   SCIPclockStart(stat->redcoststrtime, set);

   /* get temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &cstat, lp->ncols) );

   /* get basis status for columns and LP objective value */
   CHECK_OKAY( SCIPlpGetBase(lp, cstat, NULL) );
   lpobjval = SCIPlpGetObjval(lp);

   /* check reduced costs for non-basic columns */
   for( c = 0; c < lp->ncols; ++c )
   {
      switch( cstat[c] )
      {
      case SCIP_BASESTAT_LOWER:
         redcost = SCIPcolGetRedcost(lp->cols[c], stat);
         assert(!SCIPsetIsFeasNegative(set, redcost));
         if( SCIPsetIsFeasPositive(set, redcost) )
         {
            var = SCIPcolGetVar(lp->cols[c]);
            oldlb = SCIPvarGetLbLocal(var);
            oldub = SCIPvarGetUbLocal(var);
            if( SCIPsetIsFeasLT(set, oldlb, oldub) )
            {
               /* calculate reduced cost based bound */
               newbd = (primal->upperbound - lpobjval) / redcost + oldlb;
               SCIPvarAdjustUb(var, set, &newbd);

               /* continuous variables should be tightened at least by 10 percent */
               if( (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && newbd < oldub - 0.5)
                  || newbd < oldub - 0.1*(oldub - oldlb) )
               {
                  /* strengthen upper bound */
                  debugMessage("redcost strengthening upper bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
                     SCIPvarGetName(var), oldlb, oldub, oldlb, newbd, primal->upperbound, lpobjval, redcost);
                  CHECK_OKAY( SCIPnodeAddBoundchg(tree->actnode, memhdr, set, stat, tree, lp, branchcand, eventqueue,
                                 var, newbd, SCIP_BOUNDTYPE_UPPER, NULL) );
                  stat->nredcoststrfound++;
               }
            }
         }
         break;

      case SCIP_BASESTAT_BASIC:
         break;

      case SCIP_BASESTAT_UPPER:
         redcost = SCIPcolGetRedcost(lp->cols[c], stat);
         assert(!SCIPsetIsFeasPositive(set, redcost));
         if( SCIPsetIsFeasNegative(set, redcost) )
         {
            var = SCIPcolGetVar(lp->cols[c]);
            oldlb = SCIPvarGetLbLocal(var);
            oldub = SCIPvarGetUbLocal(var);
            if( SCIPsetIsFeasLT(set, oldlb, oldub) )
            {
               /* calculate reduced cost based bound */
               newbd = (primal->upperbound - lpobjval) / redcost + oldub;
               SCIPvarAdjustLb(var, set, &newbd);

               /* continuous variables should be tightened at least by 10 percent */
               if( (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && newbd > oldlb + 0.5)
                  || newbd > oldlb + 0.1*(oldub - oldlb) )
               {
                  /* strengthen lower bound */
                  debugMessage("redcost strengthening lower bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
                     SCIPvarGetName(var), oldlb, oldub, newbd, oldub, primal->upperbound, lpobjval, redcost);
                  CHECK_OKAY( SCIPnodeAddBoundchg(tree->actnode, memhdr, set, stat, tree, lp, branchcand, eventqueue,
                                 var, newbd, SCIP_BOUNDTYPE_LOWER, NULL) );
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

   /* resolve the LP (the optimal solution should stay the same) */
   CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, FALSE) );
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(SCIPsetIsFeasEQ(set, SCIPlpGetObjval(lp), lpobjval));

   /* stop redcost strengthening activation timer */
   SCIPclockStop(stat->redcoststrtime, set);

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
   Bool*            infeasible,         /**< pointer to store whether LP/pseudo solution is infeasible */
   Bool*            solveagain          /**< pointer to store whether node has to be solved again */
   )
{
   RESULT result;
   Real actpseudoobjval;
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
   assert(infeasible != NULL);
   assert(solveagain != NULL);

   /* enforce constraints: branching, separating, reducing domains */
   assert(!(*infeasible));
   assert(!(*solveagain));

   /**@todo avoid checking the same pseudosolution twice */

   /* enforce constraints by branching, applying additional cutting planes (if LP is being processed),
    * introducing new constraints, or tighten the domains
    */
   debugMessage("enforcing constraints\n");

   /* reset the constraint handler's enforcement calls */
   for( h = 0; h < set->nconshdlrs; ++h )
      SCIPconshdlrResetEnfo(set->conshdlrs[h]);

   /* check, if the solution is infeasible anyway due to it's objective value */
   if( tree->actnodehaslp )
      objinfeasible = FALSE;
   else
   {
      actpseudoobjval = SCIPtreeGetActPseudoobjval(tree, set);
      objinfeasible = SCIPsetIsLT(set, actpseudoobjval, tree->actnode->lowerbound);
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
            CHECK_OKAY( SCIPconshdlrEnforceLPSol(conshdlrs_enfo[h], memhdr, set, prob, sepastore, &result) );
         }
         else
         {
            CHECK_OKAY( SCIPconshdlrEnforcePseudoSol(conshdlrs_enfo[h], memhdr, set, prob, objinfeasible, &result) );
            if( SCIPsepastoreGetNCuts(sepastore) != 0 )
            {
               errorMessage("pseudo enforcing method of constraint handler <%s> separated cuts\n",
                  SCIPconshdlrGetName(conshdlrs_enfo[h]));
               return SCIP_INVALIDRESULT;
            }
         }
         assert(set->buffer->firstfree == 0);
         debugMessage("enforcing of <%s> returned result %d\n", SCIPconshdlrGetName(conshdlrs_enfo[h]), result);

         switch( result )
         {  /*lint --e{788}*/
         case SCIP_DIDNOTRUN:
            assert(tree->nchildren == 0);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            assert(objinfeasible);
            *infeasible = TRUE;
            break;

         case SCIP_FEASIBLE:
            assert(tree->nchildren == 0);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            break;

         case SCIP_INFEASIBLE:
            assert(tree->nchildren == 0);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            *infeasible = TRUE;
            break;

         case SCIP_CUTOFF:
            assert(tree->nchildren == 0);

            /* the found cuts are of no use, because the node is infeasible anyway */
            CHECK_OKAY( SCIPsepastoreClearCuts(sepastore, memhdr, set, lp) );

            *infeasible = TRUE;
            resolved = TRUE;
            break;

         case SCIP_SEPARATED:
            assert(tree->nchildren == 0);
            assert(SCIPsepastoreGetNCuts(sepastore) > 0);

            /* apply found cuts */
            CHECK_OKAY( SCIPsepastoreApplyCuts(sepastore, memhdr, set, stat, tree, lp, branchcand, eventqueue) );

            *infeasible = TRUE;
            *solveagain = TRUE;
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
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            *infeasible = TRUE;
            resolved = TRUE;
            enforceagain = TRUE; /* the newly added constraints have to be enforced themselves */
            break;

         case SCIP_BRANCHED:
            assert(tree->nchildren >= 1);
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
         }
         assert(!(*solveagain) || (resolved && *infeasible));
      }
      assert(!objinfeasible || *infeasible);
      debugMessage(" -> enforcing result: infeasible=%d, solveagain=%d, resolved=%d, enforceagain=%d\n",
         *infeasible, *solveagain, resolved, enforceagain);
   }
   while( enforceagain );
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);
   assert(set->buffer->firstfree == 0);

   if( *infeasible && !resolved )
   {
      /* the node is infeasible, but no constraint handler could resolve the infeasibility
       * -> branch on the pseudo solution
       * -> e.g. select non-fixed binary or integer variable x with value x', create three
       *    sons: x <= x'-1, x = x', and x >= x'+1.
       *    In the left and right branch, the actual solution is cutted off. In the middle
       *    branch, the constraints can hopefully reduce domains of other variables to cut
       *    off the actual solution.
       */

      assert(!(*solveagain));
   
      /* branch on pseudo solution */
      CHECK_OKAY( SCIPbranchPseudo(set->scip, &result) );
      assert(set->buffer->firstfree == 0);

      switch( result )
      {  /*lint --e{788}*/
      case SCIP_CUTOFF:
         assert(tree->nchildren == 0);
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

         if( prob->ncont == 0 && set->npricers == 0 )
            resolved = TRUE;
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
      }
   }
   assert(*infeasible || !resolved);

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
   int actdepth;
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
      actdepth = SCIPnodeGetDepth(tree->actnode);
      lpforkdepth = tree->actlpfork != NULL ? SCIPnodeGetDepth(tree->actlpfork) : -1;
      for( h = 0; h < set->nheurs; ++h )
      {
         CHECK_OKAY( SCIPheurExec(set->heurs[h], set, primal, actdepth, lpforkdepth, tree->actnodehaslp, &result) );
         *foundsol = *foundsol || (result == SCIP_FOUNDSOL);
      }
   }

   return SCIP_OKAY;
}

/** main solving loop */
RETCODE SCIPsolveCIP(
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
   Bool solveagain;
   Bool initiallpsolved;
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

   /**@todo every variable, where zero is not the best bound (w.r.t. objective function) has to be in the problem */

   /* sort constraint handlers by priorities */
   ALLOC_OKAY( duplicateMemoryArray(&conshdlrs_sepa, set->conshdlrs, set->nconshdlrs) );
   ALLOC_OKAY( duplicateMemoryArray(&conshdlrs_enfo, set->conshdlrs, set->nconshdlrs) );
   SCIPbsortPtr((void**)conshdlrs_sepa, set->nconshdlrs, SCIPconshdlrCompSepa);
   SCIPbsortPtr((void**)conshdlrs_enfo, set->nconshdlrs, SCIPconshdlrCompEnfo);

   while( !SCIPsolveIsStopped(set, stat) )
   {
      assert(set->buffer->firstfree == 0);

      /* update the memory saving flag, switch algorithms respectively */
      SCIPstatUpdateMemsaveMode(stat, set);

      /* get the current node selector */
      nodesel = SCIPsetGetActNodesel(set, stat);

      /* inform tree about the current node selector */
      CHECK_OKAY( SCIPtreeSetNodesel(tree, set, stat, nodesel) );

      /* select next node to process */
      CHECK_OKAY( SCIPnodeselSelect(nodesel, set, &actnode) );
      assert(set->buffer->firstfree == 0);

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
      assert(set->buffer->firstfree == 0);

      /* activate selected node */
      CHECK_OKAY( SCIPnodeActivate(actnode, memhdr, set, stat, tree, lp, branchcand, eventqueue) );
      assert(tree->actnode == actnode);
      assert(set->buffer->firstfree == 0);

      /* process the delayed events */
      CHECK_OKAY( SCIPeventqueueProcess(eventqueue, memhdr, set, tree, lp, branchcand, eventfilter) );
      assert(set->buffer->firstfree == 0);

      /* issue NODEACTIVATED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEACTIVATED) );
      CHECK_OKAY( SCIPeventChgNode(&event, actnode) );
      CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );

      /* stop node activation timer */
      SCIPclockStop(stat->nodeactivationtime, set);

      /* if no more node was selected, we finished optimization */
      if( actnode == NULL )
         break;

      debugMessage("Processing node %lld in depth %d, %d siblings\n", 
         stat->nnodes, SCIPnodeGetDepth(actnode), tree->nsiblings);
      
      /* domain propagation */
      CHECK_OKAY( propagateDomains(memhdr, set, prob, tree, &cutoff) );

      if( cutoff )
      {
         debugMessage("domain propagation determined that node is infeasible\n");
         tree->actnode->lowerbound = primal->upperbound;
         infeasible = TRUE;
      }
      else
      {
         debugMessage("actual pseudosolution: obj=%g", SCIPtreeGetActPseudoobjval(tree, set));
         debug(SCIPprobPrintPseudoSol(prob, set));
         
         /* check, if we want to solve the LP at the selected node:
          * - solve the LP, if the lp solve depth and frequency demand solving
          * - solve the root LP, if the LP solve frequency is set to 0
          * - solve the root LP, if there are continuous variables present
          * - don't solve the node if its cut off by the pseudo objective value anyway
          */
         tree->actnodehaslp = (set->lpsolvedepth == -1 || SCIPnodeGetDepth(actnode) <= set->lpsolvedepth);
         tree->actnodehaslp = tree->actnodehaslp
            && (set->lpsolvefreq >= 1 && SCIPnodeGetDepth(actnode) % set->lpsolvefreq == 0);
         tree->actnodehaslp = tree->actnodehaslp || (SCIPnodeGetDepth(actnode) == 0 && set->lpsolvefreq == 0);
         tree->actnodehaslp = tree->actnodehaslp || (SCIPnodeGetDepth(actnode) == 0 && prob->ncont > 0);
         tree->actnodehaslp = tree->actnodehaslp
            && SCIPsetIsLT(set, SCIPtreeGetActPseudoobjval(tree, set), primal->upperbound);
         
         /* external node solving loop */
         initiallpsolved = FALSE;
         do
         {
            Real actpsobjval;

            solveagain = FALSE;
            infeasible = FALSE;
            
            stat->npricerounds = 0;
            stat->nseparounds = 0;

            /* check, if we want to solve the LP at this node */
            if( tree->actnodehaslp )
            {
               if( !initiallpsolved )
               {
                  /* load and solve the initial LP of the node */
                  CHECK_OKAY( solveNodeInitialLP(memhdr, set, stat, prob, tree, lp, sepastore, 
                                 branchcand, eventfilter, eventqueue) );
                  assert(lp->solved);
                  initiallpsolved = TRUE;
               }
               assert(SCIPsepastoreGetNCuts(sepastore) == 0);

               /* solve-redcost-propagate loop */
               do
               {
                  /* continue solving the LP with price and cut */
                  CHECK_OKAY( solveNodeLP(memhdr, set, stat, prob, tree, lp, pricestore, sepastore, cutpool, lpconflict,
                                 primal, branchcand, eventfilter, eventqueue, conshdlrs_sepa, &cutoff) );
                  assert(cutoff || lp->solved);
                  
                  /* reduced cost bound strengthening */
                  if( !cutoff && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL
                     && !SCIPsetIsInfinity(set, primal->upperbound)
                     && SCIPsetIsLT(set, SCIPlpGetObjval(lp), primal->upperbound) )
                  {
                     Longint oldnboundchanges;
                     
                     oldnboundchanges = stat->nboundchanges;

                     /* apply reduced cost bound strengthening */
                     CHECK_OKAY( redcostStrengthening(memhdr, set, stat, prob, tree, lp, branchcand, primal, eventqueue) );
                     assert(lp->solved);

                     /* apply further domain propagation, if addional bounds have been changed */
                     if( stat->nboundchanges != oldnboundchanges )
                     {
                        CHECK_OKAY( propagateDomains(memhdr, set, prob, tree, &cutoff) );
                        if( cutoff )
                        {
                           debugMessage("additional domain propagation determined that node is infeasible\n");
                           tree->actnode->lowerbound = primal->upperbound;
                        }
                     }
                  }
               }
               while( !cutoff && !lp->solved );
            }
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            assert(set->buffer->firstfree == 0);
            
            /* update lower bound w.r.t. the pseudo solution */
            actpsobjval = SCIPtreeGetActPseudoobjval(tree, set);
            tree->actnode->lowerbound = MAX(tree->actnode->lowerbound, actpsobjval);
            tree->actnode->lowerbound = MIN(tree->actnode->lowerbound, primal->upperbound);
            debugMessage(" -> new lower bound: %g (pseudoobj: %g)\n", tree->actnode->lowerbound, actpsobjval);
            
            /* check for infeasible node by bounding */
            if( cutoff || SCIPsetIsGE(set, tree->actnode->lowerbound, primal->upperbound) )
            {
               debugMessage("node is infeasible (lower=%g, upper=%g)\n", tree->actnode->lowerbound, primal->upperbound);
               infeasible = TRUE;
            }
            else
            {
               Longint oldnboundchanges;

               assert(!tree->actnodehaslp || lp->solved);

               oldnboundchanges = stat->nboundchanges;

               /* enforce constraints */
               CHECK_OKAY( enforceConstraints(memhdr, set, stat, prob, tree, lp, sepastore, branchcand, primal, eventqueue,
                              conshdlrs_enfo, &infeasible, &solveagain) );

               /* apply further domain propagation, if addional bounds have been changed and the node is not yet finished */
               if( solveagain && stat->nboundchanges != oldnboundchanges )
               {
                  CHECK_OKAY( propagateDomains(memhdr, set, prob, tree, &cutoff) );
                  solveagain = !cutoff;
               }
            }
            assert(!solveagain || !cutoff);
         }
         while( solveagain );
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);

         if( !infeasible )
         {
            SOL* sol;
            Longint oldnsolsfound;

            assert(!tree->actnodehaslp || lp->solved);

            /* issue NODEFEASIBLE event */
            CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEFEASIBLE) );
            CHECK_OKAY( SCIPeventChgNode(&event, actnode) );
            CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );
            assert(set->buffer->firstfree == 0);
               
            /* found a feasible solution */
            oldnsolsfound = primal->nsolsfound;
            if( tree->actnodehaslp )
            {
               /* start clock for LP solutions */
               SCIPclockStart(stat->lpsoltime, set);

               /* add solution to storage */
               CHECK_OKAY( SCIPsolCreateLPSol(&sol, memhdr, stat, tree, lp, NULL) );
               CHECK_OKAY( SCIPprimalAddSolFree(primal, memhdr, set, stat, prob, tree, lp, eventfilter, &sol) );
               stat->nlpsolsfound += primal->nsolsfound - oldnsolsfound;

               /* stop clock for LP solutions */
               SCIPclockStop(stat->lpsoltime, set);
            }
            else
            {
               /* start clock for pseudo solutions */
               SCIPclockStart(stat->pseudosoltime, set);

               /* add solution to storage */
               CHECK_OKAY( SCIPsolCreatePseudoSol(&sol, memhdr, set, stat, tree, NULL) );
               CHECK_OKAY( SCIPprimalAddSolFree(primal, memhdr, set, stat, prob, tree, lp, eventfilter, &sol) );
               stat->npssolsfound += primal->nsolsfound - oldnsolsfound;

               /* stop clock for pseudo solutions */
               SCIPclockStop(stat->pseudosoltime, set);
            }
         }
      }

      /* if node solution is not feasible, issue NODEINFEASIBLE or NODEBRANCHED event */
      if( infeasible )
      {
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
         CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );
         assert(set->buffer->firstfree == 0);
      }

      /* call primal heuristics */
      CHECK_OKAY( primalHeuristics(set, tree, primal, &foundsol) );
      assert(set->buffer->firstfree == 0);
      
      /* display node information line */
      CHECK_OKAY( SCIPdispPrintLine(set, stat, (SCIPnodeGetDepth(actnode) == 0) && infeasible && !foundsol) );
      assert(set->buffer->firstfree == 0);
      
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

