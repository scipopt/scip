/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   solve.c
 * @brief  main solving loop and node processing
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "cons.h"
#include "var.h"
#include "disp.h"
#include "solve.h"


static
RETCODE initRootLP(                     /**< constructs the LP of the root node */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   LP*              lp,                 /**< LP data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   VAR* var;
   COL* col;
   int v;

   debugMessage("init root LP\n");
   if( !set->usepricing )
   {
      /* if we do not use pricing, add the all variables to LP */
      for( v = 0; v < prob->nvars; ++v )
      {
         var = prob->vars[v];
         assert(var->probindex >= 0);
         
         /* transform variable into column variable, if needed */
         if( var->varstatus == SCIP_VARSTATUS_LOOSE )
         {
            CHECK_OKAY( SCIPvarColumn(var, memhdr, set, lp, stat) );
         }
         assert(var->varstatus == SCIP_VARSTATUS_COLUMN);
         
         col = var->data.col;
         assert(col != NULL);
         assert(!col->inlp);
         assert(col->lpipos == -1);
         debugMessage("adding initial variable <%s>\n", var->name);
         CHECK_OKAY( SCIPlpAddCol(lp, set, col) );
      }
   }

   return SCIP_OKAY;
}

static
RETCODE solveLP(                        /**< solves the LP with simplex algorithm */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   STAT*            stat                /**< problem statistics */
   )
{
   assert(lp != NULL);

   debugMessage("solving LP: primalfeasible=%d, dualfeasible=%d, solved=%d\n", 
      lp->primalfeasible, lp->dualfeasible, lp->solved);
   if( !lp->solved )
   {
      if( lp->dualfeasible || !lp->primalfeasible )
      {
         debugMessage("solving dual LP\n");
         CHECK_OKAY( SCIPlpSolveDual(lp, memhdr, set, stat) );
      }
      else
      {
         debugMessage("solving primal LP\n");
         CHECK_OKAY( SCIPlpSolvePrimal(lp, memhdr, set, stat) );
      }

      switch( SCIPlpGetSolstat(lp) )
      {
      case SCIP_LPSOLSTAT_OPTIMAL:
         CHECK_OKAY( SCIPlpGetSol(lp, memhdr, set, stat) );
         break;
      case SCIP_LPSOLSTAT_INFEASIBLE:
         if( set->usepricing )
         {
            CHECK_OKAY( SCIPlpGetDualfarkas(lp, memhdr, set) );
         }
         break;
      case SCIP_LPSOLSTAT_UNBOUNDED:
         CHECK_OKAY( SCIPlpGetUnboundedSol(lp, memhdr, set, stat) );
         break;
      case SCIP_LPSOLSTAT_OBJLIMIT:
         if( set->usepricing )
         {
            CHECK_OKAY( SCIPlpGetSol(lp, memhdr, set, stat) );
         }
         break;
      case SCIP_LPSOLSTAT_ITERLIMIT:
         errorMessage("LP solver reached iteration limit -- this should not happen!");
         return SCIP_ERROR;
      case SCIP_LPSOLSTAT_TIMELIMIT:
         todoMessage("time limit exceeded processing");
         return SCIP_ERROR;
      case SCIP_LPSOLSTAT_ERROR:
         errorMessage("Error in LP solver");
         return SCIP_LPERROR;
      default:
         errorMessage("Unknown LP solution status");
         return SCIP_ERROR;
      }
   }

   return SCIP_OKAY;
}

static
RETCODE solveNodeLP(                    /**< solve a single node with price and cut */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICE*           price,              /**< pricing storage */
   SEPA*            sepa,               /**< separation storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   PRIMAL*          primal,             /**< primal data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   CONSHDLR**       conshdlrs_sepa      /**< constraint handlers sorted by separation priority */
   )
{
   RESULT result;
   Bool mustprice;
   Bool mustsepar;
   Bool root;
   int h;

   assert(set != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(lp != NULL);
   assert(price != NULL);
   assert(sepa != NULL);
   assert(cutpool != NULL);
   assert(primal != NULL);
   assert(set->nconshdlrs == 0 || conshdlrs_sepa != NULL);

   /* load the LP into the solver and load the LP state */
   debugMessage("loading LP\n");
   CHECK_OKAY( SCIPtreeLoadLP(tree, memhdr, set, lp) );
   
   /* init root node LP */
   root = (tree->actnode->depth == 0);
   if( root )
   {
      assert(stat->nlp == 0);
      CHECK_OKAY( initRootLP(memhdr, set, stat, prob, lp, tree) );
   }

   /* only the bounds and the constraint list may have been changed since the LP state node
    * -> dual simplex is applicable as first solver
    */
   debugMessage("node: solve LP\n");
   CHECK_OKAY( solveLP(memhdr, set, lp, stat) );
   assert(lp->solved);

   /* price-and-cut loop */
   mustprice = TRUE;
   mustsepar = TRUE;
   while( (set->usepricing && mustprice) || mustsepar )
   {
      debugMessage("-------- node solving loop --------\n");

      /* if the LP is unbounded, we don't need to price */
      mustprice &= (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDED);

      /* pricing (has to be done completely to get a valid lower bound) */
      while( set->usepricing && mustprice )
      {
         assert(lp->solved);
         assert(lp->lpsolstat != SCIP_LPSOLSTAT_UNBOUNDED);

         debugMessage("pricing:\n");
         CHECK_OKAY( SCIPpriceVars(price, memhdr, set, stat, prob, lp, tree, branchcand, eventqueue) );
         mustprice = !lp->solved;
         mustsepar |= !lp->solved;
         
         /* after adding columns, the LP should be primal feasible such that primal simplex is applicable;
          * if LP was infeasible, we have to use dual simplex
          */
         debugMessage("pricing: solve LP\n");
         CHECK_OKAY( solveLP(memhdr, set, lp, stat) );
         assert(lp->solved);

         /* reset bounds temporarily set by pricer to their original values */
         debugMessage("reset bounds\n");
         CHECK_OKAY( SCIPpriceResetBounds(price, memhdr, set, stat, lp, tree, branchcand, eventqueue) );
         mustprice |= !lp->solved;
         mustsepar |= !lp->solved;

         /* solve LP (with dual simplex) */
         debugMessage("reset bounds: solve LP\n");
         CHECK_OKAY( solveLP(memhdr, set, lp, stat) );
         assert(lp->solved);

         /* if the LP is unbounded, we can stop pricing */
         mustprice &= (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDED);
      }
      assert(lp->solved);

      /* update lower bound w.r.t. the the LP solution */
      tree->actnode->lowerbound = SCIPlpGetObjval(lp);
      tree->actnode->lowerbound = MIN(tree->actnode->lowerbound, primal->upperbound);

      /* if the LP is infeasible, we don't need to separate cuts */
      mustsepar &= (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_INFEASIBLE);

      /* separation (has not to be done completely, because we just want to increase the lower bound) */
      if( mustsepar )
      {
         assert(lp->solved);

         /* global cut pool separation */
         debugMessage("global cut pool separation\n");
         assert(SCIPsepaGetNCuts(sepa) == 0);
         CHECK_OKAY( SCIPcutpoolSeparate(cutpool, memhdr, set, stat, lp, sepa, root) );

         /* constraint separation */
         debugMessage("constraint separation\n");
         for( h = 0; h < set->nconshdlrs && SCIPsepaGetNCuts(sepa) == 0; ++h )
         {
            CHECK_OKAY( SCIPconshdlrSeparate(conshdlrs_sepa[h], set, &result) );
         }
         
         /* cut separation */
         if( SCIPsepaGetNCuts(sepa) == 0 )
         {
            todoMessage("cut separation");
         }

         /* apply found cuts */
         CHECK_OKAY( SCIPsepaApplyCuts(sepa, memhdr, set, tree, lp) );
         mustprice |= !lp->solved;
         mustsepar = !lp->solved;
         
         /* solve LP (with dual simplex) */
         debugMessage("separation: solve LP\n");
         CHECK_OKAY( solveLP(memhdr, set, lp, stat) );
         assert(lp->solved);
      }
   }

   /* update lower bound w.r.t. the the LP solution */
   tree->actnode->lowerbound = SCIPlpGetObjval(lp);
   tree->actnode->lowerbound = MIN(tree->actnode->lowerbound, primal->upperbound);

   return SCIP_OKAY;
}

static
RETCODE pseudoBranch(                   /**< branches on infeasible pseudo solution */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   VAR** pseudocands;
   int npseudocands;
   
   assert(result != NULL);

   /* get pseudo branching candidates (the non-fixed binary/integer/implicit integer variables) */
   CHECK_OKAY( SCIPbranchcandGetPseudoCands(branchcand, set, prob, &pseudocands, &npseudocands) );

   if( npseudocands == 0 )
   {
      /* all integer variables in the infeasible pseudo solution are fixed,
       * - if no continous variables exist, the infeasible pseudo solution is completely fixed, and the node is infeasible
       * - if at least one continous variable exist, we cannot resolve the infeasibility by branching -> solve LP
       */
      if( prob->ncont == 0 )
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTRUN;

   /* call external pseudo branching methods */
   todoMessage("external pseudo branching");
   
   /* if no pseudo branching method succeeded in choosing a branching, just branch on the first non-fixed variable */
   if( *result == SCIP_DIDNOTRUN )
   {
      CHECK_OKAY( SCIPtreeBranchVar(tree, memhdr, set, stat, lp, branchcand, eventqueue, pseudocands[0]) );
      *result = SCIP_BRANCHED;
   }

   assert(*result != SCIP_DIDNOTRUN);

   return SCIP_OKAY;
}

static
RETCODE enforceConstraints(             /**< enforces constraints by branching, separation, or domain reduction */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   PRIMAL*          primal,             /**< primal data */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   CONSHDLR**       conshdlrs_enfo,     /**< constraint handlers for enforcing constraints, sorted by priority */
   Bool*            infeasible,         /**< pointer to store whether LP/pseudo solution is infeasible */
   Bool*            solveagain          /**< pointer to store whether node has to be solved again */
   )
{
   RESULT result;
   Bool resolved;
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
   resolved = FALSE;

   todoMessage("avoid checking the same pseudosolution twice");

   /* enforce constraints by branching, applying additional cutting planes (if LP is being processed),
    * introducing new constraints, or tighten the domains
    */
   debugMessage("enforcing constraints\n");

   for( h = 0; h < set->nconshdlrs && !resolved; ++h )
   {
      if( tree->actnodehaslp )
      {
         CHECK_OKAY( SCIPconshdlrEnforceLPSol(conshdlrs_enfo[h], set, &result) );
      }
      else
      {
         CHECK_OKAY( SCIPconshdlrEnforcePseudoSol(conshdlrs_enfo[h], set, &result) );
      }
      switch( result )
      {
      case SCIP_CUTOFF:
         assert(tree->nchildren == 0);
         *infeasible = TRUE;
         resolved = TRUE;
         break;
      case SCIP_BRANCHED:
         assert(tree->nchildren >= 1);
         *infeasible = TRUE;
         resolved = TRUE;
         break;
      case SCIP_REDUCEDDOM:
         assert(tree->nchildren == 0);
         *infeasible = TRUE;
         *solveagain = TRUE;
         resolved = TRUE;
         break;
      case SCIP_SEPARATED:
         assert(tree->nchildren == 0);
         if( !tree->actnodehaslp )
         {
            char s[255];
            sprintf(s, "enforcing method of constraint handler <%s> separated cuts, but LP is not processed",
               SCIPconshdlrGetName(conshdlrs_enfo[h]));
            errorMessage(s);
            return SCIP_INVALIDRESULT;
         }
         *infeasible = TRUE;
         *solveagain = TRUE;
         resolved = TRUE;
         break;
      case SCIP_INFEASIBLE:
         assert(tree->nchildren == 0);
         *infeasible = TRUE;
         break;
      case SCIP_FEASIBLE:
         assert(tree->nchildren == 0);
         break;
      default:
         {
            char s[255];
            sprintf(s, "invalid result code <%d> from enforcing method of constraint handler <%s>",
               result, SCIPconshdlrGetName(conshdlrs_enfo[h]));
            errorMessage(s);
            return SCIP_INVALIDRESULT;
         }
         break;
      }
      assert(!(*solveagain) || (resolved && *infeasible));
   }
   debugMessage(" -> enforcing result: infeasible=%d, solveagain=%d, resolved=%d\n",
      *infeasible, *solveagain, resolved);

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
   
      CHECK_OKAY( pseudoBranch(memhdr, set, stat, prob, tree, lp, branchcand, eventqueue, &result) );

      switch( result )
      {
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
         /* we didn't resolve the infeasibility, but the node may be feasible due to contious variables
          * (this can only happen, if the LP was not solved)
          *  -> we have no other possibility than solving the LP
          */
         assert(tree->nchildren == 0);
         assert(prob->ncont > 0);
         assert(!tree->actnodehaslp);
               
         /* solve the LP in the next loop */
         tree->actnodehaslp = TRUE;
         *solveagain = TRUE;
         break;

      default:
         {
            char s[255];
            sprintf(s, "invalid result code <%d> from pseudoBranch()", result);
            errorMessage(s);
            abort();
         }
      }
   }
   assert(*infeasible || !resolved);

   return SCIP_OKAY;
}

RETCODE SCIPsolveCIP(                   /**< main solving loop */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICE*           price,              /**< pricing storage */
   SEPA*            sepa,               /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   PRIMAL*          primal,             /**< primal data */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   CONSHDLR** conshdlrs_sepa;
   CONSHDLR** conshdlrs_enfo;
   CONSHDLR** conshdlrs_chck;
   NODE* actnode;
   RESULT result;
   Bool cutoff;
   Bool infeasible;
   Bool solveagain;
   Bool propagain;
   int h;

   assert(set != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(price != NULL);
   assert(sepa != NULL);
   assert(branchcand != NULL);
   assert(cutpool != NULL);
   assert(primal != NULL);
   assert(eventqueue != NULL);

   todoMessage("every variable, where zero is not the best bound (w.r.t. objective function) has to be in the problem");

   /* sort constraint handlers by priorities */
   duplicateMemoryArray(conshdlrs_sepa, set->conshdlrs, set->nconshdlrs);
   duplicateMemoryArray(conshdlrs_enfo, set->conshdlrs, set->nconshdlrs);
   duplicateMemoryArray(conshdlrs_chck, set->conshdlrs, set->nconshdlrs);
   SCIPbsortPtr((void**)conshdlrs_sepa, set->nconshdlrs, SCIPconshdlrCompSepa);
   SCIPbsortPtr((void**)conshdlrs_enfo, set->nconshdlrs, SCIPconshdlrCompEnfo);
   SCIPbsortPtr((void**)conshdlrs_chck, set->nconshdlrs, SCIPconshdlrCompChck);

   while( stat->nnodes < set->nodelimit )
   {
      /* select next node to process */
      CHECK_OKAY( SCIPnodeselSelect(set->nodesel, set->scip, &actnode) );

      if( actnode != NULL )
      {
         /* update statistics */
         stat->nnodes++;
         stat->maxdepth = MAX(stat->maxdepth, actnode->depth);
         if( actnode->nodetype == SCIP_NODETYPE_LEAF )
         {
            stat->plungedepth = 0;
            debugMessage("selected leaf node, lowerbound=%g, plungedepth=%d\n", actnode->lowerbound, stat->plungedepth);
         }
         else if( actnode->nodetype == SCIP_NODETYPE_CHILD )
         {
            stat->plungedepth++;
            debugMessage("selected child node, lowerbound=%g, plungedepth=%d\n", actnode->lowerbound, stat->plungedepth);
         }
         else
         {
            assert(actnode->nodetype == SCIP_NODETYPE_SIBLING);
            debugMessage("selected sibling node, lowerbound=%g, plungedepth=%d\n", actnode->lowerbound, stat->plungedepth);
         }
      }
      
      /* delay events in node activation */
      CHECK_OKAY( SCIPeventqueueDelay(eventqueue) );

      /* activate selected node */
      CHECK_OKAY( SCIPnodeActivate(actnode, memhdr, set, stat, lp, tree, branchcand, eventqueue) );
      assert(tree->actnode == actnode);

      /* process the delayed events */
      CHECK_OKAY( SCIPeventqueueProcess(eventqueue, memhdr, set, tree, lp, branchcand) );

      /* if no more node was selected, we finished optimisation */
      if( actnode == NULL )
         break;

      debugMessage("Processing node %lld in depth %d, %d siblings\n", stat->nnodes, actnode->depth, tree->nsiblings);
      
      /* presolve node */
      todoMessage("node presolving");

      /* domain propagation */
      todoMessage("move domain propagation into an own method; limit number of loops");
      cutoff = FALSE;
      do
      {
         propagain = FALSE;
         for( h = 0; h < set->nconshdlrs && !cutoff; ++h )
         {
            CHECK_OKAY( SCIPconshdlrPropagate(conshdlrs_enfo[h], set, actnode->depth, &result) );
            propagain |= (result == SCIP_REDUCEDDOM);
            cutoff |= (result == SCIP_CUTOFF);
         }
      }
      while( propagain && !cutoff );

      if( cutoff )
      {
         debugMessage("node preprocessing determined that node is infeasible\n");
         tree->actnode->lowerbound = primal->upperbound;
      }
      else
      {
         debugMessage("actual pseudosolution: obj=%g", tree->actpseudoobjval);
         debug(SCIPprobPrintPseudoSol(prob, set));
         
         /* check, if we want to solve the LP at the selected node */
         tree->actnodehaslp = (set->lpsolvefreq >= 1 && actnode->depth % set->lpsolvefreq == 0);
         /* solve at least the root LP, if there are continous variables present */
         tree->actnodehaslp |= (actnode->depth == 0 && prob->ncont > 0);
         /* don't solve the node if its cut off by the pseudo objective value anyway */
         tree->actnodehaslp &= SCIPsetIsL(set, tree->actpseudoobjval, primal->upperbound);
         
         /* external node solving loop */
         do
         {
            solveagain = FALSE;
            infeasible = FALSE;
            
            /* check, if we want to solve the LP at this node */
            if( tree->actnodehaslp )
            {
               /* solve the LP */
               CHECK_OKAY( solveNodeLP(memhdr, set, stat, prob, tree, lp, price, sepa, cutpool, primal, branchcand, 
                              eventqueue, conshdlrs_sepa) );
               assert(lp->solved);
            }
            assert(SCIPsepaGetNCuts(sepa) == 0);
            
            /* update lower bound w.r.t. the pseudo solution */
            tree->actnode->lowerbound = MAX(tree->actnode->lowerbound, tree->actpseudoobjval);
            tree->actnode->lowerbound = MIN(tree->actnode->lowerbound, primal->upperbound);
            
            /* check for infeasible node by bounding */
            if( SCIPsetIsGE(set, tree->actnode->lowerbound, primal->upperbound) )
            {
               debugMessage("node is infeasible (lower=%g, upper=%g)\n", tree->actnode->lowerbound, primal->upperbound);
               infeasible = TRUE;
            }
            else
            {
               /* enforce constraints */
               CHECK_OKAY( enforceConstraints(memhdr, set, stat, prob, tree, lp, branchcand, primal, eventqueue,
                              conshdlrs_enfo, &infeasible, &solveagain) );
            }
         }
         while( solveagain );
         
         if( !infeasible )
         {
            SOL* sol;
            
            /* found a feasible solution */
            if( tree->actnodehaslp )
            {
               CHECK_OKAY( SCIPsolCreateLPSol(&sol, memhdr, stat, lp, NULL) );
            }
            else
            {
               CHECK_OKAY( SCIPsolCreatePseudoSol(&sol, memhdr, stat, tree, NULL) );
            }
            CHECK_OKAY( SCIPprimalAddSolMove(primal, memhdr, set, stat, prob, tree, lp, &sol) );
         }
      }

      /* call primal heuristics */
      todoMessage("primal heuristics");
      
      /* display node information line */
      CHECK_OKAY( SCIPdispPrintLine(set, stat, FALSE) );
      
      debugMessage("Processing of node in depth %d finished. %d siblings, %d children, %d leaves left\n", 
         SCIPnodeGetDepth(actnode), tree->nsiblings, tree->nchildren, SCIPtreeGetNLeaves(tree));
      debugMessage("**********************************************************************\n");
   }

   debugMessage("Problem solving finished\n");

   /* free sorted constraint handler arrays */
   freeMemoryArrayNull(conshdlrs_sepa);
   freeMemoryArrayNull(conshdlrs_enfo);
   freeMemoryArrayNull(conshdlrs_chck);

   return SCIP_OKAY;
}

