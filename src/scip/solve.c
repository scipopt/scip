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
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory buffers */
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
         assert(var->inprob);
         
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
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory buffers */
   LP*              lp,                 /**< LP data */
   STAT*            stat                /**< problem statistics */
   )
{
   assert(lp != NULL);

   if( !lp->solved )
   {
      if( lp->dualfeasible || !lp->primalfeasible )
      {
         debugMessage("solving dual LP\n");
         CHECK_OKAY( SCIPlpSolveDual(lp, set, memhdr, stat) );
      }
      else
      {
         debugMessage("solving primal LP\n");
         CHECK_OKAY( SCIPlpSolvePrimal(lp, set, memhdr, stat) );
      }
      switch( SCIPlpGetSolstat(lp) )
      {
      case SCIP_LPSOLSTAT_OPTIMAL:
         CHECK_OKAY( SCIPlpGetSol(lp, set, memhdr, stat) );
         return SCIP_FEASIBLE;
      case SCIP_LPSOLSTAT_INFEASIBLE:
         if( set->usepricing )
         {
            CHECK_OKAY( SCIPlpGetDualfarkas(lp, set, memhdr) );
         }
         return SCIP_INFEASIBLE;
      case SCIP_LPSOLSTAT_UNBOUNDED:
         CHECK_OKAY( SCIPlpGetUnboundedSol(lp, set, memhdr, stat) );
         return SCIP_UNBOUNDED;
      case SCIP_LPSOLSTAT_OBJLIMIT:
         if( set->usepricing )
         {
            CHECK_OKAY( SCIPlpGetSol(lp, set, memhdr, stat) );
         }
         return SCIP_INFEASIBLE;         
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
   else
   {
      switch( SCIPlpGetSolstat(lp) )
      {
      case SCIP_LPSOLSTAT_OPTIMAL:
         return SCIP_FEASIBLE;
      case SCIP_LPSOLSTAT_INFEASIBLE:
         return SCIP_INFEASIBLE;
      case SCIP_LPSOLSTAT_UNBOUNDED:
         return SCIP_UNBOUNDED;
      case SCIP_LPSOLSTAT_OBJLIMIT:
         return SCIP_INFEASIBLE;         
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
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory buffers */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICE*           price,              /**< pricing storage */
   SEPA*            sepa,               /**< separation storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   PRIMAL*          primal,             /**< primal data */
   CONSHDLR**       conshdlrs_sepa      /**< constraint handlers sorted by separation priority */
   )
{
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
      CHECK_OKAY( initRootLP(set, memhdr, stat, prob, lp, tree) );
   }

   /* only the bounds and the constraint list may have been changed since the LP state node
    * -> dual simplex is applicable as first solver
    */
   debugMessage("node: solve LP\n");
   CHECK_OKAY( solveLP(set, memhdr, lp, stat) );
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
         CHECK_OKAY( SCIPpriceVars(price, set, memhdr, stat, prob, lp, tree) );
         mustprice = !lp->solved;
         mustsepar |= !lp->solved;
         
         /* after adding columns, the LP should be primal feasible such that primal simplex is applicable;
          * if LP was infeasible, we have to use dual simplex
          */
         debugMessage("pricing: solve LP\n");
         CHECK_OKAY( solveLP(set, memhdr, lp, stat) );
         assert(lp->solved);

         /* reset bounds temporarily set by pricer to their original values */
         debugMessage("reset bounds\n");
         CHECK_OKAY( SCIPpriceResetBounds(price, memhdr, set, lp, tree) );
         mustprice |= !lp->solved;
         mustsepar |= !lp->solved;

         /* solve LP (with dual simplex) */
         debugMessage("reset bounds: solve LP\n");
         CHECK_OKAY( solveLP(set, memhdr, lp, stat) );
         assert(lp->solved);

         /* if the LP is unbounded, we can stop pricing */
         mustprice &= (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDED);
      }
      assert(lp->solved);

      /* update lower bound */
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
            CHECK_OKAY( SCIPconshdlrSeparate(conshdlrs_sepa[h], set) );
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
         CHECK_OKAY( solveLP(set, memhdr, lp, stat) );
         assert(lp->solved);
      }
   }

   /* remember that this node is solved correctly */
   tree->correctlpdepth = tree->actnode->depth;
   tree->actnode->lowerbound = SCIPlpGetObjval(lp);
   tree->actnode->lowerbound = MIN(tree->actnode->lowerbound, primal->upperbound);

   return SCIP_OKAY;
}

RETCODE SCIPsolveCIP(                   /**< main solving loop */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory buffers */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICE*           price,              /**< pricing storage */
   SEPA*            sepa,               /**< separation storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   PRIMAL*          primal              /**< primal data */
   )
{
   CONSHDLR** conshdlrs_sepa;
   CONSHDLR** conshdlrs_enfo;
   CONSHDLR** conshdlrs_chck;
   NODE* actnode;
   RETCODE retcode;
   Bool feasible;
   Bool resolved;
   int h;

   assert(set != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(price != NULL);
   assert(sepa != NULL);
   assert(cutpool != NULL);
   assert(primal != NULL);

   /* sort constraint handlers by priorities */
   duplicateMemoryArray(conshdlrs_sepa, set->conshdlrs, set->nconshdlrs);
   duplicateMemoryArray(conshdlrs_enfo, set->conshdlrs, set->nconshdlrs);
   duplicateMemoryArray(conshdlrs_chck, set->conshdlrs, set->nconshdlrs);
   SCIPbsortPtr((void**)conshdlrs_sepa, set->nconshdlrs, SCIPconshdlrCompSepa);
   SCIPbsortPtr((void**)conshdlrs_enfo, set->nconshdlrs, SCIPconshdlrCompEnfo);
   SCIPbsortPtr((void**)conshdlrs_chck, set->nconshdlrs, SCIPconshdlrCompChck);

   actnode = NULL;
   do
   {
      /* select next node to process */
      CHECK_OKAY( SCIPnodeselSelect(set->nodesel, set->scip, &actnode) );

      /* activate selected node */
      CHECK_OKAY( SCIPnodeActivate(actnode, memhdr, set, lp, tree) );

      if( actnode != NULL )
      {
         stat->nnodes++;
         stat->maxdepth = MAX(stat->maxdepth, actnode->depth);

         debugMessage("Processing node %d in depth %d\n", stat->nnodes, SCIPnodeGetDepth(actnode));
         
         /* presolve node */
         todoMessage("node presolving + domain propagation");
         
#ifndef NDEBUG
         debugMessage("actual pseudosolution: ");
         debug( SCIPsolPrint(tree->actpseudosol, set, NULL) );
#endif

         /* solve node */
         todoMessage("decide whether to solve the LP or not");
#if 0
         if( ... )
         {
            CHECK_OKAY( solveNodeLP(set, memhdr, stat, prob, tree, lp, price, sepa, cutpool, primal, conshdlrs_sepa) );
         }
         else
         {
            CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
         }
#else
         CHECK_OKAY( solveNodeLP(set, memhdr, stat, prob, tree, lp, price, sepa, cutpool, primal, conshdlrs_sepa) );
#endif

         /* check for infeasible node */
         if( SCIPsetIsGE(set, actnode->lowerbound, primal->upperbound) )
         {
            debugMessage("node is infeasible (lower=%g, upper=%g)\n", actnode->lowerbound, primal->upperbound);
            feasible = FALSE;
            resolved = TRUE;
         }
         else
         {
            /* enforce constraints by branching, applying additional cutting planes, or introducing new constraints */
            assert(SCIPsepaGetNCuts(sepa) == 0);
            feasible = TRUE;
            resolved = FALSE;
            for( h = 0; h < set->nconshdlrs && !resolved; ++h )
            {
               CHECK_OKAY( retcode = SCIPconshdlrEnforce(conshdlrs_enfo[h], set) );
               switch( retcode )
               {
               case SCIP_BRANCHED:
               case SCIP_REDUCEDDOM:
               case SCIP_SEPARATED:
                  resolved = TRUE;
                  feasible = FALSE;
                  break;
               case SCIP_INFEASIBLE:
                  feasible = FALSE;
                  break;
               case SCIP_FEASIBLE:
                  break;
               default:
                  {
                     char s[255];
                     sprintf(s, "invalid return code <%d> from enforcing method of constraint handler <%s>",
                        retcode, SCIPconshdlrGetName(conshdlrs_enfo[h]));
                     errorMessage(s);
                  }
                  break;
               }
            }
         }

         if( !feasible && !resolved )
         {
            /* the node is infeasible, but no constraint handler could resolve the infeasibility
             * -> select non-fixed binary or integer variable x with value x', create three sons: 
             *    x <= x'-1, x = x', and x >= x'+1.
             *    In the left and right branch, the actual solution is cutted off. In the middle
             *    branch, the constraints can hopefully reduce domains of other variables to cut
             *    off the actual solution.
             */
            todoMessage("standard resolve for infeasible nodes -> 3-way branching on single variable");
         }

         if( feasible )
         {
            SOL* sol;

            /* found a feasible solution */
            assert(!resolved);
            CHECK_OKAY( SCIPsolCreateLPSol(&sol, memhdr, set, stat, lp) );
            CHECK_OKAY( SCIPprimalAddSol(primal, memhdr, set, stat, tree, lp, &sol) );
            CHECK_OKAY( SCIPsolRelease(&sol, memhdr, set, lp) );
         }

         /* call primal heuristics */
         todoMessage("primal heuristics");

         /* display node information line */
         CHECK_OKAY( SCIPdispPrintLine(set, stat, FALSE) );

         debugMessage("Processing of node in depth %d finished. %d siblings, %d children, %d leaves left\n", 
            SCIPnodeGetDepth(actnode), tree->nsiblings, tree->nchildren, SCIPtreeGetNLeaves(tree));
         debugMessage("**********************************************************************\n");
      }
   }
   while( actnode != NULL && stat->nnodes < set->nodelimit );

   debugMessage("Problem solving finished\n");

   /* free sorted constraint handler arrays */
   freeMemoryArrayNull(conshdlrs_sepa);
   freeMemoryArrayNull(conshdlrs_enfo);
   freeMemoryArrayNull(conshdlrs_chck);

   return SCIP_OKAY;
}

