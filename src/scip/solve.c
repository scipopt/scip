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
#include "solve.h"


static
RETCODE solveDualLP(                    /**< solves the LP with dual simplex algorithm */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory buffers */
   LP*              lp                  /**< LP data */
   )
{
   assert(lp != NULL);

   if( !lp->solved )
   {
      LPSOLSTAT lpsolstat;
      
      CHECK_OKAY( SCIPlpSolveDual(lp, set, memhdr, &lpsolstat) );
      switch( lpsolstat )
      {
      case SCIP_LPSOLSTAT_OPTIMAL:
         CHECK_OKAY( SCIPlpGetSol(lp, set, memhdr) );
         break;
      case SCIP_LPSOLSTAT_INFEASIBLE:
         todoMessage("LP infeasibility processing");
         return SCIP_ERROR;
      case SCIP_LPSOLSTAT_UNBOUNDED:
         todoMessage("LP unboundness processing");
         return SCIP_ERROR;
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
RETCODE solvePrimalLP(                  /**< solves the LP with primal simplex algorithm */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory buffers */
   LP*              lp                  /**< LP data */
   )
{
   assert(lp != NULL);

   if( !lp->solved )
   {
      LPSOLSTAT lpsolstat;
      
      CHECK_OKAY( SCIPlpSolvePrimal(lp, set, memhdr, &lpsolstat) );
      switch( lpsolstat )
      {
      case SCIP_LPSOLSTAT_OPTIMAL:
         CHECK_OKAY( SCIPlpGetSol(lp, set, memhdr) );
         break;
      case SCIP_LPSOLSTAT_INFEASIBLE:
         todoMessage("LP infeasibility processing");
         return SCIP_ERROR;
      case SCIP_LPSOLSTAT_UNBOUNDED:
         todoMessage("LP unboundness processing");
         return SCIP_ERROR;
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
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory buffers */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICE*           price,              /**< pricing storage */
   SEPA*            sepa,               /**< separation storage */
   CONSHDLR**       conshdlrs_sepa      /**< constraint handlers sorted by separation priority */
   )
{
   Bool mustprice;
   Bool mustsepar;
   int h;

   assert(set != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(lp != NULL);
   assert(price != NULL);
   assert(set->nconshdlrs == 0 || conshdlrs_sepa != NULL);

   /* load the LP into the solver and load the LP state */
   debugMessage("loading LP\n");
   CHECK_OKAY( SCIPtreeLoadLP(tree, memhdr, set, lp) );
   
   /* only the bounds and the constraint list may have been changed since the LP state node
    * -> dual simplex is applicable as first solver
    */
   debugMessage("dual solve\n");
   CHECK_OKAY( solveDualLP(set, memhdr, lp) );
   assert(lp->solved);

   /* price-and-cut loop */
   mustprice = TRUE;
   mustsepar = TRUE;
   while( mustprice || mustsepar )
   {
      debugMessage("-------- node solving loop --------\n");

      /* pricing (has to be done completely to get a valid lower bound) */
      while( mustprice )
      {
         assert(lp->solved);

         debugMessage("pricing\n");
         CHECK_OKAY( SCIPpriceVars(price, set, memhdr, stat, prob, lp) );
         mustprice = !lp->solved;
         mustsepar |= !lp->solved;
            
         /* solve LP with primal simplex */
         debugMessage("pricing: primal solve\n");
         CHECK_OKAY( solvePrimalLP(set, memhdr, lp) );
         assert(lp->solved);
         
         /* reset bounds temporarily set by pricer to their original values */
         debugMessage("reset bounds\n");
         CHECK_OKAY( SCIPpriceResetBounds(price, memhdr, set, lp) );
         mustprice = !lp->solved;
         mustsepar |= !lp->solved;

         /* solve LP with dual simplex */
         debugMessage("reset bounds: dual solve\n");
         CHECK_OKAY( solveDualLP(set, memhdr, lp) );
         assert(lp->solved);
      }

      /* separation (has not to be done completely, because we just want to increase the lower bound) */
      if( mustsepar )
      {
         assert(lp->solved);

         /* constraint separation */
         debugMessage("constraint separation\n");
         assert(SCIPsepaGetNCuts(sepa) == 0);
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
         
         /* solve LP with dual simplex */
         debugMessage("separation: dual solve\n");
         CHECK_OKAY( solveDualLP(set, memhdr, lp) );
         assert(lp->solved);
      }
   }

   /* remember that this node is solved correctly */
   tree->correctlpdepth = tree->actnode->depth;
   /* tree->pathnlpcols[tree->actnode->depth] = lp->ncols;  ???
    * tree->pathnlprows[tree->actnode->depth] = lp->nrows;  ???
    */
   return SCIP_OKAY;
}

RETCODE SCIPsolveCIP(                   /**< main solving loop */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory buffers */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICE*           price,              /**< pricing storage */
   SEPA*            sepa                /**< separation storage */
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
         debugMessage("Processing node in depth %d\n", SCIPnodeGetDepth(actnode));
         
         /* presolve node */
         todoMessage("node presolving + domain propagation");
         
         /* solve node */
         todoMessage("decide whether to solve the LP or not");
         if( TRUE )
         {
            CHECK_OKAY( solveNodeLP(set, memhdr, stat, prob, tree, lp, price, sepa, conshdlrs_sepa) );
         }
         else
         {
            CHECK_OKAY( SCIPnodeReleaseLPIState(tree->actlpfork, memhdr, lp) );
         }

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
            /* found a feasible solution */
            assert(!resolved);
            todoMessage("feasible solution found");
         }

         /* call primal heuristics */
         todoMessage("primal heuristics");
         
         debugMessage("Processing of node in depth %d finished. %d siblings, %d children, %d leaves left\n", 
            SCIPnodeGetDepth(actnode), tree->nsiblings, tree->nchildren, SCIPtreeGetNLeaves(tree));
      }
   }
   while( actnode != NULL );

   debugMessage("Problem solving finished\n");

   /* free sorted constraint handler arrays */
   freeMemoryArrayNull(conshdlrs_sepa);
   freeMemoryArrayNull(conshdlrs_enfo);
   freeMemoryArrayNull(conshdlrs_chck);

   return SCIP_OKAY;
}

