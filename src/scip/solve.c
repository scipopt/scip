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
         CHECK_OKAY( SCIPpriceVars(price, memhdr, set, stat, prob, lp, tree) );
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
         CHECK_OKAY( SCIPpriceResetBounds(price, memhdr, set, stat, lp, tree) );
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

RETCODE SCIPsolveCIP(                   /**< main solving loop */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
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
   RESULT result;
   Bool infeasible;
   Bool resolved;
   Bool solveagain;
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
         if( actnode->nodetype == SCIP_NODETYPE_CHILD )
         {
            stat->plungedepth++;
            debugMessage("selected child node, lowerbound=%g, plungedepth=%d\n", actnode->lowerbound, stat->plungedepth);
         }
         else if( actnode->nodetype == SCIP_NODETYPE_SIBLING )
         {
            debugMessage("selected sibling node, lowerbound=%g, plungedepth=%d\n", actnode->lowerbound, stat->plungedepth);
         }
         else
         {
            assert(actnode->nodetype == SCIP_NODETYPE_LEAF);            
            stat->plungedepth = 0;
            debugMessage("selected leaf node, lowerbound=%g, plungedepth=%d\n", actnode->lowerbound, stat->plungedepth);
         }
      }
      
      /* activate selected node */
      CHECK_OKAY( SCIPnodeActivate(actnode, memhdr, set, stat, lp, tree) );

      /* if no more node was selected, we finished optimisation */
      if( actnode == NULL )
         break;

      debugMessage("Processing node %d in depth %d\n", stat->nnodes, actnode->depth);
      
      /* presolve node */
      todoMessage("node presolving + domain propagation");
      
      debugMessage("actual pseudosolution: obj=%g", tree->actpseudoobjval);
      debug(SCIPprobPrintPseudoSol(prob, set));
      
      /* check, if we want to solve the LP at the selected node */
      tree->actnodehaslp = (set->lpsolvefreq >= 1 && actnode->depth % set->lpsolvefreq == 0);
      /* solve at least the root LP, if there are continous variables present */
      tree->actnodehaslp |= (actnode->depth == 0 && prob->ncont > 0);

      /* external node solving loop */
      do
      {
         /* check, if we want to solve the LP at this node */
         if( tree->actnodehaslp )
         {
            /* solve the LP */
            CHECK_OKAY( solveNodeLP(memhdr, set, stat, prob, tree, lp, price, sepa, cutpool, primal, conshdlrs_sepa) );
            assert(lp->solved);
         }

         /* update lower bound w.r.t. the pseudo solution */
         tree->actnode->lowerbound = MAX(tree->actnode->lowerbound, tree->actpseudoobjval);
         tree->actnode->lowerbound = MIN(tree->actnode->lowerbound, primal->upperbound);
         
         /* check for infeasible node */
         if( SCIPsetIsGE(set, actnode->lowerbound, primal->upperbound) )
         {
            debugMessage("node is infeasible (lower=%g, upper=%g)\n", actnode->lowerbound, primal->upperbound);
            infeasible = TRUE;
            resolved = TRUE;
            solveagain = FALSE;
         }
         else
         {
            todoMessage("avoid checking the same pseudosolution twice");

            /* enforce constraints by branching, applying additional cutting planes (if LP is being processed),
             * introducing new constraints, or tighten the domains
             */
            assert(SCIPsepaGetNCuts(sepa) == 0);

            debugMessage("enforcing constraints\n");

            infeasible = FALSE;
            resolved = FALSE;
            solveagain = FALSE;
            for( h = 0; h < set->nconshdlrs && !resolved; ++h )
            {
               CHECK_OKAY( SCIPconshdlrEnforce(conshdlrs_enfo[h], set, tree->actnodehaslp, &result) );
               switch( result )
               {
               case SCIP_BRANCHED:
                  infeasible = TRUE;
                  resolved = TRUE;
                  break;
               case SCIP_REDUCEDDOM:
                  infeasible = TRUE;
                  resolved = TRUE;
                  solveagain = TRUE;
                  break;
               case SCIP_SEPARATED:
                  if( !tree->actnodehaslp )
                  {
                     char s[255];
                     sprintf(s, "enforcing method of constraint handler <%s> separated cuts, but LP is not processed",
                        SCIPconshdlrGetName(conshdlrs_enfo[h]));
                     errorMessage(s);
                     return SCIP_INVALIDRESULT;
                  }
                  infeasible = TRUE;
                  resolved = TRUE;
                  solveagain = TRUE;
                  break;
               case SCIP_INFEASIBLE:
                  infeasible = TRUE;
                  break;
               case SCIP_FEASIBLE:
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
               assert(!solveagain || (resolved && infeasible));
            }
            debugMessage(" -> enforcing result: infeasible=%d, resolved=%d, solveagain=%d\n",
               infeasible, resolved, solveagain);
         }

         if( infeasible && !resolved )
         {
            VAR* var;
            Real solval;
            int v;
            
            assert(!solveagain);

            /* the node is infeasible, but no constraint handler could resolve the infeasibility
             * -> select non-fixed binary or integer variable x with value x', create three sons: 
             *    x <= x'-1, x = x', and x >= x'+1.
             *    In the left and right branch, the actual solution is cutted off. In the middle
             *    branch, the constraints can hopefully reduce domains of other variables to cut
             *    off the actual solution.
             */
            todoMessage("select variable fixing");
            
            /* the following is just a temporary implementation of variable fixing */
            
            /* scan the integer variables for a non-fixed variable */
            for( v = 0; v < prob->nbin + prob->nint + prob->nimpl && !resolved; ++v )
            {
               var = prob->vars[v];
               /* check, if variable is already fixed */
               if( !SCIPsetIsEQ(set, var->dom.lb, var->dom.ub) )
               {
                  NODE* node;
                  Real fixval;
                  
                  /* create child nodes with x <= x'-1, x = x', and x >= x'+1 */
                  solval = SCIPvarGetSol(var, tree);
                  fixval = SCIPsetCeil(set, solval);
                  assert(SCIPsetIsEQ(set, SCIPsetCeil(set, solval), SCIPsetFloor(set, solval)));
                  assert(SCIPsetIsIntegral(set, solval));
                  assert(SCIPsetIsIntegral(set, var->dom.lb));
                  assert(SCIPsetIsIntegral(set, var->dom.ub));
                  assert(SCIPsetIsGE(set, solval, var->dom.lb));
                  assert(SCIPsetIsLE(set, solval, var->dom.ub));
                  
                  debugMessage("selected variable <%s> with value %g to fix to %g\n", var->name, solval, fixval);
                  
                  /* create child node with x = x' */
                  debugMessage(" -> creating child: <%s> == %g\n", var->name, SCIPsetCeil(set, solval));
                  /*CHECK_OKAY( SCIPcreateChild(scip, &node) );
                    CHECK_OKAY( SCIPchgNodeUb(scip, node, var, SCIPceil(scip, primsol)) );
                    CHECK_OKAY( SCIPchgNodeLb(scip, node, var, SCIPfloor(scip, primsol)) );*/
                  CHECK_OKAY( SCIPnodeCreate(&node, memhdr, set, tree) );
                  if( !SCIPsetIsEQ(set, var->dom.lb, fixval) )
                  {
                     CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, lp, tree, var, fixval,
                                    SCIP_BOUNDTYPE_LOWER) );
                  }
                  if( !SCIPsetIsEQ(set, var->dom.ub, fixval) )
                  {
                     CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, lp, tree, var, fixval, 
                                    SCIP_BOUNDTYPE_UPPER) );
                  }
                  
                  /* create child node with x <= x'-1, if this would be feasible */
                  if( SCIPsetIsGE(set, fixval-1, var->dom.lb) )
                  {
                     debugMessage(" -> creating child: <%s> <= %g\n", var->name, fixval-1);
                     /*CHECK_OKAY( SCIPcreateChild(scip, &node) );
                       CHECK_OKAY( SCIPchgNodeUb(scip, node, var, SCIPceil(scip, primsol)-1) );*/
                     CHECK_OKAY( SCIPnodeCreate(&node, memhdr, set, tree) );
                     CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, lp, tree, var, fixval-1,
                                    SCIP_BOUNDTYPE_UPPER) );
                  }
                  
                  /* create child node with x >= x'+1, if this would be feasible */
                  if( SCIPsetIsLE(set, fixval+1, var->dom.ub) )
                  {
                     debugMessage(" -> creating child: <%s> >= %g\n", var->name, fixval+1);
                     /*CHECK_OKAY( SCIPcreateChild(scip, &node) );
                       CHECK_OKAY( SCIPchgNodeLb(scip, node, var, SCIPfloor(scip, primsol)+1) );*/
                     CHECK_OKAY( SCIPnodeCreate(&node, memhdr, set, tree) );
                     CHECK_OKAY( SCIPnodeAddBoundchg(node, memhdr, set, stat, lp, tree, var, fixval+1,
                                    SCIP_BOUNDTYPE_LOWER) );
                  }
                  
                  resolved = TRUE;
               }
            }
            
            if( !resolved && prob->ncont == 0 )
            {
               /* node is infeasible, because all variables are fixed, and the pseudo solution (the only point in the
                * given bounds) is infeasible
                */
               resolved = TRUE;
            }
            
            if( !resolved )
            {
               /* we didn't resolve the infeasibility, but the node may be feasible due to contious variables
                * (this can only happen, if the LP was not solved)
                *  -> we have no other possibility than solving the LP
                */
               assert(prob->ncont > 0);
               assert(!tree->actnodehaslp);
               
               /* solve the LP in the next loop */
               tree->actnodehaslp = TRUE;
               solveagain = TRUE;
            }
         }
      }
      while( solveagain );
      
      if( !infeasible )
      {
         SOL* sol;
         
         /* found a feasible solution */
         assert(!resolved);
         if( tree->actnodehaslp )
         {
            CHECK_OKAY( SCIPsolCreateLPSol(&sol, memhdr, set, stat, lp) );
         }
         else
         {
            CHECK_OKAY( SCIPsolCreatePseudoSol(&sol, memhdr, set, stat, prob) );
         }
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

   debugMessage("Problem solving finished\n");

   /* free sorted constraint handler arrays */
   freeMemoryArrayNull(conshdlrs_sepa);
   freeMemoryArrayNull(conshdlrs_enfo);
   freeMemoryArrayNull(conshdlrs_chck);

   return SCIP_OKAY;
}

