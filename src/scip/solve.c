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
   PROB*            transprob,          /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICE*           price               /**< pricing storage */
   )
{
   Bool nodesolved;
   Bool pricingsolved;

   /* load the LP into the solver and load the LP state */
   CHECK_OKAY( SCIPtreeLoadLP(tree, memhdr, set, lp) );
   
   /* only the bounds and the constraint list may have been changed since the LP state node
    * -> dual simplex is applicable as first solver
    */

   /* price-and-cut loop */
   do
   {
      debugMessage("node solving loop\n");

      nodesolved = TRUE;

      /* solve LP with dual simplex */
      CHECK_OKAY( solveDualLP(set, memhdr, lp) );
      assert(lp->solved);

      /* pricing */
      do
      {
         CHECK_OKAY( SCIPpriceVars(price, set, memhdr, stat, transprob, lp) );
         pricingsolved = lp->solved;

         /* solve LP with primal simplex */
         CHECK_OKAY( solvePrimalLP(set, memhdr, lp) );
         assert(lp->solved);
      }
      while( !pricingsolved );

      /* reset bounds temporarily set by pricer to their original values */
      CHECK_OKAY( SCIPpriceResetBounds(price, set, lp) );
      nodesolved &= lp->solved;

      /* solve LP with dual simplex */
      CHECK_OKAY( solveDualLP(set, memhdr, lp) );
      assert(lp->solved);

      /* constraint separation */
      todoMessage("constraint separation");
      nodesolved &= lp->solved;

      /* solve LP with dual simplex */
      CHECK_OKAY( solveDualLP(set, memhdr, lp) );
      assert(lp->solved);

      /* cut separation */
      todoMessage("cut separation");
      nodesolved &= lp->solved;

      /* solve LP with dual simplex */
      CHECK_OKAY( solveDualLP(set, memhdr, lp) );
      assert(lp->solved);
   }
   while( !nodesolved );

   return SCIP_OKAY;
}

RETCODE SCIPsolveCIP(                   /**< main solving loop */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory buffers */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            transprob,          /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICE*           price               /**< pricing storage */
   )
{
   NODE* actnode;

   assert(set != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);

   /* select next node to process */
   CHECK_OKAY( SCIPnodeselSelect(set->nodesel, set->scip, &actnode) );

   while( actnode != NULL )
   {
      debugMessage("Processing node in depth %d (type %d)\n", SCIPnodeGetDepth(actnode), SCIPnodeGetType(actnode));

      /* activate selected node */
      CHECK_OKAY( SCIPnodeActivate(actnode, memhdr, set, lp, tree) );
      debugMessage("%d siblings, %d children, %d leaves left\n", 
         SCIPtreeGetNSiblings(tree), SCIPtreeGetNChildren(tree), SCIPtreeGetNLeaves(tree));

      /* presolve node */
      todoMessage("node presolving");

      /* solve node */
      CHECK_OKAY( solveNodeLP(set, memhdr, stat, transprob, tree, lp, price) );
      
      /* branch, if necessary */
      todoMessage("branching");

      /* select next node to process */
      CHECK_OKAY( SCIPnodeselSelect(set->nodesel, set->scip, &actnode) );
   }

   return SCIP_OKAY;
}

