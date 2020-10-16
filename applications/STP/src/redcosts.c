/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   redcosts.c
 * @brief  Reduced cost based routines for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file encompasses various routines for reduced costs based computations
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


//#define SCIP_DEBUG
#include "redcosts.h"
#include "portab.h"


/**@name Local methods
 *
 * @{
 */



/**@} */

/**@name Interface methods
 *
 * @{
 */


/** reduced costs available? */
SCIP_Bool redcosts_forLPareAvailable(
   SCIP*                 scip                /**< SCIP structure */
)
{
   /* only execute if current node has an LP */
   if( !SCIPhasCurrentNodeLP(scip) )
   {
      SCIPdebugMessage("!SCIPhasCurrentNodeLP \n");
      return FALSE;
   }

   /* only execute dualcostVarfixing if optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL (%d) \n", SCIPgetLPSolstat(scip));
      return FALSE;
   }

   /* only execute if current LP is valid relaxation */
   if( !SCIPisLPRelax(scip) )
   {
      SCIPdebugMessage("!SCIPisLPRelax \n");
      return FALSE;
   }

   /* we cannot apply reduced cost strengthening if no simplex basis is available */
   if( !SCIPisLPSolBasic(scip) )
   {
      SCIPdebugMessage("!SCIPisLPSolBasic \n");
      return FALSE;
   }

   /* reduced cost strengthening can only be applied if cutoff is finite */
   if( SCIPisInfinity(scip, SCIPgetCutoffbound(scip)) )
   {
      SCIPdebugMessage("!SCIPgetCutoffbound \n");
      return FALSE;
   }

   SCIPdebugMessage("reduced costs are available! \n");


   return TRUE;
}


/** initialize reduced costs*/
void redcosts_forLPget(
   SCIP*                 scip,               /**< SCIP structure */
   SCIP_VAR**            vars,               /**< variables (in) */
   const GRAPH*          graph,              /**< graph data */
   SCIP_Real*            redcosts            /**< reduced costs (out) */
   )
{
   const int nedges = graph_get_nEdges(graph);

   assert(nedges >= 0);
   assert(vars && redcosts && scip);

   for( int e = 0; e < nedges; e++ )
   {
      assert(SCIPvarIsBinary(vars[e]));

      /* variable is already fixed, we must not trust the reduced cost */
      if( SCIPvarGetLbLocal(vars[e]) + 0.5 > SCIPvarGetUbLocal(vars[e]) )
      {
         if( SCIPvarGetLbLocal(vars[e]) > 0.5 )
            redcosts[e] = 0.0;
         else
         {
            assert(SCIPvarGetUbLocal(vars[e]) < 0.5);
            redcosts[e] = FARAWAY;
         }
      }
      else
      {
         if( SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, vars[e])) )
         {
            assert(!SCIPisDualfeasNegative(scip, SCIPgetVarRedcost(scip, vars[e])));
            redcosts[e] = SCIPgetVarRedcost(scip, vars[e]);
         }
         else
         {
            assert(!SCIPisDualfeasPositive(scip, SCIPgetVarRedcost(scip, vars[e])));
            assert(SCIPisFeasEQ(scip, SCIPgetSolVal(scip, NULL, vars[e]), 1.0) || SCIPisDualfeasZero(scip, SCIPgetVarRedcost(scip, vars[e])));
            redcosts[e] = 0.0;
         }
      }

      if( redcosts[e] < 0.0 )
         redcosts[e] = 0.0;
   }
#ifdef SCIP_DISABLED_CODE
   if( graph_pc_isPcMw(graph) )
   {
      /* we do some clean-up */
      const int nnodes = graph_get_nNodes(graph);
      const int root = graph->source;

      assert(graph->term2edge);

      for( int i = 0; i < nnodes; i++ )
      {
         if( graph_pc_knotIsDummyTerm(graph, i) && i != root )
         {
            const int edge2dummy = flipedge(graph->term2edge[i]);

            assert(edge2dummy >= 0 && graph->head[edge2dummy] == i);
            assert(EQ(graph->cost[edge2dummy], 0.0));
            assert(EQ(redcosts[edge2dummy], 0.0) || EQ(redcosts[edge2dummy], FARAWAY));

            if( EQ(redcosts[edge2dummy], 0.0) )
               redcosts[edge2dummy] = 0.0;
         }
      }
   }
#endif
}
