/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   HeurFrats.cpp
 * @brief  fractional travelling salesman heuristic - Rounding heuristic for TSP
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "HeurFrats.h"
#include "ProbDataTSP.h"

using namespace tsp;
using namespace std;


/*
 * Local methods
 */


/*
 * Callback methods of primal heuristic
 */


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_DECL_HEURFREE(HeurFrats::scip_free)
{
   return SCIP_OKAY;
}
   
/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_DECL_HEURINIT(HeurFrats::scip_init)
{
   ProbDataTSP*   probdata;  

   /* create heuristic data */
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
   
   /* load the problem specific data */
   probdata = dynamic_cast<ProbDataTSP*>(SCIPgetObjProbData(scip));
   assert(probdata != NULL);

   graph = probdata->getGraph();
   assert(graph != NULL);

   capture_graph(graph);

   return SCIP_OKAY;
}
   
/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_DECL_HEUREXIT(HeurFrats::scip_exit)
{
   /* free everything which was created in scip_init */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   release_graph(&graph);

   return SCIP_OKAY;
}
   
/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 */
SCIP_DECL_HEURINITSOL(HeurFrats::scip_initsol)
{
   return SCIP_OKAY;
}
   
/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_DECL_HEUREXITSOL(HeurFrats::scip_exitsol)
{
   return SCIP_OKAY;
}
   
/** execution method of primal heuristic */
SCIP_DECL_HEUREXEC(HeurFrats::scip_exec)
{  /*lint --e{715}*/

   SCIP_SOL* newsol;
   GRAPHNODE* currnode;   
   SCIP_Bool* visited;   
   int nnodes;
   int i;
   SCIP_Bool success;

   assert(result != NULL);
   /* since the timing is SCIP_HEURTIMING_AFTERLPNODE, the current node should have an LP */
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* get the working solution from heuristic's local data */
   assert(sol != NULL);

   /* copy the current LP solution to the working solution */
   SCIP_CALL( SCIPlinkLPSol(scip, sol) );

   *result = SCIP_DIDNOTFIND;

   /* choose the first node as starting point*/
   currnode = &graph->nodes[0];   
   nnodes = graph->nnodes;
   success = TRUE;
   
   /* allocate local memory */
   SCIP_CALL( SCIPcreateSol (scip, &newsol, heur) );      
   SCIP_CALL( SCIPallocBufferArray(scip, &visited, nnodes) ); 
   BMSclearMemoryArray(visited, nnodes);
   
   assert(currnode->id == 0);
   visited[0] = TRUE;

   /*exactly nnodes edges have to be inserted into the tour */
   for( i = 0; i < nnodes; i++ )
   {
      GRAPHEDGE* edge;
      SCIP_Real bestval; 
      GRAPHEDGE* bestedge;

      /* initialization */
      bestedge = NULL;
      bestval = -1;

      /* the graph works with adjacency lists */
      edge = currnode->first_edge; 
      
      /* the last edge is treated separately */
      if( i != nnodes-1 )
      {
         while( edge != NULL )
         {
            /* update, if an edge has a better LP value AND was not visited yet AND was not globally fixed to zero */
            if( SCIPgetSolVal(scip, sol, edge->var) > bestval && !visited[edge->adjac->id] 
               && SCIPvarGetUbGlobal(edge->var) == 1.0 )
            {
               bestval = SCIPgetSolVal(scip, sol, edge->var);
               bestedge = edge;
            }
            edge = edge->next;
         }
      }
      else
      {
         GRAPHNODE* finalnode;
         finalnode = &graph->nodes[0]; 

         /* find the last edge which closes the tour */
         while( edge != NULL )
         {
            if( edge->adjac == finalnode )
            {
               if( SCIPvarGetUbGlobal(edge->var) == 1.0 )
               {
                  bestval =  SCIPgetSolVal(scip, sol, edge->var);
                  bestedge = edge;
               }
               break;
            }
            edge = edge->next;
         }
      }

      /* it may happen that we were not able to build a complete tour */
      if( bestval == -1 )
      {
         success = FALSE;
         break;
      }
      /* assert that the data is not corrupted */
      assert(bestedge != NULL);
      assert(SCIPisFeasLE(scip, 0.0, bestval) && SCIPisFeasLE(scip, bestval, 1.0));
      assert(bestval == SCIPgetSolVal(scip, sol, bestedge->var));

      /* fix the variable which represents the best edge to one in the new solution and proceed to next node */
      SCIP_CALL( SCIPsetSolVal(scip, newsol, bestedge->var, 1.0) );
      currnode = bestedge->adjac;
      assert(currnode != NULL);
      assert(0 <= currnode->id && currnode->id <= nnodes-1);
      if( i != nnodes-1 )
         assert(!visited[currnode->id]);
      visited[currnode->id] = TRUE;
   }
   /* if we were able to construct a complete tour, try to add the solution to SCIP */
   if( success )
   {
      for( i = 0; i < nnodes; i++ )
         assert(visited[graph->nodes[i].id]);
      
      success = FALSE;
      /* due to construction we already know, that the solution will be feasible */
      SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, FALSE, FALSE, FALSE, &success) );
      if( success )
         *result = SCIP_FOUNDSOL;  
   }
   /* free all local memory */
   SCIP_CALL( SCIPfreeSol(scip, &newsol) );      
   SCIPfreeBufferArray(scip, &visited);
   
   return SCIP_OKAY;
}

/** clone method which will be used to copy a objective plugin */
SCIP_DECL_HEURCLONE(scip::ObjCloneable* HeurFrats::clone)
{
   return new HeurFrats(scip);
}
