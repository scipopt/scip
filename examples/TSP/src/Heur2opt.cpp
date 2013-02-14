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
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   Heur2opt.cpp
 * @brief  2-Optimum - combinatorial improvement heuristic for TSP
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>
#include <iostream>

#include "objscip/objscip.h"
#include "GomoryHuTree.h"
#include "Heur2opt.h"
#include "ProbDataTSP.h"

using namespace tsp;
using namespace std;


/** method finding the edge going from the node with id index1 to the node with id index2 */
static
GRAPHEDGE* findEdge(
   GRAPHNODE*         nodes,              /**< all nodes of the graph */
   GRAPHNODE*         node1,              /**< id of the node where the searched edge starts */
   GRAPHNODE*         node2               /**< id of the node where the searched edge ends */
   )
{
   GRAPHEDGE* edge =  node1->first_edge;

   // regard every outgoing edge of node index1 and stop if adjacent to node index2
   while( edge != NULL )
   {
      if( edge->adjac == node2 )
         break;
      edge = edge->next;
   }
   return edge;   
}


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_DECL_HEURFREE(Heur2opt::scip_free)
{
   return SCIP_OKAY;
}
   

/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_DECL_HEURINIT(Heur2opt::scip_init)
{
   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_DECL_HEUREXIT(Heur2opt::scip_exit)
{
   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 */
SCIP_DECL_HEURINITSOL(Heur2opt::scip_initsol)
{
   ProbDataTSP* probdata = dynamic_cast<ProbDataTSP*>(SCIPgetObjProbData(scip));
   graph_ = probdata->getGraph();
   capture_graph(graph_);

   ncalls_ = 0;
   sol_ = NULL;
   SCIP_CALL( SCIPallocMemoryArray(scip, &tour_, graph_->nnodes) );
   
   return SCIP_OKAY;
}

   
/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_DECL_HEUREXITSOL(Heur2opt::scip_exitsol)
{
   release_graph(&graph_);
   SCIPfreeMemoryArray(scip, &tour_);

   return SCIP_OKAY;
}


/** execution method of primal heuristic 2-Opt */
SCIP_DECL_HEUREXEC(Heur2opt::scip_exec)
{  
   assert( heur != NULL );
   SCIP_SOL* sol = SCIPgetBestSol( scip );
   bool newsol;

   // check whether a new solution was found meanwhile
   if( sol != sol_ )
   {
      sol_ = sol;
      ncalls_ = 0;
      newsol = true;
   }
   else
      newsol = false;

   ncalls_++;

   int nnodes = graph_->nnodes;

   // some cases need not to be handled
   if( nnodes < 4 || sol == NULL || ncalls_ >= nnodes )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   *result= SCIP_DIDNOTFIND;

   GRAPHNODE* nodes = graph_->nodes; 

   // get tour from sol and sort edges by length, if new solution was found
   if( newsol )
   {
      GRAPHEDGE* edge;
      GRAPHEDGE* lastedge = NULL;
      GRAPHNODE* node = &nodes[0];
      int i = 0; 

      do
      {
         edge = node->first_edge;      
         while( edge != NULL )
         {
            // find the next edge of the tour 
            if( edge->back != lastedge && SCIPgetSolVal(scip, sol, edge->var) > 0.5 )
            {
               node = edge->adjac;
               lastedge = edge;

               int j;
               // shift edge through the (sorted) array 
               for(j = i; j > 0 && tour_[j-1]->length < edge->length; j-- )
                  tour_[j] = tour_[j-1];
                               
               // and insert the edge at the right position
               tour_[j] = edge; 

               i++;
               break;
            }
            edge = edge->next;

         }
      }while ( node != &nodes[0] );
      assert( i == nnodes );

   }

   GRAPHEDGE** edges2test;
   SCIP_CALL( SCIPallocBufferArray(scip, &edges2test, 4) ); 

   // test current edge with all 'longer' edges for improvement 
   // if swapping with crossing edges (though do 2Opt for one edge)
   for( int i = 0; i < ncalls_ && *result != SCIP_FOUNDSOL; i++ )
   {
      edges2test[0] = tour_[ncalls_];
      edges2test[1] = tour_[i];
      edges2test[2] = findEdge( nodes, edges2test[0]->back->adjac, edges2test[1]->back->adjac );  
      edges2test[3] = findEdge( nodes, edges2test[0]->adjac, edges2test[1]->adjac );
      assert( edges2test[2] != NULL );
      assert( edges2test[3] != NULL );
             
      // if the new solution is better and variables are not fixed, update and end
      if( edges2test[0]->length + edges2test[1]->length > edges2test[2]->length + edges2test[3]->length 
         &&  SCIPvarGetLbGlobal(edges2test[0]->var) == 0.0
         &&  SCIPvarGetLbGlobal(edges2test[1]->var) == 0.0
         &&  SCIPvarGetUbGlobal(edges2test[2]->var) == 1.0
         &&  SCIPvarGetUbGlobal(edges2test[3]->var) == 1.0 )
      {

         SCIP_Bool success;
         SCIP_SOL* swapsol; // copy of sol with 4 edges swapped 

         SCIP_CALL( SCIPcreateSol(scip, &swapsol, heur) );

         // copy the old solution
         for( int j = 0; j < nnodes; j++)
         {
            SCIP_CALL( SCIPsetSolVal(scip, swapsol, tour_[j]->var, 1.0) );
         }

         // and replace two edges
         SCIP_CALL( SCIPsetSolVal(scip, swapsol, edges2test[0]->var, 0.0) );
         SCIP_CALL( SCIPsetSolVal(scip, swapsol, edges2test[1]->var, 0.0) );
         SCIP_CALL( SCIPsetSolVal(scip, swapsol, edges2test[2]->var, 1.0) );
         SCIP_CALL( SCIPsetSolVal(scip, swapsol, edges2test[3]->var, 1.0) );
         SCIP_CALL( SCIPaddSolFree(scip, &swapsol, &success) );

         assert(success);
         *result = SCIP_FOUNDSOL;                   
         ncalls_ = 0;

      }
   }
   SCIPfreeBufferArray(scip, &edges2test);

   return SCIP_OKAY;
}

/** clone method which will be used to copy a objective plugin */
SCIP_DECL_HEURCLONE(scip::ObjCloneable* Heur2opt::clone)
{
   return new Heur2opt(scip);
}
