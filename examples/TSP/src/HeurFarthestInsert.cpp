/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   HeurFarthestInsert.cpp
 * @brief  farthest insert - combinatorial heuristic for TSP
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <iostream>
#include <cassert>

#include "objscip/objscip.h"
#include "GomoryHuTree.h"
#include "HeurFarthestInsert.h"
#include "ProbDataTSP.h"

using namespace tsp;
using namespace std;

/** primal heuristic data */
struct SCIP_HeurData
{  
   scip::ObjHeur*     objheur;            /**< primal heuristic object */
   SCIP_Bool          deleteobject;       /**< should the primal heuristic object be deleted when heuristic is freed? */
};

/** method finding the edge going from the node with id index1 to the node with id index2 */
GRAPHEDGE* findEdge(
   GRAPHNODE*         nodes,              /**< all nodes of the graph */
   int                index1,             /**< id of the node where the searched edge starts */
   int                index2              /**< id of the node where the searched edge ends */
   )
{
   GRAPHEDGE* edge = nodes[index1].first_edge;

   // regard every outgoing edge of node index1 and stop if adjacent to node index2
   while( edge != NULL )
   {
      if( edge->adjac->id == index2 )
         break; 
      edge = edge->next;
   }
   
   return edge;
}

/** method updating the distances of the nodes to the tour after having inserted one node with id index */
void updateDistances(
   GRAPHNODE*         nodes,              /**< all nodes of the graph */
   double*            dist,               /**< array with current distances of all nodes to the subtour */
   int                index               /**< id of the inserted node */
   )
{ 
   GRAPHEDGE* edge = nodes[index].first_edge;
   
   // regard all outgoing edges of the node and update if the length and therefore the distance of the adjacent is smaller
   while( edge != NULL )
   {
      if( dist[edge->adjac->id] > edge->length  )
         dist[edge->adjac->id] = edge->length;
      edge = edge->next;
   }
}


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_RETCODE HeurFarthestInsert::scip_free(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_RETCODE HeurFarthestInsert::scip_init(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   ProbDataTSP* probdata = dynamic_cast<ProbDataTSP*>(SCIPgetObjProbData(scip));
   graph_ = probdata->getGraph();
   capture_graph(graph_);
   return SCIP_OKAY;
}
   
/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_RETCODE HeurFarthestInsert::scip_exit(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   release_graph(&graph_);

   return SCIP_OKAY;
}
   
/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 */
SCIP_RETCODE HeurFarthestInsert::scip_initsol(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_RETCODE HeurFarthestInsert::scip_exitsol(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
/** execution method of primal heuristic
 *
 *  Searches for feasible primal solutions. The method is called in the node processing loop.
 *
 *  possible return values for *result:
 *  - SCIP_FOUNDSOL   : at least one feasible primal solution was found
 *  - SCIP_DIDNOTFIND : the heuristic searched, but did not find a feasible solution
 *  - SCIP_DIDNOTRUN  : the heuristic was skipped
 *  - SCIP_DELAYED    : the heuristic was skipped, but should be called again as soon as possible, disregarding
 *                      its frequency
 */
SCIP_RETCODE HeurFarthestInsert::scip_exec(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur,               /**< the primal heuristic itself */
   SCIP_HEURTIMING    heurtiming,         /**< current point in the node solving loop */
   SCIP_RESULT*       result              /**< pointer to store the result of the heuristic call */
   )
{   
   int nnodes = graph_->nnodes;

   if( nnodes < 3 || SCIPgetNRuns(scip) > 1 )
      *result = SCIP_DIDNOTRUN;
   else
   {   
      bool* subtour;
      int i;
      double*       dist;

      GRAPHNODE* startnode;   
      GRAPHNODE* node;
      GRAPHEDGE* edge; 

      GRAPHEDGE** bestedges;         // will contain the best insertion of a given node into a subtour
      GRAPHEDGE** edges;             // will contain some insertion of a given node into a subtour
      GRAPHEDGE** successor;         // stores the successor of a node in the current subtour    
      GRAPHNODE* nodes = graph_->nodes;
    
      for( i = 0; i < nnodes; i++ )
         assert( i == nodes[i].id );

      // memory allocation
      SCIP_CALL( SCIPallocBufferArray(scip, &subtour, nnodes) ); 
      SCIP_CALL( SCIPallocBufferArray(scip, &dist, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &successor, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &edges, 3) );
      SCIP_CALL( SCIPallocBufferArray(scip, &bestedges, 3) );

      BMSclearMemoryArray(subtour, nnodes);
      for( i = 0; i < nnodes; i++ )
         dist[i] = DBL_MAX;
      
      // building up a 3-circle
      subtour[0] = true;
      dist[0] = 0.0;
      updateDistances(nodes, dist, 0);
      subtour[1] = true;
      dist[1] = 0.0; 
      updateDistances(nodes, dist, 1);
      subtour[2] = true;
      dist[2] = 0.0;
      updateDistances(nodes, dist, 2);
      edge = findEdge(nodes,0,1);
      assert(edge != NULL);
      successor[0] = edge;
      edge = findEdge(nodes,1,2);
      assert(edge != NULL);
      successor[1] = edge;
      edge = findEdge(nodes,2,0);
      assert(edge != NULL);
      successor[2] = edge;

      double maxmin;
      double min;
      int newnodeindex;

      // widen the subtour by one node each step until you have a complete tour, actually the farthest insert heuritic
      for( int subtourlength = 3; subtourlength < nnodes; subtourlength++ )
      {   
         // find the node with the maximal distance to the tour
         maxmin = 0.0;
         newnodeindex = -1;
         for( i = 0; i < nnodes; i++)
            
            if( maxmin < dist[i] || (maxmin == dist[i] && !subtour[i]) )
            {
               maxmin = dist[i];
               newnodeindex = i;
            }
         assert(newnodeindex >= 0);

         // find connection to one node in the tour 
         BMSclearMemoryArray(bestedges, 3);
         edge = nodes[newnodeindex].first_edge;
         startnode = NULL;
        
         while( edge != NULL )
         {
            if( subtour[edge->adjac->id] )
               break;
            edge = edge->next;
         }
       
         assert(subtour[edge->adjac->id]);

         // find best insertion of the new node by trying to replace any edge connecting  by the two edges connecting
         // its end node with the new node
         min = DBL_MAX;
         edges[0] = edge;
         startnode = edge->adjac;
         node = startnode;

         // succeed to the next edge in the subtour 
         do
         {
            edges[1] = successor[node->id];
            edges[2] = findEdge(nodes, edges[1]->adjac->id, newnodeindex);

            // check, whether you have find a better insertion
            if( edges[0]->back->length - edges[1]->length + edges[2]->back->length < min)
            {
               min = edges[0]->back->length - edges[1]->length + edges[2]->back->length;
               for( i = 0; i < 3; i++ )
                  bestedges[i] = edges[i];
            } 
            node = edges[1]->adjac;
            edges[0] = edges[2]->back;
         }
         while( node != startnode);         

         // bestedges should contain a 3-cycle (modulo orientation) connecting new node with two incident ones of the tour
         assert(bestedges[0]->adjac->id == bestedges[1]->back->adjac->id);
         assert(bestedges[1]->adjac->id == bestedges[2]->back->adjac->id);
         assert(bestedges[2]->adjac->id == bestedges[0]->back->adjac->id);
         assert(subtour[bestedges[0]->adjac->id]);
         assert(subtour[bestedges[1]->adjac->id]);
         assert(bestedges[2]->adjac->id == newnodeindex);
         assert(!subtour[newnodeindex]);

         // now officially insert the new node into the tour
         successor[bestedges[0]->adjac->id] = bestedges[0]->back;
         successor[bestedges[2]->adjac->id] = bestedges[2]->back;
         dist[newnodeindex] = 0.0;
         subtour[newnodeindex] = true;
         updateDistances(nodes, dist, newnodeindex);
      }

      SCIP_SOL* sol;
      SCIP_Bool success;

      // now create a solution out of the edges stored in successor and try to add it to SCIP
      SCIP_CALL( SCIPcreateSol (scip, &sol, heur) );      
      for( i = 0; i < nnodes; i++ )
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, successor[i]->var, 1.0) );
      }
      SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, &success) );
      if( success )
         *result = SCIP_FOUNDSOL;  
      SCIP_CALL( SCIPfreeSol(scip, &sol) );

      // free all local memory
      SCIPfreeBufferArray(scip, &bestedges);
      SCIPfreeBufferArray(scip, &edges);
      SCIPfreeBufferArray(scip, &successor);
      SCIPfreeBufferArray(scip, &dist);
      SCIPfreeBufferArray(scip, &subtour);
   }

   return SCIP_OKAY;
}



