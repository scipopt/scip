/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: HeurFarthestInsert.cpp,v 1.1 2005/03/03 16:43:34 bzfberth Exp $"

/**@file   HeurFarthestInsert.cpp
 * @brief  farthestinsert  heuristic
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <iostream>
#include <cassert>
#include "gminucut.h"
#include "HeurFarthestInsert.h"
#include "TSPProbData.h"

using namespace tsp;
using namespace std;

/** primal heuristic data */
struct HeurData
{  
   scip::ObjHeur*   objheur;         /**< primal heuristic object */
   Bool             deleteobject;    /**< should the primal heuristic object be deleted when heuristic is freed? */
};

/** method finding the edge going from the node with id index1 to the node with id index2 */
GRAPHEDGE* findEdge(
   GRAPHNODE*    nodes,              /**< all nodes of the graph */
   int           index1,             /**< id of the node where the searched edge starts */
   int           index2              /**< id of the node where the searched edge ends */
   )
{
   GRAPHEDGE* startedge;
   GRAPHEDGE* edge;

   startedge = nodes[index1].first_edge;
   assert(startedge != NULL);
   edge = startedge;

   // regard every outgoing edge of node index1 and stop if adjacent to node index2
   do
   {
      if( edge->adjac->id == index2 )
         return edge;
      edge = edge->next;
   }
   while(startedge != edge);

   return NULL;
}

/** method updating the distances of the nodes to the tour after having inserted one node with id index */
void updateDistances(
   GRAPHNODE*    nodes,              /**< all nodes of the graph */
   double*       dist,               /**< array with current distances of all nodes to the subtour */
   int           index               /**< id of the inserted node */
   )
{ 
   GRAPHEDGE* startedge;
   GRAPHEDGE* edge;
   
   startedge = nodes[index].first_edge;
   assert(startedge != NULL);
   edge = startedge;
   // regard all outgoing edges of the node and update if the length and therefore the distance of the adjacent is smaller
   do
   {
      if( dist[edge->adjac->id] > edge->length  )
         dist[edge->adjac->id] = edge->length;
      edge = edge->next;
   }
   while(startedge != edge);
}


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
RETCODE HeurFarthestInsert::scip_free(
   SCIP*         scip,               /**< SCIP data structure */
   HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
/** initialization method of primal heuristic (called after problem was transformed) */
RETCODE HeurFarthestInsert::scip_init(
   SCIP*         scip,               /**< SCIP data structure */
   HEUR*         heur                /**< the primal heuristic itself */
   )
{
   TSPProbData* probdata = dynamic_cast<TSPProbData*>(SCIPgetObjProbData(scip));
   graph_ = probdata->getGraph();
   capture_graph(graph_);
   return SCIP_OKAY;
}
   
/** deinitialization method of primal heuristic (called before transformed problem is freed) */
RETCODE HeurFarthestInsert::scip_exit(
   SCIP*         scip,               /**< SCIP data structure */
   HEUR*         heur                /**< the primal heuristic itself */
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
RETCODE HeurFarthestInsert::scip_initsol(
   SCIP*         scip,               /**< SCIP data structure */
   HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
RETCODE HeurFarthestInsert::scip_exitsol(
   SCIP*         scip,               /**< SCIP data structure */
   HEUR*         heur                /**< the primal heuristic itself */
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
RETCODE HeurFarthestInsert::scip_exec(
   SCIP*         scip,               /**< SCIP data structure */
   HEUR*         heur,               /**< the primal heuristic itself */
   RESULT*       result              /**< pointer to store the result of the heuristic call */
   )
{   
   int nnodes = graph_->nnodes;

   if( nnodes < 3 )
      *result = SCIP_DIDNOTRUN;
   else
   {   
      bool* subtour;
      int i;
      double*  dist;

      GRAPHNODE* startnode;   
      GRAPHNODE* node;
      GRAPHEDGE* edge; 
      GRAPHEDGE* startedge;
      GRAPHEDGE** bestedges;         // will contain the best insertion of a given node into a subtour
      GRAPHEDGE** edges;             // will contain some insertion of a given node into a subtour
      GRAPHEDGE** successor;         // stores the successor of a node in the current subtour    
      GRAPHNODE* nodes = graph_->nodes;
    
      for( i = 0; i < nnodes; i++ )
         assert( i == nodes[i].id );

      //memory allociation
      CHECK_OKAY( SCIPallocBufferArray(scip, &subtour, nnodes) ); 
      CHECK_OKAY( SCIPallocBufferArray(scip, &dist, nnodes) );
      CHECK_OKAY( SCIPallocBufferArray(scip, &successor, nnodes) );
      CHECK_OKAY( SCIPallocBufferArray(scip, &edges, 3) );
      CHECK_OKAY( SCIPallocBufferArray(scip, &bestedges, 3) );

      clearMemoryArray(subtour, nnodes);
      for( i = 0; i < nnodes; i++ )
         dist[i] = DBL_MAX;
      
      //building up a 3-circle
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
      for(int subtourlength = 3; subtourlength < nnodes; subtourlength++)
      {   
         //find the node with the maximal distance to the tour
         maxmin = 0.0;
         newnodeindex = -1;
         for( i = 0; i < nnodes; i++)
            
            if( maxmin < dist[i] || ( maxmin == dist[i] && !subtour[i]) )
            {
               maxmin = dist[i];
               newnodeindex = i;
            }
         assert( newnodeindex >= 0 );

         // find connection to one node in the tour 
         clearMemoryArray(bestedges, 3);
         startedge = nodes[newnodeindex].first_edge;
         startnode = NULL;
         assert(startedge != NULL);
         edge = startedge;
         do
         {
            if( subtour[edge->adjac->id] )
               break;
            edge = edge->next;
         }
         while(startedge != edge); 
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
               for( i = 0; i < 3; i++)
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

      SOL* sol;
      Bool success;

      // now create a solution out of the edges stored in successor and try to add it to SCIP
      CHECK_OKAY( SCIPcreateSol (scip, &sol, heur) );      
      for( i = 0; i < nnodes; i++)
      {
         CHECK_OKAY( SCIPsetSolVal(scip, sol, successor[i]->var, 1.0) );
      }
      CHECK_OKAY( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, &success) );
      if( success )
         *result = SCIP_FOUNDSOL;  
      CHECK_OKAY( SCIPfreeSol(scip, &sol) );

      // free all local memory
      SCIPfreeBufferArray(scip, &bestedges);
      SCIPfreeBufferArray(scip, &edges);
      SCIPfreeBufferArray(scip, &successor);
      SCIPfreeBufferArray(scip, &dist);
      SCIPfreeBufferArray(scip, &subtour);
   }

   return SCIP_OKAY;
}



