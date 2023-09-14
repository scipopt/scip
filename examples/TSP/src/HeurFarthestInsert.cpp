/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
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

/** method finding the edge going from the node with id index1 to the node with id index2 */
static
GRAPHEDGE* findEdge(
   GRAPHNODE*            nodes,              /**< all nodes of the graph */
   int                   index1,             /**< id of the node where the searched edge starts */
   int                   index2              /**< id of the node where the searched edge ends */
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
static
void updateDistances(
   GRAPHNODE*            nodes,              /**< all nodes of the graph */
   double*               dist,               /**< array with current distances of all nodes to the subtour */
   int                   idx                 /**< id of the inserted node */
   )
{
   GRAPHEDGE* edge = nodes[idx].first_edge;

   /* Regard all outgoing edges of the node and update, if the length and therefore the distance of the adjacent is
    * smaller and the edge is not fixed to 0. */
   while( edge != NULL )
   {
      if( dist[edge->adjac->id] > edge->length && SCIPvarGetUbGlobal(edge->var) != 0.0 )
         dist[edge->adjac->id] = edge->length;
      edge = edge->next;
   }
}


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_DECL_HEURFREE(HeurFarthestInsert::scip_free)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_DECL_HEURINIT(HeurFarthestInsert::scip_init)
{  /*lint --e{715}*/
   ProbDataTSP* probdata = dynamic_cast<ProbDataTSP*>(SCIPgetObjProbData(scip));
   graph_ = probdata->getGraph(); /*lint !e613*/
   capture_graph(graph_);
   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_DECL_HEUREXIT(HeurFarthestInsert::scip_exit)
{  /*lint --e{715}*/
   release_graph(&graph_);

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 */
SCIP_DECL_HEURINITSOL(HeurFarthestInsert::scip_initsol)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_DECL_HEUREXITSOL(HeurFarthestInsert::scip_exitsol)
{  /*lint --e{715}*/
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
SCIP_DECL_HEUREXEC(HeurFarthestInsert::scip_exec)
{  /*lint --e{715}*/
   int nnodes = graph_->nnodes; /*lint !e613*/
   int nedges = graph_->nedges; /*lint !e613*/

   SCIP_Bool hasFixedEdges = FALSE;
   for(int e = 0; e < nedges; ++e)
   {
      GRAPHEDGE* edge = &(graph_->edges[e]); /*lint !e613*/
      if( SCIPvarGetLbGlobal(edge->var) == 1.0 )
      {
         hasFixedEdges = TRUE;
         break;
      }
   }

   // no longer need "SCIPgetNRuns(scip) > 1" since we now respect fixed variables after restart
   if( nnodes < 3 || hasFixedEdges )
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
      GRAPHNODE* nodes = graph_->nodes; /*lint !e613*/

      for( i = 0; i < nnodes; i++ )
         assert( i == nodes[i].id );

      // memory allocation
      SCIP_CALL( SCIPallocBufferArray(scip, &subtour, nnodes) ); /*lint !e530*/
      SCIP_CALL( SCIPallocBufferArray(scip, &dist, nnodes) ); /*lint !e530*/
      SCIP_CALL( SCIPallocBufferArray(scip, &successor, nnodes) ); /*lint !e530*/
      SCIP_CALL( SCIPallocBufferArray(scip, &edges, 3) ); /*lint !e530*/
      SCIP_CALL( SCIPallocBufferArray(scip, &bestedges, 3) ); /*lint !e530*/

      BMSclearMemoryArray(subtour, nnodes);
      for( i = 0; i < nnodes; i++ )
         dist[i] = DBL_MAX;

      // building up a 3-circle, only using edges not fixed to 0
      SCIP_Bool foundThreeCircle = FALSE;
      for(int u = 0; u < nnodes - 2 && !foundThreeCircle; ++u)
      {
         for(int v = u + 1; v < nnodes - 1 && !foundThreeCircle; ++v)
         {
            GRAPHEDGE * uv = findEdge(nodes, u, v);
            assert(uv != NULL);
            if( SCIPvarGetUbGlobal(uv->var) == 0.0 )
               continue;

            for(int w = v + 1; (w < nnodes) && !foundThreeCircle; ++w) /*lint !e845*/
            {
               GRAPHEDGE * vw = findEdge(nodes, v, w);
               assert(vw != NULL);
               GRAPHEDGE * wu = findEdge(nodes, w, u);
               assert(wu != NULL);

               if( SCIPvarGetUbGlobal(vw->var) == 0.0 || SCIPvarGetUbGlobal(wu->var) == 0.0 )
                  continue;
               else
               {
                  foundThreeCircle = TRUE;

                  subtour[u] = true;
                  dist[u] = 0.0;
                  updateDistances(nodes, dist, u);
                  subtour[v] = true;
                  dist[v] = 0.0;
                  updateDistances(nodes, dist, v);
                  subtour[w] = true;
                  dist[w] = 0.0;
                  updateDistances(nodes, dist, w);
                  successor[u] = uv;
                  successor[v] = vw;
                  successor[w] = wu;
               } // foundThreeCircle with no fixed variables
            } // for w
         } // for v
      } // for u

      if( !foundThreeCircle )
      {
         *result = SCIP_DIDNOTFIND;
      }
      else
      {
         double maxmin;
         double minval;
         int newnodeindex;

         SCIP_Bool couldNotInsert = FALSE;

         // widen the subtour by one node each step until you have a complete tour, actually the farthest insert heuritic
         int subtourlength = 3;
         for(; subtourlength < nnodes; subtourlength++ )
         {
            // find the node with the maximal distance to the tour
            maxmin = 0.0;
            newnodeindex = -1;
            for( i = 0; i < nnodes; i++)
            {
               if( (maxmin < dist[i] && dist[i] != DBL_MAX) || (maxmin == dist[i] && !subtour[i]) ) /*lint !e777*/
               {
                  maxmin = dist[i];
                  newnodeindex = i;
               }
            }
            if(newnodeindex == -1)
            {
               couldNotInsert = TRUE;
               break;
            }

            // find connection to one node in the tour
            BMSclearMemoryArray(bestedges, 3);
            edge = nodes[newnodeindex].first_edge;
            startnode = NULL;

            while( edge != NULL )
            {
               if( subtour[edge->adjac->id] && SCIPvarGetUbGlobal(edge->var) != 0.0 )
                  break;
               edge = edge->next;
            }

            assert(edge != NULL);
            assert(subtour[edge->adjac->id]);

            /* find best insertion of the new node by trying to replace any edge connecting by the two edges connecting
             * its end node with the new node */
            minval = DBL_MAX;
            edges[0] = edge;
            startnode = edge->adjac;
            node = startnode;

            // succeed to the next edge in the subtour
            do
            {
               edges[1] = successor[node->id];
               edges[2] = findEdge(nodes, edges[1]->adjac->id, newnodeindex);
               assert( edges[2] != NULL );

               // check, whether you have found a better (feasible) insertion
               if( edges[0]->back->length - edges[1]->length + edges[2]->back->length < minval
                   && SCIPvarGetUbGlobal(edges[0]->var) != 0.0
                   && SCIPvarGetUbGlobal(edges[2]->var) != 0.0 )
               {
                  minval = edges[0]->back->length - edges[1]->length + edges[2]->back->length;
                  for( i = 0; i < 3; i++ )
                     bestedges[i] = edges[i];
               }
               node = edges[1]->adjac;
               edges[0] = edges[2]->back;
            }
            while( node != startnode);

            if( minval == DBL_MAX ) /*lint !e777*/
            {
               couldNotInsert = TRUE;
               break;
            }
            else
            {
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
            } // minval < DBL_MAX
         } // for subtourlength

         if(couldNotInsert)
         {
            *result = SCIP_DIDNOTFIND;
         }
         else
         {
            SCIP_SOL* sol;
            SCIP_Bool success;

            // now create a solution out of the edges stored in successor and try to add it to SCIP
            SCIP_CALL( SCIPcreateSol (scip, &sol, heur) );
            for( i = 0; i < nnodes; i++ )
            {
               SCIP_CALL( SCIPsetSolVal(scip, sol, successor[i]->var, 1.0) );
            }

            SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

            if( success )
               *result = SCIP_FOUNDSOL;
            else
               *result = SCIP_DIDNOTFIND;

            SCIP_CALL( SCIPfreeSol(scip, &sol) );
         } // couldNotInsert == FALSE
      } // foundThreeCircle == TRUE

      // free all local memory
      SCIPfreeBufferArray(scip, &bestedges);
      SCIPfreeBufferArray(scip, &edges);
      SCIPfreeBufferArray(scip, &successor);
      SCIPfreeBufferArray(scip, &dist);
      SCIPfreeBufferArray(scip, &subtour);
   }

   return SCIP_OKAY;
}

/** clone method which will be used to copy a objective plugin */
SCIP_DECL_HEURCLONE(scip::ObjCloneable* HeurFarthestInsert::clone) /*lint !e665*/
{
   return new HeurFarthestInsert(scip);
}
