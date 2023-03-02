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

/**@file   ProbDataTSP.cpp
 * @brief  C++ problem data for TSP
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "objscip/objscip.h"
#include "ProbDataTSP.h"
#include "GomoryHuTree.h"

using namespace tsp;
using namespace scip;

/** copies given graph */
static
SCIP_RETCODE copy_graph(
   GRAPH**               graph,              /**< pointer to store the copied graph */
   GRAPH*                sourcegraph         /**< graph to be copied */
   )
{
   assert( graph != NULL );
   assert( sourcegraph != NULL );

   // copy graph the way it is created in the file reader
   int n = sourcegraph->nnodes;
   int m = sourcegraph->nedges;

   // create_graphs allocates memory for 2 anti-parallel arcs for each edge
   if( ! create_graph(n, 2*m, graph))
      return SCIP_NOMEMORY;

   // copy nodes
   for(int i = 0; i < n; ++i)
   {
      GRAPHNODE* node       = &((*graph)->nodes[i]);
      GRAPHNODE* sourcenode = &(sourcegraph->nodes[i]);

      assert(sourcenode->id == i);

      node->x          = sourcenode->x;
      node->y          = sourcenode->y;
      node->id         = sourcenode->id;
      node->first_edge = NULL;
   }

   // copy edges
   int e = 0;
   for(int i = 0; i < n - 1; ++i)
   {
      GRAPHNODE* nodestart = &((*graph)->nodes[i]);
      for(int j = i + 1; j < n; ++j)
      {
         GRAPHNODE* nodeend = &((*graph)->nodes[j]);
         GRAPHEDGE* edgeforw  = &((*graph)->edges[e]);
         GRAPHEDGE* edgebackw = &((*graph)->edges[e + m]);

         // construct two 'parallel' halfedges
         edgeforw->adjac  = nodeend;
         edgebackw->adjac = nodestart;
         edgeforw->back   = edgebackw;
         edgebackw->back  = edgeforw;

         // copy length
         edgeforw->length  = sourcegraph->edges[e].length;
         edgebackw->length = edgeforw->length;

         // insert one of the halfedges into the edge list of the node
         if (nodestart->first_edge == NULL)
         {
            nodestart->first_edge = edgeforw;
            nodestart->first_edge->next = NULL;
         }
         else
         {
            edgeforw->next = nodestart->first_edge;
            nodestart->first_edge = edgeforw;
         }

         // ditto
         if (nodeend->first_edge == NULL)
         {
            nodeend->first_edge = edgebackw;
            nodeend->first_edge->next = NULL;
         }
         else
         {
            edgebackw->next = nodeend->first_edge;
            nodeend->first_edge = edgebackw;
         }

         ++e;
      } // for j
   } // for i

   return SCIP_OKAY;
}

/** copies user data if you want to copy it to a subscip */
SCIP_RETCODE ProbDataTSP::scip_copy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 sourcescip,         /**< source SCIP main data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap which stores the mapping of source variables to
                                              * corresponding target variables */
   SCIP_HASHMAP*         consmap,            /**< a hashmap which stores the mapping of source contraints to
                                              * corresponding target constraints */
   ObjProbData**         objprobdata,        /**< pointer to store the copied problem data object */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_RESULT*          result              /**< pointer to store the result of the call */
   )
{
   // get source prob data and its graph
   ProbDataTSP* sourceprobdatatsp;
   sourceprobdatatsp = dynamic_cast<ProbDataTSP*>(SCIPgetObjProbData(sourcescip));
   assert( sourceprobdatatsp != NULL );

   GRAPH* sourcegraph = sourceprobdatatsp->graph_;
   assert( sourcegraph != NULL );

   // copy graph
   GRAPH* graph = NULL;
   SCIP_CALL( copy_graph(&graph, sourcegraph) );

   // copy and link variables
   int m = graph->nedges;
   for(int e = 0; e < m; ++e)
   {
      SCIP_Bool success;
      GRAPHEDGE* edgeforw  = &(graph->edges[e]);
      GRAPHEDGE* edgebackw = &(graph->edges[e + m]);
      assert( sourcegraph->edges[e].var != NULL );

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcegraph->edges[e].var, &(edgeforw->var), varmap, consmap, global, &success) );
      SCIP_CALL( SCIPcaptureVar(scip, edgeforw->var) );
      assert(success);
      assert(edgeforw->var != NULL);

      // anti-parallel arcs share variable
      edgebackw->var = edgeforw->var;
      SCIP_CALL( SCIPcaptureVar(scip, edgebackw->var) );
   }

   // allocate memory for target prob data
   ProbDataTSP* probdatatsp = new ProbDataTSP(graph);
   assert( probdatatsp != NULL );

   // save data pointer
   assert( objprobdata != NULL );
   *objprobdata = probdatatsp;

   // graph is captured by ProbDataTSP(graph)
   release_graph(&graph);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** destructor of user problem data to free original user data (called when original problem is freed) */
SCIP_RETCODE ProbDataTSP::scip_delorig(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   for( int i = 0; i < graph_->nedges; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &graph_->edges[i].back->var) );
      SCIP_CALL( SCIPreleaseVar(scip, &graph_->edges[i].var) );
   }
   release_graph(&graph_);

   return SCIP_OKAY;
}

/** destructor of user problem data to free original user data (called when original problem is freed) */
SCIP_RETCODE ProbDataTSP::scip_deltrans(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   for( int i = 0; i < graph_->nedges; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &graph_->edges[i].back->var) );
      SCIP_CALL( SCIPreleaseVar(scip, &graph_->edges[i].var) );
   }
   release_graph(&graph_);

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
SCIP_RETCODE ProbDataTSP::scip_trans(
   SCIP*                 scip,               /**< SCIP data structure */
   ObjProbData**         objprobdata,        /**< pointer to store the transformed problem data object */
   SCIP_Bool*            deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
   )
{  /*lint --e{715}*/
   assert( objprobdata != NULL );
   assert( deleteobject != NULL );

   assert( graph_ != NULL );

   // copy graph
   GRAPH* transgraph = NULL;
   SCIP_CALL( copy_graph(&transgraph, graph_) );

   // copy and link variables
   int m = transgraph->nedges;
   for(int e = 0; e < m; ++e)
   {
      GRAPHEDGE* edgeforw  = &(transgraph->edges[e]);
      GRAPHEDGE* edgebackw = &(transgraph->edges[e + m]);
      assert( graph_->edges[e].var != NULL );

      SCIP_CALL( SCIPgetTransformedVar(scip, graph_->edges[e].var, &(edgeforw->var)) );
      SCIP_CALL( SCIPcaptureVar(scip, edgeforw->var) );

      edgebackw->var = edgeforw->var; // anti-parallel arcs share variable
      assert( edgebackw->var != NULL );

      SCIP_CALL( SCIPcaptureVar(scip, edgebackw->var) );
   }

   // allocate memory for target prob data
   ProbDataTSP* transprobdatatsp = new ProbDataTSP(transgraph);
   assert( transprobdatatsp != NULL );

   // save data pointer
   assert( objprobdata != NULL );
   *objprobdata = transprobdatatsp;

   // graph is captured by ProbDataTSP(graph)
   release_graph(&transgraph);

   *deleteobject = TRUE;

   return SCIP_OKAY;
}
