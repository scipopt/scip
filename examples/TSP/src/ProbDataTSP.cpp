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
   GRAPH** graph,                        /**< pointer to store the copied graph */
   GRAPH* sourcegraph                    /**< graph to be copied */
   )
{
   assert( graph != NULL );
   assert( sourcegraph != NULL );

   // copy graph the way it is created in the file reader
   int n = sourcegraph->nnodes;
   int m = sourcegraph->nedges;

   // create_graphs allocates memory for 2 anti-parallel arcs for each edge
   if(!create_graph(n, 2*m, graph))
      return SCIP_NOMEMORY;

   // copy nodes
   for(int i = 0; i < n; ++i)
   {
      GRAPHNODE * node       = &((*graph)->nodes[i]);
      GRAPHNODE * sourcenode = &(sourcegraph->nodes[i]);

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
      GRAPHNODE * nodestart = &((*graph)->nodes[i]);
      for(int j = i + 1; j < n; ++j)
      {
	 GRAPHNODE * nodeend = &((*graph)->nodes[j]);
	 
	 GRAPHEDGE * edgeforw  = &((*graph)->edges[e]);
	 GRAPHEDGE * edgebackw = &((*graph)->edges[e + m]);

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
                   
	 // dito
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
   SCIP*           scip,         /**< SCIP data structure */
   SCIP*           sourcescip,   /**< source SCIP main data structure */
   SCIP_HASHMAP*   varmap,       /**< a hashmap which stores the mapping of source variables to
				  * corresponding target variables */  
   SCIP_HASHMAP*   consmap,      /**< a hashmap which stores the mapping of source contraints to
				  * corresponding target constraints */ 
   ObjProbData**   objprobdata,  /**< pointer to store the copied problem data object */
   SCIP_Bool       global,       /**< create a global or a local copy? */
   SCIP_RESULT*    result        /**< pointer to store the result of the call */
   )
{
   // get source prob data and its graph
   ProbDataTSP * sourceprobdatatsp = NULL;
   sourceprobdatatsp = dynamic_cast<ProbDataTSP *>(SCIPgetObjProbData(sourcescip));
   assert( sourceprobdatatsp != NULL );
   GRAPH * sourcegraph = sourceprobdatatsp->graph_;
   assert( sourcegraph != NULL );

   // copy graph
   GRAPH * graph = NULL;
   SCIP_CALL( copy_graph(&graph, sourcegraph) );

   // copy and link variables
   int m = graph->nedges;
   for(int e = 0; e < m; ++e)
   {
      SCIP_Bool success;
      GRAPHEDGE * edgeforw  = &(graph->edges[e]);
      GRAPHEDGE * edgebackw = &(graph->edges[e + m]);

      assert( sourcegraph->edges[e].var != NULL );
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcegraph->edges[e].var, &(edgeforw->var), varmap, consmap, global, &success) );
      assert(success);
      assert(edgeforw->var != NULL);

      // anti-parallel arcs share variable
      edgebackw->var = edgeforw->var;
      SCIP_CALL( SCIPcaptureVar(scip, edgebackw->var) );
   }

   // allocate memory for target prob data
   ProbDataTSP * probdatatsp = new ProbDataTSP(graph);
   assert( probdatatsp != NULL );

   // save data pointer
   assert( objprobdata != NULL );
   *objprobdata = probdatatsp;
   
   // graph is captured by ProbDataTSP(graph)
   release_graph(&graph);
   
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** destructor of user problem data to free original user data (called when original problem is freed)
 *
 *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
 *  because all the work to delete the user problem data can be done in the destructor of the user problem
 *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
 *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
 *  longer needed.
 */
SCIP_RETCODE ProbDataTSP::scip_delorig(
   SCIP*              scip                /**< SCIP data structure */
   )
{
   for( int i = 0; i < graph_->nedges; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &graph_->edges[i].var) );
   }
   release_graph(&graph_);

   return SCIP_OKAY;
}

/** destructor of user problem data to free original user data (called when original problem is freed)
 *
 *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
 *  because all the work to delete the user problem data can be done in the destructor of the user problem
 *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
 *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
 *  longer needed.
 */
SCIP_RETCODE ProbDataTSP::scip_deltrans(
   SCIP*              scip                /**< SCIP data structure */
   )
{
   for( int i = 0; i < graph_->nedges; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &graph_->edges[i].var) );
   }
   release_graph(&graph_);
   
   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed)
 *
 *  The user has two possibilities to implement this method:
 *   1. Return the pointer to the original problem data object (this) as pointer to the transformed problem data
 *      object. The user may modify some internal attributes, but he has to make sure, that these modifications are
 *      reversed in the scip_deltrans() method, such that the original problem data is restored. In this case,
 *      he should set *deleteobject to FALSE, because the problem data must not be destructed by SCIP after the
 *      solving process is terminated.
 *   2. Call the copy constructor of the problem data object and return the created copy as transformed problem
 *      data object. In this case, he probably wants to set *deleteobject to TRUE, thus letting SCIP call the
 *      destructor of the object if the transformed problem data is no longer needed.
 */
SCIP_RETCODE ProbDataTSP::scip_trans(
   SCIP*              scip,               /**< SCIP data structure */
   ObjProbData**      objprobdata,        /**< pointer to store the transformed problem data object */
   SCIP_Bool*         deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
   )
{  /*lint --e{715}*/
   assert( objprobdata != NULL );
   assert( deleteobject != NULL );

   assert( graph_ != NULL );

   // copy graph
   GRAPH * transgraph = NULL;
   SCIP_CALL( copy_graph(&transgraph, graph_) );

   // copy and link variables
   int m = transgraph->nedges;
   for(int e = 0; e < m; ++e)
   {
      GRAPHEDGE * edgeforw  = &(transgraph->edges[e]);
      GRAPHEDGE * edgebackw = &(transgraph->edges[e + m]);

      assert( graph_->edges[e].var != NULL );
      SCIP_CALL( SCIPgetTransformedVar(scip, graph_->edges[e].var, &(edgeforw->var)) );
      edgebackw->var = edgeforw->var; // anti-parallel arcs share variable
      assert( edgebackw->var != NULL );
      SCIP_CALL( SCIPcaptureVar(scip, edgebackw->var) );
   }

   // allocate memory for target prob data
   ProbDataTSP * transprobdatatsp = new ProbDataTSP(transgraph);
   assert( transprobdatatsp != NULL );

   // save data pointer
   assert( objprobdata != NULL );
   *objprobdata = transprobdatatsp;
   
   // graph is captured by ProbDataTSP(graph)
   release_graph(&transgraph);
   
   *deleteobject = TRUE;

   return SCIP_OKAY;
}      
