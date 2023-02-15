/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2023 Zuse Institute Berlin                                */
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

/**@file   symmetry_graph.c
 * @ingroup PUBLICCOREAPI
 * @brief  methods for dealing with symmetry detection graphs
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/symmetry_graph.h"
#include "scip/scip.h"
#include "scip/misc.h"
#include <symmetry/struct_symmetry.h>
#include <symmetry/type_symmetry.h>


/** creates and initializes a symmetry detection graph with memory for the given number of nodes and edges
 *
 *  @note at some point, the graph needs to be freed!
 */
SCIP_RETCODE SCIPcreateSymgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH**           graph,              /**< pointer to hold symmetry detection graph */
   SCIP_VAR**            symvars,            /**< variables used in symmetry detection */
   int                   nsymvars,           /**< number of variables used in symmetry detection */
   int                   nopnodes,           /**< number of operator nodes */
   int                   nvalnodes,          /**< number of value nodes */
   int                   nconsnodes,         /**< number of constraint nodes */
   int                   nedges              /**< number of edges */
   )
{
   int nnodes;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(symvars != NULL);
   assert(nopnodes >= 0);
   assert(nvalnodes >= 0);
   assert(nconsnodes >= 0);
   assert(nedges >= 0);

   nnodes = nopnodes + nvalnodes + nconsnodes;

   SCIP_CALL( SCIPallocBlockMemory(scip, graph) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*graph)->nodetypes, nnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*graph)->nodeinfopos, nnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*graph)->ops, nopnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*graph)->vals, nvalnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*graph)->conss, nconsnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*graph)->lhs, nconsnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*graph)->rhs, nconsnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*graph)->edgefirst, nedges) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*graph)->edgesecond, nedges) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*graph)->edgevals, nedges) );

   (*graph)->nnodes = 0;
   (*graph)->maxnnodes = nnodes;
   (*graph)->nopnodes = 0;
   (*graph)->maxnopnodes = nopnodes;
   (*graph)->nvalnodes = 0;
   (*graph)->maxnvalnodes = nvalnodes;
   (*graph)->nconsnodes = 0;
   (*graph)->maxnconsnodes = nconsnodes;
   (*graph)->islocked = FALSE;
   (*graph)->nedges = 0;
   (*graph)->maxnedges = nedges;
   (*graph)->symvars = symvars;
   (*graph)->nsymvars = nsymvars;
   (*graph)->nvarcolors = -1;
   (*graph)->uniqueedgetype = FALSE;

   /* to avoid reallocation, allocate memory for colors later */
   (*graph)->varcolors = NULL;
   (*graph)->opcolors = NULL;
   (*graph)->valcolors = NULL;
   (*graph)->conscolors = NULL;
   (*graph)->edgecolors = NULL;

   return SCIP_OKAY;
}

/** frees a symmetry detection graph */
SCIP_RETCODE SCIPfreeSymgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH**           graph               /**< pointer to symmetry detection graph */
   )
{
   assert(scip != NULL);
   assert(graph != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*graph)->edgecolors, (*graph)->nedges);
   SCIPfreeBlockMemoryArrayNull(scip, &(*graph)->conscolors, (*graph)->nconsnodes);
   SCIPfreeBlockMemoryArrayNull(scip, &(*graph)->valcolors, (*graph)->nvalnodes);
   SCIPfreeBlockMemoryArrayNull(scip, &(*graph)->opcolors, (*graph)->nopnodes);
   SCIPfreeBlockMemoryArrayNull(scip, &(*graph)->varcolors, (*graph)->nsymvars);

   SCIPfreeBlockMemoryArray(scip, &(*graph)->edgevals, (*graph)->maxnedges);
   SCIPfreeBlockMemoryArray(scip, &(*graph)->edgesecond, (*graph)->maxnedges);
   SCIPfreeBlockMemoryArray(scip, &(*graph)->edgefirst, (*graph)->maxnedges);
   SCIPfreeBlockMemoryArray(scip, &(*graph)->rhs, (*graph)->maxnconsnodes);
   SCIPfreeBlockMemoryArray(scip, &(*graph)->lhs, (*graph)->maxnconsnodes);
   SCIPfreeBlockMemoryArray(scip, &(*graph)->conss, (*graph)->maxnconsnodes);
   SCIPfreeBlockMemoryArray(scip, &(*graph)->vals, (*graph)->maxnvalnodes);
   SCIPfreeBlockMemoryArray(scip, &(*graph)->ops, (*graph)->maxnopnodes);
   SCIPfreeBlockMemoryArray(scip, &(*graph)->nodeinfopos, (*graph)->maxnnodes);
   SCIPfreeBlockMemoryArray(scip, &(*graph)->nodetypes, (*graph)->maxnnodes);
   SCIPfreeBlockMemory(scip, graph);

   return SCIP_OKAY;
}

/** copies an existing graph and changes variable nodes according to a permutation */
SCIP_RETCODE SCIPcopySymgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH**           graph,              /**< pointer to hold copy of symmetry detection graph */
   SYM_GRAPH*            origgraph,          /**< graph to be copied */
   int*                  perm,               /**< permutation of variables */
   SYM_SPEC              fixedtype           /**< variable types that must be fixed by symmetries */
   )
{
   SYM_NODETYPE nodetype;
   int* nodeinfopos;
   int nnodes;
   int first;
   int second;
   int i;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(origgraph != NULL);
   assert(perm != NULL);

   SCIP_CALL( SCIPcreateSymgraph(scip, graph, origgraph->symvars, origgraph->nsymvars,
         origgraph->nopnodes, origgraph->nvalnodes, origgraph->nconsnodes, origgraph->nedges) );

   /* copy nodes */
   nnodes = origgraph->nnodes;
   nodeinfopos = origgraph->nodeinfopos;
   for( i = 0; i < nnodes; ++i )
   {
      nodetype = origgraph->nodetypes[i];

      switch( nodetype )
      {
      case SYM_NODETYPE_OPERATOR:
         (void) SCIPaddSymgraphOpnode(scip, *graph, origgraph->ops[nodeinfopos[i]]);
         break;
      case SYM_NODETYPE_VAL:
         (void) SCIPaddSymgraphValnode(scip, *graph, origgraph->vals[nodeinfopos[i]]);
         break;
      default:
         assert(nodetype == SYM_NODETYPE_CONS);
         (void) SCIPaddSymgraphConsnode(scip, *graph, origgraph->conss[nodeinfopos[i]],
               origgraph->lhs[nodeinfopos[i]], origgraph->rhs[nodeinfopos[i]]);
      }
   }

   /* copy edges */
   for( i = 0; i < origgraph->nedges; ++i )
   {
      first = SCIPgetSymgraphEdgeFirst(origgraph, i);
      second = SCIPgetSymgraphEdgeSecond(origgraph, i);

      /* possibly adapt edges due to permutation of variables */
      if( first < 0 )
         first = -perm[-first - 1] - 1;
      if( second < 0 )
         second = -perm[-second - 1] - 1;

      SCIP_CALL( SCIPaddSymgraphEdge(scip, *graph, first, second,
            ! SCIPisEQ(scip, SCIPinfinity(scip), origgraph->edgevals[i]), origgraph->edgevals[i]) );
   }

   SCIP_CALL( SCIPcomputeSymgraphColors(scip, *graph, fixedtype) );

   return SCIP_OKAY;
}

/** adds a symmetry detection graph for a linear constraint to existing graph */
SCIP_RETCODE SCIPextendPermsymDetectionGraphLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_VAR**            vars,               /**< variable array of linear constraint */
   SCIP_Real*            vals,               /**< coefficients of linear constraint */
   int                   nvars,              /**< number of variables in linear constraint */
   SCIP_CONS*            cons,               /**< constraint for which we encode symmetries */
   SCIP_Real             lhs,                /**< left-hand side of constraint */
   SCIP_Real             rhs,                /**< right-hand side of constraint */
   SCIP_Bool*            success             /**< pointer to store whether graph could be built */
   )
{
   int rhsnodeidx;
   int varnodeidx;
   int i;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(nvars >= 0);
   assert(cons != NULL);
   assert(success != NULL);
   assert(! graph->islocked);

   *success = TRUE;

#ifndef NDEBUG
   /* check whether variable nodes exist in the graph */
   for( i = 0; i < nvars; ++i )
   {
      assert(0 <= SCIPvarGetProbindex(vars[i]) && SCIPvarGetProbindex(vars[i]) < graph->nsymvars);
   }
#endif

   /* create node corresponding to right-hand side */
   rhsnodeidx = SCIPaddSymgraphConsnode(scip, graph, cons, lhs, rhs);

   /* create edges */
   for( i = 0; i < nvars; ++i )
   {
      varnodeidx = SCIPgetSymgraphVarnodeidx(scip, graph, vars[i]);
      SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, rhsnodeidx, varnodeidx, TRUE, vals[i]) );
   }

   return SCIP_OKAY;
}

/** adds nodes and edges corresponding to the aggregation of a variable to a symmetry detection graph */
SCIP_RETCODE SCIPaddSymgraphVarAggegration(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   rootidx,            /**< index of root node of the aggegration */
   SCIP_VAR**            vars,               /**< array of variables in aggregation */
   SCIP_Real*            vals,               /**< coefficients of variables */
   int                   nvars,              /**< number of variables in aggregation */
   SCIP_Real             constant            /**< constant of aggregation */
   )
{
   int nodeidx;
   int j;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(rootidx >= 0);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(nvars >= 0);

#ifndef NDEBUG
   /* check whether variable nodes exist in the graph */
   for( j = 0; j < nvars; ++j )
   {
      assert(0 <= SCIPvarGetProbindex(vars[j]) && SCIPvarGetProbindex(vars[j]) < graph->nsymvars);
   }
#endif

   /* add edges incident to variables in aggregation */
   for( j = 0; j < nvars; ++j )
   {
      nodeidx = SCIPgetSymgraphVarnodeidx(scip, graph, vars[j]);
      SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, rootidx, nodeidx, TRUE, vals[j]) );
   }

   /* possibly add node for constant */
   if( ! SCIPisZero(scip, constant) )
   {
      nodeidx = SCIPaddSymgraphValnode(scip, graph, constant);

      SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, rootidx, nodeidx, FALSE, 0.0) );
   }

   return SCIP_OKAY;
}

/*
 * methods for adding and accessing nodes
 */

/** ensures that the node-based arrays in symmetry graph are sufficiently long */
static
SCIP_RETCODE ensureNodeArraysSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   addsize             /**< required additional size of node-based arrays */
   )
{
   assert(scip != NULL);
   assert(graph != NULL);
   assert(addsize > 0);

   if( graph->nnodes + addsize > graph->maxnnodes )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, graph->nnodes + addsize);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &graph->nodetypes, graph->maxnnodes, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &graph->nodeinfopos, graph->maxnnodes, newsize) );
      graph->maxnnodes = newsize;
   }

   return SCIP_OKAY;
}

/** adds an operator node to a symmetry detection graph and returns its node index */
int SCIPaddSymgraphOpnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   op                  /**< int associated with operator of node */
   )
{
   int nodeidx;

   assert(scip != NULL);
   assert(graph != NULL);

   /* we can only add nodes if symmetry colors have not been computed yet */
   if( graph->islocked )
   {
      SCIPerrorMessage("Cannot add nodes to a graph for which colors have already been computed.\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( ensureNodeArraysSize(scip, graph, 1) );

   if( graph->nopnodes >= graph->maxnopnodes )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, graph->nopnodes + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &graph->ops, graph->maxnopnodes, newsize) );
      graph->maxnopnodes = newsize;
   }

   graph->nodetypes[graph->nnodes] = SYM_NODETYPE_OPERATOR;
   graph->nodeinfopos[graph->nnodes] = graph->nopnodes;
   graph->ops[graph->nopnodes] = op;

   nodeidx = graph->nnodes;
   ++graph->nnodes;
   ++graph->nopnodes;

   return nodeidx;
}

/** adds a value node to a symmetry detection graph and returns its node index */
int SCIPaddSymgraphValnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_Real             val                 /**< value of node */
   )
{
   int nodeidx;

   assert(scip != NULL);
   assert(graph != NULL);

   /* we can only add nodes if symmetry colors have not been computed yet */
   if( graph->islocked )
   {
      SCIPerrorMessage("Cannot add nodes to a graph for which colors have already been computed.\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( ensureNodeArraysSize(scip, graph, 1) );

   if( graph->nvalnodes >= graph->maxnvalnodes )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, graph->nvalnodes + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &graph->vals, graph->maxnvalnodes, newsize) );
      graph->maxnvalnodes = newsize;
   }

   graph->nodetypes[graph->nnodes] = SYM_NODETYPE_VAL;
   graph->nodeinfopos[graph->nnodes] = graph->nvalnodes;
   graph->vals[graph->nvalnodes] = val;

   nodeidx = graph->nnodes;
   ++graph->nnodes;
   ++graph->nvalnodes;

   return nodeidx;
}

/** adds a constraint node to a symmetry detection graph and returns its node index */
int SCIPaddSymgraphConsnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_CONS*            cons,               /**< constraint of node */
   SCIP_Real             lhs,                /**< left-hand side of node */
   SCIP_Real             rhs                 /**< right-hand side of node */
   )
{
   int nodeidx;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(cons != NULL);

   /* we can only add nodes if symmetry colors have not been computed yet */
   if( graph->islocked )
   {
      SCIPerrorMessage("Cannot add nodes to a graph for which colors have already been computed.\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( ensureNodeArraysSize(scip, graph, 1) );

   if( graph->nconsnodes >= graph->maxnconsnodes )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, graph->nconsnodes + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &graph->conss, graph->maxnconsnodes, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &graph->lhs, graph->maxnconsnodes, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &graph->rhs, graph->maxnconsnodes, newsize) );
      graph->maxnconsnodes = newsize;
   }

   graph->nodetypes[graph->nnodes] = SYM_NODETYPE_CONS;
   graph->nodeinfopos[graph->nnodes] = graph->nconsnodes;
   graph->conss[graph->nconsnodes] = cons;
   graph->lhs[graph->nconsnodes] = lhs;
   graph->rhs[graph->nconsnodes] = rhs;

   nodeidx = graph->nnodes;
   ++graph->nnodes;
   ++graph->nconsnodes;

   return nodeidx;
}

/** returns the (hypothetical) node index of a variable */
int SCIPgetSymgraphVarnodeidx(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_VAR*             var                 /**< variable */
   )
{
   int nodeidx;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(var != NULL);

   nodeidx = -SCIPvarGetProbindex(var) - 1;
   assert(nodeidx != INT_MAX);

   return nodeidx;
}

/** updates the lhs of a constraint node */
SCIP_RETCODE SCIPupdateSymgraphLhs(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   nodeidx,            /**< index of constraint node in graph */
   SCIP_Real             newlhs              /**< new left-hand side of node */
   )
{
   assert(graph != NULL);
   assert(nodeidx >= 0);
   assert(graph->nodetypes[nodeidx] == SYM_NODETYPE_CONS);

   graph->lhs[graph->nodeinfopos[nodeidx]] = newlhs;

   return SCIP_OKAY;
}

/** updates the rhs of a constraint node */
SCIP_RETCODE SCIPupdateSymgraphRhs(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   nodeidx,            /**< index of constraint node in graph */
   SCIP_Real             newrhs              /**< new reft-hand side of node */
   )
{
   assert(graph != NULL);
   assert(nodeidx >= 0);
   assert(graph->nodetypes[nodeidx] == SYM_NODETYPE_CONS);

   graph->rhs[graph->nodeinfopos[nodeidx]] = newrhs;

   return SCIP_OKAY;
}

/*
 * methods for adding edges
 */

/** ensures that the edge-based arrays in symmetry graph are sufficiently long */
static
SCIP_RETCODE ensureEdgeArraysSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   addsize             /**< required additional size of edge-based arrays */
   )
{
   assert(scip != NULL);
   assert(graph != NULL);
   assert(addsize > 0);

   if( graph->nedges + addsize > graph->maxnedges )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, graph->nedges + addsize);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &graph->edgefirst, graph->maxnedges, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &graph->edgesecond, graph->maxnedges, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &graph->edgevals, graph->maxnedges, newsize) );
      graph->maxnedges = newsize;
   }

   return SCIP_OKAY;
}

/** adds an edge to a symmetry detection graph */
SCIP_RETCODE SCIPaddSymgraphEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   first,              /**< first node index of edge */
   int                   second,             /**< second node index of edge */
   SCIP_Bool             hasval,             /**< whether the edge has a value */
   SCIP_Real             val                 /**< value of the edge (is it has a value) */
   )
{
   assert(scip != NULL);
   assert(graph != NULL);

   /* we can only add edges if symmetry colors have not been computed yet */
   if( graph->islocked )
   {
      SCIPerrorMessage("Cannot add edges to a graph for which colors have already been computed.\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( ensureEdgeArraysSize(scip, graph, 1) );

   graph->edgefirst[graph->nedges] = first;
   graph->edgesecond[graph->nedges] = second;
   if( hasval )
      graph->edgevals[graph->nedges] = val;
   else
      graph->edgevals[graph->nedges] = SCIPinfinity(scip);

   graph->nedges += 1;

   return SCIP_OKAY;
}

/*
 * methods to compute colors
 */

/** compares two variables for symmetry detection
 *
 *  Variables are sorted first by their type, then by their objective coefficient,
 *  then by their lower bound, and then by their upper bound.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
int compareVars(
   SCIP_VAR*             var1,               /**< first variable for comparison */
   SCIP_VAR*             var2                /**< second variable for comparison */
   )
{
   assert(var1 != NULL);
   assert(var2 != NULL);

   if( SCIPvarGetType(var1) < SCIPvarGetType(var2) )
      return -1;
   if( SCIPvarGetType(var1) > SCIPvarGetType(var2) )
      return 1;

   if( SCIPvarGetObj(var1) < SCIPvarGetObj(var2) )
      return -1;
   if( SCIPvarGetObj(var1) > SCIPvarGetObj(var2) )
      return 1;

   if( SCIPvarGetLbGlobal(var1) < SCIPvarGetLbGlobal(var2) )
      return -1;
   if( SCIPvarGetLbGlobal(var1) > SCIPvarGetLbGlobal(var2) )
      return 1;

   if( SCIPvarGetUbGlobal(var1) < SCIPvarGetUbGlobal(var2) )
      return -1;
   if( SCIPvarGetUbGlobal(var1) > SCIPvarGetUbGlobal(var2) )
      return 1;

   return 0;
}

/** sorts nodes of a symmetry detection graph
 *
 *  Nodes are sorted first by their nodetype, then by the value of the corresponding
 *  type information, and then by their consinfo (first lhs, then rhs, then conshdlr).
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortVarnodes)
{
   SCIP_VAR** vars;

   vars = (SCIP_VAR**) dataptr;

   return compareVars(vars[ind1], vars[ind2]);
}

/** compares two operators
 *
 *  Operators are sorted by their int values.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
int compareOps(
   int                   op1,                /**< first operator in comparison */
   int                   op2                 /**< second operator in comparison */
   )
{
   if( op1 < op2 )
      return -1;
   else if( op1 > op2 )
      return 1;

   return 0;
}

/** sorts operators corresponding to SCIP_EXPRHDLR*
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortOpnodes)
{
   int* vals;

   vals = (int*) dataptr;

   return compareOps(vals[ind1], vals[ind2]);
}

/** sorts real values
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortReals)
{
   SCIP_Real* vals;

   vals = (SCIP_Real*) dataptr;

   if( vals[ind1] < vals[ind2] )
      return -1;
   if( vals[ind1] > vals[ind2] )
      return 1;

   return 0;
}

/** compares constraint nodes
 *
 *  Nodes are sorted by their type of constraint, then by the lhs, and then by the rhs.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
int compareConsnodes(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< underlying symmetry detection graph */
   int                   ind1,               /**< index of first constraint node */
   int                   ind2                /**< index of second constraint node */
   )
{
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;

   assert(graph != NULL);
   assert(0 <= ind1 && ind1 < graph->nconsnodes);
   assert(0 <= ind2 && ind2 < graph->nconsnodes);

   cons1 = graph->conss[ind1];
   cons2 = graph->conss[ind2];

   if( SCIPconsGetHdlr(cons1) < SCIPconsGetHdlr(cons2) )
      return -1;
   if( SCIPconsGetHdlr(cons1) > SCIPconsGetHdlr(cons2) )
      return 1;

   /* use SCIP's comparison functions if available */
   if( scip != NULL )
   {
      if( SCIPisLT(scip, graph->lhs[ind1], graph->lhs[ind2]) )
         return -1;
      if( SCIPisGT(scip, graph->lhs[ind1], graph->lhs[ind2]) )
         return 1;

      if( SCIPisLT(scip, graph->rhs[ind1], graph->rhs[ind2]) )
         return -1;
      if( SCIPisGT(scip, graph->rhs[ind1], graph->rhs[ind2]) )
         return 1;
   }
   else
   {
      if( graph->lhs[ind1] < graph->lhs[ind2] )
         return -1;
      if( graph->lhs[ind1] > graph->lhs[ind2] )
         return 1;

      if( graph->rhs[ind1] < graph->rhs[ind2] )
         return -1;
      if( graph->rhs[ind1] > graph->rhs[ind2] )
         return 1;
   }

   return 0;
}

/** sorts constraint nodes
 *
 *  Nodes are sorted by their type of constraint, then by the lhs, and then by the rhs.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortConsnodes)
{
   return compareConsnodes(NULL, (SYM_GRAPH*) dataptr, ind1, ind2);
}

/** sorts edges
 *
 *  Edges are sorted by their weights.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortEdges)
{
   SYM_GRAPH* G;

   G = (SYM_GRAPH*) dataptr;

   if( G->edgevals[ind1] < G->edgevals[ind2] )
      return -1;
   if( G->edgevals[ind1] > G->edgevals[ind2] )
      return 1;

   return 0;
}

/** returns whether a node of the symmetry detection graph needs to be fixed */
static
SCIP_Bool isFixedVar(
   SCIP_VAR*             var,                /**< active problem variable */
   SYM_SPEC              fixedtype           /**< variable types that must be fixed by symmetries */
   )
{
   assert(var != NULL);

   if ( (fixedtype & SYM_SPEC_INTEGER) && SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER )
      return TRUE;
   if ( (fixedtype & SYM_SPEC_BINARY) && SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      return TRUE;
   if ( (fixedtype & SYM_SPEC_REAL) &&
      (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT) )
      return TRUE;
   return FALSE;
}

/** computes colors of nodes and edges
 *
 * Colors are detected by sorting different types of nodes (variables, operators, values, and constraint) and edges.
 * If two consecutive nodes of the same type differ (e.g., different variable type), they are assigned a new color.
 */
SCIP_RETCODE SCIPcomputeSymgraphColors(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SYM_SPEC              fixedtype           /**< variable types that must be fixed by symmetries */
   )
{
   SCIP_VAR* prevvar;
   SCIP_VAR* thisvar;
   SCIP_Real prevval;
   SCIP_Real thisval;
   int* perm;
   int len;
   int i;
   int color = 0;

   assert(scip != NULL);
   assert(graph != NULL);

   /* terminate early if colors have already been computed */
   if( graph->islocked )
      return SCIP_OKAY;

   /* lock graph to be extended */
   graph->islocked = TRUE;

   /* allocate memory for colors */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &graph->varcolors, graph->nsymvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &graph->opcolors, graph->nopnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &graph->valcolors, graph->nvalnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &graph->conscolors, graph->nconsnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &graph->edgecolors, graph->nedges) );

   /* allocate permutation of arrays, will be initialized by SCIPsort() */
   len = MAX(MAX3(graph->nsymvars, graph->nopnodes, graph->nvalnodes), MAX(graph->nconsnodes, graph->nedges));
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, len) );

   /* find colors of variable nodes */
   assert(graph->nsymvars > 0);
   SCIPsort(perm, SYMsortVarnodes, (void*) graph->symvars, graph->nsymvars);

   graph->varcolors[perm[0]] = color;
   prevvar = graph->symvars[perm[0]];

   for( i = 1; i < graph->nsymvars; ++i )
   {
      thisvar = graph->symvars[perm[i]];

      if( compareVars(prevvar, thisvar) != 0 || isFixedVar(prevvar, fixedtype) )
         ++color;

      graph->varcolors[perm[i]] = color + 1;
      prevvar = thisvar;
   }
   graph->nvarcolors = color;

   /* find colors of operator nodes */
   if( graph->nopnodes > 0 )
   {
      int prevop;
      int thisop;

      SCIPsort(perm, SYMsortOpnodes, (void*) graph->ops, graph->nopnodes);

      graph->opcolors[perm[0]] = ++color;
      prevop = graph->ops[perm[0]];

      for( i = 1; i < graph->nopnodes; ++i )
      {
         thisop = graph->ops[perm[i]];

         if( compareOps(prevop, thisop) != 0 )
            ++color;

         graph->opcolors[perm[i]] = color;
         prevop = thisop;
      }
   }

   /* find colors of value nodes */
   if( graph->nvalnodes > 0 )
   {
      SCIPsort(perm, SYMsortReals, (void*) graph->vals, graph->nvalnodes);

      graph->valcolors[perm[0]] = ++color;
      prevval = graph->vals[perm[0]];

      for( i = 1; i < graph->nvalnodes; ++i )
      {
         thisval = graph->vals[perm[i]];

         if( ! SCIPisEQ(scip, prevval, thisval) )
         {
            ++color;
            prevval = thisval;
         }

         graph->valcolors[perm[i]] = color;
      }
   }

   /* find colors of constraint nodes */
   if( graph->nconsnodes > 0 )
   {
      SCIPsort(perm, SYMsortConsnodes, (void*) graph, graph->nconsnodes);

      graph->conscolors[perm[0]] = ++color;

      for( i = 1; i < graph->nconsnodes; ++i )
      {
         if( compareConsnodes(scip, graph, perm[i-1], perm[i]) != 0 )
            ++color;

         graph->conscolors[perm[i]] = color;
      }
   }

   /* find colors of edges */
   if( graph->nedges > 0 )
   {
      SCIPsort(perm, SYMsortEdges, (void*) graph, graph->nedges);

      graph->edgecolors[perm[0]] = ++color;
      prevval = graph->edgevals[perm[0]];

      for( i = 1; i < graph->nedges; ++i )
      {
         thisval = graph->edgevals[perm[i]];

         /* terminate if edges are not colored anymore */
         if( SCIPisInfinity(scip, thisval) )
            break;

         if( ! SCIPisEQ(scip, prevval, thisval) )
         {
            ++color;
            prevval = thisval;
         }

         graph->edgecolors[perm[i]] = color;
      }

      /* check whether all edges are equivalent */
      if( i == graph->nedges && graph->edgecolors[perm[0]] == graph->edgecolors[perm[i-1]] )
         graph->uniqueedgetype = TRUE;

      /* assign uncolored edges color -1 */
      for( ; i < graph->nedges; ++i )
         graph->edgecolors[perm[i]] = -1;
   }

   SCIPfreeBufferArray(scip, &perm);

   return SCIP_OKAY;
}


/* general methods */

/** returns the number of variables in a symmetry detection graph */
int SCIPgetSymgraphNVars(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   )
{
   assert(graph != NULL);

   return graph->nsymvars;
}

/** returns the number of constraint nodes in a symmetry detection graph */
int SCIPgetSymgraphNConsnodes(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   )
{
   assert(graph != NULL);

   return graph->nconsnodes;
}

/** returns the number of non-variable nodes in a graph */
int SCIPgetSymgraphNNodes(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   )
{
   assert(graph != NULL);

   return graph->nnodes;
}

/** returns the number of edges in a graph */
int SCIPgetSymgraphNEdges(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   )
{
   assert(graph != NULL);

   return graph->nedges;
}

/** return the first node index of an edge */
int SCIPgetSymgraphEdgeFirst(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   edgeidx             /**< index of edge */
   )
{
   assert(graph != NULL);
   assert(0 <= edgeidx && edgeidx < graph->nedges);

   return graph->edgefirst[edgeidx];
}

/** return the second node index of an edge */
int SCIPgetSymgraphEdgeSecond(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   edgeidx             /**< index of edge */
   )
{
   assert(graph != NULL);
   assert(0 <= edgeidx && edgeidx < graph->nedges);

   return graph->edgesecond[edgeidx];
}

/** returns the color of a variable node */
int SCIPgetSymgraphVarnodeColor(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   nodeidx             /**< index of variable node */
   )
{
   assert(graph != NULL);
   assert(0 <= nodeidx && nodeidx < graph->nsymvars);
   assert(graph->islocked);

   return graph->varcolors[nodeidx];
}

/** returns the type of a node */
SYM_NODETYPE SCIPgetSymgraphNodeType(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   nodeidx             /**< index of node */
   )
{
   assert(graph != NULL);
   assert(0 <= nodeidx && nodeidx < graph->nnodes);

   return graph->nodetypes[nodeidx];
}

/** returns the color of a non-variable node */
int SCIPgetSymgraphNodeColor(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   nodeidx             /**< index of node */
   )
{
   assert(graph != NULL);
   assert(0 <= nodeidx && nodeidx < graph->nnodes);
   assert(graph->islocked);

   switch( graph->nodetypes[nodeidx] )
   {
   case SYM_NODETYPE_OPERATOR:
      return graph->opcolors[graph->nodeinfopos[nodeidx]];
   case SYM_NODETYPE_VAL:
      return graph->valcolors[graph->nodeinfopos[nodeidx]];
   default:
      assert(graph->nodetypes[nodeidx] == SYM_NODETYPE_CONS);
   }

   return graph->conscolors[graph->nodeinfopos[nodeidx]];
}

/** returns whether an edge is colored */
SCIP_Bool SCIPisSymgraphEdgeColored(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   edgeidx             /**< index of edge */
   )
{
   assert(graph != NULL);
   assert(0 <= edgeidx && edgeidx < graph->nedges);

   if( ! graph->islocked || graph->edgecolors[edgeidx] == -1 )
      return FALSE;

   return TRUE;
}

/** returns color of an edge */
SCIP_Bool SCIPgetSymgraphEdgeColor(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   edgeidx             /**< index of edge */
   )
{
   assert(graph != NULL);
   assert(0 <= edgeidx && edgeidx < graph->nedges);
   assert(SCIPisSymgraphEdgeColored(graph, edgeidx));

   return graph->edgecolors[edgeidx];
}

/** returns the number of unique variable colors in the graph */
int SCIPgetSymgraphNVarcolors(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   )
{
   assert(graph != NULL);

   if( graph->nvarcolors < 0 )
      return graph->nsymvars;

   return graph->nvarcolors;
}

/** returns whether the graph has a unique edge type */
SCIP_Bool SCIPhasGraphUniqueEdgetype(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   )
{
   assert(graph != NULL);

   return graph->uniqueedgetype;
}

/** returns the color of a symmetry */

/** Transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant.
 *
 *  @note @p constant needs to be initialized!
 */
SCIP_RETCODE SCIPgetActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to vars array to get active variables for */
   SCIP_Real**           scalars,            /**< pointer to scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(scalars != NULL);
   assert(*vars != NULL);
   assert(*scalars != NULL);
   assert(nvars != NULL);
   assert(constant != NULL);

   if( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

      if( requiredsize > *nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, scalars, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, requiredsize, constant, &requiredsize, TRUE) );
         assert(requiredsize <= *nvars);
      }
   }
   else
   {
      for( v = 0; v < *nvars; ++v )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&(*vars)[v], &(*scalars)[v], constant) );
      }
   }
   return SCIP_OKAY;
}

/** frees symmetry information of an expression */
SCIP_RETCODE SCIPfreeSymDataExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_EXPRDATA**        symdata             /**< symmetry information of an expression */
   )
{
   assert(scip != NULL);
   assert(symdata != NULL);

   if( (*symdata)->nconstants > 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(*symdata)->constants, (*symdata)->nconstants);
   }
   if( (*symdata)->ncoefficients > 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(*symdata)->coefficients, (*symdata)->ncoefficients);
      SCIPfreeBlockMemoryArrayNull(scip, &(*symdata)->children, (*symdata)->ncoefficients);
   }
   SCIPfreeBlockMemory(scip, symdata);

   return SCIP_OKAY;
}

/** returns number of coefficients from exprdata */
int SCIPgetSymExprdataNConstants(
   SYM_EXPRDATA*         symdata             /**< symmetry information of an expression */
   )
{
   assert(symdata != NULL);

   return symdata->nconstants;
}

/** returns number of coefficients form exprdata */
SCIP_Real* SCIPgetSymExprdataConstants(
   SYM_EXPRDATA*         symdata             /**< symmetry information of an expression */
   )
{
   assert(symdata != NULL);

   return symdata->constants;
}

/** gets coefficient of expression from parent expression */
SCIP_RETCODE SCIPgetCoefSymData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression for which coefficient needs to be found */
   SCIP_EXPR*            parentexpr,         /**< parent of expr */
   SCIP_Real*            coef,               /**< buffer to store coefficient */
   SCIP_Bool*            success             /**< whether a coefficient is found */
   )
{
   SYM_EXPRDATA* symdata;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(parentexpr != NULL);
   assert(coef != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* parent does not assign coefficients to its children */
   if( ! SCIPexprhdlrHasGetSymData(SCIPexprGetHdlr(parentexpr)) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetSymDataExpr(scip, parentexpr, &symdata) );

   /* parent does not assign coefficients to its children */
   if( symdata->ncoefficients < 1 )
   {
      SCIP_CALL( SCIPfreeSymDataExpr(scip, &symdata) );

      return SCIP_OKAY;
   }

   /* search for expr in the children of parentexpr */
   for( i = 0; i < symdata->ncoefficients; ++i )
   {
      if( symdata->children[i] == expr )
      {
         *coef = symdata->coefficients[i];
         *success = TRUE;
         break;
      }
   }

   SCIP_CALL( SCIPfreeSymDataExpr(scip, &symdata) );

   return SCIP_OKAY;
}
