/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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
 *  @note At some point, the graph needs to be freed!
 */
SCIP_RETCODE SCIPcreateSymgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SYMTYPE           symtype,            /**< type of symmetries encoded in graph */
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
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(*graph)->isfixedvar, nsymvars) );

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
   (*graph)->symtype = symtype;
   (*graph)->infinity = SCIPinfinity(scip);
   (*graph)->consnodeperm = NULL;

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
   switch( (*graph)->symtype )
   {
   case SYM_SYMTYPE_PERM:
      SCIPfreeBlockMemoryArrayNull(scip, &(*graph)->varcolors, (*graph)->nsymvars);
      break;
   default:
      assert((*graph)->symtype == SYM_SYMTYPE_SIGNPERM);
      SCIPfreeBlockMemoryArrayNull(scip, &(*graph)->varcolors, 2 * (*graph)->nsymvars);
   } /*lint !e788*/

   SCIPfreeBlockMemoryArrayNull(scip, &(*graph)->consnodeperm, (*graph)->nconsnodes);
   SCIPfreeBlockMemoryArray(scip, &(*graph)->isfixedvar, (*graph)->nsymvars);
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
{  /*lint --e{788}*/
   SYM_NODETYPE nodetype;
   int* nodeinfopos;
   int nnodes;
   int first;
   int second;
   int nodeidx;
   int i;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(origgraph != NULL);
   assert(perm != NULL);

   SCIP_CALL( SCIPcreateSymgraph(scip, origgraph->symtype, graph, origgraph->symvars, origgraph->nsymvars,
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
         SCIP_CALL( SCIPaddSymgraphOpnode(scip, *graph, origgraph->ops[nodeinfopos[i]], &nodeidx) );
         break;
      case SYM_NODETYPE_VAL:
         SCIP_CALL( SCIPaddSymgraphValnode(scip, *graph, origgraph->vals[nodeinfopos[i]], &nodeidx) );
         break;
      default:
         assert(nodetype == SYM_NODETYPE_CONS);
         SCIP_CALL( SCIPaddSymgraphConsnode(scip, *graph, origgraph->conss[nodeinfopos[i]],
               origgraph->lhs[nodeinfopos[i]], origgraph->rhs[nodeinfopos[i]], &nodeidx) );
      }
      assert(0 <= nodeidx && nodeidx < nnodes);
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
            ! SCIPisInfinity(scip, origgraph->edgevals[i]), origgraph->edgevals[i]) );
   }

   SCIP_CALL( SCIPcomputeSymgraphColors(scip, *graph, fixedtype) );

   return SCIP_OKAY;
}

/** adds a symmetry detection graph for a linear constraint to existing graph
 *
 *  For permutation symmetries, a constraint node is added that is connected to all
 *  variable nodes in the constraint. Edges are colored according to the variable coefficients.
 *  For signed permutation symmetries, also edges connecting the constraint node and
 *  the negated variable nodes are added, these edges are colored by the negative coefficients.
 */
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
   SCIP_CALL( SCIPaddSymgraphConsnode(scip, graph, cons, lhs, rhs, &rhsnodeidx) );

   /* create edges */
   for( i = 0; i < nvars; ++i )
   {
      if( SCIPgetSymgraphSymtype(graph) == SYM_SYMTYPE_SIGNPERM )
      {
         varnodeidx = SCIPgetSymgraphNegatedVarnodeidx(scip, graph, vars[i]);
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, rhsnodeidx, varnodeidx, TRUE, -vals[i]) );

         varnodeidx = SCIPgetSymgraphVarnodeidx(scip, graph, vars[i]);
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, rhsnodeidx, varnodeidx, TRUE, vals[i]) );
      }
      else
      {
         assert(SCIPgetSymgraphSymtype(graph) == SYM_SYMTYPE_PERM);

         varnodeidx = SCIPgetSymgraphVarnodeidx(scip, graph, vars[i]);
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, rhsnodeidx, varnodeidx, TRUE, vals[i]) );
      }
   }

   return SCIP_OKAY;
}

/** adds nodes and edges corresponding to the aggregation of a variable to a symmetry detection graph
 *
 *  For permutation symmetries, the root node is connected with all variable nodes in the aggregation.
 *  Edges are colored according to the variable coefficients.
 *  For signed permutation symmetries, also edges connecting the root node and the negated variable
 *  nodes are added, these edges are colored by the negative coefficients.
 */
SCIP_RETCODE SCIPaddSymgraphVarAggregation(
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
      if( SCIPgetSymgraphSymtype(graph) == SYM_SYMTYPE_SIGNPERM )
      {
         nodeidx = SCIPgetSymgraphNegatedVarnodeidx(scip, graph, vars[j]);
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, rootidx, nodeidx, TRUE, -vals[j]) );

         nodeidx = SCIPgetSymgraphVarnodeidx(scip, graph, vars[j]);
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, rootidx, nodeidx, TRUE, vals[j]) );
      }
      else
      {
         assert(SCIPgetSymgraphSymtype(graph) == SYM_SYMTYPE_PERM);

         nodeidx = SCIPgetSymgraphVarnodeidx(scip, graph, vars[j]);
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, rootidx, nodeidx, TRUE, vals[j]) );
      }
   }

   /* possibly add node for constant */
   if( ! SCIPisZero(scip, constant) )
   {
      SCIP_CALL( SCIPaddSymgraphValnode(scip, graph, constant, &nodeidx) );
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
SCIP_RETCODE SCIPaddSymgraphOpnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   op,                 /**< int associated with operator of node */
   int*                  nodeidx             /**< pointer to hold index of created node */
   )
{
   assert(scip != NULL);
   assert(graph != NULL);
   assert(nodeidx != NULL);

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

   *nodeidx = graph->nnodes;
   ++graph->nnodes;
   ++graph->nopnodes;

   return SCIP_OKAY;
}

/** adds a value node to a symmetry detection graph and returns its node index */
SCIP_RETCODE SCIPaddSymgraphValnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_Real             val,                /**< value of node */
   int*                  nodeidx             /**< pointer to hold index of created node */
   )
{
   assert(scip != NULL);
   assert(graph != NULL);
   assert(nodeidx != NULL);

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

   *nodeidx = graph->nnodes;
   ++graph->nnodes;
   ++graph->nvalnodes;

   return SCIP_OKAY;
}

/** adds a constraint node to a symmetry detection graph and returns its node index */
SCIP_RETCODE SCIPaddSymgraphConsnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_CONS*            cons,               /**< constraint of node */
   SCIP_Real             lhs,                /**< left-hand side of node */
   SCIP_Real             rhs,                /**< right-hand side of node */
   int*                  nodeidx             /**< pointer to hold index of created node */
   )
{
   assert(scip != NULL);
   assert(graph != NULL);
   assert(cons != NULL);
   assert(nodeidx != NULL);

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
   graph->lhs[graph->nconsnodes] = MAX(lhs, -graph->infinity);
   graph->rhs[graph->nconsnodes] = MIN(rhs, graph->infinity);

   *nodeidx = graph->nnodes;
   ++graph->nnodes;
   ++graph->nconsnodes;

   return SCIP_OKAY;
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
   assert(graph->symtype == SYM_SYMTYPE_PERM || graph->symtype == SYM_SYMTYPE_SIGNPERM);

   nodeidx = -SCIPvarGetProbindex(var) - 1;
   assert(nodeidx != INT_MAX);

   return nodeidx;
}

/** returns the (hypothetical) node index of a negated variable */
int SCIPgetSymgraphNegatedVarnodeidx(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_VAR*             var                 /**< variable */
   )
{
   int nodeidx;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(var != NULL);
   assert(graph->symtype == SYM_SYMTYPE_SIGNPERM);
   assert(SCIPgetSymgraphVarnodeidx(scip, graph, var) < 0 );

   nodeidx = SCIPgetSymgraphVarnodeidx(scip, graph, var) - graph->nsymvars;
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

/** registers a variable node (corresponding to active variable) to be fixed by symmetry */
SCIP_RETCODE SCIPfixSymgraphVarnode(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_VAR*             var                 /**< active variable that needs to be fixed */
   )
{
   int varidx;

   assert(graph != NULL);
   assert(var != NULL);

   varidx = SCIPvarGetProbindex(var);
   assert(0 <= varidx && varidx < graph->nsymvars);

   graph->isfixedvar[varidx] = TRUE;

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

/** compares two variables for permutation symmetry detection
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
   SCIP*                 scip,               /**< SCIP pointer (or NULL for exact comparison) */
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

   /* use SCIP's comparison functions if available */
   if( scip == NULL )
   {
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
   }
   else
   {
      if( SCIPisLT(scip, SCIPvarGetObj(var1), SCIPvarGetObj(var2)) )
         return -1;
      if( SCIPisGT(scip, SCIPvarGetObj(var1), SCIPvarGetObj(var2)) )
         return 1;

      if( SCIPisLT(scip, SCIPvarGetLbGlobal(var1), SCIPvarGetLbGlobal(var2)) )
         return -1;
      if( SCIPisGT(scip, SCIPvarGetLbGlobal(var1), SCIPvarGetLbGlobal(var2)) )
         return 1;

      if( SCIPisLT(scip, SCIPvarGetUbGlobal(var1), SCIPvarGetUbGlobal(var2)) )
         return -1;
      if( SCIPisGT(scip, SCIPvarGetUbGlobal(var1), SCIPvarGetUbGlobal(var2)) )
         return 1;
   }

   return 0;
}

/** compares two variables for permutation symmetry detection
 *
 *  Variables are sorted first by whether they are fixed, then by their type, then by
 *  their objective coefficient, then by their lower bound, and then by their upper bound.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
int compareVarsFixed(
   SCIP*                 scip,               /**< SCIP pointer (or NULL for exact comparison) */
   SCIP_VAR*             var1,               /**< first variable for comparison */
   SCIP_VAR*             var2,               /**< second variable for comparison */
   SCIP_Bool             isfixed1,           /**< whether var1 needs to be fixed */
   SCIP_Bool             isfixed2            /**< whether var2 needs to be fixed */
   )
{
   assert(var1 != NULL);
   assert(var2 != NULL);

   if( (! isfixed1) && isfixed2 )
      return -1;
   if( isfixed1 && (! isfixed2) )
      return 1;

   return compareVars(scip, var1, var2);
}

/** sorts nodes of a permutation symmetry detection graph
 *
 *  Variables are sorted first by whether they are fixed, then by their type, then by
 *  their objective coefficient, then by their lower bound, and then by their upper bound.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortVarnodesPermsym)
{
   SYM_GRAPH* graph;
   SCIP_VAR** vars;
   SCIP_Bool* isfixedvar;

   graph = (SYM_GRAPH*) dataptr;
   assert(graph != NULL);

   vars = graph->symvars;
   assert(vars != NULL);

   isfixedvar = graph->isfixedvar;
   assert(isfixedvar != NULL);

   return compareVarsFixed(NULL, vars[ind1], vars[ind2], isfixedvar[ind1], isfixedvar[ind2]);
}

/** compares two variables for signed permutation symmetry detection
 *
 *  Variables are sorted first by their type, then by their objective coefficient,
 *  then by their lower bound, and then by their upper bound.
 *  To take signed permutations into account, variable domains are centered at origin
 *  if the domain is finite.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
int compareVarsSignedPerm(
   SCIP*                 scip,               /**< SCIP pointer (or NULL for exact comparison) */
   SCIP_VAR*             var1,               /**< first variable for comparison */
   SCIP_VAR*             var2,               /**< second variable for comparison */
   SCIP_Bool             isneg1,             /**< whether var1 needs to be negated */
   SCIP_Bool             isneg2,             /**< whether var2 needs to be negated */
   SCIP_Real             infinity            /**< values as least as large as this are regarded as infinite */
   )
{
   SCIP_Real ub1;
   SCIP_Real ub2;
   SCIP_Real lb1;
   SCIP_Real lb2;
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real mid;

   assert(var1 != NULL);
   assert(var2 != NULL);

   /* use SCIP's comparison functions if available */
   if( scip == NULL )
   {
      if( SCIPvarGetType(var1) < SCIPvarGetType(var2) )
         return -1;
      if( SCIPvarGetType(var1) > SCIPvarGetType(var2) )
         return 1;
   }
   else
   {
      if( SCIPvarGetType(var1) < SCIPvarGetType(var2) )
         return -1;
      if( SCIPvarGetType(var1) > SCIPvarGetType(var2) )
         return 1;
   }

   obj1 = isneg1 ? -SCIPvarGetObj(var1) : SCIPvarGetObj(var1);
   obj2 = isneg2 ? -SCIPvarGetObj(var2) : SCIPvarGetObj(var2);

   /* use SCIP's comparison functions if available */
   if( scip == NULL )
   {
      if( obj1 < obj2 )
         return -1;
      if( obj1 > obj2 )
         return 1;
   }
   else
   {
      if( SCIPisLT(scip, obj1, obj2) )
         return -1;
      if( SCIPisGT(scip, obj1, obj2) )
         return 1;
   }

   /* adapt lower and upper bounds if domain is finite */
   lb1 = SCIPvarGetLbGlobal(var1);
   lb2 = SCIPvarGetLbGlobal(var2);
   ub1 = SCIPvarGetUbGlobal(var1);
   ub2 = SCIPvarGetUbGlobal(var2);
   if( ub1 < infinity && -lb1 < infinity )
   {
      mid = (lb1 + ub1) / 2;
      lb1 -= mid;
      ub1 -= mid;
   }
   if( ub2 < infinity && -lb2 < infinity )
   {
      mid = (lb2 + ub2) / 2;
      lb2 -= mid;
      ub2 -= mid;
   }

   /* for negated variables, flip upper and lower bounds */
   if( isneg1 )
   {
      mid = lb1;
      lb1 = -ub1;
      ub1 = -mid;
   }
   if( isneg2 )
   {
      mid = lb2;
      lb2 = -ub2;
      ub2 = -mid;
   }

   /* use SCIP's comparison functions if available */
   if( scip == NULL )
   {
      if( lb1 < lb2 )
         return -1;
      if( lb1 > lb2 )
         return 1;

      if( ub1 < ub2 )
         return -1;
      if( ub1 > ub2 )
         return 1;
   }
   else
   {
      if( SCIPisLT(scip, lb1, lb2) )
         return -1;
      if( SCIPisGT(scip, lb1, lb2) )
         return 1;

      if( SCIPisLT(scip, ub1, ub2) )
         return -1;
      if( SCIPisGT(scip, ub1, ub2) )
         return 1;
   }

   return 0;
}

/** compares two variables for signed permutation symmetry detection
 *
 *  Variables are sorted first by whether they are fixed, then by their type, then
 *  by their objective coefficient, then by their lower bound and then by their upper bound.
 *  To take signed permutations into account, variable domains are centered at origin
 *  if the domain is finite.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
int compareVarsFixedSignedPerm(
   SCIP*                 scip,               /**< SCIP pointer (or NULL for exact comparison) */
   SCIP_VAR*             var1,               /**< first variable for comparison */
   SCIP_VAR*             var2,               /**< second variable for comparison */
   SCIP_Bool             isfixed1,           /**< whether var1 needs to be fixed */
   SCIP_Bool             isfixed2,           /**< whether var2 needs to be fixed */
   SCIP_Bool             isneg1,             /**< whether var1 needs to be negated */
   SCIP_Bool             isneg2,             /**< whether var2 needs to be negated */
   SCIP_Real             infinity            /**< values as least as large as this are regarded as infinite */
   )
{
   assert(var1 != NULL);
   assert(var2 != NULL);

   if( (! isfixed1) && isfixed2 )
      return -1;
   if( isfixed1 && (! isfixed2) )
      return 1;

   return compareVarsSignedPerm(scip, var1, var2, isneg1, isneg2, infinity);
}

/** sorts nodes of a signed permutation symmetry detection graph
 *
 *  Variables are sorted first by whether they are fixed, then by their type, then
 *  by their objective coefficient, then by their lower bound and then by their upper bound.
 *  To take signed permutations into account, variable domains are centered at origin
 *  if the domain is finite.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortVarnodesSignedPermsym)
{
   SYM_GRAPH* graph;
   SCIP_VAR** vars;
   SCIP_Bool* isfixedvar;
   int nsymvars;
   int locind1;
   int locind2;
   SCIP_Bool isneg1 = FALSE;
   SCIP_Bool isneg2 = FALSE;

   graph = (SYM_GRAPH*) dataptr;
   assert(graph != NULL);

   nsymvars = graph->nsymvars;
   vars = graph->symvars;
   assert(nsymvars > 0);
   assert(vars != NULL);

   isfixedvar = graph->isfixedvar;
   assert(isfixedvar != NULL);

   locind1 = ind1;
   if( locind1 >= nsymvars )
   {
      isneg1 = TRUE;
      locind1 -= nsymvars;
   }
   locind2 = ind2;
   if( locind2 >= nsymvars )
   {
      isneg2 = TRUE;
      locind2 -= nsymvars;
   }

   return compareVarsFixedSignedPerm(NULL, vars[locind1], vars[locind2], isfixedvar[locind1], isfixedvar[locind2],
      isneg1, isneg2, graph->infinity);
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
   SCIP_Bool previsneg;
   SCIP_Bool thisisneg;
   int* perm;
   int nusedvars;
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

   /* possibly fix variables */
   for( i = 0; i < graph->nsymvars; ++i )
   {
      if( isFixedVar(graph->symvars[i], fixedtype) )
         graph->isfixedvar[i] = TRUE;
   }

   /* get number of variables used in symmetry detection graph */
   switch( graph->symtype )
   {
   case SYM_SYMTYPE_PERM:
      nusedvars = graph->nsymvars;
      break;
   default:
      assert(graph->symtype == SYM_SYMTYPE_SIGNPERM);
      nusedvars = 2 * graph->nsymvars;
   } /*lint !e788*/

   /* allocate memory for colors */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &graph->varcolors, nusedvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &graph->opcolors, graph->nopnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &graph->valcolors, graph->nvalnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &graph->conscolors, graph->nconsnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &graph->edgecolors, graph->nedges) );

   /* allocate permutation of arrays, will be initialized by SCIPsort() */
   len = graph->nedges;
   if ( graph->nopnodes > len )
      len = graph->nopnodes;
   if ( graph->nvalnodes > len )
      len = graph->nvalnodes;
   if ( graph->nconsnodes > len )
      len = graph->nconsnodes;
   if ( nusedvars > len )
      len = nusedvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &perm, len) );

   /* find colors of variable nodes */
   assert(graph->nsymvars > 0);
   switch( graph->symtype )
   {
   case SYM_SYMTYPE_PERM:
      SCIPsort(perm, SYMsortVarnodesPermsym, (void*) graph, nusedvars);

      graph->varcolors[perm[0]] = color;
      prevvar = graph->symvars[perm[0]];

      for( i = 1; i < nusedvars; ++i )
      {
         thisvar = graph->symvars[perm[i]];

         if( graph->isfixedvar[i] || compareVars(scip, prevvar, thisvar) != 0 )
            ++color;

         graph->varcolors[perm[i]] = color;
         prevvar = thisvar;
      }
      graph->nvarcolors = color;
      break;
   default:
      assert(graph->symtype == SYM_SYMTYPE_SIGNPERM);

      SCIPsort(perm, SYMsortVarnodesSignedPermsym, (void*) graph, nusedvars);

      graph->varcolors[perm[0]] = color;

      /* store information about first variable */
      if( perm[0] < graph->nsymvars )
      {
         previsneg = FALSE;
         prevvar = graph->symvars[perm[0]];
      }
      else
      {
         previsneg = TRUE;
         prevvar = graph->symvars[perm[0] - graph->nsymvars];
      }

      /* compute colors of remaining variables */
      for( i = 1; i < nusedvars; ++i )
      {
         if( perm[i] < graph->nsymvars )
         {
            thisisneg = FALSE;
            thisvar = graph->symvars[perm[i]];
         }
         else
         {
            thisisneg = TRUE;
            thisvar = graph->symvars[perm[i] - graph->nsymvars];
         }

         if( graph->isfixedvar[i % graph->nsymvars]
            || compareVarsSignedPerm(scip, prevvar, thisvar, previsneg, thisisneg, graph->infinity) != 0 )
            ++color;

         graph->varcolors[perm[i]] = color;
         prevvar = thisvar;
         previsneg = thisisneg;
      }
      graph->nvarcolors = color;
   } /*lint !e788*/

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
            ++color;

         graph->valcolors[perm[i]] = color;
         prevval = thisval;
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

      /* check whether edges are colored; due to sorting, only check first edge */
      if( SCIPisInfinity(scip, graph->edgevals[perm[0]]) )
      {
         /* all edges are uncolored */
         for( i = 0; i < graph->nedges; ++i )
            graph->edgecolors[perm[i]] = -1;
      }
      else
      {
         /* first edge is colored */
         graph->edgecolors[perm[0]] = ++color;
         prevval = graph->edgevals[perm[0]];

         for( i = 1; i < graph->nedges; ++i )
         {
            thisval = graph->edgevals[perm[i]];

            /* terminate if edges are not colored anymore */
            if( SCIPisInfinity(scip, thisval) )
               break;

            if( ! SCIPisEQ(scip, prevval, thisval) )
               ++color;

            graph->edgecolors[perm[i]] = color;
            prevval = thisval;
         }

         /* check whether all edges are equivalent */
         if( i == graph->nedges && graph->edgecolors[perm[0]] == graph->edgecolors[perm[i-1]] )
            graph->uniqueedgetype = TRUE;

         /* assign uncolored edges color -1 */
         for( ; i < graph->nedges; ++i )
            graph->edgecolors[perm[i]] = -1;
      }
   }

   SCIPfreeBufferArray(scip, &perm);

   return SCIP_OKAY;
}


/* general methods */

/** returns the type of symmetries encoded in graph */
SYM_SYMTYPE SCIPgetSymgraphSymtype(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   )
{
   assert(graph != NULL);

   return graph->symtype;
}

/** returns the variables in a symmetry detection graph */
SCIP_VAR** SCIPgetSymgraphVars(
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   )
{
   assert(graph != NULL);

   return graph->symvars;
}

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
   assert(graph->islocked);

   switch( graph->symtype )
   {
   case SYM_SYMTYPE_PERM:
      assert(0 <= nodeidx && nodeidx < graph->nsymvars);
      break;
   default:
      assert(graph->symtype == SYM_SYMTYPE_SIGNPERM);
      assert(0 <= nodeidx && nodeidx < 2 * graph->nsymvars);
   } /*lint !e788*/

   return graph->varcolors[nodeidx];
}

/** returns the type of a node */
SYM_NODETYPE SCIPgetSymgraphNodeType(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   nodeidx             /**< index of node */
   )
{
   assert(graph != NULL);
   assert(nodeidx < graph->nnodes);

   if( nodeidx < 0 )
      return SYM_NODETYPE_VAR;

   return graph->nodetypes[nodeidx];
}

/** returns the color of a non-variable node */
int SCIPgetSymgraphNodeColor(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   nodeidx             /**< index of node */
   )
{  /*lint --e{788}*/
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
int SCIPgetSymgraphEdgeColor(
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

/** creates consnodeperm array for symmetry detection graph
 *
 *  @note @p colors of symmetry detection graph must have been computed
 */
SCIP_RETCODE SCIPcreateSymgraphConsnodeperm(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   )
{
   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->islocked);

   /* either there are no constraint nodes or we have already computed the permutation */
   if( graph->nconsnodes <= 0 || graph->consnodeperm != NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &graph->consnodeperm, graph->nconsnodes) );
   SCIPsort(graph->consnodeperm, SYMsortConsnodes, (void*) graph, graph->nconsnodes);

   return SCIP_OKAY;
}

/** frees consnodeperm array for symmetry detection graph */
SCIP_RETCODE SCIPfreeSymgraphConsnodeperm(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   )
{
   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->islocked);

   SCIPfreeBlockMemoryArrayNull(scip, &graph->consnodeperm, graph->nconsnodes);

   return SCIP_OKAY;
}

/** returns consnodeperm array for symmetry detection graph
 *
 *  @note @p colors of symmetry detection graph must have been computed
 */
int* SCIPgetSymgraphConsnodeperm(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_GRAPH*            graph               /**< symmetry detection graph */
   )
{
   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->islocked);

   SCIP_CALL_ABORT( SCIPcreateSymgraphConsnodeperm(scip, graph) );

   return graph->consnodeperm;
}

/** Transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant.
 *
 *  For permutation symmetries, active variables as encoded in SCIP are used. For signed permutation symmetries,
 *  active variables are shifted such that their domain is centered at 0 (if both their upper and lower bounds
 *  are finite).
 *
 *  @note @p constant needs to be initialized!
 */
SCIP_RETCODE SCIPgetSymActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SYMTYPE           symtype,            /**< type of symmetries for which variables are required */
   SCIP_VAR***           vars,               /**< pointer to vars array to get active variables for */
   SCIP_Real**           scalars,            /**< pointer to scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   SCIP_Real ub;
   SCIP_Real lb;
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

   /* possibly post-process active variables */
   if( symtype == SYM_SYMTYPE_SIGNPERM )
   {
      /* center variables at origin if their domain is finite */
      for (v = 0; v < *nvars; ++v)
      {
         ub = SCIPvarGetUbGlobal((*vars)[v]);
         lb = SCIPvarGetLbGlobal((*vars)[v]);

         if ( SCIPisInfinity(scip, ub) || SCIPisInfinity(scip, -lb) )
            continue;

         *constant += (*scalars)[v] * (ub + lb) / 2;
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
