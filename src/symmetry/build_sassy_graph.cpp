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

/**@file   build_sassy_graph.c
 * @brief  methods to build sassy graph for symmetry detection
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "build_sassy_graph.h"
#include "scip/expr_var.h"
#include "scip/expr_sum.h"
#include "scip/expr_pow.h"
#include "scip/expr.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_linear.h"
#include "scip/scip_mem.h"
#include "scip/symmetry_graph.h"


/* ------------------- auxiliary functions ------------------- */

/** returns whether an edge is considered in grouping process */
static
SCIP_Bool isEdgeGroupable(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   edgeidx,            /**< index of edge to be checked */
   SCIP_Bool             groupbycons         /**< whether edges are grouped by constraints */
   )
{
   int first;
   int second;

   assert(graph != NULL);

   first = SCIPgetSymgraphEdgeFirst(graph, edgeidx);
   second = SCIPgetSymgraphEdgeSecond(graph, edgeidx);

   /* uncolored edges are not grouped */
   if ( ! SCIPisSymgraphEdgeColored(graph, edgeidx) )
      return FALSE;

   /* two variable nodes are connected */
   if ( first < 0 && second < 0 )
      return FALSE;

   if ( ! groupbycons )
   {
      /* grouping by variables requires one variable node */
      if ( first < 0 || second < 0 )
         return TRUE;
   }
   else
   {
      /* check whether there is exactly one constraint node in edge */
      if ( first >= 0 && second >= 0 )
      {
         if ( (SCIPgetSymgraphNodeType(graph, first) == SYM_NODETYPE_CONS
               && SCIPgetSymgraphNodeType(graph, second) != SYM_NODETYPE_CONS)
            || (SCIPgetSymgraphNodeType(graph, first) != SYM_NODETYPE_CONS
               && SCIPgetSymgraphNodeType(graph, second) == SYM_NODETYPE_CONS) )
            return TRUE;
      }
      else if ( first >= 0 )
      {
         if ( SCIPgetSymgraphNodeType(graph, first) == SYM_NODETYPE_CONS )
            return TRUE;
      }
      else
      {
         if ( SCIPgetSymgraphNodeType(graph, second) == SYM_NODETYPE_CONS )
            return TRUE;
      }
   }

   return FALSE;
}

/** adds grouped edges all of which have one common endpoint to a graph
 *
 * The grouping mechanism works by sorting the edges according to their color. If two
 * edges have the same color, they share the same intermediate node, which is connected
 * to the common node and the other endpoints of equivalent edges.
 */
static
SCIP_RETCODE addOrDetermineEffectOfGroupedEdges(
   SCIP*                 scip,               /**< SCIP pointer */
   sassy::static_graph*  G,                  /**< graph which gets extended */
   SCIP_Bool             determinesize,      /**< whether only the effect of grouping on the graph shall be checked */
   int*                  internodeid,        /**< (initialized) pointer to store the ID of the next intermediate node */
   int**                 degrees,            /**< pointer to array of degrees for nodes in G */
   int*                  maxdegrees,         /**< (initialized) pointer to maximum number of entries degrees can hold */
   int*                  nnodes,             /**< (initialized) pointer to store the number of */
   int*                  nedges,             /**< (initialized) pointer to store the number of */
   int                   commonnodeidx,      /**< index of common node in G */
   int*                  neighbors,          /**< neighbors of common node */
   int*                  colors,             /**< colors of edges to neighbors */
   int                   nneighbors,         /**< number of neighbors */
   int*                  naddednodes,        /**< pointer to store number of nodes added to G */
   int*                  naddededges         /**< pointer to store number of edges added to G */
   )
{
   int curcolor;
   int curstart;
   int e;
   int f;

   assert( G != NULL || determinesize );
   assert( internodeid != NULL );
   assert( *internodeid >= 0 );
   assert( degrees != NULL );
   assert( maxdegrees != NULL );
   assert( *maxdegrees > 0 );
   assert( nnodes != NULL );
   assert( nedges != NULL );
   assert( neighbors != NULL );
   assert( colors != NULL );
   assert( naddednodes != NULL );
   assert( naddededges != NULL );
   assert( commonnodeidx <= *internodeid );

   *naddednodes = 0;
   *naddededges = 0;

   /* sort edges according to color */
   SCIPsortIntInt(colors, neighbors, nneighbors);

   /* iterate over groups of identical edges and group them, ignoring the last group */
   curcolor = colors[0];
   curstart = 0;
   for (e = 1; e < nneighbors; ++e)
   {
      /* if a new group started, add edges for previous group */
      if ( colors[e] != curcolor )
      {
         if ( determinesize )
         {
            SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *internodeid + 1) );
            (*degrees)[*internodeid] = 1;
            ++(*degrees)[commonnodeidx];

            for (f = curstart; f < e; ++f)
            {
               assert( neighbors[f] <= *internodeid );
               ++(*degrees)[*internodeid];
               ++(*degrees)[neighbors[f]];
            }
         }
         else
         {
            G->add_vertex((unsigned) curcolor, (*degrees)[*internodeid]);
            G->add_edge((unsigned) commonnodeidx, (unsigned) *internodeid);

            for (f = curstart; f < e; ++f)
               (*G).add_edge((unsigned) neighbors[f], (unsigned) *internodeid);
         }
         *naddednodes += 1;
         *naddededges += e - curstart + 1;
         ++(*internodeid);

         curcolor = colors[e];
         curstart = e;
      }
   }

   /* add edges of last group */
   if ( determinesize )
   {
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *internodeid + 1) );

      (*degrees)[*internodeid] = 1;
      ++(*degrees)[commonnodeidx];

      for (f = curstart; f < nneighbors; ++f)
      {
         assert( neighbors[f] <= *internodeid );
         ++(*degrees)[*internodeid];
         ++(*degrees)[neighbors[f]];
      }
   }
   else
   {
      G->add_vertex((unsigned) curcolor, (unsigned) (*degrees)[*internodeid]);
      G->add_edge((unsigned) commonnodeidx, (unsigned) *internodeid);

      for (f = curstart; f < nneighbors; ++f)
         G->add_edge((unsigned) neighbors[f], (unsigned) *internodeid);
   }
   *naddednodes += 1;
   *naddededges += nneighbors - curstart + 1;
   ++(*internodeid);

   return SCIP_OKAY;
}

/** either creates a graph or determines its size */
static
SCIP_RETCODE createOrDetermineSizeGraph(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_Bool             determinesize,      /**< whether only the size of the graph shall be determined */
   sassy::static_graph*  G,                  /**< graph to be constructed */
   int*                  nnodes,             /**< pointer to store the total number of nodes in graph */
   int*                  nedges,             /**< pointer to store the total number of edges in graph */
   int**                 degrees,            /**< pointer to store the degrees of the nodes */
   int*                  maxdegrees,         /**< pointer to store the maximal size of the degree array */
   SCIP_Bool*            success             /**< pointer to store whether the construction was successful */
   )
{
   SYM_SYMTYPE symtype;
   SYM_NODETYPE comparetype;
   SCIP_Bool groupByConstraints;
   int* groupfirsts = NULL;
   int* groupseconds = NULL;
   int* groupcolors = NULL;
   int ngroupedges = 0;
   int nvarnodestoadd;
   int internodeid;
   int nsymvars;
   int nsymedges;
   int first;
   int second;
   int color;
   int e;
   int j;

   assert( scip != NULL );
   assert( graph != NULL );
   assert( G != NULL || determinesize );
   assert( nnodes != NULL );
   assert( nedges != NULL );
   assert( degrees != NULL );
   assert( maxdegrees != NULL );
   assert( success != NULL );

   *success = TRUE;

   /* collect basic information from symmetry detection graph */
   nsymvars = SCIPgetSymgraphNVars(graph);
   symtype = SCIPgetSymgraphSymtype(graph);
   switch ( symtype )
   {
   case SYM_SYMTYPE_PERM:
      nvarnodestoadd = nsymvars;
      break;
   default:
      assert( symtype == SYM_SYMTYPE_SIGNPERM );
      nvarnodestoadd = 2 * nsymvars;
   }

   /* possibly find number of nodes in sassy graph */
   if ( determinesize )
   {
      /* initialize counters */
      *nnodes = SCIPgetSymgraphNNodes(graph) + nvarnodestoadd;
      *nedges = 0;

      /* possibly allocate memory for degrees (will grow dynamically) */
      *degrees = NULL;
      *maxdegrees = 0;
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *nnodes + 100) );
      for (j = 0; j < *nnodes; ++j)
         (*degrees)[j] = 0;
   }
   else
   {
      /* add nodes for variables */
      for (j = 0; j < nvarnodestoadd; ++j)
         G->add_vertex((unsigned) SCIPgetSymgraphVarnodeColor(graph, j), (*degrees)[j]);

      /* add nodes for remaining nodes of graph */
      for (j = 0; j < SCIPgetSymgraphNNodes(graph); ++j)
         G->add_vertex((unsigned) SCIPgetSymgraphNodeColor(graph, j), (*degrees)[nvarnodestoadd + j]);
   }

   /* determine grouping depending on the number of rhs vs. variables */
   groupByConstraints = SCIPgetSymgraphNConsnodes(graph) < SCIPgetSymgraphNVars(graph);

   /* allocate arrays to collect edges to be grouped */
   nsymedges = SCIPgetSymgraphNEdges(graph);
   SCIP_CALL( SCIPallocBufferArray(scip, &groupfirsts, nsymedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &groupseconds, nsymedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &groupcolors, nsymedges) );

   /* loop through all edges of the symmetry detection graph and either get degrees of nodes or add edges */
   internodeid = SCIPgetSymgraphNNodes(graph) + nvarnodestoadd;
   for (e = 0; e < SCIPgetSymgraphNEdges(graph); ++e)
   {
      first = SCIPgetSymgraphEdgeFirst(graph, e);
      second = SCIPgetSymgraphEdgeSecond(graph, e);

      /* get the first and second node in edge (corrected by variable shift) */
      if ( first < 0 )
         first = -first - 1;
      else
         first += nvarnodestoadd;
      if ( second < 0 )
         second = -second - 1;
      else
         second += nvarnodestoadd;
      assert(first >= 0);
      assert(second >= 0);

      /* check whether edge is used for grouping */
      if ( ! SCIPhasGraphUniqueEdgetype(graph) && isEdgeGroupable(graph, e, groupByConstraints) )
      {
         /* store edge, first becomes the cons or var node */
         comparetype = groupByConstraints ? SYM_NODETYPE_CONS : SYM_NODETYPE_VAR;

         if ( SCIPgetSymgraphNodeType(graph, SCIPgetSymgraphEdgeFirst(graph, e)) == comparetype )
         {
            groupfirsts[ngroupedges] = first;
            groupseconds[ngroupedges] = second;
         }
         else
         {
            groupfirsts[ngroupedges] = second;
            groupseconds[ngroupedges] = first;
         }
         groupcolors[ngroupedges++] = SCIPgetSymgraphEdgeColor(graph, e);
      }
      else
      {
         /* immediately add edge or increase degrees */
         assert(0 <= first && first < *nnodes);
         assert(0 <= second && second < *nnodes);

         /* possibly split edge if it is colored */
         if ( ! SCIPhasGraphUniqueEdgetype(graph) && SCIPisSymgraphEdgeColored(graph, e) )
         {
            if ( determinesize )
            {
               SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, internodeid + 1) );

               ++(*degrees)[first];
               ++(*degrees)[second];
               (*degrees)[internodeid] = 2;
               ++(*nnodes);
               *nedges += 2;
               ++internodeid;
            }
            else
            {
               assert( internodeid < *nnodes );

               color = SCIPgetSymgraphEdgeColor(graph, e);
               G->add_vertex((unsigned) color, (unsigned) (*degrees)[internodeid]);

               G->add_edge((unsigned) first, (unsigned) internodeid);
               G->add_edge((unsigned) second, (unsigned) internodeid++);
            }
         }
         else
         {
            if ( determinesize )
            {
               ++(*degrees)[first];
               ++(*degrees)[second];
               ++(*nedges);
            }
            else
            {
               if ( first < second )
                  G->add_edge((unsigned) first, (unsigned) second);
               else
                  G->add_edge((unsigned) second, (unsigned) first);
            }
         }
      }
   }

   /* possibly add groupable edges */
   if ( ngroupedges > 0 )
   {
      int firstidx = 0;
      int firstnodeidx;
      int naddednodes;
      int naddededges;

      /* sort edges according to their first nodes */
      SCIPsortIntIntInt(groupfirsts, groupseconds, groupcolors, ngroupedges);
      firstnodeidx = groupfirsts[0];

      for (j = 1; j < ngroupedges; ++j)
      {
         /* if a new first node has been found, group the edges of the previous first node; ignoring the last group */
         if ( groupfirsts[j] != firstnodeidx )
         {
            SCIP_CALL( addOrDetermineEffectOfGroupedEdges(scip, G, determinesize, &internodeid, degrees, maxdegrees,
                  nnodes, nedges, firstnodeidx, &groupseconds[firstidx], &groupcolors[firstidx], j - firstidx,
                  &naddednodes, &naddededges) );

            firstidx = j;
            firstnodeidx = groupfirsts[j];

            if ( determinesize )
            {
               *nnodes += naddednodes;
               *nedges += naddededges;
            }
         }
      }

      /* process the last group */
      SCIP_CALL( addOrDetermineEffectOfGroupedEdges(scip, G, determinesize, &internodeid, degrees, maxdegrees,
            nnodes, nedges, firstnodeidx, &groupseconds[firstidx], &groupcolors[firstidx], ngroupedges - firstidx,
            &naddednodes, &naddededges) );

      if ( determinesize )
      {
         *nnodes += naddednodes;
         *nedges += naddededges;
      }
   }

   SCIPfreeBufferArray(scip, &groupcolors);
   SCIPfreeBufferArray(scip, &groupseconds);
   SCIPfreeBufferArray(scip, &groupfirsts);

   /* for signed permutation, also add edges connecting a variable and its negation */
   switch ( SCIPgetSymgraphSymtype(graph) )
   {
   case SYM_SYMTYPE_SIGNPERM:
      if ( determinesize )
      {
         for (j = 0; j < nvarnodestoadd; ++j)
            ++(*degrees)[j];
         *nedges += nsymvars;
      }
      else
      {
         for (j = 0; j < nsymvars; ++j)
            G->add_edge((unsigned) j, (unsigned) j + nsymvars);
      }
      break;
   default:
      assert( SCIPgetSymgraphSymtype(graph) == SYM_SYMTYPE_PERM );
   }

   if ( determinesize )
   {
      SCIPdebugMsg(scip, "#nodes: %d\n", *nnodes);
      SCIPdebugMsg(scip, "#nodes for variables: %d\n", nvarnodestoadd);
      SCIPdebugMsg(scip, "#nodes for rhs: %d\n", SCIPgetSymgraphNConsnodes(graph));
      SCIPdebugMsg(scip, "#edges: %d\n", *nedges);
   }

   return SCIP_OKAY;
}

/** either creates a graph for checking symmetries or determines its size
 *
 *  The input are two graphs and the graph to be constructed consists of copies
 *  of the two input graphs, in which non-variable nodes are colored according
 *  to the colors used in symmetry detection. Each variable gets a unique color.
 */
static
SCIP_RETCODE createOrDetermineSizeGraphCheck(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_GRAPH*            graph1,             /**< first symmetry detection graph */
   SYM_GRAPH*            graph2,             /**< second symmetry detection graph */
   SCIP_Bool             determinesize,      /**< whether only the size of the graph shall be determined */
   sassy::static_graph*  G,                  /**< graph to be constructed */
   int*                  nnodes,             /**< pointer to store the total number of nodes in graph */
   int*                  nedges,             /**< pointer to store the total number of edges in graph */
   int**                 degrees,            /**< pointer to store the degrees of the nodes */
   int*                  maxdegrees,         /**< pointer to store the maximal size of the degree array */
   int*                  nnodesfromG1,       /**< pointer to store number of nodes in sassy graph arising from G1 (or NULL) */
   SCIP_Bool*            success             /**< pointer to store whether the graph could be built */
   )
{
   SYM_SYMTYPE symtype;
   SYM_NODETYPE comparetype;
   SCIP_Bool groupByConstraints;
   SYM_GRAPH* graph;
   int* nvarused1 = NULL;
   int* nvarused2 = NULL;
   int* varlabel = NULL;
   int* groupfirsts = NULL;
   int* groupseconds = NULL;
   int* groupcolors = NULL;
   int ngroupedges = 0;
   int nusedvars = 0;
   int nodeshift;
   int curnnodes;
   int nvarnodestoadd;
   int internodeid;
   int nsymvars;
   int nsymedges;
   int first;
   int second;
   int color;
   int e;
   int i;
   int j;

   assert( scip != NULL );
   assert( graph1 != NULL );
   assert( graph2 != NULL );
   assert( G != NULL || determinesize );
   assert( nnodes != NULL );
   assert( nedges != NULL );
   assert( degrees != NULL );
   assert( maxdegrees != NULL );
   assert( success != NULL );

   *success = FALSE;

   /* graphs cannot be symmetric */
   if ( SCIPgetSymgraphNEdges(graph1) != SCIPgetSymgraphNEdges(graph2)
      || SCIPgetSymgraphNVars(graph1) != SCIPgetSymgraphNVars(graph2) )
      return SCIP_OKAY;

   /* collect basic information from symmetry detection graph */
   nsymvars = SCIPgetSymgraphNVars(graph1);
   nsymedges = SCIPgetSymgraphNEdges(graph1);
   symtype = SCIPgetSymgraphSymtype(graph1);
   switch ( symtype )
   {
   case SYM_SYMTYPE_PERM:
      nvarnodestoadd = nsymvars;
      break;
   default:
      assert( symtype == SYM_SYMTYPE_SIGNPERM );
      nvarnodestoadd = 2 * nsymvars;
   }

   /* find the variables that are contained in an edge */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nvarused1, nvarnodestoadd) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nvarused2, nvarnodestoadd) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varlabel, nvarnodestoadd) );

   for (e = 0; e < nsymedges; ++e)
   {
      first = SCIPgetSymgraphEdgeFirst(graph1, e);
      second = SCIPgetSymgraphEdgeSecond(graph1, e);
      if ( first < 0 )
         nvarused1[-first - 1] += 1;
      if ( second < 0 )
         nvarused1[-second - 1] += 1;

      first = SCIPgetSymgraphEdgeFirst(graph2, e);
      second = SCIPgetSymgraphEdgeSecond(graph2, e);
      if ( first < 0 )
         nvarused2[-first - 1] += 1;
      if ( second < 0 )
         nvarused2[-second - 1] += 1;
   }

   for (j = 0; j < nvarnodestoadd; ++j)
   {
      /* graphs cannot be identical */
      if ( nvarused1[j] != nvarused2[j] )
      {
         SCIPfreeBufferArray(scip, &varlabel);
         SCIPfreeBufferArray(scip, &nvarused2);
         SCIPfreeBufferArray(scip, &nvarused1);

         return SCIP_OKAY;
      }

      /* relabel variables by restricting to variables used in constraint (or their negation) */
      if ( nvarused1[j] > 0 || nvarused1[j % SCIPgetSymgraphNVars(graph1)] > 0 )
         varlabel[j] = nusedvars++;
      else
         varlabel[j] = -1;
   }

   /* possibly find number of nodes in sassy graph and allocate memory for dynamic array */
   if ( determinesize )
   {
      *degrees = NULL;
      *maxdegrees = 0;
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees,
            SCIPgetSymgraphNNodes(graph1) + SCIPgetSymgraphNNodes(graph2) + 2 * nusedvars + 100) );

      *nnodes = 0;
      *nedges = 0;
   }

   /* determine grouping depending on the number of rhs vs. variables */
   groupByConstraints = SCIPgetSymgraphNConsnodes(graph1) < SCIPgetSymgraphNVars(graph1);

   /* allocate arrays to collect edges to be grouped */
   SCIP_CALL( SCIPallocBufferArray(scip, &groupfirsts, nsymedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &groupseconds, nsymedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &groupcolors, nsymedges) );

   /* collect information or generate graphs, we shift the node indices of the second graph when adding them to G */
   nodeshift = 0;
   for (i = 0; i < 2; ++i)
   {
      curnnodes = 0;
      graph = i == 0 ? graph1 : graph2;
      ngroupedges = 0;

      /* possibly add nodes for variables and remaining nodes, each variable gets a unique color */
      if ( determinesize )
      {
         /* add nodes for variables */
         for (j = 0; j < nvarnodestoadd; ++j)
         {
            if ( varlabel[j] >= 0 )
            {
               SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *nnodes + 1) );
               (*degrees)[nodeshift + varlabel[j]] = 0;
               ++(*nnodes);
               ++curnnodes;
            }
         }

         /* add nodes for remaining nodes of graph */
         for (j = 0; j < SCIPgetSymgraphNNodes(graph); ++j)
         {
            SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *nnodes + 1) );
            (*degrees)[nodeshift + nusedvars + j] = 0;
            ++(*nnodes);
            ++curnnodes;
         }
      }
      else
      {
         /* add nodes for variables, each variable gets a unique color */
         for (j = 0; j < nvarnodestoadd; ++j)
         {
            if ( varlabel[j] >= 0 )
            {
               G->add_vertex((unsigned) j, (*degrees)[nodeshift + varlabel[j]]);
               ++curnnodes;
            }
         }

         /* add nodes for remaining nodes of graph, ensure that colors do not conflict with variable colors */
         for (j = 0; j < SCIPgetSymgraphNNodes(graph); ++j)
         {
            G->add_vertex((unsigned) nusedvars + SCIPgetSymgraphNodeColor(graph, j),
               (*degrees)[nodeshift + nusedvars + j]);
            ++curnnodes;
         }
      }

      /* loop through all edges of the symmetry detection graph and either get degrees of nodes or add edges */
      internodeid = nodeshift + curnnodes;
      for (e = 0; e < nsymedges; ++e)
      {
         first = SCIPgetSymgraphEdgeFirst(graph, e);
         second = SCIPgetSymgraphEdgeSecond(graph, e);

         /* get the first and second node in edge (corrected by variable shift) */
         if ( first < 0 )
            first = varlabel[-first - 1];
         else
            first = nusedvars + first;
         if ( second < 0 )
            second = varlabel[-second - 1];
         else
            second = nusedvars + second;

         /* check whether edge is used for grouping */
         if ( ! SCIPhasGraphUniqueEdgetype(graph) && isEdgeGroupable(graph, e, groupByConstraints) )
         {
            /* store edge, first becomes the cons or var node */
            comparetype = groupByConstraints ? SYM_NODETYPE_CONS : SYM_NODETYPE_VAR;

            if ( SCIPgetSymgraphNodeType(graph, SCIPgetSymgraphEdgeFirst(graph, e)) == comparetype )
            {
               groupfirsts[ngroupedges] = nodeshift + first;
               groupseconds[ngroupedges] = nodeshift + second;
            }
            else
            {
               groupfirsts[ngroupedges] = nodeshift + second;
               groupseconds[ngroupedges] = nodeshift + first;
            }
            groupcolors[ngroupedges++] = nusedvars + SCIPgetSymgraphEdgeColor(graph, e);
         }
         else
         {
            /* immediately add edge or increase degrees */
            assert(0 <= first && first < *nnodes);
            assert(0 <= second && second < *nnodes);

            /* possibly split edge if it is colored */
            if ( ! SCIPhasGraphUniqueEdgetype(graph) && SCIPisSymgraphEdgeColored(graph, e) )
            {
               if ( determinesize )
               {
                  SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, internodeid + 1) );

                  ++(*degrees)[nodeshift + first];
                  ++(*degrees)[nodeshift + second];
                  (*degrees)[internodeid] = 2;
                  ++(*nnodes);
                  *nedges += 2;
               }
               else
               {
                  assert( internodeid < *nnodes );

                  color = SCIPgetSymgraphEdgeColor(graph, e);
                  G->add_vertex((unsigned) nusedvars + color, (unsigned) (*degrees)[internodeid]);
                  G->add_edge((unsigned) nodeshift + first, (unsigned) internodeid);
                  G->add_edge((unsigned) nodeshift + second, (unsigned) internodeid);
               }
               ++internodeid;
               ++curnnodes;
            }
            else
            {
               if ( determinesize )
               {
                  ++(*degrees)[nodeshift + first];
                  ++(*degrees)[nodeshift + second];
                  ++(*nedges);
               }
               else
               {
                  if ( first < second )
                     G->add_edge((unsigned) nodeshift + first, (unsigned) nodeshift + second);
                  else
                     G->add_edge((unsigned) nodeshift + second, (unsigned) nodeshift + first);
               }
            }
         }
      }

      /* possibly add groupable edges */
      if ( ngroupedges > 0 )
      {
         int firstidx = 0;
         int firstnodeidx;
         int naddednodes;
         int naddededges;

         /* sort edges according to their first nodes */
         SCIPsortIntIntInt(groupfirsts, groupseconds, groupcolors, ngroupedges);
         firstnodeidx = groupfirsts[0];

         for (j = 1; j < ngroupedges; ++j)
         {
            /* if a new first node has been found, group the edges of the previous first node; ignoring the last group */
            if ( groupfirsts[j] != firstnodeidx )
            {
               SCIP_CALL( addOrDetermineEffectOfGroupedEdges(scip, G, determinesize, &internodeid, degrees, maxdegrees,
                     nnodes, nedges, firstnodeidx, &groupseconds[firstidx], &groupcolors[firstidx], j - firstidx,
                     &naddednodes, &naddededges) );

               firstidx = j;
               firstnodeidx = groupfirsts[j];

               if ( determinesize )
               {
                  *nnodes += naddednodes;
                  *nedges += naddededges;
               }
               curnnodes += naddednodes;
            }
         }

         /* process the last group */
         SCIP_CALL( addOrDetermineEffectOfGroupedEdges(scip, G, determinesize, &internodeid, degrees, maxdegrees,
               nnodes, nedges, firstnodeidx, &groupseconds[firstidx], &groupcolors[firstidx], ngroupedges - firstidx,
               &naddednodes, &naddededges) );

         if ( determinesize )
         {
            *nnodes += naddednodes;
            *nedges += naddededges;
         }
         curnnodes += naddednodes;
      }

      /* for signed permutation, also add edges connecting a variable and its negation */
      if ( SCIPgetSymgraphSymtype(graph1) == SYM_SYMTYPE_SIGNPERM )
      {
         if ( determinesize )
         {
            for (j = 0; j < nusedvars; ++j)
               ++(*degrees)[nodeshift + j];
            (*nedges) += nusedvars / 2;
         }
         else
         {
            for (j = 0; j < nusedvars/2; ++j)
               G->add_edge((unsigned) nodeshift + j, (unsigned) nodeshift + j + nusedvars/2);
         }
      }
      nodeshift = curnnodes;

      if ( i == 0 && nnodesfromG1 != NULL )
         *nnodesfromG1 = curnnodes;
   }

   SCIPfreeBufferArray(scip, &groupcolors);
   SCIPfreeBufferArray(scip, &groupseconds);
   SCIPfreeBufferArray(scip, &groupfirsts);

   SCIPfreeBufferArray(scip, &varlabel);
   SCIPfreeBufferArray(scip, &nvarused2);
   SCIPfreeBufferArray(scip, &nvarused1);

   *success = TRUE;

   return SCIP_OKAY;
}


/** compute generators of symmetry group */
SCIP_RETCODE SYMbuildSassyGraph(
   SCIP*                 scip,               /**< SCIP pointer */
   sassy::static_graph*  sassygraph,         /**< pointer to hold sassy graph being created */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_Bool*            success             /**< pointer to store whether sassygraph could be built */
   )
{
   int* degrees;
   int maxdegrees;
   int nnodes;
   int nedges;

   assert( scip != NULL );
   assert( sassygraph != NULL );
   assert( graph != NULL );

   *success = FALSE;

   /* determine number of nodes and edges */
   SCIP_CALL( createOrDetermineSizeGraph(scip, graph, TRUE, NULL, &nnodes, &nedges, &degrees, &maxdegrees, success) );

   if ( ! *success )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0,
         "Stopped symmetry computation: Symmetry graph would become too large.\n");
      return SCIP_OKAY;
   }

   /* init graph */
   (*sassygraph).initialize_graph((unsigned) nnodes, (unsigned) nedges);

   /* add the nodes for linear and nonlinear constraints to the graph */
   SCIP_CALL( createOrDetermineSizeGraph(scip, graph, FALSE, sassygraph,
         &nnodes, &nedges, &degrees, &maxdegrees, success) );

   SCIPfreeBlockMemoryArray(scip, &degrees, maxdegrees);

   SCIPdebugMsg(scip, "Symmetry detection graph has %d nodes.\n", nnodes);

   return SCIP_OKAY;
}

/** returns whether two given graphs are identical */
SCIP_RETCODE SYMbuildSassyGraphCheck(
   SCIP*                 scip,               /**< SCIP pointer */
   sassy::static_graph*  sassygraph,         /**< pointer to hold sassy graph being created */
   SYM_GRAPH*            G1,                 /**< first graph */
   SYM_GRAPH*            G2,                 /**< second graph */
   int*                  nnodes,             /**< pointer to store number of nodes in sassy graph */
   int*                  nnodesfromG1,       /**< pointer to store number of nodes in sassy graph arising from G1 */
   SCIP_Bool*            success             /**< pointer to store whether sassygraph could be built */
   )
{
   int* degrees = NULL;
   int maxdegrees = 0;
   int nedges;

   assert( scip != NULL );
   assert( sassygraph != NULL );
   assert( G1 != NULL );
   assert( G2 != NULL );
   assert( nnodes != NULL );
   assert( nnodesfromG1 != NULL );
   assert( success != NULL );

   *success = FALSE;
   *nnodes = 0;
   *nnodesfromG1 = 0;

   /* some simple checks */
   if ( G1->nnodes != G2->nnodes ||  G1->nopnodes != G2->nopnodes || G1->nvalnodes != G2->nvalnodes
      || G1->nconsnodes != G2->nconsnodes || G1->nedges != G2->nedges )
      return SCIP_OKAY;

   /* determine number of nodes and edges */
   SCIP_CALL_ABORT( createOrDetermineSizeGraphCheck(scip, G1, G2, TRUE, NULL,
         nnodes, &nedges, &degrees, &maxdegrees, nnodesfromG1, success) );

   if ( ! *success )
   {
      assert( degrees == NULL );
      assert( maxdegrees == 0 );
      return SCIP_OKAY;
   }
   if ( *nnodes % 2 != 0 )
   {
      assert( degrees != NULL );
      assert( maxdegrees > 0 );

      SCIPfreeBlockMemoryArray(scip, &degrees, maxdegrees);
      return SCIP_OKAY;
   }

   /* init graph */
   (*sassygraph).initialize_graph((unsigned) *nnodes, (unsigned) nedges);

   /* add the nodes for linear and nonlinear constraints to the graph */
   SCIP_CALL_ABORT( createOrDetermineSizeGraphCheck(scip, G1, G2, FALSE, sassygraph,
         nnodes, &nedges, &degrees, &maxdegrees, NULL, success) );
   assert( *success );

   SCIPfreeBlockMemoryArray(scip, &degrees, maxdegrees);

   return SCIP_OKAY;
}
