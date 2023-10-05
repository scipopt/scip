/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
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

/**@file   compute_symmetry_sassy_bliss.c
 * @brief  interface for symmetry computations to sassy as a preprocessor to bliss
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "compute_symmetry.h"

/* include bliss */
#include <bliss/defs.hh>
#include <bliss/graph.hh>

/* include sassy */
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4189)  // local variable is initialized but not referenced
# pragma warning(disable: 4388)  // compare signed and unsigned expression
# pragma warning(disable: 4456)  // shadowed variable
# pragma warning(disable: 4430)  // missing type specifier
#endif

/* the actual include */
#include <sassy/preprocessor.h>

#ifdef __GNUC__
#pragma GCC diagnostic warning "-Wunused-but-set-variable"
#pragma GCC diagnostic warning "-Wsign-compare"
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wshadow"
#endif

#ifdef _MSC_VER
# pragma warning(pop)
#endif

#include <sassy/tools/bliss_converter.h>

#include "scip/expr_var.h"
#include "scip/expr_sum.h"
#include "scip/expr_pow.h"
#include "scip/expr.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_linear.h"
#include "scip/scip_mem.h"
#include "scip/symmetry_graph.h"


/** struct for symmetry callback */
struct SYMMETRY_Data
{
   SCIP*                 scip;               /**< SCIP pointer */
   SYM_SYMTYPE           symtype;            /**< type of symmetries that need to be computed */
   int                   npermvars;          /**< number of variables for permutations */
   int                   nperms;             /**< number of permutations */
   int**                 perms;              /**< permutation generators as (nperms x npermvars) matrix */
   int                   nmaxperms;          /**< maximal number of permutations */
   int                   maxgenerators;      /**< maximal number of generators constructed (= 0 if unlimited) */
   SCIP_Bool             restricttovars;     /**< whether permutations shall be restricted to variables */
};


/* ------------------- hook functions ------------------- */

/** callback function for sassy */  /*lint -e{715}*/
static
void sassyhook(
   void*                 user_param,         /**< parameter supplied at call to sassy */
   int                   n,                  /**< dimension of permutations */
   const int*            aut,                /**< permutation */
   int                   nsupp,              /**< support size */
   const int*            suppa               /**< support list */
   )
{
   assert( aut != NULL );
   assert( user_param != NULL );

   SYMMETRY_Data* data = static_cast<SYMMETRY_Data*>(user_param);
   assert( data->scip != NULL );
   assert( data->maxgenerators >= 0 );

   /* make sure we do not generate more that maxgenerators many permutations */
   if ( data->maxgenerators != 0 && data->nperms >= data->maxgenerators )
      return;

   /* copy first part of automorphism */
   bool isIdentity = true;
   int* p = 0;
   int permlen;
   if ( data->restricttovars )
   {
      switch ( data->symtype )
      {
      case SYM_SYMTYPE_PERM:
         permlen = data->npermvars;
         break;
      default:
         assert( data->symtype == SYM_SYMTYPE_SIGNPERM );
         permlen = 2 * data->npermvars;
      }
   }
   else
      permlen = n;

   /* check whether permutation is identity */
   for (int j = 0; j < permlen; ++j)
   {
      if ( (int) aut[j] != j )
         isIdentity = false;
   }

   /* don't store identity permutations */
   if ( isIdentity )
      return;

   if ( SCIPallocBlockMemoryArray(data->scip, &p, permlen) != SCIP_OKAY )
      return;

   /* store symmetry */
   for (int j = 0; j < permlen; ++j)
      p[j] = (int) aut[j];

   /* check whether we should allocate space for perms */
   if ( data->nmaxperms <= 0 )
   {
      if ( data->maxgenerators == 0 )
         data->nmaxperms = 100;   /* seems to cover many cases */
      else
         data->nmaxperms = data->maxgenerators;

      if ( SCIPallocBlockMemoryArray(data->scip, &data->perms, data->nmaxperms) != SCIP_OKAY )
         return;
   }
   else if ( data->nperms >= data->nmaxperms )    /* check whether we need to resize */
   {
      int newsize = SCIPcalcMemGrowSize(data->scip, data->nperms + 1);
      assert( newsize >= data->nperms );
      assert( data->maxgenerators == 0 );

      if ( SCIPreallocBlockMemoryArray(data->scip, &data->perms, data->nmaxperms, newsize) != SCIP_OKAY )
         return;

      data->nmaxperms = newsize;
   }

   data->perms[data->nperms++] = p;
}


/* ------------------- other functions ------------------- */

/** returns whether an edge is considered in grouping process */
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
         }
         else
         {
            G->add_vertex((unsigned) curcolor, (*degrees)[*internodeid]);
            G->add_edge((unsigned) commonnodeidx, (unsigned) *internodeid);
         }
         *naddednodes += 1;

         if ( determinesize )
         {
            for (f = curstart; f < e; ++f)
            {
               ++(*degrees)[*internodeid];
               ++(*degrees)[neighbors[f]];
            }
         }
         else
         {
            for (f = curstart; f < e; ++f)
               (*G).add_edge((unsigned) neighbors[f], (unsigned) *internodeid);
         }
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
      *nnodes = SCIPgetSymgraphNNodes(graph) + nvarnodestoadd;
   else
   {
      /* add nodes for variables */
      for (j = 0; j < nvarnodestoadd; ++j)
         G->add_vertex((unsigned) SCIPgetSymgraphVarnodeColor(graph, j), (*degrees)[j]);

      /* add nodes for remaining nodes of graph */
      for (j = 0; j < SCIPgetSymgraphNNodes(graph); ++j)
         G->add_vertex((unsigned) SCIPgetSymgraphNodeColor(graph, j), (*degrees)[nvarnodestoadd + j]);
   }

   /* possibly allocate memory for degrees (will grow dynamically) */
   if ( determinesize )
   {
      *degrees = NULL;
      *maxdegrees = 0;
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *nnodes + 100) );
      for (j = 0; j < *nnodes; ++j)
         (*degrees)[j] = 0;
   }

   /* possibly initialize counters */
   if ( determinesize )
      *nedges = 0;

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
   curnnodes = 0;
   for (i = 0; i < 2; ++i)
   {
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
      internodeid = curnnodes;
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
            (*nedges) += nsymvars;
         }
         else
         {
            for (j = 0; j < nusedvars; ++j)
               G->add_edge((unsigned) nodeshift + j, (unsigned) nodeshift + j + nusedvars);
         }
      }
      nodeshift = curnnodes;
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

/** return whether symmetry can be computed */
SCIP_Bool SYMcanComputeSymmetry(void)
{
   return TRUE;
}

/** return name of external program used to compute generators */
char*
initStaticSymmetryName( )
{
   char* blissname = new char[100];
#ifdef BLISS_PATCH_PRESENT
   (void) SCIPsnprintf(blissname, 100, "bliss %sp", bliss::version);
#else
   (void) SCIPsnprintf(blissname, 100, "bliss %s", bliss::version);
#endif
   return blissname;
}

/** return name of external program used to compute generators */
char*
initStaticSymmetryAddName( )
{
   char* sassyname = new char[100];
   (void) SCIPsnprintf(sassyname, 100, "sassy %d.%d", SASSY_VERSION_MAJOR, SASSY_VERSION_MINOR);
   return sassyname;
}

static const char* symmetryname = initStaticSymmetryName();
static const char* symmetryaddname = initStaticSymmetryAddName();

/** return name of external program used to compute generators */
const char* SYMsymmetryGetName(void)
{
   return symmetryname;
}

/** return description of external program used to compute generators */
const char* SYMsymmetryGetDesc(void)
{
   return "Computing Graph Automorphisms by T. Junttila and P. Kaski (users.aalto.fi/~tjunttil/bliss/)";
}

/** return name of additional external program used for computing symmetries */
const char* SYMsymmetryGetAddName(void)
{
   return symmetryaddname;
}

/** return description of additional external program used to compute symmetries */
const char* SYMsymmetryGetAddDesc(void)
{
   return "Symmetry preprocessor by Markus Anders (github.com/markusa4/sassy)";
}

/** computes autormorphisms of a graph */
static
SCIP_RETCODE computeAutomorphisms(
   SCIP*                 scip,               /**< SCIP pointer */
   SYM_SYMTYPE           symtype,            /**< type of symmetries that need to be computed */
   sassy::static_graph*  G,                  /**< pointer to graph for that automorphisms are computed */
   int                   nsymvars,           /**< number of variables encoded in graph */
   int                   maxgenerators,      /**< maximum number of generators to be constructed (=0 if unlimited) */
   int***                perms,              /**< pointer to store generators as (nperms x npermvars) matrix */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations
                                              *   (needed for freeing storage) */
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of size of group */
   SCIP_Bool             restricttovars,     /**< whether permutations shall be restricted to variables */
   SCIP_Real*            symcodetime         /**< pointer to store the time for symmetry code */
   )
{
   SCIP_Real oldtime;

   assert( scip != NULL );
   assert( G != NULL );
   assert( maxgenerators >= 0 );
   assert( perms != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( log10groupsize != NULL );
   assert( symcodetime != NULL );

   /* init */
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *log10groupsize = 0;
   *symcodetime = 0.0;

   /* init data */
   struct SYMMETRY_Data data;
   data.scip = scip;
   data.symtype = symtype;
   data.npermvars = nsymvars;
   data.nperms = 0;
   data.nmaxperms = 0;
   data.maxgenerators = maxgenerators;
   data.perms = NULL;
   data.restricttovars = restricttovars;

   oldtime = SCIPgetSolvingTime(scip);

   /* set up sassy preprocessor */
   sassy::preprocessor sassy;

   /* turn off some preprocessing that generates redudant permutations */
   sassy::configstruct sconfig;
   sconfig.CONFIG_PREP_DEACT_PROBE = true;
   sassy.configure(&sconfig);

   /* lambda function to have access to data and pass it to sassyhook above */
   sassy::sassy_hook sassyglue = [&](int n, const int* p, int nsupp, const int* suppa) {
      sassyhook((void*)&data, n, p, nsupp, suppa);
   };

   /* call sassy to reduce graph */
   sassy.reduce(G, &sassyglue);

   /* create bliss graph */
   bliss::Graph blissgraph(0);

   /* convert sassy to bliss graph */
   convert_sassy_to_bliss(G, &blissgraph);

#ifdef SCIP_OUTPUT
   blissgraph.write_dot("debug.dot");
#endif

#ifdef SCIP_DISABLED_CODE
   char filename[SCIP_MAXSTRLEN];
   (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s.dimacs", SCIPgetProbName(scip));
   FILE* fp = fopen(filename, "w");
   if ( fp )
   {
      blissgraph.write_dimacs(fp);
      fclose(fp);
   }
#endif


   /* compute automorphisms */
   bliss::Stats stats;

   /* Prefer splitting partition cells corresponding to variables over those corresponding
    * to inequalities. This is because we are only interested in the action
    * of the automorphism group on the variables, and we need a base for this action */
   blissgraph.set_splitting_heuristic(bliss::Graph::shs_f);
   /* disable component recursion as advised by Tommi Junttila from bliss */
   blissgraph.set_component_recursion(false);

   /* do not use a node limit, but set generator limit */
#ifdef BLISS_PATCH_PRESENT
   blissgraph.set_search_limits(0, (unsigned) maxgenerators);
#endif

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   /* lambda function to have access to stats and terminate the search if maxgenerators are reached */
   long unsigned int terminatesearch = INT_MAX;
   if ( maxgenerators != 0 )
      terminatesearch = (long unsigned int) maxgenerators;
   auto term = [&]() {
      return (stats.get_nof_generators() >= terminatesearch);
   };

   auto hook = [&](unsigned int n, const unsigned int* aut) {
      sassy.bliss_hook(n, aut);
   };

   /* start search */
   blissgraph.find_automorphisms(stats, hook, term);
#else
   /* start search */
   blissgraph.find_automorphisms(stats, sassy::preprocessor::bliss_hook, (void*) &sassy);
#endif
   *symcodetime = SCIPgetSolvingTime(scip) - oldtime;

#ifdef SCIP_OUTPUT
   (void) stats.print(stdout);
#endif

   /* prepare return values */
   if ( data.nperms > 0 )
   {
      *perms = data.perms;
      *nperms = data.nperms;
      *nmaxperms = data.nmaxperms;
   }
   else
   {
      assert( data.perms == NULL );
      assert( data.nmaxperms == 0 );

      *perms = NULL;
      *nperms = 0;
      *nmaxperms = 0;
   }

   /* determine log10 of symmetry group size */
   *log10groupsize = (SCIP_Real) log10l(stats.get_group_size_approx());

   return SCIP_OKAY;
}

/** compute generators of symmetry group */
SCIP_RETCODE SYMcomputeSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of size of group */
   SCIP_Real*            symcodetime         /**< pointer to store the time for symmetry code */
   )
{
   SCIP_Bool success = FALSE;
   int* degrees;
   int maxdegrees;
   int nnodes;
   int nedges;

   assert( scip != NULL );
   assert( maxgenerators >= 0 );
   assert( graph != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( perms != NULL );
   assert( log10groupsize != NULL );
   assert( symcodetime != NULL );

   /* init */
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *log10groupsize = 0;
   *symcodetime = 0.0;

   /* determine number of nodes and edges */
   SCIP_CALL( createOrDetermineSizeGraph(scip, graph, TRUE, NULL, &nnodes, &nedges, &degrees, &maxdegrees, &success) );

   if ( ! success )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0,
         "Stopped symmetry computation: Symmetry graph would become too large.\n");
      return SCIP_OKAY;
   }

   /* create sassy graph */
   sassy::static_graph sassygraph;

   /* init graph */
   sassygraph.initialize_graph((unsigned) nnodes, (unsigned) nedges);

   /* add the nodes for linear and nonlinear constraints to the graph */
   SCIP_CALL( createOrDetermineSizeGraph(scip, graph, FALSE, &sassygraph,
         &nnodes, &nedges, &degrees, &maxdegrees, &success) );

   SCIPfreeBlockMemoryArray(scip, &degrees, maxdegrees);

   SCIPdebugMsg(scip, "Symmetry detection graph has %d nodes.\n", nnodes);

   /* compute symmetries */
   SCIP_CALL( computeAutomorphisms(scip, SCIPgetSymgraphSymtype(graph), &sassygraph, SCIPgetSymgraphNVars(graph),
         maxgenerators, perms, nperms, nmaxperms, log10groupsize, TRUE, symcodetime) );

   return SCIP_OKAY;
}

/** returns whether two given graphs are identical */
SCIP_Bool SYMcheckGraphsAreIdentical(
   SCIP*                 scip,               /**< SCIP pointer */
   SYM_SYMTYPE           symtype,            /**< type of symmetries to be checked */
   SYM_GRAPH*            G1,                 /**< first graph */
   SYM_GRAPH*            G2                  /**< second graph */
   )
{
   int** perms;
   int* degrees = NULL;
   int maxdegrees = 0;
   int nnodes;
   int nedges;
   int nperms;
   int nmaxperms;
   int nnodesfromG1;
   SCIP_Real symcodetime = 0.0;
   SCIP_Real log10groupsize;
   SCIP_Bool success;

   /* some simple checks */
   if ( G1->nnodes != G2->nnodes ||  G1->nopnodes != G2->nopnodes || G1->nvalnodes != G2->nvalnodes
      || G1->nconsnodes != G2->nconsnodes || G1->nedges != G2->nedges )
      return FALSE;

   /* determine number of nodes and edges */
   SCIP_CALL_ABORT( createOrDetermineSizeGraphCheck(scip, G1, G2, TRUE, NULL,
         &nnodes, &nedges, &degrees, &maxdegrees, &success) );

   if ( ! success )
   {
      assert( degrees == NULL );
      assert( maxdegrees == 0 );
      return FALSE;
   }
   if ( nnodes % 2 != 0 )
   {
      assert( degrees != NULL );
      assert( maxdegrees > 0 );

      SCIPfreeBlockMemoryArray(scip, &degrees, maxdegrees);
      return FALSE;
   }

   /* create sassy graph */
   sassy::static_graph sassygraph;

   /* init graph */
   sassygraph.initialize_graph((unsigned) nnodes, (unsigned) nedges);

   /* add the nodes for linear and nonlinear constraints to the graph */
   SCIP_CALL_ABORT( createOrDetermineSizeGraphCheck(scip, G1, G2, FALSE, &sassygraph,
         &nnodes, &nedges, &degrees, &maxdegrees, &success) );
   assert( success );

   SCIPfreeBlockMemoryArray(scip, &degrees, maxdegrees);

   /* compute symmetries */
   SCIP_CALL_ABORT( computeAutomorphisms(scip, SCIPgetSymgraphSymtype(G1), &sassygraph, nnodes, 0,
         &perms, &nperms, &nmaxperms, &log10groupsize, FALSE, &symcodetime) );

   /* since G1 and G2 are connected and disjoint, they are isomorphic iff there is a permutation
    * mapping a node from G1 to a node of G2
    */
   success = FALSE;
   nnodesfromG1 = nnodes / 2;
   for (int p = 0; p < nperms && ! success; ++p)
   {
      for (int i = 0; i < nnodesfromG1; ++i)
      {
         if ( perms[p][i] >= nnodesfromG1 )
         {
            success = TRUE;
            break;
         }
      }
   }

   for (int p = 0; p < nperms; ++p)
   {
      SCIPfreeBlockMemoryArray(scip, &perms[p], nnodes);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &perms, nmaxperms);

   return success;
}
