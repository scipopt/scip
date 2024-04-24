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

/**@file   compute_symmetry_bliss.cpp
 * @brief  interface for symmetry computations to bliss
 * @author Marc Pfetsch
 * @author Thomas Rehn
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "compute_symmetry.h"

/* include bliss graph */
#include <bliss/defs.hh>
#include <bliss/graph.hh>

#include <string.h>
#include <vector>
#include <list>
#include <math.h>

#include "scip/expr_var.h"
#include "scip/expr_sum.h"
#include "scip/expr_pow.h"
#include "scip/expr.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_linear.h"
#include "scip/symmetry.h"
#include "scip/symmetry_graph.h"

using std::vector;


/** struct for bliss callback */
struct BLISS_Data
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

/** callback function for bliss */
static
void blisshook(
   void*                 user_param,         /**< parameter supplied at call to bliss */
   unsigned int          n,                  /**< size of aut vector */
   const unsigned int*   aut                 /**< automorphism */
   )
{
   assert( aut != NULL );
   assert( user_param != NULL );

   BLISS_Data* data = static_cast<BLISS_Data*>(user_param);
   assert( data->scip != NULL );
   assert( data->maxgenerators >= 0);

   /* make sure we do not generate more that maxgenerators many permutations, if the limit in bliss is not available */
   if ( data->maxgenerators != 0 && data->nperms >= data->maxgenerators )
      return;

   /* copy first part of automorphism */
   bool isIdentity = true;
   int* p = 0;
   int permlen;
   if ( data->restricttovars )
   {
      if ( data->symtype == SYM_SYMTYPE_PERM )
         permlen = data->npermvars;
      else
         permlen = 2 * data->npermvars;
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

/** return whether symmetry can be computed */
SCIP_Bool SYMcanComputeSymmetry(void)
{
   return TRUE;
}

/** return name of external program used to compute generators */
const char* SYMsymmetryGetName(void)
{
#ifdef BLISS_PATCH_PRESENT
   return "bliss " BLISS_VERSION "p";
#else
   return "bliss " BLISS_VERSION;
#endif
}

/** return description of external program used to compute generators */
const char* SYMsymmetryGetDesc(void)
{
   return "Computing Graph Automorphism Groups by T. Junttila and P. Kaski (users.aalto.fi/~tjunttil/bliss)";
}

/** return name of additional external program used for computing symmetries */
const char* SYMsymmetryGetAddName(void)
{
   return NULL;
}

/** return description of additional external program used to compute symmetries */
const char* SYMsymmetryGetAddDesc(void)
{
   return NULL;
}

/** returns whether an edge is considered in grouping process */
SCIP_Bool isEdgeGroupable(
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int                   edgeidx,            /**< index of edge to be checked */
   SCIP_Bool             groupbycons         /**< whether edges are grouped by constraints */
   )
{
   assert(graph != NULL);

   int first = SCIPgetSymgraphEdgeFirst(graph, edgeidx);
   int second = SCIPgetSymgraphEdgeSecond(graph, edgeidx);

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
SCIP_RETCODE addGroupedEdges(
   bliss::Graph*         G,                  /**< pointer to graph which gets extended */
   int                   commonnodeidx,      /**< index of common node in G */
   int*                  neighbors,          /**< neighbors of common node */
   int*                  colors,             /**< colors of edges to neighbors */
   int                   nneighbors,         /**< number of neighbors */
   int*                  naddednodes,        /**< pointer to store number of nodes added to G */
   int*                  naddededges         /**< pointer to store number of edges added to G */
   )
{
   assert( G != NULL );
   assert( neighbors != NULL );
   assert( colors != NULL );
   assert( naddednodes != NULL );
   assert( naddededges != NULL );

   *naddednodes = 0;
   *naddededges = 0;

   /* sort edges according to color */
   SCIPsortIntInt(colors, neighbors, nneighbors);

   /* iterate over groups of identical edges and group them, ignoring the last group */
   int curcolor = colors[0];
   int curstart = 0;
   for (int e = 1; e < nneighbors; ++e)
   {
      /* if a new group started, add edges for previous group */
      if ( colors[e] != curcolor )
      {
         int internode = (*G).add_vertex(curcolor);
         (*G).add_edge(commonnodeidx, internode);
         *naddednodes += 1;

         for (int f = curstart; f < e; ++f)
            (*G).add_edge(internode, neighbors[f]);
         *naddededges += e - curstart + 1;

         curcolor = colors[e];
         curstart = e;
      }
   }

   /* add edges of last group */
   int internode = (*G).add_vertex(curcolor);
   (*G).add_edge(commonnodeidx, internode);
   *naddednodes += 1;

   for (int f = curstart; f < nneighbors; ++f)
      (*G).add_edge(internode, neighbors[f]);
   *naddededges += nneighbors - curstart + 1;

   return SCIP_OKAY;
}

/** computes autormorphisms of a graph */
static
SCIP_RETCODE computeAutomorphisms(
   SCIP*                 scip,               /**< SCIP pointer */
   SYM_SYMTYPE           symtype,            /**< type of symmetries that need to be computed */
   bliss::Graph*         G,                  /**< pointer to graph for that automorphisms are computed */
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
   assert( perms != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( log10groupsize != NULL );
   assert( maxgenerators >= 0 );
   assert( symcodetime != NULL );

   /* init */
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *log10groupsize = 0;
   *symcodetime = 0.0;

   bliss::Stats stats;
   BLISS_Data data;
   data.scip = scip;
   data.symtype = symtype;
   data.npermvars = nsymvars;
   data.nperms = 0;
   data.nmaxperms = 0;
   data.maxgenerators = maxgenerators;
   data.perms = NULL;
   data.restricttovars = restricttovars;

   /* Prefer splitting partition cells corresponding to variables over those corresponding
    * to inequalities. This is because we are only interested in the action
    * of the automorphism group on the variables, and we need a base for this action */
   G->set_splitting_heuristic(bliss::Graph::shs_f);
   /* disable component recursion as advised by Tommi Junttila from bliss */
   G->set_component_recursion(false);

   oldtime = SCIPgetSolvingTime(scip);
#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   /* lambda function to have access to data and pass it to the blisshook above */
   auto reportglue = [&](unsigned int n, const unsigned int* aut) {
      blisshook((void*)&data, n, aut);
   };

   /* lambda function to have access to data and terminate the search if maxgenerators are reached */
   auto term = [&]() {
      return (maxgenerators != 0 && data.nperms >= maxgenerators); /* check the number of generators that we have created so far */
   };

   /* start search */
   G->find_automorphisms(stats, reportglue, term);
#else

   /* Older bliss versions do not allow to terminate with a limit on the number of generators unless patched. */
#ifdef BLISS_PATCH_PRESENT
   /* If patch is present, do not use a node limit, but set generator limit. This approach is not very accurate, since
    * it limits the generators considered in bliss and not the ones that we generate (the ones that work on the variable
    * set). */
   G->set_search_limits(0, (unsigned) maxgenerators);
#endif

   /* start search */
   G->find_automorphisms(stats, blisshook, (void*) &data);
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
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations
                                              *   (needed for freeing storage) */
   int***                perms,              /**< pointer to store generators as (nperms x npermvars) matrix */
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of size of group */
   SCIP_Real*            symcodetime         /**< pointer to store the time for symmetry code */
   )
{
   int nvarnodestoadd;
   int first;
   int second;
   int nnodes = 0;
   int nedges = 0;

   assert( scip != NULL );
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

   /* create bliss graph */
   bliss::Graph G(0);

   /* add nodes for variables
    *
    * Variable nodes come first to easily extract the variable permutation.
    * For signed permutations, the first nsymvars nodes correspond to original
    * variables, and the second nsymvars nodes to their negations.
    */
   int nsymvars = SCIPgetSymgraphNVars(graph);
   SYM_SYMTYPE symtype = SCIPgetSymgraphSymtype(graph);
   switch ( symtype )
   {
   case SYM_SYMTYPE_PERM:
      nvarnodestoadd = nsymvars;
      break;
   default:
      assert( symtype == SYM_SYMTYPE_SIGNPERM );
      nvarnodestoadd = 2 * nsymvars;
   }

   for (int v = 0; v < nvarnodestoadd; ++v)
   {
      const int color = SCIPgetSymgraphVarnodeColor(graph, v);

#ifndef NDEBUG
      int node = (int) G.add_vertex((unsigned) color);
      assert( node == v );
#else
      (void) G.add_vertex((unsigned) color);
#endif

      ++nnodes;
   }

   /* add nodes for non-variable nodes */
   int nsymnodes = SCIPgetSymgraphNNodes(graph);
   for (int v = 0; v < nsymnodes; ++v)
   {
      const int color = SCIPgetSymgraphNodeColor(graph, v);

#ifndef NDEBUG
      int node = (int) G.add_vertex((unsigned) color);
      assert( node == nvarnodestoadd + v );
#else
      (void) G.add_vertex((unsigned) color);
#endif

      ++nnodes;
   }

   /* add edges to bliss graph
    *
    * Edges containing neither a variable or constraint node are added immediately.
    * Remaining edges are collected and we group these edges based on their weight.
    */
   const bool groupByConstraints = SCIPgetSymgraphNConsnodes(graph) < SCIPgetSymgraphNVars(graph);
   int nsymedges = SCIPgetSymgraphNEdges(graph);
   int* groupfirsts = NULL;
   int* groupseconds = NULL;
   int* groupcolors = NULL;
   int ngroupedges = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &groupfirsts, nsymedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &groupseconds, nsymedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &groupcolors, nsymedges) );

   for (int e = 0; e < nsymedges; ++e)
   {
      first = SCIPgetSymgraphEdgeFirst(graph, e);
      second = SCIPgetSymgraphEdgeSecond(graph, e);

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
         SYM_NODETYPE comparetype = groupByConstraints ? SYM_NODETYPE_CONS : SYM_NODETYPE_VAR;

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
         /* immediately add edge */
         assert(0 <= first && first < nnodes);
         assert(0 <= second && second < nnodes);

         /* possibly split edge if it is colored */
         if ( ! SCIPhasGraphUniqueEdgetype(graph) && SCIPisSymgraphEdgeColored(graph, e) )
         {
            const int color = SCIPgetSymgraphEdgeColor(graph, e);

            int inter = G.add_vertex((unsigned) color);

            G.add_edge(first, inter);
            G.add_edge(second, inter);

            ++nnodes;
            ++nedges;
         }
         else
            G.add_edge(first, second);
         ++nedges;
      }
   }

   /* possibly add groupable edges */
   if ( ngroupedges > 0 )
   {
      /* sort edges according to their first nodes */
      SCIPsortIntIntInt(groupfirsts, groupseconds, groupcolors, ngroupedges);

      int firstidx = 0;
      int firstnodeidx = groupfirsts[0];
      int naddednodes;
      int naddededges;

      for (int i = 1; i < ngroupedges; ++i)
      {
         /* if a new first node has been found, group the edges of the previous first node; ignoring the last group */
         if ( groupfirsts[i] != firstnodeidx )
         {
            SCIP_CALL( addGroupedEdges(&G, firstnodeidx, &groupseconds[firstidx],
                  &groupcolors[firstidx], i - firstidx, &naddednodes, &naddededges) );

            firstidx = i;
            firstnodeidx = groupfirsts[i];

            nnodes += naddednodes;
            nedges += naddededges;
         }
      }

      /* process the last group */
      SCIP_CALL( addGroupedEdges(&G, firstnodeidx, &groupseconds[firstidx],
            &groupcolors[firstidx], ngroupedges - firstidx, &naddednodes, &naddededges) );

      nnodes += naddednodes;
      nedges += naddededges;
   }

   SCIPfreeBufferArray(scip, &groupcolors);
   SCIPfreeBufferArray(scip, &groupseconds);
   SCIPfreeBufferArray(scip, &groupfirsts);

   /* for signed permutation, also add edges connecting a variable and its negation */
   switch ( SCIPgetSymgraphSymtype(graph) )
   {
   case SYM_SYMTYPE_SIGNPERM:
      for (int v = 0; v < nsymvars; ++v)
         G.add_edge(v, v + nsymvars);
      nedges += nsymvars;
      break;
   default:
      assert( SCIPgetSymgraphSymtype(graph) == SYM_SYMTYPE_PERM );
   }

   assert( (int) G.get_nof_vertices() == nnodes );
   assert( nedges >= SCIPgetSymgraphNEdges(graph) );
   SCIPdebugMsg(scip, "Symmetry detection graph has %d nodes and %d edges.\n", nnodes, nedges);

   /* compute automorphisms */
   SCIP_CALL( computeAutomorphisms(scip, symtype, &G, nsymvars, maxgenerators,
         perms, nperms, nmaxperms, log10groupsize, TRUE, symcodetime) );

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
   int* nvarused1 = NULL;
   int* nvarused2 = NULL;
   int* varlabel = NULL;
   int nusedvars = 0;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( G1 != NULL );
   assert( G2 != NULL );

   /* some simple checks */
   if ( G1->nnodes != G2->nnodes ||  G1->nopnodes != G2->nopnodes || G1->nvalnodes != G2->nvalnodes
      || G1->nconsnodes != G2->nconsnodes || G1->nedges != G2->nedges )
      return FALSE;

   /* check whether the variables used in G1 are the same as in G2 */
   switch ( symtype )
   {
   case SYM_SYMTYPE_PERM:
      nvars = G1->nsymvars;
      break;
   default:
      assert( symtype == SYM_SYMTYPE_SIGNPERM );
      nvars = 2 * G1->nsymvars;
   }
   SCIP_CALL_ABORT( SCIPallocClearBufferArray(scip, &nvarused1, nvars) );
   SCIP_CALL_ABORT( SCIPallocClearBufferArray(scip, &nvarused2, nvars) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &varlabel, nvars) );

   for (i = 0; i < G1->nedges; ++i)
   {
      if ( G1->edgefirst[i] < 0 )
         nvarused1[-G1->edgefirst[i] - 1] += 1;
      if ( G2->edgefirst[i] < 0 )
         nvarused2[-G2->edgefirst[i] - 1] += 1;
      if ( G1->edgesecond[i] < 0 )
         nvarused1[-G1->edgesecond[i] - 1] += 1;
      if ( G2->edgesecond[i] < 0 )
         nvarused2[-G2->edgesecond[i] - 1] += 1;
   }

   for (i = 0; i < nvars; ++i)
   {
      if ( nvarused1[i] != nvarused2[i] )
      {
         SCIPfreeBufferArray(scip, &varlabel);
         SCIPfreeBufferArray(scip, &nvarused2);
         SCIPfreeBufferArray(scip, &nvarused1);

         return FALSE;
      }

      /* relabel variables by restricting to variables used in constraint (or their negation) */
      if ( nvarused1[i] > 0 || nvarused1[i % G1->nsymvars] > 0 )
         varlabel[i] = nusedvars++;
      else
         varlabel[i] = -1;
   }

   /* construct bliss graph containing the (disjoint) union of the two graphs */
   bliss::Graph G(0);

   /* copy of G1 */
   for (i = 0; i < nusedvars; ++i)
      G.add_vertex(i);

   for (i = 0; i < G1->nnodes; ++i)
      G.add_vertex(nusedvars + SCIPgetSymgraphNodeColor(G1, i));

   for (i = 0; i < G1->nedges; ++i)
   {
      int first = G1->edgefirst[i];
      int second = G1->edgesecond[i];

      if ( first < 0 )
         first = varlabel[-first - 1];
      else
         first = nusedvars + first;
      assert( first >= 0 );

      if ( second < 0 )
         second = varlabel[-second - 1];
      else
         second = nusedvars + second;
      assert( second >= 0 );

      if ( SCIPisSymgraphEdgeColored(G1, i) )
      {
         int inter = G.add_vertex(nusedvars + SCIPgetSymgraphEdgeColor(G1, i));
         G.add_edge(first, inter);
         G.add_edge(second, inter);
      }
      else
         G.add_edge(first, second);
   }

   /* in case of signed permutations, also connect variables with their negation */
   switch ( symtype )
   {
   case SYM_SYMTYPE_SIGNPERM:
      for (i = 0; i < G1->nsymvars; ++i)
      {
         if ( nvarused1[i] > 0 || nvarused1[i + G1->nsymvars])
            G.add_edge(varlabel[i], varlabel[i + G1->nsymvars]);
      }
      break;
   default:
      assert( symtype == SYM_SYMTYPE_PERM );
   }

   /* copy of G2 */
   int nodeshift = G.get_nof_vertices();
   for (i = 0; i < nusedvars; ++i)
      G.add_vertex(i);

   for (i = 0; i < G2->nnodes; ++i)
      G.add_vertex(nusedvars + SCIPgetSymgraphNodeColor(G2, i));

   for (i = 0; i < G2->nedges; ++i)
   {
      int first = G2->edgefirst[i];
      int second = G2->edgesecond[i];

      if ( first < 0 )
         first = nodeshift + varlabel[-first - 1];
      else
         first = nodeshift + nusedvars + first;
      assert( first >= 0 );

      if ( second < 0 )
         second = nodeshift + varlabel[-second - 1];
      else
         second = nodeshift + nusedvars + second;
      assert( second >= 0 );

      if ( SCIPisSymgraphEdgeColored(G2, i) )
      {
         int inter = G.add_vertex(nusedvars + SCIPgetSymgraphEdgeColor(G2, i));
         G.add_edge(first, inter);
         G.add_edge(second, inter);
      }
      else
         G.add_edge(first, second);
   }

   /* in case of signed permutations, also connect variables with their negation */
   switch ( symtype )
   {
   case SYM_SYMTYPE_SIGNPERM:
      for (i = 0; i < G2->nsymvars; ++i)
      {
         if ( nvarused2[i] > 0 || nvarused2[i + G2->nsymvars])
            G.add_edge(nodeshift + varlabel[i], nodeshift + varlabel[i + G2->nsymvars]);
      }
      break;
   default:
      assert( symtype == SYM_SYMTYPE_PERM );
   }

   /* compute automorphisms */
   int** perms;
   int nperms;
   int nmaxperms;
   SCIP_Real log10groupsize;
   int n = G.get_nof_vertices();
   int nnodesfromG1 = nusedvars + G1->nnodes;
   SCIP_Real symcodetime = 0.0;

   SCIP_CALL_ABORT( computeAutomorphisms(scip, symtype, &G, n, 0,
         &perms, &nperms, &nmaxperms, &log10groupsize, FALSE, &symcodetime) );

   /* since G1 and G2 are connected and disjoint, they are isomorphic iff there is a permutation
    * mapping a node from G1 to a node of G2
    */
   SCIP_Bool success = FALSE;
   for (int p = 0; p < nperms && ! success; ++p)
   {
      for (i = 0; i < nnodesfromG1; ++i)
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
      SCIPfreeBlockMemoryArray(scip, &perms[p], n);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &perms, nmaxperms);

   SCIPfreeBufferArray(scip, &varlabel);
   SCIPfreeBufferArray(scip, &nvarused2);
   SCIPfreeBufferArray(scip, &nvarused1);

   return success;
}
