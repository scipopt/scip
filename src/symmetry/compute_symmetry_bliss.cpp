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
   int                   npermvars;          /**< number of variables for permutations */
   int                   nperms;             /**< number of permutations */
   int**                 perms;              /**< permutation generators as (nperms x npermvars) matrix */
   int                   nmaxperms;          /**< maximal number of permutations */
   int                   maxgenerators;      /**< maximal number of generators constructed (= 0 if unlimited) */
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
   if ( SCIPallocBlockMemoryArray(data->scip, &p, data->npermvars) != SCIP_OKAY )
      return;

   for (int j = 0; j < data->npermvars; ++j)
   {
      /* convert index of variable-level 0-nodes to variable indices */
      p[j] = (int) aut[j];
      if ( p[j] != j )
         isIdentity = false;
   }

   /* ignore trivial generators, i.e. generators that only permute the constraints */
   if ( isIdentity )
   {
      SCIPfreeBlockMemoryArray(data->scip, &p, data->npermvars);
      return;
   }

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

/** callback function for bliss in case of reflection symmetries */
static
void blisshookReflSym(
   void*                 user_param,         /**< parameter supplied at call to bliss */
   unsigned int          n,                  /**< size of aut vector */
   const unsigned int*   aut                 /**< automorphism */
   )
{
   assert( aut != NULL );
   assert( user_param != NULL );

   BLISS_Data* data = static_cast<BLISS_Data*>(user_param);
   assert( data->scip != NULL );
   assert( 2 * data->npermvars < (int) n );
   assert( data->maxgenerators >= 0);

   /* make sure we do not generate more that maxgenerators many permutations, if the limit in bliss is not available */
   if ( data->maxgenerators != 0 && data->nperms >= data->maxgenerators )
      return;

   /* copy first part of automorphism */
   bool isIdentity = true;
   int* p = 0;
   if ( SCIPallocBlockMemoryArray(data->scip, &p, 2 * data->npermvars) != SCIP_OKAY )
      return;

   for (int j = 0; j < 2 * data->npermvars; ++j)
   {
      /* convert index of variable-level 0-nodes to variable indices */
      p[j] = (int) aut[j];
      if ( p[j] != j )
         isIdentity = false;
   }

   /* ignore trivial generators, i.e. generators that only permute the constraints */
   if ( isIdentity )
   {
      SCIPfreeBlockMemoryArray(data->scip, &p, 2 * data->npermvars);
      return;
   }

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

/** returns the color of a tree node in its respective class */
int colorInClass(
   int                   nodeidx,            /**< index of node in tree */
   SYM_REFLSYMDATA*      reflsymdata,        /**< data for reflection symmetry detection */
   SCIP_Bool             inverse             /**< whether the color of the negated input shall be returned */
   )
{
   int color;

   assert( reflsymdata != NULL );

   switch ( reflsymdata->trees[nodeidx] )
   {
   case SYM_NODETYPE_OPERATOR :
      color = reflsymdata->opscolors[reflsymdata->treemap[nodeidx]];
      break;
   case SYM_NODETYPE_COEF :
      if ( inverse )
         color = reflsymdata->invcoefcolors[reflsymdata->treemap[nodeidx]];
      else
         color = reflsymdata->coefcolors[reflsymdata->treemap[nodeidx]];
      break;
   case SYM_NODETYPE_VAR :
      if ( inverse )
         color = reflsymdata->invvarcolors[reflsymdata->treemap[nodeidx]];
      else
         color = reflsymdata->varcolors[reflsymdata->treemap[nodeidx]];
      break;
   default :
      assert( reflsymdata->trees[nodeidx] == SYM_NODETYPE_VAL );
      color = reflsymdata->valcolors[reflsymdata->treemap[nodeidx]];
   }

   return color;
}

/** returns the types of flips allowed for variable arguments of an operator */
static
SYM_FLIPTYPE getOperatorFliptype(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_REFLSYMDATA*      reflsymdata,        /**< data of CIP */
   int                   opidx               /**< index of operator */
   )
{
   SYM_NODETYPE* trees;
   int* treemap;
   SCIP_Real* treevals;
   SCIP_EXPRHDLR** treeops;
   SCIP_EXPRHDLR* op;

   assert( reflsymdata != NULL );
   assert( 0 <= opidx && opidx < reflsymdata->ntrees );
   assert( reflsymdata->trees[opidx] == SYM_NODETYPE_OPERATOR );

   trees = reflsymdata->trees;
   treemap = reflsymdata->treemap;
   treevals = reflsymdata->treevals;
   treeops = reflsymdata->treeops;

   /* operators without proper arguments do not allow any flips */
   if ( opidx >= reflsymdata->ntrees - 2 )
      return SYM_FLIPTYPE_NONE;

   /* check which types of flips are allowed by an operator */
   op = treeops[treemap[opidx]];
   if ( op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_OR || op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_XOR ||
      op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_BDDISJ || op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_GEQ ||
      op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_OBJ || op == SCIPfindExprhdlr(scip, "sum") )
      return SYM_FLIPTYPE_SHIFT_ODD;
   else if ( treeops[treemap[opidx]] == SCIPfindExprhdlr(scip, "prod") ||
      treeops[treemap[opidx]] == SCIPfindExprhdlr(scip, "sin") )
      return SYM_FLIPTYPE_ODD;
   else if ( op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_SOS1 || op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_SOS2 ||
      op ==  SCIPfindExprhdlr(scip, "abs") || op ==  SCIPfindExprhdlr(scip, "cos") )
      return SYM_FLIPTYPE_EVEN;
   else if ( op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_AND || op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_INDICATOR ||
      op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_EQ )
      return SYM_FLIPTYPE_NONE;
   else if ( treeops[treemap[opidx]] == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_SIGNPOWER )
   {
      if ( trees[opidx + 1] == SYM_NODETYPE_VAL )
      {
         if ( trees[opidx + 2] != SYM_NODETYPE_COEF )
            return SYM_FLIPTYPE_NONE;
         else
            return SYM_FLIPTYPE_EVEN;
      }
      else
      {
         if ( trees[opidx + 1] != SYM_NODETYPE_COEF )
            return SYM_FLIPTYPE_NONE;
         else
            return SYM_FLIPTYPE_EVEN;
      }
   }
   else if ( treeops[treemap[opidx]] == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_POWER )
   {
      /* power 1 yields odd function */
      if ( trees[opidx + 1] != SYM_NODETYPE_VAL )
      {
         /* has not a single variable as argument */
         if ( trees[opidx + 2] != SYM_NODETYPE_COEF )
            return SYM_FLIPTYPE_NONE;
         else
            return SYM_FLIPTYPE_ODD;
      }

      /* check whether exponent leads to even/odd function */
      if ( ! SCIPisIntegral(scip, treevals[treemap[opidx + 1]]) )
         return SYM_FLIPTYPE_NONE;
      else if ( (int) treevals[treemap[opidx + 1]] % 2 == 0 )
         return SYM_FLIPTYPE_EVEN;
      else
         return SYM_FLIPTYPE_ODD;
   }
   else if ( treeops[treemap[opidx]] == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_BIPROD )
      return SYM_FLIPTYPE_BIPROD;

   return SYM_FLIPTYPE_NONE;
}

/** creates the graph for detecting reflection symmetries */
SCIP_RETCODE SCIPcreateReflectionSymmetryDetectionGraph(
   SCIP*                 scip,               /**< SCIP pointer */
   bliss::Graph*         G,                  /**< graph to be constructed */
   SYM_REFLSYMDATA*      reflsymdata,        /**< data of CIP */
   SCIP_Bool&            success             /**< whether the construction was successful */
   )
{
   int nnodes = 0;
   int nusedcolors = 0;

   assert( scip != NULL );
   assert( G != NULL );
   assert( reflsymdata != NULL );

   success = TRUE;

   SCIPdebugMsg(scip, "Add variable nodes to symmetry detection graph.\n");

   /* add nodes for variables */
   for (int v = 0; v < reflsymdata->ntreevars; ++v)
   {
      const int color = reflsymdata->varcolors[v];
      assert( 0 <= color && color < reflsymdata->nuniquevars );

#ifndef NDEBUG
      int node = (int) G->add_vertex((unsigned) color);
      assert( node == v );
#else
      (void) G->add_vertex((unsigned) color);
#endif

      ++nnodes;
   }

   /* add nodes for negated variables and add edges to the original variable*/
   for (int v = 0; v < reflsymdata->ntreevars; ++v)
   {
      const int color = reflsymdata->invvarcolors[v];
      assert( 0 <= color && color < reflsymdata->nuniquevars );

      int node = (int) G->add_vertex((unsigned) color);
      assert( node == reflsymdata->ntreevars + v );

      G->add_edge(node, node - reflsymdata->ntreevars);

      ++nnodes;
   }
   nusedcolors = reflsymdata->nuniquevars;

   /* add nodes for rhs coefficients of constraints */
   for (int c = 0; c < reflsymdata->nuniquerhs; ++c)
   {
#ifndef NDEBUG
      int node = (int) G->add_vertex((unsigned) (nusedcolors + c));
      assert( node == 2 * reflsymdata->ntreevars + c );
#else
      (void) G->add_vertex((unsigned) (nusedcolors + c));
#endif
   }
   assert( (int) G->get_nof_vertices() == 2 * reflsymdata->ntreevars + reflsymdata->nuniquerhs );
   nusedcolors += reflsymdata->nuniquerhs;

   /* starting positions of different types of nodes */
   int rhsstart = 2 * reflsymdata->ntreevars;

   /* starting positions of colors of different type of nodes */
   int rhscolstart = reflsymdata->nuniquevars;
   int opcolstart = rhscolstart + reflsymdata->nuniquerhs;
   int coefcolstart = opcolstart + reflsymdata->nuniqueops;
   int valcolstart = coefcolstart + reflsymdata->nuniquecoefs;

   /* add nodes and edges for constraints */
   for (int c = 0; c < reflsymdata->ntreerhs; ++c)
   {
      /* keep track of path to root node, possible since tree is stored in DFS fashion */
      std::vector<unsigned> pathtoroot;
      std::vector<unsigned> pathtorootidx;
      std::vector<int> fliptypes;

      /* iterate through tree and create nodes and edges for constraint */
      for (int i = reflsymdata->treebegins[c]; i < reflsymdata->treebegins[c + 1]; ++i)
      {
         /* remove elements from the path that are not predecessors of the current node */
         while ( pathtoroot.size() > 0 && (int) pathtorootidx.back() > reflsymdata->treeparentidx[i] )
         {
            pathtorootidx.pop_back();
            pathtoroot.pop_back();
            fliptypes.pop_back();
         }

         if ( reflsymdata->treeparentidx[i] == -1 )
         {
            /* add root node and connect it with corresponding rhs node */
            assert( reflsymdata->trees[i] == SYM_NODETYPE_OPERATOR );

            unsigned node = G->add_vertex((unsigned) opcolstart + colorInClass(i, reflsymdata, FALSE));
            G->add_edge(node, (unsigned) (rhsstart + reflsymdata->rhscolors[c]));

            pathtoroot.push_back(node);
            pathtorootidx.push_back(i);
            fliptypes.push_back(getOperatorFliptype(scip, reflsymdata, i));
         }
         else if ( reflsymdata->trees[i] == SYM_NODETYPE_OPERATOR )
         {
            /* create node for the operator and connect it with its parent node */
            unsigned node = G->add_vertex((unsigned) opcolstart + colorInClass(i, reflsymdata, FALSE));
            G->add_edge(node, pathtoroot.back());

            SYM_FLIPTYPE fliptype = getOperatorFliptype(scip, reflsymdata, i);
            if ( fliptype != SYM_FLIPTYPE_BIPROD )
            {
               pathtoroot.push_back(node);
               pathtorootidx.push_back(i);
               fliptypes.push_back(fliptype);
            }
            else
            {
               /* already create nodes/edges for values, coefficients, and variables */
               assert( i + 1 < reflsymdata->treebegins[c + 1] - 1 );
               ++i;
               if ( reflsymdata->trees[i] == SYM_NODETYPE_VAL )
               {
                  unsigned valnode = G->add_vertex((unsigned) valcolstart + colorInClass(i, reflsymdata, FALSE));
                  G->add_edge(valnode, node);
                  ++i;
               }

               assert( i < reflsymdata->treebegins[c + 1] - 3 );
               assert( reflsymdata->trees[i] == SYM_NODETYPE_COEF );
               assert( reflsymdata->trees[i + 1] == SYM_NODETYPE_VAR );
               assert( reflsymdata->trees[i + 2] == SYM_NODETYPE_COEF );
               assert( reflsymdata->trees[i + 3] == SYM_NODETYPE_VAR );
               assert( colorInClass(i, reflsymdata, FALSE) == colorInClass(i + 2, reflsymdata, FALSE) );

               unsigned coef1node = G->add_vertex((unsigned) coefcolstart + colorInClass(i, reflsymdata, FALSE));
               unsigned var1node = reflsymdata->treevaridx[reflsymdata->treemap[i + 1]];
               unsigned invvar1node = var1node + reflsymdata->ntreevars;
               unsigned coef2node = G->add_vertex((unsigned) coefcolstart + colorInClass(i + 2, reflsymdata, FALSE));
               unsigned var2node = reflsymdata->treevaridx[reflsymdata->treemap[i + 3]];
               unsigned invvar2node = var2node + reflsymdata->ntreevars;
               G->add_edge(coef1node, node);
               G->add_edge(coef2node, node);
               G->add_edge(var1node, coef1node);
               G->add_edge(var2node, coef1node);
               G->add_edge(invvar1node, coef2node);
               G->add_edge(invvar2node, coef2node);
               i += 3;
            }
         }
         else if ( reflsymdata->trees[i] == SYM_NODETYPE_COEF )
         {
            assert( i < reflsymdata->treebegins[c + 1] - 1 );
            assert( reflsymdata->trees[i + 1] == SYM_NODETYPE_VAR );

            /* create node for coefficient (and its negation if flips are allowed) */
            unsigned coefnode = G->add_vertex((unsigned) coefcolstart + colorInClass(i, reflsymdata, FALSE));
            unsigned varnode = reflsymdata->treevaridx[reflsymdata->treemap[i + 1]];

            G->add_edge(coefnode, pathtoroot.back());
            G->add_edge(coefnode, varnode);

            switch ( fliptypes.back() )
            {
            case SYM_FLIPTYPE_EVEN :
            {
               unsigned invvarnode = varnode + reflsymdata->ntreevars;
               unsigned invcoefnode = G->add_vertex((unsigned) coefcolstart + colorInClass(i, reflsymdata, TRUE));
               G->add_edge(coefnode, invcoefnode);
               G->add_edge(invcoefnode, pathtoroot.back());
               G->add_edge(invcoefnode, invvarnode);
               G->add_edge(coefnode, invvarnode);
               G->add_edge(invcoefnode, varnode);
               break;
            }
            case SYM_FLIPTYPE_SHIFT_ODD :
            case SYM_FLIPTYPE_ODD :
            {
               unsigned invcoefnode = G->add_vertex((unsigned) coefcolstart + colorInClass(i, reflsymdata, TRUE));
               unsigned invvarnode = varnode + reflsymdata->ntreevars;
               G->add_edge(coefnode, invcoefnode);
               G->add_edge(invcoefnode, pathtoroot.back());
               G->add_edge(invcoefnode, invvarnode);
               break;
            }
            default :
               assert ( fliptypes.back() == SYM_FLIPTYPE_NONE );
            }

            /* we have already taken the variable into account */
            ++i;
         }
         else if ( reflsymdata->trees[i] == SYM_NODETYPE_VAR )
         {
            SCIPinfoMessage(scip, NULL, "ERROR: terminate construction of symmetry detection graph, each variable should have a coefficient.\n");

            success = FALSE;
         }
         else
         {
            assert( reflsymdata->trees[i] == SYM_NODETYPE_VAL );

            /* create node for value and connect it with its parent node */
            unsigned node = G->add_vertex((unsigned) valcolstart + colorInClass(i, reflsymdata, FALSE));

            G->add_edge(node, pathtoroot.back());
         }
      }
   }

   return SCIP_OKAY;
}

/** return whether symmetry can be computed */
SCIP_Bool SYMcanComputeSymmetry(void)
{
   return TRUE;
}

char*
initStaticBlissName( );

static char* blissname = initStaticBlissName();

char*
initStaticBlissName( )
{
   blissname = new char[100];
#ifdef BLISS_PATCH_PRESENT
   (void) snprintf(blissname, 100, "bliss %sp", bliss::version);
#else
   (void) snprintf(blissname, 100, "bliss %s", bliss::version);
#endif
   return blissname;
}


/** return name of external program used to compute generators */
const char* SYMsymmetryGetName(void)
{
   return blissname;
}

/** return description of external program used to compute generators */
const char* SYMsymmetryGetDesc(void)
{
   return "Computing Graph Automorphism Groups by T. Junttila and P. Kaski (www.tcs.hut.fi/Software/bliss/)";
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
   int*                  naddednodes,        /**< buffer to hold number of nodes added to G */
   int*                  naddededges         /**< buffer to hold number of edges added to G */
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
   bliss::Graph*         G,                  /**< pointer to graph for that automorphisms are computed */
   int                   nsymvars,           /**< number of variables encoded in graph */
   int                   maxgenerators,      /**< maximum number of generators to be constructed (=0 if unlimited) */
   int***                perms,              /**< pointer to store generators as (nperms x npermvars) matrix */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations
                                              *   (needed for freeing storage) */
   SCIP_Real*            log10groupsize      /**< pointer to store log10 of size of group */
   )
{
   assert( scip != NULL );
   assert( G != NULL );
   assert( perms != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( log10groupsize != NULL );

   bliss::Stats stats;
   BLISS_Data data;
   data.scip = scip;
   data.npermvars = nsymvars;
   data.nperms = 0;
   data.nmaxperms = 0;
   data.maxgenerators = maxgenerators;
   data.perms = NULL;

   /* Prefer splitting partition cells corresponding to variables over those corresponding
    * to inequalities. This is because we are only interested in the action
    * of the automorphism group on the variables, and we need a base for this action */
   G->set_splitting_heuristic(bliss::Graph::shs_f);
   /* disable component recursion as advised by Tommi Junttila from bliss */
   G->set_component_recursion(false);

   /* do not use a node limit, but set generator limit */
#ifdef BLISS_PATCH_PRESENT
   G->set_search_limits(0, (unsigned) maxgenerators);
#endif

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   /* lambda function to have access to data and pass it to the blisshook above */
   auto reportglue = [&](unsigned int n, const unsigned int* aut) {
      blisshook((void*)&data, n, aut);
   };

   /* lambda function to have access to stats and terminate the search if maxgenerators are reached */
   long unsigned int terminatesearch = INT_MAX;
   if ( maxgenerators != 0 )
      terminatesearch = (long unsigned int) maxgenerators;
   auto term = [&]() {
      return (stats.get_nof_generators() >= terminatesearch);
   };

   /* start search */
   G->find_automorphisms(stats, reportglue, term);
#else
   /* start search */
   G->find_automorphisms(stats, blisshook, (void*) &data);
#endif


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
   SCIP_Real*            log10groupsize      /**< pointer to store log10 of size of group */
   )
{
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

   /* init */
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *log10groupsize = 0;

   /* create bliss graph */
   bliss::Graph G(0);

   /* add nodes for variables
    *
    * variable nodes come first to easily extract the variable permutation */
   int nsymvars = SCIPgetSymgraphNVars(graph);
   for (int v = 0; v < nsymvars; ++v)
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
      assert( node == nsymvars + v );
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
   int* groupfirsts;
   int* groupseconds;
   int* groupcolors;
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
         first += nsymvars;
      if ( second < 0 )
         second = -second - 1;
      else
         second += nsymvars;

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

   assert( (int) G.get_nof_vertices() == nnodes );
   SCIPdebugMsg(scip, "Symmetry detection graph has %d nodes and %d edges.\n", nnodes, nedges);

   /* compute automorphisms */
   SCIP_CALL( computeAutomorphisms(scip, &G, nsymvars, maxgenerators, perms, nperms, nmaxperms, log10groupsize) );

   return SCIP_OKAY;
}

/** returns whether two given graphs are identical */
SCIP_Bool SYMcheckGraphsAreIdentical(
   SCIP*                 scip,               /**< SCIP pointer */
   SYM_GRAPH*            G1,                 /**< first graph */
   SYM_GRAPH*            G2                  /**< second graph */
   )
{
   int* nvarused1;
   int* nvarused2;
   int* varlabel;
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
   nvars = G1->nsymvars;
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nvarused1, nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nvarused2, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varlabel, nvars) );

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

      /* relabel variables by restricting to variables used in constraint */
      if ( nvarused1[i] > 0 )
         varlabel[i] = nusedvars;
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
         int inter = G.add_vertex(SCIPgetSymgraphEdgeColor(G1, i));
         G.add_edge(first, inter);
         G.add_edge(second, inter);
      }
      else
         G.add_edge(first, second);

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
         int inter = G.add_vertex(SCIPgetSymgraphEdgeColor(G2, i));
         G.add_edge(first, inter);
         G.add_edge(second, inter);
      }
      else
         G.add_edge(first, second);
   }

   /* compute automorphisms */
   int** perms;
   int nperms;
   int nmaxperms;
   SCIP_Real log10groupsize;
   int n = G.get_nof_vertices();
   int nnodesfromG1 = nusedvars + G1->nnodes;

   SCIP_CALL( computeAutomorphisms(scip, &G, n, 0,
         &perms, &nperms, &nmaxperms, &log10groupsize) );

   /* since G1 and G2 are connected and disjoint, they are isomorphic iff there is a permutation
    * mapping a node from G1 to a node of G2
    */
   SCIP_Bool success = FALSE;
   for (int p = 0; p < nperms && ! success; ++p)
   {
      for (i = 0; i < nnodesfromG1; ++i)
      {
         if ( perms[p][i] >= nnodesfromG1 )
            success = TRUE;
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

/** compute generators of reflection symmetry group */
SCIP_RETCODE SYMcomputeReflectionSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP pointer */
   SYM_REFLSYMDATA*      reflsymdata,        /**< data for computing reflection symmetries */
   int                   maxgenerators,      /**< maximum number of generators to be computed */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   SCIP_Real*            log10groupsize      /**< pointer to store size of group */
   )
{
   SCIP_Bool success;

   assert( scip != NULL );
   assert( reflsymdata != NULL );

   /* create bliss graph */
   bliss::Graph G(0);

   /* create nodes corresponding to variables */
   SCIP_CALL( SCIPcreateReflectionSymmetryDetectionGraph(scip, &G, reflsymdata, success) );

   if ( !success )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Graph construction failed.\n");
      return SCIP_OKAY;
   }

#ifdef SCIP_OUTPUT
   G.write_dot("debug.dot");
#endif

   SCIPdebugMsg(scip, "Symmetry detection graph has %u nodes.\n", G.get_nof_vertices());

   /* compute automorphisms */
   bliss::Stats stats;
   BLISS_Data data;
   data.scip = scip;
   data.npermvars = reflsymdata->ntreevars;
   data.nperms = 0;
   data.nmaxperms = 0;
   data.maxgenerators = maxgenerators;
   data.perms = NULL;

   /* Prefer splitting partition cells corresponding to variables over those corresponding
    * to inequalities. This is because we are only interested in the action
    * of the automorphism group on the variables, and we need a base for this action */
   G.set_splitting_heuristic(bliss::Graph::shs_f);
   /* disable component recursion as advised by Tommi Junttila from bliss */
   G.set_component_recursion(false);

   /* do not use a node limit, but set generator limit */
#ifdef BLISS_PATCH_PRESENT
   G.set_search_limits(0, (unsigned) maxgenerators);
#endif

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   /* lambda function to have access to data and pass it to the blisshook above */
   auto reportglue = [&](unsigned int n, const unsigned int* aut) {
      blisshookReflSym((void*)&data, n, aut);
   };

   /* lambda function to have access to stats and terminate the search if maxgenerators are reached */
   auto term = [&]() {
      return (stats.get_nof_generators() >= (long unsigned int) maxgenerators);
   };

   /* start search */
   G.find_automorphisms(stats, reportglue, term);
#else
   /* start search */
   G.find_automorphisms(stats, blisshook, (void*) &data);
#endif

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
   }

   /* determine log10 of symmetry group size */
   *log10groupsize = (SCIP_Real) log10l(stats.get_group_size_approx());


   printf("Found %d generators, log10 group size %.1f\n", data.nperms, (SCIP_Real) log10l(stats.get_group_size_approx()));
   for (int p = 0; p < data.nperms; ++p)
   {
      SCIP_Bool hasreflection = FALSE;

      for (int v = 0; v < data.npermvars; ++v)
      {
         if ( data.perms[p][v] >= data.npermvars )
         {
            hasreflection = TRUE;
            break;
         }
      }

      if ( hasreflection )
      {
         // print permutation in cycle representation
         SCIP_Bool* covered;

         SCIP_CALL( SCIPallocClearBufferArray(scip, &covered, data.npermvars) );

         printf("perm %d: ", p);

         for (int v = 0; v < data.npermvars; ++v)
         {
            if ( covered[v] )
               continue;

            if ( data.perms[p][v] == v )
               continue;

            printf("(%d", v);
            covered[v] = TRUE;

            int w = data.perms[p][v];
            while ( w != v )
            {
               if ( w >= data.npermvars )
               {
                  printf(" -%d", w);
                  covered[w - data.npermvars] = TRUE;
               }
               else
               {
                  printf(" %d", w);
                  covered[w] = TRUE;
               }

               w = data.perms[p][w];
            }
            printf(")");
         }

         printf("\n");
         for (int v = 0; v < data.npermvars; ++v)
            covered[v] = FALSE;

         for (int v = 0; v < data.npermvars; ++v)
         {
            if ( covered[v] )
               continue;

            if ( data.perms[p][v] == v )
               continue;

            printf("(%s", SCIPvarGetName(SCIPgetVars(scip)[v]));
            covered[v] = TRUE;

            int w = data.perms[p][v];
            while ( w != v )
            {
               if ( w >= data.npermvars )
               {
                  printf(" -%s", SCIPvarGetName(SCIPgetVars(scip)[w - data.npermvars]));
                  covered[w - data.npermvars] = TRUE;
               }
               else
               {
                  printf(" %s", SCIPvarGetName(SCIPgetVars(scip)[w]));
                  covered[w] = TRUE;
               }

               w = data.perms[p][w];
            }
            printf(")");
         }

         SCIPfreeBufferArray(scip, &covered);
         printf("\n\n");
      }
   }

   return SCIP_OKAY;
}
