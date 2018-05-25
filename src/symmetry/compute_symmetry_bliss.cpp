/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   compute_symmetry_bliss.cpp
 * @brief  interface for symmetry computations to bliss
 * @author Marc Pfetsch
 * @author Thomas Rehn
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "compute_symmetry.h"

/* include bliss graph */
#include <bliss/defs.hh>
#include <bliss/graph.hh>

#include <vector>
#include <list>
#include <math.h>

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

   BLISS_Data* data = (BLISS_Data*) user_param;
   assert( data->scip != NULL );
   assert( data->perms != NULL );
   assert( data->npermvars < (int) n );
   assert( data->maxgenerators >= 0);

   /* make sure we do not generate more that maxgenerators many permutations, if the limit is bliss is not available */
   if ( data->maxgenerators != 0 && data->nperms >= data->maxgenerators )
      return;

   /* check whether we need to resize */
   if ( data->nperms >= data->nmaxperms )
   {
      int newsize = SCIPcalcMemGrowSize(data->scip, data->nperms);
      SCIP_RETCODE retcode = SCIPreallocBlockMemoryArray(data->scip, &data->perms, data->nmaxperms, newsize);
      if ( retcode != SCIP_OKAY )
         return;
      data->nmaxperms = newsize;
   }

   /* copy first part of automorphism */
   bool isIdentity = true;
   int* p = 0;
   SCIP_RETCODE retcode = SCIPallocBlockMemoryArray(data->scip, &p, data->npermvars);
   if ( retcode != SCIP_OKAY )
      return;
   for (int j = 0; j < data->npermvars; ++j)
   {
      /* convert index of variable-level 0-nodes to variable indices */
      p[j] = (int) aut[j];
      if ( p[j] != j )
         isIdentity = false;
   }

   /*  ignore trivial generators, i.e. generators that only permute the constraints */
   if ( isIdentity )
   {
      SCIPfreeBlockMemoryArray(data->scip, &p, data->npermvars);
   }
   else
      data->perms[data->nperms++] = p;
}


/** Construct colored graph for symmetry computations
 *
 *  Construct bipartite graph:
 *  - Each variable gets a different node.
 *  - Each constraint gets a different node.
 *  - Each matrix coefficient gets a different node that is conntected to the two nodes
 *    corresponding to the constraint and variable.
 *
 *  Each different variable, rhs, and matrix coefficient type gets a different color that is
 *  attached to the corresponding entries.
 */
static
SCIP_RETCODE fillGraphByColoredCoefficients(
   SCIP*                 scip,               /**< SCIP instance */
   bliss::Graph*         G,                  /**< Graph to be constructed */
   SYM_MATRIXDATA*       matrixdata,         /**< data for MIP matrix */
   int&                  nnodes,             /**< number of nodes in graph */
   int&                  nedges,             /**< number of edges in graph */
   SCIP_Bool&            success             /**< whether the construction was successful */
   )
{
   SCIPdebugMsg(scip, "Building graph with colored coefficient nodes.\n");

   nnodes = 0;
   nedges = 0;
   success = FALSE;

   /* add nodes for variables */
   for (int v = 0; v < matrixdata->npermvars; ++v)
   {
      const int color = matrixdata->permvarcolors[v];
      assert( 0 <= color && color < matrixdata->nuniquevars );

#ifndef NDEBUG
      int node = (int) G->add_vertex((unsigned) color);
      assert( node == v );
#else
      (void) G->add_vertex((unsigned) color);
#endif

      ++nnodes;
   }
   assert( (int) G->get_nof_vertices() == matrixdata->npermvars );

   /* add nodes for rhs of constraints */
   for (int c = 0; c < matrixdata->nrhscoef; ++c)
   {
      const int color = matrixdata->rhscoefcolors[c];
      assert( 0 <= color && color < matrixdata->nuniquerhs );

#ifndef NDEBUG
      int node = (int) G->add_vertex((unsigned) (matrixdata->nuniquevars + color));
      assert( node == matrixdata->npermvars + c );
#else
      (void) G->add_vertex((unsigned) (matrixdata->nuniquevars + color));
#endif

      ++nnodes;
   }
   assert( (int) G->get_nof_vertices() == matrixdata->npermvars + matrixdata->nrhscoef );

   /* Grouping of nodes depends on the number of nodes in the bipartite graph class.
    * If there are more variables than constraints, we group by constraints.
    * That is, given several variable nodes which are incident to one constraint node by the same color,
    * we join these variable nodes to the constaint node by only one intermediate node.
    */
   const bool groupByConstraints = matrixdata->nrhscoef < matrixdata->npermvars;
   if ( groupByConstraints )
      SCIPdebugMsg(scip, "Group intermediate nodes by constraints.\n");
   else
      SCIPdebugMsg(scip, "Group intermediate nodes by variables.\n");

   /* "colored" edges based on all matrix coefficients - loop through ordered matrix coefficients */
   int nusedcolors = matrixdata->nuniquevars + matrixdata->nuniquerhs;

   int ninternodes;
   if ( groupByConstraints )
      ninternodes = matrixdata->nrhscoef;
   else
      ninternodes = matrixdata->npermvars;

   int* internodes = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &internodes, ninternodes) ); /*lint !e530*/
   for (int l = 0; l < ninternodes; ++l)
      internodes[l] = -1;

   /* We pass through the matrix coeficients, grouped by color, i.e., different coefficients. If the coeffients appear
    * in the same row or column, it suffices to only generate a single node (depending on groupByConstraints). We store
    * this node in the array internodes. In order to avoid reinitialization, we store the node number with increasing
    * numbers for each color. The smallest number for the current color is stored in firstcolornodenumber. */
   int oldcolor = -1;
#ifndef NDEBUG
   SCIP_Real oldcoef = SCIP_INVALID;
#endif
   int firstcolornodenumber = -1;
   for (int j = 0; j < matrixdata->nmatcoef; ++j)
   {
      int idx = matrixdata->matidx[j];
      assert( 0 <= idx && idx < matrixdata->nmatcoef );

      /* find color corresponding to matrix coefficient */
      const int color = matrixdata->matcoefcolors[idx];
      assert( 0 <= color && color < matrixdata->nuniquemat );

      assert( 0 <= matrixdata->matrhsidx[idx] && matrixdata->matrhsidx[idx] < matrixdata->nrhscoef );
      assert( 0 <= matrixdata->matvaridx[idx] && matrixdata->matvaridx[idx] < matrixdata->npermvars );

      const int rhsnode = matrixdata->npermvars + matrixdata->matrhsidx[idx];
      const int varnode = matrixdata->matvaridx[idx];
      assert( matrixdata->npermvars <= rhsnode && rhsnode < matrixdata->npermvars + matrixdata->nrhscoef );
      assert( rhsnode < (int) G->get_nof_vertices() );
      assert( varnode < (int) G->get_nof_vertices() );

      /* if we have only one color, we do not need intermediate nodes */
      if ( matrixdata->nuniquemat == 1 )
      {
         G->add_edge((unsigned) varnode, (unsigned) rhsnode);
         ++nedges;
      }
      else
      {
         /* if new group of coefficients has been reached */
         if ( color != oldcolor )
         {
            assert( ! SCIPisEQ(scip, oldcoef, matrixdata->matcoef[idx]) );
            oldcolor = color;
            firstcolornodenumber = nnodes;
#ifndef NDEBUG
            oldcoef = matrixdata->matcoef[idx];
#endif
         }
         else
            assert( SCIPisEQ(scip, oldcoef, matrixdata->matcoef[idx]) );

         int varrhsidx;
         if ( groupByConstraints )
            varrhsidx = matrixdata->matrhsidx[idx];
         else
            varrhsidx = matrixdata->matvaridx[idx];
         assert( 0 <= varrhsidx && varrhsidx < ninternodes );

         if ( internodes[varrhsidx] < firstcolornodenumber )
         {
            internodes[varrhsidx] = (int) G->add_vertex((unsigned) (nusedcolors + color));
            ++nnodes;
         }
         assert( internodes[varrhsidx] >= matrixdata->npermvars + matrixdata->nrhscoef );
         assert( internodes[varrhsidx] >= firstcolornodenumber );

         /* determine whether graph would be too large for bliss (can only handle int) */
         if ( nnodes >= INT_MAX/2 )
         {
            SCIPfreeBufferArray(scip, &internodes);
            return SCIP_OKAY;
         }

         G->add_edge((unsigned) varnode, (unsigned) internodes[varrhsidx]);
         G->add_edge((unsigned) rhsnode, (unsigned) internodes[varrhsidx]);
         nedges += 2;
      }
   }
   SCIPfreeBufferArray(scip, &internodes);

   success = TRUE; /*lint !e838*/

   return SCIP_OKAY;
}


/** return whether symmetry can be computed */
SCIP_Bool SYMcanComputeSymmetry(void)
{
   return TRUE;
}

/** static variable for holding the name of bliss */
static char blissname[100];

/** return name of external program used to compute generators */
const char* SYMsymmetryGetName(void)
{
#ifdef BLISS_PATCH_PRESENT
   sprintf(blissname, "bliss %sp", bliss::version);
#else
   sprintf(blissname, "bliss %s", bliss::version);
#endif
   return blissname;
}

/** return description of external program used to compute generators */
const char* SYMsymmetryGetDesc(void)
{
   return "Computing Graph Automorphism Groups by T. Junttila and P. Kaski (http://www.tcs.hut.fi/Software/bliss/)";
}

/** compute generators of symmetry group */
SCIP_RETCODE SYMcomputeSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   SYM_MATRIXDATA*       matrixdata,         /**< data for MIP matrix */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   SCIP_Real*            log10groupsize      /**< pointer to store size of group */
   )
{
   assert( scip != NULL );
   assert( matrixdata != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( perms != NULL );
   assert( log10groupsize != NULL );
   assert( maxgenerators >= 0 );

   /* init */
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *log10groupsize = 0;

   int nnodes = 0;
   int nedges = 0;

   /* create bliss graph */
   bliss::Graph G(0);

   SCIP_Bool success = FALSE;
   SCIP_CALL( fillGraphByColoredCoefficients(scip, &G, matrixdata, nnodes, nedges, success) );
   if ( ! success )
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
   data.npermvars = matrixdata->npermvars;
   data.nperms = 0;
   data.nmaxperms = 100 * matrixdata->npermvars;
   data.maxgenerators = maxgenerators;
   data.perms = NULL;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &data.perms, data.nmaxperms) );

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

   /* start search */
   G.find_automorphisms(stats, blisshook, (void*) &data);
#ifdef SCIP_OUTPUT
   (void) stats.print(stdout);
#endif

   /* prepare return values */
   *perms = data.perms;
   *nperms = data.nperms;
   *nmaxperms = data.nmaxperms;

   /* determine log10 of symmetry group size */
   *log10groupsize = (SCIP_Real) log10l(stats.get_group_size_approx());

   return SCIP_OKAY;
}
