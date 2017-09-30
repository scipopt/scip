/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define SCIP_OUTPUT
/**@file   compute_symmetry_bliss.cpp
 * @brief  interface for symmetry computations to bliss
 * @author Marc Pfetsch
 * @author Thomas Rehn
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "compute_symmetry.h"

/* include bliss graph */
#include <graph.hh>

#include <vector>
#include <list>
#include <map>

/** struct for bliss callback */
struct BLISS_Data
{
   SCIP*                 scip;               /**< SCIP pointer */
   int                   npermvars;          /**< number of variables for permutations */
   int                   nperms;             /**< number of permutations */
   int**                 perms;              /**< permutation generators as (nperms x npermvars) matrix */
   int                   nmaxperms;          /**< maximal number of permutations */
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
   SCIP_Bool             local,              /**< Use local variable bounds? */
   int                   npermvars,          /**< number of variables for permutations */
   SCIP_VAR**            permvars,           /**< variables on which permutations act */
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
   for (int v = 0; v < npermvars; ++v)
   {
      SCIP_VAR* var = permvars[v];
      assert( var != 0 );

      SYM_VARTYPE vt;
      vt.obj = SCIPvarGetObj(var);
      if ( local )
      {
         vt.lb = SCIPvarGetLbLocal(var);
         vt.ub = SCIPvarGetUbLocal(var);
      }
      else
      {
         vt.lb = SCIPvarGetLbGlobal(var);
         vt.ub = SCIPvarGetUbGlobal(var);
      }
      vt.type = SCIPvarGetType(var);

      assert( SCIPhashtableExists(matrixdata->vartypemap, (void*) &vt) );
      SYM_VARTYPE* vtr = (SYM_VARTYPE*) SCIPhashtableRetrieve(matrixdata->vartypemap, (void*) &vt);
      const int color = vtr->color;
      assert( color < matrixdata->nuniquevars );

      (void) G->add_vertex((unsigned) color);
      ++nnodes;
   }
   assert( (int) G->get_nof_vertices() == npermvars );

   /* store nodes corresponding to rhs (we need the original order of rows) */
   int* rhsnodemap = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rhsnodemap, matrixdata->nrhscoef) );
#ifndef NDEBUG
   for (int j = 0; j < matrixdata->nrhscoef; ++j)
      rhsnodemap[j] = matrixdata->nrhscoef + 1;
#endif

   /* add nodes for rhs of constraints */
   for (int c = 0; c < matrixdata->nrhscoef; ++c)
   {
      int idx = matrixdata->rhsidx[c];
      assert( 0 <= idx && idx < matrixdata->nrhscoef );

      SYM_RHSTYPE rt;
      rt.val = matrixdata->rhscoef[idx];
      rt.sense = matrixdata->rhssense[idx];
      assert( SCIPhashtableExists(matrixdata->rhstypemap, (void*) &rt) );

      SYM_RHSTYPE* rtr = (SYM_RHSTYPE*) SCIPhashtableRetrieve(matrixdata->rhstypemap, (void*) &rt);
      const int color = rtr->color;
      assert( color < matrixdata->nuniquerhs );

      int node = G->add_vertex((unsigned) (matrixdata->nuniquevars + color));
      assert( node == npermvars + c );
      rhsnodemap[idx] = node;
      ++nnodes;
   }
   assert( (int) G->get_nof_vertices() == npermvars + matrixdata->nrhscoef );

   typedef std::pair<int, int> InterPair;
   typedef std::map<InterPair, int> IntermediatesMap;
   IntermediatesMap groupedIntermediateNodes;

   /* Grouping of nodes depends on the number of nodes in the bipartite graph class.
    * If there are more variables than constraints, we group by constraints.
    * That is, given several variable nodes which are incident to one constraint node by the same color,
    * we join these variable nodes to the constaint node by only one intermediate node.
    */
   const bool groupByConstraints = matrixdata->nrhscoef < static_cast<int>(npermvars);
   if (groupByConstraints)
      SCIPdebugMsg(scip, "Group intermediate nodes by constraints.\n");
   else
      SCIPdebugMsg(scip, "Group intermediate nodes by variables.\n");

   /* "colored" edges based on all matrix coefficients */
   int nusedcolors = matrixdata->nuniquevars + matrixdata->nuniquerhs;
   for (int j = 0; j < matrixdata->nmatcoef; ++j)
   {
      int idx = matrixdata->matidx[j];
      assert( 0 <= idx && idx < matrixdata->nmatcoef );

      /* find color corresponding to matrix coefficient */
      SCIP_Real val = matrixdata->matcoef[j];

      SYM_MATTYPE mt;
      mt.val = val;
      assert( SCIPhashtableExists(matrixdata->mattypemap, (void*) &mt) );

      SYM_MATTYPE* mtr = (SYM_MATTYPE*) SCIPhashtableRetrieve(matrixdata->mattypemap, (void*) &mt);
      const int color = mtr->color;
      assert( color < matrixdata->nuniquemat );

      assert( matrixdata->matrhsidx[idx] < matrixdata->nrhscoef );
      assert( matrixdata->matvaridx[idx] < npermvars );

      const int rhsnode = rhsnodemap[matrixdata->matrhsidx[idx]];
      const int varnode = matrixdata->matvaridx[idx];
      assert( rhsnode < (int) G->get_nof_vertices() );
      assert( varnode < (int) G->get_nof_vertices() );
      assert( npermvars <= rhsnode && rhsnode < npermvars + matrixdata->nrhscoef );

      /* if we have only one color, we do not need intermediate nodes */
      if ( matrixdata->nuniquemat == 1)
      {
         G->add_edge((unsigned) varnode, (unsigned) rhsnode);
         ++nedges;
      }
      else
      {
         InterPair key;
         if (groupByConstraints)
            key = std::make_pair(matrixdata->matrhsidx[idx], color);
         else
            key = std::make_pair(matrixdata->matvaridx[idx], color);

         IntermediatesMap::const_iterator keyIt = groupedIntermediateNodes.find(key);
         int intermediatenode = 0;
         if (keyIt == groupedIntermediateNodes.end())
         {
            intermediatenode = G->add_vertex((unsigned) (nusedcolors + color));
            groupedIntermediateNodes[ key ] = intermediatenode;
            ++nnodes;
         }
         else
         {
            intermediatenode = (*keyIt).second;
         }
         assert( intermediatenode >= static_cast<int>(npermvars + matrixdata->nrhscoef) );

         /* determine whether graph would be too large for bliss (can only handle int) */
         if ( intermediatenode >= INT_MAX )
         {
            SCIPfreeBlockMemoryArray(scip, &rhsnodemap, matrixdata->nrhscoef);
            return SCIP_OKAY;
         }
         G->add_edge((unsigned) varnode, (unsigned) intermediatenode);
         G->add_edge((unsigned) rhsnode, (unsigned) intermediatenode);
         nedges += 2;
      }
   }

   SCIPfreeBlockMemoryArray(scip, &rhsnodemap, matrixdata->nrhscoef);

   success = TRUE;

   return SCIP_OKAY;
}


/** return whether symmetry can be computed */
SCIP_Bool SYMcanComputeSymmetry(void)
{
   return TRUE;
}

/** return name of external program used to compute generators */
const char* SYMsymmetryGetName(void)
{
   return "bliss";
}

/** return description of external program used to compute generators */
const char* SYMsymmetryGetDesc(void)
{
   return "computing graph automorphism groups by T. Junttila and P. Kaski (http://www.tcs.hut.fi/Software/bliss/)";
}

/** compute generators of symmetry group */
SCIP_RETCODE SYMcomputeSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   int                   fixedtype,          /**< variable types that must be fixed by symmetries */
   SCIP_Bool             local,              /**< Use local variable bounds? */
   int                   npermvars,          /**< number of variables for permutations */
   SCIP_VAR**            permvars,           /**< variables on which permutations act */
   SYM_MATRIXDATA*       matrixdata,         /**< data for MIP matrix */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms               /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   )
{
   assert( scip != NULL );
   assert( matrixdata != NULL );
   assert( permvars != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( perms != NULL );
   assert( maxgenerators >= 0 );

   /* init */
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;

   int nnodes = 0;
   int nedges = 0;

   /* create bliss graph */
   bliss::Graph* G = new bliss::Graph();

   SCIP_Bool success = FALSE;
   SCIP_CALL( fillGraphByColoredCoefficients(scip, G, local, npermvars, permvars, matrixdata, nnodes, nedges, success) );
   if ( ! success )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Graph construction failed.\n");
      delete G;
      return SCIP_OKAY;
   }

#ifdef SCIP_OUTPUT
   G->write_dot("debug.dot");
#endif

   SCIPdebugMsg(scip, "Symmetry detection graph has %u nodes.\n", G->get_nof_vertices());

   /* compute automorphisms */
   bliss::Stats stats;
   BLISS_Data data;
   data.scip = scip;
   data.npermvars = npermvars;
   data.nperms = 0;
   data.nmaxperms = 100 * npermvars;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &data.perms, data.nmaxperms) );

   /* Prefer splitting partition cells corresponding to variables over those corresponding
    * to inequalities. This is because we are only interested in the action
    * of the automorphism group on the variables, and we need a base for this action */
   G->set_splitting_heuristic(bliss::Graph::shs_f);
   /* disable component recursion as advised by Tommi Junttila from bliss */
   G->set_component_recursion(false);

   /* do not use a node limit, but set generator limit */
   G->set_search_limits(0, (unsigned) maxgenerators);

   /* start search */
   G->find_automorphisms(stats, blisshook, (void*) &data);
#ifdef SCIP_OUTPUT
   (void) stats.print(stdout);
#endif

   /* free graph */
   delete G;

   /* prepare return values */
   *perms = data.perms;
   *nperms = data.nperms;
   *nmaxperms = data.nmaxperms;

   return SCIP_OKAY;
}
