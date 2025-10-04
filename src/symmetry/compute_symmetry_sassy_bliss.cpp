/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   compute_symmetry_sassy_bliss.cpp
 * @brief  interface for symmetry computations to sassy as a preprocessor to bliss
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "compute_symmetry.h"

/* include bliss */
#include <bliss/defs.hh>
#include <bliss/graph.hh>

/* include sassy (as part of dejavu) */
#include "build_dejavu_graph.h"
#include <dejavu/tools/bliss_converter.h>

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
   int                   nnodessdg;          /**< number of non-variable nodes in symmetry detection graph */
   SYM_GROUPTYPE         symgrouptype;       /**< type of symmetry group for which generators are computed */
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
   int nsymvars;
   switch ( data->symtype )
   {
   case SYM_SYMTYPE_PERM:
      nsymvars = data->npermvars;
      break;
   default:
      assert( data->symtype == SYM_SYMTYPE_SIGNPERM );
      nsymvars = 2 * data->npermvars;
   }

   int permlen;
   switch( data->symgrouptype )
   {
   case SYM_GROUPTYPE_VAR:
      permlen = nsymvars;
      break;
   case SYM_GROUPTYPE_NODE:
      permlen = nsymvars + data->nnodessdg;
      break;
   default:
      assert( data->symgrouptype == SYM_GROUPTYPE_SDG );
      permlen = n;
   }

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
   if ( data->symgrouptype != SYM_GROUPTYPE_NODE )
   {
      for (int j = 0; j < permlen; ++j)
         p[j] = (int) aut[j];
   }
   else
   {
      int cnt = 0;
      for (int j = nsymvars; j < permlen; ++j, ++cnt)
         p[cnt] = (int) aut[j] - nsymvars;
      for (int j = 0; j < nsymvars; ++j, ++cnt)
         p[cnt] = (int) aut[j] + data->nnodessdg;
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
   return "Computing Graph Automorphisms by T. Junttila and P. Kaski (users.aalto.fi/~tjunttil/bliss)";
}

#define STR(x) #x
#define XSTR(x) STR(x)

/** return name of additional external program used for computing symmetries */
const char* SYMsymmetryGetAddName(void)
{
   return "sassy " XSTR(SASSY_VERSION_MAJOR) "." XSTR(SASSY_VERSION_MINOR);
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
   dejavu::static_graph* G,                  /**< pointer to graph for that automorphisms are computed */
   int                   nnodessdg,          /**< number of non-variable nodes in symmetry detection graph */
   int                   nsymvars,           /**< number of variables encoded in graph */
   int                   maxgenerators,      /**< maximum number of generators to be constructed (=0 if unlimited) */
   int***                perms,              /**< pointer to store generators as (nperms x npermvars) matrix */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations
                                              *   (needed for freeing storage) */
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of size of group */
   SYM_GROUPTYPE         symgrouptype,       /**< type of symmetry group for which generators are computed */
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
   data.nnodessdg = nnodessdg;
   data.symgrouptype = symgrouptype;

   oldtime = SCIPgetSolvingTime(scip);

   /* set up sassy preprocessor */
   dejavu::preprocessor sassy;

   /* lambda function to have access to data and pass it to sassyhook above */
   dejavu_hook sassyglue = [&](int n, const int* p, int nsupp, const int* suppa) {
      sassyhook((void*)&data, n, p, nsupp, suppa);
   };

   /* call sassy to reduce graph */
   sassy.reduce(G, &sassyglue);

   /* create bliss graph */
   bliss::Graph blissgraph(0);

   /* convert sassy to bliss graph */
   convert_dejavu_to_bliss(G, &blissgraph);

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

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   /* lambda function to have access to data and terminate the search if maxgenerators are reached */
   auto term = [&]() {
      return (maxgenerators != 0 && data.nperms >= maxgenerators); /* check the number of generators that we have created so far */
   };

   auto hook = [&](unsigned int n, const unsigned int* aut) {
      sassy.bliss_hook(n, aut);
   };

   /* start search */
   blissgraph.find_automorphisms(stats, hook, term);
#else

   /* Older bliss versions do not allow to terminate with a limit on the number of generators unless patched. */
#ifdef BLISS_PATCH_PRESENT
   /* If patch is present, do not use a node limit, but set generator limit. This approach is not very accurate, since
    * it limits the generators considered in bliss and not the ones that we generate (the ones that work on the variable
    * set). */
   G->set_search_limits(0, (unsigned) maxgenerators);
#endif

   /* start search */
   blissgraph.find_automorphisms(stats, dejavu::preprocessor::bliss_hook, (void*) &sassy);
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

/** compute generators of variable symmetry group */
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

   /* create sassy graph */
   dejavu::static_graph sassygraph;

   SCIP_CALL( SYMbuildDejavuGraph(scip, &sassygraph, graph, &success) );

#ifdef WRITE_GRAPH
   std::string filename = std::string(SCIPgetProbName(scip)) + std::string(".dimacs");
   sassygraph.dump_dimacs(filename);
#endif

   /* compute symmetries */
   SCIP_CALL( computeAutomorphisms(scip, SCIPgetSymgraphSymtype(graph), &sassygraph,
         SCIPgetSymgraphNNodes(graph), SCIPgetSymgraphNVars(graph), maxgenerators,
         perms, nperms, nmaxperms, log10groupsize, SYM_GROUPTYPE_VAR, symcodetime) );

   return SCIP_OKAY;
}

/** compute generators of symmetry group of symmetry detection graph
 *
 *  If the symmetry detection graph (SDG) has k nodes, the first k entries of a generator correspond to the nodes
 *  of the SDG. The remaining entries of the generator correspond to the variables (and possibly their negation).
 */
SCIP_RETCODE SYMcomputeSymmetryGeneratorsNode(
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

   /* create sassy graph */
   dejavu::static_graph sassygraph;

   SCIP_CALL( SYMbuildDejavuGraph(scip, &sassygraph, graph, &success) );

#ifdef WRITE_GRAPH
   std::string filename = std::string(SCIPgetProbName(scip)) + std::string(".dimacs");
   sassygraph.dump_dimacs(filename);
#endif

   /* compute symmetries */
   SCIP_CALL( computeAutomorphisms(scip, SCIPgetSymgraphSymtype(graph), &sassygraph,
         SCIPgetSymgraphNNodes(graph), SCIPgetSymgraphNVars(graph), maxgenerators,
         perms, nperms, nmaxperms, log10groupsize, SYM_GROUPTYPE_NODE, symcodetime) );

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
   int nnodes;
   int nperms;
   int nmaxperms;
   int nnodesfromG1;
   SCIP_Real symcodetime = 0.0;
   SCIP_Real log10groupsize;
   SCIP_Bool success;

   /* create sassy graph */
   dejavu::static_graph sassygraph;

   SCIP_CALL( SYMbuildDejavuGraphCheck(scip, &sassygraph, G1, G2, &nnodes, &nnodesfromG1, &success) );

   if ( ! success )
      return FALSE;

   /* compute symmetries */
   SCIP_CALL_ABORT( computeAutomorphisms(scip, SCIPgetSymgraphSymtype(G1), &sassygraph, -1, -1, 0,
         &perms, &nperms, &nmaxperms, &log10groupsize, SYM_GROUPTYPE_SDG, &symcodetime) );

   /* since G1 and G2 are connected and disjoint, they are isomorphic iff there is a permutation
    * mapping a node from G1 to a node of G2
    */
   success = FALSE;
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
