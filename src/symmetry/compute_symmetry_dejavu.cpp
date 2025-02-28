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

/**@file   compute_symmetry_dejavu.cpp
 * @brief  interface for symmetry computations to dejavu
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "compute_symmetry.h"

/* possibly turn on the dejavu debug mode */
#ifdef SCIP_DEBUG
#define DEJDEBUG
#endif

#include <string>                 /* for dejvu */
#include "build_dejavu_graph.h"   /* also includes dejavu */

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

/** callback function for dejavu */  /*lint -e{715}*/
static
void dejavuhook(
   void*                 user_param,         /**< parameter supplied at call to dejavu */
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
   case SYM_GROUPTYPE_SDG:
      permlen = nsymvars + data->nnodessdg;
      break;
   default:
      assert( data->symgrouptype == SYM_GROUPTYPE_FULL );
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
   if ( data->symgrouptype != SYM_GROUPTYPE_SDG )
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

#define STR(x) #x
#define XSTR(x) STR(x)

/** return name of external program used to compute generators */
const char* SYMsymmetryGetName(void)
{
   return "dejavu " XSTR(DEJAVU_VERSION_MAJOR) "." XSTR(DEJAVU_VERSION_MINOR);
}

/** return description of external program used to compute generators */
const char* SYMsymmetryGetDesc(void)
{
   return "Monte Carlo symmetry computation library by M. Anders (github.com/markusa4/dejavu)";
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

   /* set up dejavu */
   dejavu::solver dejavu;

   dejavu.set_prefer_dfs(true);

#ifndef SCIP_DEBUG
   dejavu.set_print(false);
#endif

   /* lambda function to have access to data and pass it to dejavuhook above */
   dejavu_hook dejavuhookglue = [&](int n, const int* p, int nsupp, const int* suppa) {
      dejavuhook((void*)&data, n, p, nsupp, suppa);
   };

   /* call dejavu */
   dejavu.automorphisms(G, &dejavuhookglue);
   *symcodetime = SCIPgetSolvingTime(scip) - oldtime;

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
   dejavu::big_number size = dejavu.get_automorphism_group_size();
   *log10groupsize = (SCIP_Real) size.exponent + log10(size.mantissa);

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

   /* create dejavu graph */
   dejavu::static_graph dejavugraph;

   SCIP_CALL( SYMbuildDejavuGraph(scip, &dejavugraph, graph, &success) );

#ifdef WRITE_GRAPH
   std::string filename = std::string(SCIPgetProbName(scip)) + std::string(".dimacs");
   dejavugraph.dump_dimacs(filename);
#endif

   /* compute symmetries */
   SCIP_CALL( computeAutomorphisms(scip, SCIPgetSymgraphSymtype(graph), &dejavugraph,
         SCIPgetSymgraphNNodes(graph), SCIPgetSymgraphNVars(graph), maxgenerators,
         perms, nperms, nmaxperms, log10groupsize, SYM_GROUPTYPE_VAR, symcodetime) );

   return SCIP_OKAY;
}

/** compute generators of symmetry group of symmetry detection graph
 *
 *  If the symmetry detection graph (SDG) has k nodes, the first k entries of a generator correspond to the nodes
 *  of the SDG. The remaining entries of the generator correspond to the variables (and possibly their negation).
 */
SCIP_RETCODE SYMcomputeSymmetryGeneratorsSDG(
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

   /* create dejavu graph */
   dejavu::static_graph dejavugraph;

   SCIP_CALL( SYMbuildDejavuGraph(scip, &dejavugraph, graph, &success) );

#ifdef WRITE_GRAPH
   std::string filename = std::string(SCIPgetProbName(scip)) + std::string(".dimacs");
   dejavugraph.dump_dimacs(filename);
#endif

   /* compute symmetries */
   SCIP_CALL( computeAutomorphisms(scip, SCIPgetSymgraphSymtype(graph), &dejavugraph,
         SCIPgetSymgraphNNodes(graph), SCIPgetSymgraphNVars(graph), maxgenerators,
         perms, nperms, nmaxperms, log10groupsize, SYM_GROUPTYPE_SDG, symcodetime) );

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

   /* create dejavu graph */
   dejavu::static_graph dejavugraph;

   SCIP_CALL( SYMbuildDejavuGraphCheck(scip, &dejavugraph, G1, G2, &nnodes, &nnodesfromG1, &success) );

   if ( ! success )
      return FALSE;

   /* compute symmetries */
   SCIP_CALL_ABORT( computeAutomorphisms(scip, SCIPgetSymgraphSymtype(G1), &dejavugraph, -1, -1, 0,
         &perms, &nperms, &nmaxperms, &log10groupsize, SYM_GROUPTYPE_FULL, &symcodetime) );

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
