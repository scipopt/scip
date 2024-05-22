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

/**@file   compute_symmetry_sassy_nauty.c
 * @brief  interface for symmetry computations to sassy as a preprocessor to nauty
 * @author Marc Pfetsch
 * @author Gioni Mexi
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "build_sassy_graph.h"
#include "compute_symmetry.h"

/* the following determines whether nauty or traces is used: */
#define NAUTY

#ifdef NAUTY
#include "nauty/nauty.h"
#include "nauty/nausparse.h"
#else
#include "nauty/traces.h"
#endif

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

#include <sassy/preprocessor.h>
#ifdef NAUTY
#include "sassy/tools/nauty_converter.h"
#else
#include "sassy/tools/traces_converter.h"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic warning "-Wunused-but-set-variable"
#pragma GCC diagnostic warning "-Wsign-compare"
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wshadow"
#endif

#ifdef _MSC_VER
# pragma warning(pop)
#endif

#include "build_sassy_graph.h"

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


/** return whether symmetry can be computed */
SCIP_Bool SYMcanComputeSymmetry(void)
{
   return TRUE;
}


/** return name of external program used to compute generators */
const char* SYMsymmetryGetName(void)
{
#ifdef NAUTY
   return "Nauty " NAUTYVERSION;
#else
   return "Traces " NAUTYVERSION;
#endif
}

/** return description of external program used to compute generators */
const char* SYMsymmetryGetDesc(void)
{
#ifdef NAUTY
   return "Computing Graph Automorphism Groups by Brendan D. McKay (users.cecs.anu.edu.au/~bdm/nauty)";
#else
   return "Computing Graph Automorphism Groups by Adolfo Piperno (pallini.di.uniroma1.it)";
#endif
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

   /* first, convert the graph */
   sparsegraph sg;
   DYNALLSTAT(int, lab, lab_sz);
   DYNALLSTAT(int, ptn, ptn_sz);

#ifdef NAUTY
   convert_sassy_to_nauty(G, &sg, &lab, &lab_sz, &ptn, &ptn_sz);
   statsblk stats;
   DYNALLSTAT(int, orbits, orbits_sz);
   DYNALLOC1(int, orbits, orbits_sz, sg.nv, "malloc");
   DEFAULTOPTIONS_SPARSEGRAPH(options);
   /* init callback functions for nauty (accumulate the group generators found by nauty) */
   options.writeautoms = FALSE;
   options.userautomproc = sassy::preprocessor::nauty_hook;
   options.defaultptn = FALSE; /* use color classes */
   *log10groupsize = 0.0;
   if(sg.nv > 0) {
      sparsenauty(&sg, lab, ptn, orbits, &options, &stats, NULL);
      *log10groupsize = (SCIP_Real) stats.grpsize2;
   }
#else
   convert_sassy_to_traces(&sassygraph, &sg, &lab, &lab_sz, &ptn, &ptn_sz);
   TracesStats stats;
   DYNALLSTAT(int, orbits, orbits_sz);
   DYNALLOC1(int, orbits, orbits_sz, sg.nv, "malloc");
   DEFAULTOPTIONS_TRACES(options);
   /* init callback functions for traces (accumulate the group generators found by traces) */
   options.writeautoms = FALSE;
   options.userautomproc = sassy::preprocessor::traces_hook;
   options.defaultptn = FALSE; /* use color classes */
   if(sg.nv > 0) {
      Traces(&sg, lab, ptn, orbits, &options, &stats, NULL);
   }
#endif

   /* clean up */
   DYNFREE(lab, lab_sz);
   DYNFREE(ptn, ptn_sz);
   SG_FREE(sg);

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

   return SCIP_OKAY;
}

/** compute generators of symmetry group */
SCIP_RETCODE SYMcomputeSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   SYM_GRAPH*            symgraph,           /**< symmetry detection graph */
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
   assert( symgraph != NULL );
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
   sassy::static_graph sassygraph;

   SCIP_CALL( SYMbuildSassyGraph(scip, &sassygraph, symgraph, &success) );

   /* compute symmetries */
   SCIP_CALL( computeAutomorphisms(scip, SCIPgetSymgraphSymtype(symgraph), &sassygraph, SCIPgetSymgraphNVars(symgraph),
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
   int nnodes;
   int nperms;
   int nmaxperms;
   int nnodesfromG1;
   SCIP_Real symcodetime = 0.0;
   SCIP_Real log10groupsize;
   SCIP_Bool success;

   /* create sassy graph */
   sassy::static_graph sassygraph;

   SCIP_CALL( SYMbuildSassyGraphCheck(scip, &sassygraph, G1, G2, &nnodes, &nnodesfromG1, &success) );

   if ( ! success )
      return FALSE;

   /* compute symmetries */
   SCIP_CALL_ABORT( computeAutomorphisms(scip, SCIPgetSymgraphSymtype(G1), &sassygraph, nnodes, 0,
         &perms, &nperms, &nmaxperms, &log10groupsize, FALSE, &symcodetime) );

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
