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

/**@file   sym_linorbred.c
 * @brief  symmetry handler for lexicographic reduction and orbital reduction
 * @author Christopher Hojny
 * @author Jasper van Doornmalen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_setppc.h"
#include "scip/cons_orbitope.h"
#include "scip/event_shadowtree.h"
#include "scip/pub_cons.h"
#include "scip/pub_implics.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sym.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_event.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sym.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include "scip/sym_lexorbred.h"
#include "scip/symmetry.h"
#include "scip/symmetry_orbital.h"
#include "scip/symmetry_lexred.h"
#include "scip/type_implics.h"
#include "scip/type_set.h"

/* symmetry handler properties */
#define SYM_NAME            "sym_lexorbred"
#define SYM_DESC            "symmetry handler for lexicographic reduction and orbital reduction"
#define SYM_PRIORITY            -10000       /**< priority of try-add function*/
#define SYM_PROPPRIORITY       -100000       /**< priority of propagation method */
#define SYM_PROPTIMING SCIP_PROPTIMING_BEFORELP /**< timing of propagator method */
#define SYM_PROPFREQ                 1       /**< frequence of propagator method */
#define SYM_DELAYPROP            FALSE       /**< Should propagation method be delayed, if other propagators found reductions? */

/* default value of parameters */
#define DEFAULT_USELEXRED         TRUE       /**< Shall lexicographic reduction be used? */
#define DEFAULT_USEORBRED         TRUE       /**< Shall orbital reduction be used? */
#define DEFAULT_COMPUTENEWPERMS   TRUE       /**< Shall additional permutations of symmetry component be computed? */
#define DEFAULT_MAXNNEWPERMS       100       /**< maximum number of additional permutations of symmetry component
                                              *   that is computed (used to have better approximation of stabilizer);
                                              *   -1: unbounded */
#define DEFAULT_USEPPPUPGRADE     TRUE       /**< Shall static constraints be used if at least 50% of perms can be handled by pp-orbisacks? */

/*
 * Data structures
 */

/** symmetry component data */
struct SCIP_SymCompData
{
   SCIP_CONS**           conss;              /**< static symmetry handling constraints */
   int                   nconss;             /**< number of static symmetry handling constraints */
   int                   maxnconss;          /**< maximum number of constraints conss can hold */
   SCIP_LEXREDDATA*      lexreddata;         /**< container for lexicographic reduction propagation; */
   SCIP_Bool             active;             /**< whether lexicographic reduction is active on this component */
};

/** symmetry handler data */
struct SCIP_SymhdlrData
{
   SCIP_EVENTHDLR*       shadowtreeeventhdlr; /**< shadow tree event handler */
   SCIP_ORBITALREDDATA*  orbitalreddata;     /**< container for orbital reduction data */
   SCIP_Bool             uselexred;          /**< Shall lexicographic reduction be used? */
   SCIP_Bool             useorbred;          /**< Shall orbital reduction be used? */
   SCIP_Bool             computenewperms;    /**< Shall additional permutations of symmetry component be computed? */
   int                   maxnnewperms;       /**< maximum number of additional permutations of symmetry component
                                              *   that is computed (used to have better approximation of stabilizer);
                                              *   -1: unbounded */
   SCIP_Bool             useppupgrade;       /**< Shall static constraints be used if at least 50% of perms can be
                                              *   handled by pp-orbisacks? */

};

/** compare function for sorting an array by the addresses of its members  */
static
SCIP_DECL_SORTPTRCOMP(sortByPointerValue)
{
   /* @todo move to misc.c? */
   if ( elem1 < elem2 )
      return -1;
   else if ( elem1 > elem2 )
      return +1;
   return 0;
}

/*
 * Local methods
 */

/** checks whether two arrays that are sorted with the same comparator have a common element */
static
SCIP_Bool checkSortedArraysHaveOverlappingEntry(
   void**                arr1,               /**< first array */
   int                   narr1,              /**< number of elements in first array */
   void**                arr2,               /**< second array */
   int                   narr2,              /**< number of elements in second array */
   SCIP_DECL_SORTPTRCOMP((*compfunc))        /**< comparator function that was used to sort arr1 and arr2; must define a total ordering */
   )
{
   /* @todo move to misc.c? */
   int it1;
   int it2;
   int cmp;

   assert(arr1 != NULL || narr1 == 0);
   assert(narr1 >= 0);
   assert(arr2 != NULL || narr2 == 0);
   assert(narr2 >= 0);
   assert(compfunc != NULL);

   /* there is no overlap if one of the two arrays is empty */
   if( narr1 <= 0 )
      return FALSE;
   if( narr2 <= 0 )
      return FALSE;

   it1 = 0;
   it2 = 0;

   while( TRUE )  /*lint !e716*/
   {
      cmp = compfunc(arr1[it1], arr2[it2]);
      if( cmp < 0 )
      {
         /* comparison function determines arr1[it1] < arr2[it2]
          * increase iterator for arr1
          */
         if( ++it1 >= narr1 )
            break;
         continue;
      }
      else if( cmp > 0 )
      {
         /* comparison function determines arr1[it1] > arr2[it2]
          * increase iterator for arr2
          */
         if( ++it2 >= narr2 )
            break;
         continue;
      }
      else
      {
         /* the entries arr1[it1] and arr2[it2] are the same with respect to the comparison function */
         assert(cmp == 0);
         return TRUE;
      }
   }

   /* no overlap detected */
   assert(it1 >= narr1 || it2 >= narr2);
   return FALSE;
}

/** returns whether a permutation is an involution */
static
SCIP_Bool isInvolution(
   int*                  perm,               /**< permutation */
   int                   lenperm,            /**< length of permutation */
   SCIP_Bool*            istransposition     /**< pointer to store whether permutation is a transposition */
   )
{
   int lensupport = 0;
   int i;

   assert(perm != NULL);
   assert(lenperm > 0);
   assert(istransposition != NULL);
   *istransposition = FALSE;

   for( i = 0; i < lenperm; ++i )
   {
      if( perm[perm[i]] != i )
         return FALSE;
      if( perm[i] != i )
         ++lensupport;
   }

   if( lensupport == 2 )
      *istransposition = TRUE;
   return TRUE;
}

/** tries to generate involutions based on permutations in symmetry component
 *
 *  An involution is a permutation whose 2-fold application is the identity. Involutions
 *  are of particular interest, because their support is (in practice) often very small.
 *  Propagations based on involutions thus can be executed rather quickly. Moreover,
 *  when computing stabilizer subgroups by a filtering mechanism, it is more likely that
 *  a (sparse) involution is not filtered in contrast to dense non-involutions.
 *
 *  To create new involutions, we are given a list of involutions that have been found
 *  during symmetry detection. We then iterate through all pairs (p,q) of involutions in
 *  this list and generate the new involutions p * q (if p and q commute) and p * q * p
 *  (if p and q do not commute).
 */
static
SCIP_RETCODE tryGenerateInvolutions(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< (signed) permutations matrix */
   int                   nperms,             /**< number of permutations */
   int                   npermvars,          /**< number of variables the permutations act on */
   int                   maxnnewinvolus,     /**< maximum number of involutions to be computed */
   int***                newperms,           /**< pointer to store new (signed) permutations */
   int*                  nnewperms,          /**< pointer to store number of new permutations */
   int*                  lennewperms         /**< pointer to store length of newperms */
   )
{
   int* tmpperm;
   int* perm1;
   int* perm2;
   int permlen;
   int p;
   int q;
   int i;
   SCIP_Bool commute;
   SCIP_Bool istransposition;
   SCIP_Bool dynamicmemsize;

   /* check whether we shall run */
   if( maxnnewinvolus <= 0 )
      return SCIP_OKAY;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms > 0);
   assert(npermvars > 0);
   assert(newperms != NULL);
   assert(nnewperms != NULL);
   assert(lennewperms != NULL);

   dynamicmemsize = maxnnewinvolus == -1;
   *lennewperms = dynamicmemsize ? 100 : maxnnewinvolus;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, newperms, *lennewperms) );
   *nnewperms = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpperm, npermvars) );

   permlen = symtype == (int)SYM_SYMTYPE_PERM ? npermvars : 2 * npermvars;

   /* try to generate new involutions by combining two involutions p and q
    *
    * If p and q commute, create involution p*q. Otherwise, create involutions p*q*p and q*p*q.
    */
   for( p = 0; p < nperms; ++p )
   {
      perm1 = perms[p];
      if( !isInvolution(perm1, npermvars, &istransposition) )
         continue;

      /* it seems promising to only combine involutions with at least two cycles each */
      if ( istransposition )
         continue;

      for( q = p + 1; q < nperms; ++q )
      {
         perm2 = perms[q];
         if( !isInvolution(perm2, npermvars, &istransposition) )
            continue;
         if( istransposition )
            continue;

         /* check whether perm1 and perm2 commute */
         commute = TRUE;
         for( i = 0; i < npermvars && commute; ++i )
         {
            if( perm1[perm2[i]] != perm2[perm1[i]] )
               commute = FALSE;
         }
         /* only consider involutions that have non-disjoint support */
         if( commute )
         {
            /* permutations commute, store perm1 * perm2 if we do not know it yet */
            for( i = 0; i < npermvars; ++i )
               tmpperm[i] = perm1[perm2[i]];

            if( isPermKnown(tmpperm, npermvars, perms, nperms) )
               continue;
            if ( isPermKnown(tmpperm, npermvars, *newperms, *nnewperms) )
               continue;
            if( dynamicmemsize && *nnewperms == *lennewperms )
            {
               SCIP_CALL( SCIPensureBlockMemoryArray(scip, newperms, lennewperms, *lennewperms + 1) );
            }
            assert( *nnewperms < *lennewperms );

            /* recompute permutation, because we possibly also need the entries for negated variables */
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*newperms)[*nnewperms], permlen) );
            for( i = 0; i < permlen; ++i )
               (*newperms)[*nnewperms][i] = perm1[perm2[i]];
            ++(*nnewperms);
         }
         else
         {
            /* permutations do not commute, compute perm1 * perm2 * perm1 */
            for( i = 0; i < npermvars; ++i )
               tmpperm[i] = perm1[perm2[perm1[i]]];

            /* do not store the permutation if it is already known
             * (also for signed permutations, it is sufficient to iterate over the first npermvars
             *  entries, because the permutation on the negated variables can be derived from these entries)
             **/
            if( isPermKnown(tmpperm, npermvars, perms, nperms) )
               continue;
            if( isPermKnown(tmpperm, npermvars, *newperms, *nnewperms) )
               continue;

            /* we do not know the permutation yet, store it */
            if( dynamicmemsize && *nnewperms == *lennewperms )
            {
               SCIP_CALL( SCIPensureBlockMemoryArray(scip, newperms, lennewperms, *lennewperms + 1) );
            }
            assert(*nnewperms < *lennewperms);

            /* recompute permutation, because we possibly also need the entries for negated variables */
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*newperms)[*nnewperms], permlen) );
            for( i = 0; i < permlen; ++i )
               (*newperms)[*nnewperms][i] = perm1[perm2[perm1[i]]];
            ++(*nnewperms);

            if( !dynamicmemsize && *nnewperms == *lennewperms )
               break;

            /* compute perm2 * perm1 * perm2 */
            for( i = 0; i < npermvars; ++i )
               tmpperm[i] = perm2[perm1[perm2[i]]];

            /* do not store the permutation if it is already known */
            if( isPermKnown(tmpperm, npermvars, perms, nperms) )
               continue;
            if ( isPermKnown(tmpperm, npermvars, *newperms, *nnewperms) )
               continue;

            /* recompute permutation, because we possibly also need the entries for negated variables */
            if( dynamicmemsize && *nnewperms == *lennewperms )
            {
               SCIP_CALL( SCIPensureBlockMemoryArray(scip, newperms, lennewperms, *lennewperms + 1) );
            }
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*newperms)[*nnewperms], permlen) );
            for( i = 0; i < permlen; ++i )
               (*newperms)[*nnewperms][i] = perm2[perm1[perm2[i]]];
            ++(*nnewperms);
         }

         if( !dynamicmemsize && *nnewperms >= *lennewperms )
            break;
      }
      if( !dynamicmemsize && *nnewperms >= *lennewperms )
         break;
   }

   SCIPfreeBufferArray(scip, &tmpperm);

   if( *nnewperms == 0 )
   {
      SCIPfreeBlockMemoryArray(scip, newperms, *lennewperms);
   }

   return SCIP_OKAY;
}

/** add lexicographic reduction for given permutations */
static
SCIP_RETCODE addLexRed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMCOMPDATA*     symcompdata,        /**< data of symmetry component */
   SCIP_EVENTHDLR*       shadowtreeeventhdlr,/**< shadow tree event handler */
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< permutations */
   int                   nperms,             /**< number of permutations */
   int**                 moreperms,          /**< additional permutations */
   int                   nmoreperms,         /**< number of additional permutations */
   SCIP_VAR**            permvars,           /**< variables permutations act on */
   int                   npermvars,          /**< numnber of variables */
   SCIP_Real*            permvardomaincenter,/**< array of centers of variable domains */
   SCIP_Bool*            success             /**< pointer to store whether lexred could be added */
   )
{
   SCIP_Bool locsuccess;
   int p;

   assert(scip != NULL);
   assert(symcompdata != NULL);
   assert(shadowtreeeventhdlr != NULL);
   assert(perms != NULL);
   assert(nperms > 0);
   assert(moreperms != NULL || nmoreperms == 0);
   assert(permvars != NULL);
   assert(npermvars > 0);
   assert(permvardomaincenter != NULL || symtype != SYM_SYMTYPE_SIGNPERM);
   assert(success != NULL);

   *success = FALSE;

   SCIP_CALL( SCIPincludeLexicographicReduction(scip, &symcompdata->lexreddata, shadowtreeeventhdlr) );
   assert(symcompdata->lexreddata != NULL);

   /* add all permutations to lexicograph reduction */
   for( p = 0; p < nperms; ++p )
   {
      SCIP_CALL( SCIPlexicographicReductionAddPermutation(scip, symcompdata->lexreddata,
            permvars, npermvars, perms[p], (SYM_SYMTYPE) symtype, permvardomaincenter, TRUE, &locsuccess) );
      *success |= locsuccess;
   }
   for( p = 0; p < nmoreperms; ++p )
   {
      SCIP_CALL( SCIPlexicographicReductionAddPermutation(scip, symcompdata->lexreddata,
            permvars, npermvars, moreperms[p], (SYM_SYMTYPE) symtype, permvardomaincenter, TRUE, &locsuccess) );
      *success |= locsuccess;
   }
   symcompdata->active = *success;

   return SCIP_OKAY;
}

/** adds orbital reduction for given permutations */
static
SCIP_RETCODE addOrbRed(
   SCIP*                 scip,               /**< SCIP datat structure */
   SCIP_ORBITALREDDATA*  orbitalreddata,     /**< data for orbital reduction */
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< permutations */
   int                   nperms,             /**< number of permutations */
   int**                 moreperms,          /**< additional permutations */
   int                   nmoreperms,         /**< number of additional permutations */
   SCIP_VAR**            permvars,           /**< variables the permutations act on */
   int                   npermvars,          /**< number of variables */
   SCIP_Bool*            success             /**< pointer to store whether orbital reduction could be added */
   )
{
   SCIP_Bool freeproperperms = FALSE;
   int** properperms;
   int nproperperms;
   int p;

   assert(scip != NULL);
   assert(orbitalreddata != NULL);
   assert(perms != NULL);
   assert(nperms > 0);
   assert(moreperms != NULL || nmoreperms == 0);
   assert(permvars != NULL);
   assert(npermvars > 0);
   assert(success != NULL);

   /* store proper permutations in single data structure */
   if( symtype != SYM_SYMTYPE_PERM || nmoreperms > 0 )
   {
      int i;

      nproperperms = 0;
      freeproperperms = TRUE;

      SCIP_CALL( SCIPallocBufferArray(scip, &properperms, nperms + nmoreperms) );
      for( p = 0; p < nperms; ++p )
      {
         if( isProperPerm(symtype, perms[p], npermvars) )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &properperms[nproperperms], npermvars) );
            for( i = 0; i < npermvars; ++i )
               properperms[nproperperms][i] = perms[p][i];
            ++nproperperms;
         }
      }
      for( p = 0; p < nmoreperms; ++p )
      {
         if( isProperPerm(symtype, moreperms[p], npermvars) )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &properperms[nproperperms], npermvars) );
            for( i = 0; i < npermvars; ++i )
               properperms[nproperperms][i] = moreperms[p][i];
            ++nproperperms;
         }
      }
   }
   else
   {
      properperms = perms;
      nproperperms = nperms;
   }

   SCIP_CALL( SCIPorbitalReductionAddComponent(scip, orbitalreddata, permvars, npermvars,
         properperms, nproperperms, success) );

   if( freeproperperms )
   {
      for( p = nproperperms - 1; p >= 0; --p )
      {
         SCIPfreeBufferArray(scip, &properperms[p]);
      }
      SCIPfreeBufferArray(scip, &properperms);
   }

   return SCIP_OKAY;
}

/** tries to apply pp-orbitope upgrade if at least 50% of the permutations to pp-orbisacks */
static
SCIP_RETCODE tryPackingPartitioningOrbisackUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 perms,              /**< permutations */
   int                   nperms,             /**< number of permutations */
   SCIP_Bool             hassignedperm,      /**< whether perms contains a signed permutation */
   SCIP_VAR**            permvars,           /**< variables the permutations act on */
   int                   npermvars,          /**< number of variables */
   SCIP_HASHMAP*         permvarmap,         /**< map of variables to indices in permvars array */
   SCIP_CONS***          addedconss,         /**< pointer to store created constraints */
   int*                  naddedconss,        /**< pointer to store number of added constraints */
   int*                  maxnaddedconss,     /**< pointer to store maximum number addedconss can hold */
   int                   id,                 /**< identifier of symmetry component */
   SCIP_Bool*            success             /**< pointer to store if the packing partitioning upgrade succeeded */
   )
{
   int c;
   int i;
   int j;
   int p;
   int* perm;
   SCIP_CONSHDLR* setppcconshdlr;
   SCIP_CONS** setppcconss;
   SCIP_CONS* cons;
   SCIP_CONS** setppconsssort;
   int nsetppconss;
   int nsetppcvars;
   SCIP_VAR** setppcvars;
   int nsetppcconss;
   int** pporbisackperms;
   int npporbisackperms;
   SCIP_VAR* var;
   int varid;
   SCIP_CONS*** permvarssetppcconss;
   int* npermvarssetppcconss;
   int* maxnpermvarssetppcconss;
   int maxntwocycles;
   int ntwocycles;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms > 0);
   assert(success != NULL);

   /* we did not upgrade yet */
   *success = FALSE;
   *naddedconss = 0;

   /* currently, we cannot handle signed permutations */
   if( hassignedperm )
      return SCIP_OKAY;

   setppcconshdlr = SCIPfindConshdlr(scip, "setppc");
   if( setppcconshdlr == NULL )
      return SCIP_OKAY;

   nsetppcconss = SCIPconshdlrGetNConss(setppcconshdlr);
   if( nsetppcconss == 0 )
      return SCIP_OKAY;

   setppcconss = SCIPconshdlrGetConss(setppcconshdlr);
   assert(setppcconss != NULL);

   /* collect non-covering constraints and sort by pointer for easy intersection finding */
   SCIP_CALL( SCIPallocBufferArray(scip, &setppconsssort, nsetppcconss) );
   nsetppconss = 0;
   for( c = 0; c < nsetppcconss; ++c )
   {
      cons = setppcconss[c];

      /* only packing or partitioning constraints, no covering types */
      if( SCIPgetTypeSetppc(scip, cons) == SCIP_SETPPCTYPE_COVERING )
         continue;

      setppconsssort[nsetppconss++] = cons;
   }
   SCIPsortPtr((void**) setppconsssort, sortByPointerValue, nsetppconss);

   /* For each permvar, introduce an array of setppc constraints (initially NULL) for each variable,
    * and populate it with the setppc constraints that it contains. This array follows the ordering by cons ptr address.
    */
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &permvarssetppcconss, npermvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &npermvarssetppcconss, npermvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &maxnpermvarssetppcconss, npermvars) );
   for( c = 0; c < nsetppconss; ++c )
   {
      cons = setppconsssort[c];
      assert(cons != NULL);

      setppcvars = SCIPgetVarsSetppc(scip, cons);
      nsetppcvars = SCIPgetNVarsSetppc(scip, cons);

      for( i = 0; i < nsetppcvars; ++i )
      {
         var = setppcvars[i];
         assert(var != NULL);
         varid = SCIPhashmapGetImageInt(permvarmap, (void*) var);
         assert(varid == INT_MAX || varid < npermvars);
         assert(varid >= 0);
         if( varid < npermvars )
         {
            SCIP_CALL( SCIPensureBlockMemoryArray(scip, &(permvarssetppcconss[varid]),
                  &maxnpermvarssetppcconss[varid], npermvarssetppcconss[varid] + 1) );
            assert(npermvarssetppcconss[varid] < maxnpermvarssetppcconss[varid]);
            permvarssetppcconss[varid][npermvarssetppcconss[varid]++] = cons;
         }
      }
   }

   /* for all permutations, test involutions on binary variables and test if they are captured by setppc conss */
   SCIP_CALL( SCIPallocBufferArray(scip, &pporbisackperms, nperms) );
   maxntwocycles = 0;
   npporbisackperms = 0;
   for( p = 0; p < nperms; ++p )
   {
      perm = perms[p];
      ntwocycles = 0;

      /* check if the binary orbits are involutions */
      for( i = 0; i < npermvars; ++i )
      {
         j = perm[i];

         /* ignore fixed points in permutation */
         if( i == j )
            continue;
         /* only check for situations where i and j are binary variables */
         assert(SCIPgetSymInferredVarType(permvars[i]) == SCIPgetSymInferredVarType(permvars[j]));
         if( SCIPgetSymInferredVarType(permvars[i]) != SCIP_VARTYPE_BINARY )
            continue;
         /* the permutation must be an involution on binary variables */
         if( perm[j] != i )
            goto NEXTPERMITER;
         /* i and j are a two-cycle, so we find this once for i and once for j. Only handle this once for i < j. */
         if( i > j )
            continue;
         /* disqualify permutation if i and j are not in a common set packing constraint */
         if( !checkSortedArraysHaveOverlappingEntry((void**) permvarssetppcconss[i], npermvarssetppcconss[i],
             (void**) permvarssetppcconss[j], npermvarssetppcconss[j], sortByPointerValue) )
            goto NEXTPERMITER;
         ++ntwocycles;
      }

      /* The permutation qualifies if all binary variables are either a reflection or in a 2-cycle. There must be at
       * least one binary 2-cycle, because otherwise the permutation is the identity, or it permutes
       * nonbinary variables.
       */
      if( ntwocycles > 0 )
      {
         pporbisackperms[npporbisackperms++] = perm;
         if( ntwocycles > maxntwocycles )
            maxntwocycles = ntwocycles;
      }

   NEXTPERMITER:
      ;
   }

   /* if at least 50% of such permutations are packing-partitioning type, apply packing upgrade */
   if( npporbisackperms * 2 >= nperms )
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_VAR** ppvarsblock;
      SCIP_VAR*** ppvarsmatrix;
      SCIP_VAR** row;
      int nrows;

      assert(npporbisackperms > 0);
      assert(maxntwocycles > 0);

      /* instead of allocating and re-allocating multiple times, recycle the ppvars array */
      SCIP_CALL( SCIPallocBufferArray(scip, &ppvarsblock, 2 * maxntwocycles) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ppvarsmatrix, maxntwocycles) );
      for( i = 0; i < maxntwocycles; ++i )
         ppvarsmatrix[i] = &(ppvarsblock[2 * i]);

      /* for each of these perms, create the packing orbitope matrix and add constraint*/
      for( p = 0; p < npporbisackperms; ++p )
      {
         perm = pporbisackperms[p];

         /* populate ppvarsmatrix */
         nrows = 0;
         for( i = 0; i < npermvars; ++i )
         {
            j = perm[i];

            /* ignore fixed points in permutation, and only consider rows with i < j */
            if( i >= j )
               continue;
            /* only check for situations where i and j are binary variables */
            assert(SCIPgetSymInferredVarType(permvars[i]) == SCIPgetSymInferredVarType(permvars[j]));
            if( SCIPgetSymInferredVarType(permvars[i]) != SCIP_VARTYPE_BINARY )
               continue;
            assert(perm[j] == i);
            assert(checkSortedArraysHaveOverlappingEntry((void**) permvarssetppcconss[i], npermvarssetppcconss[i],
               (void**) permvarssetppcconss[j], npermvarssetppcconss[j], sortByPointerValue));

            assert(nrows < maxntwocycles);
            row = ppvarsmatrix[nrows++];
            row[0] = permvars[i];
            row[1] = permvars[j];
            assert(row[0] != row[1]);
         }
         assert(nrows > 0);

         /* create constraint, use same parameterization as in orbitope packing partitioning checker */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "orbitope_pp_upgrade_lexred_%d_%d", id, p);
         SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, name, ppvarsmatrix, SCIP_ORBITOPETYPE_PACKING, nrows, 2,
               FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPensureBlockMemoryArray(scip, addedconss, maxnaddedconss, *naddedconss + 1) );
         (*addedconss)[(*naddedconss)++] = cons;
         SCIP_CALL( SCIPaddCons(scip, cons) );
      }

      SCIPfreeBufferArray(scip, &ppvarsmatrix);
      SCIPfreeBufferArray(scip, &ppvarsblock);

      *success = TRUE;
   }

   /* free pp orbisack array */
   SCIPfreeBufferArray(scip, &pporbisackperms);

   /* clean the non-clean arrays */
   for( varid = 0; varid < npermvars; ++varid )
   {
      assert((permvarssetppcconss[varid] == NULL) == (maxnpermvarssetppcconss[varid] == 0));
      assert(npermvarssetppcconss[varid] >= 0);
      assert(maxnpermvarssetppcconss[varid] >= 0);
      assert(npermvarssetppcconss[varid] <= maxnpermvarssetppcconss[varid]);
      if( npermvarssetppcconss[varid] == 0 )
         continue;
      SCIPfreeBlockMemoryArray(scip, &permvarssetppcconss[varid], maxnpermvarssetppcconss[varid]);
      permvarssetppcconss[varid] = NULL;
      npermvarssetppcconss[varid] = 0;
      maxnpermvarssetppcconss[varid] = 0;
   }
   SCIPfreeCleanBufferArray(scip, &maxnpermvarssetppcconss);
   SCIPfreeCleanBufferArray(scip, &npermvarssetppcconss);
   SCIPfreeCleanBufferArray(scip, &permvarssetppcconss);
   SCIPfreeBufferArray(scip, &setppconsssort);

   return SCIP_OKAY;
}

/** returns whether a list of permutations contains a signed permutation */
static
SCIP_Bool hasSignedPerm(
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< permutations */
   int                   nperms,             /**< number of permutations */
   int                   npermvars           /**< number of variables the permutations act on */
   )
{
   int p;
   int i;

   assert(perms != NULL);
   assert(nperms > 0);
   assert(npermvars > 0);

   if( symtype == SYM_SYMTYPE_PERM )
      return FALSE;

   for( p = 0; p < nperms; ++p )
   {
      for( i = 0; i < npermvars; ++i )
      {
         if( perms[p][i] >= npermvars )
            return TRUE;
      }
   }

   return FALSE;
}

/*
 * Callback methods of symmetry handler
 */

/** addition method for symmetry method handler plugins (tries to add symmetry handling method for given symmetries) */
static
SCIP_DECL_SYMHDLRTRYADD(symhdlrTryaddLexOrbRed)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;
   int** newperms = NULL;
   int nnewperms = 0;
   int maxnnewperms = 0;
   int p;

   assert(success != NULL);
   assert(symhdlr != NULL);
   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms >= 0);
   assert(permvars != NULL || npermvars == 0);
   assert(naddedconss != NULL);

   *success = FALSE;

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   if( !symhdlrdata->uselexred && !symhdlrdata->useorbred )
      return SCIP_OKAY;

   /* possibly create new permutations */
   if( symhdlrdata->computenewperms )
   {
      SCIP_CALL( tryGenerateInvolutions(scip, symtype, perms, nperms, npermvars, symhdlrdata->maxnnewperms,
            &newperms, &nnewperms, &maxnnewperms) );
   }
   assert(symhdlrdata->shadowtreeeventhdlr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, symcompdata) );
   (*symcompdata)->lexreddata = NULL;
   (*symcompdata)->active = FALSE;
   (*symcompdata)->conss = NULL;
   (*symcompdata)->nconss = 0;
   (*symcompdata)->maxnconss = 0;

   /* check whether static constraints shall be used in case there are many set packing/partitioning constraints */
   if( symhdlrdata->useppupgrade )
   {
      /* it is sufficient to check for original permutations, because newperms has signed perm iff perms has */
      if( !hasSignedPerm(symtype, perms, nperms, npermvars) )
      {
         int** allperms;
         int nallperms;
         int cnt = 0;
         int i;

         /* store permutations in a single array */
         nallperms = nperms + nnewperms;
         if( nnewperms > 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &allperms, nallperms) );
            for( p = 0; p < nperms; ++p )
            {
               SCIP_CALL( SCIPallocBufferArray(scip, &allperms[p], npermvars) );
               for( i = 0; i < npermvars; ++i )
                  allperms[p][i] = perms[p][i];
            }
            for( p = 0, cnt = nperms; p < nnewperms; ++p, ++cnt )
            {
               SCIP_CALL( SCIPallocBufferArray(scip, &allperms[cnt], npermvars) );
               for( i = 0; i < npermvars; ++i )
                  allperms[cnt][i] = newperms[p][i];
            }
         }
         else
            allperms = perms;

         SCIP_CALL( tryPackingPartitioningOrbisackUpgrade(scip, allperms, nallperms, FALSE,
               permvars, npermvars, permvarmap, &(*symcompdata)->conss, &(*symcompdata)->nconss,
               &(*symcompdata)->maxnconss, id, success) );
         *naddedconss = (*symcompdata)->nconss;

         if( nnewperms > 0 )
         {
            for( p = nallperms - 1; p >= 0; --p )
            {
               SCIPfreeBufferArray(scip, &allperms[p]);
            }
            SCIPfreeBufferArray(scip, &allperms);
         }

         if( *success )
            goto FREEMEMORY;
      }
   }
   if( symhdlrdata->uselexred )
   {
      SCIP_CALL( addLexRed(scip, *symcompdata, symhdlrdata->shadowtreeeventhdlr,
            symtype, perms, nperms, newperms, nnewperms, permvars, npermvars, permvardomcenter, success) );
   }
   if( symhdlrdata->useorbred )
   {
      SCIP_Bool locsuccess = FALSE;

      SCIP_CALL( addOrbRed(scip, symhdlrdata->orbitalreddata, symtype, perms, nperms, newperms, nnewperms,
            permvars, npermvars, &locsuccess) );
      *success |= locsuccess;
   }

 FREEMEMORY:
   /* free new permutations */
   for( p = nnewperms - 1; p >= 0; --p )
   {
      SCIPfreeBlockMemoryArray(scip, &newperms[p], symtype == SYM_SYMTYPE_PERM ? npermvars : 2 * npermvars);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &newperms, maxnnewperms);

   return SCIP_OKAY;
}

/** solving process deinitialization method of symmetry handler (called before branch and bound process data is freed) */
static
SCIP_DECL_SYMHDLREXITSOL(symhdlrExitsolLexOrbRed)
{ /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;
   SCIP_SYMCOMPDATA* symdata;
   int s;

   assert(symcomps != NULL || nsymcomps == 0);
   assert(symhdlr != NULL);

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   for( s = 0; s < nsymcomps; ++s )
   {
      assert(symcomps[s] != NULL);

      symdata = SCIPsymcompGetData(symcomps[s]);
      assert(symdata != NULL);
      assert(symdata->lexreddata != NULL || !symdata->active);

      SCIP_CALL( SCIPlexicographicReductionReset(scip, symdata->lexreddata) );
   }
   assert(symhdlrdata->orbitalreddata != NULL);

   SCIP_CALL( SCIPorbitalReductionReset(scip, symhdlrdata->orbitalreddata) );

   return SCIP_OKAY;
}

/** deinitialization method of symmetry handler (called before transformed problem is freed) */
static
SCIP_DECL_SYMHDLREXIT(symhdlrExitLexOrbRed)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;
   SCIP_SYMCOMPDATA* symdata;
   int s;
   int c;

   assert(symcomps != NULL || nsymcomps == 0);
   assert(symhdlr != NULL);

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   for( s = 0; s < nsymcomps; ++s )
   {
      assert(symcomps[s] != NULL);

      symdata = SCIPsymcompGetData(symcomps[s]);
      assert(symdata != NULL);
      assert(symdata->lexreddata != NULL || !symdata->active);

      if( symdata->active )
      {
         SCIP_CALL( SCIPlexicographicReductionReset(scip, symdata->lexreddata) );
         SCIP_CALL( SCIPlexicographicReductionFree(scip, &symdata->lexreddata) );
      }

      for( c = 0; c < symdata->nconss; ++c )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &symdata->conss[c]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &symdata->conss, symdata->maxnconss);

      SCIPfreeBlockMemory(scip, &symdata);
   }
   assert(symhdlrdata->orbitalreddata != NULL);

   SCIP_CALL( SCIPorbitalReductionReset(scip, symhdlrdata->orbitalreddata) );

   return SCIP_OKAY;
}

/** destructor of symmetry handler to free symmetry handler data (called when SCIP is exiting) */
static
SCIP_DECL_SYMHDLRFREE(symhdlrFreeLexOrbRed)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;

   assert(scip != NULL);
   assert(symhdlr != NULL);

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);
   assert(symhdlrdata->orbitalreddata != NULL);

   SCIP_CALL( SCIPorbitalReductionFree(scip, &symhdlrdata->orbitalreddata) );
   SCIPfreeBlockMemory(scip, &symhdlrdata);

   return SCIP_OKAY;
}

/** domain propagation method of symmetry handler */
static
SCIP_DECL_SYMHDLRPROP(symhdlrPropLexOrbRed)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;
   SCIP_SYMCOMPDATA* symcompdata;
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool didrun = FALSE;
   SCIP_Bool didrunlocal;
   int nredlocal;
   int nreds = 0;
   int s;

   assert(symhdlr != NULL);
   assert(symcomps != NULL || nsymcomps == 0);
   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   /* do not run if we are in the root or not yet solving */
   if ( SCIPgetDepth(scip) <= 0 || SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   /* run orbital reduction */
   SCIP_CALL( SCIPorbitalReductionPropagate(scip, symhdlrdata->orbitalreddata,
         &infeasible, &nreds, &didrun) );
   if ( infeasible )
      return SCIP_OKAY;

   /* run lexicographic reduction */
   for( s = 0; s < nsymcomps; ++s )
   {
      assert(symcomps[s] != NULL);

      symcompdata = SCIPsymcompGetData(symcomps[s]);
      assert(symcompdata != NULL);

      if( !symcompdata->active )
         continue;

      SCIP_CALL( SCIPlexicographicReductionPropagate(scip, symcompdata->lexreddata,
            &infeasible, &nredlocal, &didrunlocal) );
      nreds += nredlocal;
      didrun |= didrunlocal;
      if ( infeasible )
         return SCIP_OKAY;
   }

   if( infeasible )
      *result = SCIP_CUTOFF;
   else if( nreds > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( didrun )
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** presolving method of symmetry handler */
static
SCIP_DECL_SYMHDLRPRESOL(symhdlrPresolLexOrbRed)
{  /*lint --e{715}*/
   SCIP_SYMCOMPDATA* symdata;
   int s;
   int c;

   assert(result != NULL);
   assert(symcomps != NULL || nsymcomps == 0);

   *result = nsymcomps > 0 ? SCIP_DIDNOTFIND : SCIP_DIDNOTRUN;

   for( s = 0; s < nsymcomps; ++s )
   {
      symdata = SCIPsymcompGetData(symcomps[s]);

      for( c = 0; c < symdata->nconss; ++c )
      {
         SCIP_CALL( SCIPpresolCons(scip, symdata->conss[c], nrounds, presoltiming, nnewfixedvars, nnewaggrvars,
               nnewchgvartypes, nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs,
               nnewchgsides, nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes, ndelconss, naddconss,
               nupgdconss, nchgcoefs, nchgsides, result) );

         /* exit if cutoff or unboundedness has been detected */
         if( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
         {
            SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n",
               SCIPconsGetName(symdata->conss[c]));
            return SCIP_OKAY;
         }
      }

   }

   return SCIP_OKAY;
}

/** include symmetry handler for lexicographic reduction and orbital reduction */
SCIP_RETCODE SCIPincludeSymhdlrLexOrbRed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SYMHDLRDATA* symhdlrdata = NULL;
   SCIP_EVENTHDLR* eventhdlr;

   assert(scip != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &symhdlrdata) );

   SCIP_CALL( SCIPincludeSymhdlrBasic(scip, SYM_NAME, SYM_DESC, SYM_PRIORITY, SYM_PROPPRIORITY, 0, -1,
         SYM_PROPFREQ, -1, SYM_DELAYPROP, FALSE, 1.0, 1, SYM_PROPTIMING, SCIP_PRESOLTIMING_FAST,
         symhdlrTryaddLexOrbRed, NULL, symhdlrFreeLexOrbRed, NULL, symhdlrExitLexOrbRed,
         NULL, symhdlrExitsolLexOrbRed, NULL, NULL, NULL, NULL, symhdlrPropLexOrbRed,
         NULL, symhdlrPresolLexOrbRed, symhdlrdata) );

   /* include shadow tree event handler if it is not included yet */
   eventhdlr = SCIPfindEventhdlr(scip, "event_shadowtree");
   if( eventhdlr == NULL )
   {
      SCIP_CALL( SCIPincludeEventHdlrShadowTree(scip, &symhdlrdata->shadowtreeeventhdlr) );
      assert(symhdlrdata->shadowtreeeventhdlr != NULL);
   }
   else
      symhdlrdata->shadowtreeeventhdlr = eventhdlr;

   SCIP_CALL( SCIPincludeOrbitalReduction(scip, &symhdlrdata->orbitalreddata,
         symhdlrdata->shadowtreeeventhdlr) );
   assert(symhdlrdata->orbitalreddata != NULL);

   /* add parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/uselexicographicreduction",
         "Shall the lexicographic reduction algorithm be used?",
         &symhdlrdata->uselexred, TRUE, DEFAULT_USELEXRED, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/useorbitalreduction",
         "Shall the orbital reduction algorithm be used?",
         &symhdlrdata->useorbred, TRUE, DEFAULT_USEORBRED, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/computenewperms",
         "Shall additional permutations of symmetry component be computed?",
         &symhdlrdata->computenewperms, TRUE, DEFAULT_COMPUTENEWPERMS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "symmetries/" SYM_NAME "/maxnnewperms",
         "maximum number of additional permutations of symmetry component that is computed " \
         "(used to have better approximation of stabilizer); -1: unbounded",
         &symhdlrdata->maxnnewperms, TRUE, DEFAULT_MAXNNEWPERMS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/checkppupgrade",
         "Shall static constraints be used if at least 50 percent of perms can be handled by pp-orbisacks?",
         &symhdlrdata->useppupgrade, TRUE, DEFAULT_USEPPPUPGRADE, NULL, NULL) );

   return SCIP_OKAY;
}
