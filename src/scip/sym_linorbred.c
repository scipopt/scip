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
 * @brief  symmetry handler for linear reduction and orbital reduction
 * @author Christopher Hojny
 * @author Jasper van Doornmalen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_shadowtree.h"
#include "scip/pub_implics.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sym.h"
#include "scip/pub_var.h"
#include "scip/scip_event.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sym.h"
#include "scip/scip_var.h"
#include "scip/sym_linorbred.h"
#include "scip/symmetry.h"
#include "scip/symmetry_orbital.h"
#include "scip/symmetry_lexred.h"
#include "scip/type_implics.h"

/* symmetry handler properties */
#define SYM_NAME            "sym_linorbred"
#define SYM_DESC            "symmetry handler for linear reduction and orbital reduction"
#define SYM_PRIORITY            -10000       /**< priority of try-add function*/
#define SYM_PROPPRIORITY       -100000       /**< priority of propagation method */
#define SYM_PROPTIMING SCIP_PROPTIMING_BEFORELP /**< timing of propagator method */
#define SYM_PROPFREQ                 1       /**< frequence of propagator method */
#define SYM_DELAYPROP            FALSE       /**< Should propagation method be delayed, if other propagators found reductions? */

/* default value of parameters */
#define DEFAULT_USELINRED         TRUE       /**< Shall linear reduction be used? */
#define DEFAULT_USEORBRED         TRUE       /**< Shall orbital reduction be used? */
#define DEFAULT_COMPUTENEWPERMS   TRUE       /**< Shall additional permutations of symmetry component be computed? */
#define DEFAULT_MAXNNEWPERMS       100       /**< maximum number of additional permutations of symmetry component
                                              *   that is computed (used to have better approximation of stabilizer);
                                              *   -1: unbounded */

/*
 * Data structures
 */

/** symmetry component data */
struct SCIP_SymCompData
{
   SCIP_ORBITALREDDATA*  orbitalreddata;     /**< container for orbital reduction data */
   SCIP_LEXREDDATA*      lexreddata;         /**< container for lexicographic reduction propagation; */
};

/** symmetry handler data */
struct SCIP_SymhdlrData
{
   SCIP_EVENTHDLR*       shadowtreeeventhdlr; /**< shadow tree event handler */
   SCIP_Bool             uselinred;          /**< Shall linear reduction be used? */
   SCIP_Bool             useorbred;          /**< Shall orbital reduction be used? */
   SCIP_Bool             computenewperms;    /**< Shall additional permutations of symmetry component be computed? */
   int                   maxnnewperms;       /**< maximum number of additional permutations of symmetry component
                                              *   that is computed (used to have better approximation of stabilizer);
                                              *   -1: unbounded */
};

/*
 * Local methods
 */

/* @symtodo avoid code duplication with sym_sst */
/** returns whether a permutation is already contained in a list of permutations */
static
SCIP_Bool isPermKnown(
   int*                  perm,               /**< permutation to be checked */
   int                   permlen,            /**< length of permutation */
   int**                 knownperms,         /**< list of known permutations (possibly longer than nknownperms) */
   int                   nknownperms         /**< number of known permutations to be checked */
   )
{
   int p;
   int i;

   assert(perm != NULL);
   assert(permlen >= 0);
   assert(knownperms != NULL);
   assert(nknownperms >= 0);

   for( p = 0; p < nknownperms; ++p )
   {
      for( i = 0; i < permlen; ++i )
      {
         /* knownperms[p] and perm differ */
         if( perm[i] != knownperms[p][i] )
            break;
      }
      /* loop did not terminate early, knownperms[p] and perm coincide */
      if( i == permlen )
         return TRUE;
   }

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

/*
 * Callback methods of symmetry handler
 */

/** addition method for symmetry method handler plugins (tries to add symmetry handling method for given symmetries) */
static
SCIP_DECL_SYMHDLRTRYADD(symhdlrTryaddLinOrbRed)
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

   *success = FALSE;

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   if( !symhdlrdata->uselinred && !symhdlrdata->useorbred )
      return SCIP_OKAY;

   /* possibly create new permutations */
   if( symhdlrdata->computenewperms )
   {
      SCIP_CALL( tryGenerateInvolutions(scip, symtype, perms, nperms, npermvars, symhdlrdata->maxnnewperms,
            &newperms, &nnewperms, &maxnnewperms) );
   }
   assert(symhdlrdata->shadowtreeeventhdlr != NULL);


   SCIP_CALL( SCIPallocBlockMemory(scip, symcompdata) );

   if( symhdlrdata->uselinred )
   {
      SCIP_Bool locsuccess;
      SCIP_Real* permvardomaincenter;

      SCIP_CALL( SCIPallocBufferArray(scip, &permvardomaincenter, npermvars) );
      for( p = 0; p < npermvars; ++p )
         permvardomaincenter[p] = (SCIPvarGetUbLocal(permvars[p]) + SCIPvarGetLbLocal(permvars[p])) / 2;

      SCIP_CALL( SCIPincludeLexicographicReduction(scip, &(*symcompdata)->lexreddata,
            symhdlrdata->shadowtreeeventhdlr) );
      assert((*symcompdata)->lexreddata != NULL);

      /* add all permutations to lexicograph reduction */
      for( p = 0; p < nperms; ++p )
      {
         SCIP_CALL( SCIPlexicographicReductionAddPermutation(scip, (*symcompdata)->lexreddata,
               permvars, npermvars, perms[p], (SYM_SYMTYPE) symtype, permvardomaincenter, TRUE, &locsuccess) );
         *success |= locsuccess;
      }
      for( p = 0; p < nnewperms; ++p )
      {
         SCIP_CALL( SCIPlexicographicReductionAddPermutation(scip, (*symcompdata)->lexreddata,
               permvars, npermvars, newperms[p], (SYM_SYMTYPE) symtype, permvardomaincenter, TRUE, &locsuccess) );
         *success |= locsuccess;
      }

      SCIPfreeBufferArray(scip, &permvardomaincenter);
   }


   /* free new permutations */
   for( p = nnewperms - 1; p >= 0; --p )
   {
      SCIPfreeBlockMemoryArray(scip, &newperms[p], symtype == SYM_SYMTYPE_PERM ? npermvars : 2 * npermvars);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &newperms, maxnnewperms);

   return SCIP_OKAY;
}

/** deinitialization method of symmetry handler (called before transformed problem is freed) */
static
SCIP_DECL_SYMHDLREXIT(symhdlrExitLinOrbRed)
{  /*lint --e{715}*/
   SCIP_SYMCOMPDATA* symdata;
   int s;

   for( s = 0; s < nsymcomps; ++s )
   {
      symdata = SCIPsymcompGetData(symcomps[s]);

      if( symdata->lexreddata != NULL )
      {
         SCIP_CALL( SCIPlexicographicReductionReset(scip, symdata->lexreddata) );
         SCIP_CALL( SCIPlexicographicReductionFree(scip, &symdata->lexreddata) );
      }

      SCIPfreeBlockMemory(scip, &symdata);
   }

   return SCIP_OKAY;
}

/** destructor of symmetry handler to free symmetry handler data (called when SCIP is exiting) */
static
SCIP_DECL_SYMHDLRFREE(symhdlrFreeLinOrbRed)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;

   assert(scip != NULL);
   assert(symhdlr != NULL);

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &symhdlrdata);

   return SCIP_OKAY;
}

/** domain propagation method of symmetry handler */
static
SCIP_DECL_SYMHDLRPROP(symhdlrPropLinOrbRed)
{  /*lint --e{715}*/
   SCIP_SYMCOMPDATA* symcompdata;
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool didrun = FALSE;
   SCIP_Bool didrunlocal;
   int nredlocal;
   int nreds = 0;
   int s;

   assert(symcomps != NULL || nsymcomps == 0);
   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   for( s = 0; s < nsymcomps; ++s )
   {
      assert(symcomps[s] != NULL);

      symcompdata = SCIPsymcompGetData(symcomps[s]);
      assert(symcompdata != NULL);

      SCIP_CALL( SCIPlexicographicReductionPropagate(scip, symcompdata->lexreddata,
            &infeasible, &nredlocal, &didrunlocal) );
      nreds += nredlocal;
      didrun |= didrunlocal;
      if ( infeasible )
         return SCIP_OKAY;
   }
   if( nreds > 0 || infeasible )
      printf("nreds %d, inf %d\n", nreds, infeasible);

   if( infeasible )
      *result = SCIP_CUTOFF;
   else if( nreds > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( didrun )
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** include symmetry handler for linear reduction and orbital reduction */
SCIP_RETCODE SCIPincludeSymhdlrLinOrbRed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SYMHDLRDATA* symhdlrdata = NULL;
   SCIP_EVENTHDLR* eventhdlr;

   assert(scip != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &symhdlrdata) );

   SCIP_CALL( SCIPincludeSymhdlrBasic(scip, SYM_NAME, SYM_DESC, SYM_PRIORITY, SYM_PROPPRIORITY, 0, -1,
         SYM_PROPFREQ, -1, SYM_DELAYPROP, FALSE, 1.0, 0, SYM_PROPTIMING, SCIP_PRESOLTIMING_FAST,
         symhdlrTryaddLinOrbRed, NULL, symhdlrFreeLinOrbRed, NULL, symhdlrExitLinOrbRed,
         NULL, NULL, NULL, NULL, NULL, NULL, symhdlrPropLinOrbRed,
         NULL, NULL, symhdlrdata) );

   /* include shadow tree event handler if it is not included yet */
   eventhdlr = SCIPfindEventhdlr(scip, "event_shadowtree");
   if( eventhdlr == NULL )
   {
      SCIP_CALL( SCIPincludeEventHdlrShadowTree(scip, &symhdlrdata->shadowtreeeventhdlr) );
      assert(symhdlrdata->shadowtreeeventhdlr != NULL);
   }
   else
      symhdlrdata->shadowtreeeventhdlr = eventhdlr;

   /* add parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/uselinearreduction",
         "Shall the linear reduction algorithm be used?",
         &symhdlrdata->uselinred, TRUE, DEFAULT_USELINRED, NULL, NULL) );

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

   return SCIP_OKAY;
}
