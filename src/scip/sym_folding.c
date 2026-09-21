/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   sym_folding.c
 * @ingroup DEFPLUGINS_SYM
 * @brief  symmetry handler for folding-based reductions
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#include "scip/sym_folding.h"
#include "scip/symmetry.h"


/* fundamental symmetry handler properties */
#define SYM_NAME                 "folding"
#define SYM_DESC                 "symmetry handler for folding based reductions"
#define SYM_PRIORITY             50000       /**< priority of try-add function*/

/* optional symmetry handler properties */
#define SYM_PRESOLPRIORITY       0 /**< priority of the symmetry handler for presolving */
#define SYM_PRESOLTIMING    SCIP_PRESOLTIMING_FAST /**< presolving timing of the symmetry handler (fast, medium, or exhaustive) */
#define SYM_MAXPREROUNDS         1 /**< maximal number of presolving rounds the symmetry handler participates in (-1: no limit) */

/*
 * Data structures
 */

/** symmetry component data for components handled by folding symmetry handler */
struct SCIP_SymCompData
{
   SCIP_CONS**           conss;              /**< list of generated constraints */
   int                   nconss;             /**< number of generated constraints */
};

/*
 * Local methods
 */

/** checks whether all affected variables are continuous and possibly returns those */
static
SCIP_RETCODE collectSymvars(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 perms,              /**< (proper) permutations */
   int                   nperms,             /**< number of permutations */
   SCIP_VAR**            permvars,           /**< variables the permutations act on */
   int                   npermvars,          /**< number of variables the permutations act on */
   SCIP_VAR***           affectedvars,       /**< pointer to store affected variables */
   int*                  naffectedvars,      /**< pointer to store number of affected variables */
   SCIP_Bool*            success             /**< pointer to store whether method ran successfully */
   )
{
   int v;
   int p;
   int i;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms > 0);
   assert(permvars != NULL);
   assert(npermvars > 0);
   assert(affectedvars != NULL);
   assert(naffectedvars != NULL);
   assert(success != NULL);

   *success = FALSE;
   *naffectedvars = 0;

   /* terminate early if there are not enough continuous variables */
   if( SCIPgetNContVars(scip) <= 1 )
      return SCIP_OKAY;

   /* check whether all affected variables are continuous and count them */
   for( v = 0; v < npermvars; ++v )
   {
      for( p = 0; p < nperms; ++p )
      {
         if( perms[p][v] != v )
         {
            if( SCIPvarGetType(permvars[v]) != SCIP_VARTYPE_CONTINUOUS )
               return SCIP_OKAY;
            else
            {
               ++(*naffectedvars);
               break;
            }
         }
      }
   }
   assert(*naffectedvars > 0);

   /* collect affected variables
    * Don't collect this already above, because in many cases the check will fail and memory allocation is expensive. */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, affectedvars, *naffectedvars) );
   for( v = 0, i = 0; v < npermvars; ++v )
   {
      for( p = 0; p < nperms; ++p )
      {
         if( perms[p][v] != v )
         {
            (*affectedvars)[i++] = permvars[v];
            break;
         }
      }
   }
   assert(*naffectedvars == i);

   *success = TRUE;

   return SCIP_OKAY;
}

/** returns whether a constraint is classified as convex */
static
SCIP_Bool isConsConvex(
   SCIP_CONS*            cons                /**< constraint to be checked for convexity */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(cons != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);

   /* treat nonlinear constraints differently, because their convexity might not be known at creation */
   if( strcmp(SCIPconshdlrGetName(conshdlr), "nonlinear") == 0
      && (SCIPgetCurvatureNonlinear(cons) == SCIP_EXPRCURV_CONVEX
         || SCIPgetCurvatureNonlinear(cons) == SCIP_EXPRCURV_LINEAR ) )
      return TRUE;

   return SCIPconsIsRelaxConvex(cons);
}

/** checks whether all constraints containing an affected variable are convex */
static
SCIP_RETCODE checkRelevantConssAreConvex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            affectedvars,       /**< affected variables */
   int                   naffectedvars,      /**< number of affected variables */
   SCIP_Bool*            success             /**< pointer to store whether all relevant constraints are convex */
   )
{
   SCIP_HASHSET* affvarsset = NULL;
   SCIP_CONS** conss;
   SCIP_VAR** locvars;
   int nlocvars;
   int nvars;
   int nconss;
   int c;
   int v;

   assert(scip != NULL);
   assert(affectedvars != NULL);
   assert(naffectedvars > 0);
   assert(success != NULL);

   *success = FALSE;

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   if( nconss <= 0 )
      return SCIP_OKAY;
   assert(conss != NULL);

   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &locvars, nvars) );

   *success = TRUE;
   for( c = 0; c < nconss && *success; ++c )
   {
      /* check whether non-convex constraints do not contain affected variables */
      if( isConsConvex(conss[c]) )
         continue;

      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nlocvars, success) );

      if( !(*success) )
         break;

      SCIP_CALL( SCIPgetConsVars(scip, conss[c], locvars, nvars, success) );
      assert(*success);

      /* store the affected variables in a hashset for checking overlap with locvars */
      if( affvarsset == NULL )
      {
         SCIP_CALL( SCIPhashsetCreate(&affvarsset, SCIPblkmem(scip), naffectedvars) );
         for( v = 0; v < naffectedvars; ++v )
         {
            SCIP_CALL( SCIPhashsetInsert(affvarsset, SCIPblkmem(scip), affectedvars[v]) );
         }
      }

      /* check whether a variable from the other array is contained in the hashset */
      for( v = 0; v < nlocvars; ++v )
      {
         if( SCIPhashsetExists(affvarsset, locvars[v]) )
         {
            *success = FALSE;
            break;
         }
      }
   }

   if( affvarsset != NULL )
   {
      SCIPhashsetFree(&affvarsset, SCIPblkmem(scip));
   }

   SCIPfreeBufferArray(scip, &locvars);

   return SCIP_OKAY;
}

/** applies the folding step, i.e., imposes that symmetric variables take same value */
static
SCIP_RETCODE applyFolding(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 perms,              /**< permutations */
   int                   nperms,             /**< number of permutations */
   SCIP_VAR**            permvars,           /**< variables the permutations act on */
   int                   npermvars,          /**< number of variables the permutations act on */
   SCIP_SYMCOMPDATA**    symcompdata,        /**< pointer to store data of symmetry component */
   int                   id,                 /**< ID of symmetry component */
   SCIP_Bool             addconss,           /**< Shall constraints be added to fold vars? (otherwise: aggregation) */
   int*                  naddedconss,        /**< pointer to store number of added constraints */
   int*                  naggrvars           /**< pointer to store number of aggregated variables */
   )
{
   int* orbits;
   int* orbitbegins;
   int norbits;
   int first;
   int nvars;
   int oidx;
   int i;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms > 0);
   assert(permvars != NULL);
   assert(npermvars > 0);
   assert(symcompdata != NULL);
   assert(naddedconss != NULL);
   assert(naggrvars != NULL);

   *naddedconss = 0;
   *naggrvars = 0;
   nvars = SCIPgetNVars(scip);
   assert(nvars >= npermvars);

   /* compute orbits of symmetry component */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, nvars) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, FALSE, permvars, npermvars, perms, nperms, orbits, orbitbegins, &norbits) );
   assert(norbits > 0);

   if( addconss )
   {
      /* add constraints to apply folding */
      char name[SCIP_MAXSTRLEN];
      SCIP_VAR* tmpvars[2] = { NULL, NULL };
      SCIP_Real coefs[2] = { 1.0, -1.0 };
      int nconss;
      int c = 0;

      /* compute number of conss to be added (per orbit, the representative does not need to be folded) */
      nconss = orbitbegins[norbits] - norbits;

      /* allocate memory for component data and create constraints  */
      SCIP_CALL( SCIPallocBlockMemory(scip, symcompdata) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*symcompdata)->conss, nconss) );
      (*symcompdata)->nconss = nconss;

      for( oidx = 0; oidx < norbits; ++oidx )
      {
         first = orbits[orbitbegins[oidx]];
         tmpvars[0] = permvars[first];

         for( i = orbitbegins[oidx] + 1; i < orbitbegins[oidx + 1]; ++i )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "folding_symcomp_%d_%d_%d_%d", id, oidx, first, orbits[i]);

            tmpvars[1] = permvars[orbits[i]];
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &(*symcompdata)->conss[c], name, 2, tmpvars, coefs, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, (*symcompdata)->conss[c++]) );
            assert(c <= nconss);
         }
      }
      assert(c == nconss);

      *naddedconss = nconss;
   }
   else
   {
      /* aggregate variables to apply folding */
      SCIP_Bool infeasible;
      SCIP_Bool redundant;
      SCIP_Bool aggregated;

      for( oidx = 0; oidx < norbits; ++oidx )
      {
         first = orbits[orbitbegins[oidx]];

         for( i = orbitbegins[oidx] + 1; i < orbitbegins[oidx + 1]; ++i )
         {
            SCIP_CALL( SCIPaggregateVars(scip, permvars[first], permvars[orbits[i]], 1.0, -1.0, 0.0,
                  &infeasible, &redundant, &aggregated) );
            if( aggregated )
               ++(*naggrvars);
         }
      }
      assert(*naggrvars <= orbitbegins[norbits] - norbits);
   }

   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   return SCIP_OKAY;
}

/*
 * Callback methods of symmetry handler
 */

/** addition method for symmetry method handler plugins (tries to add symmetry handling method for given symmetries) */
static
SCIP_DECL_SYMHDLRTRYADD(symTryaddFolding)
{  /*lint --e{715}*/
   SCIP_VAR** symvars;
   int nsymvars;
   SCIP_Bool addconss;

   assert(success != NULL);
   assert(naddedconss != NULL);
   assert(naggrvars != NULL);

   *success = FALSE;

   /* only apply reductions for permutation symmetries */
   if( symtype != SYM_SYMTYPE_PERM )
      return SCIP_OKAY;

   /* check whether all affected variables are continuous and possibly collect these variables */
   SCIP_CALL( collectSymvars(scip, perms, nperms, permvars, npermvars, &symvars, &nsymvars, success) );
   if( !(*success) )
      return SCIP_OKAY;

   /* check whether all constraints that contain an affected variable are convex */
   SCIP_CALL( checkRelevantConssAreConvex(scip, symvars, nsymvars, success) );
   if( !(*success) )
      goto FREEMEM;

   /* apply folding */
   addconss = (SCIPgetStage(scip) <= SCIP_STAGE_INITPRESOLVE);
   SCIP_CALL( applyFolding(scip, perms, nperms, permvars, npermvars, symcompdata,
         id, addconss, naddedconss, naggrvars) );
   *success = TRUE;

 FREEMEM:
   SCIPfreeBlockMemoryArray(scip, &symvars, nsymvars);

   return SCIP_OKAY;
}

/** copy method for symmetry handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SYMHDLRCOPY(symCopyFolding)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(symhdlr != NULL);
   assert(strcmp(SCIPsymhdlrGetName(symhdlr), SYM_NAME) == 0);

   /* call inclusion method of symmetry handler */
   SCIP_CALL( SCIPincludeSymhdlrFolding(scip) );

   return SCIP_OKAY;
}

/** deinitialization method of symmetry handler (called before transformed problem is freed) */
static
SCIP_DECL_SYMHDLREXIT(symExitFolding)
{  /*lint --e{715}*/
   SCIP_SYMCOMPDATA* symdata;
   int s;
   int c;

   for( s = 0; s < nsymcomps; ++s )
   {
      symdata = SCIPsymcompGetData(symcomps[s]);
      assert(symdata->conss == NULL || symdata->nconss > 0);

      for( c = 0; c < symdata->nconss; ++c )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &symdata->conss[c]) );
      }
      SCIPfreeBlockMemoryArray(scip, &symdata->conss, symdata->nconss);
      SCIPfreeBlockMemory(scip, &symdata);
   }

   return SCIP_OKAY;
}

/** presolving method of symmetry handler */
static
SCIP_DECL_SYMHDLRPRESOL(symPresolFolding)
{  /*lint --e{715}*/
   SCIP_SYMCOMPDATA* symdata;
   int s;
   int c;

   assert(result != NULL);

   *result = nsymcomps > 0 ? SCIP_DIDNOTFIND : SCIP_DIDNOTRUN;

   for( s = 0; s < nsymcomps; ++s )
   {
      symdata = SCIPsymcompGetData(symcomps[s]);

      for( c = 0; c < symdata->nconss; ++c )
      {
         if( !SCIPconsIsActive(symdata->conss[c]) )
            continue;

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

/*
 * constraint specific interface methods
 */

/** creates the handler for folding symmetry handler and includes it in SCIP */
SCIP_RETCODE SCIPincludeSymhdlrFolding(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SYMHDLRDATA* symhdlrdata = NULL;
   SCIP_SYMHDLR* symhdlr;

   assert(scip != NULL);

   SCIP_CALL( SCIPincludeSymhdlrBasic(scip, &symhdlr, SYM_NAME, SYM_DESC, SYM_PRIORITY,
         symTryaddFolding, symhdlrdata) );

   /* include non-fundamental callback methods */
   SCIP_CALL( SCIPsetSymhdlrCopy(scip, symhdlr, symCopyFolding) );
   SCIP_CALL( SCIPsetSymhdlrExit(scip, symhdlr, symExitFolding) );
   SCIP_CALL( SCIPsetSymhdlrPresol(scip, symhdlr, symPresolFolding,
         SYM_MAXPREROUNDS, SYM_PRESOLPRIORITY, SYM_PRESOLTIMING) );

   return SCIP_OKAY;
}
