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

/**@file   heur_indtwoopt.c
 * @ingroup PRIMALHEURISTICS
 * @brief  modified twoopt primal heuristic using indicator variables
 * @author Timo Berthold
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_indtwoopt.h"

#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"


#define HEUR_NAME             "indtwoopt"
#define HEUR_DESC             "modified 2-opt heuristic which tries to improve setting of single integer variables"
#define HEUR_DISPCHAR         '7'
#define HEUR_PRIORITY         -200000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_LARGESTWEIGHT FALSE
#define MAXCANDS              200

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   lastsolindex;       /**< index of the last solution for which indtwoopt was performed */
   SCIP_Bool             largestweight;      /**< should we start with the variables of largest weight? */
   SCIP_Bool             triedstart;         /**< whether we have tried to generate a starting solution */
   SCIP_CONSHDLR*        indicatorconshdlr;  /**< constraint handler for indicator constraints (if available) */
   int                   nindconss;          /**< number of indicator constraints */
   int                   nvars;              /**< number of variables */
   int                   nnonz;              /**< number of entries in indidx array */
   int                   maxnnonz;           /**< maximal number of entries in indidx array */
   SCIP_CONS**           indcons;            /**< indicator constraints that contain variables */
   int*                  indbeg;             /**< start of list of indicator constraints indices */
   SCIP_Bool*            isindvar;           /**< whether variable is an binary/slack variable of an indicator constraint */
};


/*
 * Local methods
 */

/** perform one try of changing two variables */
static
SCIP_RETCODE tryValue(
   SCIP*                 scip,               /**< scip pointer */
   SCIP_HEURDATA*        heurdata,           /**< data of heuristic */
   SCIP_Real             bestobj,            /**< best objective so far */
   SCIP_SOL*             sol,                /**< solution to be changed */
   SCIP_VAR*             var1,               /**< variable 1 to be changed */
   SCIP_Real             val1,               /**< old value of variable 1 */
   SCIP_Real             newval1,            /**< new value of variable 1 */
   SCIP_VAR*             var2,               /**< variable 2 to be changed */
   SCIP_Real             val2,               /**< old value of variable 2 */
   SCIP_Real             newval2,            /**< new value of variable 2 */
   SCIP_Bool*            success             /**< whether objective could be improved */
   )
{
   SCIP_Bool changed;
   SCIP_Real obj;
   int idx;
   int j;

   /* tentatively set variables to new values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, var1, newval1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, var2, newval2) );

   /* try to fix indicators for variable 1 */
   idx = SCIPvarGetProbindex(var1);
   assert( 0 <= idx && idx < heurdata->nvars );
   assert( heurdata->indbeg[idx] <= heurdata->indbeg[idx+1] );
   assert( 0 <= heurdata->indbeg[idx] && heurdata->indbeg[idx] <= heurdata->nnonz );
   for (j = heurdata->indbeg[idx]; j < heurdata->indbeg[idx+1]; ++j)
   {
      assert( heurdata->indcons[j] != NULL );
      SCIP_CALL( SCIPmakeIndicatorFeasible(scip, heurdata->indcons[j], sol, &changed) );
   }

   /* try to fix indicators for variable 2 */
   idx = SCIPvarGetProbindex(var2);
   assert( 0 <= idx && idx < heurdata->nvars );
   assert( heurdata->indbeg[idx] <= heurdata->indbeg[idx+1] );
   assert( 0 <= heurdata->indbeg[idx] && heurdata->indbeg[idx] <= heurdata->nnonz );
   for (j = heurdata->indbeg[idx]; j < heurdata->indbeg[idx+1]; ++j)
   {
      assert( heurdata->indcons[j] != NULL );
      SCIP_CALL( SCIPmakeIndicatorFeasible(scip, heurdata->indcons[j], sol, &changed) );
   }

   obj = SCIPgetSolOrigObj(scip, sol);

   /* check solution */
   *success = FALSE;
   if ( SCIPisLT(scip, obj, bestobj) )
   {
      SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );
   }

   if ( ! *success )
   {
      /* reset value */
      SCIP_CALL( SCIPsetSolVal(scip, sol, var2, val2) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, var1, val1) );

      idx = SCIPvarGetProbindex(var1);
      assert( heurdata->indbeg[idx] <= heurdata->indbeg[idx+1] );
      assert( 0 <= heurdata->indbeg[idx] && heurdata->indbeg[idx] <= heurdata->nnonz );
      for (j = heurdata->indbeg[idx]; j < heurdata->indbeg[idx+1]; ++j)
      {
         assert( heurdata->indcons[j] != NULL );
         SCIP_CALL( SCIPmakeIndicatorFeasible(scip, heurdata->indcons[j], sol, &changed) );
      }

      idx = SCIPvarGetProbindex(var2);
      assert( heurdata->indbeg[idx] <= heurdata->indbeg[idx+1] );
      assert( 0 <= heurdata->indbeg[idx] && heurdata->indbeg[idx] <= heurdata->nnonz );
      for (j = heurdata->indbeg[idx]; j < heurdata->indbeg[idx+1]; ++j)
      {
         assert( heurdata->indcons[j] != NULL );
         SCIP_CALL( SCIPmakeIndicatorFeasible(scip, heurdata->indcons[j], sol, &changed) );
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyIndtwoopt)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurIndtwoopt(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeIndtwoopt)
{   /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );
   assert( heurdata->indcons == NULL );
   assert( heurdata->indbeg == NULL );
   assert( heurdata->isindvar == NULL );

   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
#define heurInitIndtwoopt NULL

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitIndtwoopt NULL

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolIndtwoopt)
{
   SCIP_HEURDATA* heurdata;
   SCIP_CONS** conss;
   SCIP_Bool* isindvar;
   SCIP_CONS** indcons;
   int* indbeg;
   int* varidx;
   int nindconss;
   int maxnnonz;
   int nnonz;
   int nconss;
   int nvars;
   int c;
   int v;

   assert( scip != NULL );
   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   SCIPdebugMessage("Initialize solution process ...\n");

   /* set up heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert( heurdata->indcons == NULL );
   assert( heurdata->indbeg == NULL );
   assert( heurdata->isindvar == NULL );

   heurdata->lastsolindex = -1;
   heurdata->triedstart = FALSE;

   /* get indicator constraint handler */
   heurdata->indicatorconshdlr = SCIPfindConshdlr(scip, "indicator");

   /* exit if indicator constraint handler is not present */
   if ( heurdata->indicatorconshdlr == NULL )
      return SCIP_OKAY;

   /* get data */
   nconss = SCIPconshdlrGetNConss(heurdata->indicatorconshdlr);
   conss = SCIPconshdlrGetConss(heurdata->indicatorconshdlr);
   nvars = SCIPgetNVars(scip);

   /* init data */
   heurdata->nindconss = nconss;  /* this number might be corrected below */
   heurdata->nvars = nvars;
   heurdata->nnonz = 0;
   heurdata->maxnnonz = 0;

   /* exit if there are no indicator constraints or no variables */
   if ( nconss == 0 || nvars == 0 )
      return SCIP_OKAY;

   /* set up storage */
   maxnnonz = 10 * nconss;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &varidx, maxnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &indcons, maxnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &isindvar, nvars) );

   /* init array */
   for (v = 0; v < nvars; ++v)
      isindvar[v] = FALSE;

   nnonz = 0;
   nindconss = 0;
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONS* cons;
      SCIP_CONS* lincons;
      SCIP_VAR** linvars;
      int nlinvars;
      SCIP_VAR* var;
      int idx;

      cons = conss[c];
      assert( cons != NULL );

      /* avoid non-active indicator constraints */
      if ( ! SCIPconsIsActive(cons) )
         continue;

      ++nindconss;
      var = SCIPgetSlackVarIndicator(cons);
      idx = SCIPvarGetProbindex(var);
      /* avoid removed slack variables (may happen if presolving/propagation) was incomplete */
      if ( idx < 0 )
         continue;
      assert( idx < nvars );
      isindvar[idx] = TRUE;

      var = SCIPgetBinaryVarIndicator(cons);
      if ( SCIPvarIsNegated(var) )
         var = SCIPvarGetNegatedVar(var);
      idx = SCIPvarGetProbindex(var);
      /* avoid removed indicator variables (may happen if presolving/propagation) was incomplete */
      if ( idx < 0 )
         continue;
      assert( idx < nvars );
      isindvar[idx] = TRUE;

      /* fill matrix data */
      lincons = SCIPgetLinearConsIndicator(cons);
      nlinvars = SCIPgetNVarsLinear(scip, lincons);
      linvars = SCIPgetVarsLinear(scip, lincons);

      for (v = 0; v < nlinvars; ++v)
      {
         if ( nnonz >= maxnnonz )
         {
            int newsize = SCIPcalcMemGrowSize(scip, maxnnonz+1);
            assert( newsize > maxnnonz );
            SCIP_CALL_ABORT( SCIPreallocBlockMemoryArray(scip, &varidx, maxnnonz, newsize) );
            SCIP_CALL_ABORT( SCIPreallocBlockMemoryArray(scip, &indcons, maxnnonz, newsize) );
            maxnnonz = newsize;
         }

         idx = SCIPvarGetProbindex(linvars[v]);
         /* skip non active variables */
         if ( idx >= 0 )
         {
            varidx[nnonz] = idx;
            assert( varidx[nnonz] >= 0 );
            indcons[nnonz++] = cons;

            /* capture cons, since we have to access it later */
            SCIP_CALL( SCIPcaptureCons(scip, cons) );
         }
      }
   }
   assert( nnonz <= maxnnonz );
   assert( nindconss <= nconss );

   /* sort matrix */
   SCIPsortIntPtr(varidx, (void*) indcons, nnonz);

   /* fill in beg-array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &indbeg, nvars+1) );

   c = 0;
   for (v = 0; v < nvars; ++v)
   {
      indbeg[v] = c;
      while ( c < nnonz && varidx[c] == v )
         ++c;
   }
   assert( indbeg[0] == 0 );
   assert( c <= nnonz );
   indbeg[nvars] = nnonz;

   /* store data */
   heurdata->nindconss = nindconss;
   heurdata->nnonz = nnonz;
   heurdata->maxnnonz = maxnnonz;
   heurdata->indcons = indcons;
   heurdata->indbeg = indbeg;
   heurdata->isindvar = isindvar;

   /* free storage */
   SCIPfreeBlockMemoryArray(scip, &varidx, maxnnonz);

   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolIndtwoopt)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( scip != NULL );
   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   SCIPdebugMessage("Exit solution process ...\n");

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   if ( heurdata->indcons != NULL )
   {
      int j;

      for (j = 0; j < heurdata->nnonz; ++j)
      {
         assert( heurdata->indcons[j] != NULL );
         SCIP_CALL( SCIPreleaseCons(scip, &heurdata->indcons[j]) );
      }

      SCIPfreeBlockMemoryArray(scip, &heurdata->indbeg, heurdata->nvars + 1);
      SCIPfreeBlockMemoryArray(scip, &heurdata->indcons, heurdata->maxnnonz);
      SCIPfreeBlockMemoryArray(scip, &heurdata->isindvar, heurdata->nvars);
      heurdata->maxnnonz = 0;
      heurdata->nnonz = 0;
   }
   assert( heurdata->indbeg == NULL );

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecIndtwoopt)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;
   SCIP_SOL* bestsol;                        /* incumbent solution */
   SCIP_SOL* worksol;                        /* heuristic's working solution */
   SCIP_VAR** vars;                          /* SCIP variables                */
   SCIP_VAR** candvars;
   SCIP_Real* candvals;
   SCIP_Bool success;
   SCIP_Bool changed;
   SCIP_Real lb;
   SCIP_Real ub;
   int ncands;

   int nchgbound;
   int nbinvars;
   int nintvars;
   int nimplvars;
   int nvars;
   int i;
   int k;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   *result = SCIP_DIDNOTRUN;

   /* exit if no indicator constraints are present */
   if ( heurdata->indicatorconshdlr == NULL )
      return SCIP_OKAY;

   /* exit if there are no indicator constraints */
   if ( heurdata->indicatorconshdlr == NULL || heurdata->nindconss == 0 || heurdata->nvars == 0 )
      return SCIP_OKAY;

   /* get problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, &nimplvars, NULL) );
   nintvars += nbinvars + nimplvars;

   /* do not run if there are no discrete variables */
   if ( nintvars == 0 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* check sizes */
   assert( nvars == heurdata->nvars );

   *result = SCIP_DELAYED;
   worksol = NULL;

   /* we can only work on solutions valid in the transformed space */
   bestsol = SCIPgetBestSol(scip);
   if ( bestsol != NULL && SCIPsolGetOrigin(bestsol) != SCIP_SOLORIGIN_ORIGINAL )
   {
      /* we only want to process each solution once */
      if ( heurdata->lastsolindex == SCIPsolGetIndex(bestsol) )
         return SCIP_OKAY;

      SCIPdebugMessage("Starting indicator 2-opt heuristic ...\n");

      /* initialize data */
      heurdata->lastsolindex = SCIPsolGetIndex(bestsol);
      SCIP_CALL( SCIPcreateSolCopy(scip, &worksol, bestsol) );
      SCIPsolSetHeur(worksol, heur);
      SCIPdebugMessage("Using solution of value %g.\n", SCIPgetSolOrigObj(scip, worksol) );
   }
   else
   {
      if ( heurdata->triedstart )
         return SCIP_OKAY;

      SCIPdebugMessage("Starting indicator 2-opt heuristic ...\n");

      /* create 0-solution */
      SCIP_CALL( SCIPcreateSol(scip, &worksol, heur) );
      heurdata->triedstart = TRUE;
      /* try to fix indicators */
      SCIP_CALL( SCIPmakeIndicatorsFeasible(scip, heurdata->indicatorconshdlr, worksol, &changed) );
      SCIPdebugMessage("Trying starting solution of value %g.\n", SCIPgetSolOrigObj(scip, worksol));
      SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, TRUE, TRUE, FALSE, &success) );
   }
   assert( worksol != NULL );

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("Starting bound adjustment in 2-opt heuristic (indicator version) ... \n");

   /* maybe change solution values due to global bound changes first */
   nchgbound = 0;
   for (i = nvars - 1; i >= 0; --i)
   {
      SCIP_VAR* var;
      SCIP_Real solval;

      var = vars[i];
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      solval = SCIPgetSolVal(scip, worksol, var);

      /* old solution value is smaller than the actual lower bound */
      if ( SCIPisFeasLT(scip, solval, lb) )
      {
         /* set the solution value to the global lower bound */
         SCIP_CALL( SCIPsetSolVal(scip, worksol, var, lb) );
         ++nchgbound;
         SCIPdebugMessage("var <%s> type %d, old solval %g now fixed to lb %g\n", SCIPvarGetName(var), SCIPvarGetType(var), solval, lb);
      }
      /* old solution value is greater than the actual upper bound */
      else if ( SCIPisFeasGT(scip, solval, SCIPvarGetUbGlobal(var)) )
      {
         /* set the solution value to the global upper bound */
         SCIP_CALL( SCIPsetSolVal(scip, worksol, var, ub) );
         ++nchgbound;
         SCIPdebugMessage("var <%s> type %d, old solval %g now fixed to ub %g\n", SCIPvarGetName(var), SCIPvarGetType(var), solval, ub);
      }
   }

   SCIPdebugMessage("number of bound changes (due to global bounds) = %d\n", nchgbound);

   SCIP_CALL( SCIPallocBufferArray(scip, &candvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candvals, nvars) );

   /* select candidate integer variables and determine rounding locks */
   ncands = 0;
   for (i = 0; i < nintvars; i++)
   {
      SCIP_VAR* var;

      var = vars[i];

      /* if it is not a slack or indicator variable */
      assert( 0 <= SCIPvarGetProbindex(var) && SCIPvarGetProbindex(var) < heurdata->nvars );
      if ( ! heurdata->isindvar[SCIPvarGetProbindex(var)] )
      {
         SCIP_Real objsum;
         int j;

         /* collect sum of objective values of indicator constraints containing variable */
         objsum = 0.0;
         assert( 0 <= heurdata->indbeg[i] && heurdata->indbeg[i] <= heurdata->nnonz );
         assert( heurdata->indbeg[i] <= heurdata->indbeg[i+1] && heurdata->indbeg[i+1] <= heurdata->nnonz );
         for (j = heurdata->indbeg[i]; j < heurdata->indbeg[i+1]; ++j)
         {
            SCIP_VAR* binvar;
            SCIP_CONS* cons;

            cons = heurdata->indcons[j];
            assert( cons != NULL );
            if ( SCIPconsIsActive(cons) )
            {
               binvar = SCIPgetBinaryVarIndicator(cons);
               if ( SCIPvarIsNegated(binvar) )
                  binvar = SCIPvarGetNegatedVar(binvar);
               objsum += SCIPvarGetObj(binvar);
            }
         }

         candvars[ncands] = var;
         /* @todo: use weighted locks - take objective of indicator constraints into account */
         candvals[ncands++] = (SCIP_Real) (SCIPvarGetNLocksDown(var) + SCIPvarGetNLocksUp(var)) + objsum;
      }
   }

   /* in degenerate cases we have only slack variables */
   if ( ncands > 0 )
   {
      SCIP_Real bestobj;

      /* sort arrays with respect to the values */
      if ( heurdata->largestweight )
      {
         SCIPsortDownRealPtr(candvals, (void**)candvars, ncands);
      }
      else
      {
         SCIPsortRealPtr(candvals, (void**)candvars, ncands);
      }

      bestobj = SCIPgetSolOrigObj(scip, worksol);
      SCIPdebugMessage("Starting with solution of value %g and %d candidate variables.\n", bestobj, ncands);

      /* loop through integer variables */
      for (i = 0; i < ncands && i < MAXCANDS; ++i)
      {
         SCIP_Real newval1;
         SCIP_VAR* var1;
         SCIP_Real val1;
         SCIP_Real lb1;
         SCIP_Real ub1;

         var1 = candvars[i];
         assert( SCIPvarGetType(var1) != SCIP_VARTYPE_CONTINUOUS );

         val1 = SCIPgetSolVal(scip, worksol, var1);
         lb1 = SCIPvarGetLbGlobal(var1);
         ub1 = SCIPvarGetUbGlobal(var1);
         assert( SCIPisIntegral(scip, lb1) );
         assert( SCIPisIntegral(scip, ub1) );

         /* check if we can shift down variable 1 */
         if ( SCIPisGT(scip, val1, lb1) )
         {
            if ( SCIPisIntegral(scip, val1) )
               newval1 = val1 - 1.0;
            else
               newval1 = SCIPfloor(scip, val1);

            /* loop through integer variables */
            for (k = i+1; k < ncands && k < MAXCANDS; ++k)
            {
               SCIP_Real newval2;
               SCIP_VAR* var2;
               SCIP_Real val2;
               SCIP_Real lb2;
               SCIP_Real ub2;

               var2 = candvars[k];
               assert( SCIPvarGetType(var2) != SCIP_VARTYPE_CONTINUOUS );
               assert( var1 != var2 );

               val2 = SCIPgetSolVal(scip, worksol, var2);
               lb2 = SCIPvarGetLbGlobal(var2);
               ub2 = SCIPvarGetUbGlobal(var2);
               assert( SCIPisIntegral(scip, lb2) );
               assert( SCIPisIntegral(scip, ub2) );

               /* check if we can shift down variable 2 */
               if ( SCIPisGT(scip, val2, lb2) )
               {
                  if ( SCIPisIntegral(scip, val2) )
                     newval2 = val2 - 1.0;
                  else
                     newval2 = SCIPfloor(scip, val2);

                  /* try whether setting variables to new values improves the objective */
                  SCIP_CALL( tryValue(scip, heurdata, bestobj, worksol, var1, val1, newval1, var2, val2, newval2, &success) );

                  if ( success )
                  {
                     SCIPdebugMessage("Changed <%s> [%d <= %f <= %d] to %d and <%s> [%d <= %f <= %d] to %d (solvalue: %f).\n",
                        SCIPvarGetName(var1), (int) lb1, val1, (int) ub1, (int) newval1,
                        SCIPvarGetName(var2), (int) lb2, val2, (int) ub2, (int) newval2,
                        SCIPgetSolOrigObj(scip, worksol));

                     bestobj = SCIPgetSolOrigObj(scip, worksol);
                     *result = SCIP_FOUNDSOL;
                     val1 = newval1;
                     val2 = newval2;
                     assert( SCIPisEQ(scip, SCIPgetSolVal(scip, worksol, var1), val1) );
                     assert( SCIPisEQ(scip, SCIPgetSolVal(scip, worksol, var2), val2) );
                  }
               }

               /* check if we can shift up variable 1 */
               if ( SCIPisLT(scip, val2, ub2) )
               {
                  if ( SCIPisIntegral(scip, val2) )
                     newval2 = val2 + 1.0;
                  else
                     newval2 = SCIPceil(scip, val2);

                  /* try whether setting variables to new values improves the objective */
                  SCIP_CALL( tryValue(scip, heurdata, bestobj, worksol, var1, val1, newval1, var2, val2, newval2, &success) );

                  if ( success )
                  {
                     SCIPdebugMessage("Changed <%s> [%d <= %f <= %d] to %d and <%s> [%d <= %f <= %d] to %d (solvalue: %f).\n",
                        SCIPvarGetName(var1), (int) lb1, val1, (int) ub1, (int) newval1,
                        SCIPvarGetName(var2), (int) lb2, val2, (int) ub2, (int) newval2,
                        SCIPgetSolOrigObj(scip, worksol));

                     bestobj = SCIPgetSolOrigObj(scip, worksol);
                     *result = SCIP_FOUNDSOL;
                     val1 = newval1;
                     val2 = newval2;
                     assert( SCIPisEQ(scip, SCIPgetSolVal(scip, worksol, var1), val1) );
                     assert( SCIPisEQ(scip, SCIPgetSolVal(scip, worksol, var2), val2) );
                  }
               }
#ifdef SCIP_DEBUG
               if ( *result == SCIP_FOUNDSOL )
               {
                  SCIP_Bool feasible;
                  SCIP_CALL( SCIPcheckSol(scip, worksol, TRUE, TRUE, TRUE, TRUE, &feasible) );
                  assert( feasible );
               }
#endif
            }
         }

         /* check if we can shift up variable 1 */
         if ( SCIPisLT(scip, val1, ub1) )
         {
            if ( SCIPisIntegral(scip, val1) )
               newval1 = val1 + 1.0;
            else
               newval1 = SCIPceil(scip, val1);

            /* loop through integer variables */
            for (k = i+1; k < ncands && k < MAXCANDS; k++)
            {
               SCIP_Real newval2;
               SCIP_VAR* var2;
               SCIP_Real val2;
               SCIP_Real lb2;
               SCIP_Real ub2;

               var2 = candvars[k];
               assert( SCIPvarGetType(var2) != SCIP_VARTYPE_CONTINUOUS );
               assert( var1 != var2 );

               val2 = SCIPgetSolVal(scip, worksol, var2);
               lb2 = SCIPvarGetLbGlobal(var2);
               ub2 = SCIPvarGetUbGlobal(var2);
               assert( SCIPisIntegral(scip, lb2) );
               assert( SCIPisIntegral(scip, ub2) );

               /* check if we can shift down variable 2 */
               if ( SCIPisGT(scip, val2, lb2) )
               {
                  if ( SCIPisIntegral(scip, val2) )
                     newval2 = val2 - 1.0;
                  else
                     newval2 = SCIPfloor(scip, val2);

                  /* try whether setting variables to new values improves the objective */
                  SCIP_CALL( tryValue(scip, heurdata, bestobj, worksol, var1, val1, newval1, var2, val2, newval2, &success) );

                  if ( success )
                  {
                     SCIPdebugMessage("Changed <%s> [%d <= %f <= %d] to %d and <%s> [%d <= %f <= %d] to %d (solvalue: %f).\n",
                        SCIPvarGetName(var1), (int) lb1, val1, (int) ub1, (int) newval1,
                        SCIPvarGetName(var2), (int) lb2, val2, (int) ub2, (int) newval2,
                        SCIPgetSolOrigObj(scip, worksol));

                     bestobj = SCIPgetSolOrigObj(scip, worksol);
                     *result = SCIP_FOUNDSOL;
                     val1 = newval1;
                     val2 = newval2;
                     assert( SCIPisEQ(scip, SCIPgetSolVal(scip, worksol, var1), val1) );
                     assert( SCIPisEQ(scip, SCIPgetSolVal(scip, worksol, var2), val2) );
                  }
               }

               /* check if we can shift up variable 1 */
               if ( SCIPisLT(scip, val2, ub2) )
               {
                  if ( SCIPisIntegral(scip, val2) )
                     newval2 = val2 + 1.0;
                  else
                     newval2 = SCIPceil(scip, val2);

                  /* try whether setting variables to new values improves the objective */
                  SCIP_CALL( tryValue(scip, heurdata, bestobj, worksol, var1, val1, newval1, var2, val2, newval2, &success) );

                  if ( success )
                  {
                     SCIPdebugMessage("Changed <%s> [%d <= %f <= %d] to %d and <%s> [%d <= %f <= %d] to %d (solvalue: %f).\n",
                        SCIPvarGetName(var1), (int) lb1, val1, (int) ub1, (int) newval1,
                        SCIPvarGetName(var2), (int) lb2, val2, (int) ub2, (int) newval2,
                        SCIPgetSolOrigObj(scip, worksol));

                     bestobj = SCIPgetSolOrigObj(scip, worksol);
                     *result = SCIP_FOUNDSOL;
                     val1 = newval1;
                     val2 = newval2;
                     assert( SCIPisEQ(scip, SCIPgetSolVal(scip, worksol, var1), val1) );
                     assert( SCIPisEQ(scip, SCIPgetSolVal(scip, worksol, var2), val2) );
                  }
               }
#ifdef SCIP_DEBUG
               if ( *result == SCIP_FOUNDSOL )
               {
                  SCIP_Bool feasible;
                  SCIP_CALL( SCIPcheckSol(scip, worksol, TRUE, TRUE, TRUE, TRUE, &feasible) );
                  assert( feasible );
               }
#endif
            }
         }
      }
   }
   SCIPdebugMessage("Finished 2-opt heuristic.\n");

   if ( *result == SCIP_FOUNDSOL )
   {
      /* SCIPdebug( SCIPprintSol(scip, worksol, NULL, FALSE) ); */
      SCIP_CALL( SCIPtrySol(scip, worksol, TRUE, TRUE, TRUE, TRUE, TRUE, &success) );
   }

   SCIPfreeBufferArray(scip, &candvals);
   SCIPfreeBufferArray(scip, &candvars);

   SCIP_CALL( SCIPfreeSol(scip, &worksol) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the twoopt primal heuristic for indicators and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurIndtwoopt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heurdata->lastsolindex = -1;
   heurdata->triedstart = FALSE;
   heurdata->nindconss = 0;
   heurdata->nvars = 0;
   heurdata->nnonz = 0;
   heurdata->maxnnonz = 0;
   heurdata->indcons = NULL;
   heurdata->indbeg = NULL;
   heurdata->isindvar = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyIndtwoopt, heurFreeIndtwoopt, heurInitIndtwoopt, heurExitIndtwoopt,
         heurInitsolIndtwoopt, heurExitsolIndtwoopt, heurExecIndtwoopt,
         heurdata) );

   /* add indtwoopt primal heuristic parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/indtwoopt/largestweight",
         "should we start with the variables of largest weight?",
         &heurdata->largestweight, TRUE, DEFAULT_LARGESTWEIGHT, NULL, NULL) );

   return SCIP_OKAY;
}
