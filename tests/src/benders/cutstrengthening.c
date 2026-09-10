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

/**@file   cutstrengthening.c
 * @brief  unit test confirming that SCIPextractBendersIISMasterSolution finds an IIS of the master solution
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scip_benders.h"
#include "scip/pub_benders.h"
#include "scip/cons_linear.h"

/* struct_benders.h gives access to benders->dfbsdata (the DFBS gating fields aren't exposed via any public
 * getter/setter or SCIP parameter), and struct_scip.h/struct_stat.h are needed for the stat->status workaround
 * in setup() below.
 */
#include "scip/scip_message.h"
#include "scip/struct_benders.h"
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/type_benders.h"

#define NMASTERVARS   12

/* UNIT TEST BENDERS */

#define BENDERS_NAME                "cutstrengthtest"
#define BENDERS_DESC                "Benders' decomposition for testing non-convex cut strengthening"
#define BENDERS_PRIORITY            0
#define BENDERS_CUTLP            TRUE   /**< should Benders' cut be generated for LP solutions */
#define BENDERS_CUTPSEUDO        TRUE   /**< should Benders' cut be generated for pseudo solutions */
#define BENDERS_CUTRELAX         TRUE   /**< should Benders' cut be generated for relaxation solutions */
#define BENDERS_SHAREAUXVARS    FALSE   /**< should this Benders' share the highest priority Benders' aux vars */

/*
 * Data structures
 */

/** Benders' decomposition data.
 *
 * iisvars/niisvars describe the "true IIS" for the current test scenario: the subset of master
 * variable indices (into SCIPgetVars(scip)) that the mock subproblem solve considers necessary and
 * sufficient to reproduce the target result. infeasmode selects whether that target result is
 * infeasibility (mirroring SCIP_BENDERSSUBSTATUS_INFEAS) or a specific objective value (mirroring
 * SCIP_BENDERSSUBSTATUS_AUXVIOL).
 */
struct SCIP_BendersData
{
   int                   nsubproblems;       /**< the number of subproblems in the Benders' decomposition */
   SCIP_SOL*             origsol;            /**< the original solution */
   int*                  iisvars;            /**< indices of the variables forming the true IIS, not owned */
   int                   niisvars;           /**< the number of variables in the true IIS */
   SCIP_Bool             infeasmode;         /**< TRUE: oracle reports infeasibility; FALSE: oracle reports targetobj */
   SCIP_Real             targetobj;          /**< the objective value to report when the IIS is present and !infeasmode */
   SCIP_Bool             firstsolve;         /**< flag to indicate whether the first subproblem solve has occurred */
};

/*
 * Diagnostic printing helpers
 */

/** prints the names of the master variables (of the first NMASTERVARS variables of scip) that are set to 1.0
 * in sol, prefixed with the given label
 */
static
void printActiveVars(
   SCIP*                 scip,
   const char*           label,
   SCIP_SOL*             sol
   )
{
   SCIP_VAR** vars;
   int i;

   vars = SCIPgetVars(scip);

   SCIPinfoMessage(scip, NULL, "%s:", label);
   for( i = 0; i < NMASTERVARS; i++ )
   {
      if( SCIPisEQ(scip, SCIPgetSolVal(scip, sol, vars[i]), 1.0) )
         SCIPinfoMessage(scip, NULL, " %s", SCIPvarGetName(vars[i]));
   }
   SCIPinfoMessage(scip, NULL, "\n");
}

/** prints the master variable indices in idx, prefixed with the given label */
static
void printIndices(
   SCIP*                 scip,
   const char*           label,
   int*                  idx,
   int                   n
   )
{
   int i;

   SCIPinfoMessage(scip, NULL, "%s:", label);
   for( i = 0; i < n; i++ )
      SCIPinfoMessage(scip, NULL, " x%d", idx[i]);
   SCIPinfoMessage(scip, NULL, "\n");
}

/*
 * Callback methods for Benders' decomposition
 */

#define bendersCopyTest NULL

static
SCIP_DECL_BENDERSFREE(bendersFreeTest)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);

   SCIPfreeBlockMemory(scip, &bendersdata);

   return SCIP_OKAY;
}

#define bendersInitTest NULL
#define bendersExitTest NULL
#define bendersInitpreTest NULL
#define bendersExitpreTest NULL
#define bendersInitsolTest NULL
#define bendersExitsolTest NULL

/** mapping method between the master problem variables and the subproblem variables of Benders' decomposition.
 * Never invoked in these tests, since the subproblem is registered as NULL and generateBendersCuts is never called.
 */
static
SCIP_DECL_BENDERSGETVAR(bendersGetvarTest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cutstrengthtest Benders' decomposition not implemented\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** registers a NULL subproblem, since the mock solve callbacks below never build/solve a real SCIP subproblem */
static
SCIP_DECL_BENDERSCREATESUB(bendersCreatesubTest)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPaddBendersSubproblem(scip, benders, NULL) );

   SCIPbendersSetSubproblemType(benders, probnumber, SCIP_BENDERSSUBTYPE_NONCONVEXDIS);

   return SCIP_OKAY;
}

#define bendersPresubsolveTest NULL

/** the convex-relaxation solve callback: a constant pass-through. It always reports a low, finite
 * objective well below any target used in the tests, so that the internal cut-strengthening solve loop
 * never short-circuits based on the convex solve; the oracle logic in bendersSolvesubTest below always
 * decides the outcome.
 */
static
SCIP_DECL_BENDERSSOLVESUBCONVEX(bendersSolvesubconvexTest)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);

   /* printing the original solution and the expected IIS at the first call of this method. */
   if( !bendersdata->firstsolve )
   {
      SCIPinfoMessage(scip, NULL, "\n");
      printIndices(scip, "expected IIS", bendersdata->iisvars, bendersdata->niisvars);
      printActiveVars(scip, "original solution", sol);
      bendersdata->firstsolve = TRUE;
   }

   (*objective) = 0.0;
   (*result) = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** the IIS oracle: reports the configured target result (infeasible, or a target objective) if and only if
 * every variable in the configured "true IIS" has solution value 1.0 in sol; otherwise reports a
 * feasible result with an objective that cannot be confused with the target.
 */
static
SCIP_DECL_BENDERSSOLVESUB(bendersSolvesubTest)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;
   SCIP_VAR** vars;
   SCIP_Bool iiscomplete;
   int i;

   bendersdata = SCIPbendersGetData(benders);
   vars = SCIPgetVars(scip);

   iiscomplete = TRUE;
   for( i = 0; i < bendersdata->niisvars; i++ )
   {
      if( !SCIPisEQ(scip, SCIPgetSolVal(scip, sol, vars[bendersdata->iisvars[i]]), 1.0) )
      {
         iiscomplete = FALSE;
         break;
      }
   }

   printActiveVars(scip, "  [DFBS oracle] probing", sol);
   SCIPinfoMessage(scip, NULL, "     -> %s\n", iiscomplete ? "target reproduced" : "target NOT reproduced");

   if( iiscomplete && bendersdata->infeasmode )
   {
      (*objective) = SCIPinfinity(scip);
      (*result) = SCIP_INFEASIBLE;
   }
   else if( iiscomplete && !bendersdata->infeasmode )
   {
      (*objective) = bendersdata->targetobj;
      (*result) = SCIP_FEASIBLE;
   }
   else
   {
      /* in the feasible case, if the iis is not completed, then we return feasible but with an incorrect objective */
      (*objective) = 0.0;
      (*result) = SCIP_FEASIBLE;
   }

   return SCIP_OKAY;
}

#define bendersPostsolveTest NULL

static
SCIP_DECL_BENDERSFREESUB(bendersFreesubTest)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/*
 * Benders' decomposition specific interface methods
 */

static
SCIP_RETCODE SCIPcreateBendersTest(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   )
{
   SCIP_BENDERS* benders;
   SCIP_BENDERSDATA* bendersdata;

   assert(scip != NULL);

   benders = SCIPfindBenders(scip, BENDERS_NAME);
   bendersdata = SCIPbendersGetData(benders);

   bendersdata->nsubproblems = nsubproblems;

   SCIP_CALL( SCIPactivateBenders(scip, benders, nsubproblems) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE SCIPincludeBendersTest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BENDERSDATA* bendersdata;
   SCIP_BENDERS* benders;

   SCIP_CALL( SCIPallocBlockMemory(scip, &bendersdata) );
   bendersdata->nsubproblems = 0;
   bendersdata->origsol = NULL;
   bendersdata->iisvars = NULL;
   bendersdata->niisvars = 0;
   bendersdata->infeasmode = TRUE;
   bendersdata->targetobj = 0.0;
   bendersdata->firstsolve = FALSE;

   benders = NULL;

   SCIP_CALL( SCIPincludeBendersBasic(scip, &benders, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, BENDERS_CUTLP,
         BENDERS_CUTPSEUDO, BENDERS_CUTRELAX, BENDERS_SHAREAUXVARS, bendersGetvarTest, bendersCreatesubTest,
         bendersdata) );
   assert(benders != NULL);

   SCIP_CALL( SCIPsetBendersCopy(scip, benders, bendersCopyTest) );
   SCIP_CALL( SCIPsetBendersFree(scip, benders, bendersFreeTest) );
   SCIP_CALL( SCIPsetBendersInit(scip, benders, bendersInitTest) );
   SCIP_CALL( SCIPsetBendersExit(scip, benders, bendersExitTest) );
   SCIP_CALL( SCIPsetBendersInitpre(scip, benders, bendersInitpreTest) );
   SCIP_CALL( SCIPsetBendersExitpre(scip, benders, bendersExitpreTest) );
   SCIP_CALL( SCIPsetBendersInitsol(scip, benders, bendersInitsolTest) );
   SCIP_CALL( SCIPsetBendersExitsol(scip, benders, bendersExitsolTest) );
   SCIP_CALL( SCIPsetBendersPresubsolve(scip, benders, bendersPresubsolveTest) );
   SCIP_CALL( SCIPsetBendersSolveAndFreesub(scip, benders, bendersSolvesubconvexTest, bendersSolvesubTest,
         bendersFreesubTest) );
   SCIP_CALL( SCIPsetBendersPostsolve(scip, benders, bendersPostsolveTest) );

   return SCIP_OKAY;
}
/* END UNIT TEST BENDERS */

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_BENDERS* benders;
static int nsubproblems = 1;

static
void setup(void)
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR* var;
   int i;

   SCIP_CALL( SCIPcreate(&scip) );

   SCIP_CALL( SCIPincludeConshdlrLinear(scip) );

   SCIP_CALL( SCIPincludeBendersTest(scip) );
   benders = SCIPfindBenders(scip, BENDERS_NAME);

   SCIP_CALL( SCIPcreateProbBasic(scip, "test") );

   /* presolving is disabled so that the master variables are not fixed/aggregated away, keeping their
    * indices in SCIPgetVars() stable and predictable across the transformation to the solving stage
    */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   for( i = 0; i < NMASTERVARS; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(scip, var) );
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   SCIP_CALL( SCIPcreateBendersTest(scip, nsubproblems) );

   /* pushing SCIP through presolving/transformation. This triggers SCIPbendersInitpre, which creates the
    * auxiliary variable for the (single) subproblem.
    */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* TESTscipSetStage reaches SCIP_STAGE_SOLVING by registering a heuristic that calls SCIPinterruptSolve()
    * on the first call, so that SCIPsolve() halts immediately. This leaves stat->status permanently at
    * SCIP_STATUS_USERINTERRUPT, which makes SCIPisStopped() return TRUE forever - the very first guard in
    * the internal cut strengthening routine. Reset it so the function under test can actually execute.
    */
   scip->stat->status = SCIP_STATUS_UNKNOWN;
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );
}

/** configures the mock subproblem's IIS oracle for the current test */
static
void setIIS(
   SCIP_SOL*             origsol,
   int*                  iisvars,
   int                   niisvars,
   SCIP_Bool             infeasmode,
   SCIP_Real             targetobj
   )
{
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   bendersdata->origsol = origsol;
   bendersdata->iisvars = iisvars;
   bendersdata->niisvars = niisvars;
   bendersdata->infeasmode = infeasmode;
   bendersdata->targetobj = targetobj;
}

/** resets the DFBS gating fields so that SCIPextractBendersIISMasterSolution always executes its search,
 * bypassing the offset/frequency gating that is not exposed as a public parameter.
 */
static
void resetDfbsGating(void)
{
   benders->dfbsdata->offset = 0;
   benders->dfbsdata->freq = 1;
   benders->dfbsdata->numcalls = 0;
   benders->dfbsdata->totalnodes = 0;
   benders->dfbsdata->maxnodes = 5000;
}

/** counts the number of variables set to 1.0 in sol among the first n variables of scip */
static
int countActive(
   SCIP_SOL*             sol,
   int                   n
   )
{
   SCIP_VAR** vars;
   int count;
   int i;

   vars = SCIPgetVars(scip);
   count = 0;
   for( i = 0; i < n; i++ )
   {
      if( SCIPisEQ(scip, SCIPgetSolVal(scip, sol, vars[i]), 1.0) )
         count++;
   }

   return count;
}

/** returns TRUE iff every index in idx has solution value 1.0 in sol */
static
SCIP_Bool allActive(
   SCIP_SOL*             sol,
   int*                  idx,
   int                   n
   )
{
   SCIP_VAR** vars;
   int i;

   vars = SCIPgetVars(scip);
   for( i = 0; i < n; i++ )
   {
      if( !SCIPisEQ(scip, SCIPgetSolVal(scip, sol, vars[idx[i]]), 1.0) )
         return FALSE;
   }

   return TRUE;
}

/*
 * IIS-correctness tests
 */

Test(cutstrengthening, singleVariableIIS, .init = setup, .fini = teardown,
   .description = "DFBS reduces an 8-of-12 active master solution down to a single-variable IIS (infeasible substatus)"
   )
{
   SCIP_SOL* sol;
   SCIP_SOL* newsol;
   int active[8] = {0, 1, 2, 3, 4, 5, 6, 7};
   int iis[1] = {3};
   int origactive;
   SCIP_Bool notnull;
   SCIP_Bool alliispresent = FALSE;
   int reducedactive = -1;
   int i;

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   for( i = 0; i < 8; i++ )
      SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetVars(scip)[active[i]], 1.0) );

   setIIS(sol, iis, 1, TRUE, 0.0);
   resetDfbsGating();

   newsol = NULL;
   SCIP_CALL( SCIPextractBendersIISMasterSolution(scip, benders, sol, &newsol, 0, SCIP_BENDERSENFOTYPE_CHECK,
         SCIP_BENDERSSUBSTATUS_INFEAS) );

   if( newsol != NULL )
      printActiveVars(scip, "extracted solution", newsol);
   else
      SCIPinfoMessage(scip, NULL, "extracted solution: (NULL)\n");
   SCIPinfoMessage(scip, NULL, "\n");

   origactive = countActive(sol, NMASTERVARS);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   notnull = (newsol != NULL);
   if( notnull )
   {
      alliispresent = allActive(newsol, iis, 1);
      reducedactive = countActive(newsol, NMASTERVARS);
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
   }

   cr_assert( notnull );
   cr_assert( alliispresent, "expected IIS member x%d to be active in the reduced solution", iis[0] );
   cr_assert( reducedactive < origactive );
   cr_assert_eq( reducedactive, 1,
      "expected DFBS to isolate the single-variable IIS exactly, found %d active variables", reducedactive );
}

Test(cutstrengthening, twoVariableIIS, .init = setup, .fini = teardown,
   .description = "DFBS reduces an 8-of-12 active master solution down to a two-variable IIS (infeasible substatus)"
   )
{
   SCIP_SOL* sol;
   SCIP_SOL* newsol;
   int active[8] = {0, 1, 2, 3, 4, 5, 6, 7};
   int iis[2] = {2, 7};
   int origactive;
   SCIP_Bool notnull;
   SCIP_Bool alliispresent = FALSE;
   int reducedactive = -1;
   int i;

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   for( i = 0; i < 8; i++ )
      SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetVars(scip)[active[i]], 1.0) );

   setIIS(sol, iis, 2, TRUE, 0.0);
   resetDfbsGating();

   newsol = NULL;
   SCIP_CALL( SCIPextractBendersIISMasterSolution(scip, benders, sol, &newsol, 0, SCIP_BENDERSENFOTYPE_CHECK,
         SCIP_BENDERSSUBSTATUS_INFEAS) );

   if( newsol != NULL )
      printActiveVars(scip, "extracted solution", newsol);
   else
      SCIPinfoMessage(scip, NULL, "extracted solution: (NULL)\n");
   SCIPinfoMessage(scip, NULL, "\n");

   origactive = countActive(sol, NMASTERVARS);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   notnull = (newsol != NULL);
   if( notnull )
   {
      alliispresent = allActive(newsol, iis, 2);
      reducedactive = countActive(newsol, NMASTERVARS);
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
   }

   cr_assert( notnull );
   cr_assert( alliispresent, "expected all IIS members to be active in the reduced solution" );
   cr_assert( reducedactive < origactive );
   cr_assert_eq( reducedactive, 2,
      "expected DFBS to isolate the 2-variable IIS exactly, found %d active variables", reducedactive );
}

Test(cutstrengthening, multiVariableIIS, .init = setup, .fini = teardown,
   .description = "DFBS reduces an 8-of-12 active master solution down to a 3-variable IIS (infeasible substatus)"
   )
{
   SCIP_SOL* sol;
   SCIP_SOL* newsol;
   int active[8] = {0, 1, 2, 3, 4, 5, 6, 7};
   int iis[3] = {1, 4, 6};
   int origactive;
   SCIP_Bool notnull;
   SCIP_Bool alliispresent = FALSE;
   int reducedactive = -1;
   int i;

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   for( i = 0; i < 8; i++ )
      SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetVars(scip)[active[i]], 1.0) );

   setIIS(sol, iis, 3, TRUE, 0.0);
   resetDfbsGating();

   newsol = NULL;
   SCIP_CALL( SCIPextractBendersIISMasterSolution(scip, benders, sol, &newsol, 0, SCIP_BENDERSENFOTYPE_CHECK,
         SCIP_BENDERSSUBSTATUS_INFEAS) );

   if( newsol != NULL )
      printActiveVars(scip, "extracted solution", newsol);
   else
      SCIPinfoMessage(scip, NULL, "extracted solution: (NULL)\n");
   SCIPinfoMessage(scip, NULL, "\n");

   origactive = countActive(sol, NMASTERVARS);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   notnull = (newsol != NULL);
   if( notnull )
   {
      alliispresent = allActive(newsol, iis, 3);
      reducedactive = countActive(newsol, NMASTERVARS);
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
   }

   cr_assert( notnull );
   cr_assert( alliispresent, "expected all IIS members to be active in the reduced solution" );
   cr_assert( reducedactive < origactive );
   cr_assert_eq( reducedactive, 3,
      "expected DFBS to isolate the 3-variable IIS exactly, found %d active variables", reducedactive );
}

Test(cutstrengthening, largeActiveSetFourVariableIIS, .init = setup, .fini = teardown,
   .description = "DFBS reduces an 11-of-12 active master solution towards a four-variable IIS (infeasible substatus)"
   )
{
   SCIP_SOL* sol;
   SCIP_SOL* newsol;
   int active[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
   int iis[4] = {2, 5, 8, 10};
   int origactive;
   SCIP_Bool notnull;
   SCIP_Bool alliispresent = FALSE;
   int reducedactive = -1;
   int i;

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   for( i = 0; i < 11; i++ )
      SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetVars(scip)[active[i]], 1.0) );

   setIIS(sol, iis, 4, TRUE, 0.0);
   resetDfbsGating();

   newsol = NULL;
   SCIP_CALL( SCIPextractBendersIISMasterSolution(scip, benders, sol, &newsol, 0, SCIP_BENDERSENFOTYPE_CHECK,
         SCIP_BENDERSSUBSTATUS_INFEAS) );

   if( newsol != NULL )
      printActiveVars(scip, "extracted solution", newsol);
   else
      SCIPinfoMessage(scip, NULL, "extracted solution: (NULL)\n");
   SCIPinfoMessage(scip, NULL, "\n");

   origactive = countActive(sol, NMASTERVARS);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   notnull = (newsol != NULL);
   if( notnull )
   {
      alliispresent = allActive(newsol, iis, 4);
      reducedactive = countActive(newsol, NMASTERVARS);
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
   }

   cr_assert( notnull );
   cr_assert( alliispresent, "expected all IIS members to be active in the reduced solution" );
   cr_assert( reducedactive < origactive,
      "expected some reduction from the original %d active variables, found %d", origactive, reducedactive );
   cr_assert_eq( reducedactive, 4,
      "expected DFBS to isolate the 4-variable IIS exactly, found %d active variables", reducedactive );
}

Test(cutstrengthening, wholeActiveSetIsIIS, .init = setup, .fini = teardown,
   .description = "DFBS cannot reduce the active set when every active variable is required"
   )
{
   SCIP_SOL* sol;
   SCIP_SOL* newsol;
   int active[8] = {0, 1, 2, 3, 4, 5, 6, 7};
   int prevfreq;
   SCIP_Bool notnull;
   int reducedactive = -1;
   int i;

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   for( i = 0; i < 8; i++ )
      SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetVars(scip)[active[i]], 1.0) );

   setIIS(sol, active, 8, TRUE, 0.0);
   resetDfbsGating();
   prevfreq = benders->dfbsdata->freq;

   newsol = NULL;
   SCIP_CALL( SCIPextractBendersIISMasterSolution(scip, benders, sol, &newsol, 0, SCIP_BENDERSENFOTYPE_CHECK,
         SCIP_BENDERSSUBSTATUS_INFEAS) );

   if( newsol != NULL )
      printActiveVars(scip, "extracted solution", newsol);
   else
      SCIPinfoMessage(scip, NULL, "extracted solution: (NULL)\n");
   SCIPinfoMessage(scip, NULL, "\n");

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   notnull = (newsol != NULL);
   if( notnull )
   {
      reducedactive = countActive(newsol, NMASTERVARS);
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
   }

   cr_assert( notnull );
   cr_assert_eq( reducedactive, 8,
      "expected no reduction since the whole active set is required, found %d active variables", reducedactive );
   cr_assert_eq( benders->dfbsdata->freq, prevfreq + 5,
      "expected the adaptive frequency to increase by 5 after a call with no reduction" );
}

Test(cutstrengthening, auxviolObjectiveIIS, .init = setup, .fini = teardown,
   .description = "DFBS reduces the active set based on reproducing the subproblem objective (AUXVIOL substatus)"
   )
{
   SCIP_SOL* sol;
   SCIP_SOL* newsol;
   SCIP_VAR* auxvar;
   int active[8] = {0, 1, 2, 3, 4, 5, 6, 7};
   int iis[3] = {2, 5, 7};
   SCIP_Real targetobj = 100.0;
   int origactive;
   SCIP_Bool notnull;
   SCIP_Bool alliispresent = FALSE;
   int reducedactive = -1;
   int i;

   auxvar = SCIPbendersGetAuxiliaryVar(benders, 0);
   cr_assert( auxvar != NULL );

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   for( i = 0; i < 8; i++ )
      SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetVars(scip)[active[i]], 1.0) );
   /* the auxiliary variable value must differ from the target objective by more than a relative 0.1 gap,
    * otherwise the internal cut strengthening routine exits immediately (see the AUXVIOL gap check in
    * src/scip/benders.c)
    */
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 0.0) );

   SCIPbendersSetSubproblemObjval(benders, 0, targetobj);
   setIIS(sol, iis, 3, FALSE, targetobj);
   resetDfbsGating();

   newsol = NULL;
   SCIP_CALL( SCIPextractBendersIISMasterSolution(scip, benders, sol, &newsol, 0, SCIP_BENDERSENFOTYPE_CHECK,
         SCIP_BENDERSSUBSTATUS_AUXVIOL) );

   if( newsol != NULL )
      printActiveVars(scip, "extracted solution", newsol);
   else
      SCIPinfoMessage(scip, NULL, "extracted solution: (NULL)\n");
   SCIPinfoMessage(scip, NULL, "\n");

   origactive = countActive(sol, NMASTERVARS);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   notnull = (newsol != NULL);
   if( notnull )
   {
      alliispresent = allActive(newsol, iis, 3);
      reducedactive = countActive(newsol, NMASTERVARS);
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
   }

   cr_assert( notnull );
   cr_assert( alliispresent, "expected all IIS members to be active in the reduced solution" );
   cr_assert( reducedactive < origactive );
   cr_assert_eq( reducedactive, 3,
      "expected DFBS to isolate the 3-variable IIS exactly, found %d active variables", reducedactive );
}

Test(cutstrengthening, auxviolSingleVariableIIS, .init = setup, .fini = teardown,
   .description = "DFBS reduces an 8-of-12 active master solution down to a single-variable IIS (AUXVIOL substatus)"
   )
{
   SCIP_SOL* sol;
   SCIP_SOL* newsol;
   SCIP_VAR* auxvar;
   int active[8] = {0, 1, 2, 3, 4, 5, 6, 7};
   int iis[1] = {4};
   SCIP_Real targetobj = 100.0;
   int origactive;
   SCIP_Bool notnull;
   SCIP_Bool alliispresent = FALSE;
   int reducedactive = -1;
   int i;

   auxvar = SCIPbendersGetAuxiliaryVar(benders, 0);
   cr_assert( auxvar != NULL );

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   for( i = 0; i < 8; i++ )
      SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetVars(scip)[active[i]], 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 0.0) );

   SCIPbendersSetSubproblemObjval(benders, 0, targetobj);
   setIIS(sol, iis, 1, FALSE, targetobj);
   resetDfbsGating();

   newsol = NULL;
   SCIP_CALL( SCIPextractBendersIISMasterSolution(scip, benders, sol, &newsol, 0, SCIP_BENDERSENFOTYPE_CHECK,
         SCIP_BENDERSSUBSTATUS_AUXVIOL) );

   if( newsol != NULL )
      printActiveVars(scip, "extracted solution", newsol);
   else
      SCIPinfoMessage(scip, NULL, "extracted solution: (NULL)\n");
   SCIPinfoMessage(scip, NULL, "\n");

   origactive = countActive(sol, NMASTERVARS);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   notnull = (newsol != NULL);
   if( notnull )
   {
      alliispresent = allActive(newsol, iis, 1);
      reducedactive = countActive(newsol, NMASTERVARS);
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
   }

   cr_assert( notnull );
   cr_assert( alliispresent, "expected IIS member x%d to be active in the reduced solution", iis[0] );
   cr_assert( reducedactive < origactive );
   cr_assert_eq( reducedactive, 1,
      "expected DFBS to isolate the single-variable IIS exactly, found %d active variables", reducedactive );
}

Test(cutstrengthening, auxviolWholeActiveSetIsIIS, .init = setup, .fini = teardown,
   .description = "DFBS cannot reduce the active set when every active variable is required (AUXVIOL substatus)"
   )
{
   SCIP_SOL* sol;
   SCIP_SOL* newsol;
   SCIP_VAR* auxvar;
   int active[8] = {0, 1, 2, 3, 4, 5, 6, 7};
   SCIP_Real targetobj = 100.0;
   int prevfreq;
   SCIP_Bool notnull;
   int reducedactive = -1;
   int i;

   auxvar = SCIPbendersGetAuxiliaryVar(benders, 0);
   cr_assert( auxvar != NULL );

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   for( i = 0; i < 8; i++ )
      SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetVars(scip)[active[i]], 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 0.0) );

   SCIPbendersSetSubproblemObjval(benders, 0, targetobj);
   setIIS(sol, active, 8, FALSE, targetobj);
   resetDfbsGating();
   prevfreq = benders->dfbsdata->freq;

   newsol = NULL;
   SCIP_CALL( SCIPextractBendersIISMasterSolution(scip, benders, sol, &newsol, 0, SCIP_BENDERSENFOTYPE_CHECK,
         SCIP_BENDERSSUBSTATUS_AUXVIOL) );

   if( newsol != NULL )
      printActiveVars(scip, "extracted solution", newsol);
   else
      SCIPinfoMessage(scip, NULL, "extracted solution: (NULL)\n");
   SCIPinfoMessage(scip, NULL, "\n");

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   notnull = (newsol != NULL);
   if( notnull )
   {
      reducedactive = countActive(newsol, NMASTERVARS);
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
   }

   cr_assert( notnull );
   cr_assert_eq( reducedactive, 8,
      "expected no reduction since the whole active set is required, found %d active variables", reducedactive );
   cr_assert_eq( benders->dfbsdata->freq, prevfreq + 5,
      "expected the adaptive frequency to increase by 5 after a call with no reduction" );
}
