/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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

/**@file   nullsubproblem.c
 * @brief  unit test for supplying NULL subproblem to SCIPaddBendersSubproblem
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scip/benders.h"
#include "scip/pub_benders.h"
#include "scip/scip.h"


#define SUBPROBOBJ   404040.40

/* UNIT TEST BENDERS */

#define BENDERS_NAME                "test"
#define BENDERS_DESC                "Benders' decomposition template"
#define BENDERS_PRIORITY            0
#define BENDERS_CUTLP            TRUE   /**< should Benders' cut be generated for LP solutions */
#define BENDERS_CUTPSEUDO        TRUE   /**< should Benders' cut be generated for pseudo solutions */
#define BENDERS_CUTRELAX         TRUE   /**< should Benders' cut be generated for relaxation solutions */
#define BENDERS_SHAREAUXVARS    FALSE   /**< should this Benders' share the highest priority Benders' aux vars */

/*
 * Data structures
 */

/* TODO: fill in the necessary Benders' decomposition data */

/** Benders' decomposition data */
struct SCIP_BendersData
{
   int                   nsubproblems;       /**< the number of subproblems in the Benders' decomposition */
};

/*
 * Local methods
 */

/*
 * Callback methods for Benders' decomposition
 */

/** copy method for benders plugins (called when SCIP copies plugins) */
#define bendersCopyTest NULL

/** destructor of Benders' decomposition to free user data (called when SCIP is exiting) */
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


/** initialization method of Benders' decomposition (called after problem was transformed) */
#define bendersInitTest NULL


/** deinitialization method of Benders' decomposition (called before transformed problem is freed) */
#define bendersExitTest NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin)
 *
 *  This function is called immediately after the auxiliary variables are created in the master problem. The callback
 *  provides the user an opportunity to add variable data to the auxiliary variables.
 */
#define bendersInitpreTest NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define bendersExitpreTest NULL


/** solving process initialization method of Benders' decomposition (called when branch and bound process is about to begin) */
#define bendersInitsolTest NULL


/** solving process deinitialization method of Benders' decomposition (called before branch and bound process data is freed) */
#define bendersExitsolTest NULL


/** mapping method between the master problem variables and the subproblem variables of Benders' decomposition */
static
SCIP_DECL_BENDERSGETVAR(bendersGetvarTest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of test Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** the method for creating the Benders' decomposition subproblem. This method is called during the initialisation stage
 *  (after the master problem was transformed)
 *
 *  This method must create the SCIP instance for the subproblem and add the required variables and constraints. In
 *  addition, the settings required for the solving the problem must be set here. However, some settings will be
 *  overridden by the standard solving method included in the Benders' decomposition framework. If a special solving
 *  method is desired, the user can implement the bendersSolvesubTest callback.
 */
static
SCIP_DECL_BENDERSCREATESUB(bendersCreatesubTest)
{  /*lint --e{715}*/

   /* adding NULL subproblems */
   SCIP_CALL( SCIPaddBendersSubproblem(scip, benders, NULL) );

   /* specifying the subproblem type */
   SCIPbendersSetSubproblemType(benders, probnumber, SCIP_BENDERSSUBTYPE_CONVEXCONT);

   return SCIP_OKAY;
}

/** called before the subproblem solve for Benders' decomposition */
#define bendersPresubsolveTest NULL

/** the solving method for a convex subproblem for Benders' decomposition. In this method the subproblem is setup with
 *  the given solution and then solved.
 *  NOTE: if either bendersSolvesubconvexTest or bendersSolvesubTest callbacks are implemented then the bendersFreesubTest
 *  callback must be implemented
 */
static
SCIP_DECL_BENDERSSOLVESUBCONVEX(bendersSolvesubconvexTest)
{  /*lint --e{715}*/
   (*objective) = SUBPROBOBJ;
   (*result) = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** the subproblem solving method for Benders' decomposition. In this method the subproblem is setup with the given
 *  solution and then solved.
 *  NOTE: if either bendersSolvesubconvexTest or bendersSolvesubTest callbacks are implemented then the bendersFreesubTest
 *  callback must be implemented
 */
static
SCIP_DECL_BENDERSSOLVESUB(bendersSolvesubTest)
{  /*lint --e{715}*/
   (*objective) = SUBPROBOBJ;
   (*result) = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

#define bendersPostsolveTest NULL


/** the subproblem freeing method for Benders' decomposition. This is called between subproblem solves to clear the
 *  solving data. Generally this will only require a call to SCIPfreeTransform. However, depending on the problem it
 *  could additional freeing methods.
 *  NOTE: the bendersFreesubTest callback must be implemented if the bendersSolvesubTest is implemented */
static
SCIP_DECL_BENDERSFREESUB(bendersFreesubTest)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}




/*
 * Benders' decomposition specific interface methods
 */

/** Creates a test Benders' decomposition algorithm and activates it in SCIP */
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

/** creates the test Benders' decomposition and includes it in SCIP */
static
SCIP_RETCODE SCIPincludeBendersTest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BENDERSDATA* bendersdata;
   SCIP_BENDERS* benders;

   /* create test Benders' decomposition data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &bendersdata) );
   bendersdata->nsubproblems = 0;

   benders = NULL;

   /* include Benders' decomposition */
   SCIP_CALL( SCIPincludeBendersBasic(scip, &benders, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, BENDERS_CUTLP,
         BENDERS_CUTPSEUDO, BENDERS_CUTRELAX, BENDERS_SHAREAUXVARS, bendersGetvarTest, bendersCreatesubTest,
         bendersdata) );
   assert(benders != NULL);

   /* set non fundamental callbacks via setter functions */
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
static int nsubproblems = 2;


#define EPS    1e-5

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* including the test Benders' decomposition */
   SCIP_CALL( SCIPincludeBendersTest(scip) );

   benders = SCIPfindBenders(scip, BENDERS_NAME);

   /* creating the problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test") );

   /* creating the test Benders' decomposition */
   SCIP_CALL( SCIPcreateBendersTest(scip, nsubproblems) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );
}


Test(bd, SCIPbendersSetSubproblemObjval, .init = setup, .fini = teardown,
   .description = "check SCIPbendersSetSubproblemObjval() subroutine of the Benders' decomposition core"
   )
{
   int i;

   for( i = 0; i < nsubproblems; i++ )
   {
      SCIPbendersSetSubproblemObjval(benders, i, SCIPinfinity(scip));
      cr_assert( SCIPisGE(scip, SCIPbendersGetSubproblemObjval(benders, i), SUBPROBOBJ) );
   }
}

Test(bd, SCIPsolveBendersSubproblem, .init = setup, .fini = teardown,
   .description = "check SCIPsolveBendersSubproblem() subroutine of the Benders' decomposition core"
   )
{
   int i;

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* calling the solve method for each subproblem */
   for( i = 0; i < nsubproblems; i++ )
   {
      SCIP_Bool infeasible;
      SCIP_Real objective;

      infeasible = FALSE;
      objective = SCIPinfinity(scip);

      /* informing the Benders' decomposition core that the subproblem is setup */
      SCIPbendersSubproblemIsSetup(benders, i);

      SCIP_CALL( SCIPsolveBendersSubproblem(scip, benders, NULL, i, &infeasible, FALSE, &objective) );

      cr_assert( !infeasible );
      cr_assert( SCIPisEQ(scip, objective, SUBPROBOBJ) );
   }
}

Test(bd, SCIPcomputeBendersSubproblemLowerbound, .init = setup, .fini = teardown,
   .description = "check SCIPcomputeBendersSubproblemLowerbound() subroutine of the Benders' decomposition core"
   )
{
   int i;

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* calling the solve method for each subproblem */
   for( i = 0; i < nsubproblems; i++ )
   {
      SCIP_Bool infeasible;
      SCIP_Real lowerbound;

      infeasible = FALSE;
      lowerbound = SCIPinfinity(scip);

      SCIP_CALL( SCIPcomputeBendersSubproblemLowerbound(scip, benders, i, &lowerbound, &infeasible) );

      cr_assert( !infeasible );
      cr_assert( SCIPisLE(scip, lowerbound, -SCIPinfinity(scip)) );
   }
}

Test(bd, SCIPbendersSubproblem, .init = setup, .fini = teardown,
   .description = "check SCIPbendersSubproblem() subroutine of the Benders' decomposition core"
   )
{
   int i;

   /* calling the solve method for each subproblem */
   for( i = 0; i < nsubproblems; i++ )
   {
      cr_assert( SCIPbendersSubproblem(benders, i) == NULL );
   }
}
