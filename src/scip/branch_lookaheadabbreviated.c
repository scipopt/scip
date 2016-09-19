/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_LookaheadAbbreviated.c
 * @brief  LookaheadAbbreviated branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/branch_lookaheadabbreviated.h"


#define BRANCHRULE_NAME            "LookaheadAbbreviated"
#define BRANCHRULE_DESC            "branching rule template"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0


/*
 * Data structures
 */

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of branching rule
 */

/* TODO: Implement all necessary branching rule methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for branchrule plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BRANCHCOPY(branchCopyLookaheadAbbreviated)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of LookaheadAbbreviated branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchCopyLookaheadAbbreviated NULL
#endif

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_BRANCHFREE(branchFreeLookaheadAbbreviated)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of LookaheadAbbreviated branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchFreeLookaheadAbbreviated NULL
#endif


/** initialization method of branching rule (called after problem was transformed) */
#if 0
static
SCIP_DECL_BRANCHINIT(branchInitLookaheadAbbreviated)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of LookaheadAbbreviated branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitLookaheadAbbreviated NULL
#endif


/** deinitialization method of branching rule (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BRANCHEXIT(branchExitLookaheadAbbreviated)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of LookaheadAbbreviated branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitLookaheadAbbreviated NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BRANCHINITSOL(branchInitsolLookaheadAbbreviated)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of LookaheadAbbreviated branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolLookaheadAbbreviated NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolLookaheadAbbreviated)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of LookaheadAbbreviated branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolLookaheadAbbreviated NULL
#endif


/** branching execution method for fractional LP solutions */
#if 0
static
SCIP_DECL_BRANCHEXECLP(branchExeclpLookaheadAbbreviated)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of LookaheadAbbreviated branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExeclpLookaheadAbbreviated NULL
#endif


/** branching execution method for external candidates */
#if 0
static
SCIP_DECL_BRANCHEXECEXT(branchExecextLookaheadAbbreviated)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of LookaheadAbbreviated branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecextLookaheadAbbreviated NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
#if 0
static
SCIP_DECL_BRANCHEXECPS(branchExecpsLookaheadAbbreviated)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of LookaheadAbbreviated branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecpsLookaheadAbbreviated NULL
#endif


/*
 * branching rule specific interface methods
 */

/** creates the LookaheadAbbreviated branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleLookaheadAbbreviated(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create LookaheadAbbreviated branching rule data */
   branchruledata = NULL;
   /* TODO: (optional) create branching rule specific data here */

   branchrule = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyLookaheadAbbreviated) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeLookaheadAbbreviated) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitLookaheadAbbreviated) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitLookaheadAbbreviated) );
   SCIP_CALL( SCIPsetBranchruleInitsol(scip, branchrule, branchInitsolLookaheadAbbreviated) );
   SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolLookaheadAbbreviated) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpLookaheadAbbreviated) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextLookaheadAbbreviated) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsLookaheadAbbreviated) );

   /* add LookaheadAbbreviated branching rule parameters */
   /* TODO: (optional) add branching rule specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
