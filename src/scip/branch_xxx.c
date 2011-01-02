/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch_xxx.c,v 1.28 2011/01/02 11:10:51 bzfheinz Exp $"

/**@file   branch_xxx.c
 * @ingroup BRANCHINGRULES
 * @brief  xxx branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/branch_xxx.h"


#define BRANCHRULE_NAME            "xxx"
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
SCIP_DECL_BRANCHCOPY(branchCopyXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchCopyXxx NULL
#endif

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_BRANCHFREE(branchFreeXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchFreeXxx NULL
#endif


/** initialization method of branching rule (called after problem was transformed) */
#if 0
static
SCIP_DECL_BRANCHINIT(branchInitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitXxx NULL
#endif


/** deinitialization method of branching rule (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BRANCHEXIT(branchExitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitXxx NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BRANCHINITSOL(branchInitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolXxx NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolXxx NULL
#endif


/** branching execution method for fractional LP solutions */
#if 0
static
SCIP_DECL_BRANCHEXECLP(branchExeclpXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExeclpXxx NULL
#endif


/** branching execution method for external candidates */
#if 0
static
SCIP_DECL_BRANCHEXECEXT(branchExecextXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecextXxx NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
#if 0
static
SCIP_DECL_BRANCHEXECPS(branchExecpsXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecpsXxx NULL
#endif





/*
 * branching rule specific interface methods
 */

/** creates the xxx branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleXxx(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create xxx branching rule data */
   branchruledata = NULL;
   /* TODO: (optional) create branching rule specific data here */

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH, 
	 BRANCHRULE_MAXBOUNDDIST,
         branchCopyXxx,
         branchFreeXxx, branchInitXxx, branchExitXxx,
         branchInitsolXxx, branchExitsolXxx,
         branchExeclpXxx, branchExecextXxx, branchExecpsXxx,
         branchruledata) );

   /* add xxx branching rule parameters */
   /* TODO: (optional) add branching rule specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
