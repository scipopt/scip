/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch_xxx.c,v 1.11 2005/02/14 13:35:39 bzfpfend Exp $"

/**@file   branch_xxx.c
 * @brief  xxx branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/branch_xxx.h"


#define BRANCHRULE_NAME            "xxx"
#define BRANCHRULE_DESC            "branching rule template"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH      -1




/*
 * Data structures
 */

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct BranchData
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

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#if 0
static
DECL_BRANCHFREE(branchFreeXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx branching rule not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchFreeXxx NULL
#endif


/** initialization method of branching rule (called after problem was transformed) */
#if 0
static
DECL_BRANCHINIT(branchInitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx branching rule not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitXxx NULL
#endif


/** deinitialization method of branching rule (called before transformed problem is freed) */
#if 0
static
DECL_BRANCHEXIT(branchExitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx branching rule not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitXxx NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#if 0
static
DECL_BRANCHINITSOL(branchInitsolXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx branching rule not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolXxx NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#if 0
static
DECL_BRANCHEXITSOL(branchExitsolXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx branching rule not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolXxx NULL
#endif


/** branching execution method for fractional LP solutions */
#if 0
static
DECL_BRANCHEXECLP(branchExeclpXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx branching rule not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExeclpXxx NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
#if 0
static
DECL_BRANCHEXECPS(branchExecpsXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx branching rule not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecpsXxx NULL
#endif





/*
 * branching rule specific interface methods
 */

/** creates the xxx branching rule and includes it in SCIP */
RETCODE SCIPincludeBranchruleXxx(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   BRANCHRULEDATA* branchruledata;

   /* create xxx branching rule data */
   branchruledata = NULL;
   /* TODO: (optional) create branching rule specific data here */

   /* include branching rule */
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
         branchFreeXxx, branchInitXxx, branchExitXxx,
         branchInitsolXxx, branchExitsolXxx,
         branchExeclpXxx, branchExecpsXxx,
         branchruledata) );

   /* add xxx branching rule parameters */
   /* TODO: (optional) add branching rule specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
