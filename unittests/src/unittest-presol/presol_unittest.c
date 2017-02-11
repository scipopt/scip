/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_unittest.c
 * @brief  unittest presolver
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "presol_unittest.h"


#define PRESOL_NAME            "unittest"
#define PRESOL_DESC            "presolver template"
#define PRESOL_PRIORITY        20010001 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS             20 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_ALWAYS /* timing of the presolver (fast, medium, or exhaustive) */


/*
 * Data structures
 */

/* TODO: fill in the necessary presolver data */

/** presolver data */
/*struct SCIP_PresolData
{
};
*/

/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of presolver
 */

/* TODO: Implement all necessary presolver methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PRESOLCOPY(presolCopyUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolCopyUnittest NULL
#endif


/** destructor of presolver to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_PRESOLFREE(presolFreeUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolFreeUnittest NULL
#endif


/** initialization method of presolver (called after problem was transformed) */
#if 0
static
SCIP_DECL_PRESOLINIT(presolInitUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitUnittest NULL
#endif


/** deinitialization method of presolver (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PRESOLEXIT(presolExitUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitUnittest NULL
#endif


/** presolving initialization method of presolver (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_PRESOLINITPRE(presolInitpreUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitpreUnittest NULL
#endif


/** presolving deinitialization method of presolver (called after presolving has been finished) */
#if 0
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitpreUnittest NULL
#endif


/** execution method of presolver */
/* TODO: complete */
static
SCIP_DECL_PRESOLEXEC(presolExecUnittest)
{
   SCIP_VAR** vars;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;

   assert(scip != NULL);

   vars = SCIPgetVars(scip);
#if 0

   int i;
   int nvars;
   SCIP_Bool tightened;
   nvars = SCIPgetNVars(scip);
   /* tighten lower bound of the first variable to 1.0 */
   SCIP_CALL( SCIPtightenVarLb( scip, vars[0], 1.0, TRUE, &infeasible, &tightened) );

   /* tighten upper bounds to 2.0 */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPtightenVarUb( scip, vars[i], 2.0, TRUE, &infeasible, &tightened) );
   }

   /* fix first variable to 1 */
   SCIP_CALL( SCIPfixVar(scip, vars[0], 1.0, &infeasible, &fixed) );
#endif


   SCIP_CALL( SCIPfixVar(scip, vars[0], 2.0, &infeasible, &fixed) );
   SCIPdebugMessage("FIXED: %d INFEASIBLE: %d\n" , fixed, infeasible);

   if( fixed )
   {
      (*nfixedvars)++;
   }


   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the unittest presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create unittest presolver data */
   presoldata = NULL;
   /* TODO: (optional) create presolver specific data here */

   presol = NULL;

   /* include presolver */
#if 0
   /* use SCIPincludePresol() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_DELAY,
         presolCopyUnittest, presolFreeUnittest, presolInitUnittest, presolExitUnittest, presolInitpreUnittest, presolExitpreUnittest, presolExecUnittest,
         presoldata) );
#else
   /* use SCIPincludePresolBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecUnittest,
         presoldata) );

   assert(presol != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyUnittest) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeUnittest) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitUnittest) );
   SCIP_CALL( SCIPsetPresolExit(scip, presol, presolExitUnittest) );
   SCIP_CALL( SCIPsetPresolInitpre(scip, presol, presolInitpreUnittest) );
   SCIP_CALL( SCIPsetPresolExitpre(scip, presol, presolExitpreUnittest) );
#endif

   /* add unittest presolver parameters */
   /* TODO: (optional) add presolver specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
