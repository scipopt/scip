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

/**@file   relax_unittest.c
 * @brief  unittest relaxator
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "relax_unittest.h"


#define RELAX_NAME             "relax-unittest"
#define RELAX_DESC             "relaxator template"
#define RELAX_PRIORITY         101
#define RELAX_FREQ             2


/*
 * Data structures
 */

/* TODO: fill in the necessary relaxator data */

/** relaxator data */
struct SCIP_RelaxData
{
   int ncalls;
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of relaxator
 */

/* TODO: Implement all necessary relaxator methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for relaxator plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_RELAXCOPY(relaxCopyUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxCopyUnittest NULL
#endif

/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_RELAXFREE(relaxFreeUnittest)
{  /*lint --e{715}*/
   /* call destructor of relaxation handler */

   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);

   assert(relaxdata != NULL);

   SCIPfreeMemory(scip, &relaxdata);
   SCIPrelaxSetData(relax, NULL);

   return SCIP_OKAY;
}


/** initialization method of relaxator (called after problem was transformed) */
#if 0
static
SCIP_DECL_RELAXINIT(relaxInitUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxInitUnittest NULL
#endif


/** deinitialization method of relaxator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_RELAXEXIT(relaxExitUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxExitUnittest NULL
#endif


/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_RELAXINITSOL(relaxInitsolUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxInitsolUnittest NULL
#endif


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_RELAXEXITSOL(relaxExitsolUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxExitsolUnittest NULL
#endif


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecUnittest)
{
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);

   assert(relaxdata != NULL);

   relaxdata->ncalls++;

   return SCIP_OKAY;
}


/*
 * relaxator specific interface methods
 */

/** creates the unittest relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax;

   /* create unittest relaxator data */
   SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );
   relaxdata->ncalls = 0;

   relax = NULL;

   /* include relaxator */
#if 0
   /* use SCIPincludeRelax() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeRelax(scip, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxCopyUnittest, relaxFreeUnittest, relaxInitUnittest, relaxExitUnittest, relaxInitsolUnittest, relaxExitsolUnittest, relaxExecUnittest,
         relaxdata) );
#else
   /* use SCIPincludeRelaxBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecUnittest, relaxdata) );

   assert(relax != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetRelaxCopy(scip, relax, relaxCopyUnittest) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeUnittest) );
   SCIP_CALL( SCIPsetRelaxInit(scip, relax, relaxInitUnittest) );
   SCIP_CALL( SCIPsetRelaxExit(scip, relax, relaxExitUnittest) );
   SCIP_CALL( SCIPsetRelaxInitsol(scip, relax, relaxInitsolUnittest) );
   SCIP_CALL( SCIPsetRelaxExitsol(scip, relax, relaxExitsolUnittest) );
#endif

   /* add unittest relaxator parameters */
   /* TODO: (optional) add relaxator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** get the number of calls of the relaxator */
int SCIPgetNcallsUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = SCIPfindRelax(scip, RELAX_NAME);

   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);

   assert(relaxdata != NULL);

   return relaxdata->ncalls;
}
