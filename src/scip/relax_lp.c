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

/**@file   relax_lp.c
 * @brief  lp relaxator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/relax_lp.h"


#define RELAX_NAME             "lp"
#define RELAX_DESC             "relaxator template"
#define RELAX_PRIORITY         0
#define RELAX_FREQ             1
#define RELAX_FULLLPINFO       TRUE




/*
 * Data structures
 */

/* TODO: fill in the necessary relaxator data */

/** relaxator data */
struct SCIP_RelaxData
{

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
SCIP_DECL_RELAXCOPY(relaxCopyLp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lp relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxCopyLp NULL
#endif

/** destructor of relaxator to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_RELAXFREE(relaxFreeLp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lp relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxFreeLp NULL
#endif


/** initialization method of relaxator (called after problem was transformed) */
#if 0
static
SCIP_DECL_RELAXINIT(relaxInitLp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lp relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxInitLp NULL
#endif


/** deinitialization method of relaxator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_RELAXEXIT(relaxExitLp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lp relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxExitLp NULL
#endif


/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_RELAXINITSOL(relaxInitsolLp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lp relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxInitsolLp NULL
#endif


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_RELAXEXITSOL(relaxExitsolLp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lp relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxExitsolLp NULL
#endif


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecLp)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}





/*
 * relaxator specific interface methods
 */

/** creates the lp relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxLp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax;

   /* create lp relaxator data */
   relaxdata = NULL;
   /* TODO: (optional) create relaxator specific data here */

   relax = NULL;

   /* include relaxator */
#if 0
   /* use SCIPincludeRelax() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeRelax(scip, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ, RELAX_FULLLPINFO,
         relaxCopyLp, relaxFreeLp, relaxInitLp, relaxExitLp, relaxInitsolLp, relaxExitsolLp, relaxExecLp,
         relaxdata) );
#else
   /* use SCIPincludeRelaxBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ, RELAX_FULLLPINFO,
         relaxExecLp, relaxdata) );

   assert(relax != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetRelaxCopy(scip, relax, relaxCopyLp) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeLp) );
   SCIP_CALL( SCIPsetRelaxInit(scip, relax, relaxInitLp) );
   SCIP_CALL( SCIPsetRelaxExit(scip, relax, relaxExitLp) );
   SCIP_CALL( SCIPsetRelaxInitsol(scip, relax, relaxInitsolLp) );
   SCIP_CALL( SCIPsetRelaxExitsol(scip, relax, relaxExitsolLp) );
#endif

   /* add lp relaxator parameters */
   /* TODO: (optional) add relaxator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
