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

/**@file   benders_xyz.c
 * @brief  xyz Benders' decomposition algorithm
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/benders_xyz.h"
#include "scip/pub_benders.h"
#include "scip/bendersdefcuts.h"


#define BENDERS_NAME            "xyz"
#define BENDERS_DESC            "Benders' decomposition template"
#define BENDERS_PRIORITY        0
#define BENDERS_CUTLP        TRUE   /**< should Benders' cut be generated for LP solutions */
#define BENDERS_CUTPSEUDO    TRUE   /**< should Benders' cut be generated for pseudo solutions */
#define BENDERS_CUTRELAX     TRUE   /**< should Benders' cut be generated for relaxation solutions */




/*
 * Data structures
 */

/* TODO: fill in the necessary Benders' decomposition data */

/** Benders' decomposition data */
struct SCIP_BendersData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods for Benders' decomposition
 */

/* TODO: Implement all necessary Benders' decomposition methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for benders plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BENDERSCOPY(bendersCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decompostion not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersCopyXyz NULL
#endif

/** destructor of Benders' decomposition to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_BENDERSFREE(bendersFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersFreeXyz NULL
#endif


/** initialization method of Benders' decomposition (called after problem was transformed) */
#if 0
static
SCIP_DECL_BENDERSINIT(bendersInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitXyz NULL
#endif


/** deinitialization method of Benders' decomposition (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BENDERSEXIT(bendersExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitXyz NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin)
 *
 *  This function is called immediately after the auxiliary variables are created in the master problem. The callback
 *  provides the user an opportunity to add variable data to the auxiliary variables.
 */
#if 0
static
SCIP_DECL_BENDERSINITPRE(bendersInitpreXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitpreXyz NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_BENDERSEXITPRE(bendersExitpreXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitpreXyz NULL
#endif


/** solving process initialization method of Benders' decomposition (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BENDERSINITSOL(bendersInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitsolXyz NULL
#endif


/** solving process deinitialization method of Benders' decomposition (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BENDERSEXITSOL(bendersExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitsolXyz NULL
#endif


/** mapping method between the master problem variables and the subproblem variables of Benders' decomposition */
static
SCIP_DECL_BENDERSGETMASTERVAR(bendersGetmastervarXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** the execution method for Benders' decomposition */
static
SCIP_DECL_BENDERSEXEC(bendersExecXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** the subproblem solving method for Benders' decomposition */
static
SCIP_DECL_BENDERSSOLVESUB(bendersSolvesubXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

#if 0
/** the post-solve method for Benders' decomposition */
static
SCIP_DECL_BENDERSPOSTSOLVE(bendersPostsolveXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersPostsolveXyz NULL
#endif


/** the subproblem freeing method for Benders' decomposition */
static
SCIP_DECL_BENDERSFREESUB(bendersFreesubXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}




/*
 * Benders' decomposition specific interface methods
 */

/** creates the xyz Benders' decomposition and includes it in SCIP */
SCIP_RETCODE SCIPincludeBendersXyz(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   )
{
   SCIP_BENDERSDATA* bendersdata;
   SCIP_BENDERS* benders;

   /* create xyz Benders' decomposition data */
   bendersdata = NULL;
   /* TODO: (optional) create Benders' decomposition specific data here */

   benders = NULL;

   /* include Benders' decomposition */
#if 0
   /* use SCIPincludeBenders() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBenders(scip, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, nsubproblems,
         bendersCopyXyz, bendersFreeXyz, bendersInitXyz, bendersExitXyz, bendersInitpreXyz, bendersExitpreXyz,
         bendersInitsolXyz, bendersExitsolXyz, bendersGetmastervarXyz, bendersExecXyz,
         bendersPostsolveXyz, bendersFreesubXyz, bendersdata) );
#else
   /* use SCIPincludeBendersBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBendersBasic(scip, &benders, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, nsubproblems,
         BENDERS_CUTLP, BENDERS_CUTPSEUDO, BENDERS_CUTRELAX, bendersGetmastervarXyz, bendersExecXyz, bendersSolvesubXyz,
         bendersFreesubXyz, bendersdata) );
   assert(benders != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBendersCopy(scip, benders, bendersCopyXyz) );
   SCIP_CALL( SCIPsetBendersFree(scip, benders, bendersFreeXyz) );
   SCIP_CALL( SCIPsetBendersInit(scip, benders, bendersInitXyz) );
   SCIP_CALL( SCIPsetBendersExit(scip, benders, bendersExitXyz) );
   SCIP_CALL( SCIPsetBendersInitpre(scip, benders, bendersInitpreXyz) );
   SCIP_CALL( SCIPsetBendersExitpre(scip, benders, bendersExitpreXyz) );
   SCIP_CALL( SCIPsetBendersInitsol(scip, benders, bendersInitsolXyz) );
   SCIP_CALL( SCIPsetBendersExitsol(scip, benders, bendersExitsolXyz) );
   SCIP_CALL( SCIPsetBendersPostsolve(scip, benders, bendersPostsolveXyz) );
#endif

   /* OPTIONAL: including the default cuts for Benders' decomposition */
#if 0
   SCIP_CALL( SCIPincludeBendersDefaultCuts(scip, benders) );
#endif

   /* add xyz Benders' decomposition parameters */
   /* TODO: (optional) add Benders' decomposition specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
