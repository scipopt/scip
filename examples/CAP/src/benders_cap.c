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

/**@file   benders_cap.c
 * @brief  Benders' decomposition algorithm for the capacitated facility location problem
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scipdefplugins.h"
#include "scip/pub_benders.h"
#include "scip/bendersdefcuts.h"
#include "scip/cons_linear.h"

#include "benders_cap.h"
#include "vardata_cap.h"
#include "probdata_cap.h"



#define BENDERS_NAME            "cap"
#define BENDERS_DESC            "Benders' decomposition template"
#define BENDERS_PRIORITY        0
#define BENDERS_CUTLP         TRUE   /**< should Benders' cut be generated for LP solutions */
#define BENDERS_CUTPSEUDO     TRUE   /**< should Benders' cut be generated for pseudo solutions */
#define BENDERS_CUTRELAX      TRUE   /**< should Benders' cut be generated for relaxation solutions */


/*
 * Data structures
 */

/* TODO: fill in the necessary Benders' decomposition data */

/** Benders' decomposition data */
struct SCIP_BendersData
{
   SCIP**                subproblems;        /**< the Benders' decomposition subproblems */
};




/*
 * Local methods
 */

/** creates the Benders' decomposition subproblem */
static
SCIP_RETCODE createSubproblem(
   SCIP*                 subproblem,         /**< the Benders' decomposition subproblem */
   SCIP_PROBDATA*        probdata,           /**< the probdata for the original problem */
   int                   probnumber          /**< the subproblem number that is being created */
   )
{
   SCIP_CONS* demandcons;
   SCIP_CONS** capconss;
   SCIP_CONS* cons;
   SCIP_VAR** facilityvars;
   SCIP_Real** costs;
   SCIP_Real* demands;
   SCIP_Real capacity;
   int nfacilities;
   SCIP_VAR* var;
   SCIP_VARDATA* vardata;
   int i;
   int j;
   char name[SCIP_MAXSTRLEN];

   assert(subproblem != NULL);

   facilityvars = SCIPprobdataGetFacilityVars(probdata);
   costs = SCIPprobdataGetCosts(probdata);
   demands = SCIPprobdataGetDemands(probdata);
   capacity = SCIPprobdataGetCapacity(probdata);
   nfacilities = SCIPprobdataGetNFacilities(probdata);

   SCIP_CALL( SCIPallocBufferArray(subproblem, &capconss, nfacilities) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(subproblem) );

   /* adds the capacity constraints to the subproblem */
   for( i = 0; i < nfacilities; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "capacity_%d", i);
      SCIP_CALL( SCIPcreateConsBasicLinear(subproblem, &cons, name, 0, NULL, NULL, -SCIPinfinity(subproblem), 0.0) );

      SCIP_CALL( SCIPaddCons(subproblem, cons) );

      capconss[i] = cons;
   }


   /* adds the demand constraints to the subproblem */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "demand_%d", probnumber);
   SCIP_CALL( SCIPcreateConsBasicLinear(subproblem, &cons, name, 0, NULL, NULL, demands[probnumber],
         SCIPinfinity(subproblem)) );

   SCIP_CALL( SCIPaddCons(subproblem, cons) );

   demandcons = cons;


   for( i = 0; i < nfacilities; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "subprob_%d_facility_%d", probnumber, i);
      SCIP_CALL( SCIPcreateVarCAP(subproblem, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );

      SCIP_CALL( SCIPaddVar(subproblem, var) );

      /* creates the variable data */
      SCIP_CALL( SCIPvardataCreateCAP(subproblem, &vardata, SUBPROB, 1) );

      /* adds the master variable to the variable mapping */
      SCIPvardataAddVarMapping(vardata, facilityvars[i], -1);

      /* add the variable data to the variable */
      SCIPvarSetData(var, vardata);

      /* getting the variable data for the master variable */
      vardata = SCIPvarGetData(facilityvars[i]);

      /* adds the subproblem variable to the variable mapping */
      SCIPvardataAddVarMapping(vardata, var, probnumber);

      /* adding the variable to the capacity constriants */
      SCIP_CALL( SCIPaddCoefLinear(subproblem, capconss[i], var, -capacity) );

      /* releases the variable */
      SCIP_CALL( SCIPreleaseVar(subproblem, &var) );
   }

   /* adding the customer variables to the subproblem */
   for( j = 0; j < nfacilities; j++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "subprob_%d_customer(%d,%d)", probnumber, probnumber, j);
      SCIP_CALL( SCIPcreateVarCAP(subproblem, &var, name, 0, SCIPinfinity(subproblem), costs[j][probnumber],
            SCIP_VARTYPE_CONTINUOUS) );

      SCIP_CALL( SCIPaddVar(subproblem, var) );

      /* creates the variable data */
      SCIP_CALL( SCIPvardataCreateCAP(subproblem, &vardata, SUBPROB, 1) );

      /* add the variable data to the variable */
      SCIPvarSetData(var, vardata);

      if( costs[j][probnumber] > 0 )
      {
         /* adding the variable to the capacity constriants */
         SCIP_CALL( SCIPaddCoefLinear(subproblem, capconss[j], var, 1.0) );

         /* adding the variable to the demand constraints */
         SCIP_CALL( SCIPaddCoefLinear(subproblem, demandcons, var, 1.0) );
      }

      /* releases the variable */
      SCIP_CALL( SCIPreleaseVar(subproblem, &var) );
   }

   SCIPfreeBufferArray(subproblem, &capconss);

   return SCIP_OKAY;
}


/*
 * Callback methods for Benders' decomposition
 */

/* TODO: Implement all necessary Benders' decomposition methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for benders plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BENDERSCOPY(bendersCopyCap)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cap Benders' decompostion not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersCopyCap NULL
#endif

/** destructor of Benders' decomposition to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BENDERSFREE(bendersFreeCap)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;
   int nsubproblems;
   int i;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);
   nsubproblems = SCIPbendersGetNSubproblems(benders);

   for( i = 0; i < nsubproblems; i++ )
   {
      SCIP_CALL( SCIPfree(&bendersdata->subproblems[i]) );
   }

   SCIPfreeBlockMemoryArray(scip, &bendersdata->subproblems, SCIPbendersGetNSubproblems(benders));

   return SCIP_OKAY;
}


/** initialization method of Benders' decomposition (called after problem was transformed) */
#if 0
static
SCIP_DECL_BENDERSINIT(bendersInitCap)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cap Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitCap NULL
#endif



/** deinitialization method of Benders' decomposition (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BENDERSEXIT(bendersExitCap)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cap Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitCap NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin)
 *
 *  This function is called immediately after the auxiliary variables are created in the master problem. The callback
 *  provides the user an opportunity to add variable data to the auxiliary variables.
 */
#if 0
static
SCIP_DECL_BENDERSINITPRE(bendersInitpreCap)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cap Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitpreCap NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_BENDERSEXITPRE(bendersExitpreCap)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cap Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitpreCap NULL
#endif


/** solving process initialization method of Benders' decomposition (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BENDERSINITSOL(bendersInitsolCap)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}
#else
#define bendersInitsolCap NULL
#endif


/** solving process deinitialization method of Benders' decomposition (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BENDERSEXITSOL(bendersExitsolCap)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cap Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitsolCap NULL
#endif


/** mapping method between the master problem variables and the subproblem variables of Benders' decomposition */
static
SCIP_DECL_BENDERSGETVAR(bendersGetvarCap)
{  /*lint --e{715}*/
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(benders != NULL);
   assert(var != NULL);

   vardata = SCIPvarGetData(var);

   assert((SCIPvardataGetProb(vardata) == MASTER && probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders))
      || (SCIPvardataGetProb(vardata) == SUBPROB && probnumber == -1));

   return SCIPvardataGetMappedVar(vardata, probnumber);
}

/** the execution method for Benders' decomposition */
static
SCIP_DECL_BENDERSEXEC(bendersExecCap)
{  /*lint --e{715}*/

   //SCIP_CALL( SCIPprintTransProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** the method for creating the Benders' decomposition subproblem. This method is called during the initialisation stage
 *  (after the master problem was transformed)
 *
 *  This method must create the SCIP instance for the subproblem and add the required variables and constraints. In
 *  addition, the settings required for the solving the problem must be set here. However, some settings will be
 *  overridden by the standard solving method included in the Benders' decomposition framework. If a special solving
 *  method is desired, the user can implement the bendersSolvesubXyz callback.
 */
static
SCIP_DECL_BENDERSCREATESUB(bendersCreatesubCap)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;
   SCIP_PROBDATA* probdata;
   char name[SCIP_MAXSTRLEN];


   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);
   probdata = SCIPgetProbData(scip);

   /* creating the subproblem */
   SCIP_CALL( SCIPcreate(&bendersdata->subproblems[probnumber]) );

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sub_%d", SCIPgetProbName(scip), probnumber);

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(bendersdata->subproblems[probnumber], name) );

   SCIP_CALL( createSubproblem(bendersdata->subproblems[probnumber], probdata, probnumber) );

   /* adding the subproblem to the Benders' decomposition structure */
   SCIP_CALL( SCIPaddBendersSubproblem(scip, benders, bendersdata->subproblems[probnumber]) );

   return SCIP_OKAY;
}



/** the subproblem solving method for Benders' decomposition */
#if 0
static
SCIP_DECL_BENDERSSOLVESUB(bendersSolvesubCap)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cap Benders' decomposition not implemented yet\n");

   return SCIP_OKAY;
}
#else
#define bendersSolvesubCap NULL
#endif


#if 0
/** the post-solve method for Benders' decomposition */
static
SCIP_DECL_BENDERSPOSTSOLVE(bendersPostsolveCap)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cap Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersPostsolveCap NULL
#endif


/** the subproblem freeing method for Benders' decomposition */
#if 0
static
SCIP_DECL_BENDERSFREESUB(bendersFreesubCap)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cap Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersFreesubCap NULL
#endif




/*
 * Benders' decomposition specific interface methods
 */

/** creates the cap Benders' decomposition and includes it in SCIP */
SCIP_RETCODE SCIPincludeBendersCap(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   )
{
   SCIP_BENDERS* benders;
   SCIP_BENDERSDATA* bendersdata;

   /* create cap Benders' decomposition data */
   benders = NULL;
   bendersdata = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &bendersdata) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &bendersdata->subproblems, nsubproblems) );

   /* include Benders' decomposition */
#if 0
   /* use SCIPincludeBenders() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBenders(scip, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, nsubproblems,
         bendersCopyCap, bendersFreeCap, bendersInitCap, bendersExitCap, bendersInitpreCap, bendersExitpreCap,
         bendersInitsolCap, bendersExitsolCap, bendersGetvarCap, bendersExecCap, bendersCreatesubCap,
         bendersSolvesubCap, bendersPostsolveCap, bendersFreesubCap, bendersdata) );
#else
   /* use SCIPincludeBendersBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBendersBasic(scip, &benders, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, nsubproblems,
         BENDERS_CUTLP, BENDERS_CUTPSEUDO, BENDERS_CUTRELAX, bendersGetvarCap, bendersExecCap, bendersCreatesubCap,
         bendersdata) );
   assert(benders != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBendersCopy(scip, benders, bendersCopyCap) );
   SCIP_CALL( SCIPsetBendersFree(scip, benders, bendersFreeCap) );
   SCIP_CALL( SCIPsetBendersInit(scip, benders, bendersInitCap) );
   SCIP_CALL( SCIPsetBendersExit(scip, benders, bendersExitCap) );
   SCIP_CALL( SCIPsetBendersInitpre(scip, benders, bendersInitpreCap) );
   SCIP_CALL( SCIPsetBendersExitpre(scip, benders, bendersExitpreCap) );
   SCIP_CALL( SCIPsetBendersInitsol(scip, benders, bendersInitsolCap) );
   SCIP_CALL( SCIPsetBendersExitsol(scip, benders, bendersExitsolCap) );
   SCIP_CALL( SCIPsetBendersSolvesub(scip, benders, bendersSolvesubCap) );
   SCIP_CALL( SCIPsetBendersPostsolve(scip, benders, bendersPostsolveCap) );
   SCIP_CALL( SCIPsetBendersFreesub(scip, benders, bendersFreesubCap) );
#endif

   /* including the default cuts for Benders' decomposition */
   SCIP_CALL( SCIPincludeBendersDefaultCuts(scip, benders) );

   /* add cap Benders' decomposition parameters */
   /* TODO: (optional) add Benders' decomposition specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
