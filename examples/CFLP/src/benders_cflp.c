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

/**@file   benders_cflp.c
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

#include "benders_cflp.h"
#include "vardata_cflp.h"
#include "probdata_cflp.h"



#define BENDERS_NAME            "cflp"
#define BENDERS_DESC            "Benders' decomposition template"
#define BENDERS_PRIORITY        0
#define BENDERS_CUTLP         TRUE   /**< should Benders' cut be generated for LP solutions */
#define BENDERS_CUTPSEUDO     TRUE   /**< should Benders' cut be generated for pseudo solutions */
#define BENDERS_CUTRELAX      TRUE   /**< should Benders' cut be generated for relaxation solutions */


/* event handler properties */
#define EVENTHDLR_NAME         "benders_cflp"
#define EVENTHDLR_DESC         "node focus event handler for Benders' decomposition"

/* ---------------- Callback methods of event handler ---------------- */

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecCflp)
{  /*lint --e{715}*/

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* sending an interrupt solve signal to return the control back to the Benders' decomposition plugin.
    * This will ensure the SCIP stage is SCIP_STAGE_SOLVING, allowing the use of probing mode. */
   SCIP_CALL( SCIPinterruptSolve(scip) );

   SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1));

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolCflp)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL));

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolCflp)
{
   assert(scip != NULL);

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1));

   return SCIP_OKAY;
}



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
      SCIP_CALL( SCIPcreateVarCFLP(subproblem, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );

      SCIP_CALL( SCIPaddVar(subproblem, var) );

      /* creates the variable data */
      SCIP_CALL( SCIPvardataCreateCFLP(subproblem, &vardata, SUBPROB, 1) );

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
      SCIP_CALL( SCIPcreateVarCFLP(subproblem, &var, name, 0, SCIPinfinity(subproblem), costs[j][probnumber],
            SCIP_VARTYPE_CONTINUOUS) );

      SCIP_CALL( SCIPaddVar(subproblem, var) );

      /* creates the variable data */
      SCIP_CALL( SCIPvardataCreateCFLP(subproblem, &vardata, SUBPROB, 1) );

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


/** sets up the subproblem with the solution to the master problem */
/*  to set up the subproblem, we need to fix all facility variables to the master problem solution */
static
SCIP_RETCODE setupSubproblem(
   SCIP*                 masterprob,         /**< the master problem */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SOL*             sol,                /**< the current solution to the master problem */
   int                   probnumber          /**< the subproblem number */
   )
{
   SCIP* subproblem;
   SCIP_PROBDATA* probdata;
   SCIP_VAR** facilityvars;
   SCIP_VARDATA* vardata;
   SCIP_Real solval;
   int nfacilities;
   int i;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;

   assert(masterprob != NULL);
   assert(benders != NULL);

   subproblem = SCIPbendersSubproblem(benders, probnumber);

   /* ending probing mode to reset the current node */
   SCIP_CALL( SCIPendProbing(subproblem) );

   /* starting probing mode to fix variables for the subproblem */
   SCIP_CALL( SCIPstartProbing(subproblem) );

   probdata = SCIPgetProbData(masterprob);
   nfacilities = SCIPprobdataGetNFacilities(probdata);
   facilityvars = SCIPprobdataGetFacilityVars(probdata);

   /* looping over all facility variables and fixing the variables in the subproblem */
   for( i = 0; i < nfacilities; i++ )
   {
      solval = SCIPgetSolVal(masterprob, sol, facilityvars[i]);
      vardata = SCIPvarGetData(facilityvars[i]);

      /* fixing the variable in the subproblem */
      SCIP_CALL( SCIPfixVar(subproblem, SCIPvardataGetMappedVar(vardata, probnumber), solval, &infeasible, &fixed) );

      assert(fixed);
      assert(!infeasible);
   }

   return SCIP_OKAY;
}


/* solving the lp of the subproblems to generated Benders' cuts */
static
SCIP_RETCODE solveSubproblemLP(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   SCIP_Real*            objval,             /**< the sum of the objective function values over all subproblems */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible          /**< a flag to indicate whether all subproblems are feasible */
   )
{
   SCIP* subproblem;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;

   /* previous parameter settings */
   int prevCutoffParam;
   char prevInitAlgParam;
   char prevResolveAlgParam;
   SCIP_Bool prevDualParam;

   assert(masterprob != NULL);
   assert(benders != NULL);
   assert(objval != NULL);
   assert(infeasible != NULL);

   (*objval) = 0;
   (*infeasible) = FALSE;


   /* TODO: This should be solved just as an LP, so as a MIP. There is too much overhead with the MIP.
    * Need to change status check for checking the LP. */
   subproblem = SCIPbendersSubproblem(benders, probnumber);

   assert(SCIPisLPConstructed(subproblem));
   /* modifying all of the parameters */

   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/disablecutoff", &prevCutoffParam) );
   SCIPsetIntParam(subproblem, "lp/disablecutoff", 1);

   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/initalgorithm", &prevInitAlgParam) );
   SCIPsetCharParam(subproblem, "lp/initalgorithm", 'd');
   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/resolvealgorithm", &prevResolveAlgParam) );
   SCIPsetCharParam(subproblem, "lp/resolvealgorithm", 'd');

   SCIP_CALL( SCIPgetBoolParam(subproblem, "misc/alwaysgetduals", &prevDualParam) );
   SCIPsetBoolParam(subproblem, "misc/alwaysgetduals", TRUE);

   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );

   SCIP_CALL( SCIPsolveProbingLP(subproblem, -1, &lperror, &cutoff) );

   assert(!lperror);

   if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      (*objval) += SCIPgetLPObjval(subproblem);
   }
   else if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_INFEASIBLE )
   {
      (*infeasible) = TRUE;
      (*objval) = SCIPinfinity(masterprob);
   }
   else
      assert(FALSE);

   //SCIP_CALL( SCIPprintStatistics(subprob, NULL) );

   SCIPsetIntParam(subproblem, "lp/disablecutoff", prevCutoffParam);
   SCIPsetCharParam(subproblem, "lp/initalgorithm", prevInitAlgParam);
   SCIPsetCharParam(subproblem, "lp/resolvealgorithm", prevResolveAlgParam);
   SCIPsetBoolParam(subproblem, "misc/alwaysgetduals", prevDualParam);

   return SCIP_OKAY;
}





/* solving the subproblems to generated Benders' cuts */
static
SCIP_RETCODE solveSubproblem(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   SCIP_Real*            objval,             /**< the sum of the objective function values over all subproblems */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible          /**< a flag to indicate whether all subproblems are feasible */
   )
{
   SCIP* subproblem;
   SCIP_SOL* bestsol;

   /* previous parameter settings */
   int prevCutoffParam;
   int prevPropMaxroundsParam;
   int prevPropMaxroundsRootParam;
   char prevInitAlgParam;
   char prevResolveAlgParam;
   SCIP_Bool prevConfParam;
   SCIP_Bool prevDualParam;

   assert(masterprob != NULL);
   assert(benders != NULL);
   assert(objval != NULL);
   assert(infeasible != NULL);

   (*objval) = 0;
   (*infeasible) = FALSE;


   /* TODO: This should be solved just as an LP, so as a MIP. There is too much overhead with the MIP.
    * Need to change status check for checking the LP. */
   subproblem = SCIPbendersSubproblem(benders, probnumber);

   /* modifying all of the parameters */

   /* Do we have to disable presolving? If yes, we have to store all presolving parameters. */
   SCIPsetPresolving(subproblem, SCIP_PARAMSETTING_OFF, TRUE);

   /* Disabling heuristics so that the problem is not trivially solved */
   SCIPsetHeuristics(subproblem, SCIP_PARAMSETTING_OFF, TRUE);

   /* store parameters that are changed for the generation of the subproblem cuts */
   SCIP_CALL( SCIPgetBoolParam(subproblem, "conflict/enable", &prevConfParam) );
   SCIPsetParam(subproblem, "conflict/enable", FALSE);

   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/disablecutoff", &prevCutoffParam) );
   SCIPsetIntParam(subproblem, "lp/disablecutoff", 1);

   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/initalgorithm", &prevInitAlgParam) );
   SCIPsetCharParam(subproblem, "lp/initalgorithm", 'd');
   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/resolvealgorithm", &prevResolveAlgParam) );
   SCIPsetCharParam(subproblem, "lp/resolvealgorithm", 'd');

   SCIP_CALL( SCIPgetBoolParam(subproblem, "misc/alwaysgetduals", &prevDualParam) );
   SCIPsetBoolParam(subproblem, "misc/alwaysgetduals", TRUE);

   //SCIPinfoMessage(subproblem, NULL, "Pricing problem %d\n", probnumber);
   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
   //SCIP_CALL( SCIPsetBoolParam(subproblem, "display/lpinfo", TRUE) );

   SCIP_CALL( SCIPgetIntParam(subproblem, "propagating/maxrounds", &prevPropMaxroundsParam) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxrounds", 0) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "propagating/maxroundsroot", &prevPropMaxroundsRootParam) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxroundsroot", 0) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "constraints/linear/propfreq", -1) );

   SCIP_CALL( SCIPsolve(subproblem) );

   bestsol = SCIPgetBestSol(subproblem);

   if( SCIPgetStatus(subproblem) == SCIP_STATUS_OPTIMAL )
   {
      (*objval) += SCIPgetSolTransObj(subproblem, bestsol);
   }
   else if( SCIPgetStatus(subproblem) == SCIP_STATUS_INFEASIBLE )
   {
      (*infeasible) = TRUE;
      (*objval) = SCIPinfinity(masterprob);
   }
   else if( SCIPgetStatus(subproblem) != SCIP_STATUS_USERINTERRUPT )
      assert(FALSE);

   //SCIP_CALL( SCIPprintStatistics(subprob, NULL) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
   /* resetting the parameter settings to the previous state */
   SCIPsetPresolving(subproblem, SCIP_PARAMSETTING_DEFAULT, TRUE);
   SCIPsetHeuristics(subproblem, SCIP_PARAMSETTING_DEFAULT, TRUE);
   SCIPsetBoolParam(subproblem, "conflict/enable", prevConfParam);
   SCIPsetIntParam(subproblem, "lp/disablecutoff", prevCutoffParam);
   SCIPsetCharParam(subproblem, "lp/initalgorithm", prevInitAlgParam);
   SCIPsetCharParam(subproblem, "lp/resolvealgorithm", prevResolveAlgParam);
   SCIPsetBoolParam(subproblem, "misc/alwaysgetduals", prevDualParam);
   SCIPsetIntParam(subproblem, "propagating/maxrounds", prevPropMaxroundsParam);
   SCIPsetIntParam(subproblem, "propagating/maxroundsroot", prevPropMaxroundsRootParam);

   return SCIP_OKAY;
}





/*
 * Callback methods for Benders' decomposition
 */

/* TODO: Implement all necessary Benders' decomposition methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for benders plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BENDERSCOPY(bendersCopyCflp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cflp Benders' decompostion not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersCopyCflp NULL
#endif

/** destructor of Benders' decomposition to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BENDERSFREE(bendersFreeCflp)
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
      SCIP_CALL( SCIPendProbing(scip) );
      SCIP_CALL( SCIPfree(&bendersdata->subproblems[i]) );
   }

   SCIPfreeBlockMemoryArray(scip, &bendersdata->subproblems, SCIPbendersGetNSubproblems(benders));

   return SCIP_OKAY;
}


/** initialization method of Benders' decomposition (called after problem was transformed) */
static
SCIP_DECL_BENDERSINIT(bendersInitCflp)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;
   SCIP_PROBDATA* probdata;
   SCIP_Real objval;
   int nsubproblems;
   char name[SCIP_MAXSTRLEN];
   SCIP_Bool infeasible;
   int i;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);
   probdata = SCIPgetProbData(scip);
   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* creating the subproblems */
   for( i = 0; i < nsubproblems; i++ )
   {
      SCIP_EVENTHDLR* eventhdlr;
      SCIP_Bool cutoff;

      SCIP_CALL( SCIPcreate(&bendersdata->subproblems[i]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sub_%d", SCIPgetProbName(scip), i);

      /* create problem in SCIP and add non-NULL callbacks via setter functions */
      SCIP_CALL( SCIPcreateProbBasic(bendersdata->subproblems[i], name) );

      /* include event handler into SCIP */
      SCIP_CALL( SCIPincludeEventhdlrBasic(bendersdata->subproblems[i], &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
            eventExecCflp, NULL) );
      SCIP_CALL( SCIPsetEventhdlrInitsol(bendersdata->subproblems[i], eventhdlr, eventInitsolCflp) );
      SCIP_CALL( SCIPsetEventhdlrExitsol(bendersdata->subproblems[i], eventhdlr, eventExitsolCflp) );
      assert(eventhdlr != NULL);

      SCIP_CALL( createSubproblem(bendersdata->subproblems[i], probdata, i) );

      /* adding the subproblem to the Benders' decomposition structure */
      SCIP_CALL( SCIPaddBendersSubproblem(scip, benders, bendersdata->subproblems[i]) );

      /* Getting the problem into the right SCIP stage for solving */
      SCIP_CALL( solveSubproblem(scip, benders, &objval, i, &infeasible) );

      /* Constructing the LP that can be solved in later iterations */
      SCIP_CALL( SCIPconstructLP(bendersdata->subproblems[i], &cutoff) );

      /* starting probing mode */
      SCIP_CALL( SCIPstartProbing(bendersdata->subproblems[i]) );
   }


   return SCIP_OKAY;
}



/** deinitialization method of Benders' decomposition (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BENDERSEXIT(bendersExitCflp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cflp Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitCflp NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin)
 *
 *  This function is called immediately after the auxiliary variables are created in the master problem. The callback
 *  provides the user an opportunity to add variable data to the auxiliary variables.
 */
#if 0
static
SCIP_DECL_BENDERSINITPRE(bendersInitpreCflp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cflp Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitpreCflp NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_BENDERSEXITPRE(bendersExitpreCflp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cflp Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitpreCflp NULL
#endif


/** solving process initialization method of Benders' decomposition (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BENDERSINITSOL(bendersInitsolCflp)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}
#else
#define bendersInitsolCflp NULL
#endif


/** solving process deinitialization method of Benders' decomposition (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BENDERSEXITSOL(bendersExitsolCflp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cflp Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitsolCflp NULL
#endif


/** mapping method between the master problem variables and the subproblem variables of Benders' decomposition */
static
SCIP_DECL_BENDERSGETMASTERVAR(bendersGetmastervarCflp)
{  /*lint --e{715}*/
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(benders != NULL);
   assert(var != NULL);

   vardata = SCIPvarGetData(var);

   assert(SCIPvardataGetProb(vardata) == SUBPROB);

   return SCIPvardataGetMappedVar(vardata, -1);
}

/** the execution method for Benders' decomposition */
static
SCIP_DECL_BENDERSEXEC(bendersExecCflp)
{  /*lint --e{715}*/

   //SCIP_CALL( SCIPprintTransProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** the subproblem solving method for Benders' decomposition */
static
SCIP_DECL_BENDERSSOLVESUB(bendersSolvesubCflp)
{  /*lint --e{715}*/
   SCIP_Real objval;

   assert(scip != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));
   assert(SCIPbendersSubproblem(benders, probnumber) != NULL);

   SCIP_CALL( setupSubproblem(scip, benders, sol, probnumber) );

   /* solves the subproblem */
   SCIP_CALL( solveSubproblemLP(scip, benders, &objval, probnumber, infeasible) );

   return SCIP_OKAY;
}

#if 0
/** the post-solve method for Benders' decomposition */
static
SCIP_DECL_BENDERSPOSTSOLVE(bendersPostsolveCflp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cflp Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersPostsolveCflp NULL
#endif


/** the subproblem freeing method for Benders' decomposition */
static
SCIP_DECL_BENDERSFREESUB(bendersFreesubCflp)
{  /*lint --e{715}*/

   /* freeing the transform of the subproblems so that it can be updated in the next iteration. */
   //SCIP_CALL( SCIPfreeTransform(SCIPbendersSubproblem(benders, probnumber)) );

   return SCIP_OKAY;
}




/*
 * Benders' decomposition specific interface methods
 */

/** creates the cflp Benders' decomposition and includes it in SCIP */
SCIP_RETCODE SCIPincludeBendersCflp(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   )
{
   SCIP_BENDERS* benders;
   SCIP_BENDERSDATA* bendersdata;

   /* create cflp Benders' decomposition data */
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
         bendersCopyCflp, bendersFreeCflp, bendersInitCflp, bendersExitCflp, bendersInitpreCflp, bendersExitpreCflp,
         bendersInitsolCflp, bendersExitsolCflp, bendersGetmastervarCflp, bendersExecCflp,
         bendersPostsolveCflp, bendersFreesubCflp, bendersdata) );
#else
   /* use SCIPincludeBendersBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBendersBasic(scip, &benders, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, nsubproblems,
         BENDERS_CUTLP, BENDERS_CUTPSEUDO, BENDERS_CUTRELAX, bendersGetmastervarCflp, bendersExecCflp, bendersSolvesubCflp,
         bendersFreesubCflp, bendersdata) );
   assert(benders != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBendersCopy(scip, benders, bendersCopyCflp) );
   SCIP_CALL( SCIPsetBendersFree(scip, benders, bendersFreeCflp) );
   SCIP_CALL( SCIPsetBendersInit(scip, benders, bendersInitCflp) );
   SCIP_CALL( SCIPsetBendersExit(scip, benders, bendersExitCflp) );
   SCIP_CALL( SCIPsetBendersInitpre(scip, benders, bendersInitpreCflp) );
   SCIP_CALL( SCIPsetBendersExitpre(scip, benders, bendersExitpreCflp) );
   SCIP_CALL( SCIPsetBendersInitsol(scip, benders, bendersInitsolCflp) );
   SCIP_CALL( SCIPsetBendersExitsol(scip, benders, bendersExitsolCflp) );
   SCIP_CALL( SCIPsetBendersPostsolve(scip, benders, bendersPostsolveCflp) );
#endif

   /* including the default cuts for Benders' decomposition */
   SCIP_CALL( SCIPincludeBendersDefaultCuts(scip, benders) );

   /* OPTIONAL: including the default cuts for Benders' decomposition */
#if 0
   SCIP_CALL( SCIPincludeBendersDefaultCuts(scip, benders) );
#endif

   /* add cflp Benders' decomposition parameters */
   /* TODO: (optional) add Benders' decomposition specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
