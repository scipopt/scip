/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define SCIP_STATISTIC
/* #define SCIP_DEBUG */

/**@file   branch_cloud.c
 * @brief  cloud branching rule
 * @author Timo Berthold
 * @author Domenico Salvagnin
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_cloud.h"
#include "scip/branch_fullstrong.h"


#define BRANCHRULE_NAME            "cloud"
#define BRANCHRULE_DESC            "branching rule that considers several alternative LP optima"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define DEFAULT_USECLOUD           TRUE      /**< should a cloud of points be used? */

/*
 * Data structures
 */

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
   int                   lastcand;           /**< last evaluated candidate of last branching rule execution */
   SCIP_Bool             usecloud;           /**< should a cloud of points be used? */
   SCIP_CLOCK*           cloudclock;         /**< clock for cloud diving */
   int                   ntried;             /**< number of times the cloud was tried */
   int                   nuseful;            /**< number of times the cloud was useful */
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
SCIP_DECL_BRANCHCOPY(branchCopyCloud)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cloud branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchCopyCloud NULL
#endif

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeCloud)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPstatisticMessage("time spent diving in cloud branching: %g\n", SCIPgetClockTime(scip, branchruledata->cloudclock));
   SCIPstatisticMessage("success rate of cloud branching: %g\n", (SCIP_Real)branchruledata->nuseful / branchruledata->ntried);
   SCIP_CALL( SCIPfreeClock(scip, &(branchruledata->cloudclock)) );
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitCloud)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* initialize branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   branchruledata->lastcand = 0;
   branchruledata->nuseful = 0;
   branchruledata->ntried = 0;
   SCIP_CALL( SCIPcreateClock(scip, &(branchruledata->cloudclock)) );

   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BRANCHEXIT(branchExitCloud)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cloud branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitCloud NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BRANCHINITSOL(branchInitsolCloud)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cloud branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolCloud NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolCloud)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cloud branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolCloud NULL
#endif


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpCloud)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   SCIP_VAR** lpcands;
   SCIP_VAR** lpcandscopy;

   SCIP_VAR** vars;                          /* SCIP variables                */
   SCIP_ROW** lprows;
   SCIP_Real* lpcandsfrac;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfraccopy;
   SCIP_Real* lpcandssolcopy;
   SCIP_Real* lpcandsmin;
   SCIP_Real* lpcandsmax;

   SCIP_Real bestdown;
   SCIP_Real bestup;
   SCIP_Real bestscore;
   SCIP_Real provedbound;

   SCIP_Bool bestdownvalid;
   SCIP_Bool bestupvalid;
   SCIP_Bool newpoint;
   SCIP_Bool lperror;

   int nlpcands;
   int npriolpcands;
   int nvars;
   int bestcand;
   int nlprows;
   int i;
   int counter;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   SCIPdebugMessage("Execlp method of "BRANCHRULE_NAME" branching\n");

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, &npriolpcands) );
   nlpcands = SCIPgetNLPBranchCands(scip);
   assert(nlpcands > 0);

   /* get problem variables and LP row data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   nlprows = SCIPgetNLPRows(scip);
   lprows = SCIPgetLPRows(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &lpcandsmin, nlpcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lpcandsmax, nlpcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lpcandscopy, nlpcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lpcandsfraccopy, nlpcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lpcandssolcopy, nlpcands) );
   BMScopyMemoryArray(lpcandsmin, lpcandssol, nlpcands);
   BMScopyMemoryArray(lpcandsmax, lpcandssol, nlpcands);
   BMScopyMemoryArray(lpcandssolcopy, lpcandssol, nlpcands);
   BMScopyMemoryArray(lpcandsfraccopy, lpcandsfrac, nlpcands);
   BMScopyMemoryArray(lpcandscopy, lpcands, nlpcands);

   SCIP_CALL(  SCIPstartClock(scip, branchruledata->cloudclock) );
   branchruledata->ntried++;

   /* start diving to calculate the solution cloud */
   SCIP_CALL( SCIPstartDive(scip) );

   /* fix variables with nonzero reduced costs to reduce LP to the optimal face */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real solval;
      solval = SCIPgetSolVal(scip, NULL, vars[i]);

      if( !SCIPisFeasZero(scip, SCIPgetVarRedcost(scip, vars[i])) )
      {
         /* printf("var: %s [%g,%g] : redcost %g solval %g\n",SCIPvarGetName(vars[i]), SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i]),SCIPgetVarRedcost(scip, vars[i]), solval ); */
         SCIP_CALL( SCIPchgVarLbDive(scip, vars[i], solval) );
         SCIP_CALL( SCIPchgVarUbDive(scip, vars[i], solval) );
      }
      else if( SCIPvarGetType(vars[i]) == SCIP_VARTYPE_INTEGER && !SCIPisIntegral(scip, solval) )
      {
         /* printf("var: %s [%g,%g] : redcost %g solval %g\n",SCIPvarGetName(vars[i]), SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i]),SCIPgetVarRedcost(scip, vars[i]), solval ); */
         SCIP_CALL( SCIPchgVarLbDive(scip, vars[i], SCIPfloor(scip, solval)) );
         SCIP_CALL( SCIPchgVarUbDive(scip, vars[i], SCIPceil(scip, solval)) );
      }

      SCIP_CALL( SCIPchgVarObjDive(scip, vars[i], 0.0) );
   }

   /* fix LP rows with nonzero dual solution to reduce LP to the optimal face */
   for( i = 0; i < nlprows; ++i )
   {
      SCIP_Real dualsol;
      dualsol = SCIProwGetDualsol(lprows[i]);
      if( !SCIPisZero(scip, dualsol) )
      {
         if( dualsol > 0 && SCIPisFeasEQ(scip,SCIProwGetLhs(lprows[i]), SCIPgetRowActivity(scip,lprows[i])) )
         {
            /* printf("row %s lhs: %g = activity: %g   rhs: %g dualsol: %g \n", SCIProwGetName(lprows[i]),SCIProwGetLhs(lprows[i]), SCIPgetRowActivity(scip,lprows[i]),SCIProwGetRhs(lprows[i]), dualsol); */
            SCIP_CALL( SCIPchgRowRhsDive(scip, lprows[i], SCIProwGetLhs(lprows[i])) );
         }
         else if( dualsol < 0 && SCIPisFeasEQ(scip,SCIProwGetRhs(lprows[i]), SCIPgetRowActivity(scip,lprows[i])) )
         {
            /* printf("row %s lhs: %g   activity: %g = rhs: %g dualsol: %g \n", SCIProwGetName(lprows[i]),SCIProwGetLhs(lprows[i]), SCIPgetRowActivity(scip,lprows[i]),SCIProwGetRhs(lprows[i]), dualsol); */
            SCIP_CALL( SCIPchgRowLhsDive(scip, lprows[i], SCIProwGetRhs(lprows[i])) );
         }
      }
   }

   newpoint = TRUE;
   counter = 0;

   /* loop that generates new cloud point */
   while( newpoint && branchruledata->usecloud )
   {
#ifdef NDEBUG
      SCIP_RETCODE retcode;
#endif

      /* apply feasibility pump objective function to fractional variables */
      for( i = 0; i < nlpcands; ++i)
      {
         SCIP_Real frac;
         frac = SCIPfrac(scip, SCIPgetSolVal(scip, NULL, lpcandscopy[i]));

         if( !SCIPisZero(scip, frac) && !SCIPisIntegral(scip, lpcandsmin[i]) && !SCIPisIntegral(scip, lpcandsmax[i]) )
         {
            if( frac < 0.5 )
            {
               SCIP_CALL( SCIPchgVarObjDive(scip, vars[i], 1.0) );
            }
            else
            {
               SCIP_CALL( SCIPchgVarObjDive(scip, vars[i], -1.0) );
            }
         }
      }

      /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
       * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
#ifdef NDEBUG
      retcode = SCIPsolveDiveLP(scip, -1, &lperror);
      if( retcode != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Error while solving LP in "BRANCHRULE_NAME"; LP solve terminated with code <%d>\n",retcode);
      }
#else
      SCIP_CALL( SCIPsolveDiveLP(scip, -1, &lperror) );
#endif

      newpoint = FALSE;
      for( i = 0; i < nlpcands; ++i)
      {
         SCIP_Real solval;
         solval = SCIPgetSolVal(scip, NULL, lpcandscopy[i]);

         if( SCIPisFeasIntegral(scip,solval) && !SCIPisFeasIntegral(scip, lpcandsmin[i]) && !SCIPisFeasIntegral(scip, lpcandsmax[i]) )
            newpoint = TRUE;

         lpcandsmin[i] = MIN(lpcandsmin[i], solval);
         lpcandsmax[i] = MAX(lpcandsmax[i], solval);
      }

      if( newpoint )
         counter++;
   }
   SCIPdebugMessage("considered %d additional points in the cloud\n",counter);

   /* terminate the diving */
   SCIP_CALL( SCIPendDive(scip) );

   SCIP_CALL(  SCIPstopClock(scip, branchruledata->cloudclock) );

   if( counter > 0 )
   {
      counter = 0;

      for( i = 0; i < nlpcands; ++i)
      {
         if( !SCIPisFeasIntegral(scip, lpcandsmin[i]) && !SCIPisFeasIntegral(scip, lpcandsmax[i]) )
         {
            assert(counter <= i);
            lpcandscopy[counter] = lpcandscopy[i];
            lpcandssolcopy[counter] = lpcandssolcopy[i];
            lpcandsfraccopy[counter] = lpcandsfraccopy[i];
            counter++;
         }
      }
      SCIPdebugMessage("skipped %d/%d strong branching candidates\n", nlpcands - counter, nlpcands);
      if( nlpcands - counter > 0 )
         branchruledata->nuseful++;
   }
   else
      counter = nlpcands;

   SCIP_CALL( SCIPselectVarStrongBranching(scip, lpcandscopy, lpcandssolcopy, lpcandsfraccopy, counter, counter, /* replace second counter ??????????? */
      &branchruledata->lastcand, allowaddcons,
      &bestcand, &bestdown, &bestup, &bestscore, &bestdownvalid, &bestupvalid, &provedbound, result) );

   /* perform the branching */
   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM && *result != SCIP_CONSADDED && counter > 0 ) /* ??????? */
   {
      SCIP_NODE* downchild;
      SCIP_NODE* upchild;
      SCIP_VAR* var;
      SCIP_Bool allcolsinlp;
      SCIP_Bool exactsolve;

      assert(*result == SCIP_DIDNOTRUN);
      assert(0 <= bestcand && bestcand < nlpcands);
      assert(SCIPisLT(scip, provedbound, SCIPgetCutoffbound(scip)));

      var = lpcandscopy[bestcand];

      /* perform the branching */
      SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%g, down=%g, up=%g, score=%g)\n",
         counter, bestcand, SCIPvarGetName(var), lpcandssolcopy[bestcand], bestdown, bestup, bestscore);
      SCIP_CALL( SCIPbranchVar(scip, var, &downchild, NULL, &upchild) );
      assert(downchild != NULL);
      assert(upchild != NULL);

      /* check, if we want to solve the problem exactly, meaning that strong branching information is not useful
       * for cutting off sub problems and improving lower bounds of children
       */
      exactsolve = SCIPisExactSolve(scip);

      /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
      allcolsinlp = SCIPallColsInLP(scip);

      /* update the lower bounds in the children */
      if( allcolsinlp && !exactsolve )
      {
         SCIP_CALL( SCIPupdateNodeLowerbound(scip, downchild, bestdownvalid ? MAX(bestdown, provedbound) : provedbound) );
         SCIP_CALL( SCIPupdateNodeLowerbound(scip, upchild, bestupvalid ? MAX(bestup, provedbound) : provedbound) );
      }
      SCIPdebugMessage(" -> down child's lowerbound: %g\n", SCIPnodeGetLowerbound(downchild));
      SCIPdebugMessage(" -> up child's lowerbound: %g\n", SCIPnodeGetLowerbound(upchild));

      *result = SCIP_BRANCHED;
   }

   SCIPfreeBufferArray(scip, &lpcandscopy);
   SCIPfreeBufferArray(scip, &lpcandssolcopy);
   SCIPfreeBufferArray(scip, &lpcandsfraccopy);
   SCIPfreeBufferArray(scip, &lpcandsmax);
   SCIPfreeBufferArray(scip, &lpcandsmin);

   return SCIP_OKAY;
}


/** branching execution method for external candidates */
#if 0
static
SCIP_DECL_BRANCHEXECEXT(branchExecextCloud)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cloud branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecextCloud NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
#if 0
static
SCIP_DECL_BRANCHEXECPS(branchExecpsCloud)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cloud branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecpsCloud NULL
#endif


/*
 * branching rule specific interface methods
 */

/** creates the cloud branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleCloud(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create cloud branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
   branchruledata->lastcand = 0;

   /* include branching rule */
   branchrule = NULL;
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyCloud) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeCloud) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitCloud) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitCloud) );
   SCIP_CALL( SCIPsetBranchruleInitsol(scip, branchrule, branchInitsolCloud) );
   SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolCloud) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpCloud) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextCloud) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsCloud) );

   /* add cloud branching rule parameters */

   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/"BRANCHRULE_NAME"/usecloud",
         "should a cloud of points be used? ",
         &branchruledata->usecloud, FALSE, DEFAULT_USECLOUD, NULL, NULL) );

   return SCIP_OKAY;
}
