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
#include "scip/branch_allfullstrong.h"


#define BRANCHRULE_NAME            "cloud"
#define BRANCHRULE_DESC            "branching rule that considers several alternative LP optima"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define DEFAULT_USECLOUD           TRUE      /**< should a cloud of points be used? */
#define DEFAULT_USEUNION           FALSE     /**< should the union of candidates be used? */
#define DEFAULT_MAXPOINTS          -1        /**< maximum number of points for the cloud (-1 means no limit) */
#define DEFAULT_MINSUCCESSRATE     0.0       /**< minimum success rate for the cloud */
#define DEFAULT_MINSUCCESSUNION    0.0       /**< minimum success rate for the union */
#define DEFAULT_MAXDEPTHUNION      65000     /**< maximum depth for the union */
#define DEFAULT_ONLYF2             FALSE     /**< should only F2 be used? */

/*
 * Data structures
 */

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
   int                   lastcand;           /**< last evaluated candidate of last branching rule execution */
   SCIP_Bool             usecloud;           /**< should a cloud of points be used? */
   SCIP_Bool             useunion;           /**< should the union of candidates be used? */
   SCIP_Bool             onlyF2;             /**< should only F2 be used? */
   int                   maxpoints;          /**< maximum number of points for the cloud (-1 means no limit) */
   SCIP_Real             minsuccessrate;     /**< minimum success rate for the cloud */
   SCIP_Real             minsuccessunion;    /**< minimum success rate for the union */
   SCIP_CLOCK*           cloudclock;         /**< clock for cloud diving */
   SCIP_Bool*            skipdown;
   SCIP_Bool*            skipup;
   int                   ntried;             /**< number of times the cloud was tried */
   int                   ntriedunions;       /**< number of times the cloud was tried */
   int                   nuseful;            /**< number of times the cloud was useful (at least one LP skipped) */
   int                   nusefulunions;      /**< number of times the union was useful (took candidate from new list) */
   int                   ncloudpoints;       /**< sum of cloud points taken over all nodes with at least two poitns in cloud */
   int                   nsavedlps;          /**< sum of saved LPs taken over all nodes with at least two points in cloud */
   int                   maxdepthunion;      /**< maximum depth for the union */
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/** hash key retrieval function for variables */
static
SCIP_DECL_HASHGETKEY(hashGetKeyVar)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff the indices of both variables are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqVar)
{  /*lint --e{715}*/
   return (key1 == key2);
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValVar)
{  /*lint --e{715}*/
   assert(SCIPvarGetIndex((SCIP_VAR*) key) >= 0);
   return (unsigned int) SCIPvarGetIndex((SCIP_VAR*) key);
}

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

   if( branchruledata->cloudclock != NULL)
   {
      int ntried = branchruledata->ntried;
      int nuseful = branchruledata->nuseful;
      int ncloudpoints = branchruledata->ncloudpoints;
      int nsavedlps = branchruledata->nsavedlps;

      SCIPstatisticMessage("time spent diving in cloud branching: %g\n", SCIPgetClockTime(scip, branchruledata->cloudclock));
      SCIPstatisticMessage("cloud branching tried: %6d      found cloud: %6d \n", ntried, nuseful);
      SCIPstatisticMessage("cloud used points: %6d      saved LPs: %6d \n", ncloudpoints, nsavedlps);
      SCIPstatisticMessage("cloud success rates useful/tried: %8.6g points/useful: %8.6g  saved/useful: %8.6g \n",
         ntried == 0 ? -1 : (SCIP_Real)nuseful / ntried,  nuseful == 0 ? -1 : (SCIP_Real)ncloudpoints / nuseful, nuseful == 0 ? -1 :  (SCIP_Real)nsavedlps / nuseful);
      SCIP_CALL( SCIPfreeClock(scip, &(branchruledata->cloudclock)) );
   }

   SCIPfreeMemoryArrayNull(scip, &branchruledata->skipdown);
   SCIPfreeMemoryArrayNull(scip, &branchruledata->skipup);

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
   branchruledata->nusefulunions = 0;
   branchruledata->ntried = 0;
   branchruledata->ntriedunions = 0;
   branchruledata->ncloudpoints = 0;
   branchruledata->nsavedlps = 0;

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
   SCIP_HASHTABLE* candtable;

   SCIP_VAR** lpcands;
   SCIP_VAR** lpcandscopy;
   SCIP_VAR** newlpcands;

   SCIP_VAR** vars;                          /* SCIP variables                */
   SCIP_ROW** lprows;
   SCIP_Real* lpcandsfrac;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfraccopy;
   SCIP_Real* lpcandssolcopy;
   SCIP_Real* lpcandsmin;
   SCIP_Real* lpcandsmax;
   SCIP_Real* newlpcandsmin;
   SCIP_Real* newlpcandsmax;

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
   int ncomplete;

   int newlpcandssize;
   int nnewlpcands;

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

   if( branchruledata->skipdown == NULL )
   {
      assert(branchruledata->skipup == NULL);

      SCIPallocMemoryArray(scip, &branchruledata->skipdown, SCIPgetNVars(scip));
      SCIPallocMemoryArray(scip, &branchruledata->skipup, SCIPgetNVars(scip));
   }

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, &npriolpcands) );
   nlpcands = SCIPgetNLPBranchCands(scip);
   assert(nlpcands > 0);

   /* reset skipping arrays to zero */
   BMSclearMemoryArray(branchruledata->skipdown, nlpcands);
   BMSclearMemoryArray(branchruledata->skipup, nlpcands);

   /* get problem variables and LP row data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   nlprows = SCIPgetNLPRows(scip);
   lprows = SCIPgetLPRows(scip);

   newlpcandssize = MAX(3*nlpcands, 100);
   nnewlpcands = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &lpcandsmin, nlpcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lpcandsmax, nlpcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lpcandscopy, nlpcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lpcandsfraccopy, nlpcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lpcandssolcopy, nlpcands) );
   if( branchruledata->useunion && SCIPgetDepth(scip) < branchruledata->maxdepthunion && !branchruledata->onlyF2)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &newlpcandsmin, newlpcandssize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newlpcandsmax, newlpcandssize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newlpcands, newlpcandssize) );
      SCIP_CALL( SCIPhashtableCreate(&candtable, SCIPblkmem(scip), SCIPcalcHashtableSize(newlpcandssize+nlpcands),
            hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );
   }
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

   /* store varibles from original fractional solution in hashtable; those are the only ones with a negative image */
   for( i = 0; i < nlpcands && branchruledata->useunion && !branchruledata->onlyF2 && SCIPgetDepth(scip) < branchruledata->maxdepthunion; ++i)
   {
      SCIP_CALL( SCIPhashtableInsert(candtable, lpcands[i]) );
   }

   newpoint = TRUE;
   counter = 0;

   /* loop that generates new cloud point */
   while( newpoint && branchruledata->usecloud )
   {
#ifdef NDEBUG
      SCIP_RETCODE retcode;
#endif

      SCIP_VAR** currbranchcands;

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

      if( lperror || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
         break;

      /* check if a solution has been found */
      if( SCIPgetNLPBranchCands(scip) == 0 )
      {
         SCIP_Bool success;
         SCIP_SOL* sol;

         /* create solution from diving LP */
         SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
         SCIP_CALL( SCIPlinkLPSol(scip, sol) );
         SCIPdebugMessage("cloud branching found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, sol));

         /* try to add solution to SCIP */
         SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, FALSE, FALSE, &success) );

         /* check, if solution was feasible and good enough */
         if( success )
         {
            SCIPdebugMessage(" -> solution was feasible and good enough\n");
            SCIP_CALL( SCIPendDive(scip) );
            *result = SCIP_CUTOFF;
            goto TERMINATE;
         }
      }

      /* update cloud intervals for candidates that have been fractional in original LP */
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

      if( branchruledata->useunion && !branchruledata->onlyF2 && SCIPgetDepth(scip) < branchruledata->maxdepthunion )
      {
         /* update cloud intervals for candidates that have been integral in original LP, but have been fractional in previous cloud points */
         for( i = 0; i < nnewlpcands; ++i)
         {
            SCIP_Real solval;
            solval = SCIPgetSolVal(scip, NULL, newlpcands[i]);
            newlpcandsmin[i] = MIN(newlpcandsmin[i], solval);
            newlpcandsmax[i] = MAX(newlpcandsmax[i], solval);
         }

         /* update cloud intervals for candidates that have been integral so far, but are fractional in current LP */
         SCIP_CALL( SCIPgetLPBranchCands(scip, &currbranchcands, NULL, NULL, NULL, NULL) );
         for( i = 0; i < SCIPgetNLPBranchCands(scip); ++i)
         {
            /* insert new variable in hashtable and newlpcands array */
            if( !SCIPhashtableExists(candtable, currbranchcands[i]) )
            {
               SCIP_Real solval;
               solval = SCIPgetSolVal(scip, NULL, currbranchcands[i]);
               assert(!SCIPisFeasIntegral(scip,solval));

               /* realloc array for new LP candidates, if necessary */
               if( nnewlpcands == newlpcandssize)
               {
                  newlpcandssize *= 3;
                  SCIP_CALL( SCIPreallocBufferArray(scip, &newlpcands, newlpcandssize) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &newlpcandsmin, newlpcandssize) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &newlpcandsmax, newlpcandssize) );
               }

               SCIP_CALL( SCIPhashtableInsert(candtable, currbranchcands[i]) );
               newlpcands[nnewlpcands] = currbranchcands[i];

               /* initialize the cloud interval; incorporation of the integral value in original LP solution will be done later */
               newlpcandsmin[nnewlpcands] = solval;
               newlpcandsmax[nnewlpcands] = solval;
               nnewlpcands++;
            }
         }
      }

      if( newpoint )
         counter++;

      if( branchruledata->maxpoints != -1 && counter >= branchruledata->maxpoints )
         break;
   }

   SCIPdebugMessage("considered %d additional points in the cloud\n",counter);

   if( nnewlpcands > 0 )
      SCIPdebugMessage("  --> found %d new LP cands\n",nnewlpcands);

   /* terminate the diving */
   SCIP_CALL( SCIPendDive(scip) );

   /* update cloud intervals for candidates that have been fractional in previous cloud points with their integral bound from original LP */
   for( i = 0; i < nnewlpcands && branchruledata->useunion && !branchruledata->onlyF2 && SCIPgetDepth(scip) < branchruledata->maxdepthunion ; ++i)
   {
      SCIP_Real solval;
      solval = SCIPgetSolVal(scip, NULL, newlpcands[i]);
      assert(SCIPisFeasIntegral(scip, solval));
      newlpcandsmin[i] = MIN(newlpcandsmin[i], solval);
      newlpcandsmax[i] = MAX(newlpcandsmax[i], solval);
   }

   SCIP_CALL(  SCIPstopClock(scip, branchruledata->cloudclock) );
   ncomplete = nlpcands;

   if( counter > 0 )
   {
      branchruledata->ncloudpoints += (counter+1);
      branchruledata->nuseful++;

      counter = 0;

      /* sort all variables for which both bounds are fractional to the front */
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

      /* should only be in that if condition when at least one bound could be made integral */
      assert(nlpcands - counter > 0);

      ncomplete = counter;

      /* filter all variables for which exactly one interval bound is fractional */
      for( i = 0; i < nlpcands && !branchruledata->onlyF2; ++i)
      {
         if( SCIPisFeasIntegral(scip, lpcandsmin[i]) != SCIPisFeasIntegral(scip, lpcandsmax[i]) )
         {
            assert(counter < nlpcands);
            lpcandscopy[counter] = lpcandscopy[i];
            lpcandssolcopy[counter] = lpcandssolcopy[i];
            lpcandsfraccopy[counter] = lpcandsfraccopy[i];

            if( SCIPisFeasIntegral(scip, lpcandsmin[i]) )
               branchruledata->skipdown[counter] = TRUE;
            if( SCIPisFeasIntegral(scip, lpcandsmax[i]) )
               branchruledata->skipup[counter] = TRUE;
            assert(branchruledata->skipdown[counter] != branchruledata->skipup[counter]);

            counter++;
         }
      }

      SCIPdebugMessage("can fully skip %d/%d strong branching candidates\n", nlpcands - counter, nlpcands);
      SCIPdebugMessage("can half  skip %d/%d strong branching candidates\n", counter - ncomplete, nlpcands);
   }
   else
      counter = nlpcands;

   /* if cloud sampling was not successful enough, disable it */
   if( branchruledata->usecloud &&
      branchruledata->ntried > 100 &&
      (SCIP_Real)branchruledata->nuseful / branchruledata->ntried < branchruledata->minsuccessrate )
   {
      SCIPdebugMessage("Disabling cloud branching (not effective)\n");
      branchruledata->usecloud = FALSE;
   }


   /* if union usage was not successful enough, disable it */
   if( branchruledata->useunion &&
      branchruledata->ntriedunions > 10 &&
      (SCIP_Real)branchruledata->nusefulunions / branchruledata->ntriedunions < branchruledata->minsuccessunion )
   {
      SCIPdebugMessage("Disabling union usage (not effective)\n");
      branchruledata->useunion = FALSE;
   }

   if( branchruledata->onlyF2 )
      counter = MAX(counter,1);

   SCIP_CALL( SCIPselectVarStrongBranching(scip, lpcandscopy, lpcandssolcopy, lpcandsfraccopy, branchruledata->skipdown, branchruledata->skipup, counter, counter, /* replace second counter ??????????? */
         ncomplete, &branchruledata->lastcand, allowaddcons,
         &bestcand, &bestdown, &bestup, &bestscore, &bestdownvalid, &bestupvalid, &provedbound, result) );

   if( branchruledata->lastcand <= ncomplete )
   {
      SCIPdebugMessage("saved %d of %d LPs\n", 2*(nlpcands - ncomplete), 2*nlpcands);
      branchruledata->nsavedlps += 2*(nlpcands - ncomplete);
   }
   else
   {
      SCIPdebugMessage("saved %d of %d LPs\n", 2*(nlpcands - counter)+counter - ncomplete, 2*nlpcands);
      branchruledata->nsavedlps += 2*(nlpcands - counter)+counter - ncomplete;
   }

   /* perform the branching */
   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM && *result != SCIP_CONSADDED && counter > 0 ) /* ??????? */
   {
      SCIP_NODE* downchild;
      SCIP_NODE* upchild;
      SCIP_VAR* var;
      SCIP_Bool allcolsinlp;
      SCIP_Bool exactsolve;
      SCIP_Bool newselected;

      assert(*result == SCIP_DIDNOTRUN);
      assert(0 <= bestcand && bestcand < nlpcands);
      assert(SCIPisLT(scip, provedbound, SCIPgetCutoffbound(scip)));

      var = lpcandscopy[bestcand];
      newselected = FALSE;

      /* if there are new candidates variables, also try them */
      if( branchruledata->useunion && !branchruledata->onlyF2 && SCIPgetDepth(scip) < branchruledata->maxdepthunion && branchruledata->lastcand > ncomplete && nnewlpcands > 0 )
      {
         counter = 0;
         /* reset skipping arrays to zero */
         BMSclearMemoryArray(branchruledata->skipdown, nnewlpcands);
         BMSclearMemoryArray(branchruledata->skipup, nnewlpcands);

         /* get new LP candidates with one fractional bound */
         for( i = 0; i < nnewlpcands; ++i)
         {
            if( SCIPisFeasIntegral(scip, newlpcandsmin[i]) != SCIPisFeasIntegral(scip, newlpcandsmax[i]) )
            {
               newlpcands[counter] = newlpcands[i];

               if( SCIPisFeasIntegral(scip, newlpcandsmin[i]) )
                  branchruledata->skipdown[counter] = TRUE;
               if( SCIPisFeasIntegral(scip, newlpcandsmax[i]) )
                  branchruledata->skipup[counter] = TRUE;
               assert(branchruledata->skipdown[counter] != branchruledata->skipup[counter]);

               counter++;
            }
         }

         /* when there are new candidates, also try these */
         if( counter > 0 )
         {
            SCIP_Real newdown;
            SCIP_Real newup;
            SCIP_Real newscore;
            int newcand;
            SCIP_Bool newdownvalid;
            SCIP_Bool newupvalid;
            SCIP_Real newbound;
            int lastcand;

            branchruledata->ntriedunions++;
            lastcand = 0;
            newscore = -SCIPinfinity(scip);
            SCIP_CALL( SCIPselectVarPseudoStrongBranching(scip, newlpcands, branchruledata->skipdown, branchruledata->skipup, counter, counter, /* replace second counter ??????????? */
                  &lastcand, allowaddcons,
                  &newcand, &newdown, &newup, &newscore, &newdownvalid, &newupvalid, &newbound, result) );

            if( *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM || *result == SCIP_CONSADDED  ) /* ??????? */
               goto TERMINATE;

            if( newscore > bestscore )
            {
               bestcand = newcand;
               var = newlpcands[newcand];
               bestdown = newdown;
               bestup = newup;
               bestdownvalid = newdownvalid;
               bestupvalid = newupvalid;
               bestscore = newscore;
               newselected = TRUE;
               branchruledata->nusefulunions++;
            }
         }
      }

      /* perform the branching */
      if( !newselected )
      {
         SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%g, down=%g, up=%g, score=%g)\n",
            counter, bestcand, SCIPvarGetName(var), lpcandssolcopy[bestcand], bestdown, bestup, bestscore);
      }
      else
      {
         SCIPdebugMessage(" -> selected from %d new candidates,  candidate %d: variable <%s> (down=%g, up=%g, score=%g)\n",
            counter, bestcand, SCIPvarGetName(var), bestdown, bestup, bestscore);
      }

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

 TERMINATE:
   if( branchruledata->useunion && !branchruledata->onlyF2 && SCIPgetDepth(scip) < branchruledata->maxdepthunion )
   {
      SCIPhashtableFree(&candtable);
      SCIPfreeBufferArray(scip, &newlpcands);
      SCIPfreeBufferArray(scip, &newlpcandsmax);
      SCIPfreeBufferArray(scip, &newlpcandsmin);
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
   branchruledata->skipup = NULL;
   branchruledata->skipdown = NULL;

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
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/"BRANCHRULE_NAME"/onlyF2",
         "should only F2 be used? ",
         &branchruledata->onlyF2, FALSE, DEFAULT_ONLYF2, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/"BRANCHRULE_NAME"/useunion",
         "should the union of candidates be used?",
         &branchruledata->useunion, FALSE, DEFAULT_USEUNION, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/"BRANCHRULE_NAME"/maxpoints",
         "maximum number of points for the cloud (-1 means no limit)",
         &branchruledata->maxpoints, FALSE, DEFAULT_MAXPOINTS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/"BRANCHRULE_NAME"/minsuccessrate",
         "minimum success rate for the cloud",
         &branchruledata->minsuccessrate, FALSE, DEFAULT_MINSUCCESSRATE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/"BRANCHRULE_NAME"/minsuccessunion",
         "minimum success rate for the union",
         &branchruledata->minsuccessunion, FALSE, DEFAULT_MINSUCCESSUNION, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/"BRANCHRULE_NAME"/maxdepthunion",
         "maximum depth for the union",
         &branchruledata->maxdepthunion, FALSE, DEFAULT_MAXDEPTHUNION, 0, 65000, NULL, NULL) );

   return SCIP_OKAY;
}
