/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_solvingphase.c
 * @brief  event handler for solving phase dependent parameter adjustment
 * @author Gregor Hendel
 *
 * this event handler provides methods to support parameter adjustment at every new of the three solving phases:
 *   - Feasibility phase - before the first solution is found
 *   - Improvement phase - after the first solution was found until an optimal solution is found or believed to be found
 *   - Proof phase - the remaining time of the solution process after an optimal or believed-to-be optimal incumbent has been found.
 *
 * Of course, this event handler cannot detect by itself whether a given incumbent is optimal prior to termination of the
 * solution process. It rather uses heuristic transitions based on properties of the search tree in order to
 * determine the appropriate stage. Settings files can be passed to this event handler for each of the three phases.
 *
 * This approach of phase-based parameter adjustment was first presented in
 *
 * Gregor Hendel
 * Empirical Analysis of Solving Phases in Mixed-Integer Programming
 * Master thesis, Technical University Berlin (2014)
 *
 * with the main results also available from
 *
 * Gregor Hendel
 * Exploiting solving phases in mixed-integer programs (2015)
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_solvingphase.h"
#include "string.h"
#include "scip/event_treeinfos.h"
#include "scip/event_logregression.h"

#define EVENTHDLR_NAME         "solvingphase"
#define EVENTHDLR_DESC         "event handler to adjust settings depending on current stage"

#define MAXLINELEN 1024           /**< limit for line length */
#define DEFAULT_SOLUFILENAME "myshort.solu"/**< default filename to read solution value from */
#define DEFAULT_SOLUFILEPATH "//nfs//optimi//kombadon//bzfhende//projects//scip-git//check//testset//" /**< default solution file path */
#define DEFAULT_SETTINGFILEPATH "/nfs/optimi/kombadon/bzfhende/projects/scip-git/settings/%s"/**< default settings file path */
#define DEFAULT_SETNAME "default.set" /**< default settings file name */
#define EVENTHDLR_EVENT SCIP_EVENTTYPE_BESTSOLFOUND | SCIP_EVENTTYPE_NODESOLVED /**< the actual event to be caught */
#define TRANSITIONMETHODS "elor" /**< which heuristic transition method: (e)stimate based, (l)ogarithmic regression based, (o)ptimal value based (cheat!),
                                      (r)ank-1 node based? */
#define DEFAULT_TRANSITIONMETHOD 'r' /**< the default transition method */
#define DEFAULT_NODEOFFSET 50          /**< default node offset before transition to proof phase is active */
#define DEFAULT_FALLBACK FALSE         /**< should the phase transition fall back to suboptimal phase? */
#define DEFAULT_INTERRUPTOPTIMAL FALSE /**< should solving process be interrupted if optimal solution was found? */
#define DEFAULT_ADJUSTRELPSWEIGHTS FALSE /**< should the scoring weights of the hybrid reliability pseudo cost branching rule be adjusted? */
#define DEFAULT_USEFILEWEIGHTS   FALSE /**< should weights from a weight file be used to adjust branching score weights? */
#define DEFAULT_USEWEIGHTEDQUOTIENTS TRUE /**< should weighted quotients be used to adjust branching score weights? */

#define DEFAULT_ENABLED FALSE     /**< should the event handler be executed? */
#define DEFAULT_TESTMODE FALSE    /**< should the event handler test the criteria? */
/*
 * Data structures
 */
/** enumerator to represent the current solving phase */
enum SolvingPhase
{
   SOLVINGPHASE_UNINITIALIZED = -1,          /**< solving phase has not been initialized yet */
   SOLVINGPHASE_FEASIBILITY   = 0,           /**< no solution was found until now */
   SOLVINGPHASE_IMPROVEMENT   = 1,           /**< current incumbent solution is suboptimal */
   SOLVINGPHASE_PROOF         = 2            /**< current incumbent is optimal */
};
typedef enum SolvingPhase SOLVINGPHASE;

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_Bool            enabled;             /**< should the event handler be executed? */
   char*                solufilename;        /**< file to parse solution information from */
   char*                feassetparam;        /**< settings file parameter for the feasibility phase */
   char*                improvesetparam;     /**< settings file parameter for the improvement phase */
   char*                proofsetparam;       /**< settings file parameter for the proof phase */
   SCIP_Real            optimalvalue;        /**< value of optimal solution of the problem */
   SCIP_Real            lambda12;            /**< gap to reach for phase transition feasibility phase -> improvement phase */
   SCIP_Real            lambda23;            /**< gap to reach for phase transition improvement phase -> proof phase */
   SOLVINGPHASE         solvingphase;        /**< the current solving phase */
   char                 transitionmethod;    /**< how to do the transition from phase2 -> phase3? (c)orrected estimate based,
                                                  (e)stimate based, (l)ogarithmic regression based, (o)ptimal value based (cheat!),
                                                  (r)ank-1 node based */
   SCIP_Longint         nodeoffset;           /**< node offset for triggering rank-1 node based phased transition */
   SCIP_Bool            fallback;             /**< should the phase transition fall back to improvement phase? */
   SCIP_Bool            interruptoptimal;     /**< interrupt after optimal solution was found */
   SCIP_Bool            adjustrelpsweights;   /**< should the relpscost cutoff weights be adjusted? */
   SCIP_Bool            useweightedquotients; /**< should weighted quotients between infeasible and pruned leaf nodes be considered? */
   SCIP_EVENTHDLR*      linregeventhdlr;      /**< pointer to the linear regression event handler */
   SCIP_Bool            testmode;             /**< should transitions be tested only, but not triggered? */
   SCIP_Bool            rank1reached;         /**< has the rank-1 transition into proof phase been reached? */
   SCIP_Bool            estimatereached;      /**< has the best-estimate transition been reached? */
   SCIP_Bool            optimalreached;       /**< is the incumbent already optimal? */
   SCIP_Bool            logreached;           /**< has a logarithmic phase transition been reached? */
};

/*
 * Local methods
 */

/** returns the optimal value for this instance (as passed to the event handler) */
SCIP_Real SCIPgetOptimalSolutionValue(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   return eventhdlrdata->optimalvalue;
}

/** checks if rank-1 transition has been reached, that is, when all open nodes have a best-estimate higher than the best
 *  previously checked node at this depth
 */
static
SCIP_Bool checkRankOneTransition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   /* at least one solution is required for the transition */
   if( SCIPgetNSols(scip) > 0 )
      return (SCIPgetNNodes(scip) > eventhdlrdata->nodeoffset && SCIPgetNRank1Nodes(scip) == 0);
   else
      return FALSE;
}

/** check if Best-Estimate criterion was reached, that is, when the active estimate is not better than the current incumbent solution */
static
SCIP_Bool checkEstimateCriterion(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*   eventhdlrdata
   )
{
   if( SCIPgetNSols(scip) > 0 )
      return SCIPgetNNodes(scip) > eventhdlrdata->nodeoffset && SCIPgetNNodesBelowIncumbent(scip) == 0;
   else
      return FALSE;
}

/** check if logarithmic phase transition has been reached */
static
SCIP_Bool checkLogCriterion(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*   eventhdlrdata
   )
{
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real axisintercept = SCIPgetCurrentTangentAxisIntercept(scip, eventhdlrdata->linregeventhdlr);
      if( !SCIPisInfinity(scip, axisintercept) )
      {
         SCIP_Real primalbound;
         SCIP_Real lambda;
         SCIP_Real firstprimalbound = SCIPgetFirstPrimalBound(scip);
         primalbound = SCIPgetPrimalbound(scip);
         lambda = (axisintercept - primalbound) / (firstprimalbound - primalbound);
         if( SCIPisNegative(scip, lambda) )
            return TRUE;
      }
   }
   return FALSE;
}

/** check if incumbent solution is optimal with respect to a slightly modified gap definition */
static
SCIP_Bool checkOptimalSolution(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*   eventhdlrdata
   )
{
   SCIP_Real referencevalue;
   SCIP_Real primalbound;

   referencevalue = eventhdlrdata->optimalvalue;
   primalbound = SCIPgetPrimalbound(scip);
   if(!SCIPisInfinity(scip, REALABS(primalbound)) && !SCIPisInfinity(scip, referencevalue) )
   {
      SCIP_Real max = MAX3(1.0, REALABS(primalbound), REALABS(referencevalue));
      if( EPSZ((primalbound - referencevalue)/max, 1e-9) )
         return TRUE;
   }
   return FALSE;
}

/** check if we are in the proof phase */
static
SCIP_Bool transitionPhase3(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata
   )
{
   if( eventhdlrdata->solvingphase == SOLVINGPHASE_PROOF && !eventhdlrdata->fallback )
      return TRUE;

   /* check criterion based on selected transition method */
   switch( eventhdlrdata->transitionmethod )
   {
      case 'r':

         /* check rank-1 transition */
         if( checkRankOneTransition(scip, eventhdlrdata) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "reached rank-1 transition: nodes: %lld, rank-1: %d bound: %9.5g time: %.2f\n",
                  SCIPgetNNodes(scip), SCIPgetNRank1Nodes(scip), SCIPgetPrimalbound(scip), SCIPgetSolvingTime(scip));
            return TRUE;
         }
         break;
      case 'o':

         /* cheat and use knowledge about optimal solution */
         if( checkOptimalSolution(scip, eventhdlrdata) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "optimal solution found: %lld, bound: %9.5g time: %.2f\n",
                  SCIPgetNNodes(scip), SCIPgetPrimalbound(scip), SCIPgetSolvingTime(scip));
            return TRUE;
         }
         break;
      case 'e':

         /* check best-estimate transition */
         if( checkEstimateCriterion(scip, eventhdlrdata) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "reached best-estimate transition: nodes: %lld, estimate: %d bound: %9.5g time: %.2f\n",
                  SCIPgetNNodes(scip), SCIPgetNNodesBelowIncumbent(scip), SCIPgetPrimalbound(scip), SCIPgetSolvingTime(scip));
            return TRUE;
         }
         return FALSE;
         break;
      case 'l':

         /* check logarithmic transition */
         if( checkLogCriterion(scip, eventhdlrdata) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "reached a logarithmic phase transition: %.2f", SCIPgetSolvingTime(scip));
            return TRUE;
         }
         break;
      default:
         return FALSE;
      break;
   }

   return FALSE;
}

/* determine the solving phase: feasibility phase if no solution was found yet, otherwise improvement phase or proof phase
 * depending on whether selected transition criterion was already reached */
static
void determineSolvingPhase(
   SCIP* scip,
   SCIP_EVENTHDLRDATA* eventhdlrdata
   )
{
   /* without solution, we are in the feasibility phase */
   if( SCIPgetNSols(scip) == 0 )
      eventhdlrdata->solvingphase = SOLVINGPHASE_FEASIBILITY;
   else if( eventhdlrdata->solvingphase != SOLVINGPHASE_PROOF || eventhdlrdata->fallback )
      eventhdlrdata->solvingphase = SOLVINGPHASE_IMPROVEMENT;

   if( eventhdlrdata->solvingphase == SOLVINGPHASE_IMPROVEMENT && transitionPhase3(scip, eventhdlrdata) )
      eventhdlrdata->solvingphase = SOLVINGPHASE_PROOF;
}

/* apply the phase based settings: A phase transition invokes completely new parameters */
static
SCIP_RETCODE applySolvingPhase(
   SCIP* scip,
   SCIP_EVENTHDLRDATA* eventhdlrdata
   )
{
   FILE* file;
   SOLVINGPHASE phasebefore;
   char paramfilename[256];

   if( eventhdlrdata->solvingphase == SOLVINGPHASE_PROOF && !eventhdlrdata->fallback )
      return SCIP_OKAY;

   phasebefore = eventhdlrdata->solvingphase;
   determineSolvingPhase(scip, eventhdlrdata);

   if( eventhdlrdata->solvingphase == SOLVINGPHASE_PROOF && eventhdlrdata->transitionmethod == 'o' &&
         eventhdlrdata->interruptoptimal )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Solution is optimal. Calling user interruption\n");
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   if( phasebefore == eventhdlrdata->solvingphase )
      return SCIP_OKAY;


   switch (eventhdlrdata->solvingphase)
   {
   case SOLVINGPHASE_FEASIBILITY:
      sprintf(paramfilename, DEFAULT_SETTINGFILEPATH, eventhdlrdata->nosolsetname);
      break;
   case SOLVINGPHASE_IMPROVEMENT:
      sprintf(paramfilename, DEFAULT_SETTINGFILEPATH, eventhdlrdata->suboptsetname);
      break;
   case SOLVINGPHASE_PROOF:
      sprintf(paramfilename, DEFAULT_SETTINGFILEPATH, eventhdlrdata->optsetname);
      break;
   default:
      SCIPdebugMessage("Unknown solving phase: %d -> ABORT!\n ", eventhdlrdata->solvingphase);
      SCIPABORT();
      break;
   }
   assert(paramfilename != NULL);
   file = fopen(paramfilename, "r");


   if( file == NULL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,"Changed solving phase to %d \n", eventhdlrdata->solvingphase);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,"Parameter file %s not found--keeping settings as before \n", paramfilename);
   }
   else
   {
      char transitionmethod;
      SCIP_Bool interruptoptimal;

      fclose(file);
      SCIPduplicateMemoryArray(scip, &eventhdlrdata->nosolsetname, eventhdlrdata->feassetparam, strlen(eventhdlrdata->feassetparam) + 1);
      SCIPduplicateMemoryArray(scip, &eventhdlrdata->suboptsetname, eventhdlrdata->improvesetparam, strlen(eventhdlrdata->improvesetparam) + 1);
      SCIPduplicateMemoryArray(scip, &eventhdlrdata->optsetname, eventhdlrdata->proofsetparam, strlen(eventhdlrdata->proofsetparam) + 1);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,"Changed solving phase to %d -- \n", eventhdlrdata->solvingphase);
      //SCIP_CALL( SCIPresetParams(scip) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Reading parameters from file %s\n", paramfilename);
      interruptoptimal = eventhdlrdata->interruptoptimal;
      transitionmethod = eventhdlrdata->transitionmethod;

      SCIP_CALL( SCIPreadParams(scip, paramfilename) );

      eventhdlrdata->enabled = TRUE;
      eventhdlrdata->transitionmethod = transitionmethod;
      eventhdlrdata->interruptoptimal = interruptoptimal;
   }

   if( eventhdlrdata->solvingphase == SOLVINGPHASE_PROOF && eventhdlrdata->adjustrelpsweights && SCIPgetPhase(scip) == SCIP_STAGE_SOLVING )
   {
      SCIP_Longint objleaves;
      SCIP_Real quotient;
      SCIP_Real newcutoffweight;
      SCIP_Real newconflictweight;
      SCIP_Real cutoffweight;
      SCIP_Real conflictweight;
      SCIP_Longint cutoffleaves;

      objleaves = SCIPgetNObjLeaves(scip);
      cutoffleaves = SCIPgetNInfeasLeaves(scip);
      objleaves = MAX(objleaves, 1);
      cutoffleaves = MAX(cutoffleaves, 1);

      quotient = cutoffleaves / (SCIP_Real)objleaves;

      newcutoffweight = quotient;
      newconflictweight = quotient;

      SCIP_CALL( SCIPgetRealParam(scip, "branching/relpscost/conflictweight", &conflictweight) );
      SCIP_CALL( SCIPgetRealParam(scip, "branching/relpscost/cutoffweight", &cutoffweight) );

      if( eventhdlrdata->useweightedquotients )
      {
         newcutoffweight *= cutoffweight;
         newconflictweight *= conflictweight;
      }

      SCIP_CALL( SCIPsetRealParam(scip, "branching/relpscost/conflictweight", newconflictweight) );
      SCIP_CALL( SCIPsetRealParam(scip, "branching/relpscost/cutoffweight", newcutoffweight) );

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
            "  Adjusting relpscost weights, (quot = %.3g): cutoffweight %.4g --> %.4g, confweight: %.4g --> %.4g \n", quotient,
            cutoffweight, newcutoffweight, conflictweight, newconflictweight);
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of event handler
 */

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_EVENTCOPY(eventCopySolvingphase)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* call inclusion method of node selector */
   SCIP_CALL( SCIPincludeEventHdlrSolvingphase(scip) );

   return SCIP_OKAY;
}
/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeSolvingphase)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   dataFreeStrings(scip, eventhdlrdata);
   SCIPfreeMemory(scip, &eventhdlrdata);
   eventhdlrdata = NULL;

   return SCIP_OKAY;

}


/** initialization method of event handler (called after problem was transformed) */

static
SCIP_DECL_EVENTINIT(eventInitSolvingphase)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   const char* probname;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   /* read solufile information for problem */

   /* free memory from previous run first (can happen if several problems are consecutively read in one session) */
   dataFreeStrings(scip, eventhdlrdata);

   eventhdlrdata->solvingphase = SOLVINGPHASE_UNINITIALIZED;


   probname = SCIPstringGetBasename(SCIPgetProbName(scip));
   eventhdlrdata->optimalvalue = SCIPinfinity(scip);
   SCIP_CALL( searchSolufileForProbname(scip, probname, eventhdlrdata) );
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Optimal value for problem %s from solufile %s: %16.9g\n", probname, eventhdlrdata->solufilename, eventhdlrdata->optimalvalue);

   if( eventhdlrdata->enabled )
   {
      SCIP_CALL( applySolvingPhase(scip, eventhdlrdata) );
   }

   if( eventhdlrdata->enabled || eventhdlrdata->testmode )
   {
      SCIP_CALL( SCIPcatchEvent(scip, EVENTHDLR_EVENT, eventhdlr, NULL, NULL) );
   }
   eventhdlrdata->optimalreached = FALSE;
   eventhdlrdata->logreached = FALSE;
   eventhdlrdata->rank1reached = FALSE;
   eventhdlrdata->estimatereached = FALSE;

   eventhdlrdata->linregeventhdlr = SCIPfindEventhdlr(scip, "logregression");
   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolSolvingphase)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecSolvingphase)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   assert(SCIPeventGetType(event) & (EVENTHDLR_EVENT));

   if( eventhdlrdata->enabled )
   {
      SCIP_CALL( applySolvingPhase(scip, eventhdlrdata) );
   }

   if( eventhdlrdata->testmode )
   {
      if( !eventhdlrdata->rank1reached && checkRankOneTransition(scip, eventhdlrdata) )
      {
         eventhdlrdata->rank1reached = TRUE;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Rank 1 criterion reached after %lld nodes, %.2f sec.\n", SCIPgetNNodes(scip), SCIPgetSolvingTime(scip));
      }

      if( !eventhdlrdata->estimatereached && checkEstimateCriterion(scip, eventhdlrdata) )
      {
         eventhdlrdata->estimatereached = TRUE;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Estimate criterion reached after %lld nodes, %.2f sec.\n", SCIPgetNNodes(scip), SCIPgetSolvingTime(scip));
      }

      if( !eventhdlrdata->optimalreached && checkOptimalSolution(scip, eventhdlrdata) )
      {
         eventhdlrdata->optimalreached = TRUE;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Optimum reached after %lld nodes, %.2f sec.\n", SCIPgetNNodes(scip), SCIPgetSolvingTime(scip));
      }

      if( !eventhdlrdata->logreached && checkLogCriterion(scip, eventhdlrdata) )
      {
         eventhdlrdata->logreached = TRUE;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  Log criterion reached after %lld nodes, %.2f sec.\n", SCIPgetNNodes(scip), SCIPgetSolvingTime(scip));
      }
   }

   return SCIP_OKAY;
}

/** creates event handler for Solvingphase event */
SCIP_RETCODE SCIPincludeEventHdlrSolvingphase(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create Solvingphase event handler data */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   assert(eventhdlrdata != NULL);

   eventhdlrdata->solufilename = NULL;
   eventhdlrdata->feassetparam = NULL;
   eventhdlrdata->improvesetparam = NULL;
   eventhdlrdata->proofsetparam = NULL;

   eventhdlrdata->linregeventhdlr = NULL;

   eventhdlr = NULL;


   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecSolvingphase, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopySolvingphase) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeSolvingphase) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitSolvingphase) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolSolvingphase) );

   /* add Solvingphase event handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/"EVENTHDLR_NAME"/enabled", "should the event handler be executed?",
         &eventhdlrdata->enabled, FALSE, DEFAULT_ENABLED, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/"EVENTHDLR_NAME"/testmode", "should the event handler test for phase transition?",
         &eventhdlrdata->testmode, FALSE, DEFAULT_TESTMODE, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "eventhdlr/"EVENTHDLR_NAME"/solufilename","file to parse solution information from",
         &eventhdlrdata->solufilename, FALSE, DEFAULT_SOLUFILENAME, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "eventhdlr/"EVENTHDLR_NAME"/nosolsetname", "bla",
            &eventhdlrdata->feassetparam, FALSE, DEFAULT_SETNAME, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "eventhdlr/"EVENTHDLR_NAME"/suboptsetname", "settings file for suboptimal solving phase",
               &eventhdlrdata->improvesetparam, FALSE, DEFAULT_SETNAME, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "eventhdlr/"EVENTHDLR_NAME"/optsetname", "settings file for optimal solving phase",
               &eventhdlrdata->proofsetparam, FALSE, DEFAULT_SETNAME, NULL, NULL) );



   SCIP_CALL( SCIPaddLongintParam(scip, "eventhdlr/"EVENTHDLR_NAME"/nodeoffset", "node offset", &eventhdlrdata->nodeoffset,
         FALSE, DEFAULT_NODEOFFSET, 1, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/"EVENTHDLR_NAME"/fallback", "should the event handler fall back from optimal phase?",
            &eventhdlrdata->fallback, FALSE, DEFAULT_FALLBACK, NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip ,"eventhdlr/"EVENTHDLR_NAME"/transitionmethod", "transition method 'e','l','o','r'",
         &eventhdlrdata->transitionmethod, FALSE, DEFAULT_TRANSITIONMETHOD, TRANSITIONMETHODS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/"EVENTHDLR_NAME"/interruptoptimal", "should the event handler interrupt after optimal solution was found?",
               &eventhdlrdata->interruptoptimal, FALSE, DEFAULT_INTERRUPTOPTIMAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/"EVENTHDLR_NAME"/adjustrelpsweights", "should the branching score weights for cutoffs and "
               "conflicts be adjusted after optimal solution was found?", &eventhdlrdata->adjustrelpsweights, FALSE,
               DEFAULT_ADJUSTRELPSWEIGHTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/"EVENTHDLR_NAME"/useweightedquotients", "use weighted quotients?", &eventhdlrdata->useweightedquotients,
           FALSE, DEFAULT_USEWEIGHTEDQUOTIENTS, NULL, NULL) );

   return SCIP_OKAY;
}
