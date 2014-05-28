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

/**@file   event_solvingstage.c
 * @brief  eventhdlr for solving stage dependent parameter adjustment
 * @author Gregor Hendel
 *
 * this event handler is used to apply dynamic parameter adjustment depending on the
 * progress of the solving process.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_solvingstage.h"
#include "string.h"
#include "scip/event_treeinfos.h"

#define EVENTHDLR_NAME         "solvingstage"
#define EVENTHDLR_DESC         "event handler for events which influence the solving stage"

/*
 * Data structures
 */
#define MAXLINELEN 1024   /**< limit for line length */
#define DEFAULT_ENABLED FALSE /**< should the event handler be executed? */
#define DEFAULT_SOLUFILENAME "myshort.solu" /**< default filename to read solution value from */
#define DEFAULT_SOLUFILEPATH "//nfs//optimi//kombadon//bzfhende//projects//scip-git//check//testset//"
#define DEFAULT_SETTINGFILEPATH "/nfs/optimi/kombadon/bzfhende/projects/scip-git/settings/%s"
#define DEFAULT_SETNAME "default.set"
#define EVENTHDLR_EVENT SCIP_EVENTTYPE_BESTSOLFOUND | SCIP_EVENTTYPE_NODESOLVED /**< the actual event to be caught */
#define DEFAULT_LAMBA12 0.6 /* gap to reach for phase transition phase 1 -> phase 2 */
#define DEFAULT_LAMBA23 0.03 /* gap to reach for phase transition phase 2 -> phase 3 */
#define MAXGAP 1
#define TRANSITIONMETHODS "celor" /**< how to do the transition from phase2 -> phase3? (c)orrected estimate based,
                                      (e)stimate beased, (l)ogarithmic regression based, (o)ptimal value based (cheat!),
                                      (r)ank1 node based */
#define DEFAULT_TRANSITIONMETHOD 'r' /**< the default transition method */
#define DEFAULT_NODEOFFSET 50        /**< default node offset */
#define DEFAULT_FALLBACK FALSE       /**< should the phase transition fall back to suboptimal stage? */

/** enumerator to represent the event handler solving stage */
enum SolvingStage
{
   SOLVINGSTAGE_UNINITIALIZED = -1,         /**< solving stage has not been initialized yet */
   SOLVINGSTAGE_NOSOLUTION = 0,             /**< no solution was found until now */
   SOLVINGSTAGE_SUBOPTIMAL = 1,             /**< current incumbent solution is suboptimal */
   SOLVINGSTAGE_OPTIMAL    = 2,             /**< current incumbent is optimal */
   SOLVINGSTAGE_INFEASIBLE = 3              /**< the problem is allegedly infeasible */
};
typedef enum SolvingStage SOLVINGSTAGE;

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_Bool            enabled;             /**< should the event handler be executed? */
   char*                solufilename;        /**< file to parse solution information from */
   char*                nosolsetnameparam;
   char*                suboptsetnameparam;
   char*                optsetnameparam;
   char*                infeasiblesetnameparam;
   char*                nosolsetname;
   char*                suboptsetname;
   char*                optsetname;
   char*                infeasiblesetname;
   SCIP_Real            optimalvalue;        /**< value of optimal solution of the problem */
   SCIP_Real            lambda12;            /**< gap to reach for phase transition phase 1 -> phase 2 */
   SCIP_Real            lambda23;            /**< gap to reach for phase transition phase 2 -> phase 3 */
   SOLVINGSTAGE         solvingstage;        /**< the current solving stage */
   char                 transitionmethod;    /**< how to do the transition from phase2 -> phase3? (c)orrected estimate based,
                                                  (e)stimate beased, (l)ogarithmic regression based, (o)ptimal value based (cheat!),
                                                  (r)ank1 node based */
   SCIP_Longint         nodeoffset;          /**< node offset for triggering rank1 node based phased transition */
   SCIP_Bool            fallback;             /**< should the phase transition fall back to suboptimal stage? */
};

/*
 * Local methods
 */

SCIP_Real SCIPgetOptimalSolutionValue(
   SCIP*                 scip
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   return eventhdlrdata->optimalvalue;
}
static
SCIP_RETCODE searchSolufileForProbname(
   SCIP*                scip,
   const char*          probname,
   SCIP_EVENTHDLRDATA*  eventhdlrdata
   )
{
   SCIP_FILE* file;
   char linebuf[MAXLINELEN];
   char solufilename[256];
   sprintf(solufilename, "/optimi/kombadon/bzfhende/projects/scip-git/check/testset/%s", eventhdlrdata->solufilename);
   file = NULL;
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Trying to open %s\n", solufilename);
   file = SCIPfopen(solufilename, "r");

   if( file == NULL )
   {
      return SCIP_NOFILE;
   }
   while( SCIPfgets(linebuf, (int)sizeof(linebuf), file) != NULL )
   {
      char* status;
      char* name;
      char* obj;

      status = strtok(linebuf, " ");
      name = strtok(NULL, " ");

      if( strcasecmp(name, probname) != 0 )
         continue;

      if( strcmp(status, "=opt=") != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "No optimal solution value known for this problem!\n");
         break;
      }
      else
      {
         char* endofptr;
         double val;
         obj = strtok(NULL, " ");
         val = strtod(obj, &endofptr);
         eventhdlrdata->optimalvalue = val;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Optimal value for problem: %16.9g\n", eventhdlrdata->optimalvalue);
         break;
      }
   }
   SCIPfclose(file);

   return SCIP_OKAY;
}

/* gap function to compare values, returns a value between 0.0 and MAXGAP */
static
SCIP_Real getGap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first for gap */
   SCIP_Real             val2                /**< second value for gap */
   )
{
   assert(scip != NULL);

   /* if one value is infinite return the maximum gap */
   if( SCIPisInfinity(scip, REALABS(val1)) || SCIPisInfinity(scip, REALABS(val2)) )
      return MAXGAP;
   /* feasibly equal values have zero gap */
   if( SCIPisFeasEQ(scip, val1,val2) )
      return 0.0;
   else
   {
      /* calculate the gap as |val1 - val2|/max(|val1|,|val2|), and bound it by MAXGAP (if val1, val2 have different signs)*/
      SCIP_Real absdiff;
      SCIP_Real maxabs;
      SCIP_Real gap;
      absdiff = REALABS(val1 - val2);
      maxabs = MAX(REALABS(val1), REALABS(val2));

      gap = absdiff / maxabs;

      return MIN(gap, MAXGAP);
   }
}

static
SCIP_Bool transitionPhase3(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata
   )
{
   SCIP_Real primalbound;
   SCIP_Real referencevalue;

   if( eventhdlrdata->solvingstage == SOLVINGSTAGE_OPTIMAL && !eventhdlrdata->fallback )
      return TRUE;

   switch( eventhdlrdata->transitionmethod )
   {
      case 'r':
         if( SCIPgetNNodes(scip) > eventhdlrdata->nodeoffset && SCIPgetNRank1Nodes(scip) == 0 )
         {

            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "triggering a rank 1 transition: nodes: %lld, rank1: %d bound: %9.5g\n",
                  SCIPgetNNodes(scip), SCIPgetNRank1Nodes(scip), SCIPgetPrimalbound(scip));
            return TRUE;
         }
         break;
      case 'o':
         referencevalue = eventhdlrdata->optimalvalue;
         primalbound = SCIPgetPrimalbound(scip);
         if(!SCIPisInfinity(scip, REALABS(primalbound)) && !SCIPisInfinity(scip, referencevalue) )
            return SCIPisFeasZero(scip, getGap(scip, primalbound, referencevalue));
         break;
      case 'c':
         return FALSE;
      case 'e':
         return FALSE;
      case 'l':
         return FALSE;
      default:
         return FALSE;
      break;
   }

   return FALSE;
}

/* determine the solving phase -> phase 1: gap > lambda12, phase2: lambda12 > gap > lambda23, phase 3: gap <= lambda23 */
static
void determineSolvingStage(
   SCIP* scip,
   SCIP_EVENTHDLRDATA* eventhdlrdata
   )
{
   /* use the optimal value as reference value unless it is not available in which case we take the dual bound instead */
   if( SCIPgetNSols(scip) == 0 )
      eventhdlrdata->solvingstage = SOLVINGSTAGE_NOSOLUTION;
   else if( eventhdlrdata->solvingstage != SOLVINGSTAGE_OPTIMAL || eventhdlrdata->fallback )
      eventhdlrdata->solvingstage = SOLVINGSTAGE_SUBOPTIMAL;

   if( eventhdlrdata->solvingstage == SOLVINGSTAGE_SUBOPTIMAL && transitionPhase3(scip, eventhdlrdata) )
      eventhdlrdata->solvingstage = SOLVINGSTAGE_OPTIMAL;
}

/* apply the phase based settings: A phase transition invokes completely new parameters */
static
SCIP_RETCODE applySolvingStage(
   SCIP* scip,
   SCIP_EVENTHDLRDATA* eventhdlrdata
   )
{
   FILE* file;
   SOLVINGSTAGE stagebefore;
   char paramfilename[256];

   stagebefore = eventhdlrdata->solvingstage;
   determineSolvingStage(scip, eventhdlrdata);

   if( stagebefore != SOLVINGSTAGE_NOSOLUTION && stagebefore == eventhdlrdata->solvingstage )
      return SCIP_OKAY;

   if( eventhdlrdata->solvingstage == SOLVINGSTAGE_OPTIMAL && !eventhdlrdata->fallback )
      return SCIP_OKAY;

   switch (eventhdlrdata->solvingstage)
   {
   case SOLVINGSTAGE_NOSOLUTION:
      sprintf(paramfilename, DEFAULT_SETTINGFILEPATH, eventhdlrdata->nosolsetname);
      break;
   case SOLVINGSTAGE_SUBOPTIMAL:
      sprintf(paramfilename, DEFAULT_SETTINGFILEPATH, eventhdlrdata->suboptsetname);
      break;
   case SOLVINGSTAGE_OPTIMAL:
      sprintf(paramfilename, DEFAULT_SETTINGFILEPATH, eventhdlrdata->optsetname);
      break;
   case SOLVINGSTAGE_INFEASIBLE:
      sprintf(paramfilename, DEFAULT_SETTINGFILEPATH, eventhdlrdata->infeasiblesetname);
      break;
   default:
      SCIPdebugMessage("Unknown solving stage: %d -> ABORT!\n ", eventhdlrdata->solvingstage);
      SCIPABORT();
      break;
   }
   assert(paramfilename != NULL);
   file = fopen(paramfilename, "r");
   if( file == NULL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,"Changed solving stage to %d \n", eventhdlrdata->solvingstage);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,"Parameter file %s not found--keeping settings as before \n", paramfilename);
   }
   else
   {
      fclose(file);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,"Changed solving stage to %d -- \n", eventhdlrdata->solvingstage);
      //SCIP_CALL( SCIPresetParams(scip) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Reading parameters from file %s\n", paramfilename);
      SCIP_CALL( SCIPreadParams(scip, paramfilename) );

      eventhdlrdata->enabled = TRUE;
   }

   return SCIP_OKAY;
}

/* frees strings from the eventhandler data */
static
void dataFreeStrings(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*   eventhdlrdata
   )
{
   /* read solufile information for problem */
   /* free memory from previous run first (can happen if several problems are consecutively read in one session) */
   if (eventhdlrdata->nosolsetname != NULL )
   {
      assert(eventhdlrdata->suboptsetname != NULL);
      assert(eventhdlrdata->optsetname != NULL);
      assert(eventhdlrdata->infeasiblesetname != NULL);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->nosolsetname);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->suboptsetname);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->optsetname);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->infeasiblesetname);
   }
}
/*
 * Callback methods of event handler
 */

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_EVENTCOPY(eventCopySolvingstage)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* call inclusion method of node selector */
   SCIP_CALL( SCIPincludeEventHdlrSolvingstage(scip) );

   return SCIP_OKAY;
}
/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeSolvingstage)
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
SCIP_DECL_EVENTINIT(eventInitSolvingstage)
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

   eventhdlrdata->solvingstage = SOLVINGSTAGE_UNINITIALIZED;

   SCIPduplicateMemoryArray(scip, &eventhdlrdata->nosolsetname, eventhdlrdata->nosolsetnameparam, strlen(eventhdlrdata->nosolsetnameparam) + 1);
   SCIPduplicateMemoryArray(scip, &eventhdlrdata->suboptsetname, eventhdlrdata->suboptsetnameparam, strlen(eventhdlrdata->suboptsetnameparam) + 1);
   SCIPduplicateMemoryArray(scip, &eventhdlrdata->optsetname, eventhdlrdata->optsetnameparam, strlen(eventhdlrdata->optsetnameparam) + 1);
   SCIPduplicateMemoryArray(scip, &eventhdlrdata->infeasiblesetname, eventhdlrdata->infeasiblesetnameparam, strlen(eventhdlrdata->infeasiblesetnameparam) + 1);

   probname = SCIPstringGetBasename(SCIPgetProbName(scip));
   eventhdlrdata->optimalvalue = SCIPinfinity(scip);
   SCIP_CALL( searchSolufileForProbname(scip, probname, eventhdlrdata) );
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Optimal value for problem %s from solufile %s: %16.9g\n", probname, eventhdlrdata->solufilename, eventhdlrdata->optimalvalue);

   if( eventhdlrdata->enabled )
   {
      SCIP_CALL( applySolvingStage(scip, eventhdlrdata) );
      SCIP_CALL( SCIPcatchEvent(scip, EVENTHDLR_EVENT, eventhdlr, NULL, NULL) );
   }
   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolSolvingstage)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecSolvingstage)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   assert(SCIPeventGetType(event) & (EVENTHDLR_EVENT));

   if( eventhdlrdata->enabled )
   {
      SCIP_CALL( applySolvingStage(scip, eventhdlrdata) );
   }

   return SCIP_OKAY;
}

/** creates event handler for Solvingstage event */
SCIP_RETCODE SCIPincludeEventHdlrSolvingstage(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create Solvingstage event handler data */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   assert(eventhdlrdata != NULL);

   eventhdlrdata->solufilename = NULL;
   eventhdlrdata->nosolsetnameparam = NULL;
   eventhdlrdata->suboptsetnameparam = NULL;
   eventhdlrdata->optsetnameparam = NULL;
   eventhdlrdata->infeasiblesetnameparam = NULL;

   eventhdlrdata->nosolsetname = NULL;
   eventhdlrdata->suboptsetname= NULL;
   eventhdlrdata->optsetname= NULL;
   eventhdlrdata->infeasiblesetname= NULL;

   eventhdlr = NULL;


   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecSolvingstage, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopySolvingstage) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeSolvingstage) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitSolvingstage) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolSolvingstage) );

   /* add Solvingstage event handler parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "eventhdlr/"EVENTHDLR_NAME"/lambda12", "gap for phase transition 1->2",
            &eventhdlrdata->lambda12, FALSE, DEFAULT_LAMBA12, 0,MAXGAP, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "eventhdlr/"EVENTHDLR_NAME"/lambda23", "gap for phase transition 2->3",
               &eventhdlrdata->lambda23, FALSE, DEFAULT_LAMBA23, 0,MAXGAP, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/"EVENTHDLR_NAME"/enabled", "should the event handler be executed?",
         &eventhdlrdata->enabled, FALSE, DEFAULT_ENABLED, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "eventhdlr/"EVENTHDLR_NAME"/solufilename","file to parse solution information from",
         &eventhdlrdata->solufilename, FALSE, DEFAULT_SOLUFILENAME, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "eventhdlr/"EVENTHDLR_NAME"/nosolsetname", "bla",
            &eventhdlrdata->nosolsetnameparam, FALSE, DEFAULT_SETNAME, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "eventhdlr/"EVENTHDLR_NAME"/suboptsetname", "settings file for suboptimal solving stage",
               &eventhdlrdata->suboptsetnameparam, FALSE, DEFAULT_SETNAME, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "eventhdlr/"EVENTHDLR_NAME"/optsetname", "settings file for optimal solving stage",
               &eventhdlrdata->optsetnameparam, FALSE, DEFAULT_SETNAME, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "eventhdlr/"EVENTHDLR_NAME"/infeasiblesetname", "settings file for infeasible solving stage",
               &eventhdlrdata->infeasiblesetnameparam, FALSE, DEFAULT_SETNAME, NULL, NULL) );


   SCIP_CALL( SCIPaddLongintParam(scip, "eventhdlr/"EVENTHDLR_NAME"/nodeoffset", "node offset", &eventhdlrdata->nodeoffset,
         FALSE, DEFAULT_NODEOFFSET, 1, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/"EVENTHDLR_NAME"/fallback", "should the event handler fall back from optimal stage?",
            &eventhdlrdata->fallback, FALSE, DEFAULT_FALLBACK, NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip ,"eventhdlr/"EVENTHDLR_NAME"/transitionmethod", "transition method 'c','e','l','o','r'",
         &eventhdlrdata->transitionmethod, FALSE, DEFAULT_TRANSITIONMETHOD, TRANSITIONMETHODS, NULL, NULL) );

   return SCIP_OKAY;
}
