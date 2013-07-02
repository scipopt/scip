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
#define SCIP_DEBUG
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

#define EVENTHDLR_NAME         "solvingstage"
#define EVENTHDLR_DESC         "event handler for events which influence the solving stage"

/*
 * Data structures
 */
#define MAXLINELEN 1024   /**< limit for line length */
#define DEFAULT_ENABLED TRUE /**< should the event handler be executed? */
#define DEFAULT_SOLUFILENAME "short.solu" /**< default filename to read solution value from */
#define DEFAULT_SOLUFILEPATH "//nfs//optimi//kombadon//bzfhende//projects//scip-git//check//results//"
#define EVENTHDLR_EVENT SCIP_EVENTTYPE_BESTSOLFOUND /**< the actual event to be caught */
/* TODO: fill in the necessary event handler data */

/** enumerator to represent the event handler solving stage
 */
enum SolvingStage
{
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
   SCIP_Real            optimalvalue;        /**< value of optimal solution of the problem */
   SOLVINGSTAGE         solvingstage;        /**< the current solving stage */
};

/*
 * Local methods
 */

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
   SCIPdebugMessage("Trying to open %s\n", solufilename);
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
      obj = strtok(NULL, " ");

      /* SCIPdebugMessage("%s %s %s", status, name, obj); */
      if( strcasecmp(name, probname) != 0 )
         continue;

      if( strcmp(status, "=opt=") != 0 )
      {
         SCIPdebugMessage("No optimal solution value known for this problem!\n");
         break;
      }
      else
      {
         char* endofptr;
         double val;

         val = strtod(obj, &endofptr);
         eventhdlrdata->optimalvalue = val;
         SCIPdebugMessage("Optimal value for problem: %16.9g\n", eventhdlrdata->optimalvalue);
         break;
      }
   }
   SCIPfclose(file);

   return SCIP_OKAY;
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
   probname = SCIPgetProbName(scip);
   /* read solufile information for problem */
   eventhdlrdata->optimalvalue = SCIPinfinity(scip);
   SCIP_CALL( searchSolufileForProbname(scip, probname, eventhdlrdata) );

   SCIPdebugMessage("Optimal value for problem %s from solufile %s: %16.9g\n", probname, eventhdlrdata->solufilename, eventhdlrdata->optimalvalue);
   eventhdlrdata->solvingstage = SOLVINGSTAGE_NOSOLUTION;

   if( eventhdlrdata->enabled )
      SCIP_CALL( SCIPcatchEvent(scip, EVENTHDLR_EVENT, eventhdlr, NULL, NULL) );
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

   assert(SCIPeventGetType(event) & EVENTHDLR_EVENT);

   if( eventhdlrdata->solvingstage == SOLVINGSTAGE_NOSOLUTION )
   {
      eventhdlrdata->solvingstage = SOLVINGSTAGE_SUBOPTIMAL;
      SCIPdebugMessage("Changed solving stage of event handler to SUBOPTIMAL!\n");
   }

   if( !SCIPisInfinity(scip, eventhdlrdata->optimalvalue) && SCIPisFeasEQ(scip, eventhdlrdata->optimalvalue, SCIPgetPrimalbound(scip)) )
   {
      SCIP_NODESEL* nodeselector;

      eventhdlrdata->solvingstage = SOLVINGSTAGE_OPTIMAL;


      nodeselector = SCIPfindNodesel(scip, "dfs");
      assert(nodeselector != NULL);

      SCIPdebugMessage("Changed solving stage to OPTIMAL\n");
      SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, FALSE);

      SCIP_CALL( SCIPsetNodeselStdPriority(scip, nodeselector, 1000000) );
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_PARAMCHGD(paramChgdEnabled)
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);

   if( SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE || SCIPgetStage(scip) > SCIP_STAGE_EXITSOLVE )
      return SCIP_OKAY;

   eventhdlr = NULL;
   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   if( eventhdlrdata->enabled )
      SCIP_CALL( SCIPcatchEvent(scip, EVENTHDLR_EVENT, eventhdlr, NULL, NULL) );
   else
      SCIP_CALL( SCIPdropEvent(scip, EVENTHDLR_EVENT, eventhdlr, NULL, -1) );

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
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/"EVENTHDLR_NAME"/enabled", "should the event handler be executed?",
         &eventhdlrdata->enabled, FALSE, DEFAULT_ENABLED, paramChgdEnabled, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "eventhdlr/"EVENTHDLR_NAME"/solufilename","file to parse solution information from",
         &eventhdlrdata->solufilename, FALSE, DEFAULT_SOLUFILENAME, NULL, NULL) );

   return SCIP_OKAY;
}
