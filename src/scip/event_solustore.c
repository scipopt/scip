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

/**@file   event_solustore.c
 * @brief  eventhdlr for solustore event
 * @author Jakob Witzig
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_solustore.h"

#define EVENTHDLR_NAME         "solustore"
#define EVENTHDLR_DESC         "event handler to increase the size of the solution storage"

#define DEFAULT_MAXSOL         1000000
#define DEFAULT_MAXORIGSOL     1000000

/*
 * Data structures
 */

/** event handler data */
struct SCIP_EventhdlrData
{
   int                   maxorigsol;              /** maximal number of solutions candidates to store in the solution storage of the original problem */
   int                   maxsol;                  /** maximal number of solutions to store in the solution storage */
   SCIP_Bool             reopt;                   /** is reoptimization enabled? */
   SCIP_Bool             init;                    /** is initialized */
};

void SCIPeventhdlrSolustoreDisable(
   SCIP*                 scip
)
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->reopt = FALSE;
   return;
}

/*
 * Callback methods of event handler
 */

/** copy method for event handler plugins (called when SCIP copies plugins) */
#define eventCopySolustore NULL;

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#define eventExitSolustore NULL;

/** frees specific event data */
#define eventDeleteSolustore NULL;

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINIT(eventInitSolustore)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /** HACK */
   /** check if all variable are binary, if not, disable reoptimization */
   if( !eventhdlrdata->init )
   {
      SCIP_CALL( SCIPgetBoolParam(scip, "reoptimization/enable", &eventhdlrdata->reopt) );
      if( eventhdlrdata->reopt && SCIPgetNVars(scip) - SCIPgetNBinVars(scip) > 0 )
         eventhdlrdata->reopt = FALSE;

      eventhdlrdata->init = TRUE;
   }

   if( eventhdlrdata->init && eventhdlrdata->reopt )
   {
      eventhdlrdata->maxorigsol = DEFAULT_MAXORIGSOL;
      eventhdlrdata->maxsol = DEFAULT_MAXSOL;

      /* add solustore event handler parameters */
      SCIP_CALL(SCIPunfixParam(scip, "limits/maxorigsol"));
      SCIP_CALL(SCIPunfixParam(scip, "limits/maxsol"));
      SCIPsetIntParam(scip, "limits/maxorigsol", eventhdlrdata->maxorigsol);
      SCIPsetIntParam(scip, "limits/maxsol", eventhdlrdata->maxsol);
   }

   return SCIP_OKAY;
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeSolustore)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolSolustore)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   int savesols;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIP_CALL( SCIPgetIntParam(scip, "reoptimization/savesols", &savesols) );

   if(eventhdlrdata->reopt && savesols != 0)
   {
      SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, NULL, NULL));
   }

return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolSolustore)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   int savesols;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIP_CALL( SCIPgetIntParam(scip, "reoptimization/savesols", &savesols) );

   if(eventhdlrdata->reopt && savesols != 0)
   {
      SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, NULL, -1));
   }

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecSolustore)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   int nsols;
   int param_val;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   if(!eventhdlrdata->reopt)
      return SCIP_OKAY;

   nsols = SCIPgetNSols(scip);

   if(nsols >= 0.95*eventhdlrdata->maxsol || nsols >= 0.95*eventhdlrdata->maxorigsol)
   {
      eventhdlrdata->maxsol = MIN(2147483640,2*eventhdlrdata->maxsol);
      eventhdlrdata->maxorigsol = MIN(2147483640,2*eventhdlrdata->maxorigsol);

      SCIPchgIntParam(scip, SCIPgetParam(scip, "limits/maxorigsol"), eventhdlrdata->maxsol);
      SCIPchgIntParam(scip, SCIPgetParam(scip, "limits/maxsol"), eventhdlrdata->maxsol);
   }

   if(eventhdlrdata->maxsol == INT_MAX || eventhdlrdata->maxorigsol == INT_MAX)
      printf("WARNING: maximum size of solutions storage reached ! ! ! !\n");

   SCIPgetIntParam(scip, "limits/maxsol", &param_val);
   if(nsols > param_val)
   {
      printf(" *** %d saved solutions and maxsol storage of size %d\n", nsols, param_val);
      assert(nsols <= param_val);
   }

   SCIPgetIntParam(scip, "limits/maxorigsol", &param_val);
   if(nsols > param_val)
   {
      printf(" *** %d saved orig solutions and maxorigsol storage of size %d\n", nsols, param_val);
      assert(nsols <= param_val);
   }

   return SCIP_OKAY;
}

/** creates event handler for solustore event */
SCIP_RETCODE SCIPincludeEventHdlrSolustore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   eventhdlrdata = NULL;

   /** create nodereopt event handler data */
   SCIP_CALL(SCIPallocMemory(scip, &eventhdlrdata));

   eventhdlr = NULL;

   SCIP_CALL(SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecSolustore, eventhdlrdata));

   eventhdlrdata->init = FALSE;

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL(SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeSolustore));
   SCIP_CALL(SCIPsetEventhdlrInit(scip, eventhdlr, eventInitSolustore));
   SCIP_CALL(SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolSolustore));
   SCIP_CALL(SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolSolustore));

   return SCIP_OKAY;
}
