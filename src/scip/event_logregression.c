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

/**@file   event_logregression.c
 * @author Gregor Hendel
 *
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_logregression.h"
#include "scip/event_treeinfos.h"
#include "scip/event_estimation.h"
#include "string.h"
#include "math.h"

#define EVENTHDLR_NAME         "logregression"
#define EVENTHDLR_DESC         "event handler for log regression"

/*
 * Data structures
 */
#define DEFAULT_ENABLED FALSE /**< should the event handler be executed? */
#define EVENTHDLR_EVENT SCIP_EVENTTYPE_BESTSOLFOUND /**< the actual event to be caught */
#define DEFAULT_XTYPE    'n' /**< default type to use for log regression - (t)ime, (n)odes, (l)p iterations */
#define XTYPES           "lnt"
#define SQUARED(x) (x * x)
/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_Bool            enabled;             /**< should the event handler be executed? */
   char                 xtype;               /**< type to use for log regression - (t)ime, (n)odes, (l)p iterations */
   SCIP_Real            lastx;
   SCIP_Real            lasty;
   int                  n;
   SCIP_Real            c0;
   SCIP_Real            c1;
   SCIP_Real            sumx;
   SCIP_Real            sumy;
   SCIP_Real            sumxy;
   SCIP_Real            sumx2;
   SCIP_Real            sumy2;
   SCIP_Real            corrcoef;
};

/*
 * Local methods
 */

static
SCIP_RETCODE updateRegression(
   SCIP* scip,
   SCIP_EVENTHDLRDATA* eventhdlrdata
   )
{
   SCIP_Real x;
   SCIP_Real y;

   y = SCIPgetPrimalbound(scip);

   assert(!SCIPisInfinity(scip, y) && !SCIPisInfinity(scip, -y));

   switch( eventhdlrdata->xtype )
   {
      case 'l':
         if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVED )
            x = (SCIP_Real)SCIPgetNLPIterations(scip);
         else
            x = 1;
         break;
      case 'n':
         if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVED )
            x = (SCIP_Real)SCIPgetNTotalNodes(scip);
         else
            x = 1;
         break;
      case 't':
         x = SCIPgetSolvingTime(scip);
         break;
      default:
         return SCIP_INVALIDDATA;
         break;
   }

   /* prevent the calculation of logarithm too close to zero */
   x = MAX(x, .1);
   x = log(x);

   /* replace last observation if too close */
   if( eventhdlrdata->n > 0 && SCIPisEQ(scip, eventhdlrdata->lastx, x) )
   {
      eventhdlrdata->sumx2 -= SQUARED(eventhdlrdata->lastx);
      eventhdlrdata->sumy2 -= SQUARED(eventhdlrdata->lasty);
      eventhdlrdata->sumy -= eventhdlrdata->lasty;
      eventhdlrdata->sumx -= eventhdlrdata->lastx;
      eventhdlrdata->sumxy -= eventhdlrdata->lastx * eventhdlrdata->lasty;
   }
   else
      ++(eventhdlrdata->n);

   eventhdlrdata->lastx = x;
   eventhdlrdata->lasty = y;
   eventhdlrdata->sumx += x;
   eventhdlrdata->sumx2 += SQUARED(x);
   eventhdlrdata->sumxy += x * y;
   eventhdlrdata->sumy += y;
   eventhdlrdata->sumy2 += SQUARED(y);

   /* abort if there are not enough data points available */
   if( eventhdlrdata->n <= 2 )
      return SCIP_OKAY;

   /* compute slope                 */
   eventhdlrdata->c1 = (eventhdlrdata->n * eventhdlrdata->sumxy  -  eventhdlrdata->sumx * eventhdlrdata->sumy) /
       (eventhdlrdata->n * eventhdlrdata->sumx2 - SQUARED(eventhdlrdata->sumx));

   /* compute y-intercept           */
   eventhdlrdata->c0 = (eventhdlrdata->sumy * eventhdlrdata->sumx2  -  eventhdlrdata->sumx * eventhdlrdata->sumxy) /
       (eventhdlrdata->n * eventhdlrdata->sumx2  -  SQUARED(eventhdlrdata->sumx));

   /* compute correlation coeff     */
   eventhdlrdata->corrcoef = (eventhdlrdata->sumxy - eventhdlrdata->sumx * eventhdlrdata->sumy / eventhdlrdata->n) /
            sqrt((eventhdlrdata->sumx2 - SQUARED(eventhdlrdata->sumx)/eventhdlrdata->n) *
            (eventhdlrdata->sumy2 - SQUARED(eventhdlrdata->sumy)/eventhdlrdata->n));

   return SCIP_OKAY;
}

/*
 * Callback methods of event handler
 */

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_EVENTCOPY(eventCopyLogregression)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* call inclusion method of node selector */
   SCIP_CALL( SCIPincludeEventHdlrLogregression(scip) );

   return SCIP_OKAY;
}
/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeLogregression)
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
SCIP_DECL_EVENTINIT(eventInitLogregression)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   /* read solufile information for problem */

   if( eventhdlrdata->enabled )
   {
      SCIP_CALL( SCIPcatchEvent(scip, EVENTHDLR_EVENT, eventhdlr, NULL, NULL) );
   }
   eventhdlrdata->c0 = SCIP_INVALID;
   eventhdlrdata->c1 = SCIP_INVALID;
   eventhdlrdata->sumx = 0;
   eventhdlrdata->sumx2 = 0;
   eventhdlrdata->sumxy = 0;
   eventhdlrdata->sumy = 0;
   eventhdlrdata->sumy2 = 0;

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecLogregression)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   int nleavesbelowincumbent, nleavesbelowinccorrected;
   SCIP_Real minestimate, mincorrectedestimate;
   SCIP_Real rootcorrestimate;
   int nnodesleft;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   assert(SCIPeventGetType(event) & EVENTHDLR_EVENT);

   if( eventhdlrdata->enabled )
   {
      SCIP_CALL( updateRegression(scip, eventhdlrdata) );
   }

   nnodesleft = SCIPgetStage(scip) == SCIP_STAGE_SOLVING ? SCIPgetNNodesLeft(scip) : 0;

   mincorrectedestimate = rootcorrestimate = minestimate = SCIP_INVALID;
   nleavesbelowincumbent = nleavesbelowinccorrected = 0;

   SCIP_CALL( SCIPgetCorrectedEstimateData(scip, &mincorrectedestimate, &rootcorrestimate, &minestimate, &nleavesbelowinccorrected, &nleavesbelowincumbent, TRUE) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "    Currentestimates (root corr, min, mincorr, nbelow, ncorrbelow, n): %9.5f %9.5f %9.5f %d %d %d\n",
         SCIPretransformObj(scip, rootcorrestimate),
         SCIPretransformObj(scip, minestimate),
         SCIPretransformObj(scip, mincorrectedestimate),
         nleavesbelowincumbent, nleavesbelowinccorrected, nnodesleft);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " N rank1 nodes: %d (%d)", SCIPgetNRank1Nodes(scip), nnodesleft);
   SCIPverbMessage(scip ,SCIP_VERBLEVEL_NORMAL, NULL, " Time, nodes, LP iters: %.2f %"SCIP_LONGINT_FORMAT" %"SCIP_LONGINT_FORMAT"\n",
         SCIPgetSolvingTime(scip), SCIPgetNNodes(scip), SCIPgetNLPIterations(scip))
   ;
   if( eventhdlrdata->n >= 3 )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "   regression updated: y = %.2f log(x) + %.2f, corr=%.2f\n", eventhdlrdata->c1, eventhdlrdata->c0,eventhdlrdata->corrcoef);

   return SCIP_OKAY;
}

/** creates event handler for Logregression event */
SCIP_RETCODE SCIPincludeEventHdlrLogregression(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create Logregression event handler data */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   assert(eventhdlrdata != NULL);

   eventhdlr = NULL;
   eventhdlrdata->lastx = SCIP_INVALID;
   eventhdlrdata->lasty = SCIP_INVALID;
   eventhdlrdata->n = 0;



   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecLogregression, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopyLogregression) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeLogregression) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitLogregression) );

   /* add Logregression event handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/"EVENTHDLR_NAME"/enabled", "should the event handler be executed?",
         &eventhdlrdata->enabled, FALSE, DEFAULT_ENABLED, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "eventhdlr/"EVENTHDLR_NAME"/xtype", "x type for log regression - (t)ime, (n)odes, (l)p iterations",
                 &eventhdlrdata->xtype, FALSE, DEFAULT_XTYPE, XTYPES, NULL, NULL) );

   return SCIP_OKAY;
}
