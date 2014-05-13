/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_estimation.c
 * @brief  eventhdlr for estimation event
 * @author Gregor Hendel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "scip/event_estimation.h"
#include "scip/scip.h"
#include <string.h>

#define EVENTHDLR_NAME         "estimation"
#define EVENTHDLR_DESC         "event handler for estimation event"
#define EVENT_TYPE_ESTIMATION SCIP_EVENTTYPE_NODEBRANCHED | SCIP_EVENTTYPE_NODEFEASIBLE | SCIP_EVENTTYPE_NODEINFEASIBLE
#define DEFAULT_ENABLED        FALSE
#define HEADER "NUMBER  PRIMAL LOWER ESTIM ROOTLPPSESTIM DEPTH NCANDS NODELOWER LPSTAT"
#define DEFAULT_FILEOUTPUT FALSE
/*
 * Data structures
 */

/* TODO: fill in the necessary event handler data */

/** depth information structure */
struct DepthInfo
{
   int nnodes;
   SCIP_Real minestimate;
};
typedef struct DepthInfo DEPTHINFO;

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_VAR**           rootlpcands;
   int                  nrootlpcands;
   SCIP_Real*           rootlpcandsfrac;
   SCIP_Real            rootestim;
   SCIP_Real            correctionfactor;
   SCIP_Bool            enabled;            /**< is the estimation event handler enabled? */
   int                  eventfilterpos;     /**< the event filter position, or -1, if event has not (yet) been caught */
   char*                filename;           /**< file to keep node results */
   SCIP_Bool            fileoutput;         /**< should file output be generated? */
   DEPTHINFO**          depthinfos;
   int                  maxdepth;
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */
static
SCIP_Real calcCorrectionFactor(
   SCIP* scip,
   SCIP_Real rootlpsol,
   SCIP_Real rootestim,
   SCIP_Real incumbentsol
   )
{
   assert(SCIPisFeasGE(scip, incumbentsol, rootlpsol));
   assert(SCIPisFeasGE(scip, rootestim, rootlpsol));

   if( SCIPisInfinity(scip, incumbentsol) )
      return 1.0;
   if( rootestim <= incumbentsol )
      return 1.0;
   if( SCIPisFeasEQ(scip, rootestim, rootlpsol) )
      return 1.0;

   return (incumbentsol - rootlpsol) / (rootestim - rootlpsol);
}
/** free root LP solution from event handler data */
static
void dataFreeRootLPSol(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*   eventhdlrdata
   )
{
   assert(scip != NULL);
   assert(eventhdlrdata != NULL);

   assert(eventhdlrdata->nrootlpcands == 0 || eventhdlrdata->rootlpcands != NULL);

   if( eventhdlrdata->rootlpcands != NULL )
   {
      SCIPfreeMemoryArray(scip, &eventhdlrdata->rootlpcands);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->rootlpcandsfrac);
   }
   eventhdlrdata->nrootlpcands = 0;
}

/** stores root LP solution (after separation) in eventhdlrdata */
static
SCIP_RETCODE storeRootLPSol(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*     eventhdlrdata
   )
{
   SCIP_VAR** lpcandsroot;
   SCIP_Real* lpcandsfrac;

   assert(scip != NULL);
   assert(eventhdlrdata != NULL);
   assert(SCIPgetDepth(scip) == 0);

   /* free previously stored root LP solutions */
   if( eventhdlrdata->rootlpcands != NULL )
   {
      dataFreeRootLPSol(scip, eventhdlrdata);
   }
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OBJLIMIT )
   {
      SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcandsroot, NULL, &lpcandsfrac, &eventhdlrdata->nrootlpcands, NULL, NULL) );
   }

   if( eventhdlrdata->nrootlpcands > 0 )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &eventhdlrdata->rootlpcands, lpcandsroot, eventhdlrdata->nrootlpcands) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &eventhdlrdata->rootlpcandsfrac, lpcandsfrac, eventhdlrdata->nrootlpcands) );
   }
   return SCIP_OKAY;
}

static
SCIP_Real recalcClassicEstimate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_Real estimate;
   int i;
   SCIP_Real pscostcontribution;
   SCIP_Real incumbentsol;
   estimate = SCIPgetLowerboundRoot(scip);

   if( eventhdlrdata->nrootlpcands == 0 || eventhdlrdata->rootlpcands == NULL )
      return SCIP_INVALID;

   pscostcontribution = 0.0;
   for( i = 0; i < eventhdlrdata->nrootlpcands; ++i )
   {
      SCIP_Real upscore;
      SCIP_Real downscore;

      upscore = SCIPgetVarPseudocost(scip, eventhdlrdata->rootlpcands[i], SCIP_BRANCHDIR_UPWARDS);
      downscore = SCIPgetVarPseudocost(scip, eventhdlrdata->rootlpcands[i], SCIP_BRANCHDIR_DOWNWARDS);

      upscore *= (1.0 - eventhdlrdata->rootlpcandsfrac[i]);
      downscore *= eventhdlrdata->rootlpcandsfrac[i];
      pscostcontribution += MIN(upscore, downscore);
   }
   incumbentsol = SCIPgetUpperbound(scip);
   eventhdlrdata->correctionfactor = calcCorrectionFactor(scip, estimate, estimate + pscostcontribution, incumbentsol);

   eventhdlrdata->rootestim = estimate + eventhdlrdata->correctionfactor * pscostcontribution;

   return eventhdlrdata->rootestim;
}

SCIP_Real SCIPgetRootLPSolPscostEstimate(
   SCIP*                 scip
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   return recalcClassicEstimate(scip, eventhdlrdata);
}

static
SCIP_Real nodeGetCorrectedEstimate(
   SCIP_NODE*            node,
   SCIP_Real             correctionfactor
   )
{
   SCIP_Real nodelower;
   SCIP_Real estimdelta;
   nodelower = SCIPnodeGetLowerbound(node);
   estimdelta = SCIPnodeGetEstimate(node) - nodelower;
   assert(estimdelta >= 0);

   return nodelower + correctionfactor * estimdelta;
}

static
void nodesUpdateEstimateData(
   SCIP_NODE** nodes,
   int         nnodes,
   SCIP_Real   incumbentsol,
   SCIP_Real*  minestimate,
   SCIP_Real*  mincorrectedestimate,
   int*        nnodesbelowincumbent,
   int*        nnodesbelowincumbentcorrected,
   SCIP_Real   correctionfactor
   )
{
   int n;

   for( n = 0; n < nnodes; ++n )
   {
      SCIP_Real estim;
      SCIP_Real correctedestim;

      estim = SCIPnodeGetEstimate(nodes[n]);
      correctedestim = nodeGetCorrectedEstimate(nodes[n], correctionfactor);

      if( estim < incumbentsol )
         (*nnodesbelowincumbent)++;
      if( correctedestim < incumbentsol )
         (*nnodesbelowincumbentcorrected)++;

      if( estim < *minestimate )
         *minestimate = estim;
      if( correctedestim < *mincorrectedestimate )
         *mincorrectedestimate = correctedestim;
   }
}

SCIP_RETCODE SCIPgetCorrectedEstimateData(
   SCIP*                 scip,
   SCIP_Real*            mincorrectedestimate,
   SCIP_Real*            rootcorrectedestim,
   SCIP_Real*            minestimate,
   int*                  nnodesbelowincumbentcorrected,
   int*                  nnodesbelowincumbent,
   SCIP_Bool             recalcestim
)
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   SCIP_NODE** leaves;
   SCIP_NODE** children;
   SCIP_NODE** siblings;
   SCIP_Real incumbentsol;

   int nleaves;
   int nchildren;
   int nsiblings;

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   recalcestim = (recalcestim || eventhdlrdata->rootestim == SCIP_INVALID );

   if( recalcestim )
      *rootcorrectedestim = recalcClassicEstimate(scip , eventhdlrdata);
   else
      *rootcorrectedestim = eventhdlrdata->rootestim;

   incumbentsol = SCIPgetUpperbound(scip);

   *nnodesbelowincumbent = 0;
   *nnodesbelowincumbentcorrected = 0;
   *mincorrectedestimate = SCIPinfinity(scip);
   *minestimate = SCIPinfinity(scip);

   nleaves = nchildren = nsiblings = 0;
   SCIP_CALL( SCIPgetOpenNodesData(scip, &leaves, &children, &siblings, &nleaves, &nchildren, &nsiblings) );
   if( nchildren > 0 )
      nodesUpdateEstimateData(children, nchildren, incumbentsol, minestimate, mincorrectedestimate, nnodesbelowincumbent, nnodesbelowincumbentcorrected,
            eventhdlrdata->correctionfactor);
   if( nsiblings > 0 )
      nodesUpdateEstimateData(siblings, nsiblings, incumbentsol, minestimate, mincorrectedestimate, nnodesbelowincumbent, nnodesbelowincumbentcorrected,
            eventhdlrdata->correctionfactor);
   if( nleaves > 0 )
         nodesUpdateEstimateData(leaves, nleaves, incumbentsol, minestimate, mincorrectedestimate, nnodesbelowincumbent, nnodesbelowincumbentcorrected,
               eventhdlrdata->correctionfactor);


   return SCIP_OKAY;

}

static
SCIP_RETCODE createDepthinfo(
   SCIP* scip,
   DEPTHINFO** depthinfo
   )
{
   SCIP_CALL( SCIPallocMemory(scip, depthinfo) );

   (*depthinfo)->minestimate = SCIPinfinity(scip);
   (*depthinfo)->nnodes = 0;

   return SCIP_OKAY;
}

static
SCIP_RETCODE freeDepthinfo(
   SCIP* scip,
   DEPTHINFO** depthinfo
   )
{
   SCIPfreeMemory(scip, depthinfo);

   return SCIP_OKAY;
}

static
void updateDepthinfo(
   DEPTHINFO* depthinfo,
   SCIP_NODE* node
   )
{
   if( SCIPnodeGetEstimate(node) < depthinfo->minestimate )
      depthinfo->minestimate = SCIPnodeGetEstimate(node);

   ++(depthinfo->nnodes);
}
static
SCIP_RETCODE storeDepthInfo(
   SCIP* scip,
   SCIP_EVENTHDLRDATA* eventhdlrdata,
   SCIP_NODE* node
   )
{
   int nodedepth;
   int newsize;
   int oldsize;

   nodedepth = SCIPnodeGetDepth(node);
   oldsize = eventhdlrdata->maxdepth;
   newsize = oldsize;
   if( oldsize == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->depthinfos, 10) );
      newsize = 10;
   }
   else if(nodedepth >= eventhdlrdata->maxdepth)
   {
      assert(nodedepth > 0);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &eventhdlrdata->depthinfos, 2 * nodedepth) );
      newsize = 2 * nodedepth;
   }

   if( newsize > oldsize )
   {
      int c;
      for( c = oldsize; c < newsize; ++c )
      {
         SCIP_CALL( createDepthinfo(scip, &(eventhdlrdata->depthinfos[c])) );

      }

      eventhdlrdata->maxdepth = newsize;
   }

   assert(newsize > nodedepth);

   updateDepthinfo(eventhdlrdata->depthinfos[nodedepth], node);

   return SCIP_OKAY;
}


/*
 * Callback methods of event handler
 */

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolEstimation)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata->eventfilterpos == -1);
   eventhdlrdata->correctionfactor = 1.0;
   eventhdlrdata->rootestim = SCIP_INVALID;
   eventhdlrdata->depthinfos = NULL;
   eventhdlrdata->maxdepth = 0;

   if( eventhdlrdata->enabled )
   {
      SCIP_CALL( SCIPcatchEvent(scip, EVENT_TYPE_ESTIMATION, eventhdlr, NULL, &eventhdlrdata->eventfilterpos) );
      assert(eventhdlrdata->eventfilterpos >= 0);
   }
   return SCIP_OKAY;
}

static
SCIP_DECL_EVENTFREE(eventFreeEstimation)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata->rootlpcands == NULL && eventhdlrdata->rootlpcandsfrac == NULL);
   assert(eventhdlrdata != NULL);

   if( eventhdlrdata->filename != NULL )
      SCIPfreeMemoryArray(scip, &eventhdlrdata->filename);

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

static
SCIP_DECL_EVENTINIT(eventInitEstimation)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   FILE* file;
   char filename[SCIP_MAXSTRLEN];
   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   if( eventhdlrdata->filename != NULL )
      SCIPfreeMemory(scip, &eventhdlrdata->filename);

   sprintf(filename, "/OPTI/bzfhende/nodefiles/%s.npf",SCIPstringGetBasename(SCIPgetProbName(scip)));

   SCIPduplicateMemoryArray(scip, &eventhdlrdata->filename, filename, strlen(filename) + 1);

   /* open the file and overwrite its content */
   if( eventhdlrdata->fileoutput && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_NORMAL )
   {
      file = fopen(eventhdlrdata->filename, "w+");
      assert(file != NULL);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, file, HEADER"\n");
      fclose(file);
   }

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitEstimation)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   dataFreeRootLPSol(scip, eventhdlrdata);

   assert(eventhdlrdata->eventfilterpos == -1 || eventhdlrdata->enabled );
   if( eventhdlrdata->eventfilterpos >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, EVENT_TYPE_ESTIMATION, eventhdlr, NULL, eventhdlrdata->eventfilterpos) );
      eventhdlrdata->eventfilterpos = -1;
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolEstimation)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   dataFreeRootLPSol(scip, eventhdlrdata);

   assert(eventhdlrdata->eventfilterpos == -1 || eventhdlrdata->enabled );
   if( eventhdlrdata->eventfilterpos >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, EVENT_TYPE_ESTIMATION, eventhdlr, NULL, eventhdlrdata->eventfilterpos) );
      eventhdlrdata->eventfilterpos = -1;
   }

   if( eventhdlrdata->maxdepth > 0 )
   {
      int c;
      for( c = 0; c < eventhdlrdata->maxdepth; ++c )
      {
         SCIP_CALL( freeDepthinfo(scip, &(eventhdlrdata->depthinfos[c])) );
      }
      SCIPfreeMemoryArray(scip, &eventhdlrdata->depthinfos);
      eventhdlrdata->maxdepth = 0;
   }

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecEstimation)
{  /*lint --e{715}*/
   SCIP_NODE* focusnode;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata->eventfilterpos >= 0);

   focusnode = SCIPgetCurrentNode(scip);

   if( focusnode == NULL )
      return SCIP_OKAY;

   /* store the last root LP solution, if it is not already integer feasible */
   if( SCIPnodeGetDepth(focusnode) == 0 && SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEBRANCHED )
   {
      SCIP_CALL( storeRootLPSol(scip, eventhdlrdata) );
   }

   SCIP_CALL( storeDepthInfo(scip, eventhdlrdata, focusnode) );


   if( eventhdlrdata->fileoutput && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_NORMAL )
   {
      FILE* file;
      SCIP_Real lowerbound;
      SCIP_Real estimate;
      SCIP_Real upperbound;
      SCIP_NODE* lowerboundnode;
      int depth;
      int nlpcands;

      nlpcands = -1;

      if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
         nlpcands = SCIPgetNLPBranchCands(scip);


      lowerboundnode = SCIPgetBestboundNode(scip);

      lowerbound = lowerboundnode != NULL ? SCIPnodeGetLowerbound(lowerboundnode) : SCIPinfinity(scip);
      estimate = SCIPnodeGetEstimate(focusnode);
      upperbound = SCIPgetUpperbound(scip);
      depth = SCIPnodeGetDepth(focusnode);
      file = fopen(eventhdlrdata->filename, "a");
      assert(file != NULL);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, file, "%lld %15g %15g %15g %15g %d %d %15g %d\n",
            SCIPnodeGetNumber(focusnode),
            SCIPgetExternalValue(scip, upperbound),
            SCIPgetExternalValue(scip, lowerbound),
            SCIPgetExternalValue(scip, estimate),
            SCIPgetExternalValue(scip, SCIPgetRootLPSolPscostEstimate(scip)),
            depth, nlpcands, SCIPgetExternalValue(scip, SCIPnodeGetLowerbound(focusnode)), SCIPgetLPSolstat(scip));
      fclose(file);
   }

   return SCIP_OKAY;
}

/** creates event handler for estimation event */
SCIP_RETCODE SCIPincludeEventHdlrEstimation(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create estimation event handler data */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   assert(eventhdlrdata != NULL);
   eventhdlrdata->rootlpcands = NULL;
   eventhdlrdata->rootlpcandsfrac = NULL;
   eventhdlrdata->eventfilterpos = -1;
   eventhdlrdata->filename = NULL;
   eventhdlrdata->depthinfos = NULL;
   eventhdlrdata->maxdepth = 0;

   eventhdlr = NULL;
   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecEstimation, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitEstimation) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolEstimation) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolEstimation) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeEstimation) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitEstimation) );

   /* add estimation event handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/estimation/enabled","enable event handler to perform estimation",
               &eventhdlrdata->enabled, FALSE, DEFAULT_ENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/estimation/fileoutput","should file output be generated? ",
         &eventhdlrdata->fileoutput, FALSE, DEFAULT_FILEOUTPUT, NULL, NULL) );

   return SCIP_OKAY;
}
