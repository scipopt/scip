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

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_VAR**           rootlpcands;
   int                  nrootlpcands;
   SCIP_Real*           rootlpcandsfrac;
   SCIP_Real            rootcorrestim;
   SCIP_Real            pscostcontribution;
   SCIP_Real            correctionfactor;
   SCIP_Bool            enabled;            /**< is the estimation event handler enabled? */
   int                  eventfilterpos;     /**< the event filter position, or -1, if event has not (yet) been caught */
   char*                filename;           /**< file to keep node results */
   SCIP_Bool            fileoutput;         /**< should file output be generated? */
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
   if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
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
void recalcClassicEstimate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_Real estimate;
   int i;
   SCIP_Real incumbentsol;
   estimate = SCIPgetLowerboundRoot(scip);

   if( eventhdlrdata->nrootlpcands == 0 || eventhdlrdata->rootlpcands == NULL )
      return;

   eventhdlrdata->pscostcontribution = 0.0;
   for( i = 0; i < eventhdlrdata->nrootlpcands; ++i )
   {
      SCIP_Real upscore;
      SCIP_Real downscore;

      upscore = SCIPgetVarPseudocost(scip, eventhdlrdata->rootlpcands[i], SCIP_BRANCHDIR_UPWARDS);
      downscore = SCIPgetVarPseudocost(scip, eventhdlrdata->rootlpcands[i], SCIP_BRANCHDIR_DOWNWARDS);

      upscore *= (1.0 - eventhdlrdata->rootlpcandsfrac[i]);
      downscore *= eventhdlrdata->rootlpcandsfrac[i];
      eventhdlrdata->pscostcontribution += MIN(upscore, downscore);
   }
   incumbentsol = SCIPgetUpperbound(scip);
   eventhdlrdata->correctionfactor = calcCorrectionFactor(scip, estimate, estimate + eventhdlrdata->pscostcontribution, incumbentsol);

   eventhdlrdata->rootcorrestim = estimate + eventhdlrdata->correctionfactor * eventhdlrdata->pscostcontribution;
}

SCIP_RETCODE SCIPupdatePsCostContributionRootSolEstimate(
   SCIP* scip,
   SCIP_VAR* var,
   SCIP_Real solvaldelta,
   SCIP_Real oldpscostval
   )
{
   SCIP_Real upscore, downscore, oldupscore, olddownscore;
   SCIP_Real downfrac, upfrac;
   SCIP_VARSTATUS status;
   SCIP_Real rootlpsol;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   status = SCIPvarGetStatus(var);
   if( !(status == SCIP_VARSTATUS_COLUMN || status == SCIP_VARSTATUS_FIXED) )
      return SCIP_OKAY;
   rootlpsol = SCIPvarGetRootSol(var);

   if( SCIPisFeasIntegral(scip, rootlpsol) )
      return SCIP_OKAY;

   if( SCIPisPositive(scip, solvaldelta) )
   {
      oldupscore = oldpscostval;
      olddownscore = SCIPgetVarPseudocost(scip, var, SCIP_BRANCHDIR_DOWNWARDS);
   }
   else if( SCIPisFeasNegative(scip, solvaldelta) )
   {
      oldupscore = SCIPgetVarPseudocost(scip, var, SCIP_BRANCHDIR_UPWARDS);
      olddownscore = oldpscostval;
   }
   else
      return SCIP_OKAY;

   eventhdlrdata = SCIPeventhdlrGetData(SCIPfindEventhdlr(scip, EVENTHDLR_NAME));

   if( eventhdlrdata->nrootlpcands == 0 || eventhdlrdata->rootlpcands == NULL )
      return SCIP_OKAY;

   upscore = SCIPgetVarPseudocost(scip, var, SCIP_BRANCHDIR_UPWARDS);
   downscore = SCIPgetVarPseudocost(scip, var, SCIP_BRANCHDIR_DOWNWARDS);
   assert(oldupscore >= 0 && upscore >= 0 && downscore >= 0 && olddownscore >= 0);
   downfrac = (rootlpsol - SCIPfloor(scip, rootlpsol));
   upfrac = SCIPceil(scip, rootlpsol) - rootlpsol;
   olddownscore *= downfrac;
   oldupscore *= upfrac;
   downscore *= downfrac;
   upscore *= upfrac;

   eventhdlrdata->pscostcontribution -= MIN(oldupscore, olddownscore);
   eventhdlrdata->pscostcontribution += MIN(upscore, downscore);

   eventhdlrdata->correctionfactor = calcCorrectionFactor(scip, SCIPgetLowerboundRoot(scip), SCIPgetLowerboundRoot(scip) + eventhdlrdata->pscostcontribution, SCIPgetUpperbound(scip));

   return SCIP_OKAY;
}

SCIP_Real SCIPgetRootLPSolPscostEstimate(
   SCIP*                 scip,
   SCIP_Bool             corrected
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   if( eventhdlrdata->rootlpcands == NULL )
      return SCIP_INVALID;

   if( corrected )
      return SCIPgetLowerboundRoot(scip) + eventhdlrdata->correctionfactor * eventhdlrdata->pscostcontribution;
   else
      return SCIPgetLowerboundRoot(scip) + eventhdlrdata->pscostcontribution;
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
   int*        rootestimumbent,
   int*        rootestimumbentcorrected,
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
         (*rootestimumbent)++;
      if( correctedestim < incumbentsol )
         (*rootestimumbentcorrected)++;

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
   int*                  rootestimumbentcorrected,
   int*                  rootestimumbent,
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

   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   recalcestim = (recalcestim || eventhdlrdata->rootcorrestim == SCIP_INVALID );

   if( recalcestim )
      recalcClassicEstimate(scip , eventhdlrdata);

   *rootcorrectedestim = eventhdlrdata->rootcorrestim;

   incumbentsol = SCIPgetUpperbound(scip);

   *rootestimumbent = 0;
   *rootestimumbentcorrected = 0;
   *mincorrectedestimate = SCIPinfinity(scip);
   *minestimate = SCIPinfinity(scip);

   nleaves = nchildren = nsiblings = 0;
   SCIP_CALL( SCIPgetOpenNodesData(scip, &leaves, &children, &siblings, &nleaves, &nchildren, &nsiblings) );
   if( nchildren > 0 )
      nodesUpdateEstimateData(children, nchildren, incumbentsol, minestimate, mincorrectedestimate, rootestimumbent, rootestimumbentcorrected,
            eventhdlrdata->correctionfactor);
   if( nsiblings > 0 )
      nodesUpdateEstimateData(siblings, nsiblings, incumbentsol, minestimate, mincorrectedestimate, rootestimumbent, rootestimumbentcorrected,
            eventhdlrdata->correctionfactor);
   if( nleaves > 0 )
         nodesUpdateEstimateData(leaves, nleaves, incumbentsol, minestimate, mincorrectedestimate, rootestimumbent, rootestimumbentcorrected,
               eventhdlrdata->correctionfactor);


   return SCIP_OKAY;

}

/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTCOPY(eventCopyEstimation)
{
   SCIP_CALL( SCIPincludeEventHdlrEstimation(scip) );

   return SCIP_OKAY;
}
/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolEstimation)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata->eventfilterpos == -1);
   eventhdlrdata->correctionfactor = 1.0;
   eventhdlrdata->rootcorrestim = SCIP_INVALID;
   eventhdlrdata->pscostcontribution = 0;

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
            SCIPgetExternalValue(scip, eventhdlrdata->rootcorrestim),
            depth, nlpcands, SCIPgetExternalValue(scip, SCIPnodeGetLowerbound(focusnode)), SCIPgetLPSolstat(scip));
      fclose(file);
   }

   return SCIP_OKAY;
}

#define DISP_NAME_ROOTESTIM         "rootestim"
#define DISP_DESC_ROOTESTIM         "current number of encountered infeasible leaves"
#define DISP_HEAD_ROOTESTIM         "rootestim"
#define DISP_WIDT_ROOTESTIM         13
#define DISP_PRIO_ROOTESTIM         40000
#define DISP_POSI_ROOTESTIM         1000
#define DISP_STRI_ROOTESTIM         TRUE

static
SCIP_DECL_DISPOUTPUT(dispOutputRootestim)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_Real rootestim;
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_ROOTESTIM) == 0);
   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(SCIPfindEventhdlr(scip, EVENTHDLR_NAME));
   recalcClassicEstimate(scip, eventhdlrdata);

   if( eventhdlrdata->nrootlpcands == 0 || eventhdlrdata->rootlpcands == NULL )
      rootestim = SCIP_INVALID;
   else
      rootestim = SCIPgetLowerboundRoot(scip) + eventhdlrdata->pscostcontribution;

   if( rootestim == SCIP_INVALID )
      SCIPinfoMessage(scip, file, "     --      ");
   else
      SCIPinfoMessage(scip, file, "%13.6e", SCIPretransformObj(scip, rootestim));
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
   eventhdlrdata->nrootlpcands = 0;
   eventhdlrdata->rootlpcands = NULL;
   eventhdlrdata->rootlpcandsfrac = NULL;
   eventhdlrdata->eventfilterpos = -1;
   eventhdlrdata->filename = NULL;

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
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopyEstimation) );

   /* add estimation event handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/estimation/enabled","enable event handler to perform estimation",
               &eventhdlrdata->enabled, FALSE, DEFAULT_ENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/estimation/fileoutput","should file output be generated? ",
         &eventhdlrdata->fileoutput, FALSE, DEFAULT_FILEOUTPUT, NULL, NULL) );

   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_ROOTESTIM, DISP_DESC_ROOTESTIM, DISP_HEAD_ROOTESTIM, SCIP_DISPSTATUS_ON,
         NULL, NULL, NULL, NULL, NULL, NULL, dispOutputRootestim, NULL, DISP_WIDT_ROOTESTIM, DISP_PRIO_ROOTESTIM, DISP_POSI_ROOTESTIM,
         DISP_STRI_ROOTESTIM) );

   return SCIP_OKAY;
}
