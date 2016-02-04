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

/**@file   nodesel_linprojection.c
 * @brief  linprojection node selector
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/nodesel_linprojection.h"
#include "scip/event_solvingphase.h"

#define NODESEL_NAME            "linprojection"
#define NODESEL_DESC            "uses a linear projection method to estimate the best MIP solution in the nodes subtree"
#define NODESEL_STDPRIORITY     -30000
#define NODESEL_MEMSAVEPRIORITY -30000

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "linprojection"
#define CONSHDLR_DESC          "constraint handler for linear projection node selector"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          -1 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY   INT_MIN /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP/**< propagation timing mask of the constraint handler*/

#define CONSHDLR_MAXPREROUNDS         0 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */

#define EVENTHDLR_NAME         "linprojection"
#define EVENTHDLR_DESC         "event handler for linprojection node events"
#define EVENT_TYPE_LINPROJECTION SCIP_EVENTTYPE_NODEBRANCHED | SCIP_EVENTTYPE_NODEFEASIBLE | SCIP_EVENTTYPE_NODEINFEASIBLE

#define DEFAULT_ADDNODEDATA FALSE

/* structure to collect descriptive statistics */
struct SCIP_DescStat
{
   char*                name;
   SCIP_Real            min;
   SCIP_Real            max;
   SCIP_Real            sumofobservations;
   SCIP_Real            sumofsquaredobservations;
   SCIP_Longint         nobservations;
};
typedef struct SCIP_DescStat SCIP_DESCSTAT;

/*
 * Data structures
 */
struct SCIP_ConsData
{
   int                   nodenlpcands;       /**< number of LP branch candidates at a node */
};


struct SCIP_EventhdlrData
{
   SCIP_NODESELDATA*    nodeseldata;         /**< data of linear projection node selector */
};
/* TODO: fill in the necessary node selector data */

/** node selector data */
struct SCIP_NodeselData
{
   SCIP_CONSHDLR*        conshdlr;           /**< the constraint handler for the node data */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for linprojection events */
   SCIP_CONS*            rootnodecons;       /**< pointer to hold the constraint for the (effective) root node */
   int                   rootconsdepth;      /**< depth of the root constraint (can be > 0 in case of bounded nodes!)*/
   int                   minnlpcands;        /**< minimum number of LP candidates found so far */
   SCIP_Real             minlowerbound;      /**< best lower bound of a node with minnlpcands branching candidates */
   SCIP_Real             projectionslope;    /**< slope of linear estimate of best MIP-feasible solution attainable in a subtree */
   SCIP_Real             cutoffbound;        /**< estimation of best MIP objective attainable */
   SCIP_DESCSTAT*        newestimationstat;  /**< descriptive statistics of new linear projection estimate */
   SCIP_DESCSTAT*        oldestimationstat;  /**< descriptive statistics of classic estimation method */
   SCIP_Real             optsolvalue;        /**< optimal solution value according to solu file */
   int                   nrootlpcands;       /**< number of root LP candidates */
   SCIP_Bool             addnodedata;        /**< should data be added as node constraints? */
   int                   eventfilterpos;     /**< event filter position for node events */
};

static
void SCIPdescstatReset(
   SCIP*                scip,
   SCIP_DESCSTAT*       descstat
   )
{
   assert(descstat != NULL);

   descstat->min = SCIPinfinity(scip);
   descstat->max = -descstat->min;
   descstat->nobservations = 0l;
   descstat->sumofobservations = 0.0;
   descstat->sumofsquaredobservations = 0.0;
}
static
SCIP_RETCODE SCIPdescstatCreate(
   SCIP*                scip,
   SCIP_DESCSTAT**      descstat,
   const char*          name
   )
{
   assert(scip != NULL);
   assert(descstat != NULL);
   assert(name != NULL);

   SCIP_CALL( SCIPallocMemory(scip, descstat) );

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*descstat)->name, name, strlen(name) + 1) );

   SCIPdescstatReset(scip, *descstat);

   return SCIP_OKAY;
}

static
void SCIPdescstatFree(
   SCIP*               scip,
   SCIP_DESCSTAT**     descstat
   )
{
   assert(scip != NULL);
   assert(descstat != NULL);
   assert(*descstat != NULL);
   assert((*descstat)->name != NULL);

   SCIPfreeMemoryArray(scip, &(*descstat)->name);

   SCIPfreeMemory(scip, descstat);
}

static
void SCIPdescstatAddObservation(
   SCIP_DESCSTAT*        descstat,
   SCIP_Real             newobservation
   )
{
   assert(descstat != NULL);

   ++descstat->nobservations;
   descstat->sumofobservations += newobservation;
   descstat->sumofsquaredobservations += newobservation * newobservation;
   if( newobservation < descstat->min )
      descstat->min = newobservation;
   if( newobservation > descstat->max )
      descstat->max = newobservation;
}


static
void SCIPprintDescstat(
   SCIP*                 scip,
   SCIP_DESCSTAT*        descstat,
   FILE*                 file)
{
   SCIP_Real meanval;
   SCIP_Real variance;

   assert(scip != NULL);
   assert(descstat != NULL);
   meanval = descstat->nobservations == 0 ? 0.0 : descstat->sumofobservations / (SCIP_Real)descstat->nobservations;
   variance = descstat->nobservations == 0 ? 0.0 : descstat->sumofsquaredobservations / (SCIP_Real)descstat->nobservations;
   variance -= meanval * meanval;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, file, "%s: range [%15.4g, %15.4g], mean %15.4g, variance %15.4g (%"SCIP_LONGINT_FORMAT" observations)\n",
         descstat->name, descstat->min, descstat->max, meanval, variance, descstat->nobservations);
}

/** creates and captures a linprojection constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
static
SCIP_RETCODE createConsBasicLinprojection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESELDATA*     nodeseldata,
   SCIP_NODE*            node,
   int                   nodenlpcands
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   char name[64];
   /* find the xyz constraint handler */
   assert(node != NULL);
   assert(nodeseldata != NULL);
   assert(nodenlpcands >= 0);

   conshdlr = nodeseldata->conshdlr;
   assert(conshdlr != NULL);

   sprintf(name, "node_%"SCIP_LONGINT_FORMAT"", SCIPnodeGetNumber(node));
   /* create constraint data */
   SCIP_CALL( SCIPallocMemory(scip, &consdata) );

   /* TODO: create and store constraint specific data here */
   consdata->nodenlpcands = nodenlpcands;
   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, &cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
         TRUE, FALSE, FALSE, TRUE, TRUE) );

   SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );

   /* if the node is effectively the root node, the constraint might be updated to a global one. We store this as the
    * unique root node cons data and release it later
    */
   if( SCIPconsIsLocal(cons) )
   {
      /* local constraints can be immediately released because they are kept in the node structure */
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   else
   {
      /* update the root node constraint: release a previously allocated root cons and keep the depth of the node to
       * to later identify this node when the constraint data is retrieved
       */
      if( nodeseldata->rootnodecons != NULL )
      {
         /* a previously allocated root node constraint can be released at this point */
         SCIP_CALL( SCIPreleaseCons(scip, &nodeseldata->rootnodecons) );
      }
      /* store the new root node constraint */
      nodeseldata->rootnodecons = cons;
      nodeseldata->rootconsdepth = SCIPnodeGetDepth(node);
   }
   return SCIP_OKAY;
}


static
SCIP_CONSDATA* nodeGetLinprojectionData(
   SCIP_NODE*           node,
   SCIP_NODESELDATA*    nodeseldata
   )
{
   SCIP_CONSSETCHG* conssetchg;
   SCIP_CONS** nodeconss;
   int         nnodeconss;
   int c;

   assert(node != NULL);
   assert( SCIPnodeGetDepth(node) >= nodeseldata->rootconsdepth );
   conssetchg = SCIPnodeGetConssetchg(node);

   /* this node is effectively the root node */
   if( SCIPnodeGetDepth(node) == nodeseldata->rootconsdepth )
   {
      assert(nodeseldata->rootnodecons != NULL);
      return SCIPconsGetData(nodeseldata->rootnodecons);
   }
   else if( conssetchg == NULL )
      return NULL;

   SCIPconssetchgGetAddedConsData(conssetchg, &nodeconss, &nnodeconss);
   assert(nnodeconss == 0 || nodeconss != NULL);

   /* loop over node constraints to find the unique one belonging to the linprojection conshdlr */
   for( c = nnodeconss - 1; c >= 0; --c )
   {
      SCIP_CONS* cons;

      cons = nodeconss[c];

      if(SCIPconsGetHdlr(cons) != nodeseldata->conshdlr)
         continue;

      return SCIPconsGetData(cons);
   }
   return NULL;

}

static
SCIP_RETCODE storeNodeInformation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESELDATA*     nodeseldata,
   SCIP_NODE*            node,
   int                   nodenlpcands
   )
{

   SCIP_CONSDATA* linprojectiondata;
   linprojectiondata = nodeGetLinprojectionData(node, nodeseldata);
   if( linprojectiondata == NULL )
      createConsBasicLinprojection(scip, nodeseldata, node, nodenlpcands);
   else
      linprojectiondata->nodenlpcands = nodenlpcands;

   return SCIP_OKAY;
}

static
int nodeGetNLPCands(
   SCIP_NODE*           node,
   SCIP_NODESELDATA*    nodeseldata
   )
{

   SCIP_CONSDATA* linprojectiondata;

   assert(node != NULL);
   assert(nodeseldata != NULL);

   linprojectiondata = nodeGetLinprojectionData(node, nodeseldata);
   if( linprojectiondata == NULL )
      return -1;
   else
      return linprojectiondata->nodenlpcands;
}

static
SCIP_Real nodeGetLinprojectionEstimate(
   SCIP_NODE*           node,
   SCIP_NODESELDATA*    nodeseldata
   )
{
   SCIP_Real nodelowerbound;
   int       nodenlpcands;

   nodelowerbound = SCIPnodeGetLowerbound(node);
   nodenlpcands = nodeGetNLPCands(node, nodeseldata);

   /* use parent number of LP cands which is usually an upper bound on the exspected LP candidate number for this node */
   if( nodenlpcands == -1 )
   {
      SCIP_NODE* parent;
      parent = SCIPnodeGetParent(node);
      if( parent != NULL )
         nodenlpcands = nodeGetNLPCands(parent, nodeseldata);
      assert( parent == NULL || nodenlpcands >= 0 );
   }
   return nodelowerbound + nodeseldata->projectionslope * nodenlpcands;
}

/** select a node with lowest projection estimate from the nodes array */
static
void selectBestNode(
   SCIP*                scip,
   SCIP_NODESELDATA*    nodeseldata,
   SCIP_NODE**          nodes,
   SCIP_NODE**          selnode,
   SCIP_Real*           bestnodeestimate,
   SCIP_Real*           bestnodelowerbound,
   int                  nnodes
   )
{
   int c;

   for( c = 0; c < nnodes; ++c )
   {
      SCIP_Real nodelinprojectionestimate;
      assert(nodes[c] != NULL);
      nodelinprojectionestimate = nodeGetLinprojectionEstimate(nodes[c], nodeseldata);
      if( *selnode == NULL || SCIPisLT(scip, nodelinprojectionestimate, *bestnodeestimate) )
      {
         *bestnodeestimate = nodelinprojectionestimate;
         *selnode = nodes[c];
         *bestnodelowerbound = SCIPnodeGetLowerbound(nodes[c]);
      }
   }
}

/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of constraint handler
 */
static
SCIP_DECL_CONSDELETE(consDeleteLinprojection)
{
   assert(scip != NULL);
   assert(*consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   SCIPfreeMemory(scip, consdata);

   return SCIP_OKAY;
}
static
SCIP_DECL_CONSCHECK(consCheckLinprojection)
{
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSLOCK(consLockLinprojection)
{
   return SCIP_OKAY;
}

static
SCIP_DECL_CONSENFOLP(consEnfolpLinprojection)
{
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}
static
SCIP_DECL_CONSENFOPS(consEnfopsLinprojection)
{
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
void updateLinprojectionslope(
   SCIP*                 scip,
   SCIP_NODESELDATA*     nodeseldata,
   SCIP_Real             newminobjective,
   int                   newmincands
   )
{
   SCIP_NODE* root;
   SCIP_Real rootlowerbound;
   SCIP_CONSDATA* rootconsdata;
   int nrootlpcands;

   assert(newmincands < nodeseldata->minnlpcands || SCIPisLT(scip, newminobjective, nodeseldata->minlowerbound));
   root = SCIPgetRootNode(scip);
   assert(root != NULL);

   rootlowerbound = SCIPgetLowerboundRoot(scip);
   assert(nodeseldata->rootnodecons != NULL);
   rootconsdata = SCIPconsGetData(nodeseldata->rootnodecons);
   nrootlpcands = rootconsdata->nodenlpcands;

   nodeseldata->minlowerbound = newminobjective;
   nodeseldata->minnlpcands = newmincands;
   assert(newmincands <= nrootlpcands);
   assert(SCIPisSumGE(scip, newminobjective, rootlowerbound));

   if( newmincands < nrootlpcands )
   {
      SCIP_Real objdiff;

      objdiff = newminobjective - rootlowerbound;
      objdiff = MAX(objdiff, 0.0);
      nodeseldata->projectionslope = objdiff / (SCIP_Real)(nrootlpcands - newmincands);
   }
   else
      nodeseldata->projectionslope = 0.0;

   nodeseldata->cutoffbound = SCIPgetLowerbound(scip) + nodeseldata->projectionslope * nrootlpcands;
   nodeseldata->cutoffbound = MIN(SCIPgetUpperbound(scip), nodeseldata->cutoffbound);
}


/*
 * Callback methods of event handler
 */
/* event execution method of linprojection event handler; method keeps track of node LP branch cands */
static
SCIP_DECL_EVENTEXEC(eventExecLinprojection)
{
   SCIP_NODE* focusnode;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_NODESELDATA* nodeseldata;

   SCIP_Real nodelowerbound;
   SCIP_Real currentupperbound;
   int nlpcands;
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(EVENTHDLR_NAME, SCIPeventhdlrGetName(eventhdlr)) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   nodeseldata = eventhdlrdata->nodeseldata;
   assert(nodeseldata != NULL);
   focusnode = SCIPgetCurrentNode(scip);
   assert(focusnode != NULL);

   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEBRANCHED )
      nlpcands = SCIPgetNLPBranchCands(scip);
   else
      nlpcands = 0;
   assert(nlpcands == 0 || SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEBRANCHED);
   /* init or update the node data */
   SCIP_CALL( storeNodeInformation(scip, nodeseldata, focusnode, nlpcands) );

   nodelowerbound = SCIPnodeGetLowerbound(focusnode);
   currentupperbound = SCIPgetUpperbound(scip);

   if( SCIPgetNLimSolsFound(scip) > 0 && SCIPisLT(scip, currentupperbound, nodeseldata->minlowerbound) )
   {
      assert(!SCIPisInfinity(scip, currentupperbound));
      updateLinprojectionslope(scip, nodeseldata, currentupperbound, 0);
   }
   else if( nlpcands < nodeseldata->minnlpcands ||
         (nlpcands == nodeseldata->minnlpcands && SCIPisLT(scip, nodelowerbound, nodeseldata->minlowerbound)) )
   {
      updateLinprojectionslope(scip, nodeseldata, nodelowerbound, nlpcands);
   }
   return SCIP_OKAY;

}


/*
 * Callback methods of node selector
 */

/* TODO: Implement all necessary node selector methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_NODESELCOPY(nodeselCopyLinprojection)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeNodeselLinprojection(scip) );

   return SCIP_OKAY;
}

static
SCIP_DECL_NODESELINIT(nodeselInitLinprojection)
{
   return SCIP_OKAY;
}

static
SCIP_DECL_NODESELINITSOL(nodeselInitsolLinprojection)
{
   SCIP_NODESELDATA* nodeseldata;

   assert(scip != NULL);
   assert(nodesel != NULL);
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   nodeseldata->minnlpcands = INT_MAX;
   nodeseldata->minlowerbound = SCIP_REAL_MAX;
   nodeseldata->projectionslope = 0.0;
   nodeseldata->rootconsdepth = -1;
   nodeseldata->rootnodecons = NULL;
   nodeseldata->cutoffbound = SCIPinfinity(scip);

   SCIPdescstatReset(scip, nodeseldata->newestimationstat);
   SCIPdescstatReset(scip, nodeseldata->oldestimationstat);

   nodeseldata->optsolvalue = SCIPgetOptimalSolutionValue(scip);

   if( SCIPisInfinity(scip, REALABS(nodeseldata->optsolvalue)) && SCIPgetNRuns(scip) <= 1 )
   {
      SCIPwarningMessage(scip, "Optimal solution value is infinite!\n");
   }

   if( nodeseldata->addnodedata && nodeseldata->eventfilterpos == -1 )
   {
      SCIP_CALL( SCIPcatchEvent(scip, EVENT_TYPE_LINPROJECTION, nodeseldata->eventhdlr, NULL, &nodeseldata->eventfilterpos) );
   }
   return SCIP_OKAY;
}

static
SCIP_DECL_NODESELEXITSOL(nodeselExitsolLinprojection)
{
   SCIP_NODESELDATA* nodeseldata;

   assert(scip != NULL);
   assert(nodesel != NULL);
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->eventfilterpos == -1 || nodeseldata->addnodedata );

   if( nodeseldata->rootnodecons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &nodeseldata->rootnodecons) );
   }

   if( nodeseldata->eventfilterpos >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, EVENT_TYPE_LINPROJECTION, nodeseldata->eventhdlr, NULL, nodeseldata->eventfilterpos) );
      nodeseldata->eventfilterpos = -1;
   }
   nodeseldata->rootnodecons = NULL;
   nodeseldata->rootconsdepth = -1;

   return SCIP_OKAY;
}

static
SCIP_DECL_NODESELFREE(nodeselFreeLinprojection)
{
   SCIP_NODESELDATA* nodeseldata;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(nodeseldata->eventhdlr);
   assert(eventhdlrdata != NULL);
   SCIPfreeMemory(scip, &eventhdlrdata);
   SCIPeventhdlrSetData(nodeseldata->eventhdlr, NULL);

   SCIPstatistic(
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Old and new estimation statistics:\n");
      SCIPprintDescstat(scip, nodeseldata->oldestimationstat, NULL);
      SCIPprintDescstat(scip, nodeseldata->newestimationstat, NULL);
      )

   SCIPdescstatFree(scip, &nodeseldata->oldestimationstat);
   SCIPdescstatFree(scip, &nodeseldata->newestimationstat);

   SCIPfreeMemory(scip, &nodeseldata);
   SCIPnodeselSetData(nodesel, NULL);

   return SCIP_OKAY;
}
/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectLinprojection)
{
   SCIP_NODE** leaves;
   SCIP_NODE** children;
   SCIP_NODE** siblings;
   SCIP_NODESELDATA* nodeseldata;
   SCIP_Real bestnodeestimate;
   SCIP_Real cutoffbound;
   SCIP_Real bestlowerbound;
   int nleaves;
   int nchildren;
   int nsiblings;

   assert(nodesel != NULL);
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   SCIP_CALL( SCIPgetOpenNodesData(scip, &leaves, &children, &siblings, &nleaves, &nchildren, &nsiblings) );

   if( nsiblings + nchildren + nleaves == 0 )
      return SCIP_OKAY;

   *selnode = NULL;
   bestnodeestimate = SCIPinfinity(scip);
   bestlowerbound = SCIPinfinity(scip);
   cutoffbound = nodeseldata->cutoffbound;
   selectBestNode(scip, nodeseldata, children, selnode, &bestnodeestimate, &bestlowerbound, nchildren);

   if( *selnode == NULL || SCIPisGT(scip, bestlowerbound, cutoffbound) )
   {
      selectBestNode(scip, nodeseldata, siblings, selnode, &bestnodeestimate, &bestlowerbound, nsiblings);
      if( *selnode == NULL || SCIPisGT(scip, bestlowerbound, cutoffbound) )
      {
         SCIP_Real bestclassicestimate;
         SCIP_NODE* bestleaf;
         SCIP_NODE* bestchild;
         SCIP_NODE* bestsibling;

         selectBestNode(scip, nodeseldata, leaves, selnode, &bestnodeestimate, &bestlowerbound, nleaves);

         bestleaf = SCIPgetBestLeaf(scip);
         bestclassicestimate = SCIPinfinity(scip);
         if( bestleaf != NULL )
            bestclassicestimate = SCIPnodeGetEstimate(bestleaf);
         bestsibling = SCIPgetBestSibling(scip);
         if( bestsibling != NULL && SCIPisLT(scip, SCIPnodeGetEstimate(bestsibling), bestclassicestimate) )
            bestclassicestimate = SCIPnodeGetEstimate(bestsibling);
         bestchild = SCIPgetBestChild(scip);
         if( bestchild != NULL && SCIPisLT(scip, SCIPnodeGetEstimate(bestchild), bestclassicestimate) )
            bestclassicestimate = SCIPnodeGetEstimate(bestchild);


         if( !(SCIPisInfinity(scip, REALABS(bestclassicestimate)) || SCIPisInfinity(scip, REALABS(bestnodeestimate)) ))
         {
            SCIP_Real extclassicestimate;
            SCIP_Real extnewestimate;

            extclassicestimate = SCIPgetExternalValue(scip, bestclassicestimate);
            extnewestimate = SCIPgetExternalValue(scip, bestnodeestimate);

            if( !SCIPisInfinity(scip, REALABS(nodeseldata->optsolvalue)) )
            {
               SCIPdescstatAddObservation(nodeseldata->oldestimationstat, extclassicestimate - nodeseldata->optsolvalue);
               SCIPdescstatAddObservation(nodeseldata->newestimationstat, extnewestimate - nodeseldata->optsolvalue);
            }
            else
            {
               SCIPdescstatAddObservation(nodeseldata->oldestimationstat, extclassicestimate);
               SCIPdescstatAddObservation(nodeseldata->newestimationstat, extnewestimate);
            }
         }

      }
   }
   assert(*selnode != NULL);
   SCIPdebugMessage("Selection of node %"SCIP_LONGINT_FORMAT" with estimate %g\n", SCIPnodeGetNumber(*selnode), bestnodeestimate);

   return SCIP_OKAY;
}

/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompLinprojection)
{
   SCIP_Real nodelb1;
   SCIP_Real nodelb2;

   nodelb1 = SCIPnodeGetEstimate(node1);
   nodelb2 = SCIPnodeGetEstimate(node2);

   if( SCIPisGT(scip, nodelb1, nodelb2) )
      return 1;
   else if( SCIPisLT(scip, nodelb1, nodelb2) )
      return -1;
   return 0;
}

/*
 * node selector specific interface methods
 */

/** creates the linprojection node selector and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselLinprojection(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESELDATA* nodeseldata;
   SCIP_NODESEL* nodesel;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   /* create linprojection node selector data */
   nodeseldata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &nodeseldata) );
   nodeseldata->conshdlr = NULL;
   nodeseldata->nrootlpcands = 0;
   nodeseldata->eventfilterpos = -1;

   SCIP_CALL( SCIPdescstatCreate(scip, &nodeseldata->oldestimationstat, "classic estimate") );
   SCIP_CALL( SCIPdescstatCreate(scip, &nodeseldata->newestimationstat, "new estimate    ") );

   /* include dummy constraint handler to represent linprojection node data */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &nodeseldata->conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
               CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS, consEnfolpLinprojection, consEnfopsLinprojection, consCheckLinprojection,
               consLockLinprojection, NULL) );
   assert(nodeseldata->conshdlr != NULL);

   SCIP_CALL( SCIPsetConshdlrDelete(scip, nodeseldata->conshdlr, consDeleteLinprojection) );

   /* include event handler for linprojection events */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   eventhdlrdata->nodeseldata = nodeseldata;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &nodeseldata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecLinprojection, eventhdlrdata) );
   assert(nodeseldata->eventhdlr != NULL);

   /* include node selector */
   nodesel = NULL;
   SCIP_CALL( SCIPincludeNodeselBasic(scip, &nodesel, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
          nodeselSelectLinprojection, nodeselCompLinprojection, nodeseldata) );
   assert(nodesel != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetNodeselInit(scip, nodesel, nodeselInitLinprojection) );
   SCIP_CALL( SCIPsetNodeselCopy(scip, nodesel, nodeselCopyLinprojection) );
   SCIP_CALL( SCIPsetNodeselInitsol(scip, nodesel, nodeselInitsolLinprojection) );
   SCIP_CALL( SCIPsetNodeselFree(scip, nodesel, nodeselFreeLinprojection) );
   SCIP_CALL( SCIPsetNodeselExitsol(scip, nodesel, nodeselExitsolLinprojection) );

   /* add linprojection node selector parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "nodeselection/"NODESEL_NAME"/addnodedata", "should data be added as node constraints?",
         &nodeseldata->addnodedata, TRUE, DEFAULT_ADDNODEDATA, NULL, NULL) );
   return SCIP_OKAY;
}
