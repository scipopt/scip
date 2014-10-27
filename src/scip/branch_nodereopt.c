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

/**@file   branch_nodereopt.c
 * @brief  nodereopt branching rule
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <assert.h>
#include <string.h>
#include <math.h>

#include "scip/branch_pseudo.h"
#include "scip/branch_nodereopt.h"
#include "scip/branch_relpscost.h"
#include "scip/event_nodereopt.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/scip.h"
#include "scip/clock.h"
#include "scip/struct_scip.h"
#include "scip/struct_lp.h"
#include "scip/struct_stat.h"
#include "scip/nodesel_estimate.h"
#include "scip/tree.h"
#include "scip/reopt.h"

#define BRANCHRULE_NAME            "nodereopt"
#define BRANCHRULE_DESC            "branching rule for node reoptimization"
#define BRANCHRULE_PRIORITY        536870911
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

/*
 * Data structures
 */

/** Data for logic-or constraints to split nodes into two disjoint parts */

struct LogicOrData
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   REOPT_CONSTYPE constype;
   int nvars;
};

struct SolData
{
   SCIP_Real**           vars;           /** solution values of the variables order by the SCIPvarGetIndex() */
   int                   nsols;
   int                   size;
};

/** Data for ancestor branching paths */
struct NodeData
{
   SCIP_VAR**            vars;
   SCIP_Real*            varbounds;
   SCIP_BOUNDTYPE*       varboundtypes;

   SCIP_VAR**            afterdualvars;
   SCIP_Real*            afterdualvarbounds;
   SCIP_BOUNDTYPE*       afterdualvarboundtypes;

   SCIP_Bool             pseudobranched;
   SCIP_QUEUE*           nodechilds;
   SCIP_LPISTATE*        lpistate;
   SCIP_QUEUE*           conss;
   LOGICORDATA**         dualcons;
   LOGICORDATA*          consandprop;
   SOLDATA*              soldata;

   int                   nvars;
   int                   nafterdualvars;

   int                   parentID;
   SCIP_REOPTTYPE        reopttype;

   int                   allocmem;
   int                   afterdualallocmem;
};
typedef struct NodeData NODEDATA;

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Bool             init;                    /** is the data structure initialized? */

   /** Savings for reoptimization */
   NODEDATA**            nodedata;
   SCIP_QUEUE*           openIDs;                 /** queue of empty IDs in node-data array; ID 0 is always reserved for the root */
   int*                  currentleafIDs;          /** array of all IDs corresponding to nodes which are revived, key: node number */
   int*                  nodetoid;                /** array of all IDs corresponding to saved nodes; we need this to find the IDs of the
                                                      parent node. key: node number */
   int                   allocmemsizenodedata;    /** size of allocated memory for nodes. */
   int                   allocmemsizeNodeID;      /** size of array which links IDs and node numbers */
   int                   nsavednodes;             /** number of saved nodes */

   /** global logic-or constraints to separate binary solutions */
   LOGICORDATA**         globalcons;
   SCIP_Bool             sepaglbsols;             /** save global constraints to separate solution found so far */
   SCIP_Bool             sepaglbsubtree;          /** save global constraints to separate infeasible subtrees */
   int                   nconsglobal;             /** length of the two lists above */
   int                   nconsglobaladded;        /** count number of added constraints (only for statistic) */
   int                   allocmemglobalcons;      /** allocated memory in the two lists above */

   /** global infeasibility */
   SCIP_Bool             infeasibleproblem;       /** reoptimization proves infeasibilty */
   SCIP_Bool             infeasibleconsadded;     /** the infeasibility constraint is added */

   /** Statistic stuff */
   int                   nfeasnodes;              /** number of feasible nodes during the solving process */
   int                   nfeasnodesmax;
   int                   nfeasnodesround;         /** number of feasible nodes found in this run */

   int                   ninfeasnodes;            /** number of infeasible nodes during the solving process */
   int                   ninfeasnodesmax;
   int                   ninfeasnodesround;       /** number of infeasible nodes in this run. */

   int                   nprunednodes;            /** number of pruned nodes during the solving process */
   int                   nprunednodesmax;
   int                   nprunednodesround;       /** number of pruned nodes in this run. */

   int                   npbnodes;                /** number of pb nodes during the solving process */
   int                   npbnodesmax;
   int                   npbnodesround;           /** number of pb nodes in this run. */

   int                   nrediednodes;            /** number of redied nodes during the solving process */
   int                   nrediednodesmax;         /** maximal number of redied nodes */
   int                   nrediednodesround;       /** number of redied nodes in this run. */

   int                   nrestarts;               /** number of restarts */
   int                   nrestartsround;          /** number of restarts ins this run */
   int                   firstrestart;            /** first round a restart was processed */
   int                   lastrestart;             /** last round a restart was processed */
   int                   nruns;                   /** number of restarts */

   int                   ninfsubtrees;            /** number of infesible subtrees found */
   int                   nrevivednodes;

   /** clock */
   SCIP_CLOCK*           inittime;                /** time needed to initialize */
   SCIP_CLOCK*           updatesolutime;          /** time needed to update solutions */
   SCIP_CLOCK*           savetime;                /** time needed to save nodes */

   /** time for hypothetical multi-threading with 8 cores */
   SCIP_CLOCK**          threadclock;
   SCIP_CLOCK*           lpclock;
   int                   nlp;

   /** Parameter stuff */
   SCIP_Bool             reopt;                   /** boolean flag if the reoptimization should enabled. */
   SCIP_Bool             restart;                 /** should the problem be restarted? */
   SCIP_Bool             savelpbasis;             /** save the LP basis for feasible and strong-branched nodes */
   SCIP_Bool             dynamicdiffofnodes;      /** should the maximal number of bound changes in two ancestor nodes
                                                      calculated depending on number of variables? */
   SCIP_Bool             sepabestsol;             /** should the best solution be separated and forbidden, i.e, constraint shortest path? */
   SCIP_Bool             sepalocsols;             /** save local constraints */

   SCIP_Real             objsimrootLP;            /** similarity of two sequential objective function to disable solving the root LP. */

   int                   maxsavednodes;           /** maximal number of saved nodes for the reoptimization */
   int                   maxdiffofnodes;          /** maximal number of bound changes in two ancestor nodes such that
                                                      the path get not shrunk */
   int                   solvelp;                 /** at which reopttype should the LP be solved? (1: transit,
                                                      2: w/ added logicor, 3: strong branched, 4: only leafs). */
   int                   solvelpdiff;             /** number of bound changes between two ancestor nodes to solve
                                                      the LP, only if solvelp = 0 */

   /** heuristic */
   int                   cpressnodes;             /** type of nodes which should compressed 0: all, 1: pruned, 2: feasible */
   SCIP_Bool             reducetofrontier;        /** delete nodes which were not revived in the last iteration */
   SCIP_Bool             shrink;                  /** shrink transit nodes */
   SCIP_Bool             strongbranchinginit;     /** strong branching at the root nodes */

   SCIP_CLOCK*           lrtime;                  /** time spend in LR heuristic */
   SCIP_Real             lrloss;                  /** loss of the last LR */
   SCIP_Bool             lrenable;                /** enable LR heuristic */
   int                   lrcalls;                 /** number calls of LR heuristic */
   int                   lrsuccess;               /** number successfully calls of LR heuristic */
   int                   lrdepth;                 /** depth of the last LR */
   int                   lrk;                     /** size of the compression */

   SCIP_CLOCK*           wctime;                  /** time spend in WC heuristic */
   SCIP_Bool             wcenable;                /** enable WC heuristic */
   int                   wccalls;                 /** number calls of WC heuristic */
   int                   wcsuccess;               /** number successfully calls of WC heuristic */
   int                   wck;                     /** size of the compression */
   int                   wclastnodeID;            /** ID of the last compressed node */
   int                   wcdepth;                 /** number of fixed variables of the last successful run */
   int                   lastprunedID;
   SCIP_Longint          lastprunednr;
   SCIP_Real             lastpruneddualbound;

   SCIP_Real             localdelay;
   SCIP_Bool             dynamiclocaldelay;
   int                   ndelsubtrees;
};

/*
 * check allocated memory and if neccessary reallocate
 */
static
SCIP_RETCODE checkMemory(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   int                   mem
)
{
   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(nodeID >= 0);
   assert(nodeID < branchruledata->allocmemsizenodedata);
   assert(branchruledata->nodedata[nodeID] != NULL);
   assert(mem >= 0);

   if( branchruledata->nodedata[nodeID]->allocmem < mem )
   {
      branchruledata->nodedata[nodeID]->allocmem = mem + 1;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &branchruledata->nodedata[nodeID]->vars, branchruledata->nodedata[nodeID]->allocmem) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &branchruledata->nodedata[nodeID]->varbounds, branchruledata->nodedata[nodeID]->allocmem) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &branchruledata->nodedata[nodeID]->varboundtypes, branchruledata->nodedata[nodeID]->allocmem) );
   }

   return SCIP_OKAY;
}

/*
 *  static methods
 */
/** TODO: Braucht man das ueberhaupt, wenn Knoten mit PBCons und SepaCons immer gespeichert bleiben?! */
static
SCIP_RETCODE saveLocalConssData(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_NODE*            node,
   int                   nodeID
)
{
   SCIP_CONS** addedcons;
   SCIP_Real constant;
   SCIP_Real scalar;
   int var;
   int consnr;
   int naddedcons;

   assert(scip != NULL );
   assert(node != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(branchruledata->nodedata[nodeID]->conss != NULL );

   /** save the added pseudo-constraint */
   if(SCIPnodeGetNAddedcons(scip, node) > 0)
   {
      naddedcons = SCIPnodeGetNAddedcons(scip, node);

      /** get memory */
      SCIP_CALL( SCIPallocMemoryArray(scip, &addedcons, naddedcons) );
      SCIP_CALL( SCIPnodeGetAddedcons(scip, node, addedcons) );

      for(consnr = 0; consnr < naddedcons; consnr++)
      {
         LOGICORDATA* consdata;
         SCIPallocMemory(scip, &consdata);

         consdata->nvars = SCIPgetNVarsLogicor(scip, addedcons[consnr]);

         if( strcmp("sepasolution", SCIPconsGetName(addedcons[consnr])) == 0 )
            consdata->constype = REOPT_CONSTYPE_SEPASOLUTION;
         else if( strcmp("infsubtree", SCIPconsGetName(addedcons[consnr])) == 0 )
            consdata->constype = REOPT_CONSTYPE_INFSUBTREE;
         else if( strcmp("strbranch", SCIPconsGetName(addedcons[consnr])) == 0 )
            consdata->constype = REOPT_CONSTYPE_STRBRANCHED;

         assert(consdata->constype == REOPT_CONSTYPE_SEPASOLUTION
             || consdata->constype == REOPT_CONSTYPE_INFSUBTREE
             || consdata->constype == REOPT_CONSTYPE_STRBRANCHED);

         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vars, SCIPgetVarsLogicor(scip, addedcons[consnr]), consdata->nvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->vals, consdata->nvars) );

         for(var = 0; var < consdata->nvars; var++)
         {
            constant = 0;
            scalar = 1;

            if(!SCIPvarIsOriginal(consdata->vars[var]))
            {
               if(SCIPvarIsNegated(consdata->vars[var]))
               {
                  SCIP_CALL(SCIPvarGetOrigvarSum(&consdata->vars[var], &scalar, &constant));
                  consdata->vals[var] = 1;
               }
               else
               {
                  SCIP_CALL(SCIPvarGetOrigvarSum(&consdata->vars[var], &scalar, &constant));
                  consdata->vals[var] = 0;
               }
               assert(consdata->vars[var] != NULL );
            }
            assert(SCIPvarIsOriginal(consdata->vars[var]));
         }

         SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->conss, (void* ) (size_t ) consdata) );
      }

      assert(SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss) == naddedcons);
      SCIPfreeMemoryArray(scip, &addedcons);
   }

   return SCIP_OKAY;
}

/*
 * Reallocate memory for pruned nodes which have to be stored.
 */
static
SCIP_RETCODE reallocNodedata(
   SCIP* scip,
   SCIP_BRANCHRULEDATA* branchruledata)
{
   int oldsize;
   int pos;

   assert(scip != NULL );
   assert(branchruledata != NULL );

   oldsize = branchruledata->allocmemsizenodedata;
   branchruledata->allocmemsizenodedata = 2 * branchruledata->allocmemsizenodedata;

   /** data of node information */
   SCIP_CALL( SCIPreallocMemoryArray(scip, &branchruledata->nodedata, branchruledata->allocmemsizenodedata) );

   /** write 0 in empty slots and fill the queue */
   for(pos = oldsize; pos < branchruledata->allocmemsizenodedata; pos++)
   {
      SCIP_CALL( SCIPqueueInsert(branchruledata->openIDs, (void* ) (size_t ) pos) );
      branchruledata->nodedata[pos] = NULL;
   }
   assert(SCIPqueueNElems(branchruledata->openIDs) == branchruledata->allocmemsizenodedata / 2);

   return SCIP_OKAY;
}

#ifdef CCHECK
/*
 * Check if the reoptimization is consistent.
 */
static
SCIP_RETCODE checkConsistency(
      SCIP* scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   SCIP_Bool consistent;

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL);

   consistent = SCIPeventhdlrNodereoptCheckConsistency(scip, scip->stat->ncreatednodesrun-1);

   if(!consistent)
   {
      printf("\n** Reoptnodes are not consistent. **\n\n");
      assert(consistent);
   }

   return SCIP_OKAY;
}
#endif

static
SCIP_RETCODE resizeGlobalsConsData(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata
)
{
   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   branchruledata->allocmemglobalcons = 2 * branchruledata->allocmemglobalcons;

   SCIP_CALL( SCIPreallocMemoryArray(scip, &(branchruledata->globalcons), branchruledata->allocmemglobalcons) );

   return SCIP_OKAY;
}

/*
 * allocate memory for nodedata and initialize all node parameters
 */
static
SCIP_RETCODE initNode(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID
)
{
   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   if(branchruledata->nodedata[nodeID] == NULL )
   {
      SCIP_CALL(SCIPallocMemory(scip, &branchruledata->nodedata[nodeID]));
      branchruledata->nodedata[nodeID]->conss = NULL;
      branchruledata->nodedata[nodeID]->lpistate = NULL;
      branchruledata->nodedata[nodeID]->nodechilds = NULL;
      branchruledata->nodedata[nodeID]->nvars = 0;
      branchruledata->nodedata[nodeID]->nafterdualvars = 0;
      branchruledata->nodedata[nodeID]->parentID = -1;
      branchruledata->nodedata[nodeID]->pseudobranched = FALSE;
      branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_NONE;
      branchruledata->nodedata[nodeID]->allocmem = 0;
      branchruledata->nodedata[nodeID]->afterdualallocmem = 0;
      branchruledata->nodedata[nodeID]->vars = NULL;
      branchruledata->nodedata[nodeID]->varbounds = NULL;
      branchruledata->nodedata[nodeID]->varboundtypes = NULL;
      branchruledata->nodedata[nodeID]->afterdualvars = NULL;
      branchruledata->nodedata[nodeID]->afterdualvarbounds = NULL;
      branchruledata->nodedata[nodeID]->afterdualvarboundtypes = NULL;
      branchruledata->nodedata[nodeID]->soldata = NULL;

      /* [0] : constraint of the last iteration (need to add in this iteration);
       * [1] : constraint of the current iteration (need to add in the next iteration) */
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeID]->dualcons, 2) );
      branchruledata->nodedata[nodeID]->dualcons[0] = NULL;
      branchruledata->nodedata[nodeID]->dualcons[1] = NULL;

      branchruledata->nodedata[nodeID]->consandprop = NULL;
   }
   else
   {
      assert(branchruledata->nodedata[nodeID]->nvars == 0);
      branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_NONE;
   }

   return SCIP_OKAY;
}

/**
 * transform variable and bounds back to the originals
 */
static
SCIP_RETCODE transformIntoOrig(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID
)
{
   int varnr;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(nodeID >= 1);
   assert(branchruledata->nodedata[nodeID] != NULL );

   /* transform branching variables and bound changes applied before the first dual reduction */
   for(varnr = 0; varnr < branchruledata->nodedata[nodeID]->nvars; varnr++)
   {
      SCIP_Real constant;
      SCIP_Real scalar;

      scalar = 1;
      constant = 0;

      if(!SCIPvarIsOriginal(branchruledata->nodedata[nodeID]->vars[varnr]))
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&branchruledata->nodedata[nodeID]->vars[varnr], &scalar, &constant)) ;
         branchruledata->nodedata[nodeID]->varbounds[varnr] = (branchruledata->nodedata[nodeID]->varbounds[varnr] - constant) / scalar;
      }
      assert(SCIPvarIsOriginal(branchruledata->nodedata[nodeID]->vars[varnr]));
   }

   /* transform bound changes applied after the first dual reduction */
   for(varnr = 0; varnr < branchruledata->nodedata[nodeID]->nafterdualvars; varnr++)
   {
      SCIP_Real constant;
      SCIP_Real scalar;

      scalar = 1;
      constant = 0;

      if(!SCIPvarIsOriginal(branchruledata->nodedata[nodeID]->afterdualvars[varnr]))
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&branchruledata->nodedata[nodeID]->afterdualvars[varnr], &scalar, &constant)) ;
         branchruledata->nodedata[nodeID]->afterdualvarbounds[varnr] = (branchruledata->nodedata[nodeID]->afterdualvarbounds[varnr] - constant) / scalar;
      }
      assert(SCIPvarIsOriginal(branchruledata->nodedata[nodeID]->afterdualvars[varnr]));
   }

   return SCIP_OKAY;
}

/**
 * save the LPI state
 */
static
SCIP_RETCODE saveLPIstate(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_NODE*            node,
   int                   nodeID
)
{
   SCIP_LPI* lpi;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(branchruledata->savelpbasis);
   assert(nodeID >= 0);
   assert(branchruledata->nodedata[nodeID] != NULL );

   if( SCIPgetNCutsApplied(scip) > 0 )
      return SCIP_OKAY;

   /** save lpistate */
   SCIP_CALL( SCIPgetLPI(scip, &lpi) );

   /** check if the lpistate is already set */
   if(branchruledata->nodedata[nodeID]->lpistate != NULL )
   {
      SCIPlpiFreeState(lpi, SCIPblkmem(scip), &branchruledata->nodedata[nodeID]->lpistate);
   }

   SCIP_CALL( SCIPlpiGetState(lpi, SCIPblkmem(scip), &branchruledata->nodedata[nodeID]->lpistate) );

   assert(branchruledata->nodedata[nodeID]->lpistate != NULL );
   assert(SCIPisLPSolBasic(scip) && scip->lp->primalfeasible);

#ifdef DEBUG_MODE
   {
      int ncols, nrows;

      SCIP_CALL(SCIPlpiGetNRows(lpi, &nrows));
      SCIP_CALL(SCIPlpiGetNCols(lpi, &ncols));

      printf(
            ">>>> saved LPI state %p in SCIP %p, id %d, #%"SCIP_LONGINT_FORMAT", solstat %d, LPI nrows %d, ncols %d\n",
            (void*) (branchruledata->nodedata[nodeID]->lpistate), (void*) scip, nodeID,
            SCIPnodeGetNumber(node), SCIPgetLPSolstat(scip), nrows, ncols);
   }
#endif

   return SCIP_OKAY;
}

/**
 * save ancestor branching information
 */
static
SCIP_RETCODE saveAncestorBranchings(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_NODE*            node,
   SCIP_NODE*            parent,
   int                   nodeID,
   int                   parentID
)
{
   int nbranchvars;

   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(node != NULL );
   assert(parent != NULL );
   assert(nodeID >= 1);
   assert(branchruledata->nodedata[nodeID] != NULL );
   assert(parentID == 0 || branchruledata->nodedata[parentID] != NULL ); /* if the root is the next saved node, the nodedata can be NULL */

   /* allocate memory */
   if (branchruledata->nodedata[nodeID]->allocmem == 0)
   {
      assert(branchruledata->nodedata[nodeID]->vars == NULL );
      assert(branchruledata->nodedata[nodeID]->varbounds == NULL );
      assert(branchruledata->nodedata[nodeID]->varboundtypes == NULL );

      /** Allocate memory for node information */
      branchruledata->nodedata[nodeID]->allocmem = 10;
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[nodeID]->vars), branchruledata->nodedata[nodeID]->allocmem) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[nodeID]->varbounds), branchruledata->nodedata[nodeID]->allocmem) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[nodeID]->varboundtypes), branchruledata->nodedata[nodeID]->allocmem) );
   }

   assert(branchruledata->nodedata[nodeID]->allocmem > 0);
   assert(branchruledata->nodedata[nodeID]->nvars == 0);

   SCIPnodeGetAncestorBranchingsReopt(node, parent, branchruledata->nodedata[nodeID]->vars,
         branchruledata->nodedata[nodeID]->varbounds,
         branchruledata->nodedata[nodeID]->varboundtypes, &nbranchvars,
         branchruledata->nodedata[nodeID]->allocmem);

   if( nbranchvars >  branchruledata->nodedata[nodeID]->allocmem )
   {
      SCIP_CALL( checkMemory(scip, branchruledata, nodeID, nbranchvars) );

      SCIPnodeGetAncestorBranchingsReopt(node, parent, branchruledata->nodedata[nodeID]->vars,
               branchruledata->nodedata[nodeID]->varbounds,
               branchruledata->nodedata[nodeID]->varboundtypes, &nbranchvars,
               branchruledata->nodedata[nodeID]->allocmem);
   }

   assert(nbranchvars <= branchruledata->nodedata[nodeID]->allocmem); /* this should be the case */

   branchruledata->nodedata[nodeID]->nvars = nbranchvars;

   assert(nbranchvars <= branchruledata->nodedata[nodeID]->allocmem);
   assert(branchruledata->nodedata[nodeID]->vars != NULL );

   return SCIP_OKAY;
}

/*
 * save bound changes made after dual operations, e.g. strong branching.
 */
static
SCIP_RETCODE saveAfterDualBranchings(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_NODE*            node,
   int                   nodeID
)
{
   int nbranchvars;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(node != NULL);
   assert(0 < nodeID && nodeID < branchruledata->allocmemsizenodedata);
   assert(branchruledata->nodedata[nodeID] != NULL );

   /* allocate memory */
   if (branchruledata->nodedata[nodeID]->afterdualallocmem == 0)
   {
      assert(branchruledata->nodedata[nodeID]->afterdualvars == NULL );
      assert(branchruledata->nodedata[nodeID]->afterdualvarbounds == NULL );
      assert(branchruledata->nodedata[nodeID]->afterdualvarboundtypes == NULL );

      /** Allocate memory for node information */
      branchruledata->nodedata[nodeID]->afterdualallocmem = 10;
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[nodeID]->afterdualvars), branchruledata->nodedata[nodeID]->afterdualallocmem) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[nodeID]->afterdualvarbounds), branchruledata->nodedata[nodeID]->afterdualallocmem) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[nodeID]->afterdualvarboundtypes), branchruledata->nodedata[nodeID]->afterdualallocmem) );
   }

   assert(branchruledata->nodedata[nodeID]->afterdualallocmem > 0);
   assert(branchruledata->nodedata[nodeID]->nafterdualvars >= 0);

   SCIPnodeGetAfterDualBranchingsReopt(node,
         branchruledata->nodedata[nodeID]->afterdualvars,
         branchruledata->nodedata[nodeID]->afterdualvarbounds,
         branchruledata->nodedata[nodeID]->afterdualvarboundtypes,
         branchruledata->nodedata[nodeID]->nafterdualvars,
         &nbranchvars,
         branchruledata->nodedata[nodeID]->afterdualallocmem - branchruledata->nodedata[nodeID]->nafterdualvars);

   if( nbranchvars > branchruledata->nodedata[nodeID]->afterdualallocmem - branchruledata->nodedata[nodeID]->nafterdualvars)
   {
      int newsize;
      newsize = branchruledata->nodedata[nodeID]->afterdualallocmem - branchruledata->nodedata[nodeID]->nafterdualvars + nbranchvars + 1;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(branchruledata->nodedata[nodeID]->afterdualvars), newsize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(branchruledata->nodedata[nodeID]->afterdualvarbounds), newsize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(branchruledata->nodedata[nodeID]->afterdualvarboundtypes), newsize) );
      branchruledata->nodedata[nodeID]->afterdualallocmem = newsize;

      SCIPnodeGetAfterDualBranchingsReopt(node,
            branchruledata->nodedata[nodeID]->afterdualvars,
            branchruledata->nodedata[nodeID]->afterdualvarbounds,
            branchruledata->nodedata[nodeID]->afterdualvarboundtypes,
            branchruledata->nodedata[nodeID]->nafterdualvars,
            &nbranchvars,
            branchruledata->nodedata[nodeID]->afterdualallocmem - branchruledata->nodedata[nodeID]->nafterdualvars);
   }

   assert(nbranchvars <= branchruledata->nodedata[nodeID]->afterdualallocmem); /* this should be the case */

   branchruledata->nodedata[nodeID]->nafterdualvars = nbranchvars;

   assert(nbranchvars <= branchruledata->nodedata[nodeID]->afterdualallocmem);
   assert(branchruledata->nodedata[nodeID]->afterdualvars != NULL );

   return SCIP_OKAY;
}

/*
 * delete the LPI state and all LPI states below in the tree
 */
static
SCIP_RETCODE deleteLPIstate(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   SCIP_Bool             recursive
)
{
   SCIP_LPI* lpi;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(branchruledata->savelpbasis);
   assert(branchruledata->nodedata[nodeID] != NULL );
   assert(nodeID > 0);

   /** save lpistate */
   SCIP_CALL( SCIPgetLPI(scip, &lpi) );

   /** check if the lpistate is already set */
   if(branchruledata->nodedata[nodeID]->lpistate != NULL )
   {
      SCIPlpiFreeState(lpi, SCIPblkmem(scip), &branchruledata->nodedata[nodeID]->lpistate);
   }

   /** delete all LPI states below */
   if(recursive && branchruledata->nodedata[nodeID]->nodechilds != NULL )
   {
      int nchilds;
      int curchild;
      int childID;

      nchilds = SCIPqueueNElems(branchruledata->nodedata[nodeID]->nodechilds);
      curchild = 1;

      while(curchild <= nchilds)
      {
         childID = (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[nodeID]->nodechilds);
         assert(branchruledata->nodedata[childID] != NULL );

         /** delete LP state recursive */
         SCIP_CALL( deleteLPIstate(scip, branchruledata, childID, recursive) );

         /** insert childID to 'nodechilds' of nodeID */
         SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->nodechilds, (void* ) (size_t ) childID) );

         curchild++;
      }

      assert(SCIPqueueNElems(branchruledata->nodedata[nodeID]->nodechilds) == nchilds);
   }

   return SCIP_OKAY;
}

/**
 * find redundant constraints and delete the LPI state, if some
 * redundant constraints are found
 */
static
SCIP_RETCODE findRedundantConss(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   SCIP_Bool             deleteLPI
)
{
   SCIP_Bool redundant;
   int nredundantconss;
   int noldconss;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(nodeID > 0);
   assert(branchruledata->nodedata[nodeID] != NULL );

   redundant = FALSE;
   noldconss = SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss);
   nredundantconss = 0;

   if(branchruledata->nodedata[nodeID]->conss != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss))
   {
      int nconss;
      int consnr;

      nconss = SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss);
      for(consnr = 0; consnr < nconss; consnr++)
      {
         LOGICORDATA* consdata;
         consdata = (LOGICORDATA*) SCIPqueueRemove(branchruledata->nodedata[nodeID]->conss);

         assert(consdata->nvars > 1);

         /** if the number of variables inside the constraint is greater than nvars + 2, then the
          *  constraint can not be redundant
          */
         if(consdata->nvars - branchruledata->nodedata[nodeID]->nvars < 2)
         {
            int nfixedvars;
            int ndiffvars;
            int varnr;
            int consvarnr;
            nfixedvars = 0;
            ndiffvars = 0;

            for(varnr = 0; varnr < (int) branchruledata->nodedata[nodeID]->nvars; varnr++)
            {
               for(consvarnr = 0; consvarnr < (int) consdata->nvars; consvarnr++)
               {
                  if(strcmp(SCIPvarGetName(branchruledata->nodedata[nodeID]->vars[varnr]), SCIPvarGetName(consdata->vars[consvarnr])) == 0)
                  {
                     if(SCIPisFeasEQ(scip, branchruledata->nodedata[nodeID]->varbounds[varnr], consdata->vals[consvarnr]))
                        nfixedvars++;
                     else
                     {
                        ndiffvars++;
                        goto REDUNDANT;
                     }
                  }

               }
            }

            /** the constraint is not redundant */
            if(consdata->nvars - nfixedvars > 1)
            {
               SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->conss, (void* ) (size_t ) consdata) );
            }
            else
            {
               REDUNDANT:

#ifdef DEBUG_MODE
               {
                  printf("... redundant constraint found at ID %d\n", nodeID);
               }
#endif

               assert(consdata->nvars - nfixedvars <= 1 || ndiffvars >= 1);
               redundant = TRUE;
               nredundantconss++;

               /** free redundant data */
               SCIPfreeMemoryArray(scip, &consdata->vars);
               SCIPfreeMemoryArray(scip, &consdata->vals);
               SCIPfreeMemory(scip, &consdata);
            }
         }
         else
         {
            SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->conss, (void* ) (size_t ) consdata) );
         }
      }

      assert(noldconss == nredundantconss + SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss));

      if(redundant)
      {
         SCIP_CALL( deleteLPIstate(scip, branchruledata, nodeID, TRUE) );
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE saveGlobalCons(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_NODE*            node,
   REOPT_CONSTYPE        constype
)
{
   SCIP_Real constant;
   SCIP_Real scalar;
   int nbinvars;
   int varnr;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(branchruledata->sepaglbsols || branchruledata->sepaglbsubtree);
   assert(SCIPgetCurrentNode(scip) == node);

   /** realloc memory if necessary */
   if(branchruledata->nconsglobal == branchruledata->allocmemglobalcons)
   {
      SCIP_CALL( resizeGlobalsConsData(scip, branchruledata) );
   }

   nbinvars = SCIPgetNBinVars(scip);

   /** ensure that all variables are binary */
   if( nbinvars > 0 )
   {
      /** allocate memory at position nconsglobal and increase the counter */
      SCIP_CALL( SCIPallocMemory(scip, &branchruledata->globalcons[branchruledata->nconsglobal]) );
      branchruledata->globalcons[branchruledata->nconsglobal]->constype = constype;
      branchruledata->globalcons[branchruledata->nconsglobal]->nvars = 0;

      if( constype == REOPT_CONSTYPE_SEPASOLUTION && branchruledata->sepaglbsols)
      {
         SCIP_VAR** vars;

         SCIP_CALL( SCIPduplicateMemoryArray(scip, &vars, SCIPgetOrigVars(scip), SCIPgetNOrigVars(scip)) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->globalcons[branchruledata->nconsglobal]->vars, nbinvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->globalcons[branchruledata->nconsglobal]->vals, nbinvars) );

         for(varnr = 0; varnr < SCIPgetNOrigVars(scip); varnr++)
         {
            if( SCIPvarGetType(vars[varnr]) == SCIP_VARTYPE_BINARY )
            {
               int pos;

               scalar = 1;
               constant = 0;

               pos = branchruledata->globalcons[branchruledata->nconsglobal]->nvars;
               branchruledata->globalcons[branchruledata->nconsglobal]->vars[pos] = vars[varnr];

               if(!SCIPvarIsOriginal(branchruledata->globalcons[branchruledata->nconsglobal]->vars[pos]) )
               {
                  SCIP_CALL( SCIPvarGetOrigvarSum( &branchruledata->globalcons[branchruledata->nconsglobal]->vars[pos], &scalar, &constant) );
               }
               assert(SCIPvarGetStatus(branchruledata->globalcons[branchruledata->nconsglobal]->vars[pos]) == SCIP_VARSTATUS_ORIGINAL);

               branchruledata->globalcons[branchruledata->nconsglobal]->vals[pos] = (branchruledata->globalcons[branchruledata->nconsglobal]->vals[pos] - constant) * scalar;

               branchruledata->globalcons[branchruledata->nconsglobal]->nvars++;
            }
         }

         SCIPfreeMemoryArray(scip, &vars);
      }
      else if( constype == REOPT_CONSTYPE_INFSUBTREE && branchruledata->sepaglbsubtree )
      {
         SCIP_BOUNDTYPE* boundtypes;

         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->globalcons[branchruledata->nconsglobal]->vars, nbinvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->globalcons[branchruledata->nconsglobal]->vals, nbinvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &boundtypes, nbinvars) );

         SCIPnodeGetAncestorBranchingsReopt(node, SCIPgetRootNode(scip),
               branchruledata->globalcons[branchruledata->nconsglobal]->vars,
               branchruledata->globalcons[branchruledata->nconsglobal]->vals,
               boundtypes,
               &branchruledata->globalcons[branchruledata->nconsglobal]->nvars,
               nbinvars);

         for(varnr = 0; varnr < branchruledata->globalcons[branchruledata->nconsglobal]->nvars; varnr++)
         {
            scalar = 1;
            constant = 0;

            if(!SCIPvarIsOriginal(branchruledata->globalcons[branchruledata->nconsglobal]->vars[varnr]) )
            {
               SCIP_CALL( SCIPvarGetOrigvarSum( &branchruledata->globalcons[branchruledata->nconsglobal]->vars[varnr], &scalar, &constant) );
            }
            assert(SCIPvarGetStatus(branchruledata->globalcons[branchruledata->nconsglobal]->vars[varnr]) == SCIP_VARSTATUS_ORIGINAL);

            branchruledata->globalcons[branchruledata->nconsglobal]->vals[varnr] = (branchruledata->globalcons[branchruledata->nconsglobal]->vals[varnr] - constant) / scalar;

            assert(SCIPisFeasEQ(scip, branchruledata->globalcons[branchruledata->nconsglobal]->vals[varnr], 0)
                || SCIPisFeasEQ(scip, branchruledata->globalcons[branchruledata->nconsglobal]->vals[varnr], 1));
         }

         if( boundtypes != NULL )
            SCIPfreeMemoryArray(scip, &boundtypes);

         branchruledata->ninfsubtrees++;

      }

      assert(branchruledata->globalcons[branchruledata->nconsglobal]->vars != NULL );
      assert(branchruledata->globalcons[branchruledata->nconsglobal]->vals != NULL );

      branchruledata->nconsglobal++;
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE saveLocalConsData(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_NODE*            node,
   int                   nodeID,
   REOPT_CONSTYPE        constype
)
{
   SCIP_Real constant;
   SCIP_Real scalar;
   int varnr;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(SCIPgetCurrentNode(scip) == node);
   assert(SCIPnodeGetReoptID(node) > 0);

   /** ensure that all variables are binary */
   if( SCIPgetNIntVars(scip) == 0 )
   {
      SCIP_VAR** vars;
      LOGICORDATA* consdata;

      /* allocate memory */
      SCIP_CALL( SCIPallocMemory(scip, &consdata) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->vars, SCIPgetNOrigBinVars(scip)) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->vals, SCIPgetNOrigBinVars(scip)) );
      consdata->nvars = 0;
      consdata->constype = constype;

      vars = SCIPgetOrigVars(scip);

      for(varnr = 0; varnr < SCIPgetNOrigVars(scip); varnr++)
      {
         scalar = 1;
         constant = 0;

         if( SCIPvarGetType(vars[varnr]) == SCIP_VARTYPE_BINARY )
         {
            SCIP_Real lb;
            SCIP_Real ub;

            lb = SCIPvarGetLbLocal(SCIPvarGetTransVar(vars[varnr]));
            ub = SCIPvarGetUbLocal(SCIPvarGetTransVar(vars[varnr]));

            if( SCIPisFeasLT(scip, lb, ub) )
            {
               consdata->vars[consdata->nvars] = vars[varnr];

               if(!SCIPvarIsOriginal(consdata->vars[consdata->nvars]) )
               {
                  SCIP_CALL( SCIPvarGetOrigvarSum( &consdata->vars[consdata->nvars], &scalar, &constant) );
               }

               consdata->vals[consdata->nvars] = (consdata->vals[consdata->nvars] - constant) * scalar;

               consdata->nvars++;
            }
         }
      }

      if( consdata->nvars > 0 )
      {
         if( branchruledata->nodedata[nodeID]->conss == NULL )
         {
            SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[nodeID]->conss, 2, 2) );
         }

         SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->conss, (void*) consdata) );
      }
   }

   return SCIP_OKAY;
}

static
int lengthBranchPath(
   SCIP_NODE*            node,
   SCIP_BRANCHRULEDATA*  branchruledata
)
{
   int length;

   assert(node != 0);
   assert(branchruledata != NULL);
   assert(branchruledata->nodedata != NULL);

   if( SCIPnodeGetDepth(node) == 0 )
      return 0;
   else
   {
      SCIP_NODE* parent;

      assert(SCIPnodeGetReoptID(node) >= 0);
      assert(branchruledata->nodedata[SCIPnodeGetReoptID(node)] != NULL);

      parent = SCIPnodeGetParent(node);
      length = branchruledata->nodedata[SCIPnodeGetReoptID(node)]->nvars;
      while(SCIPnodeGetDepth(parent) != 0)
      {
         if( SCIPnodeGetReopttype(parent) >= SCIP_REOPTTYPE_TRANSIT )
         {
            assert(branchruledata->nodedata[SCIPnodeGetReoptID(parent)] != NULL);
            length += branchruledata->nodedata[SCIPnodeGetReoptID(parent)]->nvars;
         }
         parent = SCIPnodeGetParent(parent);
      }
   }

   return length;
}

static
int lengthBranchPathByID(
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID
)
{
   int length;
   int parentID;

   assert(branchruledata != NULL);
   assert(branchruledata->nodedata != NULL);
   assert(nodeID >= 0);
   assert(nodeID < branchruledata->allocmemsizenodedata);

   length = 0;

   if( nodeID == 0 )
      return length;
   else
   {
      parentID = nodeID;

      assert(branchruledata->nodedata[parentID] != NULL);

      while( parentID != 0)
      {
         length += branchruledata->nodedata[parentID]->nvars;
         parentID = branchruledata->nodedata[parentID]->parentID;
      }
   }

   return length;
}

/**
 * add a infeasible constraint, because the problem is proven to be infeasible
 */
static
SCIP_RETCODE addInfeasibleConstraint(
   SCIP*                 scip
)
{
   SCIP_CONS* cons;

   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "infeasibleproblem", 0, NULL, NULL, 1, 0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}


/**
 * add all saved global valid constraints
 */
SCIP_RETCODE SCIPbranchruleNodereoptAddGlobalCons(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int consnr;
   int varnr;

   assert(scip != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL );

   if(!branchruledata->init || !branchruledata->reopt )
      return SCIP_OKAY;

   assert(branchruledata->init);
   assert(branchruledata->reopt);

   /** add global constraints if same exists */
   if( branchruledata->nconsglobal > 0 && !branchruledata->infeasibleconsadded )
      printf("add %d global cons.\n", branchruledata->nconsglobal);

   for(consnr = 0; consnr < branchruledata->nconsglobal; consnr++)
   {
      SCIP_CONS* cons;
      const char* consname;
      assert(branchruledata->globalcons[consnr]->vars != NULL );

      if( !branchruledata->infeasibleproblem )
      {
         for(varnr = 0; varnr < branchruledata->globalcons[consnr]->nvars; varnr++)
         {
            if( SCIPisFeasEQ(scip, branchruledata->globalcons[consnr]->vals[varnr], 1) )
            {
               assert(SCIPvarIsOriginal(branchruledata->globalcons[consnr]->vars[varnr]));
               SCIP_CALL( SCIPgetNegatedVar(scip, branchruledata->globalcons[consnr]->vars[varnr], &branchruledata->globalcons[consnr]->vars[varnr]) );
               assert(SCIPvarIsNegated(branchruledata->globalcons[consnr]->vars[varnr]));
            }
         }

         if( branchruledata->globalcons[consnr]->constype == REOPT_CONSTYPE_SEPASOLUTION )
            consname = "globalsepasolution";
         else if( branchruledata->globalcons[consnr]->constype == REOPT_CONSTYPE_INFSUBTREE )
            consname = "globalinfsubtree";
         else
         {
            assert( branchruledata->globalcons[consnr]->constype == REOPT_CONSTYPE_STRBRANCHED );
            consname = "globalstrbranch";
         }

         SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, consname,
               branchruledata->globalcons[consnr]->nvars, branchruledata->globalcons[consnr]->vars,
               FALSE, FALSE, TRUE, branchruledata->sepabestsol, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE));

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         branchruledata->nconsglobaladded++;
      }

      SCIPfreeMemoryArray(scip, &branchruledata->globalcons[consnr]->vars);
      SCIPfreeMemoryArray(scip, &branchruledata->globalcons[consnr]->vals);
      SCIPfreeMemory(scip, &branchruledata->globalcons[consnr]);
   }

   branchruledata->nconsglobal = 0;

   return SCIP_OKAY;
}

/**
 * fix all variables from a logic-or constraint
 */
static
SCIP_RETCODE fixVarsInNode(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_NODE*            node,
   LOGICORDATA*          consdata,
   SCIP_Bool             saveinnodedata,
   SCIP_Bool             negate,
   int                   nodeID,
   SCIP_Bool             checkredundancy
)
{
   SCIP_BOUNDTYPE boundtype;
   SCIP_VAR* var;
   SCIP_Real val;
   int varnr;
   int varpos;

   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(consdata != NULL );
   assert(consdata->nvars > 0);
   assert(node != NULL );
   assert(!saveinnodedata || nodeID > 0);

   SCIPdebugMessage(" -> fix variables (save in nodedata: %d)\n", saveinnodedata);

   /** allocate memory if necessary */
   if( saveinnodedata )
   {
      if(branchruledata->nodedata[nodeID]->allocmem == 0)
      {
         assert(nodeID > 0);
         assert(branchruledata->nodedata[nodeID]->vars == NULL );
         assert(branchruledata->nodedata[nodeID]->varbounds == NULL );
         assert(branchruledata->nodedata[nodeID]->varboundtypes == NULL );

         branchruledata->nodedata[nodeID]->allocmem = consdata->nvars;
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeID]->vars, branchruledata->nodedata[nodeID]->allocmem) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeID]->varbounds, branchruledata->nodedata[nodeID]->allocmem) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeID]->varboundtypes, branchruledata->nodedata[nodeID]->allocmem) );
      }

      /* check memory */
      SCIP_CALL( checkMemory(scip, branchruledata, nodeID, branchruledata->nodedata[nodeID]->nvars+consdata->nvars) );
   }

   for(varnr = 0; varnr < consdata->nvars; varnr++)
   {
      SCIP_Bool bndchgd;

      var = consdata->vars[varnr];
      val = consdata->vals[varnr];
      boundtype = SCIPisFeasEQ(scip, val, 1) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;
      SCIP_CALL(SCIPvarGetProbvarBound(&var, &val, &boundtype));
      assert(SCIPvarIsTransformedOrigvar(var));

      if(negate)
      {
         val = 1 - val;
         boundtype = (SCIP_BOUNDTYPE) (SCIP_BOUNDTYPE_UPPER - boundtype);
      }

      bndchgd = FALSE;

      if(boundtype == SCIP_BOUNDTYPE_LOWER
      && SCIPisGT(scip, val, SCIPvarGetLbLocal(var))
      && SCIPisFeasLE(scip, val, SCIPvarGetUbLocal(var)))
      {
         SCIP_CALL(SCIPchgVarLbNode(scip, node, var, val));
         bndchgd = TRUE;
      }
      else if(boundtype == SCIP_BOUNDTYPE_UPPER
           && SCIPisLT(scip, val, SCIPvarGetUbLocal(var))
           && SCIPisFeasGE(scip, val, SCIPvarGetLbLocal(var)))
      {
         SCIP_CALL(SCIPchgVarUbNode(scip, node, var, val));
         bndchgd = TRUE;
      }
      else if(boundtype != SCIP_BOUNDTYPE_LOWER && boundtype != SCIP_BOUNDTYPE_UPPER)
      {
         printf("** Unknown boundtype: %d **\n", boundtype);
         assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);
      }

      SCIPdebugMessage("  <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", val);

      /** add variable and bound to branching path information, because we don't want to delete this data */
      if( saveinnodedata && bndchgd )
      {
         SCIP_Real constant;
         SCIP_Real scalar;

         constant = 0;
         scalar = 1;
         varpos = branchruledata->nodedata[nodeID]->nvars;

         branchruledata->nodedata[nodeID]->vars[varpos] = var;

         if (!SCIPvarIsOriginal(branchruledata->nodedata[nodeID]->vars[varpos]))
         {
            SCIP_CALL( SCIPvarGetOrigvarSum(&branchruledata->nodedata[nodeID]->vars[varpos], &scalar, &constant) );
         }
         assert(SCIPvarIsOriginal(branchruledata->nodedata[nodeID]->vars[varpos]));
         branchruledata->nodedata[nodeID]->varbounds[varpos] = (val - constant) * scalar;
         branchruledata->nodedata[nodeID]->varboundtypes[varpos] = (SCIPisFeasEQ(scip, branchruledata->nodedata[nodeID]->varbounds[varpos], 0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
         branchruledata->nodedata[nodeID]->nvars++;
      }
   }

   if(nodeID > 0 && checkredundancy && branchruledata->nodedata[nodeID]->conss != NULL
   && SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss))
   {
      SCIP_CALL( findRedundantConss(scip, branchruledata, nodeID, branchruledata->savelpbasis) );
   }

   return SCIP_OKAY;
}

/**
 * add a given logic-or constraint or fix the corresponding variable
 */
static
SCIP_RETCODE addLocalConsToNode(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_NODE*            node,
   LOGICORDATA*          consdata,
   SCIP_Bool             checkredundancy,
   int                   nodeID
)
{
   assert(branchruledata != NULL );
   assert(node != NULL );
   assert(branchruledata->reopt);
   assert(consdata != NULL );
   assert(consdata->nvars > 0);
   assert(consdata->vars != NULL );
   assert(consdata->vals != NULL );

   /** if the constraints consists of only one variable, we can fix the variable */
   if(consdata->nvars > 1)
   {
      SCIP_CONS* cons;
      SCIP_VAR** vars;
      SCIP_Bool initial;
      SCIP_Bool sepa;
      const char* consname;
      int varnr;

      /** copy variables and negate them if necessary */
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &vars, consdata->vars, consdata->nvars) );
      for(varnr = 0; varnr < consdata->nvars; varnr++)
      {
         SCIP_Real val;
         SCIP_BOUNDTYPE boundtype;

         assert(SCIPvarIsOriginal(vars[varnr]));

         val = consdata->vals[varnr];
         boundtype = SCIPisFeasEQ(scip, val, 1) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;
         SCIP_CALL(SCIPvarGetProbvarBound(&vars[varnr], &val, &boundtype));
         assert(SCIPvarIsTransformed(vars[varnr]));

         if( SCIPisFeasEQ(scip, val, 1) )
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, vars[varnr], &vars[varnr]) );
            assert(SCIPvarIsNegated(vars[varnr]));
         }
      }

      /* if we want to save LP basis, the constraints need to not initial */
      if( consdata->constype == REOPT_CONSTYPE_SEPASOLUTION )
      {
         consname = "sepasolution";
         initial = FALSE;
         sepa = FALSE;
      }
      else if( consdata->constype == REOPT_CONSTYPE_INFSUBTREE )
      {
         consname = "infsubtree";
         initial = TRUE;
         sepa = FALSE;
      }
      else
      {
         assert( consdata->constype == REOPT_CONSTYPE_STRBRANCHED);
         consname = "strbranch";
         initial = TRUE;
         sepa = FALSE;
      }

      /** add old constraints to child1 */
      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, consname, consdata->nvars, vars,
            initial, sepa, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE));

      SCIPdebugMessage(" -> added constraint\n");
      SCIPdebugPrintCons(scip, cons, NULL);

      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIPfreeMemoryArray(scip, &vars);
   }
   else
   {
      SCIP_CALL( fixVarsInNode(scip, branchruledata, node, consdata, nodeID > 0 ? TRUE : FALSE, TRUE, nodeID, checkredundancy) );
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE getLastSavedNode(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_NODE*            node,
   SCIP_NODE**           parent,
   int*                  parentID,
   int*                  nbndchgs
)
{
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   (*nbndchgs) = 0;
   (*parent) = node;

   /** look for a saved parent along the root-path */
   while((*parent) != SCIPgetRootNode(scip))
   {
      (*nbndchgs) += SCIPnodeGetNDomchg(*parent, TRUE, TRUE, FALSE);
      (*parent) = SCIPnodeGetParent(*parent);

      if((*parent) == SCIPgetRootNode(scip))
      {
         (*parentID) = 0;
         break;
      }
      else if(SCIPnodeGetReopttype((*parent)) >= SCIP_REOPTTYPE_TRANSIT)
      {
         assert(SCIPnodeGetReoptID((*parent)) < branchruledata->allocmemsizenodedata);
         (*parentID) = SCIPnodeGetReoptID((*parent));
         break;
      }
   }

   return SCIP_OKAY;
}

/**
 * delete the data for node nodeID in nodedata
 */
static
SCIP_RETCODE deleteNodeData(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   SCIP_Bool             exitsolve
)
{
   assert(branchruledata != NULL );
   assert(branchruledata->nodedata[nodeID] != NULL );

   SCIPdebugMessage("delete node at ID %d, hard remove: %u\n", nodeID, exitsolve);

   if (exitsolve)
   {
      /** delete data for constraints */
      if (branchruledata->nodedata[nodeID]->conss != NULL )
      {
         while (!SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss))
         {
            LOGICORDATA* consdata = SCIPqueueRemove(branchruledata->nodedata[nodeID]->conss);
            SCIPfreeMemoryArray(scip, &consdata->vals);
            SCIPfreeMemoryArray(scip, &consdata->vars);
            SCIPfreeMemory(scip, &consdata);
         }
         SCIPqueueFree(&branchruledata->nodedata[nodeID]->conss);
      }

      /* free list of children */
      if (branchruledata->nodedata[nodeID]->nodechilds != NULL )
      {
         SCIPqueueFree(&branchruledata->nodedata[nodeID]->nodechilds);
      }

      /* delete dual constraint */
      if( branchruledata->nodedata[nodeID]->dualcons != NULL )
      {
         if( branchruledata->nodedata[nodeID]->dualcons[0] != NULL )
         {
            SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->dualcons[0]->vals);
            SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->dualcons[0]->vars);
            SCIPfreeMemory(scip, &branchruledata->nodedata[nodeID]->dualcons[0]);
         }

         if( branchruledata->nodedata[nodeID]->dualcons[1] != NULL )
         {
            SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->dualcons[1]->vals);
            SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->dualcons[1]->vars);
            SCIPfreeMemory(scip, &branchruledata->nodedata[nodeID]->dualcons[1]);
         }
         SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->dualcons);
      }

      /* free boundtypes */
      if (branchruledata->nodedata[nodeID]->varboundtypes != NULL )
      {
         SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->varboundtypes);
      }

      /* free bounds */
      if (branchruledata->nodedata[nodeID]->varbounds != NULL )
      {
         SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->varbounds);
      }

      /* free variables */
      if (branchruledata->nodedata[nodeID]->vars != NULL )
      {
         SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->vars);
      }

      /* free afterdual-boundtypes */
      if (branchruledata->nodedata[nodeID]->afterdualvarboundtypes != NULL )
      {
         SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->afterdualvarboundtypes);
      }

      /* free afterdual-bounds */
      if (branchruledata->nodedata[nodeID]->afterdualvarbounds != NULL )
      {
         SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->afterdualvarbounds);
      }

      /* free afterdual-variables */
      if (branchruledata->nodedata[nodeID]->afterdualvars != NULL )
      {
         SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->afterdualvars);
      }

      SCIPfreeMemory(scip, &branchruledata->nodedata[nodeID]);
   }
   else
   {
      /** remove and delete all constraints */
      if (branchruledata->nodedata[nodeID]->conss != NULL )
      {
         while (!SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss))
         {
            LOGICORDATA* consdata = SCIPqueueRemove(branchruledata->nodedata[nodeID]->conss);
            SCIPfreeMemoryArray(scip, &consdata->vals);
            SCIPfreeMemoryArray(scip, &consdata->vars);
            SCIPfreeMemory(scip, &consdata);
         }
      }

      /* remove all children */
      if (branchruledata->nodedata[nodeID]->nodechilds != NULL )
      {
         SCIPqueueClear(branchruledata->nodedata[nodeID]->nodechilds);
      }

      branchruledata->nodedata[nodeID]->nvars = 0;
      branchruledata->nodedata[nodeID]->pseudobranched = FALSE;
      branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_NONE;
   }

   assert(branchruledata->nodedata[nodeID] == NULL
        || branchruledata->nodedata[nodeID]->conss == NULL
        || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss));
   assert(branchruledata->nodedata[nodeID] == NULL
       || branchruledata->nodedata[nodeID]->nodechilds == NULL
       || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds));

   branchruledata->nsavednodes--;

   return SCIP_OKAY;
}

/*
 * Clear the data structure after branching under the root node.
 * This will be delete all all information of pruned nodes and reset the statistics for the next pricing round.
 */
static
SCIP_RETCODE clearNodes(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_Bool             exitsolve
)
{
   int nodeID;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   /** clear queue with open IDs */
   SCIPqueueClear(branchruledata->openIDs);

   /** delete all data about nodes */
   for(nodeID = 0; nodeID < branchruledata->allocmemsizenodedata; nodeID++)
   {
      if(branchruledata->nodedata[nodeID] != NULL )
      {
         SCIP_CALL( deleteNodeData(scip, branchruledata, nodeID, TRUE) );
         assert(branchruledata->nodedata[nodeID] == NULL || branchruledata->nodedata[nodeID]->nvars == 0);
         assert(branchruledata->nodedata[nodeID] == NULL || branchruledata->nodedata[nodeID]->allocmem <= SCIPgetNOrigVars(scip));
      }

      if(nodeID > 0 && !exitsolve)
      {
         SCIP_CALL( SCIPqueueInsert(branchruledata->openIDs, (void* ) (size_t ) nodeID) );
      }
   }
   assert(exitsolve || SCIPqueueNElems(branchruledata->openIDs) == branchruledata->allocmemsizenodedata - 1);

   return SCIP_OKAY;
}

/**
 * change bound along the root path
 */
static
SCIP_RETCODE changeAncestorBranchings(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   SCIP_NODE*            child1,
   SCIP_NODE*            child2
)
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_BOUNDTYPE* boundtypes;
   int varnr;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(branchruledata->nodedata[nodeID] != NULL );
   assert(branchruledata->nodedata[nodeID]->conss != NULL || branchruledata->nodedata[nodeID]->vars != NULL );
   assert(branchruledata->nodedata[nodeID]->conss != NULL || branchruledata->nodedata[nodeID]->varbounds != NULL );
   assert(branchruledata->nodedata[nodeID]->conss != NULL || branchruledata->nodedata[nodeID]->varboundtypes != NULL );
   assert(child1 != NULL );

   /** copy memory to ensure that only original variables are saved */
   if( branchruledata->nodedata[nodeID]->nvars == 0
    && branchruledata->nodedata[nodeID]->nafterdualvars == 0)
      return SCIP_OKAY;

   /**
    * Chance bounds due to all ancestors.
    * Fix all decisions that were applied before the first dual reduction.
    */
   if( branchruledata->nodedata[nodeID]->nvars > 0 )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &vars, branchruledata->nodedata[nodeID]->vars, branchruledata->nodedata[nodeID]->nvars) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &vals, branchruledata->nodedata[nodeID]->varbounds, branchruledata->nodedata[nodeID]->nvars) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &boundtypes, branchruledata->nodedata[nodeID]->varboundtypes, branchruledata->nodedata[nodeID]->nvars) );

      SCIPdebugMessage(" -> change ancestor branchings:\n");

      for(varnr = 0; varnr < branchruledata->nodedata[nodeID]->nvars; varnr++)
      {
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Real newbound;

         assert(SCIPvarIsOriginal(vars[varnr]));

         SCIP_CALL( SCIPvarGetProbvarBound(&vars[varnr], &vals[varnr], &boundtypes[varnr]) );

         assert(SCIPvarIsTransformed(vars[varnr]));

         oldlb = SCIPvarGetLbLocal(vars[varnr]);
         oldub = SCIPvarGetUbLocal(vars[varnr]);
         newbound = vals[varnr];

         if(boundtypes[varnr] == SCIP_BOUNDTYPE_LOWER
         && SCIPisGT(scip, newbound, oldlb)
         && SCIPisFeasLE(scip, newbound, oldub))
         {
            SCIP_CALL( SCIPchgVarLbNode(scip, child1, vars[varnr], newbound) );

            if(child2 != NULL )
            {
               SCIP_CALL( SCIPchgVarLbNode(scip, child2, vars[varnr], newbound) );
            }
         }
         else if(boundtypes[varnr] == SCIP_BOUNDTYPE_UPPER
              && SCIPisLT(scip, newbound, oldub)
              && SCIPisFeasGE(scip, newbound, oldlb))
         {
            SCIP_CALL( SCIPchgVarUbNode(scip, child1, vars[varnr], newbound) );

            if(child2 != NULL )
            {
               SCIP_CALL( SCIPchgVarUbNode(scip, child2, vars[varnr], newbound) );
            }
         }
         else if(boundtypes[varnr] != SCIP_BOUNDTYPE_LOWER
              && boundtypes[varnr] != SCIP_BOUNDTYPE_UPPER)
         {
            printf("** Unknown boundtype: %d **\n", boundtypes[varnr]);
            assert(boundtypes[varnr] == SCIP_BOUNDTYPE_LOWER || boundtypes[varnr] == SCIP_BOUNDTYPE_UPPER);
         }

         SCIPdebugMessage("  <%s> %s %g\n", SCIPvarGetName(vars[varnr]), boundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=", newbound);
      }

      /** free memory */
      SCIPfreeMemoryArray(scip, &vars);
      SCIPfreeMemoryArray(scip, &vals);
      SCIPfreeMemoryArray(scip, &boundtypes);
   }

   /**
    * Chance bounds due to all ancestors.
    * Fix all decisions that were applied after the first dual reduction.
    *
    * Additionally, all these bound changes can be transfered to "normal" bound changes
    * because the node is split.
    */
   if( branchruledata->nodedata[nodeID]->nafterdualvars > 0 )
   {
      assert(branchruledata->nodedata[nodeID]->reopttype != SCIP_REOPTTYPE_INFSUBTREE);

      /* reallocate memory to transfer the after-dual bound changes */
      SCIP_CALL( checkMemory(scip, branchruledata, nodeID, branchruledata->nodedata[nodeID]->nvars + branchruledata->nodedata[nodeID]->nafterdualvars) );

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &vars, branchruledata->nodedata[nodeID]->afterdualvars, branchruledata->nodedata[nodeID]->nafterdualvars) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &vals, branchruledata->nodedata[nodeID]->afterdualvarbounds, branchruledata->nodedata[nodeID]->nafterdualvars) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &boundtypes, branchruledata->nodedata[nodeID]->afterdualvarboundtypes, branchruledata->nodedata[nodeID]->nafterdualvars) );

      SCIPdebugMessage(" -> change ancestor branchings (after dual):\n");

      for(varnr = 0; varnr < branchruledata->nodedata[nodeID]->nafterdualvars; varnr++)
      {
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Real newbound;

         /* transfer the bound change */
         assert(SCIPvarIsOriginal(branchruledata->nodedata[nodeID]->afterdualvars[varnr]));
         branchruledata->nodedata[nodeID]->vars[branchruledata->nodedata[nodeID]->nvars] = branchruledata->nodedata[nodeID]->afterdualvars[varnr];
         branchruledata->nodedata[nodeID]->varbounds[branchruledata->nodedata[nodeID]->nvars] = branchruledata->nodedata[nodeID]->afterdualvarbounds[varnr];
         branchruledata->nodedata[nodeID]->varboundtypes[branchruledata->nodedata[nodeID]->nvars] = branchruledata->nodedata[nodeID]->afterdualvarboundtypes[varnr];
         branchruledata->nodedata[nodeID]->nvars++;

         assert(SCIPvarIsOriginal(vars[varnr]));

         SCIP_CALL( SCIPvarGetProbvarBound(&vars[varnr], &vals[varnr], &boundtypes[varnr]) );

         assert(SCIPvarIsTransformed(vars[varnr]));

         oldlb = SCIPvarGetLbLocal(vars[varnr]);
         oldub = SCIPvarGetUbLocal(vars[varnr]);
         newbound = vals[varnr];

         if(boundtypes[varnr] == SCIP_BOUNDTYPE_LOWER
         && SCIPisGT(scip, newbound, oldlb)
         && SCIPisFeasLE(scip, newbound, oldub))
         {
            SCIP_CALL( SCIPchgVarLbNode(scip, child1, vars[varnr], newbound) );
         }
         else if(boundtypes[varnr] == SCIP_BOUNDTYPE_UPPER
              && SCIPisLT(scip, newbound, oldub)
              && SCIPisFeasGE(scip, newbound, oldlb))
         {
            SCIP_CALL( SCIPchgVarUbNode(scip, child1, vars[varnr], newbound) );
         }
         else if(boundtypes[varnr] != SCIP_BOUNDTYPE_LOWER
              && boundtypes[varnr] != SCIP_BOUNDTYPE_UPPER)
         {
            printf("** Unknown boundtype: %d **\n", boundtypes[varnr]);
            assert(boundtypes[varnr] == SCIP_BOUNDTYPE_LOWER || boundtypes[varnr] == SCIP_BOUNDTYPE_UPPER);
         }

         SCIPdebugMessage("  <%s> %s %g\n", SCIPvarGetName(vars[varnr]), boundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=", newbound);
      }

      /** free memory */
      SCIPfreeMemoryArray(scip, &vars);
      SCIPfreeMemoryArray(scip, &vals);
      SCIPfreeMemoryArray(scip, &boundtypes);

      /* mark the after-dual variables deleted */
      branchruledata->nodedata[nodeID]->nafterdualvars = 0;
   }

   return SCIP_OKAY;
}

/*
 * move all children to the next saved node
 *
 * if a node get branched and has neither an added constraints nor
 * exceed the number of bound changes the compared to the next saved
 * node the parameter 'maxdiffofnodes'.
 */
static
SCIP_RETCODE moveChildrenUp(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_NODE*            parentnode,
   int                   nodeID,
   int                   parentID
)
{
   int childID;
   int varnr;
   int nvars;

   assert(branchruledata != NULL);
   assert(nodeID >= 1);
   assert(parentID >= 0);
   assert(branchruledata->nodedata[nodeID]->nodechilds != NULL);

   while(!SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds))
   {
      childID = (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[nodeID]->nodechilds);

      /* check the memory */
      SCIP_CALL( checkMemory(scip, branchruledata, childID, branchruledata->nodedata[childID]->nvars + branchruledata->nodedata[nodeID]->nvars) );
      assert(branchruledata->nodedata[childID]->allocmem >= branchruledata->nodedata[childID]->nvars + branchruledata->nodedata[nodeID]->nvars);

      /** save branching information */
      for(varnr = 0; varnr < branchruledata->nodedata[nodeID]->nvars; varnr++)
      {
         nvars = branchruledata->nodedata[childID]->nvars;
         branchruledata->nodedata[childID]->vars[nvars] = branchruledata->nodedata[nodeID]->vars[varnr];
         branchruledata->nodedata[childID]->varbounds[nvars] = branchruledata->nodedata[nodeID]->varbounds[varnr];
         branchruledata->nodedata[childID]->varboundtypes[nvars] = branchruledata->nodedata[nodeID]->varboundtypes[varnr];
         branchruledata->nodedata[childID]->nvars++;
      }

      /* update the ID of the parent node */
      branchruledata->nodedata[childID]->parentID = parentID;

      /* insert the node as a child */
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->nodechilds, (void* ) (size_t ) childID));
   }

   return SCIP_OKAY;
}

/*
 * generate local constraints and add them
 */
static
SCIP_RETCODE genLocalCons(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   SCIP_NODE*            child1,
   SCIP_NODE*            child2,
   SCIP_Bool             lastalreadyadded
)
{
   int ncons;
   int curcons;
   int varnr;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(child1 != NULL );

   assert(SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss) > 0);
   ncons = SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss);
   curcons = 1;

   /** logic-or cons */
   while(curcons <= ncons - (int) lastalreadyadded)
   {
      LOGICORDATA* consdata;
      SCIP_CONS* cons1;
      SCIP_CONS* cons2;
      SCIP_VAR** vars;
      const char* consname;

      consdata = (LOGICORDATA*) (size_t) SCIPqueueRemove(branchruledata->nodedata[nodeID]->conss);
      assert(consdata->nvars > 1);

      /** copy variables and negate them if necessary */
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &vars, consdata->vars, consdata->nvars) );
      for(varnr = 0; varnr < consdata->nvars; varnr++)
      {
         SCIP_Real val;
         SCIP_BOUNDTYPE boundtype;

         assert(SCIPvarIsOriginal(vars[varnr]));

         val = consdata->vals[varnr];
         boundtype = SCIPisFeasEQ(scip, val, 1) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;
         SCIP_CALL( SCIPvarGetProbvarBound(&vars[varnr], &val, &boundtype) );
         assert(SCIPvarIsTransformed(vars[varnr]));

         if ( SCIPisFeasEQ(scip, val, 1) )
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, vars[varnr], &vars[varnr]) );
            assert(SCIPvarIsNegated(vars[varnr]));
         }
      }

      if( consdata->constype == REOPT_CONSTYPE_SEPASOLUTION )
         consname = "sepasolution";
      else if( consdata->constype == REOPT_CONSTYPE_INFSUBTREE )
      {
         consname = "infsubtree";
      }
      else
      {
         assert( consdata->constype == REOPT_CONSTYPE_STRBRANCHED);
         consname = "strbranch";
      }

      /** add old constraints to child1 */
      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons1, consname, consdata->nvars, vars,
            FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE));

      SCIP_CALL( SCIPaddConsNode(scip, child1, cons1, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons1) );

      /** add old constraints to child2, if child was PB */
      if(child2 != NULL )
      {
         SCIP_CALL( SCIPcreateConsLogicor(scip, &cons2, consname, consdata->nvars, vars,
               FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE));

         SCIP_CALL( SCIPaddConsNode(scip, child2, cons2, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
      }

      SCIPfreeMemoryArray(scip, &vars);

      curcons++;
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->conss, (void* ) (size_t ) consdata) );
   }

   SCIPdebugMessage(" -> add %d constraints to %lld\n", ncons - (int) lastalreadyadded, SCIPnodeGetNumber(child1));
   if( child2 != NULL )
   {
      SCIPdebugMessage(" -> add %d constraints to %lld\n", ncons - (int) lastalreadyadded, SCIPnodeGetNumber(child2));
   }

   return SCIP_OKAY;
}

/* delete all nodes in the by nodeID induced subtree */
static
SCIP_RETCODE deleteChildrenBelow(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   SCIP_Bool             delnodeitself
)
{
   assert(branchruledata != NULL );
   assert(nodeID >= 0);
   assert(branchruledata->init && branchruledata->reopt);
   assert(branchruledata->nodedata[nodeID] != NULL );

   /** delete all children below */
   if (branchruledata->nodedata[nodeID]->nodechilds != NULL )
   {
      while (!SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds))
      {
         SCIP_CALL( deleteChildrenBelow(scip, branchruledata, (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[nodeID]->nodechilds), TRUE) );
      }
   }

   /** delete node data*/
   if( delnodeitself )
   {
      SCIP_CALL(deleteNodeData(scip, branchruledata, nodeID, TRUE));
      SCIP_CALL(SCIPqueueInsert(branchruledata->openIDs, (void*) (size_t) nodeID));
   }

   return SCIP_OKAY;
}

/*
 * heuristic that calculates a largest induced subtree (LR)
 */
static
SCIP_RETCODE findLR(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_VAR**            varsLR1,
   SCIP_VAR**            varsLR2,
   SCIP_Real*            valsLR1,
   SCIP_Real*            valsLR2,
   SCIP_BOUNDTYPE*       boundsLR1,
   SCIP_BOUNDTYPE*       boundsLR2,
   int                   nallocs,
   int*                  nvarsLR1,
   int*                  nvarsLR2,
   SCIP_Real*            lossLR1,
   SCIP_Real*            lossLR2,
   int*                  SetToCompress,
   int                   nnodesToCompress,
   int*                  LR1,
   int*                  LR2,
   int*                  nnodesLR1,
   int*                  nnodesLR2
)
{
   SCIP_Bool* commonvars1;
   SCIP_Bool* commonvars2;
   SCIP_Bool* commonvars1_best;
   SCIP_Bool* commonvars2_best;
   SCIP_Bool* nodevars;
   const char** varnames;
   int* nvars;
   int node;
   int parentID;
   int var;

   SCIP_Real tmp_lossLR1;
   SCIP_Real tmp_lossLR2;
   SCIP_Real tmp_minloss;
   SCIP_Real minloss;
   int startnode;

   int tmp_nnodesLR1;
   int tmp_nnodesLR2;
   int tmp_nvarsLR1;
   int tmp_nvarsLR2;

   assert(scip != NULL);
   assert(branchruledata != NULL);

   /* start time */
   SCIP_CALL( SCIPstartClock(scip, branchruledata->lrtime) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &commonvars1, 2*SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &commonvars2, 2*SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &commonvars1_best, 2*SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &commonvars2_best, 2*SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nodevars, 2*SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &varnames, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nvars, nnodesToCompress) );

   *nnodesLR1 = 0;
   *nnodesLR2 = 0;
   *nvarsLR1 = 0;
   *nvarsLR2 = 0;
   *lossLR1 = SCIPinfinity(scip);
   *lossLR2 = SCIPinfinity(scip);
   minloss = SCIPinfinity(scip);

   for(startnode = 0; startnode < nnodesToCompress; startnode++)
   {
      tmp_nnodesLR1 = 0;
      tmp_nnodesLR2 = 0;
      tmp_nvarsLR1 = 0;
      tmp_nvarsLR2 = 0;

      /*
       * clear bool arrays
       * the array has structure: [x1=0, x1=1, x2=0, x2=1, ... , xn=0, xn=1]
       * */
      for(var = 0; var < 2*SCIPgetNVars(scip); var++)
      {
         commonvars1[var] = FALSE;
         commonvars2[var] = FALSE;
      }

      for(node = 0; node < nnodesToCompress; node++)
         nvars[node] = 0;

      /* initialize array with the first node */
      parentID = SetToCompress[startnode];
      assert(branchruledata->nodedata[parentID] != NULL);
      assert(branchruledata->nodedata[parentID]->nvars >= 1);

      while( parentID != 0 )
      {
         assert( branchruledata->nodedata[parentID] != NULL);
         for(var = 0; var < branchruledata->nodedata[parentID]->nvars; var++)
         {
            if( SCIPisFeasEQ(scip, branchruledata->nodedata[parentID]->varbounds[var], 0))
            {
               commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])] = TRUE;
               (tmp_nvarsLR1)++;
            }
            else
            {
               commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])+1] = TRUE;
               (tmp_nvarsLR1)++;
            }
            varnames[SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])] = SCIPvarGetName(branchruledata->nodedata[parentID]->vars[var]);
         }

         /* get strong branched vars */
         if( branchruledata->nodedata[parentID]->pseudobranched )
         {
            LOGICORDATA* consdata;
            int nvarscons;

            assert(branchruledata->nodedata[parentID]->dualcons[0] != NULL);
            assert(branchruledata->nodedata[parentID]->dualcons[0]->nvars > 0);

            for(var = 0; var < branchruledata->nodedata[parentID]->dualcons[0]->nvars; var++)
            {
               if( SCIPisFeasEQ(scip, branchruledata->nodedata[parentID]->dualcons[0]->vals[var], 0) )
               {
                  commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualcons[0]->vars[var])] = TRUE;
                  (tmp_nvarsLR1)++;
               }
               else
               {
                  commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualcons[0]->vars[var])+1] = TRUE;
                  (tmp_nvarsLR1)++;
               }
               varnames[SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualcons[0]->vars[var])] = SCIPvarGetName(branchruledata->nodedata[parentID]->dualcons[0]->vars[var]);
            }
         }

         /* go to the next saved node above */
         parentID = branchruledata->nodedata[parentID]->parentID;
      }

      /* add the first node to LR1 */
      LR1[tmp_nnodesLR1] = SetToCompress[startnode];
      (tmp_nnodesLR1)++;

      /* set length of root path */
      nvars[startnode] = (tmp_nvarsLR1);

      tmp_lossLR1 = nvars[startnode];
      tmp_lossLR2 = nvars[startnode];

      /* iterate over all other nodes */
      for(node = 0; node < nnodesToCompress; node++)
      {
         if( node != startnode )
         {
            int ncommonvars;

            parentID = SetToCompress[node];
            ncommonvars = 0;

            /* set length of root path */
            nvars[node] = 0;

            /* clear bool array */
            for(var = 0; var < 2*SCIPgetNVars(scip); var++)
               nodevars[var] = FALSE;

            /* check if the intersection is empty and delete variables which are set to a different bound in the current node */
            while( parentID != 0 )
            {
               assert( branchruledata->nodedata[parentID] != NULL);
               for(var = 0; var < branchruledata->nodedata[parentID]->nvars; var++)
               {
                  assert(branchruledata->nodedata[parentID]->vars != NULL);
                  assert(branchruledata->nodedata[parentID]->varbounds != NULL);

                  if( !nodevars[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])]
                   && !nodevars[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])+1] )
                  {
                     if( SCIPisFeasEQ(scip, branchruledata->nodedata[parentID]->varbounds[var], 0) )
                     {
                        /* add the variable to the var array */
                        nodevars[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])] = TRUE;

                        if( commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])] )
                           ncommonvars++;
                     }
                     else
                     {
                        /* add the variable to the var array */
                        nodevars[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])+1] = TRUE;

                        if( commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])+1] )
                           ncommonvars++;
                     }

                     /* add variables name that are not part of the branching path of the first seen node */
                     if( varnames[SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])] == NULL )
                        varnames[SCIPvarGetIndex(branchruledata->nodedata[parentID]->vars[var])] = SCIPvarGetName(branchruledata->nodedata[parentID]->vars[var]);

                     /* set length of root path */
                     nvars[node]++;
                  }
               }

               /* get strong branched vars */
               if( branchruledata->nodedata[parentID]->pseudobranched )
               {
                  LOGICORDATA* consdata;
                  int nvarscons;

                  assert(branchruledata->nodedata[parentID]->dualcons[0] != NULL);
                  assert(branchruledata->nodedata[parentID]->dualcons[0]->nvars > 0);

                  for(var = 0; var < branchruledata->nodedata[parentID]->dualcons[0]->nvars; var++)
                  {
                     if( SCIPisFeasEQ(scip, branchruledata->nodedata[parentID]->dualcons[0]->vals[var], 0) )
                     {
                        /* add the variable to the var array */
                        nodevars[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualcons[0]->vars[var])] = TRUE;

                        if( commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualcons[0]->vars[var])] )
                           ncommonvars++;
                     }
                     else
                     {
                        /* add the variable to the var array */
                        nodevars[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualcons[0]->vars[var])+1] = TRUE;

                        if( commonvars1[2*SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualcons[0]->vars[var])+1] )
                           ncommonvars++;
                     }
                     if( varnames[SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualcons[0]->vars[var])] == NULL )
                        varnames[SCIPvarGetIndex(branchruledata->nodedata[parentID]->dualcons[0]->vars[var])] = SCIPvarGetName(branchruledata->nodedata[parentID]->dualcons[0]->vars[var]);

                     /* set length of root path */
                     nvars[node]++;
                  }
               }

               /* go to the next saved node above */
               parentID = branchruledata->nodedata[parentID]->parentID;
            }

            /* delete all variable in commonvars1 and commonvars2 that are not part of vars, respectively. */
            for(var = 0; var < 2*SCIPgetNVars(scip); var++)
            {
               if( ncommonvars > 0 && commonvars1[var] && !nodevars[var] )
               {
                  commonvars1[var] = FALSE;
                  tmp_nvarsLR1--;
                  assert(ncommonvars <= tmp_nvarsLR1);
               }
               /* this is the first node in LR2 */
               else if( ncommonvars == 0 && tmp_nnodesLR2 == 0 && nodevars[var])
               {
                  commonvars2[var] = TRUE;
                  (tmp_nvarsLR2)++;
               }
               else if( ncommonvars == 0 && commonvars2[var] && !nodevars[var] )
               {
                  commonvars2[var] = FALSE;
                  tmp_nvarsLR2--;
                  assert(tmp_nnodesLR2 >= 0);
               }
            }

            assert(ncommonvars == 0 || ncommonvars == tmp_nvarsLR1);

            /* the intersection is empty */
            if( ncommonvars == 0 )
            {
               LR2[tmp_nnodesLR2] = SetToCompress[node];
               (tmp_nnodesLR2)++;
            }
            else
            {
               LR1[tmp_nnodesLR1] = SetToCompress[node];
               (tmp_nnodesLR1)++;
            }

            tmp_lossLR1 += nvars[node];
            tmp_lossLR2 += nvars[node];
         }
      }

      assert(tmp_nnodesLR1 + tmp_nnodesLR2 == nnodesToCompress);

      /* calculate loss of LR1 and LR2 */
      for(node = 0; node < tmp_nnodesLR1; node++)
         (tmp_lossLR1) -= (tmp_nvarsLR1);

      for(node = 0; node < *nnodesLR2; node++)
         (tmp_lossLR2) -= (tmp_nvarsLR2);

      if( tmp_nnodesLR1 > 0 && tmp_nnodesLR2 > 0 )
         tmp_minloss = MIN(tmp_lossLR1/tmp_nnodesLR1, tmp_lossLR2/tmp_nnodesLR2);
      else if( tmp_nnodesLR1 > 0 )
         tmp_minloss = tmp_lossLR1/tmp_nnodesLR1;
      else
         tmp_minloss = tmp_lossLR2/tmp_nnodesLR2;

      if( tmp_minloss < minloss )
      {
         *lossLR1 = tmp_lossLR1;
         *lossLR2 = tmp_lossLR2;

         *nnodesLR1 = tmp_nnodesLR1;
         *nnodesLR2 = tmp_nnodesLR2;

         /* update minloss */
         minloss = tmp_minloss;

         /* collect data for LR1 */
         *nvarsLR1 = 0;
         *nvarsLR2 = 0;
         for(var = 0; var < 2*SCIPgetNVars(scip); var++)
         {
            if( commonvars1[var] )
            {
               varsLR1[*nvarsLR1] = SCIPfindVar(scip, varnames[(int)((SCIP_Real)var/2)]);
               valsLR1[*nvarsLR1] =  var%2 == 0 ? 0 : 1;
               boundsLR1[*nvarsLR1] = (SCIP_BOUNDTYPE) (SCIPisFeasEQ(scip, valsLR1[*nvarsLR1], 0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);

               assert(varsLR1[*nvarsLR1] != NULL);
               (*nvarsLR1)++;
            }

            if( varsLR2 != NULL && commonvars2[var] )
            {
               varsLR2[*nvarsLR2] = SCIPfindVar(scip, varnames[(int)((SCIP_Real)var/2)]);
               valsLR2[*nvarsLR2] =  var%2 == 0 ? 0 : 1;
               boundsLR2[*nvarsLR2] = (SCIP_BOUNDTYPE) (SCIPisFeasEQ(scip, valsLR2[*nvarsLR2], 0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);

               assert(varsLR2[*nvarsLR2] != NULL);
               (*nvarsLR2)++;
            }
         }
      }
   }

   /* free memory */
   SCIPfreeBlockMemoryArray(scip, &commonvars1, 2*SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &commonvars2, 2*SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &nodevars, 2*SCIPgetNVars(scip));
   SCIPfreeMemoryArray(scip, &varnames);
   SCIPfreeBlockMemoryArray(scip, &nvars, nnodesToCompress);

   branchruledata->lrcalls++;

   /* stop time */
   SCIP_CALL( SCIPstopClock(scip, branchruledata->lrtime) );

   return SCIP_OKAY;
}

/*
 * generate weak compression
 */
static
SCIP_RETCODE genWC(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_VAR**            vars,
   SCIP_Real*            vals,
   SCIP_BOUNDTYPE*       bounds,
   LOGICORDATA**         conss,
   int                   nvars,
   int                   nconss
)
{
   LOGICORDATA* consdata;
   int nodeIDtobranch;
   int nodeIDleaf;
   int var;
   int c;
   int parentID;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(bounds != NULL);
   assert(conss != NULL);
   assert(nvars > 0);
   assert(nconss >= 0);

   /*
    * fix variables and create node with a corresponding logic-or constraint
    */
   assert(!SCIPqueueIsEmpty(branchruledata->openIDs));

   nodeIDtobranch = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
   assert(branchruledata->nodedata[nodeIDtobranch] == NULL || branchruledata->nodedata[nodeIDtobranch]->nvars == 0);

   /* initialize node */
   SCIP_CALL( initNode(scip, branchruledata, nodeIDtobranch) );
   branchruledata->nodedata[nodeIDtobranch]->reopttype = SCIP_REOPTTYPE_TRANSIT;
   branchruledata->nodedata[nodeIDtobranch]->parentID = 0;

   /* allocate memory for branching decisions */
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &branchruledata->nodedata[nodeIDtobranch]->vars, vars, nvars) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &branchruledata->nodedata[nodeIDtobranch]->varbounds, vals, nvars) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &branchruledata->nodedata[nodeIDtobranch]->varboundtypes, bounds, nvars) );
   branchruledata->nodedata[nodeIDtobranch]->allocmem = nvars;
   branchruledata->nodedata[nodeIDtobranch]->nvars = nvars;

   /* allocate queue for child nodes */
   SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[nodeIDtobranch]->nodechilds, 2, 2) );

   nodeIDleaf = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
   assert(branchruledata->nodedata[nodeIDleaf] == NULL || branchruledata->nodedata[nodeIDleaf]->nvars == 0);

   /* initialize node */
   SCIP_CALL( initNode(scip, branchruledata, nodeIDleaf) );
   branchruledata->nodedata[nodeIDleaf]->reopttype = nvars > 1 ? SCIP_REOPTTYPE_LOGICORNODE : SCIP_REOPTTYPE_LEAF;
   branchruledata->nodedata[nodeIDleaf]->parentID = 0;

   if( nvars > 1 )
   {
      /* create queue for local constaints */
      SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[nodeIDleaf]->conss, 1, 1) );

      /* create and add the constraint */
      SCIP_CALL( SCIPallocMemory(scip, &consdata) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vars, vars, nvars) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vals, vals, nvars) );
      consdata->nvars = nvars;
      consdata->constype = REOPT_CONSTYPE_STRBRANCHED;
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeIDleaf]->conss, (void*) consdata) );
   }
   else
   {
      assert(nvars == 1);

      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeIDleaf]->vars, 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeIDleaf]->varbounds, 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeIDleaf]->varboundtypes, 1) );
      branchruledata->nodedata[nodeIDleaf]->allocmem = 1;
      branchruledata->nodedata[nodeIDleaf]->nvars = 1;

      branchruledata->nodedata[nodeIDleaf]->vars[0] = vars[0];
      branchruledata->nodedata[nodeIDleaf]->varbounds[0] = 1 - vals[0];
      branchruledata->nodedata[nodeIDleaf]->varboundtypes[0] = (SCIP_BOUNDTYPE) ( SCIP_BOUNDTYPE_UPPER - bounds[0]);
   }

   /* add both nodes as child nodes to the root node */
   if( branchruledata->nodedata[0]->nodechilds == NULL )
   {
      SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[0]->nodechilds, 2, 2) );
   }
   SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[0]->nodechilds, (void*) (size_t) nodeIDtobranch) );
   SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[0]->nodechilds, (void*) (size_t) nodeIDleaf) );

   parentID = nodeIDtobranch;

   /*
    * create child nodes with constraints and corresponding fixings
    */
   for(c = 0; c < nconss; c++)
   {
      LOGICORDATA* consdatabranch;

      /*
       * create node fixings (will be leaf node)
       */
      nodeIDleaf = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
      assert(branchruledata->nodedata[nodeIDleaf] == NULL || branchruledata->nodedata[nodeIDleaf]->nvars == 0);

      /* initialize the node */
      SCIP_CALL( initNode(scip, branchruledata, nodeIDleaf) );
      branchruledata->nodedata[nodeIDleaf]->reopttype = SCIP_REOPTTYPE_LEAF;
      branchruledata->nodedata[nodeIDleaf]->parentID = parentID;

      /* allocate memory for the fixings */
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeIDleaf]->vars, conss[c]->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeIDleaf]->varbounds, conss[c]->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeIDleaf]->varboundtypes, conss[c]->nvars) );
      branchruledata->nodedata[nodeIDleaf]->allocmem = conss[c]->nvars;
      branchruledata->nodedata[nodeIDleaf]->nvars = conss[c]->nvars;

      /* copy the variables, bound, and boundtypes */
      for(var = 0; var < conss[c]->nvars; var++)
      {
         branchruledata->nodedata[nodeIDleaf]->vars[var] = conss[c]->vars[var];
         branchruledata->nodedata[nodeIDleaf]->varbounds[var] = conss[c]->vals[var];
         branchruledata->nodedata[nodeIDleaf]->varboundtypes[var] = ( SCIPisFeasEQ(scip, conss[c]->vals[var], 1) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER );
      }

      /* add the node as a child of the last branched node */
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->nodechilds, (void*) (size_t) nodeIDleaf) );

      /*
       * create node with constraint (need to be branched)
       */
      nodeIDtobranch = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
      assert(branchruledata->nodedata[nodeIDtobranch] == NULL || branchruledata->nodedata[nodeIDtobranch]->nvars == 0);

      /* initialize node */
      SCIP_CALL( initNode(scip, branchruledata, nodeIDtobranch) );
      if( branchruledata->nodedata[nodeIDtobranch]->nodechilds == NULL )
      {
         SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[nodeIDtobranch]->nodechilds, 2, 2) );
      }
      branchruledata->nodedata[nodeIDtobranch]->reopttype = SCIP_REOPTTYPE_LOGICORNODE;
      branchruledata->nodedata[nodeIDtobranch]->parentID = parentID;

      /* create queue for local constraints */
      SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[nodeIDtobranch]->conss, 1, 1) );

      /* add the constraint */
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeIDtobranch]->conss, (void*) conss[c]) );

      /* add the node as a child of the last branched node */
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->nodechilds, (void*) (size_t) nodeIDtobranch) );

      /* set the current node as the next node to branch */
      parentID = nodeIDtobranch;
   }

   branchruledata->wclastnodeID = nodeIDtobranch;

   return SCIP_OKAY;
}



/*
 * reduce the tree to feasibility
 */
static
SCIP_RETCODE genLC(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_VAR***           vars,
   SCIP_Real**           bounds,
   SCIP_BOUNDTYPE**      boundtypes,
   int                   nnodes,
   int*                  nvars,
   SCIP_QUEUE**          conss,
   int                   newparentID,
   int*                  nodeID_cons,
   SCIP_Real*            avgdepth
)
{
   int consID;
   int nodeID;
   int var;
   int newID;

   assert(scip != NULL);
   assert(branchruledata != NULL);

   /* start time */
   SCIP_CALL( SCIPstartClock(scip, branchruledata->wctime) );

   branchruledata->wccalls++;
   (*avgdepth) = 0;

   /* ensure that parent node is allocated */
   assert(branchruledata->nodedata[newparentID] != NULL);
   if( branchruledata->nodedata[newparentID]->nodechilds == NULL )
   {
      SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[newparentID]->nodechilds, nnodes+1, 2) );
   }

   /* save node with added constraints */
   if (SCIPqueueIsEmpty(branchruledata->openIDs))
   {
      SCIP_CALL(reallocNodedata(scip, branchruledata));
   }

   newID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
   assert(newID >= 1);
   assert(newID < branchruledata->allocmemsizenodedata);
   assert(branchruledata->nodedata[newID] == NULL
      || (branchruledata->nodedata[newID]->nvars == 0
          && (branchruledata->nodedata[newID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[newID]->conss))) );

   if( branchruledata->nodedata[newID] == NULL )
   {
      SCIP_CALL( initNode(scip, branchruledata, newID) );
   }

   /* set the root as parent node*/
   branchruledata->nodedata[newID]->parentID = newparentID;

   /* set the reopttype */
   branchruledata->nodedata[newID]->reopttype = SCIP_REOPTTYPE_LOGICORNODE;

   *nodeID_cons = newID;
   consID = newID;

   /* increase number of saved nodes */
   branchruledata->nsavednodes++;

   if( branchruledata->nodedata[newID]->conss == NULL )
   {
      SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[newID]->conss, nnodes, 2) );
   }

   /* add the node as a child of the root node */
   assert(branchruledata->nodedata[newparentID]->nodechilds != NULL);
   SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[newparentID]->nodechilds, (void*) (size_t) newID) );

   /** save the new generated nodes */

   /* store nodes */
   for(nodeID = 0; nodeID < nnodes; nodeID++)
   {
      if (SCIPqueueIsEmpty(branchruledata->openIDs))
      {
         SCIP_CALL(reallocNodedata(scip, branchruledata));
      }

      newID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
      assert(newID >= 1);
      assert(newID < branchruledata->allocmemsizenodedata);
      assert(nvars[nodeID] > 0);
      assert(vars[nodeID] != NULL);
      assert(bounds[nodeID] != NULL);
      assert(boundtypes[nodeID] != NULL);
      assert(branchruledata->nodedata[newID] == NULL
         || (branchruledata->nodedata[newID]->nvars == 0
             && (branchruledata->nodedata[newID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[newID]->conss))) );

      if( branchruledata->nodedata[newID] == NULL )
      {
         SCIP_CALL( initNode(scip, branchruledata, newID) );
      }

      /* set the parent node*/
      branchruledata->nodedata[newID]->parentID = newparentID;

      /* set the reopttype */
      branchruledata->nodedata[newID]->reopttype = SCIP_REOPTTYPE_LEAF;

      /* increase number of saved nodes */
      branchruledata->nsavednodes++;

      /* allocate memory */
      if( branchruledata->nodedata[newID]->allocmem == 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[newID]->vars, nvars[nodeID]) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[newID]->varbounds, nvars[nodeID]) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[newID]->varboundtypes, nvars[nodeID]) );
         branchruledata->nodedata[newID]->allocmem = nvars[nodeID];
      }

      assert(branchruledata->nodedata[newID]->vars != NULL);
      assert(branchruledata->nodedata[newID]->varbounds != NULL);
      assert(branchruledata->nodedata[newID]->varboundtypes != NULL);

      /* copy variable data */
      for(var = 0; var < nvars[nodeID]; var++)
      {
         branchruledata->nodedata[newID]->vars[var] = vars[nodeID][var];
         branchruledata->nodedata[newID]->varbounds[var] = bounds[nodeID][var];
         branchruledata->nodedata[newID]->varboundtypes[var] = boundtypes[nodeID][var];
         branchruledata->nodedata[newID]->nvars++;
      }

      (*avgdepth) += nvars[nodeID];

      /* copy constraint data */
      if( branchruledata->nodedata[newID]->conss == NULL && !SCIPqueueIsEmpty(conss[nodeID]) )
      {
         SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[newID]->conss, SCIPqueueNElems(conss[nodeID])+1, 2) );
      }
      assert(branchruledata->nodedata[newID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[newID]->conss));

      while( !SCIPqueueIsEmpty(conss[nodeID]) )
      {
         SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[newID]->conss, SCIPqueueRemove(conss[nodeID])) );
      }

      /* create a constraint and add them to the node at consID */
      if( nvars[nodeID] == 1 )
      {
         /* allocate memory */
         if( branchruledata->nodedata[consID]->allocmem == 0 )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[consID]->vars, nvars[nodeID]) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[consID]->varbounds, nvars[nodeID]) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[consID]->varboundtypes, nvars[nodeID]) );
            branchruledata->nodedata[consID]->allocmem = SCIPgetNOrigBinVars(scip);
         }

         /* fix the variable to the negated values */
         branchruledata->nodedata[consID]->vars[branchruledata->nodedata[consID]->nvars] = vars[nodeID][0];
         branchruledata->nodedata[consID]->varbounds[branchruledata->nodedata[consID]->nvars] = 1 - bounds[nodeID][0];
         branchruledata->nodedata[consID]->varboundtypes[branchruledata->nodedata[consID]->nvars] = ( 1 - boundtypes[nodeID][0] == 0 ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER );
         branchruledata->nodedata[consID]->nvars++;
      }
      else
      {
         LOGICORDATA* consdata;

         SCIP_CALL( SCIPallocMemory(scip, &consdata) );
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vars, vars[nodeID], nvars[nodeID]) );
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vals, bounds[nodeID], nvars[nodeID]) );
         consdata->nvars = nvars[nodeID];
         consdata->constype = REOPT_CONSTYPE_INFSUBTREE;

         /* add the constraint to the node */
         SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[consID]->conss, (void*) consdata) );
      }

      /* add the node as a child node */
      assert(branchruledata->nodedata[newparentID]->nodechilds != NULL);
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[newparentID]->nodechilds, (void*) (size_t) newID) );
   }

   (*avgdepth) /= nnodes;

   /* for each constraint that consists of exactly one variable, we fixed the variable to the negated value,
    * thus, the number of constraints + number of variables has to be nfeasnodeIDs */
   assert(nnodes == SCIPqueueNElems(branchruledata->nodedata[consID]->conss) + branchruledata->nodedata[consID]->nvars);

   /* stop time */
   SCIP_CALL( SCIPstopClock(scip, branchruledata->wctime) );

   return SCIP_OKAY;
}

/*
 * find weak compression
 */
static
SCIP_RETCODE findWC(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   SCIP_VAR**            vars,
   SCIP_Real*            vals,
   SCIP_BOUNDTYPE*       bounds,
   int*                  nvars,
   LOGICORDATA**         conss,
   int*                  nconss
)
{
   int parentID;
   int var;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(bounds != NULL);
   assert(nvars != NULL);


   parentID = nodeID;
   *nvars = 0;
   *nconss = 0;

   /* go up to the root */
   while( parentID != 0 )
   {
      assert(branchruledata->nodedata[parentID] != NULL);

      /* copy branching information */
      for(var = 0; var < branchruledata->nodedata[parentID]->nvars; var++)
      {
         vars[*nvars] = branchruledata->nodedata[parentID]->vars[var];
         vals[*nvars] = branchruledata->nodedata[parentID]->varbounds[var];
         bounds[*nvars] = branchruledata->nodedata[parentID]->varboundtypes[var];
         (*nvars) += 1;
      }

      /* collect all added constaints along the root path */
      if( branchruledata->nodedata[parentID]->conss != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[parentID]->conss) )
      {
         int ncons;
         int cons;

         ncons = SCIPqueueNElems(branchruledata->nodedata[parentID]->conss);
         cons = 0;
         while( cons < ncons )
         {
            LOGICORDATA* consdata;
            LOGICORDATA* consdataCopy;

            consdata = (LOGICORDATA*) SCIPqueueRemove(branchruledata->nodedata[parentID]->conss);

            if( consdata->constype == REOPT_CONSTYPE_STRBRANCHED ) //|| consdata->constype == REOPT_CONSTYPE_INFSUBTREE )
            {
               SCIP_CALL( SCIPallocMemory(scip, &consdataCopy) );
               SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdataCopy->vars, consdata->vars, consdata->nvars) );
               SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdataCopy->vals, consdata->vals, consdata->nvars) );
               consdataCopy->nvars = consdata->nvars;
               consdataCopy->constype = consdata->constype;

               conss[*nconss] = consdataCopy;
               (*nconss) += 1;
            }
            cons++;
            SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->conss, (void*) consdata) );
         }
      }

      parentID = branchruledata->nodedata[parentID]->parentID;
   }

   return SCIP_OKAY;
}


/*
 * find a lazy compression
 */
static
SCIP_RETCODE findLC(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int*                  SetToCompress,
   int                   nnodesToCompress,
   SCIP_VAR***           vars,
   SCIP_Real**           vals,
   SCIP_BOUNDTYPE**      bounds,
   int*                  nvars,
   SCIP_QUEUE**          conss
)
{
   int nodeID;
   int parentID;
   int var;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(nnodesToCompress >= 0);
   assert(SetToCompress != NULL || nnodesToCompress == 0);
   assert(vars != NULL || nnodesToCompress == 0);
   assert(vals != NULL || nnodesToCompress == 0);
   assert(bounds != NULL || nnodesToCompress == 0);
   assert(nvars != NULL || nnodesToCompress == 0);

   /* collect the root path */
   for(nodeID = 0; nodeID < nnodesToCompress; nodeID++)
   {
      assert(vars[nodeID] != NULL);
      assert(vals[nodeID] != NULL);
      assert(bounds[nodeID] != NULL);

      parentID = SetToCompress[nodeID];
      nvars[nodeID] = 0;

      /* go up to the root */
      while( parentID != 0 )
      {
         assert(branchruledata->nodedata[parentID] != NULL);

         /* copy branching information */
         for(var = 0; var < branchruledata->nodedata[parentID]->nvars; var++)
         {
            vars[nodeID][nvars[nodeID]] = branchruledata->nodedata[parentID]->vars[var];
            vals[nodeID][nvars[nodeID]] = branchruledata->nodedata[parentID]->varbounds[var];
            bounds[nodeID][nvars[nodeID]] = branchruledata->nodedata[parentID]->varboundtypes[var];
            nvars[nodeID]++;
         }

         /* collect all added constaints along the root path */
         if( branchruledata->nodedata[parentID]->conss != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[parentID]->conss) )
         {
            int ncons;
            int cons;

            ncons = SCIPqueueNElems(branchruledata->nodedata[parentID]->conss);
            cons = 0;
            while( cons < ncons )
            {
               LOGICORDATA* consdata;
               LOGICORDATA* consdataCopy;

               consdata = (LOGICORDATA*) SCIPqueueRemove(branchruledata->nodedata[parentID]->conss);

               if( consdata->constype == REOPT_CONSTYPE_STRBRANCHED )
               {
                  SCIP_CALL( SCIPallocMemory(scip, &consdataCopy) );
                  SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdataCopy->vars, consdata->vars, consdata->nvars) );
                  SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdataCopy->vals, consdata->vals, consdata->nvars) );
                  consdataCopy->nvars = consdata->nvars;
                  consdataCopy->constype = consdata->constype;
               }
               cons++;
               SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->conss, (void*) consdata) );
            }
         }

         /* collect all dual branching information */
         if( branchruledata->nodedata[parentID]->pseudobranched )
         {
            LOGICORDATA* consdata;
            int nvarscons;

            nvarscons = SCIPbranchrulePseudoGetNPseudoVars(scip, parentID);
            SCIP_CALL( SCIPallocMemory(scip, &consdata) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->vars, nvarscons) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->vals, nvarscons) );
            consdata->nvars = -1;
            SCIP_CALL( SCIPbranchrulePseudoGenerateCons(scip, consdata, &(consdata->nvars), nvarscons, parentID, FALSE, FALSE) );

            assert(consdata->nvars <= nvarscons);

            for(var = 0; var < consdata->nvars; var++)
            {
               vars[nodeID][nvars[nodeID]] = consdata->vars[var];
               vals[nodeID][nvars[nodeID]] = consdata->vals[var];
               bounds[nodeID][nvars[nodeID]] = ( SCIPisFeasEQ(scip, consdata->vals[var], 0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER );
               nvars[nodeID]++;
            }

            SCIPfreeMemoryArray(scip, &consdata->vals);
            SCIPfreeMemoryArray(scip, &consdata->vars);
            SCIPfreeMemory(scip, &consdata);
         }

         parentID = branchruledata->nodedata[parentID]->parentID;
      }

      assert(nvars[nodeID] > 0);
   }

   return SCIP_OKAY;
}

/*
 * generate the LR
 */
static
SCIP_RETCODE genLR(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_VAR**            vars,
   SCIP_Real*            vals,
   SCIP_BOUNDTYPE*       bounds,
   int                   nvars,
   int                   parentID,
   int*                  consnodeID
)
{
   LOGICORDATA* consdata;
   int LRnodeID;
   int notLRnodeID;
   int var;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);
   assert(nvars == 0 || bounds != NULL);
   assert(!SCIPqueueIsEmpty(branchruledata->openIDs));

   /* start time */
   SCIP_CALL( SCIPstartClock(scip, branchruledata->lrtime) );

   /** ensure that the parent node is allocated */
   assert(0 <= parentID && parentID < branchruledata->allocmemsizenodedata);
   assert(branchruledata->nodedata[parentID] != NULL);

   /** generate LR child */
   if (SCIPqueueIsEmpty(branchruledata->openIDs))
   {
      SCIP_CALL(reallocNodedata(scip, branchruledata));
   }

   LRnodeID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
   assert(branchruledata->nodedata[LRnodeID] == NULL
      || (branchruledata->nodedata[LRnodeID]->nvars == 0
          && (branchruledata->nodedata[LRnodeID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[LRnodeID]->conss))) );
   assert(LRnodeID > 0);
   assert(branchruledata->nodedata[LRnodeID] == NULL || branchruledata->nodedata[LRnodeID]->nvars == 0);

   if( branchruledata->nodedata[LRnodeID] == NULL )
   {
      SCIP_CALL( initNode(scip, branchruledata, LRnodeID) );
   }

   if(branchruledata->nodedata[LRnodeID]->allocmem == 0)
   {
      assert(branchruledata->nodedata[LRnodeID]->vars == NULL );
      assert(branchruledata->nodedata[LRnodeID]->varbounds == NULL );
      assert(branchruledata->nodedata[LRnodeID]->varboundtypes == NULL );

      /** Allocate memory for node information */
      branchruledata->nodedata[LRnodeID]->allocmem = SCIPgetNOrigVars(scip);
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[LRnodeID]->vars), branchruledata->nodedata[LRnodeID]->allocmem) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[LRnodeID]->varbounds), branchruledata->nodedata[LRnodeID]->allocmem) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[LRnodeID]->varboundtypes), branchruledata->nodedata[LRnodeID]->allocmem) );
   }

   /* check memory */
   SCIP_CALL( checkMemory(scip, branchruledata, LRnodeID, nvars) );
   assert(branchruledata->nodedata[LRnodeID]->allocmem >= nvars);

   /* copy bounds */
   for(var = 0; var < nvars; var++)
   {
      assert(vars[var] != NULL);
      branchruledata->nodedata[LRnodeID]->vars[var] = vars[var];
      branchruledata->nodedata[LRnodeID]->varbounds[var] = vals[var];
      branchruledata->nodedata[LRnodeID]->varboundtypes[var] = bounds[var];
      branchruledata->nodedata[LRnodeID]->nvars++;
   }

   /* set reopttype */
   branchruledata->nodedata[LRnodeID]->reopttype = SCIP_REOPTTYPE_LEAF;

   /* set parentID */
   branchruledata->nodedata[LRnodeID]->parentID = parentID;

   /* increase number of saved nodes */
   branchruledata->nsavednodes++;

   /** generate child with constraint */
   if (SCIPqueueIsEmpty(branchruledata->openIDs))
   {
      SCIP_CALL(reallocNodedata(scip, branchruledata));
   }

   notLRnodeID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
   assert(branchruledata->nodedata[notLRnodeID] == NULL
      || (branchruledata->nodedata[notLRnodeID]->nvars == 0
          && (branchruledata->nodedata[notLRnodeID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[notLRnodeID]->conss))) );
   assert(branchruledata->nodedata[notLRnodeID] == NULL || branchruledata->nodedata[notLRnodeID]->nvars == 0);
   assert(notLRnodeID > 0);

   if( branchruledata->nodedata[notLRnodeID] == NULL )
   {
      SCIP_CALL( initNode(scip, branchruledata, notLRnodeID) );
   }

   if( nvars > 1 )
   {
      if (branchruledata->nodedata[notLRnodeID]->conss == NULL)
      {
         SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[notLRnodeID]->conss, 2, 2) );
      }

      SCIP_CALL( SCIPallocMemory(scip, &consdata) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vars, vars, nvars) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vals, vals, nvars) );
      consdata->nvars = nvars;
      consdata->constype = REOPT_CONSTYPE_INFSUBTREE;

      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[notLRnodeID]->conss, (void*) consdata) );
   }
   else
   {
      if(branchruledata->nodedata[notLRnodeID]->allocmem == 0)
      {
         assert(branchruledata->nodedata[notLRnodeID]->vars == NULL );
         assert(branchruledata->nodedata[notLRnodeID]->varbounds == NULL );
         assert(branchruledata->nodedata[notLRnodeID]->varboundtypes == NULL );

         /** Allocate memory for node information */
         branchruledata->nodedata[notLRnodeID]->allocmem = SCIPgetNOrigVars(scip);
         SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[notLRnodeID]->vars), branchruledata->nodedata[notLRnodeID]->allocmem) );
         SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[notLRnodeID]->varbounds), branchruledata->nodedata[notLRnodeID]->allocmem) );
         SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[notLRnodeID]->varboundtypes), branchruledata->nodedata[notLRnodeID]->allocmem) );
      }

      /* check memory */
      SCIP_CALL( checkMemory(scip, branchruledata, notLRnodeID, 1) );

      assert(vars[0] != NULL);
      branchruledata->nodedata[notLRnodeID]->vars[0] = vars[0];
      branchruledata->nodedata[notLRnodeID]->varbounds[0] = 1 - vals[0];
      branchruledata->nodedata[notLRnodeID]->varboundtypes[0] = (SCIP_BOUNDTYPE) (1 - bounds[0]);
      branchruledata->nodedata[notLRnodeID]->nvars++;
   }

   /* set reopttype */
   branchruledata->nodedata[notLRnodeID]->reopttype = SCIP_REOPTTYPE_LOGICORNODE;

   /* set parentID */
   branchruledata->nodedata[notLRnodeID]->parentID = parentID;

   /* set the ID of the node with the added constraint */
   (*consnodeID) = notLRnodeID;

   /* add this two nodes as child nodes below the root */
   if( branchruledata->nodedata[parentID]->nodechilds == NULL )
   {
      SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[parentID]->nodechilds, 2, 2) );
   }

   SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->nodechilds, (void*) (size_t) LRnodeID) );
   SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->nodechilds, (void*) (size_t) notLRnodeID) );

   /* stop time */
   SCIP_CALL( SCIPstopClock(scip, branchruledata->lrtime) );

   return SCIP_OKAY;
}

/*
 * run heuristics to compress the search tree
 */
static
SCIP_RETCODE runHeuristics(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_Bool*            success
)
{
   SCIP_VAR** varsLR1;
   SCIP_VAR** varsLR2;
   SCIP_Real* valsLR1;
   SCIP_Real* valsLR2;
   SCIP_BOUNDTYPE* boundsLR1;
   SCIP_BOUNDTYPE* boundsLR2;
   SCIP_Real lossLR1;
   SCIP_Real lossLR2;
   SCIP_Real minLoss;
   SCIP_Bool successLR1;
   SCIP_Bool successLR2;
   int depth;
   int nnodesLR1;
   int nnodesLR2;
   int nodeID_cons;
   int nodeID;
   int nnodesToCompress;
   int nallocvars;
   int nvarsLR1;
   int nvarsLR2;
   int nrepresentatives;
   int allocsize;
   int* SetToCompress;
   int* LR1;
   int* LR2;
   int nodeWC;
   int fixedvarsWC;

   SCIP_VAR** varsWC;
   SCIP_Real* valsWC;
   SCIP_BOUNDTYPE* boundsWC;
   int nvarsWC;
   LOGICORDATA** conssWC;
   int nconssWC;

   assert(scip != NULL);
   assert(branchruledata != NULL);

   depth = 0;
   nnodesLR1 = 0;
   nnodesLR2 = 0;
   nodeID_cons = 0;
   nodeID = 0;
   nnodesToCompress = 0;
   nallocvars = 0;
   nvarsLR1 = 0;
   nvarsLR2 = 0;
   nrepresentatives = 0;
   allocsize = 0;

   lossLR1 = SCIPinfinity(scip);
   lossLR2 = SCIPinfinity(scip);

   nodeWC = -1;
   nvarsWC = 0;
   nconssWC = 0;
   fixedvarsWC = 1;

   /* allocate general memory */
   allocsize = branchruledata->allocmemsizenodedata - SCIPqueueNElems(branchruledata->openIDs);
   nallocvars = SCIPgetNBinVars(scip);
   SCIP_CALL( SCIPallocMemoryArray(scip, &SetToCompress, allocsize) );

   /* allocate memory for LR*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &varsLR1, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &valsLR1, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &boundsLR1, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &boundsLR2, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &varsLR2, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &valsLR2, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &LR1, allocsize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &LR2, allocsize) );

   /* allocate memory for WC */
   SCIP_CALL( SCIPallocMemoryArray(scip, &varsWC, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &valsWC, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &boundsWC, nallocvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &conssWC, 20) );

   /* collect nodeIDs */
   nnodesToCompress= 0;
   for(nodeID = 1; nodeID < branchruledata->allocmemsizenodedata; nodeID++)
   {
      if( branchruledata->nodedata[nodeID] != NULL
       && lengthBranchPathByID(branchruledata, nodeID) >= 1
       && (branchruledata->nodedata[nodeID]->nodechilds == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds)) )
      {
         if( branchruledata->cpressnodes == 0 )
         {
            SetToCompress[nnodesToCompress] = nodeID;
            nnodesToCompress++;
         }
         else if( branchruledata->cpressnodes == 1 && branchruledata->nodedata[nodeID]->reopttype != SCIP_REOPTTYPE_FEASIBLE )
         {
            SetToCompress[nnodesToCompress] = nodeID;
            nnodesToCompress++;
         }
         else if( branchruledata->cpressnodes == 2 && branchruledata->nodedata[nodeID]->reopttype == SCIP_REOPTTYPE_FEASIBLE )
         {
            SetToCompress[nnodesToCompress] = nodeID;
            nnodesToCompress++;
         }

//         if( branchruledata->wcenable && (branchruledata->nodedata[nodeID]->nodechilds == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds))
//          && branchruledata->nodedata[nodeID]->reopttype == SCIP_REOPTTYPE_PRUNED )
//         {
//            int nfixedvars;
//            nfixedvars = lengthBranchPathByID(branchruledata, nodeID);
//
//            if( nfixedvars > fixedvarsWC && nfixedvars < SCIPgetNBinVars(scip) && nfixedvars > branchruledata->wcdepth )
//            {
//               fixedvarsWC = nfixedvars;
//               nodeWC = nodeID;
//            }
//         }
      }
   }

   nrepresentatives = 0;
   nodeID_cons = 0;
   *success = FALSE;
   successLR1 = FALSE;
   successLR2 = FALSE;
   depth = 0;
   minLoss = SCIPinfinity(scip);

   /* find LR nodes */
   if( nnodesToCompress > 0 && branchruledata->lrenable )
   {
      SCIP_CALL( findLR(scip, branchruledata, varsLR1, varsLR2, valsLR1, valsLR2, boundsLR1, boundsLR2,
            nallocvars, &nvarsLR1, &nvarsLR2, &lossLR1, &lossLR2, SetToCompress, nnodesToCompress, LR1, LR2,
            &nnodesLR1, &nnodesLR2) );

      minLoss = MIN(lossLR1/nnodesLR1, lossLR2/nnodesLR2);
      depth = lossLR1/nnodesLR1 <= lossLR2/nnodesLR2 ? nvarsLR1 : nvarsLR2;

      if( branchruledata->lrenable && (minLoss < branchruledata->lrloss || depth > branchruledata->lrdepth) )
      {
         if( nnodesLR1 >= 1 && nvarsLR1 > 1 ) /* if nnnodesLR1feas = 1 then LR is equivalent to LC */
            successLR1 = TRUE;
         if( nnodesLR2 >= 1 && nvarsLR2 > 1 ) /* if nnnodesLR1feas = 1 then LR is equivalent to LC */
            successLR2 = TRUE;
      }
   }

   if( !successLR1 && !successLR2 && !branchruledata->wcenable )
      goto SKIP;

   assert(!branchruledata->lrenable || nnodesToCompress == nnodesLR1 + nnodesLR2);

   if( !successLR1 && !successLR2 && branchruledata->lastprunedID > 0 && branchruledata->wcenable)
   {
      int nfixedvars;
      nfixedvars = lengthBranchPathByID(branchruledata, branchruledata->lastprunedID);

      branchruledata->wccalls++;

      if( nfixedvars > branchruledata->wcdepth && nfixedvars < SCIPgetNBinVars(scip) )
      {
         fixedvarsWC = nfixedvars;
         nodeWC = branchruledata->lastprunedID;
      }
   }

   if( nodeWC == branchruledata->wclastnodeID )
      nodeWC = -1;


   /* find weak compression of nodes */
   if( !successLR1 && !successLR2 && branchruledata->wcenable && nodeWC > 0 )
   {
      SCIP_CALL( findWC(scip, branchruledata, nodeWC, varsWC, valsWC, boundsWC, &nvarsWC, conssWC, &nconssWC) );

//      /* ensure that the representation is much smaller than the current search front */
//      if( 3 + 2*nconssWC >= (SCIP_Real)(branchruledata->allocmemsizenodedata - SCIPqueueNElems(branchruledata->openIDs))/2 )
//         nodeWC = -1;
   }

   assert(successLR1 + (nodeWC > 0) <= 1);
   assert(successLR2 + (nodeWC > 0) <= 1);

   if( !successLR1 && !successLR2 && nodeWC <= 0 )
      goto SKIP;

   /* if at least one heuristic was successful we clear the node data */
   if( successLR1 || successLR2 || nodeWC > 0 )
   {
      /* reset the saved data */
      SCIP_CALL( clearNodes(scip, branchruledata, FALSE) );
      SCIP_CALL( SCIPbranchrulePseudoReset(scip, TRUE, FALSE) );

      /* initialize the root data */
      SCIP_CALL( initNode(scip, branchruledata, 0) );
   }

   /* compress the nodes */
   if( successLR1 || successLR2 || nodeWC > 0 )
   {
      if( successLR1 || successLR2 )
      {
         assert(branchruledata->lrenable);

         /* do the compression */
         if( successLR1 && lossLR1 <= lossLR2 )
         {
            /* generate LR1 */
            assert(branchruledata->nodedata[0] != NULL);
            assert(successLR1);

            SCIP_CALL( genLR(scip, branchruledata, varsLR1, valsLR1, boundsLR1, nvarsLR1, 0, &nodeID_cons) );

            assert(0 < nodeID_cons && nodeID_cons < branchruledata->allocmemsizenodedata);
            nrepresentatives = 2;

            /* generate LR2 */
            if( successLR2 )
            {
               assert(branchruledata->nodedata[nodeID_cons] != NULL);
               assert(successLR2);

               SCIP_CALL( genLR(scip, branchruledata, varsLR2, valsLR2, boundsLR2, nvarsLR2, nodeID_cons, &nodeID_cons) );

               assert(0 < nodeID_cons && nodeID_cons < branchruledata->allocmemsizenodedata);
               nrepresentatives++;
            }

            branchruledata->lrsuccess++;
            branchruledata->lrloss = minLoss;
            branchruledata->lrdepth = MAX(branchruledata->lrdepth, depth);

            *success = TRUE;
         }
         else if( successLR2 && lossLR1 > lossLR2 )
         {
            /* generate LR2 */
            assert(branchruledata->nodedata[0] != NULL);
            SCIP_CALL( genLR(scip, branchruledata, varsLR2, valsLR2, boundsLR2, nvarsLR2, nodeID_cons, &nodeID_cons) );
            assert(0 < nodeID_cons && nodeID_cons < branchruledata->allocmemsizenodedata);
            nrepresentatives = 2;

            /* generate LR1 if nnodeLR1 > 0 */
            if( successLR2 )
            {
               assert(branchruledata->nodedata[nodeID_cons] != NULL);
               SCIP_CALL( genLR(scip, branchruledata, varsLR1, valsLR1, boundsLR1, nvarsLR1, nodeID_cons, &nodeID_cons) );
               assert(0 < nodeID_cons && nodeID_cons < branchruledata->allocmemsizenodedata);
               nrepresentatives++;
            }

            branchruledata->lrsuccess++;
            branchruledata->lrloss = minLoss;
            branchruledata->lrdepth = MAX(branchruledata->lrdepth, depth);

            *success = TRUE;
         }

         branchruledata->lrk += nrepresentatives;
      }
      else
      {
         if( nodeWC > 0)
         {
            SCIP_CALL( genWC(scip, branchruledata, varsWC, valsWC, boundsWC, conssWC, nvarsWC, nconssWC) );

            branchruledata->wcsuccess++;
            nrepresentatives = nconssWC+2;
            branchruledata->wck += nrepresentatives;
            branchruledata->wcdepth = nvarsWC + 1 >= 0.75*SCIPgetNBinVars(scip) ? 1 : nvarsWC;
            *success = TRUE;
         }
      }
   }

   if( *success )
   {
      printf("** reoptimization ** heuristic compression of the search frontier:\n");
      printf("*               nnodes       loss loss/nodes      depth\n");

      printf("* nodes    %10d %10s %10s %10s\n", nnodesToCompress, "", "", "");

      if( successLR1 )
         printf("*    LR1 : %10d %10.2f %10.2f %10d\n", nnodesLR1, lossLR1, lossLR1/nnodesLR1, nvarsLR1);
      else
         printf("*    LR1 : %10d %10s %10s %10s\n", 0, "-", "-", "-");

      if( successLR2 )
         printf("*    LR2 : %10d %10.2f %10.2f %10d\n", nnodesLR2, lossLR2, lossLR2/nnodesLR2, nvarsLR2);
      else
         printf("*    LR2 : %10d %10s %10s %10s\n", 0, "-", "-", "-");

      if( nodeWC > 0 )
         printf("*     WC : %10d %10s %10s %10d (nvars)\n", 1, "-", "-", nvarsWC);
      else
         printf("*     WC : %10d %10s %10s %10s\n", 0, "-", "-", "-");

      printf("* \n");
      printf("* compression of size: %5d\n", nrepresentatives);
      printf("************************************************************\n");
   }

   SKIP:

   /* allocate memory for WC */
   SCIPfreeMemoryArray(scip, &varsWC);
   SCIPfreeMemoryArray(scip, &valsWC);
   SCIPfreeMemoryArray(scip, &boundsWC);
   SCIPfreeMemoryArray(scip, &conssWC);

   /* free memory for LR */
   SCIPfreeMemoryArray(scip, &LR2);
   SCIPfreeMemoryArray(scip, &LR1);
   SCIPfreeMemoryArray(scip, &valsLR2);
   SCIPfreeMemoryArray(scip, &varsLR2);
   SCIPfreeMemoryArray(scip, &boundsLR2);
   SCIPfreeMemoryArray(scip, &boundsLR1);
   SCIPfreeMemoryArray(scip, &valsLR1);
   SCIPfreeMemoryArray(scip, &varsLR1);

   /* free general memory */
   SCIPfreeMemoryArray(scip, &SetToCompress);

   return SCIP_OKAY;
}

static
SCIP_Real getCutoffbound(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID
)
{
   SCIP_Real cutoffbound;
   SCIP_VAR** vars;
   int sol;
   int var;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(0 <= nodeID && nodeID < branchruledata->allocmemsizenodedata);
   assert(branchruledata->nodedata[nodeID] != NULL);

   cutoffbound = SCIPinfinity(scip);
   vars = SCIPgetOrigVars(scip);

   if( branchruledata->nodedata[nodeID]->soldata == NULL || branchruledata->nodedata[nodeID]->soldata->nsols == 0 )
      return -cutoffbound;

   for(sol = 0; sol < branchruledata->nodedata[nodeID]->soldata->nsols; sol++)
   {
      SCIP_Real objval;
      objval = 0;

      for(var = 0; var < SCIPgetNOrigVars(scip); var++)
      {
         int idx;
         idx = SCIPvarGetIndex(vars[var]);

         objval += (SCIPvarGetObj(vars[var]) * branchruledata->nodedata[nodeID]->soldata->vars[sol][idx]);
      }

      cutoffbound = MIN(cutoffbound, objval);
   }

   return cutoffbound;
}

/*
 * check if the node is infeasible or redundant due to strong branching
 */
static
SCIP_RETCODE dryBranch(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   SCIP_Bool*            redundant,
   SCIP_Bool*            cutoff
)
{
   int var;
   int nredundantvars;
   int* redundantvars;
   SCIP_Real ub;
   SCIP_Real lb;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(1 <= nodeID && nodeID < branchruledata->allocmemsizenodedata);
   assert(branchruledata->nodedata[nodeID] != NULL);

   *cutoff = FALSE;
   *redundant = FALSE;
   nredundantvars = 0;

   SCIP_CALL( SCIPallocMemoryArray(scip, &redundantvars, branchruledata->nodedata[nodeID]->nvars) );

   SCIPdebugMessage("start dry branching at ID %d:\n", nodeID);

   if( branchruledata->nodedata[nodeID]->nvars > 0 )
   {
      for(var = 0; var < branchruledata->nodedata[nodeID]->nvars && !(*cutoff); var++)
      {
         assert(SCIPvarIsOriginal(branchruledata->nodedata[nodeID]->vars[var]));
         SCIP_CALL( SCIPvarGetProbvarBound(&branchruledata->nodedata[nodeID]->vars[var],
               &branchruledata->nodedata[nodeID]->varbounds[var],
               &branchruledata->nodedata[nodeID]->varboundtypes[var]) );
         assert(SCIPvarIsTransformed(branchruledata->nodedata[nodeID]->vars[var]));

         lb = SCIPvarGetLbLocal(branchruledata->nodedata[nodeID]->vars[var]);
         ub = SCIPvarGetUbLocal(branchruledata->nodedata[nodeID]->vars[var]);

         /* check for infeasibility */
         if( SCIPisFeasEQ(scip, lb, ub) && !SCIPisFeasEQ(scip, lb, branchruledata->nodedata[nodeID]->varbounds[var]) )
         {
            SCIPdebugMessage(" -> <%s> is fixed to %g, can not change bound to %g -> cutoff\n",
               SCIPvarGetName(branchruledata->nodedata[nodeID]->vars[var]), lb, branchruledata->nodedata[nodeID]->varbounds[var]);

            *cutoff = TRUE;
         }

         /* check for redundancy */
         if( SCIPisFeasEQ(scip, lb, ub) && SCIPisFeasEQ(scip, lb, branchruledata->nodedata[nodeID]->varbounds[var]) )
         {
            SCIPdebugMessage(" -> <%s> is already fixed to %g -> redundant bound change\n",
               SCIPvarGetName(branchruledata->nodedata[nodeID]->vars[var]), lb);

            redundantvars[nredundantvars] = var;
            nredundantvars++;
         }
      }

      /* delete redundant variables */
      if( nredundantvars == branchruledata->nodedata[nodeID]->nvars
       && (branchruledata->nodedata[nodeID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss)) )
      {
         *redundant = TRUE;
         SCIPdebugMessage(" -> node at ID %d is redundant\n", nodeID);
      }
      else if( nredundantvars > 0 )
      {
         /* all bound changes are redundant but there are added constraints */
         SCIPdebugMessage(" -> node at ID %d has %d/%d redundant bound changes, %d local constraints -> reduce to %d\n",
               nodeID, nredundantvars, branchruledata->nodedata[nodeID]->nvars,
               branchruledata->nodedata[nodeID]->nodechilds == NULL ? 0 : SCIPqueueNElems(branchruledata->nodedata[nodeID]->nodechilds),
               branchruledata->nodedata[nodeID]->nvars-nredundantvars);
         nredundantvars--;
         while( nredundantvars >= 0 )
         {
            var = redundantvars[nredundantvars];
            if( var < branchruledata->nodedata[nodeID]->nvars-1 )
            {
               branchruledata->nodedata[nodeID]->vars[var] = branchruledata->nodedata[nodeID]->vars[branchruledata->nodedata[nodeID]->nvars-1];
               branchruledata->nodedata[nodeID]->varbounds[var] = branchruledata->nodedata[nodeID]->varbounds[branchruledata->nodedata[nodeID]->nvars-1];
               branchruledata->nodedata[nodeID]->varboundtypes[var] = branchruledata->nodedata[nodeID]->varboundtypes[branchruledata->nodedata[nodeID]->nvars-1];
            }
            else
            {
               branchruledata->nodedata[nodeID]->vars[var] = NULL;
               branchruledata->nodedata[nodeID]->varbounds[var] = 0;
               branchruledata->nodedata[nodeID]->varboundtypes[var] = 0;
            }
            branchruledata->nodedata[nodeID]->nvars--;
            nredundantvars--;
         }
      }

      /* retransform into original space */
      SCIP_CALL( transformIntoOrig(scip, branchruledata, nodeID) );
   }
   else if( branchruledata->nodedata[nodeID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss) )
   {
      *redundant = TRUE;
      SCIPdebugMessage(" -> node at ID %d is redundant\n", nodeID);
   }

   /* free memory */
   SCIPfreeMemoryArray(scip, &redundantvars);

   return SCIP_OKAY;
}


/*
 * Execute the branching of nodes with additional constraints.
 */
static
SCIP_RETCODE Exec(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_RESULT*          result
)
{
   SCIP_Bool redundant;
   SCIP_Bool cutoff;
   int nodeID;
   int childID;
   int curChild;
   int depth;
   int nChilds;
   int ncreatedchilds;
   int nredundant;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   nodeID = SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) ? 0 : SCIPnodeGetReoptID(SCIPgetCurrentNode(scip));

   SCIPdebugMessage("branch reoptimized node %llu, ID %d:\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), nodeID);
   if( branchruledata->nodedata[nodeID]->dualcons[0] != NULL )
   {
      SCIPdebugMessage(" -> dual information: %u (constype: %d, nvars: %d)\n",
            branchruledata->nodedata[nodeID]->pseudobranched,
            branchruledata->nodedata[nodeID]->dualcons[0]->constype,
            branchruledata->nodedata[nodeID]->dualcons[0]->nvars);
   }
   else
   {
      SCIPdebugMessage(" -> dual information: %u (constype: -, nvars: -)\n", branchruledata->nodedata[nodeID]->pseudobranched);
   }

   SCIPdebugMessage(" -> nchildren: %d\n", branchruledata->nodedata[nodeID]->nodechilds == NULL ? 0 : SCIPqueueNElems(branchruledata->nodedata[nodeID]->nodechilds));

   /* calculate local similarity */
   if( branchruledata->localdelay > -1 )
   {
      /* reset the counter at the beginning of the bandb procedure */
      if( nodeID == 0 )
         branchruledata->ndelsubtrees = 0;

      SCIP_Real localsim;
      SCIP_Real localscalar;
      SCIP_Real localoldnorm;
      SCIP_Real localnewnorm;
      SCIP_Real locallb;
      SCIP_Real localub;
      SCIP_Real oldcoef;
      SCIP_Real newcoef;
      int var;
      int idx;

      localsim = 0.0;
      localscalar = 0.0;
      localoldnorm = 0.0;
      localnewnorm = 0.0;

      for(var = 0; var < SCIPgetNOrigBinVars(scip); var++)
      {
         locallb = SCIPvarGetLbLocal(SCIPgetOrigVars(scip)[var]);
         localub = SCIPvarGetUbLocal(SCIPgetOrigVars(scip)[var]);

         if( SCIPisFeasLT(scip, locallb, localub) )
         {
            idx = SCIPvarGetIndex(SCIPgetOrigVars(scip)[var]);
            oldcoef = SCIPreoptGetObjCoef(scip->reopt, scip->stat->reopt_nruns-1, idx);
            newcoef = SCIPreoptGetObjCoef(scip->reopt, scip->stat->reopt_nruns, idx);

            localscalar += (oldcoef * newcoef);
            localoldnorm += pow(oldcoef, 2);
            localnewnorm += pow(newcoef, 2);
         }
      }

      if(localnewnorm == 0 || localoldnorm == 0 || localscalar == 0)
         localsim = 0.0;
      else
         localsim = localscalar/(sqrt(localoldnorm)*sqrt(localnewnorm));

      if( SCIPisLT(scip, localsim, branchruledata->localdelay) )
      {
         SCIP_CALL( deleteChildrenBelow(scip, branchruledata, nodeID, FALSE) );

         if( branchruledata->nodedata[nodeID]->pseudobranched )
         {
            if( branchruledata->nodedata[nodeID]->dualcons[0] != NULL )
            {
               SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->dualcons[0]->vars);
               SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->dualcons[0]->vals);
               SCIPfreeMemory(scip, &branchruledata->nodedata[nodeID]->dualcons[0]);
               branchruledata->nodedata[nodeID]->dualcons[0] = NULL;
            }

            if( branchruledata->nodedata[nodeID]->dualcons[1] != NULL )
            {
               SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->dualcons[1]->vars);
               SCIPfreeMemoryArray(scip, &branchruledata->nodedata[nodeID]->dualcons[1]->vals);
               SCIPfreeMemory(scip, &branchruledata->nodedata[nodeID]->dualcons[1]);
               branchruledata->nodedata[nodeID]->dualcons[1] = NULL;
            }

            branchruledata->nodedata[nodeID]->pseudobranched = FALSE;
            branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_LEAF;
         }

         branchruledata->ndelsubtrees++;

         SCIPdebugMessage(" -> local similarity: %.4f (solve (sub)problem from scratch)\n", localsim);
         printf(" -> local similarity: %.4f (solve (sub)problem from scratch)\n", localsim);

         *result = SCIP_DIDNOTRUN;
         goto TERMINATE;
      }
      else
         SCIPdebugMessage(" -> local similarity: %.4f\n", localsim);
   }

   /**
    * current node is equal to the root and the root was pseudo-branched
    * we have to create two child nodes; one with the pseudo-constraint and
    * one with the negated fixings.
    */
   if(SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) && branchruledata->nodedata[0]->dualcons[0] != NULL)
   {
      LOGICORDATA* consdata;
      int child1_ID;
      int child2_ID;
      int v;

      assert(branchruledata->nodedata[0]->dualcons[0] != NULL);
      assert(branchruledata->nodedata[0]->dualcons[0]->vars != NULL);
      assert(branchruledata->nodedata[0]->dualcons[0]->vals != NULL);
      assert(branchruledata->nodedata[0]->dualcons[0]->nvars > 0);

      /***********************************************************************
       *
       * create two child nodes of the root node but do not create them,
       * add them to the data structure only.
       *
       * A: the added constraint is a strong branching constraint
       *   1. get a free slot to save all information, initialize the data
       *   2. copy all decisions from strong branching
       *   3. move all children of root node to child1 (if some exists)
       *
       * B: the added constraints is an infeasible subtree constraint
       *   DO NOTHING
       ***********************************************************************/

      switch (branchruledata->nodedata[0]->dualcons[0]->constype) {
         case REOPT_CONSTYPE_STRBRANCHED:

            /** add this node (the branching path is empty, thats way this method will on
             * allocate memory and create a link between node number and ID) */
            if (SCIPqueueIsEmpty(branchruledata->openIDs))
            {
               SCIP_CALL(reallocNodedata(scip, branchruledata));
            }

            child1_ID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
            assert(branchruledata->nodedata[child1_ID] == NULL
               || (branchruledata->nodedata[child1_ID]->nvars == 0
                   && (branchruledata->nodedata[child1_ID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[child1_ID]->conss))) );
            SCIP_CALL( initNode(scip, branchruledata, child1_ID) );
            branchruledata->nsavednodes++;

            SCIPdebugMessage("add a child at ID %d\n", child1_ID);

            /* tag node as TRANSIT */
            branchruledata->nodedata[child1_ID]->reopttype = SCIP_REOPTTYPE_TRANSIT;
            SCIPdebugMessage(" -> set reopttype: %d\n", SCIP_REOPTTYPE_TRANSIT);

            /** move all child nodes (w/o the current) of the pseudo-branched root node to child1 with child1_ID */
            if( branchruledata->nodedata[0]->nodechilds != NULL )
            {
               assert(branchruledata->nodedata[child1_ID]->nodechilds == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[child1_ID]->nodechilds));

               if(branchruledata->nodedata[child1_ID]->nodechilds == NULL )
               {
                  SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[child1_ID]->nodechilds, SCIPqueueNElems(branchruledata->nodedata[0]->nodechilds), 2) );
               }

               while(!SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds))
               {
                  SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[child1_ID]->nodechilds, SCIPqueueRemove(branchruledata->nodedata[0]->nodechilds)) );
               }
            }
            assert(branchruledata->nodedata[0]->nodechilds == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds));

            /* set parentID */
            branchruledata->nodedata[child1_ID]->parentID = 0;

            /* insert child1_ID as a root child, we have to do this because this node will be deleted after branching
             * by RemoveNode() */
            if( branchruledata->nodedata[0]->nodechilds == NULL )
            {
               SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[0]->nodechilds, 2, 2) );
            }
            SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[0]->nodechilds, (void* ) (size_t ) child1_ID) );

            /** fix the corresponding variables in child1 */
            if( branchruledata->nodedata[child1_ID]->allocmem == 0 )
            {
               SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[child1_ID]->vars, branchruledata->nodedata[0]->dualcons[0]->nvars) );
               SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[child1_ID]->varbounds, branchruledata->nodedata[0]->dualcons[0]->nvars) );
               SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[child1_ID]->varboundtypes, branchruledata->nodedata[0]->dualcons[0]->nvars) );
               branchruledata->nodedata[child1_ID]->allocmem = branchruledata->nodedata[0]->dualcons[0]->nvars;
            }
            else
            {
               SCIP_CALL( checkMemory(scip, branchruledata, child1_ID, branchruledata->nodedata[0]->dualcons[0]->nvars) );
            }
            assert(branchruledata->nodedata[child1_ID]->nvars == 0);

            for(v = 0; v < branchruledata->nodedata[0]->dualcons[0]->nvars; v++)
            {
               branchruledata->nodedata[child1_ID]->vars[v] = branchruledata->nodedata[0]->dualcons[0]->vars[v];
               branchruledata->nodedata[child1_ID]->varbounds[v] = branchruledata->nodedata[0]->dualcons[0]->vals[v];
               branchruledata->nodedata[child1_ID]->varboundtypes[v] = (SCIPisFeasEQ(scip, branchruledata->nodedata[0]->dualcons[0]->vals[v], 0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
               branchruledata->nodedata[child1_ID]->nvars++;
            }

         break;

         case REOPT_CONSTYPE_INFSUBTREE:
            break;

         default:
            break;
      }

      /***********************************************************************
       *
       * create copy of root and add all decisions from strong branching
       *
       * 1. add and release the logic-or constraint
       * 2. if the added constraint includes more than one variable, save the
       *    node
       *
       ***********************************************************************/

      /* add the logic-or constraint 'cons' to child2 */
      if (SCIPqueueIsEmpty(branchruledata->openIDs))
      {
         SCIP_CALL(reallocNodedata(scip, branchruledata));
      }

      child2_ID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
      assert(branchruledata->nodedata[child2_ID] == NULL
         || (branchruledata->nodedata[child2_ID]->nvars == 0
             && (branchruledata->nodedata[child2_ID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[child2_ID]->conss))) );
      SCIP_CALL( initNode(scip, branchruledata, child2_ID) );
      branchruledata->nsavednodes++;

      SCIPdebugMessage("add a child at ID %d\n", child2_ID);

      if( branchruledata->nodedata[0]->dualcons[0]->nvars > 1 )
      {
         /* tag node as LOGICOR */
         branchruledata->nodedata[child2_ID]->reopttype = SCIP_REOPTTYPE_LOGICORNODE;
         SCIPdebugMessage(" -> set reopttype: %d\n", SCIP_REOPTTYPE_LOGICORNODE);

         if( branchruledata->nodedata[child2_ID]->conss == NULL )
         {
            SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[child2_ID]->conss, 1, 2) );
         }
         assert(SCIPqueueIsEmpty(branchruledata->nodedata[child2_ID]->conss));

         SCIP_CALL( SCIPallocMemory(scip, &consdata) );
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vars, branchruledata->nodedata[0]->dualcons[0]->vars, branchruledata->nodedata[0]->dualcons[0]->nvars) );
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vals, branchruledata->nodedata[0]->dualcons[0]->vals, branchruledata->nodedata[0]->dualcons[0]->nvars) );
         consdata->nvars = branchruledata->nodedata[0]->dualcons[0]->nvars;
         consdata->constype = branchruledata->nodedata[0]->dualcons[0]->constype;

         SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[child2_ID]->conss, (void*) consdata) );
      }
      else
      {
         /* tag node as LEAF */
         branchruledata->nodedata[child2_ID]->reopttype = SCIP_REOPTTYPE_LEAF;
         SCIPdebugMessage(" -> set reopttype: %d\n", SCIP_REOPTTYPE_LEAF);

         if( branchruledata->nodedata[child2_ID]->allocmem == 0 )
         {
            branchruledata->nodedata[child2_ID]->allocmem = 1;
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &branchruledata->nodedata[child2_ID]->vars, branchruledata->nodedata[0]->dualcons[0]->vars, 1) );
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &branchruledata->nodedata[child2_ID]->varbounds, branchruledata->nodedata[0]->dualcons[0]->vals, 1) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[child2_ID]->varboundtypes, 1) );
         }
         else
         {
            assert(branchruledata->nodedata[child2_ID]->nvars == 0);
            SCIP_CALL( checkMemory(scip, branchruledata, child2_ID, 1) );
            branchruledata->nodedata[child2_ID]->vars[0] = branchruledata->nodedata[0]->dualcons[0]->vars[0];
            branchruledata->nodedata[child2_ID]->varbounds[0] = branchruledata->nodedata[0]->dualcons[0]->vals[0];
         }
         branchruledata->nodedata[child2_ID]->varboundtypes[0] = (SCIPisFeasEQ(scip, branchruledata->nodedata[0]->dualcons[0]->vals[0], 0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
         branchruledata->nodedata[child2_ID]->nvars = 1;
      }

      /* set the root as parent */
      branchruledata->nodedata[child2_ID]->parentID = 0;

      /* ad the ID as an child node */
      if( branchruledata->nodedata[0]->nodechilds == NULL )
      {
         SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[0]->nodechilds, 2, 2) );
      }
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[0]->nodechilds, (void* ) (size_t ) child2_ID) );

      /** free consdata and move */
      SCIPfreeMemoryArray(scip, &branchruledata->nodedata[0]->dualcons[0]->vars);
      SCIPfreeMemoryArray(scip, &branchruledata->nodedata[0]->dualcons[0]->vals);
      SCIPfreeMemory(scip, &branchruledata->nodedata[0]->dualcons[0]);
      branchruledata->nodedata[0]->dualcons[0] = NULL;

      if( branchruledata->nodedata[0]->dualcons[1] != NULL )
      {
         branchruledata->nodedata[0]->dualcons[0] = branchruledata->nodedata[0]->dualcons[1];
         branchruledata->nodedata[0]->dualcons[1] = NULL;
      }

      /** remove pseudobranched flag */
      branchruledata->nodedata[0]->pseudobranched = (branchruledata->nodedata[0]->dualcons[0] != NULL);
      branchruledata->nsavednodes--;

      nodeID = 0;
      goto REVIVE;
   }
   else if(SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip)
         && !branchruledata->nodedata[0]->pseudobranched)
   {
      assert(branchruledata->nodedata[0]->dualcons[0] == NULL);
      nodeID = 0;
      goto REVIVE;
   }

   /**
    * we have to generate all nodes in currentleafIDs and one node extra with the fixed bounds if the node was pseudo-branched.
    * all nodes with a non empty queue in 'nodechilds' will be add to currentleafIDs.
    */
   assert(nodeID >= 1);

   REVIVE:

   assert(branchruledata->nodedata[nodeID]->nodechilds != NULL );

   /** iterate over all children */
   nChilds = SCIPqueueNElems(branchruledata->nodedata[nodeID]->nodechilds);
   assert(nChilds > 0);

   ncreatedchilds = 0;
   nredundant = 0;
   curChild = 1;
   depth = lengthBranchPathByID(branchruledata, nodeID);

   cutoff = FALSE;
   redundant = FALSE;

   curChild = 1;

   /* check if the child node is infeasible or redundant due to strong branching.
    * infeasible nodes can be deleted. if a node is redundant, we move all children
    * to the current node, i.e., we shrink the node w/o copying the bound changes
    * because these are redundant.
    * */
   while( curChild <= nChilds && branchruledata->strongbranchinginit )
   {
      childID = (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[nodeID]->nodechilds);

      SCIP_CALL( dryBranch(scip, branchruledata, childID, &redundant, &cutoff) );

      if( redundant )
      {
         /* shrink the node */
         SCIPdebugMessage(" -> shrink redundant node: nchild = %d\n",
               branchruledata->nodedata[childID]->nodechilds == NULL ? 0 : SCIPqueueNElems(branchruledata->nodedata[childID]->nodechilds));

         if( branchruledata->nodedata[childID]->nodechilds != NULL )
         {
            SCIP_CALL( moveChildrenUp(scip, branchruledata, SCIPgetCurrentNode(scip), childID, nodeID) );
         }

         /* delete redundant node */
         SCIP_CALL( deleteNodeData(scip, branchruledata, childID, TRUE) );
         SCIP_CALL( SCIPqueueInsert(branchruledata->openIDs, (void*) (size_t) childID) );

         nredundant++;
      }
      else if( cutoff )
      {
         /* delete the induced subtree and the node itself */
         SCIP_CALL( deleteChildrenBelow(scip, branchruledata, childID, TRUE) );
      }
      else
      {
         /* the node is still a child node and need to be reoptimized */
         SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->nodechilds, (void*) (size_t) childID) );
      }

      curChild++;
   }

   nChilds = SCIPqueueNElems(branchruledata->nodedata[nodeID]->nodechilds);
   curChild = 1;

   while (curChild <= nChilds)
   {
      SCIP_NODE* child1;
      SCIP_NODE* child2;
      SCIP_Bool savedconsdata;

      savedconsdata = FALSE;

      child1 = NULL;
      child2 = NULL;

      childID = (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[nodeID]->nodechilds);
      /** childID has to be greater or equal than 1, because the ID 0 is
          reserved for the root and the root could never be a child */
      assert(childID >= 1);
      assert(childID < branchruledata->allocmemsizenodedata);

      SCIPdebugMessage("reoptimize node from ID %d\n", childID);

      ncreatedchilds++;

      /** the node at position childID was strongbranched, we have to split this node */
      if (branchruledata->nodedata[childID]->dualcons[0] != NULL )
      {
         assert(branchruledata->nodedata[childID]->reopttype == SCIP_REOPTTYPE_STRBRANCHED
             || branchruledata->nodedata[childID]->reopttype == SCIP_REOPTTYPE_INFSUBTREE);

         assert(branchruledata->nodedata[childID]->dualcons[0] != NULL);
         assert(branchruledata->nodedata[childID]->dualcons[0]->vars != NULL);
         assert(branchruledata->nodedata[childID]->dualcons[0]->vals != NULL);
         assert(branchruledata->nodedata[childID]->dualcons[0]->nvars > 0);

         SCIPdebugMessage(" -> dual information: %d (constype: %d, nvars %d)\n",
               branchruledata->nodedata[childID]->pseudobranched,
               branchruledata->nodedata[childID]->dualcons[0]->constype,
               branchruledata->nodedata[childID]->dualcons[0]->nvars);

         /* either the generated constraint split up the node into two new nodes
          * or the constraints separates a already found solution or cutoff an
          * infeasible subtree, respectively. */
         switch (branchruledata->nodedata[childID]->dualcons[0]->constype) {
            case REOPT_CONSTYPE_STRBRANCHED:
               assert(branchruledata->nodedata[childID]->reopttype == SCIP_REOPTTYPE_STRBRANCHED);

               /** the constraint split the node into two new nodes */
               SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
               SCIP_CALL( SCIPcreateChild(scip, &child2, 0.0, SCIPgetLocalTransEstimate(scip)) );
               assert(child1 != NULL );
               assert(child2 != NULL );

               SCIPdebugMessage(" -> reoptimize node from ID, depth %d\n", depth + branchruledata->nodedata[childID]->nvars + branchruledata->nodedata[childID]->dualcons[0]->nvars);

               /** change the bound in child1 and child2 along the branching path at position childID*/
               SCIP_CALL( changeAncestorBranchings(scip, branchruledata, childID, child1, child2) );

               /** fix all decisions from strong branching in child1 iff the constraint is of type STRBRNACHED */
               SCIP_CALL( fixVarsInNode(scip, branchruledata, child1, branchruledata->nodedata[childID]->dualcons[0], TRUE, FALSE, childID, FALSE) );

               SCIPdebugMessage(" -> new node from ID, depth %d\n", depth + branchruledata->nodedata[childID]->nvars);

               /** add 'pseudocons' to child2; the added constraint can not be redundant */
               SCIP_CALL( addLocalConsToNode(scip, branchruledata, child2, branchruledata->nodedata[childID]->dualcons[0], FALSE, -1) );

               break;

            case REOPT_CONSTYPE_INFSUBTREE:
               assert(branchruledata->nodedata[childID]->reopttype == SCIP_REOPTTYPE_INFSUBTREE);

               SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
               assert(child1 != NULL );

               SCIPdebugMessage(" -> new node %lld, depth %d\n", SCIPnodeGetNumber(child1), depth + branchruledata->nodedata[childID]->nvars);

               /** change the bound in child1 along the branching path at position childID*/
               SCIP_CALL( changeAncestorBranchings(scip, branchruledata, childID, child1, NULL) );

               /** add the constraint to child1; the added constraint can not be redundant */
               SCIP_CALL( addLocalConsToNode(scip, branchruledata, child1, branchruledata->nodedata[childID]->dualcons[0], FALSE, childID) );

               /** add the constraint to the nodedata at childID if the number of variables
                *  is greater than 1, else the variable was already fixed by addLocalConsToNode */
               if( branchruledata->nodedata[childID]->dualcons[0]->nvars > 1 )
               {
                  LOGICORDATA* consdata;

                  if( branchruledata->nodedata[childID]->conss == NULL )
                  {
                     SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[childID]->conss, 1, 2) );
                  }

                  SCIP_CALL( SCIPallocMemory(scip, &consdata) );
                  SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vars, branchruledata->nodedata[childID]->dualcons[0]->vars, branchruledata->nodedata[childID]->dualcons[0]->nvars) );
                  SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->vals, branchruledata->nodedata[childID]->dualcons[0]->vals, branchruledata->nodedata[childID]->dualcons[0]->nvars) );
                  consdata->nvars = branchruledata->nodedata[childID]->dualcons[0]->nvars;
                  consdata->constype = branchruledata->nodedata[childID]->dualcons[0]->constype;

                  SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[childID]->conss, (void*) consdata) );
               }

               if( branchruledata->nodedata[childID]->dualcons[0]->constype == REOPT_CONSTYPE_INFSUBTREE )
                  branchruledata->ninfsubtrees++;

               break;

            default:
               assert(branchruledata->nodedata[childID]->dualcons[0]->constype == REOPT_CONSTYPE_INFSUBTREE
                   || branchruledata->nodedata[childID]->dualcons[0]->constype == REOPT_CONSTYPE_STRBRANCHED);
               break;
         }

         /** free consdata and move */
         SCIPfreeMemoryArray(scip, &branchruledata->nodedata[childID]->dualcons[0]->vars);
         SCIPfreeMemoryArray(scip, &branchruledata->nodedata[childID]->dualcons[0]->vals);
         SCIPfreeMemory(scip, &branchruledata->nodedata[childID]->dualcons[0]);
         branchruledata->nodedata[childID]->dualcons[0] = NULL;

         if( branchruledata->nodedata[childID]->dualcons[1] != NULL )
         {
            branchruledata->nodedata[childID]->dualcons[0] = branchruledata->nodedata[childID]->dualcons[1];
            branchruledata->nodedata[childID]->dualcons[1] = NULL;
         }

      }
      else
      {
         /**
          * node at position childID was not pseudobranched, we generate only the node itself
          */
         assert(branchruledata->nodedata[childID]->dualcons[0] == NULL);

         SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
         assert(child1 != NULL );
         assert(child2 == NULL );

         SCIPdebugMessage(" -> dual information: 0\n");
         SCIPdebugMessage(" -> create new node %lld, depth %d\n", SCIPnodeGetNumber(child1), depth + branchruledata->nodedata[childID]->nvars);

         if( branchruledata->nodedata[childID]->vars != NULL )
         {
            /** Chance bounds due to all ancestors */
            SCIP_CALL( changeAncestorBranchings(scip, branchruledata, childID, child1, NULL) );
         }
         else
         {
            assert(branchruledata->nodedata[childID]->conss != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[childID]->conss));
         }
      }

      /** set LPI state */
      if(branchruledata->savelpbasis
      && branchruledata->nodedata[childID]->lpistate != NULL
      && branchruledata->nodedata[childID]->conss == NULL
      && SCIPreoptGetSimToPrevious(scip->reopt) > 0.9)
      {
         SCIP_CALL( SCIPchildSetLpistate(child1, branchruledata->nodedata[childID]->lpistate) );
         SCIPdebugMessage("use LPI from previous round in node %lld (sim = %.4f)\n",SCIPnodeGetNumber(child1), SCIPreoptGetSimToPrevious(scip->reopt));
         printf("use LPI from previous round in node %lld (sim = %.4f)\n",SCIPnodeGetNumber(child1), SCIPreoptGetSimToPrevious(scip->reopt));
      }

      /** add local constraint from an iteration before (if some exists) to child1 */
      if(branchruledata->nodedata[childID]->conss != NULL && SCIPqueueNElems(branchruledata->nodedata[childID]->conss) - savedconsdata > 0 )
      {
         /** generate all local constraints and add them to child1 and child2 (if exists) */
         SCIP_CALL(genLocalCons(scip, branchruledata, childID, child1, child2, savedconsdata));
      }

      /** remove flag 'pseudobranched' */
      branchruledata->nodedata[childID]->pseudobranched = (branchruledata->nodedata[childID]->dualcons[0] != NULL);

      /** insert child1 to the child list of nodeID */
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->nodechilds, (void* ) (size_t ) childID) );
      SCIPnodeSetReoptID(child1, childID);

      /* set the REOPTTYPE */
      assert(branchruledata->nodedata[childID]->reopttype >= SCIP_REOPTTYPE_TRANSIT);
      if( branchruledata->nodedata[childID]->reopttype == SCIP_REOPTTYPE_STRBRANCHED )
         branchruledata->nodedata[childID]->reopttype = SCIP_REOPTTYPE_TRANSIT;
      SCIPnodeSetReopttype(child1, branchruledata->nodedata[childID]->reopttype);
      SCIPdebugMessage(" -> set reopttype: %d\n", branchruledata->nodedata[childID]->reopttype);

      /** check if child2 includes some added constraints, save this node and set REOPT to true */
      if(child2 != NULL && SCIPnodeGetNAddedcons(scip, child2) > 0)
      {
         assert(branchruledata->nodedata[childID] != NULL);
         SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, child2, SCIP_REOPTTYPE_LEAF, TRUE) );
      }

      curChild++;
   }

   if( ncreatedchilds == 0 )
      *result = SCIP_DIDNOTRUN;
   else
      *result = SCIP_BRANCHED;

   TERMINATE:

   SCIPdebugMessage("*******************************************\n");

   return SCIP_OKAY;
}

/*
 * Check if a restart because of (in)feasible constraints should be done.
 * This depends on the parameter-setting of the user.
 */
static
void restartCheck(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata
)
{
   SCIP_Real sim;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   if( branchruledata->restart )
   {
      printf("** reoptimization ** forced restart.\n");
   }
   /* check if maximal number if saved nodes is reached */
   else if (branchruledata->nsavednodes > branchruledata->maxsavednodes)
   {
      branchruledata->restart = TRUE;
      printf("** reoptimization ** reset search tree: nodelimit exceeded %d / %d\n", branchruledata->nsavednodes, branchruledata->maxsavednodes);
   }
   /* the current and the previous objective are not similar enough */
   else if( SCIPcheckRestartReopt(scip, &sim) )
   {
      branchruledata->restart = TRUE;
      printf("** reoptimization ** reset search tree: new objective is not similar enough (%.4f).\n", sim);
   }

   if( branchruledata->restart )
   {
      if( branchruledata->firstrestart == -1 )
         branchruledata->firstrestart = branchruledata->nruns;

      branchruledata->lastrestart = branchruledata->nruns;
   }

   return;
}

/*
 * change the reopttype of the subtree induced by nodeID
 */
static
SCIP_RETCODE changeReopttypeOfSubtree(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   SCIP_REOPTTYPE        reopttype
)
{
   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(branchruledata->nodedata[nodeID] != NULL);

   if( branchruledata->nodedata[nodeID]->nodechilds != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds) )
   {
      int childID;
      int nchildIDs;
      int seenIDs;

      nchildIDs = SCIPqueueNElems(branchruledata->nodedata[nodeID]->nodechilds);
      seenIDs = 0;

      while( seenIDs < nchildIDs )
      {
         /* get childID */
         childID = (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[nodeID]->nodechilds);
         assert(branchruledata->nodedata[childID] != NULL);

         if( branchruledata->nodedata[childID]->reopttype != SCIP_REOPTTYPE_STRBRANCHED
          && branchruledata->nodedata[childID]->reopttype != SCIP_REOPTTYPE_INFSUBTREE
          && branchruledata->nodedata[childID]->reopttype != SCIP_REOPTTYPE_INFEASIBLE )
            branchruledata->nodedata[childID]->reopttype = reopttype;

         /* change reopttype of subtree */
         SCIP_CALL( changeReopttypeOfSubtree(scip, branchruledata, childID, reopttype) );

         /* insert the childID */
         SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->nodechilds, (void*) (size_t) childID) );
         seenIDs++;
      }
   }

   return SCIP_OKAY;
}

/*
 * update statictics
 */
static
void updateStatistics(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata
)
{
   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(branchruledata->reopt);

   if (branchruledata->nfeasnodesround > branchruledata->nfeasnodesmax)
      branchruledata->nfeasnodesmax = branchruledata->nfeasnodesround;
   branchruledata->nfeasnodes += branchruledata->nfeasnodesround;
   branchruledata->nfeasnodesround = 0;

   if (branchruledata->ninfeasnodesround > branchruledata->ninfeasnodesmax)
      branchruledata->ninfeasnodesmax = branchruledata->ninfeasnodesround;
   branchruledata->ninfeasnodes += branchruledata->ninfeasnodesround;
   branchruledata->ninfeasnodesround = 0;

   if (branchruledata->nprunednodesround > branchruledata->nprunednodesmax)
      branchruledata->nprunednodesmax = branchruledata->nprunednodesround;
   branchruledata->nprunednodes += branchruledata->nprunednodesround;
   branchruledata->nprunednodesround = 0;

   if (branchruledata->npbnodesround > branchruledata->npbnodesmax)
      branchruledata->npbnodesmax = branchruledata->npbnodesround;
   branchruledata->npbnodes += branchruledata->npbnodesround;
   branchruledata->npbnodesround = 0;

   if (branchruledata->nrediednodesround > branchruledata->nrediednodesmax)
      branchruledata->nrediednodesmax = branchruledata->nrediednodesround;
   branchruledata->nrediednodes += branchruledata->nrediednodesround;
   branchruledata->nrediednodesround = 0;

   branchruledata->nruns++;
   branchruledata->nrestartsround = 0;
   branchruledata->nrevivednodes = 0;

   return;
}


/***********************
 * non-static methods
 ***********************/

/*
 * checks if reoptimization detects infeasibility
 */
SCIP_RETCODE SCIPbranchnodereoptCheckFeasibility(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL );
   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   if( !branchruledata->init || !branchruledata->reopt)
      return SCIP_OKAY;

   if( (branchruledata->nodedata[0]->nodechilds == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds))
     && branchruledata->nodedata[0]->reopttype == SCIP_REOPTTYPE_TRANSIT)
   {
      branchruledata->infeasibleproblem = TRUE;
      SCIP_CALL( addInfeasibleConstraint(scip) );
      branchruledata->infeasibleconsadded = TRUE;
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE shrinkNode(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_NODE*            node,
   int                   nodeID,
   SCIP_Bool*            shrank
)
{
   SCIP_NODE* parent;
   int ndomchgs;
   int parentID;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(node != NULL);
   assert(0 < nodeID && nodeID < branchruledata->allocmemsizenodedata);
   assert(branchruledata->nodedata[nodeID] != NULL);

   if( branchruledata->nodedata[nodeID]->nodechilds != NULL
    && !SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds) )
   {
      ndomchgs = 0;
      parentID = -1;
      parent = NULL;

      SCIP_CALL( getLastSavedNode(scip, branchruledata, node, &parent, &parentID, &ndomchgs) );

      assert(parentID != nodeID);
      assert(branchruledata->nodedata[parentID] != NULL );
      assert(branchruledata->nodedata[parentID]->nodechilds != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[parentID]->nodechilds));

      /* check if we want move all children to the next saved node above */
      if( (!branchruledata->dynamicdiffofnodes && ndomchgs <= branchruledata->maxdiffofnodes)
          ||(branchruledata->dynamicdiffofnodes && ndomchgs <= ceil(log10((SCIP_Real) (SCIPgetNOrigBinVars(scip) - MIN(SCIPgetNOrigBinVars(scip)-1, lengthBranchPath(parent, branchruledata))))/log10(2.0))) )
      {
         SCIPdebugMessage("shrink node %lld at ID %d, nchild %d\n", SCIPnodeGetNumber(node), nodeID,
               branchruledata->nodedata[nodeID]->nodechilds == NULL ? 0 : SCIPqueueNElems(branchruledata->nodedata[nodeID]->nodechilds));

         /* copy the references of child nodes  to the parent*/
         SCIP_CALL( moveChildrenUp(scip, branchruledata, parent, nodeID, parentID) );

         /* delete the current node */
         while( (int) (size_t) SCIPqueueFirst(branchruledata->nodedata[parentID]->nodechilds) != nodeID )
         {
            SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[parentID]->nodechilds, SCIPqueueRemove(branchruledata->nodedata[parentID]->nodechilds)) );
         }
         assert((int) (size_t) SCIPqueueFirst(branchruledata->nodedata[parentID]->nodechilds) == nodeID);
         SCIPqueueRemove(branchruledata->nodedata[parentID]->nodechilds);

         SCIP_CALL( deleteNodeData(scip, branchruledata, nodeID, TRUE) );
         SCIP_CALL( SCIPqueueInsert(branchruledata->openIDs, (void*) (size_t) nodeID) );

         *shrank = TRUE;
      }
   }

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPbranchruleNodereoptSetRootLPI(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL );
   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(branchruledata->savelpbasis);
   assert(SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip));

   /** set the LPI state for the root node (if some exists) */
   if( branchruledata->nodedata[0] != NULL
    && branchruledata->nodedata[0]->lpistate != NULL
    && SCIPreoptGetSimToPrevious(scip->reopt) > 0.5)
   {
#ifdef DEBUG_MODE
      printf("** reoptimization ** use LPI from previous round in root node (sim = %.4f)\n", SCIPreoptGetSimToPrevious(scip->reopt));
#endif
      scip->tree->focuslpistate = branchruledata->nodedata[0]->lpistate;
   }

   return SCIP_OKAY;
}

/*
 * Add a pruned node the data structure.
 */
SCIP_RETCODE SCIPbranchruleNodereoptAddNode(
   SCIP*                 scip,
   SCIP_NODE*            node,                    /** current node */
   SCIP_REOPTTYPE        reopttype,
   SCIP_Bool             saveafterdual
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_NODE* parent;
   SCIP_Bool shrank;
   int nodeID;
   int parentID;
   int ndualvars;
   int pos;

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   parentID = -1;
   parent = NULL;
   shrank = FALSE;

   if( branchruledata->maxsavednodes == 0 )
      return SCIP_OKAY;

   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(reopttype == SCIP_REOPTTYPE_TRANSIT
       || reopttype == SCIP_REOPTTYPE_INFSUBTREE
       || reopttype == SCIP_REOPTTYPE_STRBRANCHED
       || reopttype == SCIP_REOPTTYPE_LOGICORNODE
       || reopttype == SCIP_REOPTTYPE_LEAF
       || reopttype == SCIP_REOPTTYPE_PRUNED
       || reopttype == SCIP_REOPTTYPE_FEASIBLE);

   /** start clock */
   SCIP_CALL (SCIPstartSaveTime(scip) );

   if (SCIPnodeGetReoptID(node) >= 1)
   {
      assert(reopttype != SCIP_REOPTTYPE_LEAF);

      nodeID = SCIPnodeGetReoptID(node);
      assert(nodeID < branchruledata->allocmemsizenodedata);
      assert(branchruledata->nodedata[nodeID] != NULL);

      SCIPdebugMessage("update node %lld at ID %u:\n", SCIPnodeGetNumber(node), nodeID);

      if( saveafterdual )
      {
         SCIPdebugMessage(" -> update after-dual bound changes.\n");
         SCIP_CALL( saveAfterDualBranchings(scip, branchruledata, node, nodeID) );

         /* ensure that all variables are original */
         SCIP_CALL( transformIntoOrig(scip, branchruledata, nodeID) );
      }
#ifdef SCIP_DEBUG
         int varnr;

         SCIPdebugMessage(" -> nvars: %d, ncons: %d, parentID: %d, reopttype: %d\n",
               branchruledata->nodedata[nodeID]->nvars,
               branchruledata->nodedata[nodeID]->conss != NULL ? SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss) : 0,
               branchruledata->nodedata[nodeID]->parentID, reopttype);
         SCIPdebugMessage(" -> saved variables:\n");

         for (varnr = 0; varnr < branchruledata->nodedata[nodeID]->nvars; varnr++)
            SCIPdebugMessage("  <%s> %s %f\n", SCIPvarGetName(branchruledata->nodedata[nodeID]->vars[varnr]),
                  branchruledata->nodedata[nodeID]->varboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                        "=>" : "<=", branchruledata->nodedata[nodeID]->varbounds[varnr]);

         for (varnr = 0; varnr < branchruledata->nodedata[nodeID]->nafterdualvars; varnr++)
            SCIPdebugMessage("  <%s> %s %f (after dual red.)\n", SCIPvarGetName(branchruledata->nodedata[nodeID]->afterdualvars[varnr]),
                  branchruledata->nodedata[nodeID]->afterdualvarboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                        "=>" : "<=", branchruledata->nodedata[nodeID]->afterdualvarbounds[varnr]);
#endif

      /** update LPI state if node is pseudobranched or feasible */
      switch (reopttype) {
         case SCIP_REOPTTYPE_TRANSIT:
            assert(branchruledata->nodedata[nodeID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss));

            if( branchruledata->shrink )
            {
               SCIP_CALL( shrinkNode(scip, branchruledata, node, nodeID, &shrank) );
            }

            goto TRANSIT;

            break;

         case SCIP_REOPTTYPE_LOGICORNODE:
         case SCIP_REOPTTYPE_LEAF:
            goto TRANSIT;
            break;

         case SCIP_REOPTTYPE_INFSUBTREE:
            /* delete the whole subtree induced be the current node */
            SCIP_CALL( deleteChildrenBelow(scip, branchruledata, nodeID, FALSE) );
            goto PSEUDO;
            break;

         case SCIP_REOPTTYPE_STRBRANCHED:
            /* dive through all children and change the reopttype to LEAF */
            SCIP_CALL( changeReopttypeOfSubtree(scip, branchruledata, nodeID, SCIP_REOPTTYPE_PRUNED) );
            goto PSEUDO;
            break;

         case SCIP_REOPTTYPE_FEASIBLE:
            /* delete the subtree */
            if( branchruledata->reducetofrontier )
            {
               SCIP_CALL( deleteChildrenBelow(scip, branchruledata, nodeID, FALSE) );
            }
            /* dive through all children and change the reopttype to LEAF */
            else
            {
               SCIP_CALL( changeReopttypeOfSubtree(scip, branchruledata, nodeID, SCIP_REOPTTYPE_PRUNED) );
            }
            goto FEASIBLE;
            break;

         case SCIP_REOPTTYPE_PRUNED:
            branchruledata->nrediednodesround++;
            /* delete the subtree */
            if( branchruledata->reducetofrontier )
            {
               SCIP_CALL( deleteChildrenBelow(scip, branchruledata, nodeID, FALSE) );
            }
            /* dive through all children and change the reopttype to LEAF */
            else
            {
               SCIP_CALL( changeReopttypeOfSubtree(scip, branchruledata, nodeID, SCIP_REOPTTYPE_PRUNED) );
            }
            goto PRUNED;

         default:
            break;
      }

      /** stop clock */
      SCIP_CALL(SCIPstopSaveTime(scip));

      return SCIP_OKAY;
   }

   /* get new IDs */
   if (SCIPqueueIsEmpty(branchruledata->openIDs))
   {
      SCIP_CALL(reallocNodedata(scip, branchruledata));
   }

   /** the current node is the root node and it can be only */
   if (SCIPnodeGetDepth(node) == 0)
   {
      nodeID = 0;

      switch (reopttype) {
         case SCIP_REOPTTYPE_TRANSIT:
            goto TRANSIT;
            break;

         case SCIP_REOPTTYPE_INFSUBTREE:
         case SCIP_REOPTTYPE_STRBRANCHED:
            branchruledata->nodedata[nodeID]->reopttype = reopttype;
            branchruledata->nodedata[nodeID]->pseudobranched = TRUE;
            branchruledata->nodedata[nodeID]->nvars = 0;

            if( reopttype == SCIP_REOPTTYPE_INFSUBTREE )
            {
               /* delete the whole subtree induced be the current node */
               SCIP_CALL( deleteChildrenBelow(scip, branchruledata, nodeID, FALSE) );
            }

            SCIPdebugMessage("update node %lld at ID %u:\n", SCIPnodeGetNumber(node), nodeID);
            SCIPdebugMessage(" -> nvars: 0, ncons: 0, parentID: -, reopttype: %d\n", reopttype);

            goto PSEUDO;
            break;

         case SCIP_REOPTTYPE_FEASIBLE:
            branchruledata->nfeasnodesround++;

            if( !branchruledata->nodedata[nodeID]->pseudobranched )
               branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_FEASIBLE;

            if( branchruledata->nodedata[0]->nodechilds != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds) )
            {
              /* delete the subtree */
               if( branchruledata->reducetofrontier )
               {
                  SCIP_CALL( deleteChildrenBelow(scip, branchruledata, nodeID, FALSE) );
                  SCIP_CALL( SCIPbranchrulePseudoReset(scip, TRUE, FALSE) );

               }
               /* dive through all children and change the reopttype to LEAF */
               else
               {
                  SCIP_CALL( changeReopttypeOfSubtree(scip, branchruledata, nodeID, SCIP_REOPTTYPE_PRUNED) );
               }
            }

            /* save the basis if the node */
            if( branchruledata->savelpbasis )
            {
               SCIP_CALL(saveLPIstate(scip, branchruledata, node, nodeID));
            }

            SCIPdebugMessage("update node %lld at ID %u:\n", SCIPnodeGetNumber(node), nodeID);
            SCIPdebugMessage(" -> nvars: 0, ncons: 0, parentID: -, reopttype: %d\n", reopttype);

            break;

         case SCIP_REOPTTYPE_PRUNED:
            branchruledata->nprunednodesround++;

            if( !branchruledata->nodedata[nodeID]->pseudobranched )
               branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_PRUNED;

            if( branchruledata->nodedata[0]->nodechilds != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds) )
            {
               /* delete the subtree */
               if( branchruledata->reducetofrontier )
               {
                  SCIP_CALL( deleteChildrenBelow(scip, branchruledata, nodeID, FALSE) );
                  SCIP_CALL( SCIPbranchrulePseudoReset(scip, TRUE, FALSE) );
               }
               /* dive through all children and change the reopttype to LEAF */
               else
               {
                  SCIP_CALL( changeReopttypeOfSubtree(scip, branchruledata, nodeID, SCIP_REOPTTYPE_PRUNED) );
               }
            }

            SCIPdebugMessage("update node %lld at ID %u:\n", SCIPnodeGetNumber(node), nodeID);
            SCIPdebugMessage(" -> nvars: 0, ncons: 0, parentID: -, reopttype: %d\n", reopttype);

            break;

         default:
            assert(reopttype == SCIP_REOPTTYPE_TRANSIT
                || reopttype == SCIP_REOPTTYPE_INFSUBTREE
                || reopttype == SCIP_REOPTTYPE_STRBRANCHED
                || reopttype == SCIP_REOPTTYPE_PRUNED
                || reopttype == SCIP_REOPTTYPE_FEASIBLE);
            break;
      }

      /** stop clock */
      SCIP_CALL(SCIPstopSaveTime(scip));

      return SCIP_OKAY;
   }
   else
   {
      int nbndchgdiff;

      /** check if we realy want to save this node:
       *  1. save the node if reopttype is at least LOGICORNODE
       *  2. save the node if the number of bound changes of this node
       *     and the last saved node is at least a given number n */

      /* get the ID of the last saved node or 0 for the root */
      SCIP_CALL( getLastSavedNode(scip, branchruledata, node, &parent, &parentID, &nbndchgdiff) );

      if( reopttype < SCIP_REOPTTYPE_INFSUBTREE
        && ((!branchruledata->dynamicdiffofnodes && nbndchgdiff <= branchruledata->maxdiffofnodes)
            ||(branchruledata->dynamicdiffofnodes && nbndchgdiff <= ceil(log10((SCIP_Real)(SCIPgetNOrigBinVars(scip) - MIN(SCIPgetNOrigBinVars(scip)-1,lengthBranchPath(parent, branchruledata))))/log10(2.0))) ) )
      {
         /** stop clock */
         SCIP_CALL(SCIPstopSaveTime(scip));

         return SCIP_OKAY;
      }


      /** get an empty slot (to ensure that ID != NULL, the queue should contains only IDs >= 1);
       * the 0 is always reserved for the root */
      if (SCIPqueueIsEmpty(branchruledata->openIDs))
      {
         SCIP_CALL(reallocNodedata(scip, branchruledata));
      }

      nodeID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
      assert(branchruledata->nodedata[nodeID] == NULL
         || (branchruledata->nodedata[nodeID]->nvars == 0
             && (branchruledata->nodedata[nodeID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss))) );
      assert(nodeID >= 1 && nodeID < branchruledata->allocmemsizenodedata);
      assert(SCIPgetRootNode(scip) != node);

      /** get memory for nodedata */
      assert(branchruledata->nodedata[nodeID] == NULL || branchruledata->nodedata[nodeID]->nvars == 0);
      SCIP_CALL(initNode(scip, branchruledata, nodeID));
      branchruledata->nodedata[nodeID]->parentID = parentID;

      assert(parent != NULL );
      assert((parent == SCIPgetRootNode(scip) && parentID == 0) || (parent != SCIPgetRootNode(scip) && parentID > 0));
      assert(nodeID >= 1);

      /** create the queue of "child nodes" if they not exist */
      if (branchruledata->nodedata[parentID]->nodechilds == NULL )
      {
         SCIP_CALL(SCIPqueueCreate(&branchruledata->nodedata[parentID]->nodechilds, 10, 2));
      }

      /** add the "child node" */
      assert(branchruledata->nodedata[parentID]->nodechilds != NULL );
      SCIP_CALL(SCIPqueueInsert(branchruledata->nodedata[parentID]->nodechilds, (void* ) (size_t ) nodeID));

      /* save branching path */
      SCIP_CALL( saveAncestorBranchings(scip, branchruledata, node, parent, nodeID, parentID) );

      /* save bound changes after some dual reduction */
      if( saveafterdual )
      {
         SCIP_CALL( saveAfterDualBranchings(scip, branchruledata, node, nodeID) );
      }
      else
      {
         SCIPdebugMessage(" -> skip saving bound changes after dual reductions.\n");
      }

      /** transform all bounds of branched variables and ensure that they are original. */
      SCIP_CALL( transformIntoOrig(scip, branchruledata, nodeID) );

      /** save pseudo-constraints (if one exists) */
      if (SCIPnodeGetNAddedcons(scip, node) >= 1)
      {
         assert(branchruledata->nodedata[nodeID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss));
         if (branchruledata->nodedata[nodeID]->conss == NULL )
         {
            SCIP_CALL(SCIPqueueCreate(&branchruledata->nodedata[nodeID]->conss, SCIPnodeGetNAddedcons(scip, node), 2));
         }
         SCIP_CALL(saveLocalConssData(scip, branchruledata, node, nodeID));
      }

      /* set ID */
      SCIPnodeSetReoptID(node, nodeID);

      /* set the REOPTTYPE */
      SCIPnodeSetReopttype(node, reopttype);

      branchruledata->nsavednodes++;

#ifdef SCIP_DEBUG
      int varnr;
      SCIPdebugMessage("save node #%lld\n", SCIPnodeGetNumber(node));
      SCIPdebugMessage(" -> ID %d, nvars %d, ncons %d, reopttype %d\n",
            nodeID, branchruledata->nodedata[nodeID]->nvars + branchruledata->nodedata[nodeID]->nafterdualvars,
            branchruledata->nodedata[nodeID]->conss == NULL ? 0 : SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss),
            reopttype);
      for (varnr = 0; varnr < branchruledata->nodedata[nodeID]->nvars; varnr++)
      {
         SCIPdebugMessage("  <%s> %s %f\n", SCIPvarGetName(branchruledata->nodedata[nodeID]->vars[varnr]),
               branchruledata->nodedata[nodeID]->varboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                     "=>" : "<=", branchruledata->nodedata[nodeID]->varbounds[varnr]);
      }
      for (varnr = 0; varnr < branchruledata->nodedata[nodeID]->nafterdualvars; varnr++)
      {
         SCIPdebugMessage("  <%s> %s %f (after dual red.)\n", SCIPvarGetName(branchruledata->nodedata[nodeID]->afterdualvars[varnr]),
               branchruledata->nodedata[nodeID]->afterdualvarboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                     "=>" : "<=", branchruledata->nodedata[nodeID]->afterdualvarbounds[varnr]);
      }
#endif
   }

   switch (reopttype) {
      case SCIP_REOPTTYPE_TRANSIT:
      case SCIP_REOPTTYPE_LOGICORNODE:
      case SCIP_REOPTTYPE_LEAF:
         TRANSIT:

         if( !shrank )
         {
            branchruledata->nodedata[nodeID]->reopttype = reopttype;

            if( branchruledata->savelpbasis
             && reopttype != SCIP_REOPTTYPE_LOGICORNODE
             && SCIPgetCurrentNode(scip) == node
             && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
            {
               SCIP_CALL( saveLPIstate(scip, branchruledata, node, nodeID) );
               branchruledata->nsavednodes++;
            }
         }
         else
         {
            SCIPnodeSetReoptID(node, -1);
            SCIPnodeSetReopttype(node, SCIP_REOPTTYPE_NONE);
         }
         break;

      case SCIP_REOPTTYPE_INFSUBTREE:
      case SCIP_REOPTTYPE_STRBRANCHED:
         PSEUDO:

         branchruledata->nodedata[nodeID]->reopttype = reopttype;
         branchruledata->nodedata[nodeID]->pseudobranched = TRUE;
         branchruledata->npbnodesround++;

         if (nodeID > 0)
         {
            assert(node != SCIPgetRootNode(scip));
            SCIP_CALL(SCIPbranchrulePseudoLinkIDs(scip, nodeID));
         }

         /* save the basis if the node */
         if( branchruledata->savelpbasis
          && branchruledata->savelpbasis
          && reopttype == SCIP_REOPTTYPE_STRBRANCHED
          && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
         {
            SCIP_CALL(saveLPIstate(scip, branchruledata, node, nodeID));
         }

         /* finish the collection of dual information */
         SCIP_CALL( SCIPbranchrulePseudoNodeFinished(scip, node, (reopttype == SCIP_REOPTTYPE_STRBRANCHED ? REOPT_CONSTYPE_STRBRANCHED : REOPT_CONSTYPE_INFSUBTREE)) );

         /* get the dual information */
         ndualvars = SCIPbranchrulePseudoGetNPseudoVars(scip, nodeID);
         pos = 0;

         if( branchruledata->nodedata[nodeID]->dualcons[0] != NULL )
         {
            assert(branchruledata->nodedata[nodeID]->dualcons[1] == NULL);
            pos = 1;
         }
         assert(branchruledata->nodedata[nodeID]->dualcons[pos] == NULL);

         SCIP_CALL( SCIPallocMemory(scip, &branchruledata->nodedata[nodeID]->dualcons[pos]) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeID]->dualcons[pos]->vars, ndualvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeID]->dualcons[pos]->vals, ndualvars) );
         branchruledata->nodedata[nodeID]->dualcons[pos]->nvars = -1;

         SCIP_CALL( SCIPbranchrulePseudoGenerateCons(scip, branchruledata->nodedata[nodeID]->dualcons[pos],
              &(branchruledata->nodedata[nodeID]->dualcons[pos]->nvars), ndualvars, nodeID, FALSE, TRUE) );
         assert(branchruledata->nodedata[nodeID]->dualcons[pos]->nvars > 0);
         assert(branchruledata->nodedata[nodeID]->dualcons[pos]->vars != NULL );
         assert(branchruledata->nodedata[nodeID]->dualcons[pos]->vals != NULL );
         assert(branchruledata->nodedata[nodeID]->dualcons[pos]->constype == REOPT_CONSTYPE_INFSUBTREE
             || branchruledata->nodedata[nodeID]->dualcons[pos]->constype == REOPT_CONSTYPE_STRBRANCHED);

         SCIPdebugMessage(" -> save dual information: node %lld, nvars %d, constype %d\n",
               SCIPnodeGetNumber(node), branchruledata->nodedata[nodeID]->dualcons[pos]->nvars,
               branchruledata->nodedata[nodeID]->dualcons[pos]->constype);

         break;

      case SCIP_REOPTTYPE_FEASIBLE:
         FEASIBLE:
         branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_FEASIBLE;
         branchruledata->nodedata[nodeID]->pseudobranched = FALSE;
         branchruledata->nfeasnodesround++;

         /**
          * save all information of the current feasible solution to separate this
          * solution in a following round (but only if all variablea are binary)
          * TODO: Verbesserungswuerdig
          */
         if( branchruledata->sepalocsols && nodeID > 0 )
         {
            SCIP_CALL( saveLocalConsData(scip, branchruledata, node, nodeID, REOPT_CONSTYPE_SEPASOLUTION) );
         }

         /* save the basis if the node */
         if( branchruledata->savelpbasis )
         {
            SCIP_CALL(saveLPIstate(scip, branchruledata, node, nodeID));
         }

         break;

      case SCIP_REOPTTYPE_PRUNED:
         PRUNED:

         branchruledata->lastprunedID = nodeID;

         branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_PRUNED;
         branchruledata->nodedata[nodeID]->pseudobranched = FALSE;
         branchruledata->nprunednodesround++;

         break;

      default:
         assert(reopttype == SCIP_REOPTTYPE_TRANSIT
             || reopttype == SCIP_REOPTTYPE_LOGICORNODE
             || reopttype == SCIP_REOPTTYPE_LEAF
             || reopttype == SCIP_REOPTTYPE_INFSUBTREE
             || reopttype == SCIP_REOPTTYPE_STRBRANCHED
             || reopttype == SCIP_REOPTTYPE_FEASIBLE
             || reopttype == SCIP_REOPTTYPE_PRUNED);
         break;
   }

   /** stop clock */
   SCIP_CALL(SCIPstopSaveTime(scip));

   return SCIP_OKAY;
}

/*
 * Returns the number of added constraints
 */
int SCIPbranchruleNodereoptGetNAddedConss(
   SCIP*                 scip,
   SCIP_NODE*            node
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int nodeID;

   assert(scip != NULL);
   assert(node != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   nodeID = SCIPnodeGetReoptID(node);

   if( nodeID >= 1 && branchruledata->nodedata[nodeID]->conss != NULL)
      return MAX( SCIPnodeGetNAddedcons(scip, node), SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss));
   else
      return SCIPnodeGetNAddedcons(scip, node);
}

/*
 * set the candidate for the WC heuristic
 */
void SCIPbranchruleNodereoptSetWCCand(
   SCIP*                 scip,
   SCIP_NODE*            node
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int nodeID;

   assert(scip != NULL);
   assert(node != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( SCIPnodeGetLowerbound(node) > branchruledata->lastpruneddualbound
    && SCIPisFeasLT(scip, fabs(SCIPnodeGetLowerbound(node)), SCIPinfinity(scip)) )
   {
      branchruledata->lastpruneddualbound = SCIPnodeGetLowerbound(node);
      branchruledata->lastprunednr = SCIPnodeGetNumber(node);
      branchruledata->lastprunedID = SCIPnodeGetReoptID(node);
   }
}

/*
 * save unexplored nodes
 */
SCIP_RETCODE SCIPbranchruleNodereoptSaveOpenNodes(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_NODE** opennodes;
   int nodeID;
   int nopennodes;
   int opennode;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPdebugMessage("-> save %d open nodes.\n", SCIPgetNNodesLeft(scip));

   /* save open leaves */
   SCIP_CALL( SCIPgetLeaves(scip, &opennodes, &nopennodes) );
   for(opennode = 0; opennode < nopennodes; opennode++)
   {
      SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, opennodes[opennode], SCIP_REOPTTYPE_PRUNED, TRUE) );
   }

   /* save open children */
   SCIP_CALL( SCIPgetChildren(scip, &opennodes, &nopennodes) );
   for(opennode = 0; opennode < nopennodes; opennode++)
   {
      SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, opennodes[opennode], SCIP_REOPTTYPE_PRUNED, TRUE) );
   }

   /* save open siblings */
   SCIP_CALL( SCIPgetSiblings(scip, &opennodes, &nopennodes) );
   for(opennode = 0; opennode < nopennodes; opennode++)
   {
      SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, opennodes[opennode], SCIP_REOPTTYPE_PRUNED, TRUE) );
   }

   return SCIP_OKAY;
}

/*
 * Save global constraint to separate solution
 */
SCIP_RETCODE SCIPbranchruleNodereoptSaveGlobaleCons(
   SCIP*                 scip,
   SCIP_SOL*             sol,
   SCIP_SET*             set,
   SCIP_STAT*            stat
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int nbinvars;
   int varnr;

   assert(scip != NULL);
   assert(sol != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(branchruledata->init);

   /** realloc memory if necessary */
   if(branchruledata->nconsglobal == branchruledata->allocmemglobalcons)
   {
      SCIP_CALL( resizeGlobalsConsData(scip, branchruledata) );
   }

   nbinvars = SCIPgetNOrigBinVars(scip);

   /** ensure that there exists binary variables */
   if( nbinvars > 0 )
   {
      SCIP_VAR** vars;

      /** allocate memory at position nconsglobal and increase the counter */
      SCIP_CALL( SCIPallocMemory(scip, &branchruledata->globalcons[branchruledata->nconsglobal]) );
      branchruledata->globalcons[branchruledata->nconsglobal]->constype = REOPT_CONSTYPE_SEPASOLUTION;
      branchruledata->globalcons[branchruledata->nconsglobal]->nvars = 0;

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &vars, SCIPgetOrigVars(scip), SCIPgetNOrigVars(scip)) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->globalcons[branchruledata->nconsglobal]->vars, nbinvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->globalcons[branchruledata->nconsglobal]->vals, nbinvars) );

      for(varnr = 0; varnr < SCIPgetNOrigVars(scip); varnr++)
      {
         if( SCIPvarGetType(vars[varnr]) == SCIP_VARTYPE_BINARY )
         {
            SCIP_Real scalar;
            SCIP_Real constant;
            SCIP_Real solval;
            int pos;

            scalar = 1;
            constant = 0;

            pos = branchruledata->globalcons[branchruledata->nconsglobal]->nvars;
            branchruledata->globalcons[branchruledata->nconsglobal]->vars[pos] = vars[varnr];

            solval = SCIPgetSolVal(scip, sol, branchruledata->globalcons[branchruledata->nconsglobal]->vars[pos]);

            if(!SCIPvarIsOriginal(branchruledata->globalcons[branchruledata->nconsglobal]->vars[pos]) )
            {
               SCIP_CALL( SCIPvarGetOrigvarSum( &branchruledata->globalcons[branchruledata->nconsglobal]->vars[pos], &scalar, &constant) );
            }
            assert(SCIPvarGetStatus(branchruledata->globalcons[branchruledata->nconsglobal]->vars[pos]) == SCIP_VARSTATUS_ORIGINAL);

            branchruledata->globalcons[branchruledata->nconsglobal]->vals[pos] = (solval - constant) * scalar;

            branchruledata->globalcons[branchruledata->nconsglobal]->nvars++;
         }
      }

      SCIPfreeMemoryArray(scip, &vars);
   }
   branchruledata->nconsglobal++;

   return SCIP_OKAY;
}

/*
 * Remove node from data structure
 */
SCIP_RETCODE SCIPbranchruleNodereoptRemoveNode(
   SCIP*                 scip,
   SCIP_NODE*            node,
   SCIP_Bool             branched,
   SCIP_Bool             infeasible
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_NODE* parent;
   int parentID;
   int firstID;
   int childID;
   int nodeID;

   assert(SCIPgetRootNode(scip) != node);
   assert(infeasible + branched == 1);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   /* get the nodeID and ensure that everything was hashed correctly */
   nodeID = SCIPnodeGetReoptID(node);
   assert(branchruledata->nodedata[nodeID] != NULL );
   assert((branched && (branchruledata->nodedata[nodeID]->nodechilds == NULL
         || SCIPqueueNElems(branchruledata->nodedata[nodeID]->nodechilds) == 0))
         || infeasible);

   /** find the next saved node above, ensure that the current node is in the child-list, delete the current node from this list */
   parent = SCIPnodeGetParent(node);
   while(parent != SCIPgetRootNode(scip) && SCIPnodeGetReopttype(parent) < SCIP_REOPTTYPE_TRANSIT)
   {
      parent = SCIPnodeGetParent(parent);
   }

   if(parent == SCIPgetRootNode(scip))
      parentID = 0;
   else
   {
      assert(1 <= SCIPnodeGetReoptID(parent) && SCIPnodeGetReoptID(parent) < branchruledata->allocmemsizenodedata);
      parentID = SCIPnodeGetReoptID(parent);
   }
   assert(branchruledata->nodedata[parentID] != NULL );
   assert(branchruledata->nodedata[parentID]->nodechilds != NULL );
   assert(!SCIPqueueIsEmpty(branchruledata->nodedata[parentID]->nodechilds));

   childID = (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[parentID]->nodechilds);
   firstID = childID;
   if (childID != nodeID)
   {
      do
      {
         SCIP_CALL(SCIPqueueInsert(branchruledata->nodedata[parentID]->nodechilds, (void* ) (size_t ) childID));
         childID = (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[parentID]->nodechilds);
         assert(childID != firstID);
      }
      while (childID != nodeID);
   }

   /* if the current node was the only child node of the next saved node above,
    * then this node is infeasible, too.  */
   if( SCIPqueueIsEmpty(branchruledata->nodedata[parentID]->nodechilds)
    && branchruledata->nodedata[parentID]->reopttype != SCIP_REOPTTYPE_STRBRANCHED)
   {
      printf("*** inf node %lld -> remove parent %lld\n", SCIPnodeGetNumber(node), SCIPnodeGetNumber(parent));
      /* if the parent node is the root node, we have to add a global constraint
       * which implies infeasibility */
      if( parent != SCIPgetRootNode(scip) )
      {
         if( infeasible )
            branchruledata->ninfeasnodes++;

         SCIP_CALL( SCIPbranchruleNodereoptRemoveNode(scip, parent, branched, infeasible) );
      }
   }
   else
   {
      /* this current node is the highest infeasible node in this subtree, we can add
       * a global constraint */
      if( infeasible )
         SCIPbranchruleNodereoptInfNode(scip, node);
   }

   /** delete all children below if the node is infeasible */
   if (infeasible && branchruledata->nodedata[nodeID]->nodechilds != NULL )
   {
      while (!SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds))
      {
         SCIP_CALL( deleteChildrenBelow(scip, branchruledata, (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[nodeID]->nodechilds), TRUE) );
      }
   }

   /* delete node data */
   assert(branchruledata->nodedata[nodeID]->nodechilds == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds));
   SCIP_CALL(deleteNodeData(scip, branchruledata, nodeID, TRUE));
   SCIPnodeSetReoptID(node, -1);
   SCIPnodeSetReopttype(node, SCIP_REOPTTYPE_NONE);

   /* insert the empty slot the queue of empty slots */
   SCIP_CALL(SCIPqueueInsert(branchruledata->openIDs, (void*) (size_t ) nodeID));

   return SCIP_OKAY;
}

/*
 *  Save information if infeasible nodes
 */
SCIP_RETCODE SCIPbranchruleNodereoptInfNode(
   SCIP*                 scip,
   SCIP_NODE*            node
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(SCIPgetRootNode(scip) != node);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   if( branchruledata->sepaglbsubtree )
   {
      SCIP_CALL( saveGlobalCons(scip, branchruledata, node, REOPT_CONSTYPE_INFSUBTREE) );
   }

   branchruledata->ninfeasnodesround++;

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPbranchruleNodereoptSolveLP(
   SCIP*                 scip,
   SCIP_NODE*            node,
   SCIP_Bool*            solvelp
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int nodeID;

   assert(scip != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL );

   /* get the ID */
   nodeID = (SCIPgetRootNode(scip) == node ? 0 : SCIPnodeGetReoptID(node));

   (*solvelp) = TRUE;

   if( nodeID == 0 )
   {
      if( branchruledata->nodedata[0]->nodechilds != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds) )
      {
         if( SCIPreoptGetSimToPrevious(scip->reopt) >= branchruledata->objsimrootLP )
            (*solvelp) = FALSE;
      }
   }
   else
      switch (branchruledata->solvelp) {
         /* solve all LPs */
         case 0:
            if( SCIPnodeGetReopttype(node) < SCIP_REOPTTYPE_LEAF )
            {
               if( branchruledata->nodedata[nodeID]->nvars < branchruledata->solvelpdiff)
                  (*solvelp) = FALSE;
            }
            break;

         default:
            if( branchruledata->nodedata[nodeID]->nodechilds != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds) )
            {
               if( branchruledata->nodedata[nodeID]->nvars < branchruledata->solvelpdiff && (int) SCIPnodeGetReopttype(node) < branchruledata->solvelp )
                  (*solvelp) = FALSE;
            }
            break;
      }

   assert(*solvelp || (branchruledata->nodedata[nodeID]->nodechilds != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds)));

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPbranchruleNodereoptRestartCheck(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_Bool run;

   assert(scip != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   /* try heuristics and check for restart */
   if (branchruledata->init && !branchruledata->infeasibleproblem)
   {
      /* tmp for Benjamin */
//      int nodeID;
//      int nleafs;
//      int parentID;
//      int v;
//      SCIP_Real d;
//
//      nleafs = 0;
//      d = 0.0;
//      for(nodeID = 0; nodeID < branchruledata->allocmemsizenodedata; nodeID++)
//      {
//         if( branchruledata->nodedata[nodeID] != NULL && (branchruledata->nodedata[nodeID]->nodechilds == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds)) )
//         {
//            nleafs++;
//            printf("node #%d (%d fixed vars).\n", nodeID, lengthBranchPathByID(branchruledata, nodeID));
//
//            parentID = nodeID;
//            while( parentID != 0 )
//            {
//               for(v = 0; v < branchruledata->nodedata[parentID]->nvars; v++)
//               {
//                  printf("%s %s\n", SCIPvarGetName(branchruledata->nodedata[parentID]->vars[v]), branchruledata->nodedata[parentID]->varboundtypes[v] == SCIP_BOUNDTYPE_UPPER ? "<= 0" : ">= 1");
//                  d++;
//               }
//               parentID = branchruledata->nodedata[parentID]->parentID;
//            }
//            printf("\n");
//         }
//      }
//      printf("%d leaf nodes with avg depth %.2f\n", nleafs, d/nleafs);
      /********************/

      run = FALSE;

      if( branchruledata->lrenable || branchruledata->wcenable )
      {
         SCIP_CALL( runHeuristics(scip, branchruledata, &run) );
      }

      /* update statics
       * we have to update the statistics after running the heuristics
       * because we use this statistic to decide of a heuristic should run
       * or not */
      updateStatistics(scip, branchruledata);

      /* if no heuristic ran, we check if need to restart due to
       * exceeding the limit of saved nodes */
      if( !run )
      {
         /* check if a restart should be forced */
         restartCheck(scip, branchruledata);

         if (branchruledata->restart)
         {
            /* clear nodedate */
            SCIP_CALL(clearNodes(scip, branchruledata, FALSE));

            /* clear saved information in branch_pseudo */
            SCIP_CALL(SCIPbranchrulePseudoReset(scip, TRUE, FALSE));

            /* init data for the root node */
            SCIP_CALL( initNode(scip, branchruledata, 0) );

            /* reset flag */
            branchruledata->restart = FALSE;

            /* decrease number of restarts */
            branchruledata->nrestarts++;

            /* set the number of saved nodes to 0 */
            branchruledata->nsavednodes = 0;
         }
         /* increase the localdelay by 0.5% if node subproblem was solved from scratch */
         else if( branchruledata->ndelsubtrees == 0 && branchruledata->dynamiclocaldelay)
         {
            SCIPdebugMessage("increase local delay: %.4f -> %.4f\n", branchruledata->localdelay, MIN(1.005*branchruledata->localdelay, 0.99));
            branchruledata->localdelay = MIN(1.005*branchruledata->localdelay, 0.99);
            branchruledata->ndelsubtrees = -1;
         }
      }
   }

   return SCIP_OKAY;
}

/* get all statistic information */
SCIP_RETCODE SCIPbranchruleNodereoptGetStatistic(
   SCIP*                 scip,
   int*                  nfeasnodes,
   int*                  nfeasnodesround,
   int*                  ninfeasnodes,
   int*                  ninfeasnodesround,
   int*                  nprunednodes,
   int*                  nprunednodesround,
   int*                  nrediednodes,
   int*                  nrediednodesround,
   int*                  nruns,
   int*                  nrestarts,
   int*                  firstrestart,
   int*                  lastrestart,
   int*                  nrestartsround,
   int*                  ninfsubtrees,
   int*                  liscalls,
   int*                  lissucces,
   int*                  lisk,
   SCIP_Real*            listime,
   int*                  lccalls,
   int*                  lcsucces,
   int*                  lck,
   SCIP_Real*            lctime,
   SCIP_Real*            lptime,
   SCIP_Real*            lptime_maxthread
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->reopt);

   /* feasible nodes */
   if( nfeasnodes != NULL)
      (*nfeasnodes) = branchruledata->nfeasnodes + branchruledata->nfeasnodesround;

   if( nfeasnodesround != NULL)
      (*nfeasnodesround) = branchruledata->nfeasnodesround;

   /* infeasible nodes */
   if( ninfeasnodes != NULL)
      (*ninfeasnodes) = branchruledata->ninfeasnodes + branchruledata->ninfeasnodesround;

   if( ninfeasnodesround != NULL)
      (*ninfeasnodesround) = branchruledata->ninfeasnodesround;

   /* pruned nodes */
   if( nprunednodes != NULL)
      (*nprunednodes) = branchruledata->nprunednodes + branchruledata->nprunednodesround;

   if( nprunednodesround != NULL)
      (*nprunednodesround) = branchruledata->nprunednodesround;

   /* redied nodes */
   if( nrediednodes != NULL)
      (*nrediednodes) = branchruledata->nrediednodes + branchruledata->nrediednodesround;

   if( nrediednodesround != NULL)
      (*nrediednodesround) = branchruledata->nrediednodesround;

   /* other stuff */
   if( nruns != NULL)
      (*nruns) = branchruledata->nruns;

   if( nrestarts != NULL)
      (*nrestarts) = branchruledata->nrestarts;

   if( firstrestart != NULL)
      (*firstrestart) = branchruledata->firstrestart;

   if( lastrestart != NULL)
      (*lastrestart) = branchruledata->lastrestart;

   if( nrestartsround != NULL)
      (*nrestartsround) = branchruledata->nrestartsround;

   if( ninfsubtrees != NULL )
      (*ninfsubtrees) = branchruledata->ninfsubtrees;

   if( liscalls != NULL )
      (*liscalls) = branchruledata->lrcalls;

   if( lissucces != NULL )
      (*lissucces) = branchruledata->lrsuccess;

   if( lisk != NULL )
      (*lisk) = branchruledata->lrk;

   if( listime != NULL )
      (*listime) = SCIPgetClockTime(scip, branchruledata->lrtime);

   if( lccalls != NULL )
      (*lccalls) = branchruledata->wccalls;

   if( lcsucces != NULL )
      (*lcsucces) = branchruledata->wcsuccess;

   if( lck != NULL )
      (*lck) = branchruledata->wck;

   if( lctime != NULL )
      (*lctime) = SCIPgetClockTime(scip, branchruledata->wctime);

   if( lptime != NULL )
      (*lptime) = SCIPgetClockTime(scip, branchruledata->lpclock);

   if( lptime_maxthread != NULL )
   {
      int thread;
      (*lptime_maxthread) = 0;
      for(thread = 0; thread < 8; thread++)
      {
         (*lptime_maxthread) = MAX(*lptime_maxthread, SCIPgetClockTime(scip, branchruledata->threadclock[thread]));
      }
   }

   return SCIP_OKAY;
}

/* start clock for update solutions */
SCIP_RETCODE SCIPbranchruleNodereoptStartUpdatesoluTime(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->reopt);

   SCIP_CALL( SCIPstartClock(scip, branchruledata->updatesolutime) );

   return SCIP_OKAY;
}

/* stop clock for update solutions */
SCIP_RETCODE SCIPbranchruleNodereoptStopUpdatesoluTime(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->reopt);

   SCIP_CALL( SCIPstopClock(scip, branchruledata->updatesolutime) );

   return SCIP_OKAY;
}

SCIP_Real SCIPbranchruleNodereoptGetCutoffbound(
   SCIP*                 scip,
   int                   nodeID
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->reopt);

   return getCutoffbound(scip, branchruledata, nodeID);
}

SCIP_RETCODE SCIPbranchruledataNodereoptLPTimeStart(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->reopt);

   SCIP_CALL( SCIPstartClock(scip, branchruledata->lpclock) );
   SCIP_CALL( SCIPstartClock(scip, branchruledata->threadclock[branchruledata->nlp%8]) );

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPbranchruledataNodereoptLPTimeStop(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->reopt);

   SCIP_CALL( SCIPstopClock(scip, branchruledata->lpclock) );
   SCIP_CALL( SCIPstopClock(scip, branchruledata->threadclock[branchruledata->nlp%8]) );
   branchruledata->nlp++;

   return SCIP_OKAY;
}

/* remove LP branching candidates that are already strong branched  */
static
SCIP_RETCODE checkLPBranchCands(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   SCIP_VAR**            branchcands,
   SCIP_Real*            branchcandssol,
   SCIP_Real*            branchcandsfrac,
   int*                  nbranchcands
)
{
   SCIP_VAR* transvar;
   int cand;
   int var;
   int ndelcands;
   int* delcands;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(0 <= nodeID && nodeID < branchruledata->allocmemsizenodedata);
   assert(*nbranchcands == 0 || branchcands != NULL);
   assert(*nbranchcands == 0 || branchcandsfrac != NULL);
   assert(*nbranchcands == 0 || nbranchcands != NULL);

   if( *nbranchcands == 0 || branchruledata->nodedata[nodeID]->dualcons[0] == NULL || branchruledata->nodedata[nodeID]->dualcons[0]->nvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &delcands, *nbranchcands) );

   ndelcands = 0;
   for(cand = 0; cand < *nbranchcands; cand++)
   {
      SCIP_VAR* origvar;
      SCIP_Real constant;
      SCIP_Real scalar;
      SCIP_Bool entering;
      SCIP_Bool leaving;
      int id;

      origvar = branchcands[cand];
      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );
      assert(SCIPvarIsOriginal(origvar));

      id = SCIPvarGetIndex(origvar);
      SCIPreoptEnterOrLeaveObj(scip->reopt, id, &entering, &leaving);

      if( SCIPreoptIsObjCoefNegated(scip->reopt, id) )
      {
         for(var = 0; var < branchruledata->nodedata[nodeID]->dualcons[0]->nvars; var++)
         {
            transvar = SCIPvarGetTransVar(branchruledata->nodedata[nodeID]->dualcons[0]->vars[var]);
            if( SCIPvarGetIndex(branchcands[cand]) == SCIPvarGetIndex(transvar) )
            {
               delcands[ndelcands] = cand;
               ndelcands++;
            }
         }
      }
      else
      {
         delcands[ndelcands] = cand;
         ndelcands++;
      }
   }

   if( ndelcands > 0 )
   {
      ndelcands--;
      while( ndelcands >= 0 )
      {
         cand = delcands[ndelcands];
         if( cand < *nbranchcands-1 )
         {
            branchcands[cand] = branchcands[*nbranchcands-1];
            branchcandssol[cand] = branchcandssol[*nbranchcands-1];
            branchcandsfrac[cand] = branchcandsfrac[*nbranchcands-1];
         }
         *nbranchcands -= 1;
         ndelcands--;
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */

#define branchCopynodereopt NULL;

static
SCIP_DECL_BRANCHEXIT(branchExitnodereopt)
{
   int var;

   /* release and unlock all variables */
   for(var = 0; var < SCIPgetNOrigVars(scip); var++)
   {
      SCIP_VAR* transvar;
      transvar =  SCIPvarGetTransVar(SCIPgetOrigVars(scip)[var]);
      SCIP_CALL( SCIPaddVarLocks(scip, transvar, -1, -1) );
      SCIP_CALL( SCIPreleaseVar(scip, &transvar) );
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_BRANCHINIT(branchInitnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   int var;

   assert(scip != NULL );
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL );

   /** HACK */
   /** check if all variable are binary, if not, disable reoptimization */
   if (!branchruledata->init)
   {
      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/enable", &branchruledata->reopt));
      if (branchruledata->reopt && SCIPgetNOrigIntVars(scip) > 0)
      {
         branchruledata->reopt = FALSE;
         printf(">> reoptimization is disabled, because the problem has non binary variables.\n");
      }
   }

   /** initialize the data and change parameters */
   if (!branchruledata->init && branchruledata->reopt)
   {
      int slot;

      /** parameter stuff */
      branchruledata->restart = FALSE;
      SCIP_CALL(SCIPgetIntParam(scip, "reoptimization/maxsavednodes", &branchruledata->maxsavednodes));

      if( branchruledata->maxsavednodes == -1 )
      {
         branchruledata->maxsavednodes = floor(0.5319 * SCIPgetNBinVars(scip));
      }

      SCIP_CALL(SCIPgetIntParam(scip, "reoptimization/maxdiffofnodes", &branchruledata->maxdiffofnodes));

      if( branchruledata->maxdiffofnodes == -1 )
      {
         branchruledata->maxdiffofnodes = floor(log10((SCIP_Real)SCIPgetNBinVars(scip))/log10(2))-1;
      }

      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/dynamicdiffofnodes", &branchruledata->dynamicdiffofnodes));
      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/savelpbasis", &branchruledata->savelpbasis));
      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/globalcons/sepainfsubtrees", &branchruledata->sepaglbsubtree));
      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/globalcons/sepasols", &branchruledata->sepaglbsols));
      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/localcons/sepasols", &branchruledata->sepalocsols));
      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/sepabestsol", &branchruledata->sepabestsol));

      SCIP_CALL(SCIPgetRealParam(scip, "reoptimization/objsimrootLP", &branchruledata->objsimrootLP));

      SCIP_CALL(SCIPgetIntParam(scip, "reoptimization/solvelp", &branchruledata->solvelp));
      SCIP_CALL(SCIPgetIntParam(scip, "reoptimization/solvelpdiff", &branchruledata->solvelpdiff));

      SCIP_CALL(SCIPgetRealParam(scip, "reoptimization/delay", &branchruledata->localdelay));

      /** Savings for reoptimization */
      branchruledata->allocmemsizenodedata = MIN(SCIPgetNOrigBinVars(scip), branchruledata->maxsavednodes) + 1;
      branchruledata->nsavednodes = 0;

      branchruledata->nfeasnodes = 0;
      branchruledata->nfeasnodesround = 0;
      branchruledata->nfeasnodesmax = 0;

      branchruledata->ninfeasnodes = 0;
      branchruledata->ninfeasnodesround = 0;
      branchruledata->ninfeasnodesmax = 0;

      branchruledata->nprunednodes = 0;
      branchruledata->nprunednodesround = 0;
      branchruledata->nprunednodesmax = 0;

      branchruledata->npbnodes = 0;
      branchruledata->npbnodesround = 0;
      branchruledata->npbnodesmax = 0;

      branchruledata->nrediednodes = 0;
      branchruledata->nrediednodesround = 0;
      branchruledata->nrediednodesmax = 0;

      branchruledata->nrestarts = 0;
      branchruledata->firstrestart = -1;
      branchruledata->lastrestart = -1;
      branchruledata->nrestartsround = 0;
      branchruledata->nruns = 0;

      branchruledata->ninfsubtrees = 0;
      branchruledata->nrevivednodes = 0;

      branchruledata->infeasibleproblem = FALSE;
      branchruledata->infeasibleconsadded = FALSE;

      /** globals constraints */
      branchruledata->nconsglobal = 0;
      branchruledata->nconsglobaladded = 0;
      branchruledata->allocmemglobalcons = 50;

      /** heuristics */
      SCIP_CALL(SCIPclockCreate(&branchruledata->lrtime, SCIP_CLOCKTYPE_DEFAULT));
      branchruledata->lrsuccess = 0;
      branchruledata->lrcalls = 0;
      branchruledata->lrdepth = 0;
      branchruledata->lrk = 0;
      branchruledata->lrloss = SCIPinfinity(scip);

      SCIP_CALL(SCIPclockCreate(&branchruledata->wctime, SCIP_CLOCKTYPE_DEFAULT));
      branchruledata->wcsuccess= 0;
      branchruledata->wccalls = 0;
      branchruledata->wck = 0;
      branchruledata->wcdepth = 1;
      branchruledata->wclastnodeID = -1;
      branchruledata->lastprunedID = -1;
      branchruledata->lastprunednr = -1;
      branchruledata->lastpruneddualbound = -SCIPinfinity(scip);

      /** clock */
      SCIP_CALL(SCIPclockCreate(&branchruledata->inittime, SCIP_CLOCKTYPE_DEFAULT));
      SCIP_CALL(SCIPclockCreate(&branchruledata->updatesolutime, SCIP_CLOCKTYPE_DEFAULT));
      SCIP_CALL(SCIPclockCreate(&branchruledata->savetime, SCIP_CLOCKTYPE_DEFAULT));

      /** hypothetical multi-threading */
      SCIP_CALL(SCIPallocMemoryArray(scip, &branchruledata->threadclock, 8));
      for(slot = 0; slot < 8; slot++)
      {
         SCIP_CALL(SCIPclockCreate(&branchruledata->threadclock[slot], SCIP_CLOCKTYPE_DEFAULT));
      }
      SCIP_CALL(SCIPclockCreate(&branchruledata->lpclock, SCIP_CLOCKTYPE_DEFAULT));
      branchruledata->nlp = 0;

      /* allocate data if nodes shall saved */
      if( branchruledata->maxsavednodes > 0 )
      {
         /** data structure for nodes, init data for root node at position 0 */
         SCIP_CALL(SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata), branchruledata->allocmemsizenodedata));

         /** allocate memory for the root */
         SCIP_CALL( initNode(scip, branchruledata, 0) );

         /** data structure to handle IDs */
         SCIP_CALL(SCIPqueueCreate(&(branchruledata->openIDs), branchruledata->allocmemsizenodedata - 1, 2));

         /** fill the queue with free IDs 1,...,allocmemsizenodedata-1; the ID 0 is always reserved for the root node */
         slot = 1;
         while (slot < branchruledata->allocmemsizenodedata)
         {
            SCIP_CALL(SCIPqueueInsert(branchruledata->openIDs, (void* ) (size_t ) slot));
            slot++;
         }
         assert(SCIPqueueNElems(branchruledata->openIDs) == branchruledata->allocmemsizenodedata - 1);
      }

      /** data for global logic-or constraints to separate binary solutions */
      SCIP_CALL(SCIPallocClearMemoryArray(scip, &(branchruledata->globalcons), branchruledata->allocmemglobalcons));

      /** mark data structure initialized */
      branchruledata->init = TRUE;
   }

   branchruledata->lastprunedID = -1;
   branchruledata->lastpruneddualbound = -SCIPinfinity(scip);

   branchruledata->ndelsubtrees = -1;

   /** change parameters */
   /* propagator */
   SCIP_CALL(SCIPsetIntParam(scip, "propagating/dualfix/freq", -1));
   SCIP_CALL(SCIPsetIntParam(scip, "propagating/dualfix/maxprerounds", 0));

   SCIP_CALL(SCIPsetIntParam(scip, "propagating/genvbounds/freq", -1));
   SCIP_CALL(SCIPsetIntParam(scip, "propagating/genvbounds/maxprerounds", 0));

   SCIP_CALL(SCIPsetIntParam(scip, "propagating/obbt/freq", -1));
   SCIP_CALL(SCIPsetIntParam(scip, "propagating/obbt/maxprerounds", 0));

   SCIP_CALL(SCIPsetIntParam(scip, "propagating/pseudoobj/freq", -1));
   SCIP_CALL(SCIPsetIntParam(scip, "propagating/pseudoobj/maxprerounds", 0));

   SCIP_CALL(SCIPsetIntParam(scip, "propagating/redcost/freq", -1));
   SCIP_CALL(SCIPsetIntParam(scip, "propagating/redcost/maxprerounds", 0));

   SCIP_CALL(SCIPsetIntParam(scip, "propagating/rootredcost/freq", -1));
   SCIP_CALL(SCIPsetIntParam(scip, "propagating/rootredcost/maxprerounds", 0));

   SCIP_CALL(SCIPsetIntParam(scip, "propagating/vbounds/freq", 1));
   SCIP_CALL(SCIPsetIntParam(scip, "propagating/vbounds/maxprerounds", -1));

   /* constraint handler */
   SCIP_CALL(SCIPsetBoolParam(scip, "constraints/abspower/dualpresolve", FALSE));
   SCIP_CALL(SCIPsetBoolParam(scip, "constraints/and/dualpresolving", FALSE));
   SCIP_CALL(SCIPsetBoolParam(scip, "constraints/indicator/dualreductions", FALSE));
   SCIP_CALL(SCIPsetBoolParam(scip, "constraints/knapsack/dualpresolving", FALSE));
   SCIP_CALL(SCIPsetBoolParam(scip, "constraints/linear/dualpresolving", FALSE));
   SCIP_CALL(SCIPsetBoolParam(scip, "constraints/logicor/dualpresolving", FALSE));
   SCIP_CALL(SCIPsetBoolParam(scip, "constraints/setppc/dualpresolving", FALSE));

   SCIP_CALL(SCIPsetBoolParam(scip, "conflict/enable", FALSE));

   /* capture and lock all variables */
   for(var = 0; var < SCIPgetNOrigVars(scip); var++)
   {
      SCIP_VAR* transvar;
      transvar =  SCIPvarGetTransVar(SCIPgetOrigVars(scip)[var]);
      SCIP_CALL( SCIPcaptureVar(scip, transvar) );
      SCIP_CALL( SCIPaddVarLocks(scip, transvar, 1, 1) );
   }

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreenodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL );
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL );

   /** free data structure only if reoptimization is enabeld */
   if (branchruledata->init )
   {
      int nodeID;
      int t;

      assert(branchruledata->reopt);

      if( branchruledata->maxsavednodes > 0 )
      {
         SCIP_CALL( clearNodes(scip, branchruledata, TRUE) );

         /** free data for handle saved nodes */
         SCIPqueueFree(&branchruledata->openIDs);
         SCIPfreeMemoryArray(scip, &branchruledata->nodedata);
      }

      /** free clock */
      SCIPfreeClock(scip, &branchruledata->inittime);
      SCIPfreeClock(scip, &branchruledata->updatesolutime);
      SCIPfreeClock(scip, &branchruledata->savetime);
      SCIPfreeClock(scip, &branchruledata->lrtime);
      SCIPfreeClock(scip, &branchruledata->wctime);

      SCIPfreeClock(scip, &branchruledata->lpclock);
      for(t = 0; t < 8; t++)
         SCIPfreeClock(scip, &branchruledata->threadclock[t]);
      SCIPfreeMemoryArray(scip, &branchruledata->threadclock);

      /** free data for global constraints */
      for (nodeID = 0; nodeID < branchruledata->nconsglobal; nodeID++)
      {
         if (branchruledata->globalcons != NULL )
         {
            SCIPfreeMemoryArray(scip, &branchruledata->globalcons[nodeID]->vars);
            SCIPfreeMemoryArray(scip, &branchruledata->globalcons[nodeID]->vals);
            SCIPfreeMemory(scip, &branchruledata->globalcons[nodeID]);
         }
      }
      SCIPfreeMemoryArray(scip, &branchruledata->globalcons);

      branchruledata->init = FALSE;
      branchruledata->reopt = FALSE;
   }
   assert(!branchruledata->init);

   SCIPfreeMemory(scip, &branchruledata);

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   SCIP_Real* branchcandsfrac;
   int nbranchcands;

   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   if ( branchruledata->reopt )
   {
      if ((SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip) && branchruledata->nodedata[0]->pseudobranched)
       || (SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip)
         && branchruledata->nodedata[0]->nodechilds != NULL
         && !SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds))
       || (SCIPgetCurrentNode(scip) != SCIPgetRootNode(scip)
         && SCIPnodeGetReopttype(SCIPgetCurrentNode(scip)) >= SCIP_REOPTTYPE_TRANSIT
         && branchruledata->nodedata[SCIPnodeGetReoptID(SCIPgetCurrentNode(scip))]->nodechilds != NULL
         && !SCIPqueueIsEmpty(branchruledata->nodedata[SCIPnodeGetReoptID(SCIPgetCurrentNode(scip))]->nodechilds)))
      {
         /** start time */
         SCIP_CALL( SCIPstartInitTime(scip) );
         SCIP_CALL( SCIPstartClock(scip, branchruledata->inittime));

         if( branchruledata->strongbranchinginit
          && SCIPreoptGetSimToPrevious(scip->reopt) < branchruledata->objsimrootLP
          && SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip))
         {
            /* get branching candidates */
            SCIP_CALL( SCIPgetLPBranchCands(scip, &branchcands, &branchcandssol, &branchcandsfrac, NULL, &nbranchcands, NULL) );

            /* run strong branching initialization */
            if( nbranchcands > 0 )
            {
               /* select only some 'good' candidates */
               SCIP_CALL( checkLPBranchCands(scip, branchruledata, 0, branchcands, branchcandssol, branchcandsfrac, &nbranchcands) );

               if( nbranchcands > 0 )
               {
                  SCIP_CALL( SCIPexecRelpscostBranching(scip, TRUE, branchcands, branchcandssol, branchcandsfrac, nbranchcands, FALSE, result) );
               }
               assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM);
            }
         }

         if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM)
         {
            if (SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip))
            {
               SCIP_CALL(Exec(scip, branchruledata, result));
               if( *result == SCIP_BRANCHED )
               {
                  if( SCIPbranchrulePseudoIsPseudoBranched(scip, SCIPgetCurrentNode(scip)) )
                  {
                     SCIP_CALL(SCIPbranchrulePseudoReset(scip, FALSE, TRUE));
                  }
                  else
                  {
                     SCIP_CALL(SCIPbranchrulePseudoReset(scip, FALSE, FALSE));
                  }
               }
            }
            else
            {
               assert(1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));
               assert(SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)) < branchruledata->allocmemsizenodedata);
               SCIP_CALL(Exec(scip, branchruledata, result));
            }
         }

         /** stop clock */
         SCIP_CALL( SCIPstopInitTime(scip) );
         SCIP_CALL( SCIPstopClock(scip, branchruledata->inittime));
      }
   }

   return SCIP_OKAY;
}

/** branching execution method for external candidates */
static SCIP_DECL_BRANCHEXECEXT(branchExecextnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   SCIP_Real* branchcandsfrac;
   int nbranchcands;

   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   if (branchruledata->reopt)
   {
      if ((SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip) && branchruledata->nodedata[0]->pseudobranched)
       || (SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip)
         && branchruledata->nodedata[0]->nodechilds != NULL
         && !SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds))
       || (SCIPgetCurrentNode(scip) != SCIPgetRootNode(scip)
         && SCIPnodeGetReopttype(SCIPgetCurrentNode(scip)) >= SCIP_REOPTTYPE_TRANSIT
         && branchruledata->nodedata[SCIPnodeGetReoptID(SCIPgetCurrentNode(scip))]->nodechilds != NULL
         && !SCIPqueueIsEmpty(branchruledata->nodedata[SCIPnodeGetReoptID(SCIPgetCurrentNode(scip))]->nodechilds)))
      {
         /** start time */
         SCIP_CALL( SCIPstartInitTime(scip) );
         SCIP_CALL( SCIPstartClock(scip, branchruledata->inittime));

         if( branchruledata->strongbranchinginit
          && SCIPreoptGetSimToPrevious(scip->reopt) < branchruledata->objsimrootLP
          && SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip) )
         {
            /* get branching candidates */
            SCIP_CALL( SCIPgetLPBranchCands(scip, &branchcands, &branchcandssol, &branchcandsfrac, NULL, &nbranchcands, NULL) );

            /* run strong branching initialization */
            if( nbranchcands > 0 )
            {
               /* select only some 'good' candidates */
               SCIP_CALL( checkLPBranchCands(scip, branchruledata, 0, branchcands, branchcandssol, branchcandsfrac, &nbranchcands) );

               if( nbranchcands > 0 )
               {
                  SCIP_CALL( SCIPexecRelpscostBranching(scip, TRUE, branchcands, branchcandssol, branchcandsfrac, nbranchcands, FALSE, result) );
               }
               assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM);
            }
         }

         if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM)
         {
            if (SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip))
            {
               SCIP_CALL(Exec(scip, branchruledata, result));
               if( *result == SCIP_BRANCHED )
               {
                  if( SCIPbranchrulePseudoIsPseudoBranched(scip, SCIPgetCurrentNode(scip)) )
                  {
                     SCIP_CALL(SCIPbranchrulePseudoReset(scip, FALSE, TRUE));
                  }
                  else
                  {
                     SCIP_CALL(SCIPbranchrulePseudoReset(scip, FALSE, FALSE));
                  }
               }
            }
            else
            {
               assert(1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));
               assert(SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)) < branchruledata->allocmemsizenodedata);
               SCIP_CALL(Exec(scip, branchruledata, result));
            }
         }

         /** stop clock */
         SCIP_CALL( SCIPstopInitTime(scip) );
         SCIP_CALL( SCIPstopClock(scip, branchruledata->inittime));
      }
   }

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static SCIP_DECL_BRANCHEXECPS(branchExecpsnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   SCIP_Real* branchcandsfrac;
   int nbranchcands;

   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   if (branchruledata->reopt)
   {
      if ((SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip) && branchruledata->nodedata[0]->pseudobranched)
       || (SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip)
         && branchruledata->nodedata[0]->nodechilds != NULL
         && !SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds))
       || (SCIPgetCurrentNode(scip) != SCIPgetRootNode(scip)
         && SCIPnodeGetReopttype(SCIPgetCurrentNode(scip)) >= SCIP_REOPTTYPE_TRANSIT
         && branchruledata->nodedata[SCIPnodeGetReoptID(SCIPgetCurrentNode(scip))]->nodechilds != NULL
         && !SCIPqueueIsEmpty(branchruledata->nodedata[SCIPnodeGetReoptID(SCIPgetCurrentNode(scip))]->nodechilds)))
      {
         /** start time */
         SCIP_CALL( SCIPstartInitTime(scip) );
         SCIP_CALL( SCIPstartClock(scip, branchruledata->inittime));

         if( branchruledata->strongbranchinginit
          && SCIPreoptGetSimToPrevious(scip->reopt) < branchruledata->objsimrootLP
          && SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip) )
         {
            /* get branching candidates */
            SCIP_CALL( SCIPgetLPBranchCands(scip, &branchcands, &branchcandssol, &branchcandsfrac, NULL, &nbranchcands, NULL) );

            /* run strong branching initialization */
            if( nbranchcands > 0 )
            {
               /* select only some 'good' candidates */
               SCIP_CALL( checkLPBranchCands(scip, branchruledata, 0, branchcands, branchcandssol, branchcandsfrac, &nbranchcands) );

               if( nbranchcands > 0 )
               {
                  SCIP_CALL( SCIPexecRelpscostBranching(scip, TRUE, branchcands, branchcandssol, branchcandsfrac, nbranchcands, FALSE, result) );
               }
               assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM);
            }
         }

         if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM)
         {
            if (SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip))
            {
               SCIP_CALL(Exec(scip, branchruledata, result));
               if( *result == SCIP_BRANCHED )
               {
                  if( SCIPbranchrulePseudoIsPseudoBranched(scip, SCIPgetCurrentNode(scip)) )
                  {
                     SCIP_CALL(SCIPbranchrulePseudoReset(scip, FALSE, TRUE));
                  }
                  else
                  {
                     SCIP_CALL(SCIPbranchrulePseudoReset(scip, FALSE, FALSE));
                  }
               }
            }
            else
            {
               assert(1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));
               assert(SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)) < branchruledata->allocmemsizenodedata);
               SCIP_CALL(Exec(scip, branchruledata, result));
            }
         }

         /** stop clock */
         SCIP_CALL( SCIPstopInitTime(scip) );
         SCIP_CALL( SCIPstopClock(scip, branchruledata->inittime));
      }
   }

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the nodereopt branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleNodereopt(
   SCIP*                 scip                     /**< SCIP data structure */
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL );

   /* create nodereopt branching rule data */
   SCIP_CALL(SCIPallocMemory(scip, &branchruledata));
   branchruledata->init = FALSE;

   /* include nodereopt branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC,
         BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata));

   assert(branchrule != NULL );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL(SCIPsetBranchruleFree(scip, branchrule, branchFreenodereopt));
   SCIP_CALL(SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpnodereopt));
   SCIP_CALL(SCIPsetBranchruleExecExt(scip, branchrule, branchExecextnodereopt));
   SCIP_CALL(SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsnodereopt));
   SCIP_CALL(SCIPsetBranchruleInit(scip, branchrule, branchInitnodereopt));
   SCIP_CALL(SCIPsetBranchruleExit(scip, branchrule, branchExitnodereopt));

   /* parameter */

   SCIP_CALL( SCIPaddBoolParam(scip, "reoptimization/shrinktransit", "upgrade the child nodes of transit nodes to the next save node above.",
         &branchruledata->shrink, TRUE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "reoptimization/strongbranchinginit", "use strong branching at the root node if we reoptimize the root?",
         &branchruledata->strongbranchinginit, TRUE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "reoptimization/reducetofrontier", "delete stored nodes which were not revived.",
         &branchruledata->reducetofrontier, TRUE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "reoptimization/LargestSubtree", "enable largest subtree heuristic.",
         &branchruledata->lrenable, TRUE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "reoptimization/LazyCompr", "enable weak compression heuristic.",
         &branchruledata->wcenable, TRUE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "reoptimization/compressnodes", "compress 0: all, 1: pruned, 2: feasible nodes.",
         &branchruledata->cpressnodes, TRUE, 2, 0, 2, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "reoptimization/dynlocdelay", "increase the local delay by 0.5% if no subproblem was restarted.",
         &branchruledata->dynamiclocaldelay, TRUE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
