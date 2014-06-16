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
//#define DEBUG_MODE

#include <assert.h>
#include "scip/branch_pseudo.h"
#include "scip/branch_nodereopt.h"
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

#define BRANCHRULE_NAME            "nodereopt"
#define BRANCHRULE_DESC            "branching rule for node reoptimization"
#define BRANCHRULE_PRIORITY        536870911
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

/*
 * Data structures
 */

/** Data for logic-or constraints to split nodes into two disjoint parts */

enum Reopt_ConsType
{
   REOPT_CONSTYPE_SEPA = 0,
   REOPT_CONSTYPE_PSEUDO = 1
};
typedef enum Reopt_ConsType REOPT_CONSTYPE;

struct LogicOrData
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   REOPT_CONSTYPE constype;
   int nvars;
};

/** Data for ancestor branching paths */
struct NodeData
{
   SCIP_VAR**            vars;
   SCIP_Real*            varbounds;
   SCIP_BOUNDTYPE*       varboundtypes;
   SCIP_Bool             pseudobranched;
   SCIP_QUEUE*           nodechilds;
   SCIP_LPISTATE*        lpistate;
   SCIP_QUEUE*           conss;
   int                   nvars;
   SCIP_Longint          oldnumber;
   SCIP_REOPTTYPE        reopttype;
   int                   redies;
   int                   allocmem;
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
   SCIP_Bool             saveglbcons;             /** save global constraints */
   int                   nconsglobal;             /** length of the two lists above */
   int                   nconsglobaladded;        /** count number of added constraints (only for statistic) */
   int                   allocmemglobalcons;      /** allocated memory in the two lists above */

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
   int                   nruns;                   /** number of restarts */

   /** clock */
   SCIP_CLOCK*           inittime;                /** time needed to initialize */
   SCIP_CLOCK*           updatesolutime;          /** time needed to update solutions */
   SCIP_CLOCK*           savetime;                /** time needed to save nodes */

   /** Parameter stuff */
   SCIP_Bool             reopt;                   /** boolean flag if the reoptimization should enabled. */
   SCIP_Bool             restart;                 /** should the problem be restarted? */
   SCIP_Bool             savelpbasis;             /** save the LP basis for feasible and strong-branched nodes */
   SCIP_Bool             dynamicdiffofnodes;      /** should the maximal number of bound changes in two ancestor nodes
                                                      calculated depending on number of variables? */
   int                   maxsavednodes;           /** maximal number of saved nodes for the reoptimization */
   int                   maxdiffofnodes;          /** maximal number of bound changes in two ancestor nodes such that
                                                      the path get not shrunk */
   int                   solvelp;                 /** strategy for solving LPs (0: only at leafs, 1: at each node,
                                                      2: only below root, 3: difference previous node on the root path) */
   int                   solvelpdiff;             /** difference is solvelp = 3 */
};

/*
 *  static methods
 */
/** TODO: Braucht man das überhaupt, wenn Knoten mit PBCons und SepaCons immer gespeichert bleiben?! */
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
         consdata->constype = REOPT_CONSTYPE_PSEUDO;

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

/*
 * Reallocate memory for links between IDs and node numbers.
 */
static
SCIP_RETCODE reallocNodeID(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata)
{
   int oldsize;
   int pos;

   assert(scip != NULL );
   assert(branchruledata != NULL );

   oldsize = branchruledata->allocmemsizeNodeID;
   branchruledata->allocmemsizeNodeID *= 2;

   /** data of node information */
   SCIP_CALL( SCIPreallocMemoryArray(scip, &branchruledata->nodetoid, branchruledata->allocmemsizeNodeID) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &branchruledata->currentleafIDs, branchruledata->allocmemsizeNodeID) );

   /** write 0 in empty slots and fill the queue */
   for(pos = oldsize; pos < branchruledata->allocmemsizeNodeID; pos++)
   {
      branchruledata->nodetoid[pos] = 0;
      branchruledata->currentleafIDs[pos] = 0;
   }

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
      branchruledata->nodedata[nodeID]->oldnumber = 0;
      branchruledata->nodedata[nodeID]->pseudobranched = FALSE;
      branchruledata->nodedata[nodeID]->redies = 0;
      branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_NONE;
      branchruledata->nodedata[nodeID]->allocmem = 0;
      branchruledata->nodedata[nodeID]->vars = NULL;
      branchruledata->nodedata[nodeID]->varbounds = NULL;
      branchruledata->nodedata[nodeID]->varboundtypes = NULL;
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
   assert(nodeID >= 1);
   assert(branchruledata->nodedata[nodeID] != NULL );

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
      branchruledata->nodedata[nodeID]->allocmem = SCIPgetNVars(scip);
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[nodeID]->vars), branchruledata->nodedata[nodeID]->allocmem) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[nodeID]->varbounds), branchruledata->nodedata[nodeID]->allocmem) );
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata[nodeID]->varboundtypes), branchruledata->nodedata[nodeID]->allocmem) );
   }

   SCIPnodeGetAncestorBranchingsReopt(node, parent, branchruledata->nodedata[nodeID]->vars,
         branchruledata->nodedata[nodeID]->varbounds,
         branchruledata->nodedata[nodeID]->varboundtypes, &nbranchvars,
         branchruledata->nodedata[nodeID]->allocmem);

   assert(nbranchvars <= branchruledata->nodedata[nodeID]->allocmem); /* this should be the case */

   branchruledata->nodedata[nodeID]->nvars = nbranchvars;

   assert(nbranchvars <= branchruledata->nodedata[nodeID]->allocmem);
   assert(branchruledata->nodedata[nodeID]->vars != NULL );

   /** transform the bounds of branched variables and ensure that they are original. */
   SCIP_CALL( transformIntoOrig(scip, branchruledata, nodeID) );

   return SCIP_OKAY;
}

/*
 * checks if the last added logic-or constraints dominates
 * a previous added constraint.
 */
static
SCIP_RETCODE findMinimalConss(
   SCIP*                 scip,
   int                   nodeID
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   LOGICORDATA** consdatas;
   int nconss;
   int cons1;
   int cons2;
   int varnr1;
   int varnr2;

   assert(scip != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(nodeID > 0);
   assert(branchruledata->nodedata[nodeID] != NULL );
   assert(branchruledata->nodedata[nodeID]->conss != NULL );
   assert(!SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss));

   SCIP_CALL( SCIPallocMemoryArray(scip, &consdatas, SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss)) );
   nconss = 0;

   while(!SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss))
   {
      consdatas[nconss] = (LOGICORDATA*) SCIPqueueRemove(branchruledata->nodedata[nodeID]->conss);
      nconss++;
   }

   /** if the variables of cons2 are a subset the variables of cons1, then cons2 dominates cons1 */
   cons2 = nconss - 1;
   for(cons1 = 0; cons1 < cons2; cons1++)
   {
      if(consdatas[cons1]->nvars > consdatas[cons2]->nvars)
      {
         SCIP_Bool subset;
         int nsubsetvars;
         subset = FALSE;
         nsubsetvars = 0;

#ifdef DEBUG_MODE
         {
            printf("... comparing constraint %d with %d:\n", cons1, cons2);
         }
#endif

         for(varnr2 = 0; varnr2 < (int) consdatas[cons2]->nvars; varnr2++)
         {
            for(varnr1 = 0; varnr1 < (int) consdatas[cons1]->nvars; varnr1++)
            {
               if(strcmp(SCIPvarGetName(consdatas[cons1]->vars[varnr1]), SCIPvarGetName(consdatas[cons2]->vars[varnr2])) == 0)
               {
                  if(consdatas[cons1]->vals[cons1] != consdatas[cons2]->vals[cons2])
                  {
#ifdef DEBUG_MODE
                     {
                        printf("... variable %s %s %f is not equal %s %s %f\n",
                              SCIPvarGetName(consdatas[cons1]->vars[varnr1]),
                              consdatas[cons1]->vals[cons1] == 1 ? ">=" : "<=",
                              consdatas[cons1]->vals[cons1],
                              SCIPvarGetName(consdatas[cons2]->vars[varnr2]),
                              consdatas[cons1]->vals[cons2] == 1 ? ">=" : "<=",
                              consdatas[cons1]->vals[cons2]);
                     }
#endif
                     goto NOTDOMINATED;
                  }
                  else
                  {
                     nsubsetvars++;
                  }
               }
            }
         }

         NOTDOMINATED:

         if (nsubsetvars == (int) consdatas[cons2]->nvars)
            subset = TRUE;

         if (!subset)
         {
            SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->conss, (void* ) (size_t ) consdatas[cons1]) );
         }
         else
         {
            /** free data of dominated constraint */
            SCIPfreeMemoryArray(scip, &consdatas[cons1]->vars);
            SCIPfreeMemoryArray(scip, &consdatas[cons1]->vals);
            SCIPfreeMemory(scip, &consdatas[cons1]);
         }

#ifdef DEBUG_MODE
         {
            if(subset)
               printf("... constraint %d is dominated by %d:\n", cons1, cons2);
         }
#endif
      }
      else
      {
         SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->conss, (void* ) (size_t ) consdatas[cons1]) );
      }
   }
   SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->conss, (void* ) (size_t ) consdatas[cons2]) );

   SCIPfreeMemoryArray(scip, &consdatas);

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
   SCIP_NODE*            node
)
{
   SCIP_Real constant;
   SCIP_Real scalar;
   int nbinvars;
   int varnr;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(branchruledata->saveglbcons);
   assert(SCIPgetCurrentNode(scip) == node);

   /** realloc memory if necessary */
   if(branchruledata->nconsglobal == branchruledata->allocmemglobalcons)
   {
      SCIP_CALL( resizeGlobalsConsData(scip, branchruledata) );
   }

   nbinvars = SCIPgetNOrigBinVars(scip);

   /** ensure that all variables are binary */
   if(nbinvars == SCIPgetNOrigVars(scip))
   {
      /** allocate memory at position nconsglobal and increase the counter */
      SCIP_CALL( SCIPallocMemory(scip, &branchruledata->globalcons[branchruledata->nconsglobal]) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &branchruledata->globalcons[branchruledata->nconsglobal]->vars, SCIPgetOrigVars(scip), nbinvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->globalcons[branchruledata->nconsglobal]->vals, nbinvars) );

      branchruledata->globalcons[branchruledata->nconsglobal]->constype = REOPT_CONSTYPE_SEPA;

      for(varnr = 0; varnr < nbinvars; varnr++)
      {
         scalar = 1;
         constant = 0;

         if(!SCIPvarIsOriginal(branchruledata->globalcons[branchruledata->nconsglobal]->vars[varnr]) )
         {
            SCIP_CALL( SCIPvarGetOrigvarSum( &branchruledata->globalcons[branchruledata->nconsglobal]->vars[varnr], &scalar, &constant) );
         }
         assert(SCIPvarGetStatus(branchruledata->globalcons[branchruledata->nconsglobal]->vars[varnr]) == SCIP_VARSTATUS_ORIGINAL);

         if(SCIPisFeasEQ(scip, abs(SCIPvarGetSol(branchruledata->globalcons[branchruledata->nconsglobal]->vars[varnr], TRUE)), 1))
         {
            branchruledata->globalcons[branchruledata->nconsglobal]->vals[varnr] = 1;
         }
         else
         {
            assert(SCIPisFeasEQ(scip, abs(SCIPvarGetSol(branchruledata->globalcons[branchruledata->nconsglobal]->vars[varnr], TRUE)), 0));
            branchruledata->globalcons[branchruledata->nconsglobal]->vals[varnr] = 0;
         }
      }

      assert(branchruledata->globalcons[branchruledata->nconsglobal]->vars != NULL );
      assert(branchruledata->globalcons[branchruledata->nconsglobal]->vals != NULL );

      branchruledata->nconsglobal++;
   }

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPbranchruleNodereoptAddGlobalCons(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int consnr;
   int varnr;
   int nodenr;
   int nbinvars;

   assert(scip != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL );

   if(!branchruledata->init || !branchruledata->reopt)
      return SCIP_OKAY;

   assert(branchruledata->init);
   assert(branchruledata->reopt);

   nbinvars = SCIPgetNBinVars(scip);

   /** add global constraints if same exists */
   for(consnr = 0; consnr < branchruledata->nconsglobal; consnr++)
   {
      SCIP_CONS* cons;
      assert(branchruledata->globalcons[consnr]->vars != NULL );

      for(varnr = 0; varnr < nbinvars; varnr++)
      {
         if(branchruledata->globalcons[consnr]->vals[varnr] == 1)
         {
            assert(SCIPvarIsOriginal(branchruledata->globalcons[consnr]->vars[varnr]));
            SCIP_CALL( SCIPgetNegatedVar(scip, branchruledata->globalcons[consnr]->vars[varnr], &branchruledata->globalcons[consnr]->vars[varnr]) );
            assert(SCIPvarIsNegated(branchruledata->globalcons[consnr]->vars[varnr]));
         }
      }

      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons,
            (branchruledata->globalcons[consnr]->constype == REOPT_CONSTYPE_SEPA ? "globalsepalogicor" : "globalpseudologicor" ),
            nbinvars, branchruledata->globalcons[consnr]->vars,
            TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE));

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIPfreeMemoryArray(scip, &branchruledata->globalcons[consnr]->vars);
      SCIPfreeMemoryArray(scip, &branchruledata->globalcons[consnr]->vals);
      SCIPfreeMemory(scip, &branchruledata->globalcons[consnr]);
      branchruledata->nconsglobaladded++;
   }

   branchruledata->nconsglobal = 0;

   return SCIP_OKAY;
}

static
SCIP_RETCODE sortConsData(
   SCIP*                 scip,
   SCIP_VAR***           vars,
   SCIP_Real**           vals,
   int                   nvars
)
{
   SCIP_VAR* tmp_var;
   SCIP_Real tmp_val;
   int firstzero;
   int lastone;

   firstzero = 0;
   lastone = nvars - 1;

   while(firstzero < lastone)
   {
      while(*vals[firstzero] == 1)
      {
         firstzero++;
      }

      while(*vals[lastone] == 0)
      {
         lastone--;
      }

      assert(*vals[firstzero] == 0);
      assert(*vals[lastone] == 1);

      tmp_var = (*vars[firstzero]);
      tmp_val = (*vals[firstzero]);

      (*vars[firstzero]) = (*vars[lastone]);
      (*vals[firstzero]) = (*vals[lastone]);

      (*vars[lastone]) = tmp_var;
      (*vals[lastone]) = tmp_val;
   }

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

   /** allocate memory if necessary */
   if(saveinnodedata && branchruledata->nodedata[nodeID]->allocmem == 0)
   {
      assert(nodeID > 0);
      assert(branchruledata->nodedata[nodeID]->vars == NULL );
      assert(branchruledata->nodedata[nodeID]->varbounds == NULL );
      assert(branchruledata->nodedata[nodeID]->varboundtypes == NULL );

      branchruledata->nodedata[nodeID]->allocmem = SCIPgetNVars(scip);
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeID]->vars, consdata->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeID]->varbounds, consdata->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &branchruledata->nodedata[nodeID]->varboundtypes, consdata->nvars) );
   }

   for(varnr = 0; varnr < consdata->nvars; varnr++)
   {
      var = consdata->vars[varnr];
      val = consdata->vals[varnr];
      boundtype = val == 1 ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;
      SCIP_CALL(SCIPvarGetProbvarBound(&var, &val, &boundtype));
      assert(SCIPvarIsTransformedOrigvar(var));

      if(negate)
      {
         val = 1 - val;
         boundtype = (SCIP_BOUNDTYPE) (SCIP_BOUNDTYPE_UPPER - boundtype);
      }

      if(boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisGT(scip, val, SCIPvarGetLbLocal(var))
          && SCIPisFeasLE(scip, val, SCIPvarGetUbLocal(var)))
      {
         SCIP_CALL(SCIPchgVarLbNode(scip, node, var, val));
      }
      else if(boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisLT(scip, val, SCIPvarGetUbLocal(var))
               && SCIPisFeasGE(scip, val, SCIPvarGetLbLocal(var)))
      {
         SCIP_CALL(SCIPchgVarUbNode(scip, node, var, val));
      }
      else if(boundtype != SCIP_BOUNDTYPE_LOWER && boundtype != SCIP_BOUNDTYPE_UPPER)
      {
         printf("** Unknown boundtype: %d **\n", boundtype);
         assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);
      }
#ifdef DEBUG_MODE
      {
         printf("\tfix variable %s %s %f in node %llu\n", SCIPvarGetName(var),
               boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", val, SCIPnodeGetNumber(node));
      }
#endif

      /** add variable and bound to branching path information, because we don't want to delete this data */
      if(saveinnodedata)
      {
         SCIP_Real constant;
         SCIP_Real scalar;

         constant = 0;
         scalar = 1;
         varpos = branchruledata->nodedata[nodeID]->nvars;

         branchruledata->nodedata[nodeID]->vars[varpos] = var;

         if (!SCIPvarIsOriginal(branchruledata->nodedata[nodeID]->vars[varpos]))
         {
            if (val == 1)
            {
               SCIP_CALL( SCIPvarGetOrigvarSum(&branchruledata->nodedata[nodeID]->vars[varpos], &scalar, &constant) );
               branchruledata->nodedata[nodeID]->varbounds[varpos] = 1;
               branchruledata->nodedata[nodeID]->varboundtypes[varpos] = SCIP_BOUNDTYPE_LOWER;
            }
            else
            {
               SCIP_CALL( SCIPvarGetOrigvarSum(&branchruledata->nodedata[nodeID]->vars[varpos], &scalar, &constant) );
               branchruledata->nodedata[nodeID]->varbounds[varpos] = 0;
               branchruledata->nodedata[nodeID]->varboundtypes[varpos] = SCIP_BOUNDTYPE_UPPER;
            }
         }
         assert(SCIPvarIsOriginal(branchruledata->nodedata[nodeID]->vars[varpos]));
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
   SCIP_Bool             checkredundancy
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
      int varnr;

      /** copy variables and negate them if necessary */
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &vars, consdata->vars, consdata->nvars) );
      for(varnr = 0; varnr < consdata->nvars; varnr++)
      {
         SCIP_Real val;
         SCIP_BOUNDTYPE boundtype;

         assert(SCIPvarIsOriginal(vars[varnr]));

         val = consdata->vals[varnr];
         boundtype = val == 1 ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;
         SCIP_CALL(SCIPvarGetProbvarBound(&vars[varnr], &val, &boundtype));
         assert(SCIPvarIsTransformed(vars[varnr]));

         if(val == 1)
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, vars[varnr], &vars[varnr]) );
            assert(SCIPvarIsNegated(vars[varnr]));
         }
      }

      /** add old constraints to child1 */
      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons,
            (consdata->constype == REOPT_CONSTYPE_SEPA ? "localsepalogicor" : "localpseudologicor" ),
            consdata->nvars, vars,
            TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE));

#ifdef DEBUG_MODE
      {
         printf("\tadded constraint in node %llu:\n\t>> ", SCIPnodeGetNumber(node));
         SCIP_CALL(SCIPprintCons(scip, cons, NULL));
      }
#endif
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIPfreeMemoryArray(scip, &vars);
   }
   else
   {
      SCIP_CALL( fixVarsInNode(scip, branchruledata, node, consdata, FALSE, TRUE, 0, checkredundancy) );
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
      (*nbndchgs) += SCIPnodeGetNDomchg(*parent);
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

   assert(branchruledata->nodedata[(*parentID)]->oldnumber == SCIPnodeGetNumber((*parent)));

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

      SCIPfreeMemory(scip, &branchruledata->nodedata[nodeID]);
      branchruledata->nodedata[nodeID] = NULL;
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
      branchruledata->nodedata[nodeID]->oldnumber = nodeID == 0 ? 1 : 0;
      branchruledata->nodedata[nodeID]->pseudobranched = FALSE;
      branchruledata->nodedata[nodeID]->redies = 0;
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
   SCIP_BRANCHRULEDATA*  branchruledata
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
         SCIP_CALL( deleteNodeData(scip, branchruledata, nodeID, FALSE) );
         assert(branchruledata->nodedata[nodeID]->nvars == 0);
         assert(branchruledata->nodedata[nodeID]->allocmem <= SCIPgetNVars(scip));
      }

      if(nodeID > 0)
      {
         SCIP_CALL( SCIPqueueInsert(branchruledata->openIDs, (void* ) (size_t ) nodeID) );
      }
   }
   assert(SCIPqueueNElems(branchruledata->openIDs) == branchruledata->allocmemsizenodedata - 1);

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
   assert(branchruledata->nodedata[nodeID]->vars != NULL );
   assert(branchruledata->nodedata[nodeID]->varbounds != NULL );
   assert(branchruledata->nodedata[nodeID]->varboundtypes != NULL );
   assert(child1 != NULL );

   /** copy memory to ensure that only original variables are saved */
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &vars, branchruledata->nodedata[nodeID]->vars, branchruledata->nodedata[nodeID]->nvars) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &vals, branchruledata->nodedata[nodeID]->varbounds, branchruledata->nodedata[nodeID]->nvars) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &boundtypes, branchruledata->nodedata[nodeID]->varboundtypes, branchruledata->nodedata[nodeID]->nvars) );

   /**
    * Chance bounds due to all ancestors.
    * Fix all decisions from branching path child1 and child2
    */
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

#ifdef DEBUG_MODE
      {
         printf("\tvar %s %s %f\n", SCIPvarGetName(vars[varnr]),
               boundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=", newbound);
      }
#endif

      if(boundtypes[varnr] == SCIP_BOUNDTYPE_LOWER && SCIPisGT(scip, newbound, oldlb)
          && SCIPisFeasLE(scip, newbound, oldub))
      {
         SCIP_CALL( SCIPchgVarLbNode(scip, child1, vars[varnr], newbound) );

         if(child2 != NULL )
         {
            SCIP_CALL( SCIPchgVarLbNode(scip, child2, vars[varnr], newbound) );
         }
      }
      else if(boundtypes[varnr] == SCIP_BOUNDTYPE_UPPER && SCIPisLT(scip, newbound, oldub)
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
   }

   /** free memory */
   SCIPfreeMemoryArray(scip, &vars);
   SCIPfreeMemoryArray(scip, &vals);
   SCIPfreeMemoryArray(scip, &boundtypes);

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
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID,
   int                   destnodeID
)
{
   assert(branchruledata != NULL);
   assert(nodeID >= 1);
   assert(destnodeID >= 0);
   assert(branchruledata->nodedata[nodeID]->nodechilds != NULL);

   while(!SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds))
   {
      SCIP_VAR** vars;
      int childID;
      int varnr;

      childID = (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[nodeID]->nodechilds);
      SCIP_CALL(SCIPduplicateMemoryArray(scip, &vars, branchruledata->nodedata[nodeID]->vars, branchruledata->nodedata[nodeID]->nvars));

      assert(branchruledata->nodedata[childID]->allocmem >= branchruledata->nodedata[childID]->nvars + branchruledata->nodedata[nodeID]->nvars);

      /** save branching information */
      for(varnr = 0; varnr < branchruledata->nodedata[nodeID]->nvars; varnr++)
      {
         branchruledata->nodedata[childID]->vars[branchruledata->nodedata[childID]->nvars] = vars[varnr];
         branchruledata->nodedata[childID]->varbounds[branchruledata->nodedata[childID]->nvars] = branchruledata->nodedata[nodeID]->varbounds[varnr];
         branchruledata->nodedata[childID]->varboundtypes[branchruledata->nodedata[childID]->nvars] = branchruledata->nodedata[nodeID]->varboundtypes[varnr];
         branchruledata->nodedata[childID]->nvars++;
      }

      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[destnodeID]->nodechilds, (void* ) (size_t ) childID));
      SCIPfreeMemoryArray(scip, &vars);
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
   SCIP_NODE*            child2
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
   while(curcons <= ncons)
   {
      LOGICORDATA* consdata;
      SCIP_CONS* cons1;
      SCIP_CONS* cons2;
      SCIP_VAR** vars;

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
         boundtype = val == 1 ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;
         SCIP_CALL( SCIPvarGetProbvarBound(&vars[varnr], &val, &boundtype) );
         assert(SCIPvarIsTransformed(vars[varnr]));

         if (val == 1)
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, vars[varnr], &vars[varnr]) );
            assert(SCIPvarIsNegated(vars[varnr]));
         }
      }

      /** add old constraints to child1 */
      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons1,
            (consdata->constype == REOPT_CONSTYPE_SEPA ? "localsepalogicor" : "localpseudologicor" ),
            consdata->nvars, vars,
            TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE));

#ifdef DEBUG_MODE
      {
         printf("\tadd constraint to node %lld from some round before.\n\t>> ",
               SCIPnodeGetNumber(child1));
         SCIP_CALL(SCIPprintCons(scip, cons1, NULL));
      }
#endif

      SCIP_CALL( SCIPaddConsNode(scip, child1, cons1, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons1) );

      /** add old constraints to child2, if child was PB */
      if(child2 != NULL )
      {
         SCIP_CALL( SCIPcreateConsLogicor(scip, &cons2,
               (consdata->constype == REOPT_CONSTYPE_SEPA ? "localsepalogicor" : "localpseudologicor" ),
               consdata->nvars, vars,
               TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE));

#ifdef DEBUG_MODE
         {
            printf("\n\tadd constraint to node %lld from some round before.\n\t>> ", SCIPnodeGetNumber(child2));
            SCIP_CALL( SCIPprintCons(scip, cons2, NULL) );
         }
#endif

         SCIP_CALL( SCIPaddConsNode(scip, child2, cons2, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
      }

      SCIPfreeMemoryArray(scip, &vars);

      curcons++;
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->conss, (void* ) (size_t ) consdata) );
   }

   return SCIP_OKAY;
}

/*
 * Execute the branching of nodes with additional constraints.
 */
static
SCIP_RETCODE Exec(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata
)
{
   int nodeID;
   int curChild;
   int nChilds;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

#ifdef DEBUG_MODE
   {
      printf("** BRANCH NODE %llu **\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
   }
#endif

   /**
    * current node is equal to the root and the root was pseudo-branched
    * we have to create two child nodes; one with the pseudo-constraint and
    * one with the negated fixings.
    */
   if(SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) && branchruledata->nodedata[0]->pseudobranched)
   {
      SCIP_NODE* child1;
      SCIP_NODE* child2;
      LOGICORDATA* consdata;
      int child1_ID;
      int nvars;

      child1 = NULL;
      child2 = NULL;

      SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
      SCIP_CALL( SCIPcreateChild(scip, &child2, 0.0, SCIPgetLocalTransEstimate(scip)) );
      assert(child1 != NULL );
      assert(child2 != NULL );

      /** get a local logic-or constraint from brach_pseudo.c */
      nvars = SCIPbranchrulePseudoGetNPseudoVars(scip, 0);
      SCIP_CALL( SCIPallocMemory(scip, &consdata) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->vars, nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->vals, nvars) );
      consdata->nvars = -1;
      consdata->constype = REOPT_CONSTYPE_PSEUDO;
      SCIP_CALL( SCIPbranchrulePseudoGenerateCons(scip, consdata, &(consdata->nvars), nvars, 0, FALSE) );
      assert(consdata->nvars > 0);
      assert(consdata->vars != NULL );
      assert(consdata->vals != NULL );

      /***********************************************************************
       *
       * create copy of root and add all decisions from strong branching
       *
       * 1. get a free slot to save all information, initialize the data
       * 2. fix all decisions from strong branching
       * 3. move all children of root node to child1
       * 4. add child1 to the list of current leafs from re-optimization
       *
       ***********************************************************************/
#ifdef DEBUG_MODE
      {
         printf("NEW NODE #%llu:\n", SCIPnodeGetNumber(child1));
      }
#endif

      /** add this node (the branching path is empty, thats way this method will on
       * allocate memory and create a link between node number and ID) */
      child1_ID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
      SCIP_CALL( initNode(scip, branchruledata, child1_ID) );
      SCIP_CALL( SCIPnodeSetReoptID(child1, child1_ID) );

      /* tag node as LEAF */
      SCIP_CALL( SCIPnodeSetReopttype(child1, SCIP_REOPTTYPE_LEAF) );

      branchruledata->nodedata[child1_ID]->oldnumber = SCIPnodeGetNumber(child1);

      /** move all child nodes (w/o the current) of the pseudo-branched root node to child1 with child1_ID */
      assert(branchruledata->nodedata[child1_ID]->nodechilds == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[child1_ID]->nodechilds));

      if(branchruledata->nodedata[child1_ID]->nodechilds == NULL )
      {
         SCIP_CALL( SCIPqueueCreate(&branchruledata->nodedata[child1_ID]->nodechilds, SCIPqueueNElems(branchruledata->nodedata[0]->nodechilds), 2) );
      }

      while(!SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds))
      {
         SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[child1_ID]->nodechilds, SCIPqueueRemove(branchruledata->nodedata[0]->nodechilds)) );
      }
      assert(SCIPqueueIsEmpty(branchruledata->nodedata[0]->nodechilds));

      /* insert child1_ID as a root child, we have to do this because this node will be deleted after branching
       * by RemoveNode() */
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[0]->nodechilds, (void* ) (size_t ) child1_ID) );

      /** fix the corresponding variables from 'cons' in child1 */
      SCIP_CALL( fixVarsInNode(scip, branchruledata, child1, consdata, TRUE, FALSE, child1_ID, FALSE) );

      /***********************************************************************
       *
       * create copy of root and add all decisions from strong branching
       *
       * 1. add and release the logic-or constraint
       * 2. if the added constraint includes more than one variable, save the
       *    node
       *
       ***********************************************************************/
#ifdef DEBUG_MODE
      {
         printf("NEW NODE #%llu:\n", SCIPnodeGetNumber(child2));
      }
#endif

      /* add the logic-or constraint 'cons' to child2 */
      SCIP_CALL( addLocalConsToNode(scip, branchruledata, child2, consdata, FALSE) );

      /** free consdata */
      SCIPfreeMemoryArray(scip, &consdata->vars);
      SCIPfreeMemoryArray(scip, &consdata->vals);
      SCIPfreeMemory(scip, &consdata);

      /** check if child2 includes some added constraints and save this node*/
      if(SCIPnodeGetNAddedcons(scip, child2) > 0)
      {
         SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, child2, SCIP_REOPTTYPE_LEAF) );

      }

      /** remove pseudobranched flag */
      branchruledata->nodedata[0]->pseudobranched = FALSE;
      branchruledata->nsavednodes--;

      return SCIP_OKAY;
   }
   else if(SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip)
         && !branchruledata->nodedata[0]->pseudobranched)
   {
      nodeID = 0;
      goto REVIVE;
   }

   /**
    * we have to generate all nodes in currentleafIDs and one node extra with the fixed bounds if the node was pseudo-branched.
    * all nodes with a non empty queue in 'nodechilds' will be add to currentleafIDs.
    */
   nodeID = SCIPnodeGetReoptID(SCIPgetCurrentNode(scip));
   assert(nodeID >= 1);

   REVIVE:

   assert(branchruledata->nodedata[nodeID]->nodechilds != NULL );

   /** iterate over all children */
   nChilds = SCIPqueueNElems(branchruledata->nodedata[nodeID]->nodechilds);
   curChild = 1;
   while (curChild <= nChilds)
   {
      SCIP_NODE* child1;
      SCIP_NODE* child2;
      int childID;

      child1 = NULL;
      child2 = NULL;

      childID = (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[nodeID]->nodechilds);
      assert(childID >= 1); /** childID has to be greater or equal than 1, because the ID 0 is
       reserved for the root and the root could never be a child */
      assert(childID < branchruledata->allocmemsizenodedata);

#ifdef DEBUG_MODE
      {
         printf("CHILDID: %u\n", childID);
      }
#endif

      /** if all children were infeasible in the LP, we can remove the 'pseudobracnhed' flag
       * and have to generate only the node itself.
       */
      if(branchruledata->nodedata[childID]->pseudobranched
      && (branchruledata->nodedata[childID]->nodechilds == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[childID]->nodechilds)))
      {
         SCIP_CALL( SCIPbranchrulePseudoDelInformation(scip, childID) );
         branchruledata->nodedata[childID]->pseudobranched = FALSE;
      }

      /** the node at position childID was pseudobranched, we have to split this node */
      if (branchruledata->nodedata[childID]->pseudobranched)
      {
         LOGICORDATA* consdata;
         int nvars;

         assert(branchruledata->nodedata[childID]->reopttype == SCIP_REOPTTYPE_STRBRANCHED);

         /** generate two children */
         SCIP_CALL(SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)));
         SCIP_CALL(SCIPcreateChild(scip, &child2, 0.0, SCIPgetLocalTransEstimate(scip)));
         assert(child1 != NULL );
         assert(child2 != NULL );

#ifdef DEBUG_MODE
         {
            if (branchruledata->nodedata[childID]->nodechilds != NULL && !SCIPqueueIsEmpty(
                      branchruledata->nodedata[childID]->nodechilds))
            {
               printf("NEW NODEs #%llu and #%llu from #%llu:\n", SCIPnodeGetNumber(child1),
                     SCIPnodeGetNumber(child2), SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
            }
            else
            {
               printf("NEW NODE #%llu from #%llu:\n", SCIPnodeGetNumber(child2),
                     SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
            }
         }
#endif

         /** change the bound in child1 and child2 along the branching path at position childID*/
         SCIP_CALL(changeAncestorBranchings(scip, branchruledata, childID, child1, child2));

         /** get logic-or constraint from branch_pseudo.c */
         nvars = SCIPbranchrulePseudoGetNPseudoVars(scip, childID);
         SCIP_CALL( SCIPallocMemory(scip, &consdata) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->vars, nvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->vals, nvars) );
         SCIP_CALL( SCIPbranchrulePseudoGenerateCons(scip, consdata, &(consdata->nvars), nvars, childID, FALSE) );
         assert(consdata->nvars > 0);
         assert(consdata->vars != NULL );
         assert(consdata->vals != NULL );

         /** fix all decisions from strong branching in child1 */
         SCIP_CALL( fixVarsInNode(scip, branchruledata, child1, consdata, TRUE, FALSE, childID, branchruledata->savelpbasis) );

         /** update node number at position childID */
         branchruledata->nodedata[childID]->oldnumber = SCIPnodeGetNumber(child1);

         /** add 'pseudocons' to child2; the added constraint can not be redundant */
         SCIP_CALL( addLocalConsToNode(scip, branchruledata, child2, consdata, FALSE) );

         if (consdata->nvars == 1)
            branchruledata->nodedata[childID]->pseudobranched = FALSE;

         /** free consdata */
         SCIPfreeMemoryArray(scip, &consdata->vars);
         SCIPfreeMemoryArray(scip, &consdata->vals);
         SCIPfreeMemory(scip, &consdata);
      }
      else
      {
         /**
          * node at position childID was not pseudobranched, we generate only the node itself
          */

         SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
         assert(child1 != NULL );
         assert(child2 == NULL );

#ifdef DEBUG_MODE
         {
            printf("NEW NODE #%llu from #%llu:\n", SCIPnodeGetNumber(child1),
                  SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
         }
#endif

         /** Chance bounds due to all ancestors */
         SCIP_CALL( changeAncestorBranchings(scip, branchruledata, childID, child1, NULL) );

         /** update node number at position childID */
         branchruledata->nodedata[childID]->oldnumber = SCIPnodeGetNumber(child1);
      }

      /** set LPI state */
      if(branchruledata->savelpbasis && branchruledata->nodedata[childID]->lpistate != NULL )
      {
         SCIP_CALL( SCIPchildSetLpistate(child1, branchruledata->nodedata[childID]->lpistate) );
      }

      /** add local constraint from an iteration before (if some exists) to child1 */
      if(branchruledata->nodedata[childID]->conss != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[childID]->conss) )
      {
         int nconss;

         nconss = SCIPqueueNElems(branchruledata->nodedata[childID]->conss);

         /** generate all local constraints and add them to child1 and child2 (if exists) */
         SCIP_CALL(genLocalCons(scip, branchruledata, childID, child1, child2));

         assert(SCIPqueueNElems(branchruledata->nodedata[childID]->conss) == nconss);
         assert(child1 == NULL || (child1 != NULL && SCIPnodeGetNAddedcons(scip, child1) == nconss));
         assert(child2 == NULL || (child2 != NULL && SCIPnodeGetNAddedcons(scip, child2) == (branchruledata->nodedata[childID]->pseudobranched ? nconss + 1 : nconss)));
      }

      /** remove flag 'pseudobranched' */
      branchruledata->nodedata[childID]->pseudobranched = FALSE;

      /** insert child1 to the child list of nodeID */
      SCIP_CALL( SCIPqueueInsert(branchruledata->nodedata[nodeID]->nodechilds, (void* ) (size_t ) childID) );
      SCIP_CALL( SCIPnodeSetReoptID(child1, childID) );

      /* set the REOPTTYPE */
      assert(branchruledata->nodedata[childID]->reopttype >= SCIP_REOPTTYPE_TRANSIT);
      SCIPnodeSetReopttype(child1, branchruledata->nodedata[childID]->reopttype);

      /** check if child2 includes some added constraints, save this node and set REOPT to true */
      if(child2 != NULL && SCIPnodeGetNAddedcons(scip, child2) > 0)
      {
         assert(branchruledata->nodedata[childID] != NULL);
         SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, child2, SCIP_REOPTTYPE_LEAF) );
      }

      curChild++;
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE deleteChildrenBelow(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID
)
{
   assert(branchruledata != NULL );
   assert(nodeID >= 1);
   assert(branchruledata->init && branchruledata->reopt);
   assert(branchruledata->nodedata[nodeID] != NULL );

   /** delete all children below */
   if (branchruledata->nodedata[nodeID]->nodechilds != NULL )
   {
      while (!SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds))
      {
         SCIP_CALL( deleteChildrenBelow(scip, branchruledata, (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[nodeID]->nodechilds)) );
      }
   }

   /** delete node data*/
   SCIP_CALL(deleteNodeData(scip, branchruledata, nodeID, TRUE));

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
   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   /* check if maximal number if saved nodes is reached */
   if (branchruledata->nsavednodes > branchruledata->maxsavednodes)
      branchruledata->restart = TRUE;

   return;
}

/*
 * non-static methods
 */
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
   if (branchruledata->nodedata[0] != NULL && branchruledata->nodedata[0]->lpistate != NULL )
   {
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
   SCIP_REOPTTYPE        reopttype
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_NODE* parent;
   int nodeID;
   int parentID;

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL );
   assert(branchruledata->reopt);
   assert(reopttype == SCIP_REOPTTYPE_TRANSIT
       || reopttype == SCIP_REOPTTYPE_LOGICORNODE
       || reopttype == SCIP_REOPTTYPE_STRBRANCHED
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

#ifdef DEBUG_MODE
      {
         int varnr;
         printf(
               "UPDATE NODE #%lld AT ID %u WITH %u BRANCHINGS AND %u ADDED CONS, REOPTTYPE %d\n",
               SCIPnodeGetNumber(node), nodeID, branchruledata->nodedata[nodeID]->nvars,
               SCIPnodeGetNAddedcons(scip, node), reopttype);

         for (varnr = 0; varnr < branchruledata->nodedata[nodeID]->nvars; varnr++)
         {
            printf("\t%s %s %f\n", SCIPvarGetName(branchruledata->nodedata[nodeID]->vars[varnr]),
                  branchruledata->nodedata[nodeID]->varboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                        "=>" : "<=", branchruledata->nodedata[nodeID]->varbounds[varnr]);
         }
      }
#endif
      /** update the saved node number */
      assert(nodeID < branchruledata->allocmemsizenodedata);
      branchruledata->nodedata[nodeID]->oldnumber = SCIPnodeGetNumber(node);

      /** update LPI state if node is pseudobranched or feasible */
      switch (reopttype) {
         case SCIP_REOPTTYPE_TRANSIT:
            assert(branchruledata->nodedata[nodeID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss));

            if( branchruledata->nodedata[nodeID]->nodechilds != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds) )
            {
               SCIP_NODE* parent;
               int parentID;
               int nbndchgdiff;

               SCIP_CALL( getLastSavedNode(scip, branchruledata, node, &parent, &parentID, &nbndchgdiff) );
               assert(parentID != nodeID);
               assert(branchruledata->nodedata[parentID] != NULL );
               assert(branchruledata->nodedata[parentID]->nodechilds != NULL && !SCIPqueueIsEmpty(branchruledata->nodedata[parentID]->nodechilds));

               /* check if we can move all children to the next saved node above */
               if( nbndchgdiff <= branchruledata->maxdiffofnodes )
               {
                  SCIP_CALL( moveChildrenUp(branchruledata, nodeID, parentID) );

                  /* delete this node */
                  SCIP_CALL( SCIPbranchruleNodereoptRemoveNode(scip, node, TRUE, FALSE) );

                  /** stop clock */
                  SCIP_CALL(SCIPstopSaveTime(scip));

                  return SCIP_OKAY;
               }
            }
            goto TRANSIT;

            break;

         case SCIP_REOPTTYPE_LOGICORNODE:
         case SCIP_REOPTTYPE_LEAF:
            goto TRANSIT;
            break;

         case SCIP_REOPTTYPE_STRBRANCHED:
            if (branchruledata->savelpbasis)
            {
               SCIP_CALL(saveLPIstate(scip, branchruledata, node, nodeID));
            }

            goto PSEUDO;
            break;

         case SCIP_REOPTTYPE_FEASIBLE:
            if (branchruledata->savelpbasis)
            {
               SCIP_CALL(saveLPIstate(scip, branchruledata, node, nodeID));
            }
            goto FEASIBLE;
            break;

         case SCIP_REOPTTYPE_PRUNED:
            branchruledata->nrediednodesround++;
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

   /** the current node is the root node and it can be only
    * TRANSIT, BRANCHED or FEASIBLE */
   if (SCIPnodeGetDepth(node) == 0)
   {
      nodeID = 0;

      switch (reopttype) {
         case SCIP_REOPTTYPE_TRANSIT:
            branchruledata->nodedata[nodeID]->reopttype = reopttype;
            goto TRANSIT;
            break;

         case SCIP_REOPTTYPE_STRBRANCHED:
            branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_STRBRANCHED;
            branchruledata->nodedata[nodeID]->pseudobranched = TRUE;
            branchruledata->nodedata[nodeID]->nvars = 0;
#ifdef DEBUG_MODE
            {
               printf("SAVE NODE #%llu AT ID %d, reopttype %d\n",
               SCIPnodeGetNumber(node), nodeID, reopttype);
            }
#endif
            goto PSEUDO;
            break;

         case SCIP_REOPTTYPE_FEASIBLE:
            branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_FEASIBLE;
            if( branchruledata->saveglbcons )
            {
               SCIP_CALL(saveGlobalCons(scip, branchruledata, node));
            }
            goto FEASIBLE;
            break;

         default:
            assert(reopttype == SCIP_REOPTTYPE_TRANSIT
                || reopttype == SCIP_REOPTTYPE_STRBRANCHED
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

      if( reopttype < SCIP_REOPTTYPE_LOGICORNODE
        && nbndchgdiff <= branchruledata->maxdiffofnodes )
      {
         /** stop clock */
         SCIP_CALL(SCIPstopSaveTime(scip));

         return SCIP_OKAY;
      }

      /** get an empty slot (to ensure that ID != NULL, the queue should contains only IDs >= 1);
       * the 0 is always reserved for the root */
      nodeID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
      assert(nodeID >= 1 && nodeID < branchruledata->allocmemsizenodedata);
      assert(SCIPgetRootNode(scip) != node);

      /** get memory for nodedata */
      assert(branchruledata->nodedata[nodeID] == NULL || branchruledata->nodedata[nodeID]->nvars == 0);
      SCIP_CALL(initNode(scip, branchruledata, nodeID));
      branchruledata->nodedata[nodeID]->oldnumber = SCIPnodeGetNumber(node);

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
      SCIP_CALL(saveAncestorBranchings(scip, branchruledata, node, parent, nodeID, parentID));

      /** save pseudo-constraints (if one exists) */
      if (SCIPnodeGetNAddedcons(scip, node) >= 1)
      {
         assert(branchruledata->nodedata[nodeID]->conss == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->conss));
         if (branchruledata->nodedata[nodeID]->conss == NULL )
         {
            SCIP_CALL(SCIPqueueCreate(&branchruledata->nodedata[nodeID]->conss, SCIPnodeGetNAddedcons(scip, node), 2));
         }
         SCIP_CALL(saveLocalConssData(scip, branchruledata, node, nodeID));

//         if( SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss) >= 2 )
//         {
//            SCIP_CALL( findMinimalConss(scip, nodeID) );
//         }
      }

      /** save the basis if the node is feasible or branched by strong branching,
       in other cases the LP is possibly not solved */
      if (branchruledata->savelpbasis
       && (reopttype == SCIP_REOPTTYPE_STRBRANCHED || reopttype == SCIP_REOPTTYPE_FEASIBLE))
      {
         SCIP_CALL(saveLPIstate(scip, branchruledata, node, nodeID));
      }

      /* set ID */
      SCIP_CALL(SCIPnodeSetReoptID(node, nodeID));

      /* set the REOPTTYPE */
      SCIPnodeSetReopttype(node, reopttype);

      branchruledata->nsavednodes++;

#ifdef DEBUG_MODE
      {
         int varnr;
         printf(
               "SAVE NODE #%llu AT ID %d WITH %d BRANCHINGS AND %d ADDED CONS, REOPTTYPE %d \n",
               SCIPnodeGetNumber(node), nodeID, branchruledata->nodedata[nodeID]->nvars,
               branchruledata->nodedata[nodeID]->conss == NULL ?
                     0 : SCIPqueueNElems(branchruledata->nodedata[nodeID]->conss), reopttype);
         for (varnr = 0; varnr < branchruledata->nodedata[nodeID]->nvars; varnr++)
         {
            printf("\t%s %s %f\n", SCIPvarGetName(branchruledata->nodedata[nodeID]->vars[varnr]),
                  branchruledata->nodedata[nodeID]->varboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                        "=>" : "<=", branchruledata->nodedata[nodeID]->varbounds[varnr]);
         }
      }
#endif
   }

   switch (reopttype) {
      case SCIP_REOPTTYPE_TRANSIT:
      case SCIP_REOPTTYPE_LOGICORNODE:
      case SCIP_REOPTTYPE_LEAF:
         TRANSIT:
         branchruledata->nodedata[nodeID]->reopttype = reopttype;
         break;

      case SCIP_REOPTTYPE_STRBRANCHED:
         PSEUDO:
         branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_STRBRANCHED;
         branchruledata->nodedata[nodeID]->pseudobranched = TRUE;
         branchruledata->npbnodesround++;

         if (nodeID > 0)
         {
            assert(node != SCIPgetRootNode(scip));
            SCIP_CALL(SCIPbranchrulePseudoLinkIDs(scip, nodeID));
         }
         break;

      case SCIP_REOPTTYPE_FEASIBLE:
         FEASIBLE:
         branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_FEASIBLE;
         branchruledata->nodedata[nodeID]->pseudobranched = FALSE;
         branchruledata->nfeasnodesround++;

         if (branchruledata->saveglbcons && SCIPgetRootNode(scip) == node)
         {
            SCIP_CALL(saveGlobalCons(scip, branchruledata, node));
         }

         /**
          * save all information of the current feasible solution to separate this
          * solution in a following round (but only if all variablea are binary)
          * TODO: Verbesserungswürdig
          */
   //      if( SCIPgetNVars(scip) == SCIPgetNBinVars(scip) )
   //      {
   //         SCIP_CALL( saveLocalConsData(scip, node, nodeID) );
   //         SCIP_CALL( SCIPbranchruleNodereoptSaveGlobalCons(scip, node) );
   //      }
         break;

      case SCIP_REOPTTYPE_PRUNED:
         PRUNED:
         branchruledata->nodedata[nodeID]->reopttype = SCIP_REOPTTYPE_PRUNED;
         branchruledata->nodedata[nodeID]->pseudobranched = FALSE;
         branchruledata->nprunednodesround++;
         break;

      default:
         assert(reopttype == SCIP_REOPTTYPE_TRANSIT
             || reopttype == SCIP_REOPTTYPE_LOGICORNODE
             || reopttype == SCIP_REOPTTYPE_LEAF
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

   /** delete all children below if the node is infeasible */
   if (infeasible && branchruledata->nodedata[nodeID]->nodechilds != NULL )
   {
      while (!SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds))
      {
         SCIP_CALL( deleteChildrenBelow(scip, branchruledata, (int) (size_t) SCIPqueueRemove(branchruledata->nodedata[nodeID]->nodechilds)) );
      }
   }

   /* delete node data */
   assert(branchruledata->nodedata[nodeID]->nodechilds == NULL || SCIPqueueIsEmpty(branchruledata->nodedata[nodeID]->nodechilds));
   SCIP_CALL(deleteNodeData(scip, branchruledata, nodeID, TRUE));
   SCIP_CALL(SCIPnodeSetReoptID(node, -1));
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
   SCIP_NODE*            node,
   SCIP_EVENT*           event
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

//   SCIP_CALL( SCIPbranchruleNodereoptSaveGlobalCons(scip, node) );

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
   nodeID = SCIPnodeGetReoptID(node);

   (*solvelp) = TRUE;

   /* the root LP should be always solved */
   if( nodeID == 0 )
      return SCIP_OKAY;
   else
      switch (branchruledata->solvelp) {
         /* solve all LPs */
         case 0:
            if( SCIPnodeGetReopttype(node) < SCIP_REOPTTYPE_LEAF )
            {
               if( SCIPnodeGetDepth(node) % branchruledata->solvelpdiff != 0
                && branchruledata->nodedata[nodeID]->nvars < branchruledata->solvelpdiff)
               (*solvelp) = FALSE;
            }
            break;

         default:
            if( (int) SCIPnodeGetReopttype(node) < branchruledata->solvelp )
               (*solvelp) = FALSE;
            break;
      }

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPbranchruleNodereoptRestartCheck(
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

   if (branchruledata->init)
   {
      /* check if a restart should be forced */
      restartCheck(scip, branchruledata);

      if (branchruledata->restart)
      {
         /* clear nodedate */
         SCIP_CALL(clearNodes(scip, branchruledata));

         /* clear saved information in branch_pseudo */
         SCIP_CALL(SCIPbranchrulePseudoReset(scip, TRUE));

         /* reset flag */
         branchruledata->restart = FALSE;

         /* decrease number of restarts */
         branchruledata->nrestarts++;

         /* set the number of saved nodes to 0 */
         branchruledata->nsavednodes = 0;

         printf(">> restart reoptimization ...\n");
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
   int*                  nrestartsround
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
      (*nfeasnodes) = branchruledata->nfeasnodes;

   if( nfeasnodesround != NULL)
      (*nfeasnodesround) = branchruledata->nfeasnodesround;

   /* infeasible nodes */
   if( ninfeasnodes != NULL)
      (*ninfeasnodes) = branchruledata->ninfeasnodes;

   if( ninfeasnodesround != NULL)
      (*ninfeasnodesround) = branchruledata->ninfeasnodesround;

   /* pruned nodes */
   if( nprunednodes != NULL)
      (*nprunednodes) = branchruledata->nprunednodes;

   if( nprunednodesround != NULL)
      (*nprunednodesround) = branchruledata->nprunednodesround;

   /* redied nodes */
   if( nrediednodes != NULL)
      (*nrediednodes) = branchruledata->nrediednodes;

   if( nrediednodesround != NULL)
      (*nrediednodesround) = branchruledata->nrediednodesround;

   /* other stuff */
   if( nruns != NULL)
      (*nruns) = branchruledata->nruns;

   if( nrestarts != NULL)
      (*nrestarts) = branchruledata->nrestarts;

   if( nrestartsround != NULL)
      (*nrestartsround) = branchruledata->nrestartsround;

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

/*
 * Callback methods of branching rule
 */

#define branchCopynodereopt NULL;
#define branchExitnodereopt NULL;
#define branchInitsolnodereopt NULL;

static
SCIP_DECL_BRANCHINIT(branchInitnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL );
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL );

   /** HACK */
   /** check if all variable are binary, if not, disable reoptimization */
   if (!branchruledata->init)
   {
      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/enable", &branchruledata->reopt));
      if (branchruledata->reopt && SCIPgetNVars(scip) - SCIPgetNBinVars(scip) > 0)
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
      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/dynamicdiffofnodes", &branchruledata->dynamicdiffofnodes));

      if( !branchruledata->dynamicdiffofnodes )
      {
         SCIP_CALL(SCIPgetIntParam(scip, "reoptimization/maxdiffofnodes", &branchruledata->maxdiffofnodes));
      }
      else
         branchruledata->maxdiffofnodes = MAX(0, floor(log10(SCIPgetNBinVars(scip))/log10(2))-2);

      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/savelpbasis", &branchruledata->savelpbasis));
      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/saveglbcons", &branchruledata->saveglbcons));

      SCIP_CALL(SCIPgetIntParam(scip, "reoptimization/solvelp", &branchruledata->solvelp));
      SCIP_CALL(SCIPgetIntParam(scip, "reoptimization/solvelpdiff", &branchruledata->solvelpdiff));

      /** Savings for reoptimization */
      branchruledata->allocmemsizenodedata = MIN(500, branchruledata->maxsavednodes) + 1;
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
      branchruledata->nrestartsround = 0;
      branchruledata->nruns = 0;

      /** globals constraints */
      branchruledata->nconsglobal = 0;
      branchruledata->nconsglobaladded = 0;
      branchruledata->allocmemglobalcons = 50;

      /** clock */
      SCIP_CALL(SCIPclockCreate(&branchruledata->inittime, SCIP_CLOCKTYPE_DEFAULT));
      SCIP_CALL(SCIPclockCreate(&branchruledata->updatesolutime, SCIP_CLOCKTYPE_DEFAULT));
      SCIP_CALL(SCIPclockCreate(&branchruledata->savetime, SCIP_CLOCKTYPE_DEFAULT));

      /** data structure for nodes, init data for root node at position 0 */
      SCIP_CALL(SCIPallocClearMemoryArray(scip, &(branchruledata->nodedata), branchruledata->allocmemsizenodedata));

      /** allocate memory for the root */
      SCIP_CALL(SCIPallocMemory(scip, &branchruledata->nodedata[0]));
      branchruledata->nodedata[0]->conss = NULL;
      branchruledata->nodedata[0]->lpistate = NULL;
      branchruledata->nodedata[0]->nodechilds = NULL;
      branchruledata->nodedata[0]->nvars = 0;
      branchruledata->nodedata[0]->oldnumber = 1;
      branchruledata->nodedata[0]->pseudobranched = FALSE;
      branchruledata->nodedata[0]->redies = 0;
      branchruledata->nodedata[0]->allocmem = 0;
      branchruledata->nodedata[0]->vars = NULL;
      branchruledata->nodedata[0]->varbounds = NULL;
      branchruledata->nodedata[0]->varboundtypes = NULL;

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

      /** data for global logic-or constraints to separate binary solutions */
      SCIP_CALL(SCIPallocClearMemoryArray(scip, &(branchruledata->globalcons), branchruledata->allocmemglobalcons));

      /** change parameters */
//      SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrounds", 0));
//      SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrestarts", 0));
      SCIP_CALL(SCIPsetIntParam(scip, "propagating/maxrounds", 0));
      SCIP_CALL(SCIPsetIntParam(scip, "propagating/maxroundsroot", 0));

      if (branchruledata->savelpbasis)
      {
         SCIP_CALL(SCIPsetIntParam(scip, "separating/maxrounds", 0));
         SCIP_CALL(SCIPsetIntParam(scip, "separating/maxroundsroot", 0));
      }

      SCIP_CALL(SCIPsetBoolParam(scip, "conflict/enable", FALSE));
      SCIP_CALL(SCIPsetBoolParam(scip, "lp/disablecutoff", FALSE));
      SCIP_CALL(SCIPsetBoolParam(scip, "misc/calcintegral", FALSE));

      /** mark data structure initialized */
      branchruledata->init = TRUE;
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_BRANCHEXITSOL(branchExitsolnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL );
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL );

   if (branchruledata->init && branchruledata->reopt)
   {
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

      if (branchruledata->reopt)
      {
         branchruledata->nruns++;
         branchruledata->nrestartsround = 0;
      }
#ifdef CCHECK
      if( SCIPgetStage(scip) == SCIP_STAGE_EXITSOLVE &&
            SCIPgetStage(scip) != SCIP_STAGE_PROBLEM &&
            SCIPgetNTotalNodes(scip) > 1 &&
            SCIPhashmapGetNEntries(branchruledata->nodetoid) >= 1)
      {
         SCIP_CALL( checkConsistency(scip) );
      }
#endif
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
   if (branchruledata->init)
   {
      int nodeID;

      assert(branchruledata->reopt);

      /** free data of unrevived nodes */
      for (nodeID = 0; nodeID < branchruledata->allocmemsizenodedata; nodeID++)
      {
         if (branchruledata->nodedata[nodeID] != NULL )
         {
            SCIP_CALL(deleteNodeData(scip, branchruledata, nodeID, TRUE));
         }
      }

      /** free data for handle saved nodes */
      SCIPqueueFree(&branchruledata->openIDs);
      SCIPfreeMemoryArray(scip, &branchruledata->nodedata);

      /** free clock */
      SCIPfreeClock(scip, &branchruledata->inittime);
      SCIPfreeClock(scip, &branchruledata->updatesolutime);
      SCIPfreeClock(scip, &branchruledata->savetime);

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

   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   if (branchruledata->reopt)
   {
      if (branchruledata->nsavednodes > 0)
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

            if (SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip))
            {
               SCIP_CALL(Exec(scip, branchruledata));
               SCIP_CALL(SCIPbranchrulePseudoReset(scip, FALSE));
            }
            else
            {
               assert(1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));
               assert(SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)) < branchruledata->allocmemsizenodedata);
               SCIP_CALL(Exec(scip, branchruledata));
            }
            *result = SCIP_BRANCHED;

            /** stop clock */
            SCIP_CALL( SCIPstopInitTime(scip) );
            SCIP_CALL( SCIPstopClock(scip, branchruledata->inittime));
         }
      }
   }

   return SCIP_OKAY;
}

/** branching execution method for external candidates */
static SCIP_DECL_BRANCHEXECEXT(branchExecextnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   if (branchruledata->reopt)
   {
      if (branchruledata->nsavednodes > 0)
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

            if (SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip))
            {
               SCIP_CALL(Exec(scip, branchruledata));
               SCIP_CALL(SCIPbranchrulePseudoReset(scip, FALSE));
            }
            else
            {
               assert(1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));
               assert(SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)) < branchruledata->allocmemsizenodedata);
               SCIP_CALL(Exec(scip, branchruledata));
            }
            *result = SCIP_BRANCHED;

            /** stop clock */
            SCIP_CALL( SCIPstopInitTime(scip) );
            SCIP_CALL( SCIPstopClock(scip, branchruledata->inittime));
         }
      }
   }

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static SCIP_DECL_BRANCHEXECPS(branchExecpsnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   if (branchruledata->reopt)
   {
      if (branchruledata->nsavednodes > 0)
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

            if (SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip))
            {
               SCIP_CALL(Exec(scip, branchruledata));
               SCIP_CALL(SCIPbranchrulePseudoReset(scip, FALSE));
            }
            else
            {
               assert(1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));
               assert(SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)) < branchruledata->allocmemsizenodedata);
               SCIP_CALL(Exec(scip, branchruledata));
            }
            *result = SCIP_BRANCHED;

            /** stop clock */
            SCIP_CALL( SCIPstopInitTime(scip) );
            SCIP_CALL( SCIPstopClock(scip, branchruledata->inittime));
         }
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
   SCIP_CALL(SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolnodereopt));

   return SCIP_OKAY;
}
