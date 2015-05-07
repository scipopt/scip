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
 * @brief  branching rule to reconstruct the search tree
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>
#include <math.h>

#include "scip/branch_nodereopt.h"
#include "scip/branch_relpscost.h"
#include "scip/cons_logicor.h"
#include "scip/scip.h"
#include "scip/tree.h"
#include "scip/reopt.h"

#define BRANCHRULE_NAME            "nodereopt"
#define BRANCHRULE_DESC            "branching rule for node reoptimization"
#define BRANCHRULE_PRIORITY        100000
#define BRANCHRULE_MAXDEPTH            -1
#define BRANCHRULE_MAXBOUNDDIST         1.0

#define DEFAULT_USESPLITCONS        TRUE

/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Real             simsolverootlp;          /**< threshold to skip solving the root LP */
   SCIP_Bool             initialized;             /**< is the data structure initialized? */
   SCIP_Bool             reopt;                   /**< is reoptimization enabled?  */
   SCIP_Bool             strongbranchinginit;     /**< run a strong branching initialization? */
   SCIP_Bool             usesplitcons;            /**< use a constraint to handle dual bound changes */

   /** Statistic stuff */
   int                   nsplits;                 /**< number of nodes split by the branching rule */
   int                   nrevivednodes;           /**< number of nodes reoptimized by the branching rule */
};

/*
 *  static methods
 */


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
   SCIP_Bool localrestart;
   SCIP_NODE* node;
   SCIP_CONS* cons;

   int* childids;
   int nchilds;
   int nodeID;
   int childID;
   int c;
   int ncreatedchilds;

   assert(scip != NULL );
   assert(branchruledata != NULL );
   assert(branchruledata->reopt);

   node = SCIPgetCurrentNode(scip);

   nodeID = SCIPgetRootNode(scip) == node ? 0 : SCIPnodeGetReoptID(node);

   /* calculate local similarity and delete the induced subtree if
    * the similarity is to low */
   localrestart = FALSE;
   SCIP_CALL( SCIPcheckLocalRestart(scip, node, &localrestart) );

   if( localrestart )
   {
      *result = SCIP_DIDNOTRUN;
      goto TERMINATE;
   }

   SCIPdebugMessage("current node is %lld, ID %d:\n", SCIPnodeGetNumber(node), nodeID);

   /**
    * current node is equal to the root and the root was pseudo-branched
    * we have to create two child nodes; one with the pseudo-constraint and
    * one with the negated fixings.
    */
   if(SCIPgetRootNode(scip) == node && SCIPnodeSplit(scip, node) )
   {
      LOGICORDATA* consdata;
      int v;

      consdata = NULL;

      /* allocate buffer memory */
      SCIP_CALL( SCIPallocBuffer(scip, &consdata) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consdata->vars, SCIPgetNBinVars(scip)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consdata->vals, SCIPgetNBinVars(scip)) );
      consdata->varssize = SCIPgetNBinVars(scip);
      consdata->nvars = 0;

      /* get the constraint which splits the node */
      SCIP_CALL( SCIPnodeGetSplitCons(scip, nodeID, consdata) );

      /**
       * if the bound changes based on dual information induces an
       * infeasible subtree, we add this constraint as a global one.
       * otherwise, we split the root into two dummy nodes and
       * proceed as like the root is a transit node.
       */
      switch (consdata->constype) {
         case REOPT_CONSTYPE_STRBRANCHED:

            /* split the node */
            if( branchruledata->usesplitcons )
            {
               /* split the root into two dummy nodes */
               SCIP_CALL( SCIPsplitReoptRoot(scip) );
            }
            else
            {
               /* generate consdata->vars+1 nodes */
               SCIP_CALL( SCIPinterdictReoptRoot(scip, consdata->vars, consdata->vals, consdata->nvars) );
            }
         break;

         case REOPT_CONSTYPE_INFSUBTREE:
            assert( SCIPgetNReoptChildrenIDs(scip, node) == 0 );

            /* the logic-or constraint induces an infeasible subtree
             * and at the beginning of the next iteration the constraint can
             * be added globally.
             */
            SCIP_CALL( SCIPaddReoptGlbCons(scip, consdata) );

            /* create a logic-or constraint and add them to the current root */
            for(v = 0; v < consdata->nvars; v++)
            {
               if( SCIPisFeasEQ(scip, consdata->vals[v], 1) )
               {
                  assert(SCIPvarIsOriginal(consdata->vars[v]));
                  SCIP_CALL( SCIPgetNegatedVar(scip, consdata->vars[v], &consdata->vars[v]) );
                  assert(SCIPvarIsNegated(consdata->vars[v]));
               }
            }

            SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, "glbinfsub", consdata->nvars, consdata->vars,
                  TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE) );

            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );

            /* remove all dual information leading to split the node */
            SCIP_CALL( SCIPnodeReoptResetDualcons(scip, node) );

            /* we can return because nothing is to do */
            *result = SCIP_DIDNOTRUN;

            /* free buffer memory */
            SCIPfreeBufferArray(scip, &consdata->vals);
            SCIPfreeBufferArray(scip, &consdata->vars);
            SCIPfreeBuffer(scip, &consdata);

            goto TERMINATE;

            break;

         default:
            break;
      }

      /* free buffer memory */
      SCIPfreeBufferArray(scip, &consdata->vals);
      SCIPfreeBufferArray(scip, &consdata->vars);
      SCIPfreeBuffer(scip, &consdata);

      nodeID = 0;
      goto REVIVE;
   }
   else if(SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) && !SCIPnodeSplit(scip, node))
   {
      nodeID = 0;
      goto REVIVE;
   }

   assert(nodeID >= 1);

   REVIVE:

   /* get the IDs of all child nodes */
   SCIP_CALL( SCIPallocBufferArray(scip, &childids, SCIPgetNReoptChildrenIDs(scip, node)) );
   SCIP_CALL( SCIPgetReoptChildrenIDs(scip, node, childids, SCIPgetNReoptChildrenIDs(scip, node), &nchilds) );
   assert(SCIPgetNReoptChildrenIDs(scip, node) == nchilds);

   ncreatedchilds = 0;

   for(c = 0; c < nchilds; c++)
   {
      SCIP_NODE* child1;
      SCIP_NODE* child2;
      SCIP_REOPTTYPE reopttype;

      child1 = NULL;
      child2 = NULL;

      childID = childids[c];
      assert(childID >= 1);

      reopttype = SCIPreoptGetNodeType(scip, childID);

      SCIPdebugMessage("process child at ID %d\n", childID);

      /** check weather the constraints contains variable and if so, check the type of the constraint */
      if( reopttype == SCIP_REOPTTYPE_STRBRANCHED || reopttype == SCIP_REOPTTYPE_INFSUBTREE )
      {
        SCIPdebugMessage(" -> node %s\n", reopttype == SCIP_REOPTTYPE_STRBRANCHED ? "need to split" : "includes infeasible subtree");

         /* split the node into two disjoint node or cut off an infeasible subset */
         switch( reopttype ) {
            case SCIP_REOPTTYPE_STRBRANCHED:

               if( branchruledata->usesplitcons )
               {
                  /** the constraint split the node into two new nodes */
                  SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
                  SCIP_CALL( SCIPcreateChild(scip, &child2, 0.0, SCIPgetLocalTransEstimate(scip)) );
                  assert(child1 != NULL );
                  assert(child2 != NULL );

                  SCIPdebugMessage(" -> create 2 nodes: #%lld and #%lld\n", SCIPnodeGetNumber(child1),
                        SCIPnodeGetNumber(child2));

                  /* reoptimize the node:
                   *  - apply the bound changes along the stored branching path in both nodes
                   *  - optional: apply the all bound changes based on primal information caught between the
                   *    first and second bound change based on dual information to child1
                   *  - fix all bound changes based on dual information in child1 and add the
                   *    corresponding constraint to child 2
                   *  - add all local constraints to both nodes
                   */
                  SCIP_CALL( SCIPapplyReopt(scip, child1, child2, childID) );

                  ncreatedchilds += 2;
               }
               else
               {
                  SCIP_NODE** childs;
                  int* permutation;
                  int nchilds;
                  int n;

                  nchilds = SCIPgetReoptnodeNDualBoundChgs(scip, childID) + 1;

                  /* allocate buffer */
                  SCIP_CALL( SCIPallocBufferArray(scip, &childs, nchilds) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &permutation, nchilds) );

                  /* calculate the variable ordering
                   * @todo find a non-trivial ordering
                   */
                  for( n = 0; n < nchilds; n++ )
                  {
                     permutation[n] = n;
                  }

                  for( n = 0; n < nchilds; n++ )
                  {
                     SCIP_CALL( SCIPcreateChild(scip, &childs[n], 0.0, SCIPgetLocalTransEstimate(scip)) );
                  }

                  SCIPdebugMessage(" -> create %d nodes: %lld -- %lld\n", nchilds, SCIPnodeGetNumber(childs[0]), SCIPnodeGetNumber(childs[nchilds-1]));

                  /* reoptimize the node:
                   *  - apply the bound changes along the stored branching path in all nodes
                   *  - optional: apply the all bound changes based on primal information caught between the
                   *    first and second bound change based on dual information to child1
                   *  - fix all bound changes based on dual information in childs[0]
                   *  - fix the first c bound changes in child[c] as in childs[0] and negate the (c+1)th bound change in child[c]
                   *  - add all local constraints to all nodes
                   */
                  SCIP_CALL( SCIPapplyReoptInterdiction(scip, childID, childs, nchilds, permutation) );

                  ncreatedchilds += nchilds;

                  /* remember the reconstructed original node*/
                  child1 = childs[0];

                  /* free buffer */
                  SCIPfreeBufferArray(scip, &permutation);
                  SCIPfreeBufferArray(scip, &childs);

                  assert(child1 != NULL);
               }

               break;

            case SCIP_REOPTTYPE_INFSUBTREE:

               /* create only one child node */
               SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
               assert(child1 != NULL);
               assert(child2 == NULL);

               SCIPdebugMessage(" -> create 1 node: #%lld\n", SCIPnodeGetNumber(child1));

               /* reoptimize the node:
                *  - apply the bound changes along the stored branching path to child1
                *  - add constraint to child 1 to cut off the infeasible subtree
                *  - add all local constraints to child1 */
               SCIP_CALL( SCIPapplyReopt(scip, NULL, child1, childID) );

               ncreatedchilds++;

               break;

            default:
               assert(reopttype == SCIP_REOPTTYPE_STRBRANCHED || reopttype == SCIP_REOPTTYPE_STRBRANCHED);
               break;
         }
      }
      else
      {
         /**
          * node at position childID includes no bound changes based on dual information
          */
         SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
         assert(child1 != NULL );
         assert(child2 == NULL );

         /* reoptimize the node:
          *  - apply the bound changes along the stored branching path to child1
          *  - add all local constraints to child1 */
         SCIP_CALL( SCIPapplyReopt(scip, child1, NULL, childID) );

         ncreatedchilds++;
      }

#ifdef SCIP_DISABLED_CODE
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
      branchruledata->nodedata[childID]->dualfixing = (branchruledata->nodedata[childID]->dualconscur != NULL);
#endif

      /* update the reoptid */
      SCIPnodeSetReoptID(child1, childID);

      /* set the REOPTTYPE */
      assert(reopttype >= SCIP_REOPTTYPE_TRANSIT);
      if( reopttype == SCIP_REOPTTYPE_STRBRANCHED )
         reopttype = SCIP_REOPTTYPE_TRANSIT;
      SCIPnodeSetReopttype(child1, reopttype);
      SCIPdebugMessage(" -> set reopttype: %d\n", reopttype);

      /** check if child2 includes some added constraints and save the node */
      if(child2 != NULL && SCIPnodeGetNAddedConss(child2) > 0)
      {
         SCIP_CALL( SCIPaddReoptnode(scip, child2, SCIP_REOPTTYPE_LEAF, TRUE, SCIPgetLocalTransEstimate(scip)) );
      }
   }

   if( ncreatedchilds == 0 )
      *result = SCIP_DIDNOTRUN;
   else
   {
      int ncreatedchilds_tmp;
      ncreatedchilds_tmp = ncreatedchilds;
      /* increase the counter for reoptimized nodes */
      while( ncreatedchilds_tmp > 0 )
      {
         SCIPenforceNReoptnodes(scip);
         --ncreatedchilds_tmp;
      }
      *result = SCIP_BRANCHED;
   }

   /* free the buffer memory */
   SCIPfreeBufferArray(scip, &childids);

  TERMINATE:

   SCIPdebugMessage("**** finish reoptimizing %d child nodes of node %lld ****\n", ncreatedchilds, SCIPnodeGetNumber(node));

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */

#define branchCopynodereopt NULL;
#define branchExitnodereopt NULL;

static
SCIP_DECL_BRANCHINIT(branchInitnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL );
   assert(branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL );

   /** check if reoptimization is enabled */
   if (!branchruledata->initialized)
   {
      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/enable", &branchruledata->reopt));
   }

   /** initialize the data and change parameters */
   if (!branchruledata->initialized && branchruledata->reopt)
   {
      /** statistic */
      branchruledata->nrevivednodes = 0;
      branchruledata->nsplits = 0;

      /** mark data structure initialized */
      branchruledata->initialized = TRUE;

      /** get parameters  */
      SCIP_CALL( SCIPgetBoolParam(scip, "reoptimization/strongbranchinginit", &branchruledata->strongbranchinginit) );
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

   /** free data structure only if reoptimization is enabled */
   if (branchruledata->initialized )
   {
      assert(branchruledata->reopt);

      branchruledata->initialized = FALSE;
      branchruledata->reopt = FALSE;
   }
   assert(!branchruledata->initialized);

   SCIPfreeMemory(scip, &branchruledata);

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;
#ifdef SCIP_DISABLED_CODE
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   SCIP_Real* branchcandsfrac;
   int nbranchcands;
#endif

   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   if ( branchruledata->reopt && SCIPreoptimizeNode(scip, SCIPgetCurrentNode(scip)) )
   {

#ifdef SCIP_DISABLED_CODE
      if( branchruledata->strongbranchinginit
       && SCIPreoptGetSimToPrevious(scip->reopt) < branchruledata->simsolverootlp
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
#endif

      if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM)
      {
         if( SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip) )
         {
            SCIP_CALL( Exec(scip, branchruledata, result) );
         }
         else
         {
            assert(1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));

            SCIP_CALL( Exec(scip, branchruledata, result) );
         }
      }
   }

   return SCIP_OKAY;
}

/** branching execution method for external candidates */
static SCIP_DECL_BRANCHEXECEXT(branchExecextnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;
#ifdef SCIP_DISABALED_CODE
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   SCIP_Real* branchcandsfrac;
   int nbranchcands;
#endif

   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   if ( branchruledata->reopt && SCIPreoptimizeNode(scip, SCIPgetCurrentNode(scip)) )
   {
#ifdef SCIP_DISABLED_CODE
      if( branchruledata->strongbranchinginit
       && SCIPreoptGetSimToPrevious(scip->reopt) < branchruledata->simsolverootlp
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
#endif

      if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM)
      {
         if( SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip) )
         {
            SCIP_CALL( Exec(scip, branchruledata, result) );
         }
         else
         {
            assert(1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));

            SCIP_CALL( Exec(scip, branchruledata, result) );
         }
      }
   }

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static SCIP_DECL_BRANCHEXECPS(branchExecpsnodereopt)
{
   SCIP_BRANCHRULEDATA* branchruledata;
#ifdef SCIP_DISABLED_CODE
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   SCIP_Real* branchcandsfrac;
   int nbranchcands;
#endif

   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   branchruledata = SCIPbranchruleGetData(branchrule);

   *result = SCIP_DIDNOTRUN;

   if ( branchruledata->reopt && SCIPreoptimizeNode(scip, SCIPgetCurrentNode(scip)) )
   {
#ifdef SCIP_DISABLED_CODE
      if( branchruledata->strongbranchinginit
       && SCIPreoptGetSimToPrevious(scip->reopt) < branchruledata->simsolverootlp
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
#endif

      if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM )
      {
         if( SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip) )
         {
            SCIP_CALL( Exec(scip, branchruledata, result) );
         }
         else
         {
            assert(1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));

            SCIP_CALL( Exec(scip, branchruledata, result) );
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
   branchruledata->initialized = FALSE;

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

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/"BRANCHRULE_NAME"/usesplitcons", "use a constraint to handle dual bound change, otherwise use some kind of interdiction branching",
         &branchruledata->usesplitcons, TRUE, DEFAULT_USESPLITCONS, NULL, NULL) );

   return SCIP_OKAY;
}
