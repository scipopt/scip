/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_ryanfoster.c
 * @ingroup BRANCHINGRULES
 * @brief  Ryan/Foster branching rule
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 *
 * @page BRANCHING Ryan/Foster branching
 *
 * Ryan/Foster branching is one of the branching rules which are very useful for the used integer programs model. A
 * standard variable branching has the disadvantage that the zero branch is more or less useless. This is the case since
 * we only forbid one packing out of exponential many. The one branch on the other side reduces the problem since
 * certain items are packed. This leads to an very unbalanced search tree.
 *
 * The idea of Ryan/Foster is to branch in a way that we say that on the one branch a certain pair of items are always
 * together and an the other branch they are never together. Note that in both case it is allowed that packings are
 * used which contain none of the two items.
 *
 * There are two issue to be taken care off:
 * -# How do we select the pair of itmes?
 * -# How do we realize such a branching within \SCIP?
 *
 * @section SELECTION How do we select the pair of items?
 *
 * To select a pair of items, we have to know for each packing to items which are contained. Since every packing is a
 * variable and each item is a set covering constraint, we have to know for each variable in which set covering
 * constraints it appears (this means, has a coefficient of 1.0). Since \SCIP is constraint based, it is not possible to
 * get this information in general. To overcome this issue we use to functionality to add variable data to each
 * variable. This variable data contains the required information. This means, the constraints in which this variable
 * appears (see vardata_binpacking.c for more details). Having this variable data, it is now possible to get the
 * information which items belong to which packing. Therefore, we can use the Ryan/Foster idea to select a pair of
 * items.
 *
 * @section SAMEDIFFBRANCHING How do we realize such a branching within SCIP?
 *
 * After we selected a pair of items to branch on, the questions how to realize that with \SCIP. Since \SCIP is
 * constraint based, it is really easy to do that. We implement a constraint handler which handles these
 * informations. Therefore, see cons_samediff.c. This constraint handler does not only stores the branching
 * decisions. It also takes care of the fact that in each node all packing which are not feasible for that node a fixed
 * locally to zero. For more details we refer to the source code of the constraint handler.
 * 
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_ryanfoster.h"
#include "cons_samediff.h"
#include "probdata_binpacking.h"
#include "vardata_binpacking.h"


#define BRANCHRULE_NAME            "RyanFoster"
#define BRANCHRULE_DESC            "Ryan/Foster branching rule"
#define BRANCHRULE_PRIORITY        50000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0


/*
 * Data structures
 */

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
};

/*
 * Callback methods of branching rule
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BRANCHCOPY(branchCopyRyanFoster)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleRyanFoster(scip) );
 
   return SCIP_OKAY;
}
#else
#define branchCopyRyanFoster NULL
#endif

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreeRyanFoster NULL

/** initialization method of branching rule (called after problem was transformed) */
#define branchInitRyanFoster NULL

/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitRyanFoster NULL

/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolRyanFoster NULL

/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolRyanFoster NULL

/** TS: whats this ?*/
#define branchExecrelRyanFoster NULL

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpRyanFoster)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
      
   SCIP_Real** pairweights;
   int npairweights;
   
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   SCIP_Real bestvalue;
   SCIP_Real value;

   SCIP_NODE* childsame;
   SCIP_NODE* childdiffer;
   SCIP_CONS* conssame;
   SCIP_CONS* consdiffer;

   SCIP_VARDATA* vardata;
   int* consids;
   int nconsids;
   int nitems;
   
   int id1;
   int id2;

   int i;
   int j;
   int v;

   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(result != NULL);
   
   SCIPdebugMessage("start branching at node %"SCIP_LONGINT_FORMAT", depth %d\n", SCIPgetNNodes(scip), SCIPgetDepth(scip));

   *result = SCIP_DIDNOTRUN;
   
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   
   nitems = SCIPprobdataGetNItems(probdata);
   npairweights = (nitems*nitems - 3*nitems +2) / 2;

   /* allocate memory for triangle matrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &pairweights, nitems) );
   for( i = 0; i < nitems; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &pairweights[i], nitems) );
      
      for( j = 0; j < nitems; ++j )
         pairweights[i][j] = 0.0;
   }
      
   /* get fractional LP candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands) );
   assert(nlpcands > 0);

   /* compute weigthts for each order pair */
   for( v = 0; v < nlpcands; ++v )
   {
      assert(lpcands[v] != NULL);

      /* get variable data which contains the information to which constraints/items the variable belongs */
      vardata = SCIPvarGetData(lpcands[v]);

      consids = SCIPvardataGetConsids(vardata);
      nconsids = SCIPvardataGetNConsids(vardata);
      assert(nconsids > 0);

      /* loop over all constraints/itmes the variable belongs to */
      for( i = 0; i < nconsids; ++i )
      {
         id1 = consids[i];
         for( j = i+1; j < nconsids; ++j )
         {
            id2 = consids[j];
            assert(id1 < id2);

            pairweights[id2][id1] += lpcandsfrac[v];

            assert( SCIPisFeasLE(scip, pairweights[id2][id1], 1.0) );
            assert( SCIPisFeasGE(scip, pairweights[id2][id1], 0.0) );
         }
      }   
   }
   
   /* select branching */
   bestvalue = 0.0;
   id1 = -1;
   id2 = -1;

   for( i = 0; i < nitems; ++i )
   {
      for( j = 0; j < i+1; ++j )
      {
         value = MIN(pairweights[i][j], 1-pairweights[i][j]);

         if( bestvalue < value )
         {
            bestvalue = value;
            id1 = j;
            id2 = i;
         }
      }
   }

   assert( bestvalue > 0.0 );
   assert( id1 >= 0 && id1 < nitems);
   assert( id2 >= 0 && id2 < nitems);

   /* free memory for triangle matrix */
   for( i = 0; i < nitems; ++i )
   {
      SCIPfreeBufferArray(scip, &pairweights[i]);
   }
   SCIPfreeBufferArray(scip, &pairweights);
   
   SCIPdebugMessage("branch on order pair <%d,%d> with weight <%g>\n", 
      SCIPprobdataGetIds(probdata)[id1], SCIPprobdataGetIds(probdata)[id2], bestvalue);

   /* create the branch-and-bound tree child nodes of the current node */
   SCIP_CALL( SCIPcreateChild(scip, &childsame, 0.0, SCIPgetLocalTransEstimate(scip)) );
   SCIP_CALL( SCIPcreateChild(scip, &childdiffer, 0.0, SCIPgetLocalTransEstimate(scip)) );

   /* create corresponding constraints */
   SCIP_CALL( SCIPcreateConsSamediff(scip, &conssame, "same", id1, id2, SAME, childsame, TRUE) );
   SCIP_CALL( SCIPcreateConsSamediff(scip, &consdiffer, "differ", id1, id2, DIFFER, childdiffer, TRUE) );

  /* add constraints to nodes */
   SCIP_CALL( SCIPaddConsNode(scip, childsame, conssame, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childdiffer, consdiffer, NULL) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &conssame) );
   SCIP_CALL( SCIPreleaseCons(scip, &consdiffer) );

  *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** branching execution method for relaxation solutions */
#define branchExecrelRyanFoster NULL

/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsRyanFoster NULL

/*
 * branching rule specific interface methods
 */

/** creates the ryan foster branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleRyanFoster(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create ryan foster branching rule data */
   branchruledata = NULL;
   
   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH, 
	 BRANCHRULE_MAXBOUNDDIST,
         branchCopyRyanFoster,
         branchFreeRyanFoster, branchInitRyanFoster, branchExitRyanFoster,
         branchInitsolRyanFoster, branchExitsolRyanFoster,
         branchExeclpRyanFoster, branchExecrelRyanFoster, branchExecpsRyanFoster,
         branchruledata) );

   /* add ryan foster branching rule parameters */
   /* TODO: (optional) add branching rule specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
