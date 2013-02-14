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

/**@file   branch_ryanfoster.c
 * @ingroup BRANCHINGRULES
 * @brief  Ryan/Foster branching rule
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file implements the Ryan/Foster branching rule. For more details see \ref BRANCHING page.
 *
 * @page BRANCHING Ryan/Foster branching
 *
 * Ryan/Foster branching is a very useful branching rule  for the integer program model in use. A
 * standard variable branching has the disadvantage that the zero branch is more or less useless because
 * we only forbid one packing out of exponential many. On the other hand, the branch fixing a packing reduces the problem since
 * certain items are packed. This leads to a very unbalanced search tree.
 *
 * The branching idea of Ryan/Foster is to select a pair of items which is either a) forced to be packed together or b)
 * not allowed to be packed together. Note that in both cases, it is allowed to use packings
 * which contain none of the two items.
 *
 * There are two issues to be taken care off:
 * -# How do we select the pair of items?
 * -# How do we realize such a branching within \SCIP?
 *
 * @section SELECTION How do we select the pair of items?
 *
 * To select a pair of items, we have to know for each packing the items which are contained. Since every packing is a
 * variable and each item is a set covering constraint, we have to know for each variable in which set covering
 * constraints it appears (this means, has a coefficient of 1.0). Since \SCIP is constraint based, it is in general
 * not possible to get this information directly. To overcome this issue, we use the functionality to add
 * \ref vardata_binpacking.c "variable data" to every
 * variable. This variable data contains the constraints in which this variable appears (see vardata_binpacking.c for more details).
 * With the help of the variable data, it is now possible to get the
 * information which items belong to which packing. Therefore, we can use the Ryan/Foster idea to select a pair of
 * items.
 *
 * @section SAMEDIFFBRANCHING How do we realize such a branching within SCIP?
 *
 * After having selected a pair of items to branch on, the question now is how to realize such a branching with \SCIP.
 * Since \SCIP is
 * constraint based, it is really easy to do that. We implement a constraint handler which handles the
 * information, see cons_samediff.c. This constraint handler does not only store the branching
 * decisions. Furthermore, it also ensures that all packing which are not feasible at a particular node are
 * locally fixed to zero. For more details, we refer to the \ref cons_samediff.c "source code of the constraint handler".
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_ryanfoster.h"
#include "cons_samediff.h"
#include "probdata_binpacking.h"
#include "vardata_binpacking.h"

/**@name Branching rule properties
 *
 * @{
 */

#define BRANCHRULE_NAME            "RyanFoster"
#define BRANCHRULE_DESC            "Ryan/Foster branching rule"
#define BRANCHRULE_PRIORITY        50000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

/**@} */

/**@name Callback methods
 *
 * @{
 */

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpRyanFoster)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;

   SCIP_Real** pairweights;

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

      /* loop over all constraints/items the variable belongs to */
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

/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the ryan foster branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleRyanFoster(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create ryan foster branching rule data */
   branchruledata = NULL;
   branchrule = NULL;
   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
         BRANCHRULE_MAXBOUNDDIST, branchruledata) );
   assert(branchrule != NULL);

   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpRyanFoster) );

   return SCIP_OKAY;
}

/**@} */
