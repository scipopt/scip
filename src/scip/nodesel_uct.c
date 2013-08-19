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

/**@file   nodesel_uct.c
 * @brief  uct node selector
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/nodesel_uct.h"


#define NODESEL_NAME            "uct"
#define NODESEL_DESC            "node selector template"
#define NODESEL_STDPRIORITY     10
#define NODESEL_MEMSAVEPRIORITY 0
#define DEFAULT_WEIGHT          0.1
#define DEFAULT_NODELIMIT       31
#define EVENTHDLR_NAME          "event_uct"
#define EVENTHDLR_EVENTTYPE     SCIP_EVENTTYPE_NODEFOCUSED

/*
 * Data structures
 */

/* TODO: fill in the necessary node selector data */

/** node selector data */
/* todo: extra int for array length because nodelimit is a user parameter */
struct SCIP_NodeselData
{
   int*  nodevisits;
   int   nodelimit;
   int   nselections;
   SCIP_Real weight;
   SCIP_NODE* lastfocusnode;
};

#if 0
struct SCIP_EventhdlrData
{
   SCIP_NODESEL* nodesel;
};
#endif
/*
 * Local methods
 */

/* put your local methods here, and declare them static */
static
int nodeGetVisits(
   SCIP_NODESEL* nodesel,
   SCIP_NODE*   node
   )
{
   SCIP_NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   if( SCIPnodeGetNumber(node) >= nodeseldata->nodelimit )
      return 0;
   else
      return nodeseldata->nodevisits[(int)SCIPnodeGetNumber(node)];
}

static
void backpropagateVisits(
   SCIP_NODESELDATA* nodeseldata
   )
{
   SCIP_NODE* pathnode;
   int* visits;

   assert(nodeseldata != NULL);
   visits = nodeseldata->nodevisits;
   assert(visits != NULL);

   pathnode = nodeseldata->lastfocusnode;
   assert(pathnode != NULL);
   do
   {
      SCIP_Longint nodenumber;

      nodenumber = SCIPnodeGetNumber(pathnode);
      if( nodenumber < nodeseldata->nodelimit )
         ++(visits[(int)nodenumber]);

      pathnode = SCIPnodeGetParent(pathnode);
   }
   while( pathnode != NULL );
}

/*
 * Callback methods of node selector
 */

/* TODO: Implement all necessary node selector methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_NODESELCOPY(nodeselCopyUct)
{  /*lint --e{715}*/
   assert(scip != NULL);
   SCIP_CALL( SCIPincludeNodeselUct(scip) );

   return SCIP_OKAY;
}

static
SCIP_DECL_NODESELINITSOL(nodeselInitsolUct)
{
   SCIP_NODESELDATA* nodeseldata;
   assert(scip != NULL);
   assert(nodesel != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);

   assert(nodeseldata != NULL);
   if( nodeseldata->nodevisits == NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(nodeseldata->nodevisits), 2 * nodeseldata->nodelimit) );
   }
   BMSclearMemoryArray(nodeseldata->nodevisits, 2 * nodeseldata->nodelimit);
   nodeseldata->lastfocusnode = NULL;
   nodeseldata->nselections = 0;

   return SCIP_OKAY;
}

static
SCIP_DECL_NODESELFREE(nodeselFreeUct)
{
   SCIP_NODESELDATA* nodeseldata;
   assert(scip != NULL);
   assert(nodesel != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);
   if( nodeseldata->nodevisits != NULL )
   {
      SCIPfreeMemoryArray(scip, &nodeseldata->nodevisits);
   }
   SCIPfreeBlockMemory(scip, &nodeseldata);

   SCIPnodeselSetData(nodesel, NULL);

   return SCIP_OKAY;
}



/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectUct)
{

   SCIP_NODESELDATA* nodeseldata;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   /* select next node as best node with respect to UCT-based comparison method */
   *selnode = SCIPgetBestNode(scip);

   if(*selnode == NULL)
      return SCIP_OKAY;

   /* count the number of selections */
   ++nodeseldata->nselections;

   /* switch back to default node selection rule if the node limit is exceeded */
   if( nodeseldata->nselections == nodeseldata->nodelimit )
   {
      SCIP_NODESEL** nodesels;
      int nnodesels;
      int newpriority;
      int n;

      nodesels = SCIPgetNodesels(scip);
      nnodesels = SCIPgetNNodesels(scip);
      newpriority = SCIPnodeselGetStdPriority(nodesel);
      for( n = 0; n < nnodesels; ++n )
      {
         newpriority = MIN(newpriority, SCIPnodeselGetStdPriority(nodesels[n]));
      }
      SCIPdebugMessage("Reached node limit of UCT node selection rule -> switching to default\n");
      SCIP_CALL( SCIPsetNodeselStdPriority(scip, nodesel, newpriority - 1) );
   }
   else if( SCIPnodeGetType(*selnode) != SCIP_NODETYPE_CHILD && nodeseldata->lastfocusnode != NULL )
   {
      /* trigger backpropagation of visits upwards from last focus node if the newly selected node is not a child node
       * of the focus node
       */
      SCIPdebugMessage("Backpropagating node visits from node number %"SCIP_LONGINT_FORMAT"\n", SCIPnodeGetNumber(nodeseldata->lastfocusnode));
      backpropagateVisits(nodeseldata);
   }

   /* keep selected node as last focus node for next time */
   nodeseldata->lastfocusnode = *selnode;

   return SCIP_OKAY;
}

/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompUct)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;
   SCIP_NODE* parent1;
   SCIP_NODE* parent2;
   SCIP_Real score1;
   SCIP_Real score2;
   int nodevisits1;
   int nodevisits2;
   int parentvisits1;
   int parentvisits2;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   score1 = SCIPnodeGetLowerbound(node1) / SCIPgetLowerboundRoot(scip);
   score2 = SCIPnodeGetLowerbound(node2) / SCIPgetLowerboundRoot(scip);

   nodevisits1 = nodeGetVisits(nodesel, node1);
   nodevisits2 = nodeGetVisits(nodesel, node2);

   parent1 = SCIPnodeGetParent(node1);
   parent2 = SCIPnodeGetParent(node2);

   assert(parent1 != NULL && parent2 != NULL);

   parentvisits1 = nodeGetVisits(nodesel, parent1);
   parentvisits2 = nodeGetVisits(nodesel, parent2);
   SCIPdebugMessage("Comparing nodes: %10.2g %d %d -- %10.2g %d %d (lower bound, visits, parent visits)\n",
         SCIPnodeGetLowerbound(node1), nodevisits1, parentvisits1, SCIPnodeGetLowerbound(node2), nodevisits2, parentvisits2 );
   if( parentvisits1 > 0 || parentvisits2 > 0 )
   {
      score1 -= nodeseldata->weight * parentvisits1 / (SCIP_Real)(1 + nodevisits1);
      score2 -= nodeseldata->weight * parentvisits2 /(SCIP_Real)(1 + nodevisits2);
   }
   if( (SCIPisInfinity(scip, score1) && SCIPisInfinity(scip, score2)) ||
         (SCIPisInfinity(scip, -score1) && SCIPisInfinity(scip, -score2))
         || SCIPisEQ(scip, score1, score2) )
   {
      return 0;
   }
   if( SCIPisLT(scip, score1, score2) )
      return -1;
   else
   {
      assert(SCIPisGT(scip, score1, score2));
      return 1;
   }
}

/*
 * node selector specific interface methods
 */

/** creates the uct node selector and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselUct(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESELDATA* nodeseldata;
   SCIP_NODESEL* nodesel;

   /* create uct node selector data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nodeseldata) );

   nodesel = NULL;
   nodeseldata->nodevisits = NULL;
   /* use SCIPincludeNodeselBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeNodeselBasic(scip, &nodesel, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
          nodeselSelectUct, nodeselCompUct, nodeseldata) );

   assert(nodesel != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetNodeselCopy(scip, nodesel, nodeselCopyUct) );
   SCIP_CALL( SCIPsetNodeselInitsol(scip, nodesel, nodeselInitsolUct) );
   SCIP_CALL( SCIPsetNodeselFree(scip, nodesel, nodeselFreeUct) );

   /* add uct node selector parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "nodeselection/"NODESEL_NAME"/nodelimit", "maximum number of nodes before switching to default rule",
         &nodeseldata->nodelimit, TRUE, DEFAULT_NODELIMIT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "nodeselection/"NODESEL_NAME"/weight", "weight for visit quotient of node selection rule",
         &nodeseldata->weight, TRUE, DEFAULT_WEIGHT, 0.0, 1.0, NULL, NULL) );
   return SCIP_OKAY;
}
