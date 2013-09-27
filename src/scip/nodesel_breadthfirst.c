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

/**@file nodesel_breadthfirst.h
 * @ingroup NODESELECTORS
 * @brief node selector for breadth-first search
 * @author Stefan Heinz
 * @author Gregor Hendel
 *
 * This node selector performs breadth-first search, i.e., it completely evaluates an entire level of the search tree before
 * proceeding to the next level. At one level, nodes are processed in the order they were created by the branching rule.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/nodesel_breadthfirst.h"


#define NODESEL_NAME             "breadthfirst"
#define NODESEL_DESC             "breadth first search"
#define NODESEL_STDPRIORITY        -10000
#define NODESEL_MEMSAVEPRIORITY  -1000000

/** node selector data */
struct SCIP_NodeselData
{
   SCIP_Longint          nodelimit;          /**< limit of node selections after which node selection is turned off */
};


/** switches to a different node selection rule by assigning the lowest priority of all node selectors to bfs */
static
SCIP_RETCODE turnoffNodeSelector(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel             /**< the node selector to be turned off */
   )
{
   SCIP_NODESEL** nodesels;
   int nnodesels;
   int newpriority;
   int prio;
   int n;

   nodesels = SCIPgetNodesels(scip);
   nnodesels = SCIPgetNNodesels(scip);
   newpriority = SCIPnodeselGetStdPriority(nodesel);

   /* loop over node selectors to find minimum priority */
   for( n = 0; n < nnodesels; ++n )
   {
      prio = SCIPnodeselGetStdPriority(nodesels[n]);
      newpriority = MIN(newpriority, prio);
   }
   newpriority = MAX(newpriority, INT_MIN + 1);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Reached node limit of %s node selection rule -> switching to default\n",
         SCIPnodeselGetName(nodesel));
   SCIP_CALL( SCIPsetNodeselStdPriority(scip, nodesel, newpriority - 1) );

   return SCIP_OKAY;
}
/*
 * Callback methods
 */

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_NODESELCOPY(nodeselCopyBreadthfirst)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);

   /* call inclusion method of node selector */
   SCIP_CALL( SCIPincludeNodeselBreadthfirst(scip) );

   return SCIP_OKAY;
}

/** destructor of node selector to free user data (called when SCIP is exiting) */
static
SCIP_DECL_NODESELFREE(nodeselFreeBreadthfirst)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   /* free user data of node selector */
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   SCIPfreeMemory(scip, &nodeseldata);
   SCIPnodeselSetData(nodesel, nodeseldata);

   return SCIP_OKAY;
}

/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectBreadthfirst)
{  /*lint --e{715}*/

   SCIP_NODESELDATA* nodeseldata;
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);
   /* siblings come before leaves at the same level. Sometimes it can occur that no leaves are left except for children */
   *selnode = SCIPgetBestSibling(scip);
   if( *selnode == NULL )
   {
      *selnode = SCIPgetBestLeaf(scip);
      if( *selnode == NULL )
         *selnode=SCIPgetBestChild(scip);
   }
   if( *selnode != NULL )
   {
      SCIPdebugMessage("Selecting next node number %"SCIP_LONGINT_FORMAT" at depth %d\n", SCIPnodeGetNumber(*selnode), SCIPnodeGetDepth(*selnode));
   }

   if( nodeseldata->nodelimit > -1 && SCIPgetNNodes(scip) >= nodeseldata->nodelimit )
   {
      SCIP_CALL( turnoffNodeSelector(scip, nodesel) );
   }

   return SCIP_OKAY;
}


/** node comparison method of breadth first search: nodes with lower depth are preferred; in case of a tie, the node
 *  which was created earlier (and therefore has a smaller node number) is preferred */
static
SCIP_DECL_NODESELCOMP(nodeselCompBreadthfirst)
{  /*lint --e{715}*/
   int depth1;
   int depth2;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   depth1 = SCIPnodeGetDepth(node1);
   depth2 = SCIPnodeGetDepth(node2);

   /* if depths differ, prefer node with smaller depth */
   if( depth1 < depth2 )
      return -1;
   else if( depth1 > depth2 )
      return +1;
   else
   {
      /* depths are equal; prefer node with smaller number */
      SCIP_Longint number1;
      SCIP_Longint number2;

      number1 = SCIPnodeGetNumber(node1);
      number2 = SCIPnodeGetNumber(node2);
      assert(number1 != number2);

      if( number1 < number2 )
         return -1;
      else
         return +1;
   }
}

/*
 * breadth first specific interface methods
 */

/** creates the node selector for breadth first search and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselBreadthfirst(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESELDATA* nodeseldata;
   SCIP_NODESEL* nodesel;

   /* create breadthfirst node selector data */
   nodeseldata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &nodeseldata) );
   /* include node selector */
   SCIP_CALL( SCIPincludeNodeselBasic(scip, &nodesel, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
         nodeselSelectBreadthfirst, nodeselCompBreadthfirst, nodeseldata) );

   assert(nodesel != NULL);

   /* set non-fundamental callback functions via setter functions */
   SCIP_CALL ( SCIPsetNodeselCopy(scip, nodesel, nodeselCopyBreadthfirst) );
   SCIP_CALL ( SCIPsetNodeselFree(scip, nodesel, nodeselFreeBreadthfirst) );

   SCIP_CALL( SCIPaddLongintParam(scip, "nodeselection/"NODESEL_NAME"/nodelimit",
         "maximum number of nodes before switching to default rule (-1 for no limit)",
         &nodeseldata->nodelimit, TRUE, -1, -1, SCIP_LONGINT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
