/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nodesel_bfs.c,v 1.23 2003/12/15 17:45:33 bzfpfend Exp $"

/**@file   nodesel_bfs.c
 * @brief  node selector for best first search
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "nodesel_bfs.h"


#define NODESEL_NAME             "bfs"
#define NODESEL_DESC             "best first search"
#define NODESEL_STDPRIORITY      100000
#define NODESEL_MEMSAVEPRIORITY       0
#define NODESEL_LOWESTFIRST        TRUE   /**< are the nodes sorted such that the lowest bound node comes first? */




/*
 * Default parameter settings
 */

#define MINPLUNGEDEPTH               10 /**< minimal plunging depth, before new best node may be selected */
#define MAXPLUNGEDEPTH               25 /**< maximal plunging depth, before new best node is forced to be selected */
#define MAXPLUNGEQUOT               1.5 /**< maximal quotient (actlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where plunging is performed */



/** node selector data for best first search node selection */
struct NodeselData
{
   Real             maxplungequot;      /**< maximal quotient (actlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where plunging is performed */
   int              maxplungedepth;     /**< maximal plunging depth, before new best node is forced to be selected */
   int              minplungedepth;     /**< minimal plunging depth, before new best node may be selected */
};



/*
 * Callback methods
 */

static
DECL_NODESELFREE(nodeselFreeBfs)
{  /*lint --e{715}*/
   NODESELDATA* nodeseldata;

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

static
DECL_NODESELSELECT(nodeselSelectBfs)
{  /*lint --e{715}*/
   NODESELDATA* nodeseldata;
   Real actlowerbound;
   Real avglowerbound;
   Real lowerbound;
   int plungedepth;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   *selnode = NULL;

   /* get node selector user data */
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   /* check, if we want to plunge once more */
   actlowerbound = SCIPgetLocalLowerbound(scip);
   avglowerbound = SCIPgetAvgLowerbound(scip);
   lowerbound = SCIPgetLowerbound(scip);
   plungedepth = SCIPgetPlungeDepth(scip);
   if( plungedepth < nodeseldata->minplungedepth
      || (plungedepth < nodeseldata->maxplungedepth
         && (SCIPisEQ(scip, avglowerbound, lowerbound)
            || (actlowerbound - lowerbound)/(avglowerbound - lowerbound) < nodeseldata->maxplungequot)) )
   {
      /* we want to plunge again: prefer children over siblings, and siblings over leaves */
      *selnode = SCIPgetBestChild(scip);
      if( *selnode == NULL )
      {
         *selnode = SCIPgetBestSibling(scip);
         if( *selnode == NULL )
         {
            *selnode = SCIPgetBestLeaf(scip);
         }
      }
   }
   else
   {
      /* we don't want to plunge again: select best node from the tree */
      *selnode = SCIPgetBestNode(scip);
   }

   return SCIP_OKAY;
}

static
DECL_NODESELCOMP(nodeselCompBfs)
{  /*lint --e{715}*/
   Real lowerbound1;
   Real lowerbound2;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   lowerbound1 = SCIPnodeGetLowerbound(node1);
   lowerbound2 = SCIPnodeGetLowerbound(node2);
   if( lowerbound1 < lowerbound2 )
      return -1;
   else if( lowerbound1 > lowerbound2 )
      return +1;
   else
   {
      NODETYPE nodetype1;
      NODETYPE nodetype2;

      nodetype1 = SCIPnodeGetType(node1);
      nodetype2 = SCIPnodeGetType(node2);
      if( nodetype1 == SCIP_NODETYPE_CHILD && nodetype2 != SCIP_NODETYPE_CHILD )
         return -1;
      else if( nodetype1 != SCIP_NODETYPE_CHILD && nodetype2 == SCIP_NODETYPE_CHILD )
         return +1;
      else if( nodetype1 == SCIP_NODETYPE_SIBLING && nodetype2 != SCIP_NODETYPE_SIBLING )
         return -1;
      else if( nodetype1 != SCIP_NODETYPE_SIBLING && nodetype2 == SCIP_NODETYPE_SIBLING )
         return +1;
      else
      {
         int depth1;
         int depth2;
         
         depth1 = SCIPnodeGetDepth(node1);
         depth2 = SCIPnodeGetDepth(node2);
         if( depth1 > depth2 )
            return -1;
         else if( depth1 < depth2 )
            return +1;
         else
            return 0;
      }
   }
}





/*
 * bfs specific interface methods
 */

/** creates the node selector for best first search and includes it in SCIP */
RETCODE SCIPincludeNodeselBfs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   NODESELDATA* nodeseldata;

   /* allocate and initialize node selector data; this has to be freed in the destructor */
   CHECK_OKAY( SCIPallocMemory(scip, &nodeseldata) );
   nodeseldata->maxplungequot = MAXPLUNGEQUOT;
   nodeseldata->maxplungedepth = MAXPLUNGEDEPTH;
   nodeseldata->minplungedepth = MINPLUNGEDEPTH;

   /* include node selector */
   CHECK_OKAY( SCIPincludeNodesel(scip, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
                  NODESEL_LOWESTFIRST,
                  nodeselFreeBfs, NULL, NULL, nodeselSelectBfs, nodeselCompBfs,
                  nodeseldata) );

   /* add node selector parameters */
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "nodeselection/bfs/minplungedepth",
                  "minimal plunging depth, before new best node may be selected",
                  &nodeseldata->minplungedepth, MINPLUNGEDEPTH, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "nodeselection/bfs/maxplungedepth",
                  "maximal plunging depth, before new best node is forced to be selected",
                  &nodeseldata->maxplungedepth, MAXPLUNGEDEPTH, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "nodeselection/bfs/maxplungequot",
                  "maximal quotient (actlowerbound - lowerbound)/(avglowerbound - lowerbound) where plunging is performed",
                  &nodeseldata->maxplungequot, MAXPLUNGEQUOT, 0.0, REAL_MAX, NULL, NULL) );
   
   return SCIP_OKAY;
}

