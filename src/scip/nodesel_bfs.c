/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nodesel_bfs.c
 * @brief  node selector for best first search
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "nodesel_bfs.h"


#define NODESEL_NAME "bfs"
#define NODESEL_DESC "Best First Search"


/*
 * Default parameter settings
 */

#define SCIP_DEFAULT_MAXPLUNGEDEPTH  10 /**< maximal plunging depth, before new best node is forced to be selected */



/** node selector data for best first search node selection */
struct NodeselData
{
   int              maxplungedepth;     /**< maximal plunging depth, before new best node is forced to be selected */
};



/*
 * Callback methods
 */

static
DECL_NODESELFREE(SCIPnodeselFreeBfs)
{
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
DECL_NODESELSELECT(SCIPnodeselSelectBfs)
{
   NODESELDATA* nodeseldata;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   *selnode = NULL;

   /* get node selector user data */
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   /* check, if we want to plunge once more */
   if( SCIPgetPlungeDepth(scip) < nodeseldata->maxplungedepth )
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
DECL_NODESELCOMP(SCIPnodeselCompBfs)
{
   Real lowerbound1;
   Real lowerbound2;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   lowerbound1 = SCIPnodeGetLowerBound(node1);
   lowerbound2 = SCIPnodeGetLowerBound(node2);
   if( SCIPisLT(scip, lowerbound1, lowerbound2) )
      return -1;
   else if( SCIPisGT(scip, lowerbound1, lowerbound2) )
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

   /* allocate and initialise node selector data; this has to be freed in the destructor */
   CHECK_OKAY( SCIPallocMemory(scip, &nodeseldata) );
   nodeseldata->maxplungedepth = SCIP_DEFAULT_MAXPLUNGEDEPTH;

   /* include node selector */
   CHECK_OKAY( SCIPincludeNodesel(scip, NODESEL_NAME, NODESEL_DESC,
                  SCIPnodeselFreeBfs, NULL, NULL, SCIPnodeselSelectBfs, SCIPnodeselCompBfs,
                  nodeseldata, FALSE) );

   return SCIP_OKAY;
}

