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

/**@file   nodesel_dfs.c
 * @brief  node selector for depth first search
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "nodesel_dfs.h"


#define NODESEL_NAME "dfs"
#define NODESEL_DESC "Depth First Search"




/*
 * Callback methods
 */

static
DECL_NODESELSELECT(SCIPnodeselSelectDfs)
{
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   CHECK_OKAY( SCIPgetBestNode(scip, selnode) );

   return SCIP_OKAY;
}

static
DECL_NODESELCOMP(SCIPnodeselCompDfs)
{
   int depth1;
   int depth2;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   depth1 = SCIPnodeGetDepth(node1);
   depth2 = SCIPnodeGetDepth(node2);
   if( depth1 > depth2 )
      return -1;
   else if( depth1 < depth2 )
      return +1;
   else
   {
      Real lowerbound1;
      Real lowerbound2;

      lowerbound1 = SCIPnodeGetLowerBound(node1);
      lowerbound2 = SCIPnodeGetLowerBound(node2);
      if( lowerbound1 < lowerbound2 )
         return -1;
      else if( lowerbound1 > lowerbound2 )
         return +1;
      else
         return 0;
   }
}





/*
 * dfs specific interface methods
 */

/** creates the node selector for depth first search and includes it in SCIP */
RETCODE SCIPincludeNodeselDfs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( SCIPincludeNodesel(scip, NODESEL_NAME, NODESEL_DESC,
                  NULL, NULL, NULL, SCIPnodeselSelectDfs, SCIPnodeselCompDfs,
                  NULL, FALSE) );

   return SCIP_OKAY;
}

