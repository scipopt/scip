/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nodesel_dfs.h
 * @brief  node selector for depth first search
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __NODESEL_DFS_H__
#define __NODESEL_DFS_H__


#include "scip.h"


DECL_NODESELINIT(SCIPnodeselInitDfs);
DECL_NODESELEXIT(SCIPnodeselExitDfs);
DECL_NODESELSLCT(SCIPnodeselSlctDfs);
DECL_NODESELCOMP(SCIPnodeselCompDfs);


extern
RETCODE SCIPincludeNodeselDfs(          /**< creates the node selector for depth first search and includes it in SCIP */
   SCIP*            scip                /**< SCIP data structure */
   );

#endif
