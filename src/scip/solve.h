/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: solve.h,v 1.32 2005/01/31 12:21:02 bzfpfend Exp $"

/**@file   solve.h
 * @brief  internal methods for main solving loop and node processing
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SOLVE_H__
#define __SOLVE_H__


#include <stdio.h>

#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_event.h"
#include "type_lp.h"
#include "type_prob.h"
#include "type_primal.h"
#include "type_tree.h"
#include "type_pricestore.h"
#include "type_sepastore.h"
#include "type_cutpool.h"
#include "type_conflict.h"



/** returns whether the solving process will be / was stopped before proving optimality;
 *  if the solving process was stopped, stores the reason as status in stat
 */
extern
Bool SCIPsolveIsStopped(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** applies domain propagation on current node */
extern
RETCODE SCIPpropagateDomains(
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   int              depth,              /**< depth level to use for propagator frequency checks */
   int              maxrounds,          /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   );

/** main solving loop */
extern
RETCODE SCIPsolveCIP(
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   MEM*             mem,                /**< block memory pools */
   PROB*            prob,               /**< transformed problem after presolve */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICESTORE*      pricestore,         /**< pricing storage */
   SEPASTORE*       sepastore,          /**< separation storage */
   CUTPOOL*         cutpool,            /**< global cut pool */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   CONFLICT*        conflict,           /**< conflict analysis data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            restart             /**< should solving process be started again with presolving? */
   );


#endif
