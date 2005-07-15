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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepastore.h,v 1.22 2005/07/15 17:20:18 bzfpfend Exp $"

/**@file   sepastore.h
 * @brief  internal methods for storing separated cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPASTORE_H__
#define __SCIP_SEPASTORE_H__


#include "scip/def.h"
#include "scip/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_tree.h"
#include "scip/type_sepastore.h"
#include "scip/type_branch.h"



/** creates separation storage */
extern
RETCODE SCIPsepastoreCreate(
   SEPASTORE**      sepastore           /**< pointer to store separation storage */
   );

/** frees separation storage */
extern
RETCODE SCIPsepastoreFree(
   SEPASTORE**      sepastore           /**< pointer to store separation storage */
   );

/** informs separation storage, that the setup of the initial LP starts now */
extern
void SCIPsepastoreStartInitialLP(
   SEPASTORE*       sepastore           /**< separation storage */
   );

/** informs separation storage, that the setup of the initial LP is now finished */
extern
void SCIPsepastoreEndInitialLP(
   SEPASTORE*       sepastore           /**< separation storage */
   );

/** informs separation storage, that the following cuts should be used in any case */
extern
void SCIPsepastoreStartForceCuts(
   SEPASTORE*       sepastore           /**< separation storage */
   );

/** informs separation storage, that the following cuts should no longer be used in any case */
extern
void SCIPsepastoreEndForceCuts(
   SEPASTORE*       sepastore           /**< separation storage */
   );

/** adds cut to separation storage and captures it;
 *  if the cut should be forced to enter the LP, an infinite score has to be used
 */
extern
RETCODE SCIPsepastoreAddCut(
   SEPASTORE*       sepastore,          /**< separation storage */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< LP data */
   ROW*             cut,                /**< separated cut */
   Bool             forcecut,           /**< should the cut be forced to enter the LP? */
   Bool             root                /**< are we at the root node? */
   );

/** adds cuts to the LP and clears separation storage */
extern
RETCODE SCIPsepastoreApplyCuts(
   SEPASTORE*       sepastore,          /**< separation storage */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            cutoff              /**< pointer to store whether an empty domain was created */
   );

/** clears the separation storage without adding the cuts to the LP */
extern
RETCODE SCIPsepastoreClearCuts(
   SEPASTORE*       sepastore,          /**< separation storage */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< LP data */
   );

/** get number of cuts in the separation storage */
extern
int SCIPsepastoreGetNCuts(
   SEPASTORE*       sepastore           /**< separation storage */
   );

/** get total number of cuts found so far */
extern
int SCIPsepastoreGetNCutsFound(
   SEPASTORE*       sepastore           /**< separation storage */
   );

/** get number of cuts found so far in current separation round */
extern
int SCIPsepastoreGetNCutsFoundRound(
   SEPASTORE*            sepastore                /**< separation storage */
   );

/** get total number of cuts stored (and possibly removed again) in current separation round */
extern
int SCIPsepastoreGetNCutsStored(
   SEPASTORE*            sepastore                /**< separation storage */
   );

/** get total number of cuts applied to the LPs */
extern
int SCIPsepastoreGetNCutsApplied(
   SEPASTORE*            sepastore                /**< separation storage */
   );


#endif
