/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepastore.h,v 1.11 2004/04/29 15:20:40 bzfpfend Exp $"

/**@file   sepastore.h
 * @brief  internal methods for storing separated cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SEPASTORE_H__
#define __SEPASTORE_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_event.h"
#include "type_lp.h"
#include "type_tree.h"
#include "type_sepastore.h"
#include "type_branch.h"



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

/** adds cut to separation storage and captures it;
 *  if the cut should be forced to enter the LP, an infinite score has to be used
 */
extern
RETCODE SCIPsepastoreAddCut(
   SEPASTORE*       sepastore,          /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   ROW*             cut,                /**< separated cut */
   Real             score,              /**< separation score of cut (the larger, the better the cut) */
   Bool             root                /**< are we at the root node? */
   );

/** adds cuts to the LP and clears separation storage */
extern
RETCODE SCIPsepastoreApplyCuts(
   SEPASTORE*       sepastore,          /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
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
   MEMHDR*          memhdr,             /**< block memory */
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

/** get total number of cuts applied to the LPs */
extern
int SCIPsepastoreGetNCutsApplied(
   SEPASTORE*            sepastore                /**< separation storage */
   );


#endif
