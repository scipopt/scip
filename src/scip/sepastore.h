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

/**@file   sepastore.h
 * @brief  methods and datastructures for storing separated cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SEPASTORE_H__
#define __SEPASTORE_H__


typedef struct SepaStore SEPASTORE;     /**< storage for separated variables */


#include "def.h"
#include "retcode.h"
#include "set.h"
#include "mem.h"
#include "lp.h"
#include "tree.h"


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

/** adds cut to separation storage and captures it */
extern
RETCODE SCIPsepastoreAddCut(
   SEPASTORE*       sepastore,          /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** clears the separation storage without adding the cuts to the LP */
extern
RETCODE SCIPsepastoreClearCuts(
   SEPASTORE*       sepastore,          /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
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
