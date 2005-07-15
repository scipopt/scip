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
#pragma ident "@(#) $Id: cutpool.h,v 1.21 2005/07/15 17:20:07 bzfpfend Exp $"

/**@file   cutpool.h
 * @brief  internal methods for storing cuts in a cut pool
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CUTPOOL_H__
#define __SCIP_CUTPOOL_H__


#include "scip/def.h"
#include "scip/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_sepastore.h"
#include "scip/type_cutpool.h"
#include "scip/pub_cutpool.h"



/** creates cut pool */
extern
RETCODE SCIPcutpoolCreate(
   CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   int              agelimit            /**< maximum age a cut can reach before it is deleted from the pool */
   );

/** frees cut pool */
extern
RETCODE SCIPcutpoolFree(
   CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   );

/** removes all rows from the cut pool */
extern
RETCODE SCIPcutpoolClear(
   CUTPOOL*         cutpool,            /**< cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   );

/** if not already existing, adds row to cut pool and captures it */
extern
RETCODE SCIPcutpoolAddRow(
   CUTPOOL*         cutpool,            /**< cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   ROW*             row                 /**< cutting plane to add */
   );

/** adds row to cut pool and captures it; doesn't check for multiple cuts */
extern
RETCODE SCIPcutpoolAddNewRow(
   CUTPOOL*         cutpool,            /**< cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   ROW*             row                 /**< cutting plane to add */
   );

/** removes the LP row from the cut pool */
extern
RETCODE SCIPcutpoolDelRow(
   CUTPOOL*         cutpool,            /**< cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< current LP data */
   ROW*             row                 /**< row to remove */
   );

/** separates cuts of the cut pool */
extern
RETCODE SCIPcutpoolSeparate(
   CUTPOOL*         cutpool,            /**< cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< current LP data */
   SEPASTORE*       sepastore,          /**< separation storage */
   Bool             root,               /**< are we at the root node? */
   RESULT*          result              /**< pointer to store the result of the separation call */
   );


#endif
