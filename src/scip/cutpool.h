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
#pragma ident "@(#) $Id: cutpool.h,v 1.12 2004/02/05 14:12:35 bzfpfend Exp $"

/**@file   cutpool.h
 * @brief  internal methods for storing cuts in a cut pool
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CUTPOOL_H__
#define __CUTPOOL_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_result.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_lp.h"
#include "type_sepastore.h"
#include "type_cutpool.h"
#include "pub_cutpool.h"



/** creates cut pool */
extern
RETCODE SCIPcutpoolCreate(
   CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   int              agelimit            /**< maximum age a cut can reach before it is deleted from the pool */
   );

/** frees cut pool */
extern
RETCODE SCIPcutpoolFree(
   CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   );

/** if not already existing, adds row to cut pool and captures it */
extern
RETCODE SCIPcutpoolAddRow(
   CUTPOOL*         cutpool,            /**< cut pool */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row                 /**< cutting plane to add */
   );

/** adds row to cut pool and captures it; doesn't check for multiple cuts */
extern
RETCODE SCIPcutpoolAddNewRow(
   CUTPOOL*         cutpool,            /**< cut pool */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row                 /**< cutting plane to add */
   );

/** removes the LP row from the cut pool */
extern
RETCODE SCIPcutpoolDelRow(
   CUTPOOL*         cutpool,            /**< cut pool */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< current LP data */
   ROW*             row                 /**< row to remove */
   );

/** separates cuts of the cut pool */
extern
RETCODE SCIPcutpoolSeparate(
   CUTPOOL*         cutpool,            /**< cut pool */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< current LP data */
   SEPASTORE*       sepastore,          /**< separation storage */
   Bool             root,               /**< are we at the root node? */
   RESULT*          result              /**< pointer to store the result of the separation call */
   );


#endif
