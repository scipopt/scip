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

/**@file   cutpool.h
 * @brief  methods and datastructures for storing cuts in a cut pool
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CUTPOOL_H__
#define __CUTPOOL_H__


typedef struct Cutpool CUTPOOL;         /**< storage for pooled cuts */


#include "def.h"
#include "retcode.h"
#include "set.h"
#include "mem.h"
#include "lp.h"
#include "sepa.h"


extern
RETCODE SCIPcutpoolCreate(              /**< creates cut pool */
   CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   int              agelimit            /**< maximum age a cut can reach before it is deleted from the pool */
   );

extern
RETCODE SCIPcutpoolFree(                /**< frees cut pool */
   CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPcutpoolAddRow(              /**< if not already existing, adds row to cut pool and captures it */
   CUTPOOL*         cutpool,            /**< cut pool */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row                 /**< cutting plane to add */
   );

extern
RETCODE SCIPcutpoolAddNewRow(           /**< adds row to cut pool and captures it; doesn't check for multiple cuts */
   CUTPOOL*         cutpool,            /**< cut pool */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row                 /**< cutting plane to add */
   );

extern
RETCODE SCIPcutpoolDelRow(              /**< removes the LP row from the cut pool */
   CUTPOOL*         cutpool,            /**< cut pool */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< actual LP data */
   ROW*             row                 /**< row to remove */
   );

extern
RETCODE SCIPcutpoolSeparate(            /**< separates cuts of the cut pool */
   CUTPOOL*         cutpool,            /**< cut pool */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< actual LP data */
   SEPA*            sepa,               /**< separation storage */
   Bool             root                /**< are we at the root node? */
   );

extern
int SCIPcutpoolGetNCuts(                /**< get number of cuts in the cut pool */
   CUTPOOL*         cutpool             /**< cut pool */
   );



#endif
