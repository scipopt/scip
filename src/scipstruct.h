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

/**@file   scipstruct.h
 * @brief  SCIP main data structure
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIPSTRUCT_H__
#define __SCIPSTRUCT_H__


#include "set.h"
#include "mem.h"
#include "interrupt.h"
#include "clock.h"
#include "prob.h"
#include "stat.h"
#include "tree.h"
#include "lp.h"
#include "price.h"
#include "sepa.h"
#include "branch.h"
#include "cutpool.h"
#include "conflict.h"
#include "primal.h"
#include "event.h"
#include "dialog.h"


/** SCIP main data structure */
struct Scip
{
   STAGE            stage;              /**< SCIP operation stage */
   SET*             set;                /**< global SCIP settings */
   MEM*             mem;                /**< block memory buffers */
   INTERRUPT*       interrupt;          /**< CTRL-C interrupt data */
   CLOCK*           totaltime;          /**< total SCIP running time */
   PROB*            origprob;           /**< original problem data */
   PROB*            transprob;          /**< transformed problem after presolve */
   STAT*            stat;               /**< dynamic problem statistics */
   TREE*            tree;               /**< branch and bound tree */
   LP*              lp;                 /**< LP data */
   PRICE*           price;              /**< storage for priced variables */
   SEPASTORE*       sepastore;          /**< storage for separated cuts */
   BRANCHCAND*      branchcand;         /**< storage for branching candidates */
   CUTPOOL*         cutpool;            /**< global cut pool */
   CONFLICT*        conflict;           /**< conflict analysis data for propagation conflicts */
   LPCONFLICT*      lpconflict;         /**< conflict analysis data for infeasible LP conflicts */
   PRIMAL*          primal;             /**< primal data and solution storage */
   EVENTFILTER*     eventfilter;        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue;         /**< event queue to cache events and process them later (bound change events) */
   DIALOGHDLR*      dialoghdlr;         /**< dialog handler for user interface */
};


#endif
