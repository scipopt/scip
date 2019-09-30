/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   primalex.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for collecting exact primal CIP solutions and primal informations
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRIMALEX_H__
#define __SCIP_PRIMALEX_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_lpex.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_sol.h"
#include "scip/type_primal.h"
#include "scip/type_tree.h"
#include "scip/type_reopt.h"
#include "scip/type_heur.h"

#include "scip/struct_primal.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * exact methods
 */

/** creates exact primal data */
extern
SCIP_RETCODE SCIPprimalexCreate(
   SCIP_PRIMALEX**       primal,             /**< pointer to exact primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PRIMAL*          fpstorage           /**< fp storage */
   );


/** frees exact primal data */
extern
SCIP_RETCODE SCIPprimalexFree(
   SCIP_PRIMALEX**       primal,             /**< pointer to exact primal data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** adds exact primal solution to solution storage, frees the solution afterwards */
extern
SCIP_RETCODE SCIPprimalexTrySolFree(
   SCIP_PRIMALEX*        primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_SOLEX**          sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool             printreason,        /**< Should all the reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** clears primal data */
SCIP_RETCODE SCIPprimalexClear(
   SCIP_PRIMALEX**       primal,             /**< pointer to primal data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** inserts solution into the global array of all existing primal solutions */
extern
SCIP_RETCODE SCIPprimalexSolexCreated(
   SCIP_PRIMALEX*        primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOLEX*           sol                 /**< primal CIP solution */
   );

/** removes solution from the global array of all existing primal solutions */
extern
void SCIPprimalexSolexFreed(
   SCIP_PRIMALEX*          primal,             /**< primal data */
   SCIP_SOLEX*             sol                 /**< primal CIP solution */
   );

#ifdef __cplusplus
}
#endif

#endif