/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   var.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for exact data of problem variables
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_VAREX_H__
#define __SCIP_VAREX_H__


#include "blockmemshell/memory.h"
#include "scip/def.h"
#include "scip/rational.h"
#include "scip/type_branch.h"
#include "scip/type_cons.h"
#include "scip/type_event.h"
#include "scip/type_history.h"
#include "scip/type_implics.h"
#include "scip/type_lp.h"
#include "scip/type_message.h"
#include "scip/type_misc.h"
#include "scip/type_primal.h"
#include "scip/type_prob.h"
#include "scip/type_prop.h"
#include "scip/type_relax.h"
#include "scip/type_reopt.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_set.h"
#include "scip/type_sol.h"
#include "scip/type_stat.h"
#include "scip/type_tree.h"
#include "scip/type_var.h"

#ifndef NDEBUG
#include "scip/struct_var.h"
#else
#include "scip/event.h"
#include "scip/pub_history.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** resolves variable to columns and adds them with the coefficient to the 
 *
 */
SCIP_RETCODE SCIPvarAddToRowExact(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_ROWEX*           rowex,              /**< LP row */
   SCIP_Rational*        val                 /**< value of coefficient */
   );

/** converts transformed variable into column variable and creates LP column */
SCIP_RETCODE SCIPvarColumnExact(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LPEX*            lp                  /**< current LP data */
   );

/** changes objective value of variable */
SCIP_RETCODE SCIPvarScaleObjExact(
   SCIP_VAR*             var,                /**< variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             scale               /**< new objective value for variable */
   );

/** changes objective value of variable */
SCIP_RETCODE SCIPvarChgObjExact(
   SCIP_VAR*             var,                /**< variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Rational*        newobj              /**< new objective value for variable */
   );

/** return the status of the exact variable data */
SCIP_VARSTATUS SCIPvarGetStatusExact(
   SCIP_VAR*             var                /**< scip variabel */
   );

/** gets column of COLUMN variable */
SCIP_COLEX* SCIPvarGetColExact(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets primal LP solution value of variable */
SCIP_Rational* SCIPvarGetLPexSolex_rec(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets primal LP solution value of variable */
SCIP_Rational* SCIPvarGetLPexSolex(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets pseudo solution value of variable */
SCIP_Rational* SCIPvarGetPseudoSolex(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets current LP or pseudo solution value of variable */
SCIP_Rational* SCIPvarGetSolex(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             getlpval            /**< should the LP solution value be returned? */
   );

#ifdef __cplusplus
}
#endif

#endif