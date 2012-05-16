/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   primalex.h
 * @brief  internal methods for collecting exact primal CIP solutions and exact primal information
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRIMALEX_H__
#define __SCIP_PRIMALEX_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_sol.h"
#include "scip/type_solex.h"
#include "scip/type_primalex.h"

#include "scip/struct_primalex.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates exact primal data */
extern
SCIP_RETCODE SCIPprimalexCreate(
   SCIP_PRIMALEX**       primal              /**< pointer to exact primal data */
   );


/** frees exact primal data */
extern
SCIP_RETCODE SCIPprimalexFree(
   SCIP_PRIMALEX**       primal,             /**< pointer to exact primal data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** adds exact primal solution to solution storage, frees the solution afterwards */
extern
SCIP_RETCODE SCIPprimalexAddSolFree(
   SCIP_PRIMALEX*        primal,             /**< exact primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_SOLEX**          sol,                /**< pointer to exact primal CIP solution; is cleared in function call */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

#ifdef __cplusplus
}
#endif

#endif
