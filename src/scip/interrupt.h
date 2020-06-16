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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   interrupt.h
 * @ingroup INTERNALAPI
 * @brief  methods for catching the user CTRL-C interrupt
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_INTERRUPT_H__
#define __SCIP_INTERRUPT_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_interrupt.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates a CTRL-C interrupt data */
SCIP_RETCODE SCIPinterruptCreate(
   SCIP_INTERRUPT**      interrupt           /**< pointer to store the CTRL-C interrupt data */
   );

/** frees a CTRL-C interrupt data */
void SCIPinterruptFree(
   SCIP_INTERRUPT**      interrupt           /**< pointer to the CTRL-C interrupt data */
   );

/** captures the CTRL-C interrupt to call the SCIP's own interrupt handler */
void SCIPinterruptCapture(
   SCIP_INTERRUPT*       interrupt           /**< CTRL-C interrupt data */
   );

/** releases the CTRL-C interrupt and restores the old interrupt handler */
void SCIPinterruptRelease(
   SCIP_INTERRUPT*       interrupt           /**< CTRL-C interrupt data */
   );

/** returns whether the user interrupted by pressing CTRL-C */
SCIP_Bool SCIPinterrupted(
   void
   );

/** returns whether the process has received a SIGTERM */
SCIP_Bool SCIPterminated(
   void
   );

/** send a termination signal to the process so that SCIP tries to terminate as soon as possible */
void SCIPtryTerminate(
   void
   );

/** resets the number of interrupts to 0 */
void SCIPresetInterrupted(
   void
   );

#ifdef __cplusplus
}
#endif

#endif
