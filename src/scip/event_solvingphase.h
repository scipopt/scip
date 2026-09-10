/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_solvingphase.h
 * @ingroup DEFPLUGINS_EVENT
 * @brief  eventhdlr for solving phase dependent parameter adjustment
 * @author Gregor Hendel
 *
 * this event handler is used to apply dynamic parameter adjustment depending on the
 * progress of the solving process.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_SOLVINGPHASE_H__
#define __SCIP_EVENT_SOLVINGPHASE_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** enumerator to represent the current solving phase */
enum SCIP_SolvingPhase
{
   SCIP_SOLVINGPHASE_UNINITIALIZED = -1,     /**< solving phase has not been initialized yet */
   SCIP_SOLVINGPHASE_FEASIBILITY   =  0,     /**< no solution was found until now */
   SCIP_SOLVINGPHASE_IMPROVEMENT   =  1,     /**< current incumbent solution is suboptimal */
   SCIP_SOLVINGPHASE_PROOF         =  2      /**< current incumbent is optimal */
};
typedef enum SCIP_SolvingPhase SCIP_SOLVINGPHASE;

/** transition criteria flags reached by the solvingphase event handler */
#define SCIP_SOLVINGPHASEFLAG_NONE           UINT8_C(0x00)  /**< no transition criterion has been reached */
#define SCIP_SOLVINGPHASEFLAG_RANK1          UINT8_C(0x01)  /**< rank-1 node based transition criterion has been reached */
#define SCIP_SOLVINGPHASEFLAG_ESTIMATE       UINT8_C(0x02)  /**< best estimate transition criterion has been reached */
#define SCIP_SOLVINGPHASEFLAG_OPTIMAL        UINT8_C(0x04)  /**< optimal value transition criterion has been reached */
#define SCIP_SOLVINGPHASEFLAG_LOG            UINT8_C(0x08)  /**< logarithmic regression transition criterion has been reached */

typedef uint8_t SCIP_SOLVINGPHASEFLAG;                     /**< flag for the solving phases (bit field) */

/** creates event handler for solving phase event */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeEventHdlrSolvingphase(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the current solving phase tracked by the solvingphase event handler
 *
 *  The phase advances from SCIP_SOLVINGPHASE_FEASIBILITY to SCIP_SOLVINGPHASE_IMPROVEMENT when a solution is found,
 *  and to SCIP_SOLVINGPHASE_PROOF when a heuristic transition criterion (see parameter
 *  `solvingphases/transitionmethod`) declares that the incumbent is expected to be optimal and
 *  the solver is proving optimality.
 *
 *  The phase is tracked when the event handler is active, i.e. when `solvingphases/enabled`
 *  or `solvingphases/testmode` is set to TRUE. Otherwise, SCIP_SOLVINGPHASE_UNINITIALIZED is returned.
 */
SCIP_EXPORT
SCIP_SOLVINGPHASE SCIPgetSolvingPhase(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the bit field of reached transition criteria
 *
 *  See SCIP_SOLVINGPHASEFLAG_* for the available flags. The flags are tracked whenever the
 *  solvingphase event handler is active, i.e. when `solvingphases/enabled` or
 *  `solvingphases/testmode` is set to TRUE. Otherwise, SCIP_SOLVINGPHASEFLAG_NONE is returned.
 */
SCIP_EXPORT
SCIP_SOLVINGPHASEFLAG SCIPgetSolvingPhaseFlags(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
