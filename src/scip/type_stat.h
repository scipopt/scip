/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   type_stat.h
 * @brief  type definitions for problem statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_STAT_H__
#define __SCIP_TYPE_STAT_H__

#ifdef __cplusplus
extern "C" {
#endif

/** SCIP solving status */
enum SCIP_Status
{
   SCIP_STATUS_UNKNOWN        =  0,     /**< the solving status is not yet known */
   SCIP_STATUS_OPTIMAL        =  1,     /**< the problem was solved to optimality, an optimal solution is available */
   SCIP_STATUS_INFEASIBLE     =  2,     /**< the problem was proven to be infeasible */
   SCIP_STATUS_UNBOUNDED      =  3,     /**< the problem was proven to be unbounded */
   SCIP_STATUS_INFORUNBD      =  4,     /**< the problem was proven to be either infeasible or unbounded */
   SCIP_STATUS_USERINTERRUPT  = 10,     /**< the user interrupted the solving process (SIGINT or SCIPinterruptSolve()) */
   SCIP_STATUS_TERMINATE      = 11,     /**< the process received a SIGTERM signal */
   SCIP_STATUS_NODELIMIT      = 20,     /**< the solving process was interrupted because the node limit was reached */
   SCIP_STATUS_TOTALNODELIMIT = 21,     /**< the solving process was interrupted because the total node limit was
                                         *   reached (incl. restarts) */
   SCIP_STATUS_STALLNODELIMIT = 22,     /**< the solving process was interrupted because the stalling node limit was
                                         *   reached (no improvement w.r.t. primal bound) */
   SCIP_STATUS_TIMELIMIT      = 23,     /**< the solving process was interrupted because the time limit was reached */
   SCIP_STATUS_MEMLIMIT       = 24,     /**< the solving process was interrupted because the memory limit was reached */
   SCIP_STATUS_GAPLIMIT       = 25,     /**< the solving process was interrupted because the gap limit was reached */
   SCIP_STATUS_PRIMALLIMIT    = 26,     /**< the solving process was interrupted because the primal limit was reached */
   SCIP_STATUS_DUALLIMIT      = 27,     /**< the solving process was interrupted because the dual limit was reached */
   SCIP_STATUS_SOLLIMIT       = 28,     /**< the solving process was interrupted because the solution limit was reached */
   SCIP_STATUS_BESTSOLLIMIT   = 29,     /**< the solving process was interrupted because the solution improvement limit
                                         *   was reached */
   SCIP_STATUS_RESTARTLIMIT   = 30,     /**< the solving process was interrupted because the restart limit was reached */
};
typedef enum SCIP_Status SCIP_STATUS;

typedef struct SCIP_Stat SCIP_STAT;               /**< problem and runtime specific statistics */

#ifdef __cplusplus
}
#endif

#endif
