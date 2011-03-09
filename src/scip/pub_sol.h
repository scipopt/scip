/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_sol.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for primal CIP solutions
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_SOL_H__
#define __SCIP_PUB_SOL_H__


#include "scip/def.h"
#include "scip/type_sol.h"
#include "scip/type_heur.h"

#ifdef NDEBUG
#include "scip/struct_sol.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets origin of solution */
extern
SCIP_SOLORIGIN SCIPsolGetOrigin(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets clock time, when this solution was found */
extern
SCIP_Real SCIPsolGetTime(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets branch and bound run number, where this solution was found */
extern
int SCIPsolGetRunnum(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets node number of the specific branch and bound run, where this solution was found */
extern
SCIP_Longint SCIPsolGetNodenum(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets node's depth, where this solution was found */
extern
int SCIPsolGetDepth(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets heuristic, that found this solution (or NULL if it's from the tree) */
extern
SCIP_HEUR* SCIPsolGetHeur(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** informs the solution that it now belongs to the given primal heuristic */
extern
void SCIPsolSetHeur(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** returns unique index of given solution */
extern
int SCIPsolGetIndex(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );


#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPsolGetOrigin(sol)           ((sol)->solorigin)
#define SCIPsolGetTime(sol)             (sol)->time
#define SCIPsolGetNodenum(sol)          (sol)->nodenum
#define SCIPsolGetRunnum(sol)           (sol)->runnum
#define SCIPsolGetDepth(sol)            (sol)->depth
#define SCIPsolGetHeur(sol)             (sol)->heur
#define SCIPsolGetIndex(sol)            (sol)->index
#define SCIPsolSetHeur(sol,heur)        { (sol)->heur = heur; }
#endif

#ifdef __cplusplus
}
#endif

#endif
