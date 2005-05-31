/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_sol.h,v 1.8 2005/05/31 17:20:19 bzfpfend Exp $"

/**@file   pub_sol.h
 * @brief  public methods for primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_SOL_H__
#define __PUB_SOL_H__


#include "scip/def.h"
#include "scip/type_sol.h"
#include "scip/type_heur.h"

#ifdef NDEBUG
#include "scip/struct_sol.h"
#endif


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets origin of solution */
extern
SOLORIGIN SCIPsolGetOrigin(
   SOL*             sol                 /**< primal CIP solution */
   );

/** gets clock time, when this solution was found */
extern
Real SCIPsolGetTime(
   SOL*             sol                 /**< primal CIP solution */
   );

/** gets branch and bound run number, where this solution was found */
extern
int SCIPsolGetRunnum(
   SOL*             sol                 /**< primal CIP solution */
   );

/** gets node number of the specific branch and bound run, where this solution was found */
extern
Longint SCIPsolGetNodenum(
   SOL*             sol                 /**< primal CIP solution */
   );

/** gets node's depth, where this solution was found */
extern
int SCIPsolGetDepth(
   SOL*             sol                 /**< primal CIP solution */
   );

/** gets heuristic, that found this solution (or NULL if it's from the tree) */
extern
HEUR* SCIPsolGetHeur(
   SOL*             sol                 /**< primal CIP solution */
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

#endif

#endif
