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

/**@file   sol.h
 * @brief  datastructures and methods for storing primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SOL_H__
#define __SOL_H__


/** origin of solution: where to retrieve uncached elements */
enum SolOrigin
{
   SCIP_SOLORIGIN_ZERO      = 0,        /**< all non-cached elements in solution are equal to zero */
   SCIP_SOLORIGIN_LPSOL     = 1,        /**< all non-cached elements in solution are equal to actual LP solution */
   SCIP_SOLORIGIN_PSEUDOSOL = 2         /**< all non-cached elements in solution are equal to actual pseudo solution */
};
typedef enum SolOrigin SOLORIGIN;

typedef struct Sol SOL;                 /**< primal CIP solution */


#include "def.h"
#include "retcode.h"
#include "memory.h"
#include "misc.h"
#include "stat.h"
#include "var.h"
#include "heur.h"
#include "prob.h"
#include "tree.h"



/** creates primal CIP solution, initialized to zero */
extern
RETCODE SCIPsolCreate(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a copy of a primal CIP solution */
extern
RETCODE SCIPsolCopy(
   SOL**            sol,                /**< pointer to store the copy of the primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   SOL*             sourcesol           /**< primal CIP solution to copy */
   );

/** creates primal CIP solution, initialized to the actual LP solution */
extern
RETCODE SCIPsolCreateLPSol(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< actual LP data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to the actual pseudo solution */
extern
RETCODE SCIPsolCreatePseudoSol(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to the actual solution */
extern
RETCODE SCIPsolCreateActSol(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** frees primal CIP solution */
extern
RETCODE SCIPsolFree(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr              /**< block memory */
   );

/** copies actual LP solution into CIP solution by linking */
extern
RETCODE SCIPsolLinkLPSol(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< actual LP data */
   );

/** copies actual pseudo solution into CIP solution by linking */
extern
RETCODE SCIPsolLinkPseudoSol(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree                /**< branch-and-bound tree */
   );

/** copies actual solution (LP or pseudo solution) into CIP solution by linking */
extern
RETCODE SCIPsolLinkActSol(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   );

/** clears primal CIP solution */
extern
RETCODE SCIPsolClear(
   SOL*             sol,                /**< primal CIP solution */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree                /**< branch-and-bound tree */
   );

/** stores solution values of variables in solution's own array */
extern
RETCODE SCIPsolUnlink(
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob                /**< problem data */
   );

/** sets value of variable in primal CIP solution */
extern
RETCODE SCIPsolSetVal(
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   VAR*             var,                /**< variable to add to solution */
   Real             val                 /**< solution value of variable */
   );

/** increases value of variable in primal CIP solution */
extern
RETCODE SCIPsolIncVal(
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   VAR*             var,                /**< variable to increase solution value for */
   Real             incval              /**< increment for solution value of variable */
   );

/** returns value of variable in primal CIP solution */
extern
RETCODE SCIPsolGetVal(
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   VAR*             var,                /**< variable to get value for */
   Real*            solval              /**< pointer to store the solution value */
   );

/** checks primal CIP solution for feasibility */
extern
RETCODE SCIPsolCheck(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
   Bool*            feasible            /**< stores whether solution is feasible */
   );

/** gets objective value of primal CIP solution */
extern
Real SCIPsolGetObj(
   SOL*             sol                 /**< primal CIP solution */
   );

/** gets clock time, when this solution was found */
extern
Real SCIPsolGetTime(
   SOL*             sol                 /**< primal CIP solution */
   );

/** gets node number, where this solution was found */
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

/** outputs non-zero elements of solution to file stream */
extern
RETCODE SCIPsolPrint(
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< problem data */
   FILE*            file                /**< output file (or NULL for standard output) */
   );


#endif
