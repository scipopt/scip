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
#include "sort.h"
#include "stat.h"
#include "var.h"
#include "heur.h"
#include "prob.h"
#include "tree.h"



extern
RETCODE SCIPsolCreate(                  /**< creates primal CIP solution, initialized to zero */
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

extern
RETCODE SCIPsolCopy(                    /**< creates a copy of a primal CIP solution */
   SOL**            sol,                /**< pointer to store the copy of the primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   SOL*             sourcesol           /**< primal CIP solution to copy */
   );

extern
RETCODE SCIPsolCreateLPSol(             /**< creates primal CIP solution, initialized to the actual LP solution */
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< actual LP data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

extern
RETCODE SCIPsolCreatePseudoSol(         /**< creates primal CIP solution, initialized to the actual pseudo solution */
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

extern
RETCODE SCIPsolCreateActSol(            /**< creates primal CIP solution, initialized to the actual solution */
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

extern
RETCODE SCIPsolFree(                    /**< frees primal CIP solution */
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr              /**< block memory */
   );

extern
RETCODE SCIPsolLinkLPSol(               /**< copies actual LP solution into CIP solution by linking */
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPsolLinkPseudoSol(           /**< copies actual pseudo solution into CIP solution by linking */
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
RETCODE SCIPsolLinkActSol(              /**< copies actual solution (LP or pseudo solution) into CIP solution by linking */
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPsolClear(                   /**< clears primal CIP solution */
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat                /**< problem statistics data */
   );

extern
RETCODE SCIPsolUnlink(                  /**< stores solution values of variables in solution's own array */
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob                /**< problem data */
   );

extern
RETCODE SCIPsolSetVal(                  /**< sets value of variable in primal CIP solution */
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   VAR*             var,                /**< variable to add to solution */
   Real             val                 /**< solution value of variable */
   );

extern
RETCODE SCIPsolIncVal(                  /**< increases value of variable in primal CIP solution */
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   VAR*             var,                /**< variable to increase solution value for */
   Real             incval              /**< increment for solution value of variable */
   );

extern
RETCODE SCIPsolGetVal(                  /**< returns value of variable in primal CIP solution */
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   VAR*             var,                /**< variable to get value for */
   Real*            solval              /**< pointer to store the solution value */
   );

extern
RETCODE SCIPsolCheck(                   /**< checks primal CIP solution for feasibility */
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
   Bool*            feasible            /**< stores whether solution is feasible */
   );

extern
Real SCIPsolGetObj(                     /**< gets objective value of primal CIP solution */
   SOL*             sol                 /**< primal CIP solution */
   );

extern
Longint SCIPsolGetNodenum(              /**< gets node number, where this solution was found */
   SOL*             sol                 /**< primal CIP solution */
   );

extern
HEUR* SCIPsolGetHeur(                   /**< gets heuristic, that found this solution (or NULL if it's from the tree) */
   SOL*             sol                 /**< primal CIP solution */
   );

extern
RETCODE SCIPsolPrint(                   /**< outputs non-zero elements of solution to file stream */
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< problem data */
   FILE*            file                /**< output file (or NULL for standard output) */
   );


#endif
