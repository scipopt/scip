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

/**@file   primal.h
 * @brief  datastructures and methods for collecting primal CIP solutions and primal informations
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PRIMAL_H__
#define __PRIMAL_H__


typedef struct Primal PRIMAL;           /**< primal data */


#include "def.h"
#include "retcode.h"
#include "memory.h"
#include "set.h"
#include "var.h"
#include "lp.h"


/** primal data and solution storage */
struct Primal
{
   SOL**            sols;               /**< primal CIP solutions */
   int              solssize;           /**< size of sols array */
   int              nsols;              /**< number of primal CIP solutions stored in sols array */
   int              nsolsfound;         /**< number of primal CIP solutions found up to now */
   Real             upperbound;         /**< upper (primal) bound of CIP: objective value of best solution or user bound */
};


/** creates primal data */
extern
RETCODE SCIPprimalCreate(
   PRIMAL**         primal,             /**< pointer to primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   LP*              lp                  /**< actual LP data */
   );

/** frees primal data */
extern
RETCODE SCIPprimalFree(
   PRIMAL**         primal,             /**< pointer to primal data */
   MEMHDR*          memhdr              /**< block memory */
   );

/** adds primal solution to solution storage by moving it */
extern
RETCODE SCIPprimalAddSolMove(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   SOL**            sol                 /**< pointer to primal CIP solution; is cleared in function call */
   );

/** adds primal solution to solution storage by copying it */
extern
RETCODE SCIPprimalAddSolCopy(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   SOL*             sol                 /**< primal CIP solution */
   );

/** checks primal solution; if feasible, adds it to storage by moving it */
extern
RETCODE SCIPprimalTrySolMove(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   );

/** checks primal solution; if feasible, adds it to storage by copying it */
extern
RETCODE SCIPprimalTrySolCopy(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   SOL*             sol,                /**< primal CIP solution */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   );

#endif
