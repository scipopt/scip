/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
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


extern
RETCODE SCIPprimalCreate(               /**< creates primal data */
   PRIMAL**         primal,             /**< pointer to primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPprimalFree(                 /**< frees primal data */
   PRIMAL**         primal,             /**< pointer to primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPprimalAddSol(               /**< adds solution to primal solution storage */
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   SOL**            sol                 /**< pointer to primal CIP solution */
   );

#endif
