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
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: primal.h,v 1.26 2005/01/18 14:34:29 bzfpfend Exp $"

/**@file   primal.h
 * @brief  internal methods for collecting primal CIP solutions and primal informations
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PRIMAL_H__
#define __PRIMAL_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_event.h"
#include "type_lp.h"
#include "type_var.h"
#include "type_prob.h"
#include "type_sol.h"
#include "type_primal.h"
#include "type_tree.h"
#include "type_heur.h"

#include "struct_primal.h"



/** creates primal data */
extern
RETCODE SCIPprimalCreate(
   PRIMAL**         primal              /**< pointer to primal data */
   );

/** frees primal data */
extern
RETCODE SCIPprimalFree(
   PRIMAL**         primal,             /**< pointer to primal data */
   MEMHDR*          memhdr              /**< block memory */
   );

/** sets upper bound in primal data and in LP solver */
extern
RETCODE SCIPprimalSetUpperbound(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   Real             upperbound          /**< new upper bound */
   );

/** recalculates upper bound in primal data after a change of the problem's objective offset */
extern
RETCODE SCIPprimalUpdateUpperbound(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   );

/** returns whether the current primal bound is justified with a feasible primal solution; if not, the primal bound
 *  was set from the user as objective limit
 */
extern
Bool SCIPprimalUpperboundIsSol(
   PRIMAL*          primal              /**< primal data */
   );

/** adds primal solution to solution storage by copying it */
extern
RETCODE SCIPprimalAddSol(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL*             sol,                /**< primal CIP solution */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** adds primal solution to solution storage, frees the solution afterwards */
extern
RETCODE SCIPprimalAddSolFree(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** adds current LP/pseudo solution to solution storage */
extern
RETCODE SCIPprimalAddCurrentSol(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   HEUR*            heur,               /**< heuristic that found the solution (or NULL if it's from the tree) */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** checks primal solution; if feasible, adds it to storage by copying it */
extern
RETCODE SCIPprimalTrySol(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   );

/** checks primal solution; if feasible, adds it to storage; solution is freed afterwards */
extern
RETCODE SCIPprimalTrySolFree(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether solution was feasible and good enough to keep */
   );

/** checks current LP/pseudo solution; if feasible, adds it to storage */
extern
RETCODE SCIPprimalTryCurrentSol(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   HEUR*            heur,               /**< heuristic that found the solution (or NULL if it's from the tree) */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** inserts solution into the global array of all existing primal solutions */
extern
RETCODE SCIPprimalSolCreated(
   PRIMAL*          primal,             /**< primal data */
   SET*             set,                /**< global SCIP settings */
   SOL*             sol                 /**< primal CIP solution */
   );

/** removes solution from the global array of all existing primal solutions */
extern
void SCIPprimalSolFreed(
   PRIMAL*          primal,             /**< primal data */
   SOL*             sol                 /**< primal CIP solution */
   );

/** updates all existing primal solutions after a change in a variable's objective value */
extern
void SCIPprimalUpdateVarObj(
   PRIMAL*          primal,             /**< primal data */
   VAR*             var,                /**< problem variable */
   Real             oldobj,             /**< old objective value */
   Real             newobj              /**< new objective value */
   );


#endif
