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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: primal.h,v 1.31 2005/02/14 13:35:47 bzfpfend Exp $"

/**@file   primal.h
 * @brief  internal methods for collecting primal CIP solutions and primal informations
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PRIMAL_H__
#define __PRIMAL_H__


#include "scip/def.h"
#include "scip/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_sol.h"
#include "scip/type_primal.h"
#include "scip/type_tree.h"
#include "scip/type_heur.h"

#include "scip/struct_primal.h"



/** creates primal data */
extern
RETCODE SCIPprimalCreate(
   PRIMAL**         primal              /**< pointer to primal data */
   );

/** frees primal data */
extern
RETCODE SCIPprimalFree(
   PRIMAL**         primal,             /**< pointer to primal data */
   BLKMEM*          blkmem              /**< block memory */
   );

/** sets the cutoff bound in primal data and in LP solver */
extern
RETCODE SCIPprimalSetCutoffbound(
   PRIMAL*          primal,             /**< primal data */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   Real             cutoffbound         /**< new cutoff bound */
   );

/** sets upper bound in primal data and in LP solver */
extern
RETCODE SCIPprimalSetUpperbound(
   PRIMAL*          primal,             /**< primal data */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   Real             upperbound          /**< new upper bound */
   );

/** updates upper bound and cutoff bound in primal data after a tightening of the problem's objective limit */
extern
RETCODE SCIPprimalUpdateObjlimit(
   PRIMAL*          primal,             /**< primal data */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   );

/** recalculates upper bound and cutoff bound in primal data after a change of the problem's objective offset */
extern
RETCODE SCIPprimalUpdateObjoffset(
   PRIMAL*          primal,             /**< primal data */
   BLKMEM*          blkmem,             /**< block memory */
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
   PRIMAL*          primal,              /**< primal data */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob                /**< transformed problem after presolve */
   );

/** adds primal solution to solution storage by copying it */
extern
RETCODE SCIPprimalAddSol(
   PRIMAL*          primal,             /**< primal data */
   BLKMEM*          blkmem,             /**< block memory */
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
   BLKMEM*          blkmem,             /**< block memory */
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
   BLKMEM*          blkmem,             /**< block memory */
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
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkbounds,        /**< should the bounds of the variables be checked? */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   );

/** checks primal solution; if feasible, adds it to storage; solution is freed afterwards */
extern
RETCODE SCIPprimalTrySolFree(
   PRIMAL*          primal,             /**< primal data */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool             checkbounds,        /**< should the bounds of the variables be checked? */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether solution was feasible and good enough to keep */
   );

/** checks current LP/pseudo solution; if feasible, adds it to storage */
extern
RETCODE SCIPprimalTryCurrentSol(
   PRIMAL*          primal,             /**< primal data */
   BLKMEM*          blkmem,             /**< block memory */
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

/** retransforms all existing solutions to original problem space */
extern
RETCODE SCIPprimalRetransformSolutions(
   PRIMAL*          primal,             /**< primal data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            origprob            /**< original problem */
   );


#endif
