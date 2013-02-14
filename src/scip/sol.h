/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sol.h
 * @brief  internal methods for storing primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SOL_H__
#define __SCIP_SOL_H__


#include <stdio.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_nlp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_sol.h"
#include "scip/type_primal.h"
#include "scip/type_tree.h"
#include "scip/type_heur.h"
#include "scip/pub_sol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates primal CIP solution, initialized to zero */
extern
SCIP_RETCODE SCIPsolCreate(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution in original problem space, initialized to zero */
extern
SCIP_RETCODE SCIPsolCreateOriginal(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a copy of a primal CIP solution */
extern
SCIP_RETCODE SCIPsolCopy(
   SCIP_SOL**            sol,                /**< pointer to store the copy of the primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SOL*             sourcesol           /**< primal CIP solution to copy */
   );

/** transformes given original solution to the transformed space; a corresponding transformed solution has to be given
 *  which is copied into the existing solution and freed afterwards
 */
extern
SCIP_RETCODE SCIPsolTransform(
   SCIP_SOL*             sol,                /**< primal CIP solution to change, living in original space */
   SCIP_SOL**            transsol,           /**< pointer to corresponding transformed primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal              /**< primal data */
   );

/** creates primal CIP solution, initialized to the current LP solution */
extern
SCIP_RETCODE SCIPsolCreateLPSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to the current NLP solution */
extern
SCIP_RETCODE SCIPsolCreateNLPSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to the current relaxation solution */
extern
SCIP_RETCODE SCIPsolCreateRelaxSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to the current pseudo solution */
extern
SCIP_RETCODE SCIPsolCreatePseudoSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to the current solution */
extern
SCIP_RETCODE SCIPsolCreateCurrentSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to unknown values */
extern
SCIP_RETCODE SCIPsolCreateUnknown(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** frees primal CIP solution */
extern
SCIP_RETCODE SCIPsolFree(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PRIMAL*          primal              /**< primal data */
   );

/** copies current LP solution into CIP solution by linking */
extern
SCIP_RETCODE SCIPsolLinkLPSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** copies current NLP solution into CIP solution by linking */
extern
SCIP_RETCODE SCIPsolLinkNLPSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** copies current relaxation solution into CIP solution by linking */
extern
SCIP_RETCODE SCIPsolLinkRelaxSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_RELAXATION*      relaxation          /**< global relaxation data */
   );

/** copies current pseudo solution into CIP solution by linking */
extern
SCIP_RETCODE SCIPsolLinkPseudoSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** copies current solution (LP or pseudo solution) into CIP solution by linking */
extern
SCIP_RETCODE SCIPsolLinkCurrentSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** clears primal CIP solution */
extern
SCIP_RETCODE SCIPsolClear(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** declares all entries in the primal CIP solution to be unknown */
extern
SCIP_RETCODE SCIPsolSetUnknown(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** stores solution values of variables in solution's own array */
extern
SCIP_RETCODE SCIPsolUnlink(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< transformed problem data */
   );

/** sets value of variable in primal CIP solution */
extern
SCIP_RETCODE SCIPsolSetVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_VAR*             var,                /**< variable to add to solution */
   SCIP_Real             val                 /**< solution value of variable */
   );

/** increases value of variable in primal CIP solution */
extern
SCIP_RETCODE SCIPsolIncVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_VAR*             var,                /**< variable to increase solution value for */
   SCIP_Real             incval              /**< increment for solution value of variable */
   );

/** returns value of variable in primal CIP solution */
extern
SCIP_Real SCIPsolGetVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var                 /**< variable to get value for */
   );

/** returns value of variable in primal ray represented by primal CIP solution */
extern
SCIP_Real SCIPsolGetRayVal(
   SCIP_SOL*             sol,                /**< primal CIP solution, representing a primal ray */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var                 /**< variable to get value for */
   );


/** gets objective value of primal CIP solution in transformed problem */
extern
SCIP_Real SCIPsolGetObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< transformed problem data */
   );

/** updates primal solutions after a change in a variable's objective value */
extern
void SCIPsolUpdateVarObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             oldobj,             /**< old objective value */
   SCIP_Real             newobj              /**< new objective value */
   );

/** checks primal CIP solution for feasibility */
extern
SCIP_RETCODE SCIPsolCheck(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_Bool             printreason,        /**< should all reasons of violations be printed? */
   SCIP_Bool             checkbounds,        /**< should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool*            feasible            /**< stores whether solution is feasible */
   );

/** try to round given solution */
extern
SCIP_RETCODE SCIPsolRound(
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool*            success             /**< pointer to store whether rounding was successful */
   );

/** updates the solution value sums in variables by adding the value in the given solution */
extern
void SCIPsolUpdateVarsum(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_Real             weight              /**< weight of solution in weighted average */
   );

/** retransforms solution to original problem space */
extern
SCIP_RETCODE SCIPsolRetransform(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob            /**< original problem */
   );

/** recomputes the objective value of an original solution, e.g., when transferring solutions
 *  from the solution pool (objective coefficients might have changed in the meantime)
 */
extern
void SCIPsolRecomputeObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob            /**< original problem */
   );


/** returns whether the given solutions in transformed space are equal */
extern
SCIP_Bool SCIPsolsAreEqual(
   SCIP_SOL*             sol1,               /**< first primal CIP solution */
   SCIP_SOL*             sol2,               /**< second primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob           /**< transformed problem after presolve */
   );

/** outputs non-zero elements of solution to file stream */
extern
SCIP_RETCODE SCIPsolPrint(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< problem data (original or transformed) */
   SCIP_PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** outputs non-zero elements of solution representing a ray to file stream */
extern
SCIP_RETCODE SCIPsolPrintRay(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< problem data (original or transformed) */
   SCIP_PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** gets current position of solution in array of existing solutions of primal data */
extern
int SCIPsolGetPrimalIndex(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** sets current position of solution in array of existing solutions of primal data */
extern
void SCIPsolSetPrimalIndex(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   primalindex         /**< new primal index of solution */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPsolGetPrimalIndex(sol)      ((sol)->primalindex)
#define SCIPsolSetPrimalIndex(sol,idx)  { (sol)->primalindex = idx; }

#endif

#ifdef __cplusplus
}
#endif

#endif
