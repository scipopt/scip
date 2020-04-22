/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   solex.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for storing primal CIP solutions
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SOLEX_H__
#define __SCIP_SOLEX_H__


#include <stdio.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_lpex.h"
#include "scip/type_nlp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_sol.h"
#include "scip/type_primal.h"
#include "scip/type_tree.h"
#include "scip/type_heur.h"
#include "scip/rational.h"
#include "scip/struct_sol.h"

#ifdef __cplusplus
extern "C" {
#endif


/** creates primal CIP solution, initialized to zero */
SCIP_RETCODE SCIPsolexCreate(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to the current LP solution */
SCIP_RETCODE SCIPsolexCreateLPexSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to the current solution */
SCIP_RETCODE SCIPsolexCreateCurrentSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a copy of a primal CIP solution */
SCIP_RETCODE SCIPvalsexCopy(
   SCIP_VALSEX**         valsex,             /**< pointer to store the copy of the primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VALSEX*          sourcevals          /**< primal CIP solution to copy */
   );

/** frees primal CIP solution */
SCIP_RETCODE SCIPvalsexFree(
   SCIP_VALSEX**         valsex,             /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** clears primal CIP solution */
SCIP_RETCODE SCIPsolexClear(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** stores solution values of variables in solution's own array */
SCIP_RETCODE SCIPsolexUnlink(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< transformed problem data */
   );

/** sets value of variable in primal CIP solution */
SCIP_RETCODE SCIPsolexSetVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_VAR*             var,                /**< variable to add to solution */
   SCIP_Rational*        val                 /**< solution value of variable */
   );

/** overwrite FP solution with exact values */
SCIP_RETCODE SCIPsolexOverwriteFPSol(
   SCIP_SOL*             sol,                /**< exact primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< problem data */
   SCIP_PROB*            transprob,          /**< problem data */
   SCIP_TREE*            tree                /**< branch and bound tree, or NULL */
   );

/** returns value of variable in primal CIP solution */
void SCIPsolexGetVal(
   SCIP_Rational*        res,
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var                 /**< variable to get value for */
   );

/** gets objective value of primal CIP solution in transformed problem */
SCIP_Rational* SCIPsolexGetObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob            /**< original problem data */
   );

/** gets objective value of primal CIP solution which lives in the original problem space */
SCIP_Rational* SCIPsolexGetOrigObj(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** returns whether the given solutions are equal */
SCIP_Bool SCIPsolexsAreEqual(
   SCIP_SOL*             sol1,               /**< first primal CIP solution */
   SCIP_SOL*             sol2,               /**< second primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob           /**< transformed problem after presolve, or NULL if both solution are
                                              *   defined in the original problem space */
   );

/** outputs non-zero elements of solution to file stream */
SCIP_RETCODE SCIPsolexPrint(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< problem data (original or transformed) */
   SCIP_PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             mipstart,           /**< should only discrete variables be printed? */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** copies current exact LP solution into CIP solution by linking */
SCIP_RETCODE SCIPsolexLinkLPexSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_LPEX*            lp                  /**< current LP data */
   );

/** copies current pseudo solution into CIP solution by linking */
SCIP_RETCODE SCIPsolexLinkPseudoSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_LPEX*            lp                  /**< current LP data */
   );

/** hard-set the obj value of a solution  */
void SCIPsolSetObjVal(
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_Real             val                 /**< objective value */
   );

/** returns whether the given solution is defined on original variables */
SCIP_Bool SCIPsolexIsOriginal(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** checks primal CIP solution for feasibility
 *
 *  @note The difference between SCIPsolCheck() and SCIPcheckSolOrig() is that modifiable constraints are handled
 *        differently. There might be some variables which do not have an original counter part (e.g. in
 *        branch-and-price). Therefore, modifiable constraints can not be double-checked in the original space.
 */
SCIP_RETCODE SCIPsolexCheck(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            feasible            /**< stores whether solution is feasible */
   );

/** gets current position of solution in array of existing solutions of primal data */
int SCIPsolexGetPrimalexIndex(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** sets current position of solution in array of existing solutions of primal data */
void SCIPsolexSetPrimalexIndex(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   primalindex         /**< new primal index of solution */
   );

SCIP_Bool SCIPsolIsExactSol(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** retransforms exact part of solution to original problem space */
SCIP_RETCODE SCIPsolexRetransform(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_Bool*            hasinfval           /**< pointer to store whether the solution has infinite values */
   );

#endif