/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lpex.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for exact LP management
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_LPEX_H__
#define __SCIP_LPEX_H__


#include <stdio.h>

#include "scip/def.h"
#include "scip/rational.h"
#include "blockmemshell/memory.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_misc.h"
#include "scip/type_lpex.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_sol.h"
#include "scip/pub_lp.h"

#include "scip/struct_lpex.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * column methods
 */

/*
 * Column methods
 */

/** checks if the exact column and its fpcol are consistent */
extern
SCIP_Bool colexInSync(
   SCIP_COLEX*           colex,              /**< exact column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     msg                 /**< message handler */
   );

/** creates an LP column */
extern
SCIP_RETCODE SCIPcolexCreate(
   SCIP_COLEX**          col,                /**< pointer to column data */
   SCIP_COL*             fpcol,              /**< the corresponding fp col */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             var,                /**< variable, this column represents */
   int                   len,                /**< number of nonzeros in the column */
   SCIP_ROWEX**          rows,               /**< array with rows of column entries */
   SCIP_Rational**       vals,               /**< array with coefficients of column entries */
   SCIP_Bool             removable           /**< should the column be removed from the LP due to aging or cleanup? */
   );

/** frees an LP column */
extern
SCIP_RETCODE SCIPcolexFree(
   SCIP_COLEX**          col,                /**< pointer to LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp                  /**< current LP data */
   );

/** output column to file stream */
extern
void SCIPcolexPrint(
   SCIP_COLEX*           col,                /**< LP column */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** adds a previously non existing coefficient to an LP column */
extern
SCIP_RETCODE SCIPcolexAddCoef(
   SCIP_COLEX*           col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_Rational*        val                 /**< value of coefficient */
   );

/** deletes coefficient from column */
extern
SCIP_RETCODE SCIPcolexDelCoef(
   SCIP_COLEX*           col,                /**< column to be changed */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_ROWEX*           row                 /**< coefficient to be deleted */
   );

/** changes or adds a coefficient to an LP column */
extern
SCIP_RETCODE SCIPcolexChgCoef(
   SCIP_COLEX*           col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_Rational*        val                 /**< value of coefficient */
   );

/** increases value of an existing or nonexisting coefficient in an LP column */
extern
SCIP_RETCODE SCIPcolexIncCoef(
   SCIP_COLEX*           col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_Rational*        incval              /**< value to add to the coefficient */
   );

/** changes objective value of column */
extern
SCIP_RETCODE SCIPcolexChgObj(
   SCIP_COLEX*           col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_Rational*        newobj              /**< new objective value */
   );

/** changes lower bound of column */
extern
SCIP_RETCODE SCIPcolexChgLb(
   SCIP_COLEX*           col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_Rational*        newlb               /**< new lower bound value */
   );

/** changes upper bound of column */
extern
SCIP_RETCODE SCIPcolexChgUb(
   SCIP_COLEX*           col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_Rational*        newub               /**< new upper bound value */
   );

/*
 * row methods
 */

/** sorts row entries such that LP columns precede non-LP columns and inside both parts lower column indices precede
 *  higher ones
 */
void SCIProwexSort(
   SCIP_ROWEX*           row                 /**< row to be sorted */
   );

/** increases usage counter of LP row */
void SCIProwexCapture(
   SCIP_ROWEX*           row                 /**< LP row */
   );

/** output column to file stream */
extern
void SCIProwexPrint(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** get the index of an exact row */
int SCIProwexGetIndex(
   SCIP_ROWEX*           row                 /**< LP row */
   );

/** get the length of a row */
int SCIProwexGetNNonz(
   SCIP_ROWEX*           row                 /**< LP row */
   );

/** returns TRUE iff row is member of current LP */
extern
SCIP_Bool SCIProwexIsInLP(
   SCIP_ROWEX*           row                 /**< LP row */
   );

/** return TRUE iff row is modifiable */
SCIP_Bool SCIProwexIsModifiable(
   SCIP_ROWEX*           row                 /**< LP row */
   );

/** returns true, if an exact row for this fprow was already created */
extern
SCIP_Bool SCIProwHasExRow(
   SCIP_LPEX*            lpex,               /**< exact lp data structure */
   SCIP_ROW*             row                 /**< SCIP row */
   );

/** returns exact row corresponding to fprow, if it exists. Otherwise returns NULL */
extern
SCIP_ROWEX* SCIProwGetExRow(
   SCIP_LPEX*            lpex,               /**< exact lp data structure */
   SCIP_ROW*             row                 /**< SCIP row */
   );

/** returns exact col corresponding to fpcol, if it exists. Otherwise returns NULL */
extern
SCIP_COLEX* SCIPcolGetExCol(
   SCIP_LPEX*            lpex,               /**< exact lp data structure */
   SCIP_COL*             col                 /**< SCIP col */
   );

/** calculates the Farkas coefficient or reduced cost of a column i using the given dual Farkas vector y */
void SCIPcolexCalcFarkasRedcostCoef(
   SCIP_COLEX*           col,                /**< LP column */
   SCIP_SET*             set,                /**< SCIP settings pointer */
   SCIP_Rational*        result,             /**< rational to store the result */
   SCIP_Rational**       dual,               /**< dense dual Farkas vector, NULL to use internal row-values */
   SCIP_Bool             usefarkas           /**< should the farkas coefficient be computed ? */
   );

/** creates and captures an LP row */
extern
SCIP_RETCODE SCIProwCreateExact(
   SCIP_ROWEX**          row,                /**< pointer to LP row data */
   SCIP_ROW*             fprow,              /**< corresponding fp row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEX*            lpex,               /**< current LP data */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COLEX**          cols,               /**< array with columns of row entries */
   SCIP_Rational**       vals,               /**< array with coefficients of row entries */
   SCIP_Rational*        lhs,                /**< left hand side of row */
   SCIP_Rational*        rhs,                /**< right hand side of row */
   SCIP_ROWORIGINTYPE    origintype,         /**< type of origin of row */
   void*                 origin              /**< pointer to constraint handler or separator who created the row (NULL if unkown) */
   );

/** applies all cached changes to the LP solver */
extern
SCIP_RETCODE SCIPlpexFlush(
   SCIP_LPEX*            lp,                 /**< current exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   );

/*
 * lp methods 
 */

/** returns whether it is possible to use neumair-shcherbina bounding method */
extern
SCIP_Bool SCIPlpexBSpossible(
   SCIP_LPEX*            lp                  /**< pointer to LP data object */
   );

/** returns whether it is possible to use project and shift bounding method */
extern
SCIP_Bool SCIPlpexPSpossible(
   SCIP_LPEX*            lp                  /**< pointer to LP data object */
   );

/** checks that lp and fplp are properly synced */
extern
SCIP_Bool SCIPlpexIsSynced(
   SCIP_LPEX*            lp,                 /**< pointer to LP data object */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     msg                 /**< message handler */
   );

/** creates empty LP data object */
extern
SCIP_RETCODE SCIPlpexCreate(
   SCIP_LPEX**           lpex,               /**< pointer to LP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LP*              fplp,               /**< the normal floating point lp */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name                /**< problem name */
   );

/** frees LP data object */
SCIP_RETCODE SCIPlpexFree(
   SCIP_LPEX**           lp,                 /**< pointer to LP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** adds a column to the LP and captures the variable */
extern
SCIP_RETCODE SCIPlpexAddCol(
   SCIP_LPEX*            lpex,               /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COLEX*           colex,              /**< LP column */
   int                   depth               /**< depth in the tree where the column addition is performed */
   );

/** adds a row to the LP and captures it */
extern
SCIP_RETCODE SCIPlpexAddRow(
   SCIP_LPEX*            lpex,               /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_ROWEX*           rowex,              /**< LP row */
   int                   depth               /**< depth in the tree where the row addition is performed */
   );

/** returns the feasibility of a row for the given solution */
extern
void SCIProwexGetSolFeasibility(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Rational*        result              /**< result pointer */
   );

/** returns the activity of a row for a given solution */
extern
void SCIProwexGetSolActivity(
   SCIP_ROWEX*           rowex,              /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SOLEX*           solex,              /**< exact primal CIP solution */
   SCIP_Bool             useexact,           /**< should the exact solution be used */
   SCIP_Rational*        result              /**< resulting activity */
   );

/** decreases usage counter of LP row, and frees memory if necessary */
extern
SCIP_RETCODE SCIProwexRelease(
   SCIP_ROWEX**          row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEX*            lp                  /**< current LP data */
   );

/** frees an LP row */
extern
SCIP_RETCODE SCIProwexFree(
   SCIP_ROWEX**          row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEX*            lp                  /**< current LP data */
   );

/** ensuresr, that column array of row can store at least num entries */
extern
SCIP_RETCODE SCIProwexEnsureSize(
   SCIP_ROWEX*           row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** add constant value to a row */
extern
SCIP_RETCODE SCIProwexAddConstant(
   SCIP_ROWEX*           row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_Rational*        addval              /**< constant value to add to the row */
   );

/** increases value of an existing or nonexisting coefficient in an LP column */
extern
SCIP_RETCODE SCIProwexIncCoef(
   SCIP_ROWEX*           row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_COLEX*           col,                /**< LP column */
   SCIP_Rational*        incval              /**< valpelue to add to the coefficient */
   );

/** returns the feasibility of a row in the current LP solution: negative value means infeasibility */
extern
void SCIProwexGetLPFeasibility(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_Rational*        result              /**< rational pointer to store the result */
   );

/** returns the pseudo feasibility of a row in the current pseudo solution: negative value means infeasibility */
extern
void SCIProwexGetPseudoFeasibility(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Rational*        result              /**< rational pointer to store the result */
   );

/** returns the activity of a row in the current LP solution */
extern
SCIP_Rational* SCIProwexGetLPActivity(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEX*            lp                  /**< current LP data */
   );

/** returns the pseudo activity of a row in the current pseudo solution */
extern
SCIP_Rational* SCIProwexGetPseudoActivity(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** recalculates the current activity of a row */
extern
void SCIProwexRecalcLPActivity(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

 /** calculates the current pseudo activity of a row */
extern
void SCIProwexRecalcPseudoActivity(
   SCIP_ROWEX*           row,                /**< row data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** gets objective value of column */
extern
SCIP_Rational* SCIPcolexGetObj(
   SCIP_COLEX*           col                 /**< LP column */
   );

/** gets lower bound of column */
extern
SCIP_Rational* SCIPcolexGetLb(
   SCIP_COLEX*           col                 /**< LP column */
   );

/** gets upper bound of column */
extern
SCIP_Rational* SCIPcolexGetUb(
   SCIP_COLEX*           col                 /**< LP column */
   );

/** gets best bound of column with respect to the objective function */
extern
SCIP_Rational* SCIPcolexGetBestBound(
   SCIP_COLEX*           col                 /**< LP column */
   );

/** gets the primal LP solution of a column */
extern
SCIP_Rational* SCIPcolexGetPrimsol(
   SCIP_COLEX*           col                 /**< LP column */
   );

/** gets the minimal LP solution value, this column ever assumed */
extern
SCIP_Rational* SCIPcolexGetMinPrimsol(
   SCIP_COLEX*           col                 /**< LP column */
   );

/** gets the maximal LP solution value, this column ever assumed */
extern
SCIP_Rational* SCIPcolexGetMaxPrimsol(
   SCIP_COLEX*           col                 /**< LP column */
   );


/*
 * lp update methods
 */



/** updates current pseudo and loose objective value for a change in a variable's objective value */
extern
SCIP_RETCODE SCIPlpexUpdateVarObj(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldobj,             /**< old objective value of variable */
   SCIP_Rational*        newobj              /**< new objective value of variable */
   );

/** updates current root pseudo objective value for a global change in a variable's lower bound */
extern
SCIP_RETCODE SCIPlpexUpdateVarLbGlobal(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldlb,              /**< old lower bound of variable */
   SCIP_Rational*        newlb               /**< new lower bound of variable */
   );

/** updates current pseudo and loose objective value for a change in a variable's lower bound */
extern
SCIP_RETCODE SCIPlpexUpdateVarLb(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldlb,              /**< old lower bound of variable */
   SCIP_Rational*        newlb               /**< new lower bound of variable */
   );

/** updates current root pseudo objective value for a global change in a variable's upper bound */
extern
SCIP_RETCODE SCIPlpexUpdateVarUbGlobal(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldub,              /**< old upper bound of variable */
   SCIP_Rational*        newub               /**< new upper bound of variable */
   );

/** updates current pseudo objective value for a change in a variable's upper bound */
extern
SCIP_RETCODE SCIPlpexUpdateVarUb(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldub,              /**< old upper bound of variable */
   SCIP_Rational*        newub               /**< new upper bound of variable */
   );

/** informs LP, that given variable was added to the problem */
extern
SCIP_RETCODE SCIPlpexUpdateAddVar(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that is now a LOOSE problem variable */
   );

/** informs LP, that given variable is to be deleted from the problem */
extern
SCIP_RETCODE SCIPlpexUpdateDelVar(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that will be deleted from the problem */
   );

/** informs LP, that given formerly loose problem variable is now a column variable */
extern
SCIP_RETCODE SCIPlpexUpdateVarColumn(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   );

/** informs LP, that given formerly column problem variable is now again a loose variable */
extern
SCIP_RETCODE SCIPlpexUpdateVarLoose(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   );

/** decrease the number of loose variables by one */
extern
void SCIPlpexDecNLoosevars(
   SCIP_LPEX*            lp                  /**< current LP data */
   );

extern
SCIP_RETCODE SCIPlexGetNRows(
   SCIP_LPEX*            lp                  /**< current LP data */
   );

/** stores the LP solution in the columns and rows */
extern
SCIP_RETCODE SCIPlpexGetSol(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            dualfeasible        /**< pointer to store whether the solution is dual feasible, or NULL */
   );

/** stores LP solution with infinite objective value in the columns and rows */
extern
SCIP_RETCODE SCIPlpexGetUnboundedSol(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            rayfeasible         /**< pointer to store whether the primal ray is a feasible unboundedness proof, or NULL */
   );

/** returns primal ray proving the unboundedness of the current LP */
extern
SCIP_RETCODE SCIPlpexGetPrimalRay(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Rational**       ray                 /**< array for storing primal ray values, they are stored w.r.t. the problem index of the variables,
                                              *   so the size of this array should be at least number of active variables
                                              *   (all entries have to be initialized to 0 before) */
   );

/** stores the dual Farkas multipliers for infeasibility proof in rows. besides, the proof is checked for validity if
 *  lp/checkfarkas = TRUE.
 *
 *  @note the check will not be performed if @p valid is NULL.
 */
extern
SCIP_RETCODE SCIPlpexGetDualfarkas(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            valid               /**< pointer to store whether the Farkas proof is valid  or NULL */
   );

/** gets objective value of current LP
 *
 *  @note This method returns the objective value of the current LP solution, which might be primal or dual infeasible
 *        if a limit was hit during solving. It must not be used as a dual bound if the LP solution status is
 *        SCIP_LPSOLSTAT_ITERLIMIT or SCIP_LPSOLSTAT_TIMELIMIT.
 */
void SCIPlpexGetObjval(
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Rational*        res                 /**< result pointer to store rational */
   );

/** gets the pseudo objective value for the current search node; that is all variables set to their best (w.r.t. the
 *  objective function) local bound
 */
void SCIPlpexGetPseudoObjval(
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Rational*        res                 /**< result pointer to store rational */
   );

/** removes all columns after the given number of cols from the LP */
extern
SCIP_RETCODE SCIPlpexShrinkCols(
   SCIP_LPEX*            lpex,               /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   newncols            /**< new number of columns in the LP */
   );

/** removes and releases all rows after the given number of rows from the LP */
extern
SCIP_RETCODE SCIPlpexShrinkRows(
   SCIP_LPEX*            lpex,               /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   int                   newnrows            /**< new number of rows in the LP */
   );

/** resets the LP to the empty LP by removing all columns and rows from LP, releasing all rows, and flushing the
 *  changes to the LP solver
 */
SCIP_RETCODE SCIPlpexReset(
   SCIP_LPEX*            lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** removes all columns and rows from LP, releases all rows */
extern
SCIP_RETCODE SCIPlpexClear(
   SCIP_LPEX*            lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** checks whether primal solution satisfies all integrality restrictions exactly. 
 * This checks either the fp solution exactly or checks the exact solution, if one exists.
 */
SCIP_RETCODE SCIPlpexEnfoIntegralityExact(
   SCIP_LP* lp,
   SCIP_LPEX* lpex,
   SCIP_SET* set,
   SCIP_STAT* stat,
   SCIP_RESULT* result
   );

#ifdef __cplusplus
}
#endif

#endif
