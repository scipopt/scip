/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: lp.h,v 1.75 2004/04/29 15:20:38 bzfpfend Exp $"

/**@file   lp.h
 * @brief  internal methods for LP management
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __LP_H__
#define __LP_H__


#include <stdio.h>

#include "def.h"
#include "memory.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_lpi.h"
#include "type_misc.h"
#include "type_lp.h"
#include "type_var.h"
#include "type_prob.h"
#include "type_sol.h"
#include "pub_lp.h"

#include "struct_lp.h"



/*
 * Column methods
 */

/** creates an LP column */
extern
RETCODE SCIPcolCreate(
   COL**            col,                /**< pointer to column data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR*             var,                /**< variable, this column represents */
   int              len,                /**< number of nonzeros in the column */
   ROW**            row,                /**< array with rows of column entries */
   Real*            val,                /**< array with coefficients of column entries */
   Bool             removeable          /**< should the column be removed from the LP due to aging or cleanup? */
   );

/** frees an LP column */
extern
RETCODE SCIPcolFree(
   COL**            col,                /**< pointer to LP column */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   );

/** sorts column entries by row index */
extern
void SCIPcolSort(
   COL*             col                 /**< column to be sorted */
   );

/** adds a previously non existing coefficient to an LP column */
extern
RETCODE SCIPcolAddCoeff(
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   );

/** deletes coefficient from column */
extern
RETCODE SCIPcolDelCoeff(
   COL*             col,                /**< column to be changed */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   ROW*             row                 /**< coefficient to be deleted */
   );

/** changes or adds a coefficient to an LP column */
extern
RETCODE SCIPcolChgCoeff(
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   );

/** increases value of an existing or nonexisting coefficient in an LP column */
extern
RETCODE SCIPcolIncCoeff(
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   ROW*             row,                /**< LP row */
   Real             incval              /**< value to add to the coefficient */
   );

/** changes objective value of column */
extern
RETCODE SCIPcolChgObj(
   COL*             col,                /**< LP column to change */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             newobj              /**< new objective value */
   );

/** changes lower bound of column */
extern
RETCODE SCIPcolChgLb(
   COL*             col,                /**< LP column to change */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             newlb               /**< new lower bound value */
   );

/** changes upper bound of column */
extern
RETCODE SCIPcolChgUb(
   COL*             col,                /**< LP column to change */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             newub               /**< new upper bound value */
   );

/** gets the reduced costs of a column in last LP or after recalculation */
extern
Real SCIPcolGetRedcost(
   COL*             col,                /**< LP column */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< current LP data */
   );

/** gets the feasibility of (the dual row of) a column in last LP or after recalculation */
extern
Real SCIPcolGetFeasibility(
   COL*             col,                /**< LP column */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< current LP data */
   );

/** gets the farkas value of a column in last LP (which must be infeasible) */
extern
Real SCIPcolGetFarkas(
   COL*             col,                /**< LP column */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< current LP data */
   );

/** gets strong branching information on a column variable */
extern
RETCODE SCIPcolGetStrongbranch(
   COL*             col,                /**< LP column */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   LP*              lp,                 /**< LP data */
   int              itlim,              /**< iteration limit for strong branchings */
   Real*            down,               /**< stores dual bound after branching column down */
   Real*            up,                 /**< stores dual bound after branching column up */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   );

/** gets last strong branching information available for a column variable;
 *  returns values of SCIP_INVALID, if strong branching was not yet called on the given column;
 *  keep in mind, that the returned old values may have nothing to do with the current LP solution
 */
extern
void SCIPcolGetStrongbranchLast(
   COL*             col,                /**< LP column */
   Real*            down,               /**< stores dual bound after branching column down, or NULL */
   Real*            up,                 /**< stores dual bound after branching column up, or NULL */
   Real*            solval              /**< stores LP solution value of column at last strong branching call, or NULL */
   );




/*
 * Row methods
 */

/** creates and captures an LP row */
extern
RETCODE SCIProwCreate(
   ROW**            row,                /**< pointer to LP row data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             local,              /**< is row only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   );

/** frees an LP row */
extern
RETCODE SCIProwFree(
   ROW**            row,                /**< pointer to LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   );

/** ensures, that column array of row can store at least num entries */
extern
RETCODE SCIProwEnsureSize(
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

/** increases usage counter of LP row */
extern
void SCIProwCapture(
   ROW*             row                 /**< LP row */
   );

/** decreases usage counter of LP row, and frees memory if necessary */
extern
RETCODE SCIProwRelease(
   ROW**            row,                /**< pointer to LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   );

/** sorts row entries by column index */
extern
void SCIProwSort(
   ROW*             row                 /**< row to be sorted */
   );

/** enables delaying of row sorting */
extern
void SCIProwDelaySort(
   ROW*             row                 /**< LP row */
   );

/** disables delaying of row sorting, sorts row and merges coefficients with equal columns */
extern
void SCIProwForceSort(
   ROW*             row,                /**< LP row */
   SET*             set                 /**< global SCIP settings */
   );

/** adds a previously non existing coefficient to an LP row */
extern
RETCODE SCIProwAddCoeff(
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   COL*             col,                /**< LP column */
   Real             val                 /**< value of coefficient */
   );

/** deletes coefficient from row */
extern
RETCODE SCIProwDelCoeff(
   ROW*             row,                /**< row to be changed */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   COL*             col                 /**< coefficient to be deleted */
   );

/** changes or adds a coefficient to an LP row */
extern
RETCODE SCIProwChgCoeff(
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   COL*             col,                /**< LP column */
   Real             val                 /**< value of coefficient */
   );

/** increases value of an existing or nonexisting coefficient in an LP column */
extern
RETCODE SCIProwIncCoeff(
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   COL*             col,                /**< LP column */
   Real             incval              /**< value to add to the coefficient */
   );

/** changes constant value of a row */
extern
RETCODE SCIProwChgConstant(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   Real             constant            /**< new constant value */
   );

/** add constant value to a row */
extern
RETCODE SCIProwAddConstant(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   Real             addval              /**< constant value to add to the row */
   );

/** changes left hand side of LP row */
extern
RETCODE SCIProwChgLhs(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             lhs                 /**< new left hand side */
   );

/** changes right hand side of LP row */
extern
RETCODE SCIProwChgRhs(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             rhs                 /**< new right hand side */
   );

/** tries to find a rational representation of the row and multiplies coefficients with common denominator */
extern
RETCODE SCIProwMakeRational(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   Real             maxscale,           /**< maximal value to scale row with */
   Bool*            success             /**< stores whether row could be made rational */
   );

/** returns the activity of a row in the current LP solution */
extern
Real SCIProwGetLPActivity(
   ROW*             row,                /**< LP row */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< current LP data */
   );

/** returns the feasibility of a row in the current LP solution: negative value means infeasibility */
extern
Real SCIProwGetLPFeasibility(
   ROW*             row,                /**< LP row */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< current LP data */
   );

/** returns the pseudo activity of a row in the current pseudo solution */
extern
Real SCIProwGetPseudoActivity(
   ROW*             row,                /**< LP row */
   STAT*            stat                /**< problem statistics */
   );

/** returns the pseudo feasibility of a row in the current pseudo solution: negative value means infeasibility */
extern
Real SCIProwGetPseudoFeasibility(
   ROW*             row,                /**< LP row */
   STAT*            stat                /**< problem statistics */
   );

/** returns the activity of a row for a given solution */
extern
RETCODE SCIProwGetSolActivity(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   SOL*             sol,                /**< primal CIP solution */
   Real*            solactivity         /**< pointer to store the row's activity for the solution */
   );

/** returns the feasibility of a row for the given solution */
extern
RETCODE SCIProwGetSolFeasibility(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   SOL*             sol,                /**< primal CIP solution */
   Real*            solfeasibility      /**< pointer to store the row's feasibility for the solution */
   );

/** returns the minimal activity of a row w.r.t. the column's bounds */
extern
Real SCIProwGetMinActivity(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics data */
   );

/** returns the maximal activity of a row w.r.t. the column's bounds */
extern
Real SCIProwGetMaxActivity(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics data */
   );

/** gets maximal absolute value of row vector coefficients */
extern
Real SCIProwGetMaxval(
   ROW*             row,                /**< LP row */
   SET*             set                 /**< global SCIP settings */
   );

/** gets minimal absolute value of row vector's non-zero coefficients */
extern
Real SCIProwGetMinval(
   ROW*             row,                /**< LP row */
   SET*             set                 /**< global SCIP settings */
   );




/*
 * LP methods
 */

/** creates empty LP data object */
extern
RETCODE SCIPlpCreate(
   LP**             lp,                 /**< pointer to LP data object */
   SET*             set,                /**< global SCIP settings */
   const char*      name                /**< problem name */
   );

/** frees LP data object */
extern
RETCODE SCIPlpFree(
   LP**             lp,                 /**< pointer to LP data object */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   );

/** adds a column to the LP and captures the variable */
extern
RETCODE SCIPlpAddCol(
   LP*              lp,                 /**< LP data */
   SET*             set,                /**< global SCIP settings */
   COL*             col                 /**< LP column */
   );

/** adds a row to the LP and captures it */
extern
RETCODE SCIPlpAddRow(
   LP*              lp,                 /**< LP data */
   SET*             set,                /**< global SCIP settings */
   ROW*             row                 /**< LP row */
   );

/** removes all columns after the given number of columns from the LP */
extern
RETCODE SCIPlpShrinkCols(
   LP*              lp,                 /**< LP data */
   int              newncols            /**< new number of columns in the LP */
   );

/** removes and releases all rows after the given number of rows from the LP */
extern
RETCODE SCIPlpShrinkRows(
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              newnrows            /**< new number of rows in the LP */
   );

/** removes all columns and rows from LP, releases all rows */
extern
RETCODE SCIPlpClear(
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   );

/** remembers number of columns and rows to track the newly added ones */
extern
void SCIPlpMarkSize(
   LP*              lp                  /**< current LP data */
   );

/** gets all indices of basic columns and rows: index i >= 0 corresponds to column i, index i < 0 to row -i-1 */
extern
RETCODE SCIPlpGetBasisInd(
   LP*              lp,                 /**< LP data */
   int*             basisind            /**< pointer to store the basis indices */
   );

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
extern
RETCODE SCIPlpGetBase(
   LP*              lp,                 /**< LP data */
   int*             cstat,              /**< array to store column basis status, or NULL */
   int*             rstat               /**< array to store row basis status, or NULL */
   );

/** gets a row from the inverse basis matrix B^-1 */
extern
RETCODE SCIPlpGetBInvRow(
   LP*              lp,                 /**< LP data */
   int              r,                  /**< row number */
   Real*            coef                /**< pointer to store the coefficients of the row */
   );

/** gets a row from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A) */
extern
RETCODE SCIPlpGetBInvARow(
   LP*              lp,                 /**< LP data */
   int              r,                  /**< row number */
   Real*            binvrow,            /**< row in B^-1 from prior call to SCIPlpGetBInvRow(), or NULL */
   Real*            coef                /**< pointer to store the coefficients of the row */
   );

/** calculates a weighted sum of all LP rows; for negative weights, the left and right hand side of the corresponding
 *  LP row are swapped in the summation
 */
extern
RETCODE SCIPlpSumRows(
   LP*              lp,                 /**< LP data */
   SET*             set,                /**< global SCIP settings */
   int              nvars,              /**< number of active variables in the problem */
   Real*            weights,            /**< row weights in row summation */
   REALARRAY*       sumcoef,            /**< array to store sum coefficients indexed by variables' probindex */
   Real*            sumlhs,             /**< pointer to store the left hand side of the row summation */
   Real*            sumrhs              /**< pointer to store the right hand side of the row summation */
   );

/* calculates a MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a MIR cut.
 */
extern
RETCODE SCIPlpCalcMIR(
   LP*              lp,                 /**< LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   int              nvars,              /**< number of active variables in the problem */
   VAR**            vars,               /**< active variables in the problem */
   Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   Real             scale,              /**< additional scaling factor multiplied to all rows */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   Bool*            success             /**< pointer to store whether the returned coefficients are a valid MIR cut */
   );

/** stores LP state (like basis information) into LP state object */
extern
RETCODE SCIPlpGetState(
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   );

/** loads LP state (like basis information) into solver */
extern
RETCODE SCIPlpSetState(
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LPISTATE*        lpistate            /**< LP state information (like basis information) */
   );

/** sets the upper objective limit of the LP solver */
extern
RETCODE SCIPlpSetCutoffbound(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   Real             cutoffbound         /**< new upper objective limit */
   );

/** solves the LP with the primal or dual simplex algorithm, depending on the current basis feasibility */
extern
RETCODE SCIPlpSolve(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   Bool             fastmip,            /**< should the FASTMIP setting of the LP solver be activated? */
   Bool             fromscratch,        /**< should the LP be solved from scratch without using current basis? */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   );

/** solves the LP with simplex algorithm, and copy the solution into the column's data */
extern
RETCODE SCIPlpSolveAndEval(
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   int              itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   Bool             aging,              /**< should aging and removal of obsolete cols/rows be applied? */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   );

/** gets solution status of last solve call */
extern
LPSOLSTAT SCIPlpGetSolstat(
   LP*              lp                  /**< current LP data */
   );

/** gets objective value of last solution */
extern
Real SCIPlpGetObjval(
   LP*              lp,                 /**< current LP data */
   SET*             set                 /**< global SCIP settings */
   );

/** gets part of objective value of last solution that results from COLUMN variables only */
extern
Real SCIPlpGetColumnObjval(
   LP*              lp                  /**< current LP data */
   );

/** gets part of objective value of last solution that results from LOOSE variables only */
extern
Real SCIPlpGetLooseObjval(
   LP*              lp,                 /**< current LP data */
   SET*             set                 /**< global SCIP settings */
   );

/** gets current pseudo objective value */
extern
Real SCIPlpGetPseudoObjval(
   LP*              lp,                 /**< current LP data */
   SET*             set                 /**< global SCIP settings */
   );

/** gets pseudo objective value, if a bound of the given variable would be modified in the given way */
extern
Real SCIPlpGetModifiedPseudoObjval(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable */
   Real             oldbound,           /**< old value for bound */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** gets pseudo objective value, if a bound of the given variable would be modified in the given way;
 *  perform calculations with interval arithmetic to get an exact lower bound
 */
extern
Real SCIPlpGetModifiedProvedPseudoObjval(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable */
   Real             oldbound,           /**< old value for bound */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** updates current pseudo and loose objective value for a change in a variable's objective value */
extern
RETCODE SCIPlpUpdateVarObj(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   Real             oldobj,             /**< old objective value of variable */
   Real             newobj              /**< new objective value of variable */
   );

/** updates current pseudo and loose objective value for a change in a variable's lower bound */
extern
RETCODE SCIPlpUpdateVarLb(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   Real             oldlb,              /**< old lower bound of variable */
   Real             newlb               /**< new lower bound of variable */
   );

/** updates current pseudo objective value for a change in a variable's upper bound */
extern
RETCODE SCIPlpUpdateVarUb(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   Real             oldub,              /**< old upper bound of variable */
   Real             newub               /**< new upper bound of variable */
   );

/** informs LP, that given variable was added to the problem */
extern
RETCODE SCIPlpUpdateAddVar(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< variable that is now a LOOSE problem variable */
   );

/** informs LP, that given formerly loose problem variable is now a column variable */
extern
RETCODE SCIPlpUpdateVarColumn(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   );

/** informs LP, that given formerly column problem variable is now again a loose variable */
extern
RETCODE SCIPlpUpdateVarLoose(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   );

/** stores the LP solution in the columns and rows */
extern
RETCODE SCIPlpGetSol(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   Bool*            dualfeasible        /**< pointer to store whether the solution is dual feasible, or NULL */
   );

/** stores LP solution with infinite objective value in the columns and rows */
extern
RETCODE SCIPlpGetUnboundedSol(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   );

/** stores the dual farkas multipliers for infeasibility proof in rows */
extern
RETCODE SCIPlpGetDualfarkas(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   );

/** get number of iterations used in last LP solve */
extern
RETCODE SCIPlpGetIterations(
   LP*              lp,                 /**< current LP data */
   int*             iterations          /**< pointer to store the iteration count */
   );

/** increases age of columns with solution value 0.0 and rows with activity not at its bounds,
 *  resets age of non-zero columns and sharp rows
 */
extern
RETCODE SCIPlpUpdateAges(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   );

/** removes all columns and rows in the part of the LP created at the current node, that are too old */
extern
RETCODE SCIPlpRemoveNewObsoletes(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   );

/** removes all columns and rows in whole LP, that are too old */
extern
RETCODE SCIPlpRemoveAllObsoletes(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   );

/** removes all columns at 0.0 and rows not at their bound in the part of the LP created at the current node */
extern
RETCODE SCIPlpCleanupNew(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   );

/** removes all columns at 0.0 and rows not at their bound in the whole LP */
extern
RETCODE SCIPlpCleanupAll(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   );

/** initiates LP diving */
extern
RETCODE SCIPlpStartDive(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   );

/** quits LP diving and resets bounds and objective values of columns to the current node's values */
extern
RETCODE SCIPlpEndDive(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   VAR**            vars,               /**< array with all active variables */
   int              nvars               /**< number of active variables */
   );

/** gets proven lower (dual) bound of last LP solution */
extern
RETCODE SCIPlpGetProvedLowerbound(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   Real*            bound               /**< pointer to store proven dual bound */
   );

/** gets proven dual bound of last LP solution */
extern
RETCODE SCIPlpIsInfeasibilityProved(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   Bool*            proved              /**< pointer to store whether infeasibility is proven */
   );

/** writes LP to a file */
extern 
RETCODE SCIPlpWrite(
   LP*              lp,                 /**< current LP data */
   const char*      fname               /**< file name */
   );



#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets array with columns of the LP */
extern
COL** SCIPlpGetCols(
   LP*              lp                  /**< current LP data */
   );

/** gets current number of columns in LP */
extern
int SCIPlpGetNCols(
   LP*              lp                  /**< current LP data */
   );

/** gets array with rows of the LP */
extern
ROW** SCIPlpGetRows(
   LP*              lp                  /**< current LP data */
   );

/** gets current number of rows in LP */
extern
int SCIPlpGetNRows(
   LP*              lp                  /**< current LP data */
   );

/** gets array with newly added columns after the last mark */
extern
COL** SCIPlpGetNewcols(
   LP*              lp                  /**< current LP data */
   );

/** gets number of newly added columns after the last mark */
extern
int SCIPlpGetNNewcols(
   LP*              lp                  /**< current LP data */
   );

/** gets array with newly added rows after the last mark */
extern
ROW** SCIPlpGetNewrows(
   LP*              lp                  /**< current LP data */
   );

/** gets number of newly added rows after the last mark */
extern
int SCIPlpGetNNewrows(
   LP*              lp                  /**< current LP data */
   );

/** gets the LP solver interface */
extern
LPI* SCIPlpGetLPI(
   LP*              lp                  /**< current LP data */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPlpGetCols(lp)               ((lp)->cols)
#define SCIPlpGetNCols(lp)              ((lp)->ncols)
#define SCIPlpGetRows(lp)               ((lp)->rows)
#define SCIPlpGetNRows(lp)              ((lp)->nrows)
#define SCIPlpGetNewcols(lp)            (&((lp)->cols[(lp)->firstnewcol]))
#define SCIPlpGetNNewcols(lp)           ((lp)->ncols - (lp)->firstnewcol)
#define SCIPlpGetNewrows(lp)            (&((lp)->rows[(lp)->firstnewrow]))
#define SCIPlpGetNNewrows(lp)           ((lp)->nrows - (lp)->firstnewrow)
#define SCIPlpGetLPI(lp)                (lp)->lpi

#endif


#endif
