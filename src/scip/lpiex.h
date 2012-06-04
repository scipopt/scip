/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lpiex.h
 * @brief  interface methods for specific LP solvers
 * @author Daniel Espinoza
 * @author Kati Wolter
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_LPIEX_H__
#define __SCIP_LPIEX_H__

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_lpiex.h"

#ifdef WITH_GMP

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
extern
const char* SCIPlpiexGetSolverName(
   void
   );

/** gets pointer for LP solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the LP solver object.
 */
extern
void* SCIPlpiexGetSolverPointer(
   SCIP_LPIEX*           lpi                 /**< pointer to an LP interface structure */
   );

/**@} */




/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** calls initializator of LP solver; this is mainly needed for defining constants in extended and rational precision */
extern
void SCIPlpiexStart(
   void
   );

/** calls deinitializator of LP solver; this is needed for freeing all internal data of the solver, like constants in
 *  extended and rational precision
 */
extern
void SCIPlpiexEnd(
   void
   );

/** creates an LP problem object */
extern
SCIP_RETCODE SCIPlpiexCreate(
   SCIP_LPIEX**          lpi,                /**< pointer to an LP interface structure */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   );

/** deletes an LP problem object */
extern
SCIP_RETCODE SCIPlpiexFree(
   SCIP_LPIEX**          lpi                 /**< pointer to an LP interface structure */
   );

/**@} */




/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies LP data with column matrix into LP solver */
extern
SCIP_RETCODE SCIPlpiexLoadColLP(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   mpq_t*                obj,                /**< objective function values of columns */
   mpq_t*                lb,                 /**< lower bounds of columns */
   mpq_t*                ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   const mpq_t*          lhs,                /**< left hand sides of rows */
   const mpq_t*          rhs,                /**< right hand sides of rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each column in ind- and val-array */
   int*                  ind,                /**< row indices of constraint matrix entries */
   mpq_t*                val                 /**< values of constraint matrix entries */
   );

/** adds columns to the LP */
extern
SCIP_RETCODE SCIPlpiexAddCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   mpq_t*                obj,                /**< objective function values of new columns */
   mpq_t*                lb,                 /**< lower bounds of new columns */
   mpq_t*                ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   int*                  beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   int*                  ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   mpq_t*                val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   );

/** deletes all columns in the given range from LP */
extern
SCIP_RETCODE SCIPlpiexDelCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   );

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
extern
SCIP_RETCODE SCIPlpiexDelColset(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   );

/** adds rows to the LP */
extern
SCIP_RETCODE SCIPlpiexAddRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const mpq_t*          lhs,                /**< left hand sides of new rows */
   const mpq_t*          rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   int*                  beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            len,                /**< number of nonzeros of each row in ind- and val-array, or NULL if only nonzeros */
   int*                  ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   mpq_t*                val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   );

/** deletes all rows in the given range from LP */
extern
SCIP_RETCODE SCIPlpiexDelRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   );

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
extern
SCIP_RETCODE SCIPlpiexDelRowset(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   );

/** clears the whole LP */
extern
SCIP_RETCODE SCIPlpiexClear(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** changes lower and upper bounds of columns */
extern
SCIP_RETCODE SCIPlpiexChgBounds(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   int*                  ind,                /**< column indices */
   mpq_t*                lb,                 /**< values for the new lower bounds, or NULL */
   mpq_t*                ub                  /**< values for the new upper bounds, or NULL */
   );

/** changes left and right hand sides of rows */
extern
SCIP_RETCODE SCIPlpiexChgSides(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const mpq_t*          lhs,                /**< new values for left hand sides */
   const mpq_t*          rhs                 /**< new values for right hand sides */
   );

/** changes a single coefficient */
extern
SCIP_RETCODE SCIPlpiexChgCoef(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   mpq_t                 newval              /**< new value of coefficient */
   );

/** changes the objective sense */
extern
SCIP_RETCODE SCIPlpiexChgObjsen(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   );

/** changes objective values of columns in the LP */
extern
SCIP_RETCODE SCIPlpiexChgObj(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   mpq_t*                obj                 /**< new objective values for columns */
   );

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
extern
SCIP_RETCODE SCIPlpiexScaleRow(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   const mpq_t           scaleval            /**< scaling multiplier */
   );

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
extern
SCIP_RETCODE SCIPlpiexScaleCol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   const mpq_t           scaleval            /**< scaling multiplier */
   );

/**@} */




/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
extern
SCIP_RETCODE SCIPlpiexGetNRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   );

/** gets the number of columns in the LP */
extern
SCIP_RETCODE SCIPlpiexGetNCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   );

/** gets the number of nonzero elements in the LP constraint matrix */
extern
SCIP_RETCODE SCIPlpiexGetNNonz(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   );

/** gets columns from LP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
extern
SCIP_RETCODE SCIPlpiexGetCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   mpq_t*                lb,                 /**< buffer to store the lower bound vector, or NULL */
   mpq_t*                ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   mpq_t*                val                 /**< buffer to store values of constraint matrix entries, or NULL */
   );

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
extern
SCIP_RETCODE SCIPlpiexGetRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   mpq_t*                lhs,                /**< buffer to store left hand side vector, or NULL */
   mpq_t*                rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   mpq_t*                val                 /**< buffer to store values of constraint matrix entries, or NULL */
   );

/** gets objective coefficients from LP problem object */
extern
SCIP_RETCODE SCIPlpiexGetObj(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   mpq_t*                vals                /**< array to store objective coefficients */
   );

/** gets current bounds from LP problem object */
extern
SCIP_RETCODE SCIPlpiexGetBounds(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get bounds for */
   int                   lastcol,            /**< last column to get bounds for */
   mpq_t*                lbs,                /**< array to store lower bound values, or NULL */
   mpq_t*                ubs                 /**< array to store upper bound values, or NULL */
   );

/** gets current row sides from LP problem object */
extern
SCIP_RETCODE SCIPlpiexGetSides(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   mpq_t*                lhss,               /**< array to store left hand side values, or NULL */
   mpq_t*                rhss                /**< array to store right hand side values, or NULL */
   );

/** gets a single coefficient */
extern
SCIP_RETCODE SCIPlpiexGetCoef(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   mpq_t*                val                 /**< pointer to store the value of the coefficient */
   );

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** calls primal simplex to solve the LP */
extern
SCIP_RETCODE SCIPlpiexSolvePrimal(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** calls dual simplex to solve the LP */
extern
SCIP_RETCODE SCIPlpiexSolveDual(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
extern
SCIP_RETCODE SCIPlpiexSolveBarrier(
   SCIP_LPIEX*           lpi,                 /**< LP interface structure */
   SCIP_Bool             crossover            /**< perform crossover */
   );

/** performs strong branching iterations on all candidates */
extern
SCIP_RETCODE SCIPlpiexStrongbranch(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   const mpq_t           psol,               /**< current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   mpq_t*                down,               /**< stores dual bound after branching column down */
   mpq_t*                up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   );

/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
extern
SCIP_Bool SCIPlpiexWasSolved(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** gets information about primal and dual feasibility of the current LP solution */
extern
SCIP_RETCODE SCIPlpiexGetSolFeasibility(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   );

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
extern
SCIP_Bool SCIPlpiexExistsPrimalRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
extern
SCIP_Bool SCIPlpiexHasPrimalRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be primal unbounded */
extern
SCIP_Bool SCIPlpiexIsPrimalUnbounded(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be primal infeasible */
extern
SCIP_Bool SCIPlpiexIsPrimalInfeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be primal feasible */
extern
SCIP_Bool SCIPlpiexIsPrimalFeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
extern
SCIP_Bool SCIPlpiexExistsDualRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
extern
SCIP_Bool SCIPlpiexHasDualRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be dual unbounded */
extern
SCIP_Bool SCIPlpiexIsDualUnbounded(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be dual infeasible */
extern
SCIP_Bool SCIPlpiexIsDualInfeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be dual feasible */
extern
SCIP_Bool SCIPlpiexIsDualFeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP was solved to optimality */
extern
SCIP_Bool SCIPlpiexIsOptimal(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff current LP basis is stable */
extern
SCIP_Bool SCIPlpiexIsStable(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the objective limit was reached */
extern
SCIP_Bool SCIPlpiexIsObjlimExc(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the iteration limit was reached */
extern
SCIP_Bool SCIPlpiexIsIterlimExc(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the time limit was reached */
extern
SCIP_Bool SCIPlpiexIsTimelimExc(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns the internal solution status of the solver */
extern
int SCIPlpiexGetInternalStatus(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
extern
SCIP_RETCODE SCIPlpiexIgnoreInstability(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   );

/** gets objective value of solution */
extern
SCIP_RETCODE SCIPlpiexGetObjval(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   mpq_t*                objval              /**< stores the objective value */
   );

/** gets primal and dual solution vectors for feasible LPs */
extern
SCIP_RETCODE SCIPlpiexGetSol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   mpq_t*                objval,             /**< stores the objective value, may be NULL if not needed */
   mpq_t*                primsol,            /**< primal solution vector, may be NULL if not needed */
   mpq_t*                dualsol,            /**< dual solution vector, may be NULL if not needed */
   mpq_t*                activity,           /**< row activity vector, may be NULL if not needed */
   mpq_t*                redcost             /**< reduced cost vector, may be NULL if not needed */
   );

/** gets primal ray for unbounded LPs */
extern
SCIP_RETCODE SCIPlpiexGetPrimalRay(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   mpq_t*                ray                 /**< primal ray */
   );

/** gets dual farkas proof for infeasibility */
extern
SCIP_RETCODE SCIPlpiexGetDualfarkas(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   mpq_t*                dualfarkas          /**< dual farkas row multipliers */
   );

/** gets the number of LP iterations of the last solve call */
extern
SCIP_RETCODE SCIPlpiexGetIterations(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   );

/**@} */




/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
extern
SCIP_RETCODE SCIPlpiexGetBase(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   );

/** sets current basis status for columns and rows */
extern
SCIP_RETCODE SCIPlpiexSetBase(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   );

/** returns the indices of the basic columns and rows */
extern
SCIP_RETCODE SCIPlpiexGetBasisInd(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   );

/** get dense row of inverse basis matrix B^-1 */
extern
SCIP_RETCODE SCIPlpiexGetBInvRow(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   mpq_t*                coef                /**< pointer to store the coefficients of the row */
   );

/** get dense column of inverse basis matrix B^-1 */
extern
SCIP_RETCODE SCIPlpiexGetBInvCol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiexGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   mpq_t*                coef                /**< pointer to store the coefficients of the column */
   );

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
extern
SCIP_RETCODE SCIPlpiexGetBInvARow(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const mpq_t*          binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiexGetBInvRow(), or NULL */
   mpq_t*                coef                /**< vector to return coefficients */
   );

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
extern
SCIP_RETCODE SCIPlpiexGetBInvACol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   mpq_t*                coef                /**< vector to return coefficients */
   );

/**@} */




/*
 * LPi State Methods
 */

/**@name LPi State Methods */
/**@{ */

/** stores LPi state (like basis information) into lpistate object */
extern
SCIP_RETCODE SCIPlpiexGetState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   );

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiexGetState()
 */
extern
SCIP_RETCODE SCIPlpiexSetState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   );

/** frees LPi state information */
extern
SCIP_RETCODE SCIPlpiexFreeState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   );

/** checks, whether the given LPi state contains simplex basis information */
extern
SCIP_Bool SCIPlpiexHasStateBasis(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   );

/** reads LPi state (like basis information from a file */
extern
SCIP_RETCODE SCIPlpiexReadState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   );

/** writes LPi state (like basis information) to a file */
extern
SCIP_RETCODE SCIPlpiexWriteState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   );

/** checks whether LPi state (i.e. basis information) is dual feasbile and returns corresponding dual objective value.
 *  if wanted it will first directly test the corresponding approximate dual and primal solution
 *  (corrected via dual variables for bounds and primal variables for slacks if possible) for optimality
 *  before performing the dual feasibility test on the more expensive exact basic solution.
 */
extern
SCIP_RETCODE SCIPlpiexStateDualFeasible(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate,           /**< LPi state information (like basis information) */
   SCIP_Bool             useprestep,         /**< should approximate primal and dual solution first */
   SCIP_Real*            primalsol,          /**< approximate primal solution; or NULL to compute by exact LP solver */
   SCIP_Real*            dualsol,            /**< approximate dual solution; or NULL to compute by exact LP solver */
   SCIP_Bool*            result,             /**< pointer to store whether given LPi state is dual feasible */
   mpq_t*                dualobjval          /**< pointer to store dual objective value in case of dual feasibility */
   );

/**@} */




/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
extern
SCIP_RETCODE SCIPlpiexGetIntpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   );

/** sets integer parameter of LP */
extern
SCIP_RETCODE SCIPlpiexSetIntpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   );

/** gets floating point parameter of LP */
extern
SCIP_RETCODE SCIPlpiexGetRealpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   mpq_t*                dval                /**< buffer to store the parameter value */
   );

/** sets floating point parameter of LP */
extern
SCIP_RETCODE SCIPlpiexSetRealpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   mpq_t                 dval                /**< parameter value */
   );

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as positive infinity in the LP solver */
extern
void SCIPlpiexPosInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   mpq_t*                infval          /**< pointer to store positive infinity value of LP solver */
   );

/** checks if given value is treated as positive infinity in the LP solver */
extern
SCIP_Bool SCIPlpiexIsPosInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   const mpq_t           val
   );

/** returns value treated as negative infinity in the LP solver */
extern
void SCIPlpiexNegInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   mpq_t*                infval          /**< pointer to store negative infinity value of LP solver */
   );

/** checks if given value is treated as negative infinity in the LP solver */
extern
SCIP_Bool SCIPlpiexIsNegInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   const mpq_t           val
   );


/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
extern
SCIP_RETCODE SCIPlpiexReadLP(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   );

/** writes LP to a file */
extern
SCIP_RETCODE SCIPlpiexWriteLP(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   );

/**@} */


/*
 * Exact LU decomposition solver interface
 */

/**@name Exact LU decomposition solver interface */
/**@{ */

/** computes and stores matrix factorization within the LPIEX structure */
extern
SCIP_RETCODE SCIPlpiexCreateFactor(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   int                   dim,            /**< dimension of matrix */
   int*                  cbeg,           /**< column indices of matrix */
   int*                  clen,           /**< column lengths of matrix */
   int*                  cindx,          /**< row index of entries */
   mpq_t*                ccoef           /**< coef values of matrix */
    );


/** solves a system using the stored factorization */
extern
SCIP_RETCODE SCIPlpiexFactorSolve(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   int                   dim,            /**< dimension of matrix */
   mpq_t*                sol,            /**< solution to system */
   mpq_t*                rhs             /**< rhs of system */
   );
/**@} */

#ifdef __cplusplus
}
#endif

#endif

#endif
