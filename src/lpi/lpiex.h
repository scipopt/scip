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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lpiex.h
 * @brief  interface methods for specific exact LP solvers
 * @author Daniel Espinoza
 * @author Kati Wolter
 * @author Marc Pfetsch
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_LPIEX_H__
#define __SCIP_LPIEX_H__

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "lpi/type_lpiex.h"
#include "scip/rational.h"

#ifdef __cplusplus
extern "C" {
#endif

/** gets name and version of LP solver */
SCIP_EXPORT
const char* SCIPlpiexGetSolverName(
   void
   );

/** gets description of LP solver (developer, webpage, ...) */
SCIP_EXPORT
const char* SCIPlpiexGetSolverDesc(
   void
   );

/** gets name and version of external package required for LP solver */
SCIP_EXPORT
const char* SCIPlpiexGetExternalCodeName(
   void
   );

/** gets description of external package required for LP solver (developer, webpage, ...) */
SCIP_EXPORT
const char* SCIPlpiexGetExternalCodeDesc(
   void
   );

/** gets pointer for LP solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the LP solver object.
 */
SCIP_EXPORT
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
SCIP_EXPORT
void SCIPlpiexStart(
   void
   );

/** calls deinitializator of LP solver; this is needed for freeing all internal data of the solver, like constants in
 *  extended and rational precision
 */
SCIP_EXPORT
void SCIPlpiexEnd(
   void
   );

/** creates an LP problem object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexCreate(
   SCIP_LPIEX**          lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   );

/** deletes an LP problem object */
SCIP_EXPORT
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
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexLoadColLP(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   SCIP_Rational**       obj,                /**< objective function values of columns */
   SCIP_Rational**       lb,                 /**< lower bounds of columns */
   SCIP_Rational**       ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   SCIP_Rational**       lhs,                /**< left hand sides of rows */
   SCIP_Rational**       rhs,                /**< right hand sides of rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each column in ind- and val-array */
   int*                  ind,                /**< row indices of constraint matrix entries */
   SCIP_Rational**       val                 /**< values of constraint matrix entries */
   );

/** adds columns to the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexAddCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   SCIP_Rational**       obj,                /**< objective function values of new columns */
   SCIP_Rational**       lb,                 /**< lower bounds of new columns */
   SCIP_Rational**       ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   int*                  beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   int*                  ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   SCIP_Rational**       val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   );

/** deletes all columns in the given range from LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexDelCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   );

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexDelColset(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   );

/** adds rows to the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexAddRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   SCIP_Rational**       lhs,                /**< left hand sides of new rows */
   SCIP_Rational**       rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   int*                  beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   int*                  ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   SCIP_Rational**       val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   );

/** deletes all rows in the given range from LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexDelRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   );

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexDelRowset(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   );

/** clears the whole LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexClear(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** changes lower and upper bounds of columns */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexChgBounds(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   int*                  ind,                /**< column indices */
   SCIP_Rational**        lb,                 /**< values for the new lower bounds, or NULL */
   SCIP_Rational**        ub                  /**< values for the new upper bounds, or NULL */
   );

/** changes left and right hand sides of rows */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexChgSides(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   int*                  ind,                /**< row indices */
   SCIP_Rational**       lhs,                /**< new values for left hand sides */
   SCIP_Rational**       rhs                 /**< new values for right hand sides */
   );

/** changes a single coefficient */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexChgCoef(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Rational*        newval              /**< new value of coefficient */
   );

/** changes the objective sense */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexChgObjsen(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   );

/** changes objective values of columns in the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexChgObj(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   SCIP_Rational**       obj                 /**< new objective values for columns */
   );

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexScaleRow(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   const SCIP_Rational*  scaleval            /**< scaling multiplier */
   );

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexScaleCol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   const SCIP_Rational*  scaleval            /**< scaling multiplier */
   );

/**@} */




/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetNRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   );

/** gets the number of columns in the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetNCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   );

/** gets the objective sense of the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetObjsen(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   );

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetNNonz(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   );

/** gets columns from LP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   SCIP_Rational**       lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Rational**       ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Rational**       val                 /**< buffer to store values of constraint matrix entries, or NULL */
   );

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_Rational**       lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_Rational**       rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_Rational**       val                 /**< buffer to store values of constraint matrix entries, or NULL */
   );

/** gets column names */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetColNames(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get name from LP */
   int                   lastcol,            /**< last column to get name from LP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for col names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   );

/** gets row names */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetRowNames(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get name from LP */
   int                   lastrow,            /**< last row to get name from LP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for row names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   );

/** gets objective coefficients from LP problem object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetObj(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Rational**       vals                /**< array to store objective coefficients */
   );

/** gets current bounds from LP problem object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetBounds(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get bounds for */
   int                   lastcol,            /**< last column to get bounds for */
   SCIP_Rational**       lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Rational**       ubs                 /**< array to store upper bound values, or NULL */
   );

/** gets current row sides from LP problem object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetSides(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Rational**       lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Rational**       rhss                /**< array to store right hand side values, or NULL */
   );

/** gets a single coefficient */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetCoef(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Rational*        val                 /**< pointer to store the value of the coefficient */
   );

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** calls primal simplex to solve the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexSolvePrimal(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** calls dual simplex to solve the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexSolveDual(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexSolveBarrier(
   SCIP_LPIEX*           lpi,                 /**< LP interface structure */
   SCIP_Bool             crossover            /**< perform crossover */
   );

/** start strong branching - call before any strong branching */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexStartStrongbranch(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** end strong branching - call after any strong branching */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexEndStrongbranch(
   SCIP_LPIEX*             lpi                 /**< LP interface structure */
   );

/** performs strong branching iterations on all candidates */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexStrongbranch(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   const SCIP_Rational*  psol,               /**< current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Rational*        down,               /**< stores dual bound after branching column down */
   SCIP_Rational*        up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   );


/* TODO: Do we need the other strong branchiing methods? */

/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
SCIP_EXPORT
SCIP_Bool SCIPlpiexWasSolved(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetSolFeasibility(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   );

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_EXPORT
SCIP_Bool SCIPlpiexExistsPrimalRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_EXPORT
SCIP_Bool SCIPlpiexHasPrimalRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsPrimalUnbounded(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsPrimalInfeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsPrimalFeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_EXPORT
SCIP_Bool SCIPlpiexExistsDualRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_EXPORT
SCIP_Bool SCIPlpiexHasDualRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsDualUnbounded(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsDualInfeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsDualFeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP was solved to optimality */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsOptimal(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff current LP solution is stable
 *
 *  This function should return true if the solution is reliable, i.e., feasible and optimal (or proven
 *  infeasible/unbounded) with respect to the original problem. The optimality status might be with respect to a scaled
 *  version of the problem, but the solution might not be feasible to the unscaled original problem; in this case,
 *  SCIPlpiIsStable() should return false.
 */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsStable(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the objective limit was reached */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsObjlimExc(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the iteration limit was reached */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsIterlimExc(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the time limit was reached */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsTimelimExc(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** returns the internal solution status of the solver */
SCIP_EXPORT
int SCIPlpiexGetInternalStatus(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexIgnoreInstability(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   );

/** gets objective value of solution */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetObjval(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Rational*        objval              /**< stores the objective value */
   );

/** gets primal and dual solution vectors for feasible LPs */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetSol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Rational*        objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Rational**       primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_Rational**       dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_Rational**       activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_Rational**       redcost             /**< reduced cost vector, may be NULL if not needed */
   );

/** gets primal ray for unbounded LPs */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetPrimalRay(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Rational**       ray                 /**< primal ray */
   );

/** gets dual farkas proof for infeasibility */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetDualfarkas(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Rational**       dualfarkas          /**< dual farkas row multipliers */
   );

/** gets the number of LP iterations of the last solve call */
SCIP_EXPORT
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
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetBase(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   );

/** sets current basis status for columns and rows */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexSetBase(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   );

/** returns the indices of the basic columns and rows */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetBasisInd(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   );

/** get dense row of inverse basis matrix B^-1 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetBInvRow(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_Rational**       coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   );

/** get dense column of inverse basis matrix B^-1 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetBInvCol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiexGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   SCIP_Rational**       coef,               /**< pointer to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   );

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetBInvARow(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_Rational**       binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiexGetBInvRow(), or NULL */
   SCIP_Rational**       coef,               /**< vector to return coefficients */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   );

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetBInvACol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Rational**       coef,               /**< vector to return coefficients */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   );

/**@} */




/*
 * LPi State Methods
 */

/**@name LPi State Methods */
/**@{ */

/** stores LPi state (like basis information) into lpistate object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   );

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiexGetState()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexSetState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   );

/** clears current LPi state (like basis information) of the solver */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexClearState(
   SCIP_LPIEX*             lpi               /**< LP interface structure */
   );

/** frees LPi state information */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexFreeState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   );

/** checks, whether the given LPi state contains simplex basis information */
SCIP_EXPORT
SCIP_Bool SCIPlpiexHasStateBasis(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   );

/** reads LPi state (like basis information from a file */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexReadState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   );

/** writes LPi state (like basis information) to a file */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexWriteState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   );

/** checks whether LPi state (i.e. basis information) is dual feasbile and returns corresponding dual objective value.
 *  if wanted it will first directly test the corresponding approximate dual and primal solution
 *  (corrected via dual variables for bounds and primal variables for slacks if possible) for optimality
 *  before performing the dual feasibility test on the more expensive exact basic solution.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexStateDualFeasible(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate,           /**< LPi state information (like basis information) */
   SCIP_Bool             useprestep,         /**< should approximate primal and dual solution first */
   SCIP_Real*            primalsol,          /**< approximate primal solution; or NULL to compute by exact LP solver */
   SCIP_Real*            dualsol,            /**< approximate dual solution; or NULL to compute by exact LP solver */
   SCIP_Bool*            result,             /**< pointer to store whether given LPi state is dual feasible */
   SCIP_Rational**       dualobjval          /**< pointer to store dual objective value in case of dual feasibility */
   );

/**@} */

/*
 * LPi Pricing Norms Methods
 */

/**@name LPi Pricing Norms Methods */
/**@{ */

/** stores lpiex pricing norms into lpiexnorms object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetNorms(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information */
   );

/** loads LPi pricing norms into solver; note that the LP might have been extended with additional
 *  columns and rows since the norms were stored with SCIPlpiGetNorms()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexSetNorms(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const SCIP_LPINORMS*  lpinorms            /**< LPi pricing norms information, or NULL */
   );

/** frees LPi pricing norms information */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexFreeNorms(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information, or NULL */
   );


/**@} */

/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetIntpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   );

/** sets integer parameter of LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexSetIntpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   );

/** gets floating point parameter of LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexGetRealpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   );

/** sets floating point parameter of LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexSetRealpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   );

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as positive infinity in the LP solver */
SCIP_EXPORT
void SCIPlpiexPosInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   SCIP_Rational*        infval          /**< pointer to store positive infinity value of LP solver */
   );

/** checks if given value is treated as positive infinity in the LP solver */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsPosInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   const SCIP_Rational*  val             /**< given value */
   );

/** returns value treated as negative infinity in the LP solver */
SCIP_EXPORT
void SCIPlpiexNegInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   SCIP_Rational*        infval          /**< pointer to store negative infinity value of LP solver */
   );

/** checks if given value is treated as negative infinity in the LP solver */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsNegInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   const SCIP_Rational*  val             /**< given value */
   );

/** returns value treated as infinity in the LP solver */
SCIP_EXPORT
SCIP_Real SCIPlpiexInfinity(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/** checks if given value is treated as infinity in the LP solver */
SCIP_EXPORT
SCIP_Bool SCIPlpiexIsInfinity(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Real             val
   );

/**@} */


/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexReadLP(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   );

/** writes LP to a file */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexWriteLP(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   );

/**@} */


/** prints additional lpiex internal info */
SCIP_EXPORT
void SCIPlpiexPrintInfo(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   );

/**@} */

/*
 * Exact LU decomposition solver interface
 */

/**@name Exact LU decomposition solver interface */
/**@{ */

/** computes and stores matrix factorization within the LPIEX structure */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexCreateFactor(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   int                   dim,            /**< dimension of matrix */
   int*                  cbeg,           /**< column indices of matrix */
   int*                  clen,           /**< column lengths of matrix */
   int*                  cindx,          /**< row index of entries */
   SCIP_Rational*        ccoef           /**< coef values of matrix */
    );


/** solves a system using the stored factorization */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiexFactorSolve(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   int                   dim,            /**< dimension of matrix */
   SCIP_Rational*        sol,            /**< solution to system */
   SCIP_Rational*        rhs             /**< rhs of system */
   );
/**@} */

#ifdef __cplusplus
}
#endif

#endif
