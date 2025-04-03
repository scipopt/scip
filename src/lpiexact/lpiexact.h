/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lpiexact.h
 * @brief  interface methods for specific exact LP solvers
 * @author Daniel Espinoza
 * @author Kati Wolter
 * @author Marc Pfetsch
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_LPIEXACT_H__
#define __SCIP_LPIEXACT_H__

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "lpiexact/type_lpiexact.h"
#include "scip/rational.h"

#ifdef __cplusplus
extern "C" {
#endif

/** gets name and version of LP solver */
SCIP_EXPORT
const char* SCIPlpiExactGetSolverName(
   void
   );

/** gets description of LP solver (developer, webpage, ...) */
SCIP_EXPORT
const char* SCIPlpiExactGetSolverDesc(
   void
   );

/** gets name and version of external package required for LP solver */
SCIP_EXPORT
const char* SCIPlpiExactGetExternalCodeName(
   void
   );

/** gets description of external package required for LP solver (developer, webpage, ...) */
SCIP_EXPORT
const char* SCIPlpiExactGetExternalCodeDesc(
   void
   );

/** gets pointer for LP solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the LP solver object.
 */
SCIP_EXPORT
void* SCIPlpiExactGetSolverPointer(
   SCIP_LPIEXACT*        lpi                 /**< pointer to an LP interface structure */
   );

/**@} */




/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** calls initializator of LP solver; this is mainly needed for defining constants in extended and rational precision */
SCIP_EXPORT
void SCIPlpiExactStart(
   void
   );

/** calls deinitializator of LP solver; this is needed for freeing all internal data of the solver, like constants in
 *  extended and rational precision
 */
SCIP_EXPORT
void SCIPlpiExactEnd(
   void
   );

/** creates an LP problem object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactCreate(
   SCIP_LPIEXACT**       lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   );

/** deletes an LP problem object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactFree(
   SCIP_LPIEXACT**       lpi                 /**< pointer to an LP interface structure */
   );

/**@} */




/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies LP data with column matrix into LP solver */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactLoadColLP(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   SCIP_RATIONAL**       obj,                /**< objective function values of columns */
   SCIP_RATIONAL**       lb,                 /**< lower bounds of columns */
   SCIP_RATIONAL**       ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   SCIP_RATIONAL**       lhs,                /**< left hand sides of rows */
   SCIP_RATIONAL**       rhs,                /**< right hand sides of rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each column in ind- and val-array */
   int*                  ind,                /**< row indices of constraint matrix entries */
   SCIP_RATIONAL**       val                 /**< values of constraint matrix entries */
   );

/** adds columns to the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactAddCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   SCIP_RATIONAL**       obj,                /**< objective function values of new columns */
   SCIP_RATIONAL**       lb,                 /**< lower bounds of new columns */
   SCIP_RATIONAL**       ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   int*                  beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   int*                  ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   SCIP_RATIONAL**       val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   );

/** deletes all columns in the given range from LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactDelCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   );

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactDelColset(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   );

/** adds rows to the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactAddRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   SCIP_RATIONAL**       lhs,                /**< left hand sides of new rows */
   SCIP_RATIONAL**       rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   int*                  beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   int*                  ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   SCIP_RATIONAL**       val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   );

/** deletes all rows in the given range from LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactDelRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   );

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactDelRowset(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   );

/** clears the whole LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactClear(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** changes lower and upper bounds of columns */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactChgBounds(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   int*                  ind,                /**< column indices */
   SCIP_RATIONAL**       lb,                 /**< values for the new lower bounds, or NULL */
   SCIP_RATIONAL**       ub                  /**< values for the new upper bounds, or NULL */
   );

/** changes left and right hand sides of rows */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactChgSides(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   int*                  ind,                /**< row indices */
   SCIP_RATIONAL**       lhs,                /**< new values for left hand sides */
   SCIP_RATIONAL**       rhs                 /**< new values for right hand sides */
   );

/** changes a single coefficient */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactChgCoef(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_RATIONAL*        newval              /**< new value of coefficient */
   );

/** changes the objective sense */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactChgObjsen(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   );

/** changes objective values of columns in the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactChgObj(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   SCIP_RATIONAL**       obj                 /**< new objective values for columns */
   );

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactScaleRow(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_RATIONAL*        scaleval            /**< scaling multiplier */
   );

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactScaleCol(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_RATIONAL*        scaleval            /**< scaling multiplier */
   );

/**@} */




/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetNRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   );

/** gets the number of columns in the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetNCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   );

/** gets the objective sense of the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetObjsen(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   );

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetNNonz(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   );

/** gets columns from LP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   SCIP_RATIONAL**       lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_RATIONAL**       ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_RATIONAL**       val                 /**< buffer to store values of constraint matrix entries, or NULL */
   );

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_RATIONAL**       lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_RATIONAL**       rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_RATIONAL**       val                 /**< buffer to store values of constraint matrix entries, or NULL */
   );

/** gets column names */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetColNames(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get name from LP */
   int                   lastcol,            /**< last column to get name from LP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for col names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   );

/** gets row names */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetRowNames(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get name from LP */
   int                   lastrow,            /**< last row to get name from LP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for row names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   );

/** gets objective coefficients from LP problem object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetObj(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_RATIONAL**       vals                /**< array to store objective coefficients */
   );

/** gets current bounds from LP problem object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetBounds(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get bounds for */
   int                   lastcol,            /**< last column to get bounds for */
   SCIP_RATIONAL**       lbs,                /**< array to store lower bound values, or NULL */
   SCIP_RATIONAL**       ubs                 /**< array to store upper bound values, or NULL */
   );

/** gets current row sides from LP problem object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetSides(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_RATIONAL**       lhss,               /**< array to store left hand side values, or NULL */
   SCIP_RATIONAL**       rhss                /**< array to store right hand side values, or NULL */
   );

/** gets a single coefficient */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetCoef(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_RATIONAL*        val                 /**< pointer to store the value of the coefficient */
   );

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** calls primal simplex to solve the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactSolvePrimal(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** calls dual simplex to solve the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactSolveDual(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactSolveBarrier(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   );

/** start strong branching - call before any strong branching */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactStartStrongbranch(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** end strong branching - call after any strong branching */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactEndStrongbranch(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** performs strong branching iterations on all candidates */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactStrongbranch(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   const SCIP_RATIONAL*  psol,               /**< current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_RATIONAL*        down,               /**< stores dual bound after branching column down */
   SCIP_RATIONAL*        up,                 /**< stores dual bound after branching column up */
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
SCIP_Bool SCIPlpiExactWasSolved(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetSolFeasibility(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   );

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactExistsPrimalRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactHasPrimalRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsPrimalUnbounded(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsPrimalInfeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsPrimalFeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactExistsDualRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactHasDualRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsDualUnbounded(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsDualInfeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsDualFeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP was solved to optimality */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsOptimal(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff current LP solution is stable
 *
 *  This function should return true if the solution is reliable, i.e., feasible and optimal (or proven
 *  infeasible/unbounded) with respect to the original problem. The optimality status might be with respect to a scaled
 *  version of the problem, but the solution might not be feasible to the unscaled original problem; in this case,
 *  SCIPlpiIsStable() should return false.
 */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsStable(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the objective limit was reached */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsObjlimExc(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the iteration limit was reached */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsIterlimExc(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the time limit was reached */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsTimelimExc(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** returns the internal solution status of the solver */
SCIP_EXPORT
int SCIPlpiExactGetInternalStatus(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactIgnoreInstability(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   );

/** gets objective value of solution */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetObjval(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        objval              /**< stores the objective value */
   );

/** gets primal and dual solution vectors for feasible LPs */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetSol(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_RATIONAL**       primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_RATIONAL**       dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_RATIONAL**       activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_RATIONAL**       redcost             /**< reduced cost vector, may be NULL if not needed */
   );

/** gets primal ray for unbounded LPs */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetPrimalRay(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL**       ray                 /**< primal ray */
   );

/** gets dual farkas proof for infeasibility */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetDualfarkas(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL**       dualfarkas          /**< dual farkas row multipliers */
   );

/** gets the number of LP iterations of the last solve call */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetIterations(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
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
SCIP_RETCODE SCIPlpiExactGetBase(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   );

/** sets current basis status for columns and rows */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactSetBase(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   );

/** returns the indices of the basic columns and rows */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetBasisInd(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   );

/** get dense row of inverse basis matrix B^-1 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetBInvRow(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_RATIONAL**       coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   );

/** get dense column of inverse basis matrix B^-1 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetBInvCol(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiExactGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   SCIP_RATIONAL**       coef,               /**< pointer to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   );

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetBInvARow(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_RATIONAL**       binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiExactGetBInvRow(), or NULL */
   SCIP_RATIONAL**       coef,               /**< vector to return coefficients */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   );

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetBInvACol(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_RATIONAL**       coef,               /**< vector to return coefficients */
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
SCIP_RETCODE SCIPlpiExactGetState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   );

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiExactGetState()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactSetState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   );

/** clears current LPi state (like basis information) of the solver */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactClearState(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** frees LPi state information */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactFreeState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   );

/** checks, whether the given LPi state contains simplex basis information */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactHasStateBasis(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   );

/** reads LPi state (like basis information from a file */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactReadState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   );

/** writes LPi state (like basis information) to a file */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactWriteState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   );

/** checks whether LPi state (i.e. basis information) is dual feasbile and returns corresponding dual objective value.
 *  if wanted it will first directly test the corresponding approximate dual and primal solution
 *  (corrected via dual variables for bounds and primal variables for slacks if possible) for optimality
 *  before performing the dual feasibility test on the more expensive exact basic solution.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactStateDualFeasible(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate,           /**< LPi state information (like basis information) */
   SCIP_Bool             useprestep,         /**< should approximate primal and dual solution first */
   SCIP_Real*            primalsol,          /**< approximate primal solution; or NULL to compute by exact LP solver */
   SCIP_Real*            dualsol,            /**< approximate dual solution; or NULL to compute by exact LP solver */
   SCIP_Bool*            result,             /**< pointer to store whether given LPi state is dual feasible */
   SCIP_RATIONAL**       dualobjval          /**< pointer to store dual objective value in case of dual feasibility */
   );

/**@} */

/*
 * LPi Pricing Norms Methods
 */

/**@name LPi Pricing Norms Methods */
/**@{ */

/** stores lpiexact pricing norms into lpiexactnorms object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetNorms(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information */
   );

/** loads LPi pricing norms into solver; note that the LP might have been extended with additional
 *  columns and rows since the norms were stored with SCIPlpiGetNorms()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactSetNorms(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const SCIP_LPINORMS*  lpinorms            /**< LPi pricing norms information, or NULL */
   );

/** frees LPi pricing norms information */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactFreeNorms(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
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
SCIP_RETCODE SCIPlpiExactGetIntpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   );

/** sets integer parameter of LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactSetIntpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   );

/** gets floating point parameter of LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactGetRealpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   );

/** sets floating point parameter of LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactSetRealpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
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
void SCIPlpiExactPosInfinity(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        infval              /**< pointer to store positive infinity value of LP solver */
   );

/** checks if given value is treated as positive infinity in the LP solver */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsPosInfinity(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        val                 /**< given value */
   );

/** returns value treated as negative infinity in the LP solver */
SCIP_EXPORT
void SCIPlpiExactNegInfinity(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        infval              /**< pointer to store negative infinity value of LP solver */
   );

/** checks if given value is treated as negative infinity in the LP solver */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsNegInfinity(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        val                 /**< given value */
   );

/** returns value treated as infinity in the LP solver */
SCIP_EXPORT
SCIP_Real SCIPlpiExactInfinity(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/** checks if given value is treated as infinity in the LP solver */
SCIP_EXPORT
SCIP_Bool SCIPlpiExactIsInfinity(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< value to test */
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
SCIP_RETCODE SCIPlpiExactReadLP(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   );

/** writes LP to a file */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactWriteLP(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   );

/**@} */


/** prints additional lpiexact internal info */
SCIP_EXPORT
void SCIPlpiExactPrintInfo(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   );

/**@} */

/*
 * Exact LU decomposition solver interface
 */

/**@name Exact LU decomposition solver interface */
/**@{ */

/** computes and stores matrix factorization within the LPIEXACT structure */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactCreateFactor(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   dim,                /**< dimension of matrix */
   int*                  cbeg,               /**< column indices of matrix */
   int*                  clen,               /**< column lengths of matrix */
   int*                  cindx,              /**< row index of entries */
   SCIP_RATIONAL*        ccoef               /**< coef values of matrix */
    );


/** solves a system using the stored factorization */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactFactorSolve(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   dim,                /**< dimension of matrix */
   SCIP_RATIONAL*        sol,                /**< solution to system */
   SCIP_RATIONAL*        rhs                 /**< rhs of system */
   );
/**@} */

#ifdef __cplusplus
}
#endif

#endif
