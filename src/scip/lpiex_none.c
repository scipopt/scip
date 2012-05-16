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

/**@file   lpiex_none.c
 * @brief  dummy interface for the case no exact LP solver is needed
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#ifdef WITH_GMP
#include "gmp.h"
#endif

#include "scip/lpiex.h"
#include "scip/message.h"

#ifdef WITH_GMP

/*
 * Local Methods
 */

/** error handling method */
static
void errorMessage(
   void
   )
{
   SCIPerrorMessage("there is no exact LP solver linked to the binary (LPSEX=none)\n");
   SCIPABORT();
}

/*
 * LP Interface Methods
 */


/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiexGetSolverName(void)
{
   return "NONE";
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiexGetSolverPointer(
   SCIP_LPIEX*           lpi                 /**< pointer to an LP interface structure */
   )
{  /*lint --e{715}*/
   return (void*) NULL;
}
/**@} */


/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** calls initializator of LP solver; this is mainly needed for defining constants in extended and rational precision */
void SCIPlpiexStart(
   void
   )
{
   errorMessage();
}

/** calls deinitializator of LP solver; this is needed for freeing all internal data of the solver, like constants in
 *  extended and rational precision
 */
void SCIPlpiexEnd(
   void
   )
{
   errorMessage();
}

/** creates an LP problem object
 * @return SCIP_OK on success
 * */
SCIP_RETCODE SCIPlpiexCreate(
   SCIP_LPIEX**          lpi,                /**< pointer to an LP interface structure */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiexFree(
   SCIP_LPIEX**            lpi                 /**< pointer to an LP interface structure */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}
/**@} */


/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */


/** copies LP data with column matrix into LP solver */
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
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** adds columns to the LP */
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
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiexDelCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiexDelColset(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** adds rows to the LP */
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
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiexDelRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiexDelRowset(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiexClear(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiexChgBounds(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   int*                  ind,                /**< column indices */
   mpq_t*                lb,                 /**< values for the new lower bounds, or NULL */
   mpq_t*                ub                  /**< values for the new upper bounds, or NULL */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiexChgSides(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const mpq_t*          lhs,                /**< new values for left hand sides */
   const mpq_t*          rhs                 /**< new values for right hand sides */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiexChgCoef(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   mpq_t                 newval              /**< new value of coefficient */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiexChgObjsen(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiexChgObj(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   mpq_t*                obj                 /**< new objective values for columns */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiexScaleRow(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   const mpq_t           scaleval            /**< scaling multiplier */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPlpiexScaleCol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   const mpq_t           scaleval            /**< scaling multiplier */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/**@} */




/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiexGetNRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{  /*lint --e{715}*/
   assert(nrows != NULL);
   *nrows = 0;
   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiexGetNCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{  /*lint --e{715}*/
   assert(ncols != NULL);
   *ncols = 0;
   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiexGetNNonz(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{  /*lint --e{715}*/
   assert(nnonz != NULL);
   return SCIP_PLUGINNOTFOUND;
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
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
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
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
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiexGetObj(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   mpq_t*                vals                /**< array to store objective coefficients */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiexGetBounds(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective value for */
   int                   lastcol,            /**< last column to get objective value for */
   mpq_t*                lbs,                /**< array to store lower bound values, or NULL */
   mpq_t*                ubs                 /**< array to store upper bound values, or NULL */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiexGetSides(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   mpq_t*                lhss,               /**< array to store left hand side values, or NULL */
   mpq_t*                rhss                /**< array to store right hand side values, or NULL */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiexGetCoef(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   mpq_t*                val                 /**< pointer to store the value of the coefficient */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/**@} */

/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiexSolvePrimal(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiexSolveDual(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiexSolveBarrier(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** performs strong branching iterations on all candidates */
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
   )
{  /*lint --e{715}*/
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);
   return SCIP_PLUGINNOTFOUND;
}

/**@} */

/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiexWasSolved(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_RETCODE SCIPlpiexGetSolFeasibility(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{  /*lint --e{715}*/
   assert(primalfeasible != NULL);
   assert(dualfeasible != NULL);
   return SCIP_PLUGINNOTFOUND;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiexExistsPrimalRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiexHasPrimalRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiexIsPrimalUnbounded(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiexIsPrimalInfeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiexIsPrimalFeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiexExistsDualRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiexHasDualRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff LP is dual unbounded */
SCIP_Bool SCIPlpiexIsDualUnbounded(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff LP is dual infeasible */
SCIP_Bool SCIPlpiexIsDualInfeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiexIsDualFeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiexIsOptimal(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiexIsStable(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiexIsObjlimExc(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiexIsIterlimExc(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiexIsTimelimExc(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns the internal solution status of the solver */
int SCIPlpiexGetInternalStatus(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return 0;
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiexIgnoreInstability(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{  /*lint --e{715}*/
   assert(success != NULL);
   return SCIP_PLUGINNOTFOUND;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiexGetObjval(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   mpq_t*                objval              /**< stores the objective value */
   )
{  /*lint --e{715}*/
   assert(objval != NULL);
   return SCIP_PLUGINNOTFOUND;
}

/** gets primal and dual solution vectors */
SCIP_RETCODE SCIPlpiexGetSol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   mpq_t*                objval,             /**< stores the objective value, may be NULL if not needed */
   mpq_t*                primsol,            /**< primal solution vector, may be NULL if not needed */
   mpq_t*                dualsol,            /**< dual solution vector, may be NULL if not needed */
   mpq_t*                activity,           /**< row activity vector, may be NULL if not needed */
   mpq_t*                redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiexGetPrimalRay(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   mpq_t*                ray                 /**< primal ray */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** gets dual farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiexGetDualfarkas(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   mpq_t*                dualfarkas          /**< dual farkas row multipliers */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiexGetIterations(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/**@} */

/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiexGetBase(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiexSetBase(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** returns the indices of the basic columns and rows */
SCIP_RETCODE SCIPlpiexGetBasisInd(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  ind                 /**< basic column n gives value n, basic row m gives value -1-m */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** get dense row of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiexGetBInvRow(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   mpq_t*                coef                /**< pointer to store the coefficients of the row */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** get dense column of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiexGetBInvCol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiexGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   mpq_t*                coef                /**< pointer to store the coefficients of the column */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiexGetBInvARow(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const mpq_t*          binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiexGetBInvRow(), or NULL */
   mpq_t*                coef                /**< vector to return coefficients */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiexGetBInvACol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   mpq_t*                coef                /**< vector to return coefficients */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/**@} */

/*
 * LP State Methods
 */

/**@name LP State Methods */
/**@{ */

/** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiexGetState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiexGetState()
 */
SCIP_RETCODE SCIPlpiexSetState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiexFreeState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiexHasStateBasis(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiexReadState(
   SCIP_LPIEX*           lpi,               /**< LP interface structure */
   const char*           fname              /**< file name */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** writes LP state (like basis information) to a file */
SCIP_RETCODE SCIPlpiexWriteState(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   const char*           fname           /**< file name */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** checks whether LPi state (i.e. basis information) is dual feasible and returns corresponding dual objective value.
 *  if wanted it will first directly test the corresponding approximate dual and primal solution
 *  (corrected via dual variables for bounds and primal variables for slacks if possible) for optimality
 *  before performing the dual feasibility test on the more expensive exact basic solution.
 */
SCIP_RETCODE SCIPlpiexStateDualFeasible(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate,           /**< LPi state information (like basis information) */
   SCIP_Bool             useprestep,         /**< should approximate primal and dual solution first */
   SCIP_Real*            primalsol,          /**< approximate primal solution; or NULL to compute by exact LP solver */
   SCIP_Real*            dualsol,            /**< approximate dual solution; or NULL to compute by exact LP solver */
   SCIP_Bool*            result,             /**< pointer to store whether given LPi state is dual feasible */
   mpq_t*                dualobjval          /**< pointer to store dual objective value in case of dual feasibility */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/**@} */

/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiexGetIntpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   )
{  /*lint --e{715}*/
   return SCIP_PARAMETERUNKNOWN;
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiexSetIntpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{  /*lint --e{715}*/
   return SCIP_PARAMETERUNKNOWN;
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiexGetRealpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   mpq_t*                dval                /**< buffer to store the parameter value */
   )
{  /*lint --e{715}*/
   return SCIP_PARAMETERUNKNOWN;
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiexSetRealpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   mpq_t                 dval                /**< parameter value */
   )
{  /*lint --e{715}*/
   return SCIP_PARAMETERUNKNOWN;
}

/**@} */

/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as positive infinity in the LP solver */
void SCIPlpiexPosInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   mpq_t*                infval          /**< pointer to store positive infinity value of LP solver */
   )
{  /*lint --e{715}*/
   errorMessage();
}

/** checks if given value is treated as positive infinity in the LP solver */
SCIP_Bool SCIPlpiexIsPosInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   const mpq_t           val
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/** returns value treated as negative infinity in the LP solver */
void SCIPlpiexNegInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   mpq_t*                infval          /**< pointer to store negative infinity value of LP solver */
   )
{  /*lint --e{715}*/
   errorMessage();
}

/** checks if given value is treated as negative infinity in the LP solver */
SCIP_Bool SCIPlpiexIsNegInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   const mpq_t           val
   )
{  /*lint --e{715}*/
   errorMessage();
   return FALSE;
}

/**@} */

/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
SCIP_RETCODE SCIPlpiexReadLP(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiexWriteLP(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** computes and stores matrix factorization within the LPIEX structure */
SCIP_RETCODE SCIPlpiexCreateFactor(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   int                   dim,            /**< dimension of matrix */
   int*                  cbeg,           /**< column indices of matrix */
   int*                  clen,           /**< column lengths of matrix */
   int*                  cindx,          /**< row index of entries */
   mpq_t*                ccoef           /**< coef values of matrix */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/** solves a system using the stored factorization */
SCIP_RETCODE SCIPlpiexFactorSolve(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   int                   dim,            /**< dimension of matrix */
   mpq_t*                sol,            /**< solution to system */
   mpq_t*                rhs             /**< rhs of system */
   )
{  /*lint --e{715}*/
   return SCIP_PLUGINNOTFOUND;
}

/**@} */
#endif
