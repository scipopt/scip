/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   lpiexact_none.c
 * @ingroup LPIEXACTS
 * @brief  dummy interface for the case no LP solver is needed
 * @author Leon Eifler
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>

#include "lpiexact/lpiexact.h"
#include "scip/pub_message.h"

#define LPINAME          "NONE"              /**< name of the LPI interface */
#define LPIINFINITY       1e20               /**< infinity value */


/* globally turn off lint warnings: */
/*lint --e{715}*/

/** LP interface
 *
 *  Store several statistic values about the LP. These values are only needed in order to provide a rudimentary
 *  communication, e.g., there are asserts that check the number of rows and columns.
 */
struct SCIP_LPiExact
{
   int                   nrows;              /**< number of rows */
   int                   ncols;              /**< number of columns */
};


/*
 * Local Methods
 */

/** error handling method */
static
void errorMessageAbort(
   void
   )
{  /*lint --e{2707}*/
   SCIPerrorMessage("No exact LP solver available (LPSEXACT=none).\n");
   SCIPABORT();
}

/** error handling method */
static
void errorMessage(
   void
   )
{
   SCIPerrorMessage("No exact LP solver available (LPSEXACT=none).\n");
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
const char* SCIPlpiExactGetSolverName(
   void
   )
{
   return LPINAME;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiExactGetSolverDesc(
   void
   )
{
   return "dummy LP solver interface which solely purpose is to resolve references at linking";
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiExactGetSolverPointer(
   SCIP_LPIEXACT*        lpi                 /**< pointer to an LP interface structure */
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

/** creates an LP problem object */
SCIP_RETCODE SCIPlpiExactCreate(
   SCIP_LPIEXACT**       lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(name != NULL);
   SCIPdebugMessage("SCIPlpiExactCreate()\n");
   SCIPdebugMessage("Note that there is no exact LP solver linked to the binary\n");

   /* create empty LPI */
   SCIP_ALLOC( BMSallocMemory(lpi) );
   (*lpi)->nrows = 0;
   (*lpi)->ncols = 0;

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiExactFree(
   SCIP_LPIEXACT**       lpi                 /**< pointer to an LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   SCIPdebugMessage("SCIPlpiExactFree()\n");

   BMSfreeMemory(lpi);

   return SCIP_OKAY;
}

/**@} */


/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies LP data with column matrix into LP solver */
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
   )
{  /*lint --e{715}*/
#ifndef NDEBUG
   {
      int j;
      for( j = 0; j < nnonz; j++ )
         assert( val[j] != 0 );
   }
#endif

   assert( lpi != NULL );
   assert( lhs != NULL );
   assert( rhs != NULL );
   assert( obj != NULL );
   assert( lb != NULL );
   assert( ub != NULL );
   assert( beg != NULL );
   assert( ind != NULL );
   assert( val != NULL );

   lpi->nrows = nrows;
   lpi->ncols = ncols;
   assert( lpi->nrows >= 0 );
   assert( lpi->ncols >= 0 );

   return SCIP_OKAY;
}

/** adds columns to the LP */
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
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( lpi->ncols >= 0 );
   assert( obj != NULL );
   assert( lb != NULL );
   assert( ub != NULL) ;
   assert( nnonz == 0 || beg != NULL );
   assert( nnonz == 0 || ind != NULL );
   assert( nnonz == 0 || val != NULL );
   assert( nnonz >= 0 );
   assert( ncols >= 0 );

   lpi->ncols += ncols;

   return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiExactDelCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( lpi->ncols >= 0 );

   lpi->ncols -= lastcol - firstcol + 1;
   assert( lpi->ncols >= 0 );

   return SCIP_OKAY;
}

/** adds rows to the LP */
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
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( lpi->nrows >= 0 );
   assert( lhs != NULL );
   assert( rhs != NULL );
   assert( nnonz == 0 || beg != NULL );
   assert( nnonz == 0 || ind != NULL );
   assert( nnonz == 0 || val != NULL );

   lpi->nrows += nrows;

   return SCIP_OKAY;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiExactDelRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( lpi->nrows >= 0 );

   lpi->nrows -= lastrow - firstrow + 1;
   assert( lpi->nrows >= 0 );

   return SCIP_OKAY;
}

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiExactDelRowset(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{  /*lint --e{715}*/
   int cnt = 0;
   int i;

   assert( lpi != NULL );
   assert( dstat != NULL );
   assert( lpi->nrows >= 0 );

   for (i = 0; i < lpi->nrows; ++i)
   {
      if ( dstat[i] )
      {
         ++cnt;
         dstat[i] = -1;
      }
      else
         dstat[i] = cnt;
   }
   lpi->nrows -= cnt;
   assert( lpi->nrows >= 0 );

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiExactClear(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( lpi->nrows >= 0 );
   assert( lpi->ncols >= 0 );

   lpi->nrows = 0;
   lpi->ncols = 0;

   return SCIP_OKAY;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiExactChgBounds(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   int*                  ind,                /**< column indices or NULL if ncols is zero */
   SCIP_RATIONAL**       lb,                 /**< values for the new lower bounds or NULL if ncols is zero */
   SCIP_RATIONAL**       ub                  /**< values for the new upper bounds or NULL if ncols is zero */
   )
{  /*lint --e{715}*/
   assert( ncols == 0 || (ind != NULL && lb != NULL && ub != NULL) );
   return SCIP_OKAY;
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiExactChgSides(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   int*                  ind,                /**< row indices */
   SCIP_RATIONAL**       lhs,                /**< new values for left hand sides */
   SCIP_RATIONAL**       rhs                 /**< new values for right hand sides */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( ind != NULL );
   assert( lhs != NULL );
   assert( rhs != NULL );

   return SCIP_OKAY;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiExactChgCoef(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_RATIONAL*        newval              /**< new value of coefficient */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiExactChgObjsen(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   return SCIP_OKAY;
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiExactChgObj(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   SCIP_RATIONAL**       obj                 /**< new objective values for columns */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   return SCIP_OKAY;
}

/**@} */


/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiExactGetNRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( nrows != NULL );
   assert( lpi->nrows >= 0 );

   *nrows = lpi->nrows;

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiExactGetNCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( ncols != NULL );
   assert( lpi->ncols >= 0 );

   *ncols = lpi->ncols;

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiExactGetNNonz(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{  /*lint --e{715}*/
   assert( nnonz != NULL );
   assert( lpi != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiExactGetCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   SCIP_RATIONAL**       lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_RATIONAL**       ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_RATIONAL**       val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiExactGetRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_RATIONAL**       lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_RATIONAL**       rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_RATIONAL**       val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets column names */
SCIP_RETCODE SCIPlpiExactGetColNames(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get name from LP */
   int                   lastcol,            /**< last column to get name from LP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for col names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( colnames != NULL || namestoragesize == 0 );
   assert( namestorage != NULL || namestoragesize == 0 );
   assert( namestoragesize >= 0 );
   assert( storageleft != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets row names */
SCIP_RETCODE SCIPlpiExactGetRowNames(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get name from LP */
   int                   lastrow,            /**< last row to get name from LP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for row names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL);
   assert( rownames != NULL || namestoragesize == 0 );
   assert( namestorage != NULL || namestoragesize == 0 );
   assert( namestoragesize >= 0 );
   assert( storageleft != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets objective sense of the LP */
SCIP_RETCODE SCIPlpiExactGetObjsen(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiExactGetObj(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_RATIONAL**       vals                /**< array to store objective coefficients */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( firstcol <= lastcol );
   assert( vals != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiExactGetBounds(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective value for */
   int                   lastcol,            /**< last column to get objective value for */
   SCIP_RATIONAL**       lbs,                /**< array to store lower bound values, or NULL */
   SCIP_RATIONAL**       ubs                 /**< array to store upper bound values, or NULL */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( firstcol <= lastcol );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiExactGetSides(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_RATIONAL**       lhss,               /**< array to store left hand side values, or NULL */
   SCIP_RATIONAL**       rhss                /**< array to store right hand side values, or NULL */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( firstrow <= lastrow );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiExactGetCoef(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_RATIONAL*        val                 /**< pointer to store the value of the coefficient */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( val != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/**@} */


/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiExactSolvePrimal(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiExactSolveDual(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiExactSolveBarrier(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** start strong branching - call before any strongbranching */
SCIP_RETCODE SCIPlpiExactStartStrongbranch(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** end strong branching - call after any strongbranching */
SCIP_RETCODE SCIPlpiExactEndStrongbranch(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/**@} */


/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiExactWasSolved(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** gets information about primal and dual feasibility of the current LP solution
 *
 *  The feasibility information is with respect to the last solving call and it is only relevant if SCIPlpiWasSolved()
 *  returns true. If the LP is changed, this information might be invalidated.
 *
 *  Note that @param primalfeasible and @param dualfeasible should only return true if the solver has proved the
 *  respective LP to be feasible. Thus, the return values should be equal to the values of SCIPlpiIsPrimalFeasible() and
 *  SCIPlpiIsDualFeasible(), respectively. Note that if feasibility cannot be proved, they should return false (even if
 *  the problem might actually be feasible).
 */
SCIP_RETCODE SCIPlpiExactGetSolFeasibility(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< pointer to store primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< pointer to store dual feasibility status */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( primalfeasible != NULL );
   assert( dualfeasible != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExactExistsPrimalRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExactHasPrimalRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiExactIsPrimalUnbounded(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiExactIsPrimalInfeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiExactIsPrimalFeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExactExistsDualRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExactHasDualRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff LP is dual unbounded */
SCIP_Bool SCIPlpiExactIsDualUnbounded(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff LP is dual infeasible */
SCIP_Bool SCIPlpiExactIsDualInfeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiExactIsDualFeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiExactIsOptimal(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiExactIsObjlimExc(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiExactIsIterlimExc(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiExactIsTimelimExc(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** returns the internal solution status of the solver */
int SCIPlpiExactGetInternalStatus(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiExactIgnoreInstability(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( success != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiExactGetObjval(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        objval              /**< stores the objective value */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( objval != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}


/** gets primal and dual solution vectors for feasible LPs
 *
 *  Before calling this function, the caller must ensure that the LP has been solved to optimality, i.e., that
 *  SCIPlpiIsOptimal() returns true.
 */
SCIP_RETCODE SCIPlpiExactGetSol(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_RATIONAL**       primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_RATIONAL**       dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_RATIONAL**       activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_RATIONAL**       redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}


/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiExactGetPrimalRay(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL**       ray                 /**< primal ray */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets dual farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiExactGetDualfarkas(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL**       dualfarkas          /**< dual farkas row multipliers */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiExactGetIterations(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/**@} */


/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */


/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiExactGetBase(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiExactSetBase(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( cstat != NULL );
   assert( rstat != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
SCIP_RETCODE SCIPlpiExactGetBasisInd(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  bind                /**< pointer to store basis indices ready to keep number of rows entries */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( bind != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** get row of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiExactGetBInvRow(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_RATIONAL**       coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( coef != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** get column of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
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
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( coef != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/**@} */


/*
 * LP State Methods
 */

/**@name LP State Methods */
/**@{ */

/** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiExactGetState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( blkmem != NULL );
   assert( lpistate != NULL );
   assert( blkmem != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiExactGetState()
 */
SCIP_RETCODE SCIPlpiExactSetState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information), or NULL */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( blkmem != NULL );
   assert( lpistate != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiExactClearState(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   return SCIP_OKAY;
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiExactFreeState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( lpistate != NULL );
   assert( blkmem != NULL );
   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiExactHasStateBasis(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information), or NULL */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessageAbort();
   return FALSE;
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiExactReadState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( fname != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** writes LPi state (i.e. basis information) to a file */
SCIP_RETCODE SCIPlpiExactWriteState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( fname != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/**@} */


/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiExactGetIntpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( ival != NULL );
   return SCIP_PARAMETERUNKNOWN;
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiExactSetIntpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   return SCIP_PARAMETERUNKNOWN;
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiExactGetRealpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( dval != NULL );
   return SCIP_PARAMETERUNKNOWN;
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiExactSetRealpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   return SCIP_PARAMETERUNKNOWN;
}

/**@} */



/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiExactInfinity(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   return LPIINFINITY;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiExactIsInfinity(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< the value */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   if( val >= LPIINFINITY )
      return TRUE;
   return FALSE;
}

/**@} */


/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
SCIP_RETCODE SCIPlpiExactReadLP(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( fname != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiExactWriteLP(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/**@} */
