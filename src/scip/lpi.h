/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lpi.h
 * @brief  interface methods for specific LP solvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __LPI_H__
#define __LPI_H__


typedef struct LPi LPI;                 /**< solver dependent LP interface */
typedef struct LPiState LPISTATE;       /**< complete LP state (i.e. basis information, dual norms) */

/** objective sense */
enum ObjSen
{
   SCIP_OBJSEN_MAXIMIZE = -1,           /**< maximize objective function */
   SCIP_OBJSEN_MINIMIZE = +1            /**< minimize objective function */
};
typedef enum ObjSen OBJSEN;

/** LP solver parameters */
enum LPParam
{
   SCIP_LPPAR_FROMSCRATCH =  0,         /**< solver should start from scratch at next call */
   SCIP_LPPAR_FASTMIP     =  1,         /**< fast mip setting of LP solver */
   SCIP_LPPAR_PRICING     =  2,         /**< pricing strategy */
   SCIP_LPPAR_LPINFO      =  3,         /**< should LP solver output information to the screen? */
   SCIP_LPPAR_FEASTOL     =  4,         /**< feasibility tolerance */
   SCIP_LPPAR_LOBJLIM     =  5,         /**< lower objective limit */
   SCIP_LPPAR_UOBJLIM     =  6,         /**< upper objective limit */
   SCIP_LPPAR_LPITLIM     =  7,         /**< LP iteration limit */
   SCIP_LPPAR_LPTILIM     =  8,         /**< LP time limit */
   SCIP_LPPAR_LPITER      =  9          /**< number of simplex iterations for last LP solve */
};
typedef enum LPParam LPPARAM;

/** LP pricing strategy */
enum Pricing
{
   SCIP_PRICING_AUTO        = 0,        /**< the LP solver should use its prefered strategy */
   SCIP_PRICING_FULL        = 1,        /**< full pricing */
   SCIP_PRICING_STEEP       = 2,        /**< steepest edge pricing */
   SCIP_PRICING_STEEPQSTART = 3         /**< steepest edge pricing without initial dual norms */
};
typedef enum Pricing PRICING;

/** basis status for columns and rows */
enum BaseStat
{
   SCIP_BASESTAT_LOWER = 0,             /**< (slack) variable is at its lower bound */
   SCIP_BASESTAT_BASIC = 1,             /**< (slack) variable is basic */
   SCIP_BASESTAT_UPPER = 2              /**< (slack) variable is at its upper bound */
};
typedef enum BaseStat BASESTAT;




#include "def.h"
#include "memory.h"
#include "retcode.h"




/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
extern
const char* SCIPlpiGetSolverName(
   void
   );

/**@} */




/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** creates an LP problem object */
extern 
RETCODE SCIPlpiCreate(
   LPI**            lpi,                /**< pointer to an LP interface structure */
   const char*      name                /**< problem name */
   );

/** deletes an LP problem object */
extern
RETCODE SCIPlpiFree(
   LPI**            lpi                 /**< pointer to an LP interface structure */
   );

/**@} */




/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies LP data with column matrix into LP solver */
extern
RETCODE SCIPlpiLoadColLP(
   LPI*             lpi,                /**< LP interface structure */
   OBJSEN           objsen,             /**< objective sense */
   int              ncols,              /**< number of columns */
   const Real*      obj,                /**< objective function values of columns */
   const Real*      lb,                 /**< lower bounds of columns */
   const Real*      ub,                 /**< upper bounds of columns */
   char**           colnames,           /**< column names, or NULL */
   int              nrows,              /**< number of rows */
   const Real*      lhs,                /**< left hand sides of rows */
   const Real*      rhs,                /**< right hand sides of rows */
   char**           rownames,           /**< row names, or NULL */
   int              nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*       beg,                /**< start index of each column in ind- and val-array */
   const int*       ind,                /**< row indices of constraint matrix entries */
   const Real*      val                 /**< values of constraint matrix entries */
   );

/** adds columns to the LP */
extern 
RETCODE SCIPlpiAddCols(
   LPI*             lpi,                /**< LP interface structure */
   int              ncols,              /**< number of columns to be added */
   const Real*      obj,                /**< objective function values of new columns */
   const Real*      lb,                 /**< lower bounds of new columns */
   const Real*      ub,                 /**< upper bounds of new columns */
   char**           colnames,           /**< column names, or NULL */
   int              nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*       beg,                /**< start index of each column in ind- and val-array */
   const int*       ind,                /**< row indices of constraint matrix entries */
   const Real*      val                 /**< values of constraint matrix entries */
   );

/** deletes all columns in the given range from LP */
extern
RETCODE SCIPlpiDelCols(
   LPI*             lpi,                /**< LP interface structure */
   int              firstcol,           /**< first column to be deleted */
   int              lastcol             /**< last column to be deleted */
   );

/** deletes columns from LP; the new position of a column must not be greater that its old position */
extern 
RETCODE SCIPlpiDelColset(
   LPI*             lpi,                /**< LP interface structure */
   int*             dstat               /**< deletion status of columns
                                         *   input:  1 if column should be deleted, 0 if not
                                         *   output: new position of column, -1 if column was deleted */
   );

/** adds rows to the LP */
extern 
RETCODE SCIPlpiAddRows(
   LPI*             lpi,                /**< LP interface structure */
   int              nrows,              /**< number of rows to be added */
   const Real*      lhs,                /**< left hand sides of new rows */
   const Real*      rhs,                /**< right hand sides of new rows */
   char**           rownames,           /**< row names, or NULL */
   int              nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*       beg,                /**< start index of each row in ind- and val-array */
   const int*       ind,                /**< column indices of constraint matrix entries */
   const Real*      val                 /**< values of constraint matrix entries */
   );

/** deletes all rows in the given range from LP */
extern
RETCODE SCIPlpiDelRows(
   LPI*             lpi,                /**< LP interface structure */
   int              firstrow,           /**< first row to be deleted */
   int              lastrow             /**< last row to be deleted */
   );

/** deletes rows from LP; the new position of a row must not be greater that its old position */
extern 
RETCODE SCIPlpiDelRowset(
   LPI*             lpi,                /**< LP interface structure */
   int*             dstat               /**< deletion status of rows
                                         *   input:  1 if row should be deleted, 0 if not
                                         *   output: new position of row, -1 if row was deleted */
   );

/** clears the whole LP */
extern
RETCODE SCIPlpiClear(
   LPI*             lpi                 /**< LP interface structure */
   );

/** changes lower and upper bounds of columns */
extern 
RETCODE SCIPlpiChgBounds(
   LPI*             lpi,                /**< LP interface structure */
   int              ncols,              /**< number of columns to change bounds for */
   const int*       ind,                /**< column indices */
   const Real*      lb,                 /**< values for the new lower bounds */
   const Real*      ub                  /**< values for the new upper bounds */
   );

/** changes left and right hand sides of rows */
extern 
RETCODE SCIPlpiChgSides(
   LPI*             lpi,                /**< LP interface structure */
   int              nrows,              /**< number of rows to change sides for */
   const int*       ind,                /**< row indices */
   const Real*      lhs,                /**< new values for left hand sides */
   const Real*      rhs                 /**< new values for right hand sides */
   );

/** changes a single coefficient */
extern
RETCODE SCIPlpiChgCoef(
   LPI*             lpi,                /**< LP interface structure */
   int              row,                /**< row number of coefficient to change */
   int              col,                /**< column number of coefficient to change */
   Real             newval              /**< new value of coefficient */
   );

/** changes the objective sense */
extern 
RETCODE SCIPlpiChgObjsen(
   LPI*             lpi,                /**< LP interface structure */
   OBJSEN           objsen              /**< new objective sense */
   );

/** changes objective values of columns in the LP */
extern
RETCODE SCIPlpiChgObj(
   LPI*             lpi,                /**< LP interface structure */
   int              ncols,              /**< number of columns to change objective value for */
   int*             ind,                /**< column indices to change objective value for */
   Real*            obj                 /**< new objective values for columns */
   );

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
extern
RETCODE SCIPlpiScaleRow(
   LPI*             lpi,                /**< LP interface structure */
   int              row,                /**< row number to scale */
   Real             scaleval            /**< scaling multiplier */
   );

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
extern
RETCODE SCIPlpiScaleCol(
   LPI*             lpi,                /**< LP interface structure */
   int              col,                /**< column number to scale */
   Real             scaleval            /**< scaling multiplier */
   );

/**@} */




/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
extern
RETCODE SCIPlpiGetNRows(
   LPI*             lpi,                /**< LP interface structure */
   int*             nrows               /**< pointer to store the number of rows */
   );

/** gets the number of columns in the LP */
extern
RETCODE SCIPlpiGetNCols(
   LPI*             lpi,                /**< LP interface structure */
   int*             ncols               /**< pointer to store the number of cols */
   );

/** gets the number of nonzero elements in the LP constraint matrix */
extern
RETCODE SCIPlpiGetNNonz(
   LPI*             lpi,                /**< LP interface structure */
   int*             nnonz               /**< pointer to store the number of nonzeros */
   );

/** gets columns from LP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
extern
RETCODE SCIPlpiGetCols(
   LPI*             lpi,                /**< LP interface structure */
   int              firstcol,           /**< first column to get from LP */
   int              lastcol,            /**< last column to get from LP */
   Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*             nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*             beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*             ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   );

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
extern
RETCODE SCIPlpiGetRows(
   LPI*             lpi,                /**< LP interface structure */
   int              firstrow,           /**< first row to get from LP */
   int              lastrow,            /**< last row to get from LP */
   Real*            lhs,                /**< buffer to store left hand side vector, or NULL */
   Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*             nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*             beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*             ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   );

/** gets objective values from LP problem object */
extern
RETCODE SCIPlpiGetObj(
   LPI*             lpi,                /**< LP interface structure */
   int              firstcol,           /**< first column to get objective value for */
   int              lastcol,            /**< last column to get objective value for */
   Real*            vals                /**< array to store objective values */
   );

/** gets a single coefficient */
extern
RETCODE SCIPlpiGetCoef(
   LPI*             lpi,                /**< LP interface structure */
   int              row,                /**< row number of coefficient */
   int              col,                /**< column number of coefficient */
   Real*            val                 /**< pointer to store the value of the coefficient */
   );

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** calls primal simplex to solve the LP */
extern 
RETCODE SCIPlpiSolvePrimal(
   LPI*             lpi                 /**< LP interface structure */
   );

/** calls dual simplex to solve the LP */
extern 
RETCODE SCIPlpiSolveDual(
   LPI*             lpi                 /**< LP interface structure */
   );

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
extern 
RETCODE SCIPlpiSolveBarrier(
   LPI*             lpi                 /**< LP interface structure */
   );

/** performs strong branching iterations on all candidates */
extern 
RETCODE SCIPlpiStrongbranch(
   LPI*             lpi,                /**< LP interface structure */
   const int*       cand,               /**< candidate list */
   Real*            psol,               /**< array with current primal solution values of candidates */
   int              ncand,              /**< size of candidate list */
   int              itlim,              /**< iteration limit for strong branchings */
   Real*            down,               /**< stores dual bounds after branching candidate down */
   Real*            up,                 /**< stores dual bounds after branching candidate up */
   int*             iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   );

/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** gets information about primal and dual feasibility of the LP basis */
extern 
RETCODE SCIPlpiGetBasisFeasibility(
   LPI*             lpi,                /**< LP interface structure */
   Bool*            primalfeasible,     /**< stores primal feasibility status */
   Bool*            dualfeasible        /**< stores dual feasibility status */
   );

/** returns TRUE iff LP is primal unbounded */
extern 
Bool SCIPlpiIsPrimalUnbounded(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is primal infeasible */
extern 
Bool SCIPlpiIsPrimalInfeasible(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is dual unbounded */
extern 
Bool SCIPlpiIsDualUnbounded(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is dual infeasible */
extern 
Bool SCIPlpiIsDualInfeasible(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP was solved to optimality */
extern 
Bool SCIPlpiIsOptimal(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff actual LP basis is stable */
extern 
Bool SCIPlpiIsStable(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the objective limit was reached */
extern 
Bool SCIPlpiIsObjlimExc(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the iteration limit was reached */
extern 
Bool SCIPlpiIsIterlimExc(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the time limit was reached */
extern 
Bool SCIPlpiIsTimelimExc(
   LPI*             lpi                 /**< LP interface structure */
   );

/** gets objective value of solution */
extern
RETCODE SCIPlpiGetObjval(
   LPI*             lpi,                /**< LP interface structure */
   Real*            objval              /**< stores the objective value */
   );

/** gets primal and dual solution vectors */
extern 
RETCODE SCIPlpiGetSol(
   LPI*             lpi,                /**< LP interface structure */
   Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
   Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   Real*            activity,           /**< row activity vector, may be NULL if not needed */
   Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   );

/** gets primal ray for unbounded LPs */
extern 
RETCODE SCIPlpiGetPrimalRay(
   LPI*             lpi,                /**< LP interface structure */
   Real*            ray                 /**< primal ray */
   );

/** gets dual farkas proof for infeasibility */
extern
RETCODE SCIPlpiGetDualfarkas(
   LPI*             lpi,                /**< LP interface structure */
   Real*            dualfarkas          /**< dual farkas row multipliers */
   );

/**@} */




/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */

/** gets actual basis status for columns and rows; arrays must be large enough to store the basis status */
extern
RETCODE SCIPlpiGetBase(
   LPI*             lpi,                /**< LP interface structure */
   int*             cstat,              /**< array to store column basis status, or NULL */
   int*             rstat               /**< array to store row basis status, or NULL */
   );

/** sets actual basis status for columns and rows */
extern
RETCODE SCIPlpiSetBase(
   LPI*             lpi,                /**< LP interface structure */
   int*             cstat,              /**< array with column basis status */
   int*             rstat               /**< array with row basis status */
   );

/** returns the indices of the basic columns and rows */
extern 
RETCODE SCIPlpiGetBasisInd(
   LPI*             lpi,                /**< LP interface structure */
   int*             bind                /**< basic column n gives value n, basic row m gives value -1-m */
   );

/** get dense row of inverse basis matrix B^-1 */
extern 
RETCODE SCIPlpiGetBInvRow(
   LPI*             lpi,                /**< LP interface structure */
   int              r,                  /**< row number */
   Real*            coef                /**< pointer to store the coefficients of the row */
   );

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
extern 
RETCODE SCIPlpiGetBInvARow(
   LPI*             lpi,                /**< LP interface structure */
   int              r,                  /**< row number */
   const Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   Real*            val                 /**< vector to return coefficients */
   );

/**@} */




/*
 * LP State Methods
 */

/**@name LP State Methods */
/**@{ */

/** stores LP state (like basis information) into lpistate object */
extern
RETCODE SCIPlpiGetState(
   LPI*             lpi,                /**< LP interface structure */
   MEMHDR*          memhdr,             /**< block memory */
   LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   );

/** loads LP state (like basis information) into solver */
extern
RETCODE SCIPlpiSetState(
   LPI*             lpi,                /**< LP interface structure */
   MEMHDR*          memhdr,             /**< block memory */
   LPISTATE*        lpistate            /**< LP state information (like basis information) */
   );

/** frees LP state information */
extern
RETCODE SCIPlpiFreeState(
   LPI*             lpi,                /**< LP interface structure */
   MEMHDR*          memhdr,             /**< block memory */
   LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   );

/** reads LP state (like basis information from a file */
extern 
RETCODE SCIPlpiReadState(
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   );

/** writes LP state (like basis information) to a file */
extern 
RETCODE SCIPlpiWriteState(
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   );

/**@} */




/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
extern 
RETCODE SCIPlpiGetIntpar(
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   int*             ival                /**< buffer to store the parameter value */
   );

/** sets integer parameter of LP */
extern 
RETCODE SCIPlpiSetIntpar(
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   int              ival                /**< parameter value */
   );

/** gets floating point parameter of LP */
extern 
RETCODE SCIPlpiGetRealpar(
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   Real*            dval                /**< buffer to store the parameter value */
   );

/** sets floating point parameter of LP */
extern 
RETCODE SCIPlpiSetRealpar(
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   Real             dval                /**< parameter value */
   );

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the LP solver */
extern
Real SCIPlpiInfinity(
   LPI*             lpi                 /**< LP interface structure */
   );

/** checks if given value is treated as infinity in the LP solver */
extern
Bool SCIPlpiIsInfinity(
   LPI*             lpi,                /**< LP interface structure */
   Real             val
   );

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
extern
RETCODE SCIPlpiReadLP(
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   );

/** writes LP to a file */
extern 
RETCODE SCIPlpiWriteLP(
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   );

/**@} */



#endif
