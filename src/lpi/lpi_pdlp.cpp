#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra-semi"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "Eigen/Core"
#include "Eigen/SparseCore"
#include "ortools/base/init_google.h"
#include "ortools/pdlp/iteration_stats.h"
#include "ortools/pdlp/primal_dual_hybrid_gradient.h"
#include "ortools/pdlp/quadratic_program.h"
#include "ortools/pdlp/solve_log.pb.h"
#include "ortools/pdlp/solvers.pb.h"
#include "ortools/util/file_util.h"

#pragma GCC diagnostic pop


#include "lpi/lpi.h"
#include "scip/pub_message.h"

#include <assert.h>
#include <limits>

namespace pdlp = ::operations_research::pdlp;

struct SCIP_LPi
{
   pdlp::QuadraticProgram *linear_program;
   pdlp::PrimalDualHybridGradientParams *parameters;
   pdlp::SolverResult *result;
   SCIP_Bool lp_modified_since_last_solve;
   SCIP_OBJSEN objsen;
   // TODO: Do we need scaling
};

// there are expected to be blank since we have no crossover.
struct SCIP_LPiState
{
   
};

static const std::string pdlp_name = "PDLP";

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   return pdlp_name.c_str();
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return "PDLP Linear Solver, developed by Google, part of OR-Tools (developers.google.com/optimization)";
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   assert( lpi != NULL );

   SCIPerrorMessage("SCIPlpiGetSolverPointer() has not been implemented yet.\n");

   return NULL;
}

SCIP_RETCODE SCIPlpiSetIntegralityInformation(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   int                   ncols,              /**< length of integrality array */
   int*                  intInfo             /**< integrality array (0: continuous, 1: integer). May be NULL iff ncols is 0.  */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( ncols == 0 || ncols == lpi->linear_program->num_variables().value() );

   return SCIP_OKAY;
}

/** informs about availability of a primal simplex solving method */
SCIP_Bool SCIPlpiHasPrimalSolve(
   void
   )
{
   return FALSE;
}

/** informs about availability of a dual simplex solving method */
SCIP_Bool SCIPlpiHasDualSolve(
   void
   )
{
   return FALSE;
}

/** informs about availability of a barrier solving method */
SCIP_Bool SCIPlpiHasBarrierSolve(
   void
   )
{
   return TRUE;
}

/** creates an LP problem object */
SCIP_RETCODE SCIPlpiCreate(
   SCIP_LPI**            lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   assert( lpi != NULL );
   assert( name != NULL );

   /* Initilialize memory. */
   SCIP_ALLOC(BMSallocMemory(lpi));
   (*lpi)->linear_program = new pdlp::QuadraticProgram();
   (*lpi)->parameters = new pdlp::PrimalDualHybridGradientParams();
   (*lpi)->result = new pdlp::SolverResult();
   (*lpi)->lp_modified_since_last_solve = TRUE;
   (*lpi)->objsen = SCIP_OBJSEN_MINIMIZE;

//    /* Set problem name and objective direction. */
//    (*lpi)->linear_program->problem_name = std::make_optional(std::string(name));
//    // Note: probably use ApplyObjectiveScalingAndOffset here
//    SCIP_CALL( SCIPlpiChgObjsen(*lpi, objsen) );

//    (*lpi)->from_scratch = false;
//    (*lpi)->lp_info = false;
//    (*lpi)->pricing = SCIP_PRICING_LPIDEFAULT;
//    (*lpi)->lp_modified_since_last_solve = true;
//    (*lpi)->lp_time_limit_was_reached = false;
//    (*lpi)->conditionlimit = -1.0;
//    (*lpi)->checkcondition = false;
//    (*lpi)->niterations = 0LL;

//    (*lpi)->tmp_row = new ScatteredRow();
//    (*lpi)->tmp_column = new ScatteredColumn();

// looks like NOSCALING is an option
// #ifdef NOSCALING
//    (*lpi)->parameters->set_use_scaling(false);
// #endif

   return SCIP_OKAY;
}

SCIP_EXPORT
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   delete (*lpi)->linear_program;
   delete (*lpi)->parameters;
   delete (*lpi)->result;

   BMSfreeMemory(lpi);

   return SCIP_OKAY;
}



/** copies LP data with column matrix into LP solver */
SCIP_RETCODE SCIPlpiLoadColLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   const SCIP_Real*      obj,                /**< objective function values of columns */
   const SCIP_Real*      lb,                 /**< lower bounds of columns */
   const SCIP_Real*      ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand sides of rows */
   const SCIP_Real*      rhs,                /**< right hand sides of rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array */
   const int*            ind,                /**< row indices of constraint matrix entries */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( obj != NULL );
   assert( lb != NULL );
   assert( ub != NULL );
   assert( beg != NULL );
   assert( ind != NULL );
   assert( val != NULL );

   lpi->linear_program->ResizeAndInitialize(0, 0);
   SCIP_CALL( SCIPlpiAddRows(lpi, nrows, lhs, rhs, rownames, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddCols(lpi, ncols, obj, lb, ub, colnames, nnonz, beg, ind, val) );
   SCIP_CALL( SCIPlpiChgObjsen(lpi, objsen) );

   return SCIP_OKAY;
}

/** adds columns to the LP */
SCIP_RETCODE SCIPlpiAddCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   const SCIP_Real*      obj,                /**< objective function values of new columns */
   const SCIP_Real*      lb,                 /**< lower bounds of new columns */
   const SCIP_Real*      ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   // looks like CSC format to me
   const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( obj != NULL );
   assert( lb != NULL );
   assert( ub != NULL );
   assert( nnonz >= 0) ;
   assert( ncols >= 0) ;

   SCIPdebugMessage("adding %d columns with %d nonzeros.\n", ncols, nnonz);

   /* @todo add names */
   if ( nnonz > 0 )
   {
      assert( beg != NULL );
      assert( ind != NULL );
      assert( val != NULL );
      assert( ncols > 0 );

#ifndef NDEBUG
      /* perform check that no new rows are added */
      size_t num_rows = lpi->linear_program->constraint_matrix.nrows();
      for (int j = 0; j < nnonz; ++j)
      {
         assert( 0 <= ind[j] && ind[j] < num_rows );
         assert( val[j] != 0.0 );
      }
#endif

      int nz = 0;
      for (int i = 0; i < ncols; ++i)
      {
         const int end = (nnonz == 0 || i == ncols - 1) ? nnonz : beg[i + 1];
         while ( nz < end )
         {
            // lpi->linear_program->constraint_matrix.
            // https://stackoverflow.com/a/41777589
            // lpi->linear_program->SetCoefficient(RowIndex(ind[nz]), col, val[nz]);
            ++nz;
         }
      }
      assert( nz == nnonz );
   }

   Eigen::VectorXd new_obj = Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(obj, ncols);
   Eigen::VectorXd new_variable_lower_bounds = Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(lb, ncols); 
   Eigen::VectorXd new_variable_upper_bounds = Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(ub, ncols);
   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}


SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );

   SCIPdebugMessage("SCIPlpiClear\n");

   lpi->linear_program->ResizeAndInitialize(0, 0);
   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}

/** changes lower and upper bounds of columns */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiChgBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices or NULL if ncols is zero */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds or NULL if ncols is zero */
   const SCIP_Real*      ub                  /**< values for the new upper bounds or NULL if ncols is zero */
   )
{
   for(int i = 0; i < ncols; ++i) {
      int col_idx = ind[i];
      lpi->linear_program->variable_lower_bounds(col_idx) = lb[i];
      lpi->linear_program->variable_upper_bounds(col_idx) = ub[i];
   }

   return SCIP_OKAY;
}

/** changes left and right hand sides of rows */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiChgSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      lhs,                /**< new values for left hand sides */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   for (int i = 0; i < nrows; ++i) {
      int row_idx = ind[i];
      lpi->linear_program->constraint_lower_bounds(row_idx) = lhs[i];
      lpi->linear_program->constraint_upper_bounds(row_idx) = rhs[i];
   }

   return SCIP_OKAY;
}

/** changes a single coefficient */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiChgCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   lpi->linear_program->constraint_matrix.coeffRef(row, col) = newval;
   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   if (lpi->objsen == objsen) {
      return SCIP_OKAY;
   } else {
      for (int i = 0; i < lpi->linear_program->objective_vector.size(); ++i)
      {
         lpi->linear_program->objective_vector(i) = -lpi->linear_program->objective_vector(i);
      }
      if (lpi->objsen == SCIP_OBJSEN_MAXIMIZE) {
         assert(lpi->objsen == SCIP_OBJSEN_MINIMIZE);
         lpi->objsen = SCIP_OBJSEN_MAXIMIZE;
      } else {
         assert(lpi->objsen == SCIP_OBJSEN_MAXIMIZE);
         lpi->objsen = SCIP_OBJSEN_MINIMIZE;
      }
      return SCIP_OKAY;
   }
}

/** changes objective values of columns in the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiChgObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   const int*            ind,                /**< column indices to change objective value for */
   const SCIP_Real*      obj                 /**< new objective values for columns */
   )
{
   for(int i = 0; i < ncols; ++i) {
      int col_idx = ind[i];
      lpi->linear_program->objective_vector(col_idx) = obj[i];
   }

   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiScaleRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   if (scaleval < 0) {
      double temp = lpi->linear_program->constraint_lower_bounds(row);
      lpi->linear_program->constraint_lower_bounds(row) = lpi->linear_program->constraint_upper_bounds(row);
      lpi->linear_program->constraint_upper_bounds(row) = temp;
   }

   lpi->linear_program->constraint_lower_bounds(row) *= scaleval;
   lpi->linear_program->constraint_upper_bounds(row) *= scaleval;

   // this is quite slow because it's CSC format, not CSR formt.
   for (int col_idx = 0; col_idx < lpi->linear_program->constraint_matrix.cols(); ++col_idx) {
      if (lpi->linear_program->constraint_matrix.coeff(row, col_idx) != 0.0) {
         lpi->linear_program->constraint_matrix.coeffRef(row, col_idx) *= scaleval;
      }
   }
   return SCIP_OKAY;
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiScaleCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   lpi->linear_program->objective_vector(col) *= scaleval;

   if (scaleval < 0) {
      double temp = lpi->linear_program->variable_lower_bounds(col);
      lpi->linear_program->variable_upper_bounds(col) = lpi->linear_program->variable_lower_bounds(col);
      lpi->linear_program->variable_upper_bounds(col) = temp;
   }

   lpi->linear_program->variable_lower_bounds(col) /= scaleval;
   lpi->linear_program->variable_upper_bounds(col) /= scaleval;
   return SCIP_OKAY;
}

/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetNRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{
   *nrows = lpi->linear_program->constraint_matrix.rows();
   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   *ncols = lpi->linear_program->constraint_matrix.cols();
   return SCIP_OKAY;
}

/** gets the objective sense of the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   *objsen = lpi->objsen;
   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   ) {
   *nnonz = lpi->linear_program->constraint_matrix.nonZeros();
   return SCIP_OKAY;
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   // seems to assume that we have a CSR matrix, while Eigen uses a CSC matrix, ugh...
   // yeah, this is a CSR matrix
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   // TODO: Transpose the CSC matrix to CSR
   assert(false);
   return SCIP_NOTIMPLEMENTED;
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_Real*            lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   assert(false);
   return SCIP_NOTIMPLEMENTED;
}

/** gets column names */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetColNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get name from LP */
   int                   lastcol,            /**< last column to get name from LP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for col names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( colnames != NULL || namestoragesize == 0 );
   assert( namestorage != NULL || namestoragesize == 0 );
   assert( namestoragesize >= 0 );
   assert( storageleft != NULL );
   assert( 0 <= firstcol && firstcol <= lastcol && lastcol < lpi->linear_program->num_variables() );

   SCIPerrorMessage("SCIPlpiGetColNames() has not been implemented yet.\n");

   return SCIP_NOTIMPLEMENTED;
}


/** gets row names */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetRowNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get name from LP */
   int                   lastrow,            /**< last row to get name from LP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for row names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( rownames != NULL || namestoragesize == 0 );
   assert( namestorage != NULL || namestoragesize == 0 );
   assert( namestoragesize >= 0 );
   assert( storageleft != NULL );
   assert( 0 <= firstrow && firstrow <= lastrow && lastrow < lpi->linear_program->num_constraints() );

   SCIPerrorMessage("SCIPlpiGetRowNames() has not been implemented yet.\n");

   return SCIP_NOTIMPLEMENTED;
}


/** gets objective coefficients from LP problem object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( firstcol <= lastcol );

   SCIPdebugMessage("getting row sides %d to %d\n", firstcol, lastcol);

   int index = 0;
   for (int col = firstcol; col <= lastcol; ++col)
   {
      vals[index] = lpi->linear_program->objective_vector(col);
      ++index;
   }

   return SCIP_OKAY;
}

/** gets current bounds from LP problem object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get bounds for */
   int                   lastcol,            /**< last column to get bounds for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( firstcol <= lastcol );

   SCIPdebugMessage("getting row sides %d to %d\n", firstcol, lastcol);

   int index = 0;
   for (int col = firstcol; col <= lastcol; ++col)
   {
      if ( lbs != NULL )
         lbs[index] = lpi->linear_program->constraint_lower_bounds(col);

      if ( ubs != NULL )
         ubs[index] = lpi->linear_program->constraint_upper_bounds(col);

      ++index;
   }

   return SCIP_OKAY;
}

/** gets current row sides from LP problem object */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( firstrow <= lastrow );

   SCIPdebugMessage("getting row sides %d to %d\n", firstrow, lastrow);

   int index = 0;
   for (int row = firstrow; row <= lastrow; ++row)
   {
      if ( lhss != NULL )
         lhss[index] = lpi->linear_program->constraint_lower_bounds(row);

      if ( rhss != NULL )
         rhss[index] = lpi->linear_program->constraint_upper_bounds(row);

      ++index;
   }

   return SCIP_OKAY;
}

/** gets a single coefficient */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   *val = lpi->linear_program->constraint_matrix.coeff(row, col);
   return SCIP_OKAY;
}

/**@} */



/** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiGetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpistate != NULL);

   *lpistate = NULL;
   return SCIP_OKAY;
}

SCIP_EXPORT
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   ) {
   *lpi->result = pdlp::PrimalDualHybridGradient(*lpi->linear_program, *lpi->parameters);
   assert(!crossover);
   return SCIP_OKAY;
}

/** start strong branching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );

   /* @todo Save state and do all the branching from there. */
   return SCIP_OKAY;
}

/** end strong branching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );

   /* @todo Restore the saved state. */
   return SCIP_OKAY;
}



/** performs strong branching iterations on given @b fractional candidates */
SCIP_RETCODE SCIPlpiStrongbranchesFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< fractional current primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL) ;
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );

   SCIPerrorMessage("SCIPlpiStrongbranchesFrac - not implemented.\n");

   return SCIP_NOTIMPLEMENTED;
}

/** performs strong branching iterations on given candidates with @b integral values */
SCIP_RETCODE SCIPlpiStrongbranchesInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< current integral primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL) ;
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );

   SCIPerrorMessage("SCIPlpiStrongbranchesInt - not implemented.\n");

   return SCIP_NOTIMPLEMENTED;
}

SCIP_Bool SCIPlpiWasSolved(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );

   /* @todo Track this to avoid uneeded resolving. */
   return ( ! lpi->lp_modified_since_last_solve );
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->result != NULL );

   return lpi->result->solve_log.termination_reason() == pdlp::TERMINATION_REASON_ITERATION_LIMIT;
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );

   return static_cast<int>(lpi->result->solve_log.termination_reason());
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( success != NULL );

   *success = FALSE;

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   assert( lpi != NULL );
   assert( objval != NULL );

   // TODO: How do I know which one I want?
   *objval = lpi->result->solve_log.solution_stats().convergence_information()[0].primal_objective();

   return SCIP_OKAY;
}

/** gets the number of LP iterations of the last solve call */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   *iterations = lpi->result->solve_log.iteration_count();
   return SCIP_OKAY;
}


SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetRealSolQuality(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPSOLQUALITY     qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   *quality = SCIP_INVALID;
   return SCIP_OKAY;
}


/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */

// Need to figure out how this interacts with a barrier method.
/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   ) {
   return SCIP_ERROR;
}

/** sets current basis status for columns and rows */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiSetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const int*            cstat,              /**< array with column basis status */
   const int*            rstat               /**< array with row basis status */
   ) {
   return SCIP_ERROR;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< pointer to store basis indices ready to keep number of rows entries */
   ) {
   return SCIP_ERROR;
}

/** get row of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetBInvRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   ) {
   return SCIP_ERROR;
}

/** get column of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetBInvCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   ) {
   return SCIP_ERROR;
}

/** get row of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            coef,               /**< vector to return coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   ) {
   return SCIP_ERROR;
}

/** get column of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef,               /**< vector to return coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   ) {
   return SCIP_ERROR;
}

/**@} */

/*
 * LPi State Methods
 */

/**@name LPi State Methods */
/**@{ */

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetState()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiSetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const SCIP_LPISTATE*  lpistate            /**< LPi state information (like basis information), or NULL */
   ) {
   return SCIP_ERROR;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   ) {
   return SCIP_ERROR;
}

/** frees LPi state information */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiFreeState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   ) {
   return SCIP_ERROR;
}

/** checks, whether the given LPi state contains simplex basis information */
SCIP_EXPORT
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   ) {
   return SCIP_ERROR;
}

/** reads LPi state (like basis information from a file */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   ) {
   return SCIP_ERROR;
}

/** writes LPi state (i.e. basis information) to a file */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   ) {
   return SCIP_ERROR;
}

/**@} */


/*
 * LP Pricing Norms Methods
 */

/**@name LP Pricing Norms Methods */
/**@{ */

/* SCIP_LPiNorms stores norm information so they are not recomputed from one state to the next. */
/* @todo Implement this. */
struct SCIP_LPiNorms {};

/** stores LPi pricing norms information
 *
 *  @todo store primal norms as well?
 */
SCIP_RETCODE SCIPlpiGetNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information */
   )
{
   assert( lpi != NULL );
   assert( blkmem != NULL );
   assert( lpi->solver != NULL );
   assert( lpinorms != NULL );

   return SCIP_OKAY;
}

/** loads LPi pricing norms into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetNorms()
 */
SCIP_RETCODE SCIPlpiSetNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const SCIP_LPINORMS*  lpinorms            /**< LPi pricing norms information, or NULL */
   )
{
   assert( lpi != NULL );
   assert( blkmem != NULL );
   assert( lpi->solver != NULL );
   assert( lpinorms != NULL );

   return SCIP_OKAY;
}

/** frees pricing norms information */
SCIP_RETCODE SCIPlpiFreeNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information, or NULL */
   )
{
   assert( lpi != NULL );
   assert( blkmem != NULL );
   assert( lpi->solver != NULL );
   assert( lpinorms != NULL );

   return SCIP_OKAY;
}

/**@} */


/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the LP solver */
SCIP_EXPORT
SCIP_Real SCIPlpiInfinity(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   return std::numeric_limits<double>::infinity();
}

/** interrupts the currently ongoing lp solve or disables the interrupt */
SCIP_RETCODE SCIPlpiInterrupt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             interrupt           /**< TRUE if interrupt should be set, FALSE if it should be disabled */
   )
{
   /*lint --e{715}*/
   assert(lpi != NULL);

   return SCIP_OKAY;
}


/** checks if given value is treated as infinity in the LP solver */
SCIP_EXPORT
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{
   return val == std::numeric_limits<double>::infinity();
}

/** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( fname != NULL );

   const std::string filespec(fname);
   operations_research::MPModelProto proto;
   if ( ! ReadFileToProto(filespec, &proto) )
   {
      SCIPerrorMessage("Could not read <%s>\n", fname);
      return SCIP_READERROR;
   }
   absl::StatusOr<pdlp::QuadraticProgram> mp = pdlp::QpFromMpModelProto(proto, true, false);
   CHECK_OK(mp.status());
   *lpi->linear_program = *mp;

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( fname != NULL );

   

   absl::StatusOr<operations_research::MPModelProto> mp = pdlp::QpToMpModelProto(*lpi->linear_program);
   CHECK_OK(mp.status());
   const std::string filespec(fname);
   if ( ! WriteProtoToFile(filespec, *mp, operations_research::ProtoWriteFormat::kProtoText, true) )
   {
      SCIPerrorMessage("Could not write <%s>\n", fname);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}
