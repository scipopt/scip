#include "ortools/pdlp/quadratic_program.h"

struct SCIP_LPi
{
   pdlp::QuadraticProgram *linear_program;
   pdlp::PrimalDualHybridGradientParams *parameters;
   // TODO: Do we need scaling
};

struct SCIP_LPiState
{
   
};

struct SCIP_LPiNorms
{
};


static std::string pdlp_name = "PDLP";

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   return pdlp_name;
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
   TODO_GALVEZ;
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( ncols == 0 || ncols == lpi->linear_program->num_variables().value() );

   // /* Pass on integrality information (currently not used by Glop). */
   // for (ColIndex col(0); col < ColIndex(ncols); ++col)
   // {
   //    assert( intInfo != NULL );
   //    int info = intInfo[col.value()];
   //    assert( info == 0 || info == 1 );
   //    if ( info == 0 )
   //       lpi->linear_program->SetVariableType(col, operations_research::glop::LinearProgram::VariableType::CONTINUOUS);
   //    else
   //       lpi->linear_program->SetVariableType(col, operations_research::glop::LinearProgram::VariableType::INTEGER);
   // }

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
   (*lpi)->parameters = new pdlp::GlopParameters();

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

// #ifdef NOSCALING
//    (*lpi)->parameters->set_use_scaling(false);
// #endif

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

   lpi->linear_program->Clear();
   SCIP_CALL( SCIPlpiAddRows(lpi, nrows, lhs, rhs, rownames, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddCols(lpi, ncols, obj, lb, ub, colnames, nnonz, beg, ind, val) );
   SCIP_CALL( SCIPlpiChgObjsen(lpi, objsen) );

   return SCIP_OKAY;
}

/** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiGetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpistate != NULL);

   *lpistate = NULL;
   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(*lpi != NULL);

   SCIPdebugMessage("SCIPlpiFree()\n");

   delete (*lpi)->linear_program;
   delete (*lpi)->parameters;

   BMSfreeMemory(lpi);

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

   

   absl::StatusOr<MPModelProto> mp = pdlp::QpToMpModelProto(const QuadraticProgram& qp);
   const std::string filespec(fname);
   if ( ! WriteProtoToFile(filespec, mp, operations_research::ProtoWriteFormat::kProtoText, true) )
   {
      SCIPerrorMessage("Could not write <%s>\n", fname);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}
