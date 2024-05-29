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

/**@file   lpi_highs.cpp
 * @ingroup LPIS
 * @brief  LP interface for HiGHS 1.4 and higher
 * @author Ambros Gleixner
 * @author Julian Hall
 * @author Alexander Hoen
 * @author Gioni Mexi
 *
 * This is an implementation of SCIP's LP interface for the open-source solver HiGHS.
 *
 * The most important open todos are:
 * - tune pricing strategy
 * - tune and activate primal simplex
 * - tune and activate parallel dual simplex
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>

/* undefine CMAKE_BUILD_TYPE in case it conflicts with HiGHS */
#ifdef CMAKE_BUILD_TYPE
#define SCIP_CMAKE_BUILD_TYPE (CMAKE_BUILD_TYPE)
#undef CMAKE_BUILD_TYPE
#endif


#include <Highs.h>

#include <lp_data/HighsLpUtils.h>

/* reset CMAKE_BUILD_TYPE to its original SCIP value */
#undef CMAKE_BUILD_TYPE
#ifdef SCIP_CMAKE_BUILD_TYPE
#define CMAKE_BUILD_TYPE (SCIP_CMAKE_BUILD_TYPE)
#undef SCIP_CMAKE_BUILD_TYPE
#endif

#include "lpi/lpi.h"
#include "scip/bitencode.h"
#include "scip/pub_message.h"
#include "scip/type_lp.h"

/* #define HIGHS_DEBUGLEVEL kHighsDebugLevelExpensive */
/* #define HIGHS_LOGDEVLEVEL kHighsLogDevLevelVerbose */

/*
 * Macros, structs, etc.
 */

#define HIGHS_relDiff(val1, val2)         ( ((val1)-(val2))/(MAX3(1.0,REALABS(val1),REALABS(val2))) )

/** Macro for a single HiGHS call for which exceptions have to be caught. We make no distinction between different
 *  exception types, e.g., between memory allocation and other exceptions. Additionally, we check if HiGHS returns kOk
 *  as status and return an LP error if not.
 */
#define HIGHS_CALL(x)   do                                                                                     \
                        {                                                                                      \
                           try                                                                                 \
                           {                                                                                   \
                              HighsStatus _restat_; /*lint -e{506,774}*/                                       \
                              (_restat_ = (x));                                                                \
                              if( _restat_ == HighsStatus::kWarning )                                          \
                              {                                                                                \
                                 SCIPerrorMessage("Warning in HiGHS function call\n");                         \
                                 return SCIP_LPERROR;                                                          \
                              }                                                                                \
                              else if( _restat_ != HighsStatus::kOk )                                          \
                              {                                                                                \
                                 SCIPerrorMessage("Error in HiGHS function call\n");                           \
                                 return SCIP_LPERROR;                                                          \
                              }                                                                                \
                           }                                                                                   \
                           catch( std::exception & E )                                                         \
                           {                                                                                   \
                              std::string s = E.what();                                                        \
                              SCIPerrorMessage( "HiGHS threw an exception: %s\n", s.c_str());                  \
                              return SCIP_LPERROR;                                                             \
                           }                                                                                   \
                           catch( ... )                                                                        \
                           {                                                                                   \
                              SCIPerrorMessage("HiGHS threw an unidentified exception\n");                     \
                              return SCIP_LPERROR;                                                             \
                           }                                                                                   \
                        }                                                                                      \
                        while( FALSE )

/** A relaxed version of HIGHS_CALL that accepts status kWarning. */
#define HIGHS_CALL_WITH_WARNING(x)   do                                                                        \
                        {                                                                                      \
                           try                                                                                 \
                           {                                                                                   \
                              HighsStatus _restat_; /*lint -e{506,774}*/                                       \
                              (_restat_ = (x));                                                                \
                              if( _restat_ != HighsStatus::kOk && _restat_ != HighsStatus::kWarning )          \
                              {                                                                                \
                                 SCIPerrorMessage("Error in HiGHS in function call (returned %d)\n",           \
                                    int(_restat_));                                                            \
                                 return SCIP_LPERROR;                                                          \
                              }                                                                                \
                           }                                                                                   \
                           catch( std::exception & E )                                                         \
                           {                                                                                   \
                              std::string s = E.what();                                                        \
                              SCIPerrorMessage( "HiGHS threw an exception: %s\n", s.c_str());                  \
                              return SCIP_LPERROR;                                                             \
                           }                                                                                   \
                           catch( ... )                                                                        \
                           {                                                                                   \
                              SCIPerrorMessage("HiGHS threw an unidentified exception\n");                     \
                              return SCIP_LPERROR;                                                             \
                           }                                                                                   \
                        }                                                                                      \
                        while( FALSE )

/**@todo make thread-safe */
int nsolvecalls = 0;

/** SCIP's HiGHS class */
class HighsSCIP : public Highs
{
   bool                  _lpinfo;
   char*                 _probname;
   SCIP_MESSAGEHDLR*     _messagehdlr; /**< messagehdlr handler for printing messages, or NULL */

public:

   HighsSCIP(
      SCIP_MESSAGEHDLR*  messagehdlr = NULL, /**< message handler */
      const char*        probname = NULL     /**< name of problem */
            )
      : _lpinfo(false),
        _probname(NULL),
        _messagehdlr(messagehdlr)
   {
      /* TODO set problem name by using an internal function */
   }

   virtual ~HighsSCIP()
   {
      /* TODO free problem name */
   }
};

/** LP interface struct for HiGHS */
struct SCIP_LPi
{
   HighsSCIP*            highs;              /**< HiGHS problem class */
   int*                  cstat;              /**< array for storing column basis status */
   int*                  rstat;              /**< array for storing row basis status */
   int                   cstatsize;          /**< size of cstat array */
   int                   rstatsize;          /**< size of rstat array */
   int                   nthreads;           /**< number of threads to be used */
   SCIP_Bool             fromscratch;        /**< shall solves be performed from scratch? */
   SCIP_Bool             solved;             /**< was the current LP solved? */
   SCIP_PRICING          pricing;            /**< SCIP pricing setting  */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehdlr handler for printing messages, or NULL */
};

typedef SCIP_DUALPACKET COLPACKET;           /** each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /** each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE

/** LPi state stores basis information */
struct SCIP_LPiState
{
   int                   ncols;              /**< number of LP columns */
   int                   nrows;              /**< number of LP rows */
   COLPACKET*            packcstat;          /**< column basis status in compressed form */
   ROWPACKET*            packrstat;          /**< row basis status in compressed form */
};

/*
 * dynamic memory arrays
 */

/** resizes cstat array to have at least num entries */
static
SCIP_RETCODE ensureCstatMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   SCIPdebugMessage("calling ensureCstatMem()\n");

   assert(lpi != NULL);

   if( num > lpi->cstatsize )
   {
      int newsize;
      newsize = MAX( 2 * lpi->cstatsize, num );
      SCIP_ALLOC( BMSreallocMemoryArray( &lpi->cstat, newsize ) );
      lpi->cstatsize = newsize;
   }
   assert(num <= lpi->cstatsize);

   return SCIP_OKAY;
}

/** resizes rstat array to have at least num entries */
static
SCIP_RETCODE ensureRstatMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   SCIPdebugMessage("calling ensureRstatMem()\n");

   assert(lpi != NULL);

   if( num > lpi->rstatsize )
   {
      int newsize;

      newsize = MAX( 2 * lpi->rstatsize, num );
      SCIP_ALLOC( BMSreallocMemoryArray( &lpi->rstat, newsize ) );
      lpi->rstatsize = newsize;
   }
   assert(num <= lpi->rstatsize);

   return SCIP_OKAY;
}

/*
 * LPi state methods
 */

/** returns the number of packets needed to store column packet information */
static
int colpacketNum(
   int                   ncols               /**< number of columns to store */
   )
{
   return (ncols + (int)COLS_PER_PACKET - 1) / (int)COLS_PER_PACKET;
}

/** returns the number of packets needed to store row packet information */
static
int rowpacketNum(
   int                   nrows               /**< number of rows to store */
   )
{
   return (nrows + (int)ROWS_PER_PACKET - 1) / (int)ROWS_PER_PACKET;
}

/** store row and column basis status in a packed LPi state object */
static
void lpistatePack(
   SCIP_LPISTATE*       lpistate,            /**< pointer to LPi state data */
   const int*           cstat,               /**< basis status of columns in unpacked format */
   const int*           rstat                /**< basis status of rows in unpacked format */
   )
{
   assert(lpistate != NULL);
   assert(lpistate->packcstat != NULL);
   assert(lpistate->packrstat != NULL);

   SCIPencodeDualBit(cstat, lpistate->packcstat, lpistate->ncols);
   SCIPencodeDualBit(rstat, lpistate->packrstat, lpistate->nrows);
}

/** unpacks row and column basis status from a packed LPi state object */
static
void lpistateUnpack(
   const SCIP_LPISTATE* lpistate,            /**< pointer to LPi state data */
   int*                 cstat,               /**< buffer for storing basis status of columns in unpacked format */
   int*                 rstat                /**< buffer for storing basis status of rows in unpacked format */
)
{
   assert(lpistate != NULL);
   assert(lpistate->packcstat != NULL);
   assert(lpistate->packrstat != NULL);

   SCIPdecodeDualBit(lpistate->packcstat, cstat, lpistate->ncols);
   SCIPdecodeDualBit(lpistate->packrstat, rstat, lpistate->nrows);
}

/** creates LPi state information object */
static
SCIP_RETCODE lpistateCreate(
   SCIP_LPISTATE**       lpistate,           /**< pointer to LPi state */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   ncols,              /**< number of columns to store */
   int                   nrows               /**< number of rows to store */
   )
{
   assert(lpistate != NULL);
   assert(blkmem != NULL);
   assert(ncols >= 0);
   assert(nrows >= 0);

   int nColPackets = colpacketNum(ncols);
   int nRowPackets = rowpacketNum(nrows);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, lpistate) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packcstat, nColPackets) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packrstat, nRowPackets) );

   return SCIP_OKAY;
}

/** frees LPi state information */
static
void lpistateFree(
   SCIP_LPISTATE**       lpistate,           /**< pointer to LPi state information (like basis information) */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(blkmem != NULL);
   assert(lpistate != NULL);
   assert(*lpistate != NULL);

   int nColPackets = colpacketNum((*lpistate)->ncols);
   int nRowPackets = rowpacketNum((*lpistate)->nrows);

   BMSfreeBlockMemoryArray( blkmem, &(*lpistate)->packcstat, nColPackets );
   BMSfreeBlockMemoryArray( blkmem, &(*lpistate)->packrstat, nRowPackets );
   BMSfreeBlockMemory( blkmem, lpistate);
}


/*
 * local methods
 */

/** marks the current LP to be unsolved */
static
void invalidateSolution(
   SCIP_LPI              *lpi                /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   lpi->solved = FALSE;
}

/** converts basis statuses */
static
HighsBasisStatus basestatToHighsBasisStatus(
   const int             &stat
   )
{
   switch( stat )
   {
   case SCIP_BASESTAT_LOWER:
      return HighsBasisStatus::kLower;
   case SCIP_BASESTAT_BASIC:
      return HighsBasisStatus::kBasic;
   case SCIP_BASESTAT_UPPER:
      return HighsBasisStatus::kUpper;
   case SCIP_BASESTAT_ZERO:
      return HighsBasisStatus::kZero;
   default:
      assert( false );
      SCIPerrorMessage("Unknown Basis Status returned. Please use supported HiGHS version!\n");
      return HighsBasisStatus::kZero;
   }
}

/** returns a string representation of the simplex strategy parameter */
static
std::string simplexStrategyToString(
   const int             &strategy
   )
{
   switch( strategy )
   {
   case 0:
      return "Choose";
   case 1:
      return "Dual (serial)";
   case 2:
      return "Dual (PAMI)";
   case 3:
      return "Dual (SIP)";
   case 4:
      return "Primal";
   default:
      return "Unknown";
   }
}

/** checks that matrix values are within range defined by HiGHS parameters */
static
SCIP_RETCODE checkMatrixValue(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             value               /**< value of coefficient */
   )
{
#ifndef NDEBUG
   SCIP_Real small_matrix_value;
   SCIP_Real large_matrix_value;

   HIGHS_CALL( lpi->highs->getOptionValue("small_matrix_value", small_matrix_value) );
   HIGHS_CALL( lpi->highs->getOptionValue("large_matrix_value", large_matrix_value) );

   assert(value == 0.0 || fabs(value) > small_matrix_value);
   assert(fabs(value) < large_matrix_value);
#endif

   return SCIP_OKAY;
}

/** calls HiGHS to solve the LP with given settings */
static
SCIP_RETCODE lpiSolve(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   std::string presolvestring;

   nsolvecalls++;
   int ck_ca_n = -99999;
   const bool check_lp = nsolvecalls == ck_ca_n;

   SCIPdebugMessage("HiGHS LP solve is called for the %d time\n", nsolvecalls);

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   if( lpi->fromscratch )
   {
      HIGHS_CALL( lpi->highs->clearSolver() );
   }

   lpi->highs->zeroAllClocks();

   /* the optimization result may be reliable even if HiGHS returns a warning status, e.g., HiGHS always returns with a
    * warning status if the iteration limit was hit
    */
   HIGHS_CALL_WITH_WARNING( lpi->highs->run() );

   HighsModelStatus model_status = lpi->highs->getModelStatus();
   switch( model_status )
   {
   /* solved or resource limit reached */
   case HighsModelStatus::kModelEmpty:
   case HighsModelStatus::kOptimal:
   case HighsModelStatus::kInfeasible:
   case HighsModelStatus::kUnboundedOrInfeasible:
   case HighsModelStatus::kUnbounded:
   case HighsModelStatus::kObjectiveBound:
   case HighsModelStatus::kTimeLimit:
   case HighsModelStatus::kIterationLimit:
#ifdef SCIP_DEBUG
      {
         int simplex_strategy = -1;
         HIGHS_CALL( lpi->highs->getOptionValue("simplex_strategy", simplex_strategy) );
         SCIPdebugMessage("HiGHS terminated with model status <%s> (%d) after simplex strategy <%s> (%d)\n",
            lpi->highs->modelStatusToString(model_status).c_str(), (int)model_status,
            simplexStrategyToString(simplex_strategy).c_str(), simplex_strategy);
      }
#endif
      break;
   /* errors or cases that should not occur in this LP interface */
   case HighsModelStatus::kNotset:
   case HighsModelStatus::kLoadError:
   case HighsModelStatus::kModelError:
   case HighsModelStatus::kPresolveError:
   case HighsModelStatus::kSolveError:
   case HighsModelStatus::kPostsolveError:
   case HighsModelStatus::kSolutionLimit:
   case HighsModelStatus::kObjectiveTarget:
   case HighsModelStatus::kUnknown:
   default:
      {
         int simplex_strategy = -1;
         HIGHS_CALL( lpi->highs->getOptionValue("simplex_strategy", simplex_strategy) );
         SCIPerrorMessage("HiGHS terminated with model status <%s> (%d) after simplex strategy <%s> (%d)\n",
            lpi->highs->modelStatusToString(model_status).c_str(), (int)model_status,
            simplexStrategyToString(simplex_strategy).c_str(), simplex_strategy);
      }
      return SCIP_LPERROR;
   }

   /* if basis factorization is unavailable, this may be due to presolving; then solve again without presolve */
   HIGHS_CALL( lpi->highs->getOptionValue("presolve", presolvestring) );
   assert(presolvestring == "on" || presolvestring == "off"); /* values used in SCIPlpiSetIntpar() */
   if( !lpi->highs->hasInvert() && presolvestring == "on" )
   {
      SCIP_RETCODE retcode;

      SCIPdebugMessage("No inverse: running HiGHS again without presolve . . .\n");
      HIGHS_CALL( lpi->highs->setOptionValue("presolve", "off") );
      retcode = lpiSolve(lpi);
      if( retcode != SCIP_OKAY )
      {
         HighsModelStatus model_status = lpi->highs->getModelStatus();
         SCIPerrorMessage("HiGHS terminated with model status <%s> (%d) after trying to recover inverse\n",
            lpi->highs->modelStatusToString(model_status).c_str(), (int)model_status);
      }
      HIGHS_CALL( lpi->highs->setOptionValue("presolve", "on") );
      SCIP_CALL( retcode );
   }

   if( check_lp )
   {
      int highs_iterations;
      HIGHS_CALL( lpi->highs->getInfoValue("simplex_iteration_count", highs_iterations) );
      SCIPdebugMessage("After call %d o solve() f=%15g; Iter = %d; Status = %s\n", nsolvecalls,
         lpi->highs->getObjectiveValue(), highs_iterations,
         lpi->highs->modelStatusToString(lpi->highs->getModelStatus()).c_str());
   }

   lpi->solved = TRUE;
   return SCIP_OKAY;
}


/*
 * LP Interface Methods
 */

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

static char highsname[30];
static char highsdesc[200];

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSolverName()\n");

   snprintf(highsname, 30, "HiGHS %d.%d.%d", HIGHS_VERSION_MAJOR, HIGHS_VERSION_MINOR, HIGHS_VERSION_PATCH);
   return highsname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSolverDesc()\n");

   snprintf(highsdesc, 200, "%s [%s] [GitHash: %s]",
      "Linear optimization suite written and engineered at the University of Edinburgh",
      HIGHS_COMPILATION_DATE, HIGHS_GITHASH);
   return highsdesc;
}

/** gets pointer for LP solver - use only with great care */
void *SCIPlpiGetSolverPointer(
      SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSolverPointer()\n");
   assert(lpi != NULL);
   return (void *) lpi->highs;
}

/** pass integrality information about variables to the solver */
SCIP_RETCODE SCIPlpiSetIntegralityInformation(
   SCIP_LPI              *lpi,               /**< pointer to an LP interface structure */
   int                   ncols,              /**< length of integrality array */
   int                   *intInfo            /**< integrality array (0: continuous, 1: integer) */
   )
{
   SCIPdebugMessage("calling SCIPlpiSetIntegralityInformation()\n");

   assert( lpi != NULL );
   assert( ncols >= 0 );
   assert( ncols == 0 || intInfo != NULL );

   SCIPerrorMessage("SCIPlpiSetIntegralityInformation() has not been implemented yet\n");

   return SCIP_LPERROR;
}

/** informs about availability of a primal simplex solving method */
SCIP_Bool SCIPlpiHasPrimalSolve(
   void
   )
{
   SCIPdebugMessage("calling SCIPlpiHasPrimalSolve()\n");
   return TRUE;
}

/** informs about availability of a dual simplex solving method */
SCIP_Bool SCIPlpiHasDualSolve(
   void
)
{
   SCIPdebugMessage("calling SCIPlpiHasDualSolve()\n");
   return TRUE;
}

/** informs about availability of a barrier solving method */
SCIP_Bool SCIPlpiHasBarrierSolve(
   void
   )
{
   SCIPdebugMessage("calling SCIPlpiHasBarrierSolve()\n");
   return FALSE;
}

/**@} */

/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** creates an LP problem object */
SCIP_RETCODE SCIPlpiCreate(
   SCIP_LPI              **lpi,              /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR      *messagehdlr,       /**< message handler to use for printing messages, or NULL */
   const char            *name,              /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   SCIPdebugMessage("calling SCIPlpiCreate()\n");

   SCIP_ALLOC( BMSallocMemory(lpi) );

   (*lpi)->highs = new HighsSCIP();
   HIGHS_CALL( (*lpi)->highs->clearModel() );

   /* initialize LPI data */
   (*lpi)->cstat = NULL;
   (*lpi)->rstat = NULL;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->nthreads = 1;
   (*lpi)->fromscratch = FALSE;
   (*lpi)->solved = FALSE;
   (*lpi)->pricing = SCIP_PRICING_LPIDEFAULT;
   (*lpi)->messagehdlr = messagehdlr;

   invalidateSolution(*lpi);

   /* set objective sense */
   SCIP_CALL( SCIPlpiChgObjsen(*lpi, objsen) );

   /* set output and debug level */
   HIGHS_CALL( (*lpi)->highs->setOptionValue("output_flag", false) );
#ifdef HIGHS_LOGDEVLEVEL
   HIGHS_CALL( (*lpi)->highs->setOptionValue("log_dev_level", HIGHS_LOGDEVLEVEL) );
#endif
#ifdef HIGHS_DEBUGLEVEL
   HIGHS_CALL( (*lpi)->highs->setOptionValue("highs_debug_level", HIGHS_DEBUGLEVEL) );
#endif

   /* set default scaling */
   SCIP_CALL( SCIPlpiSetIntpar(*lpi, SCIP_LPPAR_SCALING, 1) );

   /* use presolve by default; HiGHS runs without presolving whenever a basis is available */
   SCIP_CALL( SCIPlpiSetIntpar(*lpi, SCIP_LPPAR_PRESOLVING, 1) );
   HIGHS_CALL( (*lpi)->highs->setOptionValue("lp_presolve_requires_basis_postsolve", true) );

   /* set default pricing */
   SCIP_CALL( SCIPlpiSetIntpar(*lpi, SCIP_LPPAR_PRICING, (int)(*lpi)->pricing) );

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiFree()\n");

   assert(*lpi != NULL);
   assert((*lpi)->highs != NULL);

   /* free model and solver using destructor */
   (*lpi)->highs->~HighsSCIP();

   /* free basis arrays */
   BMSfreeMemoryArrayNull(&(*lpi)->cstat);
   BMSfreeMemoryArrayNull(&(*lpi)->rstat);

   /* free LPI memory */
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
   SCIPdebugMessage("calling SCIPlpiLoadColLP()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(beg != NULL);
   assert(ind != NULL);
   assert(val != NULL);

   assert(nrows >= 0);
   assert(ncols >= 0);

   assert(nnonz == 0 || ( nrows > 0 && ncols > 0));
#ifndef NDEBUG
   for( int j = 0; j < nnonz; ++j )
   {
      assert(0 <= ind[j] && ind[j] < nrows);
      assert(val[j] != 0.0);
      SCIP_CALL( checkMatrixValue(lpi, val[j]) );
   }
#endif

   int objectiveSenseInt = objsen == SCIP_OBJSEN_MAXIMIZE ? (int)ObjSense::kMaximize : (int)ObjSense::kMinimize;
   HIGHS_CALL( lpi->highs->passModel(ncols, nrows, nnonz, 1, objectiveSenseInt, 0, obj, lb, ub, lhs, rhs, beg, ind, val) );

   assert((objsen == SCIP_OBJSEN_MAXIMIZE && lpi->highs->getLp().sense_ == ObjSense::kMaximize)
      || (objsen == SCIP_OBJSEN_MINIMIZE && lpi->highs->getLp().sense_ == ObjSense::kMinimize));

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
   const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiAddCols()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);
   assert(nnonz >= 0);
   assert(ncols >= 0);
   assert(ncols <= 0 || obj != NULL);
   assert(ncols <= 0 || lb != NULL);
   assert(ncols <= 0 || ub != NULL);

   invalidateSolution(lpi);

#ifndef NDEBUG
   if( nnonz > 0 )
   {
      /* perform check that no new rows are added - this is likely to be a mistake
       */
      int nrows = lpi->highs->getLp().num_row_;
      for( int j = 0; j < nnonz; ++j )
      {
         assert(0 <= ind[j] && ind[j] < nrows);
         assert(val[j] != 0.0);
         SCIP_CALL( checkMatrixValue(lpi, val[j]) );
      }
   }

   /* HiGHS returns with a warning if values are within the zero tolerance, but seems to continue safely simply ignoring
    * them; in debug mode we stop, in optimized mode we accept this behavior */
   HIGHS_CALL( lpi->highs->addCols(ncols, obj, lb, ub, nnonz, beg, ind, val) );
#else
   HIGHS_CALL_WITH_WARNING( lpi->highs->addCols(ncols, obj, lb, ub, nnonz, beg, ind, val) );
#endif

   return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelCols()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(lpi->highs->getLp().num_col_ >= 0);

   invalidateSolution(lpi);
   HIGHS_CALL( lpi->highs->deleteCols(firstcol, lastcol) );

   assert(lpi->highs->getLp().num_col_ >= 0);

   return SCIP_OKAY;
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              * input:  1 if column should be deleted, 0 if not
                                              * output: new position of column, -1 if column was deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelColset()\n");

   assert(lpi != NULL);
   assert(dstat != NULL);
   assert(lpi->highs->getLp().num_col_ >= 0);

   invalidateSolution(lpi);

   HIGHS_CALL( lpi->highs->deleteCols(dstat) );

   assert(lpi->highs->getLp().num_col_ >= 0);
   return SCIP_OKAY;
}

/** adds rows to the LP */
SCIP_RETCODE SCIPlpiAddRows(
      SCIP_LPI*             lpi,                /**< LP interface structure */
      int                   nrows,              /**< number of rows to be added */
      const SCIP_Real*      lhs,                /**< left hand sides of new rows */
      const SCIP_Real*      rhs,                /**< right hand sides of new rows */
      char**                rownames,           /**< row names, or NULL */
      int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
      const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
      const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
      const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiAddRows()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(nrows >= 0);
   assert(nrows <= 0 || lhs != NULL);
   assert(nrows <= 0 || rhs != NULL);
   assert(nnonz >= 0);
   assert(nnonz <= 0 || beg != NULL);
   assert(nnonz <= 0 || ind != NULL);
   assert(nnonz <= 0 || val != NULL);

   invalidateSolution(lpi);

#ifndef NDEBUG
   if( nnonz > 0 )
   {
      /* Perform check that no new columns are added - this is likely to be a mistake - and that the values are nonzero*/
      int ncols = lpi->highs->getLp().num_col_;
      for( int j = 0; j < nnonz; ++j )
      {
         assert(0 <= ind[j] && ind[j] < ncols);
         assert(val[j] != 0.0);
         SCIP_CALL( checkMatrixValue(lpi, val[j]) );
      }
   }

   /* HiGHS returns with a warning if values are within the zero tolerance, but seems to continue safely simply ignoring
    * them; in debug mode we stop, in optimized mode we accept this behavior */
   HIGHS_CALL( lpi->highs->addRows(nrows, lhs, rhs, nnonz, beg, ind, val) );
#else
   HIGHS_CALL_WITH_WARNING( lpi->highs->addRows(nrows, lhs, rhs, nnonz, beg, ind, val) );
#endif

   return SCIP_OKAY;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelRows()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(lpi->highs->getLp().num_row_ >= 0);
   assert(0 <= firstrow && firstrow <= lastrow );

   invalidateSolution(lpi);
   HIGHS_CALL( lpi->highs->deleteRows(firstrow, lastrow) );

   assert(lpi->highs->getLp().num_row_ >= 0);
   return SCIP_OKAY;
}

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelRowset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelRowset()\n");

   assert(lpi != NULL);
   assert(dstat != NULL);
   assert(lpi->highs != NULL);
   assert(lpi->highs->getLp().num_row_ >= 0);

   invalidateSolution(lpi);

   HIGHS_CALL( lpi->highs->deleteRows(dstat) );

   assert(lpi->highs->getLp().num_row_ >= 0);

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiClear()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(lpi->highs->getLp().num_row_ >= 0);
   assert(lpi->highs->getLp().num_col_ >= 0);

   invalidateSolution(lpi);

   HIGHS_CALL( lpi->highs->clearModel() );
   return SCIP_OKAY;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiChgBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices or NULL if ncols is zero */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds or NULL if ncols is zero */
   const SCIP_Real*      ub                  /**< values for the new upper bounds or NULL if ncols is zero */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgBounds()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(ind != NULL);
   assert(lb != NULL);
   assert(ub != NULL);

   invalidateSolution(lpi);

   int i;

   /* Check validity of data */
   for( i = 0; i < ncols; ++i )
   {
      assert(0 <= ind[i] && ind[i] < lpi->highs->getLp().num_col_);

      if( SCIPlpiIsInfinity(lpi, lb[i]) )
      {
         SCIPerrorMessage( "LP Error: fixing lower bound for variable %d to infinity\n", ind[i]);
         return SCIP_LPERROR;
      }
      if( SCIPlpiIsInfinity(lpi, -ub[i]) )
      {
         SCIPerrorMessage( "LP Error: fixing upper bound for variable %d to -infinity\n", ind[i]);
         return SCIP_LPERROR;
      }
   }

   HIGHS_CALL( lpi->highs->changeColsBounds(ncols, ind, lb, ub) );

   return SCIP_OKAY;
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiChgSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      lhs,                /**< new values for left hand sides */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgSides()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(ind != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);

   int i;

   invalidateSolution(lpi);

   for( i = 0; i < nrows; ++i )
      assert(0 <= ind[i] && ind[i] < lpi->highs->getLp().num_row_);

   HIGHS_CALL( lpi->highs->changeRowsBounds(nrows, ind, lhs, rhs) );

   return SCIP_OKAY;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiChgCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgCoef()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   invalidateSolution(lpi);

   SCIP_CALL( checkMatrixValue(lpi, newval) );
   HIGHS_CALL( lpi->highs->changeCoeff(row, col, newval) );

   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgObjsen()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   invalidateSolution(lpi);

   HIGHS_CALL( lpi->highs->changeObjectiveSense(objsen == SCIP_OBJSEN_MINIMIZE ? ObjSense::kMinimize : ObjSense::kMaximize) );

   return SCIP_OKAY;
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiChgObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   const int*            ind,                /**< column indices to change objective value for */
   const SCIP_Real*      obj                 /**< new objective values for columns */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgObj()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(ind != NULL);
   assert(obj != NULL);

   invalidateSolution(lpi);

   HIGHS_CALL( lpi->highs->changeColsCost(ncols, ind, obj) );

   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiScaleRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIPdebugMessage("calling SCIPlpiScaleRow()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   invalidateSolution(lpi);

   HIGHS_CALL( lpi->highs->scaleRow(row, scaleval) );

   return SCIP_OKAY;
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPlpiScaleCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIPdebugMessage("calling SCIPlpiScaleCol()\n");

   assert(lpi != NULL);
   assert(scaleval != 0.0);

   invalidateSolution(lpi);

   HIGHS_CALL( lpi->highs->scaleCol(col, scaleval) );

   return SCIP_OKAY;
}

/**@} */

/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiGetNRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNRows()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(nrows != NULL);
   *nrows = lpi->highs->getNumRow();
   assert(*nrows >= 0);

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNCols()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(ncols != NULL);
   *ncols = lpi->highs->getNumCol();
   assert(*ncols >= 0);

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNNonz()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(nnonz != NULL);
   *nnonz = lpi->highs->getNumNz();
   assert(*nnonz >= 0);

   return SCIP_OKAY;
}

/** gets columns from LP problem object; the arrays have to be large enough to
 * store all values Either both, lb and ub, have to be NULL, or both have to be
 * non-NULL, either nnonz, beg, ind, and val have to be NULL, or all of them
 * have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetCols(
      SCIP_LPI*             lpi,                /**< LP interface structure */
      int                   firstcol,           /**< first column to get from LP */
      int                   lastcol,            /**< last column to get from LP */
      SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
      SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
      int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
      int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
      int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
      SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
      )
{
   SCIPdebugMessage("calling SCIPlpiGetCols()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   int num_col;
   HIGHS_CALL( lpi->highs->getCols(firstcol, lastcol, num_col, NULL, lb, ub, *nnonz, beg, ind, val) );
   return SCIP_OKAY;
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
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
   SCIPdebugMessage("calling SCIPlpiGetRows()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   int num_row;
   HIGHS_CALL( lpi->highs->getRows(firstrow, lastrow, num_row, lhs, rhs, *nnonz, beg, ind, val) );
   return SCIP_OKAY;
}

/** gets column names */
SCIP_RETCODE SCIPlpiGetColNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get name from LP */
   int                   lastcol,            /**< last column to get name from LP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for col names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetColNames()\n");

   assert(lpi != NULL);

   SCIPerrorMessage("SCIPlpiGetColNames() has not been implemented yet\n");

   return SCIP_PLUGINNOTFOUND;
}

/** gets row names */
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
   SCIPdebugMessage("calling SCIPlpiGetRowNames()\n");

   assert(lpi != NULL);

   SCIPerrorMessage("SCIPlpiGetRowNames() has not been implemented yet\n");

   return SCIP_PLUGINNOTFOUND;
}

/** gets the objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetObjsen()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   *objsen = SCIP_OBJSEN_MINIMIZE;
   if( lpi->highs->getLp().sense_ == ObjSense::kMaximize )
      *objsen = SCIP_OBJSEN_MAXIMIZE;

   return SCIP_OKAY;
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiGetObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetObj()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->highs->getLp().num_col_);
   assert(vals != NULL);

   for( int i = firstcol; i < lastcol + 1; ++i )
      vals[i - firstcol] = lpi->highs->getLp().col_cost_[i];

   return SCIP_OKAY;
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiGetBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective value for */
   int                   lastcol,            /**< last column to get objective value for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBounds()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->highs->getLp().num_col_);

   for( int i = firstcol; i < lastcol + 1; ++i )
   {
      if( lbs != NULL )
         lbs[i - firstcol] = lpi->highs->getLp().col_lower_[i];
      if( ubs != NULL )
         ubs[i - firstcol] = lpi->highs->getLp().col_upper_[i];
   }

   return SCIP_OKAY;
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiGetSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSides()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->highs->getLp().num_row_);

   for( int i = firstrow; i < lastrow + 1; ++i )
   {
      if( lhss != NULL )
         lhss[i - firstrow] = lpi->highs->getLp().row_lower_[i];
      if( rhss != NULL )
         rhss[i - firstrow] = lpi->highs->getLp().row_upper_[i];
   }

   return SCIP_OKAY;
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiGetCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
)
{
   SCIPdebugMessage("calling SCIPlpiGetCoef()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(0 <= col && col < lpi->highs->getNumCol());
   assert(0 <= row && row < lpi->highs->getNumCol());
   assert(val != NULL);

   HIGHS_CALL( lpi->highs->getCoeff(row, col, *val) );
   return SCIP_OKAY;
}

/**@} */

/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
      SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiSolvePrimal()\n");

   assert(lpi != NULL);

   /* HiGHS' primal simplex seems to still have performance issues, so we call the dual simplex instead. */
#ifdef SCIP_WITH_HIGHSPRIMAL
   HIGHS_CALL( lpi->highs->setOptionValue("parallel", "off") );
   HIGHS_CALL( lpi->highs->setOptionValue("threads", 1) );
   HIGHS_CALL( lpi->highs->setOptionValue("simplex_strategy", 4) );
   SCIP_CALL( lpiSolve(lpi) );
#else
   SCIP_CALL( SCIPlpiSolveDual(lpi) );
#endif

   return SCIP_OKAY;
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiSolveDual()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   /* HiGHS still seems to get stuck sometimes in parallel mode, so we ignore nthreads for now. */
#ifdef SCIP_WITH_HIGHSPARALLEL
   if( lpi->nthreads == 0 || lpi->nthreads > 1 )
   {
      SCIPdebugMessage("Running HiGHS dual simplex in parallel with lpi->nthreads=%d\n", lpi->nthreads);
      HIGHS_CALL( lpi->highs->setOptionValue("parallel", "on") );
      HIGHS_CALL( lpi->highs->setOptionValue("threads", lpi->nthreads) ); /* note that also in HiGHS, 0 is the automatic setting */
      HIGHS_CALL( lpi->highs->setOptionValue("simplex_strategy", 2) ); /* PAMI */
   }
   else
#endif
   {
      SCIPdebugMessage("Running HiGHS dual simplex in serial with lpi->nthreads=%d\n", lpi->nthreads);
      HIGHS_CALL( lpi->highs->setOptionValue("parallel", "off") );
      HIGHS_CALL( lpi->highs->setOptionValue("threads", 1) );
      HIGHS_CALL( lpi->highs->setOptionValue("simplex_strategy", 1) );
   }

   SCIP_CALL( lpiSolve(lpi) );

   return SCIP_OKAY;
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{
   SCIPdebugMessage("calling SCIPlpiSolveBarrier()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   SCIPdebugMessage("HiGHS does not support Barrier - switching to dual simplex\n");
   return SCIPlpiSolveDual(lpi);
}

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiStartStrongbranch()\n");

   assert(lpi != NULL);

   /* no work necessary for current dummy implementation */
   return SCIP_OKAY;
}

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiEndStrongbranch()\n");

   assert(lpi != NULL);

   /* no work necessary for current dummy implementation */
   return SCIP_OKAY;
}

/** performs strong branching iterations on one @b fractional candidate */
SCIP_RETCODE SCIPlpiStrongbranchFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< fractional current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiStrongbranchFrac()\n");

   assert(lpi != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   /* This is a dummy implementation to satisfy the test suite. It does not perform actual strong branching. */
   SCIP_Real dualbound = (lpi->highs->getLp().sense_ == ObjSense::kMinimize
      ? -SCIPlpiInfinity(lpi) : SCIPlpiInfinity(lpi));

   if( SCIPlpiIsOptimal(lpi) )
   {
      SCIP_CALL( SCIPlpiGetObjval(lpi, &dualbound) );
   }

   *down = *up = dualbound;
   *downvalid = TRUE;
   *upvalid = TRUE;

   if( iter != NULL )
      *iter = -1;

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
   SCIPdebugMessage("calling SCIPlpiStrongbranchesFrac()\n");

   assert(lpi != NULL);
   assert(cols != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   /* This is a dummy implementation to satisfy the test suite. It does not perform actual strong branching. */
   SCIP_Real dualbound = (lpi->highs->getLp().sense_ == ObjSense::kMinimize
      ? -SCIPlpiInfinity(lpi) : SCIPlpiInfinity(lpi));

   if( SCIPlpiIsOptimal(lpi) )
   {
      SCIP_CALL( SCIPlpiGetObjval(lpi, &dualbound) );
   }

   for( int j = 0; j < ncols; ++j )
   {
      down[j] = up[j] = dualbound;
      downvalid[j] = upvalid[j] = TRUE;
   }

   if( iter != NULL )
      *iter = -1;

   return SCIP_OKAY;
}

/** performs strong branching iterations on one candidate with @b integral value */
SCIP_RETCODE SCIPlpiStrongbranchInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current integral primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiStrongbranchInt()\n");

   assert(lpi != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   /* the dummy implementation works independently of primal values. */
   SCIP_CALL( SCIPlpiStrongbranchFrac(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter) );
   return SCIP_OKAY;
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
   SCIPdebugMessage("calling SCIPlpiStrongbranchesInt()\n");

   assert(lpi != NULL);
   assert(cols != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   /* the dummy implementation works independently of primal values */
   SCIP_CALL( SCIPlpiStrongbranchesFrac(lpi, cols, ncols, psols, itlim, down, up, downvalid, upvalid, iter) );
   return SCIP_OKAY;
}

/**@} */

/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the
 * LP */
SCIP_Bool SCIPlpiWasSolved(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiWasSolved()\n");

   assert(lpi != NULL);

   return lpi->solved;
}

/** gets information about primal and dual feasibility of the current LP solution
 *
 *  The feasibility information is with respect to the last solving call and it is only relevant if SCIPlpiWasSolved()
 *  returns true. If the LP is changed, this information might be invalidated.
 *
 *  Note that @a primalfeasible and @a dualfeasible should only return true if the solver has proved the respective LP to
 *  be feasible. Thus, the return values should be equal to the values of SCIPlpiIsPrimalFeasible() and
 *  SCIPlpiIsDualFeasible(), respectively. Note that if feasibility cannot be proved, they should return false (even if
 *  the problem might actually be feasible).
 */
SCIP_RETCODE SCIPlpiGetSolFeasibility(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< pointer to store primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< pointer to store dual feasibility status */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSolFeasibility()\n");

   assert(lpi != NULL);
   assert(primalfeasible != NULL);
   assert(dualfeasible != NULL);

   *primalfeasible = SCIPlpiIsPrimalFeasible(lpi);
   *dualfeasible = SCIPlpiIsDualFeasible(lpi);

   return SCIP_OKAY;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExistsPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExistsPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   HighsModelStatus model_status = lpi->highs->getModelStatus();

   return model_status == HighsModelStatus::kUnbounded || model_status == HighsModelStatus::kUnboundedOrInfeasible;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   if( !SCIPlpiIsPrimalUnbounded(lpi) )
      return FALSE;

   /* HiGHS method does not work in this case, but we can easily construct an unbounded primal ray */
   if( lpi->highs->getNumRow() == 0 )
      return TRUE;

   bool has_primal_ray = false;
   HIGHS_CALL( lpi->highs->getPrimalRay(has_primal_ray, NULL) );
   return has_primal_ray;
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   return lpi->highs->getModelStatus() == HighsModelStatus::kUnbounded;
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   HighsModelStatus model_status = lpi->highs->getModelStatus();

   /* not sure how to query HiGHS in this case, but we can easily decide */
   if( model_status == HighsModelStatus::kModelEmpty )
   {
      int numrow = lpi->highs->getNumRow();

      assert(lpi->highs->getNumCol() == 0);

      for( int i = 0; i < numrow; i++ )
      {
         if( lpi->highs->getLp().row_lower_[i] > 0.0 || lpi->highs->getLp().row_upper_[i] < 0.0 )
            return TRUE;
      }
      return FALSE;
   }

   /* otherwise we rely on the model status */
   const bool primal_infeasible =
      model_status == HighsModelStatus::kInfeasible ||
      model_status == HighsModelStatus::kUnboundedOrInfeasible;
   return primal_infeasible;
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalFeasible()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   HighsModelStatus model_status = lpi->highs->getModelStatus();

   /* not sure how to query HiGHS in this case, but we can easily decide */
   if( model_status == HighsModelStatus::kModelEmpty )
   {
      int numrow = lpi->highs->getNumRow();

      assert(lpi->highs->getNumCol() == 0);

      for( int i = 0; i < numrow; i++ )
      {
         if( lpi->highs->getLp().row_lower_[i] > 0.0 || lpi->highs->getLp().row_upper_[i] < 0.0 )
            return FALSE;
      }
      return TRUE;
   }

   /* otherwise we rely on the model status */
   const bool primal_feasible =
      model_status == HighsModelStatus::kOptimal ||
      model_status == HighsModelStatus::kUnbounded;
   return primal_feasible;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExistsDualRay()\n");

   assert(lpi != NULL);

   return !SCIPlpiIsPrimalFeasible(lpi);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasDualRay()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   HighsModelStatus model_status = lpi->highs->getModelStatus();

   /* HiGHS does not implement this case, but we can easily decide */
   if( model_status == HighsModelStatus::kModelEmpty )
      return !SCIPlpiIsPrimalFeasible(lpi);

   /* otherwise we rely on the model status */
   bool has_dual_ray = false;
   HIGHS_CALL( lpi->highs->getDualRay(has_dual_ray, NULL) );
   return has_dual_ray;
}

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   return SCIPlpiIsDualFeasible(lpi) && !SCIPlpiIsPrimalFeasible(lpi);
}

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   HighsModelStatus model_status = lpi->highs->getModelStatus();
   const bool dual_infeasible =
      model_status == HighsModelStatus::kUnbounded ||
      model_status == HighsModelStatus::kUnboundedOrInfeasible;
   return dual_infeasible;
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualFeasible()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   HighsModelStatus model_status = lpi->highs->getModelStatus();

   if( model_status == HighsModelStatus::kOptimal || model_status == HighsModelStatus::kModelEmpty )
      return TRUE;
   else if( model_status == HighsModelStatus::kUnbounded || model_status == HighsModelStatus::kUnboundedOrInfeasible )
      return FALSE;

   int num_dual_infeasibilities = 1;
   HighsStatus status = lpi->highs->getInfoValue("num_dual_infeasibilities", num_dual_infeasibilities);
   bool has_dual_feasible_sol = (status == HighsStatus::kOk) && (num_dual_infeasibilities == 0);
   return has_dual_feasible_sol;
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsOptimal()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   HighsModelStatus model_status = lpi->highs->getModelStatus();

   if( model_status == HighsModelStatus::kModelEmpty )
      return SCIPlpiIsPrimalFeasible(lpi);
   else
   {
      assert(lpi->highs->getModelStatus() == HighsModelStatus::kOptimal || (!SCIPlpiIsPrimalFeasible(lpi) || !SCIPlpiIsDualFeasible(lpi)));
      assert(lpi->highs->getModelStatus() != HighsModelStatus::kOptimal || (SCIPlpiIsPrimalFeasible(lpi) && SCIPlpiIsDualFeasible(lpi)));
      return lpi->highs->getModelStatus() == HighsModelStatus::kOptimal;
   }
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsStable()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   /* if an objective limit is set and HiGHS claims that it is exceeded, we should check that this is indeed the case;
    * if not this points at numerical instability; note that this aligns with an assert in lp.c */
   if( SCIPlpiIsObjlimExc(lpi) )
   {
      SCIP_Real objlimit;
      SCIP_Real objvalue;

      HIGHS_CALL( lpi->highs->getOptionValue("objective_bound", objlimit) );
      HIGHS_CALL( lpi->highs->getInfoValue("objective_function_value", objvalue) );

      if( lpi->highs->getLp().sense_ == ObjSense::kMaximize )
      {
         objlimit *= -1.0;
         objvalue *= -1.0;
      }
      if( !SCIPlpiIsInfinity(lpi, objlimit) && HIGHS_relDiff(objvalue, objlimit) < -1e-9 )
         return FALSE;
   }

   return TRUE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsObjlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   return lpi->highs->getModelStatus() == HighsModelStatus::kObjectiveBound;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsIterlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   return lpi->highs->getModelStatus() == HighsModelStatus::kIterationLimit;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   return lpi->highs->getModelStatus() == HighsModelStatus::kTimeLimit;
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
      SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetInternalStatus()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   return (int) lpi->highs->getModelStatus();
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   SCIPdebugMessage("calling SCIPlpiIgnoreInstability()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   assert(success != NULL);

   *success = TRUE;
   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetObjval()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(objval != NULL);

   HIGHS_CALL( lpi->highs->getInfoValue("objective_function_value", *objval) );
   assert(lpi->highs->getModelStatus() != HighsModelStatus::kModelEmpty || *objval == 0.0);

   return SCIP_OKAY;
}

/** gets primal and dual solution vectors for feasible LPs
 *
 *  Before calling this function, the caller must ensure that the LP has been solved to optimality, i.e., that
 *  SCIPlpiIsOptimal() returns true.
 */
SCIP_RETCODE SCIPlpiGetSol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_Real*            activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSol()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   int ncols;
   int nrows;
   int i;

   if( objval != NULL )
   {
      HIGHS_CALL( lpi->highs->getInfoValue("objective_function_value", *objval) );
      assert(lpi->highs->getModelStatus() != HighsModelStatus::kModelEmpty || *objval == 0.0);
   }

   const std::vector<double> &colValue = lpi->highs->getSolution().col_value;
   const std::vector<double> &colDual = lpi->highs->getSolution().col_dual;
   const std::vector<double> &rowValue = lpi->highs->getSolution().row_value;
   const std::vector<double> &rowDual = lpi->highs->getSolution().row_dual;

   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   if( colValue.size() != (size_t) ncols || colDual.size() != (size_t) ncols
      || rowValue.size() != (size_t) nrows || rowDual.size() != (size_t) nrows )
   {
      SCIPmessagePrintWarning( lpi->messagehdlr, "In HiGHS the size of the columns values %d does not fit the number of columns %d\n", (int) colValue.size(), ncols);
      SCIPmessagePrintWarning( lpi->messagehdlr, "In HiGHS the size of the dual values %d does not fit the number of columns %d\n", (int) colDual.size(), ncols);
      SCIPmessagePrintWarning( lpi->messagehdlr, "In HiGHS the size of the rows values %d does not fit the number of rows %d\\n\"", (int) rowValue.size(), nrows);
      SCIPmessagePrintWarning( lpi->messagehdlr, "In HiGHS the size of the dual row values %d does not fit the number of rows %d\\n\"", (int) rowDual.size(), nrows);
      assert((int) rowValue.size() == nrows);
      SCIPmessagePrintWarning( lpi->messagehdlr, "HiGHS returned solution vector of inconsistent dimension\n" );
      return SCIP_LPERROR;
   }

   if( primsol != NULL )
      for( i = 0; i < ncols; i++ )
         primsol[i] = colValue[i];
   if( dualsol != NULL )
      for( i = 0; i < nrows; i++ )
         dualsol[i] = rowDual[i];
   if( activity != NULL )
      for( i = 0; i < nrows; i++ )
         activity[i] = rowValue[i];
   if( redcost != NULL )
      for( i = 0; i < ncols; i++ )
         redcost[i] = colDual[i];

   return SCIP_OKAY;
}


/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetPrimalRay()\n");

   bool success = false;

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(SCIPlpiHasPrimalRay(lpi));

   /* HiGHS does not implement this case, but we can easily construct an unbounded primal ray */
   if( lpi->highs->getNumRow() == 0 )
   {
      int numcol = lpi->highs->getNumCol();

      for( int i = 0; i < numcol; i++ )
      {
         SCIP_Real minobj = lpi->highs->getLp().col_cost_[i];
         if( lpi->highs->getLp().sense_ == ObjSense::kMaximize )
            minobj = -minobj;

         if( SCIPlpiIsInfinity(lpi, -lpi->highs->getLp().col_lower_[i]) && minobj > 0.0 )
         {
            ray[i] = -1.0;
            success = true;
         }
         else if( SCIPlpiIsInfinity(lpi, lpi->highs->getLp().col_upper_[i]) && minobj < 0.0 )
         {
            ray[i] = 1.0;
            success = true;
         }
         else
            ray[i] = 0.0;
      }
   }
   else
   {
      HIGHS_CALL( lpi->highs->getPrimalRay(success, ray) );
   }

   return success ? SCIP_OKAY : SCIP_LPERROR;
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual farkas row multipliers */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetDualfarkas()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(dualfarkas != NULL);

   HighsModelStatus model_status = lpi->highs->getModelStatus();

   /* HiGHS does not implement this case, but we can easily construct an unbounded dual ray */
   if( model_status == HighsModelStatus::kModelEmpty )
   {
      SCIP_Real dualdir = lpi->highs->getLp().sense_ == ObjSense::kMinimize ? 1.0 : -1.0;
      int numrow = lpi->highs->getNumRow();

      assert(lpi->highs->getNumCol() == 0);

      for( int i = 0; i < numrow; i++ )
      {
         if( lpi->highs->getLp().row_lower_[i] > 0.0 )
            dualfarkas[i] = dualdir;
         else if( lpi->highs->getLp().row_upper_[i] < 0.0 )
            dualfarkas[i] = -dualdir;
         else
            dualfarkas[i] = 0.0;
      }

      return SCIP_OKAY;
   }

   bool has_dual_ray = false;
   HIGHS_CALL( lpi->highs->getDualRay(has_dual_ray, dualfarkas) );

   return has_dual_ray ? SCIP_OKAY : SCIP_LPERROR;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetIterations()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(iterations != NULL);

   *iterations = 0;
   /* this may return with a warning if the last solve failed */
   HIGHS_CALL_WITH_WARNING( lpi->highs->getInfoValue("simplex_iteration_count", *iterations) );
   assert(*iterations >= 0);
   return SCIP_OKAY;
}

/** gets information about the quality of an LP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The LPI should return SCIP_INVALID for @p quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPlpiGetRealSolQuality(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPSOLQUALITY     qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetRealSolQuality()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(quality != NULL);

   *quality = SCIP_INVALID;

   return SCIP_OKAY;
}

/**@} */

/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiGetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBase()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   if( cstat != NULL )
   {
      for( int i = 0; i < lpi->highs->getLp().num_col_; ++i )
         cstat[i] = (int) lpi->highs->getBasis().col_status[i];
   }
   if( rstat != NULL )
   {
      for( int i = 0; i < lpi->highs->getLp().num_row_; ++i )
         rstat[i] = (int) lpi->highs->getBasis().row_status[i];
   }

   return SCIP_OKAY;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiSetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const int*            cstat,              /**< array with column basis status */
   const int*            rstat               /**< array with row basis status */
   )
{
   SCIPdebugMessage("calling SCIPlpiSetBase()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   HighsBasis local_highs_basis;

   local_highs_basis.col_status.resize(lpi->highs->getLp().num_col_);
   local_highs_basis.row_status.resize(lpi->highs->getLp().num_row_);

   if( cstat != NULL )
   {
      for( int i = 0; i < lpi->highs->getLp().num_col_; ++i )
         local_highs_basis.col_status[i] = basestatToHighsBasisStatus(cstat[i]);
   }
   if( rstat != NULL )
   {
      for( int i = 0; i < lpi->highs->getLp().num_row_; ++i )
         local_highs_basis.row_status[i] = basestatToHighsBasisStatus(rstat[i]);
   }
   HIGHS_CALL( lpi->highs->setBasis(local_highs_basis) );

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< pointer to store basis indices ready to keep number of rows entries */
   )
 {
    SCIPdebugMessage("calling SCIPlpiGetBasisInd()\n");

    assert(lpi != NULL);
    assert(lpi->highs != NULL);
    assert(bind != NULL);

    if( !lpi->highs->getBasis().valid )
    {
       SCIPdebugMessage( "HiGHS Basis is not valid in function call SCIPlpiGetBasisInd()\n" );
       return SCIP_ERROR;
    }
    HIGHS_CALL( lpi->highs->getBasicVariables(bind) );

    return SCIP_OKAY;
}

/** get row of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
)
{
   SCIPdebugMessage("calling SCIPlpiGetBInvRow()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   if( lpi->highs->getBasisInverseRow(r, coef, ninds, inds) != HighsStatus::kOk )
   {
      SCIP_CALL( SCIPlpiSolveDual(lpi) );
   }

   HIGHS_CALL( lpi->highs->getBasisInverseRow(r, coef, ninds, inds) );
   return SCIP_OKAY;
}

/** get column of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
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
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvCol()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   if( lpi->highs->getBasisInverseCol(c, coef, ninds, inds) != HighsStatus::kOk )
   {
      SCIP_CALL( SCIPlpiSolveDual(lpi) );
   }

   HIGHS_CALL( lpi->highs->getBasisInverseCol(c, coef, ninds, inds) );

   return SCIP_OKAY;
}

/** get row of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            coef,               /**< vector to return coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvARow()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   if( lpi->highs->getReducedRow(r, coef, ninds, inds, binvrow) != HighsStatus::kOk )
   {
      SCIP_CALL( SCIPlpiSolveDual(lpi) );
   }

   HIGHS_CALL( lpi->highs->getReducedRow(r, coef, ninds, inds, binvrow) );

   return SCIP_OKAY;
}

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This
 * means that if, internally, the LP solver uses a -1 coefficient, then rows
 * associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef,               /**< vector to return coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvACol()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   if( lpi->highs->getReducedColumn(c, coef, ninds, inds) != HighsStatus::kOk )
   {
      SCIP_CALL( SCIPlpiSolveDual(lpi) );
   }

   HIGHS_CALL( lpi->highs->getReducedColumn(c, coef, ninds, inds) );
   return SCIP_OKAY;
}

/**@} */

/*
 * LP State Methods
 */

/**@name LP State Methods */
/**@{ */

/** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiGetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetState()\n");

   assert(blkmem != NULL);

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(lpistate != NULL);

   int ncols;
   int nrows;

   ncols = lpi->highs->getLp().num_col_;
   nrows = lpi->highs->getLp().num_row_;
   assert(ncols >= 0);
   assert(nrows >= 0);

   /* allocate lpistate data */
   SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows) );

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, ncols) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   /* get unpacked basis information */
   SCIP_CALL( SCIPlpiGetBase(lpi, lpi->cstat, lpi->rstat) );

   /* pack LPi state data */
   (*lpistate)->ncols = ncols;
   (*lpistate)->nrows = nrows;
   lpistatePack(*lpistate, lpi->cstat, lpi->rstat);

   return SCIP_OKAY;
}

/** loads LPi state (like basis information) into solver; note that the LP might
 * have been extended with additional columns and rows since the state was
 * stored with SCIPlpiGetState()
 */
SCIP_RETCODE SCIPlpiSetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const SCIP_LPISTATE*  lpistate            /**< LPi state information (like basis information), or NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiSetState()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(lpistate != NULL);

   int lpncols;
   int lpnrows;
   int i;

   lpncols = lpi->highs->getLp().num_col_;
   lpnrows = lpi->highs->getLp().num_row_;
   assert(lpistate->ncols <= lpncols);
   assert(lpistate->nrows <= lpnrows);

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, lpncols) );
   SCIP_CALL( ensureRstatMem(lpi, lpnrows) );

   /* unpack LPi state data */
   lpistateUnpack(lpistate, lpi->cstat, lpi->rstat);

   /* extend the basis to the current LP beyond the previously existing columns */
   for( i = lpistate->ncols; i < lpncols; ++i )
   {
      if( !SCIPlpiIsInfinity(lpi, -lpi->highs->getLp().col_lower_[i]) )
         /* use finite lower bound */
         lpi->cstat[i] = SCIP_BASESTAT_LOWER;
      else if( !SCIPlpiIsInfinity(lpi, lpi->highs->getLp().col_upper_[i]) )
         /* use finite upper bound */
         lpi->cstat[i] = SCIP_BASESTAT_UPPER;
      else
         /* variable is free */
         lpi->cstat[i] = SCIP_BASESTAT_ZERO;
   }
   for( i = lpistate->nrows; i < lpnrows; ++i )
      lpi->rstat[i] = SCIP_BASESTAT_BASIC;

   /* load basis information */
   SCIP_CALL( SCIPlpiSetBase(lpi, lpi->cstat, lpi->rstat) );

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
)
{
   SCIPdebugMessage("calling SCIPlpiClearState()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   HIGHS_CALL( lpi->highs->clearSolver() );
   return SCIP_OKAY;
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiFreeState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
)
{
   SCIPdebugMessage("calling SCIPlpiFreeState()\n");

   assert(lpi != NULL);
   assert(lpistate != NULL);

   if( *lpistate != NULL )
      lpistateFree(lpistate, blkmem);

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information), or NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasStateBasis()\n");
   assert(lpi != NULL);
   return TRUE;
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,               /**< LP interface structure */
   const char*           fname              /**< file name */
)
{
   SCIPdebugMessage("calling SCIPlpiReadState()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   HIGHS_CALL( lpi->highs->readBasis(fname) );
   return SCIP_OKAY;
}

/** writes LPi state (i.e. basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
      SCIP_LPI*             lpi,            /**< LP interface structure */
      const char*           fname           /**< file name */
)
{
   SCIPdebugMessage("calling SCIPlpiWriteState()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   HIGHS_CALL( lpi->highs->writeBasis(fname) );
   return SCIP_OKAY;
}

/**@} */

/*
 * LP Pricing Norms Methods
 */

/**@name LP Pricing Norms Methods */
/**@{ */

/** stores LPi pricing norms information
 *  @todo Could storing norm information improve warm start performance in HiGHS?
 */
SCIP_RETCODE SCIPlpiGetNorms(
      SCIP_LPI*             lpi,                /**< LP interface structure */
      BMS_BLKMEM*           blkmem,             /**< block memory */
      SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information */
)
{
   SCIPdebugMessage("calling SCIPlpiGetNorms()\n");

   assert(lpi != NULL);
   assert(lpinorms != NULL);

   (*lpinorms) = NULL;

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
   SCIPdebugMessage("calling SCIPlpiSetNorms()\n");

   assert(lpi != NULL);
   assert(lpinorms == NULL);

   /* no work necessary */
   return SCIP_OKAY;
}

/** frees pricing norms information */
SCIP_RETCODE SCIPlpiFreeNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information, or NULL */
   )
{
   SCIPdebugMessage("calling SCIPlpiFreeNorms()\n");

   assert(lpi != NULL);
   assert(lpinorms == NULL);

   /* no work necessary */
   return SCIP_OKAY;
}

/**@} */

/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiGetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetIntpar()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(ival != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      *ival = (int) lpi->fromscratch;
      break;
   case SCIP_LPPAR_LPINFO:
      {
         bool bool_ival;
         HIGHS_CALL( lpi->highs->getOptionValue("output_flag", bool_ival) );
         *ival = bool_ival;
      }
      break;
   case SCIP_LPPAR_SCALING:
      HIGHS_CALL( lpi->highs->getOptionValue("simplex_scale_strategy", *ival) );
      assert(*ival == 0 || *ival == 2 || *ival == 4); /* values used in SCIPlpiSetIntpar() */
      if( *ival <= 0 )
         *ival = 0;
      else if( *ival <= 2 )
         *ival = 1;
      else
         *ival = 2;
      break;
   case SCIP_LPPAR_PRESOLVING:
      {
         std::string presolve;
         HIGHS_CALL( lpi->highs->getOptionValue("presolve", presolve) );
         assert(presolve == "on" || presolve == "off"); /* values used in SCIPlpiSetIntpar() */
         *ival = (presolve == "on");
      }
      break;
   case SCIP_LPPAR_PRICING:
      *ival = (int)lpi->pricing; /* store pricing method in LPI struct */
      break;
   case SCIP_LPPAR_THREADS:
      *ival = lpi->nthreads;
      break;
   case SCIP_LPPAR_LPITLIM:
      HIGHS_CALL( lpi->highs->getOptionValue("simplex_iteration_limit", *ival) );
      break;
   case SCIP_LPPAR_RANDOMSEED:
      HIGHS_CALL( lpi->highs->getOptionValue("random_seed", *ival) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
  }

   return SCIP_OKAY;
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiSetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiSetIntpar()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      assert(ival == TRUE || ival == FALSE);
      lpi->fromscratch = (SCIP_Bool) ival;
      break;
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      HIGHS_CALL( lpi->highs->setOptionValue("output_flag", (bool) ival) );
      break;
   case SCIP_LPPAR_SCALING:
      assert(ival >= 0 && ival <= 2);
      if( ival == 0 )
         /* off */
         HIGHS_CALL( lpi->highs->setOptionValue("simplex_scale_strategy", 0) );
      else if( ival == 1 )
         /* forced equilibration */
         HIGHS_CALL( lpi->highs->setOptionValue("simplex_scale_strategy", 2) );
      else
         /* max. value scaling */
         HIGHS_CALL( lpi->highs->setOptionValue("simplex_scale_strategy", 4) );
      break;
   case SCIP_LPPAR_PRESOLVING:
      assert(ival == TRUE || ival == FALSE);
      HIGHS_CALL( lpi->highs->setOptionValue("presolve", ival ? "on" : "off") );
      break;
   case SCIP_LPPAR_PRICING:
      lpi->pricing = (SCIP_PRICING)ival;
      switch( lpi->pricing )
      {
      case SCIP_PRICING_LPIDEFAULT:
      case SCIP_PRICING_PARTIAL:
      case SCIP_PRICING_AUTO:
         HIGHS_CALL( lpi->highs->setOptionValue("simplex_primal_edge_weight_strategy", -1) );
         HIGHS_CALL( lpi->highs->setOptionValue("simplex_dual_edge_weight_strategy", -1) );
         break;
      case SCIP_PRICING_DEVEX:
         HIGHS_CALL( lpi->highs->setOptionValue("simplex_primal_edge_weight_strategy", 1) );
         HIGHS_CALL( lpi->highs->setOptionValue("simplex_dual_edge_weight_strategy", 1) );
         break;
      case SCIP_PRICING_FULL:
      case SCIP_PRICING_STEEP:
      case SCIP_PRICING_STEEPQSTART:
         HIGHS_CALL( lpi->highs->setOptionValue("simplex_primal_edge_weight_strategy", 2) );
         HIGHS_CALL( lpi->highs->setOptionValue("simplex_dual_edge_weight_strategy", 2) );
         break;
      default:
         return SCIP_LPERROR;
      }
      break;
   case SCIP_LPPAR_THREADS:
      lpi->nthreads = ival;
      break;
   case SCIP_LPPAR_LPITLIM:
      HIGHS_CALL( lpi->highs->setOptionValue("simplex_iteration_limit", ival) );
      break;
   case SCIP_LPPAR_RANDOMSEED:
      HIGHS_CALL( lpi->highs->setOptionValue("random_seed", ival) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiGetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
)
{
   SCIPdebugMessage("calling SCIPlpiGetRealpar()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(dval != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      HIGHS_CALL( lpi->highs->getOptionValue("primal_feasibility_tolerance", *dval) );
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      HIGHS_CALL( lpi->highs->getOptionValue("dual_feasibility_tolerance", *dval) );
      break;
   case SCIP_LPPAR_LPTILIM:
      HIGHS_CALL( lpi->highs->getOptionValue("time_limit", *dval) );
      break;
   case SCIP_LPPAR_OBJLIM:
      HIGHS_CALL( lpi->highs->getOptionValue("objective_bound", *dval) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiSetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiSetRealpar()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      /* Primal feasibility tolerance cannot be smaller than 1e-10 */
      dval = MAX(dval, 1e-10);
      HIGHS_CALL( lpi->highs->setOptionValue("primal_feasibility_tolerance", dval) );
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      /* Dual feasibility tolerance cannot be smaller than 1e-10 */
      dval = MAX(dval, 1e-10);
      HIGHS_CALL( lpi->highs->setOptionValue("dual_feasibility_tolerance", dval) );
      break;
   case SCIP_LPPAR_LPTILIM:
      HIGHS_CALL( lpi->highs->setOptionValue("time_limit", dval) );
      break;
   case SCIP_LPPAR_OBJLIM:
      HIGHS_CALL( lpi->highs->setOptionValue("objective_bound", dval) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** interrupts the currently ongoing lp solve or disables the interrupt */
SCIP_RETCODE SCIPlpiInterrupt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             interrupt           /**< TRUE if interrupt should be set, FALSE if it should be disabled */
   )
{
   SCIPdebugMessage("calling SCIPlpiInterrupt()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   /* not implemented */
   return SCIP_OKAY;
}


/**@} */

/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiInfinity(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiInfinity()\n");

   assert(lpi != NULL);

   return kHighsInf;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val
   )
{
   SCIPdebugMessage("calling SCIPlpiIsInfinity()\n");

   assert(lpi != NULL);

   return val >= kHighsInf;
}

/**@} */

/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadLP()\n");

   assert(lpi != NULL);
   assert(lpi->highs != NULL);
   assert(fname != NULL);

   HIGHS_CALL( lpi->highs->readModel(fname) );

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiWriteLP()\n");
   assert(fname != NULL);

   assert(lpi != NULL);
   assert(lpi->highs != NULL);

   HIGHS_CALL( lpi->highs->writeModel(fname) );

   return SCIP_OKAY;
}

/**@} */
