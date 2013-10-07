/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lpi_clp.cpp
 * @ingroup LPIS
 * @brief  LP interface for Clp
 * @author Stefan Heinz
 * @author Marc Pfetsch
 * @author John Forrest
 *
 *
 * Notes on this interface:
 *
 * - Currently, Clp (Version 1.10) supports two ways of adding rows/columns from arrays: One uses a
 *   length array that for each row/column specifies the number of nonzeros to be added. The second
 *   uses the @p beg array that gives the starting index for each row/column. We use the latter
 *   variant. Since for LPI there should be no gaps in the corresponding arrays, i.e., every entry in
 *   @p val and @a ind gives a nonzero entry, one can switch between the two formats. With the current
 *   Clp implementation both formats involve an overhead:
 *
 *    - For the @p beg variant, Clp gets the end of the array from the last position in @p beg
 *      (i.e., the entry one after the last row/column) and we have to copy and extend @p beg for this
 *      purpose. In the matrix implementation a length information is then again computed.
 *
 *    - For the @p length variant, Clp computes the number of elements from this length variant and
 *      there exists no matrix implementation that uses the length information, i.e., it is recomputed
 *      again.
 *
 *    Concluding: the implementation of Clp/CoinPackeMatrix could be improved. The functions
 *    affected by this are SCIPlpiLoadColLP(), SCIPlpiAddCols(), SCIPlpiAddRows()
 *
 * - In former versions Clp used an "auxiliary model" that allows to save time when the model is
 *   scaled. This is discarded from version higher than 1.8.2.
 *
 * - Clp allows the setting of several special flags. These are now set when the FASTMIP option in
 *   SCIP is true. We tried to use the best settings, while still working correctly, see
 *   setFastmipClpParameters(). These settings probably have to be adapted to future Clp
 *   versions. Maybe more possibilities will appear.
 *
 * - At several places this interface corrects the return value of some Clp functions, e.g.,
 *   isProvenPrimalInfeasible(). Currently (version 1.10) no change in the Clp functions will be made,
 *   but this might change in the future.
 */
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <ClpSimplex.hpp>
#include <ClpPrimalColumnSteepest.hpp>
#include <ClpDualRowSteepest.hpp>
#include <CoinIndexedVector.hpp>
#include <ClpConfig.h>
#ifndef CLP_VERSION
#include <config_clp.h>
#define CLP_VERSION VERSION
#endif

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#include "scip/lpi.h"
#include "scip/bitencode.h"
#include "scip/pub_message.h"

/* do defines for windows directly her to make the lpi more independent*/
#if defined(_WIN32) || defined(_WIN64)
#define snprintf _snprintf
#endif

/* for debugging: alternatingly write files "debug_[p|d]_[0|1].mps" after each run - use with care! */
#ifdef LPI_CLP_DEBUG_WRITE_FILES
static int fileNr = 0;
#endif

/* bound for accepting primal or dual sum of infeasibilities */
#define SUMINFEASBOUND   1.0e-3

/** LP interface for Clp */
struct SCIP_LPi
{
   ClpSimplex*           clp;                        /**< Clp simiplex solver class */
   int*                  cstat;                      /**< array for storing column basis status */
   int*                  rstat;                      /**< array for storing row basis status */
   int                   cstatsize;                  /**< size of cstat array */
   int                   rstatsize;                  /**< size of rstat array */
   bool                  startscratch;               /**< start from scratch? */
   bool                  presolving;                 /**< preform preprocessing? */
   SCIP_PRICING          pricing;                    /**< SCIP pricing setting  */
   bool                  validFactorization;         /**< whether we have a valid factorization in clp */
   SCIP_Bool             solved;                     /**< was the current LP solved? */
   bool                  setFactorizationFrequency;  /**< store whether the factorization frequency is set */
   SCIP_Bool             fastmip;                    /**< are fast mip settings turned on */
};






/** Definitions for storing basis status  (copied from lpi_spx.cpp) */
typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
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
   assert(lpi != 0);

   if( num > lpi->cstatsize )
   {
      int newsize;

      newsize = MAX(2*lpi->cstatsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->cstat, newsize) );
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
   assert(lpi != 0);

   if( num > lpi->rstatsize )
   {
      int newsize;

      newsize = MAX(2*lpi->rstatsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->rstat, newsize) );
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
   return (ncols+(int)COLS_PER_PACKET-1)/(int)COLS_PER_PACKET;
}

/** returns the number of packets needed to store row packet information */
static
int rowpacketNum(
   int                   nrows               /**< number of rows to store */
   )
{
   return (nrows+(int)ROWS_PER_PACKET-1)/(int)ROWS_PER_PACKET;
}

/** store row and column basis status in a packed LPi state object */
static
void lpistatePack(
   SCIP_LPISTATE*       lpistate,            /**< pointer to LPi state data */
   const int*           cstat,               /**< basis status of columns in unpacked format */
   const int*           rstat                /**< basis status of rows in unpacked format */
   )
{
   assert(lpistate != 0);
   assert(lpistate->packcstat != 0);
   assert(lpistate->packrstat != 0);

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
   assert(lpistate != 0);
   assert(lpistate->packcstat != 0);
   assert(lpistate->packrstat != 0);

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
   assert(lpistate != 0);
   assert(blkmem != 0);
   assert(ncols >= 0);
   assert(nrows >= 0);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, lpistate) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packcstat, colpacketNum(ncols)) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packrstat, rowpacketNum(nrows)) );

   return SCIP_OKAY;
}

/** frees LPi state information */
static
void lpistateFree(
   SCIP_LPISTATE**       lpistate,           /**< pointer to LPi state information (like basis information) */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(blkmem != 0);
   assert(lpistate != 0);
   assert(*lpistate != 0);

   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packcstat, colpacketNum((*lpistate)->ncols));
   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packrstat, rowpacketNum((*lpistate)->nrows));
   BMSfreeBlockMemory(blkmem, lpistate);
}





/*
 * local methods
 */

/** marks the current LP to be unsolved */
static
void invalidateSolution(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   lpi->solved = FALSE;
}

/** set factorization frequency */
static
void setFactorizationFrequency(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   /* set the factorization frequency only once */
   if ( lpi->setFactorizationFrequency )
      return;

   lpi->clp->defaultFactorizationFrequency();
   lpi->setFactorizationFrequency = true;
}

/** this methods sets parameters of Clp */
static
void setFastmipClpParameters(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   lpi->fastmip = TRUE;

   /** Perturbation:
    *  50  - switch on perturbation
    *  100 - auto perturb if takes too long (1.0e-6 largest nonzero)
    *  101 - we are perturbed
    *  102 - don't try perturbing again
    *  - default is 100
    *  - others are for playing
    *
    * for Clp 1.8 stable: 50 seems to be 10% faster than 100
    */
   lpi->clp->setPerturbation(50);

   /** Special options description from ClpModell.hpp:
    *       1 - Don't keep changing infeasibility weight
    *       2 - Keep nonLinearCost round solves
    *       4 - Force outgoing variables to exact bound (primal)
    *       8 - Safe to use dense initial factorization
    *      16 - Just use basic variables for operation if column generation
    *      32 - Create ray even in BAB
    *      64 - Treat problem as feasible until last minute (i.e. minimize infeasibilities)
    *     128 - Switch off all matrix sanity checks
    *     256 - No row copy
    *     512 - If not in values pass, solution guaranteed, skip as much as possible
    *    1024 - In branch and bound
    *    2048 - Don't bother to re-factorize if < 20 iterations
    *    4096 - Skip some optimality checks
    *    8192 - Do Primal when cleaning up primal
    *   16384 - In fast dual (so we can switch off things)
    *   32768 - called from Osi
    *   65536 - keep arrays around as much as possible (also use maximumR/C)
    *  131072 - transposeTimes is -1.0 and can skip basic and fixed
    *  262144 - extra copy of scaled matrix
    *  524288 - Clp fast dual
    * 1048576 - don't need to finish dual (can return 3)
    *  NOTE   - many applications can call Clp but there may be some short cuts
    *           which are taken which are not guaranteed safe from all applications.
    *           Vetted applications will have a bit set and the code may test this
    *           At present I expect a few such applications - if too many I will
    *           have to re-think.  It is up to application owner to change the code
    *           if she/he needs these short cuts.  I will not debug unless in Coin
    *           repository.  See COIN_CLP_VETTED comments.
    *  0x01000000 is Cbc (and in branch and bound)
    *  0x02000000 is in a different branch and bound
    *
    *  Comments:
    *       2 - nonlinear costs are used in primal for infeasibility weight
    *       4 - in anti-degeneracy operations can move variables just off a bound
    *       8 - means dense nucleus in factorization - normally not safe in first factorization as
    *           singularity handling is not useful. Is switched on if going from dual to primal or vv.
    *      16 - Used for "real" column generation
    *      64 - Good idea, since in B&B most problems are feasible.
    *     128 - Assumes user will not create tiny or duplicate elements.
    *     256 - Normally Clp keeps a scaled row copy for speed. For very large problems you might want to turn it off.
    *     512 - Means nonbasic variables should be at bounds and basis will be reasonable.
    *    4096 - Skip some optimality checks
    *    8192 - If the primal has a perturbed problem and needs to clean up, it normally uses dual - but in some cases can be better to use primal.
    *   32768 - Just switches off some messages e.g. empty problem.
    *  131072 - used internally
    *  262144 - Normally Clp has unscaled column copy of matrix - this makes an extra scaled copy.
    *  524288 - used internally
    * 1048576 - only set by fastDual
    * 0x02000000 - main point: does allow use of disaster handler
    *
    * Cbc seems to use the following special options:
    * lpi->clp->setSpecialOptions(64|128|1024|2048|4096|32768|262144|0x01000000);
    * Sometimes 512+8192 and 8192 or 8 are used as well.
    */

   // 2048 does not seem to work
   // 65536 does not seem to work
   // 262144 does not seem to work

#ifndef NDEBUG
   // in debug mode: leave checks on
   lpi->clp->setSpecialOptions(32|64|512|1024|32768);
#else
   lpi->clp->setSpecialOptions(32|64|128|512|1024|4096|32768);
#endif

   // 8192 bit - don't even think of using primal if user asks for dual (and vv)
   lpi->clp->setMoreSpecialOptions(8192 | lpi->clp->moreSpecialOptions());

   // let memory grow only (do not shrink) - [needs specialOptions & 65536 != 0]
   // does not seem to work
   //lpi->clp->setPersistenceFlag(1);
}

/** this methods sets parameters of Clp */
static
void unsetFastmipClpParameters(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   lpi->fastmip = FALSE;

   // reset to default value:
   lpi->clp->setPerturbation(100);

   // turn off special options:
   lpi->clp->setSpecialOptions(0);

   // turn off memory enlargement
   lpi->clp->setPersistenceFlag(0);
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
const char* SCIPlpiGetSolverName(
   void
   )
{
   // Currently Clp has no function to get version, so we hard code it ...
   return "Clp "CLP_VERSION;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return "COIN-OR Linear Programming Solver developed by J. Forrest et.al. (projects.coin-or.org/Clp)";
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   return (void*) lpi->clp;
}
/**@} */




/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** creates an LP problem object */
SCIP_RETCODE SCIPlpiCreate(
   SCIP_LPI**            lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   assert(lpi != 0);

   SCIPdebugMessage("calling SCIPlpiCreate()\n");

   // create lpi object
   SCIP_ALLOC( BMSallocMemory(lpi) );
   (*lpi)->clp = new ClpSimplex();
   (*lpi)->cstat = 0;
   (*lpi)->rstat = 0;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->startscratch = true;
   (*lpi)->pricing = SCIP_PRICING_LPIDEFAULT;
   (*lpi)->validFactorization = false;
   (*lpi)->setFactorizationFrequency = false;
   (*lpi)->fastmip = FALSE;
   invalidateSolution(*lpi);

   // if you want to use saveModel()
   // (*lpi)->clp->setLengthNames(255);

   // set pricing routines

   // for primal:
   // 0 is exact devex,
   // 1 full steepest,
   // 2 is partial exact devex
   // 3 switches between 0 and 2 depending on factorization
   // 4 starts as partial dantzig/devex but then may switch between 0 and 2.
   // - currently (Clp 1.8stable) default is 3
   ClpPrimalColumnSteepest primalSteepest;
   (*lpi)->clp->setPrimalColumnPivotAlgorithm(primalSteepest);

   // for dual:
   // 0 is uninitialized,
   // 1 full,
   // 2 is partial uninitialized,
   // 3 starts as 2 but may switch to 1.
   // - currently (Clp 1.8stable) default is 3
   ClpDualRowSteepest dualSteepest;
   (*lpi)->clp->setDualRowPivotAlgorithm(dualSteepest);

   // set problem name
   (*lpi)->clp->setStrParam(ClpProbName, std::string(name) );

   // set objective sense: SCIP values are the same as the ones for Clp
   (*lpi)->clp->setOptimizationDirection(objsen);

   // turn off output by default
   (*lpi)->clp->setLogLevel(0);

   // turn off scaling by default
   (*lpi)->clp->scaling(0);

   /* set default pricing */
   SCIP_CALL( SCIPlpiSetIntpar(*lpi, SCIP_LPPAR_PRICING, (int)(*lpi)->pricing) );

   return SCIP_OKAY;
}


/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != 0);
   assert(*lpi != 0);
   assert((*lpi)->clp != 0);

   SCIPdebugMessage("calling SCIPlpiFree()\n");

   /* free LP */
   delete (*lpi)->clp;

   /* free memory */
   BMSfreeMemoryArrayNull(&(*lpi)->cstat);
   BMSfreeMemoryArrayNull(&(*lpi)->rstat);
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
   char**                colnames,           /**< column names, or 0 */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand sides of rows */
   const SCIP_Real*      rhs,                /**< right hand sides of rows */
   char**                rownames,           /**< row names, or 0 */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array */
   const int*            ind,                /**< row indices of constraint matrix entries */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
   SCIPdebugMessage("calling SCIPlpiLoadColLP()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(lhs != 0);
   assert(rhs != 0);
   assert( nnonz > beg[ncols-1] );

   invalidateSolution(lpi);

   ClpSimplex* clp = lpi->clp;

   // copy beg-array
   int* mybeg = NULL;
   SCIP_ALLOC( BMSallocMemoryArray(&mybeg, ncols + 1) );
   BMScopyMemoryArray(mybeg, beg, ncols);
   mybeg[ncols] = nnonz;   // add additional entry at end

   // load problem
   clp->loadProblem(ncols, nrows, mybeg, ind, val, lb, ub, obj, lhs, rhs);
   BMSfreeMemoryArray( &mybeg );

   // set objective sense
   clp->setOptimizationDirection(objsen);

   // copy column and rownames if necessary
   if ( colnames || rownames )
   {
      std::vector<std::string> columnNames(ncols);
      std::vector<std::string> rowNames(nrows);
      if (colnames)
      {
         for (int j = 0; j < ncols; ++j)
            columnNames[j].assign(colnames[j]);
      }
      if (rownames)
      {
         for (int i = 0; i < ncols; ++i)
            rowNames[i].assign(rownames[i]);
      }
      clp->copyNames(rowNames, columnNames);
   }

   return SCIP_OKAY;
}


/** adds columns to the LP */
SCIP_RETCODE SCIPlpiAddCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   const SCIP_Real*      obj,                /**< objective function values of new columns */
   const SCIP_Real*      lb,                 /**< lower bounds of new columns */
   const SCIP_Real*      ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or 0 */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array, or 0 if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or 0 if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or 0 if nnonz == 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiAddCols()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(obj != 0);
   assert(lb != 0);
   assert(ub != 0);
   assert(nnonz == 0 || beg != 0);
   assert(nnonz == 0 || ind != 0);
   assert(nnonz == 0 || val != 0);
   assert(nnonz >= 0);
   assert(ncols >= 0);

   invalidateSolution(lpi);

   // store number of columns for later
   int numCols = lpi->clp->getNumCols();

   // copy beg-array (if not 0)
   int* mybeg = NULL;
   SCIP_ALLOC( BMSallocMemoryArray(&mybeg, ncols + 1) );

   // if columns are not empty
   if ( nnonz != 0 )
   {
      BMScopyMemoryArray(mybeg, beg, ncols);
      mybeg[ncols] = nnonz;   // add additional entry at end

      // add columns
      lpi->clp->addColumns(ncols, lb, ub, obj, mybeg, ind, val);
   }
   else
   {
      for (int j = 0; j <= ncols; ++j)
         mybeg[j] = 0;
      // add empty columns
      lpi->clp->addColumns(ncols, lb, ub, obj, mybeg, 0, 0);
   }
   BMSfreeMemoryArray(&mybeg);

   // copy columnnames if necessary
   if ( colnames )
   {
      std::vector<std::string> columnNames(ncols);
      for (int j = 0; j < ncols; ++j)
         columnNames[j].assign(colnames[j]);
      lpi->clp->copyColumnNames(columnNames, numCols, numCols + ncols);
   }

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

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->clp->numberColumns());

   invalidateSolution(lpi);

   // Current Clp version (1.8) can't delete a range of columns; we have to use deleteColumns (see SCIPlpiDelColset)
   int num = lastcol-firstcol+1;
   int* which = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &which, num) );;

   // fill array with interval
   for (int j = firstcol; j <= lastcol; ++j)
      which[j - firstcol] = j;

   lpi->clp->deleteColumns(num, which);
   BMSfreeMemoryArray( &which );

   return SCIP_OKAY;
}


/** deletes columns from SCIP_LPI; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelColset()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(dstat != 0);

   invalidateSolution(lpi);

   // transform dstat information
   int ncols = lpi->clp->getNumCols();
   int* which = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &which, ncols) );
   int cnt = 0;
   for (int j = 0; j < ncols; ++j)
   {
      if ( dstat[j] == 1 )
         which[cnt++] = j;
   }
   lpi->clp->deleteColumns(cnt, which);
   BMSfreeMemoryArray(&which);

   // update dstat
   cnt = 0;
   for (int j = 0; j < ncols; ++j)
   {
      if ( dstat[j] == 1 )
      {
         dstat[j] = -1;
         ++cnt;
      }
      else
         dstat[j] = j - cnt;
   }

   return SCIP_OKAY;
}


/** adds rows to the LP */
SCIP_RETCODE SCIPlpiAddRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or 0 */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each row in ind- and val-array, or 0 if nnonz == 0 */
   const int*            ind,                /**< column indices of constraint matrix entries, or 0 if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or 0 if nnonz == 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiAddRows()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(lhs != 0);
   assert(rhs != 0);
   assert(nnonz == 0 || beg != 0);
   assert(nnonz == 0 || ind != 0);
   assert(nnonz == 0 || val != 0);

   invalidateSolution(lpi);

   // store number of rows for later use
   int numRows = lpi->clp->getNumRows();

   int* mybeg = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &mybeg, nrows + 1) );

   if ( nnonz != 0 )
   {
      // copy beg-array
      BMScopyMemoryArray( mybeg, beg, nrows);
      mybeg[nrows] = nnonz;   // add additional entry at end

      // add rows
      lpi->clp->addRows(nrows, lhs, rhs, mybeg, ind, val);
   }
   else
   {
      // add empty rows
      for (int i = 0; i <= nrows; ++i)
         mybeg[i] = 0;
      lpi->clp->addRows(nrows, lhs, rhs, mybeg, 0, 0);
   }
   BMSfreeMemoryArray( &mybeg );

   // copy rownames if necessary
   if ( rownames )
   {
      std::vector<std::string> rowNames(nrows);
      for (int j = 0; j < nrows; ++j)
         rowNames[j].assign(rownames[j]);
      lpi->clp->copyRowNames(rowNames, numRows, numRows + nrows);
   }

   return SCIP_OKAY;
}


/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelRows() (number: %d)\n", lastrow-firstrow+1);

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->clp->numberRows());

   invalidateSolution(lpi);

   // Current Clp version (1.8) can't delete a range of rows; we have to use deleteRows (see SCIPlpiDelRowset)
   int num = lastrow-firstrow+1;
   int* which = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &which, num) );

   // fill array with interval
   for (int i = firstrow; i <= lastrow; ++i)
      which[i - firstrow] = i;

   lpi->clp->deleteRows(num, which);

   BMSfreeMemoryArray( &which );

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

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(dstat != 0);

   invalidateSolution(lpi);

   // transform dstat information
   int nrows = lpi->clp->getNumRows();
   int* which = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &which, nrows) );
   int cnt = 0;
   for (int i = 0; i < nrows; ++i)
   {
      if ( dstat[i] == 1 )
         which[cnt++] = i;
   }
   lpi->clp->deleteRows(cnt, which);
   BMSfreeMemoryArray( &which );

   // update dstat
   cnt = 0;
   for (int i = 0; i < nrows; ++i)
   {
      if ( dstat[i] == 1 )
      {
         dstat[i] = -1;
         ++cnt;
      }
      else
         dstat[i] = i - cnt;
   }

   return SCIP_OKAY;
}


/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiClear()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   invalidateSolution(lpi);

   // We use the resize(0,0) to get rid of the model but keep all other settings
   lpi->clp->resize(0,0);

   return SCIP_OKAY;
}


/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiChgBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgBounds()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(ind != 0);
   assert(lb != 0);
   assert(ub != 0);

   invalidateSolution(lpi);

   ClpSimplex* clp = lpi->clp;

   /* We currently employ the following bug fix: the solution vector is modified to be set to the
    * corresponding bounds. This avoids one error in Clp - maybe fixed in later versions. */
   double* sol = lpi->clp->primalColumnSolution();
   const double* colLower = lpi->clp->getColLower();
   const double* colUpper = lpi->clp->getColUpper();

   for (int j = 0; j < ncols; ++j)
   {
      clp->setColumnBounds(ind[j], lb[j], ub[j]);
      if ( sol != 0 )
      {
         if( clp->statusExists() )
         {
            assert( colLower != 0 );
            assert( colUpper != 0 );
            int k = ind[j];
            switch ( clp->getColumnStatus(k) )
            {
               case ClpSimplex::isFree:
               case ClpSimplex::superBasic:
                  sol[j] = 0.0;
                  break;
               case ClpSimplex::atUpperBound:
                  sol[k] = colUpper[k];
                  assert( colUpper[k] == ub[j] );
                  break;
               case ClpSimplex::isFixed:
               case ClpSimplex::atLowerBound:
                  sol[k] = colLower[k];
                  assert( colLower[k] == lb[j] );
                  break;
               default:;
            }
         }
         else
         { /* workaround: if there is no status, we assume something */
            sol[j] = 0.0;
         }
      }
   }

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

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(ind != 0);
   assert(lhs != 0);
   assert(rhs != 0);

   invalidateSolution(lpi);

   ClpSimplex* clp = lpi->clp;

   for (int i = 0; i < nrows; ++i)
      clp->setRowBounds(ind[i], lhs[i], rhs[i]);

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

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= row && row < lpi->clp->numberRows());
   assert(0 <= col && col < lpi->clp->numberColumns());

   invalidateSolution(lpi);

   lpi->clp->matrix()->modifyCoefficient(row, col, newval);

   return SCIP_OKAY;
}


/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgObjsen()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   invalidateSolution(lpi);

   // set objective sense: SCIP values are the same as the ones for Clp
   lpi->clp->setOptimizationDirection(objsen);

   return SCIP_OKAY;
}


/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiChgObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   SCIP_Real*            obj                 /**< new objective values for columns */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgObj()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(ind != 0);
   assert(obj != 0);

   invalidateSolution(lpi);

   ClpSimplex* clp = lpi->clp;

   // updates whatsChanged in Clp (bound checking in Clp)
   for( int j = 0; j < ncols; ++j )
      clp->setObjCoeff(ind[j], obj[j]);  // inlined version of clp->setObjectiveCoefficient(ind[j], obj[j]);

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

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(scaleval != 0.0);
   assert(0 <= row && row <= lpi->clp->numberRows() );

   invalidateSolution(lpi);

   // Note: if the scaling should be performed because of numerical stability,
   // there are other more effective methods in Clp to adjust the scaling values
   // for each row.

   ClpSimplex* clp = lpi->clp;

   // adjust the sides
   double* lhs = clp->rowLower();
   double* rhs = clp->rowUpper();

   double lhsval = lhs[row];
   if( lhsval > -COIN_DBL_MAX )
      lhsval *= scaleval;
   else if( scaleval < 0.0 )
      lhsval = COIN_DBL_MAX;
   double rhsval = rhs[row];
   if( rhsval < COIN_DBL_MAX)
      rhsval *= scaleval;
   else if( scaleval < 0.0 )
      rhsval = -COIN_DBL_MAX;
   if( scaleval < 0.0 )
   {
      SCIP_Real oldlhs = lhsval;
      lhsval = rhsval;
      rhsval = oldlhs;
   }
   lhs[row] = lhsval;    // change values directly into Clp data!
   rhs[row] = rhsval;

   // apply scaling ...

   // WARNING: the following is quite expensive:
   // We have to loop over the matrix to find the row entries.
   // For columns we can do better, see @c SCIPlpiScaleCol.
   CoinPackedMatrix* M = clp->matrix();
   assert( M->getNumCols() == clp->numberColumns() );

   const CoinBigIndex* beg = M->getVectorStarts();
   const int* length = M->getVectorLengths();
   const int* ind = M->getIndices();
   double* val = M->getMutableElements();

   for (int j = 0; j < M->getNumCols(); ++j)
   {
      for (CoinBigIndex k = beg[j]; k < beg[j] + length[j]; ++k)
      {
	 if (ind[k] == row)
	    val[k] *= scaleval;
      }
   }

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

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(scaleval != 0.0);
   assert(0 <= col && col <= lpi->clp->numberColumns() );

   invalidateSolution(lpi);

   // Note: if the scaling should be performed because of numerical stability,
   // there are other more effective methods in Clp to adjust the scaling values
   // for each column.

   ClpSimplex* clp = lpi->clp;

   // adjust the objective coefficients
   double* objvec = clp->objective();          // we have direct access to the data of Clp!
   objvec[col] *= scaleval;                    // adjust the objective function value

   // adjust the bounds
   double* lb = clp->columnLower();
   double* ub = clp->columnUpper();
   double lbval = lb[col];
   double ubval = ub[col];

   if( lbval > -COIN_DBL_MAX )
      lbval /= scaleval;
   else if( scaleval < 0.0 )
      lbval = COIN_DBL_MAX;
   if( ubval < COIN_DBL_MAX )
      ubval /= scaleval;
   else if( scaleval < 0.0 )
      ubval = -COIN_DBL_MAX;
   if( scaleval < 0.0 )
   {
      SCIP_Real oldlb = lbval;
      lbval = ubval;
      ubval = oldlb;
   }
   lb[col] = lbval;        // directly adjust values into Clp data
   ub[col] = ubval;

   // apply scaling directly to matrix (adapted from ClpPackedMatrix::reallyScale)
   // See also ClpModel::gutsOfScaling ...
   CoinPackedMatrix* M = clp->matrix();
   assert( M->getNumCols() == clp->numberColumns() );

   const CoinBigIndex* beg = M->getVectorStarts();
   const int* length = M->getVectorLengths();
   double* val = M->getMutableElements();
   for (CoinBigIndex k = beg[col]; k < beg[col] + length[col]; ++k)
      val[k] *= scaleval;

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

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(nrows != 0);

   *nrows = lpi->clp->numberRows();

   return SCIP_OKAY;
}


/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNCols()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(ncols != 0);

   *ncols = lpi->clp->numberColumns();

   return SCIP_OKAY;
}


/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNNonz()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(nnonz != 0);

   *nnonz = lpi->clp->getNumElements();

   return SCIP_OKAY;
}


/** gets columns from LP problem object; the arrays have to be large enough to store all values
 *  Either both, lb and ub, have to be 0, or both have to be non-0,
 *  either nnonz, beg, ind, and val have to be 0, or all of them have to be non-0.
 */
SCIP_RETCODE SCIPlpiGetCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or 0 */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or 0 */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or 0 */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or 0 */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or 0 */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetCols()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->clp->numberColumns());

   ClpSimplex* clp = lpi->clp;

   // get lower and upper bounds for the variables
   assert( (lb != 0 && ub != 0) || (lb == 0 && ub == 0) );
   if ( lb != 0 )
   {
      const double* colLower = clp->getColLower();    // Here we can use the const versions (see SCIPchgBounds)
      const double* colUpper = clp->getColUpper();

      BMScopyMemoryArray( lb, colLower + firstcol, (lastcol - firstcol + 1));
      BMScopyMemoryArray( ub, colUpper + firstcol, (lastcol - firstcol + 1));
   }

   assert( nnonz != 0 || beg == 0);
   assert( nnonz != 0 || ind == 0);
   assert( nnonz != 0 || val == 0);

   if ( nnonz != 0 )
   {
      CoinPackedMatrix* M = clp->matrix();
      assert( M != 0 );
      assert( M->getNumCols() == clp->numberColumns() );

      const CoinBigIndex* Mbeg = M->getVectorStarts();   // can use const versions
      const int* Mlength = M->getVectorLengths();
      const int* Mind = M->getIndices();
      const double* Mval = M->getElements();

      *nnonz = 0;
      // can we use memcpy for the whole set (requires that columns are stored sequentially)
      for (int j = firstcol; j <= lastcol; ++j)
      {
         beg[j-firstcol] = *nnonz;

         BMScopyMemoryArray( (ind + (*nnonz)), Mind + Mbeg[j], Mlength[j]);
         BMScopyMemoryArray( (val + (*nnonz)), Mval + Mbeg[j], Mlength[j]);

         (*nnonz) += Mlength[j];
      }
   }

   return SCIP_OKAY;
}


/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be 0, or both have to be non-0,
 *  either nnonz, beg, ind, and val have to be 0, or all of them have to be non-0.
 */
SCIP_RETCODE SCIPlpiGetRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_Real*            lhs,                /**< buffer to store left hand side vector, or 0 */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or 0 */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or 0 */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or 0 */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or 0 */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetRows()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->clp->numberRows());

   ClpSimplex* clp = lpi->clp;
   assert( (lhs != 0 && rhs != 0) || (lhs == 0 && rhs == 0) );
   if ( lhs != 0 )
   {
      const double* rowLower = clp->getRowLower();    // Here we can use the const versions (see SCIPchgSides)
      const double* rowUpper = clp->getRowUpper();

      BMScopyMemoryArray( lhs, rowLower + firstrow, (lastrow - firstrow + 1) );
      BMScopyMemoryArray( rhs, rowUpper + firstrow, (lastrow - firstrow + 1) );
   }

   assert( nnonz != 0 || beg == 0);
   assert( nnonz != 0 || ind == 0);
   assert( nnonz != 0 || val == 0);

   if ( nnonz != 0 )
   {
      ClpMatrixBase* M = clp->rowCopy();   // get row view on matrix
      if ( M == 0 ) // can happen e.g. if no LP was solved yet ...
	 M = clp->clpMatrix()->reverseOrderedCopy();
      assert( M != 0 );
      assert( M->getNumRows() == clp->numberRows() );

      const CoinBigIndex* Mbeg = M->getVectorStarts();
      const int* Mlength = M->getVectorLengths();
      const int* Mind = M->getIndices();
      const double* Mval = M->getElements();

      *nnonz = 0;
      for( int i = firstrow; i <= lastrow; ++i )
      {
         beg[i-firstrow] = *nnonz;
         for( CoinBigIndex k = Mbeg[i]; k < Mbeg[i] + Mlength[i]; ++k )
         {
            ind[*nnonz] = Mind[k];
            val[*nnonz] = Mval[k];
            (*nnonz)++;
         }
      }
   }

   return SCIP_OKAY;
}


/** gets column names */
SCIP_RETCODE SCIPlpiGetColNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get name from LP */
   int                   lastcol,            /**< last column to get name from LP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) */
   char*                 namestorage,        /**< storage for col names */
   int                   namestoragesize,    /**< size of namestorage (if 0, storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   )
{
   SCIPerrorMessage("SCIPlpiGetColNames() has not been implemented yet.\n");
   return SCIP_LPERROR;
}


/** gets row names */
SCIP_RETCODE SCIPlpiGetRowNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get name from LP */
   int                   lastrow,            /**< last row to get name from LP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) */
   char*                 namestorage,        /**< storage for row names */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   )
{
   SCIPerrorMessage("SCIPlpiGetRowNames() has not been implemented yet.\n");
   return SCIP_LPERROR;
}


/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   SCIPdebugMessage("calling SCIPlpiIgnoreInstability()\n");

   assert(lpi != NULL);
   assert(lpi->clp != NULL);

   /* unstable situations cannot be ignored */
   *success = FALSE;

   return SCIP_OKAY;
}


/** gets the objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   assert( lpi != NULL );
   assert( lpi->clp != NULL );
   assert( objsen != NULL );

   // Clp direction of optimization (1 - minimize, -1 - maximize, 0 - ignore)
   if ( lpi->clp->getObjSense() < 0 )
      *objsen = SCIP_OBJSEN_MAXIMIZE;
   else
      *objsen = SCIP_OBJSEN_MINIMIZE;

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

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->clp->numberColumns());
   assert(vals != 0);

   const double* obj = lpi->clp->getObjCoefficients();    // Here we can use the const versions (see SCIPchgObj)

   BMScopyMemoryArray( vals, obj + firstcol, (lastcol - firstcol + 1) );

   return SCIP_OKAY;
}


/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiGetBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective value for */
   int                   lastcol,            /**< last column to get objective value for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or 0 */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBounds()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->clp->numberColumns());

   if ( lbs != 0 )
   {
      const double* colLower = lpi->clp->getColLower();    // Here we can use the const versions (see SCIPchgBounds)
      BMScopyMemoryArray( lbs, colLower + firstcol, (lastcol - firstcol + 1) );
   }

   if ( ubs != 0 )
   {
      const double* colUpper = lpi->clp->getColUpper();
      BMScopyMemoryArray( ubs, colUpper + firstcol, (lastcol - firstcol + 1) );
   }

   return SCIP_OKAY;
}


/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiGetSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or 0 */
   SCIP_Real*            rhss                /**< array to store right hand side values, or 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSides()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->clp->numberRows());

   if ( lhss != 0 )
   {
      const double* rowLower = lpi->clp->getRowLower();    // Here we can use the const versions (see SCIPchgSides)
      BMScopyMemoryArray( lhss, rowLower + firstrow, (lastrow - firstrow + 1) );
   }

   if ( rhss != 0 )
   {
      const double* rowUpper = lpi->clp->getRowUpper();
      BMScopyMemoryArray( rhss,  rowUpper + firstrow, (lastrow - firstrow + 1) );
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

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= col && col < lpi->clp->numberColumns());
   assert(0 <= row && row < lpi->clp->numberRows());
   assert(val != 0);

   *val = lpi->clp->matrix()->getCoefficient(row, col);

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
   assert(lpi != 0);
   assert(lpi->clp != 0);

   SCIPdebugMessage("calling Clp primal(): %d cols, %d rows\n", lpi->clp->numberColumns(), lpi->clp->numberRows());

#ifdef LPI_CLP_DEBUG_WRITE_FILES
   char filename[255];
   snprintf(filename, 255, "debug_p_%d.mps", fileNr);
   fileNr = fileNr % 2;
   SCIPlpiWriteLP(lpi, filename);
   SCIPdebugMessage("Wrote file <%s>\n", filename);
#endif

   invalidateSolution(lpi);

   // intialize factorization freq. depending on model size - applied only once
   setFactorizationFrequency(lpi);

   // if we want to construct a new basis
   if ( lpi->startscratch )
   {
      lpi->clp->allSlackBasis(true);   // reset basis
      lpi->validFactorization = false;
   }

   /** startFinishOptions - bits
    *  1 - do not delete work areas and factorization at end
    *  2 - use old factorization if same number of rows
    *  4 - skip as much initialization of work areas as possible (work in progress)
    *
    *  4 does not seem to work.
    */
   int startFinishOptions = 1;
   if ( lpi->validFactorization )
      startFinishOptions = startFinishOptions | 2;

   /** Primal algorithm */
   int status = lpi->clp->primal(0, startFinishOptions);

#ifdef LPI_CLP_DEBUG_WRITE_FILES
   char basisname[255];
   snprintf(basisname, 255, "debug_p_%d.bas", fileNr);
   SCIP_CALL( SCIPlpiWriteState(lpi, basisname) );
   SCIPdebugMessage("Wrote basis file <%s>\n", basisname);
   ++fileNr; /* not increased above! */
   fileNr = fileNr % 2;
#endif

   lpi->validFactorization = true;
   lpi->solved = TRUE;

   // Unfortunately the status of Clp is hard coded ...
   // -1 - did not run
   //  0 - optimal
   //  1 - primal infeasible
   //  2 - dual infeasible
   //  3 - stopped on iterations or time
   //  4 - stopped due to errors
   //  5 - stopped by event handler
   assert( status != -1 );      // did not run should not occur
   assert( status != 5 );       // begin stopped by event handler should not occur
   if ( status == 4 || status == 5 || status == -1 )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}


/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != 0);
   assert(lpi->clp != 0);

   SCIPdebugMessage("calling Clp dual(): %d cols, %d rows\n", lpi->clp->numberColumns(), lpi->clp->numberRows());

#ifdef LPI_CLP_DEBUG_WRITE_FILES
   char filename[255];
   snprintf(filename, 255, "debug_d_%d.mps", fileNr);
   SCIPlpiWriteLP(lpi, filename);
   SCIPdebugMessage("Wrote file <%s>\n", filename);
   snprintf(filename, 255, "debug_d_%d.sav", fileNr);
   // lpi->clp->saveModel(filename);
   SCIPdebugMessage("Wrote file <%s>\n", filename);
#endif

   invalidateSolution(lpi);

   // intialize factorization freq. depending on model size - applied only once
   setFactorizationFrequency(lpi);

   // if we want to construct a new basis
   if( lpi->startscratch )
   {
      lpi->clp->allSlackBasis(true);   // reset basis
      lpi->validFactorization = false;
   }

   /** startFinishOptions - bits
    *  1 - do not delete work areas and factorization at end
    *  2 - use old factorization if same number of rows
    *  4 - skip as much initialization of work areas as possible (work in progress)
    *
    *  4 does not seem to work.
    */
   int startFinishOptions = 1;
   if ( lpi->validFactorization )
      startFinishOptions = startFinishOptions | 2;

   /** Dual algorithm */
   int status = lpi->clp->dual(0, startFinishOptions);

#ifdef LPI_CLP_DEBUG_WRITE_FILES
   char basisname[255];
   snprintf(basisname, 255, "debug_d_%d.bas", fileNr);
   SCIP_CALL( SCIPlpiWriteState(lpi, basisname) );
   SCIPdebugMessage("Wrote basis file <%s>\n", basisname);
   ++fileNr; /* not increased above! */
   fileNr = fileNr % 2;
#endif

   lpi->validFactorization = true;
   lpi->solved = TRUE;

   // Unfortunately the status of Clp is hard coded ...
   // -1 - did not run
   //  0 - optimal
   //  1 - primal infeasible
   //  2 - dual infeasible
   //  3 - stopped on iterations or time
   //  4 - stopped due to errors
   //  5 - stopped by event handler
   assert( status != -1 );      // did not run should not occur
   assert( status != 5 );       // begin stopped by event handler should not occur
   if ( status == 4 || status == 5 || status == -1 )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}


/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                 /**< LP interface structure */
   SCIP_Bool             crossover            /**< perform crossover */
   )
{
   assert(lpi != 0);
   assert(lpi->clp != 0);

   SCIPdebugMessage("calling Clp barrier(): %d cols, %d rows\n", lpi->clp->numberColumns(), lpi->clp->numberRows());

   invalidateSolution(lpi);

   // Check whether we have a factorization, if yes destroy it (Clp doesn't like it ...)
   /*
   if (lpi->haveFactorization)
      lpi->clp->finish();
   */

   // call barrier
   int status = lpi->clp->barrier(crossover);
   lpi->solved = TRUE;

   // We may need to call ClpModel::status()

   // Unfortunately the status of Clp is hard coded ...
   // -1 - did not run
   //  0 - optimal
   //  1 - primal infeasible
   //  2 - dual infeasible
   //  3 - stopped on iterations or time
   //  4 - stopped due to errors
   //  5 - stopped by event handler
   assert( status != -1 );      // did not run should not occur
   assert( status != 5 );       // begin stopped by event handler should not occur
   if ( status == 4 || status == 5 || status == -1 )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}

/** start strong branching - call before any strongbranching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   // currently do nothing; in the future: use code as in OSI
   return SCIP_OKAY;
}

/** end strong branching - call after any strongbranching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   // currently do nothing; in the future: use code as in OSI
   return SCIP_OKAY;
}

/** performs strong branching iterations on one arbitrary candidate */
static
SCIP_RETCODE lpiStrongbranch(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current primal solution value of column */
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
   SCIPdebugMessage("calling SCIPlpiStrongbranch() on variable %d (%d iterations)\n", col, itlim);

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(down != 0);
   assert(up != 0);
   assert(downvalid != 0);
   assert(upvalid != 0);

   ClpSimplex* clp = lpi->clp;

   // set up output arrays
   int ncols = clp->numberColumns();
   assert( 0 <= col && col < ncols );
   double** outputSolution = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &outputSolution, 2) );
   SCIP_ALLOC( BMSallocMemoryArray( &outputSolution[0], ncols) );
   SCIP_ALLOC( BMSallocMemoryArray( &outputSolution[1], ncols) );

   int* outputStatus = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &outputStatus, 2) );

   int* outputIterations = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &outputIterations, 2) );

   // set iteration limit
   int iterlimit = clp->maximumIterations();
   clp->setMaximumIterations(itlim);

   // store objective value
   double objval = clp->objectiveValue();

   // store special options for later reset
   int specialoptions = clp->specialOptions();

   /** Clp special options:
    *       1 - Don't keep changing infeasibility weight
    *       2 - Keep nonLinearCost round solves
    *       4 - Force outgoing variables to exact bound (primal)
    *       8 - Safe to use dense initial factorization
    *      16 - Just use basic variables for operation if column generation
    *      32 - Create ray even in BAB
    *      64 - Treat problem as feasible until last minute (i.e. minimize infeasibilities)
    *     128 - Switch off all matrix sanity checks
    *     256 - No row copy
    *     512 - If not in values pass, solution guaranteed, skip as much as possible
    *    1024 - In branch and bound
    *    2048 - Don't bother to re-factorize if < 20 iterations
    *    4096 - Skip some optimality checks
    *    8192 - Do Primal when cleaning up primal
    *   16384 - In fast dual (so we can switch off things)
    *   32768 - called from Osi
    *   65536 - keep arrays around as much as possible (also use maximumR/C)
    *  131072 - transposeTimes is -1.0 and can skip basic and fixed
    *  262144 - extra copy of scaled matrix
    *  524288 - Clp fast dual
    * 1048576 - don't need to finish dual (can return 3)
    *  NOTE   - many applications can call Clp but there may be some short cuts
    *           which are taken which are not guaranteed safe from all applications.
    *           Vetted applications will have a bit set and the code may test this
    *           At present I expect a few such applications - if too many I will
    *           have to re-think.  It is up to application owner to change the code
    *           if she/he needs these short cuts.  I will not debug unless in Coin
    *           repository.  See COIN_CLP_VETTED comments.
    *  0x01000000 is Cbc (and in branch and bound)
    *  0x02000000 is in a different branch and bound
    *
    *  2048 does not seem to work
    *  262144 does not seem to work
    */
#ifndef NDEBUG
   // in debug mode: leave checks on
   clp->setSpecialOptions(64|512|1024);
#else
   clp->setSpecialOptions(64|128|512|1024|4096);
#endif

   /* 'startfinish' options for strong branching:
    *  1 - do not delete work areas and factorization at end
    *  2 - use old factorization if same number of rows
    *  4 - skip as much initialization of work areas as possible
    *      (based on whatsChanged in clpmodel.hpp) ** work in progress
    *
    *  4 does not seem to work in strong branching ...
    */
   int startFinishOptions = 1;
   if ( lpi->validFactorization )
      startFinishOptions = startFinishOptions | 2;

   // set new lower and upper bounds for variable
   *down = EPSCEIL(psol - 1.0, 1e-06);
   *up   = EPSFLOOR(psol + 1.0, 1e-06);

   /** For strong branching.  On input lower and upper are new bounds while
    *  on output they are change in objective function values (>1.0e50
    *  infeasible).  Return code is
    *   0 if nothing interesting,
    *  -1 if infeasible both ways and
    *  +1 if infeasible one way (check values to see which one(s))
    *  -2 if bad factorization
    * Solutions are filled in as well - even down, odd up - also status and number of iterations
    *
    * The bools are:
    *   bool stopOnFirstInfeasible
    *   bool alwaysFinish
    *
    * At the moment: we need alwaysFinish to get correct bounds.
    */
   //int res = clp->strongBranching(1, &col, up, down, outputSolution, outputStatus, outputIterations, false, false, startFinishOptions);
   int res = clp->strongBranching(1, &col, up, down, outputSolution, outputStatus, outputIterations, false, true, startFinishOptions);

   // reset special options
   clp->setSpecialOptions(specialoptions);

   lpi->validFactorization = true;

   *down += objval;
   *up += objval;

   // The bounds returned by CLP seem to be valid using the above options
   *downvalid = TRUE;
   *upvalid = TRUE;

   // correct iteration count
   if (iter)
      *iter = outputIterations[0] + outputIterations[1];

   // reset iteration limit
   clp->setMaximumIterations(iterlimit);

   // free local memory
   BMSfreeMemoryArray( &outputStatus );
   BMSfreeMemoryArray( &outputIterations );
   BMSfreeMemoryArray( &outputSolution[1] );
   BMSfreeMemoryArray( &outputSolution[0] );
   BMSfreeMemoryArray( &outputSolution );

   if ( res == -2 )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}

/** performs strong branching iterations on given arbitrary candidates */
static
SCIP_RETCODE lpiStrongbranches(
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
   SCIPdebugMessage("calling SCIPlpiStrongbranches() on %d variables (%d iterations)\n", ncols, itlim);

   assert( lpi != 0 );
   assert( lpi->clp != 0 );
   assert( cols != 0 );
   assert( psols != 0 );
   assert( down != 0 );
   assert( up != 0 );
   assert( downvalid != 0 );
   assert( upvalid != 0 );

   ClpSimplex* clp = lpi->clp;

   // set up output arrays
   int n = clp->numberColumns();
   assert( 0 < ncols && ncols <= n );
   double** outputSolution = NULL;
   SCIP_ALLOC( BMSallocMemoryArray( &outputSolution, 2*ncols) );
   for (int j = 0; j < 2*ncols; ++j)
   {
      SCIP_ALLOC( BMSallocMemoryArray( &(outputSolution[j]), n) );
   }

   int* outputStatus = NULL;
   SCIP_ALLOC( BMSallocMemoryArray(&outputStatus, 2*ncols) );

   int* outputIterations = NULL;
   SCIP_ALLOC( BMSallocMemoryArray(&outputIterations, 2*ncols) );

   // set iteration limit
   int iterlimit = clp->maximumIterations();
   clp->setMaximumIterations(itlim);

   // store objective value
   double objval = clp->objectiveValue();

   // store special options for later reset
   int specialoptions = clp->specialOptions();

   /** Clp special options:
    *       1 - Don't keep changing infeasibility weight
    *       2 - Keep nonLinearCost round solves
    *       4 - Force outgoing variables to exact bound (primal)
    *       8 - Safe to use dense initial factorization
    *      16 - Just use basic variables for operation if column generation
    *      32 - Create ray even in BAB
    *      64 - Treat problem as feasible until last minute (i.e. minimize infeasibilities)
    *     128 - Switch off all matrix sanity checks
    *     256 - No row copy
    *     512 - If not in values pass, solution guaranteed, skip as much as possible
    *    1024 - In branch and bound
    *    2048 - Don't bother to re-factorize if < 20 iterations
    *    4096 - Skip some optimality checks
    *    8192 - Do Primal when cleaning up primal
    *   16384 - In fast dual (so we can switch off things)
    *   32768 - called from Osi
    *   65536 - keep arrays around as much as possible (also use maximumR/C)
    *  131072 - transposeTimes is -1.0 and can skip basic and fixed
    *  262144 - extra copy of scaled matrix
    *  524288 - Clp fast dual
    * 1048576 - don't need to finish dual (can return 3)
    *  NOTE   - many applications can call Clp but there may be some short cuts
    *           which are taken which are not guaranteed safe from all applications.
    *           Vetted applications will have a bit set and the code may test this
    *           At present I expect a few such applications - if too many I will
    *           have to re-think.  It is up to application owner to change the code
    *           if she/he needs these short cuts.  I will not debug unless in Coin
    *           repository.  See COIN_CLP_VETTED comments.
    *  0x01000000 is Cbc (and in branch and bound)
    *  0x02000000 is in a different branch and bound
    *
    *  2048 does not seem to work
    *  262144 does not seem to work
    */
#ifndef NDEBUG
   // in debug mode: leave checks on
   clp->setSpecialOptions(64|512|1024);
#else
   clp->setSpecialOptions(64|128|512|1024|4096);
#endif

   /* 'startfinish' options for strong branching:
    *  1 - do not delete work areas and factorization at end
    *  2 - use old factorization if same number of rows
    *  4 - skip as much initialization of work areas as possible
    *      (based on whatsChanged in clpmodel.hpp) ** work in progress
    *
    *  4 does not seem to work in strong branching ...
    */
   int startFinishOptions = 1;
   if ( lpi->validFactorization )
      startFinishOptions = startFinishOptions | 2;

   // set new lower and upper bounds for variables
   for (int j = 0; j < ncols; ++j)
   {
      assert( 0 <= cols[j] && cols[j] < n );
      down[j] = EPSCEIL(psols[j] - 1.0, 1e-06);
      up[j]   = EPSFLOOR(psols[j] + 1.0, 1e-06);

      // The bounds returned by CLP seem to be valid using the above options
      downvalid[j] = TRUE;
      upvalid[j] = TRUE;
   }

   /** For strong branching.  On input lower and upper are new bounds while
    *  on output they are change in objective function values (>1.0e50
    *  infeasible).  Return code is
    *   0 if nothing interesting,
    *  -1 if infeasible both ways and
    *  +1 if infeasible one way (check values to see which one(s))
    *  -2 if bad factorization
    * Solutions are filled in as well - even down, odd up - also status and number of iterations
    *
    * The bools are:
    *   bool stopOnFirstInfeasible
    *   bool alwaysFinish
    *
    * At the moment: we need alwaysFinish to get correct bounds.
    */
   int res = clp->strongBranching(ncols, cols, up, down, outputSolution, outputStatus, outputIterations, false, true, startFinishOptions);

   // reset special options
   clp->setSpecialOptions(specialoptions);

   lpi->validFactorization = true;

   for (int j = 0; j < ncols; ++j)
   {
      down[j] += objval;
      up[j] += objval;

      // correct iteration count
      if (iter)
         *iter += outputIterations[2*j] + outputIterations[2*j+1];

      BMSfreeMemoryArray(&outputSolution[2*j]);
      BMSfreeMemoryArray(&outputSolution[2*j+1]);
   }

   // reset iteration limit
   clp->setMaximumIterations(iterlimit);

   // free local memory
   BMSfreeMemoryArray( &outputStatus );
   BMSfreeMemoryArray( &outputIterations );
   BMSfreeMemoryArray( &outputSolution );

   if ( res == -2 )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}

/** performs strong branching iterations on one @b fractional candidate */
SCIP_RETCODE SCIPlpiStrongbranchFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current primal solution value of column */
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
   /* pass call on to lpiStrongbranch() */
   SCIP_CALL( lpiStrongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter) );

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
   if ( iter != NULL )
      *iter = 0;

   /* pass call on to lpiStrongbranches() */
   SCIP_CALL( lpiStrongbranches(lpi, cols, ncols, psols, itlim, down, up, downvalid, upvalid, iter) );

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
   /* pass call on to lpiStrongbranch() */
   SCIP_CALL( lpiStrongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter) );

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
   if ( iter != NULL )
      *iter = 0;

   /* pass call on to lpiStrongbranches() */
   SCIP_CALL( lpiStrongbranches(lpi, cols, ncols, psols, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiWasSolved(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   return lpi->solved;
}

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_RETCODE SCIPlpiGetSolFeasibility(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSolFeasibility()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(primalfeasible != 0);
   assert(dualfeasible != 0);

   if ( lpi->clp->primalFeasible() )
      *primalfeasible = TRUE;
   else
      *primalfeasible = FALSE;

   if ( lpi->clp->dualFeasible() )
      *dualfeasible = TRUE;
   else
      *dualfeasible = FALSE;

   // say feasible if deviation is small
   if (lpi->clp->status()==0 && ( ! (*primalfeasible) || ! (*dualfeasible)) ) 
   {
      if ( !(*primalfeasible) && lpi->clp->sumPrimalInfeasibilities() < SUMINFEASBOUND ) 
      {
         lpi->clp->setNumberPrimalInfeasibilities(0);
         *primalfeasible = TRUE;
      }
      if ( !(*dualfeasible) && lpi->clp->sumDualInfeasibilities() < SUMINFEASBOUND)
      {
         lpi->clp->setNumberDualInfeasibilities(0);
         *dualfeasible = TRUE;
      }
   }

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

   assert(lpi != 0);
   assert(lpi->clp != 0);

   /* Clp seems to have a primal ray whenever it concludes "dual infeasible" (status == 2)
    * (but is not necessarily primal feasible), see ClpModel::unboundedRay(). */
   return ( lpi->clp->status() == 2 );
}


/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasPrimalRay()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   /* Clp seems to have a primal ray whenever it concludes "dual infeasible" (status == 2)
    * (but is not necessarily primal feasible), see ClpModel::unboundedRay(). */
   return ( lpi->clp->status() == 2 );
}


/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalUnbounded()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return ( lpi->clp->isProvenDualInfeasible() && lpi->clp->primalFeasible() );
}


/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalInfeasible()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   /* Should return ClpModel::isProvenPrimalInfeasible() (which returns status == 1), but the
    * following is correct (Clp will not be changed). The secondaryStatus is 1 if the dual simplex
    * detects an objective limit exceedence. The primal simplex has no such detection (will never
    * stop with objective limit exceedence). Hence we are infeasible only if status == 1 and we have
    * not stopped due to the objective limit. */
   return ( lpi->clp->status() == 1 && (lpi->clp->secondaryStatus() == 0 || lpi->clp->secondaryStatus() == 6) );
}


/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalFeasible()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return ( lpi->clp->primalFeasible() );
}


/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExistsDualRay()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   /* Clp assumes to have a dual ray whenever it concludes "primal infeasible" and the algorithm was
    * the dual simplex, (but is not necessarily dual feasible), see ClpModel::infeasibilityRay */
   return ( lpi->clp->status() == 1 && lpi->clp->secondaryStatus() == 0 && lpi->clp->algorithm() < 0 );
}


/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasDualRay()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   /* Clp assumes to have a dual ray whenever it concludes "primal infeasible" and the algorithm was
    * the dual simplex, (but is not necessarily dual feasible), see ClpModel::infeasibilityRay */
   if ( lpi->clp->rayExists() )
   {
      if ( lpi->clp->status() == 1 && lpi->clp->secondaryStatus() == 0 && lpi->clp->algorithm() < 0)
         return TRUE;
      else 
      {
         if ( lpi->clp->status() != 2 || lpi->clp->algorithm() <= 0 ) 
            lpi->clp->deleteRay();
         return FALSE;
      }
   } 
   else
      return FALSE;
}


/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualUnbounded()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   /* The dual seems to be unbounded if the status is 1 (primal unbounded), the secondaryStatus is
    * not 1 (i.e., the dual simplex has not stopped because of an objective limit exceedence), and
    * the dual is feasible. */
   return ( lpi->clp->status() == 1 && lpi->clp->secondaryStatus() == 0 && lpi->clp->dualFeasible() );
}


/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualInfeasible()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return ( lpi->clp->isProvenDualInfeasible() );
}


/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualFeasible()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return ( lpi->clp->dualFeasible() );
}


/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsOptimal()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   if ( SCIPlpiIsObjlimExc(lpi) )
      return FALSE;

   /* secondaryStatus == 6 means that the problem is empty */
   return( lpi->clp->isProvenOptimal() && (lpi->clp->secondaryStatus() == 0 || lpi->clp->secondaryStatus() == 6));
}


/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsStable()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   /*  We first check if status is ok, i.e., is one of the following:
    *   0 - optimal
    *   1 - primal infeasible
    *   2 - dual infeasible
    *   3 - stopped on iterations or time
    *   4 - stopped due to errors
    *   5 - stopped by event handler (virtual int ClpEventHandler::event())
    */

   /* Then we check the secondary status of Clp:
    *  0 - none
    *  1 - primal infeasible because dual limit reached OR (probably primal infeasible but can't prove it  - main status was 4)
    *  2 - scaled problem optimal - unscaled problem has primal infeasibilities
    *  3 - scaled problem optimal - unscaled problem has dual infeasibilities
    *  4 - scaled problem optimal - unscaled problem has primal and dual infeasibilities
    *  5 - giving up in primal with flagged variables
    *  6 - failed due to empty problem check
    *  7 - postSolve says not optimal
    *  8 - failed due to bad element check
    *  9 - status was 3 and stopped on time
    *  100 up - translation of enum from ClpEventHandler
    */
   SCIPdebugMessage("status: %d   secondary: %d\n", lpi->clp->status(), lpi->clp->secondaryStatus());
   assert( 0 <= lpi->clp->status() && lpi->clp->status() <= 5 );
   return( (lpi->clp->status() <= 3) && (lpi->clp->secondaryStatus() <= 1 || lpi->clp->secondaryStatus() == 6 || lpi->clp->secondaryStatus() == 9) );
}


/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsObjlimExc()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   /* if status == 1 (primal infeasible) and secondaryStatus == 1 then Clp hit the dual bound */
   if ( lpi->clp->status() == 1 )
   {
      if ( lpi->clp->secondaryStatus() == 1 )
	 return TRUE;
      else
	 return FALSE;
   }

   return ( lpi->clp->isObjectiveLimitTestValid() && (lpi->clp->isPrimalObjectiveLimitReached() || lpi->clp->isDualObjectiveLimitReached()) );

   /* The above code is equivalent to the following:
   if ( lpi->clp->status() == 0 || (lpi->clp->status() == 1 && lpi->clp->algorithm() < 0) || (lpi->clp->status() == 2 && lpi->clp->algorithm() > 0) )
   {
      return ( lpi->clp->isPrimalObjectiveLimitReached() || lpi->clp->isDualObjectiveLimitReached() );
   }
   */

   return FALSE;
}


/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsIterlimExc()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   /* status == 3 means that Clp stopped on time or iteration limit
    * secondary status == 9 means that status was 3 and Clp stopped on time */
   return ( lpi->clp->status() == 3 && lpi->clp->secondaryStatus() != 9 );
}


/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   /* status == 3 means that Clp stopped on time or iteration limit
    * secondary status == 9 means that status was 3 and Clp stopped on time */
   return ( lpi->clp->status() == 3 && lpi->clp->secondaryStatus() == 9 );
}


/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetInternalStatus()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return lpi->clp->status();
}


/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetObjval()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(objval != 0);

   *objval = lpi->clp->objectiveValue();

   return SCIP_OKAY;
}


/** gets primal and dual solution vectors */
SCIP_RETCODE SCIPlpiGetSol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be 0 if not needed */
   SCIP_Real*            primsol,            /**< primal solution vector, may be 0 if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be 0 if not needed */
   SCIP_Real*            activity,           /**< row activity vector, may be 0 if not needed */
   SCIP_Real*            redcost             /**< reduced cost vector, may be 0 if not needed */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSol()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   ClpSimplex* clp = lpi->clp;
   if( objval != 0 )
      *objval = clp->objectiveValue();

   if( primsol != 0 )
   {
      const double* sol = clp->getColSolution();
      BMScopyMemoryArray( primsol, sol, clp->numberColumns() );
   }
   if( dualsol != 0 )
   {
      const double* dsol = clp->getRowPrice();
      BMScopyMemoryArray( dualsol, dsol, clp->numberRows() );
   }
   if( activity != 0 )
   {
      const double* act = clp->getRowActivity();
      BMScopyMemoryArray( activity, act, clp->numberRows() );
   }
   if( redcost != 0 )
   {
      const double* red = clp->getReducedCost();
      BMScopyMemoryArray( redcost, red, clp->numberColumns() );
   }

   return SCIP_OKAY;
}


/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetPrimalRay()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(ray != 0);

   /** Unbounded ray (NULL returned if none/wrong). Up to user to use delete [] on these arrays.  */
   const double* clpray = lpi->clp->unboundedRay();

   if ( clpray == 0 )
      return SCIP_LPERROR;

   BMScopyMemoryArray( ray, clpray, lpi->clp->numberColumns() );

   delete [] clpray;

   return SCIP_OKAY;
}

/** gets dual farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual farkas row multipliers */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetDualfarkas()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(dualfarkas != 0);

   /** Infeasibility ray (NULL returned if none/wrong). Up to user to use delete [] on these arrays.  */
   const double* dualray = lpi->clp->infeasibilityRay();

   if ( dualray == 0 )
      return SCIP_LPERROR;

   BMScopyMemoryArray( dualfarkas, dualray, lpi->clp->numberRows() );

   /* convert sign - this is needed for versions <= 1.10 */
   /*
   for (int j = 0; j < lpi->clp->numberRows(); ++j)
      dualfarkas[j] = -dualfarkas[j];
   */

   delete [] dualray;

   return SCIP_OKAY;
}


/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert(lpi != 0);
   assert(iterations != 0);

   *iterations = lpi->clp->numberIterations();

   return SCIP_OKAY;
}

/** gets information about the quality of an LP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The LPI should return SCIP_INVALID for *quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPlpiGetRealSolQuality(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPSOLQUALITY     qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   assert(lpi != NULL);
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
   int*                  cstat,              /**< array to store column basis status, or 0 */
   int*                  rstat               /**< array to store row basis status, or 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBase()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   ClpSimplex* clp = lpi->clp;

   // slower but easier to understand (and portable)
   if( rstat != 0 )
   {
      for( int i = 0; i < clp->numberRows(); ++i )
      {
	 switch ( clp->getRowStatus(i) )
	 {
	 case ClpSimplex::isFree:
            rstat[i] = SCIP_BASESTAT_ZERO;
            break;
	 case ClpSimplex::basic:
            rstat[i] = SCIP_BASESTAT_BASIC;
            break;
	 case ClpSimplex::atUpperBound:
            rstat[i] = SCIP_BASESTAT_UPPER;
            break;
	 case ClpSimplex::atLowerBound:
            rstat[i] = SCIP_BASESTAT_LOWER;
            break;
	 case ClpSimplex::superBasic:
            rstat[i] = SCIP_BASESTAT_ZERO;
            break;
	 case ClpSimplex::isFixed:
	    if (clp->getRowPrice()[i] > 0.0)
	       rstat[i] = SCIP_BASESTAT_LOWER;
	    else
	       rstat[i] = SCIP_BASESTAT_UPPER;
	    break;
	 default:
            SCIPerrorMessage("invalid basis status\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
	 }
      }
   }

   if( cstat != 0 )
   {
      for( int j = 0; j < clp->numberColumns(); ++j )
      {
	 switch ( clp->getColumnStatus(j) )
	 {
	 case ClpSimplex::isFree:
            cstat[j] = SCIP_BASESTAT_ZERO;
            break;
	 case ClpSimplex::basic:
            cstat[j] = SCIP_BASESTAT_BASIC;
            break;
	 case ClpSimplex::atUpperBound:
            cstat[j] = SCIP_BASESTAT_UPPER;
            break;
	 case ClpSimplex::atLowerBound:
            cstat[j] = SCIP_BASESTAT_LOWER;
            break;
	 case ClpSimplex::superBasic:
            cstat[j] = SCIP_BASESTAT_ZERO;
            break;
	 case ClpSimplex::isFixed:
	    if (clp->getReducedCost()[j] > 0.0)
	       cstat[j] = SCIP_BASESTAT_LOWER;
	    else
	       cstat[j] = SCIP_BASESTAT_UPPER;
	    break;
	 default: SCIPerrorMessage("invalid basis status\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
	 }
      }
   }

   return SCIP_OKAY;
}


/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiSetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   )
{
   SCIPdebugMessage("calling SCIPlpiSetBase()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   invalidateSolution(lpi);

   // Adapted from OsiClpSolverInterface::setBasisStatus

   ClpSimplex* clp = lpi->clp;
   clp->createStatus();

   const double* lhs = clp->getRowLower();
   const double* rhs = clp->getRowUpper();

   assert( rstat != 0 || clp->numberRows() == 0 );
   for( int i = 0; i < clp->numberRows(); ++i )
   {
      int status = rstat[i];
      assert( 0 <= status && status <= 3 );
      assert( lhs[i] > -COIN_DBL_MAX || status != SCIP_BASESTAT_LOWER); // can't be at lower bound
      assert( rhs[i] < COIN_DBL_MAX  || status != SCIP_BASESTAT_UPPER); // can't be at upper bound

      switch ( status )
      {
      case SCIP_BASESTAT_ZERO:
	 if (lhs[i] <= -COIN_DBL_MAX && rhs[i] >= COIN_DBL_MAX)
	    clp->setRowStatus(i, ClpSimplex::isFree);
	 else
	    clp->setRowStatus(i, ClpSimplex::superBasic);
	 break;
      case SCIP_BASESTAT_BASIC:
         clp->setRowStatus(i, ClpSimplex::basic);
         break;
      case SCIP_BASESTAT_UPPER:
         clp->setRowStatus(i, ClpSimplex::atUpperBound);
         break;
      case SCIP_BASESTAT_LOWER:
	 if ( EPSEQ(rhs[i], lhs[i], 1e-6) )   // if bounds are equal
	    clp->setRowStatus(i, ClpSimplex::isFixed);
	 else
	    clp->setRowStatus(i, ClpSimplex::atLowerBound);
	 break;
      default:
         SCIPerrorMessage("invalid basis status\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   const double* lb = clp->getColLower();
   const double* ub = clp->getColUpper();

   assert( cstat != 0 || clp->numberColumns() == 0 );
   for( int j = 0; j < clp->numberColumns(); ++j )
   {
      int status = cstat[j];
      assert( 0 <= status && status <= 3 );
      assert( lb[j] > -COIN_DBL_MAX || status != SCIP_BASESTAT_LOWER); // can't be at lower bound
      assert( ub[j] < COIN_DBL_MAX  || status != SCIP_BASESTAT_UPPER); // can't be at upper bound

      switch ( status )
      {
      case SCIP_BASESTAT_ZERO:
	 if (lb[j] <= -COIN_DBL_MAX && ub[j] >= COIN_DBL_MAX)
	    clp->setColumnStatus(j, ClpSimplex::isFree);
	 else
	    clp->setColumnStatus(j, ClpSimplex::superBasic);
	 break;
      case SCIP_BASESTAT_BASIC:
         clp->setColumnStatus(j, ClpSimplex::basic);
         break;
      case SCIP_BASESTAT_UPPER:
         clp->setColumnStatus(j, ClpSimplex::atUpperBound);
         break;
      case SCIP_BASESTAT_LOWER:
	 if ( EPSEQ(ub[j], lb[j], 1e-6) )
	    clp->setColumnStatus(j, ClpSimplex::isFixed);
	 else
	    clp->setColumnStatus(j, ClpSimplex::atLowerBound);
	 break;
      default:
         SCIPerrorMessage("invalid basis status\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   /** Whats changed since last solve.
    *  Is only used when startFinishOptions used in dual or primal.
    * Bit 1 - number of rows/columns has not changed (so work arrays valid)
    *     2 - matrix has not changed
    *     4 - if matrix has changed only by adding rows
    *     8 - if matrix has changed only by adding columns
    *    16 - row lbs not changed
    *    32 - row ubs not changed
    *    64 - column objective not changed
    *   128 - column lbs not changed
    *   256 - column ubs not changed
    *	512 - basis not changed (up to user to set this to 0)
    *	      top bits may be used internally
    */
   clp->setWhatsChanged(clp->whatsChanged() & (~512));

   return SCIP_OKAY;
}


/** returns the indices of the basic columns and rows */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBasisInd()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(bind != 0);

   ClpSimplex* clp = lpi->clp;
   int nrows = clp->numberRows();
   int ncols = clp->numberColumns();

   int* idx = NULL;
   SCIP_ALLOC( BMSallocMemoryArray(&idx, nrows) );

   /* If secondaryStatus == 6, clp says the LP is empty. Mose likely this happened, because the
      matrix is empty, i.e., all rows were redundant/empty. In this case, we construct a basis
      consisting of slack variables. */
   if ( clp->secondaryStatus() == 6 )
   {
      assert( clp->getNumElements() == 0 );
      for (int i = 0; i < nrows; ++i)
	 idx[i] = ncols + i;
   }
   else
      clp->getBasics(idx);

   for (int i = 0; i < nrows; ++i)
   {
      if ( idx[i] < ncols )
         bind[i] = idx[i];
      else
         bind[i] = -1 - (idx[i] - ncols);
   }

   BMSfreeMemoryArray(&idx);

   return SCIP_OKAY;
}


/** get dense row of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiGetBInvRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvRow()\n");

   assert( lpi != 0 );
   assert( lpi->clp != 0 );
   assert( coef != 0 );
   assert( 0 <= r && r <= lpi->clp->numberRows() );

   ClpSimplex* clp = lpi->clp;
   clp->getBInvRow(r, coef);

   return SCIP_OKAY;
}


/** get dense column of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiGetBInvCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the column */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvCol()\n");

   assert( lpi != 0 );
   assert( lpi->clp != 0 );
   assert( coef != 0 );
   assert( 0 <= c && c <= lpi->clp->numberRows() ); /* basis matrix is nrows * nrows */

   ClpSimplex* clp = lpi->clp;
   clp->getBInvCol(c, coef);

   return SCIP_OKAY;
}


/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or 0 */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvARow()\n");

   assert( lpi != 0 );
   assert( lpi->clp != 0 );
   assert( coef != 0 );
   assert( 0 <= r && r <= lpi->clp->numberRows() );

   ClpSimplex* clp = lpi->clp;
   clp->getBInvARow(r, coef, 0);

   return SCIP_OKAY;
}


/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvACol()\n");

   assert( lpi != 0 );
   assert( lpi->clp != 0 );
   assert( coef != 0 );
   assert( 0 <= c && c <= lpi->clp->numberColumns() );

   ClpSimplex* clp = lpi->clp;
   clp->getBInvACol(c, coef);

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

   assert(blkmem != 0);
   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(lpistate != 0);

   int ncols = lpi->clp->numberColumns();
   int nrows = lpi->clp->numberRows();
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


/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetState()
 */
SCIP_RETCODE SCIPlpiSetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           /*blkmem*/,         /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   )
{
   int lpncols;
   int lpnrows;
   int i;

   SCIPdebugMessage("calling SCIPlpiSetState()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(lpistate != 0);

   lpncols = lpi->clp->numberColumns();
   lpnrows = lpi->clp->numberRows();
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
      SCIP_Real bnd = (lpi->clp->getColLower())[i];
      if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         bnd = (lpi->clp->getColUpper())[i];
         if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
            lpi->cstat[i] = SCIP_BASESTAT_ZERO;  /* variable is free */
         else
            lpi->cstat[i] = SCIP_BASESTAT_UPPER; /* use finite upper bound */
      }
      else
         lpi->cstat[i] = SCIP_BASESTAT_LOWER;    /* use finite lower bound */
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

   assert(lpi != 0);
   assert(lpi->clp != 0);

   lpi->clp->allSlackBasis(true);
   lpi->validFactorization = false;

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

   assert(lpi != 0);
   assert(lpistate != NULL);

   if ( *lpistate != NULL )
      lpistateFree(lpistate, blkmem);

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{
   return (lpistate != NULL);
}

/** reads LP state (like basis information) from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,            /**< LP interface structure */
   const char*           fname           /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadState()\n");

   /** Read a basis from the given filename,
    *  returns -1 on file error, 0 if no values, 1 if values
    */
   if ( lpi->clp->readBasis(fname) < 0 )
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** writes LP state (like basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,            /**< LP interface structure */
   const char*           fname           /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiWriteState()\n");

   /** Write the basis in MPS format to the specified file.
    *  If writeValues true, writes values of structurals
    *  (and adds VALUES to end of NAME card)
    *
    *  parameters:
    *  - filename
    *  - bool writeValues
    *  - int formatType  (0 - normal, 1 - extra accuracy, 2 - IEEE hex)
    */
   if ( lpi->clp->writeBasis(fname, false, 0) )
      return SCIP_WRITEERROR;

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

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(ival != 0);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      *ival = lpi->startscratch;
      break;
   case SCIP_LPPAR_SCALING:
      if( lpi->clp->scalingFlag() != 0 )     // 0 -off, 1 equilibrium, 2 geometric, 3, auto, 4 dynamic(later)
	 *ival = TRUE;
      else
	 *ival = FALSE;
      break;
   case SCIP_LPPAR_PRICING:
      *ival = (int)lpi->pricing;          // store pricing method in LPI struct
      break;
   case SCIP_LPPAR_LPINFO:
      *ival = lpi->clp->logLevel() > 0 ? TRUE : FALSE;
      break;
   case SCIP_LPPAR_LPITLIM:
      *ival = lpi->clp->maximumIterations();
      break;
   case SCIP_LPPAR_FASTMIP:
      *ival = lpi->fastmip;
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

   assert(lpi != 0);
   assert(lpi->clp != 0);

   // Handle pricing separately ...
   if( type == SCIP_LPPAR_PRICING )
   {
      // for primal:
      // 0 is exact devex,
      // 1 full steepest,
      // 2 is partial exact devex
      // 3 switches between 0 and 2 depending on factorization
      // 4 starts as partial dantzig/devex but then may switch between 0 and 2.
      // - currently (Clp 1.8) default is 3

      // for dual:
      // 0 is uninitialized,
      // 1 full,
      // 2 is partial uninitialized,
      // 3 starts as 2 but may switch to 1.
      // - currently (Clp 1.8) default is 3
      lpi->pricing = (SCIP_PRICING)ival;
      int primalmode = 0;
      int dualmode = 0;
      switch( (SCIP_PRICING)ival )
      {
      case SCIP_PRICING_AUTO:
         primalmode = 3; dualmode = 3; break;
      case SCIP_PRICING_FULL:
         primalmode = 0; dualmode = 1; break;
      case SCIP_PRICING_LPIDEFAULT:
      case SCIP_PRICING_STEEP:
         primalmode = 1; dualmode = 0; break;
      case SCIP_PRICING_STEEPQSTART:
         primalmode = 1; dualmode = 2; break;
      case SCIP_PRICING_DEVEX:
         primalmode = 2; dualmode = 3; break;
      default:
         SCIPerrorMessage("unkown pricing parameter %d!\n", ival);
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
      ClpPrimalColumnSteepest primalpivot(primalmode);
      lpi->clp->setPrimalColumnPivotAlgorithm(primalpivot);
      ClpDualRowSteepest dualpivot(dualmode);
      lpi->clp->setDualRowPivotAlgorithm(dualpivot);
      return SCIP_OKAY;
   }

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      lpi->startscratch = ival;
      break;
   case SCIP_LPPAR_SCALING:
      lpi->clp->scaling(ival == TRUE ? 3 : 0);    // 0 -off, 1 equilibrium, 2 geometric, 3, auto, 4 dynamic(later));
      break;
   case SCIP_LPPAR_PRICING:
      /* should not happen - see above */
      SCIPABORT();
      return SCIP_LPERROR; /*lint !e527*/
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      /** Amount of print out:
       *  0 - none
       *  1 - just final
       *  2 - just factorizations
       *  3 - as 2 plus a bit more
       *  4 - verbose
       *  above that 8,16,32 etc just for selective SCIPdebug
       */
      if ( ival )
	 lpi->clp->setLogLevel(2);      // lpi->clp->setLogLevel(63);
      else
         lpi->clp->setLogLevel(0);
      break;
   case SCIP_LPPAR_LPITLIM:
      lpi->clp->setMaximumIterations(ival);
      break;
   case SCIP_LPPAR_FASTMIP:
      assert(ival == TRUE || ival == FALSE);
      if( ival )
         setFastmipClpParameters(lpi);
      else
         unsetFastmipClpParameters(lpi);
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

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(dval != 0);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      *dval = lpi->clp->primalTolerance();
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      *dval = lpi->clp->dualTolerance();
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:
      /**@todo add BARRIERCONVTOL parameter */
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_LOBJLIM:
      if ( lpi->clp->optimizationDirection() > 0 )   // if minimization
	 *dval = lpi->clp->primalObjectiveLimit();
      else
	 *dval = lpi->clp->dualObjectiveLimit();
      break;
   case SCIP_LPPAR_UOBJLIM:
      if ( lpi->clp->optimizationDirection() > 0 )   // if minimization
	 *dval = lpi->clp->dualObjectiveLimit();
      else
	 *dval = lpi->clp->primalObjectiveLimit();
      break;
   case SCIP_LPPAR_LPTILIM:
      *dval = lpi->clp->maximumSeconds();
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
   SCIPdebugMessage("setting parameter %d to value %g.\n", type, dval);
   assert(lpi != 0);
   assert(lpi->clp != 0);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      lpi->clp->setPrimalTolerance(dval);
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      lpi->clp->setDualTolerance(dval);
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:
      /**@todo add BARRIERCONVTOL parameter */
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_LOBJLIM:
      if ( lpi->clp->optimizationDirection() > 0 )   // if minimization
	 lpi->clp->setPrimalObjectiveLimit(dval);
      else
	 lpi->clp->setDualObjectiveLimit(dval);
      break;
   case SCIP_LPPAR_UOBJLIM:
      if ( lpi->clp->optimizationDirection() > 0 )   // if minimization
	 lpi->clp->setDualObjectiveLimit(dval);
      else
	 lpi->clp->setPrimalObjectiveLimit(dval);
      break;
   case SCIP_LPPAR_LPTILIM:
      lpi->clp->setMaximumSeconds(dval);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

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
   SCIP_LPI*             /*lpi*/             /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiInfinity()\n");

   return COIN_DBL_MAX;
}


/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             /*lpi*/,            /**< LP interface structure */
   SCIP_Real             val
   )
{
   SCIPdebugMessage("calling SCIPlpiIsInfinity()\n");

   return (val >= COIN_DBL_MAX);
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** returns, whether the given file exists */
static
SCIP_Bool fileExists(
   const char*           filename            /**< file name */
   )
{
   FILE* f;

   f = fopen(filename, "r");
   if( f == 0 )
      return FALSE;

   fclose(f);

   return TRUE;
}

/** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadLP()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   // WARNING: can only read mps files

   if ( !fileExists(fname) )
      return SCIP_NOFILE;

   /** read file in MPS format
    * parameters:
    * filename
    * bool keepNames
    * bool ignoreErrors
    */
   if ( lpi->clp->readMps(fname, true, false) )
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiWriteLP() - %s\n", fname);

   assert(lpi != 0);
   assert(lpi->clp != 0);

   /** write file in MPS format
    *  parameters:
    *  filename
    *  int formatType  (0 - normal, 1 - extra accuracy, 2 - IEEE hex)
    *  int numberAcross (1 or 2 values should be specified on every data line in the MPS file)
    *  double objSense
    */
   if ( lpi->clp->writeMps(fname, 0, 2, lpi->clp->optimizationDirection()) )
      return SCIP_WRITEERROR;

   return SCIP_OKAY;
}

/**@} */
