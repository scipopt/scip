/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: lpi_clp.cpp,v 1.6 2004/10/28 15:26:25 bzfpfets Exp $"

/**@file   lpi_clp.cpp
 * @brief  LP interface for Clp
 * @author Marc Pfetsch
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <ClpSimplex.hpp>
#include <ClpPrimalColumnSteepest.hpp>
#include <ClpDualRowSteepest.hpp>
#include <CoinIndexedVector.hpp>
#include <ClpFactorization.hpp>
#include <ClpSimplexDual.hpp>    // ????????????????

#include <iostream>
#include <cassert>
#include <vector>
#include <string>

extern "C" 
{
#include "lpi.h"
#include "bitencode.h"
#include "message.h"
}




/** LP interface for Clp */
struct LPi
{
   ClpSimplex*      clp;                /**< Clp simiplex solver class */
   int*             cstat;              /**< array for storing column basis status */
   int*             rstat;              /**< array for storing row basis status */
   int              cstatsize;          /**< size of cstat array */
   int              rstatsize;          /**< size of rstat array */
   bool             startscratch;       /**< start from scratch? */
   bool             presolving;         /**< preform preprocessing? */
   int              pricing;            /**< scip pricing setting  */
   bool             validFactorization; /**< whether we have a valid factorization in clp */
};






/** Definitions for storing basis status  (copied from lpi_spx.cpp) */
typedef DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET DUALPACKETSIZE
typedef DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET DUALPACKETSIZE

/** LPi state stores basis information */
struct LPiState
{
   int              ncols;              /**< number of LP columns */
   int              nrows;              /**< number of LP rows */
   COLPACKET*       packcstat;          /**< column basis status in compressed form */
   ROWPACKET*       packrstat;          /**< row basis status in compressed form */
};




/*
 * dynamic memory arrays
 */

/** resizes cstat array to have at least num entries */
static
RETCODE ensureCstatMem(
   LPI*             lpi,                /**< LP interface structure */
   int              num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != 0);

   if( num > lpi->cstatsize )
   {
      int newsize;

      newsize = MAX(2*lpi->cstatsize, num);
      ALLOC_OKAY( reallocMemoryArray(&lpi->cstat, newsize) );
      lpi->cstatsize = newsize;
   }
   assert(num <= lpi->cstatsize);

   return SCIP_OKAY;
}

/** resizes rstat array to have at least num entries */
static
RETCODE ensureRstatMem(
   LPI*             lpi,                /**< LP interface structure */
   int              num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != 0);

   if( num > lpi->rstatsize )
   {
      int newsize;

      newsize = MAX(2*lpi->rstatsize, num);
      ALLOC_OKAY( reallocMemoryArray(&lpi->rstat, newsize) );
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
   int              ncols               /**< number of columns to store */
   )
{
   return (ncols+COLS_PER_PACKET-1)/COLS_PER_PACKET;
}

/** returns the number of packets needed to store row packet information */
static 
int rowpacketNum(
   int              nrows               /**< number of rows to store */
   )
{
   return (nrows+ROWS_PER_PACKET-1)/ROWS_PER_PACKET;
}

/** store row and column basis status in a packed LPi state object */
static
void lpistatePack(
   LPISTATE*       lpistate,            /**< pointer to LPi state data */
   const int*      cstat,               /**< basis status of columns in unpacked format */
   const int*      rstat                /**< basis status of rows in unpacked format */
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
   const LPISTATE* lpistate,            /**< pointer to LPi state data */
   int*            cstat,               /**< buffer for storing basis status of columns in unpacked format */
   int*            rstat                /**< buffer for storing basis status of rows in unpacked format */
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
RETCODE lpistateCreate(
   LPISTATE**       lpistate,           /**< pointer to LPi state */
   MEMHDR*          memhdr,             /**< block memory */
   int              ncols,              /**< number of columns to store */
   int              nrows               /**< number of rows to store */
   )
{
   assert(lpistate != 0);
   assert(memhdr != 0);
   assert(ncols >= 0);
   assert(nrows >= 0);

   ALLOC_OKAY( allocBlockMemory(memhdr, lpistate) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, &(*lpistate)->packcstat, colpacketNum(ncols)) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, &(*lpistate)->packrstat, rowpacketNum(nrows)) );

   return SCIP_OKAY;
}

/** frees LPi state information */
static
void lpistateFree(
   LPISTATE**       lpistate,           /**< pointer to LPi state information (like basis information) */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(memhdr != 0);
   assert(lpistate != 0);
   assert(*lpistate != 0);

   freeBlockMemoryArray(memhdr, &(*lpistate)->packcstat, colpacketNum((*lpistate)->ncols));
   freeBlockMemoryArray(memhdr, &(*lpistate)->packrstat, rowpacketNum((*lpistate)->nrows));
   freeBlockMemory(memhdr, lpistate);
}





/*
 * LP Interface Methods
 */


/*
 * Miscellaneous Methods
 */

static char clpname[MAXSTRLEN];

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   // Currently Clp has no function to get version, so we hard code it ...
   sprintf(clpname, "Clp simplex solver");
   return clpname;
}

/**@} */




/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** creates an LP problem object */
RETCODE SCIPlpiCreate(
   LPI**            lpi,                /**< pointer to an LP interface structure */
   const char*      name,               /**< problem name */
   OBJSEN           objsen              /**< objective sense */
   )
{
   assert(lpi != 0);

   // create lpi object
   ALLOC_OKAY( allocMemory(lpi) );
   (*lpi)->clp = new ClpSimplex();
   (*lpi)->cstat = 0;
   (*lpi)->rstat = 0;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->startscratch = true;
   (*lpi)->pricing = SCIP_PRICING_AUTO;
   (*lpi)->validFactorization = false;

   // set pricing routines

   // for primal:
   // 0 is exact devex, 
   // 1 full steepest, 
   // 2 is partial exact devex
   // 3 switches between 0 and 2 depending on factorization
   // 4 starts as partial dantzig/devex but then may switch between 0 and 2.
   ClpPrimalColumnSteepest primalSteepest(3);
   (*lpi)->clp->setPrimalColumnPivotAlgorithm(primalSteepest);

   // for dual:
   // 0 is uninitialized, 1 full, 2 is partial uninitialized,
   // 3 starts as 2 but may switch to 1.
   ClpDualRowSteepest steep(3);
   (*lpi)->clp->setDualRowPivotAlgorithm(steep);

   // set problem name
   (*lpi)->clp->setStrParam(ClpProbName, std::string(name) );

   // set objective sense: SCIP values are the same as the ones for Clp
   (*lpi)->clp->setOptimizationDirection(objsen);

   // turn off output by default
   (*lpi)->clp->setLogLevel(0);

   return SCIP_OKAY;
}


/** deletes an LP problem object */
RETCODE SCIPlpiFree(
   LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != 0);
   assert(*lpi != 0);
   assert((*lpi)->clp != 0);

   /* free LP */
   delete (*lpi)->clp;

   /* free memory */
   freeMemoryArrayNull(&(*lpi)->cstat);
   freeMemoryArrayNull(&(*lpi)->rstat);
   freeMemory(lpi);

   return SCIP_OKAY;
}

/**@} */




/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies LP data with column matrix into LP solver */
RETCODE SCIPlpiLoadColLP(
   LPI*             lpi,                /**< LP interface structure */
   OBJSEN           objsen,             /**< objective sense */
   int              ncols,              /**< number of columns */
   const Real*      obj,                /**< objective function values of columns */
   const Real*      lb,                 /**< lower bounds of columns */
   const Real*      ub,                 /**< upper bounds of columns */
   char**           colnames,           /**< column names, or 0 */
   int              nrows,              /**< number of rows */
   const Real*      lhs,                /**< left hand sides of rows */
   const Real*      rhs,                /**< right hand sides of rows */
   char**           rownames,           /**< row names, or 0 */
   int              nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*       beg,                /**< start index of each column in ind- and val-array */
   const int*       ind,                /**< row indices of constraint matrix entries */
   const Real*      val                 /**< values of constraint matrix entries */
   )
{
   debugMessage("calling SCIPlpiLoadColLP()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(lhs != 0);
   assert(rhs != 0);
   assert( nnonz > beg[ncols-1] );

   ClpSimplex* clp = lpi->clp;

   // copy beg-array
   int* mybeg = new int [ncols+1];
   memcpy((void *) mybeg, beg, ncols * sizeof(int));
   mybeg[ncols] = nnonz;   // add additional entry at end

   // load problem
   clp->loadProblem(ncols, nrows, mybeg, ind, val, lb, ub, obj, lhs, rhs);
   delete [] mybeg;

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
RETCODE SCIPlpiAddCols(
   LPI*             lpi,                /**< LP interface structure */
   int              ncols,              /**< number of columns to be added */
   const Real*      obj,                /**< objective function values of new columns */
   const Real*      lb,                 /**< lower bounds of new columns */
   const Real*      ub,                 /**< upper bounds of new columns */
   char**           colnames,           /**< column names, or 0 */
   int              nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*       beg,                /**< start index of each column in ind- and val-array, or 0 if nnonz == 0 */
   const int*       ind,                /**< row indices of constraint matrix entries, or 0 if nnonz == 0 */
   const Real*      val                 /**< values of constraint matrix entries, or 0 if nnonz == 0 */
   )
{
   debugMessage("calling SCIPlpiAddCols()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(obj != 0);
   assert(lb != 0);
   assert(ub != 0);
   assert(nnonz == 0 || beg != 0);
   assert(nnonz == 0 || ind != 0);
   assert(nnonz == 0 || val != 0);

   // store number of columns for later
   int numCols = lpi->clp->getNumCols(); 

   // copy beg-array
   int* mybeg = new int [ncols+1];
   memcpy((void *) mybeg, beg, ncols * sizeof(int));
   mybeg[ncols] = nnonz;   // add additional entry at end

   // add columns
   lpi->clp->addColumns(ncols, lb, ub, obj, mybeg, ind, val);
   delete [] mybeg;

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
RETCODE SCIPlpiDelCols(
   LPI*             lpi,                /**< LP interface structure */
   int              firstcol,           /**< first column to be deleted */
   int              lastcol             /**< last column to be deleted */
   )
{
   debugMessage("calling SCIPlpiDelCols()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->clp->numberColumns());

   // Clp can't delete a range of columns; we have to use deleteColumns (see SCIPlpiDelColset)
   int num = lastcol-firstcol+1;
   int* which = new int [num];

   for (int j = firstcol; j <= lastcol; ++j)
      which[j - firstcol] = j;

   lpi->clp->deleteColumns(num, which);
   delete [] which;

   return SCIP_OKAY;   
}


/** deletes columns from LP; the new position of a column must not be greater that its old position */
RETCODE SCIPlpiDelColset(
   LPI*             lpi,                /**< LP interface structure */
   int*             dstat               /**< deletion status of columns
                                         *   input:  1 if column should be deleted, 0 if not
                                         *   output: new position of column, -1 if column was deleted */
   )
{
   debugMessage("calling SCIPlpiDelColset()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(dstat != 0);

   // transform dstat information
   int ncols = lpi->clp->getNumCols(); 
   int* which = new int [ncols];  
   int cnt = 0;
   for (int j = 0; j < ncols; ++j)
   {
      if ( dstat[j] == 1 )
	 which[cnt++] = j;
   }
   lpi->clp->deleteColumns(cnt, which);
   delete [] which;

   // update dstat (is this necessary?)
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
RETCODE SCIPlpiAddRows(
   LPI*             lpi,                /**< LP interface structure */
   int              nrows,              /**< number of rows to be added */
   const Real*      lhs,                /**< left hand sides of new rows */
   const Real*      rhs,                /**< right hand sides of new rows */
   char**           rownames,           /**< row names, or 0 */
   int              nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*       beg,                /**< start index of each row in ind- and val-array, or 0 if nnonz == 0 */
   const int*       ind,                /**< column indices of constraint matrix entries, or 0 if nnonz == 0 */
   const Real*      val                 /**< values of constraint matrix entries, or 0 if nnonz == 0 */
   )
{
   debugMessage("calling SCIPlpiAddRows()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(lhs != 0);
   assert(rhs != 0);
   assert(nnonz == 0 || beg != 0);
   assert(nnonz == 0 || ind != 0);
   assert(nnonz == 0 || val != 0);

   // store number of rows for later use
   int numRows = lpi->clp->getNumRows(); 

   if ( nnonz != 0 )
   {
      // copy beg-array
      int* mybeg = new int [nrows+1];
      memcpy((void *) mybeg, beg, nrows * sizeof(int));
      mybeg[nrows] = nnonz;   // add additional entry at end
      
      // add rows
      lpi->clp->addRows(nrows, lhs, rhs, mybeg, ind, val);
      delete [] mybeg;
   }
   else
   {
      // add empty rows
      int* mybeg = new int [nrows+1];
      for (int i = 0; i <= nrows; ++i)
	 mybeg[i] = 0;
      lpi->clp->addRows(nrows, lhs, rhs, mybeg, 0, 0);
   }

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
RETCODE SCIPlpiDelRows(
   LPI*             lpi,                /**< LP interface structure */
   int              firstrow,           /**< first row to be deleted */
   int              lastrow             /**< last row to be deleted */
   )
{
   debugMessage("calling SCIPlpiDelRows()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->clp->numberRows());

   // Clp can't delete a range of rows; we have to use deleteRows (see SCIPlpiDelRowset)
   int num = lastrow-firstrow+1;
   int* which = new int [num];

   for (int i = firstrow; i <= lastrow; ++i)
      which[i - firstrow] = i;

   lpi->clp->deleteRows(num, which);
   delete [] which;

   return SCIP_OKAY;   
}


/** deletes rows from LP; the new position of a row must not be greater that its old position */
RETCODE SCIPlpiDelRowset(
   LPI*             lpi,                /**< LP interface structure */
   int*             dstat               /**< deletion status of rows
                                         *   input:  1 if row should be deleted, 0 if not
                                         *   output: new position of row, -1 if row was deleted */
   )
{
   debugMessage("calling SCIPlpiDelRowset()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(dstat != 0);

   // transform dstat information
   int nrows = lpi->clp->getNumRows(); 
   int* which = new int [nrows];
   int cnt = 0;
   for (int i = 0; i < nrows; ++i)
   {
      if ( dstat[i] == 1 )
	 which[cnt++] = i;
   }
   lpi->clp->deleteRows(cnt, which);
   delete [] which;

   // update dstat (is this necessary)
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
RETCODE SCIPlpiClear(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiClear()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   // We use the resize(0,0) to get rid of the model but keep all other settings
   lpi->clp->resize(0,0);

   return SCIP_OKAY;
}


/** changes lower and upper bounds of columns */
RETCODE SCIPlpiChgBounds(
   LPI*             lpi,                /**< LP interface structure */
   int              ncols,              /**< number of columns to change bounds for */
   const int*       ind,                /**< column indices */
   const Real*      lb,                 /**< values for the new lower bounds */
   const Real*      ub                  /**< values for the new upper bounds */
   )
{
   debugMessage("calling SCIPlpiChgBounds()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(ind != 0);
   assert(lb != 0);
   assert(ub != 0);

   ClpSimplex* clp = lpi->clp;

   /* OLD code
   double* colLower = clp->columnLower();     // we have direct access to the data of Clp!
   double* colUpper = clp->columnUpper();
   for( int j = 0; j < ncols; ++j )
   {
      assert(0 <= ind[j] && ind[j] < clp->numberColumns());
      colLower[ind[j]] = lb[j];
      colUpper[ind[j]] = ub[j];
   }
   */

   // new version (updates whatsChanged in Clp)  (bound checking in Clp)
   for( int j = 0; j < ncols; ++j )
      clp->setColBounds(ind[j], lb[j], ub[j]);
   
   return SCIP_OKAY;
}


/** changes left and right hand sides of rows */
RETCODE SCIPlpiChgSides(
   LPI*             lpi,                /**< LP interface structure */
   int              nrows,              /**< number of rows to change sides for */
   const int*       ind,                /**< row indices */
   const Real*      lhs,                /**< new values for left hand sides */
   const Real*      rhs                 /**< new values for right hand sides */
   )
{
   debugMessage("calling SCIPlpiChgSides()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(ind != 0);
   assert(lhs != 0);
   assert(rhs != 0);

   ClpSimplex* clp = lpi->clp;

   /* OLD code
   double* rowLower = clp->rowLower();     // we have direct access to the data of Clp!
   double* rowUpper = clp->rowUpper();
   for( int i = 0; i < nrows; ++i )
   {
      assert(0 <= ind[i] && ind[i] < clp->numberRows());
      rowLower[ind[i]] = lhs[i];
      rowUpper[ind[i]] = rhs[i];
   }
   */

   // new version (updates whatsChanged in Clp)  (bound checking in Clp)
   for( int i = 0; i < nrows; ++i )
      clp->setRowBounds(ind[i], lhs[i], rhs[i]);

   return SCIP_OKAY;
}


/** changes a single coefficient */
RETCODE SCIPlpiChgCoef(
   LPI*             lpi,                /**< LP interface structure */
   int              row,                /**< row number of coefficient to change */
   int              col,                /**< column number of coefficient to change */
   Real             newval              /**< new value of coefficient */
   )
{
   debugMessage("calling SCIPlpiChgCoef()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= row && row < lpi->clp->numberRows());
   assert(0 <= col && col < lpi->clp->numberColumns());

   lpi->clp->matrix()->modifyCoefficient(row, col, newval);

   return SCIP_OKAY;
}


/** changes the objective sense */
RETCODE SCIPlpiChgObjsen(
   LPI*             lpi,                /**< LP interface structure */
   OBJSEN           objsen              /**< new objective sense */
   )
{
   debugMessage("calling SCIPlpiChgObjsen()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   // set objective sense: SCIP values are the same as the ones for Clp
   lpi->clp->setOptimizationDirection(objsen);

   return SCIP_OKAY;
}


/** changes objective values of columns in the LP */
RETCODE SCIPlpiChgObj(
   LPI*             lpi,                /**< LP interface structure */
   int              ncols,              /**< number of columns to change objective value for */
   int*             ind,                /**< column indices to change objective value for */
   Real*            obj                 /**< new objective values for columns */
   )
{
   debugMessage("calling SCIPlpiChgObj()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(ind != 0);
   assert(obj != 0);

   ClpSimplex* clp = lpi->clp;

   /* OLD code
   double* objvec = clp->objective();          // we have direct access to the data of Clp!
   for( int j = 0; j < ncols; ++j )
   {
      assert(0 <= ind[j] && ind[j] < lpi->clp->numberColumns());
      objvec[ind[j]] = obj[j];
   }
   */

   // new version (updates whatsChanged in Clp)  (bound checking in Clp)
   for( int j = 0; j < ncols; ++j )
      clp->setObjCoeff(ind[j], obj[j]);

   return SCIP_OKAY;
}


/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
RETCODE SCIPlpiScaleRow(
   LPI*             lpi,                /**< LP interface structure */
   int              row,                /**< row number to scale */
   Real             scaleval            /**< scaling multiplier */
   )
{
   debugMessage("calling SCIPlpiScaleRow()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(scaleval != 0.0);
   assert(0 <= row && row <= lpi->clp->numberRows() );

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
      Real oldlhs = lhsval;
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
RETCODE SCIPlpiScaleCol(
   LPI*             lpi,                /**< LP interface structure */
   int              col,                /**< column number to scale */
   Real             scaleval            /**< scaling multiplier */
   )
{
   debugMessage("calling SCIPlpiScaleCol()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(scaleval != 0.0);
   assert(0 <= col && col <= lpi->clp->numberColumns() );

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
      Real oldlb = lbval;
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
RETCODE SCIPlpiGetNRows(
   LPI*             lpi,                /**< LP interface structure */
   int*             nrows               /**< pointer to store the number of rows */
   )
{
   debugMessage("calling SCIPlpiGetNRows()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(nrows != 0);

   *nrows = lpi->clp->numberRows();

   return SCIP_OKAY;
}


/** gets the number of columns in the LP */
RETCODE SCIPlpiGetNCols(
   LPI*             lpi,                /**< LP interface structure */
   int*             ncols               /**< pointer to store the number of cols */
   )
{
   debugMessage("calling SCIPlpiGetNCols()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(ncols != 0);

   *ncols = lpi->clp->numberColumns();

   return SCIP_OKAY;
}


/** gets the number of nonzero elements in the LP constraint matrix */
RETCODE SCIPlpiGetNNonz(
   LPI*             lpi,                /**< LP interface structure */
   int*             nnonz               /**< pointer to store the number of nonzeros */
   )
{
   debugMessage("calling SCIPlpiGetNNonz()\n");

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
RETCODE SCIPlpiGetCols(
   LPI*             lpi,                /**< LP interface structure */
   int              firstcol,           /**< first column to get from LP */
   int              lastcol,            /**< last column to get from LP */
   Real*            lb,                 /**< buffer to store the lower bound vector, or 0 */
   Real*            ub,                 /**< buffer to store the upper bound vector, or 0 */
   int*             nnonz,              /**< pointer to store the number of nonzero elements returned, or 0 */
   int*             beg,                /**< buffer to store start index of each column in ind- and val-array, or 0 */
   int*             ind,                /**< buffer to store column indices of constraint matrix entries, or 0 */
   Real*            val                 /**< buffer to store values of constraint matrix entries, or 0 */
   )
{
   debugMessage("calling SCIPlpiGetCols()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->clp->numberColumns());

   ClpSimplex* clp = lpi->clp;

   // get lower and upper bounds for the variables
   if( lb != 0 )
   {
      assert(ub != 0);

      const double* colLower = clp->getColLower();    // Here we can use the const versions (see SCIPchgBounds)
      const double* colUpper = clp->getColUpper();
      for( int j = firstcol; j <= lastcol; ++j )
      {
         lb[j-firstcol] = colLower[j];
         ub[j-firstcol] = colUpper[j];
      }
   }
   else
      assert(ub == 0);

   if( nnonz != 0 )
   {
      CoinPackedMatrix* M = clp->matrix();
      assert( M != 0 );
      assert( M->getNumCols() == clp->numberColumns() );

      const CoinBigIndex* Mbeg = M->getVectorStarts();   // can use const versions
      const int* Mlength = M->getVectorLengths();
      const int* Mind = M->getIndices();
      const double* Mval = M->getElements();

      *nnonz = 0;
      for( int j = firstcol; j <= lastcol; ++j )
      {
         beg[j-firstcol] = *nnonz;
         for( CoinBigIndex k = Mbeg[j]; k < Mbeg[j] + Mlength[j]; ++k )
         {
            ind[*nnonz] = Mind[k];
            val[*nnonz] = Mval[k];
            (*nnonz)++;
         }
      }
   }
   else
   {
      assert(beg == 0);
      assert(ind == 0);
      assert(val == 0);
   }

   return SCIP_OKAY;
}


/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be 0, or both have to be non-0,
 *  either nnonz, beg, ind, and val have to be 0, or all of them have to be non-0.
 */
RETCODE SCIPlpiGetRows(
   LPI*             lpi,                /**< LP interface structure */
   int              firstrow,           /**< first row to get from LP */
   int              lastrow,            /**< last row to get from LP */
   Real*            lhs,                /**< buffer to store left hand side vector, or 0 */
   Real*            rhs,                /**< buffer to store right hand side vector, or 0 */
   int*             nnonz,              /**< pointer to store the number of nonzero elements returned, or 0 */
   int*             beg,                /**< buffer to store start index of each row in ind- and val-array, or 0 */
   int*             ind,                /**< buffer to store row indices of constraint matrix entries, or 0 */
   Real*            val                 /**< buffer to store values of constraint matrix entries, or 0 */
   )
{
   debugMessage("calling SCIPlpiGetRows()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->clp->numberRows());

   ClpSimplex* clp = lpi->clp;
   if( lhs != 0 )
   {
      assert(rhs != 0);

      const double* rowLower = clp->getRowLower();    // Here we can use the const versions (see SCIPchgSides)
      const double* rowUpper = clp->getRowUpper();
      for( int i = firstrow; i <= lastrow; ++i )
      {
         lhs[i-firstrow] = rowLower[i];
         rhs[i-firstrow] = rowUpper[i];
      }
   }
   else
      assert(rhs == 0);

   if( nnonz != 0 )
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
   else
   {
      assert(beg == 0);
      assert(ind == 0);
      assert(val == 0);
   }

   return SCIP_OKAY;
}


/** gets objective coefficients from LP problem object */
RETCODE SCIPlpiGetObj(
   LPI*             lpi,                /**< LP interface structure */
   int              firstcol,           /**< first column to get objective coefficient for */
   int              lastcol,            /**< last column to get objective coefficient for */
   Real*            vals                /**< array to store objective coefficients */
   )
{
   debugMessage("calling SCIPlpiGetObj()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->clp->numberColumns());
   assert(vals != 0);

   const double* obj = lpi->clp->getObjCoefficients();    // Here we can use the const versions (see SCIPchgObj)
   for( int j = firstcol; j <= lastcol; ++j )
      vals[j-firstcol] = obj[j];

   return SCIP_OKAY;
}


/** gets current bounds from LP problem object */
RETCODE SCIPlpiGetBounds(
   LPI*             lpi,                /**< LP interface structure */
   int              firstcol,           /**< first column to get objective value for */
   int              lastcol,            /**< last column to get objective value for */
   Real*            lbs,                /**< array to store lower bound values, or 0 */
   Real*            ubs                 /**< array to store upper bound values, or 0 */
   )
{
   debugMessage("calling SCIPlpiGetBounds()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->clp->numberColumns());
   
   const double* colLower = lpi->clp->getColLower();    // Here we can use the const versions (see SCIPchgBounds)
   const double* colUpper = lpi->clp->getColUpper();
   for( int j = firstcol; j <= lastcol; ++j )
   {
      if (lbs != 0)
	 lbs[j-firstcol] = colLower[j];
      if (ubs != 0)
	 ubs[j-firstcol] = colUpper[j];
   }

   return SCIP_OKAY;
}


/** gets current row sides from LP problem object */
RETCODE SCIPlpiGetSides(
   LPI*             lpi,                /**< LP interface structure */
   int              firstrow,           /**< first row to get sides for */
   int              lastrow,            /**< last row to get sides for */
   Real*            lhss,               /**< array to store left hand side values, or 0 */
   Real*            rhss                /**< array to store right hand side values, or 0 */
   )
{
   // Why are the left hand sides called lhss and not lhs like everywhere ????????
   debugMessage("calling SCIPlpiGetSides()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->clp->numberRows());

   const double* rowLower = lpi->clp->getRowLower();    // Here we can use the const versions (see SCIPchgSides)
   const double* rowUpper = lpi->clp->getRowUpper();
   for( int i = firstrow; i <= lastrow; ++i )
   {
      if (lhss != 0)
	 lhss[i-firstrow] = rowLower[i];
      if (rhss != 0)
	 rhss[i-firstrow] = rowLower[i];
   }

   return SCIP_OKAY;
}


/** gets a single coefficient */
RETCODE SCIPlpiGetCoef(
   LPI*             lpi,                /**< LP interface structure */
   int              row,                /**< row number of coefficient */
   int              col,                /**< column number of coefficient */
   Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   debugMessage("calling SCIPlpiGetCoef()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(0 <= col && col < lpi->clp->numberColumns());
   assert(0 <= row && row < lpi->clp->numberRows());
   assert(val != 0);

   *val = 0.0;
   // There seems to be no direct way to get the coefficient
   CoinPackedMatrix* M = lpi->clp->matrix();
   assert( M != 0 );
   assert( M->getNumCols() == lpi->clp->numberColumns() );

   const CoinBigIndex* Mbeg = M->getVectorStarts();
   const int* Mlength = M->getVectorLengths();
   const int* Mind = M->getIndices();
   const double* Mval = M->getElements();
   for( CoinBigIndex k = Mbeg[col]; k < Mbeg[col] + Mlength[col]; ++k )
   {
      if ( Mind[k] == row )
      {
	 *val = Mval[k];
	 break;
      }
   }

   return SCIP_OKAY;
}

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */


/** calls primal simplex to solve the LP */
RETCODE SCIPlpiSolvePrimal(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling Clp primal(): %d cols, %d rows\n", lpi->clp->numberColumns(), lpi->clp->numberRows());

   assert(lpi != 0);
   assert(lpi->clp != 0);

   // Call the primal Simplex without preprocessing   
   if ( lpi->startscratch )
   {
      lpi->clp->allSlackBasis(true);   // reset basis
      lpi->validFactorization = false;
      //lpi->clp->factorize();           // ?????????????????????
   }

   // temorarily turn off scaling ???????????????????
   int scaling = lpi->clp->scalingFlag();
   lpi->clp->scaling(0);
   int status = lpi->clp->primal(0, lpi->validFactorization ? 3 : 1);
   lpi->clp->scaling(scaling);
   lpi->validFactorization = true;

   // We may need to call ClpModel::status()

   // Unfortunately the status of Clp is hard coded ...
   // 0 - optimal
   // 1 - primal infeasible
   // 2 - dual infeasible
   // 3 - stopped on iterations or time
   // 4 - stopped due to errors
   // 5 - stopped by event handler

   if ( status == 4 || status == 5 )    // begin stopped by event handler should not occur
      return SCIP_LPERROR;

   return SCIP_OKAY;
}


/** calls dual simplex to solve the LP */
RETCODE SCIPlpiSolveDual(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling Clp dual(): %d cols, %d rows\n", lpi->clp->numberColumns(), lpi->clp->numberRows());

   assert(lpi != 0);
   assert(lpi->clp != 0);

   // Call the dual Simplex without preprocessing
   if ( lpi->startscratch )
   {
      lpi->clp->allSlackBasis(true);   // reset basis
      lpi->validFactorization = false;
      //lpi->clp->factorize();  // ???????????????????????
   }

   // temorarily turn off scaling ????????????????????????
   int scaling = lpi->clp->scalingFlag();
   lpi->clp->scaling(0);
   int status = lpi->clp->dual(0, lpi->validFactorization ? 3 : 1);
   lpi->clp->scaling(scaling);
   lpi->validFactorization = true;

   //int status = lpi->clp->dual();

   //std::cout << "Clp Status: " << status << " ( " << lpi->clp->status();
   //std::cout << ")   Secondary: " << lpi->clp->secondaryStatus() << std::endl;

   // We may need to call ClpModel::status()

   // Unfortunately the status of Clp is hard coded ...
   // 0 - optimal
   // 1 - primal infeasible
   // 2 - dual infeasible
   // 3 - stopped on iterations or time
   // 4 - stopped due to errors
   // 5 - stopped by event handler

   if ( status == 4 || status == 5 )    // begin stopped by event handler should not occur
      return SCIP_LPERROR;

   return SCIP_OKAY;
}


/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
RETCODE SCIPlpiSolveBarrier(
   LPI*             lpi,                 /**< LP interface structure */
   Bool             crossover            /**< perform crossover */
   )
{
   debugMessage("calling Clp barrier(): %d cols, %d rows\n", lpi->clp->numberColumns(), lpi->clp->numberRows());

   assert(lpi != 0);
   assert(lpi->clp != 0);

   // Check whether we have a factorization, if yes destroy it (Clp doesn't like it ...)
   /*
   if (lpi->haveFactorization)
      lpi->clp->finish();
   */

   // Call barrier
   int status = lpi->clp->barrier(crossover);

   // We may need to call ClpModel::status()

   // Unfortunately the status of Clp is hard coded ...
   // 0 - optimal
   // 1 - primal infeasible
   // 2 - dual infeasible
   // 3 - stopped on iterations or time
   // 4 - stopped due to errors
   // 5 - stopped by event handler

   if ( status == 4 || status == 5 )    // begin stopped by event handler should not occur
      return SCIP_LPERROR;

   return SCIP_OKAY;
}


/** performs strong branching iterations on all candidates */
RETCODE SCIPlpiStrongbranch(
   LPI*             lpi,                /**< LP interface structure */
   int              col,                /**< column to apply strong branching on */
   Real             psol,               /**< current primal solution value of column */
   int              itlim,              /**< iteration limit for strong branchings */
   Real*            down,               /**< stores dual bound after branching column down */
   Real*            up,                 /**< stores dual bound after branching column up */
   int*             iter                /**< stores total number of strong branching iterations, or -1; may be 0 */
   )
{
   debugMessage("calling SCIPlpiStrongbranch() on variable %d (%d iterations)\n", col, itlim);

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(down != 0);
   assert(up != 0);

   ClpSimplex* clp = lpi->clp;

   // prepare call
   *down = EPSCEIL(psol - 1.0, 1e-06);
   *up   = EPSFLOOR(psol + 1.0, 1e-06);

   int ncols = clp->numberColumns();
   double** outputSolution = new double* [2];
   outputSolution[0] = new double [ncols];
   outputSolution[1] = new double [ncols];
   int* outputStatus = new int [2];
   int* outputIterations = new int [2];

   // set iteration limit
   int iterlimit = clp->maximumIterations();
   clp->setMaximumIterations(itlim);

   double objval = clp->objectiveValue();

   // the last three parameters are:
   // bool stopOnFirstInfeasible=true,
   // bool alwaysFinish=false);
   // int startFinishOptions
   int retval = ((ClpSimplexDual*) clp)->strongBranching(1, &col, up, down, outputSolution, outputStatus, outputIterations, false, true, 1);

   *down += objval;
   *up += objval;

   if (iter)
      *iter = outputIterations[0] + outputIterations[1];

   // reset iteration limit
   clp->setMaximumIterations(iterlimit);

   delete [] outputStatus;
   delete [] outputIterations;
   delete [] outputSolution[1];
   delete [] outputSolution[0];
   delete [] outputSolution;

   lpi->validFactorization = false; // ???????????

   return SCIP_OKAY;
}





/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** gets information about primal and dual feasibility of the LP basis */
RETCODE SCIPlpiGetBasisFeasibility(
   LPI*             lpi,                /**< LP interface structure */
   Bool*            primalfeasible,     /**< stores primal feasibility status */
   Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   debugMessage("calling SCIPlpiGetBasisFeasibility()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(primalfeasible != 0);
   assert(dualfeasible != 0);

   *primalfeasible = lpi->clp->primalFeasible();
   *dualfeasible   = lpi->clp->dualFeasible();

   return SCIP_OKAY;
}


/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point) */
Bool SCIPlpiHasPrimalRay(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiHasPrimalRay()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   // Clp seems to have a primal ray whenever it concludes "dual infeasible", 
   // (but is not necessarily primal feasible), see ClpModel::unboundedRay
   return ( lpi->clp->status() == 2);      // status = 2  is "dual infeasible"
}


/** returns TRUE iff LP is proven to be primal unbounded */
Bool SCIPlpiIsPrimalUnbounded(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsPrimalUnbounded()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return ( lpi->clp->isProvenDualInfeasible() && lpi->clp->primalFeasible() );
}


/** returns TRUE iff LP is proven to be primal infeasible */
Bool SCIPlpiIsPrimalInfeasible(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsPrimalInfeasible()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return ( lpi->clp->isProvenPrimalInfeasible() );
}


/** returns TRUE iff LP is proven to be primal feasible */
Bool SCIPlpiIsPrimalFeasible(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsPrimalFeasible()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return (! lpi->clp->primalFeasible() );
}


/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point) */
extern 
Bool SCIPlpiHasDualRay(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiHasDualRay()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   // Clp assumes to have a dual ray whenever it concludes "primal infeasible", 
   // (but is not necessarily dual feasible), see ClpModel::infeasibilityRay

   return ( lpi->clp->status() == 1 && ! lpi->clp->secondaryStatus() ); // status = 1  is "primal infeasible"
}


/** returns TRUE iff LP is proven to be dual unbounded */
Bool SCIPlpiIsDualUnbounded(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsDualUnbounded()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return ( lpi->clp->isProvenPrimalInfeasible() && lpi->clp->dualFeasible() );
}


/** returns TRUE iff LP is proven to be dual infeasible */
Bool SCIPlpiIsDualInfeasible(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsDualInfeasible()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return ( lpi->clp->isProvenDualInfeasible() );
}


/** returns TRUE iff LP is proven to be dual feasible */
Bool SCIPlpiIsDualFeasible(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsDualFeasible()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return (! lpi->clp->dualFeasible() );
}


/** returns TRUE iff LP was solved to optimality */
Bool SCIPlpiIsOptimal(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsOptimal()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return( lpi->clp->isProvenOptimal() && lpi->clp->secondaryStatus() == 0 );
}


/** returns TRUE iff current LP basis is stable */
Bool SCIPlpiIsStable(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsStable()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   // We first if status is ok, i.e., is one of the following:
   //    0 - optimal
   //    1 - primal infeasible
   //    2 - dual infeasible

   // Then we check the secondary status of Clp:
   /** 0 - none
       1 - primal infeasible because dual limit reached OR probably primal infeasible but can't prove it (main status 4)
       2 - scaled problem optimal - unscaled problem has primal infeasibilities
       3 - scaled problem optimal - unscaled problem has dual infeasibilities
       4 - scaled problem optimal - unscaled problem has primal and dual infeasibilities
       100 up - translation of enum from ClpEventHandler
   */
   return( (lpi->clp->status() <= 2) && (! lpi->clp->isAbandoned()) && (lpi->clp->secondaryStatus() <= 1) );
}


/** returns TRUE iff the objective limit was reached */
Bool SCIPlpiIsObjlimExc(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsObjlimExc()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return( lpi->clp->isPrimalObjectiveLimitReached() || lpi->clp->isDualObjectiveLimitReached() );
}


/** returns TRUE iff the iteration limit was reached */
Bool SCIPlpiIsIterlimExc(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsIterlimExc()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return(lpi->clp->isIterationLimitReached() );
}


/** returns TRUE iff the time limit was reached */
Bool SCIPlpiIsTimelimExc(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return(lpi->clp->status() == 3 && (! lpi->clp->isIterationLimitReached()) );
}


/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiGetInternalStatus()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   return lpi->clp->status();
}


/** gets objective value of solution */
RETCODE SCIPlpiGetObjval(
   LPI*             lpi,                /**< LP interface structure */
   Real*            objval              /**< stores the objective value */
   )
{
   debugMessage("calling SCIPlpiGetObjval()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(objval != 0);

   *objval = lpi->clp->objectiveValue();

   return SCIP_OKAY;
}


/** gets primal and dual solution vectors */
RETCODE SCIPlpiGetSol(
   LPI*             lpi,                /**< LP interface structure */
   Real*            objval,             /**< stores the objective value, may be 0 if not needed */
   Real*            primsol,            /**< primal solution vector, may be 0 if not needed */
   Real*            dualsol,            /**< dual solution vector, may be 0 if not needed */
   Real*            activity,           /**< row activity vector, may be 0 if not needed */
   Real*            redcost             /**< reduced cost vector, may be 0 if not needed */
   )
{
   debugMessage("calling SCIPlpiGetSol()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   ClpSimplex* clp = lpi->clp;
   if( objval != 0 )
      *objval = clp->objectiveValue();

   if( primsol != 0 )
   {
      const double* sol = clp->getColSolution();
      memcpy((void *) primsol, sol, clp->numberColumns() * sizeof(double));
   }
   if( dualsol != 0 )
   {
      const double* dsol = clp->getRowPrice();
      memcpy((void *) dualsol, dsol, clp->numberRows() * sizeof(double));
   }
   if( activity != 0 )
   {
      const double* act = clp->getRowActivity();
      memcpy((void *) activity, act, clp->numberRows() * sizeof(double));
   }
   if( redcost != 0 )
   {
      const double* red = clp->getReducedCost();
      memcpy((void *) redcost, red, clp->numberColumns() * sizeof(double));
   }

   return SCIP_OKAY;
}


/** gets primal ray for unbounded LPs */
RETCODE SCIPlpiGetPrimalRay(
   LPI*             lpi,                /**< LP interface structure */
   Real*            ray                 /**< primal ray */
   )
{
   debugMessage("calling SCIPlpiGetPrimalRay()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(ray != 0);

   const double* clpray = lpi->clp->unboundedRay();
   memcpy((void *) ray, clpray, lpi->clp->numberColumns() * sizeof(double));
   delete [] clpray;
   
   return SCIP_LPERROR;
}

/** gets dual farkas proof for infeasibility */
RETCODE SCIPlpiGetDualfarkas(
   LPI*             lpi,                /**< LP interface structure */
   Real*            dualfarkas          /**< dual farkas row multipliers */
   )
{
   debugMessage("calling SCIPlpiGetDualfarkas()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(dualfarkas != 0);

   const double* dualray = lpi->clp->infeasibilityRay();
   memcpy((void *) dualfarkas, dualray, lpi->clp->numberRows() * sizeof(double));
   delete [] dualray;
   
   return SCIP_LPERROR;
}


/** gets the number of LP iterations of the last solve call */
extern
RETCODE SCIPlpiGetIterations(
   LPI*             lpi,                /**< LP interface structure */
   int*             iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert(lpi != 0);
   assert(iterations != 0);

   *iterations = lpi->clp->numberIterations();

   return SCIP_OKAY;
}

/**@} */




/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
RETCODE SCIPlpiGetBase(
   LPI*             lpi,                /**< LP interface structure */
   int*             cstat,              /**< array to store column basis status, or 0 */
   int*             rstat               /**< array to store row basis status, or 0 */
   )
{
   debugMessage("calling SCIPlpiGetBase()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   ClpSimplex* clp = lpi->clp;
   
   // Clp basis status can be:
   //  isFree = 0x00,
   //  basic = 0x01,
   //  atUpperBound = 0x02,
   //  atLowerBound = 0x03,
   //  superBasic = 0x04,
   //  isFixed = 0x05

   // The conversions come from OsiClpSolverInterface

   /*
   int lookup[] = {SCIP_BASESTAT_ZERO, SCIP_BASESTAT_BASIC, SCIP_BASESTAT_UPPER, SCIP_BASESTAT_LOWER, 
		   SCIP_BASESTAT_ZERO, SCIP_BASESTAT_LOWER};   // conversion table

   if( rstat != 0 )
   {
      for( int i = 0; i < clp->numberRows(); ++i )
	 rstat[i] = lookup[clp->getRowStatus(i)];
   }

   if( cstat != 0 )
   {
      for( int j = 0; j < clp->numberColumns(); ++j )
         cstat[j] = lookup[clp->getBasisColStatus(i)];
   }
   */


   // slower but easier to understand (and portable)
   if( rstat != 0 )
   {
      for( int i = 0; i < clp->numberRows(); ++i )
      {
	 switch ( clp->getRowStatus(i) )
	 {
	 case ClpSimplex::isFree: rstat[i] = SCIP_BASESTAT_ZERO; break;
	 case ClpSimplex::basic:  rstat[i] = SCIP_BASESTAT_BASIC; break;
	 case ClpSimplex::atUpperBound: rstat[i] = SCIP_BASESTAT_UPPER; break;
	 case ClpSimplex::atLowerBound:  rstat[i] = SCIP_BASESTAT_LOWER; break;
	 case ClpSimplex::superBasic: rstat[i] = SCIP_BASESTAT_ZERO; break;
	 case ClpSimplex::isFixed:    
	    if (clp->getRowPrice()[i] > 0.0)
	       rstat[i] = SCIP_BASESTAT_LOWER; 
	    else
	       rstat[i] = SCIP_BASESTAT_UPPER; 
	    break;
	 default: errorMessage("invalid basis status\n");  abort();
	 }
      }
   }

   if( cstat != 0 )
   {
      for( int j = 0; j < clp->numberColumns(); ++j )
      {
	 switch ( clp->getColumnStatus(j) )
	 {
	 case ClpSimplex::isFree: cstat[j] = SCIP_BASESTAT_ZERO; break;
	 case ClpSimplex::basic:  cstat[j] = SCIP_BASESTAT_BASIC; break;
	 case ClpSimplex::atUpperBound: cstat[j] = SCIP_BASESTAT_UPPER; break;
	 case ClpSimplex::atLowerBound:  cstat[j] = SCIP_BASESTAT_LOWER; break;
	 case ClpSimplex::superBasic: cstat[j] = SCIP_BASESTAT_ZERO; break;
	 case ClpSimplex::isFixed: 
	    if (clp->getReducedCost()[j] > 0.0)
	       cstat[j] = SCIP_BASESTAT_LOWER; 
	    else
	       cstat[j] = SCIP_BASESTAT_UPPER; 
	    break;
	 default: errorMessage("invalid basis status\n");  abort();
	 }
      }
   }

   return SCIP_OKAY;
}


/** sets current basis status for columns and rows */
RETCODE SCIPlpiSetBase(
   LPI*             lpi,                /**< LP interface structure */
   int*             cstat,              /**< array with column basis status */
   int*             rstat               /**< array with row basis status */
   )
{
   debugMessage("calling SCIPlpiSetBase()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(cstat != 0);
   assert(rstat != 0);

   // Adapted from OsiClpSolverInterface::setBasisStatus

   ClpSimplex* clp = lpi->clp;
   clp->createStatus();
   
   const double* lhs = clp->getRowLower();
   const double* rhs = clp->getRowUpper();
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
      case SCIP_BASESTAT_BASIC: clp->setRowStatus(i, ClpSimplex::basic); break;
      case SCIP_BASESTAT_UPPER: clp->setRowStatus(i, ClpSimplex::atUpperBound); break;
      case SCIP_BASESTAT_LOWER: 
	 if ( EPSGE(rhs[i], lhs[i], 1e-6) )   // if bounds are equal (or rhs > lhs)
	    clp->setRowStatus(i, ClpSimplex::isFixed);
	 else
	    clp->setRowStatus(i, ClpSimplex::atLowerBound); 
	 break;
      default: errorMessage("invalid basis status\n");  abort();
      }
   }

   const double* lb = clp->getColLower();
   const double* ub = clp->getColUpper();
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
      case SCIP_BASESTAT_BASIC: clp->setColumnStatus(j, ClpSimplex::basic); break;
      case SCIP_BASESTAT_UPPER: clp->setColumnStatus(j, ClpSimplex::atUpperBound); break;
      case SCIP_BASESTAT_LOWER: 
	 if ( EPSGE(ub[j], lb[j], 1e-6) )
	    clp->setColumnStatus(j, ClpSimplex::isFixed);
	 else
	    clp->setColumnStatus(j, ClpSimplex::atLowerBound); 
	 break;
      default: errorMessage("invalid basis status\n");  abort();
      }
   }
   
   //clp->factorize();  // ???????????????????????????????????
   //lpi->validFactorization = false;
   clp->setWhatsChanged(clp->whatsChanged() & (~512));   // tell Clp that basis has changed

   return SCIP_OKAY;
}


/** returns the indices of the basic columns and rows */
RETCODE SCIPlpiGetBasisInd(
   LPI*             lpi,                /**< LP interface structure */
   int*             bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   debugMessage("calling SCIPlpiGetBasisInd()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(bind != 0);

   int cnt = 0;
   ClpSimplex* clp = lpi->clp;
   for( int i = 0; i < clp->numberRows(); ++i )
   {
      int status = clp->getRowStatus(i);
      if ( status == ClpSimplex::basic )
	 bind[cnt++] = -1 - i;
   }
   for( int j = 0; j < clp->numberColumns(); ++j )
   {
      int status = clp->getColumnStatus(j);
      if ( status == ClpSimplex::basic )
	 bind[cnt++] = j;
   }
   
   return SCIP_OKAY;
}


/** get dense row of inverse basis matrix B^-1 */
RETCODE SCIPlpiGetBInvRow(
   LPI*             lpi,                /**< LP interface structure */
   int              r,                  /**< row number */
   Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   debugMessage("calling SCIPlpiGetBInvRow()\n");

   assert( lpi != 0 );
   assert( lpi->clp != 0 );
   assert( coef != 0 );
   assert( 0 <= r && r <= lpi->clp->numberRows() );

   ClpSimplex* clp = lpi->clp;
   clp->getBInvRow(r, coef);

   return SCIP_OKAY;
}


/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
RETCODE SCIPlpiGetBInvARow(
   LPI*             lpi,                /**< LP interface structure */
   int              r,                  /**< row number */
   const Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or 0 */
   Real*            val                 /**< vector to return coefficients */
   )
{
   debugMessage("calling SCIPlpiGetBInvARow()\n");

   assert( lpi != 0 );
   assert( lpi->clp != 0 );
   assert( val != 0 );
   assert( 0 <= r && r <= lpi->clp->numberRows() );

   // WARNING:  binvrow is not used at the moment!  ?????????????

   ClpSimplex* clp = lpi->clp;
   clp->getBInvARow(r, val, 0);
  
   return SCIP_OKAY;
}

/**@} */




/*
 * LP State Methods
 */

/**@name LP State Methods */
/**@{ */

/** stores LPi state (like basis information) into lpistate object */
RETCODE SCIPlpiGetState(
   LPI*             lpi,                /**< LP interface structure */
   MEMHDR*          memhdr,             /**< block memory */
   LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   debugMessage("calling SCIPlpiGetState()\n");

   assert(memhdr != 0);
   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(lpistate != 0);

   int ncols = lpi->clp->numberColumns();
   int nrows = lpi->clp->numberRows();
   assert(ncols >= 0);
   assert(nrows >= 0);
   
   /* allocate lpistate data */
   CHECK_OKAY( lpistateCreate(lpistate, memhdr, ncols, nrows) );

   /* allocate enough memory for storing uncompressed basis information */
   CHECK_OKAY( ensureCstatMem(lpi, ncols) );
   CHECK_OKAY( ensureRstatMem(lpi, nrows) );

   /* get unpacked basis information */
   CHECK_OKAY( SCIPlpiGetBase(lpi, lpi->cstat, lpi->rstat) );

   /* pack LPi state data */
   (*lpistate)->ncols = ncols;
   (*lpistate)->nrows = nrows;
   lpistatePack(*lpistate, lpi->cstat, lpi->rstat);

   return SCIP_OKAY;
}


/** loads LPi state (like basis information) into solver */
RETCODE SCIPlpiSetState(
   LPI*             lpi,                /**< LP interface structure */
   MEMHDR*          /*memhdr*/,         /**< block memory */
   LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   )
{
   debugMessage("calling SCIPlpiSetState()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);
   assert(lpistate != 0);
   assert(lpistate->ncols == lpi->clp->numberColumns());
   assert(lpistate->nrows == lpi->clp->numberRows());

   /* allocate enough memory for storing uncompressed basis information */
   CHECK_OKAY( ensureCstatMem(lpi, lpistate->ncols) );
   CHECK_OKAY( ensureRstatMem(lpi, lpistate->nrows) );

   /* unpack LPi state data */
   lpistateUnpack(lpistate, lpi->cstat, lpi->rstat);

   /* load basis information */
   CHECK_OKAY( SCIPlpiSetBase(lpi, lpi->cstat, lpi->rstat) );

   return SCIP_OKAY;
}


/** frees LPi state information */
RETCODE SCIPlpiFreeState(
   LPI*             lpi,                /**< LP interface structure */
   MEMHDR*          memhdr,             /**< block memory */
   LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   debugMessage("calling SCIPlpiFreeState()\n");

   assert(lpi != 0);

   lpistateFree(lpistate, memhdr);

   return SCIP_OKAY;
}


/** reads LP state (like basis information from a file */
RETCODE SCIPlpiReadState(
   LPI*             lpi,            /**< LP interface structure */
   const char*      fname           /**< file name */
   )
{
   debugMessage("calling SCIPlpiReadState()\n");

   // Note: this might not be what is expected ...
   if ( lpi->clp->restoreModel(fname) )
      return SCIP_ERROR;
   return SCIP_OKAY;
}

/** writes LP state (like basis information) to a file */
RETCODE SCIPlpiWriteState(
   LPI*             lpi,            /**< LP interface structure */
   const char*      fname           /**< file name */
   )
{
   debugMessage("calling SCIPlpiWriteState()\n");

   // Note: this might not be what is expected ...
   if ( lpi->clp->saveModel(fname) )
      return SCIP_ERROR;
   return SCIP_OKAY;
}

/**@} */




/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
RETCODE SCIPlpiGetIntpar(
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   int*             ival                /**< buffer to store the parameter value */
   )
{
   debugMessage("calling SCIPlpiGetIntpar()\n");

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
      *ival = lpi->pricing;          // store pricing method in LPI struct
      break;
   case SCIP_LPPAR_LPINFO:
      /** Amount of print out:
	  0 - none
	  1 - just final
	  2 - just factorizations
	  3 - as 2 plus a bit more
	  4 - verbose
	  above that 8,16,32 etc just for selective debug
      */
      *ival = lpi->clp->logLevel() > 0 ? TRUE : FALSE;
      break;
   case SCIP_LPPAR_LPITLIM:
      *ival = lpi->clp->maximumIterations();
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}


/** sets integer parameter of LP */
RETCODE SCIPlpiSetIntpar(
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   int              ival                /**< parameter value */
   )
{
   debugMessage("calling SCIPlpiSetIntpar()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   // Handle pricing separately ...
   if ( type == SCIP_LPPAR_PRICING )
   {
      // primal pricing:
      // 0 is exact devex, 
      // 1 full steepest, 
      // 2 is partial exact devex
      // 3 switches between 0 and 2 depending on factorization
      // 4 starts as partial dantzig/devex but then may switch between 0 and 2.
      //
      // for dual:
      // 0 is uninitialized, 1 full, 2 is partial uninitialized,
      // 3 starts as 2 but may switch to 1.
      lpi->pricing = ival;
      int primalmode = 0;
      int dualmode = 0;
      switch (ival)
      {
      case SCIP_PRICING_AUTO: primalmode = 3; dualmode = 3; break;
      case SCIP_PRICING_FULL: primalmode = 0; dualmode = 1; break;
      case SCIP_PRICING_STEEP: primalmode = 1; dualmode = 0; break;
      case SCIP_PRICING_STEEPQSTART: primalmode = 1; dualmode = 2; break;
      default: errorMessage("unkown pricing parameter %d!\n", ival); abort();
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
      lpi->validFactorization = false; // ???????????????????????
      break;
   case SCIP_LPPAR_PRICING:
      abort();
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      if( ival )
         lpi->clp->setLogLevel(2);
      else 
         lpi->clp->setLogLevel(0);
      break;
   case SCIP_LPPAR_LPITLIM:
      lpi->clp->setMaximumIterations(ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}


/** gets floating point parameter of LP */
RETCODE SCIPlpiGetRealpar(
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   Real*            dval                /**< buffer to store the parameter value */
   )
{
   debugMessage("calling SCIPlpiGetRealpar()\n");

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
   case SCIP_LPPAR_LOBJLIM:
      if ( lpi->clp->optimizationDirection() > 0 )   // if minimization
	 *dval = lpi->clp->dualObjectiveLimit();
      else
	 *dval = lpi->clp->primalObjectiveLimit();
      break;
   case SCIP_LPPAR_UOBJLIM:
      if ( lpi->clp->optimizationDirection() > 0 )   // if minimization
	 *dval = lpi->clp->primalObjectiveLimit();
      else
	 *dval = lpi->clp->dualObjectiveLimit();
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
RETCODE SCIPlpiSetRealpar(
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   Real             dval                /**< parameter value */
   )
{
   debugMessage("calling SCIPlpiSetRealpar()\n");

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
   case SCIP_LPPAR_LOBJLIM:
      if ( lpi->clp->optimizationDirection() > 0 )   // if minimization
	 lpi->clp->setDualObjectiveLimit(dval);
      else
	 lpi->clp->setPrimalObjectiveLimit(dval);
      break;
   case SCIP_LPPAR_UOBJLIM:
      if ( lpi->clp->optimizationDirection() > 0 )   // if minimization
	 lpi->clp->setPrimalObjectiveLimit(dval);
      else
	 lpi->clp->setDualObjectiveLimit(dval);
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
Real SCIPlpiInfinity(
   LPI*             /*lpi*/             /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiInfinity()\n");

   return COIN_DBL_MAX;
}


/** checks if given value is treated as infinity in the LP solver */
Bool SCIPlpiIsInfinity(
   LPI*             /*lpi*/,            /**< LP interface structure */
   Real             val
   )
{
   debugMessage("calling SCIPlpiIsInfinity()\n");

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
Bool fileExists(
   const char*      filename            /**< file name */
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
RETCODE SCIPlpiReadLP(
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   )
{
   debugMessage("calling SCIPlpiReadLP()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   // WARNING: can only read mps files

   if( !fileExists(fname) )
      return SCIP_NOFILE;

   if( lpi->clp->readMps(fname, true, false) )
      return SCIP_READERROR;
   
   return SCIP_OKAY;
}

/** writes LP to a file */
RETCODE SCIPlpiWriteLP(
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   )
{
   debugMessage("calling SCIPlpiWriteLP()\n");

   assert(lpi != 0);
   assert(lpi->clp != 0);

   if ( lpi->clp->writeMps(fname, 0, 2, lpi->clp->optimizationDirection()) )
      return SCIP_ERROR;

   return SCIP_OKAY;
}

/**@} */
