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
#pragma ident "@(#) $Id: lpi_spx.cpp,v 1.34 2004/10/26 07:30:57 bzfpfend Exp $"

/**@file   lpi_spx.cpp
 * @brief  LP interface for SOPLEX 1.2.2 (optimized version)
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* remember the original value of the DEBUG define and undefine it */
#ifdef DEBUG
#define ___DEBUG
#undef DEBUG
#endif

#include "spxsolver.h"
#include "slufactor.h"
#include "spxsteeppr.h"
#include "spxfastrt.h"

/* reset the DEBUG define to its original SCIP value */
#undef DEBUG
#ifdef ___DEBUG
#define DEBUG
#undef ___DEBUG
#endif

#include <cassert>


/********************************************************************/
/*----------------------------- C++ --------------------------------*/
/********************************************************************/

using namespace soplex;


/** SCIP's SoPlex class */
class SPxSCIP : public SPxSolver
{
   SLUFactor        m_slu;              /**< sparse LU factorization */
   SPxSteepPR       m_price;            /**< steepest edge pricer */
   SPxFastRT        m_ratio;            /**< Harris fast ratio tester */
   char*            m_probname;         /**< problem name */
   bool             m_fromscratch;      /**< use old basis indicator */
   Real             m_objLoLimit;       /**< lower objective limit */
   Real             m_objUpLimit;       /**< upper objective limit */
   Status           m_stat;             /**< solving status */

public:
   SPxSCIP(const char* probname) 
      : SPxSolver(LEAVE, COLUMN),
        m_probname(0),
        m_fromscratch(false),
        m_objLoLimit(-soplex::infinity),
        m_objUpLimit(soplex::infinity),
        m_stat(NO_PROBLEM)
   {
      setSolver(&m_slu);
      setTester(&m_ratio);
      setPricer(&m_price);
      /* no starter, no simplifier, no scaler */

      m_slu.setUtype(SLUFactor::ETA);

      if( probname != NULL )
         setProbname(probname);
   }

   virtual ~SPxSCIP()
   {
      if( m_probname != NULL )
         spx_free(m_probname);  /*lint !e1551*/
   }

   bool getFromScratch() const
   {
      return m_fromscratch;
   }

   void setFromScratch(bool fs)
   {
      m_fromscratch = fs;
   }

   void setProbname(const char* probname)
   {
      assert(probname != NULL);
      if( m_probname != NULL )
         spx_free(m_probname);
      spx_alloc(m_probname, strlen(probname) + 1);
      strcpy(m_probname, probname);
   }

   Real getObjLoLimit() const
   {
      return m_objLoLimit;
   }

   void setObjLoLimit(Real limit)
   {
      m_objLoLimit = limit;
   }

   Real getObjUpLimit() const
   {
      return m_objUpLimit;
   }

   void setObjUpLimit(Real limit)
   {
      m_objUpLimit = limit;
   }

   virtual Status solve()
   {
      if( getFromScratch() )
         SPxSolver::reLoad();

      m_stat = SPxSolver::solve();

      assert(rep() == COLUMN);

      if( m_stat == OPTIMAL )
      {
         Real objval = value();

         if( (objval > m_objUpLimit) || (objval < m_objLoLimit) )
            m_stat = ABORT_VALUE;
      }
      return m_stat;
   }

   Status getStatus() const
   {
      return m_stat;
   }

   virtual void clear()
   {
      SPxSolver::clear();

      m_stat = NO_PROBLEM;
   }

   /* the following methods have to be reimplemented to install a workaround for a SOPLEX bug */
   virtual void addCol(const LPCol& col)
   {
      SPxSolver::addCol(col);
      if( matrixIsSetup )
         SPxBasis::loadMatrixVecs(); /* bug workaround */
   }
   virtual void addCol(SPxColId& theid, const LPCol& col)
   {
      SPxSolver::addCol(theid, col);
      if( matrixIsSetup )
         SPxBasis::loadMatrixVecs(); /* bug workaround */
   }
   virtual void addCols(const LPColSet& pset)
   {
      SPxSolver::addCols(pset);
      if( matrixIsSetup )
         SPxBasis::loadMatrixVecs(); /* bug workaround */
   }
   virtual void addCols(SPxColId theid[], const LPColSet& theset)
   {
      SPxSolver::addCols(theid, theset);
      if( matrixIsSetup )
         SPxBasis::loadMatrixVecs(); /* bug workaround */
   }
};




/********************************************************************/
/*-----------------------------  C  --------------------------------*/
/********************************************************************/

extern "C" 
{
#include "lpi.h"
#include "bitencode.h"
#include "message.h"
}


typedef DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET DUALPACKETSIZE
typedef DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET DUALPACKETSIZE



/** LP interface */
struct LPi
{
   SPxSCIP*         spx;                /**< our SPxSolver implementation */
   int*             cstat;              /**< array for storing column basis status */
   int*             rstat;              /**< array for storing row basis status */
   int              cstatsize;          /**< size of cstat array */
   int              rstatsize;          /**< size of rstat array */
};

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
   assert(lpi != NULL);

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
   assert(lpi != NULL);

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
   assert(lpistate != NULL);
   assert(lpistate->packcstat != NULL);
   assert(lpistate->packrstat != NULL);

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
   assert(lpistate != NULL);
   assert(lpistate->packcstat != NULL);
   assert(lpistate->packrstat != NULL);

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
   assert(lpistate != NULL);
   assert(memhdr != NULL);
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
   assert(memhdr != NULL);
   assert(lpistate != NULL);
   assert(*lpistate != NULL);

   freeBlockMemoryArray(memhdr, &(*lpistate)->packcstat, colpacketNum((*lpistate)->ncols));
   freeBlockMemoryArray(memhdr, &(*lpistate)->packrstat, rowpacketNum((*lpistate)->nrows));
   freeBlockMemory(memhdr, lpistate);
}




/*
 * local methods
 */

/** converts SCIP's objective sense into SOPLEX's objective sense */
static
SPxLP::SPxSense spxObjsen(
   OBJSEN           objsen              /**< SCIP's objective sense value */
   )
{
   switch( objsen )
   {
   case SCIP_OBJSEN_MAXIMIZE:
      return SPxLP::MAXIMIZE;
   case SCIP_OBJSEN_MINIMIZE:
      return SPxLP::MINIMIZE;
   default:
      errorMessage("invalid objective sense\n");
      abort();
   }
}




/*
 * LP Interface Methods
 */


/*
 * Miscellaneous Methods
 */

static char spxname[MAXSTRLEN];

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   SPxSCIP spx("tmp");
   int version;

   version = spx.version();
   sprintf(spxname, "SOPLEX %d.%d.%d", version/100, (version % 100)/10, version % 10);
   return spxname;
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
   assert(lpi != NULL);

   /* create SOPLEX object */
   ALLOC_OKAY( allocMemory(lpi) );
   ALLOC_OKAY( (*lpi)->spx = new SPxSCIP(name) );
   (*lpi)->cstat = NULL;
   (*lpi)->rstat = NULL;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;

   /* set objective sense */
   CHECK_OKAY( SCIPlpiChgObjsen(*lpi, objsen) );

   return SCIP_OKAY;
}

/** deletes an LP problem object */
RETCODE SCIPlpiFree(
   LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(*lpi != NULL);
   assert((*lpi)->spx != NULL);

   /* free LP */
   delete (*lpi)->spx;

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
   char**           colnames,           /**< column names, or NULL */
   int              nrows,              /**< number of rows */
   const Real*      lhs,                /**< left hand sides of rows */
   const Real*      rhs,                /**< right hand sides of rows */
   char**           /*rownames*/,       /**< row names, or NULL */
   int              nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*       beg,                /**< start index of each column in ind- and val-array */
   const int*       ind,                /**< row indices of constraint matrix entries */
   const Real*      val                 /**< values of constraint matrix entries */
   )
{
   debugMessage("calling SCIPlpiLoadColLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);

   SPxSCIP* spx = lpi->spx;
   LPRowSet rows(nrows);
   DSVector emptyVector(0);
   int i;

   spx->clear();

   /* set objective sense */
   spx->changeSense(spxObjsen(objsen));

   /* create empty rows with given sides */
   for( i = 0; i < nrows; ++i )
      rows.add(lhs[i], emptyVector, rhs[i]);
   spx->addRows(rows);
   
   /* create column vectors with coefficients and bounds */
   CHECK_OKAY( SCIPlpiAddCols(lpi, ncols, obj, lb, ub, colnames, nnonz, beg, ind, val) );

   return SCIP_OKAY;
}

/** adds columns to the LP */
RETCODE SCIPlpiAddCols(
   LPI*             lpi,                /**< LP interface structure */
   int              ncols,              /**< number of columns to be added */
   const Real*      obj,                /**< objective function values of new columns */
   const Real*      lb,                 /**< lower bounds of new columns */
   const Real*      ub,                 /**< upper bounds of new columns */
   char**           /*colnames*/,       /**< column names, or NULL */
   int              nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*       beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*       ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   debugMessage("calling SCIPlpiAddCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);

   SPxSCIP* spx = lpi->spx;
   LPColSet cols(ncols);
   DSVector colVector(ncols);
   int last;
   int i;
   int j;

   /* create column vectors with coefficients and bounds */
   for( i = 0; i < ncols; ++i )
   {
      colVector.clear();
      if( nnonz > 0 )
      {
         last = (i == ncols-1 ? nnonz : beg[i+1]);
         for( j = beg[i]; j < last; ++j )
            colVector.add(ind[j], val[j]);
      }
      cols.add(obj[i], lb[i], colVector, ub[i]);
   }
   spx->addCols(cols);
 
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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->nCols());

   lpi->spx->removeColRange(firstcol, lastcol);

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
   int ncols;
   int i;

   debugMessage("calling SCIPlpiDelColset()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   ncols = lpi->spx->nCols();

   /* SOPLEX' removeCols() method deletes the columns with dstat[i] < 0, so we have to negate the values */
   for( i = 0; i < ncols; ++i )
      dstat[i] *= -1;

   lpi->spx->removeCols(dstat);

   return SCIP_OKAY;   
}

/** adds rows to the LP */
RETCODE SCIPlpiAddRows(
   LPI*             lpi,                /**< LP interface structure */
   int              nrows,              /**< number of rows to be added */
   const Real*      lhs,                /**< left hand sides of new rows */
   const Real*      rhs,                /**< right hand sides of new rows */
   char**           /*rownames*/,       /**< row names, or NULL */
   int              nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*       beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*       ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   debugMessage("calling SCIPlpiAddRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);

   SPxSCIP* spx = lpi->spx;
   LPRowSet rows(nrows);
   DSVector rowVector;
   int last;
   int i;
   int j;

   /* create row vectors with given sides */
   for( i = 0; i < nrows; ++i )
   {
      rowVector.clear();
      if( nnonz > 0 )
      {
         last = (i == nrows-1 ? nnonz : beg[i+1]);
         for( j = beg[i]; j < last; ++j )
            rowVector.add(ind[j], val[j]);
      }
      rows.add(lhs[i], rowVector, rhs[i]);
   }
   spx->addRows(rows);

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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->nRows());

   lpi->spx->removeRowRange(firstrow, lastrow);

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
   int nrows;
   int i;

   debugMessage("calling SCIPlpiDelRowset()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   nrows = lpi->spx->nRows();

   /* SOPLEX' removeRows() method deletes the rows with dstat[i] < 0, so we have to negate the values */
   for( i = 0; i < nrows; ++i )
      dstat[i] *= -1;

   lpi->spx->removeRows(dstat);

   return SCIP_OKAY;   
}

/** clears the whole LP */
RETCODE SCIPlpiClear(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiClear()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   lpi->spx->clear();

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
   int i;

   debugMessage("calling SCIPlpiChgBounds()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(lb != NULL);
   assert(ub != NULL);

   for( i = 0; i < ncols; ++i )
   {
      assert(0 <= ind[i] && ind[i] < lpi->spx->nCols());
      lpi->spx->changeBounds(ind[i], lb[i], ub[i]);
   }

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
   int i;

   debugMessage("calling SCIPlpiChgSides()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);

   for( i = 0; i < nrows; ++i )
   {
      assert(0 <= ind[i] && ind[i] < lpi->spx->nRows());
      lpi->spx->changeRange(ind[i], lhs[i], rhs[i]);
   }

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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= row && row < lpi->spx->nRows());
   assert(0 <= col && col < lpi->spx->nCols());

   lpi->spx->changeElement(row, col, newval);

   return SCIP_OKAY;
}

/** changes the objective sense */
RETCODE SCIPlpiChgObjsen(
   LPI*             lpi,                /**< LP interface structure */
   OBJSEN           objsen              /**< new objective sense */
   )
{
   debugMessage("calling SCIPlpiChgObjsen()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   lpi->spx->changeSense(spxObjsen(objsen));

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
   int i;

   debugMessage("calling SCIPlpiChgObj()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(obj != NULL);

   for( i = 0; i < ncols; ++i )
   {
      assert(0 <= ind[i] && ind[i] < lpi->spx->nCols());
      lpi->spx->changeObj(ind[i], obj[i]);
   }

   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
RETCODE SCIPlpiScaleRow(
   LPI*             lpi,                /**< LP interface structure */
   int              row,                /**< row number to scale */
   Real             scaleval            /**< scaling multiplier */
   )
{
   Real lhs;
   Real rhs;

   debugMessage("calling SCIPlpiScaleRow()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(scaleval != 0.0);

   /* get the row vector and the row's sides */
   SVector rowvec = lpi->spx->rowVector(row);
   lhs = lpi->spx->lhs(row);
   rhs = lpi->spx->rhs(row);

   /* scale the row vector */
   rowvec *= scaleval;

   /* adjust the sides */
   if( lhs > -soplex::infinity )
      lhs *= scaleval;
   else if( scaleval < 0.0 )
      lhs = soplex::infinity;
   if( rhs < soplex::infinity )
      rhs *= scaleval;
   else if( scaleval < 0.0 )
      rhs = -soplex::infinity;
   if( scaleval < 0.0 )
   {
      Real oldlhs = lhs;
      lhs = rhs;
      rhs = oldlhs;
   }

   /* create the new row */
   LPRow lprow(lhs, rowvec, rhs);
   
   /* change the row in the LP */
   lpi->spx->changeRow(row, lprow);

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
   Real obj;
   Real lb;
   Real ub;

   debugMessage("calling SCIPlpiScaleCol()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(scaleval != 0.0);

   /* get the col vector and the col's bounds and objective value */
   SVector colvec = lpi->spx->colVector(col);
   obj = lpi->spx->obj(col);
   lb = lpi->spx->lower(col);
   ub = lpi->spx->upper(col);

   /* scale the col vector */
   colvec *= scaleval;

   /* scale the objective value */
   obj *= scaleval;

   /* adjust the bounds */
   if( lb > -soplex::infinity )
      lb /= scaleval;
   else if( scaleval < 0.0 )
      lb = soplex::infinity;
   if( ub < soplex::infinity )
      ub /= scaleval;
   else if( scaleval < 0.0 )
      ub = -soplex::infinity;
   if( scaleval < 0.0 )
   {
      Real oldlb = lb;
      lb = ub;
      ub = oldlb;
   }

   /* create the new col (in LPCol's constructor, the upper bound is given first!) */
   LPCol lpcol(obj, colvec, ub, lb);
   
   /* change the col in the LP */
   lpi->spx->changeCol(col, lpcol);

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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(nrows != NULL);

   *nrows = lpi->spx->nRows();

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
RETCODE SCIPlpiGetNCols(
   LPI*             lpi,                /**< LP interface structure */
   int*             ncols               /**< pointer to store the number of cols */
   )
{
   debugMessage("calling SCIPlpiGetNCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ncols != NULL);

   *ncols = lpi->spx->nCols();

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
RETCODE SCIPlpiGetNNonz(
   LPI*             lpi,                /**< LP interface structure */
   int*             nnonz               /**< pointer to store the number of nonzeros */
   )
{
   int i;

   debugMessage("calling SCIPlpiGetNNonz()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(nnonz != NULL);

   /* SOPLEX has no direct method to return the number of nonzeros, so we have to count them manually */
   *nnonz = 0;
   if( lpi->spx->nRows() < lpi->spx->nCols() )
   {
      for( i = 0; i < lpi->spx->nRows(); ++i )
         (*nnonz) += lpi->spx->rowVector(i).size();
   }
   else
   {
      for( i = 0; i < lpi->spx->nCols(); ++i )
         (*nnonz) += lpi->spx->colVector(i).size();
   }

   return SCIP_OKAY;
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
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
   )
{
   int i;
   int j;

   debugMessage("calling SCIPlpiGetCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->nCols());

   if( lb != NULL )
   {
      assert(ub != NULL);

      const Vector& lbvec = lpi->spx->lower();
      const Vector& ubvec = lpi->spx->upper();
      for( i = firstcol; i <= lastcol; ++i )
      {
         lb[i-firstcol] = lbvec[i];
         ub[i-firstcol] = ubvec[i];
      }
   }
   else
      assert(ub == NULL);

   if( nnonz != NULL )
   {
      *nnonz = 0;
      for( i = firstcol; i <= lastcol; ++i )
      {
         beg[i-firstcol] = *nnonz;
         const SVector& cvec = lpi->spx->colVector(i);
         for( j = 0; j < cvec.size(); ++j )
         {
            ind[*nnonz] = cvec.index(j);
            val[*nnonz] = cvec.value(j);
            (*nnonz)++;
         }
      }
   }
   else
   {
      assert(beg == NULL);
      assert(ind == NULL);
      assert(val == NULL);
   }

   return SCIP_OKAY;
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
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
   )
{
   int i;
   int j;

   debugMessage("calling SCIPlpiGetRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->nRows());

   if( lhs != NULL )
   {
      assert(rhs != NULL);

      const Vector& lhsvec = lpi->spx->lhs();
      const Vector& rhsvec = lpi->spx->rhs();
      for( i = firstrow; i <= lastrow; ++i )
      {
         lhs[i-firstrow] = lhsvec[i];
         rhs[i-firstrow] = rhsvec[i];
      }
   }
   else
      assert(rhs == NULL);

   if( nnonz != NULL )
   {
      *nnonz = 0;
      for( i = firstrow; i <= lastrow; ++i )
      {
         beg[i-firstrow] = *nnonz;
         const SVector& rvec = lpi->spx->rowVector(i);
         for( j = 0; j < rvec.size(); ++j )
         {
            ind[*nnonz] = rvec.index(j);
            val[*nnonz] = rvec.value(j);
            (*nnonz)++;
         }
      }
   }
   else
   {
      assert(beg == NULL);
      assert(ind == NULL);
      assert(val == NULL);
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
   int i;

   debugMessage("calling SCIPlpiGetObj()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->nCols());
   assert(vals != NULL);
   
   for( i = firstcol; i <= lastcol; ++i )
      vals[i-firstcol] = lpi->spx->obj(i);

   return SCIP_OKAY;
}

/** gets current bounds from LP problem object */
RETCODE SCIPlpiGetBounds(
   LPI*             lpi,                /**< LP interface structure */
   int              firstcol,           /**< first column to get objective value for */
   int              lastcol,            /**< last column to get objective value for */
   Real*            lbs,                /**< array to store lower bound values, or NULL */
   Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   int i;

   debugMessage("calling SCIPlpiGetBounds()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->nCols());
   
   for( i = firstcol; i <= lastcol; ++i )
   {
      if( lbs != NULL )
         lbs[i-firstcol] = lpi->spx->lower(i);
      if( ubs != NULL )
         ubs[i-firstcol] = lpi->spx->upper(i);
   }

   return SCIP_OKAY;
}

/** gets current row sides from LP problem object */
RETCODE SCIPlpiGetSides(
   LPI*             lpi,                /**< LP interface structure */
   int              firstrow,           /**< first row to get sides for */
   int              lastrow,            /**< last row to get sides for */
   Real*            lhss,               /**< array to store left hand side values, or NULL */
   Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
   int i;

   debugMessage("calling SCIPlpiGetSides()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->nRows());
   
   for( i = firstrow; i <= lastrow; ++i )
   {
      if( lhss != NULL )
         lhss[i-firstrow] = lpi->spx->lhs(i);
      if( rhss != NULL )
         rhss[i-firstrow] = lpi->spx->rhs(i);
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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= col && col < lpi->spx->nCols());
   assert(0 <= row && row < lpi->spx->nRows());
   assert(val != NULL);

   *val = lpi->spx->colVector(col)[row];

   return SCIP_OKAY;
}

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solves LP -- used for both, primal and dual simplex, because SOPLEX doesn't distinct the two cases */
static
RETCODE spxSolve(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SOPLEX solve(): %d cols, %d rows\n", lpi->spx->nCols(), lpi->spx->nRows());

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   SPxSolver::Status status = lpi->spx->solve();
   debugMessage(" -> SOPLEX status: %d, basis status: %d\n", lpi->spx->getStatus(), lpi->spx->basis().status());

   switch( status )
   {
   case SPxSolver::ABORT_TIME:
   case SPxSolver::ABORT_ITER:
   case SPxSolver::ABORT_VALUE:
   case SPxSolver::SINGULAR:
   case SPxSolver::REGULAR:
   case SPxSolver::UNKNOWN:
   case SPxSolver::OPTIMAL:
   case SPxSolver::UNBOUNDED:
   case SPxSolver::INFEASIBLE:
      return SCIP_OKAY;
   case SPxSolver::NO_PROBLEM:
   case SPxSolver::RUNNING:
   case SPxSolver::ERROR:
   default:
      return SCIP_LPERROR;
   }
}

/** calls primal simplex to solve the LP */
RETCODE SCIPlpiSolvePrimal(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiSolvePrimal()\n");

   return spxSolve(lpi);
}

/** calls dual simplex to solve the LP */
RETCODE SCIPlpiSolveDual(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiSolveDual()\n");

   return spxSolve(lpi);
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
RETCODE SCIPlpiSolveBarrier(
   LPI*             lpi,                /**< LP interface structure */
   Bool             crossover            /**< perform crossover */
   )
{
   debugMessage("calling SCIPlpiSolveBarrier()\n");

   return spxSolve(lpi);
}

/** performs strong branching iterations on all candidates */
RETCODE SCIPlpiStrongbranch(
   LPI*             lpi,                /**< LP interface structure */
   int              col,                /**< column to apply strong branching on */
   Real             psol,               /**< current primal solution value of column */
   int              itlim,              /**< iteration limit for strong branchings */
   Real*            down,               /**< stores dual bound after branching column down */
   Real*            up,                 /**< stores dual bound after branching column up */
   int*             iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SPxSCIP* spx;
   SPxSolver::VarStatus* rowstat;
   SPxSolver::VarStatus* colstat;
   SPxSolver::Status status;
   Real oldlb;
   Real oldub;
   Real newlb;
   Real newub;
   bool error;
   int oldItlim;

   debugMessage("calling SCIPlpiStrongbranch() on variable %d (%d iterations)\n", col, itlim);

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(down != NULL);
   assert(up != NULL);

   /**@todo remember, whether the last solve call was strong branching, and save/restore basis only once */
   spx = lpi->spx;
   rowstat = new SPxSolver::VarStatus[spx->nRows()];
   colstat = new SPxSolver::VarStatus[spx->nCols()]; 
   status = SPxSolver::UNKNOWN;                      
   error = false;                                 

   /* get basis and current bounds of column */
   spx->getBasis(rowstat, colstat);    
   oldlb = spx->lower(col);
   oldub = spx->upper(col);

   if( iter != NULL )
      *iter = 0;

   oldItlim = spx->terminationIter();             
   spx->setTerminationIter(itlim);

   /* down branch */
   newub = EPSCEIL(psol-1.0, 1e-06);
   if( newub >= oldlb - 0.5 )
   {
      debugMessage("strong branching down on x%d (%g) with %d iterations\n", col, psol, itlim);

      spx->changeUpper(col, newub);

      status = spx->solve();
      switch( status )
      {
      case SPxSolver::ABORT_TIME:
      case SPxSolver::ABORT_ITER:
      case SPxSolver::OPTIMAL:
         *down = spx->value();
         break;
      case SPxSolver::ABORT_VALUE:
      case SPxSolver::INFEASIBLE:
         *down = spx->terminationValue();
         break;
      default:
         error = true;
         break;
      }
      if( iter != NULL )
         (*iter) += spx->iterations();
      spx->changeUpper(col, oldub);
      spx->setBasis(rowstat, colstat);
   }
   else
      *down = spx->terminationValue();
      
   /* up branch */
   newlb = EPSFLOOR(psol+1.0, 1e-06);
   if( !error && newlb <= oldub + 0.5 )
   {
      debugMessage("strong branching  up  on x%d (%g) with %d iterations\n", col, psol, itlim);
      
      spx->changeLower(col, newlb);
      
      status = spx->solve();
      switch( status )
      {
      case SPxSolver::ABORT_TIME:
      case SPxSolver::ABORT_ITER:
      case SPxSolver::OPTIMAL:
         *up = spx->value();
         break;
      case SPxSolver::ABORT_VALUE:
      case SPxSolver::INFEASIBLE:
         *up = spx->terminationValue();
         break;
      case SPxSolver::UNBOUNDED:
      default:
         error = true;
         break;
      }
      if( iter != NULL )
         (*iter) += spx->iterations();
      spx->changeLower(col, oldlb);
      spx->setBasis(rowstat, colstat);
   }
   else
      *up = spx->terminationValue();
   
   spx->setTerminationIter(oldItlim);
   
   delete [] rowstat;
   delete [] colstat;

   if( error )
   {
      errorMessage("SOPLEX status %d returned SCIPlpiStrongbranch()\n", int(status));
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/**@} */




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
   SPxBasis::SPxStatus basestatus;

   debugMessage("calling SCIPlpiGetBasisFeasibility()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(primalfeasible != NULL);
   assert(dualfeasible != NULL);

   basestatus = lpi->spx->basis().status();
   *primalfeasible = (basestatus == SPxBasis::PRIMAL || basestatus == SPxBasis::OPTIMAL);
   *dualfeasible = (basestatus == SPxBasis::DUAL || basestatus == SPxBasis::OPTIMAL);

   return SCIP_OKAY;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point) */
Bool SCIPlpiHasPrimalRay(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiHasPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::UNBOUNDED);
}

/** returns TRUE iff LP is proven to be primal unbounded */
Bool SCIPlpiIsPrimalUnbounded(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsPrimalUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::UNBOUNDED && lpi->spx->basis().status() == SPxBasis::PRIMAL);
}

/** returns TRUE iff LP is proven to be primal infeasible */
Bool SCIPlpiIsPrimalInfeasible(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsPrimalInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::INFEASIBLE);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point) */
Bool SCIPlpiHasDualRay(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiHasDualRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::INFEASIBLE);
}

/** returns TRUE iff LP is dual unbounded */
Bool SCIPlpiIsDualUnbounded(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsDualUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::INFEASIBLE && lpi->spx->basis().status() == SPxBasis::DUAL);
}

/** returns TRUE iff LP is dual infeasible */
Bool SCIPlpiIsDualInfeasible(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsDualInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::UNBOUNDED);
}

/** returns TRUE iff LP was solved to optimality */
Bool SCIPlpiIsOptimal(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsOptimal()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::OPTIMAL);
}

/** returns TRUE iff current LP basis is stable */
Bool SCIPlpiIsStable(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsStable()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() != SPxSolver::ERROR && lpi->spx->getStatus() != SPxSolver::SINGULAR);
}

/** returns TRUE iff the objective limit was reached */
Bool SCIPlpiIsObjlimExc(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsObjlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::ABORT_VALUE);
}

/** returns TRUE iff the iteration limit was reached */
Bool SCIPlpiIsIterlimExc(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsIterlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::ABORT_ITER);
}

/** returns TRUE iff the time limit was reached */
Bool SCIPlpiIsTimelimExc(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::ABORT_TIME);
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   debugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return lpi->spx->getStatus();
}

/** gets objective value of solution */
RETCODE SCIPlpiGetObjval(
   LPI*             lpi,                /**< LP interface structure */
   Real*            objval              /**< stores the objective value */
   )
{
   debugMessage("calling SCIPlpiGetObjval()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(objval != NULL);

   *objval = lpi->spx->value();

   return SCIP_OKAY;
}

/** gets primal and dual solution vectors */
RETCODE SCIPlpiGetSol(
   LPI*             lpi,                /**< LP interface structure */
   Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
   Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   Real*            activity,           /**< row activity vector, may be NULL if not needed */
   Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{
   debugMessage("calling SCIPlpiGetSol()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   if( objval != NULL )
      *objval = lpi->spx->value();

   if( primsol != NULL )
   {
      Vector tmp(lpi->spx->nCols(), primsol);
      lpi->spx->getPrimal(tmp);
   }
   if( dualsol != NULL )
   {
      Vector tmp(lpi->spx->nRows(), dualsol);
      lpi->spx->getDual(tmp);
   }
   if( activity != NULL )
   {
      Vector tmp(lpi->spx->nRows(), activity);
      lpi->spx->getSlacks(tmp);  /* in SOPLEX, the activities are called "slacks" */
   }
   if( redcost != NULL )
   {
      Vector tmp(lpi->spx->nCols(), redcost);
      lpi->spx->getRedCost(tmp);
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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   errorMessage("SCIPlpiGetPrimalRay() not supported by SOPLEX\n");
   
   return SCIP_LPERROR;
}

/** gets dual farkas proof for infeasibility */
RETCODE SCIPlpiGetDualfarkas(
   LPI*             lpi,                /**< LP interface structure */
   Real*            dualfarkas          /**< dual farkas row multipliers */
   )
{
   debugMessage("calling SCIPlpiGetDualfarkas()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   Vector tmp(lpi->spx->nRows(), dualfarkas);
   lpi->spx->getDualfarkas(tmp);

   return SCIP_OKAY;
}

/** gets the number of LP iterations of the last solve call */
RETCODE SCIPlpiGetIterations(
   LPI*             lpi,                /**< LP interface structure */
   int*             iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   debugMessage("calling SCIPlpiGetIterations()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   *iterations = lpi->spx->iterations();

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
   int*             cstat,              /**< array to store column basis status, or NULL */
   int*             rstat               /**< array to store row basis status, or NULL */
   )
{
   int i;

   debugMessage("calling SCIPlpiGetBase()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   if( rstat != NULL )
   {
      for( i = 0; i < lpi->spx->nRows(); ++i )
      {
         switch( lpi->spx->getBasisRowStatus(i) )
         {
         case SPxSolver::BASIC:
            rstat[i] = SCIP_BASESTAT_BASIC;
            break;	  
         case SPxSolver::FIXED:
         case SPxSolver::ON_LOWER:
            rstat[i] = SCIP_BASESTAT_LOWER;
            break;
         case SPxSolver::ON_UPPER:
            rstat[i] = SCIP_BASESTAT_UPPER;
            break;
         case SPxSolver::ZERO:
            errorMessage("slack variable has basis status ZERO (should not occur)\n");
            return SCIP_LPERROR;
         default:
            errorMessage("invalid basis status\n");
            abort();
         }
      }
   }

   if( cstat != NULL )
   {
      for( i = 0; i < lpi->spx->nCols(); ++i )
      {
         switch( lpi->spx->getBasisColStatus(i) )
         {
         case SPxSolver::BASIC:
            cstat[i] = SCIP_BASESTAT_BASIC;
            break;	  
         case SPxSolver::FIXED:
            assert(lpi->spx->rep() == SPxSolver::COLUMN);
            if( lpi->spx->pVec()[i] - lpi->spx->maxObj()[i] < 0.0 )  /* reduced costs < 0 => UPPER  else => LOWER */
               cstat[i] = SCIP_BASESTAT_UPPER;
            else
               cstat[i] = SCIP_BASESTAT_LOWER;
            break;
         case SPxSolver::ON_LOWER:
            cstat[i] = SCIP_BASESTAT_LOWER;
            break;
         case SPxSolver::ON_UPPER:
            cstat[i] = SCIP_BASESTAT_UPPER;
            break;
         case SPxSolver::ZERO:
            cstat[i] = SCIP_BASESTAT_ZERO;
            break;
         default:
            errorMessage("invalid basis status\n");
            abort();
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
   int i;

   debugMessage("calling SCIPlpiSetBase()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(cstat != NULL);
   assert(rstat != NULL);

   SPxSolver::VarStatus* spxcstat = new SPxSolver::VarStatus[lpi->spx->nCols()];
   SPxSolver::VarStatus* spxrstat = new SPxSolver::VarStatus[lpi->spx->nRows()];

   for( i = 0; i < lpi->spx->nRows(); ++i )
   {
      switch( rstat[i] )
      {
      case SCIP_BASESTAT_LOWER:
         spxrstat[i] = SPxSolver::ON_LOWER;
         break;
      case SCIP_BASESTAT_BASIC:
         spxrstat[i] = SPxSolver::BASIC;
         break;
      case SCIP_BASESTAT_UPPER:
         spxrstat[i] = SPxSolver::ON_UPPER;
         break;
      case SCIP_BASESTAT_ZERO:
         errorMessage("slack variable has basis status ZERO (should not occur)\n");
         return SCIP_LPERROR;
      default:
         errorMessage("invalid basis status\n");
         abort();
      }
   }

   for( i = 0; i < lpi->spx->nCols(); ++i )
   {
      switch( cstat[i] )
      {
      case SCIP_BASESTAT_LOWER:
         spxcstat[i] = SPxSolver::ON_LOWER;
         break;
      case SCIP_BASESTAT_BASIC:
         spxcstat[i] = SPxSolver::BASIC;
         break;
      case SCIP_BASESTAT_UPPER:
         spxcstat[i] = SPxSolver::ON_UPPER;
         break;
      case SCIP_BASESTAT_ZERO:
         spxcstat[i] = SPxSolver::ZERO;
         break;
      default:
         errorMessage("invalid basis status\n");
         abort();
      }
   }
   lpi->spx->setBasis(spxrstat, spxcstat);

   delete[] spxcstat;
   delete[] spxrstat;
   
   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows */
RETCODE SCIPlpiGetBasisInd(
   LPI*             lpi,                /**< LP interface structure */
   int*             bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   int i;

   debugMessage("calling SCIPlpiGetBasisInd()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   for( i = 0; i < lpi->spx->nRows(); ++i )
   {
      SPxId id = lpi->spx->basis().baseId(i);
      if( lpi->spx->isId(id) ) /* column id? */
         bind[i] = lpi->spx->number(id);
      else                     /* row id?    */
         bind[i] = -1 - lpi->spx->number(id);
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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   Vector x(lpi->spx->nRows(), coef);
   DVector e(lpi->spx->nRows());

   /**@todo make this faster by pregenerating all dense unit vectors */
   e = lpi->spx->unitVector(r);

   /* solve system "x = e_r^T * B^-1" to get r'th row of B^-1 */
   lpi->spx->basis().coSolve(x, e);

   return SCIP_OKAY;
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
RETCODE SCIPlpiGetBInvARow(
   LPI*             lpi,                /**< LP interface structure */
   int              r,                  /**< row number */
   const Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   Real*            val                 /**< vector to return coefficients */
   )
{
   Real* buf;
   Real* binv;
   int nrows;
   int ncols;
   int c;

   debugMessage("calling SCIPlpiGetBInvARow()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   nrows = lpi->spx->nRows();
   ncols = lpi->spx->nCols();

   /* get (or calculate) the row in B^-1 */
   if( binvrow == NULL )
   {
      ALLOC_OKAY( allocMemoryArray(&buf, nrows) );
      CHECK_OKAY( SCIPlpiGetBInvRow(lpi, r, buf) );
      binv = buf;
   }
   else
   {
      buf = NULL;
      binv = const_cast<Real*>(binvrow);
   }
   assert(binv != NULL);

   /* calculate the scalar product of the row in B^-1 and A */
   soplex::Vector binvvec(nrows, binv);
   for( c = 0; c < ncols; ++c )
      val[c] = binvvec * lpi->spx->colVector(c);  /* scalar product */

   /* free memory if it was temporarily allocated */
   freeMemoryArrayNull(&buf);

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
   int ncols;
   int nrows;

   debugMessage("calling SCIPlpiGetState()\n");

   assert(memhdr != NULL);
   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpistate != NULL);

   ncols = lpi->spx->nCols();
   nrows = lpi->spx->nRows();
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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpistate != NULL);
   assert(lpistate->ncols == lpi->spx->nCols());
   assert(lpistate->nrows == lpi->spx->nRows());

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

   assert(lpi != NULL);

   lpistateFree(lpistate, memhdr);

   return SCIP_OKAY;
}

/** reads LP state (like basis information from a file */
RETCODE SCIPlpiReadState(
   LPI*             /*lpi*/,            /**< LP interface structure */
   const char*      /*fname*/           /**< file name */
   )
{
   debugMessage("calling SCIPlpiReadState()\n");

   errorMessage("SCIPlpiReadState() not implemented yet in SOPLEX interface\n");

   return SCIP_INVALIDCALL;
}

/** writes LP state (like basis information) to a file */
RETCODE SCIPlpiWriteState(
   LPI*             /*lpi*/,            /**< LP interface structure */
   const char*      /*fname*/           /**< file name */
   )
{
   debugMessage("calling SCIPlpiWriteState()\n");

   errorMessage("SCIPlpiWriteState() not implemented yet in SOPLEX interface\n");

   return SCIP_INVALIDCALL;
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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ival != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      *ival = lpi->spx->getFromScratch();
      break;
   case SCIP_LPPAR_LPINFO:
      *ival = (Param::verbose() > 0 ? TRUE : FALSE);
      break;
   case SCIP_LPPAR_LPITLIM:
      *ival = lpi->spx->terminationIter();
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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      assert(ival == TRUE || ival == FALSE);
      lpi->spx->setFromScratch(bool(ival));
      break;
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      if( ival )
         Param::setVerbose(3);
      else 
         Param::setVerbose(0);
      break;
   case SCIP_LPPAR_LPITLIM:
      lpi->spx->setTerminationIter(ival);
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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(dval != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      *dval = lpi->spx->delta();
      break;
   case SCIP_LPPAR_LOBJLIM:
      *dval = lpi->spx->getObjLoLimit();
      break;
   case SCIP_LPPAR_UOBJLIM:
      *dval = lpi->spx->getObjUpLimit();
      break;
   case SCIP_LPPAR_LPTILIM:
      *dval = lpi->spx->terminationTime();
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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      lpi->spx->setDelta(dval);
      break;
   case SCIP_LPPAR_LOBJLIM:
      lpi->spx->setObjLoLimit(dval);
      break;
   case SCIP_LPPAR_UOBJLIM:
      lpi->spx->setObjUpLimit(dval);
      break;
   case SCIP_LPPAR_LPTILIM:
      lpi->spx->setTerminationTime(dval);
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

   return soplex::infinity;
}

/** checks if given value is treated as infinity in the LP solver */
Bool SCIPlpiIsInfinity(
   LPI*             /*lpi*/,            /**< LP interface structure */
   Real             val
   )
{
   debugMessage("calling SCIPlpiIsInfinity()\n");

   return (val >= soplex::infinity);
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
   if( f == NULL )
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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   if( !fileExists(fname) )
      return SCIP_NOFILE;

   if( !lpi->spx->readFile(fname) )
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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   lpi->spx->dumpFile(fname);

   return SCIP_OKAY;
}

/**@} */

