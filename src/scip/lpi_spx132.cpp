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

/**@file   lpi_spx132.cpp
 * @ingroup LPIS
 * @brief  LP interface for SoPlex 1.3.2
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Marc Pfetsch
 *
 * This is an interface for SoPlex version 1.3.2. It is likely to also work for version 1.3.0 and
 * 1.3.1, but this is untested.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* turn off asserts; otherwise there is an assert(ptr == 0) failure at OneSpin Solutions */
#ifndef NDEBUG
#define NDEBUG
#endif

/* remember the original value of the SCIP_DEBUG define and undefine it */
#ifdef SCIP_DEBUG
#define ___DEBUG
#undef SCIP_DEBUG
#endif

#include "spxsolver.h"
#include "slufactor.h"
#include "spxsteeppr.h"
#include "spxparmultpr.h"
#include "spxdevexpr.h"
#include "spxfastrt.h"

/* check version */
#if ( (SOPLEX_VERSION > 132) || (SOPLEX_VERSION < 130) )
#error "This interface is for SoPlex version 132"
#endif


/* reset the SCIP_DEBUG define to its original SCIP value */
#undef SCIP_DEBUG
#ifdef ___DEBUG
#define SCIP_DEBUG
#undef ___DEBUG
#endif

#include <cassert>

#include "scip/message.h"

/********************************************************************/
/*----------------------------- C++ --------------------------------*/
/********************************************************************/

/* in C++ we have to use "0" instead of "(void*)0" */
#undef NULL
#define NULL 0


#include <iostream>
using namespace soplex;


/** SCIP's SoPlex class */
class SPxSCIP : public SPxSolver
{
   SLUFactor        m_slu;              /**< sparse LU factorization */
   SPxSteepPR       m_price_steep;      /**< steepest edge pricer */
   SPxParMultPR     m_price_parmult;    /**< partial multiple pricer */
   SPxDevexPR       m_price_devex;      /**< devex pricer */
   SPxFastRT        m_ratio;            /**< Harris fast ratio tester */
   char*            m_probname;         /**< problem name */
   bool             m_fromscratch;      /**< use old basis indicator */
   Real             m_objLoLimit;       /**< lower objective limit */
   Real             m_objUpLimit;       /**< upper objective limit */
   Status           m_stat;             /**< solving status */
   bool             m_lpinfo;           /**< storing whether output is turned on */

public:
   SPxSCIP(const char* probname = NULL)
      : SPxSolver(LEAVE, COLUMN),
        m_probname(0),
        m_fromscratch(false),
        m_objLoLimit(-soplex::infinity),
        m_objUpLimit(soplex::infinity),
        m_stat(NO_PROBLEM),
        m_lpinfo(false)
   {
      setSolver(&m_slu);
      setTester(&m_ratio);
      setPricer(&m_price_steep);
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

   void setAutoPricer()
   {
      setPricer(&m_price_steep);
   }

   void setSteepPricer()
   {
      setPricer(&m_price_steep);
   }

   void setParmultPricer()
   {
      setPricer(&m_price_parmult);
   }

   void setDevexPricer()
   {
      setPricer(&m_price_devex);
   }

   bool getFromScratch() const
   {
      return m_fromscratch;
   }

   void setFromScratch(bool fs)
   {
      m_fromscratch = fs;
   }

   bool getLpInfo() const
   {
      return m_lpinfo;
   }

   void setLpInfo(bool li)
   {
      m_lpinfo = li;
   }

   void setProbname(const char* probname)
   {
      int len;

      assert(probname != NULL);
      if( m_probname != NULL )
         spx_free(m_probname);
      len = (int)strlen(probname);
      spx_alloc(m_probname, len + 1);
      strncpy(m_probname, probname, len);
      m_probname[len] = '\0';
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
      /*SPxSolver::setTerminationValue(limit); this is buggy in Soplex!*/
   }

   virtual Status solve()
   {
      if( getFromScratch() )
	 SPxSolver::reLoad();

      if ( getLpInfo() )
	 Param::setVerbose(5);
      else
	 Param::setVerbose(0);

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

   /* the following methods have to be reimplemented to install a workaround for a SoPlex bug */
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
}; /*lint !e1748*/




/********************************************************************/
/*-----------------------------  C  --------------------------------*/
/********************************************************************/

#include "scip/lpi.h"
#include "scip/bitencode.h"

typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE



/** LP interface */
struct SCIP_LPi
{
   SPxSCIP*              spx;                /**< our SPxSolver implementation */
   int*                  cstat;              /**< array for storing column basis status */
   int*                  rstat;              /**< array for storing row basis status */
   int                   cstatsize;          /**< size of cstat array */
   int                   rstatsize;          /**< size of rstat array */
   SCIP_PRICING          pricing;            /**< current pricing strategy */
   SCIP_Bool             solved;             /**< was the current LP solved? */
};

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
   assert(lpi != NULL);

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
   assert(lpi != NULL);

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
   assert(blkmem != NULL);
   assert(lpistate != NULL);
   assert(*lpistate != NULL);

   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packcstat, colpacketNum((*lpistate)->ncols));
   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packrstat, rowpacketNum((*lpistate)->nrows));
   BMSfreeBlockMemory(blkmem, lpistate);
}




/*
 * local methods
 */

/** converts SCIP's objective sense into SoPlex's objective sense */
static
SPxLP::SPxSense spxObjsen(
   SCIP_OBJSEN           objsen              /**< SCIP's objective sense value */
   )
{
   switch( objsen )
   {
   case SCIP_OBJSEN_MAXIMIZE:
      return SPxLP::MAXIMIZE;
   case SCIP_OBJSEN_MINIMIZE:
      return SPxLP::MINIMIZE;
   default:
      SCIPerrorMessage("invalid objective sense\n");
      SCIPABORT();
      return SPxLP::MINIMIZE;
   }
}

/** marks the current LP to be unsolved */
static
void invalidateSolution(SCIP_LPI* lpi)
{
   assert(lpi != NULL);
   lpi->solved = FALSE;
}




/*
 * LP Interface Methods
 */


/*
 * Miscellaneous Methods
 */

static char spxname[100];

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
   sprintf(spxname, "SoPlex %d.%d.%d", version/100, (version % 100)/10, version % 10);
   return spxname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return "Linear Programming Solver developed at Zuse Institute Berlin (soplex.zib.de)";
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   return (void*) lpi->spx;
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
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   assert(lpi != NULL);

   /* create SoPlex object */
   SCIP_ALLOC( BMSallocMemory(lpi) );
   SCIP_ALLOC( (*lpi)->spx = new SPxSCIP(name) ); /*lint !e774*/
   (*lpi)->cstat = NULL;
   (*lpi)->rstat = NULL;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->pricing = SCIP_PRICING_AUTO;
   invalidateSolution(*lpi);

   /* set objective sense */
   SCIP_CALL( SCIPlpiChgObjsen(*lpi, objsen) );

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(*lpi != NULL);
   assert((*lpi)->spx != NULL);

   /* free LP */
   delete (*lpi)->spx;

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
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand sides of rows */
   const SCIP_Real*      rhs,                /**< right hand sides of rows */
   char**                /*rownames*/,       /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array */
   const int*            ind,                /**< row indices of constraint matrix entries */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
   SCIPdebugMessage("calling SCIPlpiLoadColLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);

   invalidateSolution(lpi);

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
   SCIP_CALL( SCIPlpiAddCols(lpi, ncols, obj, lb, ub, colnames, nnonz, beg, ind, val) );

   return SCIP_OKAY;
}

/** adds columns to the LP */
SCIP_RETCODE SCIPlpiAddCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   const SCIP_Real*      obj,                /**< objective function values of new columns */
   const SCIP_Real*      lb,                 /**< lower bounds of new columns */
   const SCIP_Real*      ub,                 /**< upper bounds of new columns */
   char**                /*colnames*/,       /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiAddCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);

   invalidateSolution(lpi);

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
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->nCols());

   invalidateSolution(lpi);

   lpi->spx->removeColRange(firstcol, lastcol);

   return SCIP_OKAY;
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   int ncols;
   int i;

   SCIPdebugMessage("calling SCIPlpiDelColset()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   ncols = lpi->spx->nCols();

   /* SoPlex removeCols() method deletes the columns with dstat[i] < 0, so we have to negate the values */
   for( i = 0; i < ncols; ++i )
      dstat[i] *= -1;

   lpi->spx->removeCols(dstat);

   return SCIP_OKAY;
}

/** adds rows to the LP */
SCIP_RETCODE SCIPlpiAddRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   char**                /*rownames*/,       /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiAddRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);

   invalidateSolution(lpi);

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
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->nRows());

   invalidateSolution(lpi);

   lpi->spx->removeRowRange(firstrow, lastrow);

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
   int nrows;
   int i;

   SCIPdebugMessage("calling SCIPlpiDelRowset()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   nrows = lpi->spx->nRows();

   /* SoPlex removeRows() method deletes the rows with dstat[i] < 0, so we have to negate the values */
   for( i = 0; i < nrows; ++i )
      dstat[i] *= -1;

   lpi->spx->removeRows(dstat);

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiClear()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   lpi->spx->clear();

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
   int i;

   SCIPdebugMessage("calling SCIPlpiChgBounds()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(lb != NULL);
   assert(ub != NULL);

   invalidateSolution(lpi);

   for( i = 0; i < ncols; ++i )
   {
      assert(0 <= ind[i] && ind[i] < lpi->spx->nCols());
      lpi->spx->changeBounds(ind[i], lb[i], ub[i]);
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
   int i;

   SCIPdebugMessage("calling SCIPlpiChgSides()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);

   invalidateSolution(lpi);

   for( i = 0; i < nrows; ++i )
   {
      assert(0 <= ind[i] && ind[i] < lpi->spx->nRows());
      lpi->spx->changeRange(ind[i], lhs[i], rhs[i]);
   }

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
   assert(lpi->spx != NULL);
   assert(0 <= row && row < lpi->spx->nRows());
   assert(0 <= col && col < lpi->spx->nCols());

   invalidateSolution(lpi);

   lpi->spx->changeElement(row, col, newval);

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
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   lpi->spx->changeSense(spxObjsen(objsen));

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
   int i;

   SCIPdebugMessage("calling SCIPlpiChgObj()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(obj != NULL);

   invalidateSolution(lpi);

   for( i = 0; i < ncols; ++i )
   {
      assert(0 <= ind[i] && ind[i] < lpi->spx->nCols());
      lpi->spx->changeObj(ind[i], obj[i]);
   }

   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiScaleRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIP_Real lhs;
   SCIP_Real rhs;

   SCIPdebugMessage("calling SCIPlpiScaleRow()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(scaleval != 0.0);

   invalidateSolution(lpi);

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
      SCIP_Real oldlhs = lhs;
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
SCIP_RETCODE SCIPlpiScaleCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;

   SCIPdebugMessage("calling SCIPlpiScaleCol()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(scaleval != 0.0);

   invalidateSolution(lpi);

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
      SCIP_Real oldlb = lb;
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
SCIP_RETCODE SCIPlpiGetNRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(nrows != NULL);

   *nrows = lpi->spx->nRows();

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
   assert(lpi->spx != NULL);
   assert(ncols != NULL);

   *ncols = lpi->spx->nCols();

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetNNonz()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(nnonz != NULL);

   /* SoPlex has no direct method to return the number of nonzeros, so we have to count them manually */
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
SCIP_RETCODE SCIPlpiGetCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int i;
   int j;

   SCIPdebugMessage("calling SCIPlpiGetCols()\n");

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
SCIP_RETCODE SCIPlpiGetRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_Real*            lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int i;
   int j;

   SCIPdebugMessage("calling SCIPlpiGetRows()\n");

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
   return SCIP_ERROR;
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
   return SCIP_ERROR;
}

/** gets the objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   SCIPerrorMessage("SCIPlpiGetObjsen() has not been implemented yet.\n");
   return SCIP_ERROR;
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiGetObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetObj()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->nCols());
   assert(vals != NULL);

   for( i = firstcol; i <= lastcol; ++i )
      vals[i-firstcol] = lpi->spx->obj(i);

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
   int i;

   SCIPdebugMessage("calling SCIPlpiGetBounds()\n");

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
SCIP_RETCODE SCIPlpiGetSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetSides()\n");

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
SCIP_RETCODE SCIPlpiGetCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetCoef()\n");

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

/** solves LP -- used for both, primal and dual simplex, because SoPlex doesn't distinct the two cases */
static
SCIP_RETCODE spxSolve(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SPxSolver::Type       type                /**< algorithm type */
   )
{
   SCIPdebugMessage("calling SoPlex solve(): %d cols, %d rows\n", lpi->spx->nCols(), lpi->spx->nRows());

   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( type == SPxSolver::ENTER || type == SPxSolver::LEAVE );

   invalidateSolution(lpi);

   /* set the algorithm type */
   lpi->spx->setType(type);

   SPxSolver::Status status = lpi->spx->solve();
   SCIPdebugMessage(" -> SoPlex status: %d, basis status: %d\n", lpi->spx->getStatus(), lpi->spx->basis().status());
   lpi->solved = TRUE;

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
   default:
      return SCIP_LPERROR;
   }  /*lint !e788*/
}

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiSolvePrimal()\n");

   /* SoPlex doesn't distinct between the primal and dual simplex; however
    * we can force SoPlex to start with the desired method:
    * If the representation is COLUMN:
    * - ENTER = PRIMAL
    * - LEAVE = DUAL
    */
   return spxSolve(lpi, SPxSolver::ENTER);
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiSolveDual()\n");

   /* SoPlex doesn't distinct between the primal and dual simplex; however
    * we can force SoPlex to start with the desired method:
    * If the representation is COLUMN:
    * - ENTER = PRIMAL
    * - LEAVE = DUAL
    */
   return spxSolve(lpi, SPxSolver::LEAVE);
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiSolveBarrier()\n");

   /* SoPlex doesn't distinct between the primal and dual simplex; however
    * we can force SoPlex to start with the desired method:
    * If the representation is COLUMN:
    * - ENTER = PRIMAL
    * - LEAVE = DUAL
    *
    * since SoPlex does not support barrier we switch to DUAL
    */
   return spxSolve(lpi, SPxSolver::LEAVE);
}

/** start strong branching - call before any strongbranching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   /* do nothing */
   return SCIP_OKAY;
}

/** end strong branching - call after any strongbranching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   /* do nothing */
   return SCIP_OKAY;
}

/** performs strong branching iterations on all candidates */
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
   SPxSCIP* spx;
   SPxSolver::VarStatus* rowstat;
   SPxSolver::VarStatus* colstat;
   SPxSolver::Status status;
   SCIP_Real oldlb;
   SCIP_Real oldub;
   SCIP_Real newlb;
   SCIP_Real newub;
   bool error;
   int oldItlim;

   SCIPdebugMessage("calling SCIPlpiStrongbranch() on variable %d (%d iterations)\n", col, itlim);

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   /** remember, whether the last solve call was strong branching, and save/restore basis only once */
   spx = lpi->spx;
   rowstat = new SPxSolver::VarStatus[spx->nRows()];
   colstat = new SPxSolver::VarStatus[spx->nCols()];
   status = SPxSolver::UNKNOWN;
   error = false;

   /* get basis and current bounds of column */
   (void)spx->getBasis(rowstat, colstat);
   oldlb = spx->lower(col);
   oldub = spx->upper(col);

   *downvalid = FALSE;
   *upvalid = FALSE;
   if( iter != NULL )
      *iter = 0;

   oldItlim = spx->terminationIter();
   spx->setTerminationIter(itlim);
   lpi->spx->setType(SPxSolver::LEAVE);

   /* down branch */
   newub = EPSCEIL(psol-1.0, 1e-06);
   if( newub >= oldlb - 0.5 )
   {
      SCIPdebugMessage("strong branching down on x%d (%g) with %d iterations\n", col, psol, itlim);

      spx->changeUpper(col, newub);

      status = spx->solve();
      switch( status )
      {
      case SPxSolver::OPTIMAL:
         *down = spx->value();
         *downvalid = TRUE;
         break;
      case SPxSolver::ABORT_TIME: /* Soplex does not return a proven dual bound, if it is aborted */
      case SPxSolver::ABORT_ITER:
      case SPxSolver::ABORT_CYCLING:
         *down = spx->value();
         break;
      case SPxSolver::ABORT_VALUE:
      case SPxSolver::INFEASIBLE:
         *down = spx->terminationValue();
         *downvalid = TRUE;
         break;
      default:
         error = true;
         break;
      }  /*lint !e788*/
      if( iter != NULL )
         (*iter) += spx->iterations();
      spx->changeUpper(col, oldub);
      spx->setBasis(rowstat, colstat);
   }
   else
   {
      *down = spx->terminationValue();
      *downvalid = TRUE;
   }

   /* up branch */
   if( !error )
   {
      newlb = EPSFLOOR(psol+1.0, 1e-06);
      if( newlb <= oldub + 0.5 )
      {
         SCIPdebugMessage("strong branching  up  on x%d (%g) with %d iterations\n", col, psol, itlim);

         spx->changeLower(col, newlb);

         status = spx->solve();
         switch( status )
         {
         case SPxSolver::OPTIMAL:
            *up = spx->value();
            *upvalid = TRUE;
            break;
         case SPxSolver::ABORT_TIME: /* Soplex does not return a proven dual bound, if it is aborted */
         case SPxSolver::ABORT_ITER:
         case SPxSolver::ABORT_CYCLING:
            *up = spx->value();
            break;
         case SPxSolver::ABORT_VALUE:
         case SPxSolver::INFEASIBLE:
            *up = spx->terminationValue();
            *upvalid = TRUE;
            break;
         default:
            error = true;
            break;
         }  /*lint !e788*/
         if( iter != NULL )
            (*iter) += spx->iterations();
         spx->changeLower(col, oldlb);
         spx->setBasis(rowstat, colstat);
      }
      else
      {
         *up = spx->terminationValue();
         *upvalid = TRUE;
      }
   }

   spx->setTerminationIter(oldItlim);

   delete [] rowstat;
   delete [] colstat;

   if( error )
   {
      SCIPerrorMessage("SCIPlpiStrongbranch() returned SoPlex status %d\n", int(status));
      return SCIP_LPERROR;
   }

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
   assert( iter != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if ( iter != NULL )
      *iter = 0;

   for (int j = 0; j < ncols; ++j)
   {
      /* pass call on to lpiStrongbranch() */
      SCIP_CALL( lpiStrongbranch(lpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter) );
   }
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
   assert( iter != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if ( iter != NULL )
      *iter = 0;

   for (int j = 0; j < ncols; ++j)
   {
      /* pass call on to lpiStrongbranch() */
      SCIP_CALL( lpiStrongbranch(lpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter) );
   }
   return SCIP_OKAY;
}
/**@} */




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
   SPxBasis::SPxStatus basestatus;

   SCIPdebugMessage("calling SCIPlpiGetSolFeasibility()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(primalfeasible != NULL);
   assert(dualfeasible != NULL);

   basestatus = lpi->spx->basis().status();
   *primalfeasible = (basestatus == SPxBasis::PRIMAL || basestatus == SPxBasis::OPTIMAL);
   *dualfeasible = (basestatus == SPxBasis::DUAL || basestatus == SPxBasis::OPTIMAL);

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
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::UNBOUNDED);
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
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::UNBOUNDED);
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::UNBOUNDED && lpi->spx->basis().status() == SPxBasis::PRIMAL);
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::INFEASIBLE);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SPxBasis::SPxStatus basestatus;

   SCIPdebugMessage("calling SCIPlpiIsPrimalFeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   basestatus = lpi->spx->basis().status();

   return (basestatus == SPxBasis::PRIMAL || basestatus == SPxBasis::OPTIMAL);
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
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::INFEASIBLE);
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
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::INFEASIBLE);
}

/** returns TRUE iff LP is dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::INFEASIBLE && lpi->spx->basis().status() == SPxBasis::DUAL);
}

/** returns TRUE iff LP is dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::UNBOUNDED);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SPxBasis::SPxStatus basestatus;

   SCIPdebugMessage("calling SCIPlpiIsDualFeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   basestatus = lpi->spx->basis().status();

   return (basestatus == SPxBasis::DUAL || basestatus == SPxBasis::OPTIMAL);
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsOptimal()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::OPTIMAL);
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsStable()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() != SPxSolver::ERROR && lpi->spx->getStatus() != SPxSolver::SINGULAR);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsObjlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::ABORT_VALUE);
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsIterlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::ABORT_ITER);
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::ABORT_TIME);
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return static_cast<int>(lpi->spx->getStatus());
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiIgnoreInstability()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   /* instable situations cannot be ignored */
   *success = FALSE;

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
   assert(lpi->spx != NULL);
   assert(objval != NULL);

   *objval = lpi->spx->value();

   return SCIP_OKAY;
}

/** gets primal and dual solution vectors */
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
   assert(lpi->spx != NULL);

   if( objval != NULL )
      *objval = lpi->spx->value();

   if( primsol != NULL )
   {
      Vector tmp(lpi->spx->nCols(), primsol);
      (void)lpi->spx->getPrimal(tmp);
   }
   if( dualsol != NULL )
   {
      Vector tmp(lpi->spx->nRows(), dualsol);
      (void)lpi->spx->getDual(tmp);
   }
   if( activity != NULL )
   {
      Vector tmp(lpi->spx->nRows(), activity);
      (void)lpi->spx->getSlacks(tmp);  /* in SoPlex, the activities are called "slacks" */
   }
   if( redcost != NULL )
   {
      Vector tmp(lpi->spx->nCols(), redcost);
      (void)lpi->spx->getRedCost(tmp);
   }

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiGetPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   SCIPerrorMessage("SCIPlpiGetPrimalRay() not supported by SoPlex\n");

   return SCIP_LPERROR;
}

/** gets dual farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual farkas row multipliers */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetDualfarkas()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   Vector tmp(lpi->spx->nRows(), dualfarkas);
   (void)lpi->spx->getDualfarkas(tmp);

   return SCIP_OKAY;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetIterations()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   *iterations = lpi->spx->iterations();

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
   SCIPdebugMessage("calling SCIPlpiGetRealSolQuality()\n");

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
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetBase()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   if( rstat != NULL )
   {
      for( i = 0; i < lpi->spx->nRows(); ++i )
      {
         switch( lpi->spx->getBasisRowStatus(i) )
         {
         case SPxSolver::BASIC:
            rstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/
            break;
         case SPxSolver::FIXED:
         case SPxSolver::ON_LOWER:
            rstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
            break;
         case SPxSolver::ON_UPPER:
            rstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
            break;
         case SPxSolver::ZERO:
            SCIPerrorMessage("slack variable has basis status ZERO (should not occur)\n");
            return SCIP_LPERROR;
         default:
            SCIPerrorMessage("invalid basis status\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
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
            cstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/
            break;
         case SPxSolver::FIXED:
            assert(lpi->spx->rep() == SPxSolver::COLUMN);
            if( lpi->spx->pVec()[i] - lpi->spx->maxObj()[i] < 0.0 )  /* reduced costs < 0 => UPPER  else => LOWER */
               cstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
            else
               cstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
            break;
         case SPxSolver::ON_LOWER:
            cstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
            break;
         case SPxSolver::ON_UPPER:
            cstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
            break;
         case SPxSolver::ZERO:
            cstat[i] = SCIP_BASESTAT_ZERO; /*lint !e641*/
            break;
         default:
            SCIPerrorMessage("invalid basis status\n");
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
   int i;

   SCIPdebugMessage("calling SCIPlpiSetBase()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(cstat != NULL);
   assert(rstat != NULL);

   invalidateSolution(lpi);

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
         SCIPerrorMessage("slack variable has basis status ZERO (should not occur)\n");
         delete[] spxcstat;
         delete[] spxrstat;
         return SCIP_LPERROR; /*lint !e429*/
      default:
         SCIPerrorMessage("invalid basis status\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
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
         SCIPerrorMessage("invalid basis status\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }
   lpi->spx->setBasis(spxrstat, spxcstat);

   delete[] spxcstat;
   delete[] spxrstat;

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetBasisInd()\n");

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
SCIP_RETCODE SCIPlpiGetBInvRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvRow()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   Vector x(lpi->spx->nRows(), coef); /* row of B^-1 has nrows entries */
   DVector e(lpi->spx->nRows());      /* prepare unit vector */
   e.clear();
   e[r] = 1.0;

   /* solve system "x = e_r^T * B^-1" to get r'th row of B^-1 */
   lpi->spx->basis().coSolve(x, e);

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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   Vector x(lpi->spx->nRows(), coef); /* column of B^-1 has nrows entries */
   DVector e(lpi->spx->nRows());      /* prepare unit vector */
   e.clear();
   e[c] = 1.0;

   /* solve system "x = B^-1 * e_c" to get c'th column of B^-1 */
   lpi->spx->basis().solve(x, e);

   return SCIP_OKAY;
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{
   SCIP_Real* buf;
   SCIP_Real* binv;
   int nrows;
   int ncols;
   int c;

   SCIPdebugMessage("calling SCIPlpiGetBInvARow()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   nrows = lpi->spx->nRows();
   ncols = lpi->spx->nCols();

   /* get (or calculate) the row in B^-1 */
   if( binvrow == NULL )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&buf, nrows) );
      SCIP_CALL( SCIPlpiGetBInvRow(lpi, r, buf) );
      binv = buf;
   }
   else
   {
      buf = NULL;
      binv = const_cast<SCIP_Real*>(binvrow);
   }
   assert(binv != NULL);

   /* calculate the scalar product of the row in B^-1 and A */
   soplex::Vector binvvec(nrows, binv);
   for( c = 0; c < ncols; ++c )
      coef[c] = binvvec * lpi->spx->colVector(c);  /* scalar product */ /*lint !e1702*/

   /* free memory if it was temporarily allocated */
   BMSfreeMemoryArrayNull(&buf);

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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   Vector x(lpi->spx->nRows(), coef); /* column of B^-1 has nrows entries */
   DVector col(lpi->spx->nRows());

   /* extract column c of A */
   col = lpi->spx->colVector(c);

   /* solve system "x = B^-1 * A_c" to get c'th column of B^-1 * A */
   lpi->spx->basis().solve(x, col);

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
   int ncols;
   int nrows;

   SCIPdebugMessage("calling SCIPlpiGetState()\n");

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpistate != NULL);

   ncols = lpi->spx->nCols();
   nrows = lpi->spx->nRows();
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

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpistate != NULL);

   lpncols = lpi->spx->nCols();
   lpnrows = lpi->spx->nRows();
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
      SCIP_Real bnd = lpi->spx->lower(i);
      if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         bnd = lpi->spx->lower(i);
         if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
            lpi->cstat[i] = SCIP_BASESTAT_ZERO;  /* variable is free */
         else
            lpi->cstat[i] = SCIP_BASESTAT_UPPER; /* use finite upper bound */
      }
      else
         lpi->cstat[i] = SCIP_BASESTAT_LOWER;    /* use finite lower bound */
   }
   for( i = lpistate->nrows; i < lpnrows; ++i )
      lpi->rstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/

   /* load basis information */
   SCIP_CALL( SCIPlpiSetBase(lpi, lpi->cstat, lpi->rstat) );

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiClearState()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   lpi->spx->reLoad();

   return SCIP_OKAY;
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiFreeState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiFreeState()\n");

   assert(lpi != NULL);
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
{  /*lint --e{715}*/
   return TRUE;
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             /*lpi*/,            /**< LP interface structure */
   const char*           /*fname*/           /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadState()\n");

   SCIPerrorMessage("SCIPlpiReadState() not implemented yet in SoPlex interface\n");

   return SCIP_INVALIDCALL;
}

/** writes LP state (like basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             /*lpi*/,            /**< LP interface structure */
   const char*           /*fname*/           /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiWriteState()\n");

   SCIPerrorMessage("SCIPlpiWriteState() not implemented yet in SoPlex interface\n");

   return SCIP_INVALIDCALL;
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
   assert(lpi->spx != NULL);
   assert(ival != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      *ival = lpi->spx->getFromScratch();
      break;
   case SCIP_LPPAR_LPINFO:
      *ival = lpi->spx->getLpInfo();
      break;
   case SCIP_LPPAR_LPITLIM:
      *ival = lpi->spx->terminationIter();
      break;
   case SCIP_LPPAR_PRICING:
      *ival = (int)lpi->pricing;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

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
   assert(lpi->spx != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      assert(ival == TRUE || ival == FALSE);
      lpi->spx->setFromScratch(bool(ival));
      break;
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      lpi->spx->setLpInfo(bool(ival));
      break;
   case SCIP_LPPAR_LPITLIM:
      lpi->spx->setTerminationIter(ival);
      break;
   case SCIP_LPPAR_PRICING:
      lpi->pricing = (SCIP_PRICING)ival;
      switch( lpi->pricing )
      {
      case SCIP_PRICING_LPIDEFAULT:
      case SCIP_PRICING_AUTO:
      case SCIP_PRICING_FULL:
         lpi->spx->setAutoPricer();
         break;
      case SCIP_PRICING_PARTIAL:
         lpi->spx->setParmultPricer();
         break;
      case SCIP_PRICING_STEEP:
      case SCIP_PRICING_STEEPQSTART:
         lpi->spx->setSteepPricer();
	 break;
      case SCIP_PRICING_DEVEX:
         lpi->spx->setDevexPricer();
	 break;
      default:
         return SCIP_LPERROR;
      }
      break;

   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

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
   }  /*lint !e788*/

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
   }  /*lint !e788*/

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

   return soplex::infinity;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             /*lpi*/,            /**< LP interface structure */
   SCIP_Real             val
   )
{
   SCIPdebugMessage("calling SCIPlpiIsInfinity()\n");

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
SCIP_Bool fileExists(
   const char*           filename            /**< file name */
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
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   if( !fileExists(fname) )
      return SCIP_NOFILE;

   if( !lpi->spx->readFile(fname) )
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiWriteLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   lpi->spx->dumpFile(fname);

   return SCIP_OKAY;
}

/**@} */
