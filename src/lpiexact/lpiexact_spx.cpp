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

/**@file   lpiexact_spx.cpp
 * @ingroup LPIS
 * @brief  exact LP interface for SoPlex version 2.0 and higher
 * @author Leon Eifler
 *
 * This is an implementation of SCIP's LP interface for SoPlex using the extended and improved interface of SoPlex 2.0
 *
 * For debugging purposes, the SoPlex results can be double checked with CPLEX if WITH_LPSCHECK is defined. This may
 * yield false positives, since the LP is dumped to a file for transfering it to CPLEX, hence, precision may be lost.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "lpiexact/type_lpiexact.h"
#include "lpiexact/lpiexact.h"
#include "scip/rational.h"

#define STRONGBRANCH_RESTOREBASIS            /**< if defined then in SCIPlpiStrongbranch() we restore the basis after the
                                              *   down branch and after the up branch; if false only after the end of a
                                              *   strong branching phase, which however seems to mostly increase strong
                                              *   branching time and iterations */

/* check the return value of setParam methods */
#define CHECK_SOPLEX_PARAM(x)                                                           \
   if( !x )                                                                             \
   {                                                                                    \
      SCIPmessagePrintWarning(_messagehdlr, "SoPlex: unsupported parameter value\n");   \
   }

/* remember the original value of the SCIP_DEBUG define and undefine it */
#ifdef SCIP_DEBUG
#define ___DEBUG
#undef SCIP_DEBUG
#endif

/* disable -Wclass-memaccess warnings due to dubious memcpy/realloc calls in SoPlex headers, e.g.,
 * dataarray.h:314:16: warning: ‘void* memcpy(void*, const void*, size_t)’ writing to an object of type ‘struct soplex::SPxParMultPR::SPxParMultPr_Tmp’ with no trivial copy-assignment; use copy-assignment or copy-initialization instead [-Wclass-memaccess]
 */
#ifdef __GNUC__
#if __GNUC__ >= 8
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif
#endif

/* include SoPlex solver */
#include "soplex.h"

/* define subversion for versions <= 1.5.0.1 */
#ifndef SOPLEX_SUBVERSION
#define SOPLEX_SUBVERSION 0
#endif
/* define API version for versions <= 3.0.0 */
#ifndef SOPLEX_APIVERSION
#define SOPLEX_APIVERSION 0
#endif

/* check version */
#if (SOPLEX_APIVERSION < 11)
#error "This interface requires SoPlex with API version 11 or higher, which is available from SoPlex 5.0.0."
#endif

#if (SOPLEX_APIVERSION <= 5)
#include "spxgithash.h"
#endif

#ifndef SOPLEX_WITH_GMP
#error "SOPLEX_WITH_GMP not defined: SoPlex must be build with GMP to be used with LPIEX=spx"
#endif

#ifndef SOPLEX_WITH_MPFR
#error "SOPLEX_WITH_MPFR not defined: SoPlex must be build with MPFR to be used with LPIEX=spx"
#endif

#ifndef SOPLEX_WITH_BOOST
#error "SOPLEX_WITH_BOOST not defined: SoPlex must be build with BOOST to be used with LPIEX=spx"
#endif

/* reset the SCIP_DEBUG define to its original SCIP value */
#undef SCIP_DEBUG
#ifdef ___DEBUG
#define SCIP_DEBUG
#undef ___DEBUG
#endif

/* define snprintf when using a too old MSVC version */
#if defined(_MSC_VER) && _MSC_VER < 1900
#ifndef snprintf
#define snprintf _snprintf
#endif
#endif

#define SOPLEX_VERBLEVEL                5    /**< verbosity level for LPINFO */

#include "scip/pub_message.h"
#include "scip/struct_rational.h"

/********************************************************************/
/*----------------------------- C++ --------------------------------*/
/********************************************************************/

/* MP@LE I suggest to remove the following and replace 0/NULL by nullptr */
/* in C++ we have to use "0" instead of "(void*)0" */
#undef NULL
#define NULL 0

#include <cassert>
using namespace soplex;


/** Macro for a single SoPlex call for which exceptions have to be catched - return an LP error. We
 *  make no distinction between different exception types, e.g., between memory allocation and other
 *  exceptions.
 */
#ifndef NDEBUG
#define SOPLEX_TRY(messagehdlr, x)  do                                  \
   {                                                                    \
      try                                                               \
      {                                                                 \
         (x);                                                           \
      }                                                                 \
      catch( const SPxMemoryException& E )                              \
      {                                                                 \
         std::string s = E.what();                                      \
         SCIPerrorMessage("SoPlex threw a memory exception: %s\n", s.c_str()); \
         return SCIP_ERROR;                                             \
      }                                                                 \
      catch( const SPxException& E )                                    \
      {                                                                 \
         std::string s = E.what();                                      \
         SCIPmessagePrintWarning((messagehdlr), "SoPlex threw an exception: %s\n", s.c_str()); \
         return SCIP_LPERROR;                                           \
      }                                                                 \
   }                                                                    \
   while( FALSE )

#else
#define SOPLEX_TRY(messagehdlr, x)  do                                  \
   {                                                                    \
      try                                                               \
      {                                                                 \
         (x);                                                           \
      }                                                                 \
      catch( const SPxMemoryException& E )                              \
      {                                                                 \
         std::string s = E.what();                                      \
         SCIPerrorMessage("SoPlex threw a memory exception: %s\n", s.c_str()); \
         return SCIP_ERROR;                                             \
      }                                                                 \
      catch( const SPxException& )                                      \
      {                                                                 \
         return SCIP_LPERROR;                                           \
      }                                                                 \
   }                                                                    \
   while( FALSE )
#endif

/* Macro for a single SoPlex call for which exceptions have to be catched - abort if they
 * arise. SCIP_ABORT() is not accessible here.
 */
#define SOPLEX_TRY_ABORT(x)  do                                         \
   {                                                                    \
      try                                                               \
      {                                                                 \
         (x);                                                           \
      }                                                                 \
      catch( const SPxException& E )                                    \
      {                                                                 \
         std::string s = E.what();                                      \
         SCIPerrorMessage("SoPlex threw an exception: %s\n", s.c_str()); \
         abort();                                                       \
      }                                                                 \
   }                                                                    \
   while( FALSE )

/* Set the value of a SCIP_Rational* from a SoPlex Rational */
static void RsetSpxR(
   SCIP_LPIEXACT*        lpi,                /**< exact lpi*/
   SCIP_Rational*        r,                  /**< scip rational */
   const soplex::Rational& spxr              /**< soplex rational */
   )
{
   if( SCIPlpiExactIsInfinity(lpi, double(spxr)) )
   {
      RatSetString(r, "inf");
   }
   else if( SCIPlpiExactIsInfinity(lpi, -double(spxr)) )
   {
      RatSetString(r, "-inf");
   }
   else
   {
#if defined(SCIP_WITH_BOOST) && defined(SCIP_WITH_GMP)
      r->val = spxr;
#else
      r->val = (double) spxr;
#endif
      r->isinf = false;
      r->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
   }
}

/* Set the value of a SoPlex Rational vector from a SCIP_Rational** */
static void RsetSpxVector(
   SCIP_LPIEXACT*        lpi,                /**< exact LPI */
   SCIP_Rational**       r,                  /**< SCIP_Rational array */
   VectorRational        src                 /**< SoPlex rational vector */
   )
{
   for( int i = 0; i < src.dim(); ++i )
   {
      assert(r[i] != NULL);

      RsetSpxR(lpi, r[i], src[i]);
   }
}

/** set the value of a SoPlex rational from a SCIP_Rational*
 * @todo exip: there seems to be something wrong with the = of spx rational */
static void SpxRSetRat(
   SCIP_LPIEXACT*        lpi,                /**< exact LPI */
   soplex::Rational&     spxr,               /**< SoPlex Rational*/
   SCIP_Rational*        src                 /**< SCIP_Rational */
   )
{
   if( RatIsAbsInfinity(src) )
   {
      if( RatIsPositive(src) )
         spxr = soplex::Rational(SCIPlpiExactInfinity(lpi));
      else
         spxr = soplex::Rational(-SCIPlpiExactInfinity(lpi));
   }
   else
   {
#if defined(SOPLEX_WITH_GMP) && defined(SCIP_WITH_BOOST) && defined(SCIP_WITH_GMP)
      spxr = soplex::Rational(*RatGetGMP(src));
#else
      spxr = soplex::Rational(RatApproxReal(src));
#endif
   }
}



/** SCIP's SoPlex class */
class SPxexSCIP : public SoPlex
{/*lint !e1790*/
   bool                  _lpinfo;
   bool                  _fromscratch;
   char*                 _probname;
   DataArray<SPxSolver::VarStatus> _colStat; /**< column basis status used for strong branching */
   DataArray<SPxSolver::VarStatus> _rowStat; /**< row basis status used for strong branching */
   SCIP_MESSAGEHDLR*     _messagehdlr;       /**< messagehdlr handler for printing messages, or NULL */

public:
   SPxexSCIP(
      SCIP_MESSAGEHDLR*  messagehdlr = NULL, /**< message handler */
      const char*        probname = NULL     /**< name of problem */
      )
      : _lpinfo(false),
        _fromscratch(false),
        _probname(NULL),
        _colStat(0),
        _rowStat(0),
        _messagehdlr(messagehdlr)
   {
      if ( probname != NULL )
         SOPLEX_TRY_ABORT( setProbname(probname) );

#if SOPLEX_APIVERSION >= 2
      setBoolParam(SoPlex::ENSURERAY, true);
#endif
   }

   virtual ~SPxexSCIP()
   {
      if( _probname != NULL )
         spx_free(_probname); /*lint !e1551*/

      freePreStrongbranchingBasis(); /*lint !e1551*/
   }/*lint -e1579*/

   /** get objective limit according to objective sense */
   Real getObjLimit() const
   {
      return (intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE)
         ? realParam(SoPlex::OBJLIMIT_UPPER)
         : realParam(SoPlex::OBJLIMIT_LOWER);
   }

   // @todo realize this with a member variable as before
   bool getFromScratch() const
   {
      return _fromscratch;
   }

   void setFromScratch(bool fs)
   {
      _fromscratch = fs;
   }

   // @todo member variable?
   bool getLpInfo() const
   {
      return _lpinfo;
   }

   void setLpInfo(bool lpinfo)
   {
      _lpinfo = lpinfo;
   }

   // @todo member variable?
   void setProbname(const char* probname)
   {
      size_t len;

      assert(probname != NULL);
      if( _probname != NULL )
         spx_free(_probname);

      len = strlen(probname);
      spx_alloc(_probname, len + 1);
      memcpy(_probname, probname, len + 1);
   }

   void setRep(SPxSolver::Representation p_rep)
   {
      if( p_rep == SPxSolver::COLUMN && intParam(REPRESENTATION) == REPRESENTATION_ROW )
      {
         SCIPdebugMessage("switching to column representation of the basis\n");
         CHECK_SOPLEX_PARAM(setIntParam(REPRESENTATION, REPRESENTATION_COLUMN));
      }
      else if( (p_rep == SPxSolver::ROW && intParam(REPRESENTATION) == REPRESENTATION_COLUMN) )
      {
         SCIPdebugMessage("switching to row representation of the basis\n");
         CHECK_SOPLEX_PARAM(setIntParam(REPRESENTATION, REPRESENTATION_ROW));
      }
   }

#ifndef NDEBUG
   bool checkConsistentBounds() const
   {
      for( int i = 0; i < numColsRational(); ++i )
      {
         if( lowerRational(i) > upperRational(i) )
         {
            SCIPerrorMessage("inconsistent bounds on column %d: lower=%s, upper=%s\n",
               i, lowerRational(i).str().c_str(), upperRational(i).str().c_str());
            return false;
         }
      }

      return true;
   }

   bool checkConsistentSides() const
   {
      for( int i = 0; i < numRowsRational(); ++i )
      {
         if( lhsRational(i) > rhsRational(i) )
         {
            SCIPerrorMessage("inconsistent sides on row %d: lhs=%s, rhs=%s\n",
              i, lhsRational(i).str().c_str(), rhsRational(i).str().c_str());
            return false;
         }
      }

      return true;
   }
#endif

   void trySolve(bool printwarning = true)
   {
      Real timespent;
      Real timelimit;

      try
      {
         (void) optimize();
      }
      catch(const SPxException& x)
      {
         std::string s = x.what();
         if( printwarning )
         {
            SCIPmessagePrintWarning(_messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
         }

         /* since it is not clear if the status in SoPlex are set correctly
          * we want to make sure that if an error is thrown the status is
          * not OPTIMAL anymore.
          */
         assert(status() != SPxSolver::OPTIMAL);
      }

      assert(intParam(ITERLIMIT) < 0 || numIterations() <= intParam(ITERLIMIT));

      /* update time limit */
      timespent = solveTime();
      if( timespent > 0 )
      {
         /* get current time limit */
         timelimit = realParam(TIMELIMIT);
         if( timelimit > timespent )
            timelimit -= timespent;
         else
            timelimit = 0;
         /* set new time limit */
         assert(timelimit >= 0);
         CHECK_SOPLEX_PARAM(setRealParam(TIMELIMIT, timelimit));
      }
   }

   SPxSolver::Status doSolve(bool printwarning = true)
   {
      SPxOut::Verbosity verbosity;

      SPxSolver::Status spxStatus;

      /* store and set verbosity */
      verbosity = spxout.getVerbosity();
      spxout.setVerbosity((SPxOut::Verbosity)(getLpInfo() ? SOPLEX_VERBLEVEL : 0));

      assert(checkConsistentBounds());
      assert(checkConsistentSides());

      trySolve(printwarning);
      spxStatus = status();

      /* restore verbosity */
      spxout.setVerbosity(verbosity);

      return spxStatus;
   }

   /** save the current basis */
   void savePreStrongbranchingBasis()
   {
      _rowStat.reSize(numRowsRational());
      _colStat.reSize(numColsRational());

      try
      {
         getBasis(_rowStat.get_ptr(), _colStat.get_ptr());
      }
#ifndef NDEBUG
      catch(const SPxException& x)
      {
         std::string s = x.what();
         SCIPmessagePrintWarning(_messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());

         /* since it is not clear if the status in SoPlex are set correctly
          * we want to make sure that if an error is thrown the status is
          * not OPTIMAL anymore.
          */
         assert(status() != SPxSolver::OPTIMAL);
      }
#else
      catch(const SPxException&)
      { }
#endif
   }

   /** restore basis */
   void restorePreStrongbranchingBasis()
   {
      assert(_rowStat.size() == numRowsRational());
      assert(_colStat.size() == numColsRational());

      try
      {
         setBasis(_rowStat.get_ptr(), _colStat.get_ptr());
      }
#ifndef NDEBUG
      catch(const SPxException& x)
      {
         std::string s = x.what();
         SCIPmessagePrintWarning(_messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
      catch(const SPxException&)
      {
#endif
         /* since it is not clear if the status in SoPlex are set correctly
          * we want to make sure that if an error is thrown the status is
          * not OPTIMAL anymore.
          */
         assert(status() != SPxSolver::OPTIMAL);
      }
   }

   /** if basis is in store, delete it without restoring it */
   void freePreStrongbranchingBasis()
   {
      _rowStat.clear();
      _colStat.clear();
   }

   /** is pre-strong-branching basis freed? */
   bool preStrongbranchingBasisFreed() const
   {
      return ((_rowStat.size() == 0 ) && (_colStat.size() == 0));
   }

   /** provides access for temporary storage of basis status of rows */
   DataArray<SPxSolver::VarStatus>& rowStat()
   {
      return _rowStat; /*lint !e1536*/
   }

   /** provides access for temporary storage of basis status or columns */
   DataArray<SPxSolver::VarStatus>& colStat()
   {
      return _colStat; /*lint !e1536*/
   }

}; /*lint !e1748*/




/********************************************************************/
/*-----------------------------  C  --------------------------------*/
/********************************************************************/

#include "lpiexact/lpiexact.h"
#include "scip/bitencode.h"

typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE



/** LP interface */
struct SCIP_LPiExact
{
   SPxexSCIP*            spx;                /**< our SoPlex implementation */
   int*                  cstat;              /**< array for storing column basis status */
   int*                  rstat;              /**< array for storing row basis status */
   int                   cstatsize;          /**< size of cstat array */
   int                   rstatsize;          /**< size of rstat array */
   SCIP_PRICING          pricing;            /**< current pricing strategy */
   SCIP_Bool             solved;             /**< was the current LP solved? */
   SCIP_Real             conditionlimit;     /**< maximum condition number of LP basis counted as stable (-1.0: no limit) */
   SCIP_Bool             checkcondition;     /**< should condition number of LP basis be checked for stability? */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehdlr handler to printing messages, or NULL */
};

/** LPi state stores basis information */
struct SCIP_LPiState
{
   int                   ncols;              /**< number of LP columns */
   int                   nrows;              /**< number of LP rows */
   COLPACKET*            packcstat;          /**< column basis status in compressed form */
   ROWPACKET*            packrstat;          /**< row basis status in compressed form */
};

/** LPi norms to store dual steepest edge */
struct SCIP_LPiNorms
{
   int                   nrows;              /**< number of stored norms corresponding to rows */
   int                   ncols;              /**< number of stored norms corresponding to cols */
   SCIP_Rational**       norms;              /**< norms to be (re)stored */
};



/*
 * dynamic memory arrays
 */

/** resizes cstat array to have at least num entries */
static
SCIP_RETCODE ensureCstatMem(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
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
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
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
   SCIP_LPISTATE*        lpistate,           /**< pointer to LPi state data */
   const int*            cstat,              /**< basis status of columns in unpacked format */
   const int*            rstat               /**< basis status of rows in unpacked format */
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
   const SCIP_LPISTATE*  lpistate,           /**< pointer to LPi state data */
   int*                  cstat,              /**< buffer for storing basis status of columns in unpacked format */
   int*                  rstat               /**< buffer for storing basis status of rows in unpacked format */
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

   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packcstat, nColPackets);
   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packrstat, nRowPackets);
   BMSfreeBlockMemory(blkmem, lpistate);
}




/*
 * local methods
 */


/** marks the current LP to be unsolved */
static
void invalidateSolution(SCIP_LPIEXACT* lpi)
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
static char spxdesc[200];

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiExactGetSolverName(
   void
   )
{
   SCIPdebugMessage("calling SCIPlpiExactGetSolverName()\n");

#if (SOPLEX_SUBVERSION > 0)
   snprintf(spxname, 100, "SoPlex %d.%d.%d.%d", SOPLEX_VERSION/100, (SOPLEX_VERSION % 100)/10, SOPLEX_VERSION % 10, SOPLEX_SUBVERSION); /*lint !e778 !e845*/
#else
   snprintf(spxname, 100, "SoPlex %d.%d.%d", SOPLEX_VERSION/100, (SOPLEX_VERSION % 100)/10, SOPLEX_VERSION % 10); /*lint !e778 !e845*/
#endif
   return spxname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiExactGetSolverDesc(
   void
   )
{
   snprintf(spxdesc, 200, "%s [GitHash: %s]", "Exact Linear Programming Solver developed at Zuse Institute Berlin (soplex.zib.de)", getGitHash());

   return spxdesc;
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiExactGetSolverPointer(
   SCIP_LPIEXACT*        lpi                 /**< pointer to an LP interface structure */
   )
{
   return (void*) lpi->spx;
}

/** pass integrality information about variables to the solver */
SCIP_RETCODE SCIPlpiExactSetIntegralityInformation(
   SCIP_LPIEXACT*        lpi,                /**< pointer to an LP interface structure */
   int                   ncols,              /**< length of integrality array */
   int*                  intInfo             /**< integrality array (0: continuous, 1: integer). May be NULL iff ncols is 0.  */
   )
{
   assert(ncols == lpi->spx->numColsRational() || (ncols == 0 && intInfo == NULL));
   lpi->spx->setIntegralityInformation(ncols, intInfo);
   return SCIP_OKAY;
}

/** informs about availability of a primal simplex solving method */
SCIP_Bool SCIPlpiExactHasPrimalSolve(
   void
   )
{
   return TRUE;
}

/** informs about availability of a dual simplex solving method */
SCIP_Bool SCIPlpiExactHasDualSolve(
   void
   )
{
   return TRUE;
}

/** informs about availability of a barrier solving method */
SCIP_Bool SCIPlpiExactHasBarrierSolve(
   void
   )
{
   return FALSE;
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
{
   assert(lpi != NULL);
   assert(name != NULL);

   /* create SoPlex object */
   SCIP_ALLOC( BMSallocMemory(lpi) );

   /* we use this construction to allocate the memory for the SoPlex class also via the blockmemshell */
   (*lpi)->spx = static_cast<SPxexSCIP*>(BMSallocMemoryCPP(sizeof(SPxexSCIP)));
   SOPLEX_TRY( messagehdlr, (*lpi)->spx = new ((*lpi)->spx) SPxexSCIP(messagehdlr, name) );
   (void) (*lpi)->spx->setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
   (void) (*lpi)->spx->setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
   (void) (*lpi)->spx->setRealParam(SoPlex::FEASTOL, 0.0, TRUE);
   (void) (*lpi)->spx->setRealParam(SoPlex::OPTTOL, 0.0, TRUE);
   (void) (*lpi)->spx->setBoolParam(SoPlex::FORCEBASIC, TRUE);

   (*lpi)->cstat = NULL;
   (*lpi)->rstat = NULL;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->pricing = SCIP_PRICING_LPIDEFAULT;
   (*lpi)->conditionlimit = -1.0;
   (*lpi)->checkcondition = FALSE;
   (*lpi)->messagehdlr = messagehdlr;

   invalidateSolution(*lpi);

   /* set objective sense */
   SCIP_CALL( SCIPlpiExactChgObjsen(*lpi, objsen) );

   /* set default pricing */
   SCIP_CALL( SCIPlpiExactSetIntpar(*lpi, SCIP_LPPAR_PRICING, (int)(*lpi)->pricing) );

   {
      SPxOut::Verbosity verbosity = (*lpi)->spx->spxout.getVerbosity();
      (*lpi)->spx->spxout.setVerbosity((SPxOut::Verbosity)((*lpi)->spx->getLpInfo() ? SOPLEX_VERBLEVEL : 0));
      (*lpi)->spx->printVersion();
      (*lpi)->spx->spxout.setVerbosity(verbosity);
   }

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiExactFree(
   SCIP_LPIEXACT**       lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(*lpi != NULL);
   assert((*lpi)->spx != NULL);

   /* free LP using destructor and free memory via blockmemshell */
   (*lpi)->spx->~SPxexSCIP();
   BMSfreeMemory(&((*lpi)->spx));

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
SCIP_RETCODE SCIPlpiExactLoadColLP(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   SCIP_Rational**       obj,                /**< objective function values of columns */
   SCIP_Rational**       lb,                 /**< lower bounds of columns */
   SCIP_Rational**       ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   SCIP_Rational**       lhs,                /**< left hand sides of rows */
   SCIP_Rational**       rhs,                /**< right hand sides of rows */
   char**                /*rownames*/,       /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each column in ind- and val-array */
   int*                  ind,                /**< row indices of constraint matrix entries */
   SCIP_Rational**       val                 /**< values of constraint matrix entries */
   )
{
#ifndef NDEBUG
   {
      int j;
      for( j = 0; j < nnonz; j++ )
      {
         assert(val[j] != NULL);
         assert(!RatIsZero(val[j]));
      }
   }
#endif

   SCIPdebugMessage("calling SCIPlpiExactLoadColLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(beg != NULL);
   assert(ind != NULL);
   assert(val != NULL);

   invalidateSolution(lpi);
   assert(lpi->spx->preStrongbranchingBasisFreed());

   try
   {
      SPxexSCIP* spx = lpi->spx;
      LPRowSetRational rows(nrows);
      DSVectorRational emptyVector(0);
      int i;


      spx->clearLPRational();

      /* set objective sense */
      (void) spx->setIntParam(SoPlex::OBJSENSE, (objsen == SCIP_OBJSEN_MINIMIZE ? SoPlex::OBJSENSE_MINIMIZE : SoPlex::OBJSENSE_MAXIMIZE));

      /* create empty rows with given sides */
      for( i = 0; i < nrows; ++i )
      {
         soplex::Rational spxlhs;
         soplex::Rational spxrhs;
         SpxRSetRat(lpi, spxlhs, lhs[i]);
         SpxRSetRat(lpi, spxlhs, rhs[i]);
         rows.add(spxlhs, emptyVector, spxrhs);
      }
      spx->addRowsRational(rows);

      /* create column vectors with coefficients and bounds */
      SCIP_CALL( SCIPlpiExactAddCols(lpi, ncols, obj, lb, ub, colnames, nnonz, beg, ind, val) );
      //spx->syncLPReal();
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** adds columns to the LP */
SCIP_RETCODE SCIPlpiExactAddCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   SCIP_Rational**       obj,                /**< objective function values of new columns */
   SCIP_Rational**       lb,                 /**< lower bounds of new columns */
   SCIP_Rational**       ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   int*                  beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   int*                  ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   SCIP_Rational**       val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
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
   assert(nnonz >= 0);
   assert(ncols >= 0);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

#ifndef NDEBUG
   if ( nnonz > 0 )
   {
      /* perform check that no new rows are added - this is likely to be a mistake */
      int nrows = lpi->spx->numRowsRational();
      for (int j = 0; j < nnonz; ++j)
      {
         assert( 0 <= ind[j] && ind[j] < nrows );
         assert( val[j] != NULL );
         assert( !RatIsZero(val[j]) );
      }
   }
#endif

   SPxexSCIP* spx = lpi->spx;
   try
   {
      LPColSetRational cols(ncols);
      DSVectorRational colVector(ncols);
      int start;
      int last;
      int i;

      /* create column vectors with coefficients and bounds */
      for( i = 0; i < ncols; ++i )
      {
         int j;
         soplex::Rational spxlb;
         soplex::Rational spxub;
         soplex::Rational spxobj;
         std::vector<soplex::Rational> spxval;

         SpxRSetRat(lpi, spxlb, lb[i]);
         SpxRSetRat(lpi, spxub, ub[i]);
         SpxRSetRat(lpi, spxobj, obj[i]);

         colVector.clear();
         if( nnonz > 0 )
         {
            soplex::Rational tmp;

            start = beg[i];
            last = (i == ncols-1 ? nnonz : beg[i+1]);
            spxval.reserve(last - start);
            for( j = start; j < last; ++j )
            {
               spxval.push_back(tmp);
               SpxRSetRat(lpi, spxval[j - start], val[j]);
            }
            colVector.add(last - start, ind + start, spxval.data());
         }
         cols.add(spxobj, spxlb, colVector, spxub);
      }
      spx->addColsRational(cols);
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiExactDelCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->numColsRational());

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->removeColRangeRational(firstcol, lastcol) );

   return SCIP_OKAY;
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
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
   assert(dstat != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   ncols = lpi->spx->numColsRational();

   /* SoPlex removeCols() method deletes the columns with dstat[i] < 0, so we have to negate the values */
   for( i = 0; i < ncols; ++i )
      dstat[i] *= -1;

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->removeColsRational(dstat) );

   return SCIP_OKAY;
}

/** adds rows to the LP */
SCIP_RETCODE SCIPlpiExactAddRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   SCIP_Rational**       lhs,                /**< left hand sides of new rows */
   SCIP_Rational**       rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   int*                  beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   int*                  ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   SCIP_Rational**       val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
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

   assert( lpi->spx->preStrongbranchingBasisFreed() );

#ifndef NDEBUG
   if ( nnonz > 0 )
   {
      /* perform check that no new columns are added - this is likely to be a mistake */
      int ncols = lpi->spx->numColsRational();
      for (int j = 0; j < nnonz; ++j)
      {
         assert( !RatIsZero(val[j]) );
         assert( 0 <= ind[j] && ind[j] < ncols );
      }
   }
#endif

   try
   {
      SPxexSCIP* spx = lpi->spx;
      LPRowSetRational rows(nrows);
      DSVectorRational rowVector;
      int start;
      int last;
      int i;

      /* create row vectors with given sides */
      for( i = 0; i < nrows; ++i )
      {
          soplex::Rational spxlhs;
          soplex::Rational spxrhs;
          std::vector<soplex::Rational> spxval;

          SpxRSetRat(lpi, spxlhs, lhs[i]);
          SpxRSetRat(lpi, spxrhs, rhs[i]);

         rowVector.clear();
         if( nnonz > 0 )
         {
            start = beg[i];
            last = (i == nrows-1 ? nnonz : beg[i+1]);
            spxval.reserve(last - start);

            for( int j = start; j < last; ++j )
            {
               soplex::Rational tmp;
               spxval.push_back(tmp);
               SpxRSetRat(lpi, spxval[j - start], val[j]);
            }
            rowVector.add(last - start, ind + start, spxval.data());
         }
         rows.add(spxlhs, rowVector, spxrhs);
      }
      spx->addRowsRational(rows);
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiExactDelRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->numRowsRational());

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->removeRowRangeRational(firstrow, lastrow) );

   return SCIP_OKAY;
}

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiExactDelRowset(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
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

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   nrows = lpi->spx->numRowsRational();

   /* SoPlex removeRows() method deletes the rows with dstat[i] < 0, so we have to negate the values */
   for( i = 0; i < nrows; ++i )
      dstat[i] *= -1;

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->removeRowsRational(dstat) );

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiExactClear(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactClear()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );
   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->clearLPRational() );

   return SCIP_OKAY;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiExactChgBounds(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   int*                  ind,                /**< column indices or NULL if ncols is zero */
   SCIP_Rational**       lb,                 /**< values for the new lower bounds or NULL if ncols is zero */
   SCIP_Rational**       ub                  /**< values for the new upper bounds or NULL if ncols is zero */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiChgBounds()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ncols == 0 || (ind != NULL && lb != NULL && ub != NULL));
   if( ncols <= 0 )
      return SCIP_OKAY;

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      soplex::Rational spxlb;
      soplex::Rational spxub;

      for( i = 0; i < ncols; ++i )
      {
         assert(0 <= ind[i] && ind[i] < lpi->spx->numColsRational());
         assert(lb[i] != NULL && ub[i] != NULL);

         if( RatIsInfinity(lb[i]) )
         {
            SCIPerrorMessage("LP Error: fixing lower bound for variable %d to infinity.\n", ind[i]);
            return SCIP_LPERROR;
         }
         if( RatIsNegInfinity(ub[i]) )
         {
            SCIPerrorMessage("LP Error: fixing upper bound for variable %d to -infinity.\n", ind[i]);
            return SCIP_LPERROR;
         }

         SpxRSetRat(lpi, spxlb, lb[i]);
         SpxRSetRat(lpi, spxub, ub[i]);

         lpi->spx->changeBoundsRational(ind[i], spxlb, spxub);
         assert(lpi->spx->lowerRational(ind[i]) <= lpi->spx->upperRational(ind[i]));
      }
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiExactChgSides(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   int*                  ind,                /**< row indices */
   SCIP_Rational**       lhs,                /**< new values for left hand sides */
   SCIP_Rational**       rhs                 /**< new values for right hand sides */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiChgSides()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   if( nrows <= 0 )
      return SCIP_OKAY;

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      soplex::Rational spxlhs;
      soplex::Rational spxrhs;

      for( i = 0; i < nrows; ++i )
      {
         assert(0 <= ind[i] && ind[i] < lpi->spx->numRowsRational());

         SpxRSetRat(lpi, spxlhs, lhs[i]);
         SpxRSetRat(lpi, spxrhs, rhs[i]);

         lpi->spx->changeRangeRational(ind[i], spxlhs, spxrhs);
         assert(lpi->spx->lhsRational(ind[i]) <= lpi->spx->rhsRational(ind[i]));
      }
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiExactChgCoef(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Rational*        newval              /**< new value of coefficient */
   )
{
   soplex::Rational spxval;

   SCIPdebugMessage("calling SCIPlpiChgCoef()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= row && row < lpi->spx->numRowsRational());
   assert(0 <= col && col < lpi->spx->numColsRational());

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SpxRSetRat(lpi, spxval, newval);
   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->changeElementRational(row, col, spxval) );

   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiExactChgObjsen(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgObjsen()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SOPLEX_TRY( lpi->messagehdlr, (void) lpi->spx->setIntParam(SoPlex::OBJSENSE, objsen == SCIP_OBJSEN_MINIMIZE ? SoPlex::OBJSENSE_MINIMIZE : SoPlex::OBJSENSE_MAXIMIZE ) );

   return SCIP_OKAY;
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiExactChgObj(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   SCIP_Rational**       obj                 /**< new objective values for columns */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiChgObj()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(obj != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      soplex::Rational spxobj;

      for( i = 0; i < ncols; ++i )
      {
         assert(obj[i] != NULL);
         assert(0 <= ind[i] && ind[i] < lpi->spx->numColsRational());

         SpxRSetRat(lpi, spxobj, obj[i]);
         lpi->spx->changeObjRational(ind[i], spxobj);
      }
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}


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
{
   SCIPdebugMessage("calling SCIPlpiExactGetNRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(nrows != NULL);

   *nrows = lpi->spx->numRowsRational();

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiExactGetNCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactGetNCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ncols != NULL);

   *ncols = lpi->spx->numColsRational();

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiExactGetNNonz(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiExactGetNNonz()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(nnonz != NULL);

   /* SoPlex has no direct method to return the number of nonzeros, so we have to count them manually */
   *nnonz = 0;
   if( lpi->spx->numRowsRational() < lpi->spx->numColsRational() )
   {
      for( i = 0; i < lpi->spx->numRowsRational(); ++i )
         (*nnonz) += lpi->spx->rowVectorRational(i).size();
   }
   else
   {
      for( i = 0; i < lpi->spx->numColsRational(); ++i )
         (*nnonz) += lpi->spx->colVectorRational(i).size();
   }

   return SCIP_OKAY;
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiExactGetCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   SCIP_Rational**       lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Rational**       ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_Rational**       val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int i;
   int j;

   SCIPdebugMessage("calling SCIPlpiExactGetCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->numColsRational());
   assert((lb != NULL && ub != NULL) || (lb == NULL && ub == NULL));
   assert((nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL));

   if( lb != NULL )
   {
      /** @todo exip: what about scaling? */
      if( lpi->spx->boolParam(SoPlex::PERSISTENTSCALING) )
      {
         const VectorRational& lbvec = lpi->spx->lowerRational();
         const VectorRational& ubvec = lpi->spx->upperRational();

         for( i = firstcol; i <= lastcol; ++i )
         {
            RsetSpxR(lpi, lb[i-firstcol], lbvec[i]);
            RsetSpxR(lpi, ub[i-firstcol], ubvec[i]);
         }
      }
      else
      {
         const VectorRational& lbvec = lpi->spx->lowerRational();
         const VectorRational& ubvec = lpi->spx->upperRational();
         for( i = firstcol; i <= lastcol; ++i )
         {
            RsetSpxR(lpi, lb[i-firstcol], lbvec[i]);
            RsetSpxR(lpi, ub[i-firstcol], ubvec[i]);
         }
      }
   }

   if( nnonz != NULL )
   {
      *nnonz = 0;
      for( i = firstcol; i <= lastcol; ++i )
      {
         beg[i-firstcol] = *nnonz;

         /** @todo exip: what about scaling? */
         if( lpi->spx->boolParam(SoPlex::PERSISTENTSCALING) && FALSE )
         {
         }
         else
         {
            const SVectorRational& cvec = lpi->spx->colVectorRational(i);
            for( j = 0; j < cvec.size(); ++j )
            {
               ind[*nnonz] = cvec.index(j);
               RsetSpxR(lpi, val[*nnonz], cvec.value(j));
               (*nnonz)++;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiExactGetRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_Rational**       lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_Rational**       rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Rational**       val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int i;
   int j;

   SCIPdebugMessage("calling SCIPlpiExactGetRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->numRowsRational());
   assert((lhs != NULL && rhs != NULL) || (lhs == NULL && rhs == NULL));
   assert((nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL));

   if( lhs != NULL )
   {
      /** @todo exip: what about scaling? */
      if( lpi->spx->boolParam(SoPlex::PERSISTENTSCALING) && FALSE )
      {
      }
      else
      {
         const VectorRational& lhsvec = lpi->spx->lhsRational();
         const VectorRational& rhsvec = lpi->spx->rhsRational();
         for( i = firstrow; i <= lastrow; ++i )
         {
            RsetSpxR(lpi, lhs[i-firstrow], lhsvec[i]);
            RsetSpxR(lpi, rhs[i-firstrow], rhsvec[i]);
         }
      }
   }

   if( nnonz != NULL )
   {
      *nnonz = 0;
      for( i = firstrow; i <= lastrow; ++i )
      {
         beg[i-firstrow] = *nnonz;

         /** @todo exip: what to do about scaling? */
         if( lpi->spx->boolParam(SoPlex::PERSISTENTSCALING) )
         {
         }
         else
         {
            const SVectorRational& rvec = lpi->spx->rowVectorRational(i);
            for( j = 0; j < rvec.size(); ++j )
            {
               ind[*nnonz] = rvec.index(j);
               RsetSpxR(lpi, val[*nnonz], rvec.value(j));
               (*nnonz)++;
            }
         }
      }
   }

   return SCIP_OKAY;
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
{
   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( colnames != NULL || namestoragesize == 0 );
   assert( namestorage != NULL || namestoragesize == 0 );
   assert( namestoragesize >= 0 );
   assert( storageleft != NULL );
   assert( 0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->numColsRational() );

   SCIPdebugMessage("getting column names %d to %d\n", firstcol, lastcol);

//    lpi->spx->getColNames(firstcol, lastcol, colnames, namestorage, namestoragesize, storageleft);

   return SCIP_OKAY;
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
{
   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( rownames != NULL || namestoragesize == 0 );
   assert( namestorage != NULL || namestoragesize == 0 );
   assert( namestoragesize >= 0 );
   assert( storageleft != NULL );
   assert( 0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->numRowsRational() );

   SCIPdebugMessage("getting row names %d to %d\n", firstrow, lastrow);

//    lpi->spx->getRowNames(firstrow, lastrow, rownames, namestorage, namestoragesize, storageleft);

   return SCIP_OKAY;
}

/** gets objective sense of the LP */
SCIP_RETCODE SCIPlpiExactGetObjsen(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactGetObjsen()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(objsen != NULL);

   *objsen = (lpi->spx->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE) ? SCIP_OBJSEN_MINIMIZE : SCIP_OBJSEN_MAXIMIZE;

   return SCIP_OKAY;
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiExactGetObj(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Rational**       vals                /**< array to store objective coefficients */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiExactGetObj()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->numColsRational());
   assert(vals != NULL);

   for( i = firstcol; i <= lastcol; ++i )
   {
      assert(vals[i-firstcol] != NULL);
      RsetSpxR(lpi, vals[i-firstcol], lpi->spx->objRational(i));
   }

   return SCIP_OKAY;
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiExactGetBounds(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective value for */
   int                   lastcol,            /**< last column to get objective value for */
   SCIP_Rational**       lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Rational**       ubs                 /**< array to store upper bound values, or NULL */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiExactGetBounds()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->numColsRational());

   for( i = firstcol; i <= lastcol; ++i )
   {
      if( lbs != NULL )
      {
         assert(lbs[i-firstcol] != NULL);
         RsetSpxR(lpi, lbs[i-firstcol], lpi->spx->lowerRational(i));
      }
      if( ubs != NULL )
      {
         assert(ubs[i-firstcol] != NULL);
         RsetSpxR(lpi, ubs[i-firstcol], lpi->spx->upperRational(i));
      }
   }

   return SCIP_OKAY;
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiExactGetSides(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Rational**       lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Rational**       rhss                /**< array to store right hand side values, or NULL */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiExactGetSides()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->numRowsRational());

   for( i = firstrow; i <= lastrow; ++i )
   {
      if( lhss != NULL )
      {
         assert(lhss[i-firstrow] != NULL);
         RsetSpxR(lpi, lhss[i-firstrow], lpi->spx->lhsRational(i));
      }
      if( rhss != NULL )
      {
         assert(rhss[i-firstrow] != NULL);
         RsetSpxR(lpi, rhss[i-firstrow], lpi->spx->rhsReal(i));
      }
   }

   return SCIP_OKAY;
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiExactGetCoef(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Rational*        val                 /**< pointer to store the value of the coefficient */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactGetCoef()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= col && col < lpi->spx->numColsRational());
   assert(0 <= row && row < lpi->spx->numRowsRational());
   assert(val != NULL);

   RsetSpxR(lpi, val, lpi->spx->colVectorRational(col)[row]);

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
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->spx != NULL );

   SPxOut::Verbosity verbosity;
   /* store and set verbosity */
   verbosity = lpi->spx->spxout.getVerbosity();
   lpi->spx->spxout.setVerbosity((SPxOut::Verbosity)(lpi->spx->getLpInfo() ? SOPLEX_VERBLEVEL : 0));

   SCIPdebugMessage("calling exact SoPlex solve(): %d cols, %d rows\n", lpi->spx->numColsRational(), lpi->spx->numRowsRational());

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

#ifdef WITH_LPSCHECK
   lpi->spx->setDoubleCheck(CHECK_SPXSOLVE);
#endif

   /* delete starting basis if solving from scratch */
   if( lpi->spx->getFromScratch() )
   {
      try
      {
         lpi->spx->clearBasis();
      }
#ifndef NDEBUG
      catch(const SPxException& x)
      {
         std::string s = x.what();
         SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
      catch(const SPxException&)
      {
#endif
         assert( lpi->spx->status() != SPxSolver::OPTIMAL );
         return SCIP_LPERROR;
      }
   }
   assert(!lpi->spx->getFromScratch() || lpi->spx->status() == SPxSolver::NO_PROBLEM);

   SPxSolver::Status status = lpi->spx->doSolve();
   SCIPdebugMessage(" -> SoPlex status: %d, basis status: %d\n", lpi->spx->status(), lpi->spx->basisStatus());
   lpi->solved = TRUE;

   /* restore verbosity */
   lpi->spx->spxout.setVerbosity(verbosity);

   switch( status )
   {
   case SPxSolver::ABORT_TIME:
   case SPxSolver::ABORT_ITER:
   case SPxSolver::ABORT_VALUE:
   case SPxSolver::SINGULAR:
   case SPxSolver::REGULAR:
   case SPxSolver::OPTIMAL:
#if SOPLEX_APIVERSION >= 3
   case SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS:
#endif
   case SPxSolver::UNBOUNDED:
   case SPxSolver::INFEASIBLE:
      return SCIP_OKAY;
   case SPxSolver::UNKNOWN:
   default:
      return SCIP_LPERROR;
   }  /*lint !e788*/
}

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiExactSolvePrimal(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiSolvePrimal()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   (void) lpi->spx->setIntParam(SoPlex::ALGORITHM, SoPlex::ALGORITHM_PRIMAL);
   return spxSolve(lpi);
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiExactSolveDual(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiSolveDual()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   (void) lpi->spx->setIntParam(SoPlex::ALGORITHM, SoPlex::ALGORITHM_DUAL);
   return spxSolve(lpi);
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiExactSolveBarrier(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   SCIPdebugMessage("calling SCIPlpiSolveBarrier()\n");

   /* Since SoPlex does not support barrier we switch to DUAL */
   return SCIPlpiExactSolveDual(lpi);
}

/** start strong branching - call before any strongbranching */
SCIP_RETCODE SCIPlpiExactStartStrongbranch(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );
   lpi->spx->savePreStrongbranchingBasis();

   return SCIP_OKAY;
}

/** end strong branching - call after any strongbranching */
SCIP_RETCODE SCIPlpiExactEndStrongbranch(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   assert( ! lpi->spx->preStrongbranchingBasisFreed() );
   lpi->spx->restorePreStrongbranchingBasis();
   lpi->spx->freePreStrongbranchingBasis();

   return SCIP_OKAY;
}


/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiExactWasSolved(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   return lpi->solved;
}

/** gets information about primal and dual feasibility of the current LP solution
 *
 *  The feasibility information is with respect to the last solving call and it is only relevant if SCIPlpiWasSolved()
 *  returns true. If the LP is changed, this information might be invalidated.
 *
 *  Note that @a primalfeasible and @dualfeasible should only return true if the solver has proved the respective LP to
 *  be feasible. Thus, the return values should be equal to the values of SCIPlpiIsPrimalFeasible() and
 *  SCIPlpiIsDualFeasible(), respectively. Note that if feasibility cannot be proved, they should return false (even if
 *  the problem might actually be feasible).
 */
SCIP_RETCODE SCIPlpiExactGetSolFeasibility(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< pointer to store primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< pointer to store dual feasibility status */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactGetSolFeasibility()\n");

   assert(lpi != NULL);
   assert(primalfeasible != NULL);
   assert(dualfeasible != NULL);

   *primalfeasible = SCIPlpiExactIsPrimalFeasible(lpi);
   *dualfeasible = SCIPlpiExactIsDualFeasible(lpi);

   return SCIP_OKAY;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExactExistsPrimalRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactExistsPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::UNBOUNDED);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExactHasPrimalRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return lpi->spx->hasPrimalRay();
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiExactIsPrimalUnbounded(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   assert(lpi->spx->status() != SPxSolver::UNBOUNDED || lpi->spx->basisStatus() == SPxBasis::UNBOUNDED);

   /* if SoPlex returns unbounded, this may only mean that an unbounded ray is available, not necessarily a primal
    * feasible point; hence we have to check the perturbation
    */
   return lpi->spx->status() == SPxSolver::UNBOUNDED;
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiExactIsPrimalInfeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::INFEASIBLE);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiExactIsPrimalFeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SPxBasis::SPxStatus basestatus;

   SCIPdebugMessage("calling SCIPlpiIsPrimalFeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   basestatus = lpi->spx->basisStatus();

   /* note that the solver status may be ABORT_VALUE and the basis status optimal; if we are optimal, isPerturbed() may
    * still return true as long as perturbation plus violation is within tolerances
    */
   assert(basestatus == SPxBasis::OPTIMAL || lpi->spx->status() != SPxSolver::OPTIMAL);

   return basestatus == SPxBasis::OPTIMAL || basestatus == SPxBasis::PRIMAL;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExactExistsDualRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactExistsDualRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::INFEASIBLE);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExactHasDualRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasDualRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return lpi->spx->hasDualFarkas();
}

/** returns TRUE iff LP is dual unbounded */
SCIP_Bool SCIPlpiExactIsDualUnbounded(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return lpi->spx->status() == SPxSolver::INFEASIBLE && lpi->spx->basisStatus() == SPxBasis::DUAL;
}

/** returns TRUE iff LP is dual infeasible */
SCIP_Bool SCIPlpiExactIsDualInfeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::UNBOUNDED);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiExactIsDualFeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualFeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   /* note that the solver status may be ABORT_VALUE and the basis status optimal; if we are optimal, isPerturbed() may
    * still return true as long as perturbation plus violation is within tolerances
    */
   assert(lpi->spx->basisStatus() == SPxBasis::OPTIMAL || lpi->spx->status() != SPxSolver::OPTIMAL);

   return (lpi->spx->basisStatus() == SPxBasis::OPTIMAL) || lpi->spx->basisStatus() == SPxBasis::DUAL;
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiExactIsOptimal(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsOptimal()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert((lpi->spx->basisStatus() == SPxBasis::OPTIMAL)
      == (SCIPlpiExactIsPrimalFeasible(lpi) && SCIPlpiExactIsDualFeasible(lpi)));

   /* note that the solver status may be ABORT_VALUE and the basis status optimal; if we are optimal, isPerturbed() may
    * still return true as long as perturbation plus violation is within tolerances
    */
   return (lpi->spx->basisStatus() == SPxBasis::OPTIMAL);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiExactIsObjlimExc(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsObjlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::ABORT_VALUE);
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiExactIsIterlimExc(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsIterlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::ABORT_ITER);
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiExactIsTimelimExc(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::ABORT_TIME);
}

/** returns the internal solution status of the solver */
int SCIPlpiExactGetInternalStatus(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactGetInternalStatus()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return static_cast<int>(lpi->spx->status());
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiExactIgnoreInstability(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiIgnoreInstability()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(success != NULL);

#if SOPLEX_APIVERSION >= 4
   *success = lpi->spx->ignoreUnscaledViolations();
#else
   *success = FALSE;
#endif

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiExactGetObjval(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Rational*        objval              /**< stores the objective value */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactGetObjval()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(objval != NULL);

   RsetSpxR(lpi, objval, lpi->spx->objValueRational());

   return SCIP_OKAY;
}


/** gets primal and dual solution vectors for feasible LPs
 *
 *  Before calling this function, the caller must ensure that the LP has been solved to optimality, i.e., that
 *  SCIPlpiIsOptimal() returns true.
 */
SCIP_RETCODE SCIPlpiExactGetSol(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Rational*        objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Rational**       primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_Rational**       dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_Rational**       activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_Rational**       redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{
   DVectorRational* tmpvec;
   SCIPdebugMessage("calling SCIPlpiExactGetSol()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   tmpvec = new DVectorRational(MAX(lpi->spx->numColsRational(), lpi->spx->numRowsRational()));
   tmpvec->clear();

   if( objval != NULL )
      RsetSpxR(lpi, objval, lpi->spx->objValueRational());

   try
   {
      if( primsol != NULL )
      {
         VectorRational tmp(lpi->spx->numColsRational(), tmpvec->get_ptr());
         (void)lpi->spx->getPrimalRational(tmp);
         RsetSpxVector(lpi, primsol, tmp);
      }
      if( dualsol != NULL )
      {
         VectorRational tmp(lpi->spx->numRowsRational(), tmpvec->get_ptr());
         (void)lpi->spx->getDualRational(tmp);
         RsetSpxVector(lpi, dualsol, tmp);
      }
      if( activity != NULL )
      {
         VectorRational tmp(lpi->spx->numRowsRational(), tmpvec->get_ptr());
         (void)lpi->spx->getSlacksRational(tmp);  /* in SoPlex, the activities are called "slacks" */
         RsetSpxVector(lpi, activity, tmp);
      }
      if( redcost != NULL )
      {
         VectorRational tmp(lpi->spx->numColsRational(), tmpvec->get_ptr());
         (void)lpi->spx->getRedCostRational(tmp);
         RsetSpxVector(lpi, redcost, tmp);
      }

      delete tmpvec;
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}


/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiExactGetPrimalRay(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Rational**       ray                 /**< primal ray */
   )
{  /*lint --e{715}*/
   DVectorRational* tmpvec;
   SCIPdebugMessage("calling SCIPlpiExactGetPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpi->spx->hasPrimalRay());
   assert(ray != NULL);

   tmpvec = new DVectorRational(MAX(lpi->spx->numColsRational(), lpi->spx->numRowsRational()));

   try
   {
      VectorRational tmp(lpi->spx->numColsRational(), tmpvec->get_ptr());
      (void)lpi->spx->getPrimalRayRational(tmp);
      RsetSpxVector(lpi, ray, tmp);
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   delete tmpvec;

   return SCIP_OKAY;
}

/** gets dual farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiExactGetDualfarkas(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Rational**       dualfarkas          /**< dual farkas row multipliers */
   )
{
   DVectorRational* tmpvec;
   SCIPdebugMessage("calling SCIPlpiExactGetDualfarkas()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpi->spx->hasDualFarkas());
   assert(dualfarkas != NULL);

   tmpvec = new DVectorRational(MAX(lpi->spx->numColsRational(), lpi->spx->numRowsRational()));

   try
   {
      VectorRational tmp(lpi->spx->numRowsRational(), tmpvec->get_ptr());
      (void)lpi->spx->getDualFarkasRational(tmp);
      RsetSpxVector(lpi, dualfarkas, tmp);
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }
   delete tmpvec;

   return SCIP_OKAY;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiExactGetIterations(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactGetIterations()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(iterations != NULL);

   *iterations = lpi->spx->numIterations();

   return SCIP_OKAY;
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
{
   int i;

   SCIPdebugMessage("calling SCIPlpiExactGetBase()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   if( rstat != NULL )
   {
      for( i = 0; i < lpi->spx->numRowsRational(); ++i )
      {
         switch( lpi->spx->basisRowStatus(i) )
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
         case SPxSolver::UNDEFINED:
         default:
            SCIPerrorMessage("invalid basis status\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
         }
      }
   }

   if( cstat != NULL )
   {
      for( i = 0; i < lpi->spx->numColsRational(); ++i )
      {
//         SCIP_Real val = 0.0;
         switch( lpi->spx->basisColStatus(i) )
         {
         case SPxSolver::BASIC:
            cstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/
            break;
         case SPxSolver::FIXED:
         /* Get reduced cost estimation. If the estimation is not correct this should not hurt:
         * If the basis is loaded into SoPlex again, the status is converted to FIXED again; in
         * this case there is no problem at all. If the basis is saved and/or used in some other
         * solver, it usually is very cheap to perform the pivots necessary to get an optimal
         * basis.
         * @todo implement getRedCostEst()
         * */
//          SCIP_CALL( getRedCostEst(lpi->spx, i, &val) );
//            if( val < 0.0 )  /* reduced costs < 0 => UPPER  else => LOWER */
//               cstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
//            else
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
         case SPxSolver::UNDEFINED:
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
SCIP_RETCODE SCIPlpiExactSetBase(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   )
{
   int i;
   int ncols;
   int nrows;

   SCIPdebugMessage("calling SCIPlpiExactSetBase()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   SCIP_CALL( SCIPlpiExactGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiExactGetNCols(lpi, &ncols) );

   assert(cstat != NULL || ncols == 0);
   assert(rstat != NULL || nrows == 0);

   assert( lpi->spx->preStrongbranchingBasisFreed() );
   invalidateSolution(lpi);

   DataArray<SPxSolver::VarStatus>& _colstat = lpi->spx->colStat();
   DataArray<SPxSolver::VarStatus>& _rowstat = lpi->spx->rowStat();

   _colstat.reSize(ncols);
   _rowstat.reSize(nrows);

   for( i = 0; i < nrows; ++i )
   {
      switch( rstat[i] ) /*lint !e613*/
      {
      case SCIP_BASESTAT_LOWER:
         _rowstat[i] = SPxSolver::ON_LOWER;
         break;
      case SCIP_BASESTAT_BASIC:
         _rowstat[i] = SPxSolver::BASIC;
         break;
      case SCIP_BASESTAT_UPPER:
         _rowstat[i] = SPxSolver::ON_UPPER;
         break;
      case SCIP_BASESTAT_ZERO:
         SCIPerrorMessage("slack variable has basis status ZERO (should not occur)\n");
         return SCIP_LPERROR; /*lint !e429*/
      default:
         SCIPerrorMessage("invalid basis status\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   for( i = 0; i < ncols; ++i )
   {
      switch( cstat[i] ) /*lint !e613*/
      {
      case SCIP_BASESTAT_LOWER:
         _colstat[i] = SPxSolver::ON_LOWER;
         break;
      case SCIP_BASESTAT_BASIC:
         _colstat[i] = SPxSolver::BASIC;
         break;
      case SCIP_BASESTAT_UPPER:
         _colstat[i] = SPxSolver::ON_UPPER;
         break;
      case SCIP_BASESTAT_ZERO:
         _colstat[i] = SPxSolver::ZERO;
         break;
      default:
         SCIPerrorMessage("invalid basis status\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->setBasis(_rowstat.get_ptr(), _colstat.get_ptr()) );
   lpi->spx->freePreStrongbranchingBasis();

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
SCIP_RETCODE SCIPlpiExactGetBasisInd(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  bind                /**< pointer to store basis indices ready to keep number of rows entries */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactGetBasisInd()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(bind != NULL);

   assert(lpi->spx->preStrongbranchingBasisFreed());

   lpi->spx->getBasisInd(bind);

   return SCIP_OKAY;
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
   SCIP_Rational**       coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   int i;
   SSVectorBase<soplex::Rational> tmpvec(0);
   SCIPdebugMessage("calling SCIPlpiExactGetBInvRow()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpi->spx->preStrongbranchingBasisFreed());
   assert(coef != NULL);

   assert(r >= 0);
   assert(r < lpi->spx->numRowsRational());

   if( !lpi->spx->getBasisInverseRowRational(r, tmpvec) )
      return SCIP_LPERROR;

   for( i = 0; i < tmpvec.size(); ++i )
   {
      inds[i] = tmpvec.index(i);
      RsetSpxR(lpi, coef[i], tmpvec.value(i));
   }

   *ninds = tmpvec.size();

   return SCIP_OKAY;
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
   SCIP_Rational**       coef,               /**< pointer to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   int i;
   SSVectorRational tmpvec(0);

   SCIPdebugMessage("calling SCIPlpiExactGetBInvCol()\n");

   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( lpi->spx->preStrongbranchingBasisFreed() );
   assert(coef != NULL);

   if( ! lpi->spx->getBasisInverseColRational(c, tmpvec) )
      return SCIP_LPERROR;

   for( i = 0; i < tmpvec.size(); ++i )
   {
      inds[i] = tmpvec.index(i);
      RsetSpxR(lpi, coef[i], tmpvec.value(i));
   }

   *ninds = tmpvec.size();


   return SCIP_OKAY;
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
{
   int ncols;
   int nrows;

   SCIPdebugMessage("calling SCIPlpiExactGetState()\n");

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpistate != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   ncols = lpi->spx->numColsRational();
   nrows = lpi->spx->numRowsRational();
   assert(ncols >= 0);
   assert(nrows >= 0);

   /* allocate lpistate data */
   SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows) );

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, ncols) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   /* get unpacked basis information */
   SCIP_CALL( SCIPlpiExactGetBase(lpi, lpi->cstat, lpi->rstat) );

   /* pack LPi state data */
   (*lpistate)->ncols = ncols;
   (*lpistate)->nrows = nrows;
   lpistatePack(*lpistate, lpi->cstat, lpi->rstat);

   return SCIP_OKAY;
}

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiExactGetState()
 */
SCIP_RETCODE SCIPlpiExactSetState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information), or NULL */
   )
{
   int lpncols;
   int lpnrows;
   int i;

   SCIPdebugMessage("calling SCIPlpiExactSetState()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpistate != NULL);
   /* assert(blkmem != NULL); */

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   lpncols = lpi->spx->numColsRational();
   lpnrows = lpi->spx->numRowsRational();
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
      SCIP_Real bnd = lpi->spx->lowerReal(i);
      if ( SCIPlpiExactIsInfinity(lpi, REALABS(bnd)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         bnd = lpi->spx->lowerReal(i);
         if ( SCIPlpiExactIsInfinity(lpi, REALABS(bnd)) )
            /* variable is free */
            lpi->cstat[i] = SCIP_BASESTAT_ZERO;  /*lint !e641*/
         else
            /* use finite upper bound */
            lpi->cstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
      }
      else
         /* use finite lower bound */
         lpi->cstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
   }
   for( i = lpistate->nrows; i < lpnrows; ++i )
      lpi->rstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/

   /* load basis information */
   SCIP_CALL( SCIPlpiExactSetBase(lpi, lpi->cstat, lpi->rstat) );

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiExactClearState(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiClearState()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   try
   {
      lpi->spx->clearBasis();
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      assert( lpi->spx->status() != SPxSolver::OPTIMAL );
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiExactFreeState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiFreeState()\n");

   assert(lpi != NULL);
   assert(lpistate != NULL);
   assert(blkmem != NULL);

   if ( *lpistate != NULL )
      lpistateFree(lpistate, blkmem);

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiExactHasStateBasis(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information), or NULL */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   return TRUE;
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiExactReadState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadState()\n");
   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(fname != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   bool success;
   SOPLEX_TRY( lpi->messagehdlr, success = lpi->spx->readBasisFile(fname, 0, 0) );

   return success ? SCIP_OKAY : SCIP_LPERROR;
}

/** writes LPi state (i.e. basis information) to a file */
SCIP_RETCODE SCIPlpiExactWriteState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(fname != NULL);
   SCIPdebugMessage("calling SCIPlpiWriteState()\n");

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   bool res;
   SOPLEX_TRY( lpi->messagehdlr, res = lpi->spx->writeBasisFile(fname, 0, 0) );

   if ( ! res )
      return SCIP_LPERROR;

   return SCIP_OKAY;
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
{
   int scaleparam;

   SCIPdebugMessage("calling SCIPlpiExactGetIntpar()\n");

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
      *ival = lpi->spx->intParam(SoPlex::ITERLIMIT);
      if( *ival == -1 )
          *ival = INT_MAX;
      break;
   case SCIP_LPPAR_PRESOLVING:
      *ival = lpi->spx->intParam(SoPlex::SIMPLIFIER) == SoPlex::SIMPLIFIER_AUTO;
      break;
   case SCIP_LPPAR_PRICING:
      *ival = (int) lpi->pricing;
      break;
   case SCIP_LPPAR_SCALING:
      scaleparam = lpi->spx->intParam(SoPlex::SCALER);

      if( scaleparam == SoPlex::SCALER_OFF )
         *ival = 0;
      else if( scaleparam == SoPlex::SCALER_BIEQUI )
         *ival = 1;
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 2)
      else
      {
         assert(scaleparam == SoPlex::SCALER_LEASTSQ);
         *ival = 2;
      }
#else
      else
      {
         assert(scaleparam == SoPlex::SCALER_GEO8);
         *ival = 2;
      }
#endif
      break;
#if SOPLEX_VERSION >= 201
   case SCIP_LPPAR_TIMING:
      *ival = (int) (lpi->spx->intParam(SoPlex::TIMER));
      break;
#endif
#if SOPLEX_VERSION >= 230 || (SOPLEX_VERSION == 220 && SOPLEX_SUBVERSION >= 3)
   case SCIP_LPPAR_RANDOMSEED:
      *ival = (int) lpi->spx->randomSeed();
      break;
#endif
#if SOPLEX_APIVERSION >= 1
   case SCIP_LPPAR_REFACTOR:
      *ival = (int) lpi->spx->intParam(SoPlex::FACTOR_UPDATE_MAX);
      break;
#endif
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiExactSetIntpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactSetIntpar()\n");

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
      assert( ival >= 0 );
      /* -1 <= ival, -1 meaning no time limit, 0 stopping immediately */
      if( ival >= INT_MAX )
         ival = -1;
      (void) lpi->spx->setIntParam(SoPlex::ITERLIMIT, ival);
      break;
   case SCIP_LPPAR_PRESOLVING:
      assert(ival == TRUE || ival == FALSE);
      (void) lpi->spx->setIntParam(SoPlex::SIMPLIFIER, (ival ? SoPlex::SIMPLIFIER_AUTO : SoPlex::SIMPLIFIER_OFF));
      break;
   case SCIP_LPPAR_PRICING:
      lpi->pricing = (SCIP_PRICING)ival;
      switch( lpi->pricing )
      {
      case SCIP_PRICING_LPIDEFAULT:
      case SCIP_PRICING_AUTO:
         (void) lpi->spx->setIntParam(SoPlex::PRICER, SoPlex::PRICER_AUTO);
         break;
      case SCIP_PRICING_FULL:
         (void) lpi->spx->setIntParam(SoPlex::PRICER, SoPlex::PRICER_STEEP);
         break;
      case SCIP_PRICING_PARTIAL:
         (void) lpi->spx->setIntParam(SoPlex::PRICER, SoPlex::PRICER_PARMULT);
         break;
      case SCIP_PRICING_STEEP:
         (void) lpi->spx->setIntParam(SoPlex::PRICER, SoPlex::PRICER_STEEP);
         break;
      case SCIP_PRICING_STEEPQSTART:
         (void) lpi->spx->setIntParam(SoPlex::PRICER, SoPlex::PRICER_QUICKSTEEP);
         break;
      case SCIP_PRICING_DEVEX:
         (void) lpi->spx->setIntParam(SoPlex::PRICER, SoPlex::PRICER_DEVEX);
         break;
      default:
         return SCIP_LPERROR;
      }
      break;
   case SCIP_LPPAR_SCALING:
      assert(ival >= 0 && ival <= 2);
      if( ival == 0 )
         (void) lpi->spx->setIntParam(SoPlex::SCALER, SoPlex::SCALER_OFF);
      else if( ival == 1 )
         (void) lpi->spx->setIntParam(SoPlex::SCALER, SoPlex::SCALER_BIEQUI);
      else
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 2)
         (void) lpi->spx->setIntParam(SoPlex::SCALER, SoPlex::SCALER_LEASTSQ);
#else
         (void) lpi->spx->setIntParam(SoPlex::SCALER, SoPlex::SCALER_GEO8);
#endif

      break;
#if SOPLEX_VERSION >= 201
   case SCIP_LPPAR_TIMING:
      assert(ival >= 0 && ival < 3);
      (void) lpi->spx->setIntParam(SoPlex::TIMER, ival);
      break;
#endif
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 3)
   case SCIP_LPPAR_RANDOMSEED:
      lpi->spx->setRandomSeed((unsigned long)(long)ival);
      break;
#endif
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION >= 221 && SOPLEX_SUBVERSION >= 3)
   case SCIP_LPPAR_POLISHING:
      assert(ival >= 0 && ival < 3);
      (void) lpi->spx->setIntParam(SoPlex::SOLUTION_POLISHING, ival);
      break;
#endif
#if SOPLEX_APIVERSION >= 1
   case SCIP_LPPAR_REFACTOR:
      assert(ival >= 0);
      (void) lpi->spx->setIntParam(SoPlex::FACTOR_UPDATE_MAX, ival);
      break;
#endif

   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiExactGetRealpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactGetRealpar()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(dval != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
       *dval = 0.0; //lpi->spx->feastol();
      break;
   case SCIP_LPPAR_DUALFEASTOL:
       *dval = 0.0; //lpi->spx->opttol();
      break;
   case SCIP_LPPAR_OBJLIM:
      if ( lpi->spx->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
         *dval = lpi->spx->realParam(SoPlex::OBJLIMIT_UPPER);
      else
         *dval = lpi->spx->realParam(SoPlex::OBJLIMIT_LOWER);
      break;
   case SCIP_LPPAR_LPTILIM:
      *dval = lpi->spx->realParam(SoPlex::TIMELIMIT);
      break;
   case SCIP_LPPAR_ROWREPSWITCH:
      *dval = lpi->spx->realParam(SoPlex::REPRESENTATION_SWITCH);
      if( *dval >= SCIPlpiExactInfinity(lpi) )
         *dval = -1.0;
      break;
   case SCIP_LPPAR_CONDITIONLIMIT:
      *dval = lpi->conditionlimit;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiExactSetRealpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiExactSetRealpar()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   switch( type )
   {
   case SCIP_LPPAR_OBJLIM:
      /* no restrictions on dval */
      if ( lpi->spx->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
         (void) lpi->spx->setRealParam(SoPlex::OBJLIMIT_UPPER, dval);
      else
         (void) lpi->spx->setRealParam(SoPlex::OBJLIMIT_LOWER, dval);
      break;
   case SCIP_LPPAR_LPTILIM:
      assert( dval > 0.0 );
      /* soplex requires 0 < dval < DEFAULT_INFINITY (= 1e100), -1 means unlimited */
      (void) lpi->spx->setRealParam(SoPlex::TIMELIMIT, dval);
      break;
   case SCIP_LPPAR_ROWREPSWITCH:
      /* 0 <= dval <= inf */
      assert( dval >= 0.0 || dval == -1.0 );
      if( dval == -1 )
         (void) lpi->spx->setRealParam(SoPlex::REPRESENTATION_SWITCH, SCIPlpiExactInfinity(lpi));
      else
         (void) lpi->spx->setRealParam(SoPlex::REPRESENTATION_SWITCH, dval);
      break;
   case SCIP_LPPAR_CONDITIONLIMIT:
      lpi->conditionlimit = dval;
      lpi->checkcondition = (dval >= 0.0);
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
SCIP_Real SCIPlpiExactInfinity(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   SCIPdebugMessage("calling SCIPlpiInfinity()\n");

   return lpi->spx->realParam(SoPlex::INFTY);
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiExactIsInfinity(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< the value */
   )
{
   assert(lpi != NULL);
   SCIPdebugMessage("calling SCIPlpiExactIsInfinity()\n");

   return (val >= lpi->spx->realParam(SoPlex::INFTY));
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
SCIP_RETCODE SCIPlpiExactReadLP(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(fname != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   if( !fileExists(fname) )
      return SCIP_NOFILE;

   try
   {
      assert(lpi->spx->intParam(SoPlex::READMODE) == SoPlex::READMODE_RATIONAL);
      if( !lpi->spx->readFile(fname) )
         return SCIP_READERROR;
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiExactWriteLP(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiWriteLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(fname != NULL);

   try
   {
      (void) lpi->spx->writeFileRational(fname);
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_WRITEERROR;
   }

   return SCIP_OKAY;
}

/**@} */
