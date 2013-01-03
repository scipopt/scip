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

/**@file   lpi_xprs.c
 * @ingroup LPIS
 * @brief  LP interface for Xpress-MP 16-21
 * @author Tobias Achterberg
 * @author Michael Perregaard
 * @author Livio Bertacco
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include <assert.h>

#include "xprs.h"

#if (XPVERSION < 21)
#define OLDRAYCODE 1 /* New ray functions available and public since version 21 */
#if (XPVERSION < 17) /* XPRSpostsolve public since v17 */
int XPRS_CC XPRSpostsolve( XPRSprob prob );
#endif
#if (XPVERSION >= 18) /* XPRSstrongbranch available since version 18, public since v21 */
int XPRS_CC XPRSstrongbranch( XPRSprob prob, const int _nbnd, const int *_mbndind, const char *_cbndtype,
   const double *_dbndval, const int _itrlimit, double *_dsbobjval, int *_msbstatus );
#endif
#endif
#ifndef XPRS_LPQUICKPRESOLVE
#define XPRS_LPQUICKPRESOLVE 8207
#endif

/* For SCIP we need an extra LP status which is optimal with */
/* scaled infeasibilities. */
#define XPRS_LP_OPTIMAL_SCALEDINFEAS 16

#include "scip/bitencode.h"
#include "scip/lpi.h"
#include "scip/pub_message.h"


/** output Xpress error */
static
void xprs_error(
   XPRSprob              prob,               /**< Xpress problem instance */
   int                   restat,             /**< return status */
   const char**          msg,                /**< error message on output */
   char*                 errmsg              /**< string to store last Xpress error */
   )
{
   *errmsg='\0';
   if (prob)
      XPRSgetlasterror(prob, errmsg);
   if (*errmsg)
      *msg = "LP Error: Xpress returned %d - %s\n";
   else
      *msg = "LP Error: Xpress returned %d\n";
}

#define CHECK_ZEROE(p, x) {  int restat = (x);  \
      if( restat != 0 ) {                       \
         char errmsg[512];                      \
         const char *msg;                       \
         xprs_error((p), restat, &msg, errmsg); \
         SCIPerrorMessage(msg, restat, errmsg); \
         return SCIP_LPERROR;                   \
      }                                         \
   }

#define CHECK_ZEROW(p, messagehdlr, x) {  int restat = (x);      \
      if( restat != 0 ) {                               \
         char errmsg[512];                              \
         const char *msg;                               \
         xprs_error((p), restat, &msg, errmsg);         \
         SCIPmessagePrintWarning((messagehdlr), msg, restat, errmsg);   \
      }                                                 \
   }

#define CHECK_ZEROLPIE(x) CHECK_ZEROE(lpi->xprslp, x)
#define CHECK_ZEROLPIW(x) CHECK_ZEROW(lpi->xprslp, lpi->messagehdlr, x)
#define CHECK_ZEROPLPIE(x) CHECK_ZEROE((*lpi)->xprslp, x)
#define CHECK_ZERO CHECK_ZEROLPIE


typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE

/** LP interface */
struct SCIP_LPi
{
   XPRSprob              xprslp;             /**< Xpress LP pointer */
   char                  name[64];           /**< problem name */
   int                   objsense;           /**< direction of optimization: +1 minim, -1 maxim */
   int                   solstat;            /**< solution status of last optimization call */
   int                   unbvec;             /**< primal or dual vector on which the problem is unbounded */
   char                  solmethod;          /**< method used to solve the LP */
   char*                 larray;             /**< array with 'L' entries for changing lower bounds */
   char*                 uarray;             /**< array with 'U' entries for changing upper bounds */
   char*                 senarray;           /**< array for storing row senses */
   SCIP_Real*            rhsarray;           /**< array for storing rhs values */
   SCIP_Real*            rngarray;           /**< array for storing range values */
   SCIP_Real*            valarray;           /**< array for storing coefficient values */
   int*                  cstat;              /**< array for storing column basis status */
   int*                  rstat;              /**< array for storing row basis status */
   int*                  indarray;           /**< array for storing coefficient indices */
   int                   boundchgsize;       /**< size of larray and uarray */
   int                   sidechgsize;        /**< size of senarray and rngarray */
   int                   valsize;            /**< size of valarray and indarray */
   int                   cstatsize;          /**< size of cstat array */
   int                   rstatsize;          /**< size of rstat array */
   int                   iterations;         /**< number of iterations used in the last solving call */
   SCIP_Bool             solisbasic;         /**< is current LP solution a basic solution? */
   SCIP_Real             par_lobjlim;        /**< objective lower bound */
   SCIP_Real             par_uobjlim;        /**< objective upper bound */
   int                   par_fastlp;         /**< special meta parameter for making LP reoptimize go faster */
   int                   par_presolve;       /**< need to distinguish between the users setting and the optimizer setting of presolve */
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

static int               numlp = 0;          /**< number of open LP objects */


/*
 * dynamic memory arrays
 */

/** resizes larray and uarray to have at least num entries */
static
SCIP_RETCODE ensureBoundchgMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->boundchgsize )
   {
      int newsize;
      int i;

      newsize = MAX(2*lpi->boundchgsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->larray, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->uarray, newsize) );
      for( i = lpi->boundchgsize; i < newsize; ++i )
      {
         lpi->larray[i] = 'L';
         lpi->uarray[i] = 'U';
      }
      lpi->boundchgsize = newsize;
   }
   assert(num <= lpi->boundchgsize);

   return SCIP_OKAY;
}

/** resizes senarray and rngarray to have at least num entries */
static
SCIP_RETCODE ensureSidechgMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->sidechgsize )
   {
      int newsize;

      newsize = MAX(2*lpi->sidechgsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->senarray, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->rhsarray, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->rngarray, newsize) );
      lpi->sidechgsize = newsize;
   }
   assert(num <= lpi->sidechgsize);

   return SCIP_OKAY;
}

/** resizes valarray and indarray to have at least num entries */
static
SCIP_RETCODE ensureValMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->valsize )
   {
      int newsize;

      newsize = MAX(2*lpi->valsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->valarray, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->indarray, newsize) );
      lpi->valsize = newsize;
   }
   assert(num <= lpi->valsize);

   return SCIP_OKAY;
}

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

/** stores current basis in internal arrays of LPI data structure */
static
SCIP_RETCODE getBase(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int ncols;
   int nrows;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("getBase()\n");

   XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols);
   XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows);

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, ncols) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   /* get unpacked basis information from Xpress */
   CHECK_ZERO( XPRSgetbasis(lpi->xprslp, lpi->rstat, lpi->cstat) );

   return SCIP_OKAY;
}

/** loads basis stored in internal arrays of LPI data structure into Xpress */
static
SCIP_RETCODE setBase(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   SCIPdebugMessage("setBase()\n");

   /* load basis information into Xpress */
   CHECK_ZERO( XPRSloadbasis(lpi->xprslp, lpi->rstat, lpi->cstat) );

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

/** marks the current LP to be unsolved */
static
void invalidateSolution(SCIP_LPI* lpi)
{
   assert(lpi != NULL);
   lpi->solstat = -1;
}

/** converts SCIP's objective sense into Xpress' objective sense */
static
int xprsObjsen(SCIP_OBJSEN objsen)
{
   switch( objsen )
   {
   case SCIP_OBJSEN_MAXIMIZE:
      return -1;
   case SCIP_OBJSEN_MINIMIZE:
      return +1;
   default:
      SCIPerrorMessage("invalid objective sense\n");
      SCIPABORT();
      return 0;
   }
}

/** converts SCIP's lhs/rhs pairs into Xpress' sen/rhs/rng */
static
void convertSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand side vector */
   const SCIP_Real*      rhs                 /**< right hand side vector */
   )
{
   int i;

   assert(lpi != NULL);
   assert(nrows >= 0);
   assert(lhs != NULL);
   assert(rhs != NULL);

   /* convert lhs/rhs into sen/rhs/rng */
   for( i = 0; i < nrows; ++i )
   {
      assert(lhs[i] <= rhs[i]);
      if( lhs[i] == rhs[i] ) /*lint !e777*/
      {
         assert(XPRS_MINUSINFINITY < rhs[i] && rhs[i] < XPRS_PLUSINFINITY);
         lpi->senarray[i] = 'E';
         lpi->rhsarray[i] = rhs[i];
         lpi->rngarray[i] = 0.0;
      }
      else if( lhs[i] <= XPRS_MINUSINFINITY )
      {
         assert(XPRS_MINUSINFINITY < rhs[i] && rhs[i] < XPRS_PLUSINFINITY);
         lpi->senarray[i] = 'L';
         lpi->rhsarray[i] = rhs[i];
         lpi->rngarray[i] = 0.0;
      }
      else if( rhs[i] >= XPRS_PLUSINFINITY )
      {
         assert(XPRS_MINUSINFINITY < lhs[i] && lhs[i] < XPRS_PLUSINFINITY);
         lpi->senarray[i] = 'G';
         lpi->rhsarray[i] = lhs[i];
         lpi->rngarray[i] = 0.0;
      }
      else
      {
         /* Xpress defines a ranged row to be within rhs-rng and rhs.
          */
         lpi->senarray[i] = 'R';
         lpi->rhsarray[i] = rhs[i];
         lpi->rngarray[i] = rhs[i] - lhs[i];
      }
   }
}

/** converts Xpress' sen/rhs/rng triplets into SCIP's lhs/rhs pairs */
static
void reconvertBothSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows */
   SCIP_Real*            lhs,                /**< buffer to store the left hand side vector */
   SCIP_Real*            rhs                 /**< buffer to store the right hand side vector */
   )
{
   int i;

   assert(lpi != NULL);
   assert(nrows >= 0);
   assert(lhs != NULL);
   assert(rhs != NULL);

   for( i = 0; i < nrows; ++i )
   {
      switch( lpi->senarray[i] )
      {
      case 'E':
         lhs[i] = lpi->rhsarray[i];
         rhs[i] = lpi->rhsarray[i];
         break;

      case 'L':
         lhs[i] = XPRS_MINUSINFINITY;
         rhs[i] = lpi->rhsarray[i];
         break;

      case 'G':
         lhs[i] = lpi->rhsarray[i];
         rhs[i] = XPRS_PLUSINFINITY;
         break;

      case 'R':
         assert(lpi->rngarray[i] >= 0.0);
         rhs[i] = lpi->rhsarray[i];
         lhs[i] = lpi->rhsarray[i] - lpi->rngarray[i];
         break;

      default:
         SCIPerrorMessage("invalid row sense\n");
         SCIPABORT();
      }
      assert(lhs[i] <= rhs[i]);
   }
}

/** converts Xpress' sen/rhs/rng triplets into SCIP's lhs/rhs pairs, only storing the left hand side */
static
void reconvertLhs(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows */
   SCIP_Real*            lhs                 /**< buffer to store the left hand side vector */
   )
{
   int i;

   assert(lpi != NULL);
   assert(nrows >= 0);
   assert(lhs != NULL);

   for( i = 0; i < nrows; ++i )
   {
      switch( lpi->senarray[i] )
      {
      case 'E':
         assert(lpi->rngarray[i] == 0.0);
         lhs[i] = lpi->rhsarray[i];
         break;

      case 'L':
         assert(lpi->rngarray[i] == 0.0);
         lhs[i] = XPRS_MINUSINFINITY;
         break;

      case 'G':
         assert(lpi->rngarray[i] == 0.0);
         lhs[i] = lpi->rhsarray[i];
         break;

      case 'R':
         assert(lpi->rngarray[i] >= 0.0);
         lhs[i] = lpi->rhsarray[i] - lpi->rngarray[i];
         break;

      default:
         SCIPerrorMessage("invalid row sense\n");
         SCIPABORT();
      }
   }
}

/** converts Xpress' sen/rhs/rng triplets into SCIP's lhs/rhs pairs, only storing the right hand side */
static
void reconvertRhs(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows */
   SCIP_Real*            rhs                 /**< buffer to store the right hand side vector */
   )
{
   int i;

   assert(lpi != NULL);
   assert(nrows >= 0);
   assert(rhs != NULL);

   for( i = 0; i < nrows; ++i )
   {
      switch( lpi->senarray[i] )
      {
      case 'E':
         assert(lpi->rngarray[i] == 0.0);
         rhs[i] = lpi->rhsarray[i];
         break;

      case 'L':
         assert(lpi->rngarray[i] == 0.0);
         rhs[i] = lpi->rhsarray[i];
         break;

      case 'G':
         assert(lpi->rngarray[i] == 0.0);
         rhs[i] = XPRS_PLUSINFINITY;
         break;

      case 'R':
         assert(lpi->rngarray[i] >= 0.0);
         rhs[i] = lpi->rhsarray[i];
         break;

      default:
         SCIPerrorMessage("invalid row sense\n");
         SCIPABORT();
      }
   }
}

/** converts Xpress' sen/rhs/rng triplets into SCIP's lhs/rhs pairs */
static
void reconvertSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows */
   SCIP_Real*            lhs,                /**< buffer to store the left hand side vector, or NULL */
   SCIP_Real*            rhs                 /**< buffer to store the right hand side vector, or NULL */
   )
{
   if( lhs != NULL && rhs != NULL )
      reconvertBothSides(lpi, nrows, lhs, rhs);
   else if( lhs != NULL )
      reconvertLhs(lpi, nrows, lhs);
   else if( rhs != NULL )
      reconvertRhs(lpi, nrows, rhs);
}




/*
 * LP Interface Methods
 */


/*
 * Miscellaneous Methods
 */

static char xprsname[100];

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   sprintf(xprsname, "Xpress %d", XPVERSION);

   return xprsname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return "Linear Programming Solver developed by FICO (www.fico.com/xpress)";
}

/** gets pointer for LP solver - use only with great care
 *
 *  Here we return the pointer to the LP environment.
 */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   return (void*) lpi->xprslp;
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
   int izero = 0;

   assert(sizeof(SCIP_Real) == sizeof(double)); /* Xpress only works with doubles as floating points */
   assert(sizeof(SCIP_Bool) == sizeof(int));    /* Xpress only works with ints as bools */
   assert(lpi != NULL);
   assert(numlp >= 0);

   SCIPdebugMessage("SCIPlpiCreate()\n");

   /* Initialize the Xpress library (licensing). */
   if( numlp == 0 )
   {
      CHECK_ZEROE( NULL, XPRSinit(NULL) );
   }

   /* create LP */
   SCIP_ALLOC( BMSallocMemory(lpi) );
   assert(strlen(name) < 64);
   strcpy((*lpi)->name, name);
   (*lpi)->larray = NULL;
   (*lpi)->uarray = NULL;
   (*lpi)->senarray = NULL;
   (*lpi)->rhsarray = NULL;
   (*lpi)->rngarray = NULL;
   (*lpi)->valarray = NULL;
   (*lpi)->cstat = NULL;
   (*lpi)->rstat = NULL;
   (*lpi)->indarray = NULL;
   (*lpi)->boundchgsize = 0;
   (*lpi)->sidechgsize = 0;
   (*lpi)->valsize = 0;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->iterations = 0;
   (*lpi)->solisbasic = TRUE;
   (*lpi)->solmethod = ' ';
   (*lpi)->par_lobjlim = -1e+40;
   (*lpi)->par_uobjlim = +1e+40;
   (*lpi)->par_fastlp = 1;
   (*lpi)->par_presolve = 1;
   (*lpi)->messagehdlr = messagehdlr;

   CHECK_ZEROPLPIE( XPRScreateprob(&(*lpi)->xprslp) );
   invalidateSolution(*lpi);

   /* Turn logging off until the user explicitly turns it on. This should */
   /* prevent any unwanted Xpress output from appearing in the SCIP log. */
   CHECK_ZEROPLPIE( XPRSsetintcontrol((*lpi)->xprslp, XPRS_OUTPUTLOG, 0) );

   /* Reserve some extra space for names. */
   CHECK_ZEROPLPIE( XPRSsetintcontrol((*lpi)->xprslp, XPRS_MPSNAMELENGTH, 16) );

#ifdef XPRS_SOLUTIONFILE
   /* Don't use solution files. */
   CHECK_ZEROPLPIE( XPRSsetintcontrol((*lpi)->xprslp, XPRS_SOLUTIONFILE, 0) );
#endif

   /* We need to create an empty LP in this prob since SCIP might */
   /* attempt to add rows or columns to it. */
   CHECK_ZEROPLPIE( XPRSloadlp((*lpi)->xprslp, "temp", 0, 0, NULL, NULL, NULL, NULL, &izero, NULL, NULL, NULL, NULL, NULL) );

   /* Tell Xpress to not declare a problem infeasible in presolve. */
   CHECK_ZEROPLPIE( XPRSsetintcontrol((*lpi)->xprslp, XPRS_PRESOLVE, -1) );

   numlp++;

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

   SCIPdebugMessage("SCIPlpiFree()\n");

   /* free LP */
   CHECK_ZEROPLPIE( XPRSdestroyprob(((*lpi)->xprslp)) );

   /* free memory */
   BMSfreeMemoryArrayNull(&(*lpi)->larray);
   BMSfreeMemoryArrayNull(&(*lpi)->uarray);
   BMSfreeMemoryArrayNull(&(*lpi)->senarray);
   BMSfreeMemoryArrayNull(&(*lpi)->rhsarray);
   BMSfreeMemoryArrayNull(&(*lpi)->rngarray);
   BMSfreeMemoryArrayNull(&(*lpi)->cstat);
   BMSfreeMemoryArrayNull(&(*lpi)->rstat);
   BMSfreeMemory(lpi);

   /* free environment */
   numlp--;
   if( numlp == 0 )
   {
      XPRSfree();
   }

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
   int* cnt = NULL;
   int r, c;
   int namelength;
   int cnamesize;
   int rnamesize;
   char *cnamestore = NULL;
   char *rnamestore = NULL;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("loading LP in column format into Xpress: %d cols, %d rows\n", ncols, nrows);

   invalidateSolution(lpi);

   SCIP_CALL( ensureSidechgMem(lpi, nrows) );

   /* Save objective sense for when we have to solve the LP. */
   lpi->objsense = objsen;

   /* convert lhs/rhs into sen/rhs/range tuples */
   convertSides(lpi, nrows, lhs, rhs);

   /* get the longest name since we need to ask Xpress to make enough space before loading the LP. */
   namelength = 0;
   cnamesize = 16;
   if (colnames) 
   {
      for (c = 0; c < ncols; c++)
      {
         int isize = strlen(colnames[c]);
         cnamesize += isize+1;
         if (namelength < isize)
            namelength = isize;
      }
   }
   rnamesize = 16;
   if (rownames)
   {
      for (r = 0; r < nrows; r++) 
      {
         int isize = strlen(rownames[r]);
         rnamesize += isize+1;
         if (namelength < isize)
            namelength = isize;
      }
   }
   if (namelength) 
   {
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_MPSNAMELENGTH, namelength) );
   }

   /* calculate column lengths */
   SCIP_ALLOC( BMSallocMemoryArray(&cnt, ncols) );
   for( c = 0; c < ncols-1; ++c )
   {
      cnt[c] = beg[c+1] - beg[c];
      assert(cnt[c] >= 0);
   }
   cnt[ncols-1] = nnonz - beg[ncols-1];
   assert(cnt[ncols-1] >= 0);

   /* copy data into Xpress */
   CHECK_ZERO( XPRSloadlp(lpi->xprslp, lpi->name, ncols, nrows, lpi->senarray, lpi->rhsarray,
         lpi->rngarray, obj, beg, cnt, ind, val, lb, ub) );
   if (colnames) 
   {
      /* We need all names stored consecutively in a single array. */
      int isize = 0;
      SCIP_ALLOC( BMSallocMemoryArray(&cnamestore, cnamesize) );
      for (c = 0; c < ncols; c++) 
      {
         strcpy(cnamestore+isize, colnames[c]);
         isize += strlen(colnames[c])+1;
      }
      CHECK_ZEROLPIW( XPRSaddnames(lpi->xprslp, 2, cnamestore, 0, ncols-1) );
      BMSfreeMemoryArray(&cnamestore);
   }
   if (rownames) 
   {
      /* We need all names stored consecutively in a single array. */
      int isize = 0;
      SCIP_ALLOC( BMSallocMemoryArray(&rnamestore, rnamesize) );
      for (c = 0; c < nrows; c++) 
      {
         strcpy(rnamestore+isize, rownames[c]);
         isize += strlen(rownames[c])+1;
      }
      CHECK_ZEROLPIW( XPRSaddnames(lpi->xprslp, 1, rnamestore, 0, nrows-1) );
      BMSfreeMemoryArray(&rnamestore);
   }

   /* free temporary memory */
   BMSfreeMemoryArray(&cnt);

   {
      int chk_ncols;
      int chk_nrows;
      int chk_nnonz;
      CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &chk_ncols) );
      CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &chk_nrows) );
      CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ELEMS, &chk_nnonz) );
      assert(chk_ncols == ncols);
      assert(chk_nrows == nrows);
      assert(chk_nnonz == nnonz);
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
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   int c;
   int imaxnamelength;
   int *mstart = NULL;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("adding %d columns with %d nonzeros to Xpress\n", ncols, nnonz);

   invalidateSolution(lpi);

   /* We need ncol+1 entries in the start array for Xpress. */
   SCIP_ALLOC( BMSallocMemoryArray(&mstart, ncols+1) );
   for (c = 0; c < ncols; c++)
      mstart[c] = beg[c];
   mstart[ncols] = nnonz;
   CHECK_ZERO( XPRSaddcols(lpi->xprslp, ncols, nnonz, obj, mstart, ind, val, lb, ub) );
   BMSfreeMemoryArray(&mstart);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_NAMELENGTH, &imaxnamelength) );
   imaxnamelength *= 8;
   if (colnames) 
   {
      int lp_ncols;
      char *cnamestore;
      int cnamesize = 0;
      int isize;
      CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &lp_ncols) );
      for (c = 0; c < ncols; c++)
      {
         isize = strlen(colnames[c]);
#if (XPVERSION < 19)
         /* Xpress versions older than 19 does not allow names of arbitrary length. */
         if (isize > imaxnamelength) 
            isize = imaxnamelength;
#endif
         cnamesize += isize+1;
      }
      SCIP_ALLOC( BMSallocMemoryArray(&cnamestore, cnamesize) );
      isize = 0;
      for (c = 0; c < ncols; c++) 
      {
         int i;
         for (i = 0; colnames[c][i]; i++)
         {
#if (XPVERSION < 19)
            if (i >= imaxnamelength) 
               break;
#endif
            cnamestore[isize++] = colnames[c][i];
         }
         cnamestore[isize++] = '\0';
      }
      assert(isize == cnamesize);
      CHECK_ZEROLPIW( XPRSaddnames(lpi->xprslp, 2, cnamestore, lp_ncols-ncols, lp_ncols-1) );
      BMSfreeMemoryArray(&cnamestore);
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
   int c;
   int ncols;
   int *mind = NULL;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < ncols);

   SCIPdebugMessage("deleting %d columns from Xpress\n", lastcol - firstcol + 1);

   invalidateSolution(lpi);

   SCIP_ALLOC( BMSallocMemoryArray(&mind, lastcol-firstcol+1) );
   for (c = firstcol; c <= lastcol; c++)
      mind[c-firstcol] = c;
   CHECK_ZERO( XPRSdelcols(lpi->xprslp, lastcol-firstcol+1, mind) );
   BMSfreeMemoryArray(&mind);

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
   int c_new;
   int c;
   int ndel;
   int ncols;
   int *mind = NULL;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("deleting a column set from Xpress\n");

   invalidateSolution(lpi);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );
   ndel = 0;
   for (c = 0; c < ncols; c++) 
   {
      if (dstat[c]) 
         ndel++;
   }
   SCIP_ALLOC( BMSallocMemoryArray(&mind, ndel) );
   c_new = 0;
   ndel = 0;
   for (c = 0; c < ncols; c++)
   {
      if (dstat[c])
      {
         mind[ndel++] = c;
         dstat[c] = -1;
      } 
      else 
      {
         dstat[c] = c_new++;
      }
   }
   CHECK_ZERO( XPRSdelcols(lpi->xprslp, ndel, mind) );
   BMSfreeMemoryArray(&mind);

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
   int r;
   int lp_nrows;
   int *mstart = NULL;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("adding %d rows with %d nonzeros to Xpress\n", nrows, nnonz);

   invalidateSolution(lpi);

   SCIP_CALL( ensureSidechgMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &lp_nrows) );
   convertSides(lpi, nrows, lhs, rhs);

   SCIP_ALLOC( BMSallocMemoryArray(&mstart, nrows+1) );
   for (r = 0; r < nrows; r++)
      mstart[r] = beg[r];
   mstart[nrows] = nnonz;
   CHECK_ZERO( XPRSaddrows(lpi->xprslp, nrows, nnonz, lpi->senarray, lpi->rhsarray, lpi->rngarray, mstart, ind, val) );
   BMSfreeMemoryArray(&mstart);

   if (rownames)
   {
      char *rnamestore;
      int imaxnamelength;
      int rnamesize = 0;
      int isize;

      CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_NAMELENGTH, &imaxnamelength) );
      imaxnamelength *= 8;
      for (r = 0; r < nrows; r++)
      {
         isize = strlen(rownames[r]);
#if (XPVERSION < 19)
         if (isize > imaxnamelength)
            isize = imaxnamelength;
#endif
         rnamesize += isize+1;
      }
      SCIP_ALLOC( BMSallocMemoryArray(&rnamestore, rnamesize) );
      isize = 0;
      for (r = 0; r < nrows; r++)
      {
         int i;
         for (i = 0; rownames[r][i]; i++)
         {
#if (XPVERSION < 19)
            if (i >= imaxnamelength)
               break;
#endif
            rnamestore[isize++] = rownames[r][i];
         }
         rnamestore[isize++] = '\0';
      }
      assert(isize == rnamesize);
      CHECK_ZEROLPIW( XPRSaddnames(lpi->xprslp, 1, rnamestore, lp_nrows, lp_nrows+nrows-1) );
      BMSfreeMemoryArray(&rnamestore);
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
   int r;
   int nrows;
   int *mind = NULL;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < nrows);

   SCIPdebugMessage("deleting %d rows from Xpress\n", lastrow - firstrow + 1);

   invalidateSolution(lpi);

   SCIP_ALLOC( BMSallocMemoryArray(&mind, lastrow-firstrow+1) );
   for (r = firstrow; r <= lastrow; r++)
      mind[r-firstrow] = r;
   CHECK_ZERO( XPRSdelrows(lpi->xprslp, lastrow-firstrow+1, mind) );
   BMSfreeMemoryArray(&mind);

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
   int r, r_new;
   int ndel;
   int nrows;
   int *mind = NULL;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("deleting a row set from Xpress\n");

   invalidateSolution(lpi);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   ndel = 0;
   for (r = 0; r < nrows; r++)
   {
      if (dstat[r])
         ndel++;
   }
   SCIP_ALLOC( BMSallocMemoryArray(&mind, ndel) );
   r_new = 0;
   ndel = 0;
   for (r = 0; r < nrows; r++)
   {
      if (dstat[r])
      {
         mind[ndel++] = r;
         dstat[r] = -1;
      }
      else
      {
         dstat[r] = r_new++;
      }
   }
   CHECK_ZERO( XPRSdelrows(lpi->xprslp, ndel, mind) );
   BMSfreeMemoryArray(&mind);

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int ncols;
   int nrows;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("clearing Xpress LP\n");

   invalidateSolution(lpi);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );

   if( ncols >= 1 )
   {
      SCIP_CALL( SCIPlpiDelCols(lpi, 0, ncols-1) );
   }
   if( nrows >= 1 )
   {
      SCIP_CALL( SCIPlpiDelRows(lpi, 0, nrows-1) );
   }

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("changing %d bounds in Xpress\n", ncols);
#ifdef SCIP_DEBUG
   {
      int i;
      for( i = 0; i < ncols; ++i )
         SCIPdebugPrintf("  col %d: [%g,%g]\n", ind[i], lb[i], ub[i]);
   }
#endif

   invalidateSolution(lpi);

   SCIP_CALL( ensureBoundchgMem(lpi, ncols) );

   CHECK_ZERO( XPRSchgbounds(lpi->xprslp, ncols, ind, lpi->larray, (SCIP_Real*)lb) );
   CHECK_ZERO( XPRSchgbounds(lpi->xprslp, ncols, ind, lpi->uarray, (SCIP_Real*)ub) );

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("changing %d sides in Xpress\n", nrows);

   invalidateSolution(lpi);

   SCIP_CALL( ensureSidechgMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   convertSides(lpi, nrows, lhs, rhs);

   /* change row sides */
   CHECK_ZERO( XPRSchgrowtype(lpi->xprslp, nrows, ind, lpi->senarray) );
   CHECK_ZERO( XPRSchgrhs(lpi->xprslp, nrows, ind, lpi->rhsarray) );
   CHECK_ZERO( XPRSchgrhsrange(lpi->xprslp, nrows, ind, lpi->rngarray) );

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("changing coefficient row %d, column %d in Xpress to %g\n", row, col, newval);

   invalidateSolution(lpi);

   CHECK_ZERO( XPRSchgcoef(lpi->xprslp, row, col, newval) );

   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("changing objective sense in Xpress to %d\n", objsen);

   invalidateSolution(lpi);

   lpi->objsense = xprsObjsen(objsen);

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("changing %d objective values in Xpress\n", ncols);

   CHECK_ZERO( XPRSchgobj(lpi->xprslp, ncols, ind, obj) );

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
   int nnonz;
   int ncol;
   int beg;
   int i;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(scaleval != 0.0);

   SCIPdebugMessage("scaling row %d with factor %g in Xpress\n", row, scaleval);

   invalidateSolution(lpi);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncol) );
   SCIP_CALL( ensureValMem(lpi, ncol) );

   /* get the row */
   SCIP_CALL( SCIPlpiGetRows(lpi, row, row, &lhs, &rhs, &nnonz, &beg, lpi->indarray, lpi->valarray) );

   /* scale row coefficients */
   for( i = 0; i < nnonz; ++i )
   {
      SCIP_CALL( SCIPlpiChgCoef(lpi, row, lpi->indarray[i], lpi->valarray[i] * scaleval) );
   }

   /* scale row sides */
   if( lhs > XPRS_MINUSINFINITY )
      lhs *= scaleval;
   else if( scaleval < 0.0 )
      lhs = XPRS_PLUSINFINITY;
   if( rhs < XPRS_PLUSINFINITY )
      rhs *= scaleval;
   else if( scaleval < 0.0 )
      rhs = XPRS_MINUSINFINITY;
   if( scaleval > 0.0 )
   {
      SCIP_CALL( SCIPlpiChgSides(lpi, 1, &row, &lhs, &rhs) );
   }
   else
   {
      SCIP_CALL( SCIPlpiChgSides(lpi, 1, &row, &rhs, &lhs) );
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
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   int nnonz;
   int ncol;
   int beg;
   int i;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(scaleval != 0.0);

   SCIPdebugMessage("scaling column %d with factor %g in Xpress\n", col, scaleval);

   invalidateSolution(lpi);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncol) );
   SCIP_CALL( ensureValMem(lpi, ncol) );

   /* get the column */
   SCIP_CALL( SCIPlpiGetCols(lpi, col, col, &lb, &ub, &nnonz, &beg, lpi->indarray, lpi->valarray) );

   /** get objective coefficient */
   SCIP_CALL( SCIPlpiGetObj(lpi, col, col, &obj) );

   /* scale column coefficients */
   for( i = 0; i < nnonz; ++i )
   {
      SCIP_CALL( SCIPlpiChgCoef(lpi, lpi->indarray[i], col, lpi->valarray[i] * scaleval) );
   }

   /* scale objective value */
   obj *= scaleval;
   SCIP_CALL( SCIPlpiChgObj(lpi, 1, &col, &obj) );

   /* scale column bounds */
   if( lb > XPRS_MINUSINFINITY )
      lb /= scaleval;
   else if( scaleval < 0.0 )
      lb = XPRS_PLUSINFINITY;
   if( ub < XPRS_PLUSINFINITY )
      ub /= scaleval;
   else if( scaleval < 0.0 )
      ub = XPRS_MINUSINFINITY;
   if( scaleval > 0.0 )
   {
      SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &col, &lb, &ub) );
   }
   else
   {
      SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &col, &ub, &lb) );
   }

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(nrows != NULL);

   SCIPdebugMessage("getting number of rows\n");

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, nrows) );

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(ncols != NULL);

   SCIPdebugMessage("getting number of columns\n");

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, ncols) );

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(nnonz != NULL);

   SCIPdebugMessage("getting number of non-zeros\n");

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ELEMS, nnonz) );

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
   int ncol;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncol) );
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < ncol);

   SCIPdebugMessage("getting columns %d to %d\n", firstcol, lastcol);

   if( lb != NULL )
   {
      assert(ub != NULL);

      CHECK_ZERO( XPRSgetlb(lpi->xprslp, lb, firstcol, lastcol) );
      CHECK_ZERO( XPRSgetub(lpi->xprslp, ub, firstcol, lastcol) );
   }
   else
      assert(ub == NULL);

   if( nnonz != NULL )
   {
      int c;
      int ndim;
      int *mstart = NULL;

      assert(beg != NULL);
      assert(ind != NULL);
      assert(val != NULL);

      /* get matrix entries */
      SCIP_ALLOC( BMSallocMemoryArray(&mstart, lastcol-firstcol+2) );
      SCIPlpiGetNNonz(lpi, &ndim);
      CHECK_ZERO( XPRSgetcols(lpi->xprslp, mstart, ind, val, ndim, nnonz, firstcol, lastcol) );
      assert(*nnonz <= ndim);
      assert(mstart[lastcol-firstcol+1] == *nnonz);
      for (c = 0; c < lastcol-firstcol+1; c++)
         beg[c] = mstart[c];
      BMSfreeMemoryArray(&mstart);
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
   int nrows;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < nrows);

   SCIPdebugMessage("getting rows %d to %d\n", firstrow, lastrow);

   if( lhs != NULL || rhs != NULL )
   {
      /* get row sense, rhs, and ranges */
      SCIP_CALL( ensureSidechgMem(lpi, lastrow - firstrow + 1) );
      CHECK_ZERO( XPRSgetrowtype(lpi->xprslp, lpi->senarray, firstrow, lastrow) );
      CHECK_ZERO( XPRSgetrhs(lpi->xprslp, lpi->rhsarray, firstrow, lastrow) );
      CHECK_ZERO( XPRSgetrhsrange(lpi->xprslp, lpi->rngarray, firstrow, lastrow) );

      /* convert sen/rhs/range into lhs/rhs tuples */
      reconvertSides(lpi, lastrow - firstrow + 1, lhs, rhs);
   }

   if( nnonz != NULL )
   {
      int r;
      int ndim;
      int *mstart = NULL;

      assert(beg != NULL);
      assert(ind != NULL);
      assert(val != NULL);

      /* get matrix entries */
      SCIP_ALLOC( BMSallocMemoryArray(&mstart, lastrow-firstrow+2) );
      SCIPlpiGetNNonz(lpi, &ndim);
      CHECK_ZERO( XPRSgetrows(lpi->xprslp, mstart, ind, val, ndim, nnonz, firstrow, lastrow) );
      assert(*nnonz <= ndim);
      assert(mstart[lastrow-firstrow+1] == *nnonz);
      for (r = 0; r < lastrow-firstrow+1; r++)
         beg[r] = mstart[r];
      BMSfreeMemoryArray(&mstart);
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

/** gets the objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   SCIPerrorMessage("SCIPlpiGetObjsen() has not been implemented yet.\n");
   return SCIP_LPERROR;
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiGetObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(firstcol <= lastcol);
   assert(vals != NULL);

   SCIPdebugMessage("getting objective values %d to %d\n", firstcol, lastcol);

   CHECK_ZERO( XPRSgetobj(lpi->xprslp, vals, firstcol, lastcol) );

   return SCIP_OKAY;
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiGetBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get bounds for */
   int                   lastcol,            /**< last column to get bounds for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(firstcol <= lastcol);

   SCIPdebugMessage("getting bounds %d to %d\n", firstcol, lastcol);

   if( lbs != NULL )
   {
      CHECK_ZERO( XPRSgetlb(lpi->xprslp, lbs, firstcol, lastcol) );
   }

   if( ubs != NULL )
   {
      CHECK_ZERO( XPRSgetub(lpi->xprslp, ubs, firstcol, lastcol) );
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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(firstrow <= lastrow);

   SCIPdebugMessage("getting row sides %d to %d\n", firstrow, lastrow);

   /* get row sense, rhs, and ranges */
   SCIP_CALL( ensureSidechgMem(lpi, lastrow - firstrow + 1) );
   CHECK_ZERO( XPRSgetrowtype(lpi->xprslp, lpi->senarray, firstrow, lastrow) );
   CHECK_ZERO( XPRSgetrhs(lpi->xprslp, lpi->rhsarray, firstrow, lastrow) );
   CHECK_ZERO( XPRSgetrhsrange(lpi->xprslp, lpi->rngarray, firstrow, lastrow) );

   /* convert sen/rhs/range into lhs/rhs tuples */
   reconvertSides(lpi, lastrow - firstrow + 1, lhss, rhss);

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
   int i;
   int nnonz;
   int mstart[2];
   int *mind = NULL;
   double *dval = NULL;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("getting coefficient of row %d col %d\n", row, col);

   /* To get a coefficient we need to extract the full row or column first. */
   CHECK_ZERO( XPRSgetrows(lpi->xprslp, NULL, NULL, NULL, 0, &nnonz, row, row) );
   SCIP_ALLOC( BMSallocMemoryArray(&mind, nnonz) );
   SCIP_ALLOC( BMSallocMemoryArray(&dval, nnonz) );
   CHECK_ZERO( XPRSgetrows(lpi->xprslp, mstart, mind, dval, nnonz, &nnonz, row, row) );
   *val = 0.0;
   for (i = 0; i < nnonz; i++)
   {
      if (mind[i] == col)
      {
         *val = dval[i];
         break;
      }
   }
   BMSfreeMemoryArray(&dval);
   BMSfreeMemoryArray(&mind);

   return SCIP_OKAY;
}

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solve LP */
static SCIP_RETCODE lpiSolve(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           method              /**< indicates the method to use ('p' - primal, 'd' - dual, 'b' - barrier) */
   )
{
   int ncols;
   int nrows;
   int primalinfeasible;
   int dualinfeasible;
   int state;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );
   SCIPdebugMessage("calling Xpress lp solver type %s: %d cols, %d rows\n", method, ncols, nrows);

   invalidateSolution(lpi);

   if (lpi->par_fastlp)
   {
      /* Set controls to try and speed up the lp solve. */
      int keepbasis;
      CHECK_ZERO( XPRSgetintcontrol(lpi->xprslp, XPRS_KEEPBASIS, &keepbasis) );
      if (keepbasis || !lpi->par_presolve)
      {
         /* If we are reoptimizing from a given basis then presolve might have */
         /* quite a significant overhead. */
         CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_PRESOLVE, 0) );
         CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_LPQUICKPRESOLVE, 0) );
      }
      else
      {
         /* No given basis so presolve might reduce the problem enough to speed */
         /* up the LP solve. */
         CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_PRESOLVE, -1) );
         CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_LPQUICKPRESOLVE, 1) );
      }
      /* Don't refactorize at the end of the solve. */
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_REFACTOR, 0) );
   }
   else
   {
      /* Use default settings for solving an lp (hopefully) robustly. */
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_PRESOLVE, (lpi->par_presolve) ?  -1 : 0) );
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_REFACTOR, 1) );
   }

   if (lpi->objsense > 0)
   {
      CHECK_ZERO( XPRSsetdblcontrol(lpi->xprslp, XPRS_MIPABSCUTOFF, lpi->par_uobjlim) );
      SCIPdebugMessage("calling XPRSminim()\n");
      CHECK_ZERO( XPRSminim(lpi->xprslp, method) );
   }
   else
   {
      CHECK_ZERO( XPRSsetdblcontrol(lpi->xprslp, XPRS_MIPABSCUTOFF, lpi->par_lobjlim) );
      SCIPdebugMessage("calling XPRSmaxim()\n");
      CHECK_ZERO( XPRSmaxim(lpi->xprslp, method) );
   }

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_LPSTATUS, &lpi->solstat) );
   if (lpi->solstat == XPRS_LP_UNBOUNDED || lpi->solstat == XPRS_LP_INFEAS)
   {
      CHECK_ZERO( XPRSgetunbvec(lpi->xprslp, &lpi->unbvec) );
   }
   else
      lpi->unbvec = -1;

   /* Make sure the LP is postsolved in case it was interrupted. */
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_PRESOLVESTATE, &state) );
   if (state & (2|4))
   {
      /* Problem is in a presolve state - postsolve it. */
      CHECK_ZERO( XPRSpostsolve(lpi->xprslp) );
   }

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_SIMPLEXITER, &lpi->iterations) );
   lpi->solisbasic = TRUE;
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_PRIMALINFEAS, &primalinfeasible) );
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_DUALINFEAS, &dualinfeasible) );
   SCIPdebugMessage(" -> Xpress returned solstat=%d, pinfeas=%d, dinfeas=%d (%d iterations)\n",
      lpi->solstat, primalinfeasible, dualinfeasible, lpi->iterations);

   if ((lpi->solstat == XPRS_LP_OPTIMAL) && (primalinfeasible || dualinfeasible))
      lpi->solstat = XPRS_LP_OPTIMAL_SCALEDINFEAS;

   return SCIP_OKAY;
}

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   lpi->solmethod = 'p';
   return lpiSolve(lpi, "p");
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   lpi->solmethod = 'd';
   return lpiSolve(lpi, "d");
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{
   SCIP_RETCODE retval;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   lpi->solmethod = 'b';

   CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_CROSSOVER, crossover) );

   retval = lpiSolve(lpi, "b");
   lpi->solisbasic = crossover;

   return retval;
}

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   /* currently do nothing */
   return SCIP_OKAY;
}

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   /* currently do nothing */
   return SCIP_OKAY;
}

/** performs strong branching iterations on one candidate */
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
   int objsen;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   SCIPdebugMessage("calling Xpress strong branching on variable %d (%d iterations)\n", col, itlim);

   /* results of Xpress are valid in any case */
   *downvalid = TRUE;
   *upvalid = TRUE;

   SCIPdebugMessage(" -> strong branching on integral variable\n");

   if( iter != NULL )
      *iter = 0;

   objsen = lpi->objsense;

#if (XPVERSION >= 18)
   /* From version 17.01.01 we have the dedicated strong branching function */
   /* XPRSstrongbranch(). Note: It does not work with version 17.01.00. */
   {
      int    mbndind[2];
      double dbndval[2];
      char   cbndtype[2];
      double dobjval[2];
      int    mstatus[2];

      /* Set the branching bounds (down first, up second). */
      mbndind[0]  = col;
      dbndval[0]  = EPSCEIL(psol-1.0, 1e-06);
      cbndtype[0] = 'U';
      mbndind[1]  = col;
      dbndval[1]  = EPSFLOOR(psol+1.0, 1e-06);
      cbndtype[1] = 'L';

      /* Apply strong branching to the two branches. */
      CHECK_ZERO( XPRSstrongbranch(lpi->xprslp, 2, mbndind, cbndtype, dbndval, itlim, dobjval, mstatus) );

      /* Get the objective of the down branch. */
      if ((mstatus[0] == XPRS_LP_INFEAS) || (mstatus[0] == XPRS_LP_CUTOFF_IN_DUAL))
      {
         *down = objsen == +1 ? 1e+40 : -1e+40;
      }
      else if ((mstatus[0] == XPRS_LP_OPTIMAL) || (mstatus[0] == XPRS_LP_UNFINISHED))
      {
         *down = dobjval[0];
      }
      else
      {
         /* Something weird happened. */
         *downvalid = FALSE;
      }

      /* Get the objective of the up branch. */
      if ((mstatus[1] == XPRS_LP_INFEAS) || (mstatus[1] == XPRS_LP_CUTOFF_IN_DUAL))
      {
         *up = objsen == +1 ? 1e+40 : -1e+40;
      }
      else if ((mstatus[1] == XPRS_LP_OPTIMAL) || (mstatus[1] == XPRS_LP_UNFINISHED))
      {
         *up = dobjval[1];
      }
      else
      {
         /* Something weird happened. */
         *upvalid = FALSE;
      }

      /* When using the XPRSstrongbranch function we are unable to provide */
      /* an iteration count. */
      if (iter)
         *iter = -1;
   }

#else

   {
      const char lbound = 'L';
      const char ubound = 'U';
      SCIP_Real oldlb;
      SCIP_Real oldub;
      SCIP_Real newlb;
      SCIP_Real newub;
      int olditlim;
      int it;

      /* save current LP basis and bounds*/
      SCIP_CALL( getBase(lpi) );
      CHECK_ZERO( XPRSgetlb(lpi->xprslp, &oldlb, col, col) );
      CHECK_ZERO( XPRSgetub(lpi->xprslp, &oldub, col, col) );

      /* save old iteration limit and set iteration limit to strong */
      /* branching limit */
      if( itlim > XPRS_MAXINT ) itlim = XPRS_MAXINT;
      CHECK_ZERO( XPRSgetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, &olditlim) );
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, itlim) );

      /* down branch */
      newub = EPSCEIL(psol-1.0, 1e-06);
      if( newub >= oldlb - 0.5 )
      {
         CHECK_ZERO( XPRSchgbounds(lpi->xprslp, 1, &col, &ubound, &newub) );
         SCIP_CALL( SCIPlpiSolveDual(lpi) );
         if( SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) )
            *down = objsen == +1 ? 1e+40 : -1e+40;
         else if( SCIPlpiIsOptimal(lpi) || SCIPlpiIsIterlimExc(lpi) )
         {
            SCIP_CALL( SCIPlpiGetObjval(lpi, down) );
         }
         else
            *down = objsen == +1 ? 1e+40 : -1e+40;
         if( iter != NULL )
         {
            SCIP_CALL( SCIPlpiGetIterations(lpi, &it) );
            *iter += it;
         }
         SCIPdebugMessage(" -> down (x%d <= %g): %g\n", col, newub, *down);

         CHECK_ZERO( XPRSchgbounds(lpi->xprslp, 1, &col, &ubound, &oldub) );
         SCIP_CALL( setBase(lpi) );
      }
      else
         *down = objsen == +1 ? 1e+40 : -1e+40;

      /* up branch */
      newlb = EPSFLOOR(psol+1.0, 1e-06);
      if( newlb <= oldub + 0.5 )
      {
         CHECK_ZERO( XPRSchgbounds(lpi->xprslp, 1, &col, &lbound, &newlb) );
         SCIP_CALL( SCIPlpiSolveDual(lpi) );
         if( SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) )
            *up = objsen == +1 ? 1e+40 : -1e+40;
         else if( SCIPlpiIsOptimal(lpi) || SCIPlpiIsIterlimExc(lpi) )
         {
            SCIP_CALL( SCIPlpiGetObjval(lpi, up) );
         }
         else
            *up = objsen == +1 ? 1e+40 : -1e+40;
         if( iter != NULL )
         {
            SCIP_CALL( SCIPlpiGetIterations(lpi, &it) );
            *iter += it;
         }
         SCIPdebugMessage(" -> up  (x%d >= %g): %g\n", col, newlb, *up);

         CHECK_ZERO( XPRSchgbounds(lpi->xprslp, 1, &col, &lbound, &oldlb) );
         SCIP_CALL( setBase(lpi) );
      }
      else
         *up = objsen == +1 ? 1e+40 : -1e+40;

      /* reset iteration limit */
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, olditlim) );

   }

#endif

   return SCIP_OKAY;
}

/** performs strong branching iterations on given candidates */
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
   int objsen;
   int j;

   assert( lpi != NULL );
   assert( lpi->xprslp != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );

   SCIPdebugMessage("calling Xpress strong branching on %d variables (%d iterations)\n", ncols, itlim);

   if( iter != NULL )
      *iter = 0;

   objsen = lpi->objsense;

#if (XPVERSION >= 18)
   /* From version 17.01.01 we have the dedicated strong branching function */
   /* XPRSstrongbranch(). Note: It does not work with version 17.01.00. */
   {
      int*    mbndind;
      double* dbndval;
      char*   cbndtype;
      double* dobjval;
      int*    mstatus;

      /* Set the branching bounds (down first, up second). */
      SCIP_ALLOC( BMSallocMemoryArray(&mbndind, 2*ncols) );
      SCIP_ALLOC( BMSallocMemoryArray(&dbndval, 2*ncols) );
      SCIP_ALLOC( BMSallocMemoryArray(&cbndtype, 2*ncols) );
      SCIP_ALLOC( BMSallocMemoryArray(&dobjval, 2*ncols) );
      SCIP_ALLOC( BMSallocMemoryArray(&mstatus, 2*ncols) );
      
      for (j = 0; j < ncols; ++j)
      {
         mbndind[2*j]  = cols[j];
         dbndval[2*j]  = EPSCEIL(psols[j] - 1.0, 1e-06);
         cbndtype[2*j] = 'U';

         mbndind[2*j+1]  = cols[j];
         dbndval[2*j+1]  = EPSFLOOR(psols[j] + 1.0, 1e-06);
         cbndtype[2*j+1] = 'L';
      }

      /* Apply strong branching to the 2*ncols branches. */
      CHECK_ZERO( XPRSstrongbranch(lpi->xprslp, 2*ncols, mbndind, cbndtype, dbndval, itlim, dobjval, mstatus) );

      for (j = 0; j < ncols; ++j)
      {
         upvalid[j]   = TRUE;
         downvalid[j] = TRUE;

         /* Get the objective of the down branch. */
         if ((mstatus[2*j] == XPRS_LP_INFEAS) || (mstatus[2*j] == XPRS_LP_CUTOFF_IN_DUAL))
            down[j] = objsen == +1 ? 1e+40 : -1e+40;
         else if ((mstatus[2*j] == XPRS_LP_OPTIMAL) || (mstatus[2*j] == XPRS_LP_UNFINISHED))
            down[j] = dobjval[2*j];
         else
         {
            /* Something weird happened. */
            downvalid[j] = FALSE;
         }

         /* Get the objective of the up branch. */
         if ((mstatus[2*j+1] == XPRS_LP_INFEAS) || (mstatus[2*j+1] == XPRS_LP_CUTOFF_IN_DUAL))
            up[j] = objsen == +1 ? 1e+40 : -1e+40;
         else if ((mstatus[2*j+1] == XPRS_LP_OPTIMAL) || (mstatus[2*j+1] == XPRS_LP_UNFINISHED))
            up[j] = dobjval[2*j+1];
         else
         {
            /* Something weird happened. */
            upvalid[j] = FALSE;
         }
      }

      /* When using the XPRSstrongbranch function we are unable to provide */
      /* an iteration count. */
      if (iter)
         *iter = -1;

      BMSfreeMemoryArray(&mstatus);
      BMSfreeMemoryArray(&dobjval);
      BMSfreeMemoryArray(&cbndtype);
      BMSfreeMemoryArray(&dbndval);
      BMSfreeMemoryArray(&mbndind)
   }
#else
   {
      int olditlim;

      /* save current LP basis */
      SCIP_CALL( getBase(lpi) );

      /* save old iteration limit and set iteration limit to strong */
      /* branching limit */
      if ( itlim > XPRS_MAXINT )
         itlim = XPRS_MAXINT;
      CHECK_ZERO( XPRSgetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, &olditlim) );
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, itlim) );

      for (j = 0; j < ncols; ++j)
      {
         const char lbound = 'L';
         const char ubound = 'U';
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Real newlb;
         SCIP_Real newub;
         SCIP_Real psol;
         int col;
         int it;

         upvalid[j]   = TRUE;
         downvalid[j] = TRUE;

         col = cols[j];
         psol = psols[j];

         /* save current LP bounds*/
         CHECK_ZERO( XPRSgetlb(lpi->xprslp, &oldlb, col, col) );
         CHECK_ZERO( XPRSgetub(lpi->xprslp, &oldub, col, col) );

         /* down branch */
         newub = EPSCEIL(psol-1.0, 1e-06);
         if( newub >= oldlb - 0.5 )
         {
            CHECK_ZERO( XPRSchgbounds(lpi->xprslp, 1, &col, &ubound, &newub) );
            SCIP_CALL( SCIPlpiSolveDual(lpi) );
            if( SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) )
               down[j] = objsen == +1 ? 1e+40 : -1e+40;
            else if( SCIPlpiIsOptimal(lpi) || SCIPlpiIsIterlimExc(lpi) )
            {
               SCIP_CALL( SCIPlpiGetObjval(lpi, &(down[j])) );
            }
            else
               down[j] = objsen == +1 ? 1e+40 : -1e+40;

            if( iter != NULL )
            {
               SCIP_CALL( SCIPlpiGetIterations(lpi, &it) );
               *iter += it;
            }
            SCIPdebugMessage(" -> down (x%d <= %g): %g\n", col, newub, down[j]);

            CHECK_ZERO( XPRSchgbounds(lpi->xprslp, 1, &col, &ubound, &oldub) );
            SCIP_CALL( setBase(lpi) );
         }
         else
            down[j] = objsen == +1 ? 1e+40 : -1e+40;

         /* up branch */
         newlb = EPSFLOOR(psol+1.0, 1e-06);
         if( newlb <= oldub + 0.5 )
         {
            CHECK_ZERO( XPRSchgbounds(lpi->xprslp, 1, &col, &lbound, &newlb) );
            SCIP_CALL( SCIPlpiSolveDual(lpi) );
            if( SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) )
               up[j] = objsen == +1 ? 1e+40 : -1e+40;
            else if( SCIPlpiIsOptimal(lpi) || SCIPlpiIsIterlimExc(lpi) )
            {
               SCIP_CALL( SCIPlpiGetObjval(lpi, up) );
            }
            else
               up[j] = objsen == +1 ? 1e+40 : -1e+40;
            if( iter != NULL )
            {
               SCIP_CALL( SCIPlpiGetIterations(lpi, &it) );
               *iter += it;
            }
            SCIPdebugMessage(" -> up  (x%d >= %g): %g\n", col, newlb, up[j]);

            CHECK_ZERO( XPRSchgbounds(lpi->xprslp, 1, &col, &lbound, &oldlb) );
            SCIP_CALL( setBase(lpi) );
         }
         else
            up[j] = objsen == +1 ? 1e+40 : -1e+40;
      }

      /* reset iteration limit */
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, olditlim) );
   }
#endif

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
   /* pass call on to lpiStrongbranches() */
   SCIP_CALL( lpiStrongbranches(lpi, cols, ncols, psols, itlim, down, up, downvalid, upvalid, iter) );

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

   return (lpi->solstat != -1);
}

/** gets information about primal and dual feasibility of the current LP solution */
/** here "true" should mean feasible, "false" should mean unknown                 */
SCIP_RETCODE SCIPlpiGetSolFeasibility(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(primalfeasible != NULL);
   assert(dualfeasible != NULL);

   SCIPdebugMessage("getting solution feasibility\n");

   *primalfeasible = (SCIP_Bool) (lpi->solstat==XPRS_LP_OPTIMAL || lpi->solstat == XPRS_LP_OPTIMAL_SCALEDINFEAS || (lpi->solmethod == 'p' && lpi->solstat==XPRS_LP_UNBOUNDED));
   *dualfeasible = (SCIP_Bool) (lpi->solstat==XPRS_LP_OPTIMAL || lpi->solstat == XPRS_LP_OPTIMAL_SCALEDINFEAS || (lpi->solmethod == 'd' && lpi->solstat==XPRS_LP_INFEAS));

   return SCIP_OKAY;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExistsPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == XPRS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

#if OLDRAYCODE
   if (lpi->solstat != XPRS_LP_UNBOUNDED)
      return FALSE;

   if (lpi->unbvec < 0)
      return FALSE;

   return TRUE;
#else
   {
      int hasRay;
      CHECK_ZERO( XPRSgetprimalray(lpi->xprslp, NULL, &hasRay) );
      return hasRay;
   }
#endif
}

/** returns TRUE iff LP is proven to be primal feasible and unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for primal unboundedness\n");

   /* If the solution status of Xpress is XPRS_LP_UNBOUNDED, it only means, there is an unbounded ray,
    * but not necessarily a feasible primal solution. If problem is declared LP_UNBOUNDED by dual,
    * we have no way to decide primal feasibility.
    */

   return lpi->solstat == XPRS_LP_UNBOUNDED && lpi->solmethod == 'p';
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for primal infeasibility\n");

   return (lpi->solstat == XPRS_LP_INFEAS);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for primal feasibility\n");

   /* problem is optimal or unbounded found by primal */
   return lpi->solstat == XPRS_LP_OPTIMAL || lpi->solstat == XPRS_LP_OPTIMAL_SCALEDINFEAS || (lpi->solstat == XPRS_LP_UNBOUNDED && lpi->solmethod == 'p');
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == XPRS_LP_INFEAS);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

#if OLDRAYCODE
   if (lpi->solmethod != 'p')
   {
      /* We can only get a dual ray from primal. */
      SCIP_CALL( SCIPlpiSolvePrimal(lpi) );
   }

   if (lpi->solstat != XPRS_LP_INFEAS)
      return FALSE;

   return TRUE;
#else
   {
      int hasRay;
      CHECK_ZERO( XPRSgetdualray(lpi->xprslp, NULL, &hasRay) );
      return hasRay;
   }
#endif
}

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for dual unboundedness\n");

   return ((lpi->solstat == XPRS_LP_INFEAS) && (lpi->solmethod == 'd'));
}

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for dual infeasibility\n");

   return (lpi->solstat == XPRS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for dual feasibility\n");

   /* problem is optimal or infeasible found by dual */
   return lpi->solstat == XPRS_LP_OPTIMAL || lpi->solstat == XPRS_LP_OPTIMAL_SCALEDINFEAS || (lpi->solstat == XPRS_LP_INFEAS && lpi->solmethod == 'd');
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == XPRS_LP_OPTIMAL) || (lpi->solstat == XPRS_LP_OPTIMAL_SCALEDINFEAS);
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for stability: Xpress solstat = %d\n", lpi->solstat);

   /* If the solution status of Xpress is XPRS_LP_UNBOUNDED, it only means, there is an unbounded ray,
    * but not necessarily a feasible primal solution. If primalfeasible == FALSE, we interpret this
    * result as instability, s.t. the problem is resolved from scratch
    */
   if( lpi->solstat == XPRS_LP_UNBOUNDED )
   {
      int pinfeas;

      CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_PRIMALINFEAS, &pinfeas) );

      if( pinfeas )
         return FALSE;
   }
   else if ( lpi->solstat == XPRS_LP_OPTIMAL_SCALEDINFEAS )
   {
      /* Presolved problem was solved to optimality but infeasibilities */
      /* were introduced by postsolve. */
      return FALSE;
   }

   return TRUE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == XPRS_LP_CUTOFF_IN_DUAL);
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int lpiter;
   int lpiterlimit;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_SIMPLEXITER, &lpiter) );
   CHECK_ZERO( XPRSgetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, &lpiterlimit) );

   if ( (lpi->solstat == XPRS_LP_UNFINISHED) && (lpiter >= lpiterlimit) )
      return TRUE;
   else
      return FALSE;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int lpiter;
   int lpiterlimit;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_SIMPLEXITER, &lpiter) );
   CHECK_ZERO( XPRSgetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, &lpiterlimit) );

   if ( (lpi->solstat == XPRS_LP_UNFINISHED) && (lpiter < lpiterlimit) )
      return TRUE;
   else
      return FALSE;
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   return lpi->solstat;
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   /* Nothing to do here for Xpress. */
   *success = TRUE;

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("getting solution's objective value\n");

   CHECK_ZERO( XPRSgetdblattrib(lpi->xprslp, XPRS_LPOBJVAL, objval) );

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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("getting solution\n");

   CHECK_ZERO( XPRSgetsol(lpi->xprslp, primsol, activity, dualsol, redcost) );
   if (objval)
   {
      CHECK_ZERO( XPRSgetdblattrib(lpi->xprslp, XPRS_LPOBJVAL, objval) );
   }

   if( activity != NULL )
   {
      /* Convert the slack values into activity values. */
      int r;
      int nrows;

      CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
      SCIP_CALL( ensureSidechgMem(lpi, nrows) );
      CHECK_ZERO( XPRSgetrhs(lpi->xprslp, lpi->rhsarray, 0, nrows-1) );
      for (r = 0; r < nrows; r++)
         activity[r] = lpi->rhsarray[r] - activity[r];
   }

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{
#if OLDRAYCODE
   int i;
   int irow;
   int nrows;
   int ncols;
   int cfirst;
   double dmult;

   int    *bind = NULL;
   double *bvec = NULL;
#endif

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(ray != NULL);
   assert(lpi->solstat >= 0);

#if OLDRAYCODE
   /* Check if it is possible for us to extract a primal ray. */
   if ((lpi->solstat != XPRS_LP_UNBOUNDED) || (lpi->unbvec < 0))
   {
      /* Not unbounded or the optimizer didn't return a ray index. */
      return SCIP_LPERROR;
   }

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &ncols) );
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_SPAREROWS, &cfirst) );
   cfirst += nrows;
   SCIP_CALL( getBase(lpi) );
   if (lpi->solmethod == 'd')
      return SCIP_LPERROR; /* Unboundedness found by dual - nothing we can do about it. */

   /* At this point we should be fairly confident that unboundedness */
   /* is caused by primal trying to pivot a non-basic variable into */
   /* the basis. */
   memset(ray, 0, ncols*sizeof(*ray));

   /* Get the simplex tableau column of the pivot variable. */
   SCIP_ALLOC( BMSallocMemoryArray(&bind, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&bvec, nrows) );
   memset(bvec, 0, nrows*sizeof(*bvec));
   if (lpi->unbvec >= cfirst)
   {
      /* We have a non-basic column - extract it... */
      int icol = lpi->unbvec-cfirst;
      int ifirst = 0;
      int nnonz = nrows;
      int    *mcolind = NULL;
      double *dcolval = NULL;

      assert(lpi->cstat[icol] != 1);

      dmult = lpi->cstat[icol] ? -1.0 : +1.0;
      SCIP_ALLOC( BMSallocMemoryArray(&dcolval, nrows) );
      SCIP_ALLOC( BMSallocMemoryArray(&mcolind, nrows) );
      SCIP_CALL( SCIPlpiGetCols(lpi, icol, icol, NULL, NULL, &nnonz, &ifirst, mcolind, dcolval) );

      /* ... and unpack it. */
      assert(nnonz > 0);
      for (i = 0; i < nnonz; i++)
         bvec[mcolind[i]] = dmult*dcolval[i];
      ray[icol] = dmult;
      BMSfreeMemoryArray(&dcolval);
      BMSfreeMemoryArray(&mcolind);
   }
   else
   {
      /* We have a non-basic row. */
      irow = lpi->unbvec;
      assert(lpi->rstat[irow] != 1);
      dmult = lpi->rstat[irow] ? -1.0 : +1.0;
      bvec[irow] = dmult;
   }

   /* Get the simplex tableau column and the variable basic in each row. */
   CHECK_ZERO( XPRSftran(lpi->xprslp, bvec) );
   CHECK_ZERO( XPRSgetpivotorder(lpi->xprslp, bind) );

   /* Save the ray. */
   for (irow = 0; irow < nrows; irow++)
   {
      if (bind[irow] >= cfirst)
         ray[bind[irow]-cfirst] = bvec[irow];
   }

   BMSfreeMemoryArray(&bind);
   BMSfreeMemoryArray(&bvec);
#else
   {
      int hasRay;
      CHECK_ZERO( XPRSgetprimalray(lpi->xprslp, ray, &hasRay) );
      if (!hasRay)
         return SCIP_LPERROR;
   }
#endif

   return SCIP_OKAY;
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual Farkas row multipliers */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);
   assert(dualfarkas != NULL);

#if OLDRAYCODE
   /* Check if it is possible for us to extract a dual ray. */
   if (lpi->solstat != XPRS_LP_INFEAS)
   {
      /* Not infeasible. */
      return SCIP_LPERROR;
   }
   if (lpi->solmethod != 'p')
      return SCIP_LPERROR;

   /* The required Farkas multipliers should be the duals set up by */
   /* phase I primal. */
   CHECK_ZERO( XPRSgetsol(lpi->xprslp, NULL, NULL, dualfarkas, NULL) );
#else
   {
      int hasRay;
      CHECK_ZERO( XPRSgetdualray(lpi->xprslp, dualfarkas, &hasRay) );
      if (!hasRay)
         return SCIP_LPERROR;
   }
#endif

   return SCIP_OKAY;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpi->solstat >= 0);
   assert(iterations != NULL);

   *iterations = lpi->iterations;

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
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("saving Xpress basis into %p/%p\n", (void*)rstat, (void*)cstat);

   CHECK_ZERO( XPRSgetbasis(lpi->xprslp, rstat, cstat) );

   /* Convert Xpress basis status into the SCIP format. */
#if (XPVERSION < 17)
   if (cstat)
   {
      int c;
      int ncols;

      CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );

      /* Before version 17.00.00 of Xpress, there was no superbasic status in xprs */
      for (c = 0; c < ncols; c++)
      {
         assert(cstat[c] >=0 && cstat[c] <= 2);
         if (cstat[c] == 0)
         {  /* Check if it might be super-basic. */
            double dlb;
            CHECK_ZERO( XPRSgetlb(lpi->xprslp, &dlb, c, c) );
            if (dlb <= XPRS_MINUSINFINITY)
               cstat[c] = SCIP_BASESTAT_ZERO;
         }
      }
   }
#endif

   return SCIP_OKAY;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiSetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(cstat != NULL);
   assert(rstat != NULL);

   SCIPdebugMessage("loading basis %p/%p into Xpress\n", (void*)rstat, (void*)cstat);

   invalidateSolution(lpi);

#if (XPVERSION < 17)
   {
      int c;
      int ncols;
      int *cstat_xprs = NULL;

      /* Set any super-basic variables to be at lower bound. */
      CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );
      SCIP_ALLOC( BMSallocMemoryArray(&cstat_xprs, ncols) );
      /* We can't set a super-basic status so set it at lower bound instead. */
      for (c = 0; c < ncols; c++)
         cstat_xprs[c] = (cstat[c] == 3) ? 0 : cstat[c];
      CHECK_ZERO( XPRSloadbasis(lpi->xprslp, rstat, cstat_xprs) );
      BMSfreeMemoryArray(&cstat_xprs);
   }
#else
   CHECK_ZERO( XPRSloadbasis(lpi->xprslp, rstat, cstat) );
#endif

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   int r;
   int nrows;
   int irspace;

   /* In the basis methods we assume that xprs basis flags coincide with scip, so assert it */
   assert((0 == SCIP_BASESTAT_LOWER) && (1 == SCIP_BASESTAT_BASIC) && (2 == SCIP_BASESTAT_UPPER) && (3 == SCIP_BASESTAT_ZERO));

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(bind != NULL);

   SCIPdebugMessage("getting basis information\n");

   CHECK_ZERO( XPRSgetpivotorder(lpi->xprslp, bind) );

   /* Reindex variables to match those of SCIP. */
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_SPAREROWS, &irspace) );
   irspace += nrows;

   for (r = 0; r < nrows; r++)
   {
      if (bind[r] < nrows)
         bind[r] = -bind[r]-1;
      else
      {
         assert(bind[r] >= irspace);
         bind[r] = bind[r] - irspace;
      }
   }

   return SCIP_OKAY;
}

/** get dense row of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiGetBInvRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   int nrows;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("getting binv-row %d\n", row);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   memset(coef, 0, nrows*sizeof(*coef));
   coef[row] = 1.0;
   CHECK_ZERO( XPRSbtran(lpi->xprslp, coef) );

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
   int nrows;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("getting binv-col %d\n", c);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   memset(coef, 0, nrows*sizeof(*coef));
   coef[c] = 1.0;
   CHECK_ZERO( XPRSftran(lpi->xprslp, coef) );

   return SCIP_OKAY;
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow_in,         /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            val                 /**< vector to return coefficients */
   )
{
   int c;
   int nrows, ncols;
   int nnonz;

   SCIP_Real *binvrow = NULL;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(binvrow != NULL);

   SCIPdebugMessage("getting binva-row %d\n", r);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );

   /* Get the row of the basis inverse. */
   if (binvrow_in)
      binvrow = (double *) binvrow_in;
   else
      SCIP_ALLOC( BMSallocMemoryArray(&binvrow, nrows) );

   /* We need space to extract a single column. */
   SCIP_CALL( ensureValMem(lpi, nrows) );

   for (c = 0; c < ncols; c++)
   {
      int i;
      double dsum = 0.0;

      /* Extract the column. */
      CHECK_ZERO( XPRSgetcols(lpi->xprslp, NULL, lpi->indarray, lpi->valarray, nrows, &nnonz, c, c) );

      /* Price out the column. */
      for (i = 0; i < nnonz; i++)
         dsum += binvrow[lpi->indarray[i]]*lpi->valarray[i];
      val[c] = dsum;
   }

   /* Free allocated memory. */
   if (binvrow_in == NULL)
      BMSfreeMemoryArray(&binvrow);

   return SCIP_OKAY;
}

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{
   int i;
   int nrows;
   int nnonz;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("getting binv-col %d\n", c);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );

   /* We need space to extract the column. */
   SCIP_CALL( ensureValMem(lpi, nrows) );

   /* Get the column to transform. */
   CHECK_ZERO( XPRSgetcols(lpi->xprslp, NULL, lpi->indarray, lpi->valarray, nrows, &nnonz, c, c) );

   /* Transform the column. */
   memset(coef, 0, nrows*sizeof(*coef));
   for (i = 0; i < nnonz; i++)
      coef[lpi->indarray[i]] = lpi->valarray[i];
   CHECK_ZERO( XPRSbtran(lpi->xprslp, coef) );

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

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(lpistate != NULL);

   /* if there is no basis information available (e.g. after barrier without crossover), no state can be saved */
   if( !lpi->solisbasic )
   {
      *lpistate = NULL;
      return SCIP_OKAY;
   }

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );
   assert(ncols >= 0);
   assert(nrows >= 0);

   /* allocate lpistate data */
   SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows) );

   SCIPdebugMessage("storing Xpress LPI state in %p (%d cols, %d rows)\n", (void*)*lpistate, ncols, nrows);

   /* get unpacked basis information from Xpress */
   SCIP_CALL( getBase(lpi) );

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
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   )
{
   int nrows;
   int ncols;
   int i;

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_ROWS, &nrows) );
   CHECK_ZERO( XPRSgetintattrib(lpi->xprslp, XPRS_COLS, &ncols) );

   assert(lpistate == NULL || lpistate->ncols == ncols);
   assert(lpistate == NULL || lpistate->nrows == nrows);

   /* if there was no basis information available, the LPI state was not stored */
   if( lpistate == NULL )
      return SCIP_OKAY;

   SCIPdebugMessage("loading LPI state %p (%d cols, %d rows) into Xpress\n", (void*)lpistate, lpistate->ncols, lpistate->nrows);

   if( lpistate->ncols == 0 || lpistate->nrows == 0 )
      return SCIP_OKAY;

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, lpistate->ncols) );
   SCIP_CALL( ensureRstatMem(lpi, lpistate->nrows) );

   /* unpack LPi state data */
   lpistateUnpack(lpistate, lpi->cstat, lpi->rstat);

   /* extend the basis to the current LP beyond the previously existing columns */
   for (i = lpistate->ncols; i < ncols; ++i)
   {
      SCIP_Real bnd;
      CHECK_ZERO( XPRSgetlb(lpi->xprslp, &bnd, i, i) );
      if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         CHECK_ZERO( XPRSgetub(lpi->xprslp, &bnd, i, i) );
         if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
            lpi->cstat[i] = SCIP_BASESTAT_ZERO;  /* variable is free */
         else
            lpi->cstat[i] = SCIP_BASESTAT_UPPER; /* use finite upper bound */
      }
      else
         lpi->cstat[i] = SCIP_BASESTAT_LOWER;    /* use finite lower bound */
   }
   for (i = lpistate->nrows; i < nrows; ++i)
      lpi->rstat[i] = SCIP_BASESTAT_BASIC;

   /* load basis information into Xpress */
   SCIP_CALL( setBase(lpi) );

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   /**@todo implement SCIPlpiClearState() for Xpress */
   SCIPmessagePrintWarning(lpi->messagehdlr, "Xpress interface does not implement SCIPlpiClearState()\n");

   return SCIP_OKAY;
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiFreeState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   assert(lpi != NULL);
   assert(lpistate != NULL);

   if( *lpistate != NULL )
   {
      lpistateFree(lpistate, blkmem);
   }

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   return (lpistate != NULL);
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("reading LP state from file <%s>\n", fname);

   CHECK_ZERO( XPRSreadbasis(lpi->xprslp, fname, "") );

   return SCIP_OKAY;
}

/** writes LP state (like basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("writing LP state to file <%s>\n", fname);

   CHECK_ZERO( XPRSwritebasis(lpi->xprslp, fname, "") );

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
   int ictrlval;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(ival != NULL);

   SCIPdebugMessage("getting int parameter %d\n", type);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      CHECK_ZERO( XPRSgetintcontrol(lpi->xprslp, XPRS_KEEPBASIS, &ictrlval) );
      *ival = (ictrlval == 0);
      break;
   case SCIP_LPPAR_SCALING:
      CHECK_ZERO( XPRSgetintcontrol(lpi->xprslp, XPRS_SCALING, &ictrlval) );
      *ival = (ictrlval != 0);
      break;
   case SCIP_LPPAR_PRESOLVING:
      *ival = lpi->par_presolve;
      break;
   case SCIP_LPPAR_LPINFO:
      CHECK_ZERO( XPRSgetintcontrol(lpi->xprslp, XPRS_OUTPUTLOG, &ictrlval) );
      *ival = (ictrlval != 0);
      break;
   case SCIP_LPPAR_LPITLIM:
      CHECK_ZERO( XPRSgetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, &ictrlval) );
      *ival = ictrlval;
      if( *ival >= XPRS_MAXINT )
         *ival = XPRS_MAXINT;
      break;
   case SCIP_LPPAR_FASTMIP:
      /* We treat this as a meta parameter to enable settings that make */
      /* reoptimization go faster. */
      *ival = lpi->par_fastlp;
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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("setting int parameter %d to %d\n", type, ival);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      assert(ival == TRUE || ival == FALSE);
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_KEEPBASIS, (ival == FALSE) ? 1 : 0) );
      break;
   case SCIP_LPPAR_SCALING:
      assert(ival == TRUE || ival == FALSE);
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_SCALING, (ival == TRUE) ? 35 : 0) );
      break;
   case SCIP_LPPAR_PRESOLVING:
      assert(ival == TRUE || ival == FALSE);
      lpi->par_presolve = ival;
      break;
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_OUTPUTLOG, (ival == TRUE) ? 1 : 0) );
      break;
   case SCIP_LPPAR_LPITLIM:
      ival = MIN(ival, XPRS_MAXINT);
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_LPITERLIMIT, ival) );
      break;
   case SCIP_LPPAR_FASTMIP:
      /* We treat this as a meta parameter to enable settings that make */
      /* reoptimization go faster. Nothing is set in Xpress until we solve */
      /* the problem. */
      lpi->par_fastlp = ival;
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
   int ictrlval;
   double dctrlval;

   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);
   assert(dval != NULL);

   SCIPdebugMessage("getting real parameter %d\n", type);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      CHECK_ZERO( XPRSgetdblcontrol(lpi->xprslp, XPRS_FEASTOL, &dctrlval) );
      *dval = dctrlval;
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      CHECK_ZERO( XPRSgetdblcontrol(lpi->xprslp, XPRS_OPTIMALITYTOL, &dctrlval) );
      *dval = dctrlval;
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:
      CHECK_ZERO( XPRSgetdblcontrol(lpi->xprslp, XPRS_BARGAPSTOP, &dctrlval) );
      *dval = dctrlval;
      break;
   case SCIP_LPPAR_LPTILIM:
      CHECK_ZERO( XPRSgetintcontrol(lpi->xprslp, XPRS_MAXTIME, &ictrlval) );
      *dval = (double) ictrlval;
      break;
   case SCIP_LPPAR_MARKOWITZ:
      CHECK_ZERO( XPRSgetdblcontrol(lpi->xprslp, XPRS_MARKOWITZTOL, &dctrlval) );
      *dval = dctrlval;
      break;
   case SCIP_LPPAR_LOBJLIM:
      *dval = lpi->par_lobjlim;
      break;
   case SCIP_LPPAR_UOBJLIM:
      *dval = lpi->par_uobjlim;
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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("setting real parameter %d to %g\n", type, dval);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      CHECK_ZERO( XPRSsetdblcontrol(lpi->xprslp, XPRS_FEASTOL, dval) );
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      CHECK_ZERO( XPRSsetdblcontrol(lpi->xprslp, XPRS_OPTIMALITYTOL, dval) );
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:
      CHECK_ZERO( XPRSsetdblcontrol(lpi->xprslp, XPRS_BARGAPSTOP, dval) );
      break;
   case SCIP_LPPAR_LPTILIM:
   {
      int ival = (int) dval;
      CHECK_ZERO( XPRSsetintcontrol(lpi->xprslp, XPRS_MAXTIME, ival) );
      break;
   }
   case SCIP_LPPAR_MARKOWITZ:
      CHECK_ZERO( XPRSsetdblcontrol(lpi->xprslp, XPRS_MARKOWITZTOL, dval) );
      break;
   case SCIP_LPPAR_LOBJLIM:
      lpi->par_lobjlim = dval;
      break;
   case SCIP_LPPAR_UOBJLIM:
      lpi->par_uobjlim = dval;
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
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   return XPRS_PLUSINFINITY;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   return (val >= XPRS_PLUSINFINITY);
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
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("reading LP from file <%s>\n", fname);

   CHECK_ZERO( XPRSreadprob(lpi->xprslp, fname, "") );

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->xprslp != NULL);

   SCIPdebugMessage("writing LP to file <%s>\n", fname);

   CHECK_ZERO( XPRSwriteprob(lpi->xprslp, fname, "p") );

   return SCIP_OKAY;
}

/**@} */
