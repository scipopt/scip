/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: lpi_cpx.c,v 1.47 2003/12/03 18:08:12 bzfpfend Exp $"

/**@file   lpi_cpx.c
 * @brief  LP interface for CPLEX 8.0
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "cplex.h"
#include "bitencode.h"
#include "lpi.h"
#include "message.h"


#define CHECK_ZERO(x) { int _restat_;                                               \
                        if( (_restat_ = (x)) != 0 )                                 \
                        {                                                           \
                           errorMessage("LP Error: CPLEX returned %d\n", _restat_); \
                           return SCIP_LPERROR;                                     \
                        }                                                           \
                      }

#define ABORT_ZERO(x) { int _restat_;                                               \
                        if( (_restat_ = (x)) != 0 )                                 \
                        {                                                           \
                           errorMessage("LP Error: CPLEX returned %d\n", _restat_); \
                           abort();                                                 \
                        }                                                           \
                      }

#define NOTCALLED  -1


typedef DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET DUALPACKETSIZE
typedef DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET DUALPACKETSIZE

/* CPLEX parameter lists which can be changed */
#define NUMINTPARAM  9
static const int intparam[NUMINTPARAM] = {
   CPX_PARAM_ADVIND,
   CPX_PARAM_ITLIM,
   CPX_PARAM_FASTMIP,
   CPX_PARAM_SCAIND,
   CPX_PARAM_PREIND,
   CPX_PARAM_PPRIIND,
   CPX_PARAM_DPRIIND,
   CPX_PARAM_SIMDISPLAY,
   CPX_PARAM_SCRIND
};
#define NUMDBLPARAM  5
static const int dblparam[NUMDBLPARAM] = {
   CPX_PARAM_EPRHS,
   CPX_PARAM_EPOPT,
   CPX_PARAM_OBJLLIM,
   CPX_PARAM_OBJULIM,
   CPX_PARAM_TILIM
};

/** CPLEX parameter settings */
struct CPXParam
{
   int              intparval[NUMINTPARAM]; /**< integer parameter values */
   double           dblparval[NUMDBLPARAM]; /**< double parameter values */
};
typedef struct CPXParam CPXPARAM;

/** LP interface */
struct LPi
{
   CPXLPptr         cpxlp;              /**< CPLEX LP pointer */
   int              solstat;            /**< solution status of last optimization call */
   CPXPARAM         cpxparam;           /**< actual parameter values for this LP */
   char*            larray;             /**< array with 'L' entries for changing lower bounds */
   char*            uarray;             /**< array with 'U' entries for changing upper bounds */
   char*            senarray;           /**< array for storing row senses */
   Real*            rhsarray;           /**< array for storing rhs values */
   Real*            rngarray;           /**< array for storing range values */
   Real*            valarray;           /**< array for storing coefficient values */
   int*             rngindarray;        /**< array for storing row indices with range values */
   int*             cstat;              /**< array for storing column basis status */
   int*             rstat;              /**< array for storing row basis status */
   int*             indarray;           /**< array for storing coefficient indices */
   int              boundchgsize;       /**< size of larray and uarray */
   int              sidechgsize;        /**< size of senarray, rngarray, and rngindarray */
   int              valsize;            /**< size of valarray and indarray */
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
   double*          dnorm;              /**< dual norms of variables */
};


static CPXENVptr    cpxenv = NULL;      /**< CPLEX environment */
static CPXPARAM     defparam;           /**< default CPLEX parameters */
static CPXPARAM     actparam;           /**< actual CPLEX parameters in the environment */
static int          numlp = 0;          /**< number of open LP objects */



/*
 * dynamic memory arrays
 */

/** resizes larray and uarray to have at least num entries */
static
RETCODE ensureBoundchgMem(
   LPI*             lpi,                /**< LP interface structure */
   int              num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->boundchgsize )
   {
      int newsize;
      int i;

      newsize = MAX(2*lpi->boundchgsize, num);
      ALLOC_OKAY( reallocMemoryArray(&lpi->larray, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&lpi->uarray, newsize) );
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

/** resizes senarray, rngarray, and rngindarray to have at least num entries */
static
RETCODE ensureSidechgMem(
   LPI*             lpi,                /**< LP interface structure */
   int              num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->sidechgsize )
   {
      int newsize;

      newsize = MAX(2*lpi->sidechgsize, num);
      ALLOC_OKAY( reallocMemoryArray(&lpi->senarray, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&lpi->rhsarray, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&lpi->rngarray, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&lpi->rngindarray, newsize) );
      lpi->sidechgsize = newsize;
   }
   assert(num <= lpi->sidechgsize);

   return SCIP_OKAY;
}

/** resizes valarray and indarray to have at least num entries */
static
RETCODE ensureValMem(
   LPI*             lpi,                /**< LP interface structure */
   int              num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->valsize )
   {
      int newsize;

      newsize = MAX(2*lpi->valsize, num);
      ALLOC_OKAY( reallocMemoryArray(&lpi->valarray, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&lpi->indarray, newsize) );
      lpi->valsize = newsize;
   }
   assert(num <= lpi->valsize);

   return SCIP_OKAY;
}

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
   (*lpistate)->dnorm = NULL;

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
   freeBlockMemoryArrayNull(memhdr, &(*lpistate)->dnorm, (*lpistate)->ncols);
   freeBlockMemory(memhdr, lpistate);
}



/*
 * local methods
 */

static
RETCODE getParameterValues(CPXPARAM* cpxparam)
{
   int i;
   
   assert(cpxenv != NULL);
   assert(cpxparam != NULL);

   for( i = 0; i < NUMINTPARAM; ++i )
   {
      CHECK_ZERO( CPXgetintparam(cpxenv, intparam[i], &(cpxparam->intparval[i])) );
   }
   for( i = 0; i < NUMDBLPARAM; ++i )
   {
      CHECK_ZERO( CPXgetdblparam(cpxenv, dblparam[i], &(cpxparam->dblparval[i])) );
   }

   return SCIP_OKAY;
}
   
static
void checkParameterValues(void)
{
#ifndef NDEBUG
   CPXPARAM par;
   int i;
   
   getParameterValues(&par);
   for( i = 0; i < NUMINTPARAM; ++i )
      assert(actparam.intparval[i] == par.intparval[i]);
   for( i = 0; i < NUMDBLPARAM; ++i )
      assert(actparam.dblparval[i] == par.dblparval[i]);
#endif
}

static
RETCODE setParameterValues(const CPXPARAM* cpxparam)
{
   int i;
   
   assert(cpxenv != NULL);
   assert(cpxparam != NULL);
   
   for( i = 0; i < NUMINTPARAM; ++i )
   {
      if( actparam.intparval[i] != cpxparam->intparval[i] )
      {
         actparam.intparval[i] = cpxparam->intparval[i];
         CHECK_ZERO( CPXsetintparam(cpxenv, intparam[i], actparam.intparval[i]) );
      }
   }
   for( i = 0; i < NUMDBLPARAM; ++i )
   {
      if( actparam.dblparval[i] != cpxparam->dblparval[i] )
      {
         actparam.dblparval[i] = cpxparam->dblparval[i];
         CHECK_ZERO( CPXsetdblparam(cpxenv, dblparam[i], actparam.dblparval[i]) );
      }
   }

   checkParameterValues();

   return SCIP_OKAY;
}

static
void copyParameterValues(CPXPARAM* dest, const CPXPARAM* source)
{
   int i;

   for( i = 0; i < NUMINTPARAM; ++i )
      dest->intparval[i] = source->intparval[i];
   for( i = 0; i < NUMDBLPARAM; ++i )
      dest->dblparval[i] = source->dblparval[i];
}

static
int getIntParam(LPI* lpi, const int param)
{
   int i;
   
   assert(lpi != NULL);

   for( i = 0; i < NUMINTPARAM; ++i )
      if( intparam[i] == param )
         return lpi->cpxparam.intparval[i];

   errorMessage("Unknown CPLEX integer parameter\n");
   abort();
}

static
double getDblParam(LPI* lpi, const int param)
{
   int i;

   assert(lpi != NULL);

   for( i = 0; i < NUMDBLPARAM; ++i )
      if( dblparam[i] == param )
         return lpi->cpxparam.dblparval[i];

   errorMessage("Unknown CPLEX double parameter\n");
   abort();
}

static
void setIntParam(LPI* lpi, const int param, int parval)
{
   int i;

   assert(lpi != NULL);

   for( i = 0; i < NUMINTPARAM; ++i )
      if( intparam[i] == param )
      {
         lpi->cpxparam.intparval[i] = parval;
         return;
      }

   errorMessage("Unknown CPLEX integer parameter\n");
   abort();
}

static
void setDblParam(LPI* lpi, const int param, double parval)
{
   int i;

   assert(lpi != NULL);

   for( i = 0; i < NUMDBLPARAM; ++i )
      if( dblparam[i] == param )
      {
         lpi->cpxparam.dblparval[i] = parval;
         return;
      }

   errorMessage("Unknown CPLEX double parameter\n");
   abort();
}

static
void invalidateSolution(LPI* lpi)
{
   assert(lpi != NULL);
   lpi->solstat = -1;
}

static
int cpxObjsen(OBJSEN objsen)
{
   switch( objsen )
   {
   case SCIP_OBJSEN_MAXIMIZE:
      return CPX_MAX;
   case SCIP_OBJSEN_MINIMIZE:
      return CPX_MIN;
   default:
      errorMessage("invalid objective sense\n");
      abort();
   }
}

/** converts SCIP's lhs/rhs pairs into CPLEX's sen/rhs/rng */
static
void convertSides(
   LPI*             lpi,                /**< LP interface structure */
   int              nrows,              /**< number of rows */
   const Real*      lhs,                /**< left hand side vector */
   const Real*      rhs,                /**< right hand side vector */
   int              indoffset,          /**< index of first row in LP */
   int*             rngcount            /**< pointer to store the number of range rows */
   )
{
   int i;

   assert(lpi != NULL);
   assert(nrows >= 0);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(rngcount != NULL);

   /* convert lhs/rhs into sen/rhs/rng */
   *rngcount = 0;
   for( i = 0; i < nrows; ++i )
   {
      assert(lhs[i] <= rhs[i]);
      if( lhs[i] == rhs[i] )
      {
         assert(-CPX_INFBOUND < rhs[i] && rhs[i] < CPX_INFBOUND);
         lpi->senarray[i] = 'E';
         lpi->rhsarray[i] = rhs[i];
      }
      else if( lhs[i] <= -CPX_INFBOUND )
      {
         assert(-CPX_INFBOUND < rhs[i] && rhs[i] < CPX_INFBOUND);
         lpi->senarray[i] = 'L';
         lpi->rhsarray[i] = rhs[i];
      }
      else if( rhs[i] >= CPX_INFBOUND )
      {
         assert(-CPX_INFBOUND < lhs[i] && lhs[i] < CPX_INFBOUND);
         lpi->senarray[i] = 'G';
         lpi->rhsarray[i] = lhs[i];
      }
      else
      {
         /* CPLEX defines a ranged row to be within rhs and rhs+rng.
          * -> To keep SCIP's meaning of the rhs value, we would like to use negative range values: rng := lhs - rng,
          *    but there seems to be a bug in CPLEX's presolve with negative range values:
          *    the ranged row
          *              0 <= -x <= 100000 with x >= 0 (rhs=0, rng=-100000) 
          *    would lead to the CPLEX row
          *              -x -Rg = 100000 
          *                  Rg = 0
          *    instead of the correct presolving implication  Rg = -100000.
          * -> Because of this bug, we have to use an additional rhsarray[] for the converted right hand sides and
          *    use rhsarray[i] = lhs[i] and rngarray[i] = rhs[i] - lhs[i] for ranged rows to keep the range values
          *    non-negative.
          */
         lpi->senarray[i] = 'R';
         lpi->rhsarray[i] = lhs[i];
         lpi->rngarray[*rngcount] = rhs[i] - lhs[i];
         lpi->rngindarray[*rngcount] = i + indoffset;
         (*rngcount)++;
      }
   }
}

/** converts CPLEX's sen/rhs/rng triplets into SCIP's lhs/rhs pairs */
static
void reconvertSides(
   LPI*             lpi,                /**< LP interface structure */
   int              nrows,              /**< number of rows */
   Real*            lhs,                /**< buffer to store the left hand side vector */
   Real*            rhs                 /**< buffer to store the right hand side vector */
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
         assert(lpi->rngarray[i] == 0.0);
         lhs[i] = lpi->rhsarray[i];
         rhs[i] = lpi->rhsarray[i];
         break;

      case 'L':
         assert(lpi->rngarray[i] == 0.0);
         lhs[i] = -CPX_INFBOUND;
         rhs[i] = lpi->rhsarray[i];
         break;

      case 'G':
         assert(lpi->rngarray[i] == 0.0);
         lhs[i] = lpi->rhsarray[i];
         rhs[i] = CPX_INFBOUND;
         break;

      case 'R':
         assert(lpi->rngarray[i] != 0.0);
         if( lpi->rngarray[i] > 0.0 )
         {
            lhs[i] = lpi->rhsarray[i];
            rhs[i] = lpi->rhsarray[i] + lpi->rngarray[i];
         }
         else
         {
            lhs[i] = lpi->rhsarray[i] + lpi->rngarray[i];
            rhs[i] = lpi->rhsarray[i];
         }
         break;
         
      default:
         errorMessage("invalid row sense\n");
         abort();
      }
      assert(lhs[i] <= rhs[i]);
   }
}




/*
 * LP Interface Methods
 */


/*
 * Miscellaneous Methods
 */

static char cpxname[MAXSTRLEN];

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   sprintf(cpxname, "CPLEX %.2f", (Real)CPX_VERSION/100.0);
   return cpxname;
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
   const char*      name                /**< problem name */
   )
{
   int     restat;

   assert(sizeof(Real) == sizeof(double)); /* CPLEX only works with doubles as floating points */
   assert(sizeof(Bool) == sizeof(int));    /* CPLEX only works with ints as bools */
   assert(lpi != NULL);
   assert(numlp >= 0);

   /* create environment */
   if( cpxenv == NULL )
   {
      assert(numlp == 0);
      cpxenv = CPXopenCPLEX(&restat);
      CHECK_ZERO( restat );

      /* get default parameter values */
      getParameterValues(&defparam);
      copyParameterValues(&actparam, &defparam);
   }
   assert(cpxenv != NULL);

   /* create LP */
   ALLOC_OKAY( allocMemory(lpi) );
   (*lpi)->larray = NULL;
   (*lpi)->uarray = NULL;
   (*lpi)->senarray = NULL;
   (*lpi)->rhsarray = NULL;
   (*lpi)->rngarray = NULL;
   (*lpi)->valarray = NULL;
   (*lpi)->rngindarray = NULL;
   (*lpi)->cstat = NULL;
   (*lpi)->rstat = NULL;
   (*lpi)->indarray = NULL;
   (*lpi)->boundchgsize = 0;
   (*lpi)->sidechgsize = 0;
   (*lpi)->valsize = 0;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->cpxlp = CPXcreateprob(cpxenv, &restat, name);
   CHECK_ZERO(restat);
   invalidateSolution(*lpi);
   copyParameterValues(&((*lpi)->cpxparam), &defparam);
   numlp++;

   return SCIP_OKAY;
}

/** deletes an LP problem object */
RETCODE SCIPlpiFree(
   LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(*lpi != NULL);

   /* free LP */
   CHECK_ZERO( CPXfreeprob(cpxenv, &((*lpi)->cpxlp)) );

   /* free memory */
   freeMemoryArrayNull(&(*lpi)->larray);
   freeMemoryArrayNull(&(*lpi)->uarray);
   freeMemoryArrayNull(&(*lpi)->senarray);
   freeMemoryArrayNull(&(*lpi)->rhsarray);
   freeMemoryArrayNull(&(*lpi)->rngarray);
   freeMemoryArrayNull(&(*lpi)->rngindarray);
   freeMemoryArrayNull(&(*lpi)->cstat);
   freeMemoryArrayNull(&(*lpi)->rstat);
   freeMemory(lpi);

   /* free environment */
   numlp--;
   if( numlp == 0 )
   {
      CHECK_ZERO( CPXcloseCPLEX(&cpxenv) );
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
   )
{
   int* cnt;
   int rngcount;
   int c;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_OKAY( ensureSidechgMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   convertSides(lpi, nrows, lhs, rhs, 0, &rngcount);

   /* calculate column lengths */
   ALLOC_OKAY( allocMemoryArray(&cnt, ncols) );
   for( c = 0; c < ncols-1; ++c )
   {
      cnt[c] = beg[c+1] - beg[c];
      assert(cnt[c] >= 0);
   }
   cnt[ncols-1] = nnonz - beg[ncols-1];
   assert(cnt[ncols-1] >= 0);

   /* copy data into CPLEX */
   CHECK_ZERO( CPXcopylpwnames(cpxenv, lpi->cpxlp, ncols, nrows, cpxObjsen(objsen), obj, 
                  lpi->rhsarray, lpi->senarray, beg, cnt, ind, val, lb, ub, lpi->rngarray, colnames, rownames) );

   /* free temporary memory */
   freeMemoryArray(&cnt);

   assert(CPXgetnumcols(cpxenv, lpi->cpxlp) == ncols);
   assert(CPXgetnumrows(cpxenv, lpi->cpxlp) == nrows);
   assert(CPXgetnumnz(cpxenv, lpi->cpxlp) == nnonz);

   return SCIP_OKAY;
}

/** adds columns to the LP */
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
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_ZERO( CPXaddcols(cpxenv, lpi->cpxlp, ncols, nnonz, obj, beg, ind, val, lb, ub, colnames) );

   return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
RETCODE SCIPlpiDelCols(
   LPI*             lpi,                /**< LP interface structure */
   int              firstcol,           /**< first column to be deleted */
   int              lastcol             /**< last column to be deleted */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < CPXgetnumcols(cpxenv, lpi->cpxlp));

   invalidateSolution(lpi);
   CHECK_ZERO( CPXdelcols(cpxenv, lpi->cpxlp, firstcol, lastcol) );

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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_ZERO( CPXdelsetcols(cpxenv, lpi->cpxlp, dstat) );

   return SCIP_OKAY;   
}

/** adds rows to the LP */
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
   )
{
   int rngcount;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_OKAY( ensureSidechgMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   convertSides(lpi, nrows, lhs, rhs, CPXgetnumrows(cpxenv, lpi->cpxlp), &rngcount);

   /* add rows to LP */
   CHECK_ZERO( CPXaddrows(cpxenv, lpi->cpxlp, 0, nrows, nnonz, lpi->rhsarray, lpi->senarray, beg, ind, val, NULL,
                  rownames) );
   if( rngcount > 0 )
   {
      CHECK_ZERO( CPXchgrngval(cpxenv, lpi->cpxlp, rngcount, lpi->rngindarray, lpi->rngarray) );
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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < CPXgetnumrows(cpxenv, lpi->cpxlp));

   invalidateSolution(lpi);
   CHECK_ZERO( CPXdelrows(cpxenv, lpi->cpxlp, firstrow, lastrow) );

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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_ZERO( CPXdelsetrows(cpxenv, lpi->cpxlp, dstat) );

   return SCIP_OKAY;   
}

/** clears the whole LP */
RETCODE SCIPlpiClear(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_ZERO( CPXdelcols(cpxenv, lpi->cpxlp, 0, CPXgetnumcols(cpxenv, lpi->cpxlp)-1) );
   CHECK_ZERO( CPXdelrows(cpxenv, lpi->cpxlp, 0, CPXgetnumrows(cpxenv, lpi->cpxlp)-1) );

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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_OKAY( ensureBoundchgMem(lpi, ncols) );

   CHECK_ZERO( CPXchgbds(cpxenv, lpi->cpxlp, ncols, ind, lpi->larray, (Real*)lb) );
   CHECK_ZERO( CPXchgbds(cpxenv, lpi->cpxlp, ncols, ind, lpi->uarray, (Real*)ub) );

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
   int rngcount;
   int i;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_OKAY( ensureSidechgMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   convertSides(lpi, nrows, lhs, rhs, 0, &rngcount);

   /* change row sides */
   CHECK_ZERO( CPXchgsense(cpxenv, lpi->cpxlp, nrows, ind, lpi->senarray) );
   CHECK_ZERO( CPXchgrhs(cpxenv, lpi->cpxlp, nrows, ind, lpi->rhsarray) );
   if( rngcount > 0 )
   {
      /* adjust the range count indices to the correct row indices */
      for( i = 0; i < rngcount; ++i )
      {
         assert(0 <= lpi->rngindarray[i] && lpi->rngindarray[i] < nrows);
         assert(lpi->senarray[i] == 'R');
         lpi->rngindarray[i] = ind[lpi->rngindarray[i]];
      }

      /* change the range values in CPLEX */
      CHECK_ZERO( CPXchgrngval(cpxenv, lpi->cpxlp, rngcount, lpi->rngindarray, lpi->rngarray) );
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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_ZERO( CPXchgcoef(cpxenv, lpi->cpxlp, row, col, newval) );

   return SCIP_OKAY;
}

/** changes the objective sense */
RETCODE SCIPlpiChgObjsen(
   LPI*             lpi,                /**< LP interface structure */
   OBJSEN           objsen              /**< new objective sense */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);
   
   CPXchgobjsen(cpxenv, lpi->cpxlp, cpxObjsen(objsen));

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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXchgobj(cpxenv, lpi->cpxlp, ncols, ind, obj) );

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
   int nnonz;
   int beg;
   int i;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(scaleval != 0.0);

   invalidateSolution(lpi);

   CHECK_OKAY( ensureValMem(lpi, CPXgetnumcols(cpxenv, lpi->cpxlp)) );

   /* get the row */
   CHECK_OKAY( SCIPlpiGetRows(lpi, row, row, &lhs, &rhs, &nnonz, &beg, lpi->indarray, lpi->valarray) );

   /* scale row coefficients */
   for( i = 0; i < nnonz; ++i )
   {
      CHECK_OKAY( SCIPlpiChgCoef(lpi, row, lpi->indarray[i], lpi->valarray[i] * scaleval) );
   }

   /* scale row sides */
   if( lhs > -CPX_INFBOUND )
      lhs *= scaleval;
   else if( scaleval < 0.0 )
      lhs = CPX_INFBOUND;
   if( rhs < CPX_INFBOUND )
      rhs *= scaleval;
   else if( scaleval < 0.0 )
      rhs = -CPX_INFBOUND;
   if( scaleval > 0.0 )
   {
      CHECK_OKAY( SCIPlpiChgSides(lpi, 1, &row, &lhs, &rhs) );
   }
   else
   {
      CHECK_OKAY( SCIPlpiChgSides(lpi, 1, &row, &rhs, &lhs) );
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
   Real lb;
   Real ub;
   Real obj;
   int nnonz;
   int beg;
   int i;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(scaleval != 0.0);

   invalidateSolution(lpi);

   CHECK_OKAY( ensureValMem(lpi, CPXgetnumcols(cpxenv, lpi->cpxlp)) );

   /* get the column */
   CHECK_OKAY( SCIPlpiGetCols(lpi, col, col, &lb, &ub, &nnonz, &beg, lpi->indarray, lpi->valarray) );

   /** get objective value */
   CHECK_OKAY( SCIPlpiGetObj(lpi, col, col, &obj) );

   /* scale column coefficients */
   for( i = 0; i < nnonz; ++i )
   {
      CHECK_OKAY( SCIPlpiChgCoef(lpi, lpi->indarray[i], col, lpi->valarray[i] * scaleval) );
   }

   /* scale objective value */
   obj *= scaleval;
   CHECK_OKAY( SCIPlpiChgObj(lpi, 1, &col, &obj) );

   /* scale column bounds */
   if( lb > -CPX_INFBOUND )
      lb /= scaleval;
   else if( scaleval < 0.0 )
      lb = CPX_INFBOUND;
   if( ub < CPX_INFBOUND )
      ub /= scaleval;
   else if( scaleval < 0.0 )
      ub = -CPX_INFBOUND;
   if( scaleval > 0.0 )
   {
      CHECK_OKAY( SCIPlpiChgBounds(lpi, 1, &col, &lb, &ub) );
   }
   else
   {
      CHECK_OKAY( SCIPlpiChgBounds(lpi, 1, &col, &ub, &lb) );
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
RETCODE SCIPlpiGetNRows(
   LPI*             lpi,                /**< LP interface structure */
   int*             nrows               /**< pointer to store the number of rows */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(nrows != NULL);

   *nrows = CPXgetnumrows(cpxenv, lpi->cpxlp);

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
RETCODE SCIPlpiGetNCols(
   LPI*             lpi,                /**< LP interface structure */
   int*             ncols               /**< pointer to store the number of cols */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(ncols != NULL);

   *ncols = CPXgetnumcols(cpxenv, lpi->cpxlp);

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
RETCODE SCIPlpiGetNNonz(
   LPI*             lpi,                /**< LP interface structure */
   int*             nnonz               /**< pointer to store the number of nonzeros */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(nnonz != NULL);

   *nnonz = CPXgetnumnz(cpxenv, lpi->cpxlp);

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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < CPXgetnumcols(cpxenv, lpi->cpxlp));

   if( lb != NULL )
   {
      assert(ub != NULL);

      CHECK_ZERO( CPXgetlb(cpxenv, lpi->cpxlp, lb, firstcol, lastcol) );
      CHECK_ZERO( CPXgetub(cpxenv, lpi->cpxlp, ub, firstcol, lastcol) );
   }
   else
      assert(ub == NULL);

   if( nnonz != NULL )
   {
      int surplus;

      assert(beg != NULL);
      assert(ind != NULL);
      assert(val != NULL);

      /* get matrix entries */
      CHECK_ZERO( CPXgetcols(cpxenv, lpi->cpxlp, nnonz, beg, ind, val, CPXgetnumnz(cpxenv, lpi->cpxlp), &surplus, 
                     firstcol, lastcol) );
      assert(surplus >= 0);
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
   int retcode;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < CPXgetnumrows(cpxenv, lpi->cpxlp));

   if( lhs != NULL )
   {
      assert(rhs != NULL);

      /* get row sense, rhs, and ranges */
      CHECK_OKAY( ensureSidechgMem(lpi, lastrow - firstrow + 1) );
      CHECK_ZERO( CPXgetsense(cpxenv, lpi->cpxlp, lpi->senarray, firstrow, lastrow) );
      CHECK_ZERO( CPXgetrhs(cpxenv, lpi->cpxlp, lpi->rhsarray, firstrow, lastrow) );
      retcode = CPXgetrngval(cpxenv, lpi->cpxlp, lpi->rngarray, firstrow, lastrow);
      if( retcode != CPXERR_NO_RNGVAL ) /* ignore "No range values" error */
      {
         CHECK_ZERO( retcode );
      }
      else
         clearMemoryArray(lpi->rngarray, lastrow-firstrow+1);

      /* convert sen/rhs/range into lhs/rhs tuples */
      reconvertSides(lpi, lastrow - firstrow + 1, lhs, rhs);
   }
   else
      assert(rhs == NULL);

   if( nnonz != NULL )
   {
      int surplus;

      assert(beg != NULL);
      assert(ind != NULL);
      assert(val != NULL);

      /* get matrix entries */
      CHECK_ZERO( CPXgetrows(cpxenv, lpi->cpxlp, nnonz, beg, ind, val, CPXgetnumnz(cpxenv, lpi->cpxlp), &surplus, 
                     firstrow, lastrow) );
      assert(surplus >= 0);
   }
   else
   {
      assert(beg == NULL);
      assert(ind == NULL);
      assert(val == NULL);
   }

   return SCIP_OKAY;
}

/** gets objective values from LP problem object */
RETCODE SCIPlpiGetObj(
   LPI*             lpi,                /**< LP interface structure */
   int              firstcol,           /**< first column to get objective value for */
   int              lastcol,            /**< last column to get objective value for */
   Real*            vals                /**< array to store objective values */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(firstcol <= lastcol);
   assert(vals != NULL);
   
   CHECK_ZERO( CPXgetobj(cpxenv, lpi->cpxlp, vals, firstcol, lastcol) );

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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXgetcoef(cpxenv, lpi->cpxlp, row, col, val) );

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
   int retval;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   debugMessage("calling CPLEX primal simplex: %d cols, %d rows\n",
      CPXgetnumcols(cpxenv, lpi->cpxlp), CPXgetnumrows(cpxenv, lpi->cpxlp));

   invalidateSolution(lpi);

   CHECK_OKAY( setParameterValues(&(lpi->cpxparam)) );

   retval = CPXprimopt(cpxenv, lpi->cpxlp);
   switch( retval  )
   {
   case 0:
      break;
   case CPXERR_NO_MEMORY:
      return SCIP_NOMEMORY;
   default:
      return SCIP_LPERROR;
   }

   lpi->solstat = CPXgetstat(cpxenv, lpi->cpxlp);
   debugMessage(" -> CPLEX returned solstat=%d\n", lpi->solstat);

   if( lpi->solstat == CPX_STAT_INForUNBD )
   {
      /* maybe the preprocessor solved the problem; but we need a solution, so solve again without preprocessing */
      debugMessage("CPLEX returned INForUNBD -> calling CPLEX primal simplex again without presolve\n");
      
      /* switch off preprocessing */
      setIntParam(lpi, CPX_PARAM_PREIND, CPX_OFF);
      CHECK_OKAY( setParameterValues(&(lpi->cpxparam)) );

      retval = CPXprimopt(cpxenv, lpi->cpxlp);
      switch( retval  )
      {
      case 0:
         break;
      case CPXERR_NO_MEMORY:
         return SCIP_NOMEMORY;
      default:
         return SCIP_LPERROR;
      }

      lpi->solstat = CPXgetstat(cpxenv, lpi->cpxlp);
      debugMessage(" -> CPLEX returned solstat=%d\n", lpi->solstat);

      if( lpi->solstat == CPX_STAT_INForUNBD )
      {
         /* preprocessing was not the problem; issue a warning message and treat LP as infeasible */
         errorMessage("CPLEX primal simplex returned CPX_STAT_INForUNBD after presolving was turned off\n");
      }

      /* switch on preprocessing again */
      setIntParam(lpi, CPX_PARAM_PREIND, CPX_ON);
   }

   return SCIP_OKAY;
}

/** calls dual simplex to solve the LP */
RETCODE SCIPlpiSolveDual(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   int retval;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   debugMessage("calling CPLEX dual simplex: %d cols, %d rows\n", 
      CPXgetnumcols(cpxenv, lpi->cpxlp), CPXgetnumrows(cpxenv, lpi->cpxlp));

   invalidateSolution(lpi);

   CHECK_OKAY( setParameterValues(&(lpi->cpxparam)) );

   retval = CPXdualopt(cpxenv, lpi->cpxlp);
   switch( retval  )
   {
   case 0:
      break;
   case CPXERR_NO_MEMORY:
      return SCIP_NOMEMORY;
   default:
      return SCIP_LPERROR;
   }

   lpi->solstat = CPXgetstat(cpxenv, lpi->cpxlp);
   debugMessage(" -> CPLEX returned solstat=%d\n", lpi->solstat);

   if( lpi->solstat == CPX_STAT_INForUNBD )
   {
      /* maybe the preprocessor solved the problem; but we need a solution, so solve again without preprocessing */
      debugMessage("CPLEX returned INForUNBD -> calling CPLEX dual simplex again without presolve\n");
      
      /* switch off preprocessing */
      setIntParam(lpi, CPX_PARAM_PREIND, CPX_OFF);
      CHECK_OKAY( setParameterValues(&(lpi->cpxparam)) );

      retval = CPXdualopt(cpxenv, lpi->cpxlp);
      switch( retval  )
      {
      case 0:
         break;
      case CPXERR_NO_MEMORY:
         return SCIP_NOMEMORY;
      default:
         return SCIP_LPERROR;
      }

      lpi->solstat = CPXgetstat(cpxenv, lpi->cpxlp);
      debugMessage(" -> CPLEX returned solstat=%d\n", lpi->solstat);

      if( lpi->solstat == CPX_STAT_INForUNBD )
      {
         /* preprocessing was not the problem; issue a warning message and treat LP as infeasible */
         errorMessage("CPLEX dual simplex returned CPX_STAT_INForUNBD after presolving was turned off\n");
      }

      /* switch on preprocessing again */
      setIntParam(lpi, CPX_PARAM_PREIND, CPX_ON);
   }

   return SCIP_OKAY;
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
RETCODE SCIPlpiSolveBarrier(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   int retval;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   debugMessage("calling CPLEX barrier: %d cols, %d rows\n",
      CPXgetnumcols(cpxenv, lpi->cpxlp), CPXgetnumrows(cpxenv, lpi->cpxlp));

   invalidateSolution(lpi);

   setParameterValues(&(lpi->cpxparam));

   retval = CPXhybbaropt(cpxenv, lpi->cpxlp, 0);
   switch( retval  )
   {
   case 0:
      break;
   case CPXERR_NO_MEMORY:
      return SCIP_NOMEMORY;
   default:
      return SCIP_LPERROR;
   }

   lpi->solstat = CPXgetstat(cpxenv, lpi->cpxlp);
   debugMessage(" -> CPLEX returned solstat=%d\n", lpi->solstat);

   if( lpi->solstat == CPX_STAT_INForUNBD )
   {
      /* maybe the preprocessor solved the problem; but we need a solution, so solve again without preprocessing */
      debugMessage("CPLEX returned INForUNBD -> calling CPLEX barrier again without presolve\n");
      
      CHECK_ZERO( CPXsetintparam(cpxenv, CPX_PARAM_PREIND, CPX_OFF) );
      retval = CPXhybbaropt(cpxenv, lpi->cpxlp, 0);
      switch( retval  )
      {
      case 0:
         break;
      case CPXERR_NO_MEMORY:
         return SCIP_NOMEMORY;
      default:
         return SCIP_LPERROR;
      }

      lpi->solstat = CPXgetstat(cpxenv, lpi->cpxlp);
      debugMessage(" -> CPLEX returned solstat=%d\n", lpi->solstat);

      if( lpi->solstat == CPX_STAT_INForUNBD )
      {
         /* preprocessing was not the problem; issue a warning message and treat LP as infeasible */
         errorMessage("CPLEX barrier returned CPX_STAT_INForUNBD after presolving was turned off\n");
      }

      CHECK_ZERO( CPXsetintparam(cpxenv, CPX_PARAM_PREIND, CPX_ON) );
   }

   return SCIP_OKAY;
}

/** performs strong branching iterations on all candidates */
RETCODE SCIPlpiStrongbranch(
   LPI*             lpi,                /**< LP interface structure */
   const int*       cand,               /**< candidate list */
   Real*            psol,               /**< array with current primal solution values of candidates */
   int              ncand,              /**< size of candidate list */
   int              itlim,              /**< iteration limit for strong branchings */
   Real*            down,               /**< stores dual bound after branching candidate down */
   Real*            up,                 /**< stores dual bound after branching candidate up */
   int*             iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   
   CHECK_ZERO( CPXstrongbranch(cpxenv, lpi->cpxlp, cand, ncand, down, up, itlim) );

   /* CPLEX is not able to return the iteration counts in strong branching */
   if( iter != NULL )
      *iter = -1;

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
   int pfeas;
   int dfeas;

   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(primalfeasible != NULL);
   assert(dualfeasible != NULL);

   CHECK_ZERO( CPXsolninfo(cpxenv, lpi->cpxlp, NULL, NULL, &pfeas, &dfeas) );
   *primalfeasible = (Bool)pfeas;
   *dualfeasible = (Bool)dfeas;

   return SCIP_OKAY;
}

/** returns TRUE iff LP is primal unbounded */
Bool SCIPlpiIsPrimalUnbounded(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   int primalfeasible;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpi->solstat >= 0);

   ABORT_ZERO( CPXsolninfo(cpxenv, lpi->cpxlp, NULL, NULL, &primalfeasible, NULL) );
   
   return (primalfeasible && (lpi->solstat == CPX_STAT_UNBOUNDED || lpi->solstat == CPX_STAT_INForUNBD));
}

/** returns TRUE iff LP is primal infeasible */
Bool SCIPlpiIsPrimalInfeasible(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   int primalfeasible;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpi->solstat >= 0);

   ABORT_ZERO( CPXsolninfo(cpxenv, lpi->cpxlp, NULL, NULL, &primalfeasible, NULL) );

   return (lpi->solstat == CPX_STAT_INFEASIBLE
      || (!primalfeasible && (lpi->solstat == CPX_STAT_UNBOUNDED || lpi->solstat == CPX_STAT_INForUNBD)));
}

/** returns TRUE iff LP is dual unbounded */
Bool SCIPlpiIsDualUnbounded(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   int dualfeasible;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpi->solstat >= 0);

   ABORT_ZERO( CPXsolninfo(cpxenv, lpi->cpxlp, NULL, NULL, NULL, &dualfeasible) );

   return (lpi->solstat == CPX_STAT_INFEASIBLE || (lpi->solstat == CPX_STAT_INForUNBD && dualfeasible));
}

/** returns TRUE iff LP is dual infeasible */
Bool SCIPlpiIsDualInfeasible(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   int dualfeasible;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpi->solstat >= 0);

   ABORT_ZERO( CPXsolninfo(cpxenv, lpi->cpxlp, NULL, NULL, NULL, &dualfeasible) );

   return (lpi->solstat == CPX_STAT_UNBOUNDED || (lpi->solstat == CPX_STAT_INForUNBD && !dualfeasible));
}

/** returns TRUE iff LP was solved to optimality */
Bool SCIPlpiIsOptimal(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == CPX_STAT_OPTIMAL);
}

/** returns TRUE iff actual LP basis is stable */
Bool SCIPlpiIsStable(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat != CPX_STAT_NUM_BEST && lpi->solstat != CPX_STAT_OPTIMAL_INFEAS);
}

/** returns TRUE iff the objective limit was reached */
Bool SCIPlpiIsObjlimExc(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == CPX_STAT_ABORT_OBJ_LIM);
}

/** returns TRUE iff the iteration limit was reached */
Bool SCIPlpiIsIterlimExc(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == CPX_STAT_ABORT_IT_LIM);
}

/** returns TRUE iff the time limit was reached */
Bool SCIPlpiIsTimelimExc(
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == CPX_STAT_ABORT_TIME_LIM);
}

/** gets objective value of solution */
RETCODE SCIPlpiGetObjval(
   LPI*             lpi,                /**< LP interface structure */
   Real*            objval              /**< stores the objective value */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXgetobjval(cpxenv, lpi->cpxlp, objval) );

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
   int dummy;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpi->solstat >= 0);

   CHECK_ZERO( CPXsolution(cpxenv, lpi->cpxlp, &dummy, objval, primsol, dualsol, NULL, redcost) );
   assert(dummy == lpi->solstat);

   if( activity != NULL )
   {
      CHECK_ZERO( CPXgetax(cpxenv, lpi->cpxlp, activity, 0, CPXgetnumrows(cpxenv, lpi->cpxlp)-1) );
   }

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
RETCODE SCIPlpiGetPrimalRay(
   LPI*             lpi,                /**< LP interface structure */
   Real*            ray                 /**< primal ray */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpi->solstat >= 0);

   debugMessage("calling CPLEX get primal ray: %d cols, %d rows\n",
      CPXgetnumcols(cpxenv, lpi->cpxlp), CPXgetnumrows(cpxenv, lpi->cpxlp));

   CHECK_ZERO( CPXgetray(cpxenv, lpi->cpxlp, ray) );

   return SCIP_OKAY;
}

/** gets dual farkas proof for infeasibility */
RETCODE SCIPlpiGetDualfarkas(
   LPI*             lpi,                /**< LP interface structure */
   Real*            dualfarkas          /**< dual farkas row multipliers */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpi->solstat >= 0);
   assert(dualfarkas != NULL);

   debugMessage("calling CPLEX dual farkas: %d cols, %d rows\n",
      CPXgetnumcols(cpxenv, lpi->cpxlp), CPXgetnumrows(cpxenv, lpi->cpxlp));

   CHECK_ZERO( CPXdualfarkas(cpxenv, lpi->cpxlp, dualfarkas, NULL) );

   return SCIP_OKAY;
}

/**@} */




/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */

/** gets actual basis status for columns and rows; arrays must be large enough to store the basis status */
RETCODE SCIPlpiGetBase(
   LPI*             lpi,                /**< LP interface structure */
   int*             cstat,              /**< array to store column basis status, or NULL */
   int*             rstat               /**< array to store row basis status, or NULL */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXgetbase(cpxenv, lpi->cpxlp, cstat, rstat) );

   /* because the basis status values are equally defined in SCIP and CPLEX, they don't need to be transformed */
   assert(SCIP_BASESTAT_LOWER == CPX_AT_LOWER);
   assert(SCIP_BASESTAT_BASIC == CPX_BASIC);
   assert(SCIP_BASESTAT_UPPER == CPX_AT_UPPER);

   return SCIP_OKAY;
}

/** sets actual basis status for columns and rows */
RETCODE SCIPlpiSetBase(
   LPI*             lpi,                /**< LP interface structure */
   int*             cstat,              /**< array with column basis status */
   int*             rstat               /**< array with row basis status */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(cstat != NULL);
   assert(rstat != NULL);

   invalidateSolution(lpi);

   /* because the basis status values are equally defined in SCIP and CPLEX, they don't need to be transformed */
   assert(SCIP_BASESTAT_LOWER == CPX_AT_LOWER);
   assert(SCIP_BASESTAT_BASIC == CPX_BASIC);
   assert(SCIP_BASESTAT_UPPER == CPX_AT_UPPER);

   CHECK_ZERO( CPXcopybase(cpxenv, lpi->cpxlp, cstat, rstat) );

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows */
RETCODE SCIPlpiGetBasisInd(
   LPI*             lpi,                /**< LP interface structure */
   int*             bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXgetbhead(cpxenv, lpi->cpxlp, bind, NULL) );

   return SCIP_OKAY;
}

/** get dense row of inverse basis matrix B^-1 */
RETCODE SCIPlpiGetBInvRow(
   LPI*             lpi,                /**< LP interface structure */
   int              r,                  /**< row number */
   Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXbinvrow(cpxenv, lpi->cpxlp, r, coef) );

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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXbinvarow(cpxenv, lpi->cpxlp, r, val) );

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

   assert(memhdr != NULL);
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpistate != NULL);

   ncols = CPXgetnumcols(cpxenv, lpi->cpxlp);
   nrows = CPXgetnumrows(cpxenv, lpi->cpxlp);
   assert(ncols >= 0);
   assert(nrows >= 0);
   
   /* allocate lpistate data */
   CHECK_OKAY( lpistateCreate(lpistate, memhdr, ncols, nrows) );

   /* allocate enough memory for storing uncompressed basis information */
   CHECK_OKAY( ensureCstatMem(lpi, ncols) );
   CHECK_OKAY( ensureRstatMem(lpi, nrows) );

   /* get unpacked basis information from CPLEX */
   if( getIntParam(lpi, CPX_PARAM_DPRIIND) == CPX_DPRIIND_STEEP )
   {
      ALLOC_OKAY( allocBlockMemoryArray(memhdr, &(*lpistate)->dnorm, ncols) );
      CHECK_ZERO( CPXgetbasednorms(cpxenv, lpi->cpxlp, lpi->cstat, lpi->rstat, (*lpistate)->dnorm) );
   }
   else
   {
      (*lpistate)->dnorm = NULL;
      CHECK_ZERO( CPXgetbase(cpxenv, lpi->cpxlp, lpi->cstat, lpi->rstat) );
   }

   /* pack LPi state data */
   (*lpistate)->ncols = ncols;
   (*lpistate)->nrows = nrows;
   lpistatePack(*lpistate, lpi->cstat, lpi->rstat);

   return SCIP_OKAY;
}

/** loads LPi state (like basis information) into solver */
RETCODE SCIPlpiSetState(
   LPI*             lpi,                /**< LP interface structure */
   MEMHDR*          memhdr,             /**< block memory */
   LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   )
{
   assert(memhdr != NULL);
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpistate != NULL);
   assert(lpistate->ncols == CPXgetnumcols(cpxenv, lpi->cpxlp));
   assert(lpistate->nrows == CPXgetnumrows(cpxenv, lpi->cpxlp));

   /* allocate enough memory for storing uncompressed basis information */
   CHECK_OKAY( ensureCstatMem(lpi, lpistate->ncols) );
   CHECK_OKAY( ensureRstatMem(lpi, lpistate->nrows) );

   /* unpack LPi state data */
   lpistateUnpack(lpistate, lpi->cstat, lpi->rstat);

   /* load basis information into CPLEX */
   if( lpistate->dnorm != NULL && getIntParam(lpi, CPX_PARAM_DPRIIND) == CPX_DPRIIND_STEEP )
   {
      CHECK_ZERO( CPXcopybasednorms(cpxenv, lpi->cpxlp, lpi->cstat, lpi->rstat, lpistate->dnorm) );
   }
   else
   {
      CHECK_ZERO( CPXcopybase(cpxenv, lpi->cpxlp, lpi->cstat, lpi->rstat) );
   }

   return SCIP_OKAY;
}

/** frees LPi state information */
RETCODE SCIPlpiFreeState(
   LPI*             lpi,                /**< LP interface structure */
   MEMHDR*          memhdr,             /**< block memory */
   LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   assert(lpi != NULL);

   lpistateFree(lpistate, memhdr);

   return SCIP_OKAY;
}

/** reads LP state (like basis information from a file */
RETCODE SCIPlpiReadState(
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXreadcopybase(cpxenv, lpi->cpxlp, fname) );

   return SCIP_OKAY;
}

/** writes LP state (like basis information) to a file */
RETCODE SCIPlpiWriteState(
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXmbasewrite(cpxenv, lpi->cpxlp, fname) );

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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(ival != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      *ival = (getIntParam(lpi, CPX_PARAM_ADVIND) == CPX_OFF);
      break;
   case SCIP_LPPAR_FASTMIP:
      *ival = (getIntParam(lpi, CPX_PARAM_FASTMIP) == CPX_ON);
      break;
   case SCIP_LPPAR_SCALING:
      *ival = (getIntParam(lpi, CPX_PARAM_SCAIND) == 0);
      break;
   case SCIP_LPPAR_PRICING:
      switch( getIntParam(lpi, CPX_PARAM_DPRIIND) )
      {
      case CPX_DPRIIND_FULL:
         *ival = SCIP_PRICING_FULL;
         break;
      case CPX_DPRIIND_STEEP:
         *ival = SCIP_PRICING_STEEP;
         break;
      case CPX_DPRIIND_STEEPQSTART:
         *ival = SCIP_PRICING_STEEPQSTART;
         break;
      default:
         *ival = SCIP_PRICING_AUTO;
         break;
      }
      break;
   case SCIP_LPPAR_LPINFO:
      *ival = (getIntParam(lpi, CPX_PARAM_SCRIND) == CPX_ON);
      break;
   case SCIP_LPPAR_LPITLIM:
      *ival = getIntParam(lpi, CPX_PARAM_ITLIM);
      break;
   case SCIP_LPPAR_LPITER:
      *ival = CPXgetphase1cnt(cpxenv, lpi->cpxlp) + CPXgetitcnt(cpxenv, lpi->cpxlp);
      break;
   default:
      return SCIP_LPERROR;
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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      assert(ival == TRUE || ival == FALSE);
      setIntParam(lpi, CPX_PARAM_ADVIND, (ival == FALSE) ? CPX_ON : CPX_OFF);
      break;
   case SCIP_LPPAR_FASTMIP:
      assert(ival == TRUE || ival == FALSE);
      setIntParam(lpi, CPX_PARAM_FASTMIP, (ival == TRUE) ? CPX_ON : CPX_OFF);
      break;
   case SCIP_LPPAR_SCALING:
      assert(ival == TRUE || ival == FALSE);
      setIntParam(lpi, CPX_PARAM_SCAIND, (ival == TRUE) ? 0 : -1);
      break;
   case SCIP_LPPAR_PRICING:
      switch( (PRICING)ival )
      {
      case SCIP_PRICING_AUTO:
	 setIntParam(lpi, CPX_PARAM_PPRIIND, CPX_PPRIIND_AUTO);
	 setIntParam(lpi, CPX_PARAM_DPRIIND, CPX_DPRIIND_AUTO);
         break;
      case SCIP_PRICING_FULL:
	 setIntParam(lpi, CPX_PARAM_PPRIIND, CPX_PPRIIND_FULL);
	 setIntParam(lpi, CPX_PARAM_DPRIIND, CPX_DPRIIND_FULL);
         break;
      case SCIP_PRICING_STEEP:
	 setIntParam(lpi, CPX_PARAM_PPRIIND, CPX_PPRIIND_STEEP);
	 setIntParam(lpi, CPX_PARAM_DPRIIND, CPX_DPRIIND_STEEP);
	 break;
      case SCIP_PRICING_STEEPQSTART:
	 setIntParam(lpi, CPX_PARAM_PPRIIND, CPX_PPRIIND_STEEPQSTART);
	 setIntParam(lpi, CPX_PARAM_DPRIIND, CPX_DPRIIND_STEEPQSTART);
	 break;
      default:
         return SCIP_LPERROR;
      }
      break;
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      if( ival )
	 setIntParam(lpi, CPX_PARAM_SCRIND, CPX_ON);
      else 
	 setIntParam(lpi, CPX_PARAM_SCRIND, CPX_OFF);
      break;
   case SCIP_LPPAR_LPITLIM:
      setIntParam(lpi, CPX_PARAM_ITLIM, ival);
      break;
   default:
      return SCIP_LPERROR;
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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(dval != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      *dval = getDblParam(lpi, CPX_PARAM_EPRHS);
      break;
   case SCIP_LPPAR_LOBJLIM:
      *dval = getDblParam(lpi, CPX_PARAM_OBJLLIM);
      break;
   case SCIP_LPPAR_UOBJLIM:
      *dval = getDblParam(lpi, CPX_PARAM_OBJULIM);
      break;
   case SCIP_LPPAR_LPTILIM:
      *dval = getDblParam(lpi, CPX_PARAM_TILIM);
      break;
   default:
      return SCIP_LPERROR;
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
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      setDblParam(lpi, CPX_PARAM_EPRHS, dval);
      setDblParam(lpi, CPX_PARAM_EPOPT, dval);
      break;
   case SCIP_LPPAR_LOBJLIM:
      setDblParam(lpi, CPX_PARAM_OBJLLIM, dval);
      break;
   case SCIP_LPPAR_UOBJLIM:
      setDblParam(lpi, CPX_PARAM_OBJULIM, dval);
      break;
   case SCIP_LPPAR_LPTILIM:
      setDblParam(lpi, CPX_PARAM_TILIM, dval);
      break;
   default:
      return SCIP_LPERROR;
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
   LPI*             lpi                 /**< LP interface structure */
   )
{
   return CPX_INFBOUND;
}

/** checks if given value is treated as infinity in the LP solver */
Bool SCIPlpiIsInfinity(
   LPI*             lpi,                /**< LP interface structure */
   Real             val
   )
{
   return (val >= CPX_INFBOUND);
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
RETCODE SCIPlpiReadLP(
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXreadcopyprob(cpxenv, lpi->cpxlp, fname, NULL) );

   return SCIP_OKAY;
}

/** writes LP to a file */
RETCODE SCIPlpiWriteLP(
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXwriteprob(cpxenv, lpi->cpxlp, fname, NULL) );

   return SCIP_OKAY;
}

/**@} */

