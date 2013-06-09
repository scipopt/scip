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

/**@file   lpi_grb.c
 * @ingroup LPIS
 * @brief  LP interface for Gurobi
 * @author Marc Pfetsch
 *
 * This LPI is beta!
 *
 * Several things are missing in the Gurobi interface that make this LPI relatively useless:
 *
 * - Gurobi currently does not allow to access the basis inverse.
 * - Strong branching is supported, but not documented.
 * - The support of ranged rows is complicated for the user: one has to keep track of the additional
 *   variables, which are added to generate a ranged row. Hence, one would need to adapt the count
 *   of variables and retrieve the information of ranged rows to get the correct answers.
 *
 * While the first two issues only influence the performance, the third is critical for some
 * problems, which contain ranged rows.
 *
 * @todo Check whether functions for basis inverses are correct. Which ones are the right ones?
 *
 * @todo Check whether solisbasic is correctly used.
 *
 * @todo Try quad-precision and concurrent runs.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "gurobi_c.h"
#include "scip/lpi.h"
#include "scip/pub_message.h"


#define CHECK_ZERO(messagehdlr, x) { int _restat_;                      \
      if( (_restat_ = (x)) != 0 )                                       \
      {                                                                 \
         SCIPmessagePrintWarning((messagehdlr), "Gurobi error %d: %s\n", _restat_, GRBgeterrormsg(grbenv)); \
         return SCIP_LPERROR;                                           \
      }                                                                 \
   }

#if( GRB_VERSION_MAJOR < 4 )
#define GRB_METHOD_DUAL    GRB_LPMETHOD_DUAL
#define GRB_METHOD_PRIMAL  GRB_LPMETHOD_PRIMAL
#define GRB_INT_PAR_METHOD GRB_INT_PAR_LPMETHOD
#endif

typedef unsigned int SCIP_SINGLEPACKET;                /**< storing single bits in packed format */
#define SCIP_SINGLEPACKETSIZE (sizeof(SCIP_SINGLEPACKET)*8) /**< each entry needs one bit of information */
typedef unsigned int SCIP_DUALPACKET;                  /**< storing bit pairs in packed format */
#define SCIP_DUALPACKETSIZE   (sizeof(SCIP_DUALPACKET)*4)   /**< each entry needs two bits of information */

typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE


/* Gurobi parameter lists which can be changed */
#define NUMINTPARAM 4

static const char* intparam[NUMINTPARAM] =
{
   GRB_INT_PAR_SCALEFLAG,
   GRB_INT_PAR_PRESOLVE,
   GRB_INT_PAR_SIMPLEXPRICING,
   GRB_INT_PAR_OUTPUTFLAG
};

#define NUMDBLPARAM 6

static const char* dblparam[NUMDBLPARAM] =
{
   GRB_DBL_PAR_FEASIBILITYTOL,
   GRB_DBL_PAR_OPTIMALITYTOL,
   GRB_DBL_PAR_CUTOFF,
   GRB_DBL_PAR_TIMELIMIT,
   GRB_DBL_PAR_ITERATIONLIMIT,
   GRB_DBL_PAR_MARKOWITZTOL
};

static const double dblparammin[NUMDBLPARAM] =
{
   +1e-09,               /* GRB_DBL_PAR_FEASIBILITYTOL */
   +1e-09,               /* GRB_DBL_PAR_OPTIMALITYTOL */
   -GRB_INFINITY,        /* GRB_DBL_PAR_CUTOFF */
   0,                    /* GRB_DBL_PAR_TIMELIMIT */
   0,                    /* GRB_DBL_PAR_ITERATIONLIMIT */
   1e-04                 /* GRB_DBL_PAR_MARKOWITZTOL */
};

/** Gurobi parameter settings */
struct GRBParam
{
   int                   intparval[NUMINTPARAM]; /**< integer parameter values */
   double                dblparval[NUMDBLPARAM]; /**< double parameter values */
};
typedef struct GRBParam GRBPARAM;


/** LP interface */
struct SCIP_LPi
{
   GRBmodel*             grbmodel;           /**< Gurobi model pointer */
   GRBenv*               grbenv;             /**< environment corresponding to model */
   int                   solstat;            /**< solution status of last optimization call */
   GRBPARAM              defparam;           /**< default parameter values */
   GRBPARAM              curparam;           /**< current parameter values stored in Gurobi LP */
   GRBPARAM              grbparam;           /**< current parameter values for this LP */
   char*                 senarray;           /**< array for storing row senses */
   SCIP_Real*            rhsarray;           /**< array for storing rhs values */
   SCIP_Real*            valarray;           /**< array for storing coefficient values */
   int*                  cstat;              /**< array for storing column basis status */
   int*                  rstat;              /**< array for storing row basis status */
   int*                  indarray;           /**< array for storing coefficient indices */
   int                   sidechgsize;        /**< size of senarray */
   int                   valsize;            /**< size of valarray and indarray */
   int                   cstatsize;          /**< size of cstat array */
   int                   rstatsize;          /**< size of rstat array */
   int                   iterations;         /**< number of iterations used in the last solving call */
   SCIP_Bool             solisbasic;         /**< is current LP solution a basic solution? */
   SCIP_PRICING          pricing;            /**< SCIP pricing setting  */
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


static GRBenv*           grbenv = NULL;      /**< Gurobi environment (only needed for initialization) */
static int               numlp = 0;          /**< number of open LP objects */



/*
 * dynamic memory arrays
 */

/** resizes senarray to have at least num entries */
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
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< whether basis information has successfully been obtained */
   )
{
   int ncols;
   int nrows;
   int res;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );

   SCIPdebugMessage("getBase()\n");
   if ( success != NULL )
      *success = TRUE;

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, ncols) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   /* get unpacked basis information from Gurobi */
   res = GRBgetintattrarray(lpi->grbmodel, GRB_INT_ATTR_VBASIS, 0, ncols, lpi->cstat);
   if ( res == GRB_ERROR_DATA_NOT_AVAILABLE )
   {
      /* if the model is infeasible Gurobi does not currently return basis information */
      if ( success != NULL )
         *success = FALSE;
      return SCIP_OKAY;
   }
   else if ( res != 0 )
   {
      SCIPerrorMessage("Gurobi error %d: %s\n", res, GRBgeterrormsg(lpi->grbenv));
      return SCIP_LPERROR;
   }

   res = GRBgetintattrarray(lpi->grbmodel, GRB_INT_ATTR_CBASIS, 0, nrows, lpi->rstat);
   if ( res == GRB_ERROR_DATA_NOT_AVAILABLE )
   {
      /* if the model is infeasible Gurobi does not currently return basis information */
      if ( success != NULL )
         *success = FALSE;
      return SCIP_OKAY;
   }
   else if ( res != 0 )
   {
      SCIPerrorMessage("Gurobi error %d: %s\n", res, GRBgeterrormsg(lpi->grbenv));
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** loads basis stored in internal arrays of LPI data structure into Gurobi */
static
SCIP_RETCODE setBase(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int ncols;
   int nrows;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );

   SCIPdebugMessage("setBase()\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );

   /* load basis information into Gurobi */
   CHECK_ZERO( lpi->messagehdlr, GRBsetintattrarray(lpi->grbmodel, GRB_INT_ATTR_VBASIS, 0, ncols, lpi->cstat) );
   CHECK_ZERO( lpi->messagehdlr, GRBsetintattrarray(lpi->grbmodel, GRB_INT_ATTR_CBASIS, 0, nrows, lpi->rstat) );

   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

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


/* The basis information for Gurobi is negative. So we cannot use the functions in bitencode.h/c. The functions below are a modified copy. */

/** encode a negated dual bit vector into packed format */
static
void SCIPencodeDualBitNeg(
   const int*            inp,                /**< unpacked input vector */
   SCIP_DUALPACKET*      out,                /**< buffer to store the packed vector */
   int                   count               /**< number of elements */
   )
{
   static const SCIP_DUALPACKET mask[SCIP_DUALPACKETSIZE][4] = {   /* if the packet size changes, the mask has to be updated */
      {0x00000000, 0x00000001, 0x00000002, 0x00000003},
      {0x00000000, 0x00000004, 0x00000008, 0x0000000C},
      {0x00000000, 0x00000010, 0x00000020, 0x00000030},
      {0x00000000, 0x00000040, 0x00000080, 0x000000C0},
      {0x00000000, 0x00000100, 0x00000200, 0x00000300},
      {0x00000000, 0x00000400, 0x00000800, 0x00000C00},
      {0x00000000, 0x00001000, 0x00002000, 0x00003000},
      {0x00000000, 0x00004000, 0x00008000, 0x0000C000},
      {0x00000000, 0x00010000, 0x00020000, 0x00030000},
      {0x00000000, 0x00040000, 0x00080000, 0x000C0000},
      {0x00000000, 0x00100000, 0x00200000, 0x00300000},
      {0x00000000, 0x00400000, 0x00800000, 0x00C00000},
      {0x00000000, 0x01000000, 0x02000000, 0x03000000},
      {0x00000000, 0x04000000, 0x08000000, 0x0C000000},
      {0x00000000, 0x10000000, 0x20000000, 0x30000000},
      {0x00000000, 0x40000000, 0x80000000, 0xC0000000}
   };
   int i;
   int rest;
   int nfull;

   assert(inp != NULL || count == 0);
   assert(out != NULL || count == 0);
   assert(count >= 0);
   assert(SCIP_DUALPACKETSIZE == 16);

   rest = count % (int)SCIP_DUALPACKETSIZE;
   nfull = count - rest;

   for( i = 0; i < nfull; i += (int)SCIP_DUALPACKETSIZE, inp += (int)SCIP_DUALPACKETSIZE )
   {
      assert(inp != NULL);
      assert(out != NULL);

#ifndef NDEBUG
      {
         unsigned int j;
         for( j = 0; j < SCIP_DUALPACKETSIZE; ++j )
            assert(0 <= -inp[j] && -inp[j] <= 3);
      }
#endif
      *out++ =
         mask[0][-inp[0]] | mask[1][-inp[1]] | mask[2][-inp[2]] | mask[3][inp[3]]
         | mask[4][-inp[4]] | mask[5][-inp[5]] | mask[6][-inp[6]]
         | mask[7][-inp[7]] | mask[8][-inp[8]] | mask[9][-inp[9]]
         | mask[10][-inp[10]] | mask[11][-inp[11]] | mask[12][-inp[12]]
         | mask[13][-inp[13]] | mask[14][-inp[14]] | mask[15][-inp[15]];
   }

   if( rest > 0 )
   {
      SCIP_DUALPACKET m = (SCIP_DUALPACKET) 0u;

      assert(inp != NULL);
      assert(out != NULL);

      for( i = 0; i < rest; i++ )
         m |= mask[i][-inp[i]];
      *out = m;
   }
}

/** decode a packed dual bit vector into negated unpacked format */
static
void SCIPdecodeDualBitNeg(
   const SCIP_DUALPACKET* inp,               /**< packed input vector */
   int*                  out,                /**< buffer to store unpacked vector */
   int                   count               /**< number of elements */
   )
{
   SCIP_DUALPACKET m;
   int rest;
   int nfull;
   int i;

   assert(inp != NULL || count == 0);
   assert(out != NULL || count == 0);
   assert(count >= 0);
   assert(SCIP_DUALPACKETSIZE == 16);

   rest = count % (int)SCIP_DUALPACKETSIZE;
   nfull = count - rest;

   for( i = 0; i < nfull; i += (int)SCIP_DUALPACKETSIZE )
   {
      assert(inp != NULL);
      assert(out != NULL);

      m = *inp++;

      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      m >>= 2;
      *out++ = -(int)(m & 3);
      assert(m >> 2 == 0);
   }

   if( rest > 0 )
   {
      assert(inp != NULL);
      assert(out != NULL);

      m = *inp;
      for( i = 0; i < rest; i++ )
      {
         *out++ = -(int)(m & 3);
         m >>= 2;
      }
   }
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

   SCIPencodeDualBitNeg(cstat, lpistate->packcstat, lpistate->ncols);
   SCIPencodeDualBitNeg(rstat, lpistate->packrstat, lpistate->nrows);
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

   SCIPdecodeDualBitNeg(lpistate->packcstat, cstat, lpistate->ncols);
   SCIPdecodeDualBitNeg(lpistate->packrstat, rstat, lpistate->nrows);
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

   BMSfreeBlockMemoryArrayNull(blkmem, &(*lpistate)->packcstat, colpacketNum((*lpistate)->ncols));
   BMSfreeBlockMemoryArrayNull(blkmem, &(*lpistate)->packrstat, rowpacketNum((*lpistate)->nrows));
   BMSfreeBlockMemory(blkmem, lpistate);
}



/*
 * local methods
 */

/** gets all Gurobi parameters used in LPI */
static
SCIP_RETCODE getParameterValues(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   GRBPARAM*             grbparam            /**< Gurobi parameters */
   )
{
   int i;

   assert( lpi != NULL );
   assert( lpi->grbenv != NULL );
   assert( grbparam != NULL );

   SCIPdebugMessage("getParameterValues()\n");

   for( i = 0; i < NUMINTPARAM; ++i )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetintparam(lpi->grbenv, intparam[i], &(grbparam->intparval[i])) );
   }
   for( i = 0; i < NUMDBLPARAM; ++i )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, dblparam[i], &(grbparam->dblparval[i])) );
   }

   return SCIP_OKAY;
}

/** in debug mode, checks validity of Gurobi parameters */
static
SCIP_RETCODE checkParameterValues(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
#ifndef NDEBUG
   GRBPARAM par;
   int i;

   SCIP_CALL( getParameterValues(lpi, &par) );
   for (i = 0; i < NUMINTPARAM; ++i)
      assert( lpi->curparam.intparval[i] == par.intparval[i] );
   for (i = 0; i < NUMDBLPARAM; ++i)
      assert(MAX(lpi->curparam.dblparval[i], dblparammin[i]) == par.dblparval[i]); /*lint !e777*/
#endif

   return SCIP_OKAY;
}

/** sets all Gurobi parameters used in LPI */
static
SCIP_RETCODE setParameterValues(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   GRBPARAM*             grbparam            /**< Gurobi parameters */
   )
{
   int i;

   assert( lpi != NULL );
   assert( lpi->grbenv != NULL );
   assert( grbparam != NULL );

   SCIPdebugMessage("setParameterValues()\n");

   for( i = 0; i < NUMINTPARAM; ++i )
   {
      if( lpi->curparam.intparval[i] != grbparam->intparval[i] )
      {
         SCIPdebugMessage("setting Gurobi int parameter %s from %d to %d\n",
            intparam[i], lpi->curparam.intparval[i], grbparam->intparval[i]);
         lpi->curparam.intparval[i] = grbparam->intparval[i];
         CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, intparam[i], lpi->curparam.intparval[i]) );
      }
   }
   for( i = 0; i < NUMDBLPARAM; ++i )
   {
      if( lpi->curparam.dblparval[i] != grbparam->dblparval[i] ) /*lint !e777*/
      {
         SCIPdebugMessage("setting Gurobi dbl parameter %s from %g to %g\n",
            dblparam[i], lpi->curparam.dblparval[i], MAX(grbparam->dblparval[i], dblparammin[i]));
         lpi->curparam.dblparval[i] = MAX(grbparam->dblparval[i], dblparammin[i]);
         CHECK_ZERO( lpi->messagehdlr, GRBsetdblparam(lpi->grbenv, dblparam[i], lpi->curparam.dblparval[i]) );
      }
   }

   SCIP_CALL( checkParameterValues(lpi) );

   return SCIP_OKAY;
}

/** copies Gurobi parameters from source to dest */
static
void copyParameterValues(
   GRBPARAM*             dest,               /**< destination Gurobi parameters */
   const GRBPARAM*       source              /**< original Gurobi parameters */
   )
{
   int i;

   for( i = 0; i < NUMINTPARAM; ++i )
      dest->intparval[i] = source->intparval[i];
   for( i = 0; i < NUMDBLPARAM; ++i )
      dest->dblparval[i] = source->dblparval[i];
}

/** gets a single integer parameter value */
static
SCIP_RETCODE getIntParam(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           param,              /**< parameter name */
   int*                  p                   /**< value of parameter */
   )
{
   int i;

   assert( lpi != NULL );

   for( i = 0; i < NUMINTPARAM; ++i )
   {
      if( strcmp(intparam[i], param) == 0 )
      {
         *p = lpi->grbparam.intparval[i];
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("unknown Gurobi integer parameter <%s>.\n", param);
   return SCIP_LPERROR;
}

/** sets a single integer parameter value */
static
SCIP_RETCODE setIntParam(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           param,              /**< parameter name */
   int                   parval              /**< value of parameter */
   )
{
   int i;

   assert( lpi != NULL );

   for( i = 0; i < NUMINTPARAM; ++i )
   {
      if( strcmp(intparam[i], param) == 0 )
      {
         lpi->grbparam.intparval[i] = parval;
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("unknown Gurobi integer parameter <%s>.\n", param);
   return SCIP_LPERROR;
}

/** gets a single double parameter value */
static
SCIP_RETCODE getDblParam(SCIP_LPI* lpi, const char* param, double* p)
{
   int i;

   assert(lpi != NULL);

   for( i = 0; i < NUMDBLPARAM; ++i )
   {
      if( strcmp(dblparam[i], param) == 0 )
      {
         *p = lpi->grbparam.dblparval[i];
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("unknown Gurobi double parameter <%s>.\n", param);
   return SCIP_LPERROR;
}

/** sets a single double parameter value */
static
SCIP_RETCODE setDblParam(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           param,              /**< parameter name */
   double                parval              /**< value of parameter */
   )
{
   int i;

   assert( lpi != NULL );

   for( i = 0; i < NUMDBLPARAM; ++i )
   {
      if( strcmp(dblparam[i], param) == 0 )
      {
         lpi->grbparam.dblparval[i] = parval;
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("unknown Gurobi double parameter <%s>.\n", param);
   return SCIP_LPERROR;
}

/** marks the current LP to be unsolved */
static
void invalidateSolution(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   lpi->solstat = -1;
}

/** converts SCIP's lhs/rhs pairs into Gurobi's sen/rhs */
static
void convertSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand side vector */
   const SCIP_Real*      rhs,                /**< right hand side vector */
   int*                  rngcount            /**< number of ranged rows found */
   )
{
   int i;

   assert(lpi != NULL);
   assert(nrows >= 0);
   assert(lhs != NULL);
   assert(rhs != NULL);

   /* convert lhs/rhs into sen/rhs */
   *rngcount = 0;
   for( i = 0; i < nrows; ++i )
   {
      assert(lhs[i] <= rhs[i]);

      if( lhs[i] == rhs[i] ) /*lint !e777*/
      {
         assert(-GRB_INFINITY < rhs[i] && rhs[i] < GRB_INFINITY);
         lpi->senarray[i] = GRB_EQUAL;
         lpi->rhsarray[i] = rhs[i];
      }
      else if( lhs[i] <= -GRB_INFINITY )
      {
         assert(-GRB_INFINITY < rhs[i] && rhs[i] < GRB_INFINITY);
         lpi->senarray[i] = GRB_LESS_EQUAL;
         lpi->rhsarray[i] = rhs[i];
      }
      else if( rhs[i] >= GRB_INFINITY )
      {
         assert(-GRB_INFINITY < lhs[i] && lhs[i] < GRB_INFINITY);
         lpi->senarray[i] = GRB_GREATER_EQUAL;
         lpi->rhsarray[i] = lhs[i];
      }
      else
      {
         /* Gurobi cannot handle ranged rows */
         SCIPerrorMessage("Gurobi cannot handle ranged rows.\n");
         SCIPABORT();
         (*rngcount)++;
      }
   }
}

/** converts Gurobi's sen/rhs pairs into SCIP's lhs/rhs pairs */
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

   for(  i = 0; i < nrows; ++i )
   {
      switch( lpi->senarray[i] )
      {
      case GRB_EQUAL:
         lhs[i] = lpi->rhsarray[i];
         rhs[i] = lpi->rhsarray[i];
         break;

      case GRB_LESS_EQUAL:
         lhs[i] = -GRB_INFINITY;
         rhs[i] = lpi->rhsarray[i];
         break;

      case GRB_GREATER_EQUAL:
         lhs[i] = lpi->rhsarray[i];
         rhs[i] = GRB_INFINITY;
         break;

      default:
         SCIPerrorMessage("invalid row sense\n");
         SCIPABORT();
      }
      assert(lhs[i] <= rhs[i]);
   }
}

/** converts Gurobi's sen/rhs pairs into SCIP's lhs/rhs pairs, only storing the left hand side */
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

   for(  i = 0; i < nrows; ++i )
   {
      switch( lpi->senarray[i] )
      {
      case GRB_EQUAL:
         lhs[i] = lpi->rhsarray[i];
         break;

      case GRB_LESS_EQUAL:
         lhs[i] = -GRB_INFINITY;
         break;

      case GRB_GREATER_EQUAL:
         lhs[i] = lpi->rhsarray[i];
         break;

      default:
         SCIPerrorMessage("invalid row sense\n");
         SCIPABORT();
      }
   }
}

/** converts Gurobi's sen/rhs pairs into SCIP's lhs/rhs pairs, only storing the right hand side */
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

   for(  i = 0; i < nrows; ++i )
   {
      switch( lpi->senarray[i] )
      {
      case GRB_EQUAL:
         rhs[i] = lpi->rhsarray[i];
         break;

      case GRB_LESS_EQUAL:
         rhs[i] = lpi->rhsarray[i];
         break;

      case GRB_GREATER_EQUAL:
         rhs[i] = GRB_INFINITY;
         break;

      default:
         SCIPerrorMessage("invalid row sense\n");
         SCIPABORT();
      }
   }
}

/** converts Gurobi's sen/rhs pairs into SCIP's lhs/rhs pairs */
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

static char grbname[100];

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   int major;
   int minor;
   int technical;

   GRBversion(&major, &minor, &technical);
   sprintf(grbname, "Gurobi %d.%d.%d", major, minor, technical);
   return grbname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return "Linear Programming Solver developed by Gurobi Optimization (www.gurobi.com)";
}

/** gets pointer for LP solver - use only with great care
 *
 *  Here we return the pointer to the model.
 */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   return (void*) lpi->grbmodel;
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
   assert(sizeof(SCIP_Real) == sizeof(double)); /* Gurobi only works with doubles as floating points */
   assert(sizeof(SCIP_Bool) == sizeof(int));    /* Gurobi only works with ints as bools */
   assert(lpi != NULL);
   assert(numlp >= 0);

   SCIPdebugMessage("SCIPlpiCreate()\n");

   /* create environment
    *
    * Each problem will get a copy of the original environment. Thus, grbenv is only needed once.
    */
   if ( grbenv == NULL )
   {
      /* initialize environment - no log file */
      CHECK_ZERO( messagehdlr, GRBloadenv(&grbenv, NULL) );

      /* turn off output for all models */
      CHECK_ZERO( messagehdlr, GRBsetintparam(grbenv, GRB_INT_PAR_OUTPUTFLAG, 0) );

      /* turn on that basis information for infeasible and unbounded models is available */
      CHECK_ZERO( messagehdlr, GRBsetintparam(grbenv, GRB_INT_PAR_INFUNBDINFO, 1) );
   }
   assert( grbenv != NULL );

   /* create empty LPI */
   SCIP_ALLOC( BMSallocMemory(lpi) );
   CHECK_ZERO( messagehdlr, GRBnewmodel(grbenv, &(*lpi)->grbmodel, name, 0, NULL, NULL, NULL, NULL, NULL) );

   /* get local copy of environment */
   (*lpi)->grbenv = GRBgetenv((*lpi)->grbmodel);
   (*lpi)->senarray = NULL;
   (*lpi)->rhsarray = NULL;
   (*lpi)->valarray = NULL;
   (*lpi)->cstat = NULL;
   (*lpi)->rstat = NULL;
   (*lpi)->indarray = NULL;
   (*lpi)->sidechgsize = 0;
   (*lpi)->valsize = 0;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->iterations = 0;
   (*lpi)->solisbasic = FALSE;
   (*lpi)->pricing = SCIP_PRICING_LPIDEFAULT;
   (*lpi)->messagehdlr = messagehdlr;
   invalidateSolution(*lpi);

   /* get default parameter values */
   SCIP_CALL( getParameterValues((*lpi), &((*lpi)->defparam)) );
   copyParameterValues(&((*lpi)->curparam), &((*lpi)->defparam));
   copyParameterValues(&((*lpi)->grbparam), &((*lpi)->defparam));
   ++numlp;

   /* set objective sense */
   SCIP_CALL( SCIPlpiChgObjsen(*lpi, objsen) );

   /* set default pricing */
   SCIP_CALL( SCIPlpiSetIntpar(*lpi, SCIP_LPPAR_PRICING, (*lpi)->pricing) );

   SCIPmessagePrintWarning(messagehdlr, "The Gurobi LPI is a beta version only - use with care.\n\n");

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(grbenv != NULL);
   assert(lpi != NULL);
   assert(*lpi != NULL);

   SCIPdebugMessage("SCIPlpiFree()\n");

   /* free model */
   CHECK_ZERO( (*lpi)->messagehdlr, GRBfreemodel((*lpi)->grbmodel) );

   /* free memory */
   BMSfreeMemoryArrayNull(&(*lpi)->senarray);
   BMSfreeMemoryArrayNull(&(*lpi)->rhsarray);
   BMSfreeMemoryArrayNull(&(*lpi)->cstat);
   BMSfreeMemoryArrayNull(&(*lpi)->rstat);
   BMSfreeMemory(lpi);

   /* free environment */
   --numlp;
   if( numlp == 0 )
   {
      GRBfreeenv(grbenv);
      grbenv = NULL;
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
   int* cnt;
   int rngcount;
   int c;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->grbenv != NULL);
   assert(objsen == SCIP_OBJSEN_MAXIMIZE || objsen == SCIP_OBJSEN_MINIMIZE);

   SCIPdebugMessage("loading LP in column format into Gurobi: %d cols, %d rows\n", ncols, nrows);

   invalidateSolution(lpi);

   SCIP_CALL( ensureSidechgMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   convertSides(lpi, nrows, lhs, rhs, &rngcount);
   assert( rngcount == 0 );

   /* calculate column lengths */
   SCIP_ALLOC( BMSallocMemoryArray(&cnt, ncols) );
   for( c = 0; c < ncols-1; ++c )
   {
      cnt[c] = beg[c+1] - beg[c];
      assert(cnt[c] >= 0);
   }
   cnt[ncols-1] = nnonz - beg[ncols-1];
   assert(cnt[ncols-1] >= 0);

   /* delete model */
   assert( lpi->grbmodel != NULL );
   CHECK_ZERO( lpi->messagehdlr, GRBfreemodel(lpi->grbmodel) );

   /* load model - all variables are continuous */
   CHECK_ZERO( lpi->messagehdlr, GRBloadmodel(lpi->grbenv, &(lpi->grbmodel), NULL, ncols, nrows, objsen, 0.0, (SCIP_Real*)obj,
         lpi->senarray, lpi->rhsarray, (int*)beg, cnt, (int*)ind, (SCIP_Real*)val, (SCIP_Real*)lb, (SCIP_Real*)ub, NULL, colnames, rownames) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   /* free temporary memory */
   BMSfreeMemoryArray(&cnt);

#ifndef NDEBUG
   {
      int temp;

      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &temp) );
      assert( temp == ncols);

      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &temp) );
      assert( temp == nrows);

      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMNZS, &temp) );
      assert( temp == nnonz);
   }
#endif

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
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("adding %d columns with %d nonzeros to Gurobi\n", ncols, nnonz);

   invalidateSolution(lpi);

   /* add columns - all new variables are continuous */
   CHECK_ZERO( lpi->messagehdlr, GRBaddvars(lpi->grbmodel, ncols, nnonz, (int*)beg, (int*)ind, (SCIP_Real*)val, (SCIP_Real*)obj, (SCIP_Real*)lb, (SCIP_Real*)ub, NULL, colnames) )
      CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   int j;
   int* which;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
#ifndef NDEBUG
   {
      int temp;

      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &temp) );
      assert(0 <= firstcol && firstcol <= lastcol && lastcol < temp);
   }
#endif

   SCIPdebugMessage("deleting %d columns from Gurobi\n", lastcol - firstcol + 1);

   invalidateSolution(lpi);

   /* Gurobi can't delete a range of columns, we have to set up an index array */
   SCIP_ALLOC( BMSallocMemoryArray(&which, lastcol-firstcol+1) );;
   for( j = firstcol; j <= lastcol; ++j )
      which[j - firstcol] = j;

   CHECK_ZERO( lpi->messagehdlr, GRBdelvars(lpi->grbmodel, lastcol-firstcol+1, which) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   BMSfreeMemoryArray( &which );

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
   int j, nvars, num;
   int* which;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("deleting a column set from Gurobi\n");

   invalidateSolution(lpi);

   /* Gurobi can't delete a range of columns, we have to set up an index array */
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &nvars) );

   SCIP_ALLOC( BMSallocMemoryArray(&which, nvars) );;
   num = 0;
   for( j = 0; j < nvars; ++j )
   {
      if( dstat[j] )
         which[num++] = j;
   }
   CHECK_ZERO( lpi->messagehdlr, GRBdelvars(lpi->grbmodel, num, which) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   BMSfreeMemoryArray( &which );

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
   int rngcount;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("adding %d rows with %d nonzeros to Gurobi\n", nrows, nnonz);

   invalidateSolution(lpi);

   SCIP_CALL( ensureSidechgMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   convertSides(lpi, nrows, lhs, rhs, &rngcount);
   assert( rngcount == 0 );

   /* add rows to LP */
   CHECK_ZERO( lpi->messagehdlr, GRBaddconstrs(lpi->grbmodel, nrows, nnonz, (int*)beg, (int*)ind, (SCIP_Real*)val, lpi->senarray, lpi->rhsarray, rownames) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   return SCIP_OKAY;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   int i;
   int* which;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
#ifndef NDEBUG
   {
      int nrows;
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );
      assert(0 <= firstrow && firstrow <= lastrow && lastrow < nrows);
   }
#endif

   SCIPdebugMessage("deleting %d rows from Gurobi\n", lastrow - firstrow + 1);

   invalidateSolution(lpi);

   /* Gurobi can't delete a range of rows, we have to set up an index array */
   SCIP_ALLOC( BMSallocMemoryArray(&which, lastrow-firstrow+1) );;
   for( i = firstrow; i <= lastrow; ++i )
      which[i - firstrow] = i;

   CHECK_ZERO( lpi->messagehdlr, GRBdelconstrs(lpi->grbmodel, lastrow-firstrow+1, which) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

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
   int i, num;
   int nrows;
   int* which;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("deleting a row set from Gurobi\n");

   invalidateSolution(lpi);

   /* Gurobi can't delete a range of rows, we have to set up an index array */
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&which, nrows) );;
   num = 0;
   for( i = 0; i < nrows; ++i )
   {
      if( dstat[i] )
         which[num++] = i;
   }
   CHECK_ZERO( lpi->messagehdlr, GRBdelconstrs(lpi->grbmodel, num, which) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   /* update dstat */
   num = 0;
   for( i = 0; i < nrows; ++i )
   {
      if( dstat[i] )
      {
         dstat[i] = -1;
         ++num;
      }
      else
         dstat[i] = i - num;
   }

   BMSfreeMemoryArray( &which );

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );

   SCIPdebugMessage("clearing Gurobi LP\n");

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, GRBfreemodel(lpi->grbmodel) );
   CHECK_ZERO( lpi->messagehdlr, GRBnewmodel(lpi->grbenv, &(lpi->grbmodel), "", 0, NULL, NULL, NULL, NULL, NULL) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

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
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("changing %d bounds in Gurobi\n", ncols);
#ifdef SCIP_DEBUG
   {
      int i;
      for( i = 0; i < ncols; ++i )
         SCIPdebugPrintf("  col %d: [%g,%g]\n", ind[i], lb[i], ub[i]);
   }
#endif

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrlist(lpi->grbmodel, GRB_DBL_ATTR_LB, ncols, (int*)ind, (SCIP_Real*)lb) );
   CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrlist(lpi->grbmodel, GRB_DBL_ATTR_UB, ncols, (int*)ind, (SCIP_Real*)ub) );

   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

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
   int rngcount;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("changing %d sides in Gurobi\n", nrows);

   invalidateSolution(lpi);

   /* convert lhs/rhs into sen/rhs/range tuples */
   SCIP_CALL( ensureSidechgMem(lpi, nrows) );
   convertSides(lpi, nrows, lhs, rhs, &rngcount);
   assert( rngcount == 0 );

   /* change row sides */
   CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrlist(lpi->grbmodel, GRB_DBL_ATTR_RHS, nrows, (int*)ind, lpi->rhsarray) );
   CHECK_ZERO( lpi->messagehdlr, GRBsetcharattrlist(lpi->grbmodel, GRB_CHAR_ATTR_SENSE, nrows, (int*)ind, lpi->senarray) );

   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

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
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("changing coefficient row %d, column %d in Gurobi to %g\n", row, col, newval);

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, GRBchgcoeffs(lpi->grbmodel, 1, &row, &col, &newval) );
   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(objsen == SCIP_OBJSEN_MAXIMIZE || objsen == SCIP_OBJSEN_MINIMIZE);

   SCIPdebugMessage("changing objective sense in Gurobi to %d\n", objsen);

   invalidateSolution(lpi);

   /* The objective sense of Gurobi and SCIP are equal */
   CHECK_ZERO( lpi->messagehdlr, GRBsetintattr(lpi->grbmodel, GRB_INT_ATTR_MODELSENSE, objsen) );

   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

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
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("changing %d objective values in Gurobi\n", ncols);

   CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrlist(lpi->grbmodel, GRB_DBL_ATTR_OBJ, ncols, ind, obj) );

   CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );

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
   int ncols;
   int beg;
   int i;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(scaleval != 0.0);

   SCIPdebugMessage("scaling row %d with factor %g in Gurobi\n", row, scaleval);

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
   SCIP_CALL( ensureValMem(lpi, ncols) );

   /* get the row */
   SCIP_CALL( SCIPlpiGetRows(lpi, row, row, &lhs, &rhs, &nnonz, &beg, lpi->indarray, lpi->valarray) );

   /* scale row coefficients */
   for(  i = 0; i < nnonz; ++i )
   {
      SCIP_CALL( SCIPlpiChgCoef(lpi, row, lpi->indarray[i], lpi->valarray[i] * scaleval) );
   }

   /* scale row sides */
   if( lhs > -GRB_INFINITY )
      lhs *= scaleval;
   else if( scaleval < 0.0 )
      lhs = GRB_INFINITY;
   if( rhs < GRB_INFINITY )
      rhs *= scaleval;
   else if( scaleval < 0.0 )
      rhs = -GRB_INFINITY;
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
   int ncols;
   int beg;
   int i;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(scaleval != 0.0);

   SCIPdebugMessage("scaling column %d with factor %g in Gurobi\n", col, scaleval);

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
   SCIP_CALL( ensureValMem(lpi, ncols) );

   /* get the column */
   SCIP_CALL( SCIPlpiGetCols(lpi, col, col, &lb, &ub, &nnonz, &beg, lpi->indarray, lpi->valarray) );

   /* get objective coefficient */
   SCIP_CALL( SCIPlpiGetObj(lpi, col, col, &obj) );

   /* scale column coefficients */
   for(  i = 0; i < nnonz; ++i )
   {
      SCIP_CALL( SCIPlpiChgCoef(lpi, lpi->indarray[i], col, lpi->valarray[i] * scaleval) );
   }

   /* scale objective value */
   obj *= scaleval;
   SCIP_CALL( SCIPlpiChgObj(lpi, 1, &col, &obj) );

   /* scale column bounds */
   if( lb > -GRB_INFINITY )
      lb /= scaleval;
   else if( scaleval < 0.0 )
      lb = GRB_INFINITY;
   if( ub < GRB_INFINITY )
      ub /= scaleval;
   else if( scaleval < 0.0 )
      ub = -GRB_INFINITY;
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
   assert(nrows != NULL);

   SCIPdebugMessage("getting number of rows\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, nrows) );

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   assert(lpi != NULL);
   assert(ncols != NULL);

   SCIPdebugMessage("getting number of columns\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, ncols) );

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   assert(lpi != NULL);
   assert(nnonz != NULL);

   SCIPdebugMessage("getting number of non-zeros\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMNZS, nnonz) );

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
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
#ifndef NDEBUG
   {
      int ncols;
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
      assert(0 <= firstcol && firstcol <= lastcol && lastcol < ncols);
   }
#endif

   SCIPdebugMessage("getting columns %d to %d\n", firstcol, lastcol);

   if( lb != NULL )
   {
      assert(ub != NULL);

      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_LB, firstcol, lastcol-firstcol+1, lb) );
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_UB, firstcol, lastcol-firstcol+1, ub) );
   }
   else
      assert(ub == NULL);

   if( nnonz != NULL )
   {
      assert(beg != NULL);
      assert(ind != NULL);
      assert(val != NULL);

      /* get matrix entries */
      CHECK_ZERO( lpi->messagehdlr, GRBgetvars(lpi->grbmodel, nnonz, beg, ind, val, firstcol, lastcol-firstcol+1) )
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
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
#ifndef NDEBUG
   {
      int nrows;
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );
      assert(0 <= firstrow && firstrow <= lastrow && lastrow < nrows);
   }
#endif

   SCIPdebugMessage("getting rows %d to %d\n", firstrow, lastrow);

   if( lhs != NULL || rhs != NULL )
   {
      /* get row sense and rhs */
      SCIP_CALL( ensureSidechgMem(lpi, lastrow - firstrow + 1) );
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_RHS, firstrow, lastrow-firstrow+1, lpi->rhsarray) );
      CHECK_ZERO( lpi->messagehdlr, GRBgetcharattrarray(lpi->grbmodel, GRB_CHAR_ATTR_SENSE, firstrow, lastrow-firstrow+1, lpi->senarray) );

      /* convert sen and rhs into lhs/rhs tuples */
      reconvertSides(lpi, lastrow - firstrow + 1, lhs, rhs);
   }

   if( nnonz != NULL )
   {
      assert(beg != NULL);
      assert(ind != NULL);
      assert(val != NULL);

      /* get matrix entries */
      CHECK_ZERO( lpi->messagehdlr, GRBgetconstrs(lpi->grbmodel, nnonz, beg, ind, val, firstrow, lastrow-firstrow+1) );
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
   int grbobjsen;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( objsen != NULL );

   /* note that the objective sense is define equally in SCIP (LPI) and Gurobi */
   assert( GRB_MINIMIZE == SCIP_OBJSEN_MINIMIZE );
   assert( GRB_MAXIMIZE == SCIP_OBJSEN_MAXIMIZE );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_MODELSENSE, &grbobjsen) );

   *objsen = grbobjsen;

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
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(firstcol <= lastcol);
   assert(vals != NULL);

   SCIPdebugMessage("getting objective values %d to %d\n", firstcol, lastcol);

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_OBJ, firstcol, lastcol-firstcol+1, vals) );

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
   assert(lpi->grbmodel != NULL);
#ifndef NDEBUG
   {
      int ncols;
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
      assert(0 <= firstcol && firstcol <= lastcol && lastcol < ncols);
   }
#endif

   SCIPdebugMessage("getting bounds %d to %d\n", firstcol, lastcol);

   if( lbs != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_LB, firstcol, lastcol-firstcol+1, lbs) );
   }

   if( ubs != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_UB, firstcol, lastcol-firstcol+1, ubs) );
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
   assert(lpi->grbmodel != NULL);
   assert(firstrow <= lastrow);

   SCIPdebugMessage("getting row sides %d to %d\n", firstrow, lastrow);

   /* get row sense, rhs, and ranges */
   SCIP_CALL( ensureSidechgMem(lpi, lastrow - firstrow + 1) );

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_RHS, firstrow, lastrow-firstrow+1, lpi->rhsarray) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetcharattrarray(lpi->grbmodel, GRB_CHAR_ATTR_SENSE, firstrow, lastrow-firstrow+1, lpi->senarray) );

   /* convert sen and rhs into lhs/rhs tuples */
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
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("getting coefficient of row %d col %d\n", row, col);

   CHECK_ZERO( lpi->messagehdlr, GRBgetcoeff(lpi->grbmodel, row, col, val) );

   return SCIP_OKAY;
}

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** calls primal simplex to solve the LP
 *
 *  @todo Check concurrent (GRB_METHOD_CONCURRENT or GRB_METHOD_DETERMINISTIC_CONCURRENT)
 */
SCIP_RETCODE SCIPlpiSolvePrimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int retval;
   int primalfeasible;
   int dualfeasible;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );

#ifdef SCIP_DEBUG
   {
      int ncols, nrows;
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );
      SCIPdebugMessage("calling Gurobi primal simplex: %d cols, %d rows\n", ncols, nrows);
   }
#endif

   invalidateSolution(lpi);

   SCIPdebugMessage("calling GRBoptimize() - primal\n");

   /* set primal simplex */
   SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );
   CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, GRB_METHOD_PRIMAL) );

   retval = GRBoptimize(lpi->grbmodel);
   switch( retval  )
   {
   case 0:
      break;
   case GRB_ERROR_OUT_OF_MEMORY:
      return SCIP_NOMEMORY;
   default:
      return SCIP_LPERROR;
   }

   lpi->solisbasic = TRUE;
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &lpi->solstat) );

   /*
     CHECK_ZERO( lpi->messagehdlr, CPXsolninfo(lpi->grbenv, lpi->grbmodel, NULL, NULL, &primalfeasible, &dualfeasible) );
     SCIPdebugMessage(" -> Gurobi returned solstat=%d, pfeas=%d, dfeas=%d (%d iterations)\n",
     lpi->solstat, primalfeasible, dualfeasible, lpi->iterations);
   */
   primalfeasible = FALSE;
   dualfeasible = FALSE;

   if( lpi->solstat == GRB_INF_OR_UNBD
      || (lpi->solstat == GRB_INFEASIBLE && !dualfeasible)
      || (lpi->solstat == GRB_UNBOUNDED && !primalfeasible) )
   {
      int cnt;
      int presolve;
      CHECK_ZERO( lpi->messagehdlr, GRBgetintparam(lpi->grbenv, GRB_INT_PAR_PRESOLVE, &presolve) );
      if( presolve != GRB_PRESOLVE_OFF )
      {
         /* maybe the preprocessor solved the problem; but we need a solution, so solve again without preprocessing */
         SCIPdebugMessage("presolver may have solved the problem -> calling Gurobi primal simplex again without presolve\n");

         /* switch off preprocessing */
         CHECK_ZERO( lpi->messagehdlr, GRBgetintparam(lpi->grbenv, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF) );
         SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );

         retval = GRBoptimize(lpi->grbmodel);
         switch( retval  )
         {
         case 0:
            break;
         case GRB_ERROR_OUT_OF_MEMORY:
            return SCIP_NOMEMORY;
         default:
            return SCIP_LPERROR;
         }

         CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
         lpi->iterations += cnt;
         CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &lpi->solstat) );
         SCIPdebugMessage(" -> Gurobi returned solstat=%d (%d iterations)\n", lpi->solstat, lpi->iterations);

         /* switch on preprocessing again */
         CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_AUTO) );
      }

      if( lpi->solstat == GRB_INF_OR_UNBD )
      {
         /* preprocessing was not the problem; issue a warning message and treat LP as infeasible */
         SCIPerrorMessage("Gurobi primal simplex returned GRB_INF_OR_UNBD after presolving was turned off\n");
      }
   }

   return SCIP_OKAY;
}

/** calls dual simplex to solve the LP
 *
 *  @todo Check concurrent (GRB_METHOD_CONCURRENT or GRB_METHOD_DETERMINISTIC_CONCURRENT)
 */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int retval;
   double cnt;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );

#ifdef SCIP_DEBUG
   {
      int ncols, nrows;
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );
      SCIPdebugMessage("calling Gurobi dual simplex: %d cols, %d rows\n", ncols, nrows);
   }
#endif

   invalidateSolution(lpi);

   SCIPdebugMessage("calling GRBoptimize() - dual\n");

   SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );

   /* set dual simplex */
   CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL) );

   retval = GRBoptimize(lpi->grbmodel);
   switch( retval  )
   {
   case 0:
      break;
   case GRB_ERROR_OUT_OF_MEMORY:
      return SCIP_NOMEMORY;
   default:
      return SCIP_LPERROR;
   }

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
   lpi->iterations = (int) cnt;

   lpi->solisbasic = TRUE;
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &lpi->solstat) );

   /*
     SCIPdebugMessage(" -> Gurobi returned solstat=%d, pfeas=%d, dfeas=%d (%d iterations)\n",
     lpi->solstat, primalfeasible, dualfeasible, lpi->iterations);
   */

   if( lpi->solstat == GRB_INF_OR_UNBD )
   {
      int presolve;
      CHECK_ZERO( lpi->messagehdlr, getIntParam(lpi, GRB_INT_PAR_PRESOLVE, &presolve) );

      if( presolve != GRB_PRESOLVE_OFF )
      {
         /* maybe the preprocessor solved the problem; but we need a solution, so solve again without preprocessing */
         SCIPdebugMessage("presolver may have solved the problem -> calling Gurobi dual simplex again without presolve\n");

         /* switch off preprocessing */
         CHECK_ZERO( lpi->messagehdlr, setIntParam(lpi, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF) );
         SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );

         retval = GRBoptimize(lpi->grbmodel);
         switch( retval  )
         {
         case 0:
            break;
         case GRB_ERROR_OUT_OF_MEMORY:
            return SCIP_NOMEMORY;
         default:
            return SCIP_LPERROR;
         }

         CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
         lpi->iterations += (int) cnt;
         CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &lpi->solstat) );
         SCIPdebugMessage(" -> Gurobi returned solstat=%d (%d iterations)\n", lpi->solstat, lpi->iterations);

         /* switch on preprocessing again */
         CHECK_ZERO( lpi->messagehdlr, setIntParam(lpi, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_AUTO) );
      }

      if( lpi->solstat == GRB_INF_OR_UNBD )
      {
         /* preprocessing was not the problem; issue a warning message and treat LP as infeasible */
         SCIPerrorMessage("Gurobi dual simplex returned GRB_INF_OR_UNBD after presolving was turned off\n");
      }
   }

   return SCIP_OKAY;
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{
   int retval;
   double cnt;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );

#ifdef SCIP_DEBUG
   {
      int ncols, nrows;
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );
      SCIPdebugMessage("calling Gurobi barrier: %d cols, %d rows\n", ncols, nrows);
   }
#endif

   invalidateSolution(lpi);

   SCIPdebugMessage("calling GRBoptimize() - barrier\n");

   /* set barrier */
   SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );

   if( crossover )
   {
      /* turn on crossover to automatic setting (-1) */
      CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_CROSSOVER, -1) );
   }
   else
   {
      /* turn off crossover */
      CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_CROSSOVER, 0) );
   }

   CHECK_ZERO( lpi->messagehdlr, GRBsetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER) );

   retval = GRBoptimize(lpi->grbmodel);
   switch( retval  )
   {
   case 0:
      break;
   case GRB_ERROR_OUT_OF_MEMORY:
      return SCIP_NOMEMORY;
   default:
      return SCIP_LPERROR;
   }

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
   lpi->iterations = (int) cnt;

   lpi->solisbasic = crossover;
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &lpi->solstat) );

   /*
     SCIPdebugMessage(" -> Gurobi returned solstat=%d, pfeas=%d, dfeas=%d (%d iterations)\n",
     lpi->solstat, primalfeasible, dualfeasible, lpi->iterations);
   */

   if( lpi->solstat == GRB_INF_OR_UNBD )
   {
      int presolve;
      CHECK_ZERO( lpi->messagehdlr, getIntParam(lpi, GRB_INT_PAR_PRESOLVE, &presolve) );

      if( presolve != GRB_PRESOLVE_OFF )
      {
         /* maybe the preprocessor solved the problem; but we need a solution, so solve again without preprocessing */
         SCIPdebugMessage("presolver may have solved the problem -> calling Gurobi barrier again without presolve\n");

         /* switch off preprocessing */
         CHECK_ZERO( lpi->messagehdlr, setIntParam(lpi, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF) );
         SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );

         retval = GRBoptimize(lpi->grbmodel);
         switch( retval  )
         {
         case 0:
            break;
         case GRB_ERROR_OUT_OF_MEMORY:
            return SCIP_NOMEMORY;
         default:
            return SCIP_LPERROR;
         }

         CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_ITERCOUNT, &cnt) );
         lpi->iterations += (int) cnt;
         CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_STATUS, &lpi->solstat) );
         SCIPdebugMessage(" -> Gurobi returned solstat=%d (%d iterations)\n", lpi->solstat, lpi->iterations);

         /* switch on preprocessing again */
         CHECK_ZERO( lpi->messagehdlr, setIntParam(lpi, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_AUTO) );
      }

      if( lpi->solstat == GRB_INF_OR_UNBD )
      {
         /* preprocessing was not the problem; issue a warning message and treat LP as infeasible */
         SCIPerrorMessage("Gurobi dual simplex returned GRB_INF_OR_UNBD after presolving was turned off\n");
      }
   }
   return SCIP_OKAY;
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
   SCIP_Real oldlb;
   SCIP_Real oldub;
   SCIP_Real newlb;
   SCIP_Real newub;
   SCIP_Real olditlim;
   SCIP_Bool error;
   SCIP_Bool success;
   int objsen;
   int it;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );

   SCIPdebugMessage("performing strong branching on variable %d (%d iterations)\n", col, itlim);

   SCIP_CALL( setParameterValues(lpi, &(lpi->grbparam)) );

   error = FALSE;
   *downvalid = FALSE;
   *upvalid = FALSE;
   if( iter != NULL )
      *iter = 0;

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_MODELSENSE, &objsen) );

   /* save current LP basis and bounds*/
   SCIP_CALL( getBase(lpi, &success) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_LB, col, &oldlb) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_UB, col, &oldub) );

   /* save old iteration limit and set iteration limit to strong branching limit */
   if( itlim > INT_MAX )
      itlim = INT_MAX;

   SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_ITERATIONLIMIT, &olditlim) );
   SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_ITERATIONLIMIT, (double) itlim) );

   /* down branch */
   newub = EPSCEIL(psol-1.0, 1e-06);
   if( newub >= oldlb - 0.5 )
   {
      SCIPdebugMessage("strong branching down (%g) on x%d (%g) with %d iterations\n", newub, col, psol, itlim);

      CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_UB, col, newub) );

      SCIP_CALL( SCIPlpiSolveDual(lpi) );
      if( SCIPlpiIsOptimal(lpi) || SCIPlpiIsIterlimExc(lpi) )
      {
         SCIP_CALL( SCIPlpiGetObjval(lpi, down) );
         *downvalid = TRUE;
      }
      else if( SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) )
      {
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_CUTOFF, down) );
      }
      else
         error = TRUE;

      if( iter != NULL )
      {
         SCIP_CALL( SCIPlpiGetIterations(lpi, &it) );
         *iter += it;
      }
      SCIPdebugMessage(" -> down (x%d <= %g): %g\n", col, newub, *down);

      CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_UB, col, oldub) );
      CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );
#ifdef SCIP_DEBUG
      {
         double b;
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_UB, col, &b) );
         assert( b == oldub );
      }
#endif

      if ( success )
      {
         SCIP_CALL( setBase(lpi) );
      }
   }
   else
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_CUTOFF, down) );
      *downvalid = TRUE;
   }

   /* up branch */
   if( !error )
   {
      newlb = EPSFLOOR(psol+1.0, 1e-06);
      if( newlb <= oldub + 0.5 )
      {
         SCIPdebugMessage("strong branching  up (%g) on x%d (%g) with %d iterations\n", newlb, col, psol, itlim);

         CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_LB, col, newlb) );

         SCIP_CALL( SCIPlpiSolveDual(lpi) );
         if( SCIPlpiIsOptimal(lpi) || SCIPlpiIsIterlimExc(lpi) )
         {
            SCIP_CALL( SCIPlpiGetObjval(lpi, up) );
            *upvalid = TRUE;
         }
         else if( SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) )
         {
            CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_CUTOFF, up) );
         }
         else
            error = TRUE;

         if( iter != NULL )
         {
            SCIP_CALL( SCIPlpiGetIterations(lpi, &it) );
            *iter += it;
         }
         SCIPdebugMessage(" -> up  (x%d >= %g): %g\n", col, newlb, *up);

         CHECK_ZERO( lpi->messagehdlr, GRBsetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_LB, col, oldlb) );
         CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) );
#ifdef SCIP_DEBUG
         {
            double b;
            CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrelement(lpi->grbmodel, GRB_DBL_ATTR_LB, col, &b) );
            assert( b == oldlb );
         }
#endif

         if ( success )
         {
            SCIP_CALL( setBase(lpi) );
         }
      }
      else
      {
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_CUTOFF, up) );
         *upvalid = TRUE;
      }
   }

   /* reset iteration limit */
   SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_ITERATIONLIMIT, olditlim) );
   /* CHECK_ZERO( lpi->messagehdlr, GRBupdatemodel(lpi->grbmodel) ); */

   if( error )
   {
      SCIPerrorMessage("LP error in strong branching.\n");
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
   int j;

   assert( iter != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if( iter != NULL )
      *iter = 0;

   for( j = 0; j < ncols; ++j )
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
   int j;

   assert( iter != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if( iter != NULL )
      *iter = 0;

   for( j = 0; j < ncols; ++j )
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

   return (lpi->solstat != -1);
}

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_RETCODE SCIPlpiGetSolFeasibility(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   int algo;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( lpi->solstat >= 1 );

   SCIPdebugMessage("getting solution feasibility\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, &algo) );

   if( primalfeasible != NULL )
   {
      *primalfeasible = (lpi->solstat == GRB_OPTIMAL || (lpi->solstat == GRB_UNBOUNDED && algo == GRB_METHOD_PRIMAL));
   }

   if( dualfeasible != NULL )
   {
      *dualfeasible = (lpi->solstat == GRB_OPTIMAL || (lpi->solstat == GRB_INFEASIBLE && algo == GRB_METHOD_DUAL));
   }


#if 0
   /* @todo: check whether this code is needed anymore (this was the first version) */
   SCIP_Real viol;
   SCIP_Real tol;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->solstat >= 1 );

   SCIPdebugMessage("getting solution feasibility\n");

   if( primalfeasible != NULL )
   {
      if(lpi->solstat != GRB_INF_OR_UNBD && lpi->solstat != GRB_INFEASIBLE)
      {
         /* check whether maximum scaled violation is smaller than feasibility tolerance */
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_CONSTR_SRESIDUAL, &viol) );
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_FEASIBILITYTOL, &tol) );
         *primalfeasible = (viol <= tol) ? TRUE : FALSE;
         SCIPdebugMessage("primal violation: %g  (tol: %g)\n", viol, tol);
      }
      else
         *primalfeasible = FALSE;
   }

   if( dualfeasible != NULL )
   {
      if(lpi->solstat != GRB_UNBOUNDED && lpi->solstat != GRB_INFEASIBLE)
      {
         /* check whether maximum scaled dual violation is smaller than optimality tolerance */
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_DUAL_SRESIDUAL, &viol) );
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblparam(lpi->grbenv, GRB_DBL_PAR_OPTIMALITYTOL, &tol) );
         *dualfeasible = (viol <= tol) ? TRUE : FALSE;
         SCIPdebugMessage("dual violation: %g  (tol: %g)\n", viol, tol);
      }
      else
         *dualfeasible = FALSE;
   }
#endif

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
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_UNBOUNDED);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_UNBOUNDED);
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIP_Bool primalfeasible;
   SCIP_RETCODE retcode;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for primal unboundedness\n");

   primalfeasible = FALSE; /* to fix compiler warning */
   retcode = SCIPlpiGetSolFeasibility(lpi, &primalfeasible, NULL);
   if ( retcode != SCIP_OKAY )
      return FALSE;

   /* Probably GRB_UNBOUNDED means that the problem has an unbounded ray, but not necessarily that a feasible primal solution exists. */
   return (primalfeasible && (lpi->solstat == GRB_UNBOUNDED || lpi->solstat == GRB_INF_OR_UNBD));
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for primal infeasibility\n");

   assert( lpi->solstat != GRB_INF_OR_UNBD );
   return (lpi->solstat == GRB_INFEASIBLE);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int algo;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( lpi->solstat >= 0 );

   SCIPdebugMessage("checking for primal feasibility\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, &algo) );

   return (lpi->solstat == GRB_OPTIMAL || (lpi->solstat == GRB_UNBOUNDED && algo == GRB_METHOD_PRIMAL));
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_INFEASIBLE);  /* ????????? */
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int algo;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( lpi->solstat >= 0 );

   CHECK_ZERO( lpi->messagehdlr, GRBgetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, &algo) );

   return (lpi->solstat == GRB_INFEASIBLE && algo == GRB_METHOD_DUAL);
}

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int algo;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( lpi->solstat >= 0 );

   SCIPdebugMessage("checking for dual unboundedness\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, &algo) );

   return (lpi->solstat == GRB_INFEASIBLE && algo == GRB_METHOD_DUAL);
}

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->solstat >= 0 );

   SCIPdebugMessage("checking for dual infeasibility\n");

   return (lpi->solstat == GRB_UNBOUNDED);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int algo;

   assert( lpi != NULL );
   assert( lpi->grbmodel != NULL );
   assert( lpi->grbenv != NULL );
   assert( lpi->solstat >= 0 );

   SCIPdebugMessage("checking for dual feasibility\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintparam(lpi->grbenv, GRB_INT_PAR_METHOD, &algo) );

   return (lpi->solstat == GRB_OPTIMAL || (lpi->solstat == GRB_INFEASIBLE && algo == GRB_METHOD_DUAL));
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_OPTIMAL);
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("checking for stability: Gurobi solstat = %d\n", lpi->solstat);

   return (lpi->solstat != GRB_NUMERIC);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_CUTOFF);
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_ITERATION_LIMIT);
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   return (lpi->solstat == GRB_TIME_LIMIT);
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   return lpi->solstat;
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(success != NULL);

   *success = FALSE;

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("getting solution's objective value\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_OBJVAL, objval) );

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
   int ncols;
   int nrows;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   SCIPdebugMessage("getting solution\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );
   assert( ncols >= 0 && nrows >= 0 );

   if( objval != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_OBJVAL, objval) );
   }

   if( primsol != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_X, 0, ncols, primsol) );
   }

   if( dualsol != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_PI, 0, nrows, dualsol) );
   }

   if( activity != NULL )
   {
      int i;

      /* first get the values of the slack variables */
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_SLACK, 0, nrows, activity) );

      SCIP_CALL( ensureSidechgMem(lpi, nrows) );

      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_RHS, 0, nrows, lpi->rhsarray) );
      CHECK_ZERO( lpi->messagehdlr, GRBgetcharattrarray(lpi->grbmodel, GRB_CHAR_ATTR_SENSE, 0, nrows, lpi->senarray) );

      for( i = 0; i < nrows; ++i )
      {
         switch(lpi->senarray[i])
         {
         case GRB_LESS_EQUAL:
         case GRB_EQUAL:
            activity[i] = lpi->rhsarray[i] - activity[i];
            break;
         case GRB_GREATER_EQUAL:
            activity[i] = lpi->rhsarray[i] - activity[i];
            break;
         default:
            SCIPerrorMessage("Unkown sense %c.\n", lpi->senarray[i]);
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
         }
      }
   }

   if( redcost != NULL )
   {
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_RC, 0, ncols, redcost) );
   }

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{
   int ncols;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
   assert( ncols >= 0 );

   SCIPdebugMessage("calling Gurobi get primal ray: %d cols\n", ncols);

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_UNBDRAY, 0, ncols, ray) );

   return SCIP_OKAY;
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual Farkas row multipliers */
   )
{
   int nrows;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpi->solstat >= 0);
   assert(dualfarkas != NULL);

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );
   assert( nrows >= 0 );

   SCIPdebugMessage("calling Gurobi dual Farkas: %d rows\n", nrows);

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_FARKASDUAL, 0, nrows, dualfarkas) );

   return SCIP_LPERROR;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(iterations != NULL);

   *iterations = lpi->iterations;

   return SCIP_OKAY;
}

/** gets information about the quality of an LP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The LPI should return SCIP_INVALID for quality, if the requested quantity is not available. 
 */
SCIP_RETCODE SCIPlpiGetRealSolQuality(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPSOLQUALITY     qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   assert(lpi != NULL);
   assert(quality != NULL);

   CHECK_ZERO( lpi->messagehdlr, GRBgetdblattr(lpi->grbmodel, GRB_DBL_ATTR_KAPPA, quality) );

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
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("saving Gurobi basis into %p/%p\n", (void*) cstat, (void*) rstat);

   if( rstat != 0 )
   {
      int i;
      int nrows;

      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );

      for( i = 0; i < nrows; ++i )
      {
         int stat;
         CHECK_ZERO( lpi->messagehdlr, GRBgetintattrelement(lpi->grbmodel, GRB_INT_ATTR_CBASIS, i, &stat) );

         switch( stat )
         {
         case GRB_BASIC:
            rstat[i] = SCIP_BASESTAT_BASIC;
            break;

         case GRB_NONBASIC_LOWER:
            rstat[i] = SCIP_BASESTAT_LOWER;
            break;

         case GRB_NONBASIC_UPPER:
            rstat[i] = SCIP_BASESTAT_UPPER;
            break;

         case GRB_SUPERBASIC:
            rstat[i] = SCIP_BASESTAT_ZERO;
            break;

         default:
            SCIPerrorMessage("invalid basis status %d\n", stat);
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
         }
      }
   }

   if( cstat != 0 )
   {
      int j;
      int ncols;

      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );

      for( j = 0; j < ncols; ++j )
      {
         int stat;
         CHECK_ZERO( lpi->messagehdlr, GRBgetintattrelement(lpi->grbmodel, GRB_INT_ATTR_VBASIS, j, &stat) );

         switch( stat )
         {
         case GRB_BASIC:
            cstat[j] = SCIP_BASESTAT_BASIC;
            break;

         case GRB_NONBASIC_LOWER:
            cstat[j] = SCIP_BASESTAT_LOWER;
            break;

         case GRB_NONBASIC_UPPER:
            cstat[j] = SCIP_BASESTAT_UPPER;
            break;
         case GRB_SUPERBASIC:
            cstat[j] = SCIP_BASESTAT_ZERO;
            break;

         default:
            SCIPerrorMessage("invalid basis status %d\n", stat);
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
   int i, j;
   int nrows, ncols;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(cstat != NULL);
   assert(rstat != NULL);

   SCIPdebugMessage("loading basis %p/%p into Gurobi\n", (void*) cstat, (void*) rstat);

   invalidateSolution(lpi);

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );

   for( i = 0; i < nrows; ++i )
   {
      switch( rstat[i] )
      {
      case SCIP_BASESTAT_BASIC:
         CHECK_ZERO( lpi->messagehdlr, GRBsetintattrelement(lpi->grbmodel, GRB_INT_ATTR_CBASIS, i, GRB_BASIC) );
         break;

      case SCIP_BASESTAT_LOWER:
         CHECK_ZERO( lpi->messagehdlr, GRBsetintattrelement(lpi->grbmodel, GRB_INT_ATTR_CBASIS, i, GRB_NONBASIC_LOWER) );
         break;

      case SCIP_BASESTAT_UPPER:
         CHECK_ZERO( lpi->messagehdlr, GRBsetintattrelement(lpi->grbmodel, GRB_INT_ATTR_CBASIS, i, GRB_NONBASIC_UPPER) );
         break;

      case SCIP_BASESTAT_ZERO:
         CHECK_ZERO( lpi->messagehdlr, GRBsetintattrelement(lpi->grbmodel, GRB_INT_ATTR_CBASIS, i, GRB_SUPERBASIC) );
         break;

      default:
         SCIPerrorMessage("invalid basis status %d\n", rstat[i]);
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   for( j = 0; j < ncols; ++j )
   {
      switch( cstat[j] )
      {
      case SCIP_BASESTAT_BASIC:
         CHECK_ZERO( lpi->messagehdlr, GRBsetintattrelement(lpi->grbmodel, GRB_INT_ATTR_VBASIS, j, GRB_BASIC) );
         break;

      case SCIP_BASESTAT_LOWER:
         CHECK_ZERO( lpi->messagehdlr, GRBsetintattrelement(lpi->grbmodel, GRB_INT_ATTR_VBASIS, j, GRB_NONBASIC_LOWER) );
         break;

      case SCIP_BASESTAT_UPPER:
         CHECK_ZERO( lpi->messagehdlr, GRBsetintattrelement(lpi->grbmodel, GRB_INT_ATTR_VBASIS, j, GRB_NONBASIC_UPPER) );

      case SCIP_BASESTAT_ZERO:
         CHECK_ZERO( lpi->messagehdlr, GRBsetintattrelement(lpi->grbmodel, GRB_INT_ATTR_VBASIS, j, GRB_SUPERBASIC) );
         break;

      default:
         SCIPerrorMessage("invalid basis status %d\n", cstat[j]);
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   int i, j;
   int nrows, ncols;
   int cnt;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("getting basis information\n");

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );

   cnt = 0;
   for( i = 0; i < nrows; ++i )
   {
      int stat;
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattrelement(lpi->grbmodel, GRB_INT_ATTR_CBASIS, i, &stat) );

      if( stat == GRB_BASIC )
         bind[cnt++] = -1 - i;
   }

   for( j = 0; j < ncols; ++j )
   {
      int stat;
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattrelement(lpi->grbmodel, GRB_INT_ATTR_VBASIS, j, &stat) );

      if( stat == GRB_BASIC )
         bind[cnt++] = j;
   }
   assert( cnt == nrows );

   return SCIP_OKAY;
}

/** get dense row of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiGetBInvRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   SVECTOR x;
   int nrows;
   int k;
   int j;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("getting binv-row %d\n", r);

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );

   x.len = 0;
   SCIP_ALLOC( BMSallocMemoryArray(&(x.ind), nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&(x.val), nrows) );

   /* maybe: CHECK_ZERO( lpi->messagehdlr, GRBBinvRowi(lpi->grbmodel, r, &x) ); */
   CHECK_ZERO( lpi->messagehdlr, GRBBinvi(lpi->grbmodel, r, &x) );

   /* size should be at most the number of rows */
   assert( x.len <= nrows );

   k = 0;
   for( j = 0; j < nrows; ++j )
   {
      assert( k <= x.len );
      if ( k < x.len && (x.ind)[k] == j )
         coef[j] = (x.val)[k++];
      else
         coef[j] = 0.0;
   }
   BMSfreeMemoryArray(&(x.val));
   BMSfreeMemoryArray(&(x.ind));

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
   SVECTOR x;
   int nrows;
   int k;
   int j;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("getting binv-col %d\n", c);

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );

   x.len = 0;
   SCIP_ALLOC( BMSallocMemoryArray(&(x.ind), nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&(x.val), nrows) );
   assert(0 <= c && c < nrows );

   /* maybe: CHECK_ZERO( lpi->messagehdlr, GRBBinvColj(lpi->grbmodel, c, &x) ); */
   CHECK_ZERO( lpi->messagehdlr, GRBBinvj(lpi->grbmodel, c, &x) );

   /* size should be at most the number of rows */
   assert( x.len <= nrows );

   k = 0;
   for( j = 0; j < nrows; ++j )
   {
      assert( k <= x.len );
      if ( k < x.len && (x.ind)[k] == j )
         coef[j] = (x.val)[k++];
      else
         coef[j] = 0.0;
   }
   BMSfreeMemoryArray(&(x.val));
   BMSfreeMemoryArray(&(x.ind));

   return SCIP_OKAY;
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{  /*lint --e{715}*/
   SVECTOR x;
   int ncols;
   int k;
   int j;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("getting binv-row %d\n", r);

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );

   x.len = 0;
   SCIP_ALLOC( BMSallocMemoryArray(&(x.ind), ncols) );
   SCIP_ALLOC( BMSallocMemoryArray(&(x.val), ncols) );

   CHECK_ZERO( lpi->messagehdlr, GRBBinvRowi(lpi->grbmodel, r, &x) );

   /* size should be at most the number of columns */
   assert( x.len <= ncols );

   k = 0;
   for( j = 0; j < ncols; ++j )
   {
      assert( k <= x.len );
      if ( k < x.len && (x.ind)[k] == j )
         coef[j] = (x.val)[k++];
      else
         coef[j] = 0.0;
   }
   BMSfreeMemoryArray(&(x.val));
   BMSfreeMemoryArray(&(x.ind));

   return SCIP_OKAY;
}

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{  /*lint --e{715}*/
   SVECTOR x;
   int nrows;
   int k;
   int j;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("getting binv-col %d\n", c);

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );

   x.len = 0;
   SCIP_ALLOC( BMSallocMemoryArray(&(x.ind), nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&(x.val), nrows) );

   CHECK_ZERO( lpi->messagehdlr, GRBBinvColj(lpi->grbmodel, c, &x) );

   /* size should be at most the number of rows */
   assert( x.len <= nrows );

   k = 0;
   for( j = 0; j < nrows; ++j )
   {
      assert( k <= x.len );
      if ( k < x.len && (x.ind)[k] == j )
         coef[j] = (x.val)[k++];
      else
         coef[j] = 0.0;
   }
   BMSfreeMemoryArray(&(x.val));
   BMSfreeMemoryArray(&(x.ind));

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
   SCIP_Bool success;
   int ncols;
   int nrows;

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(lpistate != NULL);

   /* if there is no basis information available, no state can be saved */
   if( !lpi->solisbasic )
   {
      *lpistate = NULL;
      return SCIP_OKAY;
   }

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
   assert(ncols >= 0);
   assert(nrows >= 0);

   /* get unpacked basis information from Gurobi */
   SCIP_CALL( getBase(lpi, &success) );

   if ( success )
   {
      SCIPdebugMessage("storing Gurobi LPI state in %p (%d cols, %d rows)\n", (void*) *lpistate, ncols, nrows);

      /* allocate lpistate data */
      SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows) );
      (*lpistate)->ncols = ncols;
      (*lpistate)->nrows = nrows;

      /* pack LPi state data */
      lpistatePack(*lpistate, lpi->cstat, lpi->rstat);
   }
   else
   {
      /* In this case no basis information is available. Since SCIP expects the information to work
         in any case, we allocate the lpistate, but do not use the packed information. This might
         happen if the model is infeasible, since Gurobi currently does not return basis information
         in this case. */
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, lpistate) );
      (*lpistate)->ncols = ncols;
      (*lpistate)->nrows = nrows;
      (*lpistate)->packrstat = NULL;
      (*lpistate)->packcstat = NULL;
   }

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
   int ncols;
   int nrows;
   int i;

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   /* if there was no basis information available, the LPI state was not stored */
   if( lpistate == NULL || lpistate->packrstat == NULL || lpistate->packcstat )
      return SCIP_OKAY;

   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMCONSTRS, &nrows) );
   CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_NUMVARS, &ncols) );
   assert(lpistate->ncols <= ncols);
   assert(lpistate->nrows <= nrows);

   SCIPdebugMessage("loading LPI state %p (%d cols, %d rows) into Gurobi LP with %d cols and %d rows\n",
      (void*) lpistate, lpistate->ncols, lpistate->nrows, ncols, nrows);

   if( lpistate->ncols == 0 || lpistate->nrows == 0 )
      return SCIP_OKAY;

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, ncols) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   /* unpack LPi state data */
   lpistateUnpack(lpistate, lpi->cstat, lpi->rstat);

   /* extend the basis to the current LP beyond the previously existing columns */
   for( i = lpistate->ncols; i < ncols; ++i )
   {
      SCIP_Real bnd;
      CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_LB, i, i, &bnd) );
      if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         CHECK_ZERO( lpi->messagehdlr, GRBgetdblattrarray(lpi->grbmodel, GRB_DBL_ATTR_UB, i, i, &bnd) );
         if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
            lpi->cstat[i] = SCIP_BASESTAT_ZERO;  /* variable is free */
         else
            lpi->cstat[i] = SCIP_BASESTAT_UPPER; /* use finite upper bound */
      }
      else
         lpi->cstat[i] = SCIP_BASESTAT_LOWER;    /* use finite lower bound */
   }
   for( i = lpistate->nrows; i < nrows; ++i )
      lpi->rstat[i] = SCIP_BASESTAT_BASIC;

   /* load basis information into Gurobi */
   SCIP_CALL( setBase(lpi) );

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   /**@todo implement SCIPlpiClearState() for Gurobi */
   SCIPmessagePrintWarning(lpi->messagehdlr, "Gurobi interface does not implement SCIPlpiClearState()\n");

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
      lpistateFree(lpistate, blkmem);

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{  /*lint --e{715}*/
   return (lpistate != NULL);
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("reading LP state from file <%s>\n", fname);

   SCIPerrorMessage("SCIPlpiReadState() not supported by Gurobi\n");

   return SCIP_LPERROR;
}

/** writes LP state (like basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("writing LP state to file <%s>\n", fname);

   SCIPerrorMessage("SCIPlpiWriteState() not supported by Gurobi\n");

   return SCIP_LPERROR;
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
   int temp;
   SCIP_Real dtemp;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(ival != NULL);

   SCIPdebugMessage("getting int parameter %d\n", type);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_FASTMIP:
      /* maybe set perturbation */
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_SCALING:
      SCIP_CALL( getIntParam(lpi, GRB_INT_PAR_SCALEFLAG, &temp) );
      assert( temp == 0 || temp == 1 );
      *ival = (temp == 1) ? TRUE : FALSE;
      break;
   case SCIP_LPPAR_PRESOLVING:
      SCIP_CALL( getIntParam(lpi, GRB_INT_PAR_PRESOLVE, &temp) );
      assert( temp == GRB_PRESOLVE_AUTO || temp == GRB_PRESOLVE_OFF || temp == GRB_PRESOLVE_CONSERVATIVE || temp == GRB_PRESOLVE_AGGRESSIVE );
      *ival = (temp == GRB_PRESOLVE_OFF) ? FALSE : TRUE;
      break;
   case SCIP_LPPAR_PRICING:
      *ival = (int) lpi->pricing;
      break;
   case SCIP_LPPAR_LPINFO:
      SCIP_CALL( getIntParam(lpi, GRB_INT_PAR_OUTPUTFLAG, &temp) );
      assert( temp == 0 || temp == 1 );
      *ival = (temp == 1) ? TRUE : FALSE;
      break;
   case SCIP_LPPAR_LPITLIM:
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_ITERATIONLIMIT, &dtemp) );
      assert( dtemp >= 0.0 );
      if( dtemp >= GRB_INFINITY )
         *ival = INT_MAX;
      else
         *ival = (int) dtemp;
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
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("setting int parameter %d to %d\n", type, ival);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      assert(ival == TRUE || ival == FALSE);
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_FASTMIP:
      assert(ival == TRUE || ival == FALSE);
      return SCIP_PARAMETERUNKNOWN;
      break;
   case SCIP_LPPAR_SCALING:
      assert(ival == TRUE || ival == FALSE);
      if( ival )
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SCALEFLAG, 1) );
      else
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SCALEFLAG, 0) );
      break;
   case SCIP_LPPAR_PRESOLVING:
      assert(ival == TRUE || ival == FALSE);
      if( ival )
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_AUTO) );
      else
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF) );
      break;
   case SCIP_LPPAR_PRICING:
      lpi->pricing = (SCIP_PRICING)ival;
      switch( (SCIP_PRICING)ival )
      {
      case SCIP_PRICING_LPIDEFAULT:
      case SCIP_PRICING_AUTO:
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SIMPLEXPRICING, GRB_SIMPLEXPRICING_AUTO) );
         break;
      case SCIP_PRICING_FULL:
         /* full does not seem to exist -> use auto */
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SIMPLEXPRICING, GRB_SIMPLEXPRICING_AUTO) );
         break;
      case SCIP_PRICING_PARTIAL:
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SIMPLEXPRICING, GRB_SIMPLEXPRICING_PARTIAL) );
         break;
      case SCIP_PRICING_STEEP:
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SIMPLEXPRICING, GRB_SIMPLEXPRICING_STEEPEST_EDGE) );
         break;
      case SCIP_PRICING_STEEPQSTART:
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SIMPLEXPRICING, GRB_SIMPLEXPRICING_STEEPEST_QUICK) );
         break;
      case SCIP_PRICING_DEVEX:
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_SIMPLEXPRICING, GRB_SIMPLEXPRICING_DEVEX) );
         break;
      default:
         return SCIP_LPERROR;
      }
      break;
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      if( ival )
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_OUTPUTFLAG, 1) );
      else
         SCIP_CALL( setIntParam(lpi, GRB_INT_PAR_OUTPUTFLAG, 0) );
      break;
   case SCIP_LPPAR_LPITLIM:
      SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_ITERATIONLIMIT, (double) ival) );
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
   int objsen;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);
   assert(dval != NULL);

   SCIPdebugMessage("getting real parameter %d\n", type);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_FEASIBILITYTOL, dval) );
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_OPTIMALITYTOL, dval) );
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:
      return SCIP_PARAMETERUNKNOWN;
      break;
   case SCIP_LPPAR_LOBJLIM:
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_MODELSENSE, &objsen) );
      if( objsen == 1 )
         SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_CUTOFF, dval) );
      else
         return SCIP_PARAMETERUNKNOWN;
      break;
   case SCIP_LPPAR_UOBJLIM:
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_MODELSENSE, &objsen) );
      if( objsen == 0 )
         SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_CUTOFF, dval) );
      else
         return SCIP_PARAMETERUNKNOWN;
      break;
   case SCIP_LPPAR_LPTILIM:
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_TIMELIMIT, dval) );
      break;
   case SCIP_LPPAR_MARKOWITZ:
      SCIP_CALL( getDblParam(lpi, GRB_DBL_PAR_MARKOWITZTOL, dval) );
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
   int objsen;

   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("setting real parameter %d to %g\n", type, dval);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_FEASIBILITYTOL, dval) );
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_OPTIMALITYTOL, dval) );
      break;
   case SCIP_LPPAR_BARRIERCONVTOL:
      return SCIP_PARAMETERUNKNOWN;
      break;
   case SCIP_LPPAR_LOBJLIM:
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_MODELSENSE, &objsen) );
      if( objsen == 1 )
         SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_CUTOFF, dval) );
      break;
   case SCIP_LPPAR_UOBJLIM:
      CHECK_ZERO( lpi->messagehdlr, GRBgetintattr(lpi->grbmodel, GRB_INT_ATTR_MODELSENSE, &objsen) );
      if( objsen == 0 )
         SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_CUTOFF, dval) );
      break;
   case SCIP_LPPAR_LPTILIM:
      SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_TIMELIMIT, dval) );
      break;
   case SCIP_LPPAR_MARKOWITZ:
      SCIP_CALL( setDblParam(lpi, GRB_DBL_PAR_MARKOWITZTOL, dval) );
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
   return GRB_INFINITY;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{  /*lint --e{715}*/
   return (val >= GRB_INFINITY);
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
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("reading LP from file <%s>\n", fname);

   CHECK_ZERO( lpi->messagehdlr, GRBread(lpi->grbmodel, fname) );

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->grbmodel != NULL);

   SCIPdebugMessage("writing LP to file <%s>\n", fname);

   CHECK_ZERO( lpi->messagehdlr, GRBwrite(lpi->grbmodel, fname) );

   return SCIP_OKAY;
}

/**@} */
