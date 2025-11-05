/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   lpiexact_qsoex.c
 * @ingroup LPIEXACTS
 * @brief  exact LP interface for QSopt_ex version >= 2.5.4 (r239)
 * @author Leon Eifler
 * @author Daniel Espinoza
 * @author Marc Pfetsch
 * @author Kati Wolter
*/

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define VERIFY_OUT */ /* uncomment to get info of QSopt_ex about verifying dual feasibility of the basis */

/* #define USEOBJLIM */  /* uncomment to pass objlimit to exact lp solver; same as in cons_exactlinear.c;  warning: QSopt_ex allows objlimits but the support is buggy; if the limit is reached, QSopt_ex does not stop but increasess the precision */

#include "scip/def.h"

#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_EGLIB)
#include "EGlib.h"
#include "QSopt_ex.h"
#endif

#include <string.h>

#include "scip/bitencode.h"
#include "scip/scip_mem.h"
#include "lpiexact/type_lpiexact.h"
#include "lpiexact/lpiexact.h"
#include "scip/pub_message.h"


/** solver name */
static char __qsstr[1024];
static char __egstr[1024];

#ifdef SCIP_WITH_GMP

typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE

/** LP interface */
struct SCIP_LPiExact
{
#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_EGLIB)
   mpq_QSprob            prob;               /**< LP struct pointer */
   mpq_factor_work*      factor;             /**< factorized matrix  */
   int                   solstat;            /**< solution status of last optimization call */
   int                   previt;             /**< previous number of simplex iterations performed */
   int                   pricing;            /**< SCIP pricing option */
#endif
   int                   rowspace;           /**< current size of internal row-related arrays */
   char*                 isen;               /**< array of length rowspace */
   mpq_t*                irhs;               /**< array of rhs rowspace */
   mpq_t*                irng;               /**< array of range rowspace */
   int*                  ircnt;              /**< array of count rowspace */
   int*                  irbeg;              /**< array of begining index rowspace */
   int                   colspace;           /**< current size of internal column-related arrays */
   int*                  iccnt;              /**< array of length colspace */
   char*                 iccha;              /**< array of type colspace */
   int                   tbsz;               /**< current size of tableau-related arrays */
   mpq_t*                itab;               /**< array of length tbsz */
   char*                 ibas;               /**< array of length tbsz */
};



/*
 * local defines
 */

/*lint --e{750} */
/** print location of the calling line */
#define __QS_PRINTLOC__ fprintf(stderr,", in (%s:%d)\n", __FILE__, __LINE__);

/** This macro is to print error messages and jump to the given point in the code, it also print the
 *  file and line where this happend */
#define QS_TESTG(A,B,...) do { { \
   if (A){ \
      fprintf(stderr,__VA_ARGS__); \
      __QS_PRINTLOC__; \
      goto B; } } } while(0)

/** This macro is to print error messages and to exit with SCIP_LPERROR */
#define QS_ERROR(A,...) do { { \
   if (A){ \
      fprintf(stderr,__VA_ARGS__); \
      __QS_PRINTLOC__; \
      return SCIP_LPERROR; } } } while(0)

/** return value macro, if the value is non-zero, write to standard error the returning code and
 *  where this happened, and return SCIP_ERROR, otherwise return normal SCIP_OKAY termination code. */
#define QS_RETURN(A) do { \
      const int __RVAL__ = (A); \
      if (__RVAL__){ \
         fprintf(stderr,"LP Error: QSopt_ex returned %d",__RVAL__); \
         __QS_PRINTLOC__; \
         return SCIP_ERROR; } \
      return SCIP_OKAY; } while(0)

/** return value macro, if the value is non-zero, write to standard error the returning code and
 *  where this happened, and return SCIP_ERROR, otherwise do nothing. */
#define QS_CONDRET(A) do { \
      const int __RVAL__ = (A); \
      if (__RVAL__){ \
         fprintf(stderr,"LP Error: QSopt_ex returned %d",__RVAL__); \
         __QS_PRINTLOC__; \
         return SCIP_LPERROR; } \
   } while(0)



/*
 * LPi state methods
 */


/** LPi state stores basis information */
struct SCIP_LPiState
{
   int                   ncols;              /**< number of LP columns */
   int                   nrows;              /**< number of LP rows */
   COLPACKET*            packcstat;          /**< column basis status in compressed form */
   ROWPACKET*            packrstat;          /**< row basis status in compressed form */
};

static
void printGMP(
   const mpq_t           val
   )
{
   char* buffer;
   buffer = (char*) malloc(mpz_sizeinbase(mpq_numref(val), 10) + mpz_sizeinbase(mpq_denref(val), 10) + 3);
   (void)mpq_get_str(buffer, 10, val);
   printf("%s \n", buffer);
   free(buffer);
}

/** returns the number of packets needed to store column packet information */
static
int colpacketNum(
   int                   ncols               /**< number of columns to store */
   )
{
   return (ncols + (int)COLS_PER_PACKET-1)/(int)COLS_PER_PACKET;
}

/** returns the number of packets needed to store row packet information */
static
int rowpacketNum(
   int                   nrows               /**< number of rows to store */
   )
{
   return (nrows + (int)ROWS_PER_PACKET-1)/(int)ROWS_PER_PACKET;
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
 * local functions
 */

/** ensure size of column-related arrays */
static inline
SCIP_RETCODE ensureTabMem(
   SCIP_LPIEXACT*        lpi,
   int const             sz
   )
{ /*lint --e{647}*/
   register int i;
   if( lpi->tbsz < sz )
   {
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->itab), sz*2) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->ibas), sz*2) );
      for( i = lpi->tbsz ; i < sz * 2 ; i++ )
         mpq_init(lpi->itab[i]);
      lpi->tbsz = sz*2;
   }
   return SCIP_OKAY;
}

/** ensure size of column-related arrays */
static inline
SCIP_RETCODE ensureColMem(
   SCIP_LPIEXACT*        lpi,
   int const             ncols
   )
{
   if( lpi->colspace < ncols )
   {
      lpi->colspace = ncols * 2;
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->iccnt), lpi->colspace) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->iccha), lpi->colspace) );
   }
   return SCIP_OKAY;
}

/** ensure size of row-related arrays */
static inline
SCIP_RETCODE ensureRowMem(
   SCIP_LPIEXACT*        lpi,
   int const             nrows
   )
{  /*lint --e{647}*/
   register int i;
   if( lpi->rowspace < nrows )
   {
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->isen), nrows * 2) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->irhs), nrows * 2) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->irng), nrows * 2) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->ircnt), nrows * 2) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->irbeg), nrows * 2) );
      for (i = lpi->rowspace ; i < nrows * 2; i++)
      {
         mpq_init(lpi->irhs[i]);
         mpq_init(lpi->irng[i]);
      }
      lpi->rowspace = nrows * 2;
   }
   return SCIP_OKAY;
}

/** transform lhs/rhs into qsopt format */
static inline
SCIP_RETCODE convertSides(
   SCIP_LPIEXACT* const  lpi,
   int const             nrows,
   SCIP_RATIONAL**       lhs,
   SCIP_RATIONAL**       rhs
   )
{  /*lint --e{663, 550, 438, 528}*/
   int state;
   register int i;

   for( i = 0; i < nrows; ++i )
   {
      mpq_t* lhsg;
      mpq_t* rhsg;
      int state1;
      int state2;

      lhsg = SCIPrationalGetGMP(lhs[i]);
      rhsg = SCIPrationalGetGMP(rhs[i]);
#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_EGLIB)
      state1 = ((mpq_cmp(*lhsg, mpq_ILL_MINDOUBLE) <= 0) ? 1U : 0U);
      state2 = ((mpq_cmp(*rhsg, mpq_ILL_MAXDOUBLE) >= 0) ? 2U : 0U);
      state = state1 | state2;
#else
      state1 = 0;
      state2 = 0;
      state = 0;
#endif
      /* state = (((mpq_cmp(*lhsg, mpq_ILL_MINDOUBLE) <= 0) ? 1U : 0U) | ((mpq_cmp(*rhsg, mpq_ILL_MAXDOUBLE) >= 0) ? 2U : 0U)); */
      lpi->ircnt[i] = 0;
      lpi->irbeg[i] = 0;
      switch( state )
      {
      case 0:
         if( SCIPrationalIsEQ(lhs[i], rhs[i]) )
         {
            lpi->isen[i] = 'E';
            mpq_set(lpi->irhs[i], *lhsg);
            mpq_set_ui(lpi->irng[i], 0UL, 1UL);
         }
         else
         {
            lpi->isen[i] = 'R';
            mpq_set(lpi->irhs[i], *lhsg);
            mpq_sub(lpi->irng[i], *rhsg, *lhsg);
            assert( mpq_sgn(lpi->irng[i]) >=0 );
         }
         break;
      case 1:
         lpi->isen[i] = 'L';
         mpq_set(lpi->irhs[i], *rhsg);
         mpq_set_ui(lpi->irng[i], 0UL, 1UL);
         break;
      case 2:
         lpi->isen[i] = 'G';
         mpq_set(lpi->irhs[i], *lhsg);
         mpq_set_ui(lpi->irng[i], 0UL, 1UL);
         break;
      default:
         SCIPerrorMessage("Error, constraint %d has no bounds!",i);
         SCIPABORT();
      }
   }
   return SCIP_OKAY;
 }  /*lint --e{528}*/


/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */


#endif
/** gets name and version of LP solver */
const char* SCIPlpiExactGetSolverName(
   void
   )
{
#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_EGLIB)
   sprintf(__qsstr, "QSopt_ex %s", string_QSopt_ex);
#else
   sprintf(__qsstr, "QSopt_ex");
#endif

   return __qsstr;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiExactGetSolverDesc(
   void
   )
{
   return "Exact Linear Programming Solver by D. Espinoza, W. Cook, S. Dash, and D. Applegate (dii.uchile.cl/~daespino/QSoptExact_doc/main.html)";
}

/** gets name and version of external package required for LP solver */
const char* SCIPlpiExactGetExternalCodeName(
   void
   )
{
#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_EGLIB)
   sprintf(__egstr, "EGlib %s", string_EGlib);
#else
   sprintf(__egstr, "EGlib");
#endif

   return __egstr;
}

/** gets description of external package required for LP solver (developer, webpage, ...) */
const char* SCIPlpiExactGetExternalCodeDesc(
   void
   )
{
   return "Library for basic structures and utilities by D. Espinoza and M. Goycoolea (dii.uchile.cl/~daespino/EGlib_doc/main.html)";
}

#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_EGLIB)

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiExactGetSolverPointer(
   SCIP_LPIEXACT*        lpi                 /**< pointer to an LP interface structure */
   )
{
   return (void*) lpi->prob;
}

/**@} */


/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** calls initializator of LP solver (if this did not already happen); this is mainly needed for defining constants in extended and rational precision */
void SCIPlpiExactStart(
   void
   )
{
   if( !__QSexact_setup )
      QSexactStart();
}

/** calls deinitializator of LP solver; this is needed for freeing all internal data of the solver, like constants in
 *  extended and rational precision
 */
void SCIPlpiExactEnd(
   void
   )
{
   if( __QSexact_setup )
      QSexactClear();
}

/** creates an LP problem object
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPlpiExactCreate(
   SCIP_LPIEXACT**       lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   register int i;

   /* QSopt_ex only works with bools as integers */
   assert(sizeof (SCIP_Bool) == sizeof (int));
   assert(lpi != NULL);

   SCIPdebugMessage("SCIPlpiExactCreate()\n");

   /* create LP */
   SCIP_ALLOC( BMSallocMemory(lpi) );
   SCIPlpiExactStart();
   memset(*lpi, 0, sizeof(struct SCIP_LPiExact));

   /* factor work is NULL unless used */
   (*lpi)->factor =  (mpq_factor_work*) NULL;

   (*lpi)->prob = mpq_QScreate_prob(name, (int) objsen);
   if( (*lpi)->prob == NULL )
   {
      SCIPerrorMessage("No memory\n");
      return SCIP_LPERROR;
   }

   (*lpi)->rowspace = 1024;
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->isen), 1024) );
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->irhs), 1024) );
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->irng), 1024) );
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->irbeg), 1024) );
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->ircnt), 1024) );

   (*lpi)->colspace = 1024;
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->iccnt), 1024) );
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->iccha), 1024) );

   (*lpi)->tbsz = 1024;
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->itab), 1024) );
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->ibas), 1024) );
   for( i = 0; i < 1024; i++ )
   {
      mpq_init((*lpi)->irhs[i]);
      mpq_init((*lpi)->irng[i]);
      mpq_init((*lpi)->itab[i]);
   }

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiExactFree(
   SCIP_LPIEXACT**       lpi                 /**< pointer to an LP interface structure */
   )
{
   register int i;

   assert(lpi != NULL);
   assert(*lpi != NULL);

   SCIPdebugMessage("SCIPlpiExactFree()\n");

   /* free factor work */
   if( (*lpi)->factor != NULL )
   {
      mpq_ILLfactor_free_factor_work((*lpi)->factor);
      BMSfreeMemoryArray( &((*lpi)->factor) );
   }

   /* free LP */
   mpq_QSfree_prob((*lpi)->prob);
   for( i = 0; i < (*lpi)->tbsz; ++i )
      mpq_clear((*lpi)->itab[i]);
   for( i = 0; i < (*lpi)->rowspace; ++i )
   {
      mpq_clear((*lpi)->irng[i]);
      mpq_clear((*lpi)->irhs[i]);
   }
   BMSfreeMemoryArray(&((*lpi)->isen));
   BMSfreeMemoryArray(&((*lpi)->irhs));
   BMSfreeMemoryArray(&((*lpi)->irng));
   BMSfreeMemoryArray(&((*lpi)->ircnt));
   BMSfreeMemoryArray(&((*lpi)->irbeg));
   BMSfreeMemoryArray(&((*lpi)->iccnt));
   BMSfreeMemoryArray(&((*lpi)->iccha));
   BMSfreeMemoryArray(&((*lpi)->itab));
   BMSfreeMemoryArray(&((*lpi)->ibas));

   /* free memory */
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
   SCIP_RATIONAL**       obj,                /**< objective function values of columns */
   SCIP_RATIONAL**       lb,                 /**< lower bounds of columns */
   SCIP_RATIONAL**       ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   SCIP_RATIONAL**       lhs,                /**< left hand sides of rows */
   SCIP_RATIONAL**       rhs,                /**< right hand sides of rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each column in ind- and val-array */
   int*                  ind,                /**< row indices of constraint matrix entries */
   SCIP_RATIONAL**       val                 /**< values of constraint matrix entries */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("loading LP in column format into QSopt_ex: %d cols, %d rows\n", ncols, nrows);

   /* delete old LP */
   SCIP_CALL( SCIPlpiExactClear(lpi) );

   /* set sense */
   if( objsen == SCIP_OBJSEN_MAXIMIZE )
   {
      rval = mpq_QSchange_objsense(lpi->prob, QS_MAX);
      QS_CONDRET(rval);
   }
   else
   {
      rval = mpq_QSchange_objsense(lpi->prob, QS_MIN);
      QS_CONDRET(rval);
   }

   /* add rows with no matrix, and then the columns, first ensure space */
   SCIP_CALL( ensureRowMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   SCIP_CALL( convertSides(lpi, nrows, lhs, rhs) );

   /* now we add the rows */
   rval = mpq_QSadd_ranged_rows(lpi->prob, nrows, lpi->ircnt, lpi->irbeg, 0, (const mpq_t*) 0, (const mpq_t*) lpi->irhs,
      lpi->isen, (const mpq_t*) lpi->irng, (const char**)rownames);
   QS_CONDRET(rval);

   SCIPlpiExactAddCols(lpi, ncols, obj, lb, ub, colnames, nnonz, beg, ind, val);

   QS_RETURN(rval);
}

/** adds columns to the LP */
SCIP_RETCODE SCIPlpiExactAddCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   SCIP_RATIONAL**       obj,                /**< objective function values of new columns */
   SCIP_RATIONAL**       lb,                 /**< lower bounds of new columns */
   SCIP_RATIONAL**       ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   int*                  beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   int*                  ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   SCIP_RATIONAL**       val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   mpq_t* valgmp;
   int rval = 0;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("adding %d columns with %d nonzeros to QSopt_ex\n", ncols, nnonz);

   lpi->solstat = 0;

   /* ensure column size */
   SCIP_CALL( ensureColMem(lpi, ncols) );

   if( nnonz > 0 )
   {
      /* compute column lengths */
      for( i = 0; i < ncols - 1; ++i )
      {
         lpi->iccnt[i] = beg[i+1] - beg[i];
         assert(lpi->iccnt[i] >= 0);
      }
      if( ncols > 0 )
      {
         lpi->iccnt[ncols-1] = nnonz - beg[ncols-1];
         assert( lpi->iccnt[ncols-1] >= 0 );
      }

      SCIP_ALLOC( BMSallocMemoryArray(&valgmp, nnonz) );
      SCIPrationalSetGMPArray(valgmp, val, nnonz);

      /* and add the columns */
      for( i = 0; i < ncols; ++i )
      {
         mpq_QSadd_col(lpi->prob, lpi->iccnt[i], &ind[beg[i]], &(valgmp[beg[i]]),
            *SCIPrationalGetGMP(obj[i]), *SCIPrationalGetGMP(lb[i]), *SCIPrationalGetGMP(ub[i]), (const char*) colnames[i]);
      }

      SCIPrationalClearArrayGMP(valgmp, nnonz);
      BMSfreeMemoryArray(&valgmp);
   }
   else
   {
      for( i = 0; i < ncols; ++i )
      {
         rval = mpq_QSnew_col(lpi->prob, *SCIPrationalGetGMP(obj[i]), *SCIPrationalGetGMP(lb[i]), *SCIPrationalGetGMP(ub[i]), (const char*) colnames[i]);
      }
   }

   QS_RETURN(rval);
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiExactDelCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   const int len = lastcol - firstcol +1;
   int rval = 0;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   assert(0 <= firstcol && len > 0 && lastcol < mpq_QSget_colcount (lpi->prob));

   SCIPdebugMessage("deleting %d columns from QSopt_ex\n", len);

   SCIP_CALL( ensureColMem(lpi, len) );
   for( i = firstcol ; i <= lastcol ; i++ )
      lpi->iccnt[i-firstcol] = i;

   rval = mpq_QSdelete_cols(lpi->prob, len, lpi->iccnt);

   QS_RETURN(rval);
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiExactDelColset(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   int rval = 0, ncols, ccnt;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   ncols = mpq_QSget_colcount(lpi->prob);
   lpi->solstat = 0;

   SCIPdebugMessage("deleting a column set from QSopt_ex\n");

   rval = mpq_QSdelete_setcols(lpi->prob,dstat);
   QS_CONDRET(rval);

   for( i = 0, ccnt = 0; i < ncols; i++ )
   {
      if( dstat[i] )
         dstat[i] = -1;
      else
         dstat[i] = ccnt++;
   }
   return SCIP_OKAY;
}

/** adds rows to the LP */
SCIP_EXPORT
SCIP_RETCODE SCIPlpiExactAddRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   SCIP_RATIONAL**       lhs,                /**< left hand sides of new rows */
   SCIP_RATIONAL**       rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   int*                  beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   int*                  ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   SCIP_RATIONAL**       val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   register int i;
   int rval = 0;
   mpq_t* valgmp;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("adding %d rows with %d nonzeros to QSopt_ex\n", nrows, nnonz);

   /* add rows with no matrix, and then the columns, first ensure space */
   SCIP_CALL( ensureRowMem (lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   SCIP_CALL( convertSides(lpi, nrows, lhs, rhs) );

   /* compute row count */
   for( i = 0; i < nrows - 1; i++ )
   {
      lpi->ircnt[i] = beg[i + 1] - beg[i];
      assert(lpi->ircnt[i] >= 0);
   }
   if( nrows > 0 )
   {
      lpi->ircnt[nrows - 1] = nnonz - beg[nrows - 1];
      assert(lpi->ircnt[nrows - 1] >= 0);
   }

   SCIP_ALLOC( BMSallocMemoryArray(&valgmp, nnonz) );
   SCIPrationalSetGMPArray(valgmp, val, nnonz);

   rval = mpq_QSadd_ranged_rows(lpi->prob, nrows, lpi->ircnt, beg, ind, (const mpq_t*) valgmp, (const mpq_t*) lpi->irhs,
      lpi->isen, (const mpq_t*) lpi->irng, (const char**) rownames);
   QS_CONDRET(rval);

   SCIPrationalClearArrayGMP(valgmp, nnonz);
   BMSfreeMemoryArray(&valgmp);

   return SCIP_OKAY;
}


/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiExactDelRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   const int len = lastrow - firstrow +1;
   int rval = 0;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   assert(0 <= firstrow && len > 0 && lastrow < mpq_QSget_rowcount (lpi->prob));

   SCIPdebugMessage("deleting %d rows from QSopt_ex\n", len);

   SCIP_CALL( ensureRowMem(lpi, len) );
   for( i = firstrow; i <= lastrow; i++ )
      lpi->ircnt[i - firstrow] = i;
   rval = mpq_QSdelete_rows(lpi->prob, len, lpi->ircnt);

   QS_RETURN(rval);
}


/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiExactDelRowset(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   int rval = 0, nrows, ccnt, ndel=0;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   nrows = mpq_QSget_rowcount(lpi->prob);
   lpi->solstat = 0;

   for( i = 0; i < nrows; ++i )
   {
      if( dstat[i] == 1 )
         ndel++;
   }

   SCIPdebugMessage("deleting a row set from QSopt_ex (%d)\n",ndel);

   rval = mpq_QSdelete_setrows(lpi->prob,dstat);
   QS_CONDRET(rval);

   for( i = 0, ccnt = 0; i < nrows; ++i )
   {
      if( dstat[i] )
         dstat[i] = -1;
      else
         dstat[i] = ccnt++;
   }
   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiExactClear(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   register int i;
   int ncols, nrows, rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("clearing QSopt_ex LP\n");
   lpi->solstat = 0;

   ncols = mpq_QSget_colcount(lpi->prob);
   nrows = mpq_QSget_rowcount(lpi->prob);
   if( ncols >=  1 )
   {
      SCIP_CALL( ensureColMem(lpi,ncols) );
      for (i = 0; i < ncols; ++i)
         lpi->iccnt[i] = i;
      rval = mpq_QSdelete_cols(lpi->prob, ncols, lpi->iccnt);
      QS_CONDRET(rval);
   }

   if( nrows >= 1 )
   {
      SCIP_CALL( ensureRowMem(lpi, nrows) );
      for (i = 0; i < nrows; ++i)
         lpi->ircnt[i] = i;
      rval = mpq_QSdelete_rows(lpi->prob, nrows, lpi->ircnt);
      QS_CONDRET(rval);
   }
   QS_RETURN(rval);
}


/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiExactChgBounds(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   int*                  ind,                /**< column indices */
   SCIP_RATIONAL**       lb,                 /**< values for the new lower bounds, or NULL */
   SCIP_RATIONAL**       ub                  /**< values for the new upper bounds, or NULL */
   )
{
   register int i;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(lb != NULL || ub != NULL);

   lpi->solstat = 0;

   SCIPdebugMessage("changing %d bounds in QSopt_ex\n", ncols);
#ifdef SCIP_DEBUG
   {
      int j;
      char s[SCIP_MAXSTRLEN];

      for( j = 0; j < ncols; ++j )
      {
         if( lb == NULL)
            gmp_snprintf(s, SCIP_MAXSTRLEN, "  col %d: [--,%Qd]\n", ind[j], *SCIPrationalGetGMP(ub[j]));
         else if( ub == NULL )
            gmp_snprintf(s, SCIP_MAXSTRLEN, "  col %d: [%Qd,--]\n", ind[j], *SCIPrationalGetGMP(lb[j]));
         else
            gmp_snprintf(s, SCIP_MAXSTRLEN, "  col %d: [%Qd,%Qd]\n", ind[j], *SCIPrationalGetGMP(lb[j]), *SCIPrationalGetGMP(ub[j]));
         SCIPdebugPrintf(s);
      }
   }
#endif

   SCIP_CALL( ensureColMem(lpi, ncols) );

   if( lb != NULL )
   {
      for( i = 0; i < ncols; i++ )
      {
         lpi->iccha[i] = 'L';
         rval = mpq_QSchange_bound(lpi->prob, ind[i], lpi->iccha[i], *SCIPrationalGetGMP(lb[i]));
         QS_CONDRET(rval);
      }
   }

   if( ub != NULL )
   {
      for( i = 0; i < ncols; i++ )
      {
         lpi->iccha[i] = 'U';
         rval = mpq_QSchange_bound(lpi->prob, ind[i], lpi->iccha[i], *SCIPrationalGetGMP(ub[i]));
         QS_CONDRET(rval);
      }
   }
   QS_RETURN(rval);
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiExactChgSides(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   int*                  ind,                /**< row indices */
   SCIP_RATIONAL**       lhs,                /**< new values for left hand sides */
   SCIP_RATIONAL**       rhs                 /**< new values for right hand sides */
   )
{
   register int i;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("changing %d sides in QSopt_ex\n", nrows);

   SCIP_CALL( ensureRowMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   SCIP_CALL( convertSides(lpi, nrows, lhs, rhs) );

   /* now we change all rows */
   for( i = 0; i < nrows; ++i )
   {
      rval = mpq_QSchange_sense(lpi->prob, ind[i], lpi->isen[i]);
      QS_CONDRET(rval);

      rval = mpq_QSchange_rhscoef(lpi->prob, ind[i], lpi->irhs[i]);
      QS_CONDRET(rval);

      if (lpi->isen[i] == 'R')
      {
         rval = mpq_QSchange_range(lpi->prob, ind[i], lpi->irng[i]);
         QS_CONDRET(rval);
      }
   }

   QS_RETURN(rval);
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiExactChgCoef(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_RATIONAL*        newval              /**< new value of coefficient */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;

   SCIPdebugMessage("changing coefficient row %d, column %d in QSopt_ex to %g\n", row, col, SCIPrationalGetReal(newval));

   rval = mpq_QSchange_coef(lpi->prob, row, col, *SCIPrationalGetGMP(newval));

   QS_RETURN(rval);
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiExactChgObjsen(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("changing objective sense in QSopt_ex to %d\n", objsen);

   /* set sense */
   if( objsen == SCIP_OBJSEN_MAXIMIZE )
   {
      rval = mpq_QSchange_objsense(lpi->prob, QS_MAX);
      QS_CONDRET(rval);
   }
   else
   {
      rval = mpq_QSchange_objsense(lpi->prob, QS_MIN);
      QS_CONDRET(rval);
   }

   QS_RETURN(rval);
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiExactChgObj(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   SCIP_RATIONAL**       obj                 /**< new objective values for columns */
   )
{
   register int i;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("changing %d objective values in QSopt_ex\n", ncols);

   for( i = 0; i < ncols; ++i )
   {
      rval = mpq_QSchange_objcoef(lpi->prob, ind[i], *SCIPrationalGetGMP(obj[i]));
      QS_CONDRET(rval);
   }
   QS_RETURN(rval);
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiExactScaleRow(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_RATIONAL*        scaleval            /**< scaling multiplier */
   )
{
   register int i;
   int rowlist[1];
   int* rowcnt = NULL, *rowbeg = NULL, *rowind = NULL;
   mpq_t* rowval = NULL, *rhs = NULL, *range = NULL;
   char* sense = NULL;
   int rval = 0;
   mpq_t svl;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   mpq_init(svl);
   lpi->solstat = 0;
   SCIPrationalDebugMessage("scaling row %d with factor %g in QSopt_ex\n", row, scaleval);

   rowlist[0] = row;
   /* get row */
   rval = mpq_QSget_ranged_rows_list(lpi->prob, 1, rowlist, &rowcnt, &rowbeg, &rowind, &rowval, &rhs, &sense, &range, 0);
   QS_TESTG(rval, CLEANUP, " ");

   /* change all coefficients in the constraint */
   for( i = 0; i < rowcnt[0]; ++i )
   {
      mpq_mul(svl, rowval[i], *SCIPrationalGetGMP(scaleval));
      rval = mpq_QSchange_coef(lpi->prob, row, rowind[i], svl);
      QS_TESTG(rval, CLEANUP, " ");
   }

   /* if we have a positive scalar, we just scale rhs and range */
   if( SCIPrationalGetSign(scaleval) >= 0 )
   {
      mpq_mul(svl, *SCIPrationalGetGMP(scaleval), rhs[0]);
      rval = mpq_QSchange_rhscoef(lpi->prob, row, svl);
      QS_TESTG(rval, CLEANUP, " ");
      if (sense[0] == 'R')
      {
         mpq_mul(svl, range[0], *SCIPrationalGetGMP(scaleval));
         rval = mpq_QSchange_range(lpi->prob, row, svl);
         QS_TESTG(rval, CLEANUP, " ");
      }
   }
   /* otherwise, we must change everything */
   else
   {
      switch( sense[0] )
      {
      case 'E':
         mpq_mul(svl, rhs[0], *SCIPrationalGetGMP(scaleval));
         rval = mpq_QSchange_rhscoef(lpi->prob, row, svl);
         QS_TESTG(rval, CLEANUP, " ");
         break;
      case 'L':
         mpq_mul(svl, rhs[0], *SCIPrationalGetGMP(scaleval));
         rval = mpq_QSchange_rhscoef(lpi->prob, row, svl);
         QS_TESTG(rval, CLEANUP, " ");
         rval = mpq_QSchange_sense(lpi->prob, row, 'G');
         QS_TESTG(rval, CLEANUP, " ");
         break;
      case 'G':
         mpq_mul(svl, rhs[0], *SCIPrationalGetGMP(scaleval));
         rval = mpq_QSchange_rhscoef(lpi->prob, row, svl);
         QS_TESTG(rval, CLEANUP, " ");
         rval = mpq_QSchange_sense(lpi->prob, row, 'L');
         QS_TESTG(rval, CLEANUP, " ");
         break;
      case 'R':
         mpq_add(svl, rhs[0], range[0]);
         mpq_mul(svl, svl, *SCIPrationalGetGMP(scaleval));
         rval = mpq_QSchange_rhscoef(lpi->prob, row, svl);
         QS_TESTG(rval, CLEANUP, " ");
         mpq_abs(svl,*SCIPrationalGetGMP(scaleval));
         mpq_mul(svl, svl, range[0]);
         rval = mpq_QSchange_range(lpi->prob, row, svl);
         QS_TESTG(rval, CLEANUP, " ");
         break;
      default:
         SCIPerrorMessage("Imposible! received sense %c (not E L G R)", sense[0]);
         rval = 1;
         goto CLEANUP;
      }
   }
   /* now we must free all received arrays */
   /* ending */
 CLEANUP:
   if (rowcnt) mpq_QSfree(rowcnt);
   if (rowbeg) mpq_QSfree(rowbeg);
   if (rowind) mpq_QSfree(rowind);
   if (rowval) mpq_EGlpNumFreeArray(rowval);
   if (rhs) mpq_EGlpNumFreeArray(rhs);
   if (sense) mpq_QSfree(sense);
   if (range) mpq_EGlpNumFreeArray(range);
   mpq_clear(svl);

   QS_RETURN(rval);
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPlpiExactScaleCol(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_RATIONAL*        scaleval            /**< scaling multiplier */
   )
{
   register int i;
   int collist[1];
   int* colcnt=0, *colbeg=0, *colind=0;
   mpq_t* colval=0, *lb=0, *ub=0, *obj=0;
   int rval = 0;
   mpq_t svl;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   mpq_init(svl);
   lpi->solstat = 0;
   SCIPdebugMessage("scaling column %d with factor %g in QSopt_ex\n",
      col, SCIPrationalGetReal(scaleval));

   /* get the column */
   collist[0] = col;
   rval = mpq_QSget_columns_list(lpi->prob, 1, collist, &colcnt, &colbeg, &colind, &colval, &obj, &lb, &ub, 0);
   QS_TESTG(rval,CLEANUP," ");

   /* scale column coefficients */
   for( i = 0; i < colcnt[0]; ++i )
   {
      mpq_mul(svl, colval[i], *SCIPrationalGetGMP(scaleval));
      rval = mpq_QSchange_coef(lpi->prob, colind[i], col, svl);
      QS_TESTG(rval,CLEANUP," ");
   }

   /* scale objective value */
   mpq_mul(svl, obj[0], *SCIPrationalGetGMP(scaleval));
   rval = mpq_QSchange_objcoef(lpi->prob, col, svl);
   QS_TESTG(rval,CLEANUP," ");

   /* scale column bounds */
   if( mpq_sgn(*SCIPrationalGetGMP(scaleval)) < 0 )
   {
      mpq_set(obj[0], lb[0]);
      mpq_neg(lb[0], ub[0]);
      mpq_neg(ub[0], obj[0]);
   }
   if( mpq_cmp(lb[0],mpq_ILL_MINDOUBLE) > 0 )
   {
      mpq_abs(svl,*SCIPrationalGetGMP(scaleval));
      mpq_mul(lb[0], lb[0], svl);
   }
   if( mpq_cmp(ub[0], mpq_ILL_MAXDOUBLE) < 0 )
   {
      mpq_abs(svl, *SCIPrationalGetGMP(scaleval));
      mpq_mul(ub[0], ub[0], svl);
   }

   if( mpq_cmp(lb[0], mpq_ILL_MINDOUBLE) < 0 )
      mpq_set(lb[0], mpq_ILL_MINDOUBLE);
   if( mpq_cmp(ub[0], mpq_ILL_MAXDOUBLE) > 0 )
      mpq_set(ub[0], mpq_ILL_MAXDOUBLE);

   rval = mpq_QSchange_bound(lpi->prob, col, 'L', lb[0]);
   QS_TESTG(rval, CLEANUP, " ");
   rval = mpq_QSchange_bound(lpi->prob, col, 'U', ub[0]);
   QS_TESTG(rval, CLEANUP, " ");

   /* ending */
 CLEANUP:
   if (colcnt) mpq_QSfree(colcnt);
   if (colbeg) mpq_QSfree(colbeg);
   if (colind) mpq_QSfree(colind);
   if (colval) mpq_EGlpNumFreeArray(colval);
   if (obj) mpq_EGlpNumFreeArray(obj);
   if (lb) mpq_EGlpNumFreeArray(lb);
   if (ub) mpq_EGlpNumFreeArray(ub);
   mpq_clear(svl);

   QS_RETURN(rval);
}
/**@} */

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
   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(nrows != NULL);

   SCIPdebugMessage("getting number of rows\n");

   *nrows = mpq_QSget_rowcount(lpi->prob);

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiExactGetNCols(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(ncols != NULL);

   SCIPdebugMessage("getting number of columns\n");

   *ncols = mpq_QSget_colcount(lpi->prob);

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiExactGetNNonz(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting number of columns\n");

   *nnonz = mpq_QSget_nzcount(lpi->prob);

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
   SCIP_RATIONAL**       lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_RATIONAL**       ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_RATIONAL**       val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int len;
   register int i;
   mpq_t* lval = NULL;
   mpq_t* llb = NULL;
   mpq_t* lub = NULL;
   int rval = 0;
   int* lcnt = NULL;
   int* lbeg = NULL;
   int* lind = NULL;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert( 0 <= firstcol && firstcol <= lastcol && lastcol < mpq_QSget_colcount (lpi->prob) );
   assert( (lb == 0 && ub == 0) || (lb != 0 && ub != 0));
   assert( (nnonz != 0 && beg != 0 && ind != 0 && val != 0) || (nnonz == 0 && beg == 0 && ind == 0 && val == 0) );

   SCIPdebugMessage("getting columns %d to %d\n", firstcol, lastcol);

   /* build col-list */
   len = lastcol - firstcol + 1;
   SCIP_CALL( ensureColMem(lpi,len) );
   for( i = 0; i < len; ++i )
      lpi->iccnt[i] = i + firstcol;

   /* get data from qsopt */
   rval = mpq_QSget_columns_list(lpi->prob, len, lpi->iccnt, nnonz ? (&lcnt) : 0, nnonz ? (&lbeg) : 0, nnonz ? (&lind) : 0,
      nnonz ? (&lval) : 0, 0, lb ? (&llb) : 0, lb ? (&lub) : 0, 0);

   QS_TESTG(rval, CLEANUP, " ");

   /* store in the user-provided data */
   if( nnonz )
   {
      assert( lbeg != NULL );
      assert( lcnt != NULL );
      assert( lind != NULL );
      assert( lval != NULL );

      *nnonz = lbeg[len-1] + lcnt[len-1];
      for( i = 0 ; i < len ; i++ )
         beg[i] = lbeg[i];  /*lint !e613*/
      for( i = 0; i < *nnonz; ++i )
      {
         ind[i] = lind[i];  /*lint !e613*/
         SCIPrationalSetGMP(val[i], lval[i]);
      }
   }
   if( lb )
   {
      assert(llb != NULL);
      assert(lub != NULL);

      for( i = 0; i < len; ++i )
      {
         SCIPrationalSetGMP(lb[i], llb[i]);
         SCIPrationalSetGMP(ub[i], lub[i]);   /*lint !e613*/
      }
   }

 CLEANUP:
   if (lval) mpq_EGlpNumFreeArray(lval);
   if (lub) mpq_EGlpNumFreeArray(lub);
   if (llb) mpq_EGlpNumFreeArray(llb);
   if (lind) mpq_QSfree(lind);
   if (lbeg) mpq_QSfree(lbeg);
   if (lcnt) mpq_QSfree(lcnt);

   QS_RETURN(rval);
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiExactGetRows(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_RATIONAL**       lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_RATIONAL**       rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_RATIONAL**       val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   const int len = lastrow - firstrow + 1;
   register int i;
   mpq_t* lval = NULL;
   mpq_t* lrhs = NULL;
   mpq_t* lrng = NULL;
   int rval = 0;
   int* lcnt = NULL;
   int* lbeg = NULL;
   int* lind = NULL;
   char* lsense = NULL;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < mpq_QSget_rowcount (lpi->prob));
   assert( (lhs == 0 && rhs == 0) || (rhs != 0 && lhs != 0));
   assert( (nnonz != 0 && beg != 0 && ind != 0 && val != 0) || (nnonz == 0 && beg == 0 && ind == 0 && val == 0));

   SCIPdebugMessage("getting rows %d to %d\n", firstrow, lastrow);

   /* build row-list */
   SCIP_CALL( ensureRowMem(lpi, len) );
   for( i = 0; i < len; ++i )
      lpi->ircnt[i] = i + firstrow;

   /* get data from qsopt */
   rval = mpq_QSget_ranged_rows_list(lpi->prob, len, lpi->ircnt, nnonz ? (&lcnt) : 0, nnonz ? (&lbeg) : 0, nnonz ? (&lind) : 0,
      nnonz ? (&lval) : 0, rhs ? (&lrhs) : 0, rhs ? (&lsense) : 0, rhs ? (&lrng) : 0, 0);
   QS_TESTG(rval, CLEANUP, " ");

   /* store in the user-provided data */
   if( nnonz )
   {
      assert(lbeg != NULL);
      assert(lcnt != NULL);
      assert(lind != NULL);
      assert(lval != NULL);

      *nnonz = lbeg[len-1] + lcnt[len-1];
      for( i = 0; i < len; i++ )
         beg[i] = lbeg[i];  /*lint !e613*/
      for( i = 0; i < *nnonz; ++i )
      {
         ind[i] = lind[i];  /*lint !e613*/
         SCIPrationalSetGMP(val[i], lval[i]);  /*lint !e613*/
      }
   }
   if( rhs )
   {
      assert(lrhs != NULL);
      assert(lrng != NULL);
      assert(lsense != NULL);

      for( i = 0; i < len; ++i )
      {
         switch( lsense[i] )
         {
         case 'R':
            SCIPrationalSetGMP(lhs[i], lrhs[i]);
            SCIPrationalSetGMP(rhs[i], lrng[i]);
            SCIPrationalAdd(rhs[i], rhs[i], lhs[i]);
            break;
         case 'E':
            SCIPrationalSetGMP(lhs[i], lrhs[i]);
            SCIPrationalSetGMP(rhs[i], lrhs[i]);
            break;
         case 'L':
            SCIPrationalSetGMP(rhs[i], lrhs[i]);
            SCIPrationalSetGMP(lhs[i], mpq_ILL_MINDOUBLE);      /*lint !e613*/
            break;
         case 'G':
            SCIPrationalSetGMP(lhs[i], lrhs[i]);
            SCIPrationalSetGMP(rhs[i], mpq_ILL_MAXDOUBLE);
            break;
         default:
            SCIPerrorMessage("Unknown sense %c from QSopt_ex", lsense[i]);
            SCIPABORT();
         }
      }
   }

 CLEANUP:
   if (lsense) mpq_QSfree(lsense);
   if (lrng) mpq_EGlpNumFreeArray(lrng);
   if (lrhs) mpq_EGlpNumFreeArray(lrhs);
   if (lval) mpq_EGlpNumFreeArray(lval);
   if (lind) mpq_QSfree(lind);
   if (lbeg) mpq_QSfree(lbeg);
   if (lcnt) mpq_QSfree(lcnt);

   QS_RETURN(rval);
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiExactGetObj(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_RATIONAL**       vals                /**< array to store objective coefficients */
   )
{
   const int len = lastcol - firstcol + 1;
   int rval = 0;
   register int i;
   mpq_t* valgmp;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < mpq_QSget_colcount (lpi->prob));

   SCIPdebugMessage("getting objective values %d to %d\n", firstcol, lastcol);

   /* build col-list */
   SCIP_CALL( ensureColMem(lpi,len) );
   for( i = 0; i < len; ++i )
      lpi->iccnt[i] = i + firstcol;

   SCIP_ALLOC( BMSallocMemoryArray(&valgmp, len) );
   SCIPrationalSetGMPArray(valgmp, vals, len);
   /* get data from qsopt */
   rval = mpq_QSget_obj_list(lpi->prob, len, lpi->iccnt, valgmp);
   SCIPrationalSetArrayGMP(vals, valgmp, len);

   SCIPrationalClearArrayGMP(valgmp, len);
   BMSfreeMemoryArray(&valgmp);

   QS_RETURN(rval);
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiExactGetBounds(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective value for */
   int                   lastcol,            /**< last column to get objective value for */
   SCIP_RATIONAL**       lbs,                /**< array to store lower bound values, or NULL */
   SCIP_RATIONAL**       ubs                 /**< array to store upper bound values, or NULL */
   )
{
   const int len = lastcol - firstcol + 1;
   int rval = 0;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(0 <= firstcol && firstcol <= lastcol&& lastcol < mpq_QSget_colcount (lpi->prob));

   SCIPdebugMessage("getting bound values %d to %d\n", firstcol, lastcol);

   /* build col-list */
   SCIP_CALL( ensureColMem(lpi,len) );
   for( i = 0; i < len; ++i )
   {
      if( lbs != NULL )
      {
         QS_CONDRET( mpq_QSget_bound(lpi->prob, i + firstcol, 'L', SCIPrationalGetGMP(lbs[i])) );
         SCIPrationalCheckInfByValue(lbs[i]);
         SCIPrationalResetFloatingPointRepresentable(lbs[i]);
      }
      if( ubs != NULL )
      {
         QS_CONDRET( mpq_QSget_bound(lpi->prob, i + firstcol, 'U', SCIPrationalGetGMP(ubs[i])) );
         SCIPrationalCheckInfByValue(ubs[i]);
         SCIPrationalResetFloatingPointRepresentable(ubs[i]);
      }
   }

   QS_RETURN(rval);
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiExactGetSides(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_RATIONAL**       lhss,               /**< array to store left hand side values, or NULL */
   SCIP_RATIONAL**       rhss                /**< array to store right hand side values, or NULL */
   )
{
   const int len = lastrow - firstrow + 1;
   register int i;
   mpq_t* lrhs = 0;
   mpq_t* lrng = 0;
   int rval = 0;
   char* lsense=0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < mpq_QSget_rowcount (lpi->prob));
   assert(rhss != 0 && lhss != 0);

   SCIPdebugMessage("getting row sides %d to %d\n", firstrow, lastrow);

   /* build row-list */
   SCIP_CALL( ensureRowMem(lpi, len) );
   for( i = 0; i < len; ++i )
      lpi->ircnt[i] = i + firstrow;

   /* get data from qsopt */
   rval = mpq_QSget_ranged_rows_list(lpi->prob, len, lpi->ircnt, 0, 0, 0, 0, &lrhs, &lsense, &lrng, 0);
   QS_TESTG(rval, CLEANUP, " ");

   /* store in the user-provided data */
   for( i = 0; i < len; ++i )
   {
      switch (lsense[i])
      {
      case 'R':
         SCIPrationalSetGMP(lhss[i], lrhs[i]);
         SCIPrationalSetGMP(rhss[i], lrng[i]);
         SCIPrationalAdd(rhss[i], rhss[i], lhss[i]);
         break;
      case 'E':
         SCIPrationalSetGMP(lhss[i], lrhs[i]);
         SCIPrationalSetGMP(rhss[i], lrhs[i]);
         break;
      case 'L':
         SCIPrationalSetGMP(rhss[i], lrhs[i]);
         SCIPrationalSetGMP(lhss[i], mpq_ILL_MINDOUBLE);
         break;
      case 'G':
         SCIPrationalSetGMP(lhss[i], lrhs[i]);
         SCIPrationalSetGMP(rhss[i], mpq_ILL_MAXDOUBLE);
         break;
      default:
         SCIPerrorMessage("Unknown sense %c from QSopt_ex", lsense[i]);
         SCIPABORT();
      }
   }

 CLEANUP:
   if (lsense)
      mpq_QSfree(lsense);
   if (lrng)
      mpq_EGlpNumFreeArray(lrng);
   if (lrhs)
      mpq_EGlpNumFreeArray(lrhs);

   QS_RETURN(rval);
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiExactGetCoef(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_RATIONAL*        val                 /**< pointer to store the value of the coefficient */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting coefficient of row %d col %d\n", row, col);

   rval = mpq_QSget_coef(lpi->prob, row, col, SCIPrationalGetGMP(val));
   SCIPrationalCheckInfByValue(val);
   SCIPrationalResetFloatingPointRepresentable(val);

   QS_RETURN(rval);
}

/**@} */

/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiExactSolvePrimal(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   int rval = 0;
   QSbasis* B;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("calling QSopt_ex primal simplex: %d cols, %d rows, %d nz\n", mpq_QSget_colcount(lpi->prob),
      mpq_QSget_rowcount(lpi->prob), mpq_QSget_nzcount(lpi->prob));

   B = mpq_QSget_basis(lpi->prob);
   rval = QSexact_solver(lpi->prob, 0, 0, B, PRIMAL_SIMPLEX, &(lpi->solstat));
   if( B )
      mpq_QSfree_basis(B);

   QS_RETURN(rval);
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiExactSolveDual(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   int rval = 0;
   QSbasis* B;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("calling QSopt_ex dual simplex: %d cols, %d rows, %d nz\n", mpq_QSget_colcount(lpi->prob),
      mpq_QSget_rowcount(lpi->prob), mpq_QSget_nzcount(lpi->prob));

   B = mpq_QSget_basis(lpi->prob);
   rval = QSexact_solver(lpi->prob, 0, 0, B, DUAL_SIMPLEX, &(lpi->solstat));
   if( B )
      mpq_QSfree_basis(B);

   QS_RETURN(rval);
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiExactSolveBarrier(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{  /*lint --e{715}*/
   return SCIPlpiExactSolveDual(lpi);
}

/**@} */

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
   assert(lpi!=NULL);
   assert(lpi->prob!=NULL);

   return (lpi->solstat != 0 && lpi->solstat != QS_LP_MODIFIED && lpi->solstat != QS_LP_CHANGE_PREC);
}

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_RETCODE SCIPlpiExactGetSolFeasibility(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting solution feasibility\n");

   *primalfeasible = FALSE;
   *dualfeasible = FALSE;

   if( lpi->solstat == QS_LP_OPTIMAL || lpi->solstat == QS_LP_UNBOUNDED )
      *primalfeasible = TRUE;

   /* @todo: check why we can conclude dual feasibility from primal infeasibility. in theory, the LP could be primal and
    * dual infeasible as well; see also SCIPlpiExactIsDualFeasible() and SCIPlpiExactIsDualInfeasible()
    */
#ifdef USEOBJLIM
   if( lpi->solstat == QS_LP_OPTIMAL || lpi->solstat == QS_LP_INFEASIBLE || lpi->solstat == QS_LP_OBJ_LIMIT )
#else
   if( lpi->solstat == QS_LP_OPTIMAL || lpi->solstat == QS_LP_INFEASIBLE )
#endif
      *dualfeasible = TRUE;

   return SCIP_OKAY;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExactExistsPrimalRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking primal ray existance\n");

   return (lpi->solstat == QS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExactHasPrimalRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for primal ray\n");

   /* the current version of QSopt_ex can not give a primal certificate of unboundness */
   return FALSE;
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiExactIsPrimalUnbounded(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for primal unboundness\n");

   return (lpi->solstat == QS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiExactIsPrimalInfeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for primal infeasibility\n");

   return (lpi->solstat == QS_LP_INFEASIBLE);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiExactIsPrimalFeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for primal feasibility\n");

   return (lpi->solstat == QS_LP_OPTIMAL || lpi->solstat == QS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExactExistsDualRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual ray availability\n");

   return (lpi->solstat == QS_LP_INFEASIBLE);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExactHasDualRay(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual ray availability\n");

   return (lpi->solstat == QS_LP_INFEASIBLE);
}

/** returns TRUE iff LP is dual unbounded */
SCIP_Bool SCIPlpiExactIsDualUnbounded(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual unboundness\n");

   return (lpi->solstat == QS_LP_INFEASIBLE);
}

/** returns TRUE iff LP is dual infeasible */
SCIP_Bool SCIPlpiExactIsDualInfeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual infeasibility\n");

   return (lpi->solstat == QS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiExactIsDualFeasible(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual feasibility\n");

#ifdef USEOBJLIM
   return (lpi->solstat == QS_LP_OPTIMAL || lpi->solstat == QS_LP_OBJ_LIMIT );
#else
   return (lpi->solstat == QS_LP_OPTIMAL);
#endif
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiExactIsOptimal(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for optimality\n");

   return (lpi->solstat == QS_LP_OPTIMAL);
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiExactIsStable(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for numerical stability\n");

   return (lpi->solstat != QS_LP_NUMERR);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiExactIsObjlimExc(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for objective limit exceeded\n");

#ifdef USEOBJLIM
   return (lpi->solstat == QS_LP_OBJ_LIMIT);
#else
   return FALSE;
#endif
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiExactIsIterlimExc(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for iteration limit exceeded\n");

   return (lpi->solstat == QS_LP_ITER_LIMIT);
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiExactIsTimelimExc(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for time limit exceeded\n");

   return (lpi->solstat == QS_LP_TIME_LIMIT);
}

/** returns the internal solution status of the solver */
int SCIPlpiExactGetInternalStatus(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting internal solution status\n");

   return lpi->solstat;
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiExactIgnoreInstability(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("ignore instability (will fail)\n");

   /* it seems that in QSopt_ex this does not make much sense */
   *success = FALSE;

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiExactGetObjval(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        objval              /**< stores the objective value */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting solution's objective value\n");

   rval = mpq_QSget_objval(lpi->prob, SCIPrationalGetGMP(objval));
   SCIPrationalCheckInfByValue(objval);
   SCIPrationalResetFloatingPointRepresentable(objval);

   QS_RETURN(rval);
}

/** gets primal and dual solution vectors */
SCIP_RETCODE SCIPlpiExactGetSol(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_RATIONAL**       primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_RATIONAL**       dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_RATIONAL**       activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_RATIONAL**       redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{
   int rval = 0, nrows, ncols;
   register int i;
   mpq_t* primsolgmp, *dualsolgmp, *redcostgmp, *objvalgmp;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting solution\n");

   nrows = mpq_QSget_rowcount(lpi->prob);
   ncols = mpq_QSget_colcount(lpi->prob);
   SCIP_CALL( ensureRowMem(lpi, nrows) );

   if( primsol == NULL )
      primsolgmp = NULL;
   else
   {
      SCIP_ALLOC( BMSallocMemoryArray(&primsolgmp, ncols) );
      SCIPrationalSetGMPArray(primsolgmp, primsol, ncols);
   }
   if( redcost == NULL )
      redcostgmp = NULL;
   else
   {
      SCIP_ALLOC( BMSallocMemoryArray(&redcostgmp, ncols) );
      SCIPrationalSetGMPArray(redcostgmp, redcost, ncols);
   }
   if( dualsol == NULL )
      dualsolgmp = NULL;
   else
   {
      SCIP_ALLOC( BMSallocMemoryArray(&dualsolgmp, nrows) );
      SCIPrationalSetGMPArray(dualsolgmp, dualsol, nrows);
   }
   if( objval != NULL )
   {
      objvalgmp = SCIPrationalGetGMP(objval);
      SCIPrationalResetFloatingPointRepresentable(objval);
   }
   else
      objvalgmp = NULL;

   rval = mpq_QSget_solution(lpi->prob, objvalgmp, primsolgmp, dualsolgmp, lpi->irng, redcostgmp);

   if( objval != NULL )
      SCIPrationalCheckInfByValue(objval);

   if( redcost != NULL )
   {
      SCIPrationalSetArrayGMP(redcost, redcostgmp, ncols);
      SCIPrationalClearArrayGMP(redcostgmp, ncols);
      BMSfreeMemoryArray(&redcostgmp);
   }
   if( primsol != NULL )
   {
      SCIPrationalSetArrayGMP(primsol, primsolgmp, ncols);
      SCIPrationalClearArrayGMP(primsolgmp, ncols);
      BMSfreeMemoryArray(&primsolgmp);
   }
   if( dualsol != NULL )
   {
      SCIPrationalSetArrayGMP(dualsol, dualsolgmp, nrows);
      SCIPrationalClearArrayGMP(dualsolgmp, nrows);
      BMSfreeMemoryArray(&dualsolgmp);
   }

   QS_CONDRET(rval);

   rval = mpq_QSget_rhs(lpi->prob, lpi->irhs);
   QS_CONDRET(rval);
   rval = mpq_QSget_senses(lpi->prob, lpi->isen);
   QS_CONDRET(rval);

   /* build back the activity */
   if( activity != NULL )
   {
      for( i = 0; i < nrows; ++i )
      {
         switch (lpi->isen[i])
         {
         case 'R':
         case 'E':
         case 'G':
            mpq_add(*SCIPrationalGetGMP(activity[i]), lpi->irhs[i], lpi->irng[i]);
            SCIPrationalResetFloatingPointRepresentable(activity[i]);
            SCIPrationalCheckInfByValue(activity[i]);
            break;
         case 'L':
            mpq_sub(*SCIPrationalGetGMP(activity[i]), lpi->irhs[i], lpi->irng[i]);
            SCIPrationalResetFloatingPointRepresentable(activity[i]);
            SCIPrationalCheckInfByValue(activity[i]);
            break;
         default:
            SCIPerrorMessage("unknown sense %c\n", lpi->isen[i]);
            SCIPABORT();
         }
      }
   }

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiExactGetPrimalRay(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL**       ray                 /**< primal ray */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPerrorMessage("SCIPlpiExactGetPrimalRay() not supported by QSopt_ex.\n");

   return SCIP_ERROR;
}

/** gets dual farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiExactGetDualfarkas(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL**       dualfarkas          /**< dual farkas row multipliers */
   )
{
   int rval = 0;
   int nrows;
   mpq_t* dualfarkasgmp;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(dualfarkas != NULL);

   SCIPdebugMessage("calling QSopt_ex dual farkas: %d cols, %d rows, %d non zeros\n", mpq_QSget_colcount (lpi->prob),
      mpq_QSget_rowcount(lpi->prob), mpq_QSget_nzcount(lpi->prob));

   nrows = mpq_QSget_rowcount(lpi->prob);\
   SCIP_ALLOC( BMSallocMemoryArray(&dualfarkasgmp, nrows) );
   SCIPrationalSetGMPArray(dualfarkasgmp, dualfarkas, nrows);

   rval = mpq_QSget_infeas_array(lpi->prob, dualfarkasgmp);

   SCIPrationalSetArrayGMP(dualfarkas, dualfarkasgmp, nrows);
   SCIPrationalClearArrayGMP(dualfarkasgmp, nrows);
   BMSfreeMemoryArray(&dualfarkasgmp);

   QS_RETURN(rval);
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiExactGetIterations(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   int rval = 0, nit;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   rval = mpq_QSget_itcnt(lpi->prob, 0, 0, 0, 0, &nit);
   QS_CONDRET(rval);

   *iterations = nit - lpi->previt;
   lpi->previt = nit;

   QS_RETURN(rval);
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
   int rval = 0, ncols, nrows;
   char* icstat = NULL;
   char* irstat = NULL;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("saving QSopt_ex basis into %p/%p\n", (void *) cstat, (void *) rstat);

   ncols = mpq_QSget_colcount(lpi->prob);
   nrows = mpq_QSget_rowcount(lpi->prob);

   SCIP_CALL( ensureTabMem(lpi, nrows + ncols) );

   icstat = lpi->ibas;
   irstat = lpi->ibas+ncols;
   rval = mpq_QSget_basis_array(lpi->prob, icstat, irstat);
   QS_CONDRET(rval);

   /* now we must transform QSopt_ex codes into SCIP codes */
   for( i = 0; i < nrows; ++i )
   {
      switch (irstat[i])
      {
      case QS_ROW_BSTAT_LOWER:
         rstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
         break;
      case QS_ROW_BSTAT_BASIC:
         rstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/
         break;
      case QS_ROW_BSTAT_UPPER:
         rstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
         break;
      default:
         SCIPerrorMessage("Unknown row basic status %c", rstat[i]);
         SCIPABORT();
      }
   }
   for( i = 0; i < ncols; ++i )
   {
      switch(icstat[i])
      {
      case QS_COL_BSTAT_LOWER:
         cstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
         break;
      case QS_COL_BSTAT_BASIC:
         cstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/
         break;
      case QS_COL_BSTAT_UPPER:
         cstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
         break;
      case QS_COL_BSTAT_FREE:
         cstat[i] = SCIP_BASESTAT_ZERO; /*lint !e641*/
         break;
      default:
         SCIPerrorMessage("Unknown column basic status %c", cstat[i]);
         SCIPABORT();
      }
   }
   QS_RETURN(rval);
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiExactSetBase(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   )
{
   int rval = 0, ncols, nrows;
   register int i;
   char* icstat=0, *irstat = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("loading basis %p/%p into QSopt_ex\n", (void *) cstat, (void *) rstat);

   ncols = mpq_QSget_colcount(lpi->prob);
   nrows = mpq_QSget_rowcount(lpi->prob);

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureColMem(lpi, ncols) );
   SCIP_CALL( ensureRowMem(lpi, nrows) );
   SCIP_CALL( ensureTabMem(lpi, nrows+ncols) );

   icstat = lpi->ibas;
   irstat = lpi->ibas + ncols;

   /* now we must transform QSopt_ex codes into SCIP codes */
   for( i = 0; i < nrows; ++i )
   {
      switch(rstat[i])
      {
      case SCIP_BASESTAT_LOWER:
         irstat[i] = QS_ROW_BSTAT_LOWER; /*lint !e641*/
         break;
      case SCIP_BASESTAT_BASIC:
         irstat[i] = QS_ROW_BSTAT_BASIC; /*lint !e641*/
         break;
      case SCIP_BASESTAT_UPPER:
         /* sense of inexact LP row is R (ranged row) since this is the only case where the basis status of the
          * slack variable is allowed to be UPPER
          */
         if( lpi->isen[i] == 'R' )
            /* sense of LPEX row is R, too */
            irstat[i] = QS_ROW_BSTAT_UPPER; /*lint !e641*/
         else
            /* sense of LPEX row is L, G or E, thus, basis status must be LOWER/BASIC. we use non-basic status LOWER
             * instead of non-basic status UPPER for slack variable in LPEX. this might happen when the inexact LP
             * is an FP relaxation of the exact LP
             */
            irstat[i] = QS_ROW_BSTAT_LOWER;
         break;
      default:
         SCIPerrorMessage("Unknown row basic status %d", rstat[i]);
         SCIPABORT();
      }
   }
   for( i = 0; i < ncols; ++i )
   {
      switch(cstat[i])
      {
      case SCIP_BASESTAT_LOWER:
         icstat[i] = QS_COL_BSTAT_LOWER; /*lint !e641*/
         break;
      case SCIP_BASESTAT_BASIC:
         icstat[i] = QS_COL_BSTAT_BASIC; /*lint !e641*/
         break;
      case SCIP_BASESTAT_UPPER:
         icstat[i] = QS_COL_BSTAT_UPPER; /*lint !e641*/
         break;
      case SCIP_BASESTAT_ZERO:
         icstat[i] = QS_COL_BSTAT_FREE; /*lint !e641*/
         break;
      default:
         SCIPerrorMessage("Unknown column basic status %d", cstat[i]);
         SCIPABORT();
      }
   }

   /* set the basis */
   rval = mpq_QSload_basis_array(lpi->prob, icstat, irstat);
   QS_RETURN(rval);
}

/** returns the indices of the basic columns and rows */
SCIP_RETCODE SCIPlpiExactGetBasisInd(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   int rval = 0, nrows, ncols;
   register int i;

   assert(lpi!=NULL);
   assert(lpi->prob!=NULL);

   SCIPdebugMessage("getting basis information\n");

   nrows = mpq_QSget_rowcount(lpi->prob);
   ncols = mpq_QSget_colcount(lpi->prob);
   rval = mpq_QSget_basis_order(lpi->prob, bind);
   QS_CONDRET(rval);

   /* transform QSopt_ex basis header into SCIP format */
   for( i = 0; i < nrows; ++i )
   {
      if( bind[i] >= ncols )
         bind[i] = -(bind[i] - ncols - 1);
   }

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

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(blkmem != NULL);
   assert(lpistate != NULL);

   ncols = mpq_QSget_colcount(lpi->prob);
   nrows = mpq_QSget_rowcount(lpi->prob);

   assert(ncols >= 0);
   assert(nrows >= 0);

   /* allocate lpistate data */
   SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows));
   SCIPdebugMessage("storing QSopt_ex LPI state in %p (%d cols, %d rows)\n", (void *) *lpistate, ncols, nrows);

   /* get unpacked basis information from QSopt_ex */
   SCIP_CALL( ensureColMem(lpi, ncols) );
   SCIP_CALL( ensureRowMem(lpi, nrows) );
   SCIP_CALL( SCIPlpiExactGetBase(lpi, lpi->iccnt, lpi->ircnt) );

   /* pack LPi state data */
   (*lpistate)->ncols = ncols;
   (*lpistate)->nrows = nrows;
   lpistatePack(*lpistate, lpi->iccnt, lpi->ircnt);

   return SCIP_OKAY;
}

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiExactGetState()
 */
SCIP_RETCODE SCIPlpiExactSetState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   )
{  /*lint --e{715} */
   register int i;
   int rval = 0;
   int ncols;
   int nrows;
   char* icstat = 0;
   char* irstat = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   /* if there was no basis information available, LPI state was not stored */
   if (lpistate == NULL)
      QS_RETURN(rval);

   /* continue test */
   ncols = mpq_QSget_colcount(lpi->prob);
   nrows = mpq_QSget_rowcount(lpi->prob);

   assert(ncols >= 0);
   assert(nrows >= 0);
   assert(lpistate->ncols <= ncols);
   assert(lpistate->nrows <= nrows);

   SCIPdebugMessage("loading LPI state %p (%d cols, %d rows) into QSopt_ex LP (%d cols and %d rows)\n", (void*) lpistate,
      lpistate->ncols, lpistate->nrows, ncols, nrows);

   if( lpistate->ncols == 0 || lpistate->nrows == 0 )
      QS_RETURN(rval);

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureColMem(lpi, ncols) );
   SCIP_CALL( ensureRowMem(lpi, nrows) );
   SCIP_CALL( ensureTabMem(lpi, nrows+ncols) );

   icstat = lpi->ibas;
   irstat = lpi->ibas + ncols;

   /* unpack LPi state data */
   lpistateUnpack(lpistate, lpi->iccnt, lpi->ircnt);

   /* extend the basis to the current LP */
   for( i = lpistate->ncols; i < ncols; ++i )
      lpi->iccnt[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/ /**@todo this has to be corrected for lb = -infinity */
   for( i = lpistate->nrows; i < nrows; ++i )
      lpi->ircnt[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/

   /* convert the loaded basis into QSopt_ex format */
   SCIPdebugMessage("basis status of SCIP lpistate rows (nrows=%d):\n", lpistate->nrows);
   for( i = 0; i < nrows; ++i )
   {
      SCIPdebugMessage("row_%d: %d (%s)\n", i, lpi->ircnt[i],
         lpi->ircnt[i] == SCIP_BASESTAT_LOWER ? "lower" : lpi->ircnt[i] == SCIP_BASESTAT_BASIC ? "basic" : "upper");

      switch( lpi->ircnt[i] )
      {
      case SCIP_BASESTAT_LOWER:
         irstat[i] = QS_ROW_BSTAT_LOWER;
         break;
      case SCIP_BASESTAT_BASIC:
         irstat[i] = QS_ROW_BSTAT_BASIC;
         break;
      case SCIP_BASESTAT_UPPER:
         /* sense of inexact LP row is R (ranged row) since this is the only case where the basis status of the
          * slack variable is allowed to be UPPER
          */
         if( lpi->isen[i] == 'R' )
            /* sense of LPEX row is R, too */
            irstat[i] = QS_ROW_BSTAT_UPPER;
         else
            /* sense of LPEX row is L, G or E, thus, basis status must be LOWER/BASIC. we use non-basic status LOWER
             * instead of non-basic status UPPER for slack variable in LPEX. this might happen when the inexact LP
             * is an FP relaxation of the exact LP
             */
            irstat[i] = QS_ROW_BSTAT_LOWER;
         break;
      default:
         SCIPerrorMessage("Unknown row basic status %d", lpi->ircnt[i]);
         SCIPABORT();
         break;
      }
   }
   for( i = 0; i < ncols; ++i )
   {
      switch(lpi->iccnt[i])
      {
      case SCIP_BASESTAT_LOWER:
         icstat[i] = QS_COL_BSTAT_LOWER;
         break;
      case SCIP_BASESTAT_BASIC:
         icstat[i] = QS_COL_BSTAT_BASIC;
         break;
      case SCIP_BASESTAT_UPPER:
         icstat[i] = QS_COL_BSTAT_UPPER;
         break;
      case SCIP_BASESTAT_ZERO:
         icstat[i] = QS_COL_BSTAT_FREE;
         break;
      default:
         SCIPerrorMessage("Unknown column basic status %d", lpi->iccnt[i]);
         SCIPABORT();
         break;
      }
   }

   /* set the basis */
   rval = mpq_QSload_basis_array(lpi->prob, icstat, irstat);
   QS_RETURN(rval);
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiExactFreeState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
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
SCIP_Bool SCIPlpiExactHasStateBasis(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{  /*lint --e{715} */
   return (lpistate != NULL);
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiExactReadState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("reading QSopt_ex LP state from file <%s>\n", fname);

   rval = mpq_QSread_and_load_basis(lpi->prob, fname);
   if( rval )
   {
      SCIPerrorMessage("Error while loading basis from file <%s>.\n", fname);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** writes LP state (like basis information) to a file */
SCIP_RETCODE SCIPlpiExactWriteState(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   QSbasis* bas = 0;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("writing QSopt_ex LP state to file <%s>\n", fname);

   bas = mpq_QSget_basis(lpi->prob);
   QS_ERROR(bas == 0, "Could not get basis from problem.");   /*lint !e820*/
   assert( bas );

   rval = mpq_QSwrite_basis(lpi->prob, bas, fname);
   mpq_QSfree(bas);
   if( rval )
   {
      SCIPerrorMessage("Could not write basis to file <%s>.\n", fname);
      return SCIP_WRITEERROR;
   }

   return SCIP_OKAY;
}

#ifdef SCIP_DISABLED_CODE
/** checks whether LPi state (i.e. basis information) is dual feasible and returns corresponding dual objective value.
 *  if wanted it will first directly test the corresponding approximate dual and primal solution
 *  (corrected via dual variables for bounds and primal variables for slacks if possible) for optimality
 *  before performing the dual feasibility test on the more expensive exact basic solution.
 */
SCIP_RETCODE SCIPlpiExactStateDualFeasible(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate,           /**< LPi state information (like basis information) */
   SCIP_Bool             useprestep,         /**< should approximate primal and dual solution first */
   SCIP_Real*            primalsol,          /**< approximate primal solution; or NULL to compute by exact LP solver */
   SCIP_Real*            dualsol,            /**< approximate dual solution; or NULL to compute by exact LP solver */
   SCIP_Bool*            result,             /**< pointer to store whether given LPi state is dual feasible */
   mpq_t*                dualobjval          /**< pointer to store dual objective value in case of dual feasibility */
   )
{  /*lint --e{715} */
   int rval = 0;
   QSbasis* B;

   *result = FALSE;

   /* loads LPi state (like basis information) into solver */
   SCIP_CALL( SCIPlpiExactSetState(lpi, blkmem, lpistate) );

   /* checks whether basis just loaded into the solver is dual feasible */
   B =  mpq_QSget_basis(lpi->prob);

#ifdef VERIFY_OUT
   rval = QSexact_verify(lpi->prob, B, (int) useprestep, primalsol, dualsol, (char*) result, dualobjval, 0);
#else
   rval = QSexact_verify(lpi->prob, B, (int) useprestep, primalsol, dualsol, (char*) result, dualobjval, 1);
#endif

   if( B )
      mpq_QSfree_basis(B);

   QS_RETURN(rval);
}
#endif

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
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(ival != NULL);

   SCIPdebugMessage("getting int parameter %d\n", type);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
   case SCIP_LPPAR_FASTMIP:
   case SCIP_LPPAR_PRESOLVING:
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_SCALING:
      rval = mpq_QSget_param(lpi->prob, QS_PARAM_SIMPLEX_SCALING, ival);
      if( *ival )
         *ival = TRUE;
      else
         *ival = FALSE;
      break;
   case SCIP_LPPAR_PRICING:
      *ival = lpi->pricing;
      break;
   case SCIP_LPPAR_LPINFO:
      rval = mpq_QSget_param(lpi->prob, QS_PARAM_SIMPLEX_DISPLAY, ival);
      if( *ival )
         *ival = TRUE;
      else
         *ival = FALSE;
      break;
   case SCIP_LPPAR_LPITLIM:
      rval = mpq_QSget_param(lpi->prob, QS_PARAM_SIMPLEX_MAX_ITERATIONS, ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   QS_RETURN(rval);
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiExactSetIntpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("setting int parameter %d to %d\n", type, ival);

   switch (type)
   {
   case SCIP_LPPAR_SCALING:
      if( ival == TRUE )
         rval = mpq_QSset_param(lpi->prob, QS_PARAM_SIMPLEX_SCALING, 1);
      else
         rval = mpq_QSset_param(lpi->prob, QS_PARAM_SIMPLEX_SCALING, 0);
      break;
   case SCIP_LPPAR_PRICING:
      lpi->pricing = ival;
      switch( ival )
      {
      case SCIP_PRICING_AUTO:
      case SCIP_PRICING_LPIDEFAULT:
      case SCIP_PRICING_FULL:
      case SCIP_PRICING_STEEP:
      case SCIP_PRICING_STEEPQSTART:
         rval = mpq_QSset_param(lpi->prob, QS_PARAM_PRIMAL_PRICING, QS_PRICE_PSTEEP);
         rval += mpq_QSset_param(lpi->prob, QS_PARAM_DUAL_PRICING, QS_PRICE_DSTEEP);
         break;
      case SCIP_PRICING_PARTIAL:
         rval = mpq_QSset_param(lpi->prob,QS_PARAM_PRIMAL_PRICING,QS_PRICE_PMULTPARTIAL);
         rval += mpq_QSset_param(lpi->prob,QS_PARAM_DUAL_PRICING,QS_PRICE_DMULTPARTIAL);
         break;
      case SCIP_PRICING_DEVEX:
         rval = mpq_QSset_param(lpi->prob,QS_PARAM_PRIMAL_PRICING,QS_PRICE_PDEVEX);
         rval += mpq_QSset_param(lpi->prob,QS_PARAM_DUAL_PRICING,QS_PRICE_DDEVEX);
         break;
      default:
         return SCIP_LPERROR;
      }
      break;
   case SCIP_LPPAR_LPINFO:
      if( ival == TRUE )
         rval = mpq_QSset_param(lpi->prob, QS_PARAM_SIMPLEX_DISPLAY, 1);
      else
         rval = mpq_QSset_param(lpi->prob, QS_PARAM_SIMPLEX_DISPLAY, 0);
      break;
   case SCIP_LPPAR_LPITLIM:
      rval = mpq_QSset_param(lpi->prob, QS_PARAM_SIMPLEX_MAX_ITERATIONS, ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   QS_RETURN(rval);
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiExactGetRealpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   int rval = 0;
   mpq_t tmpval;
   mpq_init(tmpval);

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(dval != NULL);

   SCIPdebugMessage("getting real parameter %d\n", type);

   switch( type )
   {
   case SCIP_LPPAR_OBJLIM:
      rval = mpq_QSget_param_EGlpNum(lpi->prob, QS_PARAM_OBJLLIM, &tmpval);
      break;
   case SCIP_LPPAR_LPTILIM:
      rval = mpq_QSget_param_EGlpNum(lpi->prob, QS_PARAM_SIMPLEX_MAX_TIME, &tmpval);
      break;
   default:
   case SCIP_LPPAR_MARKOWITZ:
   case SCIP_LPPAR_BARRIERCONVTOL:
   case SCIP_LPPAR_DUALFEASTOL:
   case SCIP_LPPAR_FEASTOL:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   *dval = mpq_get_d(tmpval);
   mpq_clear(tmpval);

   QS_RETURN(rval);
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiExactSetRealpar(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   int rval = 0;
   mpq_t tmpval;
   mpq_init(tmpval);
   mpq_set_d(tmpval, dval);

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("setting real parameter %d to %g\n", type, dval);

   switch( type )
   {
   case SCIP_LPPAR_LPTILIM:
      rval = mpq_QSset_param_EGlpNum(lpi->prob, QS_PARAM_SIMPLEX_MAX_TIME, tmpval);
      break;
   case SCIP_LPPAR_OBJLIM:
      rval = mpq_QSset_param_EGlpNum(lpi->prob, QS_PARAM_OBJLLIM, tmpval);
      break;
   case SCIP_LPPAR_FEASTOL:
   case SCIP_LPPAR_DUALFEASTOL:
   case SCIP_LPPAR_BARRIERCONVTOL:
   case SCIP_LPPAR_MARKOWITZ:
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   mpq_clear(tmpval);

   QS_RETURN(rval);
}

/**@} */

/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as positive infinity in the LP solver */
void SCIPlpiExactPosInfinity(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        infval              /**< pointer to store positive infinity value of LP solver */
   )
{
   assert(infval != NULL);

   SCIPrationalSetGMP(infval, mpq_ILL_MAXDOUBLE);
}

/** checks if given value is treated as positive infinity in the LP solver */
SCIP_Bool SCIPlpiExactIsPosInfinity(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        val                 /**< given value */
   )
{  /*lint --e{715} */
   return (mpq_cmp(*SCIPrationalGetGMP(val), mpq_ILL_MAXDOUBLE) >= 0);
}

/** returns value treated as negative infinity in the LP solver */
void SCIPlpiExactNegInfinity(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        infval              /**< pointer to store negative infinity value of LP solver */
   )
{  /*lint --e{715} */
   assert(infval != NULL);

   SCIPrationalSetGMP(infval, mpq_ILL_MINDOUBLE);
}

/** checks if given value is treated as negative infinity in the LP solver */
SCIP_Bool SCIPlpiExactIsNegInfinity(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_RATIONAL*        val                 /**< given value */
   )
{  /*lint --e{715} */
   return (mpq_cmp(*SCIPrationalGetGMP(val), mpq_ILL_MINDOUBLE) <= 0);
}

/** returns value treated as negative infinity in the LP solver */
SCIP_Real SCIPlpiExactInfinity(
   SCIP_LPIEXACT*        lpi                 /**< LP interface structure */
   )
{  /*lint --e{715} */
   assert(lpi != NULL);

   return mpq_get_d(mpq_ILL_MAXDOUBLE);
}

/** checks if given value is treated as negative infinity in the LP solver */
SCIP_Bool SCIPlpiExactIsInfinity(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< given value */
   )
{  /*lint --e{715} */
   return val >= mpq_get_d(mpq_ILL_MAXDOUBLE);
}

/**@} */

/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
SCIP_RETCODE SCIPlpiExactReadLP(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   int j;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("reading LP from file <%s>\n", fname);

   if( lpi->prob )
      mpq_QSfree_prob(lpi->prob);

   lpi->solstat = 0;
   lpi->previt = 0;

   /* try to extract file type */
   j = strlen(fname)-1;
   while( j >= 0 && fname[j] != '.'  )
      --j;
   if( fname[j] == '.' )
      ++j;

   /* load problem */
   lpi->prob = mpq_QSread_prob(fname, &(fname[j]));
   if( lpi->prob == 0 )
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiExactWriteLP(
   SCIP_LPIEXACT*        lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   int j;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("writing LP to file <%s>\n", fname);

   /* try to extract file type */
   j = strlen(fname) - 1;
   while( j >= 0 && fname[j] != '.' )
      --j;
   if( fname[j] == '.' )
      ++j;

   /* write problem */
   if( mpq_QSwrite_prob(lpi->prob, fname, &(fname[j])) )
      return SCIP_WRITEERROR;

   return SCIP_OKAY;
}

/** prints additional lpiex internal info */
void SCIPlpiExactPrintInfo(
   SCIP_LPIEXACT*        lpi                 /**< pointer to an LP interface structure */
   )
{
   mpq_lpinfo* lp;
   lp = lpi->prob->lp;
   SCIPerrorMessage("solstat= %d\n (solstat values: QS_LP_OPTIMAL=1, QS_LP_INFEASIBLE=2, QS_LP_UNBOUNDED=3, QS_LP_ITER_LIMIT=4, QS_LP_TIME_LIMIT=5, QS_LP_UNSOLVED=6, QS_LP_ABORTED=7, QS_LP_NUMERR=8, QS_LP_OBJ_LIMIT=9, QS_MODIFIED=100)\n", lpi->solstat );
   SCIPerrorMessage("probstat.optimal= %d\n", lp->probstat.optimal );
   SCIPerrorMessage("probstat.primal_feasible= %d\n", lp->probstat.primal_feasible );
   SCIPerrorMessage("probstat.primal_infeasible= %d\n", lp->probstat.primal_infeasible );
   SCIPerrorMessage("probstat.primal_unbounded= %d\n", lp->probstat.primal_unbounded );
   SCIPerrorMessage("probstat.dual_feasible= %d\n", lp->probstat.dual_feasible );
   SCIPerrorMessage("probstat.dual_infeasible= %d\n", lp->probstat.dual_infeasible );
   SCIPerrorMessage("probstat.dual_unbounded= %d\n", lp->probstat.dual_unbounded );
   SCIPerrorMessage("basisstat.primal_feasible= %d\n", lp->basisstat.primal_feasible );
   SCIPerrorMessage("basisstat.primal_infeasible= %d\n", lp->basisstat.primal_infeasible );
   SCIPerrorMessage("basisstat.dual_feasible= %d\n", lp->basisstat.dual_feasible );
   SCIPerrorMessage("basisstat.dual_infeasible= %d\n", lp->basisstat.dual_infeasible );
}
/**@} */
#endif
