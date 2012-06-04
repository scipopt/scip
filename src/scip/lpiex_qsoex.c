/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lpiex_qsoex.c
 * @brief  LP interface for QSopt_ex version >= 2.5.4 (r239)
 * @author Daniel Espinoza
 * @author Marc Pfetsch
 * @author Kati Wolter
*/

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define VERIFY_OUT  /** uncomment to get info of QSopt_ex about verifying dual feasibility of the basis */

//#define USEOBJLIM   /** uncomment to pass objlimit to exact lp solver; same as in cons_exactlp.c
//                     *  warning: QSopt_ex allows objlimits but the support is buggy; if the limit is reached, QSopt_ex
//                     *  does not stop but increasess the precision */

#include <string.h>
#ifdef WITH_GMP
#include "gmp.h"
#include "EGlib.h"
#include "QSopt_ex.h"
#endif

#include "scip/lpiex.h"
#include "scip/bitencode.h"
#include "scip/message.h"
#include "scip/misc.h"


#ifdef WITH_GMP

typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE

/** LP interface */
struct SCIP_LPiEx
{
   mpq_QSprob prob; /**< LP struct pointer */
   int solstat;	    /**< solution status of last optimization call */
   int previt;	    /**< previous number of simplex iterations performed */
   int rowspace;    /**< current size of internal row-related arrays */
   char* isen;	    /**< array of length rowspace */
   mpq_t* irhs;     /**< array of rhs rowspace */
   mpq_t* irng;     /**< array of range rowspace */
   int* ircnt;	    /**< array of count rowspace */
   int* irbeg;	    /**< array of begining index rowspace */
   int colspace;    /**< current size of internal column-related arrays */
   int* iccnt;	    /**< array of length colspace */
   char* iccha;	    /**< array of type colspace */
   int tbsz;	      /**< current size of tableau-related arrays */
   mpq_t* itab;     /**< array of length tbsz */
   char* ibas;	    /**< array of length tbsz */
   int pricing;	    /**< SCIP pricing option */
   mpq_factor_work* factor;  /**< factorized matrix  */
};

/** solver name */
static char __qsstr[1024];



/*
 * local defines
 */

/** print location of the calling line */
#define __QS_PRINTLOC__ fprintf(stderr,", in (%s:%d)\n", __FILE__, __LINE__);

/** This macro is to print error messages and jump to the given point in the code, it also print the
 * file and line where this happend */
#define QS_TESTG(A,B,...) do{{\
	 if (A){				\
	    fprintf(stderr,__VA_ARGS__);	\
	    __QS_PRINTLOC__;			\
	    goto B;}}}while(0)

/** This macro is to print error messages and to exit with SCIP_LPERROR */
#define QS_ERROR(A,...) do{{\
	 if (A){				\
	    fprintf(stderr,__VA_ARGS__);	\
	    __QS_PRINTLOC__;			\
	    return SCIP_LPERROR;}}}while(0)

/** return value macro, if the value is non-zero, write to standard error the returning code and
 * where this happened, and return SCIP_ERROR, otherwise return normal SCIP_OKAY termination code. */
#define QS_RETURN(A) do{\
      const int __RVAL__ = (A);						\
      if (__RVAL__){							\
	 fprintf(stderr,"LP Error: QSopt_ex returned %d",__RVAL__);	\
	 __QS_PRINTLOC__;						\
	 return SCIP_ERROR;}						\
      return SCIP_OKAY;}while(0)

/** return value macro, if the value is non-zero, write to standard error the returning code and
 * where this happened, and return SCIP_ERROR, otherwise do nothing. */
#define QS_CONDRET(A) do{\
      const int __RVAL__ = (A);						\
      if (__RVAL__){							\
	 fprintf(stderr,"LP Error: QSopt_ex returned %d",__RVAL__);	\
	 __QS_PRINTLOC__;						\
	 return SCIP_LPERROR;}						\
   }while(0)



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
 * local functions
 */

/** ensure size of column-related arrays */
static inline
SCIP_RETCODE ensureTabMem(
   SCIP_LPIEX* const lpi,
   int const sz
   )
{
   register int i;
   if (lpi->tbsz < sz)
   {
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->itab), sz*2) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->ibas), sz*2) );
      for (i = lpi->tbsz ; i < sz*2 ; i++)
	 mpq_init(lpi->itab[i]);
      lpi->tbsz = sz*2;
   }
   return SCIP_OKAY;
}

/** ensure size of column-related arrays */
static inline
SCIP_RETCODE ensureColMem(
   SCIP_LPIEX* const lpi,
   int const ncols
   )
{
   if (lpi->colspace < ncols)
   {
      lpi->colspace = ncols*2;
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->iccnt), lpi->colspace) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->iccha), lpi->colspace) );
   }
   return SCIP_OKAY;
}

/** ensure size of row-related arrays */
static inline
SCIP_RETCODE ensureRowMem(
   SCIP_LPIEX* const lpi,
   int const nrows
   )
{
   register int i;
   if (lpi->rowspace < nrows)
   {
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->isen), nrows*2) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->irhs), nrows*2) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->irng), nrows*2) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->ircnt), nrows*2) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->irbeg), nrows*2) );
      for (i = lpi->rowspace ; i < nrows*2; i++)
      {
	 mpq_init(lpi->irhs[i]);
	 mpq_init(lpi->irng[i]);
      }
      lpi->rowspace = nrows*2;
   }
   return SCIP_OKAY;
}

/** transform lhs/rhs into qsopt format */
static inline
SCIP_RETCODE convertSides(
   SCIP_LPIEX* const lpi,
   int const nrows,
   const mpq_t* const lhs,
   const mpq_t* const rhs
   )
{
   int state;
   register int i;

   for (i = 0; i < nrows; ++i)
   {
      state = (((mpq_cmp(lhs[i], mpq_ILL_MINDOUBLE) <= 0) ? 1U : 0U) | ((mpq_cmp(rhs[i], mpq_ILL_MAXDOUBLE) >= 0) ? 2U : 0U));
      lpi->ircnt[i] = 0;
      lpi->irbeg[i] = 0;
      switch (state)
      {
      case 0:
	 if ( mpq_equal(lhs[i], rhs[i]) )
	 {
	    lpi->isen[i] = 'E';
	    mpq_set(lpi->irhs[i], lhs[i]);
	    mpq_set_ui(lpi->irng[i], 0UL, 1UL);
	 }
	 else
	 {
	    lpi->isen[i] = 'R';
	    mpq_set(lpi->irhs[i], lhs[i]);
	    mpq_sub(lpi->irng[i], rhs[i], lhs[i]);
	    assert( mpq_sgn(lpi->irng[i]) >=0 );
	 }
	 break;
      case 1:
	 lpi->isen[i] = 'L';
	 mpq_set(lpi->irhs[i], rhs[i]);
	 mpq_set_ui(lpi->irng[i], 0UL, 1UL);
	 break;
      case 2:
	 lpi->isen[i] = 'G';
	 mpq_set(lpi->irhs[i], lhs[i]);
	 mpq_set_ui(lpi->irng[i], 0UL, 1UL);
	 break;
      default:
	 SCIPerrorMessage("Error, constraint %d has no bounds!",i);
	 SCIPABORT();
      }
   }
   return SCIP_OKAY;
}




/*
 * Miscellaneous Methods
 */


/**@name Miscellaneous Methods */
/**@{ */


/** gets name and version of LP solver */
const char* SCIPlpiexGetSolverName(void)
{
   char* vname = mpq_QSversion();
   snprintf (__qsstr, 1023, "%s", vname);
   mpq_QSfree(vname);
   return __qsstr;
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiexGetSolverPointer(
   SCIP_LPIEX*           lpi                 /**< pointer to an LP interface structure */
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

/** calls initializator of LP solver; this is mainly needed for defining constants in extended and rational precision */
void SCIPlpiexStart(
   void
   )
{
   assert(!__QSexact_setup);
   QSexactStart();
}

/** calls deinitializator of LP solver; this is needed for freeing all internal data of the solver, like constants in
 *  extended and rational precision
 */
void SCIPlpiexEnd(
   void
   )
{
   assert(__QSexact_setup);
   QSexactClear();
}

/** creates an LP problem object
 * @return SCIP_OK on success
 * */
SCIP_RETCODE SCIPlpiexCreate(
   SCIP_LPIEX**          lpi,                /**< pointer to an LP interface structure */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   register int i;

   /* QSopt_ex only works with bools as integers */
   assert(sizeof (SCIP_Bool) == sizeof (int));
   assert(lpi != NULL);

   SCIPdebugMessage("SCIPlpiexCreate()\n");

   /* create LP */
   SCIP_ALLOC( BMSallocMemory(lpi) );
   memset(*lpi, 0, sizeof(struct SCIP_LPiEx));

   /* factor work is NULL unless used */
   (*lpi)->factor =  (mpq_factor_work*) NULL;

   (*lpi)->prob = mpq_QScreate_prob(name, (int) objsen);
   if( (*lpi)->prob == NULL )
   {
      SCIPerrorMessage("No memory\n");
      return SCIP_LPERROR;
   }

   (*lpi)->rowspace = 1024;
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->isen),1024) );
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->irhs),1024) );
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->irng),1024) );
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->irbeg),1024) );
   SCIP_ALLOC( BMSallocMemoryArray(&((*lpi)->ircnt),1024) );

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
SCIP_RETCODE SCIPlpiexFree(
   SCIP_LPIEX**            lpi                 /**< pointer to an LP interface structure */
   )
{
   register int i;
   assert(lpi != NULL);
   assert(*lpi != NULL);

   SCIPdebugMessage("SCIPlpiexFree()\n");

   /* free factor work */
   if( (*lpi)->factor != NULL )
   {
      mpq_ILLfactor_free_factor_work((*lpi)->factor);
      BMSfreeMemoryArray( &((*lpi)->factor) );
   }

   /* free LP */
   mpq_QSfree_prob((*lpi)->prob);
   for (i = 0; i < (*lpi)->tbsz; ++i)
      mpq_clear((*lpi)->itab[i]);
   for (i = 0; i < (*lpi)->rowspace; ++i)
   {
      mpq_clear((*lpi)->irng[i]);
      mpq_clear((*lpi)->irhs[i]);
   }
   BMSfreeMemoryArray( &((*lpi)->isen) );
   BMSfreeMemoryArray( &((*lpi)->irhs) );
   BMSfreeMemoryArray( &((*lpi)->irng) );
   BMSfreeMemoryArray( &((*lpi)->ircnt) );
   BMSfreeMemoryArray( &((*lpi)->irbeg) );
   BMSfreeMemoryArray( &((*lpi)->iccnt) );
   BMSfreeMemoryArray( &((*lpi)->iccha) );
   BMSfreeMemoryArray( &((*lpi)->itab) );
   BMSfreeMemoryArray( &((*lpi)->ibas) );

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
SCIP_RETCODE SCIPlpiexLoadColLP(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   mpq_t*                obj,                /**< objective function values of columns */
   mpq_t*                lb,                 /**< lower bounds of columns */
   mpq_t*                ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   const mpq_t*          lhs,                /**< left hand sides of rows */
   const mpq_t*          rhs,                /**< right hand sides of rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each column in ind- and val-array */
   int*                  ind,                /**< row indices of constraint matrix entries */
   mpq_t*                val                 /**< values of constraint matrix entries */
   )
{
   register int i;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("loading LP in column format into QSopt_ex: %d cols, %d rows\n", ncols, nrows);

   /* delete old LP */
   SCIP_CALL( SCIPlpiexClear(lpi) );

   /* set sense */
   if (objsen == SCIP_OBJSEN_MAXIMIZE)
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

   /* ensure column size */
   SCIP_CALL( ensureColMem(lpi, ncols) );

   /* compute column lengths */
   for (i = 0; i < ncols-1; ++i)
   {
      lpi->iccnt[i] = beg[i+1] - beg[i];
      assert( lpi->iccnt[i] >= 0 );
   }
   if ( ncols > 0 )
   {
      lpi->iccnt[ncols-1] = nnonz - beg[ncols-1];
      assert( lpi->iccnt[ncols-1] >= 0 );
   }

   /* and add the columns */
   rval = mpq_QSadd_cols(lpi->prob, ncols, lpi->iccnt, beg, ind, val, obj, lb, ub, (const char**)colnames);

   QS_RETURN(rval);
}


#if 0 /* old version with some minor fix concering nnonz=0, but which does not work, beg=NULL is still not allowed */
/** adds columns to the LP */
SCIP_RETCODE SCIPlpiexAddCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   const mpq_t*          obj,                /**< objective function values of new columns */
   const mpq_t*          lb,                 /**< lower bounds of new columns */
   const mpq_t*          ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const mpq_t*          val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   int rval = 0;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("adding %d columns with %d nonzeros to QSopt_ex\n", ncols, nnonz);

   lpi->solstat = 0;

   /* ensure column size */
   SCIP_CALL(ensureColMem(lpi, ncols));

   /* compute column lengths */
   for (i = 0; i < ncols - 1; ++i)
   {
      if( nnonz > 0 )
         lpi->iccnt[i] = beg[i+1] - beg[i];
      else
         lpi->iccnt[i] = 0;
      assert(lpi->iccnt[i] >= 0);
   }
   if ( ncols > 0 )
   {
      if( nnonz > 0 )
         lpi->iccnt[ncols-1] = nnonz - beg[ncols-1];
      else
         lpi->iccnt[ncols-1] = 0;

      assert( lpi->iccnt[ncols-1] >= 0 );
   }

   /* and add the columns */
   rval = mpq_QSadd_cols(lpi->prob, ncols, lpi->iccnt, beg, ind, val, obj, lb, ub, (const char**)colnames);

   QS_RETURN(rval);
}
#else
/** @todo exip: check whether I implemented handling of case beg=ind=val=NULL correctly */
/** adds columns to the LP */
SCIP_RETCODE SCIPlpiexAddCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   mpq_t*                obj,                /**< objective function values of new columns */
   mpq_t*                lb,                 /**< lower bounds of new columns */
   mpq_t*                ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   int*                  beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   int*                  ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   mpq_t*                val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   int rval = 0;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("adding %d columns with %d nonzeros to QSopt_ex\n", ncols, nnonz);

   lpi->solstat = 0;

   /* ensure column size */
   SCIP_CALL(ensureColMem(lpi, ncols));

   if( nnonz > 0 )
   {
      /* compute column lengths */
      for (i = 0; i < ncols - 1; ++i)
      {
         lpi->iccnt[i] = beg[i+1] - beg[i];
         assert(lpi->iccnt[i] >= 0);
      }
      if ( ncols > 0 )
      {
         lpi->iccnt[ncols-1] = nnonz - beg[ncols-1];
         assert( lpi->iccnt[ncols-1] >= 0 );
      }

      /* and add the columns */
      rval = mpq_QSadd_cols(lpi->prob, ncols, lpi->iccnt, beg, ind, val, obj, lb, ub, (const char**)colnames);
   }
   else
   {
      for( i = 0; i < ncols; ++i )
      {
         rval = mpq_QSnew_col(lpi->prob, obj[i], lb[i], ub[i], (const char*) colnames[i]);
      }
   }

   QS_RETURN(rval);
}
#endif

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiexDelCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
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

   SCIP_CALL(ensureColMem(lpi, len));
   for (i = firstcol ; i <= lastcol ; i++)
      lpi->iccnt[i-firstcol] = i;

   rval = mpq_QSdelete_cols(lpi->prob, len, lpi->iccnt);

   QS_RETURN(rval);
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiexDelColset(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
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

   for (i=0, ccnt=0; i < ncols; i++)
   {
      if (dstat[i])
	 dstat[i] = -1;
      else
	 dstat[i] = ccnt++;
   }
   return SCIP_OKAY;
}


/** adds rows to the LP */
SCIP_RETCODE SCIPlpiexAddRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const mpq_t*          lhs,                /**< left hand sides of new rows */
   const mpq_t*          rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   int*                  beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            len,                /**< number of nonzeros of each row in ind- and val-array, or NULL if only nonzeros */
   int*                  ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   mpq_t*                val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   register int i;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("adding %d rows with %d nonzeros to QSopt_ex\n", nrows, nnonz);

   /* add rows with no matrix, and then the columns, first ensure space */
   SCIP_CALL( ensureRowMem (lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   SCIP_CALL( convertSides(lpi, nrows, lhs, rhs) );

   /* compute row count */
   if( len != NULL )
   {
      for( i = 0; i < nrows; i++ )
         lpi->ircnt[i] = len[i];
   }
   else
   {
      for( i = 0; i < nrows-1; i++ )
      {
         lpi->ircnt[i] = beg[i+1] - beg[i];
         assert(lpi->ircnt[i] >= 0);
      }
      if( nrows > 0 )
      {
         lpi->ircnt[nrows-1] = nnonz - beg[nrows-1];
         assert(lpi->ircnt[nrows-1] >= 0);
      }
   }

   /* now we add the rows */
   rval = mpq_QSadd_ranged_rows(lpi->prob, nrows, lpi->ircnt, beg, ind, (const mpq_t*) val, (const mpq_t*) lpi->irhs,
      lpi->isen, (const mpq_t*) lpi->irng, (const char**)rownames);
   QS_ERROR(rval, "failed adding %d rows with %d non-zeros", nrows, nnonz);

   return SCIP_OKAY;
}


/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiexDelRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
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
   for (i = firstrow; i <= lastrow; i++)
      lpi->ircnt[i-firstrow] = i;
   rval = mpq_QSdelete_rows(lpi->prob, len, lpi->ircnt);

   QS_RETURN(rval);
}


/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiexDelRowset(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
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

   for (i = 0; i < nrows; ++i)
   {
      if (dstat[i] == 1)
	 ndel++;
   }

   SCIPdebugMessage("deleting a row set from QSopt_ex (%d)\n",ndel);

   rval = mpq_QSdelete_setrows(lpi->prob,dstat);
   QS_CONDRET(rval);

   for (i=0, ccnt=0; i < nrows; ++i)
   {
      if (dstat[i])
	 dstat[i] = -1;
      else
	 dstat[i] = ccnt++;
   }
   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiexClear(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
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
   if (ncols >= 1)
   {
      SCIP_CALL( ensureColMem(lpi,ncols) );
      for (i = 0; i < ncols; ++i)
	 lpi->iccnt[i] = i;
      rval = mpq_QSdelete_cols(lpi->prob, ncols, lpi->iccnt);
      QS_CONDRET(rval);
   }

   if (nrows >= 1)
   {
      SCIP_CALL( ensureRowMem(lpi, nrows) );
      for (i = 0; i < nrows; ++i)
	 lpi->ircnt[i] = i;
      rval = mpq_QSdelete_rows(lpi->prob, nrows, lpi->ircnt);
      QS_CONDRET(rval);
   }
   return SCIP_OKAY;
}


/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiexChgBounds(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   int*                  ind,                /**< column indices */
   mpq_t*                lb,                 /**< values for the new lower bounds, or NULL */
   mpq_t*                ub                  /**< values for the new upper bounds, or NULL */
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

      for (j = 0; j < ncols; ++j)
      {
         if( lb == NULL)
            gmp_snprintf(s, SCIP_MAXSTRLEN, "  col %d: [--,%Qd]\n", ind[j], ub[j]);
         else if( ub == NULL )
            gmp_snprintf(s, SCIP_MAXSTRLEN, "  col %d: [%Qd,--]\n", ind[j], lb[j]);
         else
            gmp_snprintf(s, SCIP_MAXSTRLEN, "  col %d: [%Qd,%Qd]\n", ind[j], lb[j], ub[j]);
         SCIPdebugPrintf(s);
      }
   }
#endif

   SCIP_CALL(ensureColMem(lpi, ncols));

   if( lb != NULL )
   {
      for (i = 0; i < ncols; ++i)
         lpi->iccha[i] = 'L';

      rval = mpq_QSchange_bounds(lpi->prob, ncols, ind, lpi->iccha, (const mpq_t*) lb);
      QS_CONDRET(rval);
   }

   if( ub != NULL )
   {
      for (i = 0; i < ncols; ++i)
         lpi->iccha[i] = 'U';

      rval = mpq_QSchange_bounds(lpi->prob, ncols, ind, lpi->iccha, (const mpq_t*) ub);
   }
   QS_RETURN(rval);
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiexChgSides(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const mpq_t*          lhs,                /**< new values for left hand sides */
   const mpq_t*          rhs                 /**< new values for right hand sides */
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
   for (i = 0; i < nrows; ++i)
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

   return SCIP_OKAY;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiexChgCoef(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   mpq_t                 newval              /**< new value of coefficient */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;

   SCIPdebugMessage("changing coefficient row %d, column %d in QSopt_ex to %g\n", row, col, mpq_get_d(newval));

   rval = mpq_QSchange_coef(lpi->prob, row, col, newval);

   QS_RETURN(rval);
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiexChgObjsen(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("changing objective sense in QSopt_ex to %d\n", objsen);

   /* set sense */
   if (objsen == SCIP_OBJSEN_MAXIMIZE)
   {
      rval = mpq_QSchange_objsense(lpi->prob, QS_MAX);
      QS_CONDRET(rval);
   }
   else
   {
      rval = mpq_QSchange_objsense(lpi->prob, QS_MIN);
      QS_CONDRET(rval);
   }
   return SCIP_OKAY;
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiexChgObj(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   mpq_t*                obj                 /**< new objective values for columns */
   )
{
   register int i;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("changing %d objective values in QSopt_ex\n", ncols);

   for (i = 0; i < ncols; ++i)
   {
      rval = mpq_QSchange_objcoef(lpi->prob, ind[i], obj[i]);
      QS_CONDRET(rval);
   }
   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiexScaleRow(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   const mpq_t           scaleval            /**< scaling multiplier */
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
   SCIPdebugMessage("scaling row %d with factor %g in QSopt_ex\n", row, mpq_get_d(scaleval));

   rowlist[0] = row;
   /* get row */
   rval = mpq_QSget_ranged_rows_list(lpi->prob, 1, rowlist, &rowcnt, &rowbeg, &rowind, &rowval, &rhs, &sense, &range, 0);
   QS_TESTG(rval, CLEANUP, " ");

   /* change all coefficients in the constraint */
   for (i = 0; i < rowcnt[0]; ++i)
   {
      mpq_mul(svl, rowval[i], scaleval);
      rval = mpq_QSchange_coef(lpi->prob, row, rowind[i], svl);
      QS_TESTG(rval, CLEANUP, " ");
   }

   /* if we have a positive scalar, we just scale rhs and range */
   if (mpq_sgn(scaleval) >= 0)
   {
      mpq_mul(svl,scaleval,rhs[0]);
      rval = mpq_QSchange_rhscoef(lpi->prob, row, svl);
      QS_TESTG(rval, CLEANUP, " ");
      if (sense[0] == 'R')
      {
	 mpq_mul(svl, range[0], scaleval);
	 rval = mpq_QSchange_range(lpi->prob, row, svl);
	 QS_TESTG(rval, CLEANUP, " ");
      }
   }
   /* otherwise, we must change everything */
   else
   {
      switch(sense[0])
      {
      case 'E':
	 mpq_mul(svl,rhs[0],scaleval);
	 rval = mpq_QSchange_rhscoef(lpi->prob, row, svl);
	 QS_TESTG(rval, CLEANUP, " ");
	 break;
      case 'L':
	 mpq_mul(svl,rhs[0],scaleval);
	 rval = mpq_QSchange_rhscoef(lpi->prob, row, svl);
	 QS_TESTG(rval, CLEANUP, " ");
	 rval = mpq_QSchange_sense(lpi->prob, row, 'G');
	 QS_TESTG(rval, CLEANUP, " ");
	 break;
      case 'G':
	 mpq_mul(svl,rhs[0],scaleval);
	 rval = mpq_QSchange_rhscoef(lpi->prob, row, svl);
	 QS_TESTG(rval, CLEANUP, " ");
	 rval = mpq_QSchange_sense(lpi->prob, row, 'L');
	 QS_TESTG(rval, CLEANUP, " ");
	 break;
      case 'R':
	 mpq_add(svl,rhs[0],range[0]);
	 mpq_mul(svl,svl,scaleval);
	 rval = mpq_QSchange_rhscoef(lpi->prob, row, svl);
	 QS_TESTG(rval, CLEANUP, " ");
	 mpq_abs(svl,scaleval);
	 mpq_mul(svl,svl,range[0]);
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
SCIP_RETCODE SCIPlpiexScaleCol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   const mpq_t           scaleval            /**< scaling multiplier */
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
   SCIPdebugMessage("scaling column %d with factor %g in QSopt_ex\n", col,
	 										mpq_get_d(scaleval));

   /* get the column */
   collist[0] = col;
   rval = mpq_QSget_columns_list(lpi->prob, 1, collist, &colcnt, &colbeg, &colind, &colval, &obj, &lb, &ub, 0);
   QS_TESTG(rval,CLEANUP," ");

   /* scale column coefficients */
   for (i = 0; i < colcnt[0]; ++i)
   {
      mpq_mul(svl, colval[i], scaleval);
      rval = mpq_QSchange_coef(lpi->prob, colind[i], col, svl);
      QS_TESTG(rval,CLEANUP," ");
   }

   /* scale objective value */
   mpq_mul(svl,obj[0],scaleval);
   rval = mpq_QSchange_objcoef(lpi->prob, col, svl);
   QS_TESTG(rval,CLEANUP," ");

   /* scale column bounds */
   if (mpq_sgn(scaleval) < 0)
   {
      mpq_set(obj[0],lb[0]);
      mpq_neg(lb[0],ub[0]);
      mpq_neg(ub[0],obj[0]);
   }
   if (mpq_cmp(lb[0],mpq_ILL_MINDOUBLE)>0)
   {
      mpq_abs(svl,scaleval);
      mpq_mul(lb[0],lb[0],svl);
   }
   if (mpq_cmp(ub[0],mpq_ILL_MAXDOUBLE)<0)
   {
      mpq_abs(svl,scaleval);
      mpq_mul(ub[0],ub[0],svl);
   }

   if (mpq_cmp(lb[0],mpq_ILL_MINDOUBLE)<0)
      mpq_set(lb[0],mpq_ILL_MINDOUBLE);
   if (mpq_cmp(ub[0],mpq_ILL_MAXDOUBLE)>0)
      mpq_set(ub[0],mpq_ILL_MAXDOUBLE);

   rval = mpq_QSchange_bound(lpi->prob, col, 'L', lb[0]);
   QS_TESTG(rval,CLEANUP," ");
   rval = mpq_QSchange_bound(lpi->prob, col, 'U', ub[0]);
   QS_TESTG(rval,CLEANUP," ");

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
SCIP_RETCODE SCIPlpiexGetNRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
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
SCIP_RETCODE SCIPlpiexGetNCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
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
SCIP_RETCODE SCIPlpiexGetNNonz(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
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
SCIP_RETCODE SCIPlpiexGetCols(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   mpq_t*                lb,                 /**< buffer to store the lower bound vector, or NULL */
   mpq_t*                ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   mpq_t*                val                 /**< buffer to store values of constraint matrix entries, or NULL */
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
   for (i = 0; i < len; ++i)
      lpi->iccnt[i] = i + firstcol;

   /* get data from qsopt */
   rval = mpq_QSget_columns_list(lpi->prob, len, lpi->iccnt, nnonz ? (&lcnt) : 0, nnonz ? (&lbeg) : 0, nnonz ? (&lind) : 0,
      nnonz ? (&lval) : 0, 0, lb ? (&llb) : 0, lb ? (&lub) : 0, 0);

   QS_TESTG(rval, CLEANUP, " ");

   /* store in the user-provided data */
   if (nnonz)
   {
      assert( lbeg != NULL );
      assert( lcnt != NULL );
      assert( lind != NULL );
      assert( lval != NULL );

      *nnonz = lbeg[len-1] + lcnt[len-1];
      for (i = 0 ; i < len ; i++)
	 beg[i] = lbeg[i];  /*lint !e613*/
      for (i = 0; i < *nnonz; ++i)
      {
	 ind[i] = lind[i];  /*lint !e613*/
	 mpq_set(val[i], lval[i]);  /*lint !e613*/
      }
   }
   if (lb)
   {
      assert( llb != NULL );
      assert( lub != NULL );

      for (i = 0; i < len; ++i)
      {
	 mpq_set(lb[i], llb[i]);
	 mpq_set(ub[i], lub[i]);   /*lint !e613*/
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
SCIP_RETCODE SCIPlpiexGetRows(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   mpq_t*                lhs,                /**< buffer to store left hand side vector, or NULL */
   mpq_t*                rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   mpq_t*                val                 /**< buffer to store values of constraint matrix entries, or NULL */
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
   for (i = 0; i < len; ++i)
      lpi->ircnt[i] = i + firstrow;

   /* get data from qsopt */
   rval = mpq_QSget_ranged_rows_list(lpi->prob, len, lpi->ircnt, nnonz ? (&lcnt) : 0, nnonz ? (&lbeg) : 0, nnonz ? (&lind) : 0,
      nnonz ? (&lval) : 0, rhs ? (&lrhs) : 0, rhs ? (&lsense) : 0, rhs ? (&lrng) : 0, 0);
   QS_TESTG(rval, CLEANUP, " ");

   /* store in the user-provided data */
   if (nnonz)
   {
      assert( lbeg != NULL );
      assert( lcnt != NULL );
      assert( lind != NULL );
      assert( lval != NULL );

      *nnonz = lbeg[len-1] + lcnt[len-1];
      for (i = 0; i < len; i++)
	 beg[i] = lbeg[i];  /*lint !e613*/
      for (i = 0; i < *nnonz; ++i)
      {
	 ind[i] = lind[i];  /*lint !e613*/
	 mpq_set(val[i], lval[i]);  /*lint !e613*/
      }
   }
   if (rhs)
   {
      assert( lrhs != NULL );
      assert( lrng != NULL );
      assert( lsense != NULL );

      for (i = 0; i < len; ++i)
      {
	 switch (lsense[i])
	 {
	 case 'R':
	    mpq_set(lhs[i], lrhs[i]);            /*lint !e613*/
	    mpq_add(rhs[i], lrhs[i], lrng[i]);  /*lint !e613*/
	    break;
	 case 'E':
	    mpq_set(lhs[i], lrhs[i]);
	    mpq_set(rhs[i], lrhs[i]);   /*lint !e613*/
	    break;
	 case 'L':
	    mpq_set(rhs[i], lrhs[i]);            /*lint !e613*/
	    mpq_set(lhs[i], mpq_ILL_MINDOUBLE);      /*lint !e613*/
	    break;
	 case 'G':
	    mpq_set(lhs[i], lrhs[i]);            /*lint !e613*/
	    mpq_set(rhs[i], mpq_ILL_MAXDOUBLE);       /*lint !e613*/
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
SCIP_RETCODE SCIPlpiexGetObj(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   mpq_t*                vals                /**< array to store objective coefficients */
   )
{
   const int len = lastcol - firstcol + 1;
   int rval = 0;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < mpq_QSget_colcount (lpi->prob));

   SCIPdebugMessage("getting objective values %d to %d\n", firstcol, lastcol);

   /* build col-list */
   SCIP_CALL(ensureColMem(lpi,len));
   for (i = 0; i < len; ++i)
      lpi->iccnt[i] = i + firstcol;

   /* get data from qsopt */
   rval = mpq_QSget_obj_list(lpi->prob, len, lpi->iccnt, vals);

   QS_RETURN(rval);
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiexGetBounds(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective value for */
   int                   lastcol,            /**< last column to get objective value for */
   mpq_t*                lbs,                /**< array to store lower bound values, or NULL */
   mpq_t*                ubs                 /**< array to store upper bound values, or NULL */
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
   SCIP_CALL(ensureColMem(lpi,len));
   for (i = 0; i < len; ++i)
      lpi->iccnt[i] = i + firstcol;

   /* get data from qsopt */
   rval = mpq_QSget_bounds_list(lpi->prob, len, lpi->iccnt, lbs, ubs);

   QS_RETURN(rval);
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiexGetSides(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   mpq_t*                lhss,               /**< array to store left hand side values, or NULL */
   mpq_t*                rhss                /**< array to store right hand side values, or NULL */
   )
{
   const int len = lastrow - firstrow + 1;
   register int i;
   mpq_t* lrhs=0, *lrng=0;
   int rval = 0;
   char* lsense=0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < mpq_QSget_rowcount (lpi->prob));
   assert(rhss != 0 && lhss != 0);

   SCIPdebugMessage("getting row sides %d to %d\n", firstrow, lastrow);

   /* build row-list */
   SCIP_CALL( ensureRowMem(lpi, len) );
   for (i = 0; i < len; ++i)
      lpi->ircnt[i] = i + firstrow;

   /* get data from qsopt */
   rval = mpq_QSget_ranged_rows_list(lpi->prob, len, lpi->ircnt, 0, 0, 0, 0, &lrhs,	&lsense, &lrng, 0);
   QS_TESTG(rval, CLEANUP, " ");

   /* store in the user-provided data */
   for (i = 0; i < len; ++i)
   {
      switch (lsense[i])
      {
      case 'R':
	 mpq_set(lhss[i], lrhs[i]);
	 mpq_add(rhss[i], lrhs[i], lrng[i]);
	 break;
      case 'E':
	 mpq_set(lhss[i], lrhs[i]);
	 mpq_set(rhss[i], lrhs[i]);
	 break;
      case 'L':
	 mpq_set(rhss[i], lrhs[i]);
	 mpq_set(lhss[i], mpq_ILL_MINDOUBLE);
	 break;
      case 'G':
	 mpq_set(lhss[i], lrhs[i]);
	 mpq_set(rhss[i], mpq_ILL_MAXDOUBLE);
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
SCIP_RETCODE SCIPlpiexGetCoef(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   mpq_t*                val                 /**< pointer to store the value of the coefficient */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting coefficient of row %d col %d\n", row, col);

   rval = mpq_QSget_coef(lpi->prob, row, col, val);

   QS_RETURN(rval);
}

/**@} */

/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiexSolvePrimal(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
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
   if (B)
      mpq_QSfree_basis(B);

   QS_RETURN(rval);
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiexSolveDual(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
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
   if (B)
      mpq_QSfree_basis(B);

   QS_RETURN(rval);
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiexSolveBarrier(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{  /*lint --e{715}*/
   return SCIPlpiexSolveDual(lpi);
}

/** performs strong branching iterations on all candidates */
SCIP_RETCODE SCIPlpiexStrongbranch(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   const mpq_t           psol,               /**< current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   mpq_t*                down,               /**< stores dual bound after branching column down */
   mpq_t*                up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIPdebugMessage("Strong branching not implemented for QSopt_ex\n");
   return SCIP_LPERROR;

#if 0
   int rval = 0, nit;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   SCIPdebugMessage("calling QSopt_ex strongbranching on variable %d (%d it lim)\n", col, itlim);

   /* results of QSopt_ex are valid in any case */
   *downvalid = TRUE;
   *upvalid = TRUE;

   /* call QSopt_ex */
   rval = mpq_QSopt_strongbranch(lpi->prob, 1, &col, &psol, down, up, itlim, mpq_ILL_MAXDOUBLE);
   QS_CONDRET(rval);

   rval = mpq_QSget_itcnt(lpi->prob, 0, 0, 0, 0, &nit);
   QS_CONDRET(rval);

   *iter = nit - lpi->previt;
   lpi->previt = nit;

   return SCIP_OKAY;
#endif
}

/**@} */

/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiexWasSolved(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{
   assert(lpi!=NULL);
   assert(lpi->prob!=NULL);

   return (lpi->solstat != 0 && lpi->solstat != QS_LP_MODIFIED && lpi->solstat != QS_LP_CHANGE_PREC);
}

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_RETCODE SCIPlpiexGetSolFeasibility(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting solution feasibility\n");

   if ( lpi->solstat == QS_LP_OPTIMAL || lpi->solstat == QS_LP_UNBOUNDED)
      *primalfeasible = 1;

   /* @todo: check why we can conclude dual feasibility from primal infeasibility. in theory, the LP could be primal and
    * dual infeasible as well; see also SCIPlpiexIsDualFeasible() and SCIPlpiexIsDualInfeasible()
    */
#ifdef USEOBJLIM
   if ( lpi->solstat == QS_LP_OPTIMAL || lpi->solstat == QS_LP_INFEASIBLE || lpi->solstat == QS_LP_OBJ_LIMIT )
#else
   if ( lpi->solstat == QS_LP_OPTIMAL || lpi->solstat == QS_LP_INFEASIBLE )
#endif
      *dualfeasible = 1;

   return SCIP_OKAY;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiexExistsPrimalRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
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
SCIP_Bool SCIPlpiexHasPrimalRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for primal ray\n");

   /* the current version of QSopt_ex can not give a primal certificate of unboundness */
   return FALSE;
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiexIsPrimalUnbounded(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for primal unboundness\n");

   return (lpi->solstat == QS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiexIsPrimalInfeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for primal infeasibility\n");

   return (lpi->solstat == QS_LP_INFEASIBLE);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiexIsPrimalFeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
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
SCIP_Bool SCIPlpiexExistsDualRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
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
SCIP_Bool SCIPlpiexHasDualRay(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual ray availability\n");

   return (lpi->solstat == QS_LP_INFEASIBLE);
}

/** returns TRUE iff LP is dual unbounded */
SCIP_Bool SCIPlpiexIsDualUnbounded(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual unboundness\n");

   return (lpi->solstat == QS_LP_INFEASIBLE);
}

/** returns TRUE iff LP is dual infeasible */
SCIP_Bool SCIPlpiexIsDualInfeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual infeasibility\n");

   return (lpi->solstat == QS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiexIsDualFeasible(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
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
SCIP_Bool SCIPlpiexIsOptimal(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for optimality\n");

   return (lpi->solstat == QS_LP_OPTIMAL);
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiexIsStable(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for numerical stability\n");

   return (lpi->solstat != QS_LP_NUMERR);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiexIsObjlimExc(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
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
SCIP_Bool SCIPlpiexIsIterlimExc(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for iteration limit exceeded\n");

   return (lpi->solstat == QS_LP_ITER_LIMIT);
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiexIsTimelimExc(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for time limit exceeded\n");

   return (lpi->solstat == QS_LP_TIME_LIMIT);
}

/** returns the internal solution status of the solver */
int SCIPlpiexGetInternalStatus(
   SCIP_LPIEX*           lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting internal solution status\n");

   return lpi->solstat;
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiexIgnoreInstability(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
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
SCIP_RETCODE SCIPlpiexGetObjval(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   mpq_t*                objval              /**< stores the objective value */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting solution's objective value\n");

   rval = mpq_QSget_objval(lpi->prob, objval);

   QS_RETURN(rval);
}

/** gets primal and dual solution vectors */
SCIP_RETCODE SCIPlpiexGetSol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   mpq_t*                objval,             /**< stores the objective value, may be NULL if not needed */
   mpq_t*                primsol,            /**< primal solution vector, may be NULL if not needed */
   mpq_t*                dualsol,            /**< dual solution vector, may be NULL if not needed */
   mpq_t*                activity,           /**< row activity vector, may be NULL if not needed */
   mpq_t*                redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{
   int rval = 0, nrows;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting solution\n");

   nrows = mpq_QSget_rowcount(lpi->prob);
   SCIP_CALL( ensureRowMem(lpi, nrows) );

   rval = mpq_QSget_solution(lpi->prob, objval, primsol, dualsol, lpi->irng, redcost);
   QS_CONDRET(rval);

   rval = mpq_QSget_rhs(lpi->prob, lpi->irhs);
   QS_CONDRET(rval);
   rval = mpq_QSget_senses(lpi->prob, lpi->isen);
   QS_CONDRET(rval);

   /* build back the activity */
   if ( activity )
   {
      for (i = 0; i < nrows; ++i)
      {
	 switch (lpi->isen[i])
	 {
	 case 'R':
	 case 'E':
	 case 'G':
	    mpq_add(activity[i], lpi->irhs[i], lpi->irng[i]);
	    break;
	 case 'L':
	    mpq_sub(activity[i], lpi->irhs[i], lpi->irng[i]);
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
SCIP_RETCODE SCIPlpiexGetPrimalRay(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   mpq_t*                ray                 /**< primal ray */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPerrorMessage("SCIPlpiexGetPrimalRay() not supported by QSopt_ex.\n");

   return SCIP_ERROR;
}

/** gets dual farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiexGetDualfarkas(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   mpq_t*                dualfarkas          /**< dual farkas row multipliers */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(dualfarkas != NULL);

   SCIPdebugMessage("calling QSopt_ex dual farkas: %d cols, %d rows, %d non zeros\n", mpq_QSget_colcount (lpi->prob),
      mpq_QSget_rowcount(lpi->prob), mpq_QSget_nzcount(lpi->prob));

   rval = mpq_QSget_infeas_array(lpi->prob, dualfarkas);

   QS_RETURN(rval);
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiexGetIterations(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
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

   return SCIP_OKAY;
}

/**@} */

/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiexGetBase(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
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

   SCIPdebugMessage("saving QSopt_ex basis into %p/%p\n", cstat, rstat);

   ncols = mpq_QSget_colcount(lpi->prob);
   nrows = mpq_QSget_rowcount(lpi->prob);

   SCIP_CALL(ensureTabMem(lpi, nrows + ncols));

   icstat = lpi->ibas;
   irstat = lpi->ibas+ncols;
   rval = mpq_QSget_basis_array(lpi->prob, icstat, irstat);
   QS_CONDRET(rval);

   /* now we must transform QSopt_ex codes into SCIP codes */
   for (i = 0; i < nrows; ++i)
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
   for (i = 0; i < ncols; ++i)
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
   return SCIP_OKAY;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiexSetBase(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   )
{
   int rval = 0, ncols, nrows;
   register int i;
   char* icstat=0, *irstat = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("loading basis %p/%p into QSopt_ex\n", cstat, rstat);

   ncols = mpq_QSget_colcount(lpi->prob);
   nrows = mpq_QSget_rowcount(lpi->prob);

   SCIP_CALL(ensureTabMem(lpi, ncols));

   icstat = lpi->ibas;
   irstat = lpi->ibas + ncols;

   /* now we must transform QSopt_ex codes into SCIP codes */
   for (i = 0; i < nrows; ++i)
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
	 irstat[i] = QS_ROW_BSTAT_UPPER; /*lint !e641*/
	 break;
      default:
	 SCIPerrorMessage("Unknown row basic status %d", rstat[i]);
	 SCIPABORT();
      }
   }
   for (i = 0; i < ncols; ++i)
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
SCIP_RETCODE SCIPlpiexGetBasisInd(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int*                  ind                 /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   int rval = 0, nrows, ncols;
   register int i;

   assert(lpi!=NULL);
   assert(lpi->prob!=NULL);

   SCIPdebugMessage("getting basis information\n");

   nrows = mpq_QSget_rowcount(lpi->prob);
   ncols = mpq_QSget_colcount(lpi->prob);
   rval = mpq_QSget_basis_order(lpi->prob, ind);
   QS_CONDRET(rval);

   /* transform QSopt_ex basis header into SCIP format */
   for (i = 0; i < nrows; ++i)
   {
      if (ind[i] >= ncols)
	 ind[i] = -(ind[i] - ncols - 1);
   }

   return SCIP_OKAY;
}

/** get dense row of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiexGetBInvRow(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   mpq_t*                coef                /**< pointer to store the coefficients of the row */
   )
{
   int rval = 0;

   assert(lpi!=NULL);
   assert(lpi->prob!=NULL);
   SCIPdebugMessage("getting binv-row %d from Qsopt %d cols, %d rows, %d nonz\n", r, mpq_QSget_colcount(lpi->prob),
      mpq_QSget_rowcount(lpi->prob), mpq_QSget_nzcount(lpi->prob));

   rval = mpq_QSget_binv_row(lpi->prob, r, coef);
   QS_RETURN(rval);
}

/** get dense column of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiexGetBInvCol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiexGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   mpq_t*                coef                /**< pointer to store the coefficients of the column */
   )
{  /*lint --e{715} */
   assert(lpi!=NULL);
   assert(lpi->prob!=NULL);

   SCIPerrorMessage("SCIPlpiexGetBInvCol() not supported by QSopt_ex.\n");

   /* QSopt_ex does not provide an interface for this yet */
   return SCIP_ERROR;
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiexGetBInvARow(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const mpq_t*          binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiexGetBInvRow(), or NULL */
   mpq_t*                coef                /**< vector to return coefficients */
   )
{  /*lint --e{715} */
   int rval = 0,ncols,nrows;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting binva-row %d\n", r);

   ncols = mpq_QSget_colcount(lpi->prob);
   nrows = mpq_QSget_rowcount(lpi->prob);

   SCIP_CALL(ensureTabMem(lpi, nrows+ncols));

   rval = mpq_QSget_tableau_row(lpi->prob, r, lpi->itab);
   QS_CONDRET(rval);

   /* copy local information to the outside */
   memcpy(coef, lpi->itab, sizeof(mpq_t)*ncols);

   return SCIP_OKAY;
}

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiexGetBInvACol(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   mpq_t*                coef                /**< vector to return coefficients */
   )
{  /*lint --e{715} */
   assert(lpi!=NULL);
   assert(lpi->prob!=NULL);

   SCIPerrorMessage("SCIPlpiexGetBInvACol() not supported by QSopt_ex.\n");

   /* QSopt_ex does not provide an interface for this yet */
   return SCIP_ERROR;
}

/**@} */

/*
 * LP State Methods
 */

/**@name LP State Methods */
/**@{ */

/** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiexGetState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   int ncols, nrows;

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
   SCIPdebugMessage("storing QSopt_ex LPI state in %p (%d cols, %d rows)\n", *lpistate, ncols, nrows);

   /* get unpacked basis information from QSopt_ex */
   SCIP_CALL( ensureColMem(lpi, ncols) );
   SCIP_CALL( ensureRowMem(lpi, nrows) );
   SCIP_CALL( SCIPlpiexGetBase(lpi, lpi->iccnt, lpi->ircnt) );

   /* pack LPi state data */
   (*lpistate)->ncols = ncols;
   (*lpistate)->nrows = nrows;
   lpistatePack(*lpistate, lpi->iccnt, lpi->ircnt);

   return SCIP_OKAY;
}

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiexGetState()
 */
SCIP_RETCODE SCIPlpiexSetState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   )
{  /*lint --e{715} */
   register int i;
   int rval = 0, ncols, nrows;
   char* icstat=0, *irstat=0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   /* if there was no basis information available, LPI state was not stored */
   if (lpistate == NULL)
      return SCIP_OKAY;

   /* continue test */
   ncols = mpq_QSget_colcount(lpi->prob);
   nrows = mpq_QSget_rowcount(lpi->prob);

   assert(ncols >= 0);
   assert(nrows >= 0);
   assert(lpistate->ncols <= ncols);
   assert(lpistate->nrows <= nrows);

   SCIPdebugMessage("loading LPI state %p (%d cols, %d rows) into QSopt_ex LP (%d cols and %d rows)\n", lpistate, lpistate->ncols,
      lpistate->nrows, ncols, nrows);

   if (lpistate->ncols == 0 || lpistate->nrows == 0)
      return SCIP_OKAY;

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureColMem(lpi, ncols) );
   SCIP_CALL( ensureRowMem(lpi, nrows) );
   SCIP_CALL( ensureTabMem(lpi, nrows+ncols) );

   icstat = lpi->ibas;
   irstat = lpi->ibas + ncols;

   /* unpack LPi state data */
   lpistateUnpack(lpistate, lpi->iccnt, lpi->ircnt);

   /* extend the basis to the current LP */
   for (i = lpistate->ncols; i < ncols; ++i)
      lpi->iccnt[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/ /**@todo this has to be corrected for lb = -infinity */
   for (i = lpistate->nrows; i < nrows; ++i)
      lpi->ircnt[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/

   /* convert the loaded basis into QSopt_ex format */
   SCIPdebugMessage("basis status of SCIP lpistate rows (nrows=%d):\n", lpistate->nrows);
   for (i = 0; i < nrows; ++i)
   {
      SCIPdebugMessage("row_%d: %d (%s)\n", i, lpi->ircnt[i],
         lpi->ircnt[i] == SCIP_BASESTAT_LOWER ? "lower" : lpi->ircnt[i] == SCIP_BASESTAT_BASIC ? "basic" : "upper");

      switch(lpi->ircnt[i])
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
   for (i = 0; i < ncols; ++i)
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
SCIP_RETCODE SCIPlpiexFreeState(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   assert(lpi != NULL);
   assert(lpistate != NULL);

   if (*lpistate != NULL)
      lpistateFree(lpistate, blkmem);

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiexHasStateBasis(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{  /*lint --e{715} */
   return (lpistate != NULL);
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiexReadState(
   SCIP_LPIEX*           lpi,               /**< LP interface structure */
   const char*           fname              /**< file name */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("reading QSopt_ex LP state from file <%s>\n", fname);

   rval = mpq_QSread_and_load_basis(lpi->prob, fname);
   if ( rval )
   {
      SCIPerrorMessage("Error while loading basis from file <%s>.\n", fname);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** writes LP state (like basis information) to a file */
SCIP_RETCODE SCIPlpiexWriteState(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   const char*           fname           /**< file name */
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
   if ( rval )
   {
      SCIPerrorMessage("Could not write basis to file <%s>.\n", fname);
      return SCIP_WRITEERROR;
   }

   return SCIP_OKAY;
}

/** checks whether LPi state (i.e. basis information) is dual feasible and returns corresponding dual objective value.
 *  if wanted it will first directly test the corresponding approximate dual and primal solution
 *  (corrected via dual variables for bounds and primal variables for slacks if possible) for optimality
 *  before performing the dual feasibility test on the more expensive exact basic solution.
 */
SCIP_RETCODE SCIPlpiexStateDualFeasible(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
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

   /* loads LPi state (like basis information) into solver */
   SCIP_CALL( SCIPlpiexSetState(lpi, blkmem, lpistate) );

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

/**@} */

/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiexGetIntpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(ival != NULL);

   SCIPdebugMessage("getting int parameter %d\n", type);

   switch (type)
   {
   case SCIP_LPPAR_FROMSCRATCH:
   case SCIP_LPPAR_FASTMIP:
   case SCIP_LPPAR_PRESOLVING:
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_LPPAR_SCALING:
      rval = mpq_QSget_param(lpi->prob, QS_PARAM_SIMPLEX_SCALING, ival);
      if (*ival)
	 *ival = TRUE;
      else
	 *ival = FALSE;
      break;
   case SCIP_LPPAR_PRICING:
      *ival = lpi->pricing;
      break;
   case SCIP_LPPAR_LPINFO:
      rval = mpq_QSget_param(lpi->prob, QS_PARAM_SIMPLEX_DISPLAY, ival);
      if (*ival)
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
SCIP_RETCODE SCIPlpiexSetIntpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
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
      if (ival == TRUE)
	 rval = mpq_QSset_param(lpi->prob, QS_PARAM_SIMPLEX_SCALING, 1);
      else
	 rval = mpq_QSset_param(lpi->prob, QS_PARAM_SIMPLEX_SCALING, 0);
      break;
   case SCIP_LPPAR_PRICING:
      lpi->pricing = ival;
      switch(ival)
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
      if (ival == TRUE)
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
SCIP_RETCODE SCIPlpiexGetRealpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   mpq_t*                dval                /**< buffer to store the parameter value */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(dval != NULL);

   SCIPdebugMessage("getting real parameter %d\n", type);

   switch (type)
   {
   case SCIP_LPPAR_LOBJLIM:
      rval = mpq_QSget_param_EGlpNum(lpi->prob, QS_PARAM_OBJLLIM, dval);
      break;
   case SCIP_LPPAR_UOBJLIM:
      rval = mpq_QSget_param_EGlpNum(lpi->prob, QS_PARAM_OBJULIM, dval);
      break;
   case SCIP_LPPAR_LPTILIM:
      rval = mpq_QSget_param_EGlpNum(lpi->prob, QS_PARAM_SIMPLEX_MAX_TIME, dval);
      break;
   default:
   case SCIP_LPPAR_MARKOWITZ:
   case SCIP_LPPAR_BARRIERCONVTOL:
   case SCIP_LPPAR_DUALFEASTOL:
   case SCIP_LPPAR_FEASTOL:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   QS_RETURN(rval);
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiexSetRealpar(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   mpq_t                 dval                /**< parameter value */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("setting real parameter %d to %g\n", type, mpq_get_d(dval));

   switch (type)
   {
   case SCIP_LPPAR_LPTILIM:
      rval = mpq_QSset_param_EGlpNum(lpi->prob, QS_PARAM_SIMPLEX_MAX_TIME, dval);
      break;
   case SCIP_LPPAR_LOBJLIM:
      rval = mpq_QSset_param_EGlpNum(lpi->prob, QS_PARAM_OBJLLIM, dval);
	 break;
   case SCIP_LPPAR_UOBJLIM:
      rval = mpq_QSset_param_EGlpNum(lpi->prob, QS_PARAM_OBJULIM, dval);
      break;
   case SCIP_LPPAR_FEASTOL:
   case SCIP_LPPAR_DUALFEASTOL:
   case SCIP_LPPAR_BARRIERCONVTOL:
   case SCIP_LPPAR_MARKOWITZ:
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   QS_RETURN(rval);
}

/**@} */

/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as positive infinity in the LP solver */
void SCIPlpiexPosInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   mpq_t*                infval          /**< pointer to store positive infinity value of LP solver */
   )
{
   assert(infval != NULL);

   mpq_set(*infval, mpq_ILL_MAXDOUBLE);
}

/** checks if given value is treated as positive infinity in the LP solver */
SCIP_Bool SCIPlpiexIsPosInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   const mpq_t           val
   )
{  /*lint --e{715} */
   return (mpq_cmp(val, mpq_ILL_MAXDOUBLE) >= 0);
}

/** returns value treated as negative infinity in the LP solver */
void SCIPlpiexNegInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   mpq_t*                infval          /**< pointer to store negative infinity value of LP solver */
   )
{  /*lint --e{715} */
   assert(infval != NULL);

   mpq_set(*infval, mpq_ILL_MINDOUBLE);
}

/** checks if given value is treated as negative infinity in the LP solver */
SCIP_Bool SCIPlpiexIsNegInfinity(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   const mpq_t           val
   )
{  /*lint --e{715} */
   return (mpq_cmp(val, mpq_ILL_MINDOUBLE) <= 0);
}

/**@} */

/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
SCIP_RETCODE SCIPlpiexReadLP(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   int j;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("reading LP from file <%s>\n", fname);

   if (lpi->prob)
      mpq_QSfree_prob(lpi->prob);

   lpi->solstat = 0;
   lpi->previt = 0;

   /* try to extract file type */
   j = strlen(fname)-1;
   while (j >= 0 && fname[j] != '.' )
      --j;
   if (fname[j] == '.')
      ++j;

   /* load problem */
   lpi->prob = mpq_QSread_prob(fname, &(fname[j]));
   if ( lpi->prob == 0 )
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiexWriteLP(
   SCIP_LPIEX*           lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   int j;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("writing LP to file <%s>\n", fname);

   /* try to extract file type */
   j = strlen(fname)-1;
   while (j >= 0 && fname[j] != '.' )
      --j;
   if (fname[j] == '.')
      ++j;

   /* write problem */
   if ( mpq_QSwrite_prob(lpi->prob, fname, &(fname[j])) )
      return SCIP_WRITEERROR;

   return SCIP_OKAY;
}


/** computes and stores matrix factorization within the LPIEX structure */
SCIP_RETCODE SCIPlpiexCreateFactor(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   int                   dim,            /**< dimension of matrix */
   int*                  cbeg,           /**< column indices of matrix */
   int*                  clen,           /**< column lengths of matrix */
   int*                  cindx,          /**< row index of entries */
   mpq_t*                ccoef           /**< coef values of matrix */
   )
{
   int i;
   int rval;
   if(lpi->factor == NULL)
   {
      int nsing;
      int *singr;
      int *singc;
      int * basis;

      SCIP_ALLOC(  BMSallocMemoryArray(&lpi->factor,1) );

      SCIP_ALLOC(  BMSallocMemoryArray(&basis,dim) );
      for(i = 0; i < dim; i++)
      {
         basis[i] = i;
      }
      /* mpq_factor_work *f = (mpq_factor_work *) NULL; */
      /* f = (mpq_factor_work *) malloc (sizeof (mpq_factor_work)); */
      /* this is just a temporary fix for debugging, this should use lpi->factor!!! */
      mpq_init(lpi->factor->fzero_tol);
      mpq_init(lpi->factor->szero_tol);
      mpq_init(lpi->factor->partial_tol);
      mpq_init(lpi->factor->partial_cur);

      mpq_ILLfactor_init_factor_work (lpi->factor );/*  lpi->factor  */
      mpq_ILLfactor_create_factor_work (lpi->factor,dim);

      nsing = 0;
      singr = 0;
      singc = 0;

      rval = mpq_ILLfactor(lpi->factor,basis, cbeg,clen,cindx,ccoef,&nsing,&singr,&singc);
      assert(!rval);
      if (nsing > 0)
      {
         printf ("Matrix is nonsingular \n");
      }
      BMSfreeMemoryArray( &basis );
      if( rval )
         return SCIP_ERROR;
   }
   return SCIP_OKAY;
}

/** solves a system using the stored factorization */
SCIP_RETCODE SCIPlpiexFactorSolve(
   SCIP_LPIEX*           lpi,            /**< LP interface structure */
   int                   dim,            /**< dimension of matrix */
   mpq_t*                sol,            /**< solution to system */
   mpq_t*                rhs             /**< rhs of system */
   )
{
   int i;
   if(lpi->factor != NULL)
   {
      mpq_svector ssol;
      mpq_svector srhs;

      mpq_ILLsvector_init (& ssol);
      mpq_ILLsvector_init (& srhs);
      mpq_ILLsvector_alloc (& ssol,dim);
      mpq_ILLsvector_alloc (& srhs,dim);

      for(i = 0; i < dim; i++)
      {
	 mpq_set(srhs.coef[i],rhs[i]);
         srhs.indx[i] = i;
      }

      /* solve the system */
      mpq_ILLfactor_ftran(lpi->factor, &srhs, &ssol);

      for(i = 0; i < ssol.nzcnt ; i++)
         mpq_set_ui(sol[i],0,1);
      for(i = 0; i < ssol.nzcnt ; i++)
      {
	 mpq_set(sol[ssol.indx[i]],ssol.coef[i]);
      }

      mpq_ILLsvector_free (& ssol);
      mpq_ILLsvector_free (& srhs);
   }
   return SCIP_OKAY;
}

/**@} */
#endif
