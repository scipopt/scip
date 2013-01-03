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

/**@file   lpi_qso.c
 * @brief  LP interface for QSopt version >= 070303
 * @author Daniel Espinoza
 * @author Marc Pfetsch
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "qsopt.h"
#include "scip/bitencode.h"
#include "scip/lpi.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include <string.h>

typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE

/** LP interface */
struct SCIP_LPi
{
   QSprob                prob;               /**< LP struct pointer */
   int                   solstat;            /**< solution status of last optimization call */
   int                   previt;             /**< previous number of simplex iterations performed */
   int                   rowspace;           /**< current size of internal row-related arrays */
   char*                 isen;               /**< array of length rowspace */
   double*               irhs;               /**< array of rhs rowspace */
   double*               irng;               /**< array of range rowspace */
   int*                  ircnt;              /**< array of count rowspace */
   int*                  irbeg;              /**< array of beginning index rowspace */
   int                   colspace;           /**< current size of internal column-related arrays */
   int*                  iccnt;              /**< array of length colspace */
   char*                 iccha;              /**< array of type colspace */
   int                   tbsz;               /**< current size of tableau-related arrays */
   double*               itab;               /**< array of length tbsz */
   char*                 ibas;               /**< array of length tbsz */
   int                   pricing;            /**< SCIP pricing option */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehdlr handler to printing messages, or NULL */
};

/** solver name */
static char __qsstr[1024];



/*
 * local defines
 */

/** print location of the calling line */
#define __QS_PRINTLOC__ fprintf(stderr,", in (%s:%d)\n", __FILE__, __LINE__);

/** This macro is to print error messages and jump to the given point in the code, it also print the
 * file and line where this happened */
#define QS_TESTG(A,B,...) do{{                  \
         if (A){                                \
            fprintf(stderr,__VA_ARGS__);        \
            __QS_PRINTLOC__;                    \
            goto B;}}}while(0)

/** This macro is to print error messages and to exit with SCIP_LPERROR */
#define QS_ERROR(A,...) do{{                    \
         if (A){                                \
            fprintf(stderr,__VA_ARGS__);        \
            __QS_PRINTLOC__;                    \
            return SCIP_LPERROR;}}}while(0)

/** return value macro, if the value is non-zero, write to standard error the returning code and
 * where this happened, and return SCIP_ERROR, otherwise return normal SCIP_OKAY termination code. */
#define QS_RETURN(A) do{                                                \
      const int __RVAL__ = (A);                                         \
      if (__RVAL__){                                                    \
         fprintf(stderr,"LP Error: QSopt returned %d",__RVAL__);        \
         __QS_PRINTLOC__;                                               \
         return SCIP_ERROR;}                                            \
      return SCIP_OKAY;}while(0)

/** return value macro, if the value is non-zero, write to standard error the returning code and
 * where this happened, and return SCIP_ERROR, otherwise do nothing. */
#define QS_CONDRET(A) do{                                               \
      const int __RVAL__ = (A);                                         \
      if (__RVAL__){                                                    \
         fprintf(stderr,"LP Error: QSopt returned %d",__RVAL__);        \
         __QS_PRINTLOC__;                                               \
         return SCIP_LPERROR;}                                          \
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
static
SCIP_RETCODE ensureTabMem(
   SCIP_LPI* const lpi,
   int const sz
   )
{
   if( lpi->tbsz < sz )
   {
      lpi->tbsz = sz*2;
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->itab), lpi->tbsz) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->ibas), lpi->tbsz) );
   }
   return SCIP_OKAY;
}

/** ensure size of column-related arrays */
static
SCIP_RETCODE ensureColMem(
   SCIP_LPI* const lpi,
   int const ncols
   )
{
   if( lpi->colspace < ncols )
   {
      lpi->colspace = ncols*2;
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->iccnt), lpi->colspace) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->iccha), lpi->colspace) );
   }
   return SCIP_OKAY;
}

/** ensure size of row-related arrays */
static
SCIP_RETCODE ensureRowMem(
   SCIP_LPI* const lpi,
   int const nrows
   )
{
   if( lpi->rowspace < nrows )
   {
      lpi->rowspace = nrows*2;
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->isen), lpi->rowspace) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->irhs), lpi->rowspace) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->irng), lpi->rowspace) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->ircnt), lpi->rowspace) );
      SCIP_ALLOC( BMSreallocMemoryArray(&(lpi->irbeg), lpi->rowspace) );
   }
   return SCIP_OKAY;
}

/** transform lhs/rhs into qsopt format */
static
SCIP_RETCODE convertSides(
   SCIP_LPI* const lpi,
   int const nrows,
   const double* const lhs,
   const double* const rhs
   )
{
   int state;
   register int i;

   for( i = nrows ; i-- ; )
   {
      state = ((lhs[i] <= -QS_MAXDOUBLE ? 1U:0U) | (rhs[i] >= QS_MAXDOUBLE ? 2U:0U));
      lpi->ircnt[i] = 0;
      lpi->irbeg[i] = 0;
      switch( state )
      {
      case 0:
        /* check for equations */
         if( lhs[i] == rhs[i] )
         {
            lpi->isen[i] = 'E';
            lpi->irhs[i] = lhs[i];
            lpi->irng[i] = 0.0;
         }
         else
         {
            lpi->isen[i] = 'R';
            lpi->irhs[i] = lhs[i];
            lpi->irng[i] = rhs[i] - lhs[i];
            assert(lpi->irng[i] >= 0.0);
         }
         break;
      case 1:
         lpi->isen[i] = 'L';
         lpi->irhs[i] = rhs[i];
         lpi->irng[i] = 0;
         break;
      case 2:
         lpi->isen[i] = 'G';
         lpi->irhs[i] = lhs[i];
         lpi->irng[i] = 0;
         break;
      default:
         SCIPerrorMessage("Error, constraint %d has no bounds!",i);
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
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
const char* SCIPlpiGetSolverName(void)
{
   char* vname = QSversion();
   size_t vnamelen;
   vnamelen = strlen(vname);
   memcpy(__qsstr, vname, MIN(sizeof(__qsstr), vnamelen+1));
   __qsstr[sizeof(__qsstr)-1] = '\0';
   QSfree(vname);
   return __qsstr;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return "Linear Programming Solver developed by D. Applegate, W. Cook, S. Dash, and M. Mevenkamp (www.isye.gatech.edu/~wcook/qsopt)";
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
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

/** creates an LP problem object
 * @return SCIP_OK on success
 * */
SCIP_RETCODE SCIPlpiCreate(
   SCIP_LPI**            lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   /* QSopt only works with doubles as floating points and bool as integers */
   assert(sizeof (SCIP_Real) == sizeof (double));
   assert(sizeof (SCIP_Bool) == sizeof (int));
   assert(lpi != NULL);

   SCIPdebugMessage("SCIPlpiCreate()\n");

   /* create LP */
   SCIP_ALLOC( BMSallocMemory(lpi) );
   memset(*lpi, 0, sizeof(struct SCIP_LPi));

   (*lpi)->prob = QScreate_prob(name, (int) objsen);
   if ( (*lpi)->prob == NULL )
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

   (*lpi)->messagehdlr = messagehdlr;

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
   QSfree_prob((*lpi)->prob);
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
   register int i;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("loading LP in column format into QSopt: %d cols, %d rows\n", ncols, nrows);

   /* delete old LP */
   SCIP_CALL( SCIPlpiClear(lpi) );

   /* set sense */
   if( objsen == SCIP_OBJSEN_MAXIMIZE )
   {
      rval = QSchange_objsense(lpi->prob, QS_MAX);
      QS_CONDRET(rval);
   }
   else
   {
      rval = QSchange_objsense(lpi->prob, QS_MIN);
      QS_CONDRET(rval);
   }

   /* add rows with no matrix, and then the columns, first ensure space */
   SCIP_CALL( ensureRowMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   SCIP_CALL( convertSides(lpi, nrows, lhs, rhs) );

   /* now we add the rows */
   rval = QSadd_ranged_rows(lpi->prob, nrows, lpi->ircnt, lpi->irbeg, 0, 0, lpi->irhs, lpi->isen, lpi->irng, (const char**)rownames);
   QS_CONDRET(rval);

   /* ensure column size */
   SCIP_CALL( ensureColMem(lpi, ncols) );

   /* compute column lengths */
   for( i = 0; i < ncols-1; ++i )
   {
      lpi->iccnt[i] = beg[i+1] - beg[i];
      assert(lpi->iccnt[i] >= 0);
   }
   if( ncols > 0 )
   {
      lpi->iccnt[ncols-1] = nnonz - beg[ncols-1];
      assert(lpi->iccnt[ncols-1] >= 0);
   }

   /* and add the columns */
   rval = QSadd_cols(lpi->prob, ncols, lpi->iccnt, (int*) beg, (int*) ind, (SCIP_Real*) val, (SCIP_Real*) obj,
      (SCIP_Real*) lb, (SCIP_Real*) ub, (const char**)colnames);

   QS_RETURN(rval);
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
   int rval = 0;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("adding %d columns with %d nonzeros to QSopt\n", ncols, nnonz);

   lpi->solstat = 0;

   /* ensure column size */
   SCIP_CALL(ensureColMem(lpi, ncols));

   /* compute column lengths */
   for( i = 0; i < ncols - 1; ++i )
   {
      lpi->iccnt[i] = beg[i+1] - beg[i];
      assert(lpi->iccnt[i] >= 0);
   }
   if( ncols > 0 )
   {
      lpi->iccnt[ncols-1] = nnonz - beg[ncols-1];
      assert(lpi->iccnt[ncols-1] >= 0);
   }

   /* and add the columns */
   rval = QSadd_cols(lpi->prob, ncols, lpi->iccnt, (int*) beg, (int*) ind, (SCIP_Real*) val, (SCIP_Real*) obj,
      (SCIP_Real*) lb, (SCIP_Real*) ub, (const char**)colnames);

   QS_RETURN(rval);
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
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
   assert(0 <= firstcol && len > 0 && lastcol < QSget_colcount(lpi->prob));

   SCIPdebugMessage("deleting %d columns from QSopt\n", len);

   SCIP_CALL(ensureColMem(lpi, len));
   for( i = firstcol ; i <= lastcol ; i++ )
      lpi->iccnt[i-firstcol] = i;

   rval = QSdelete_cols(lpi->prob, len, lpi->iccnt);

   QS_RETURN(rval);
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   int rval = 0, ncols, ccnt;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   ncols = QSget_colcount(lpi->prob);
   lpi->solstat = 0;

   SCIPdebugMessage("deleting a column set from QSopt\n");

   rval = QSdelete_setcols(lpi->prob,dstat);
   QS_CONDRET(rval);

   for( i=0, ccnt=0; i < ncols; i++ )
   {
      if( dstat[i] )
         dstat[i] = -1;
      else
         dstat[i] = ccnt++;
   }
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
   register int i;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("adding %d rows with %d nonzeros to QSopt\n", nrows, nnonz);

   /* add rows with no matrix, and then the columns, first ensure space */
   SCIP_CALL( ensureRowMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   SCIP_CALL( convertSides(lpi, nrows, lhs, rhs) );

   /* compute row count */
   if( nnonz > 0 )
   {
      assert(beg != NULL);
      assert(ind != NULL);
      assert(val != NULL);

      for( i = 0 ; i < nrows -1 ; i++ )
      {
         lpi->ircnt[i] = beg[i+1] - beg[i];
         assert(lpi->ircnt[i] >= 0);
      }
      if( nrows > 0 )
      {
         lpi->ircnt[nrows-1] = nnonz - beg[nrows-1];
         assert(lpi->ircnt[nrows-1] >= 0);
      }

      /* now we add the rows */
      rval = QSadd_ranged_rows(lpi->prob, nrows, lpi->ircnt, (int*) beg, (int*) ind, (SCIP_Real*) val, lpi->irhs,
         lpi->isen, lpi->irng, (const char**)rownames);
      QS_ERROR(rval, "failed adding %d rows with %d non-zeros", nrows, nnonz);
   }
   else
   {
      for( i = 0; i < nrows -1; ++i )
      {
         lpi->ircnt[i] = 0;
         lpi->irbeg[i] = 0;
      }

      /* now we add the rows */
      rval = QSadd_ranged_rows(lpi->prob, nrows, lpi->ircnt, lpi->irbeg, (int*) ind, (SCIP_Real*) val, lpi->irhs,
         lpi->isen, lpi->irng, (const char**)rownames);
      QS_ERROR(rval, "failed adding %d rows with %d non-zeros", nrows, nnonz);
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
   char** cnames;
   char* s;
   int ncols;
   int rval;
   int j;
   int sizeleft;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(colnames != NULL || namestoragesize == 0);
   assert(namestorage != NULL || namestoragesize == 0);
   assert(namestoragesize >= 0);
   assert(storageleft != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < QSget_colcount(lpi->prob));

   SCIPdebugMessage("getting column names %d to %d\n", firstcol, lastcol);

   ncols = QSget_colcount(lpi->prob);
   SCIP_ALLOC( BMSallocMemoryArray(&cnames, ncols) );

   rval = QSget_colnames(lpi->prob, cnames);
   QS_ERROR(rval, "failed getting column names");

   /* copy column names */
   s = namestorage;
   sizeleft = namestoragesize;
   for( j = firstcol; j <= lastcol; ++j )
   {
      const char* t;
      t = cnames[j];
      if( colnames != NULL )
         colnames[j-firstcol] = s;
      while( *t != '\0' )
      {
         if( sizeleft > 0 )
            *(s++) = *(t++);
         --sizeleft;
      }
      *(s++) = '\0';
   }
   *storageleft = sizeleft;

   /* free space */
   for( j = 0; j < ncols; ++j )
      free(cnames[j]);

   return SCIP_OKAY;
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
   char** rnames;
   char* s;
   int nrows;
   int rval;
   int i;
   int sizeleft;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(rownames != NULL || namestoragesize == 0);
   assert(namestorage != NULL || namestoragesize == 0);
   assert(namestoragesize >= 0);
   assert(storageleft != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < QSget_rowcount(lpi->prob));

   SCIPdebugMessage("getting row names %d to %d\n", firstrow, lastrow);

   nrows = QSget_rowcount(lpi->prob);
   SCIP_ALLOC( BMSallocMemoryArray(&rnames, nrows) );

   rval = QSget_rownames(lpi->prob, rnames);
   QS_ERROR(rval, "failed getting row names");

   s = namestorage;
   sizeleft = namestoragesize;
   for( i = firstrow; i <= lastrow; ++i )
   {
      const char* t;
      t = rnames[i];
      if( rownames != NULL )
         rownames[i-firstrow] = s;
      while( *t != '\0' )
      {
         if( sizeleft > 0 )
            *(s++) = *(t++);
         --sizeleft;
      }
      *(s++) = '\0';
   }
   *storageleft = sizeleft;

   /* free space */
   for( i = 0; i < nrows; ++i )
      free(rnames[i]);

   return SCIP_OKAY;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
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
   assert(0 <= firstrow && len > 0 && lastrow < QSget_rowcount (lpi->prob));

   SCIPdebugMessage("deleting %d rows from QSopt\n", len);

   SCIP_CALL( ensureRowMem(lpi, len) );
   for( i = firstrow; i <= lastrow; i++ )
      lpi->ircnt[i-firstrow] = i;
   rval = QSdelete_rows(lpi->prob, len, lpi->ircnt);

   QS_RETURN(rval);
}


/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelRowset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   int rval = 0, nrows, ccnt, ndel=0;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   nrows = QSget_rowcount(lpi->prob);
   lpi->solstat = 0;

   for( i = 0; i < nrows; ++i )
   {
      if( dstat[i] == 1 )
         ndel++;
   }

   SCIPdebugMessage("deleting a row set from QSopt (%d)\n",ndel);

   rval = QSdelete_setrows(lpi->prob,dstat);
   QS_CONDRET(rval);

   for( i=0, ccnt=0; i < nrows; i++ )
   {
      if( dstat[i] )
         dstat[i] = -1;
      else
         dstat[i] = ccnt++;
   }
   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   register int i;
   int ncols, nrows, rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("clearing QSopt LP\n");
   lpi->solstat = 0;

   ncols = QSget_colcount(lpi->prob);
   nrows = QSget_rowcount(lpi->prob);
   if( ncols >= 1 )
   {
      SCIP_CALL( ensureColMem(lpi,ncols) );
      for( i = 0; i < ncols; ++i )
         lpi->iccnt[i] = i;
      rval = QSdelete_cols(lpi->prob, ncols, lpi->iccnt);
      QS_CONDRET(rval);
   }

   if( nrows >= 1 )
   {
      SCIP_CALL( ensureRowMem(lpi, nrows) );
      for( i = 0; i < nrows; ++i )
         lpi->ircnt[i] = i;
      rval = QSdelete_rows(lpi->prob, nrows, lpi->ircnt);
      QS_CONDRET(rval);
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
   register int i;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;

   SCIPdebugMessage("changing %d bounds in QSopt\n", ncols);
#ifdef SCIP_DEBUG
   {
      int j;
      for( j = 0; j < ncols; ++j )
         SCIPdebugPrintf("  col %d: [%lg,%lg]\n", ind[j], lb[j], ub[j]);
   }
#endif

   SCIP_CALL(ensureColMem(lpi, ncols));
   for( i = 0; i < ncols; ++i )
      lpi->iccha[i] = 'L';

   rval = QSchange_bounds(lpi->prob, ncols, (int*) ind, lpi->iccha, (SCIP_Real*) lb);
   QS_CONDRET(rval);

   for( i = 0; i < ncols; ++i )
      lpi->iccha[i] = 'U';

   rval = QSchange_bounds(lpi->prob, ncols, (int*) ind, lpi->iccha, (SCIP_Real*) ub);

   QS_RETURN(rval);
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
   register int i;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("changing %d sides in QSopt\n", nrows);

   SCIP_CALL( ensureRowMem(lpi, nrows) );

   /* convert lhs/rhs into sen/rhs/range tuples */
   SCIP_CALL( convertSides(lpi, nrows, lhs, rhs) );

   /* now we change all rows */
   for( i = 0; i < nrows; ++i )
   {
      rval = QSchange_sense(lpi->prob, ind[i], lpi->isen[i]);
      QS_CONDRET(rval);

      rval = QSchange_rhscoef(lpi->prob, ind[i], lpi->irhs[i]);
      QS_CONDRET(rval);

      if( lpi->isen[i] == 'R' )
      {
         rval = QSchange_range(lpi->prob, ind[i], lpi->irng[i]);
         QS_CONDRET(rval);
      }
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
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;

   SCIPdebugMessage("changing coefficient row %d, column %d in QSopt to %g\n", row, col, newval);

   rval = QSchange_coef(lpi->prob, row, col, newval);

   QS_RETURN(rval);
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("changing objective sense in QSopt to %d\n", objsen);

   /* set sense */
   if( objsen == SCIP_OBJSEN_MAXIMIZE )
   {
      rval = QSchange_objsense(lpi->prob, QS_MAX);
      QS_CONDRET(rval);
   }
   else
   {
      rval = QSchange_objsense(lpi->prob, QS_MIN);
      QS_CONDRET(rval);
   }
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
   register int i;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("changing %d objective values in QSopt\n", ncols);

   for( i = 0; i < ncols; ++i )
   {
      rval = QSchange_objcoef(lpi->prob, ind[i], obj[i]);
      QS_CONDRET(rval);
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
   register int i;
   int rowlist[1];
   int* rowcnt = NULL, *rowbeg = NULL, *rowind = NULL;
   double* rowval = NULL, *rhs = NULL, *range = NULL;
   char* sense = NULL;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("scaling row %d with factor %g in QSopt\n", row, scaleval);

   rowlist[0] = row;
   /* get row */
   rval = QSget_ranged_rows_list(lpi->prob, 1, rowlist, &rowcnt, &rowbeg, &rowind, &rowval, &rhs, &sense, &range, 0);
   QS_TESTG(rval, CLEANUP, " ");

   /* change all coefficients in the constraint */
   for( i = 0; i < rowcnt[0]; ++i )
   {
      rval = QSchange_coef(lpi->prob, row, rowind[i], rowval[i] * scaleval);
      QS_TESTG(rval, CLEANUP, " ");
   }

   /* if we have a positive scalar, we just scale rhs and range */
   if( scaleval >= 0 )
   {
      rval = QSchange_rhscoef(lpi->prob, row, rhs[0] * scaleval);
      QS_TESTG(rval, CLEANUP, " ");
      if( sense[0] == 'R' )
      {
         rval = QSchange_range(lpi->prob, row, range[0] * scaleval);
         QS_TESTG(rval, CLEANUP, " ");
      }
   }
   /* otherwise, we must change everything */
   else
   {
      switch( sense[0] )
      {
      case 'E':
         rval = QSchange_rhscoef(lpi->prob, row, rhs[0]*scaleval);
         QS_TESTG(rval, CLEANUP, " ");
         break;
      case 'L':
         rval = QSchange_rhscoef(lpi->prob, row, rhs[0]*scaleval);
         QS_TESTG(rval, CLEANUP, " ");
         rval = QSchange_sense(lpi->prob, row, 'G');
         QS_TESTG(rval, CLEANUP, " ");
         break;
      case 'G':
         rval = QSchange_rhscoef(lpi->prob, row, rhs[0]*scaleval);
         QS_TESTG(rval, CLEANUP, " ");
         rval = QSchange_sense(lpi->prob, row, 'L');
         QS_TESTG(rval, CLEANUP, " ");
         break;
      case 'R':
         rhs[0] = (rhs[0] + range[0]) * scaleval;
         range[0] = fabs(scaleval) * range[0];
         rval = QSchange_rhscoef(lpi->prob, row, rhs[0]);
         QS_TESTG(rval, CLEANUP, " ");
         rval = QSchange_range(lpi->prob, row, range[0]);
         QS_TESTG(rval, CLEANUP, " ");
         break;
      default:
         SCIPerrorMessage("Impossible! received sense %c (not E L G R)", sense[0]);
         rval = 1;
         goto CLEANUP;
      }
   }

   /* now we must free all received arrays */
   /* ending */
 CLEANUP:
   if( rowcnt != NULL )
      QSfree(rowcnt);
   if( rowbeg != NULL ) 
      QSfree(rowbeg);
   if( rowind != NULL ) 
      QSfree(rowind);
   if( rowval != NULL ) 
      QSfree(rowval);
   if( rhs != NULL ) 
      QSfree(rhs);
   if( sense != NULL ) 
      QSfree(sense);
   if( range != NULL ) 
      QSfree(range);

   QS_RETURN(rval);
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
   register int i;
   int collist[1];
   int* colcnt=0;
   int* colbeg=0;
   int* colind=0;
   double* colval=0;
   double* lb=0;
   double* ub=0;
   double* obj=0;

   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   lpi->solstat = 0;
   SCIPdebugMessage("scaling column %d with factor %g in QSopt\n", col, scaleval);

   /* get the column */
   collist[0] = col;
   rval = QSget_columns_list(lpi->prob, 1, collist, &colcnt, &colbeg, &colind, &colval, &obj, &lb, &ub, 0);
   QS_TESTG(rval,CLEANUP," ");

   /* scale column coefficients */
   for( i = 0; i < colcnt[0]; ++i )
   {
      rval = QSchange_coef(lpi->prob, colind[i], col, colval[i]*scaleval);
      QS_TESTG(rval,CLEANUP," ");
   }

   /* scale objective value */
   rval = QSchange_objcoef(lpi->prob, col, obj[0]*scaleval);
   QS_TESTG(rval,CLEANUP," ");

   /* scale column bounds */
   if( scaleval < 0 )
   {
      scaleval = -scaleval;
      obj[0] = lb[0];
      lb[0] = -ub[0];
      ub[0] = -obj[0];
   }
   if( lb[0] > -QS_MAXDOUBLE )
      lb[0] *= scaleval;
   if( ub[0] < QS_MAXDOUBLE )
      ub[0] *= scaleval;

   if( lb[0] < -QS_MAXDOUBLE )
      lb[0] = -QS_MAXDOUBLE;
   if( ub[0] > QS_MAXDOUBLE )
      ub[0] = QS_MAXDOUBLE;

   rval = QSchange_bound(lpi->prob, col, 'L', lb[0]);
   QS_TESTG(rval,CLEANUP," ");
   rval = QSchange_bound(lpi->prob, col, 'U', ub[0]);
   QS_TESTG(rval,CLEANUP," ");

   /* ending */
 CLEANUP:
   if( colcnt != NULL )
      QSfree(colcnt);
   if( colbeg != NULL ) 
      QSfree(colbeg);
   if( colind != NULL )
      QSfree(colind);
   if( colval != NULL )
      QSfree(colval);
   if( obj != NULL )
      QSfree(obj);
   if( lb != NULL )
      QSfree(lb);
   if( ub != NULL )
      QSfree(ub);

   QS_RETURN(rval);
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
   assert(lpi->prob != NULL);
   assert(nrows != NULL);

   SCIPdebugMessage("getting number of rows\n");

   *nrows = QSget_rowcount(lpi->prob);

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(ncols != NULL);

   SCIPdebugMessage("getting number of columns\n");

   *ncols = QSget_colcount(lpi->prob);

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting number of columns\n");

   *nnonz = QSget_nzcount(lpi->prob);

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
   int len;
   register int i;
   double* lval = NULL;
   double* llb = NULL;
   double* lub = NULL;
   int rval = 0;
   int* lcnt = NULL;
   int* lbeg = NULL;
   int* lind = NULL;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < QSget_colcount (lpi->prob));
   assert((lb == 0 && ub == 0) || (lb != 0 && ub != 0));
   assert((nnonz != 0 && beg != 0 && ind != 0 && val != 0) || (nnonz == 0 && beg == 0 && ind == 0 && val == 0));

   SCIPdebugMessage("getting columns %d to %d\n", firstcol, lastcol);

   /* build col-list */
   len = lastcol - firstcol + 1;
   SCIP_CALL( ensureColMem(lpi,len) );
   for( i = 0; i < len; ++i )
      lpi->iccnt[i] = i + firstcol;

   /* get data from qsopt */
   rval = QSget_columns_list(lpi->prob, len, lpi->iccnt, nnonz ? (&lcnt) : 0, nnonz ? (&lbeg) : 0, nnonz ? (&lind) : 0,
      nnonz ? (&lval) : 0, 0, lb ? (&llb) : 0, lb ? (&lub) : 0, 0);

   QS_TESTG(rval, CLEANUP, " ");

   /* store in the user-provided data */
   if( nnonz )
   {
      assert(lbeg != NULL);
      assert(lcnt != NULL);
      assert(lind != NULL);
      assert(lval != NULL);

      *nnonz = lbeg[len-1] + lcnt[len-1];
      for( i = 0 ; i < len ; i++ )
         beg[i] = lbeg[i];  /*lint !e613*/
      for( i = 0; i < *nnonz; ++i )
      {
         ind[i] = lind[i];  /*lint !e613*/
         val[i] = lval[i];  /*lint !e613*/
      }
   }
   if( lb )
   {
      assert(llb != NULL);
      assert(lub != NULL);

      for( i = 0; i < len; ++i )
      {
         lb[i] = llb[i];
         ub[i] = lub[i];   /*lint !e613*/
      }
   }

 CLEANUP:
   if( lval != NULL )
      QSfree(lval);
   if( lub != NULL )
      QSfree(lub);
   if( llb != NULL )
      QSfree(llb);
   if( lind != NULL )
      QSfree(lind);
   if( lbeg != NULL )
      QSfree(lbeg);
   if( lcnt != NULL )
      QSfree(lcnt);

   QS_RETURN(rval);
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
   const int len = lastrow - firstrow + 1;
   register int i;
   double* lval = NULL;
   double* lrhs = NULL;
   double* lrng = NULL;
   int rval = 0;
   int* lcnt = NULL;
   int* lbeg = NULL;
   int* lind = NULL;
   char* lsense = NULL;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < QSget_rowcount (lpi->prob));
   assert((lhs == 0 && rhs == 0) || (rhs != 0 && lhs != 0));
   assert((nnonz != 0 && beg != 0 && ind != 0 && val != 0) || (nnonz == 0 && beg == 0 && ind == 0 && val == 0));

   SCIPdebugMessage("getting rows %d to %d\n", firstrow, lastrow);

   /* build row-list */
   SCIP_CALL( ensureRowMem(lpi, len) );
   for( i = 0; i < len; ++i )
      lpi->ircnt[i] = i + firstrow;

   /* get data from qsopt */
   rval = QSget_ranged_rows_list(lpi->prob, len, lpi->ircnt, nnonz ? (&lcnt) : 0, nnonz ? (&lbeg) : 0, nnonz ? (&lind) : 0,
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
      for( i = 0 ; i < len; i++ )
         beg[i] = lbeg[i];  /*lint !e613*/
      for( i = 0; i < *nnonz; ++i )
      {
         ind[i] = lind[i];  /*lint !e613*/
         val[i] = lval[i];  /*lint !e613*/
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
            lhs[i] = lrhs[i];            /*lint !e613*/
            rhs[i] = lrhs[i] + lrng[i];  /*lint !e613*/
            break;
         case 'E':
            lhs[i] = rhs[i] = lrhs[i];   /*lint !e613*/
            break;
         case 'L':
            rhs[i] = lrhs[i];            /*lint !e613*/
            lhs[i] = -QS_MAXDOUBLE;      /*lint !e613*/
            break;
         case 'G':
            lhs[i] = lrhs[i];            /*lint !e613*/
            rhs[i] = QS_MAXDOUBLE;       /*lint !e613*/
            break;
         default:
            SCIPerrorMessage("Unknown sense %c from QSopt", lsense[i]);
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
         }
      }
   }

 CLEANUP:
   if( lsense != NULL )
      QSfree(lsense);
   if( lrng != NULL )
      QSfree(lrng);
   if( lrhs != NULL )
      QSfree(lrhs);
   if( lval != NULL )
      QSfree(lval);
   if( lind != NULL )
      QSfree(lind);
   if( lbeg != NULL )
      QSfree(lbeg);
   if( lcnt != NULL )
      QSfree(lcnt);

   QS_RETURN(rval);
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
   const int len = lastcol - firstcol + 1;
   int rval = 0;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < QSget_colcount (lpi->prob));

   SCIPdebugMessage("getting objective values %d to %d\n", firstcol, lastcol);

   /* build col-list */
   SCIP_CALL(ensureColMem(lpi,len));
   for( i = 0; i < len; ++i )
      lpi->iccnt[i] = i + firstcol;

   /* get data from qsopt */
   rval = QSget_obj_list(lpi->prob, len, lpi->iccnt, vals);

   QS_RETURN(rval);
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
   const int len = lastcol - firstcol + 1;
   int rval = 0;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(0 <= firstcol && firstcol <= lastcol&& lastcol < QSget_colcount (lpi->prob));

   SCIPdebugMessage("getting bound values %d to %d\n", firstcol, lastcol);

   /* build col-list */
   SCIP_CALL(ensureColMem(lpi,len));
   for( i = 0; i < len; ++i )
      lpi->iccnt[i] = i + firstcol;

   /* get data from qsopt */
   rval = QSget_bounds_list(lpi->prob, len, lpi->iccnt, lbs, ubs);

   QS_RETURN(rval);
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
   const int len = lastrow - firstrow + 1;
   register int i;
   double* lrhs=0, *lrng=0;
   int rval = 0;
   char* lsense=0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < QSget_rowcount (lpi->prob));
   assert(rhss != NULL);
   assert(lhss != NULL);

   SCIPdebugMessage("getting row sides %d to %d\n", firstrow, lastrow);

   /* build row-list */
   SCIP_CALL( ensureRowMem(lpi, len) );
   for( i = 0; i < len; ++i )
      lpi->ircnt[i] = i + firstrow;

   /* get data from qsopt */
   rval = QSget_ranged_rows_list(lpi->prob, len, lpi->ircnt, 0, 0, 0, 0, &lrhs, &lsense, &lrng, 0);
   QS_TESTG(rval, CLEANUP, " ");

   /* store in the user-provided data */
   for( i = 0; i < len; ++i )
   {
      switch( lsense[i] )
      {
      case 'R':
         lhss[i] = lrhs[i];
         rhss[i] = lrhs[i] + lrng[i];
         break;
      case 'E':
         lhss[i] = rhss[i] = lrhs[i];
         break;
      case 'L':
         rhss[i] = lrhs[i];
         lhss[i] = -QS_MAXDOUBLE;
         break;
      case 'G':
         lhss[i] = lrhs[i];
         rhss[i] = QS_MAXDOUBLE;
         break;
      default:
         SCIPerrorMessage("Unknown sense %c from QSopt", lsense[i]);
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

 CLEANUP:
   if( lsense != NULL )
      QSfree(lsense);
   if( lrng != NULL )
      QSfree(lrng);
   if( lrhs != NULL )
      QSfree(lrhs);

   QS_RETURN(rval);
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiGetCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting coefficient of row %d col %d\n", row, col);

   rval = QSget_coef(lpi->prob, row, col, val);

   QS_RETURN(rval);
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
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("calling QSopt primal simplex: %d cols, %d rows, %d nz\n", QSget_colcount(lpi->prob),
      QSget_rowcount(lpi->prob), QSget_nzcount(lpi->prob));

   rval = QSopt_primal(lpi->prob, &(lpi->solstat));

   QS_RETURN(rval);
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("calling QSopt dual simplex: %d cols, %d rows, %d nz\n", QSget_colcount(lpi->prob),
      QSget_rowcount(lpi->prob), QSget_nzcount(lpi->prob));

   rval = QSopt_dual(lpi->prob, &(lpi->solstat));

   QS_RETURN(rval);
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{  /*lint --e{715}*/
   return SCIPlpiSolveDual(lpi);
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
   int rval = 0;
   int nit;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   SCIPdebugMessage("calling QSopt strong branching on variable %d with fractional value (%d it lim)\n", col, itlim);

   /* results of QSopt are valid in any case */
   *downvalid = TRUE;
   *upvalid = TRUE;

   assert(!EPSISINT(psol, 1e-06));

   /* call QSopt */
   rval = QSopt_strongbranch(lpi->prob, 1, &col, &psol, down, up, itlim, QS_MAXDOUBLE);
   QS_CONDRET(rval);

   rval = QSget_itcnt(lpi->prob, 0, 0, 0, 0, &nit);
   QS_CONDRET(rval);

   if( iter )
      *iter = nit - lpi->previt;
   lpi->previt = nit;

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
   int rval = 0;
   int nit;
   int j;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(cols != NULL);
   assert(psols != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   SCIPdebugMessage("calling QSopt strong branching on %d variables with fractional value (%d it lim)\n", ncols, itlim);

   /* results of QSopt are valid in any case */
   for( j = 0; j < ncols; ++j )
   {
      downvalid[j] = TRUE;
      upvalid[j] = TRUE;
      assert(!EPSISINT(psols[j], 1e-06));
   }

   /* call QSopt */
   rval = QSopt_strongbranch(lpi->prob, ncols, cols, psols, down, up, itlim, QS_MAXDOUBLE);
   QS_CONDRET(rval);

   rval = QSget_itcnt(lpi->prob, 0, 0, 0, 0, &nit);
   QS_CONDRET(rval);

   if( iter )
      *iter = nit - lpi->previt;
   lpi->previt = nit;

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
   int rval = 0;
   SCIP_Real objval;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   SCIPdebugMessage("calling QSopt strong branching on variable %d with integral value (%d it lim)\n", col, itlim);

   assert(EPSISINT(psol, 1e-06));

   /* QSopt cannot directly strong branch on integral values! We thus return the current objective
    * value for both cases. Could also implement a manual search as in lpi_cpx.c
    */
   rval = QSget_objval(lpi->prob, &objval);
   QS_CONDRET(rval);

   *down = objval;
   *up = objval;
   *downvalid = TRUE;
   *upvalid = TRUE;
   
   if( iter )
      *iter = 0;

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
   int rval = 0;
   SCIP_Real objval;
   int j;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   SCIPdebugMessage("calling QSopt strong branching on %d variables with integral value (%d it lim)\n", ncols, itlim);

   /* QSopt cannot directly strong branch on integral values! We thus return the current objective
    * value for all cases. Could also implement a manual search as in lpi_cpx.c
    */
   rval = QSget_objval(lpi->prob, &objval);
   QS_CONDRET(rval);

   for( j = 0; j < ncols; ++j )
   {
      assert(EPSISINT(psols[j], 1e-06));
      down[j] = objval;
      up[j] = objval;
      downvalid[j] = TRUE;
      upvalid[j] = TRUE;
   }

   if( iter )
      *iter = 0;

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
   assert(lpi!=NULL);
   assert(lpi->prob!=NULL);

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat != 0 && lpi->solstat != QS_LP_MODIFIED);
}

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_RETCODE SCIPlpiGetSolFeasibility(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting solution feasibility\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   if( lpi->solstat == QS_LP_OPTIMAL || lpi->solstat == QS_LP_UNBOUNDED )
      *primalfeasible = 1;

   if( lpi->solstat == QS_LP_OPTIMAL || lpi->solstat == QS_LP_INFEASIBLE || lpi->solstat == QS_LP_OBJ_LIMIT )
      *dualfeasible = 1;

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
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking primal ray existence\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for primal ray\n");

   /* the current version of QSopt cannot give a primal certificate of unboundedness */
   return FALSE;
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for primal unboundedness\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for primal infeasibility\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_INFEASIBLE);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for primal feasibility\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_OPTIMAL || lpi->solstat == QS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual ray availability\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_INFEASIBLE);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual ray availability\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_INFEASIBLE);
}

/** returns TRUE iff LP is dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual unboundedness\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_INFEASIBLE);
}

/** returns TRUE iff LP is dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual infeasibility\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_UNBOUNDED);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for dual feasibility\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_OPTIMAL || lpi->solstat == QS_LP_OBJ_LIMIT);
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for optimality\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_OPTIMAL);
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for numerical stability\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat != QS_LP_NUMERR);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for objective limit exceeded\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_OBJ_LIMIT);
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for iteration limit exceeded\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_ITER_LIMIT);
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("checking for time limit exceeded\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return (lpi->solstat == QS_LP_TIME_LIMIT);
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting internal solution status\n");

   (void) QSget_status(lpi->prob, &(lpi->solstat));

   return lpi->solstat;
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("ignore instability (will fail)\n");

   /* it seems that in QSopt this does not make much sense */
   *success = FALSE;

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting solution's objective value\n");

   rval = QSget_objval(lpi->prob, objval);

   QS_RETURN(rval);
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
   int rval = 0, nrows;
   register int i;

#ifdef SCIP_DEBUG
   int stat, ncols, sense;
   char *icstat, *irstat;
#endif

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting solution\n");

   nrows = QSget_rowcount(lpi->prob);
   SCIP_CALL( ensureRowMem(lpi, nrows) );

   rval = QSget_solution(lpi->prob, objval, primsol, dualsol, lpi->irng, redcost);
   QS_CONDRET(rval);

#if 0
#ifdef SCIP_DEBUG
   QSget_status(lpi->prob, &stat);
   rval = QSget_objsense(lpi->prob, &sense);
   if( stat == QS_LP_OPTIMAL )
   {
      ncols = QSget_colcount(lpi->prob);
      QS_CONDRET(rval);

      SCIP_CALL(ensureTabMem(lpi,nrows+ncols));
      icstat = lpi->ibas;
      irstat = lpi->ibas+ncols;

      rval = QSget_basis_array(lpi->prob,icstat, irstat);
      QS_CONDRET(rval);

      for( i = ncols ; i-- ; )
      {
         switch( icstat[i] )
         {
         case QS_COL_BSTAT_BASIC:
         case QS_COL_BSTAT_FREE:
            if( fabs(redcost[i])> 1e-6 )
            {
               SCIPerrorMessage("stat col[%d] = %c, rd[%d] = %lg sense %d\n", i, icstat[i], i, redcost[i]*sense, sense);
               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }
            break;
         case QS_COL_BSTAT_UPPER:
            if( redcost[i]*sense > 1e-6 )
            {
               SCIPerrorMessage("stat col[%d] = %c, rd[%d] = %lg sense %d\n", i, icstat[i], i, redcost[i]*sense, sense);
               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }
            break;
         case QS_COL_BSTAT_LOWER:
            if( redcost[i]*sense < -1e-6 )
            {
               SCIPerrorMessage("stat col[%d] = %c, rd[%d] = %lg sense %d\n", i, icstat[i], i, redcost[i]*sense, sense);
               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }
            break;
         default:
            SCIPerrorMessage("unknown stat col[%d] = %c, rd[%d] = %lg\n", i, icstat[i], i, redcost[i]*sense);
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
         }
      }
   }
   else
   {
      SCIPerrorMessage("Getting solution with stat %d (not optimal)\n", stat);
   }
#endif
#endif

   rval = QSget_rhs(lpi->prob, lpi->irhs);
   QS_CONDRET(rval);
   rval = QSget_senses(lpi->prob, lpi->isen);
   QS_CONDRET(rval);

   /* build back the activity */
   if( activity )
   {
      for( i = 0; i < nrows; ++i )
      {
         switch( lpi->isen[i] )
         {
         case 'R':
         case 'E':
         case 'G':
            activity[i] = lpi->irhs[i] + lpi->irng[i];
            break;
         case 'L':
            activity[i] = lpi->irhs[i] - lpi->irng[i];
            break;
         default:
            SCIPerrorMessage("unknown sense %c\n", lpi->isen[i]);
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
         }
      }
   }

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPerrorMessage("SCIPlpiGetPrimalRay() not supported by QSopt.\n");

   return SCIP_LPERROR;
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual Farkas row multipliers */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(dualfarkas != NULL);

   SCIPdebugMessage("calling QSopt dual Farkas: %d cols, %d rows, %d non zeros\n", QSget_colcount (lpi->prob),
      QSget_rowcount(lpi->prob), QSget_nzcount(lpi->prob));

   rval = QSget_infeas_array(lpi->prob, dualfarkas);

   QS_RETURN(rval);
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   int rval = 0;
   int nit;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   rval = QSget_itcnt(lpi->prob, 0, 0, 0, 0, &nit);
   QS_CONDRET(rval);

   *iterations = nit - lpi->previt;
   lpi->previt = nit;

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
   int rval = 0, ncols, nrows;
   char* icstat = NULL;
   char* irstat = NULL;
   register int i;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("saving QSopt basis into %p/%p\n", (void*)cstat, (void*)rstat);

   ncols = QSget_colcount(lpi->prob);
   nrows = QSget_rowcount(lpi->prob);

   SCIP_CALL(ensureTabMem(lpi, nrows + ncols));

   icstat = lpi->ibas;
   irstat = lpi->ibas+ncols;
   rval = QSget_basis_array(lpi->prob, icstat, irstat);
   QS_CONDRET(rval);

   /* now we must transform QSopt codes into SCIP codes */
   for( i = 0; i < nrows; ++i )
   {
      switch( irstat[i] )
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
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }
   for( i = 0; i < ncols; ++i )
   {
      switch( icstat[i] )
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
         return SCIP_INVALIDDATA; /*lint !e527*/
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
   int rval = 0, ncols, nrows;
   register int i;
   char* icstat=0, *irstat = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("loading basis %p/%p into QSopt\n", (void*)cstat, (void*)rstat);

   ncols = QSget_colcount(lpi->prob);
   nrows = QSget_rowcount(lpi->prob);

   SCIP_CALL(ensureTabMem(lpi, ncols));

   icstat = lpi->ibas;
   irstat = lpi->ibas + ncols;

   /* now we must transform QSopt codes into SCIP codes */
   for( i = 0; i < nrows; ++i )
   {
      switch( rstat[i] )
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
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }
   for( i = 0; i < ncols; ++i )
   {
      switch( cstat[i] )
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
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   /* set the basis */
   rval = QSget_basis_array(lpi->prob, icstat, irstat);
   QS_RETURN(rval);
}

/** returns the indices of the basic columns and rows */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   int rval = 0, nrows, ncols;
   register int i;

   assert(lpi!=NULL);
   assert(lpi->prob!=NULL);

   SCIPdebugMessage("getting basis information\n");

   nrows = QSget_rowcount(lpi->prob);
   ncols = QSget_colcount(lpi->prob);
   rval = QSget_basis_order( lpi->prob, bind);
   QS_CONDRET(rval);

   /* transform QSopt basis header into SCIP format */
   for( i = 0; i < nrows; ++i )
   {
      if( bind[i] >= ncols )
         bind[i] = -(bind[i] - ncols - 1);
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
   int rval = 0;

   assert(lpi!=NULL);
   assert(lpi->prob!=NULL);
   SCIPdebugMessage("getting binv-row %d from Qsopt %d cols, %d rows, %d nonz\n", r, QSget_colcount(lpi->prob),
      QSget_rowcount(lpi->prob), QSget_nzcount(lpi->prob));

   rval = QSget_binv_row(lpi->prob, r, coef);
   QS_RETURN(rval);
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
{  /*lint --e{715} */
   assert(lpi!=NULL);
   assert(lpi->prob!=NULL);

   SCIPerrorMessage("SCIPlpiGetBInvCol() not supported by QSopt.\n");

   /* QSopt does not provide an interface for this yet */
   return SCIP_LPERROR;
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{  /*lint --e{715} */
   int rval = 0,ncols,nrows;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("getting binva-row %d\n", r);

   ncols = QSget_colcount(lpi->prob);
   nrows = QSget_rowcount(lpi->prob);

   SCIP_CALL(ensureTabMem(lpi, nrows+ncols));

   rval = QSget_tableau_row(lpi->prob, r, lpi->itab);
   QS_CONDRET(rval);

   /* copy local information to the outside */
   memcpy(coef, lpi->itab, sizeof(double)*ncols);

   return SCIP_OKAY;
}

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{  /*lint --e{715} */
   assert(lpi!=NULL);
   assert(lpi->prob!=NULL);

   SCIPerrorMessage("SCIPlpiGetBInvACol() not supported by QSopt.\n");

   /* QSopt does not provide an interface for this yet */
   return SCIP_LPERROR;
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
   int ncols, nrows;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(blkmem != NULL);
   assert(lpistate != NULL);

   ncols = QSget_colcount(lpi->prob);
   nrows = QSget_rowcount(lpi->prob);

   assert(ncols >= 0);
   assert(nrows >= 0);

   /* allocate lpistate data */
   SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows));
   SCIPdebugMessage("storing QSopt LPI state in %p (%d cols, %d rows)\n", (void*)*lpistate, ncols, nrows);

   /* get unpacked basis information from QSopt */
   SCIP_CALL( ensureColMem(lpi, ncols) );
   SCIP_CALL( ensureRowMem(lpi, nrows) );
   SCIP_CALL( SCIPlpiGetBase(lpi, lpi->iccnt, lpi->ircnt) );

   /* pack LPi state data */
   (*lpistate)->ncols = ncols;
   (*lpistate)->nrows = nrows;
   lpistatePack(*lpistate, lpi->iccnt, lpi->ircnt);

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
{  /*lint --e{715} */
   char* icstat = 0;
   char* irstat = 0;
   int i;
   int ncols;
   int nrows;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   /* if there was no basis information available, LPI state was not stored */
   if( lpistate == NULL )
      return SCIP_OKAY;

   /* continue test */
   ncols = QSget_colcount(lpi->prob);
   nrows = QSget_rowcount(lpi->prob);

   assert(ncols >= 0);
   assert(nrows >= 0);
   assert(lpistate->ncols <= ncols);
   assert(lpistate->nrows <= nrows);

   SCIPdebugMessage("loading LPI state %p (%d cols, %d rows) into QSopt LP with %d cols and %d rows\n", (void*)lpistate, lpistate->ncols,
      lpistate->nrows, ncols, nrows);

   if( lpistate->ncols == 0 || lpistate->nrows == 0 )
      return SCIP_OKAY;

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureColMem(lpi, ncols) );
   SCIP_CALL( ensureRowMem(lpi, nrows) );
   SCIP_CALL( ensureTabMem(lpi, nrows + ncols) );

   icstat = lpi->ibas;
   irstat = lpi->ibas + ncols;

   /* unpack LPi state data */
   lpistateUnpack(lpistate, lpi->iccnt, lpi->ircnt);

   /* extend the basis to the current LP beyond the previously existing columns */
   for( i = lpistate->ncols; i < ncols; ++i )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      /* get bounds from qsopt */
      rval = QSget_bounds_list(lpi->prob, 1, &i, &lb, &ub);
      if ( SCIPlpiIsInfinity(lpi, REALABS(lb)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         if ( SCIPlpiIsInfinity(lpi, REALABS(ub)) )
            lpi->iccnt[i] = SCIP_BASESTAT_ZERO;  /* variable is free */
         else
            lpi->iccnt[i] = SCIP_BASESTAT_UPPER; /* use finite upper bound */
      }
      else
         lpi->iccnt[i] = SCIP_BASESTAT_LOWER;    /* use finite lower bound */
   }
   for( i = lpistate->nrows; i < nrows; ++i )
      lpi->ircnt[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/

   /* convert the loaded basis into QSopt format */
   for( i = 0; i < nrows; ++i )
   {
      switch( lpi->ircnt[i] )
      {
      case SCIP_BASESTAT_LOWER:
         irstat[i] = QS_ROW_BSTAT_LOWER;
         break;
      case SCIP_BASESTAT_BASIC:
         irstat[i] = QS_ROW_BSTAT_BASIC;
         break;
      case SCIP_BASESTAT_UPPER:
         irstat[i] = QS_ROW_BSTAT_UPPER;
         break;
      default:
         SCIPerrorMessage("Unknown row basic status %d", lpi->ircnt[i]);
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }
   for( i = 0; i < ncols; ++i )
   {
      switch( lpi->iccnt[i] )
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
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   /* set the basis */
   rval = QSload_basis_array(lpi->prob, icstat, irstat);
   QS_RETURN(rval);
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   /**@todo implement SCIPlpiClearState() for QSopt */
   SCIPerrorMessage("QSopt interface does not implement SCIPlpiClearState()\n");

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
{  /*lint --e{715} */
   return (lpistate != NULL);
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("reading QSopt LP state from file <%s>\n", fname);

   rval = QSread_and_load_basis(lpi->prob, fname);
   if( rval )
   {
      SCIPerrorMessage("Error while loading basis from file <%s>.\n", fname);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** writes LP state (like basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   QSbas bas = 0;
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("writing QSopt LP state to file <%s>\n", fname);

   bas = QSget_basis(lpi->prob);
   QS_ERROR(bas == 0, "Could not get basis from problem.");   /*lint !e820*/
   assert(bas);

   rval = QSwrite_basis(lpi->prob, bas, fname);
   QSfree(bas);
   if( rval )
   {
      SCIPerrorMessage("Could not write basis to file <%s>.\n", fname);
      return SCIP_WRITEERROR;
   }

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
      rval = QSget_param(lpi->prob, QS_PARAM_SIMPLEX_SCALING,ival);
      if( *ival )
         *ival = TRUE;
      else
         *ival = FALSE;
      break;
   case SCIP_LPPAR_PRICING:
      *ival = lpi->pricing;
      break;
   case SCIP_LPPAR_LPINFO:
      rval = QSget_param(lpi->prob, QS_PARAM_SIMPLEX_DISPLAY, ival);
      if( *ival )
         *ival = TRUE;
      else
         *ival = FALSE;
      break;
   case SCIP_LPPAR_LPITLIM:
      rval = QSget_param(lpi->prob, QS_PARAM_SIMPLEX_MAX_ITERATIONS, ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   QS_RETURN(rval);
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiSetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("setting int parameter %d to %d\n", type, ival);

   switch( type )
   {
   case SCIP_LPPAR_SCALING:
      if( ival == TRUE )
         rval = QSset_param(lpi->prob, QS_PARAM_SIMPLEX_SCALING, 1);
      else
         rval = QSset_param(lpi->prob, QS_PARAM_SIMPLEX_SCALING, 0);
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
         rval = QSset_param(lpi->prob, QS_PARAM_PRIMAL_PRICING, QS_PRICE_PSTEEP);
         rval += QSset_param(lpi->prob, QS_PARAM_DUAL_PRICING, QS_PRICE_DSTEEP);
         break;
      case SCIP_PRICING_PARTIAL:
         rval = QSset_param(lpi->prob,QS_PARAM_PRIMAL_PRICING,QS_PRICE_PMULTPARTIAL);
         rval += QSset_param(lpi->prob,QS_PARAM_DUAL_PRICING,QS_PRICE_DMULTPARTIAL);
         break;
      case SCIP_PRICING_DEVEX:
         rval = QSset_param(lpi->prob,QS_PARAM_PRIMAL_PRICING,QS_PRICE_PDEVEX);
         rval += QSset_param(lpi->prob,QS_PARAM_DUAL_PRICING,QS_PRICE_DDEVEX);
         break;
      default:
         return SCIP_LPERROR;
      }
      break;
   case SCIP_LPPAR_LPINFO:
      if( ival )
         rval = QSset_param(lpi->prob, QS_PARAM_SIMPLEX_DISPLAY, 1);
      else
         rval = QSset_param(lpi->prob, QS_PARAM_SIMPLEX_DISPLAY, 0);
      break;
   case SCIP_LPPAR_LPITLIM:
      rval = QSset_param(lpi->prob, QS_PARAM_SIMPLEX_MAX_ITERATIONS, ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   QS_RETURN(rval);
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiGetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);
   assert(dval != NULL);

   SCIPdebugMessage("getting real parameter %d\n", type);

   switch( type )
   {
   case SCIP_LPPAR_LOBJLIM:
      rval = QSget_param_double(lpi->prob, QS_PARAM_OBJLLIM, dval);
      break;
   case SCIP_LPPAR_UOBJLIM:
      rval = QSget_param_double(lpi->prob, QS_PARAM_OBJULIM, dval);
      break;
   case SCIP_LPPAR_LPTILIM:
      rval = QSget_param_double(lpi->prob, QS_PARAM_SIMPLEX_MAX_TIME, dval);
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
SCIP_RETCODE SCIPlpiSetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   int rval = 0;

   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("setting real parameter %d to %g\n", type, dval);

   switch( type )
   {
   case SCIP_LPPAR_LPTILIM:
      rval = QSset_param_double(lpi->prob, QS_PARAM_SIMPLEX_MAX_TIME, dval);
      break;
   case SCIP_LPPAR_LOBJLIM:
      rval = QSset_param_double(lpi->prob, QS_PARAM_OBJLLIM, dval);
      break;
   case SCIP_LPPAR_UOBJLIM:
      rval = QSset_param_double(lpi->prob, QS_PARAM_OBJULIM, dval);
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

/** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiInfinity(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715} */
   return QS_MAXDOUBLE;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{  /*lint --e{715} */
   return (val >= QS_MAXDOUBLE);
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
   assert(lpi->prob != NULL);

   SCIPdebugMessage("reading LP from file <%s>\n", fname);

   if( lpi->prob != NULL )
      QSfree_prob(lpi->prob);

   lpi->solstat = 0;
   lpi->previt = 0;

   lpi->prob = QSread_prob(fname, "LP");
   if( lpi->prob == 0 )
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->prob != NULL);

   SCIPdebugMessage("writing LP to file <%s>\n", fname);

   if( QSwrite_prob (lpi->prob, fname, "LP") )
      return SCIP_WRITEERROR;

   return SCIP_OKAY;
}

/**@} */
