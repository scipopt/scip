/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lapack_calls.h
 * @brief  interface methods for lapack functions
 * @author Marc Pfetsch
 *
 * This file is used to call the LAPACK routine DSYEVR and DGETRF.
 *
 * LAPACK can be built with 32- or 64-bit integers, which is not visible to the outside. This interface tries to work
 * around this issue. Since the Fortran routines are called by reference, they only get a pointer. We always use 64-bit
 * integers on input, but reduce the output to 32-bit integers. We assume that all sizes can be represented in 32-bit
 * integers.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "lapack_calls.h"

#include "scip/def.h"
#include "scip/pub_message.h"                /* for debug and error message */
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/nlpi_ipopt.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/


/* we use 64 bit integers as the base type */
typedef long long int LAPACKINTTYPE;

/** Checks if a BMSallocMemory-call was successfull, otherwise returns SCIP_NOMEMORY */
#define BMS_CALL(x)   do                                                                                      \
                      {                                                                                       \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             SCIPerrorMessage("No memory in function call\n");                                \
                             return SCIP_NOMEMORY;                                                            \
                          }                                                                                   \
                      }                                                                                       \
                      while( FALSE )

/** transforms a SCIP_Real (that should be integer, but might be off by some numerical error) to an integer by adding 0.5 and rounding down */
#define SCIP_RealTOINT(x) ((LAPACKINTTYPE) (x + 0.5))

/*
 * BLAS/LAPACK Calls
 */

/**@name BLAS/LAPACK Calls */
/**@{ */

/** Define to a macro mangling the given C identifier (in lower and upper
 *  case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/** As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/** LAPACK Fortran subroutine DSYEVR */
void F77_FUNC(dsyevr, DSYEVR)(char* JOBZ, char* RANGE, char* UPLO,
   LAPACKINTTYPE* N, SCIP_Real* A, LAPACKINTTYPE* LDA,
   SCIP_Real* VL, SCIP_Real* VU,
   LAPACKINTTYPE* IL, LAPACKINTTYPE* IU,
   SCIP_Real* ABSTOL, LAPACKINTTYPE* M, SCIP_Real* W, SCIP_Real* Z,
   LAPACKINTTYPE* LDZ, LAPACKINTTYPE* ISUPPZ, SCIP_Real* WORK,
   LAPACKINTTYPE* LWORK, LAPACKINTTYPE* IWORK, LAPACKINTTYPE* LIWORK,
   LAPACKINTTYPE* INFO);

/**@} */

#ifdef SCIP_HAVE_LAPACK
/** converts a number stored in a long long int to an int, depending on big- or little endian machines
 *
 *  We assume that the number actually fits into an int. Thus, if more bits are used, we assume that the number is
 *  negative.
 */
static
int convertToInt(
   long long int         num                 /**< number to be converted */
   )
{
   long long int work;
   int checkval = 1;

   assert(sizeof(work) > sizeof(checkval)); /*lint !e506*/

   /* if we have a little-endian machine (e.g, x86), the sought value is in the bottom part */
   if ( *(int8_t*)&checkval != 0 ) /*lint !e774*/
   {
      /* if the top part is nonzero, we assume that the number is negative */
      if ( *((int8_t*)&num + 4) != 0 ) /*lint !e2662*/
      {
         work = -num;
         return -(*((int*)&work));
      }
      return *((int*)&num);
   }

   /* otherwise we have a big-endian machine (e.g., PowerPC); the sought value is in the top part */
   assert( *(int8_t*)&checkval == 0 );

   /* if the bottom part is nonzero, we assume that the number is negative */
   if ( *(int8_t*)&num != 0 ) /*lint !e774*/
   {
      work = -num;
      return -(*((int*)&work + 4)); /*lint !e2662*/
   }
   return *((int*)&num + 4);
}
#endif

/*
 * Functions
 */

/**@name Functions */
/**@{ */


/** returns whether Lapack s available, i.e., whether it has been linked in */
SCIP_Bool SCIPlapackIsAvailable(void)
{
   if ( SCIPisIpoptAvailableIpopt() )
      return TRUE;

#ifdef SCIP_HAVE_LAPACK
   return TRUE;
#endif
   return FALSE;
}

#ifdef SCIP_HAVE_LAPACK
/** computes eigenvalues of a symmetric matrix using LAPACK */
static
SCIP_RETCODE lapackComputeEigenvalues(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_Bool             geteigenvectors,    /**< Should also the eigenvectors be computed? */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix data on input (size n * n); eigenvectors on output if geteigenvectors == TRUE */
   SCIP_Real*            eigenvalues         /**< pointer to store eigenvalue */
   )
{
   LAPACKINTTYPE* IWORK;
   LAPACKINTTYPE* ISUPPZ;
   LAPACKINTTYPE N;
   LAPACKINTTYPE INFO;
   LAPACKINTTYPE LDA;
   LAPACKINTTYPE WISIZE;
   LAPACKINTTYPE IL;
   LAPACKINTTYPE IU;
   LAPACKINTTYPE M;
   LAPACKINTTYPE LDZ;
   LAPACKINTTYPE LWORK;
   LAPACKINTTYPE LIWORK;
   SCIP_Real* WORK;
   SCIP_Real* WTMP;
   SCIP_Real* Z = NULL;
   SCIP_Real ABSTOL;
   SCIP_Real WSIZE;
   SCIP_Real VL;
   SCIP_Real VU;
   char JOBZ;
   char RANGE;
   char UPLO;

   assert( bufmem != NULL );
   assert( n > 0 );
   assert( n < INT_MAX );
   assert( A != NULL );
   assert( eigenvalues != NULL );

   N = n;
   JOBZ = geteigenvectors ? 'V' : 'N';
   RANGE = 'A';
   UPLO = 'L';
   LDA  = n;
   ABSTOL = 0.0; /* we use abstol = 0, since some lapack return an error otherwise */
   VL = -1e20;
   VU = 1e20;
   IL = 0;
   IU = n;
   M = n;
   LDZ = n;
   INFO = 0LL;

   /* standard LAPACK workspace query, to get the amount of needed memory */
   LWORK = -1LL;
   LIWORK = -1LL;

   /* this computes the internally needed memory and returns this as (the first entries of [the 1x1 arrays]) WSIZE and WISIZE */
   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, NULL, &LDA,
      NULL, NULL,
      &IL, &IU,
      &ABSTOL, &M, NULL, NULL,
      &LDZ, NULL, &WSIZE,
      &LWORK, &WISIZE, &LIWORK,
      &INFO);

   if ( convertToInt(INFO) != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d.\n", convertToInt(INFO));
      return SCIP_ERROR;
   }

   /* allocate workspace */
   LWORK = SCIP_RealTOINT(WSIZE);
   LIWORK = WISIZE;

   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WORK, (int) LWORK) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &IWORK, (int) LIWORK) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WTMP, (int) N) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &ISUPPZ, 2) ); /*lint !e506*/
   if ( geteigenvectors )
   {
      BMS_CALL( BMSallocBufferMemoryArray(bufmem, &Z, n * n) );
   }

   /* call the function */
   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, A, &LDA,
      &VL, &VU,
      &IL, &IU,
      &ABSTOL, &M, WTMP, Z,
      &LDZ, ISUPPZ, WORK,
      &LWORK, IWORK, &LIWORK,
      &INFO);

   /* handle output */
   if ( convertToInt(INFO) == 0 )
   {
      int m;
      int i;
      int j;

      m = convertToInt(M);
      for (i = 0; i < m; ++i)
         eigenvalues[i] = WTMP[i];
      for (i = m; i < n; ++i)
         eigenvalues[i] = SCIP_INVALID;

      /* possibly overwrite matrix with eigenvectors */
      if ( geteigenvectors )
      {
         for (i = 0; i < m; ++i)
         {
            for (j = 0; j < n; ++j)
               A[i * n + j] = Z[i * n + j];
         }
      }
   }

   /* free memory */
   BMSfreeBufferMemoryArrayNull(bufmem, &Z);
   BMSfreeBufferMemoryArray(bufmem, &ISUPPZ);
   BMSfreeBufferMemoryArray(bufmem, &WTMP);
   BMSfreeBufferMemoryArray(bufmem, &IWORK);
   BMSfreeBufferMemoryArray(bufmem, &WORK);

   if ( convertToInt(INFO) != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d.\n", convertToInt(INFO));
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}
#endif

/** computes eigenvalues and eigenvectors of a dense symmetric matrix
 *
 *  Calls Lapack's DSYEV function.
 */
SCIP_RETCODE SCIPlapackComputeEigenvalues(
   BMS_BUFMEM*           bufmem,             /**< buffer memory (or NULL if IPOPT is used) */
   SCIP_Bool             geteigenvectors,    /**< should also eigenvectors should be computed? */
   int                   N,                  /**< dimension */
   SCIP_Real*            a,                  /**< matrix data on input (size N*N); eigenvectors on output if geteigenvectors == TRUE */
   SCIP_Real*            w                   /**< array to store eigenvalues (size N) (or NULL) */
   )
{
   /* if IPOPT is available, call its LAPACK routine */
   if ( SCIPisIpoptAvailableIpopt() )
   {
      SCIP_CALL( SCIPcallLapackDsyevIpopt(geteigenvectors, N, a, w) );
   }
   else
   {
      assert( bufmem != NULL );
#ifdef SCIP_HAVE_LAPACK
      SCIP_CALL( lapackComputeEigenvalues(bufmem, geteigenvectors, N, a, w) );
#else
      SCIPerrorMessage("Lapack not available.\n");
      return SCIP_PLUGINNOTFOUND;
#endif
   }

   return SCIP_OKAY;
}

/**@} */
