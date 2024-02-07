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

#include "scip/lapack_calls.h"

#include "scip/def.h"
#include "scip/pub_message.h"                /* for debug and error message */
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/nlpi_ipopt.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/


#ifdef SCIP_WITH_LAPACK
/* we use 64 bit integers as the base type */
typedef int64_t LAPACKINTTYPE;

/** transforms a SCIP_Real (that should be integer, but might be off by some numerical error) to an integer by adding 0.5 and rounding down */
#define SCIP_RealTOINT(x) ((LAPACKINTTYPE) (x + 0.5))

/*
 * BLAS/LAPACK Calls
 */

/**@name BLAS/LAPACK Calls */
/**@{ */

/** Define to a macro mangling the given C identifier (in lower and upper
 *  case), which must not contain underscores, for linking with Fortran. */
#ifdef FNAME_LCASE_DECOR
#define F77_FUNC(name,NAME) name ## _
#endif
#ifdef FNAME_UCASE_DECOR
#define F77_FUNC(name,NAME) NAME ## _
#endif
#ifdef FNAME_LCASE_NODECOR
#define F77_FUNC(name,NAME) name
#endif
#ifdef FNAME_UCASE_NODECOR
#define F77_FUNC(name,NAME) NAME
#endif

/* use backup ... */
#ifndef F77_FUNC
#define F77_FUNC(name,NAME) name ## _
#endif


/** LAPACK Fortran subroutine DSYEVR */
void F77_FUNC(dsyevr, DSYEVR)(char* JOBZ, char* RANGE, char* UPLO,
   LAPACKINTTYPE* N, SCIP_Real* A, LAPACKINTTYPE* LDA,
   SCIP_Real* VL, SCIP_Real* VU,
   LAPACKINTTYPE* IL, LAPACKINTTYPE* IU,
   SCIP_Real* ABSTOL, LAPACKINTTYPE* M, SCIP_Real* W, SCIP_Real* Z,
   LAPACKINTTYPE* LDZ, LAPACKINTTYPE* ISUPPZ, SCIP_Real* WORK,
   LAPACKINTTYPE* LWORK, LAPACKINTTYPE* IWORK, LAPACKINTTYPE* LIWORK,
   LAPACKINTTYPE* INFO);

/** LAPACK Fortran subroutine DGETRF */
void F77_FUNC(dgetrf, DGETRF)(LAPACKINTTYPE* M, LAPACKINTTYPE* N, SCIP_Real* A,
   LAPACKINTTYPE* LDA, LAPACKINTTYPE* IPIV, LAPACKINTTYPE* INFO);

/** LAPACK Fortran subroutine DGETRS */
void F77_FUNC(dgetrs, DGETRS)(char* TRANS, LAPACKINTTYPE* N, LAPACKINTTYPE* NRHS,
   SCIP_Real* A, LAPACKINTTYPE* LDA, LAPACKINTTYPE* IPIV, SCIP_Real* B, LAPACKINTTYPE* LDB,
   LAPACKINTTYPE* INFO);

/** LAPACK Fortran subroutine ivlayer */
void F77_FUNC(ilaver, ILAVER)(LAPACKINTTYPE* MAJOR, LAPACKINTTYPE* MINOR, LAPACKINTTYPE* PATCH);

/**@} */
#endif

/*
 * Functions
 */

/**@name Functions */
/**@{ */

/** returns whether Lapack is available, i.e., whether it has been linked in */
SCIP_Bool SCIPlapackIsAvailable(void)
{
   if ( SCIPisIpoptAvailableIpopt() )
      return TRUE;

#ifdef SCIP_WITH_LAPACK
   return TRUE;
#else
   return FALSE;
#endif
}

#ifdef SCIP_WITH_LAPACK
/** converts a number stored in a int64_t to an int, depending on big- or little endian machines
 *
 *  We assume that the number actually fits into an int. Thus, if more bits are used, we assume that the number is
 *  negative.
 */
static
int convertToInt(
   int64_t               num                 /**< number to be converted */
   )
{
   union
   {
      int64_t big;
      int     small[2];
   } work;
   int checkval = 1;

   assert(sizeof(work) > sizeof(checkval)); /*lint !e506*/

   work.big = num;

   /* if we have a little-endian machine (e.g, x86), the sought value is in the bottom part */
   if ( *(int8_t*)&checkval != 0 ) /*lint !e774*/
   {
      /* if the top part is nonzero, we assume that the number is negative */
      if ( work.small[1] != 0 ) /*lint !e2662*/
      {
         work.big = -num;
         return -work.small[0];
      }
      return work.small[0];
   }

   /* otherwise we have a big-endian machine (e.g., PowerPC); the sought value is in the top part */
   assert( *(int8_t*)&checkval == 0 );

   /* if the bottom part is nonzero, we assume that the number is negative */
   if ( work.small[0] != 0 ) /*lint !e774*/
   {
      work.big = -num;
      return -work.small[1]; /*lint !e2662*/
   }
   return work.small[1];
}
#endif

/** returns whether Lapack s available, i.e., whether it has been linked in */
void SCIPlapackVersion(
   int*                  majorver,           /**< major version number */
   int*                  minorver,           /**< minor version number */
   int*                  patchver            /**< patch version number */
   )
{
#ifdef SCIP_WITH_LAPACK
   LAPACKINTTYPE MAJOR = 0LL;
   LAPACKINTTYPE MINOR = 0LL;
   LAPACKINTTYPE PATCH = 0LL;
#endif

   assert( majorver != NULL );
   assert( minorver != NULL );
   assert( patchver != NULL );

#ifdef SCIP_WITH_LAPACK
   F77_FUNC(ilaver, ILAVER)(&MAJOR, &MINOR, &PATCH);

   *majorver = convertToInt(MAJOR);
   *minorver = convertToInt(MINOR);
   *patchver = convertToInt(PATCH);
#else
   *majorver = -1;
   *minorver = -1;
   *patchver = -1;
#endif
}

#ifdef SCIP_WITH_LAPACK
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

   SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &WORK, (int) LWORK) );
   SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &IWORK, (int) LIWORK) );
   SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &WTMP, (int) N) );
   SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &ISUPPZ, 2) ); /*lint !e506*/
   if ( geteigenvectors )
   {
      SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &Z, n * n) ); /*lint !e647*/
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
#ifdef SCIP_WITH_LAPACK
      SCIP_CALL( lapackComputeEigenvalues(bufmem, geteigenvectors, N, a, w) );
#else
      SCIPerrorMessage("Lapack not available.\n");
      return SCIP_PLUGINNOTFOUND;
#endif
   }

   return SCIP_OKAY;
}

/** solves a linear problem of the form Ax = b for a regular matrix A
 *
 *  Calls Lapacks DGETRF routine to calculate a LU factorization and uses this factorization to solve
 *  the linear problem Ax = b.
 *
 *  Code taken from nlpi_ipopt.cpp
 */
SCIP_RETCODE SCIPlapackSolveLinearEquations(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   n,                  /**< dimension */
   SCIP_Real*            A,                  /**< matrix data on input (size N*N); filled column-wise */
   SCIP_Real*            b,                  /**< right hand side vector (size N) */
   SCIP_Real*            x,                  /**< buffer to store solution (size N) */
   SCIP_Bool*            success             /**< pointer to store if the solving routine was successful */
   )
{
   assert( n > 0 );
   assert( A != NULL );
   assert( b != NULL );
   assert( x != NULL );
   assert( success != NULL );

   /* if possible, use IPOPT */
   if ( SCIPisIpoptAvailableIpopt() )
   {
      SCIP_CALL( SCIPsolveLinearEquationsIpopt(n, A, b, x, success) );
   }
   else
   {
#ifdef SCIP_WITH_LAPACK
      LAPACKINTTYPE INFO;
      LAPACKINTTYPE N;
      LAPACKINTTYPE* pivots;
      SCIP_Real* Atmp = NULL;
      SCIP_Real* btmp = NULL;

      assert( bufmem != NULL );

      SCIP_ALLOC( BMSduplicateBufferMemoryArray(bufmem, &Atmp, A, n * n) ); /*lint !e647*/
      SCIP_ALLOC( BMSduplicateBufferMemoryArray(bufmem, &btmp, b, n) );
      SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &pivots, n) );

      /* compute LU factorization */
      N = n;
      F77_FUNC(dgetrf, DGETRF)(&N, &N, Atmp, &N, pivots, &INFO);

      if ( convertToInt(INFO) != 0 )
      {
         SCIPdebugMessage("There was an error when calling DGETRF. INFO = %d\n", convertToInt(INFO));
         *success = FALSE;
      }
      else
      {
         LAPACKINTTYPE NRHS = 1LL;
         char TRANS = 'N';

         /* solve system */
         F77_FUNC(dgetrs, DGETRS)(&TRANS, &N, &NRHS, Atmp, &N, pivots, btmp, &N, &INFO);

         if ( convertToInt(INFO) != 0 )
         {
            SCIPdebugMessage("There was an error when calling DGETRF. INFO = %d\n", convertToInt(INFO));
            *success = FALSE;
         }
         else
            *success = TRUE;

         /* copy the solution */
         BMScopyMemoryArray(x, btmp, n);
      }

      BMSfreeBufferMemoryArray(bufmem, &pivots);
      BMSfreeBufferMemoryArray(bufmem, &btmp);
      BMSfreeBufferMemoryArray(bufmem, &Atmp);
#else
      SCIP_UNUSED(bufmem);

      /* call fallback solution in nlpi_ipopt_dummy */
      SCIP_CALL( SCIPsolveLinearEquationsIpopt(n, A, b, x, success) );
#endif
   }

   return SCIP_OKAY;
}

/**@} */
