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

/**@file   rectlu_factor.h
 * @brief  rectlu internal functions interface
 * @author David Applegate
 * @author Bill Cook
 * @author Sanjeeb Dash
 * @author Daniel Espinoza
 * @author Dan Steffy
 * @author Kati Wolter
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "rectlu_num.h"
#include "rectlu.h"

#ifndef __RECTLU_FACTOR_H
#define __RECTLU_FACTOR_H


#ifdef SCIP_WITH_GMP

/** initialize empty sparse qsnum vector */
void QSnum_svector_init(
   qsnum_svector *       s                   /**< sparse vector */
   );

/** free sparse qsnum vector */
void QSnum_svector_free(
   qsnum_svector *       s                   /**< sparse vector */
   );

/** initialize sparse qsnum vector with n nonzeros */
int init_sxvector(
   qsnum_svector*        v,                   /**< sparse vector */
   int                   n                    /**< nonzero count */
   );

/** free sparse qsnum vector */
int clear_sxvector (
   qsnum_svector *       v                   /**< sparse vector */
   );

/** allocate sparse qsnum vector */
int QSnum_svector_alloc(
   qsnum_svector *       s,                  /**< sparse vector */
   int                   nzcnt               /**< nonzero count */
   );


/** initializes global zero and one */
void QSnum_factor_init (void);

/** frees global zero and one */
void QSnum_factor_clear (void);

/** initializes temporary work used for factorization */
void QSnum_factor_init_factor_work(
   qsnum_factor_work *   f                   /**< factorization work */
   );

/** frees temporary work used for factorization */
void QSnum_factor_free_factor_work(
   qsnum_factor_work *   f                   /**< factorization work */
   );

/** FTRAN routine: solves Bx^t=a^t (or, x B^t = a) for x */
void QSnum_factor_ftran(
   qsnum_factor_work *   f,                  /**< factorization work data */
   qsnum_svector *       a,                  /**< rhs vector */
   qsnum_svector *       x                   /**< solution vector */
   );

/** BTRAN routine: solves B^tx^t=a^t (or, xB = a) for x */
void QSnum_factor_btran (
   qsnum_factor_work *   f,                  /**< factorization work data */
   qsnum_svector *       a,                  /**< rhs vector */
   qsnum_svector *       x                   /**< solution vector */
   );

/** allocates memory used for LU factorization */
int QSnum_factor_create_factor_work(
   qsnum_factor_work *   f,                   /**< factorization work */
   int                   dimr,                /**< row dimension of matrix */
   int                   dimc                 /**< column dimension of matrix */
   );

/**
 * build LU factorization of basis matrix or return singularity
 *
 * coefficients of matrix are ccoef in a sparse matrix column format
 * cbeg tells where each column begins in ccoef, clen tells how many nonzeros
 * are stored, and cindx stores the row indices of each coefficient
 */
int QSnum_factor(
   qsnum_factor_work *   f,                  /**< factorization work storage */
   int *                 basis,              /**< stores indices of basis columns */
   int *                 cbeg,               /**< col beginning indices for constraint matrix */
   int *                 clen,               /**< col lengths for constraint matrix */
   int *                 cindx,              /**< row indices for each cons. matrix index */
   QSnum_type *          ccoef,              /**< coef in each position of cons. matrix */
   int *                 p_nsing,            /**< returns rank deficiency (nonzero if singular) */
   int **                p_singr,            /**< indices of dependent rows */
   int **                p_singc             /**< indices of dependent columns */
   );



/*
 *  return codes for RectLU code
 */
#define E_CHECK_FAILED 6
#define E_NO_PIVOT 7
#define E_FACTOR_BLOWUP 8
#define E_UPDATE_NOSPACE 9
#define E_UPDATE_SINGULAR_ROW 10
#define E_UPDATE_SINGULAR_COL 11
#define E_SING_NO_DATA 12
#define E_SINGULAR_INTERNAL 13

/* factor used classify which vectors should be treated as sparse */
#define SPARSE_FACTOR 0.05

#endif /* __RECTLU_FACTOR_H */

#ifndef __CG_UTIL_H
#define __CG_UTIL_H

#define CG_GENERAL_ERROR -1
#define CG_NO_MEMORY 2
#define CG_NULL_PTR  3


/*
 *   Utility macros used within RectLU code
 */

#define CG_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

#define CGcheck_rval(rval,msg) {                \
      if ((rval)) {                             \
         fprintf (stderr, "%s\n", (msg));       \
         goto CLEANUP;                          \
      }                                         \
   }

#define CGcheck_NULL(item,msg) {                \
      if ((!item)) {                            \
         fprintf (stderr, "%s\n", (msg));       \
         rval = CG_NULL_PTR;                    \
         goto CLEANUP;                          \
      }                                         \
   }

#define CG_UTIL_SAFE_MALLOC(nnum,type,varname)                  \
   (type *) CGutil_allocrus (((nnum)) * sizeof (type))

#define CG_IFFREE(object,type) {                \
      if ((object)) {                           \
         CGutil_freerus ((void *) (object));    \
         object = (type *) NULL;                \
      }}

/** allocates memory */
void *CGutil_allocrus (
   size_t                size                /**< length of array to be allocated */
   );

/** frees memory */
void CGutil_freerus (
   void *                p                   /**< array to be freed */
   );

#define CG_CHECKnull(expr, msg)                 \
   { if ((expr) == NULL)  {                     \
         fprintf (stderr, msg);                 \
         rval = CG_NULL_PTR;                    \
         goto CLEANUP;                          \
      } }

#define CG_CLEANUP_IF(rval)      { if ((rval) != 0) { goto CLEANUP; } }
#define CG_CLEANUP               goto CLEANUP

#define CG_SAFE_MALLOC(lhs, n, type)            \
   { lhs = CG_UTIL_SAFE_MALLOC(n, type, lhs);   \
      if (lhs == NULL)  {                       \
         fprintf (stderr, "Out of memory\n");   \
         rval = CG_NO_MEMORY;                   \
         goto CLEANUP;                          \
      }}

#define CG_SAFE_MALLOC_no_rval(lhs, n, type)    \
   { lhs = CG_UTIL_SAFE_MALLOC(n, type, lhs);   \
      if (lhs == NULL)  {                       \
         fprintf (stderr, "Out of memory\n");   \
         goto CLEANUP;                          \
      }}
#endif

#endif /* __CG_UTIL_H */
