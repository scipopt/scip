/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                        This file is part of the program                   */
/*            RECTLU --- Algorithm for Exact Rectangular LU Factorization    */
/*                                                                           */
/*    Copyright (C) 2009-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*   RECTLU is distributed under the terms of the ZIB Academic License.      */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with RECTLU; see the file COPYING.                                 */
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


#ifdef WITH_GMP
void QSnum_svector_init (qsnum_svector * s),
     QSnum_svector_free (qsnum_svector * s);

int init_sxvector (qsnum_svector *v, int n);
int clear_sxvector (qsnum_svector *v);

int QSnum_svector_alloc (qsnum_svector * s, int nzcnt);



void
    QSnum_factor_init (void),
    QSnum_factor_clear (void),
    QSnum_factor_init_factor_work (qsnum_factor_work * f),
    QSnum_factor_free_factor_work (qsnum_factor_work * f),
    QSnum_factor_ftran (qsnum_factor_work * f, qsnum_svector * a,
        qsnum_svector * x),
    QSnum_factor_btran (qsnum_factor_work * f, qsnum_svector * a,
        qsnum_svector * x);


int QSnum_factor_create_factor_work (qsnum_factor_work * f, int dimr, int dimc),
  QSnum_factor (qsnum_factor_work * f, int *basis, int *cbeg, int *clen,
		int *cindx, QSnum_type * ccoef, int *p_nsing, int **p_singr,
		int **p_singc);
#define E_CHECK_FAILED 6
#define E_NO_PIVOT 7
#define E_FACTOR_BLOWUP 8
#define E_UPDATE_NOSPACE 9
#define E_UPDATE_SINGULAR_ROW 10
#define E_UPDATE_SINGULAR_COL 11
#define E_SING_NO_DATA 12
#define E_SINGULAR_INTERNAL 13
#define SPARSE_FACTOR 0.05

#endif /* __RECTLU_FACTOR_H */

#ifndef __CG_UTIL_H
#define __CG_UTIL_H

#define CG_GENERAL_ERROR -1
#define CG_NO_MEMORY 2
#define CG_NULL_PTR  3

#define CG_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

#define CGcheck_rval(rval,msg) {                                           \
    if ((rval)) {                                                          \
        fprintf (stderr, "%s\n", (msg));                                   \
        goto CLEANUP;                                                      \
    }                                                                      \
}

#define CGcheck_NULL(item,msg) {                                           \
    if ((!item)) {                                                         \
        fprintf (stderr, "%s\n", (msg));                                   \
        rval = CG_NULL_PTR;                                               \
        goto CLEANUP;                                                      \
    }                                                                      \
}

#define CG_UTIL_SAFE_MALLOC(nnum,type,varname)                             \
     (type *) CGutil_allocrus (((size_t) (nnum)) * sizeof (type))

#define CG_IFFREE(object,type) {                                           \
    if ((object)) {                                                         \
       CGutil_freerus ((void *) (object));                                 \
       object = (type *) NULL;                                              \
    }}

void
   *CGutil_allocrus (size_t size),
    CGutil_freerus (void *p);

#define CG_CHECKnull(expr, msg) \
    { if ((expr) == NULL)  {  \
         fprintf (stderr, msg); \
         rval = CG_NULL_PTR; \
         goto CLEANUP;  \
      } }

#define CG_CLEANUP_IF(rval)      { if ((rval) != 0) { goto CLEANUP; } }
#define CG_CLEANUP               goto CLEANUP

#define CG_SAFE_MALLOC(lhs, n, type) \
    { lhs = CG_UTIL_SAFE_MALLOC(n, type, lhs); \
      if (lhs == NULL)  {  \
         fprintf (stderr, "Out of memory\n"); \
         rval = CG_NO_MEMORY; \
         goto CLEANUP;  \
      }}

#define CG_SAFE_MALLOC_no_rval(lhs, n, type) \
    { lhs = CG_UTIL_SAFE_MALLOC(n, type, lhs); \
      if (lhs == NULL)  {  \
         fprintf (stderr, "Out of memory\n"); \
         goto CLEANUP;  \
      }}
#endif

#endif /* __CG_UTIL_H */
