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

/**@file   rectlu.h
 * @brief  rectlu user interface
 * @author David Applegate
 * @author Bill Cook
 * @author Sanjeeb Dash
 * @author Daniel Espinoza
 * @author Dan Steffy
 * @author Kati Wolter
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "rectlu_num.h"


#ifndef __RECTLU_H
#define __RECTLU_H

typedef struct qsnum_svector {
	int nzcnt;
	int *indx;
	int size;
	QSnum_type *coef;
} qsnum_svector;

typedef struct qsnum_uc_info {
    int cbeg;
    int nzcnt;
    int next;
    int prev;
    int delay;
} qsnum_uc_info;

typedef struct qsnum_ur_info {
    QSnum_type max;
    int rbeg;
    int nzcnt;
    int pivcnt;
    int next;
    int prev;
    int delay;
} qsnum_ur_info;

typedef struct qsnum_lc_info {
    int cbeg;
    int nzcnt;
    int c;
    int crank;
    int delay;
} qsnum_lc_info;

typedef struct qsnum_lr_info {
    int rbeg;
    int nzcnt;
    int r;
    int rrank;
    int delay;
} qsnum_lr_info;

typedef struct qsnum_er_info {
    int rbeg;
    int nzcnt;
    int r;
} qsnum_er_info;


typedef struct qsnum_factor_work {
  int max_k;
  QSnum_type fzero_tol;
  QSnum_type szero_tol;
  QSnum_type partial_tol;
  double ur_space_mul;
  double uc_space_mul;
  double lc_space_mul;
  double lr_space_mul;
  double er_space_mul;
  double grow_mul;
  int p;
  int etamax;
  double minmult;
  double maxmult;
  double updmaxmult;
  double dense_fract;
  int dense_min;

  QSnum_type maxelem_orig;
  int nzcnt_orig;
  QSnum_type maxelem_factor;
  int nzcnt_factor;
  QSnum_type maxelem_cur;
  int nzcnt_cur;

  QSnum_type partial_cur;

  int dimr;
  int dimc;
  int maxdim;
  int stage;
  int nstages;
  int etacnt;
  QSnum_type *work_coef;
  int *work_indx;
  qsnum_uc_info *uc_inf;
  qsnum_ur_info *ur_inf;
  qsnum_lc_info *lc_inf;
  qsnum_lr_info *lr_inf;
  qsnum_er_info *er_inf;
  int *ucindx;                         /* row index for column data */
  int *ucrind;                         /* index of column in row data */
  QSnum_type *uccoef;                       /* coefficient for column data */
  int *urindx;                         /* col index for row data */
  int *urcind;                         /* index of row in column data */
  QSnum_type *urcoef;                       /* coefficient for row data */
  int *lcindx;                         /* row index for L data */
  QSnum_type *lccoef;                       /* coefficient for L row data */
  int *lrindx;                         /* col index for L data */
  QSnum_type *lrcoef;                       /* coefficient for L col data */
  int *erindx;                         /* col index for eta data */
  QSnum_type *ercoef;                       /* coefficient for eta data */
  int *rperm;
  int *rrank;
  int *cperm;
  int *crank;
  qsnum_svector xtmp;
  int ur_freebeg;
  int ur_space;
  int uc_freebeg;
  int uc_space;
  int lc_freebeg;
  int lc_space;
  int lr_freebeg;
  int lr_space;
  int er_freebeg;
  int er_space;

  int *p_nsing;
  int **p_singr;
  int **p_singc;

  QSnum_type *dmat;
  int drows;
  int dcols;
  int dsize;
  int dense_base;
} qsnum_factor_work;


/** factorizes the matrix and points f to factor work or returns error */
int RECTLUbuildFactorization(
   qsnum_factor_work** f,                /**< pointer to store factor work*/
   int              n,                   /**< number of rows in matrix */
   int              m,                   /**< number of columns in matrix*/
   int*             basisx,              /**< basis columns of matrix */
   mpq_t*           matval,              /**< values of matrix entries */
   int*             matind,              /**< row index of matrix entries */
   int*             matbeg,              /**< start of columns in sparse matrix */
   int*             matcnt               /**< length of column in sparse matrix */
   );

/** solves the system using factor work f */
int RECTLUsolveSystem(
   qsnum_factor_work* f,                 /**< pointer to factor work*/
   int              n,                   /**< number of rows in matrix */
   int              m,                   /**< number of columns in matrix*/
   mpq_t*           rhs,                 /**< right hand side of system to solve */
   mpq_t*           sol                  /**< solution to system */
   );

/** frees factor work f */
void RECTLUfreeFactorization(
   qsnum_factor_work* f                 /**< pointer to store factor work*/
   );

#endif	/*  __RECTLU_H */
