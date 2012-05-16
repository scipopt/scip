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

/**@file   rectlu_factor.c
 * @brief  rectlu internal functions
 * @author David Applegate
 * @author Bill Cook
 * @author Sanjeeb Dash
 * @author Daniel Espinoza
 * @author Dan Steffy
 * @author Kati Wolter
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <assert.h>

#include "rectlu_factor.h"
#include "rectlu.h"

static QSnum_type qsnum_zeroLpNum;
static QSnum_type qsnum_oneLpNum;
static QSnum_type __t__;

/* From cg_util.c */
void *CGutil_allocrus (size_t size)
{
    void *mem = (void *) NULL;

    if (size == 0) { fprintf (stderr, "Warning: 0 bytes allocated\n"); }

    mem = (void *) malloc (size);
    if (mem == (void *) NULL) {
        fprintf (stderr, "Out of memory. Asked for %d bytes\n", (int) size);
    }
    return mem;
}

void CGutil_freerus (void *p)
{
    if (!p) {
        fprintf (stderr, "Warning: null pointer freed\n");
        return;
    }
    free (p);
}
/*/ From cg_util.c */







void QSnum_factor_init (void)
{
    QSnum_Init(qsnum_zeroLpNum);
    QSnum_Init(qsnum_oneLpNum);
    QSnum_Init(__t__);
    QSnum_SetZero(qsnum_zeroLpNum);
    QSnum_SetOne(qsnum_oneLpNum);
}

void QSnum_factor_clear (void)
{
    QSnum_Clear(qsnum_zeroLpNum);
    QSnum_Clear(qsnum_oneLpNum);
    QSnum_Clear(__t__);
}

void QSnum_factor_init_factor_work (qsnum_factor_work * f)
{
    QSnum_Init (f->fzero_tol);
    QSnum_Init (f->szero_tol);
    QSnum_Init (f->partial_tol);
    QSnum_Init (f->partial_cur);
    f->max_k = 1000;               /* must be less than 46340 (2^15.5) */
    QSnum_SetEpsilon (f->fzero_tol);                 /* epsilon (2^-50) */
    QSnum_SetEpsilon (f->szero_tol);                 /* epsilon (2^-50) */
    QSnum_SetDelta   (f->partial_tol);               /* delta (2^-7)    */
    f->ur_space_mul = 2.0;
    f->uc_space_mul = 1.1;
    f->lc_space_mul = 1.1;
    f->er_space_mul = 1000.0;
    f->grow_mul = 1.5;
    f->p = 4;
    f->etamax = 100;
    f->minmult = 1e3;
    f->maxmult = 1e5;
    f->updmaxmult = 1e7;
    /* f->dense_fract = .25;*/
    f->dense_fract = 1;         /*this means we only use dense for a totally dense matrix */
    /*    f->dense_min = 25;*/
    f->dense_min = 1;
    QSnum_Copy (f->partial_cur, f->partial_tol);
    f->work_coef = 0;
    f->work_indx = 0;
    f->uc_inf = 0;
    f->ur_inf = 0;
    f->lc_inf = 0;
    f->lr_inf = 0;
    f->er_inf = 0;
    f->ucindx = 0;
    f->ucrind = 0;
    f->uccoef = 0;
    f->urindx = 0;
    f->urcind = 0;
    f->urcoef = 0;
    f->lcindx = 0;
    f->lccoef = 0;
    f->lrindx = 0;
    f->lrcoef = 0;
    f->erindx = 0;
    f->ercoef = 0;
    f->rperm = 0;
    f->rrank = 0;
    f->cperm = 0;
    f->crank = 0;
    f->dmat = 0;
    f->dsize = 0;
    QSnum_svector_init (&f->xtmp);
}

void QSnum_factor_free_factor_work (qsnum_factor_work * f)
{
    QSnum_Clear (f->fzero_tol);
    QSnum_Clear (f->szero_tol);
    QSnum_Clear (f->partial_tol);
    QSnum_Clear (f->partial_cur);
    QSnum_FreeArray (f->work_coef, f->maxdim);
    CG_IFFREE (f->work_indx, int);
    CG_IFFREE (f->uc_inf, qsnum_uc_info);
    if (f->dimr + f->max_k > 0 && f->ur_inf)
    {
        unsigned int i = f->dimr + f->max_k + 1;
        while (i--)
            QSnum_Clear (f->ur_inf[i].max);
    }
    CG_IFFREE (f->ur_inf, qsnum_ur_info);
    CG_IFFREE (f->lc_inf, qsnum_lc_info);
    CG_IFFREE (f->lr_inf, qsnum_lr_info);
    CG_IFFREE (f->er_inf, qsnum_er_info);
    CG_IFFREE (f->ucindx, int);
    CG_IFFREE (f->ucrind, int);
    QSnum_FreeArray (f->uccoef, f->uc_space);
    CG_IFFREE (f->urindx, int);
    CG_IFFREE (f->urcind, int);
    QSnum_FreeArray (f->urcoef, f->ur_space);
    CG_IFFREE (f->lcindx, int);
    QSnum_FreeArray (f->lccoef, f->lc_space);
    CG_IFFREE (f->lrindx, int);
    QSnum_FreeArray (f->lrcoef, f->lr_space);
    CG_IFFREE (f->erindx, int);
    QSnum_FreeArray (f->ercoef, f->er_space);
    CG_IFFREE (f->rperm, int);
    CG_IFFREE (f->rrank, int);
    CG_IFFREE (f->cperm, int);
    CG_IFFREE (f->crank, int);
    QSnum_FreeArray (f->dmat, f->dsize);
    QSnum_svector_free (&f->xtmp);
}

int QSnum_factor_create_factor_work (qsnum_factor_work * f, int dimr, int dimc)
{
  int i=0;
  int rval;
  f->maxdim = (dimc > dimr) ? dimc : dimr; 
  int maxdim = f->maxdim;
  f->dimr       = dimr; 
  f->dimc       = dimc;
  f->etacnt = 0;
  f->work_coef = QSnum_AllocArray (maxdim);
  CG_SAFE_MALLOC (f->work_indx, maxdim, int);
  CG_SAFE_MALLOC (f->uc_inf, dimc + (f->max_k + 1), qsnum_uc_info);
  CG_SAFE_MALLOC (f->ur_inf, dimr + (f->max_k + 1), qsnum_ur_info);
  CG_SAFE_MALLOC (f->lc_inf, dimr, qsnum_lc_info);
  CG_SAFE_MALLOC (f->lr_inf, dimr, qsnum_lr_info);
  CG_SAFE_MALLOC (f->rperm, dimr, int);
  CG_SAFE_MALLOC (f->rrank, dimr, int);
  CG_SAFE_MALLOC (f->cperm, dimc, int);
  CG_SAFE_MALLOC (f->crank, dimc, int);

  for (i = dimr + f->max_k + 1; i--;)
    QSnum_Init (f->ur_inf[i].max);

  for (i = 0; i < maxdim; i++){
    QSnum_SetZero (f->work_coef[i]);
    f->work_indx[i] = 0;
  }
  for (i = 0; i < dimr; i++) {
    f->ur_inf[i].nzcnt = 0;
    f->lc_inf[i].nzcnt = 0;
    f->lr_inf[i].nzcnt = 0;
    f->rperm[i] = i;
    f->rrank[i] = i;
  }
  for (i = 0; i < dimc; i++) {
    f->uc_inf[i].nzcnt = 0;
    f->cperm[i] = i;
    f->crank[i] = i;
  }

  for (i = 0; i <= f->max_k; i++) {
    f->uc_inf[dimc + i].nzcnt = i;
    f->uc_inf[dimc + i].next = dimc + i;
    f->uc_inf[dimc + i].prev = dimc + i;
    f->ur_inf[dimr + i].nzcnt = i;
    f->ur_inf[dimr + i].next = dimr + i;
    f->ur_inf[dimr + i].prev = dimr + i;
  }

  rval = QSnum_svector_alloc (&f->xtmp, maxdim);

  CG_CLEANUP_IF (rval);

  rval = 0;

 CLEANUP:
  if (rval) {
    QSnum_factor_free_factor_work (f);
  }
  return rval;
}


static void qsnum_clear_work (qsnum_factor_work * f)
{
  int i;
  int dimr = f->dimr;
  QSnum_type *work_coef = f->work_coef;
  for (i = 0; i < dimr; i++) {
    QSnum_SetZero (work_coef[i]);
  }
}


static void qsnum_load_row (qsnum_factor_work * f, int r)
{
    QSnum_type *prow_urcoef = f->urcoef + f->ur_inf[r].rbeg;
    int *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
    int prow_nzcnt = f->ur_inf[r].nzcnt;
    QSnum_type *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    int i;
    int j;

    for (i = 0; i < prow_nzcnt; i++) {
        j = prow_urindx[i];
        QSnum_Copy (work_coef[j], prow_urcoef[i]);
        work_indx[j] = 1;
    }
}

static void qsnum_clear_row (qsnum_factor_work * f, int r)
{
    int *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
    int prow_nzcnt = f->ur_inf[r].nzcnt;
    QSnum_type *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    int i;
    int j;

    for (i = 0; i < prow_nzcnt; i++) {
        j = prow_urindx[i];
        QSnum_SetZero (work_coef[j]);
        work_indx[j] = 0;
    }
}

static int qsnum_make_ur_space (qsnum_factor_work * f, int space)
{
    QSnum_type *new_urcoef = 0;
    int *new_urindx = 0;
    int *new_urcind = 0;
    QSnum_type *urcoef = f->urcoef;
    int *urindx = f->urindx;
    int *urcind = f->urcind;
    int minspace;
    qsnum_ur_info *ur_inf = f->ur_inf;
    int dimr = f->dimr;
    int new_nzcnt = 0, old_nzcnt;
    int rbeg;
    int nzcnt;
    int i;
    int j;
    int rval;

    minspace = f->ur_space;
    nzcnt = space;
    for (i = 0; i < dimr; i++)
        nzcnt += ur_inf[i].nzcnt;
    old_nzcnt = nzcnt;
    while (nzcnt * 2 >= minspace) {
        minspace = 1+minspace*f->grow_mul;
    }

    new_urcoef = QSnum_AllocArray (minspace);
    CG_SAFE_MALLOC (new_urindx, minspace + 1, int);

    if (urcind) {
        CG_SAFE_MALLOC (new_urcind, minspace, int);
    }

    if (urcind) {
        for (j = 0; j < dimr; j++) {
            rbeg = ur_inf[j].rbeg;
            nzcnt = ur_inf[j].nzcnt;
            ur_inf[j].rbeg = new_nzcnt;
            for (i = 0; i < nzcnt; i++) {
                new_urindx[new_nzcnt] = urindx[rbeg + i];
                QSnum_Copy (new_urcoef[new_nzcnt], urcoef[rbeg + i]);
                new_urcind[new_nzcnt] = urcind[rbeg + i];
                new_nzcnt++;
            }
        }
    } else {
        for (j = 0; j < dimr; j++) {
            rbeg = ur_inf[j].rbeg;
            nzcnt = ur_inf[j].nzcnt;
            ur_inf[j].rbeg = new_nzcnt;
            for (i = 0; i < nzcnt; i++) {
                new_urindx[new_nzcnt] = urindx[rbeg + i];
                QSnum_Copy (new_urcoef[new_nzcnt], urcoef[rbeg + i]);
                new_nzcnt++;
            }
        }
    }

    for (i = new_nzcnt; i < minspace; i++) {
        new_urindx[i] = -1;
    }
    new_urindx[minspace] = 0;
    QSnum_FreeArray (f->urcoef, f->ur_space);
    f->urcoef = new_urcoef;
    new_urcoef = 0;

    CG_IFFREE (f->urindx, int);
    f->urindx = new_urindx;
    new_urindx = 0;

    CG_IFFREE (f->urcind, int);
    f->urcind = new_urcind;
    new_urcind = 0;

    f->ur_freebeg = new_nzcnt;
    f->ur_space = minspace;

    rval = 0;

CLEANUP:
    CG_IFFREE (new_urcoef, QSnum_type);
    CG_IFFREE (new_urindx, int);
    CG_IFFREE (new_urcind, int);
    return rval;
}

static int qsnum_make_uc_space (qsnum_factor_work * f, int space)
{
    QSnum_type *new_uccoef = 0;
    int *new_ucindx = 0;
    int *new_ucrind = 0;
    int uc_freebeg = f->uc_freebeg;
    QSnum_type *uccoef = f->uccoef;
    int *ucindx = f->ucindx;
    int *ucrind = f->ucrind;
    int minspace = uc_freebeg + space;
    qsnum_uc_info *uc_inf = f->uc_inf;
    int dimc = f->dimc;    
int new_nzcnt = 0;
    int cbeg;
    int nzcnt;
    int i;
    int j;
    int rval;

   if (f->uc_space * f->grow_mul > minspace) {
        minspace = f->uc_space * f->grow_mul;
    }

    CG_SAFE_MALLOC (new_ucindx, minspace + 1, int);

    if (ucrind) {
        new_uccoef = QSnum_AllocArray (minspace);
        CG_SAFE_MALLOC (new_ucrind, minspace, int);
    }

    if (ucrind) {
        for (j = 0; j < dimc; j++) {
            cbeg = uc_inf[j].cbeg;
            nzcnt = uc_inf[j].nzcnt;
            uc_inf[j].cbeg = new_nzcnt;
            for (i = 0; i < nzcnt; i++) {
                new_ucindx[new_nzcnt] = ucindx[cbeg + i];
                QSnum_Copy (new_uccoef[new_nzcnt], uccoef[cbeg + i]);
                new_ucrind[new_nzcnt] = ucrind[cbeg + i];
                new_nzcnt++;
            }
        }
    } else {
        for (j = 0; j < dimc; j++) {
            cbeg = uc_inf[j].cbeg;
            nzcnt = uc_inf[j].nzcnt;
            uc_inf[j].cbeg = new_nzcnt;
            for (i = 0; i < nzcnt; i++) {
                new_ucindx[new_nzcnt] = ucindx[cbeg + i];
                new_nzcnt++;
            }
        }
    }

    for (i = new_nzcnt; i < minspace; i++) {
        new_ucindx[i] = -1;
    }
    new_ucindx[minspace] = 0;

    QSnum_FreeArray (f->uccoef, f->uc_space);
    f->uccoef = new_uccoef;
    new_uccoef = 0;

    CG_IFFREE (f->ucindx, int);
    f->ucindx = new_ucindx;
    new_ucindx = 0;

    CG_IFFREE (f->ucrind, int);
    f->ucrind = new_ucrind;
    new_ucrind = 0;

    f->uc_freebeg = new_nzcnt;
    f->uc_space = minspace;

    rval = 0;

CLEANUP:
    CG_IFFREE (new_uccoef, QSnum_type);
    CG_IFFREE (new_ucindx, int);
    CG_IFFREE (new_ucrind, int);
    return rval; 
}

static int qsnum_make_lc_space (qsnum_factor_work * f, int space)
{
    QSnum_type *new_lccoef = 0;
    int *new_lcindx = 0;
    int lc_freebeg = f->lc_freebeg;
    QSnum_type *lccoef = f->lccoef;
    int *lcindx = f->lcindx;
    int minspace = lc_freebeg + space;
    int i;
    int rval;
    if (f->lc_space * f->grow_mul > minspace) {
        minspace = f->lc_space * f->grow_mul;
    }

    new_lccoef = QSnum_AllocArray (minspace);
    CG_SAFE_MALLOC (new_lcindx, minspace, int);

    for (i = 0; i < lc_freebeg; i++) {
        QSnum_Copy (new_lccoef[i], lccoef[i]);
        new_lcindx[i] = lcindx[i];
    }

    QSnum_FreeArray (lccoef, f->lc_space);
    f->lccoef = new_lccoef;
    new_lccoef = 0;

    CG_IFFREE (lcindx, int);
    f->lcindx = new_lcindx;
    new_lcindx = 0;

    f->lc_space = minspace;

    rval = 0;

CLEANUP:
    CG_IFFREE (new_lccoef, QSnum_type);
    CG_IFFREE (new_lcindx, int);
    return rval;
}

static void qsnum_set_col_nz (qsnum_factor_work * f, int c)
{
  qsnum_uc_info *uc_inf = f->uc_inf;
  int nzcnt = uc_inf[c].nzcnt;
  int max_k = f->max_k;
  int dimc = f->dimc;  
  if (uc_inf[c].next >= 0) {
    uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
    uc_inf[uc_inf[c].prev].next = uc_inf[c].next;

    if (nzcnt >= max_k)
      nzcnt = max_k;
    uc_inf[c].next = uc_inf[dimc + nzcnt].next;
    uc_inf[c].prev = dimc + nzcnt;
    uc_inf[dimc + nzcnt].next = c;
    uc_inf[uc_inf[c].next].prev = c;
  }
}

static void qsnum_set_row_nz (qsnum_factor_work * f, int r)
{
  qsnum_ur_info *ur_inf = f->ur_inf;
  int nzcnt = ur_inf[r].pivcnt;
  int max_k = f->max_k;
  int dimr = f->dimr;
  if (ur_inf[r].next >= 0) {
    ur_inf[ur_inf[r].next].prev = ur_inf[r].prev;
    ur_inf[ur_inf[r].prev].next = ur_inf[r].next;

    if (nzcnt >= max_k)
      nzcnt = max_k;
    ur_inf[r].next = ur_inf[dimr + nzcnt].next;
    ur_inf[r].prev = dimr + nzcnt;
    ur_inf[dimr + nzcnt].next = r;
    ur_inf[ur_inf[r].next].prev = r;
  }
}

static void qsnum_remove_col_nz (qsnum_factor_work * f, int r, int c)
{
    qsnum_uc_info *uc_inf = f->uc_inf;
    int *ucindx = f->ucindx + uc_inf[c].cbeg;
    int nzcnt = uc_inf[c].nzcnt;
    int i;

    for (i = 0; i < nzcnt; i++) {
        if (ucindx[i] == r) {
            --nzcnt;
            ucindx[i] = ucindx[nzcnt];
            ucindx[nzcnt] = -1;
            break;
        }
    }
    uc_inf[c].nzcnt = nzcnt;

    qsnum_set_col_nz (f, c);
}

static void qsnum_remove_row_nz (qsnum_factor_work * f, int r, int c)
{
    qsnum_ur_info *ur_inf = f->ur_inf;
    int *urindx = f->urindx + ur_inf[r].rbeg;
    QSnum_type *urcoef = f->urcoef + ur_inf[r].rbeg;
    int pivcnt = ur_inf[r].pivcnt;
    QSnum_type max;
    int tind;
    QSnum_type tcoef;
    int i;
    QSnum_Init (tcoef);
    QSnum_Init (max);
    QSnum_SetZero (max);

    for (i = 0; i < pivcnt; i++) {
        if (urindx[i] == c) {
            --pivcnt;
            CG_SWAP (urindx[i], urindx[pivcnt], tind);
            QSnum_SWAP (urcoef[i], urcoef[pivcnt], tcoef);
            --i;
        } else {
            QSnum_CopyMaxAbs (max, urcoef[i]);
        }
    }
    ur_inf[r].pivcnt = pivcnt;
    QSnum_Copy (ur_inf[r].max, max);
    qsnum_set_row_nz (f, r);
    QSnum_Clear (max);
    QSnum_Clear (tcoef);
}

static int qsnum_add_col_nz (qsnum_factor_work * f, int r, int c)
{
    qsnum_uc_info *uc_inf = f->uc_inf;
    int cbeg = uc_inf[c].cbeg;
    int nzcnt = uc_inf[c].nzcnt;
    int uc_freebeg = f->uc_freebeg;
    int *ucindx = f->ucindx;
    int i;
    int rval = 0;

    if (uc_inf[c].next == -1) return 0;

    if (ucindx[cbeg + nzcnt] == -1) {
        ucindx[cbeg + nzcnt] = r;
        uc_inf[c].nzcnt++;
        if (nzcnt + cbeg == uc_freebeg) {
            f->uc_freebeg = uc_freebeg + 1;
        }
    } else {
        if (uc_freebeg + nzcnt + 1 >= f->uc_space) {
            rval = qsnum_make_uc_space (f, nzcnt + 1);
            CG_CLEANUP_IF (rval);
            uc_freebeg = f->uc_freebeg;
            cbeg = uc_inf[c].cbeg;
            ucindx = f->ucindx;
        }
        for (i = 0; i < nzcnt; i++) {
            ucindx[uc_freebeg + i] = ucindx[cbeg + i];
            ucindx[cbeg + i] = -1;
        }
        ucindx[uc_freebeg + nzcnt] = r;
        uc_inf[c].cbeg = uc_freebeg;
        uc_inf[c].nzcnt++;
        f->uc_freebeg = uc_freebeg + nzcnt + 1;
    }

    qsnum_set_col_nz (f, c);
CLEANUP:
    return rval;
}

static void qsnum_disable_col (qsnum_factor_work * f, int c)
{
    qsnum_uc_info *uc_inf = f->uc_inf;

    if (uc_inf[c].next >= 0) {
        uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
        uc_inf[uc_inf[c].prev].next = uc_inf[c].next;

        uc_inf[c].next = -2;
        uc_inf[c].prev = -2;
    }
}

static void qsnum_remove_col (qsnum_factor_work * f, int c)
{
    qsnum_uc_info *uc_inf = f->uc_inf;
    int cbeg = uc_inf[c].cbeg;
    int nzcnt = uc_inf[c].nzcnt;
    int *ucindx = f->ucindx;
    int i;

    for (i = 0; i < nzcnt; i++) {
        ucindx[cbeg + i] = -1;
    }
    uc_inf[c].cbeg = 0;
    uc_inf[c].nzcnt = 0;

    if (uc_inf[c].next >= 0) {
        uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
        uc_inf[uc_inf[c].prev].next = uc_inf[c].next;

        uc_inf[c].next = -1;
        uc_inf[c].prev = -1;
    }
}

static void qsnum_remove_row (qsnum_factor_work * f, int r)
{
    qsnum_ur_info *ur_inf = f->ur_inf;

    if (ur_inf[r].next >= 0) {
        ur_inf[ur_inf[r].next].prev = ur_inf[r].prev;
        ur_inf[ur_inf[r].prev].next = ur_inf[r].next;
        ur_inf[r].next = -1;
        ur_inf[r].prev = -1;
    }
}

static void qsnum_find_coef (qsnum_factor_work * f, int r, int c, QSnum_type * coef)
{
    QSnum_type *prow_urcoef = f->urcoef + f->ur_inf[r].rbeg;
    int *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
    int i;
    int prow_nzcnt = f->ur_inf[r].nzcnt;
    QSnum_SetZero (*coef);
    for (i = 0; i < prow_nzcnt; i++) {
        if (prow_urindx[i] == c) {
            QSnum_Copy (*coef, prow_urcoef[i]);
            return;
        }
    }
    fprintf (stderr, "Coefficient not found\n");
    return;
}

static int qsnum_elim_row (qsnum_factor_work * f, int elim_r, int r, int c,
        QSnum_type * p_pivot_coef)
{
    qsnum_ur_info *ur_inf = f->ur_inf;
    QSnum_type *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    QSnum_type *urcoef = f->urcoef;
    int *urindx = f->urindx;
    int prow_beg = ur_inf[r].rbeg;
    int prow_nzcnt = ur_inf[r].nzcnt;
    int prow_pivcnt = ur_inf[r].pivcnt;
    int fill = ur_inf[elim_r].nzcnt;
    int cancel = 0;
    QSnum_type max;
    int erow_beg;
    int erow_nzcnt;
    int erow_pivcnt;
    QSnum_type x;
    int i;
    int j;
    int rval = 0;
    QSnum_type elim_coef;
    QSnum_Init (max);
    QSnum_Init (x);
    QSnum_Init (elim_coef);
    QSnum_SetZero (max);
    qsnum_find_coef (f, r, c, &elim_coef);
    QSnum_CopyDiv (elim_coef, work_coef[c]);
    QSnum_Copy (*p_pivot_coef, elim_coef);

    for (i = 0; i < prow_nzcnt; i++) {
        j = urindx[prow_beg + i];
        if (work_indx[j] == 1) {
            QSnum_Copy (x, urcoef[prow_beg + i]);
            QSnum_CopySubProd (x, elim_coef, work_coef[j]);
            if (!QSnum_NeqZeroTol (x, f->fzero_tol) || j == c) {
                cancel++;
                if (j != c) {
                    qsnum_remove_col_nz (f, r, j);
                }
                if (i < prow_pivcnt) {
                    prow_pivcnt--;
                    prow_nzcnt--;
                    urindx[prow_beg + i] = urindx[prow_beg + prow_pivcnt];
                    QSnum_Copy (urcoef[prow_beg + i], urcoef[prow_beg + prow_pivcnt]);
                    if (prow_pivcnt != prow_nzcnt) {
                        urindx[prow_beg + prow_pivcnt] = urindx[prow_beg + prow_nzcnt];
                        QSnum_Copy (urcoef[prow_beg + prow_pivcnt],
                                                 urcoef[prow_beg + prow_nzcnt]);
                    }
                } else {
                    prow_nzcnt--;
                    urindx[prow_beg + i] = urindx[prow_beg + prow_nzcnt];
                    QSnum_Copy (urcoef[prow_beg + i], urcoef[prow_beg + prow_nzcnt]);
                }
                urindx[prow_beg + prow_nzcnt] = -1;
                i--;
            } else {
                QSnum_Copy (urcoef[prow_beg + i], x);
                if (i < prow_pivcnt) {
                    QSnum_CopyMaxAbs (max, x);
                }
            }
            work_indx[j] = 0;
            fill--;
        } else {
            if (i < prow_pivcnt) {
                QSnum_CopyMaxAbs (max, urcoef[prow_beg + i]);
            }
        }
    }

    if (fill > 0) {
        ur_inf[r].nzcnt = prow_nzcnt;
        ur_inf[r].pivcnt = prow_pivcnt;
        if (fill > cancel) {
            int ur_freebeg = f->ur_freebeg;

            if (ur_freebeg + prow_nzcnt + fill >= f->ur_space) {
                rval = qsnum_make_ur_space (f, prow_nzcnt + fill);
                CG_CLEANUP_IF (rval);
                urcoef = f->urcoef;
                urindx = f->urindx;
                ur_freebeg = f->ur_freebeg;
                prow_beg = f->ur_inf[r].rbeg;
            }
            for (i = 0; i < prow_nzcnt; i++) {
                urindx[ur_freebeg + i] = urindx[prow_beg + i];
                QSnum_Copy (urcoef[ur_freebeg + i], urcoef[prow_beg + i]);
                urindx[prow_beg + i] = -1;
            }
            ur_inf[r].rbeg = ur_freebeg;
            f->ur_freebeg = ur_freebeg + prow_nzcnt + fill;
            prow_beg = ur_freebeg;
        }

        erow_beg = ur_inf[elim_r].rbeg;
        erow_nzcnt = ur_inf[elim_r].nzcnt;
        erow_pivcnt = ur_inf[elim_r].pivcnt;

        for (i = 0; i < erow_pivcnt; i++) {
            j = urindx[erow_beg + i];
            if (work_indx[j] == 1) {
                QSnum_CopyNeg (x, elim_coef);
                QSnum_CopyMult (x, urcoef[erow_beg + i]);
                if (QSnum_NeqZeroTol (x, f->fzero_tol)) {
                    rval = qsnum_add_col_nz (f, r, j);
                    CG_CLEANUP_IF (rval);
                    if (prow_pivcnt != prow_nzcnt) {
                        urindx[prow_beg + prow_nzcnt] = urindx[prow_beg + prow_pivcnt];
                        QSnum_Copy (urcoef[prow_beg + prow_nzcnt],
                                                 urcoef[prow_beg + prow_pivcnt]);
                    }
                    urindx[prow_beg + prow_pivcnt] = j;
                    QSnum_Copy (urcoef[prow_beg + prow_pivcnt], x);
                    QSnum_CopyMaxAbs (max, x);
                    prow_pivcnt++;
                    prow_nzcnt++;
                }
            } else {
                work_indx[j] = 1;
            }
        }
        for (i = erow_pivcnt; i < erow_nzcnt; i++) {
            j = urindx[erow_beg + i];
            if (work_indx[j] == 1) {
                QSnum_CopyNeg (x, elim_coef);
                QSnum_CopyMult (x, urcoef[erow_beg + i]);
                if (QSnum_NeqZeroTol (x, f->fzero_tol)) {
                    rval = qsnum_add_col_nz (f, r, j);
                    CG_CLEANUP_IF (rval);
                    urindx[prow_beg + prow_nzcnt] = j;
                    QSnum_Copy (urcoef[prow_beg + prow_nzcnt], x);
                    prow_nzcnt++;
                }
            } else {
                work_indx[j] = 1;
            }
        }
    } else {
        erow_nzcnt = ur_inf[elim_r].nzcnt;
        erow_beg = ur_inf[elim_r].rbeg;
        for (i = 0; i < erow_nzcnt; i++) {
            j = urindx[erow_beg + i];
            work_indx[j] = 1;
        }
    }

    ur_inf[r].nzcnt = prow_nzcnt;
    ur_inf[r].pivcnt = prow_pivcnt;
    QSnum_Copy (ur_inf[r].max, max);

    qsnum_set_row_nz (f, r);
CLEANUP:
    QSnum_Clear (elim_coef);
    QSnum_Clear (x);
    QSnum_Clear (max);
    return rval;
}

#define qsnum_SETPERM(f,s,r,c) {                \
        f->rperm[f->rrank[r]] = f->rperm[s];  \
        f->rrank[f->rperm[s]] = f->rrank[r];  \
        f->rperm[s] = r;                      \
        f->rrank[r] = s;                      \
                                              \
        f->cperm[f->crank[c]] = f->cperm[s];  \
        f->crank[f->cperm[s]] = f->crank[c];  \
        f->cperm[s] = c;                      \
        f->crank[c] = s;                      \
}

static int qsnum_elim (qsnum_factor_work * f, int r, int c)
{
    qsnum_uc_info *uc_inf = f->uc_inf;
    qsnum_ur_info *ur_inf = f->ur_inf;
    qsnum_lc_info *lc_inf = f->lc_inf;
    int *urindx;
    int *ucindx;
    int *lcindx;
    QSnum_type *urcoef;
    QSnum_type *lccoef;
    QSnum_type pivot_coef;
    int nzcnt;
    int lc_freebeg;
    int s = f->stage;
    int i;
    int j;
    int rval = 0;
    QSnum_Init (pivot_coef);

    if (uc_inf[c].nzcnt == 1) {
        /* col singleton */
        qsnum_SETPERM (f, s, r, c);

        lc_inf[s].cbeg = -1;
        lc_inf[s].c = r;
        lc_inf[s].nzcnt = 0;
        f->stage++;

        urindx = f->urindx + ur_inf[r].rbeg;
        urcoef = f->urcoef + ur_inf[r].rbeg;
        nzcnt = ur_inf[r].nzcnt;
        for (i = 0; i < nzcnt; i++) {
            j = urindx[i];
            qsnum_remove_col_nz (f, r, j);
            if (j == c) {
                urindx[i] = urindx[0];
                urindx[0] = c;
                QSnum_SWAP (urcoef[0], urcoef[i], pivot_coef);
            }
        }
        qsnum_remove_row (f, r);
        qsnum_remove_col (f, c);
    } else if (ur_inf[r].nzcnt == 1) {
        /* row singleton */
        --(f->nstages);
        qsnum_SETPERM (f, f->nstages, r, c);

        lc_inf[f->nstages].cbeg = -1;
        lc_inf[f->nstages].c = r;
        lc_inf[f->nstages].nzcnt = 0;

        ucindx = f->ucindx + uc_inf[c].cbeg;
        nzcnt = uc_inf[c].nzcnt;
        for (i = 0; i < nzcnt; i++) {
            j = ucindx[i];
            qsnum_remove_row_nz (f, j, c);
        }
        qsnum_remove_row (f, r);
        qsnum_remove_col (f, c);
    } else {
        qsnum_SETPERM (f, s, r, c);
        f->stage++;

        nzcnt = uc_inf[c].nzcnt;
        if (f->lc_freebeg + nzcnt >= f->lc_space)
        {
            rval = qsnum_make_lc_space (f, nzcnt);
            CG_CLEANUP_IF (rval);
        }
        lc_freebeg = f->lc_freebeg;
        lc_inf[s].cbeg = lc_freebeg;
        lc_inf[s].c = r;
        lcindx = f->lcindx;
        lccoef = f->lccoef;
        qsnum_load_row (f, r);
        ucindx = f->ucindx + uc_inf[c].cbeg;
        for (i = 0; i < nzcnt; i++) {
            j = f->ucindx[uc_inf[c].cbeg + i];
            if (j != r) {
                rval = qsnum_elim_row (f, r, j, c, &pivot_coef);
                CG_CLEANUP_IF (rval);
                lcindx[lc_freebeg] = j;
                QSnum_Copy (lccoef[lc_freebeg], pivot_coef);
                lc_freebeg++;
            }
        }
        lc_inf[s].nzcnt = lc_freebeg - lc_inf[s].cbeg;
        f->lc_freebeg = lc_freebeg;

        qsnum_clear_row (f, r);

        urindx = f->urindx + ur_inf[r].rbeg;
        urcoef = f->urcoef + ur_inf[r].rbeg;
        nzcnt = ur_inf[r].nzcnt;
        for (i = 0; i < nzcnt; i++) {
            j = urindx[i];
            qsnum_remove_col_nz (f, r, j);
            if (j == c) {
                urindx[i] = urindx[0];
                urindx[0] = c;
                QSnum_SWAP (urcoef[0], urcoef[i], pivot_coef);
            }
        }
        qsnum_remove_row (f, r);
        qsnum_remove_col (f, c);
    }
CLEANUP:
    QSnum_Clear (pivot_coef);
    return rval;
}

static void qsnum_find_pivot_column (qsnum_factor_work * f, int c, int *p_r)
{
    qsnum_uc_info *uc_inf = f->uc_inf;
    qsnum_ur_info *ur_inf = f->ur_inf;
    int *ucindx = f->ucindx;
    int nzcnt = uc_inf[c].nzcnt;
    int cbeg = uc_inf[c].cbeg;
    QSnum_type num_tmp[2];
    int bestnz = -1;
    int i;
    int r;
    QSnum_Init (num_tmp[0]);
    QSnum_Init (num_tmp[1]);
   *p_r = -1;
    for (i = 0; i < nzcnt; i++) {
        r = ucindx[cbeg + i];
        if ((bestnz == -1 || ur_inf[r].pivcnt < bestnz)) {
	    qsnum_find_coef (f, r, c, &num_tmp[0]);
            if (QSnum_Less (num_tmp[0], qsnum_zeroLpNum))
                QSnum_Sign (num_tmp[0]);
            QSnum_Copy (num_tmp[1], f->partial_cur);
            QSnum_CopyMult (num_tmp[1], ur_inf[r].max);
            if (QSnum_Leq (num_tmp[1], num_tmp[0])) {
                bestnz = ur_inf[r].pivcnt;
                *p_r = r;
            }
        }
    }
    QSnum_Clear (num_tmp[0]);
    QSnum_Clear (num_tmp[1]);
}

static void qsnum_find_pivot_row (qsnum_factor_work * f, int r, int *p_c)
{
    qsnum_uc_info *uc_inf = f->uc_inf;
    qsnum_ur_info *ur_inf = f->ur_inf;
    int *urindx = f->urindx;
    QSnum_type *urcoef = f->urcoef;
    int pivcnt = ur_inf[r].pivcnt;
    int rbeg = ur_inf[r].rbeg;
    QSnum_type thresh[2];
    int bestnz = -1;
    int i = 0;
    int c;
    QSnum_Init (thresh[0]);
    QSnum_Init (thresh[1]);
    QSnum_Copy (thresh[0], f->partial_cur);
    QSnum_CopyMult (thresh[0], ur_inf[r].max);
    QSnum_CopyAbs (thresh[1], urcoef[rbeg + i]);
    *p_c = -1;
    for (i = 0; i < pivcnt; i++) {
        c = urindx[rbeg + i];
        if ((bestnz == -1 || uc_inf[c].nzcnt < bestnz)) {
            /* DAVID FIX */
            QSnum_CopyAbs (thresh[1], urcoef[rbeg + i]);
	    if (QSnum_Leq (thresh[0], thresh[1])) {
                bestnz = uc_inf[c].nzcnt;
                *p_c = c;
            }
        }
    }
    QSnum_Clear (thresh[0]);
    QSnum_Clear (thresh[1]);
}

static int qsnum_find_pivot (qsnum_factor_work * f, int *p_r, int *p_c)
{
    qsnum_uc_info *uc_inf = f->uc_inf;
    qsnum_ur_info *ur_inf = f->ur_inf;
    int dimr = f->dimr;
    int dimc = f->dimc;
    int max_k = f->max_k;
    int p = f->p;
    int c;
    int r;
    int mm = 0;
    int n = 0;
    int m;
    int k = 2;

    if (uc_inf[dimc + 1].next != dimc + 1) {
        c = uc_inf[dimc + 1].next;
        r = f->ucindx[uc_inf[c].cbeg];
        *p_c = c;
        *p_r = r;
        return 0;
    } else if (ur_inf[dimr + 1].next != dimr + 1) {
        r = ur_inf[dimr + 1].next;
        c = f->urindx[ur_inf[r].rbeg];
        *p_c = c;
        *p_r = r;
        return 0;
    }
    *p_r = -1;
    *p_c = -1;
    for (; k <= max_k && (mm == 0 || mm > (k - 1) * (k - 1)); k++) {
        if (uc_inf[dimc + k].next != dimc + k) {
            for (c = uc_inf[dimc + k].next; c != dimc + k; c = uc_inf[c].next) {
                qsnum_find_pivot_column (f, c, &r);
                if (r >= 0) {
                    m = (uc_inf[c].nzcnt - 1) * (ur_inf[r].pivcnt - 1);
                    if (mm == 0 || m < mm) {
                        mm = m;
                        *p_c = c;
                        *p_r = r;
                        if (mm <= (k - 1) * (k - 1)) {
                            return 0;
                        }
                    }
                } else {
                    c = uc_inf[c].prev;
                    qsnum_disable_col (f, uc_inf[c].next);
                }
                n++;
                if (n >= p && mm != 0) {
                    return 0;
                }
            }
        }

        if (ur_inf[dimr + k].next != dimr + k) {
            for (r = ur_inf[dimr + k].next; r != dimr + k; r = ur_inf[r].next) {
                qsnum_find_pivot_row (f, r, &c);
                if (c >= 0) {
                    m = (uc_inf[c].nzcnt - 1) * (ur_inf[r].pivcnt - 1);
                    if (mm == 0 || m < mm) {
                        mm = m;
                        *p_c = c;
                        *p_r = r;
                        if (mm <= k * (k - 1)) {
                            return 0;
                        }
                    }
                }
                n++;
                if (n >= p && mm != 0) {
                    return 0;
                }
            }
        }
    }
    if (mm != 0) {
        return 0;
    } else {
        fprintf (stderr, "No acceptable pivot found!\n");
        return E_NO_PIVOT;
    }
}

static int qsnum_create_factor_space (qsnum_factor_work * f)
{
  qsnum_uc_info *uc_inf = f->uc_inf;
  qsnum_ur_info *ur_inf = f->ur_inf;
  int dimr = f->dimr;
  int dimc = f->dimc;
  int nzcnt;
  int i;
  int rval;

  nzcnt = 0;
  for (i = 0; i < dimr; i++) {
    nzcnt += ur_inf[i].nzcnt;
  }

  if (f->ucindx == 0) {
    f->uc_space = nzcnt * f->uc_space_mul;
    CG_SAFE_MALLOC (f->ucindx, f->uc_space + 1, int);
  }

  if (f->urindx == 0 || f->urcoef == 0) {
    CG_IFFREE (f->urindx, int);
    QSnum_FreeArray (f->urcoef, f->ur_space);
    f->ur_space = nzcnt * f->ur_space_mul;
    CG_SAFE_MALLOC (f->urindx, f->ur_space + 1, int);
    f->urcoef = QSnum_AllocArray (f->ur_space);
  }

  if (f->lcindx == 0 || f->lccoef == 0) {
    CG_IFFREE (f->lcindx, int);
    QSnum_FreeArray (f->lccoef, f->lc_space);
    f->lc_space = nzcnt * f->lc_space_mul;
    CG_SAFE_MALLOC (f->lcindx, f->lc_space, int);
    f->lccoef = QSnum_AllocArray (f->lc_space);
  }

  nzcnt = 0;
  for (i = 0; i < dimr; i++) {
    ur_inf[i].rbeg = nzcnt;
    nzcnt += ur_inf[i].nzcnt;
    ur_inf[i].nzcnt = ur_inf[i].rbeg;
  }
  f->ur_freebeg = nzcnt;

  nzcnt = 0;
  for (i = 0; i < dimc; i++) {
    uc_inf[i].cbeg = nzcnt;
    nzcnt += uc_inf[i].nzcnt;
    uc_inf[i].nzcnt = uc_inf[i].cbeg;
  }
  f->uc_freebeg = nzcnt;

  f->lc_freebeg = 0;

  rval = 0;
 CLEANUP:
  return rval;
}

static int qsnum_init_matrix (qsnum_factor_work * f, int *basis, int *cbeg,
			      int *clen, int *in_ucindx, QSnum_type * in_uccoef)
{
  qsnum_uc_info *uc_inf = f->uc_inf;
  qsnum_ur_info *ur_inf = f->ur_inf;
  int dimr = f->dimr;
  int dimc = f->dimc;
  int max_k = f->max_k;
  int *ucindx;
  int *urindx;
  QSnum_type *urcoef;
  int nzcnt;
  int beg;
  int i;
  int j;
  int r;
  int rval = 0;
  QSnum_type v;
  QSnum_type max;
  QSnum_Init (v);
  QSnum_Init (max);

  for (i = 0; i < dimr; i++) {
    ur_inf[i].nzcnt = 0;
  }
  for (i = 0; i < dimc; i++) {
    nzcnt = clen[basis[i]];
    beg = cbeg[basis[i]];
    uc_inf[i].nzcnt = nzcnt;
    for (j = 0; j < nzcnt; j++) {
      r = in_ucindx[beg + j];
      ur_inf[r].nzcnt++;
    }
  }

  rval = qsnum_create_factor_space (f);
  CG_CLEANUP_IF (rval);

  urindx = f->urindx;
  ucindx = f->ucindx;
  urcoef = f->urcoef;

  for (i = 0; i < dimc; i++) {
    nzcnt = clen[basis[i]];
    beg = cbeg[basis[i]];
    for (j = 0; j < nzcnt; j++) {
      QSnum_Copy (v, in_uccoef[beg + j]);
      if (!QSnum_NeqZeroTol (v, f->fzero_tol))
	continue;
      r = in_ucindx[beg + j];
      ucindx[uc_inf[i].nzcnt++] = r;
      urindx[ur_inf[r].nzcnt] = i;
      QSnum_Copy (urcoef[ur_inf[r].nzcnt], v);
      ur_inf[r].nzcnt++;
    }
  }

  for (i = 0; i < dimc; i++) 
    uc_inf[i].nzcnt -= uc_inf[i].cbeg;
  for (i = 0; i < dimr; i++) 
    ur_inf[i].nzcnt -= ur_inf[i].rbeg;
    

  j = f->uc_space;
  for (i = f->uc_freebeg; i < j; i++) {
    ucindx[i] = -1;
  }
  ucindx[j] = 0;

  j = f->ur_space;
  for (i = f->ur_freebeg; i < j; i++) {
    urindx[i] = -1;
  }
  urindx[j] = 0;

  for (i = 0; i < dimr; i++) {
    nzcnt = ur_inf[i].nzcnt;
    ur_inf[i].pivcnt = nzcnt;
    beg = ur_inf[i].rbeg;
    QSnum_SetZero (max);
    for (j = 0; j < nzcnt; j++) {
      QSnum_CopyMaxAbs (max, urcoef[beg + j]);
    }
    QSnum_Copy (ur_inf[i].max, max);
  }

  for (i = 0; i <= max_k; i++) {
    ur_inf[dimr + i].next = dimr + i;
    ur_inf[dimr + i].prev = dimr + i;
    uc_inf[dimc + i].next = dimc + i;
    uc_inf[dimc + i].prev = dimc + i;
  }

  for (i = 0; i < dimc; i++) {
    nzcnt = uc_inf[i].nzcnt;
    if (nzcnt >= max_k)
      nzcnt = max_k;
    uc_inf[i].next = uc_inf[dimc + nzcnt].next;
    uc_inf[i].prev = dimc + nzcnt;
    uc_inf[dimc + nzcnt].next = i;
    uc_inf[uc_inf[i].next].prev = i;
  }
  for (i = 0; i < dimr; i++) {
    nzcnt = ur_inf[i].pivcnt;
    if (nzcnt >= max_k)
      nzcnt = max_k;
    ur_inf[i].next = ur_inf[dimr + nzcnt].next;
    ur_inf[i].prev = dimr + nzcnt;
    ur_inf[dimr + nzcnt].next = i;
    ur_inf[ur_inf[i].next].prev = i;
  }

  /* sentinal for column space */
  ucindx[f->uc_space] = 0;

  qsnum_clear_work (f);

 CLEANUP:
  QSnum_Clear (max);
  QSnum_Clear (v);
  return rval;
}










static int qsnum_build_iteration_u_data (qsnum_factor_work * f)
{
    int dimr = f->dimr;
    int dimc = f->dimc;
    qsnum_ur_info *ur_inf = f->ur_inf;
    qsnum_uc_info *uc_inf = f->uc_inf;
    QSnum_type *uccoef = 0;
    int *ucindx = 0;
    int *urindx = f->urindx;
    QSnum_type *urcoef = f->urcoef;
    int *ucrind = 0;
    int *urcind = 0;
    int nzcnt;
    int beg;
    int cbeg;
    int cnzcnt;
    int uc_space = f->uc_space;
    int er_space;
    int i;
    int j;
    int k;
    int rval;
    nzcnt = 0;
    for (i = 0; i < dimr; i++) {
        nzcnt += ur_inf[i].nzcnt;
    }

    QSnum_FreeArray (f->uccoef, f->uc_space);
    uccoef = QSnum_AllocArray (nzcnt);
    f->uccoef = uccoef;

    CG_IFFREE (f->ucrind, int);
    CG_SAFE_MALLOC (ucrind, nzcnt, int);
    f->ucrind = ucrind;

    CG_IFFREE (f->urcind, int);
    CG_SAFE_MALLOC (urcind, f->ur_space, int);
    f->urcind = urcind;

    if (uc_space < nzcnt) {
        CG_IFFREE (f->ucindx, int);
        CG_SAFE_MALLOC (f->ucindx, nzcnt + 1, int);
    }
    f->uc_space = nzcnt;
    uc_space = nzcnt;
    ucindx = f->ucindx;

    for (i = 0; i < dimc; i++) {
        uc_inf[i].nzcnt = 0;
    }

    for (i = 0; i < dimr; i++) {
        nzcnt = ur_inf[i].nzcnt;
        beg = ur_inf[i].rbeg;
        for (j = 0; j < nzcnt; j++) {
            uc_inf[urindx[beg + j]].nzcnt++;
        }
        ur_inf[i].delay = 0;
    }

    nzcnt = 0;
    for (i = 0; i < dimc; i++) {
        uc_inf[i].cbeg = nzcnt;
        nzcnt += uc_inf[i].nzcnt;
        uc_inf[i].nzcnt = 0;
        uc_inf[i].delay = 0;
    }

    f->uc_freebeg = nzcnt;
    for (i = nzcnt; i < uc_space; i++) {
        ucindx[i] = -1;
    }
    ucindx[uc_space] = 0;

    for (i = 0; i < dimr; i++) {
        nzcnt = ur_inf[i].nzcnt;
        beg = ur_inf[i].rbeg;
        k = urindx[beg];
        cbeg = uc_inf[k].cbeg;
        cnzcnt = uc_inf[k].nzcnt;
        if (cnzcnt != 0) {
            ucindx[cbeg + cnzcnt] = ucindx[cbeg];
            QSnum_Copy (uccoef[cbeg + cnzcnt], uccoef[cbeg]);
            ucrind[cbeg + cnzcnt] = ucrind[cbeg];
            urcind[ur_inf[ucindx[cbeg]].rbeg + ucrind[cbeg]] = cnzcnt;
        }
        ucindx[cbeg] = i;
        QSnum_Copy (uccoef[cbeg], urcoef[beg]);
        ucrind[cbeg] = 0;
        urcind[beg] = 0;
        uc_inf[k].nzcnt = cnzcnt + 1;
        for (j = 1; j < nzcnt; j++) {
            k = urindx[beg + j];
            cbeg = uc_inf[k].cbeg;
            cnzcnt = uc_inf[k].nzcnt;
            ucindx[cbeg + cnzcnt] = i;
            QSnum_Copy (uccoef[cbeg + cnzcnt], urcoef[beg + j]);
            ucrind[cbeg + cnzcnt] = j;
            urcind[beg + j] = cnzcnt;
            uc_inf[k].nzcnt++;
        }
    }

    for (i = 0; i < dimr; i++) {
        f->rrank[f->rperm[i]] = i;
    }

    nzcnt = f->ur_space;

    for (i = f->ur_freebeg; i < nzcnt; i++) {
        urindx[i] = -1;
    }
    urindx[nzcnt] = 0;

    qsnum_clear_work (f);

    er_space = f->er_space_mul * f->etamax;
    CG_SAFE_MALLOC (f->er_inf, f->etamax, qsnum_er_info);
    CG_SAFE_MALLOC (f->erindx, er_space, int);
    f->ercoef = QSnum_AllocArray (er_space);
    f->etacnt = 0;
    f->er_freebeg = 0;
    f->er_space = er_space;

    rval = 0;

CLEANUP:
    return rval;
}

static int qsnum_build_iteration_l_data (qsnum_factor_work * f)
{
    int dimr = f->dimr;
    qsnum_lc_info *lc_inf = f->lc_inf;
    qsnum_lr_info *lr_inf = f->lr_inf;
    QSnum_type *lrcoef = 0; 
    int *lrindx = 0;
    QSnum_type *lccoef = f->lccoef;
    int *lcindx = f->lcindx;
    int nzcnt;
    int beg;
    int rnzcnt;
    int rbeg;
    int i;
    int j;
    int k;
    int c;
    int rval;

    nzcnt = 0;
    for (i = 0; i < dimr; i++) {
        nzcnt += lc_inf[i].nzcnt;
        lr_inf[i].nzcnt = 0;
        lr_inf[i].delay = 0;
        lc_inf[lc_inf[i].c].crank = i;
    }

    QSnum_FreeArray (f->lrcoef, f->lr_space);
    if (nzcnt) {
       lrcoef = QSnum_AllocArray (nzcnt);
       f->lr_space = nzcnt; /* added by dan */
       f->lrcoef = lrcoef;
    }

    CG_IFFREE (f->lrindx, int);
    CG_SAFE_MALLOC (lrindx, nzcnt + 1, int);
    f->lrindx = lrindx;

    for (i = 0; i < dimr; i++) {
        nzcnt = lc_inf[i].nzcnt;
        beg = lc_inf[i].cbeg;
        lc_inf[i].delay = 0;
        for (j = 0; j < nzcnt; j++) {
            lr_inf[lc_inf[lcindx[beg + j]].crank].nzcnt++;
        }
    }

    nzcnt = 0;
    for (i = 0; i < dimr; i++) {
        lr_inf[i].rbeg = nzcnt;
        nzcnt += lr_inf[i].nzcnt;
        lr_inf[i].nzcnt = 0;
        lr_inf[i].r = lc_inf[i].c;
        lr_inf[lr_inf[i].r].rrank = i;
    }

    for (i = 0; i < dimr; i++) {
        nzcnt = lc_inf[i].nzcnt;
        beg = lc_inf[i].cbeg;
        c = lc_inf[i].c;
        for (j = 0; j < nzcnt; j++) {
            k = lc_inf[lcindx[beg + j]].crank;
            rbeg = lr_inf[k].rbeg;
            rnzcnt = lr_inf[k].nzcnt;
            lrindx[rbeg + rnzcnt] = c;
            QSnum_Copy (lrcoef[rbeg + rnzcnt], lccoef[beg + j]);
            lr_inf[k].nzcnt++;
        }
    }

    rval = 0;

CLEANUP:
    return rval;
}

static int qsnum_handle_singularity (qsnum_factor_work * f)
{
    int rval = 0;
    int nsing;
    int *singr = 0;
    int *singc = 0;
    int i;

    if (f->p_nsing == 0 || f->p_singr == 0 || f->p_singc == 0) {
        fprintf (stderr, "singular basis, but no place for singularity data\n");
        return E_SING_NO_DATA;
    }

    nsing = f->nstages - f->stage;
    CG_SAFE_MALLOC (singr, nsing, int);
    CG_SAFE_MALLOC (singc, nsing, int);
    for (i = f->stage; i < f->nstages; i++) {
        singr[i - f->stage] = f->rperm[i];
        singc[i - f->stage] = f->cperm[i];
    }
    *f->p_nsing = nsing;
    *f->p_singr = singr;
    *f->p_singc = singc;
    singr = 0;
    singc = 0;

CLEANUP:
    CG_IFFREE (singr, int);
    CG_IFFREE (singc, int);
    return rval;
}

static int qsnum_dense_build_matrix (qsnum_factor_work * f)
{
    QSnum_type *dmat = 0;
    int stage = f->stage;
    int drows = f->nstages - stage;
    int dcols = f->dimr - stage;
    int dsize = drows * dcols;
    int *crank = f->crank;
    QSnum_type *urcoef = f->urcoef;
    int *urindx = f->urindx;
    int nzcnt;
    int beg;
    int i;
    int r;
    int j;
    int rval = 0;

    dmat = QSnum_AllocArray (dsize);
    /*  WHY no check for NULL?  -- Bico */

    for (i = 0; i < dsize; i++)
        QSnum_SetZero (dmat[i]);

    for (i = 0; i < drows; i++) {
        r = f->rperm[i + stage];
        nzcnt = f->ur_inf[r].nzcnt;
        beg = f->ur_inf[r].rbeg;
        for (j = 0; j < nzcnt; j++) {
            QSnum_Copy (dmat[i * dcols - stage + crank[urindx[beg + j]]],
                                      urcoef[beg + j]);
        }
    }

    f->drows = drows;
    f->dcols = dcols;
    f->dense_base = f->stage;
    f->dmat = dmat;
    f->dsize = dsize;
    dmat = 0;

/* CLEANUP:  */
    return rval;
}

static int qsnum_dense_find_pivot (qsnum_factor_work * f, int *p_r, int *p_c)
{
    int dcols = f->dcols;
    int drows = f->drows;
    QSnum_type *dmat = f->dmat;
    int dense_base = f->dense_base;
    int s = f->stage - dense_base;
    qsnum_ur_info *ur_inf = f->ur_inf;
    int *rperm = f->rperm;
    QSnum_type maxval;
    int max_r;
    int max_c;
    int i;
    QSnum_Init (maxval);
    QSnum_SetZero (maxval);
    max_r = -1;
    for (i = s; i < drows; i++) {
        if (QSnum_Less (maxval, ur_inf[rperm[dense_base + i]].max)) {
            QSnum_Copy (maxval, ur_inf[rperm[dense_base + i]].max);
            max_r = i;
        }
    }
    if (max_r == -1) {
        return E_NO_PIVOT;
    }

    QSnum_SetZero (maxval);
    max_c = -1;
    for (i = s; i < drows; i++) {
        QSnum_CopyMaxAbsAndDo (maxval, dmat[max_r * dcols + i], max_c = i);
    }
    if (max_c == -1) {
        return E_NO_PIVOT;
    }
    *p_r = max_r;
    *p_c = max_c;

    QSnum_Clear (maxval);
    return 0;
}

static void qsnum_dense_swap (qsnum_factor_work * f, int r, int c)
{
    int dcols = f->dcols;
    int drows = f->drows;
    QSnum_type *dmat = f->dmat;
    int dense_base = f->dense_base;
    int s = f->stage - dense_base;
    int i;
    QSnum_type v;
    QSnum_Init (v);

    if (r != s) {
        CG_SWAP (f->rperm[dense_base + s], f->rperm[dense_base + r], i);
        f->rrank[f->rperm[dense_base + s]] = dense_base + s;
        f->rrank[f->rperm[dense_base + r]] = dense_base + r;
        for (i = 0; i < dcols; i++) {
            QSnum_SWAP (dmat[s * dcols + i], dmat[r * dcols + i], v);
        }
    }
    if (c != s) {
        CG_SWAP (f->cperm[dense_base + s], f->cperm[dense_base + c], i);
        f->crank[f->cperm[dense_base + s]] = dense_base + s;
        f->crank[f->cperm[dense_base + c]] = dense_base + c;
        for (i = 0; i < drows; i++) {
            QSnum_SWAP (dmat[i * dcols + s], dmat[i * dcols + c], v);
        }
    }
    QSnum_Clear (v);
}

static void qsnum_dense_elim (qsnum_factor_work * f, int r, int c)
{
    int dcols = f->dcols;
    int drows = f->drows;
    QSnum_type *dmat = f->dmat;
    int dense_base = f->dense_base;
    int s = f->stage - dense_base;
    qsnum_ur_info *ur_inf = f->ur_inf;
    int *rperm = f->rperm;
    int i;
    int j;
    QSnum_type pivval;
    QSnum_type max;
    QSnum_type v;
    QSnum_type w;
    QSnum_Init (pivval);
    QSnum_Init (max);
    QSnum_Init (v);
    QSnum_Init (w);

    qsnum_dense_swap (f, r, c);
    f->stage++;
    QSnum_CopyFrac (pivval, qsnum_oneLpNum, dmat[s * dcols + s]);
    for (i = s + 1; i < drows; i++) {
        QSnum_Copy (v, dmat[i * dcols + s]);
        if (QSnum_NeqZero (v)) {
            QSnum_CopyMult (v, pivval);
            if (QSnum_NeqZeroTol (v, f->fzero_tol)) {
                QSnum_Copy (dmat[i * dcols + s], v);
                QSnum_SetZero (max);
                for (j = s + 1; j < drows; j++) {
                    QSnum_Copy (w, dmat[i * dcols + j]);
                    QSnum_CopySubProd (w, v, dmat[s * dcols + j]);
                    QSnum_Copy (dmat[i * dcols + j], w);
                    QSnum_CopyMaxAbs (max, w);
                }
                for (j = drows; j < dcols; j++) {
                    QSnum_Copy (w, dmat[i * dcols + j]);
                    QSnum_CopySubProd (w, v, dmat[s * dcols + j]);
                    QSnum_Copy (dmat[i * dcols + j], w);
                }
                QSnum_Copy (ur_inf[rperm[dense_base + i]].max, max);
            } else {
                QSnum_SetZero (dmat[i * dcols + s]);
            }
        }
    }
    QSnum_Clear (pivval);
    QSnum_Clear (max);
    QSnum_Clear (v);
    QSnum_Clear (w);
}

static int qsnum_dense_replace_row (qsnum_factor_work * f, int i)
{
    int dcols = f->dcols;
    int dense_base = f->dense_base;
    QSnum_type *dmat = f->dmat + i * dcols;
    QSnum_type *urcoef;
    qsnum_ur_info *ur_inf = f->ur_inf;
    int *cperm = f->cperm;
    int r = f->rperm[dense_base + i];
    int *urindx;
    int nzcnt;
    int beg;
    int j;
    int rval = 0;

    nzcnt = 0;
    for (j = i; j < dcols; j++) {
        if (QSnum_NeqZeroTol (dmat[j], f->fzero_tol)) {
            nzcnt++;
        }
    }
    if (nzcnt > ur_inf[r].nzcnt) {
        if (ur_inf[r].rbeg + ur_inf[r].nzcnt == f->ur_freebeg) {
            f->ur_freebeg = ur_inf[r].rbeg;
        }
        ur_inf[r].nzcnt = 0;
        if (f->ur_freebeg + nzcnt > f->ur_space) {
            rval = qsnum_make_ur_space (f, nzcnt);
            CG_CLEANUP_IF (rval);
        }
        ur_inf[r].rbeg = f->ur_freebeg;
        f->ur_freebeg += nzcnt;
    }
    beg = ur_inf[r].rbeg;
    urcoef = f->urcoef;
    urindx = f->urindx;
    for (j = i; j < dcols; j++) {
        if (QSnum_NeqZeroTol (dmat[j], f->fzero_tol)) {
            QSnum_Copy (urcoef[beg], dmat[j]);
            urindx[beg] = cperm[dense_base + j];
            beg++;
        }
    }
    ur_inf[r].nzcnt = beg - ur_inf[r].rbeg;
CLEANUP:
    return rval;
}

static int qsnum_dense_create_col (qsnum_factor_work * f, int i)
{
    int dcols = f->dcols;
    int drows = f->drows;
    int dense_base = f->dense_base;
    QSnum_type *dmat = f->dmat;
    QSnum_type *lccoef;
    qsnum_lc_info *lc_inf = f->lc_inf;
    int *rperm = f->rperm;
    int *lcindx;
    int nzcnt;
    int beg;
    int j;
    int rval = 0;

    nzcnt = 0;
    for (j = i + 1; j < drows; j++) {
        if (QSnum_NeqZeroTol (dmat[j * dcols + i], f->fzero_tol)) {
            nzcnt++;
        }
    }

    if (f->lc_freebeg + nzcnt >= f->lc_space) {
        rval = qsnum_make_lc_space (f, nzcnt);
        CG_CLEANUP_IF (rval);
    }
    beg = f->lc_freebeg;
    lc_inf[dense_base + i].cbeg = beg;
    lc_inf[dense_base + i].c = rperm[dense_base + i];
    lcindx = f->lcindx;
    lccoef = f->lccoef;

    for (j = i + 1; j < drows; j++) {
        if (QSnum_NeqZeroTol (dmat[j * dcols + i], f->fzero_tol)) {
            QSnum_Copy (lccoef[beg], dmat[j * dcols + i]);
            lcindx[beg] = rperm[dense_base + j];
            beg++;
        }
    }
    lc_inf[dense_base + i].nzcnt = beg - lc_inf[dense_base + i].cbeg;
    f->lc_freebeg = beg;
CLEANUP:
    return rval;
}

static int qsnum_dense_replace (qsnum_factor_work * f)
{
    int drows = f->drows;
    int rval = 0;
    int i;

    for (i = 0; i < drows; i++) {
        rval = qsnum_dense_replace_row (f, i);
        CG_CLEANUP_IF (rval);
        rval = qsnum_dense_create_col (f, i);
        CG_CLEANUP_IF (rval);
    }

    QSnum_FreeArray (f->dmat, f->dsize);
    f->dsize = 0;
    f->drows = 0;
    f->dcols = 0;
CLEANUP:
    return rval;
}

static int qsnum_dense_factor (qsnum_factor_work * f)
{
    int r;
    int c;
    int rval = 0;

    rval = qsnum_dense_build_matrix (f);
    CG_CLEANUP_IF (rval);

    while (f->stage < f->nstages) {
        r = f->stage - f->dense_base;
        rval = qsnum_dense_find_pivot (f, &r, &c);
        if (rval == E_NO_PIVOT) {
            rval = qsnum_handle_singularity (f);
            CG_CLEANUP_IF (rval);
            return E_SINGULAR_INTERNAL;
        } else {
            CG_CLEANUP_IF (rval);
        }
        qsnum_dense_elim (f, r, c);
    }

    rval = qsnum_dense_replace (f);
    CG_CLEANUP_IF (rval);

CLEANUP:
    return rval;
}

static int qsnum_factor_try (qsnum_factor_work * f, int *basis, int *cbeg,
        int *clen, int *cindx, QSnum_type * ccoef)
{
    int rval = 0;
    int r;
    int c;

    rval = qsnum_init_matrix (f, basis, cbeg, clen, cindx, ccoef);
    CG_CLEANUP_IF (rval);

    f->stage = 0;
    f->nstages = f->dimr;

    while (f->stage < f->nstages) {
        rval = qsnum_find_pivot (f, &r, &c);
        if (rval == E_NO_PIVOT) {
            printf ("NO_PIVOT\n"); fflush (stdout);
            rval = qsnum_handle_singularity (f);
            CG_CLEANUP_IF (rval);
            return 0;
        } else {
            CG_CLEANUP_IF (rval);
        }
        if (f->ur_inf[r].pivcnt > f->dense_fract * (f->nstages - f->stage) &&
                f->uc_inf[c].nzcnt > f->dense_fract * (f->nstages - f->stage) &&
                f->nstages - f->stage > f->dense_min) {
            rval = qsnum_dense_factor (f);
            if (rval == E_SINGULAR_INTERNAL) return 0;
            if (rval) return rval;
            break;
        }
        rval = qsnum_elim (f, r, c);
        CG_CLEANUP_IF (rval);
    }

    rval = qsnum_build_iteration_u_data (f);
    CG_CLEANUP_IF (rval);

    rval = qsnum_build_iteration_l_data (f);
    CG_CLEANUP_IF (rval);

CLEANUP:
    return rval;
}

int QSnum_factor (qsnum_factor_work * f, int *basis, int *cbeg, int *clen,
        int *cindx, QSnum_type * ccoef, int *p_nsing, int **p_singr,
        int **p_singc)
{
    int rval;
    f->p_nsing = p_nsing;
    f->p_singr = p_singr;
    f->p_singc = p_singc;
    *p_nsing = 0;

AGAIN:
    rval = qsnum_factor_try (f, basis, cbeg, clen, cindx, ccoef);
    if (rval == E_FACTOR_BLOWUP) {
        if (QSnum_LessDbl (f->partial_cur, 0.1)) {
            QSnum_CopyMultUi (f->partial_cur, 10);
        } else if (QSnum_LessDbl (f->partial_cur, 0.25)) {
            QSnum_CopyDbl (f->partial_cur, 0.25);
        } else if (QSnum_LessDbl (f->partial_cur, 0.5)) {
            QSnum_CopyDbl (f->partial_cur, 0.5);
        } else if (QSnum_Less (f->partial_cur, qsnum_oneLpNum)) {
            QSnum_SetOne (f->partial_cur);
        } else {
            return rval;
        }
        goto AGAIN;
    }
    return rval;
}

static void qsnum_factor_ftranl (qsnum_factor_work * f, QSnum_type * a)
{
    int *lcindx = f->lcindx;
    qsnum_lc_info *lc_inf = f->lc_inf;
    QSnum_type *lccoef = f->lccoef;
    int dimr = f->dimr;
    int beg;
    int nzcnt;
    int i;
    int j;
    QSnum_type v;
    QSnum_Init (v);

    for (i = 0; i < dimr; i++) {
        QSnum_Copy (v, a[lc_inf[i].c]);
        if (QSnum_NeqZero (v)) {
            nzcnt = lc_inf[i].nzcnt;
            beg = lc_inf[i].cbeg;
            for (j = 0; j < nzcnt; j++) {
                QSnum_CopySubProd (a[lcindx[beg + j]], v, lccoef[beg + j]);
            }
        }
    }
    QSnum_Clear (v);
}

static void qsnum_ftranl3_delay2 (qsnum_factor_work * f, int c)
{
    qsnum_lc_info *lc_inf = f->lc_inf;
    int nzcnt;
    int *indx;
    int i;
    int last;

    do {
        c = lc_inf[c].crank;
        nzcnt = lc_inf[c].nzcnt;
        indx = f->lcindx + lc_inf[c].cbeg;
        last = -1;
        for (i = 0; i < nzcnt; i++) {
            c = indx[i];
            if (lc_inf[c].delay++ == 0) {
                if (last >= 0) {
                    qsnum_ftranl3_delay2 (f, last);
                }
                last = c;
            }
        }
        c = last;
    } while (c >= 0);
}

static void qsnum_ftranl3_process2 (qsnum_factor_work * f, int c, qsnum_svector * x)
{
    qsnum_lc_info *lc_inf = f->lc_inf;
    QSnum_type *work = f->work_coef;
    int nzcnt;
    int *indx;
    QSnum_type *coef;
    int i;
    int last;
    QSnum_type v;
    QSnum_Init (v);

    do {
        QSnum_Copy (v, work[c]);
        QSnum_SetZero (work[c]);
        if (QSnum_NeqZero (v)) {
            x->indx[x->nzcnt] = c;
            QSnum_Copy (x->coef[x->nzcnt], v);
            x->nzcnt++;
        }
        c = lc_inf[c].crank;
        nzcnt = lc_inf[c].nzcnt;
        indx = f->lcindx + lc_inf[c].cbeg;
        coef = f->lccoef + lc_inf[c].cbeg;
        last = -1;
        for (i = 0; i < nzcnt; i++) {
            c = indx[i];
            QSnum_CopySubProd (work[c], v, coef[i]);
            if (--lc_inf[c].delay == 0) {
                if (last >= 0) {
                    qsnum_ftranl3_process2 (f, last, x);
                }
                last = c;
            }
        }
        c = last;
    } while (c >= 0);
    QSnum_Clear (v);
}

static void qsnum_factor_ftranl3 (qsnum_factor_work * f, qsnum_svector * a,
        qsnum_svector * x)
{
    QSnum_type *work = f->work_coef;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    QSnum_type *acoef = a->coef;
    qsnum_lc_info *lc_inf = f->lc_inf;
    int i;

    for (i = 0; i < anzcnt; i++) {
        if (lc_inf[aindx[i]].delay++ == 0) {
            qsnum_ftranl3_delay2 (f, aindx[i]);
        }
        QSnum_Copy (work[aindx[i]], acoef[i]);
    }
    x->nzcnt = 0;
    for (i = 0; i < anzcnt; i++) {
        if (--lc_inf[aindx[i]].delay == 0) {
            qsnum_ftranl3_process2 (f, aindx[i], x);
        }
    }
}

static void qsnum_factor_ftrane (qsnum_factor_work * f, QSnum_type * a)
{
    int *erindx = f->erindx;
    QSnum_type *ercoef = f->ercoef;
    qsnum_er_info *er_inf = f->er_inf;
    int etacnt = f->etacnt;
    int beg;
    int nzcnt;
    int i;
    int j;
    QSnum_type v;
    QSnum_Init (v);

    for (i = 0; i < etacnt; i++) {
        QSnum_Copy (v, a[er_inf[i].r]);
        nzcnt = er_inf[i].nzcnt;
        beg = er_inf[i].rbeg;
        for (j = 0; j < nzcnt; j++) {
            QSnum_CopySubProd (v, ercoef[beg + j], a[erindx[beg + j]]);
        }
        QSnum_Copy (a[er_inf[i].r], v);
    }
    QSnum_Clear (v);
}

static void qsnum_factor_ftrane2 (qsnum_factor_work * f, qsnum_svector * a)
{
    int *erindx = f->erindx;
    QSnum_type *ercoef = f->ercoef;
    qsnum_er_info *er_inf = f->er_inf;
    int etacnt = f->etacnt;
    int beg;
    int nzcnt;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    QSnum_type *acoef = a->coef;
    QSnum_type *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    int i;
    int j;
    int r;
    QSnum_type v;
    QSnum_Init (v);

    for (i = 0; i < anzcnt; i++) {
        QSnum_Copy (work_coef[aindx[i]], acoef[i]);
        work_indx[aindx[i]] = i + 1;
    }
    for (i = 0; i < etacnt; i++) {
        r = er_inf[i].r;
        QSnum_Copy (v, work_coef[r]);
        nzcnt = er_inf[i].nzcnt;
        beg = er_inf[i].rbeg;
        for (j = 0; j < nzcnt; j++) {
            QSnum_CopySubProd (v, ercoef[beg + j], work_coef[erindx[beg + j]]);
        }
        if (QSnum_NeqZero (v)) {
            QSnum_Copy (work_coef[r], v);
            if (work_indx[r] == 0) {
                QSnum_Copy (acoef[anzcnt], v);
                aindx[anzcnt] = r;
                work_indx[r] = anzcnt + 1;
                anzcnt++;
            } else {
                QSnum_Copy (acoef[work_indx[r] - 1], v);
            }
        } else {
            QSnum_SetZero (work_coef[r]);
            if (work_indx[r]) {
                QSnum_SetZero (acoef[work_indx[r] - 1]);
            }
        }
    }
    i = 0;
    while (i < anzcnt) {
        QSnum_SetZero (work_coef[aindx[i]]);
        work_indx[aindx[i]] = 0;
        if (QSnum_NeqZeroTol (acoef[i], f->fzero_tol)) {
            /*if (acoef[i] > fzero_tol || acoef[i] < -fzero_tol) */
            i++;
        } else {
            --anzcnt;
            QSnum_Copy (acoef[i], acoef[anzcnt]);
            aindx[i] = aindx[anzcnt];
        }
    }
    a->nzcnt = anzcnt;

    QSnum_Clear (v);
}

static void qsnum_factor_ftranu (qsnum_factor_work * f, QSnum_type * a,
        qsnum_svector * x)
{
    int *ucindx = f->ucindx;
    QSnum_type *uccoef = f->uccoef;
    qsnum_uc_info *uc_inf = f->uc_inf;
    int *cperm = f->cperm;
    int *rperm = f->rperm;
    int dimr = f->dimr;
    int xnzcnt = 0;
    int *xindx = x->indx;
    QSnum_type *xcoef = x->coef;
    int nzcnt;
    int beg;
    int i;
    int j;
    QSnum_type v;
    QSnum_Init (v);

    for (i = dimr - 1; i >= 0; i--) {
        QSnum_Copy (v, a[rperm[i]]);
        if (QSnum_NeqZero (v)) {
            /*((v = a[rperm[i]]) != 0.0) */
            j = cperm[i];
            beg = uc_inf[j].cbeg;
            QSnum_CopyDiv (v, uccoef[beg]);
            if (QSnum_NeqZeroTol (v, f->szero_tol)) {
                /*if (v > szero_tol || v < -szero_tol) */
                xindx[xnzcnt] = j;
                QSnum_Copy (xcoef[xnzcnt], v);
                xnzcnt++;
            }
            nzcnt = uc_inf[j].nzcnt;
            for (j = 1; j < nzcnt; j++)
            {
                QSnum_CopySubProd (a[ucindx[beg + j]], v, uccoef[beg + j]);
            }
            QSnum_SetZero (a[rperm[i]]);
        }
    }
    x->nzcnt = xnzcnt;
    QSnum_Clear (v);
}

static void qsnum_ftranu3_delay2 (qsnum_factor_work * f, int c)
{
    qsnum_uc_info *uc_inf = f->uc_inf;
    int nzcnt;
    int *indx;
    int i;
    int last;

    do {
        c = f->cperm[f->rrank[c]];
        nzcnt = uc_inf[c].nzcnt;
        indx = f->ucindx + uc_inf[c].cbeg;
        last = -1;
        for (i = 1; i < nzcnt; i++) {
            c = indx[i];
            if (uc_inf[c].delay++ == 0) {
                if (last >= 0) {
                    qsnum_ftranu3_delay2 (f, last);
                }
                last = c;
            }
        }
        c = last;
    } while (c >= 0);
}

static void qsnum_ftranu3_process2 (qsnum_factor_work * f, int c, qsnum_svector * x)
{
    qsnum_uc_info *uc_inf = f->uc_inf;
    QSnum_type *work = f->work_coef;
    int nzcnt;
    int *indx;
    QSnum_type *coef;
    int i;
    int last;
    QSnum_type v;
    QSnum_Init (v);

    do {
        QSnum_Copy (v, work[c]);
        QSnum_SetZero (work[c]);
        c = f->cperm[f->rrank[c]];
        nzcnt = uc_inf[c].nzcnt;
        indx = f->ucindx + uc_inf[c].cbeg;
        coef = f->uccoef + uc_inf[c].cbeg;
        QSnum_CopyDiv (v, coef[0]);
        if (QSnum_NeqZeroTol (v, f->szero_tol)) {
            /*if (v > szero_tol || v < -szero_tol) */
            x->indx[x->nzcnt] = c;
            QSnum_Copy (x->coef[x->nzcnt], v);
            x->nzcnt++;
        }
        last = -1;
        for (i = 1; i < nzcnt; i++) {
            c = indx[i];
            QSnum_CopySubProd (work[c], v, coef[i]);
            if (--uc_inf[c].delay == 0) {
                if (last >= 0) {
                    qsnum_ftranu3_process2 (f, last, x);
                }
                last = c;
            }
        }
        c = last;
    } while (c >= 0);
    QSnum_Clear (v);
}

static void qsnum_factor_ftranu3 (qsnum_factor_work * f, qsnum_svector * a,
        qsnum_svector * x)
{
    QSnum_type *work = f->work_coef;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    QSnum_type *acoef = a->coef;
    qsnum_uc_info *uc_inf = f->uc_inf;
    int i;

    for (i = 0; i < anzcnt; i++) {
        if (uc_inf[aindx[i]].delay++ == 0) {
            qsnum_ftranu3_delay2 (f, aindx[i]);
        }
        QSnum_Copy (work[aindx[i]], acoef[i]);
    }
    x->nzcnt = 0;
    for (i = 0; i < anzcnt; i++) {
        if (--uc_inf[aindx[i]].delay == 0) {
            qsnum_ftranu3_process2 (f, aindx[i], x);
        }
    }
}

void QSnum_factor_ftran (qsnum_factor_work * f, qsnum_svector * a,
        qsnum_svector * x)
{
    int i;
    int nzcnt;
    int sparse;
    int *aindx;
    QSnum_type *acoef;
    QSnum_type *work_coef = f->work_coef;

    if (a->nzcnt >= SPARSE_FACTOR * f->dimr) {
        nzcnt = a->nzcnt;
        aindx = a->indx;
        acoef = a->coef;
        for (i = 0; i < nzcnt; i++) {
            QSnum_Copy (work_coef[aindx[i]], acoef[i]);
        }
        sparse = 0;
    } else {
        sparse = 1;
    }

    if (sparse) {
        qsnum_factor_ftranl3 (f, a, &f->xtmp);
        if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dimr) {
            nzcnt = f->xtmp.nzcnt;
            aindx = f->xtmp.indx;
            acoef = f->xtmp.coef;

            for (i = 0; i < nzcnt; i++) {
                QSnum_Copy (work_coef[aindx[i]], acoef[i]);
            }
            sparse = 0;
        }
    } else {
        qsnum_factor_ftranl (f, work_coef);
    }

    if (sparse) {
        qsnum_factor_ftrane2 (f, &f->xtmp);
        if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dimr) {
            nzcnt = f->xtmp.nzcnt;
            aindx = f->xtmp.indx;
            acoef = f->xtmp.coef;

            for (i = 0; i < nzcnt; i++) {
                QSnum_Copy (work_coef[aindx[i]], acoef[i]);
            }
            sparse = 0;
        }
    } else {
        qsnum_factor_ftrane (f, work_coef);
    }

    if (sparse) {
        qsnum_factor_ftranu3 (f, &f->xtmp, x);
    } else {
        qsnum_factor_ftranu (f, work_coef, x);
    }

    return;
}

static void qsnum_factor_btranl2 (qsnum_factor_work * f, QSnum_type * x)
{
    int *lrindx = f->lrindx;
    QSnum_type *lrcoef = f->lrcoef;
    qsnum_lr_info *lr_inf = f->lr_inf;
    int dimr = f->dimr;
    int nzcnt;
    int beg;
    int i;
    int j;
    QSnum_type v;
    QSnum_Init (v);

    for (i = dimr - 1; i >= 0; i--) {
        QSnum_Copy (v, x[lr_inf[i].r]);
        if (QSnum_NeqZero (v)) {
            nzcnt = lr_inf[i].nzcnt;
            beg = lr_inf[i].rbeg;
            for (j = 0; j < nzcnt; j++) {
                QSnum_CopySubProd (x[lrindx[beg + j]], v, lrcoef[beg + j]);
            }
        }
    }
    QSnum_Clear (v);
}


static void qsnum_btranl3_delay2 (qsnum_factor_work * f, int r)
{
    qsnum_lr_info *lr_inf = f->lr_inf;
    int nzcnt;
    int *indx;
    int i;
    int last;

    do {
        r = lr_inf[r].rrank;
        nzcnt = lr_inf[r].nzcnt;
        indx = f->lrindx + lr_inf[r].rbeg;
        last = -1;
        for (i = 0; i < nzcnt; i++) {
            r = indx[i];
            if (lr_inf[r].delay++ == 0) {
                if (last >= 0) {
                    qsnum_btranl3_delay2 (f, last);
                }
                last = r;
            }
        }
        r = last;
    } while (r >= 0);
}

static void qsnum_btranl3_process2 (qsnum_factor_work * f, int r, qsnum_svector * x)
{
    qsnum_lr_info *lr_inf = f->lr_inf;
    QSnum_type *work = f->work_coef;
    int nzcnt;
    int *indx;
    QSnum_type *coef;
    int i;
    int last;
    QSnum_type v;
    QSnum_Init (v);

    do {
        QSnum_Copy (v, work[r]);
        QSnum_SetZero (work[r]);
        if (QSnum_NeqZeroTol (v, f->szero_tol)) {
            /*if (v > szero_tol || v < -szero_tol) */
            x->indx[x->nzcnt] = r;
            QSnum_Copy (x->coef[x->nzcnt], v);
            x->nzcnt++;
        }
        r = lr_inf[r].rrank;
        nzcnt = lr_inf[r].nzcnt;
        indx = f->lrindx + lr_inf[r].rbeg;
        coef = f->lrcoef + lr_inf[r].rbeg;
        last = -1;
        for (i = 0; i < nzcnt; i++) {
            r = indx[i];
            QSnum_CopySubProd (work[r], v, coef[i]);
            if (--lr_inf[r].delay == 0) {
                if (last >= 0) {
                    qsnum_btranl3_process2 (f, last, x);
                }
                last = r;
            }
        }
        r = last;
    } while (r >= 0);
    QSnum_Clear (v);
}

static void qsnum_factor_btranl3 (qsnum_factor_work * f, qsnum_svector * a,
        qsnum_svector * x)
{
    QSnum_type *work = f->work_coef;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    QSnum_type *acoef = a->coef;
    qsnum_lr_info *lr_inf = f->lr_inf;
    int i;

    for (i = 0; i < anzcnt; i++) {
        if (lr_inf[aindx[i]].delay++ == 0) {
            qsnum_btranl3_delay2 (f, aindx[i]);
        }
        QSnum_Copy (work[aindx[i]], acoef[i]);
    }
    x->nzcnt = 0;
    for (i = 0; i < anzcnt; i++) {
        if (--lr_inf[aindx[i]].delay == 0) {
            qsnum_btranl3_process2 (f, aindx[i], x);
        }
    }
}

static void qsnum_factor_btrane (qsnum_factor_work * f, QSnum_type * x)
{
    int *erindx = f->erindx;
    QSnum_type *ercoef = f->ercoef;
    qsnum_er_info *er_inf = f->er_inf;
    int etacnt = f->etacnt;
    int beg;
    int nzcnt;
    int i;
    int j;
    QSnum_type v;
    QSnum_Init (v);

    for (i = etacnt - 1; i >= 0; i--) {
        QSnum_Copy (v, x[er_inf[i].r]);
        if (QSnum_NeqZero (v)) {
            nzcnt = er_inf[i].nzcnt;
            beg = er_inf[i].rbeg;
            for (j = 0; j < nzcnt; j++) {
                QSnum_CopySubProd (x[erindx[beg + j]], v, ercoef[beg + j]);
            }
        }
    }
    QSnum_Clear (v);
}

static void qsnum_factor_btrane2 (qsnum_factor_work * f, qsnum_svector * x)
{
    int *erindx = f->erindx;
    QSnum_type *ercoef = f->ercoef;
    qsnum_er_info *er_inf = f->er_inf;
    int etacnt = f->etacnt;
    int beg;
    int nzcnt;
    int xnzcnt = x->nzcnt;
    int *xindx = x->indx;
    QSnum_type *xcoef = x->coef;
    QSnum_type *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    int i;
    int j;
    QSnum_type v;
    QSnum_Init (v);

    for (i = 0; i < xnzcnt; i++) {
        QSnum_Copy (work_coef[xindx[i]], xcoef[i]);
        work_indx[xindx[i]] = i + 1;
    }
    for (i = etacnt - 1; i >= 0; i--) {
        QSnum_Copy (v, work_coef[er_inf[i].r]);
        if (QSnum_NeqZero (v)) {
            nzcnt = er_inf[i].nzcnt;
            beg = er_inf[i].rbeg;
            for (j = 0; j < nzcnt; j++) {
                if (work_indx[erindx[beg + j]] == 0) {
                    work_indx[erindx[beg + j]] = xnzcnt;
                    xindx[xnzcnt++] = erindx[beg + j];
                }
                QSnum_CopySubProd (work_coef[erindx[beg + j]], v,
                                             ercoef[beg + j]);
            }
        }
    }

    j = 0;
    while (j < xnzcnt) {
        QSnum_Copy (xcoef[j], work_coef[xindx[j]]);
        QSnum_SetZero (work_coef[xindx[j]]);
        work_indx[xindx[j]] = 0;
        if (QSnum_Equal (xcoef[j], qsnum_zeroLpNum)) {
            --xnzcnt;
            xindx[j] = xindx[xnzcnt];
        } else {
            j++;
        }
    }
    x->nzcnt = xnzcnt;
    QSnum_Clear (v);
}

static void qsnum_factor_btranu (qsnum_factor_work * f, QSnum_type * a,
        qsnum_svector * x)
{
    int *urindx = f->urindx;
    QSnum_type *urcoef = f->urcoef;
    qsnum_ur_info *ur_inf = f->ur_inf;
    int *rperm = f->rperm;
    int *cperm = f->cperm;
    int dimc = f->dimc;
    int xnzcnt = 0;
    int *xindx = x->indx;
    QSnum_type *xcoef = x->coef;
    int nzcnt;
    int beg;
    int i;
    int j;
    QSnum_type v;
    QSnum_Init (v);

    for (i = 0; i < dimc; i++) {
        QSnum_Copy (v, a[cperm[i]]);
        if (QSnum_NeqZero (v)) {
            j = rperm[i];
            beg = ur_inf[j].rbeg;
            QSnum_CopyDiv (v, urcoef[beg]);
            if (QSnum_NeqZeroTol (v, f->szero_tol)) {
                /* if (v > szero_tol || v < -szero_tol) */
                xindx[xnzcnt] = j;
                QSnum_Copy (xcoef[xnzcnt], v);
                xnzcnt++;
            }
            nzcnt = ur_inf[j].nzcnt;
            for (j = 1; j < nzcnt; j++) {
                QSnum_CopySubProd (a[urindx[beg + j]], v, urcoef[beg + j]);
            }
            QSnum_SetZero (a[cperm[i]]);
        }
    }
    x->nzcnt = xnzcnt;
    QSnum_Clear (v);
}

static void qsnum_btranu3_delay2 (qsnum_factor_work * f, int r)
{
    qsnum_ur_info *ur_inf = f->ur_inf;
    int nzcnt;
    int *indx;
    int i;
    int last;

    do {
        r = f->rperm[f->crank[r]];
        nzcnt = ur_inf[r].nzcnt;
        indx = f->urindx + ur_inf[r].rbeg;
        last = -1;
        for (i = 1; i < nzcnt; i++) {
            r = indx[i];
            if (ur_inf[r].delay++ == 0) {
                if (last >= 0) {
                    qsnum_btranu3_delay2 (f, last);
                }
                last = r;
            }
        }
        r = last;
    } while (r >= 0);
}

static void qsnum_btranu3_process2 (qsnum_factor_work * f, int r, qsnum_svector * x)
{
    qsnum_ur_info *ur_inf = f->ur_inf;
    QSnum_type *work = f->work_coef;
    int nzcnt;
    int *indx;
    QSnum_type *coef;
    int i;
    int last;
    QSnum_type v;
    QSnum_Init (v);

    do {
        QSnum_Copy (v, work[r]);
        QSnum_SetZero (work[r]);
        r = f->rperm[f->crank[r]];
        nzcnt = ur_inf[r].nzcnt;
        indx = f->urindx + ur_inf[r].rbeg;
        coef = f->urcoef + ur_inf[r].rbeg;
        QSnum_CopyDiv (v, coef[0]);
        if (QSnum_NeqZero (v)) {
            x->indx[x->nzcnt] = r;
            QSnum_Copy (x->coef[x->nzcnt], v);
            x->nzcnt++;
        }
        last = -1;
        for (i = 1; i < nzcnt; i++) {
            r = indx[i];
            QSnum_CopySubProd (work[r], v, coef[i]);
            if (--ur_inf[r].delay == 0) {
                if (last >= 0) {
                    qsnum_btranu3_process2 (f, last, x);
                }
                last = r;
            }
        }
        r = last;
    } while (r >= 0);
    QSnum_Clear (v);
}

static void qsnum_factor_btranu3 (qsnum_factor_work * f, qsnum_svector * a,
        qsnum_svector * x)
{
    QSnum_type *work = f->work_coef;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    QSnum_type *acoef = a->coef;
    qsnum_ur_info *ur_inf = f->ur_inf;
    int i;

    for (i = 0; i < anzcnt; i++) {
        if (ur_inf[aindx[i]].delay++ == 0) {
            qsnum_btranu3_delay2 (f, aindx[i]);
        }
        QSnum_Copy (work[aindx[i]], acoef[i]);
    }
    x->nzcnt = 0;
    for (i = 0; i < anzcnt; i++) {
        if (--ur_inf[aindx[i]].delay == 0) {
            qsnum_btranu3_process2 (f, aindx[i], x);
        }
    }
}

/* QSnum_facctor_btran solves x^tB=a^t (or, B^t x = a) for x */
void QSnum_factor_btran (qsnum_factor_work * f, qsnum_svector * a,
        qsnum_svector * x)
{
    int i;
    int nzcnt;
    int sparse;
    int *aindx = a->indx;
    QSnum_type *acoef = a->coef;
    QSnum_type *work_coef = f->work_coef;
    int dimr = f->dimr;

    if (a->nzcnt >= SPARSE_FACTOR * f->dimr) {
        aindx = a->indx;
        acoef = a->coef;
        work_coef = f->work_coef;
        nzcnt = a->nzcnt;
        for (i = 0; i < nzcnt; i++) {
            QSnum_Copy (work_coef[aindx[i]], acoef[i]);
        }
        sparse = 0;
    } else {
        sparse = 1;
    }

    if (sparse) {
        qsnum_factor_btranu3 (f, a, &f->xtmp);
    } else {
        qsnum_factor_btranu (f, work_coef, &f->xtmp);
    }

    if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dimr) {
        aindx = f->xtmp.indx;
        acoef = f->xtmp.coef;
        work_coef = f->work_coef;
        nzcnt = f->xtmp.nzcnt;
        for (i = 0; i < nzcnt; i++) {
            QSnum_Copy (work_coef[aindx[i]], acoef[i]);
        }
        sparse = 0;
    } else {
        sparse = 1;
    }

    if (sparse) {
        qsnum_factor_btrane2 (f, &f->xtmp);
        if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dimr) {
            aindx = f->xtmp.indx;
            acoef = f->xtmp.coef;
            work_coef = f->work_coef;
            nzcnt = f->xtmp.nzcnt;
            for (i = 0; i < nzcnt; i++) {
                QSnum_Copy (work_coef[aindx[i]], acoef[i]);
            }
            sparse = 0;
        }
    } else {
        qsnum_factor_btrane (f, work_coef);
    }

    if (sparse) {
        qsnum_factor_btranl3 (f, &f->xtmp, x);
    } else {
        qsnum_factor_btranl2 (f, work_coef);
        dimr = f->dimr;
        nzcnt = 0;
        aindx = x->indx;
        acoef = x->coef;
        for (i = 0; i < dimr; i++) {
            if (QSnum_NeqZero (work_coef[i])) {
                if (QSnum_NeqZeroTol (work_coef[i], f->szero_tol)) {
                /*if (work_coef[i] > szero_tol || work_coef[i] < -szero_tol) */
                    aindx[nzcnt] = i;
                    QSnum_Copy (acoef[nzcnt], work_coef[i]);
                    nzcnt++;
                }
                QSnum_SetZero (work_coef[i]);
            }
        }
        x->nzcnt = nzcnt;
    }
    return;
}

void QSnum_svector_init (qsnum_svector * s)
{
    s->nzcnt = 0;
    s->indx = 0;
    s->coef = 0;
}

void QSnum_svector_free (qsnum_svector * s)
{
    CG_IFFREE (s->indx, int);
    QSnum_FreeArray (s->coef, s->nzcnt);
    s->nzcnt = 0;
}
int init_sxvector (qsnum_svector *v, int n)
{
    int rval = 0;
    int i;

    if (v) {
        v->nzcnt = n;
        v->indx = (int *) NULL;
        v->coef = (QSnum_type *) NULL;
        if (n > 0) {
            v->indx = (int *) malloc (n * sizeof (int));
            v->coef = (QSnum_type *) malloc (n * sizeof (QSnum_type));
            if (!v->indx || !v->coef) {
                fprintf (stderr, "out of memory for v\n");
                rval = 1;  goto CLEANUP;
            }
            for (i = 0; i < n; i++) QSnum_Init (v->coef[i]);
        }
    }

CLEANUP:
    return rval;
}
int clear_sxvector (qsnum_svector *v)
{
  if(v->indx) free(v->indx);
  QSnum_FreeArray(v->coef,v->nzcnt);
  return 0;
}

int QSnum_svector_alloc (qsnum_svector * s, int nzcnt)
{
    int rval = 0;

    s->nzcnt = nzcnt;
    if (nzcnt == 0) {
        s->indx = 0;
        s->coef = 0;
    } else {
        CG_SAFE_MALLOC (s->indx, nzcnt, int);
        s->coef = QSnum_AllocArray (nzcnt);
        /*  Why no check for NULL? --- Bico */
    }
    return 0;
CLEANUP:
    CG_IFFREE (s->indx, int);
    return rval;
}

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
   )
{
   int rval = 0;
   int nsing = 0;
   int *singr = 0;
   int *singc = 0;

   assert(*f == NULL);
   assert(m >= n);
   
   if(m<n)
      return 1;
   
   QSnum_factor_init ();

   *f = (qsnum_factor_work *) malloc (sizeof (qsnum_factor_work));
   CGcheck_NULL (*f, "out of memory for f");
   QSnum_factor_init_factor_work (*f);
   rval = QSnum_factor_create_factor_work (*f,n,m);
   CGcheck_rval (rval, "QSnum_factor_create_factor_work failed");
   
   rval = QSnum_factor (*f, basisx, matbeg, matcnt, matind, matval, &nsing,&singr, &singc);
   CGcheck_rval (rval, "QSsnum_factor failed");
   if (nsing > 0) {
      printf ("Matrix is NONSINGULAR!\n");
      rval = 1;  goto CLEANUP;
   }

CLEANUP:
   return rval;
}    

/** solves the system using factor work f */ 
int RECTLUsolveSystem(
   qsnum_factor_work* f,                 /**< pointer to factor work*/
   int              n,                   /**< number of rows in matrix */
   int              m,                   /**< number of columns in matrix*/ 
   mpq_t*           rhs,                 /**< right hand side of system to solve */ 
   mpq_t*           sol                  /**< solution to system, initialized to zero */
)
{
   int rval = 0;
   qsnum_svector srhs, ssol;
   int i;
   int rhsnz;
   int temp;
   
   rhsnz =0;
   for(i=0; i < n; i++)
   {
      if(mpq_sgn(rhs[i]))
         rhsnz++;
   }

   /* set up sparse vectors for solve interface */
   //rval = init_sxvector (&srhs, rhsnz); 
   QSnum_svector_init(&srhs);
   rval = QSnum_svector_alloc(&srhs, rhsnz);   
   CGcheck_rval (rval, "init_sxvector failed");

   temp = 0;
   for(i=0; i < n; i++)
   {
      if(mpq_sgn(rhs[i]))
      {
         mpq_set(srhs.coef[temp],rhs[i]);
         srhs.indx[temp] = i;
         temp++;
      }
   }

   /* sparse solution will have at most n nonzeros because it */
   /* will be a basic solution */
   //   rval = init_sxvector (&ssol, n);
   QSnum_svector_init(&ssol);
   rval = QSnum_svector_alloc(&ssol, n); 
   CGcheck_rval (rval, "init_sxvector failed");

   QSnum_factor_ftran (f, &srhs, &ssol);


   for(i = 0; i < m; i++)
      mpq_set_si(sol[i],0,1);

   /* assign the solution value */
   for(i = 0 ; i < ssol.nzcnt; i++)
      mpq_set(sol[ssol.indx[i]],ssol.coef[i]);

CLEANUP:
   clear_sxvector(&srhs);
   ssol.nzcnt = n; /* this has to be done to make sure it frees all the mpqs, nzcnt is reset to be lower in ftranu3 - dan */
   clear_sxvector(&ssol);
   return rval;
}

/** frees factor work f */ 
void RECTLUfreeFactorization(
   qsnum_factor_work* f                 /**< factor work */
)
{
   QSnum_factor_free_factor_work (f);
   QSnum_factor_clear ();
   if (f) free (f);
}
