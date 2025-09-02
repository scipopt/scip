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

/**@file   lpexact.c
 * @brief  LP management methods and data structures for exact mirror of LP
 * @author Leon Eifler
 *
 *  In LP management, we have to distinguish between the current LP and the SCIP_LP
 *  stored in the LP solver. All LP methods affect the current LP only.
 *  Before solving the current LP with the LP solver or setting an LP state,
 *  the LP solvers data has to be updated to the current LP with a call to
 *  lpExactFlush().
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "lpi/lpi.h"
#include "lpiexact/lpiexact.h"
#include "scip/clock.h"
#include "scip/cutpool.h"
#include "scip/event.h"
#include "scip/intervalarith.h"
#include "scip/lp.h"
#include "scip/lpexact.h"
#include "scip/misc.h"
#include "scip/prob.h"
#include "scip/pub_cons.h"
#include "scip/pub_lp.h"
#include "scip/pub_lpexact.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_var.h"
#include "scip/pub_tree.h"
#include "scip/rational.h"
#include "scip/scip_lp.h"
#include "scip/scip_lpexact.h"
#include "scip/scip_message.h"
#include "scip/set.h"
#include "scip/sepastoreexact.h"
#include "scip/sol.h"
#include "scip/solve.h"
#include "scip/stat.h"
#include "scip/struct_event.h"
#include "scip/struct_lpexact.h"
#include "scip/struct_prob.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_var.h"
#include "scip/struct_cutpool.h"
#include "scip/var.h"
#include <string.h>
#include <inttypes.h>

/** comparison method for sorting rows by non-decreasing index */
SCIP_DECL_SORTPTRCOMP(SCIProwExactComp)
{
   assert(elem1 != NULL);
   assert(elem2 != NULL);

   assert(((SCIP_ROWEXACT*)elem1)->fprow != NULL);
   assert(((SCIP_ROWEXACT*)elem2)->fprow != NULL);

   if( ((SCIP_ROWEXACT*)elem1)->index < ((SCIP_ROWEXACT*)elem2)->index )
      return -1;
   else if( ((SCIP_ROWEXACT*)elem1)->index > ((SCIP_ROWEXACT*)elem2)->index )
      return +1;
   else
   {
      assert(SCIProwExactGetIndex(((SCIP_ROWEXACT*)elem1))
         == SCIProwExactGetIndex(((SCIP_ROWEXACT*)elem2)));
      return 0;
   }
}

#if 0 /* enable this to check links between columns and rows in LP data structure (for debugging, very slow!) */

#ifdef NDEBUG
#define ASSERT(x) do { if( !(x) ) abort(); } while( FALSE )
#else
#define ASSERT(x) assert(x)
#endif

static SCIP_Bool msgdisp_checklinks = FALSE;


static
void checkLinks(
   SCIP_LPEXACT*         lp                  /**< current LP data */
   )
{
   SCIP_COLEXACT* col;
   SCIP_ROWEXACT* row;
   int i;
   int j;

   ASSERT(lp != NULL);

   if( !msgdisp_checklinks )
   {
      printf("LP LINK CHECKING ACTIVATED! THIS IS VERY SLOW!\n");
      msgdisp_checklinks = TRUE;
   }

   for( i = 0; i < lp->ncols; ++i )
   {
      col = lp->cols[i];
      ASSERT(col != NULL);
      ASSERT(!lp->flushed || col->lppos >= 0);
      ASSERT(!lp->flushed || col->lppos >= 0);
      ASSERT(col->nlprows <= col->len);
      ASSERT(col->lppos == -1 || col->lppos >= lp->lpifirstchgcol || col->nunlinked == 0);

      for( j = 0; j < col->len; ++j )
      {
         row = col->rows[j];
         ASSERT(row != NULL);
         ASSERT(!lp->flushed || col->lppos == -1 || col->linkpos[j] >= 0);
         ASSERT(col->linkpos[j] == -1 || row->cols[col->linkpos[j]] == col);
         ASSERT(col->linkpos[j] == -1 || SCIPrationalIsEQ(row->vals[col->linkpos[j]], col->vals[j]));
         ASSERT((j < col->nlprows) == (col->linkpos[j] >= 0 && row->lppos >= 0));
      }
   }

   for( i = 0; i < lp->nrows; ++i )
   {
      row = lp->rows[i];
      ASSERT(row != NULL);
      ASSERT(!lp->flushed || row->lppos >= 0);
      ASSERT(!lp->flushed || row->lppos >= 0);
      ASSERT(row->nlpcols <= row->len);
      ASSERT(row->lppos == -1 || row->lppos >= lp->lpifirstchgrow || row->nunlinked == 0);

      for( j = 0; j < row->len; ++j )
      {
         col = row->cols[j];
         ASSERT(col != NULL);
         ASSERT(!lp->flushed || row->lppos == -1 || row->linkpos[j] >= 0);
         ASSERT(row->linkpos[j] == -1 || col->rows[row->linkpos[j]] == row);
         ASSERT(row->linkpos[j] == -1 || SCIPrationalIsEQ(col->vals[row->linkpos[j]], row->vals[j]));
         ASSERT((j < row->nlpcols) == (row->linkpos[j] >= 0 && col->lppos >= 0));
      }
   }
}

#undef ASSERT

#else
#define checkLinks(lp) /**/
#endif

#ifndef NDEBUG
/** checks if the exact column and its fpcol are consistent */
static
SCIP_Bool colExactInSync(
   SCIP_COLEXACT*        colexact,           /**< exact column */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_COL* fpcol;

   assert(colexact != NULL);

   fpcol = colexact->fpcol;
   assert(fpcol != NULL);

   assert(colexact->var == fpcol->var);
   assert(colexact->lpipos == fpcol->lpipos);
   assert(colexact->index == fpcol->index);
   assert(colexact->len >= fpcol->nlprows);

   assert(SCIPrationalIsApproxEQReal(set, colexact->obj, fpcol->obj, SCIP_R_ROUND_NEAREST));
   assert(SCIPrationalIsApproxEQReal(set, colexact->flushedobj, fpcol->flushedobj, SCIP_R_ROUND_NEAREST));
   assert(SCIPrationalIsApproxEQReal(set, colexact->lb, fpcol->lb, SCIP_R_ROUND_DOWNWARDS) || (SCIPrationalIsNegInfinity(colexact->lb) && SCIPsetIsInfinity(set, -fpcol->lb)));
   assert(SCIPrationalIsApproxEQReal(set, colexact->ub, fpcol->ub, SCIP_R_ROUND_UPWARDS) || (SCIPrationalIsInfinity(colexact->ub) && SCIPsetIsInfinity(set, fpcol->ub)));

   return TRUE;
}

/** checks if the exact row and its fprow are consistent */
static
SCIP_Bool rowExactInSync(
   SCIP_ROWEXACT*        rowexact,           /**< exact row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     msg                 /**< message handler for debug output */ /*lint !e204*/
   ) /*lint --e{715}*/
{
   SCIP_ROW* fprow;
   SCIP_Bool synced;

   assert(rowexact != NULL);

   fprow = rowexact->fprow;

   assert(fprow != NULL);
   assert(rowexact->len >= fprow->len);
   assert(rowexact->lppos == rowexact->fprow->lppos);

   synced = SCIPrationalIsGEReal(rowexact->lhs, fprow->lhs) || (SCIPrationalIsNegInfinity(rowexact->lhs) && SCIPsetIsInfinity(set, -fprow->lhs));
   synced = synced && (SCIPrationalIsLEReal(rowexact->rhs, fprow->rhs) || (SCIPrationalIsInfinity(rowexact->rhs) && SCIPsetIsInfinity(set, fprow->rhs)));
   synced = synced && (SCIPrationalIsApproxEQReal(set, rowexact->constant, fprow->constant, SCIP_R_ROUND_NEAREST) );

   if( !synced )
   {
      SCIPdebug(SCIProwPrint(rowexact->fprow, msg, NULL));
      SCIPdebug(SCIProwExactPrint(rowexact, msg, NULL));
      SCIPABORT();
   }

   return TRUE;
} /*lint !e715*/
#endif

/** checks if the exact lp and lp are consistent (same number of rows/cols, and all cols/rows in sync) */
static
SCIP_Bool lpExactInSync(
   SCIP_LPEXACT*         lpexact,            /**< exact lp */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     msg                 /**< message handler for debug output */
   )
{
#ifndef NDEBUG
   int i;
   SCIP_LP* fplp;

   assert(lpexact != NULL);

   fplp = lpexact->fplp;
   assert(fplp != NULL);

   assert(lpexact->nrows == fplp->nrows);
   for( i = 0; i < lpexact->nrows; i++)
   {
      assert(rowExactInSync(lpexact->rows[i], set, msg));
   }

   assert(lpexact->ncols == fplp->ncols);
   for( i = 0; i < lpexact->ncols; i++)
   {
      assert(colExactInSync(lpexact->cols[i], set));
   }
#endif

   return TRUE;
}

/** ensures that rows array can store at least num entries */
static
SCIP_RETCODE ensureRowexsSize(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lpexact->nrows <= lpexact->rowssize);

   if( num > lpexact->rowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpexact->rows, newsize) );
      lpexact->rowssize = newsize;
   }
   assert(num <= lpexact->rowssize);

   return SCIP_OKAY;
}

/** sorts column entries of linked rows currently in the LP such that lower row indices precede higher ones */
static
void colExactSortLP(
   SCIP_COLEXACT*        col                 /**< column to be sorted */
   )
{
   int i;

   assert(col != NULL);

   /* check, if column is already sorted in the LP part */
   if( col->lprowssorted )
      return;

   /* sort coefficients */
   SCIPsortPtrPtrInt((void**)col->rows, (void**)col->vals, col->linkpos, SCIProwExactComp, col->nlprows );

   /* update links */
   for( i = 0; i < col->nlprows; ++i )
   {
      if( col->linkpos[i] >= 0 )
      {
         assert(col->rows[i]->cols[col->linkpos[i]] == col);
         assert(col->rows[i]->linkpos[col->linkpos[i]] >= 0);
         col->rows[i]->linkpos[col->linkpos[i]] = i;
      }
   }

   col->lprowssorted = TRUE;
}

/** sorts column entries of unlinked rows or rows currently not in the LP such that lower row indices precede higher
 *  ones
 */
static
void colExactSortNonLP(
   SCIP_COLEXACT*        col                 /**< column to be sorted */
   )
{
   int i;

   assert(col != NULL);

   /* check, if column is already sorted in the non-LP part */
   if( col->nonlprowssorted )
      return;

   /* sort coefficients */
   SCIPsortPtrPtrInt((void**)(&(col->rows[col->nlprows])),
                     (void**)(&(col->vals[col->nlprows])),
                     &(col->linkpos[col->nlprows]), SCIProwExactComp,
                     col->len - col->nlprows);

   /* update links */
   for( i = col->nlprows; i < col->len; ++i )
   {
      if( col->linkpos[i] >= 0 )
      {
         assert(col->rows[i]->cols[col->linkpos[i]] == col);
         assert(col->rows[i]->linkpos[col->linkpos[i]] >= 0);
         col->rows[i]->linkpos[col->linkpos[i]] = i;
      }
   }

   col->nonlprowssorted = TRUE;
}

/** sorts row entries of linked columns currently in the LP such that lower column indices precede higher ones */
static
void rowExactSortLP(
   SCIP_ROWEXACT*        row                 /**< row to be sorted */
   )
{
   int i;

   assert(row != NULL);

   /* check, if row is already sorted in the LP part, or if the sorting should be delayed */
   if( row->lpcolssorted || row->delaysort )
      return;

   /* sort coefficients */
   SCIPsortIntIntPtrPtrInterval(row->cols_index, row->linkpos, (void**)row->cols,
      (void**)row->vals, row->valsinterval, row->nlpcols);

   /* update links */
   for( i = 0; i < row->nlpcols; ++i )
   {
      if( row->linkpos[i] >= 0 )
      {
         assert(row->cols[i]->rows[row->linkpos[i]] == row);
         assert(row->cols[i]->linkpos[row->linkpos[i]] >= 0);
         row->cols[i]->linkpos[row->linkpos[i]] = i;
      }
   }

   row->lpcolssorted = TRUE;
}

/** sorts row entries of unlinked columns or columns currently not in the LP such that lower column indices precede
 *  higher ones
 */
static
void rowExactSortNonLP(
   SCIP_ROWEXACT*        row                 /**< row to be sorted */
   )
{
   int i;

   assert(row != NULL);

   /* check, if row is already sorted in the non-LP part, or if the sorting should be delayed */
   if( row->nonlpcolssorted || row->delaysort )
      return;

   /* sort coefficients */
   SCIPsortIntIntPtrPtrInterval(&(row->cols_index[row->nlpcols]), &(row->linkpos[row->nlpcols]),
      (void**)(&(row->cols[row->nlpcols])), (void**)&(row->vals[row->nlpcols]), &(row->valsinterval[row->nlpcols]),
      row->len - row->nlpcols);

   /* update links */
   for( i = row->nlpcols; i < row->len; ++i )
   {
      if( row->linkpos[i] >= 0 )
      {
         assert(row->cols[i]->rows[row->linkpos[i]] == row);
         assert(row->cols[i]->linkpos[row->linkpos[i]] >= 0);
         row->cols[i]->linkpos[row->linkpos[i]] = i;
      }
   }

   row->nonlpcolssorted = TRUE;
}

/** ensures, that row array of column can store at least num entries */
static
SCIP_RETCODE colExactEnsureSize(
   SCIP_COLEXACT*        col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(col != NULL);
   assert(col->len <= col->size);

   if( num > col->size )
   {
      int newsize;
      int i;

      /* realloc fpcol */
      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->rows, col->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->vals, col->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->linkpos, col->size, newsize) );

      /* realloc colexact */
      for( i = col->size; i < newsize; ++i )
      {
         SCIP_CALL( SCIPrationalCreateBlock(blkmem, &col->vals[i]) );
      }

      col->size = newsize;
   }
   assert(num <= col->size);

   return SCIP_OKAY;
}

/** ensures, that cols array can store at least num entries */
static
SCIP_RETCODE ensureColexsSize(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->ncols <= lp->colssize);

   if( num > lp->colssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->cols, newsize) );
      lp->colssize = newsize;
   }
   assert(num <= lp->colssize);

   return SCIP_OKAY;
}

/** ensures, that chgcols array can store at least num entries */
static
SCIP_RETCODE ensureChgcolsSize(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nchgcols <= lp->chgcolssize);

   if( num > lp->chgcolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->chgcols, newsize) );
      lp->chgcolssize = newsize;
   }
   assert(num <= lp->chgcolssize);

   return SCIP_OKAY;
}

/** ensures, that lpicols array can store at least num entries */
static
SCIP_RETCODE ensureLpiExactcolsSize(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nlpicols <= lp->lpicolssize);

   if( num > lp->lpicolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->lpicols, newsize) );
      lp->lpicolssize = newsize;
   }
   assert(num <= lp->lpicolssize);

   return SCIP_OKAY;
}

/** ensures, that lpirows array can store at least num entries */
static
SCIP_RETCODE ensureLpirowexactsSize(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nlpirows <= lp->lpirowssize);

   if( num > lp->lpirowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->lpirows, newsize) );
      lp->lpirowssize = newsize;
   }
   assert(num <= lp->lpirowssize);

   return SCIP_OKAY;
}

/** searches coefficient in part of the column, returns position in col vector or -1 if not found */
static
int colExactSearchCoefPart(
   SCIP_COLEXACT*        col,                /**< column to be searched in */
   const SCIP_ROWEXACT*  row,                /**< coefficient to be searched for */
   int                   minpos,             /**< first position of search range */
   int                   maxpos              /**< last position of search range */
   )
{
   int pos;
   int idx;
   int searchidx;

   assert(col != NULL);
   assert(row != NULL);

   /* binary search */
   searchidx = row->index;
   while(minpos <= maxpos)
   {
      pos = (minpos + maxpos)/2;
      assert(0 <= pos && pos < col->len);
      assert(col->rows[pos] != NULL);
      assert((pos < col->nlprows) == (col->rows[pos]->lppos >= 0 && col->linkpos[pos] >= 0));
      idx = col->rows[pos]->index;
      if( searchidx == idx )
         return pos;
      else if( searchidx < idx )
         maxpos = pos-1;
      else
         minpos = pos+1;
   }

   return -1;
}

/** searches coefficient in column, returns position in col vector or -1 if not found */
static
int colExactSearchCoef(
   SCIP_COLEXACT*        col,                /**< column to be searched in */
   const SCIP_ROWEXACT*  row                 /**< coefficient to be searched for */
   )
{
   int pos;

   assert(col != NULL);
   assert(row != NULL);

   pos = -1;

   /* search in the linked LP rows */
   if( row->lppos >= 0 )
   {
      /* column has to be sorted, such that binary search works */
      colExactSortLP(col);
      assert(col->lprowssorted);

      pos = colExactSearchCoefPart(col, row, 0, col->nlprows-1);
      if( pos >= 0 )
         return pos;
   }

   /* search in the non-LP/unlinked rows */
   if( row->lppos == -1 || col->nunlinked > 0 )
   {
      /* column has to be sorted, such that binary search works */
      colExactSortNonLP(col);
      assert(col->nonlprowssorted);

      pos = colExactSearchCoefPart(col, row, col->nlprows, col->len-1);
   }

   return pos;
}

/** searches coefficient in part of the row, returns position in col vector or -1 if not found */
static
int rowExactSearchCoefPart(
   SCIP_ROWEXACT*        row,                /**< row to be searched in */
   const SCIP_COLEXACT*  col,                /**< coefficient to be searched for */
   int                   minpos,             /**< first position of search range */
   int                   maxpos              /**< last position of search range */
   )
{
   int pos;
   int idx;
   int searchidx;

   assert(col != NULL);
   assert(row != NULL);

   /* binary search */
   searchidx = col->index;
   while(minpos <= maxpos)
   {
      pos = (minpos + maxpos)/2;
      assert(0 <= pos && pos < row->len);
      assert(row->cols[pos] != NULL);
      assert((pos < row->nlpcols) == (row->cols[pos]->lppos >= 0 && row->linkpos[pos] >= 0));
      assert(row->cols_index[pos] == row->cols[pos]->index);
      idx = row->cols_index[pos];
      if( searchidx == idx )
         return pos;
      else if( searchidx < idx )
         maxpos = pos-1;
      else
         minpos = pos+1;
   }

   return -1;
}

/** searches coefficient in row, returns position in row vector or -1 if not found;
 *  if the sorting of the row is delayed, returns -1
 */
static
int rowExactSearchCoef(
   SCIP_ROWEXACT*        row,                /**< row to be searched in */
   const SCIP_COLEXACT*  col                 /**< coefficient to be searched for */
   )
{
   int pos;

   assert(col != NULL);
   assert(row != NULL);

   if( row->delaysort )
      return -1;

   pos = -1;

   /* search in the linked LP columns */
   if( col->lppos >= 0 )
   {
      /* row has to be sorted, such that binary search works */
      rowExactSortLP(row);
      assert(row->lpcolssorted);

      pos = rowExactSearchCoefPart(row, col, 0, row->nlpcols-1);
   }

   /* search in the non-LP/unlinked columns */
   if( pos == -1 && (col->lppos == -1 || row->nunlinked > 0) )
   {
      /* row has to be sorted, such that binary search works */
      rowExactSortNonLP(row);
      assert(row->nonlpcolssorted);

      pos = rowExactSearchCoefPart(row, col, row->nlpcols, row->len-1);
   }

#ifndef NDEBUG
   /* validate result */
   assert(-1 <= pos && pos < row->len);
   if( pos >= 0 )
       assert(row->cols[pos] == col);
   else
   {
      int i;
      for( i = 0; i < row->len; ++i )
         assert(row->cols[i] != col);
   }
#endif

   return pos;
}

/*
 * Changing announcements
 */

/** announces, that the given coefficient in the constraint matrix changed */
static
void coefChangedExact(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_COLEXACT*        col,                /**< LP col */
   SCIP_LPEXACT*         lp                  /**< current LP data */
   )
{
   assert(row != NULL);
   assert(col != NULL);
   assert(lp != NULL);

   if( row->lpipos >= 0 && col->lpipos >= 0 )
   {
      assert(row->lpipos < lp->nlpirows);
      assert(col->lpipos < lp->nlpicols);

      /* we have to remember the change only in the row or in the column,
       * because the readdition of one vector would change the other automatically.
       */
      if( row->lpipos >= lp->lpifirstchgrow )
         row->coefchanged = TRUE;
      else if( col->lpipos >= lp->lpifirstchgcol )
         col->coefchanged = TRUE;
      else if( lp->lpifirstchgrow - row->lpipos <= lp->lpifirstchgcol - col->lpipos )
      {
         row->coefchanged = TRUE;
         lp->lpifirstchgrow = row->lpipos;
      }
      else
      {
         col->coefchanged = TRUE;
         lp->lpifirstchgcol = col->lpipos;
      }

      /* mark the current LP unflushed */
      lp->flushed = FALSE;
   }

   SCIPrationalSetInfinity(row->pseudoactivity);
}

/*
 * local column changing methods
 */

/** moves a coefficient in a column to a different place, and updates all corresponding data structures */
static
void colExactMoveCoef(
   SCIP_COLEXACT*        col,                /**< LP column */
   int                   oldpos,             /**< old position of coefficient */
   int                   newpos              /**< new position of coefficient */
   )
{
   assert(col != NULL);
   assert(0 <= oldpos && oldpos < col->len);
   assert(0 <= newpos && newpos < col->len);
   assert(col->rows[oldpos] != NULL);

   if( oldpos == newpos )
      return;

   SCIPrationalSetRational(col->vals[newpos], col->vals[oldpos]);
   col->rows[newpos] = col->rows[oldpos];
   SCIPrationalSetRational(col->vals[newpos], col->vals[oldpos]);
   col->linkpos[newpos] = col->linkpos[oldpos];

   /* update link position in row */
   if( col->linkpos[newpos] >= 0 )
   {
      assert(col->rows[newpos]->cols[col->linkpos[newpos]] == col);
      assert(col->rows[newpos]->linkpos[col->linkpos[newpos]] == oldpos);

      col->rows[newpos]->linkpos[col->linkpos[newpos]] = newpos;
   }

   /* update sorted flags */
   if( col->rows[newpos]->lppos >= 0 && col->linkpos[newpos] >= 0 )
      col->lprowssorted = FALSE;
   else
      col->nonlprowssorted = FALSE;
}

/** swaps two coefficients in a column, and updates all corresponding data structures */
static
SCIP_RETCODE colExactSwapCoefs(
   SCIP_COLEXACT*        col,                /**< LP column */
   BMS_BUFMEM*           buffer,             /**< buffer for temp real */
   int                   pos1,               /**< position of first coefficient */
   int                   pos2                /**< position of second coefficient */
   )
{
   SCIP_ROWEXACT* tmprow;
   SCIP_RATIONAL* tmpval;
   int tmplinkpos;

   assert(col != NULL);
   assert(0 <= pos1 && pos1 < col->len);
   assert(0 <= pos2 && pos2 < col->len);
   assert(col->rows[pos1] != NULL);

   if( pos1 == pos2 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPrationalCreateBuffer(buffer, &tmpval) );

   /* swap coefficients */
   tmprow = col->rows[pos2];
   SCIPrationalSetRational(tmpval, col->vals[pos2]);
   tmplinkpos = col->linkpos[pos2];

   col->rows[pos2] = col->rows[pos1];
   SCIPrationalSetRational(col->vals[pos2], col->vals[pos1]);
   col->linkpos[pos2] = col->linkpos[pos1];

   col->rows[pos1] = tmprow;
   SCIPrationalSetRational(col->vals[pos1], tmpval);
   col->linkpos[pos1] = tmplinkpos;

   SCIPrationalFreeBuffer(buffer, &tmpval);

   /* update link position in rows */
   if( col->linkpos[pos1] >= 0 )
   {
      assert(col->rows[pos1]->cols[col->linkpos[pos1]] == col);
      assert(col->rows[pos1]->linkpos[col->linkpos[pos1]] == pos2);

      col->rows[pos1]->linkpos[col->linkpos[pos1]] = pos1;
   }

   if( col->linkpos[pos2] >= 0 )
   {
      assert(col->rows[pos2]->cols[col->linkpos[pos2]] == col);
      assert(col->rows[pos2]->linkpos[col->linkpos[pos2]] == pos1);

      col->rows[pos2]->linkpos[col->linkpos[pos2]] = pos2;
   }

   /* update sorted flags */
   if( col->rows[pos1]->lppos >= 0 && col->linkpos[pos1] >= 0 )
      col->lprowssorted = FALSE;
   else
      col->nonlprowssorted = FALSE;
   if( col->rows[pos2]->lppos >= 0 && col->linkpos[pos2] >= 0 )
      col->lprowssorted = FALSE;
   else
      col->nonlprowssorted = FALSE;

   return SCIP_OKAY;
}

/** moves a coefficient in a row to a different place, and updates all corresponding data structures */
static
void rowExactMoveCoef(
   SCIP_ROWEXACT*        row,                /**< LP row */
   int                   oldpos,             /**< old position of coefficient */
   int                   newpos              /**< new position of coefficient */
   )
{
   assert(row != NULL);
   assert(0 <= oldpos && oldpos < row->len);
   assert(0 <= newpos && newpos < row->len);
   assert(row->cols[oldpos] != NULL);

   if( oldpos == newpos )
      return;

   row->cols[newpos] = row->cols[oldpos];
   row->cols_index[newpos] = row->cols_index[oldpos];
   SCIPrationalSetRational(row->vals[newpos], row->vals[oldpos]);
   row->valsinterval[newpos] = row->valsinterval[oldpos];
   row->linkpos[newpos] = row->linkpos[oldpos];

   /* update link position in column */
   if( row->linkpos[newpos] >= 0 )
   {
      assert(row->cols[newpos]->rows[row->linkpos[newpos]] == row);
      assert(row->cols[newpos]->linkpos[row->linkpos[newpos]] == oldpos);

      row->cols[newpos]->linkpos[row->linkpos[newpos]] = newpos;
   }

   /* update sorted flags */
   if( row->cols[newpos]->lppos >= 0 && row->linkpos[newpos] >= 0 )
      row->lpcolssorted = FALSE;
   else
      row->nonlpcolssorted = FALSE;
}

/** swaps two coefficients in a row, and updates all corresponding data structures */
static
SCIP_RETCODE rowExactSwapCoefs(
   SCIP_ROWEXACT*        row,                /**< LP row */
   BMS_BUFMEM*           buffer,             /**< buffer for temp real */
   int                   pos1,               /**< position of first coefficient */
   int                   pos2                /**< position of second coefficient */
   )
{
   SCIP_COLEXACT* tmpcol;
   SCIP_RATIONAL* tmpval;
   SCIP_INTERVAL tmp;
   int tmpindex;
   int tmplinkpos;

   assert(row != NULL);
   assert(0 <= pos1 && pos1 < row->len);
   assert(0 <= pos2 && pos2 < row->len);
   assert(row->cols[pos1] != NULL);
   assert(row->cols[pos1]->index == row->cols_index[pos1]);

   if( pos1 == pos2 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPrationalCreateBuffer(buffer, &tmpval) );

   /* swap coefficients */
   tmpcol = row->cols[pos2];
   tmpindex = row->cols_index[pos2];
   SCIPrationalSetRational(tmpval, row->vals[pos2]);
   tmp = row->valsinterval[pos2];
   tmplinkpos = row->linkpos[pos2];

   row->cols[pos2] = row->cols[pos1];
   row->cols_index[pos2] = row->cols_index[pos1];
   SCIPrationalSetRational(row->vals[pos2], row->vals[pos1]);
   row->valsinterval[pos2] = row->valsinterval[pos1];
   row->linkpos[pos2] = row->linkpos[pos1];

   row->cols[pos1] = tmpcol;
   row->cols_index[pos1] = tmpindex;
   SCIPrationalSetRational(row->vals[pos1], tmpval);
   row->valsinterval[pos1] = tmp;
   row->linkpos[pos1] = tmplinkpos;

   SCIPrationalFreeBuffer(buffer, &tmpval);

   /* update link position in columns */
   if( row->linkpos[pos1] >= 0 )
   {
      assert(row->cols[pos1]->rows[row->linkpos[pos1]] == row);
      assert(row->cols[pos1]->linkpos[row->linkpos[pos1]] == pos2);

      row->cols[pos1]->linkpos[row->linkpos[pos1]] = pos1;
   }
   if( row->linkpos[pos2] >= 0 )
   {
      assert(row->cols[pos2]->rows[row->linkpos[pos2]] == row);
      assert(row->cols[pos2]->linkpos[row->linkpos[pos2]] == pos1);

      row->cols[pos2]->linkpos[row->linkpos[pos2]] = pos2;
   }

   /* update sorted flags */
   if( row->cols[pos1]->lppos >= 0 && row->linkpos[pos1] >= 0 )
      row->lpcolssorted = FALSE;
   else
      row->nonlpcolssorted = FALSE;
   if( row->cols[pos2]->lppos >= 0 && row->linkpos[pos2] >= 0 )
      row->lpcolssorted = FALSE;
   else
      row->nonlpcolssorted = FALSE;

   return SCIP_OKAY;
}

/** forward declaration for rowExactAddCoef() */
static
SCIP_RETCODE rowExactAddCoef(
   SCIP_ROWEXACT*        row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_COLEXACT*        col,                /**< LP column */
   SCIP_RATIONAL*        val,                /**< value of coefficient */
   int                   linkpos             /**< position of row in the column's row array, or -1 */
   );


/** insert column in the chgcols list (if not already there) */
static
SCIP_RETCODE insertColChgcols(
   SCIP_COLEXACT*        col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lp                  /**< current LP data */
   )
{
   if( !col->objchanged && !col->lbchanged && !col->ubchanged )
   {
      SCIP_CALL( ensureChgcolsSize(lp, set, lp->nchgcols+1) );
      lp->chgcols[lp->nchgcols] = col;
      lp->nchgcols++;
   }

   /* mark the current LP unflushed */
   lp->flushed = FALSE;

   return SCIP_OKAY;
}

/** adds a previously non existing coefficient to an LP column */
static
SCIP_RETCODE colExactAddCoef(
   SCIP_COLEXACT*        col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_RATIONAL*        val,                /**< value of coefficient */
   int                   linkpos             /**< position of column in the row's col array, or -1 */
   )
{
   int pos;

   assert(blkmem != NULL);
   assert(col != NULL);
   assert(col->nlprows <= col->len);
   assert(col->var != NULL);
   assert(!SCIPrationalIsZero(val));

   SCIP_CALL( colExactEnsureSize(col, blkmem, set, col->len+1) );
   assert(col->rows != NULL);
   assert(col->vals != NULL);
   assert(col->linkpos != NULL);

   pos = col->len;
   col->len++;

   /* if the row is in current LP and is linked to the column, we have to insert it at the end of the linked LP rows
    * part of the column's arrays
    */
   if( row->lppos >= 0 && linkpos >= 0 )
   {
      /* move the first non-LP/not linked row to the end */
      if( col->nlprows < pos )
      {
         colExactMoveCoef(col, col->nlprows, pos);
         pos = col->nlprows;
      }
      col->nlprows++;
   }

   /* insert the row at the correct position and update the links */
   col->rows[pos] = row;

   if( col->vals[pos] != NULL )
      SCIPrationalSetRational(col->vals[pos], val);
   else
      SCIP_CALL( SCIPrationalCopyBlock(blkmem, &col->vals[pos], val) );

   col->linkpos[pos] = linkpos;
   if( linkpos == -1 )
   {
      col->nunlinked++;

      /* if the column is in current LP, we have to link it to the row, because otherwise, the primal information
       * of the row is not complete
       */
      if( col->lppos >= 0 )
      {
         /* this call might swap the current row with the first non-LP/not linked row, s.t. insertion position
          * has to be updated
          */
         SCIP_CALL( rowExactAddCoef(row, blkmem, set, eventqueue, lp, col, val, pos) );
         if( row->lppos >= 0 )
            pos = col->nlprows-1;
         linkpos = col->linkpos[pos];

         assert(0 <= linkpos && linkpos < row->len);
         assert(row->cols[linkpos] == col);
         assert(col->rows[pos] == row);
         assert(col->rows[pos]->cols[col->linkpos[pos]] == col);
         assert(col->rows[pos]->fprow->linkpos[col->linkpos[pos]] == pos);
      }
   }
   else
   {
      assert(row->linkpos[linkpos] == -1);
      assert(row->nunlinked > 0);
      row->linkpos[linkpos] = pos;
      row->nunlinked--;

      /* if the column is in current LP, now both conditions, row->cols[linkpos]->lppos >= 0 and row->linkpos[linkpos] >= 0
       * hold, so we have to move the column to the linked LP-cols part of the row's cols array
       */
      if( col->lppos >= 0 )
      {
         row->nlpcols++;
         SCIP_CALL( rowExactSwapCoefs(row, set->buffer, linkpos, row->nlpcols-1) );

         /* if no swap was necessary, mark nonlpcols to be unsorted */
         if( linkpos == row->nlpcols-1 )
            row->lpcolssorted = FALSE;
      }
   }

   /* update the sorted flags */
   if( row->lppos >= 0 && linkpos >= 0 )
   {
      assert(col->nlprows >= 1);
      assert(col->rows[col->nlprows-1] == row);
      if( col->nlprows > 1 )
      {
         col->lprowssorted = col->lprowssorted
            && (col->rows[col->nlprows-2]->index < row->index);
      }
   }
   else
   {
      assert(col->len - col->nlprows >= 1);
      assert(col->rows[col->len-1] == row);
      if( col->len - col->nlprows > 1 )
      {
         col->nonlprowssorted = col->nonlprowssorted
            && (col->rows[col->len-2]->index < row->index);
      }
   }

   coefChangedExact(row, col, lp);

   SCIPrationalDebugMessage("added coefficient %q * <%s> at position %d (%d/%d) to column <%s> (nunlinked=%d)\n",
      val, row->fprow->name, pos, col->nlprows, col->len,
      SCIPvarGetName(col->var), col->nunlinked);

   return SCIP_OKAY;
}

/** deletes coefficient at given position from column */
static
SCIP_RETCODE colExactDelCoefPos(
   SCIP_COLEXACT*        col,                /**< column to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   int                   pos                 /**< position in column vector to delete */
   )
{
   SCIP_ROWEXACT* row;

   assert(lpexact != NULL);
   assert(col != NULL);
   assert(col->var != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < col->len);
   assert(col->rows[pos] != NULL);
   assert(col->linkpos[pos] == -1 || col->rows[pos]->cols[col->linkpos[pos]] == col);
   assert((pos < col->nlprows) == (col->linkpos[pos] >= 0 && col->rows[pos]->lppos >= 0));

   row = col->rows[pos];
   assert((row->lppos >= 0) == (pos < col->nlprows));

   SCIPrationalDebugMessage("deleting coefficient %q * <%p> at position %d from column <%s>\n",
      col->vals[pos], (void*) row, pos, SCIPvarGetName(col->var));

   if( col->linkpos[pos] == -1 )
      col->nunlinked--;

   /* if row is a linked LP row, move last linked LP coefficient to position of empty slot (deleted coefficient) */
   if( pos < col->nlprows )
   {
      colExactMoveCoef(col, col->nlprows-1, pos);
      col->nlprows--;
      pos = col->nlprows;
   }

   /* move last coefficient to position of empty slot */
   colExactMoveCoef(col, col->len-1, pos);
   col->len--;

   coefChangedExact(row, col, lpexact);

   return SCIP_OKAY;
}

/** changes a coefficient at given position of an LP column */
static
SCIP_RETCODE colExactChgCoefPos(
   SCIP_COLEXACT*        col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   int                   pos,                /**< position in column vector to change */
   SCIP_RATIONAL*        val                 /**< value of coefficient */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(0 <= pos && pos < col->len);
   assert(col->rows[pos] != NULL);
   assert(col->linkpos[pos] == -1 || col->rows[pos]->cols[col->linkpos[pos]] == col);

   SCIPrationalDebugMessage("changing coefficient %q * <%s> at position %d of column <%s> to %g\n",
     col->vals[pos], col->rows[pos]->fprow->name, pos, SCIPvarGetName(col->var), val);

   if( SCIPrationalIsZero(val) )
   {
      /* delete existing coefficient */
      SCIP_CALL( colExactDelCoefPos(col, set, lp, pos) );
   }
   else if( !SCIPrationalIsEQ(col->vals[pos], val) )
   {
      /* change existing coefficient */
      SCIPrationalSetRational(col->vals[pos], val);
      coefChangedExact(col->rows[pos], col, lp);
   }

   return SCIP_OKAY;
}


/** forward declaration for rowEactAddCoef() */
static
SCIP_RETCODE rowExactAddCoef(
   SCIP_ROWEXACT*        row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_COLEXACT*        col,                /**< LP column */
   SCIP_RATIONAL*        val,                /**< value of coefficient */
   int                   linkpos             /**< position of row in the column's row array, or -1 */
   )
{
   int pos;

   assert(row != NULL);
   assert(row->nlpcols <= row->len);
   assert(blkmem != NULL);
   assert(col != NULL);
   assert(col->var != NULL);
   assert(!SCIPrationalIsZero(val));

   if( row->nlocks > 0 )
   {
      SCIPerrorMessage("cannot add a coefficient to the locked unmodifiable row <%s>\n", row->fprow->name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIProwExactEnsureSize(row, blkmem, set, row->len+1) );
   assert(row->cols != NULL);
   assert(row->vals != NULL);

   pos = row->len;
   row->len++;

   /* if the column is in current LP and is linked to the row, we have to insert it at the end of the linked LP columns
    * part of the row's arrays
    */
   if( col->lppos >= 0 && linkpos >= 0 )
   {
      /* move the first non-LP/not linked column to the end */
      if( row->nlpcols < pos )
      {
         rowExactMoveCoef(row, row->nlpcols, pos);
         pos = row->nlpcols;
      }
      row->nlpcols++;
   }

   /* insert the column at the correct position and update the links */
   row->cols[pos] = col;
   row->cols_index[pos] = col->index;
   if( row->vals[pos] == NULL )
      SCIP_CALL( SCIPrationalCopyBlock(blkmem, &row->vals[pos], val) );
   else
      SCIPrationalSetRational(row->vals[pos], val);

   SCIPintervalSetRational(&row->valsinterval[pos], row->vals[pos]);
   row->linkpos[pos] = linkpos;
   row->integral = row->integral && SCIPcolIsIntegral(col->fpcol) && SCIPrationalIsIntegral(val);
   if( linkpos == -1 )
   {
      row->nunlinked++;

      /* if the row is in current LP, we have to link it to the column, because otherwise, the dual information
       * of the column is not complete
       */
      if( row->lppos >= 0 )
      {
         /* this call might swap the current column with the first non-LP/not linked column, s.t. insertion position
          * has to be updated
          */
         SCIP_CALL( colExactAddCoef(col, blkmem, set, eventqueue, lp, row, val, pos) );
         if( col->lppos >= 0 )
            pos = row->nlpcols-1;
         linkpos = row->linkpos[pos];

         assert(0 <= linkpos && linkpos < col->len);
         assert(col->rows[linkpos] == row);
         assert(row->cols[pos] == col);
         assert(row->cols[pos]->rows[row->linkpos[pos]] == row);
         assert(row->cols[pos]->linkpos[row->linkpos[pos]] == pos);
      }
   }
   else
   {
      assert(col->linkpos[linkpos] == -1);
      assert(col->nunlinked > 0);
      col->linkpos[linkpos] = pos;
      col->nunlinked--;

      /* if the row is in current LP, now both conditions, col->rows[linkpos]->lppos >= 0 and col->linkpos[linkpos] >= 0
       * hold, so we have to move the row to the linked LP-rows part of the column's rows array
       */
      if( row->lppos >= 0 )
      {
         col->nlprows++;
         SCIP_CALL( colExactSwapCoefs(col, set->buffer, linkpos, col->nlprows-1) );

         /* if no swap was necessary, mark lprows to be unsorted */
         if( linkpos == col->nlprows-1 )
            col->lprowssorted = FALSE;
      }
   }

   /* update the sorted flags */
   if( col->lppos >= 0 && linkpos >= 0 )
   {
      assert(row->nlpcols >= 1);
      assert(row->cols[row->nlpcols-1] == col);
      if( row->nlpcols > 1 )
      {
         assert(row->cols_index[row->nlpcols-2] == row->cols[row->nlpcols-2]->index);
         row->lpcolssorted = row->lpcolssorted && (row->cols_index[row->nlpcols-2] < col->index);
      }
   }
   else
   {
      assert(row->len - row->nlpcols >= 1);
      assert(row->cols[row->len-1] == col);
      if( row->len - row->nlpcols > 1 )
      {
         assert(row->cols_index[row->len-2] == row->cols[row->len-2]->index);
         row->nonlpcolssorted = row->nonlpcolssorted && (row->cols_index[row->len-2] < col->index);
      }
   }

   coefChangedExact(row, col, lp);

   SCIPrationalDebugMessage("added coefficient %q * <%s> at position %d (%d/%d) to row <%s> (nunlinked=%d)\n",
      val, SCIPvarGetName(col->var), pos, row->nlpcols, row->len, row->fprow->name, row->nunlinked);

   return SCIP_OKAY;
}

/** deletes coefficient at given position from row */
static
SCIP_RETCODE rowExactDelCoefPos(
   SCIP_ROWEXACT*        row,                /**< row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   int                   pos                 /**< position in row vector to delete */
   )
{
   SCIP_COLEXACT* col;

   assert(row != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < row->len);
   assert(row->cols[pos] != NULL);
   assert((pos < row->nlpcols) == (row->linkpos[pos] >= 0 && row->cols[pos]->lppos >= 0));

   col = row->cols[pos];

   assert((pos < row->nlpcols) == (col->lppos >= 0 && row->linkpos[pos] >= 0));

   SCIPrationalDebugMessage("deleting coefficient %q * <%s> at position %d from row <%s>\n",
     row->vals[pos], SCIPvarGetName(col->var), pos, row->fprow->name);

   if( row->nlocks > 0 )
   {
      SCIPerrorMessage("cannot delete a coefficient from the locked unmodifiable row <%s>\n", row->fprow->name);
      return SCIP_INVALIDDATA;
   }

   if( row->linkpos[pos] == -1 )
      row->nunlinked--;

   /* if column is a linked LP column, move last linked LP coefficient to position of empty slot (deleted coefficient) */
   if( pos < row->nlpcols )
   {
      rowExactMoveCoef(row, row->nlpcols-1, pos);
      assert(!row->lpcolssorted);
      row->nlpcols--;
      pos = row->nlpcols;
   }

   /* move last coefficient to position of empty slot */
   rowExactMoveCoef(row, row->len-1, pos);
   row->len--;

   coefChangedExact(row, col, lp);

   return SCIP_OKAY;
}

/** changes a coefficient at given position of an LP row */
static
SCIP_RETCODE rowExactChgCoefPos(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   int                   pos,                /**< position in row vector to change */
   SCIP_RATIONAL*        val                 /**< value of coefficient */
   )
{
   SCIP_COLEXACT* col;

   assert(row != NULL);
   assert(lp != NULL);

   col = row->cols[pos];

   assert(col != NULL);
   assert(0 <= pos && pos < row->len);

   SCIPrationalDebugMessage("changing coefficient %q * <%s> at position %d of row <%s> to %q\n",
     row->vals[pos], SCIPvarGetName(row->cols[pos]->var), pos, row->fprow->name, val);

   if( row->nlocks > 0 )
   {
      SCIPerrorMessage("cannot change a coefficient of the locked unmodifiable row <%s>\n", row->fprow->name);
      return SCIP_INVALIDDATA;
   }

   assert(col != NULL);

   if( SCIPrationalIsZero(val) )
   {
      /* delete existing coefficient */
      SCIP_CALL( rowExactDelCoefPos(row, set, lp, pos) );
   }
   else if( !SCIPrationalIsEQ(row->vals[pos], val) )
   {
      /* change existing coefficient */
      SCIPrationalSetRational(row->vals[pos], val);
      SCIPintervalSetRational(&row->valsinterval[pos], val);
      row->integral = row->integral && SCIPcolIsIntegral(col->fpcol) && SCIPrationalIsIntegral(val);
      coefChangedExact(row, col, lp);
   }

   return SCIP_OKAY;
}

/*
 * double linked coefficient matrix methods
 */

/** insert column coefficients in corresponding rows */
static
SCIP_RETCODE colExactLink(
   SCIP_COLEXACT*        col,                /**< column data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lp                  /**< current LP data */
   )
{
   int i;

   assert(col != NULL);
   assert(col->fpcol->var != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( col->nunlinked > 0 )
   {
      SCIPsetDebugMsg(set, "linking column <%s>\n", SCIPvarGetName(col->var));

      /* unlinked rows can only be in the non-LP/unlinked rows part of the rows array */
      for( i = col->nlprows; i < col->len; ++i )
      {
         assert(!SCIPrationalIsZero(col->vals[i]));
         if( col->linkpos[i] == -1 )
         {
            /* this call might swap the current row with the first non-LP/not linked row, but this is of no harm */
            SCIP_CALL( rowExactAddCoef(col->rows[i], blkmem, set, eventqueue, lp, col, col->vals[i], i) );
         }
         assert(col->rows[i]->cols[col->linkpos[i]] == col);
         assert(col->rows[i]->linkpos[col->linkpos[i]] == i);
         assert(col->nlprows == 0 || col->rows[col->nlprows-1]->cols[col->linkpos[col->nlprows-1]] == col);
         assert(col->nlprows == 0 || col->rows[col->nlprows-1]->linkpos[col->linkpos[col->nlprows-1]] == col->nlprows-1);
      }
   }
   assert(col->nunlinked == 0);

   checkLinks(lp);

   return SCIP_OKAY;
}

/** insert row coefficients in corresponding columns */
static
SCIP_RETCODE rowExactLink(
   SCIP_ROWEXACT*        row,                /**< row data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lp                  /**< current LP data */
   )
{
   int i;

   assert(row != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( row->nunlinked > 0 )
   {
      SCIPsetDebugMsg(set, "linking row <%s>\n", row->fprow->name);

      /* unlinked columns can only be in the non-LP/unlinked columns part of the cols array */
      for( i = row->nlpcols; i < row->len; ++i )
      {
         assert(!SCIPrationalIsZero(row->vals[i]));
         if( row->linkpos[i] == -1 )
         {
            /* this call might swap the current column with the first non-LP/not linked column, but this is of no harm */
            SCIP_CALL( colExactAddCoef(row->cols[i], blkmem, set, eventqueue, lp, row, row->vals[i], i) );
         }
         assert(row->cols[i]->rows[row->linkpos[i]] == row);
         assert(row->cols[i]->linkpos[row->linkpos[i]] == i);
         assert(row->nlpcols == 0 || row->cols[row->nlpcols-1]->rows[row->linkpos[row->nlpcols-1]] == row);
         assert(row->nlpcols == 0 || row->cols[row->nlpcols-1]->linkpos[row->linkpos[row->nlpcols-1]] == row->nlpcols-1);
      }
   }
   assert(row->nunlinked == 0);

   checkLinks(lp);

   return SCIP_OKAY;
}

/** removes row coefficients from corresponding columns */
static
SCIP_RETCODE rowExactUnlink(
   SCIP_ROWEXACT*        row,                /**< row data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lp                  /**< current LP data */
   )
{
   int i;

   assert(row != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( row->nunlinked < row->len )
   {
      SCIPsetDebugMsg(set, "unlinking exact row <%p>\n", (void*) row);
      for( i = 0; i < row->len; ++i )
      {
         if( row->linkpos[i] >= 0 )
         {
            assert(row->cols[i]->rows[row->linkpos[i]] == row);
            SCIP_CALL( colExactDelCoefPos(row->cols[i], set, lp, row->linkpos[i]) );
            row->nunlinked++;
         }
      }
   }
   assert(row->nunlinked == row->len);

   return SCIP_OKAY;
}

/** updates link data after addition of column */
static
SCIP_RETCODE colExactUpdateAddLP(
   SCIP_COLEXACT*        col,                /**< LP column */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_ROWEXACT* row;
   int i;
   int pos;

   assert(col != NULL);
   assert(col->lppos >= 0);

   /* update column arrays of all linked rows */
   for( i = 0; i < col->len; ++i )
   {
      pos = col->linkpos[i];
      if( pos >= 0 )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->linkpos[pos] == i);
         assert(row->cols[pos] == col);
         assert(row->nlpcols <= pos && pos < row->len);

         row->nlpcols++;
         SCIP_CALL( rowExactSwapCoefs(row, set->buffer, pos, row->nlpcols-1) );
         assert(row->cols[row->nlpcols-1] == col);

         /* if no swap was necessary, mark lpcols to be unsorted */
         if( pos == row->nlpcols-1 )
            row->lpcolssorted = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** updates link data after addition of row */
static
SCIP_RETCODE rowExactUpdateAddLP(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_COLEXACT* col;
   int i;
   int pos;

   assert(row != NULL);
   assert(row->lppos >= 0);

   /* update row arrays of all linked columns */
   for( i = 0; i < row->len; ++i )
   {
      pos = row->linkpos[i];
      if( pos >= 0 )
      {
         col = row->cols[i];
         assert(col != NULL);
         assert(col->linkpos[pos] == i);
         assert(col->rows[pos] == row);
         assert(col->nlprows <= pos && pos < col->len);

         col->nlprows++;
         SCIP_CALL( colExactSwapCoefs(col, set->buffer, pos, col->nlprows-1) );

         /* if no swap was necessary, mark lprows to be unsorted */
         if( pos == col->nlprows-1 )
            col->lprowssorted = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** updates link data after removal of column */
static
SCIP_RETCODE colExactUpdateDelLP(
   SCIP_COLEXACT*        col,                /**< LP column */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_ROWEXACT* row;
   int i;
   int pos;

   assert(col != NULL);
   assert(col->lppos == -1);

   /* update column arrays of all linked rows */
   for( i = 0; i < col->len; ++i )
   {
      pos = col->linkpos[i];
      if( pos >= 0 )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->linkpos[pos] == i);
         assert(row->cols[pos] == col);
         assert(0 <= pos && pos < row->nlpcols);

         row->nlpcols--;
         SCIP_CALL( rowExactSwapCoefs(row, set->buffer, pos, row->nlpcols) );

         /* if no swap was necessary, mark nonlpcols to be unsorted */
         if( pos == row->nlpcols )
            row->nonlpcolssorted = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** updates link data after removal of row */
static
SCIP_RETCODE rowExactUpdateDelLP(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_COLEXACT* col;
   int i;
   int pos;

   assert(row != NULL);
   assert(row->lppos == -1);

   /* update row arrays of all linked columns */
   for( i = 0; i < row->len; ++i )
   {
      pos = row->linkpos[i];
      if( pos >= 0 )
      {
         col = row->cols[i];
         assert(col != NULL);
         assert(0 <= pos && pos < col->nlprows);
         assert(col->linkpos[pos] == i);
         assert(col->rows[pos] == row);

         col->nlprows--;
         SCIP_CALL( colExactSwapCoefs(col, set->buffer, pos, col->nlprows) );

         /* if no swap was necessary, mark lprows to be unsorted */
         if( pos == col->nlprows )
            col->nonlprowssorted = FALSE;
      }
   }

   return SCIP_OKAY;
}

/*
 * flushing methods
 */

/** resets column data to represent a column not in the LP solver */
static
void markColexDeleted(
   SCIP_COLEXACT*        col                 /**< column to be marked deleted */
   )
{
   assert(col != NULL);

   col->lpipos = -1;
   col->validredcostlp = -1;
   col->validfarkaslp = -1;
   col->basisstatus = SCIP_BASESTAT_ZERO; /*lint !e641*/
}

/** applies all cached column removals to the LP solver */
static
SCIP_RETCODE lpExactFlushDelCols(
   SCIP_LPEXACT*         lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->lpifirstchgcol <= lp->nlpicols);
   assert(lp->lpifirstchgcol <= lp->ncols);

   /* find the first column to change */
   while( lp->lpifirstchgcol < lp->nlpicols
      && lp->lpifirstchgcol < lp->ncols
      && lp->cols[lp->lpifirstchgcol]->lpipos == lp->lpifirstchgcol
      && !lp->cols[lp->lpifirstchgcol]->coefchanged )
   {
      assert(lp->cols[lp->lpifirstchgcol] == lp->lpicols[lp->lpifirstchgcol]);
      lp->lpifirstchgcol++;
   }

   /* shrink LP to the part which didn't change */
   if( lp->lpifirstchgcol < lp->nlpicols )
   {
      int i;

      assert(!lp->fplp->diving);
      SCIPdebugMessage("flushing col deletions: shrink exact LP from %d to %d columns\n", lp->nlpicols, lp->lpifirstchgcol);
      SCIP_CALL( SCIPlpiExactDelCols(lp->lpiexact, lp->lpifirstchgcol, lp->nlpicols-1) );
      for( i = lp->lpifirstchgcol; i < lp->nlpicols; ++i )
      {
         markColexDeleted(lp->lpicols[i]);
      }
      lp->nlpicols = lp->lpifirstchgcol;
      lp->flushdeletedcols = TRUE;
      lp->updateintegrality = TRUE;

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->primalchecked = FALSE;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }
   assert(lp->nlpicols == lp->lpifirstchgcol);

   return SCIP_OKAY;
}

/** applies all cached column additions to the LP solver */
static
SCIP_RETCODE lpExactFlushAddCols(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_RATIONAL** obj;
   SCIP_RATIONAL** lb;
   SCIP_RATIONAL** ub;
   int* beg;
   int* ind;
   SCIP_RATIONAL** val;
   char** name;
   SCIP_COLEXACT* col;
   int c;
   int pos;
   int nnonz;
   int naddcols;
   int naddcoefs;
   int i;
   int lpipos;

   assert(lp != NULL);
   assert(lp->lpifirstchgcol == lp->nlpicols);
   assert(blkmem != NULL);
   assert(set != NULL);

   /* if there are no columns to add, we are ready */
   if( lp->ncols == lp->nlpicols )
      return SCIP_OKAY;

   /* add the additional columns */
   assert(lp->ncols > lp->nlpicols);
   SCIP_CALL( ensureLpiExactcolsSize(lp, set, lp->ncols) );

   /* count the (maximal) number of added coefficients, calculate the number of added columns */
   naddcols = lp->ncols - lp->nlpicols;
   naddcoefs = 0;
   for( c = lp->nlpicols; c < lp->ncols; ++c )
      naddcoefs += lp->cols[c]->len;
   assert(naddcols > 0);

   /* get temporary memory for changes */
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &obj, naddcols) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &lb, naddcols) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &ub, naddcols) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &val, naddcoefs) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &beg, naddcols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ind, naddcoefs) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &name, naddcols) );

   /* fill temporary memory with column data */
   nnonz = 0;
   for( pos = 0, c = lp->nlpicols; c < lp->ncols; ++pos, ++c )
   {
      col = lp->cols[c];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatusExact(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetColExact(col->var) == col);
      assert(col->lppos == c);
      assert(nnonz + col->nlprows <= naddcoefs);

      SCIPsetDebugMsg(set, "flushing added column <%s>: ", SCIPvarGetName(col->var));

      /* Because the column becomes a member of the LP solver, it now can take values
       * different from zero. That means, we have to include the column in the corresponding
       * row vectors.
       */
      SCIP_CALL( colExactLink(col, blkmem, set, eventqueue, lp) );

      lp->lpicols[c] = col;
      col->lpipos = c;
      col->validredcostlp = -1;
      col->validfarkaslp = -1;
      col->objchanged = FALSE;
      col->lbchanged = FALSE;
      col->ubchanged = FALSE;
      col->coefchanged = FALSE;
      SCIPrationalSetRational(obj[pos], col->obj);
      SCIPrationalSetRational(lb[pos], col->lb);
      SCIPrationalSetRational(ub[pos], col->ub);

      beg[pos] = nnonz;
      name[pos] = (char*)SCIPvarGetName(col->fpcol->var);

      SCIPrationalSetRational(col->flushedobj, obj[pos]);
      SCIPrationalSetRational(col->flushedlb, lb[pos]);
      SCIPrationalSetRational(col->flushedub, ub[pos]);

      for( i = 0; i < col->nlprows; ++i )
      {
         assert(col->rows[i] != NULL);
         lpipos = col->rows[i]->lpipos;
         if( lpipos >= 0 )
         {
            assert(lpipos < lp->nrows);
            assert(nnonz < naddcoefs);
            ind[nnonz] = lpipos;
            SCIPrationalSetRational(val[nnonz], col->vals[i]);
            nnonz++;
         }
      }
#ifndef NDEBUG
      for( i = col->nlprows; i < col->len; ++i )
      {
         assert(col->rows[i] != NULL);
         assert(col->rows[i]->lpipos == -1); /* because the row deletions are already performed */
      }
#endif
   }

   /* call LP interface */
   SCIPsetDebugMsg(set, "flushing col additions: enlarge exact LP from %d to %d columns\n", lp->nlpicols, lp->ncols);
   SCIP_CALL( SCIPlpiExactAddCols(lp->lpiexact, naddcols, obj, lb, ub, name, nnonz, beg, ind, val) );
   lp->nlpicols = lp->ncols;
   lp->lpifirstchgcol = lp->nlpicols;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &name);
   SCIPsetFreeBufferArray(set, &ind);
   SCIPsetFreeBufferArray(set, &beg);
   SCIPrationalFreeBufferArray(set->buffer, &val, naddcoefs);
   SCIPrationalFreeBufferArray(set->buffer, &ub, naddcols);
   SCIPrationalFreeBufferArray(set->buffer, &lb, naddcols);
   SCIPrationalFreeBufferArray(set->buffer, &obj, naddcols);

   lp->flushaddedcols = TRUE;
   lp->updateintegrality = TRUE;

   /* mark the LP unsolved */
   lp->solved = FALSE;
   lp->dualfeasible = FALSE;
   lp->dualchecked = FALSE;
   lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

   return SCIP_OKAY;
}

/** resets row data to represent a row not in the LP solver */
static
void markRowexDeleted(
   SCIP_ROWEXACT*        row                 /**< row to be marked deleted */
   )
{
   assert(row != NULL);

   row->lpipos = -1;
   row->basisstatus = SCIP_BASESTAT_BASIC; /*lint !e641*/
   row->validactivitylp = -1;
}

/** applies all cached row removals to the LP solver */
static
SCIP_RETCODE lpExactFlushDelRows(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(lp->lpifirstchgrow <= lp->nlpirows);
   assert(lp->lpifirstchgrow <= lp->nrows);

   /* find the first row to change */
   while( lp->lpifirstchgrow < lp->nlpirows
      && lp->lpifirstchgrow < lp->nrows
      && lp->rows[lp->lpifirstchgrow]->lpipos == lp->lpifirstchgrow
      && !lp->rows[lp->lpifirstchgrow]->coefchanged )
   {
      assert(lp->rows[lp->lpifirstchgrow] == lp->lpirows[lp->lpifirstchgrow]);
      lp->lpifirstchgrow++;
   }

   /* shrink LP to the part which didn't change */
   if( lp->lpifirstchgrow < lp->nlpirows )
   {
      int i;

      SCIPsetDebugMsg(set, "flushing row deletions: shrink exact LP from %d to %d rows\n", lp->nlpirows, lp->lpifirstchgrow);
      SCIP_CALL( SCIPlpiExactDelRows(lp->lpiexact, lp->lpifirstchgrow, lp->nlpirows-1) );
      for( i = lp->lpifirstchgrow; i < lp->nlpirows; ++i )
      {
         markRowexDeleted(lp->lpirows[i]);
         SCIP_CALL( SCIProwExactRelease(&lp->lpirows[i], blkmem, set, lp) );
      }
      lp->nlpirows = lp->lpifirstchgrow;
      lp->flushdeletedrows = TRUE;

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->dualchecked = FALSE;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }
   assert(lp->nlpirows == lp->lpifirstchgrow);

   return SCIP_OKAY;
}

/** applies all cached row additions and removals to the LP solver */
static
SCIP_RETCODE lpExactFlushAddRows(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_RATIONAL** lhs;
   SCIP_RATIONAL** rhs;
   SCIP_RATIONAL** val;
   int* beg;
   int* ind;
   char** name;
   SCIP_ROWEXACT* row;
   int r;
   int pos;
   int nnonz;
   int naddrows;
   int naddcoefs;
   int i;
   int lpipos;

   assert(lp != NULL);
   assert(lp->lpifirstchgrow == lp->nlpirows);
   assert(blkmem != NULL);

   /* if there are no rows to add, we are ready */
   if( lp->nrows == lp->nlpirows )
      return SCIP_OKAY;

   /* add the additional rows */
   assert(lp->nrows > lp->nlpirows);
   SCIP_CALL( ensureLpirowexactsSize(lp, set, lp->nrows) );

   /* count the (maximal) number of added coefficients, calculate the number of added rows */
   naddrows = lp->nrows - lp->nlpirows;
   naddcoefs = 0;
   for( r = lp->nlpirows; r < lp->nrows; ++r )
      naddcoefs += lp->rows[r]->len;
   assert(naddrows > 0);

   /* get temporary memory for changes */
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &lhs, naddrows) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &rhs, naddrows) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &val, naddcoefs) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &beg, naddrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ind, naddcoefs) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &name, naddrows) );

   /* fill temporary memory with row data */
   nnonz = 0;
   for( pos = 0, r = lp->nlpirows; r < lp->nrows; ++pos, ++r )
   {
      row = lp->rows[r];
      assert(row != NULL);
      assert(row->lppos == r);
      assert(nnonz + row->nlpcols <= naddcoefs);

      SCIPsetDebugMsg(set, "flushing added exact row <%s>: ", row->fprow->name);

      /* Because the row becomes a member of the LP solver, its dual variable now can take values
       * different from zero. That means, we have to include the row in the corresponding
       * column vectors.
       */
      SCIP_CALL( rowExactLink(row, blkmem, set, eventqueue, lp) );

      SCIProwExactCapture(row);
      lp->lpirows[r] = row;
      row->lpipos = r;
      row->lhschanged = FALSE;
      row->rhschanged = FALSE;
      row->coefchanged = FALSE;

      SCIPrationalDiff(lhs[pos], row->lhs, row->constant);
      SCIPrationalDiff(rhs[pos], row->rhs, row->constant);
      beg[pos] = nnonz;
      name[pos] = row->fprow->name;

      SCIPrationalSetRational(row->flushedlhs, lhs[pos]);
      SCIPrationalSetRational(row->flushedrhs, rhs[pos]);

      SCIPrationalDebugMessage("flushing added row (SCIP_LPI): %q <=", lhs[pos]);
      for( i = 0; i < row->nlpcols; ++i )
      {
         assert(row->cols[i] != NULL);
         lpipos = row->cols[i]->lpipos;
         if( lpipos >= 0 )
         {
            assert(lpipos < lp->ncols);
            assert(nnonz < naddcoefs);
            SCIPrationalDebugMessage(" %q %d(<%s>)", row->vals[i], lpipos+1, SCIPvarGetName(row->cols[i]->fpcol->var));
            ind[nnonz] = lpipos;
            SCIPrationalSetRational(val[nnonz], row->vals[i]);
            nnonz++;
         }
      }
      SCIPrationalDebugMessage(" <= %q\n", rhs[pos]);
#ifndef NDEBUG
      for( i = row->nlpcols; i < row->len; ++i )
      {
         assert(row->cols[i] != NULL);
         assert(row->cols[i]->lpipos == -1); /* because the column deletions are already performed */
      }
#endif
   }

   /* call LP interface */
   SCIPsetDebugMsg(set, "flushing row additions: enlarge LP from %d to %d rows\n", lp->nlpirows, lp->nrows);
   SCIP_CALL( SCIPlpiExactAddRows(lp->lpiexact, naddrows, lhs, rhs, name, nnonz, beg, ind, val) );
   lp->nlpirows = lp->nrows;
   lp->lpifirstchgrow = lp->nlpirows;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &name);
   SCIPsetFreeBufferArray(set, &ind);
   SCIPsetFreeBufferArray(set, &beg);
   SCIPrationalFreeBufferArray(set->buffer, &val, naddcoefs);
   SCIPrationalFreeBufferArray(set->buffer, &rhs, naddrows);
   SCIPrationalFreeBufferArray(set->buffer, &lhs, naddrows);

   lp->flushaddedrows = TRUE;

   /* mark the LP unsolved */
   lp->solved = FALSE;
   lp->primalfeasible = FALSE;
   lp->primalchecked = FALSE;
   lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

   return SCIP_OKAY;
}

/** applies all cached column bound and objective changes to the LP */
static
SCIP_RETCODE lpExactFlushChgCols(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
#ifndef NDEBUG
   SCIP_Bool lpinone = (strcmp( SCIPlpiGetSolverName(), "NONE") == 0);
#endif
   SCIP_COLEXACT* col;
   int* objind;
   int* bdind;
   SCIP_RATIONAL** obj;
   SCIP_RATIONAL** lb;
   SCIP_RATIONAL** ub;
   int nobjchg;
   int nbdchg;
   int i;

   assert(lp != NULL);

   if( lp->nchgcols == 0 )
      return SCIP_OKAY;

   /* get temporary memory for changes */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &objind, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bdind, lp->ncols) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &obj, lp->ncols) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &lb, lp->ncols) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &ub, lp->ncols) );

   /* collect all cached bound and objective changes */
   nobjchg = 0;
   nbdchg = 0;
   for( i = 0; i < lp->nchgcols; ++i )
   {
      col = lp->chgcols[i];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatusExact(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetColExact(col->var) == col);

      if( col->lpipos >= 0 )
      {
#ifndef NDEBUG
         /* do not check consistency of data with LPI in case of LPI=none */
         if( !lpinone )
         {
            SCIP_RATIONAL* lpiobj;
            SCIP_RATIONAL* lpiub;
            SCIP_RATIONAL* lpilb;

            SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &lpiobj) );
            SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &lpilb) );
            SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &lpiub) );

            SCIP_CALL( SCIPlpiExactGetObj(lp->lpiexact, col->lpipos, col->lpipos, &lpiobj) );
            SCIP_CALL( SCIPlpiExactGetBounds(lp->lpiexact, col->lpipos, col->lpipos, &lpilb, &lpiub) );
            assert(SCIPrationalIsEQ(lpiobj, col->flushedobj));
            SCIPrationalFreeBuffer(set->buffer, &lpiub);
            SCIPrationalFreeBuffer(set->buffer, &lpilb);
            SCIPrationalFreeBuffer(set->buffer, &lpiobj);
         }
#endif

         if( col->objchanged )
         {
            if( SCIPrationalIsEQ(col->flushedobj, col->obj) ) /*lint !e777*/
            {
               assert(nobjchg < lp->ncols);
               objind[nobjchg] = col->lpipos;
               SCIPrationalSetRational(obj[nobjchg], col->obj);
               nobjchg++;
               SCIPrationalSetRational(col->flushedobj, col->obj);
            }
            col->objchanged = FALSE;
         }

         if( col->lbchanged || col->ubchanged )
         {
            if( !SCIPrationalIsEQ(col->flushedlb, col->lb) || !SCIPrationalIsEQ(col->flushedub, col->ub) ) /*lint !e777*/
            {
               assert(nbdchg < lp->ncols);
               bdind[nbdchg] = col->lpipos;
               SCIPrationalSetRational(lb[nbdchg], col->lb);
               SCIPrationalSetRational(ub[nbdchg], col->ub);
               nbdchg++;
               SCIPrationalSetRational(col->flushedlb, col->lb);
               SCIPrationalSetRational(col->flushedub, col->ub);
            }
            col->lbchanged = FALSE;
            col->ubchanged = FALSE;
         }
      }
      /* maybe lb/ub/objchanged should all be set to false when lpipos is -1 */
   }

   /* change objective values in LP */
   if( nobjchg > 0 )
   {
      SCIPsetDebugMsg(set, "flushing objective changes: change %d objective values of %d changed columns\n", nobjchg, lp->nchgcols);
      SCIP_CALL( SCIPlpiExactChgObj(lp->lpiexact, nobjchg, objind, obj) );

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->dualchecked = FALSE;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   /* change bounds in LP */
   if( nbdchg > 0 )
   {
      SCIPsetDebugMsg(set, "flushing bound changes: change %d bounds of %d changed columns\n", nbdchg, lp->nchgcols);
      SCIP_CALL( SCIPlpiExactChgBounds(lp->lpiexact, nbdchg, bdind, lb, ub) );

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->primalchecked = FALSE;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   lp->nchgcols = 0;

   /* free temporary memory */
   SCIPrationalFreeBufferArray(set->buffer, &ub, lp->ncols);
   SCIPrationalFreeBufferArray(set->buffer, &lb, lp->ncols);
   SCIPrationalFreeBufferArray(set->buffer, &obj, lp->ncols);
   SCIPsetFreeBufferArray(set, &bdind);
   SCIPsetFreeBufferArray(set, &objind);

   return SCIP_OKAY;
}

/** applies all cached row side changes to the LP */
static
SCIP_RETCODE lpExactFlushChgRows(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
#ifndef NDEBUG
   SCIP_Bool lpinone = (strcmp( SCIPlpiGetSolverName(), "NONE") == 0);
#endif
   SCIP_ROWEXACT* row;
   int* ind;
   SCIP_RATIONAL** lhs;
   SCIP_RATIONAL** rhs;
   int i;
   int nchg;

   assert(lp != NULL);

   if( lp->nchgrows == 0 )
      return SCIP_OKAY;

   /* get temporary memory for changes */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ind, lp->nrows) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &lhs, lp->nrows) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &rhs, lp->nrows) );

   /* collect all cached left and right hand side changes */
   nchg = 0;
   for( i = 0; i < lp->nchgrows; ++i )
   {
      row = lp->chgrows[i];
      assert(row != NULL);

      if( row->lpipos >= 0 )
      {
#ifndef NDEBUG
         /* do not check consistency of data with LPI in case of LPI=none */
         if( !lpinone )
         {
            SCIP_RATIONAL* lpirhs;
            SCIP_RATIONAL* lpilhs;

            SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &lpilhs) );
            SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &lpirhs) );

            SCIP_CALL( SCIPlpiExactGetSides(lp->lpiexact, row->lpipos, row->lpipos, &lpilhs, &lpirhs) );
            assert(SCIPrationalIsEQ(lpilhs, row->flushedlhs));
            assert(SCIPrationalIsEQ(lpirhs, row->flushedrhs));

            SCIPrationalFreeBuffer(set->buffer, &lpirhs);
            SCIPrationalFreeBuffer(set->buffer, &lpilhs);
         }
#endif
         if( row->lhschanged || row->rhschanged )
         {
            SCIP_RATIONAL* newlhs;
            SCIP_RATIONAL* newrhs;

            SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &newlhs) );
            SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &newrhs) );

            SCIPrationalDiff(newlhs, row->lhs, row->constant);
            SCIPrationalDiff(newrhs, row->rhs, row->constant);
            if( SCIPrationalIsEQ(row->flushedlhs, newlhs) || SCIPrationalIsEQ(row->flushedrhs, newrhs) ) /*lint !e777*/
            {
               assert(nchg < lp->nrows);
               ind[nchg] = row->lpipos;
               SCIPrationalSetRational(lhs[nchg], newlhs);
               SCIPrationalSetRational(rhs[nchg], newrhs);
               nchg++;
               SCIPrationalSetRational(row->flushedlhs, newlhs);
               SCIPrationalSetRational(row->flushedrhs, newrhs);
            }
            row->lhschanged = FALSE;
            row->rhschanged = FALSE;

            SCIPrationalFreeBuffer(set->buffer, &newrhs);
            SCIPrationalFreeBuffer(set->buffer, &newlhs);
         }
      }
   }

   /* change left and right hand sides in LP */
   if( nchg > 0 )
   {
      SCIPsetDebugMsg(set, "flushing side changes: change %d sides of %d exact rows\n", nchg, lp->nchgrows);
      SCIP_CALL( SCIPlpiExactChgSides(lp->lpiexact, nchg, ind, lhs, rhs) );

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->primalchecked = FALSE;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   lp->nchgrows = 0;

   /* free temporary memory */
   SCIPrationalFreeBufferArray(set->buffer, &rhs, lp->nrows);
   SCIPrationalFreeBufferArray(set->buffer, &lhs, lp->nrows);
   SCIPsetFreeBufferArray(set, &ind);

   return SCIP_OKAY;
}

/** gets finite part of objective value of current LP that results from LOOSE variables only.
 * returns reference, so be careful not to change!
 */
static
SCIP_RATIONAL* getFiniteLooseObjvalExact(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(lp != NULL);
   assert(set != NULL);
   assert(prob != NULL);
   assert((lp->nloosevars > 0) || (lp->looseobjvalinf == 0 && SCIPrationalIsZero(lp->looseobjval)));
   assert(lp->looseobjvalinf == 0);

   return lp->looseobjval;
}

/*
 * Column methods
 */

/** creates an LP column */
SCIP_RETCODE SCIPcolExactCreate(
   SCIP_COLEXACT**       col,                /**< pointer to column data */
   SCIP_COL*             fpcol,              /**< the corresponding fp col */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             var,                /**< variable, this column represents */
   int                   len,                /**< number of nonzeros in the column */
   SCIP_ROWEXACT**       rows,               /**< array with rows of column entries */
   SCIP_RATIONAL**       vals,               /**< array with coefficients of column entries */
   SCIP_Bool             removable           /**< should the column be removed from the LP due to aging or cleanup? */
   )
{
   int i;

   assert(col != NULL);
   assert(fpcol != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(var != NULL);
   assert(len >= 0);
   assert(len == 0 || (vals != NULL));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, col) );

   (*col)->fpcol = fpcol;

   if( len > 0 )
   {
      SCIP_CALL( SCIPrationalCopyBlockArray(blkmem, &(*col)->vals, vals, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*col)->linkpos, len) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*col)->rows, rows, len) );

      for( i = 0; i < len; ++i )
      {
         assert(!SCIPrationalIsZero(vals[i]));
         assert(rows[i] != NULL);
         (*col)->linkpos[i] = -1;
      }
   }
   else
   {
      (*col)->rows = NULL;
      (*col)->vals = NULL;
      (*col)->linkpos = NULL;
   }

   (*col)->var = var;
   SCIP_CALL( SCIPrationalCopyBlock(blkmem, &(*col)->obj, SCIPvarGetObjExact(var)) );
   SCIP_CALL( SCIPrationalCopyBlock(blkmem, &(*col)->lb, SCIPvarGetLbLocalExact(var)) );
   SCIP_CALL( SCIPrationalCopyBlock(blkmem, &(*col)->ub, SCIPvarGetUbLocalExact(var)) );
   (*col)->index = (*col)->fpcol->index;
   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &(*col)->flushedobj) );
   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &(*col)->flushedlb) );
   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &(*col)->flushedub) );
   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &(*col)->primsol) );
   SCIP_CALL( SCIPrationalCreateString(blkmem, &(*col)->redcost, "inf") );
   SCIP_CALL( SCIPrationalCreateString(blkmem, &(*col)->farkascoef, "inf") );

   (*col)->storedsolvals = NULL;
   (*col)->size = len;
   (*col)->len = len;
   (*col)->nlprows = 0;
   (*col)->lprowssorted = 0;
   (*col)->nunlinked = len;
   (*col)->lppos = -1;
   (*col)->lpipos = -1;
   (*col)->validredcostlp = -1;
   (*col)->validfarkaslp = -1;

   assert((*col)->fpcol->removable == removable);

   return SCIP_OKAY;
}

/** sets parameter of type SCIP_Real in exact LP solver, ignoring unknown parameters */
static
SCIP_RETCODE lpExactSetRealpar(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   SCIP_Real             value,              /**< value to set parameter to */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   SCIP_RETCODE retcode;

   assert(lp != NULL);
   assert(success != NULL);

   retcode = SCIPlpiExactSetRealpar(lp->lpiexact, lpparam, value);

   /* check, if parameter is unknown */
   if( retcode == SCIP_PARAMETERUNKNOWN )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   *success = TRUE;

   return retcode;
}

/** sets parameter of type SCIP_Real in exact LP solver, ignoring unknown parameters */
static
SCIP_RETCODE lpExactSetIntpar(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   int                   value,              /**< value to set parameter to */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   SCIP_RETCODE retcode;

   assert(lp != NULL);
   assert(success != NULL);

   retcode = SCIPlpiExactSetIntpar(lp->lpiexact, lpparam, value);

   /* check, if parameter is unknown */
   if( retcode == SCIP_PARAMETERUNKNOWN )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   *success = TRUE;

   return retcode;
}

/** sets the objective limit of the exact LP solver
 *
 *  Note that we are always minimizing.
 */
static
SCIP_RETCODE lpExactSetObjlim(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             objlim,             /**< new objective limit */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was actually changed */
   )
{
   assert(lp != NULL);
   assert(set != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* We disabled the objective limit in the LP solver or we want so solve exactly and thus cannot rely on the LP
    * solver's objective limit handling, so we return here and do not apply the objective limit. */
   if( set->lp_disablecutoff == 1 || (set->nactivepricers > 0 && set->lp_disablecutoff == 2) )
      return SCIP_OKAY;

   /* convert SCIP infinity value to lp-solver infinity value if necessary */
   if( SCIPsetIsInfinity(set, objlim) )
      objlim = SCIPlpiExactInfinity(lp->lpiexact);

   if( objlim != lp->lpiobjlim ) /*lint !e777*/
   {
      SCIP_CALL( lpExactSetRealpar(lp, SCIP_LPPAR_OBJLIM, objlim, success) );
      if( *success )
      {
         SCIP_Real actualobjlim;

         /* check whether the parameter was actually changed or already was at the boundary of the LP solver's parameter range */
         SCIP_CALL( SCIPlpiExactGetRealpar(lp->lpiexact, SCIP_LPPAR_OBJLIM, &actualobjlim) );
         if( actualobjlim != lp->lpiobjlim ) /*lint !e777*/
         {
            /* mark the current solution invalid */
            lp->solved = FALSE;
            lp->primalfeasible = FALSE;
            lp->primalchecked = FALSE;
            SCIPrationalSetInfinity(lp->lpobjval);
            lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         }
         lp->lpiobjlim = actualobjlim;
      }
   }

   return SCIP_OKAY;
}

/** sets the iteration limit of the LP solver */
static
SCIP_RETCODE lpExactSetIterationLimit(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   int                   itlim               /**< maximal number of LP iterations to perform, or -1 for no limit */
   )
{
   SCIP_Bool success;

   assert(lp != NULL);
   assert(itlim >= -1);

   if( itlim == -1 )
      itlim = INT_MAX;

   if( itlim != lp->lpiitlim )
   {
      SCIP_CALL( lpExactSetIntpar(lp, SCIP_LPPAR_LPITLIM, itlim, &success) );
      if( success )
      {
         if( itlim > lp->lpiitlim )
         {
            /* mark the current solution invalid */
            lp->solved = FALSE;
            SCIPrationalSetInfinity(lp->lpobjval);
            lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         }
         lp->lpiitlim = itlim;
      }
   }

   return SCIP_OKAY;
}

/** resets row data to represent a row not in the LP solver */
static
void markRowExactDeleted(
   SCIP_ROWEXACT*        row                 /**< row to be marked deleted */
   )
{
   assert(row != NULL);

   row->lpipos = -1;
   SCIPrationalSetReal(row->dualsol, 0.0);
   SCIPrationalSetInfinity(row->activity);
   SCIPrationalSetReal(row->dualfarkas, 0.0);
   row->basisstatus = SCIP_BASESTAT_BASIC; /*lint !e641*/
   row->validactivitylp = -1;
}

/** deletes the marked rows from the LP and the LP interface */
SCIP_RETCODE SCIPlpExactDelRowset(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   int*                  rowdstat            /**< deletion status of rows:  1 if row should be deleted, 0 if not */
   )
{
   SCIP_ROWEXACT* row;
   int nrows;
   int nlpirows;
   int r;
   int c;

   assert(lpexact != NULL);
   assert(rowdstat != NULL);

   nrows = lpexact->nrows;
   nlpirows = lpexact->nlpirows;
   if( nlpirows == 0 )
      return SCIP_OKAY;

   /* delete rows in LP solver */
   SCIP_CALL( SCIPlpiExactDelRowset(lpexact->lpiexact, rowdstat) );

   /* set the correct status for rows that never made it to the lpi (this is special for the exact lp) */
   c = lpexact->nlpirows - 1;
   while( rowdstat[c] == -1 )
   {
      c--;
   }

   c = rowdstat[c] + 1;
   for( r = lpexact->nlpirows; r < nrows; r++ )
   {
      if( rowdstat[r] != 1 )
      {
         rowdstat[r] = c;
         ++c;
      }
      else
         rowdstat[r] = -1;
   }

   /* update LP data respectively */
   for( r = 0; r < nrows; ++r )
   {
      row = lpexact->rows[r];
      assert(rowdstat[r] <= r);
      assert(row != NULL);
      row->lppos = rowdstat[r];
      if( rowdstat[r] == -1 )
      {
         if( row->removable )
            lpexact->nremovablerows--;

         /* mark row to be deleted from the LPI and update row arrays of all linked columns */
         markRowExactDeleted(row);
         SCIP_CALL( rowExactUpdateDelLP(row, set) );
         row->lpdepth = -1;

         /* only release lpirows if they actually exist */
         if( r < nlpirows )
         {
            assert(row == lpexact->lpirows[r]);

            SCIP_CALL( SCIProwExactRelease(&lpexact->lpirows[r], blkmem, set, lpexact) );
            lpexact->nlpirows--;
         }
         SCIP_CALL( SCIProwExactRelease(&lpexact->rows[r], blkmem, set, lpexact) );
         assert(lpexact->rows[r] == NULL);
         lpexact->nrows--;
      }
      else if( rowdstat[r] < r )
      {
         assert(lpexact->rows[rowdstat[r]] == NULL);
         assert(lpexact->lpirows[rowdstat[r]] == NULL);
         lpexact->rows[rowdstat[r]] = row;

         /* only re-order lpirows if they actually exist */
         if( r < nlpirows )
         {
            lpexact->lpirows[rowdstat[r]] = row;
            lpexact->lpirows[r] = NULL;
         }
         lpexact->rows[rowdstat[r]]->lppos = rowdstat[r];
         lpexact->rows[rowdstat[r]]->lpipos = rowdstat[r];
         lpexact->rows[r] = NULL;
      }
   }

   /* mark LP to be unsolved */
   if( lpexact->nrows < nrows )
   {
      assert(lpexact->nchgrows == 0);

      lpexact->lpifirstchgrow = lpexact->nlpirows;

      /* mark the current solution invalid */
      lpexact->solved = FALSE;
      lpexact->dualfeasible = FALSE;
      lpexact->dualchecked = FALSE;
      SCIPrationalSetInfinity(lpexact->lpobjval);
      lpexact->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   checkLinks(lpexact);

   return SCIP_OKAY;
}

/** frees an LP column */
SCIP_RETCODE SCIPcolExactFree(
   SCIP_COLEXACT**       col,                /**< pointer to LP column */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(blkmem != NULL);
   assert(col != NULL);
   assert(*col != NULL);
   assert((*col)->fpcol != NULL);

   if( (*col)->size > 0 )
   {
      SCIPrationalFreeBlockArray(blkmem, &(*col)->vals, (*col)->size);
      BMSfreeBlockMemoryArray(blkmem, &(*col)->linkpos, (*col)->size);
      BMSfreeBlockMemoryArray(blkmem, &(*col)->rows, (*col)->size);
   }
   else
      assert((*col)->vals == NULL);

   if( (*col)->storedsolvals != NULL )
   {
      SCIPrationalFreeBlock(blkmem, &(*col)->storedsolvals->primsol);
      SCIPrationalFreeBlock(blkmem, &(*col)->storedsolvals->redcost);
      BMSfreeBlockMemoryNull(blkmem, &(*col)->storedsolvals);
   }

   SCIPrationalFreeBlock(blkmem, &(*col)->obj);
   SCIPrationalFreeBlock(blkmem, &(*col)->lb);
   SCIPrationalFreeBlock(blkmem, &(*col)->ub);
   SCIPrationalFreeBlock(blkmem, &(*col)->flushedobj);
   SCIPrationalFreeBlock(blkmem, &(*col)->flushedlb);
   SCIPrationalFreeBlock(blkmem, &(*col)->flushedub);
   SCIPrationalFreeBlock(blkmem, &(*col)->primsol);
   SCIPrationalFreeBlock(blkmem, &(*col)->redcost);
   SCIPrationalFreeBlock(blkmem, &(*col)->farkascoef);

   BMSfreeBlockMemory(blkmem, col);

   return SCIP_OKAY;
}

/** output column to file stream */
void SCIPcolExactPrint(
   SCIP_COLEXACT*        col,                /**< LP column */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int r;

   assert(col != NULL);
   assert(col->fpcol != NULL);
   assert(col->fpcol->var != NULL);

   SCIPmessageFPrintInfo(messagehdlr, file, "(obj:");
   SCIPrationalMessage(messagehdlr, file, col->obj);
   SCIPmessageFPrintInfo(messagehdlr, file, ") [");
   SCIPrationalMessage(messagehdlr, file, col->lb);
   SCIPmessageFPrintInfo(messagehdlr, file, ", ");
   SCIPrationalMessage(messagehdlr, file, col->ub);
   SCIPmessageFPrintInfo(messagehdlr, file, "], ");

   /* print coefficients */
   if( col->len == 0 )
      SCIPmessageFPrintInfo(messagehdlr, file, "<empty>");
   for( r = 0; r < col->len; ++r )
   {
      assert(col->rows[r] != NULL);
      assert(col->rows[r]->fprow->name != NULL);

      if( SCIPrationalIsPositive(col->vals[r]) )
         SCIPmessageFPrintInfo(messagehdlr, file, "+");

      SCIPrationalMessage(messagehdlr, file, col->vals[r]);
      SCIPmessageFPrintInfo(messagehdlr, file, "<%s> ", col->rows[r]->fprow->name);
   }
   SCIPmessageFPrintInfo(messagehdlr, file, "\n");
}

/** adds a previously non existing coefficient to an LP column */
SCIP_RETCODE SCIPcolExactAddCoef(
   SCIP_COLEXACT*        col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_RATIONAL*        val                 /**< value of coefficient */
   )
{
   assert(lpexact != NULL);

   SCIP_CALL( colExactAddCoef(col, blkmem, set, eventqueue, lpexact, row, val, -1) );

   return SCIP_OKAY;
}

/** deletes coefficient from column */
SCIP_RETCODE SCIPcolExactDelCoef(
   SCIP_COLEXACT*        col,                /**< column to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_ROWEXACT*        row                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(col != NULL);
   assert(col->var != NULL);
   assert(lpexact != NULL);
   assert(row != NULL);

   /* search the position of the row in the column's row vector */
   pos = colExactSearchCoef(col, row);
   if( pos == -1 )
   {
      SCIPerrorMessage("coefficient for row <%s> doesn't exist in column <%s>\n", row->fprow->name, SCIPvarGetName(col->var));
      return SCIP_INVALIDDATA;
   }
   assert(0 <= pos && pos < col->fpcol->len);
   assert(col->rows[pos] == row);

   /* if row knows of the column, remove the column from the row's col vector */
   if( col->linkpos[pos] >= 0 )
   {
      assert(row->cols[col->linkpos[pos]] == col);
      assert(row->cols_index[col->linkpos[pos]] == col->index);
      assert(SCIPrationalIsEQ(row->vals[col->linkpos[pos]], col->vals[pos]));
      SCIP_CALL( rowExactDelCoefPos(row, set, lpexact, col->linkpos[pos]) );
   }

   /* delete the row from the column's row vector */
   SCIP_CALL( colExactDelCoefPos(col, set, lpexact, pos) );

   checkLinks(lpexact);

   return SCIP_OKAY;
}

/** changes or adds a coefficient to an LP column */
SCIP_RETCODE SCIPcolExactChgCoef(
   SCIP_COLEXACT*        col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_RATIONAL*        val                 /**< value of coefficient */
   )
{
   int pos;

   assert(col != NULL);
   assert(lpexact != NULL);
   assert(!lpexact->fplp->diving);
   assert(row != NULL);

   /* search the position of the row in the column's row vector */
   pos = colExactSearchCoef(col, row);

   /* check, if row already exists in the column's row vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( colExactAddCoef(col, blkmem, set, eventqueue, lpexact, row, val, -1) );
   }
   else
   {
      /* modify already existing coefficient */
      assert(0 <= pos && pos < col->len);
      assert(col->rows[pos] == row);

      /* if row knows of the column, change the corresponding coefficient in the row */
      if( col->linkpos[pos] >= 0 )
      {
         assert(row->cols[col->linkpos[pos]] == col);
         assert(row->cols_index[col->linkpos[pos]] == col->index);
         assert(SCIPrationalIsEQ(row->vals[col->linkpos[pos]], col->vals[pos]));
         SCIP_CALL( rowExactChgCoefPos(row, set, lpexact, col->linkpos[pos], val) );
      }

      /* change the coefficient in the column */
      SCIP_CALL( colExactChgCoefPos(col, set, lpexact, pos, val) );
   }

   checkLinks(lpexact);

   return SCIP_OKAY;
}

/** increases value of an existing or nonexisting coefficient in an LP column */
SCIP_RETCODE SCIPcolExactIncCoef(
   SCIP_COLEXACT*        col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_RATIONAL*        incval              /**< value to add to the coefficient */
   )
{
   int pos;

   assert(col != NULL);
   assert(lpexact != NULL);
   assert(!lpexact->fplp->diving);
   assert(row != NULL);

   if( SCIPrationalIsZero(incval) )
      return SCIP_OKAY;

   /* search the position of the row in the column's row vector */
   pos = colExactSearchCoef(col, row);

   /* check, if row already exists in the column's row vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( colExactAddCoef(col, blkmem, set, eventqueue, lpexact, row, incval, -1) );
   }
   else
   {
      /* modify already existing coefficient */
      assert(0 <= pos && pos < col->len);
      assert(col->rows[pos] == row);

      /* if row knows of the column, change the corresponding coefficient in the row */
      if( col->linkpos[pos] >= 0 )
      {
         assert(row->cols[col->linkpos[pos]] == col);
         assert(row->cols_index[col->linkpos[pos]] == col->index);
         assert(SCIPrationalIsEQ(row->vals[col->linkpos[pos]], col->vals[pos]));

         SCIPrationalAdd(incval, incval, col->vals[pos]);
         SCIP_CALL( rowExactChgCoefPos(row, set, lpexact, col->linkpos[pos], incval) );
      }

      /* change the coefficient in the column */
      SCIP_CALL( colExactChgCoefPos(col, set, lpexact, pos, incval) );
   }

   checkLinks(lpexact);

   return SCIP_OKAY;
}

/** changes objective value of column */
SCIP_RETCODE SCIPcolExactChgObj(
   SCIP_COLEXACT*        col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_RATIONAL*        newobj              /**< new objective value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatusExact(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetColExact(col->var) == col);
   assert(lpexact != NULL);

   SCIPrationalDebugMessage("changing objective value of column <%s> from %q to %q\n", SCIPvarGetName(col->var), col->obj, newobj);

   /* only add actual changes */
   if( !SCIPrationalIsEQ(col->obj, newobj) )
   {
      /* only variables with a real position in the LPI can be inserted */
      if( col->lpipos >= 0 )
      {
         /* insert column in the chgcols list (if not already there) */
         SCIP_CALL( insertColChgcols(col, set, lpexact) );

         /* mark objective value change in the column */
         col->objchanged = TRUE;

         assert(lpexact->nchgcols > 0);
      }
      /* in any case, when the sign of the objective (and thereby the best bound) changes, the variable has to enter the
       * LP and the LP has to be flushed
       */
      else if( (SCIPrationalIsNegative(col->obj) && SCIPrationalIsPositive(newobj) && SCIPrationalIsZero(col->ub))
         || (SCIPrationalIsPositive(col->obj) && SCIPrationalIsNegative(newobj) && SCIPrationalIsZero(col->lb)) )
      {
         /* mark the LP unflushed */
         lpexact->flushed = FALSE;
      }
   }

   /* store new objective function value */
   SCIPrationalSetRational(col->obj, newobj);

   return SCIP_OKAY;
}

/** changes lower bound of column */
SCIP_RETCODE SCIPcolExactChgLb(
   SCIP_COLEXACT*        col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_RATIONAL*        newlb               /**< new lower bound value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatusExact(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetColExact(col->var) == col);
   assert(lpexact != NULL);

   SCIPrationalDebugMessage("changing lower bound of column <%s> from %q to %q\n", SCIPvarGetName(col->var), col->lb, newlb);

   /* only add actual changes */
   if( !SCIPrationalIsEQ(col->lb, newlb) )
   {
      /* only variables with a real position in the LPI can be inserted */
      if( col->lpipos >= 0 )
      {
         /* insert column in the chgcols list (if not already there) */
         SCIP_CALL( insertColChgcols(col, set, lpexact) );

         /* mark bound change in the column */
         col->lbchanged = TRUE;

         assert(lpexact->nchgcols > 0);
      }
      /* in any case, when the best bound is zero and gets changed, the variable has to enter the LP and the LP has to be
       * flushed
       */
      else if( !SCIPrationalIsNegative(col->obj) && SCIPrationalIsZero(col->lb) )
      {
         /* mark the LP unflushed */
         lpexact->flushed = FALSE;
      }
   }

   SCIPrationalSetRational(col->lb, newlb);

   return SCIP_OKAY;
}

/** changes upper bound of exact column */
SCIP_RETCODE SCIPcolExactChgUb(
   SCIP_COLEXACT*        col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_RATIONAL*        newub               /**< new upper bound value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatusExact(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetColExact(col->var) == col);
   assert(lpexact != NULL);

   SCIPrationalDebugMessage("changing upper bound of column <%s> from %q to %q\n", SCIPvarGetName(col->var), col->ub, newub);

   /* only add actual changes */
   if( !SCIPrationalIsEQ(col->ub, newub) )
   {
      /* only variables with a real position in the LPI can be inserted */
      if( col->lpipos >= 0 )
      {
         /* insert column in the chgcols list (if not already there) */
         SCIP_CALL( insertColChgcols(col, set, lpexact) );

         /* mark bound change in the column */
         col->ubchanged = TRUE;

         assert(lpexact->nchgcols > 0);
      }
      /* in any case, when the best bound is zero and gets changed, the variable has to enter the LP and the LP has to be
       * flushed
       */
      else if( SCIPrationalIsNegative(col->obj) && SCIPrationalIsZero(col->ub) )
      {
         /* mark the LP unflushed */
         lpexact->flushed = FALSE;
      }
   }

   SCIPrationalSetRational(col->ub, newub);

   return SCIP_OKAY;
}

/** creates and captures an LP row */
SCIP_RETCODE SCIProwExactCreate(
   SCIP_ROWEXACT**       row,                /**< pointer to LP row data */
   SCIP_ROW*             fprow,              /**< corresponding fp row */
   SCIP_ROW*             fprowrhs,           /**< rhs-part of fp-relaxation of this row if necessary, NULL otherwise */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COLEXACT**       cols,               /**< array with columns of row entries */
   SCIP_RATIONAL**       vals,               /**< array with coefficients of row entries */
   SCIP_RATIONAL*        lhs,                /**< left hand side of row */
   SCIP_RATIONAL*        rhs,                /**< right hand side of row */
   SCIP_Bool             isfprelaxable       /**< is it possible to make fp-relaxation of this row */
   ) /*lint --e{715}*/
{
   assert(row != NULL);
   assert(fprow != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(len >= 0);
   assert(len == 0 || (cols != NULL && vals != NULL));
   assert(SCIProwGetNNonz(fprow) == len || len == 0);
   /* note, that the assert tries to avoid numerical troubles in the LP solver.
    * in case, for example, lhs > rhs but they are equal with tolerances, one could pass lhs=rhs=lhs+rhs/2 to
    * SCIProwCreate() (see cons_linear.c: detectRedundantConstraints())
    */
   assert(SCIPrationalIsLE(lhs, rhs));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, row) );

   (*row)->storedsolvals = NULL;
   (*row)->integral = TRUE;
   (*row)->fprow = fprow;
   (*row)->fprowrhs = fprowrhs;
   (*row)->nuses = 0;
   (*row)->nlocks = 0;
   fprow->rowexact = (*row);
   SCIProwExactCapture(*row);
   if( fprowrhs != NULL )
   {
      SCIProwExactCapture(*row);
      fprowrhs->rowexact = (*row);
   }

   if( len > 0 )
   {
      SCIP_VAR* var;
      int i;

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*row)->cols, cols, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*row)->cols_index, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*row)->linkpos, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*row)->valsinterval, len) );
      SCIP_CALL( SCIPrationalCopyBlockArray(blkmem, &(*row)->vals, vals, len) );

      for( i = 0; i < len; ++i )
      {
         assert(cols[i] != NULL);
         assert(!SCIPrationalIsZero(vals[i]));

         var = cols[i]->var;
         (*row)->cols_index[i] = cols[i]->index;
         (*row)->linkpos[i] = -1;
         SCIPintervalSetRational(&(*row)->valsinterval[i], vals[i]);

         if( SCIPrationalIsIntegral((*row)->vals[i]) )
            (*row)->integral = (*row)->integral && SCIPvarIsIntegral(var);
         else
         {
            (*row)->integral = FALSE;
         }
      }
   }
   else
   {
      (*row)->cols = NULL;
      (*row)->vals = NULL;
      (*row)->valsinterval = NULL;
      (*row)->linkpos = NULL;
      (*row)->cols_index = NULL;
   }

   SCIP_CALL( SCIPrationalCopyBlock(blkmem, &(*row)->lhs, lhs) );
   SCIP_CALL( SCIPrationalCopyBlock(blkmem, &(*row)->rhs, rhs) );
   SCIP_CALL( SCIPrationalCreateString(blkmem, &(*row)->flushedlhs, "-inf") );
   SCIP_CALL( SCIPrationalCreateString(blkmem, &(*row)->flushedrhs, "inf") );
   SCIP_CALL( SCIPrationalCreateString(blkmem, &(*row)->objprod, "0") );
   SCIP_CALL( SCIPrationalCreateString(blkmem, &(*row)->dualsol, "0") );
   SCIP_CALL( SCIPrationalCreateString(blkmem, &(*row)->activity, "inf") );
   SCIP_CALL( SCIPrationalCreateString(blkmem, &(*row)->dualfarkas, "0") );
   SCIP_CALL( SCIPrationalCreateString(blkmem, &(*row)->pseudoactivity, "inf") );
   SCIP_CALL( SCIPrationalCreateString(blkmem, &(*row)->constant, "0") );

   (*row)->index = stat->nrowidx;
   SCIPstatIncrement(stat, set, nrowidx);
   (*row)->size = len;
   (*row)->len = len;
   (*row)->nlpcols = 0;
   (*row)->nunlinked = len;
   (*row)->lppos = -1;
   (*row)->lpipos = -1;
   (*row)->lpdepth = -1;
   (*row)->validactivitylp = -1;
   (*row)->delaysort = FALSE;
   (*row)->lpcolssorted = TRUE;
   (*row)->nonlpcolssorted = (len <= 1);
   (*row)->delaysort = FALSE;
   (*row)->fprelaxable = isfprelaxable;
   (*row)->rhsreal = SCIPrationalRoundReal((*row)->rhs, SCIP_R_ROUND_UPWARDS);
   (*row)->lhsreal = SCIPrationalRoundReal((*row)->lhs, SCIP_R_ROUND_DOWNWARDS);
   SCIPintervalSet(&(*row)->constantreal, 0.0);
   return SCIP_OKAY;
} /*lint !e715*/

/** changes an exact row, so that all denominators are bounded by set->exact_cutmaxdenom */
static
SCIP_RETCODE rowExactCreateFromRowLimitEncodingLength(
   SCIP_ROW*             row,                /**< SCIP row */
   SCIP_ROWEXACT*        rowexact,           /**< exact row */
   SCIP_SET*             set,                /**< SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory structure */
   SCIP_EVENTQUEUE*      eventqueue,         /**< the eventqueue */
   SCIP_LPEXACT*         lpexact             /**< the exact lp */
   )
{
   int i;
   SCIP_Longint maxdenom;
   SCIP_Longint maxboundval;
   SCIP_RATIONAL* val;
   SCIP_RATIONAL* newval;
   SCIP_RATIONAL* difference;
   SCIP_Real rhschange;
   int forcegreater;

   assert(row != NULL);
   assert(set->exact_cutmaxdenom > 0);
   assert(set->exact_cutapproxmaxboundval >= 0);

   SCIPdebugMessage("approximating row ");
   SCIPdebug(SCIPprintRow(set->scip, row, NULL));

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &val) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &difference) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &newval) );

   rhschange = 0;
   maxboundval = set->exact_cutapproxmaxboundval;
   maxdenom = set->exact_cutmaxdenom;

   assert(maxdenom >= 0);
   assert(maxboundval >= 0);

   for( i = 0; i <= row->len - 1; ++i )
   {
      SCIP_VAR* var = row->cols[i]->var;
      SCIPrationalSetReal(val, row->vals[i]);

      forcegreater = 0;

      if( SCIPrationalIsNegInfinity(SCIPvarGetLbGlobalExact(var)) && SCIPrationalIsInfinity(SCIPvarGetUbGlobalExact(var)) )
         forcegreater = -2;
      else if( SCIPrationalIsNegInfinity(SCIPvarGetLbGlobalExact(var)) )
         forcegreater = 1;
      else if( SCIPrationalIsInfinity(SCIPvarGetUbGlobalExact(var)) )
         forcegreater = -1;

      if( forcegreater == -2 || SCIPrationalDenominatorIsLE(val, maxdenom) ||
            ((maxboundval > 0) && SCIPrationalIsGTReal(SCIPvarGetUbGlobalExact(var), (double) maxboundval)) ||
            SCIPrationalIsLTReal(SCIPvarGetLbGlobalExact(var), (double) -maxboundval) )
      {
         SCIPrationalSetRational(newval, val);
      }
      else
         SCIPrationalComputeApproximation(newval, val, maxdenom, forcegreater);
#ifndef NDEBUG
      if( forcegreater == 1 )
         assert(SCIPrationalIsGE(newval, val));
      else if( forcegreater == -1 )
         assert(SCIPrationalIsLE(newval, val));
#endif

      SCIPrationalDiff(difference, newval, val);
      if( SCIPrationalIsPositive(difference) )
         SCIPrationalAddProd(rowexact->rhs, difference, SCIPvarGetUbGlobalExact(var));
      else
         SCIPrationalAddProd(rowexact->rhs, difference, SCIPvarGetLbGlobalExact(var));

      if( !SCIPrationalIsZero(newval) )
      {
         SCIP_CALL( SCIProwExactAddCoef(rowexact, blkmem, set, eventqueue, lpexact, SCIPcolGetColExact(row->cols[i]), newval) );
      }

      if( SCIPrationalIsNegative(SCIPvarGetLbGlobalExact(var)) && !SCIPrationalIsZero(newval) )
      {
         rhschange += (SCIPintervalGetInf(rowexact->valsinterval[rowexact->len - 1]) - SCIPintervalGetSup(rowexact->valsinterval[rowexact->len - 1])) * SCIPvarGetLbGlobal(var);
      }
   }

   SCIPrationalComputeApproximation(newval, rowexact->rhs, maxdenom, 1);
   assert(SCIPrationalIsGE(newval, rowexact->rhs));

   SCIPrationalSetRational(rowexact->rhs, newval);

   SCIProwExactSort(rowexact);

   for( i = rowexact->fprow-> len-1; i >= 0; i-- )
   {
      SCIP_CALL( SCIProwDelCoef(rowexact->fprow, blkmem, set, eventqueue, lpexact->fplp, rowexact->fprow->cols[i]) );
   }

   for( i = 0; i < rowexact->len; i++ )
   {
      SCIP_CALL( SCIProwAddCoef(rowexact->fprow, blkmem, set, eventqueue, lpexact->fplp, rowexact->cols[i]->fpcol,
         SCIPrationalRoundReal(rowexact->vals[i], SCIP_R_ROUND_DOWNWARDS)) );
   }

   SCIP_CALL( SCIProwChgRhs(rowexact->fprow, blkmem, set, eventqueue, lpexact->fplp,
      SCIPrationalRoundReal(rowexact->rhs, SCIP_R_ROUND_UPWARDS) + rhschange) );

   SCIPrationalFreeBuffer(set->buffer, &newval);
   SCIPrationalFreeBuffer(set->buffer, &difference);
   SCIPrationalFreeBuffer(set->buffer, &val);

   SCIPdebugMessage("new row ");
   SCIPdebug(SCIPprintRowExact(set->scip, rowexact, NULL));

   return SCIP_OKAY;
}

/** creates and captures an exact LP row from a fp row
 *
 *  @note This may change the floating-point coefficients slightly if the rational representation is rounded to smaller
 *  denominators according to parameter exact/cutmaxdenom.
 */
SCIP_RETCODE SCIProwExactCreateFromRow(
   SCIP_ROW*             fprow,              /**< corresponding fp row to create from */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< the eventqueue */
   SCIP_PROB*            prob,               /**< scip prob structure */
   SCIP_LPEXACT*         lpexact             /**< current LP data */
   )
{
   SCIP_ROWEXACT** row;
   SCIP_ROWEXACT* workrow;
   int i;
   int oldnlocks;
   SCIP_RATIONAL* tmpval;
   SCIP_RATIONAL* tmplhs;
   SCIP_Real* rowvals;

   row = &(fprow->rowexact);

   /* unlock the row temporarily to be able to change it (slightly) */
   oldnlocks = (int) fprow->nlocks;
   fprow->nlocks = 0;

   assert(row != NULL);
   assert(fprow != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmpval) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmplhs) );

   if( !SCIPsetIsInfinity(set, fprow->rhs) )
      SCIPrationalSetReal(tmpval, fprow->rhs);
   else
      SCIPrationalSetInfinity(tmpval);

   if( !SCIPsetIsInfinity(set, -fprow->lhs) )
      SCIPrationalSetReal(tmplhs, fprow->lhs);
   else
      SCIPrationalSetNegInfinity(tmplhs);

   SCIP_CALL( SCIProwExactCreate(row, fprow, NULL, blkmem, set, stat, lpexact, 0, NULL, NULL, tmplhs, tmpval, TRUE) );

   workrow = *row;
   rowvals = SCIProwGetVals(fprow);
   workrow->removable = TRUE;

   SCIP_CALL( SCIProwExactEnsureSize(workrow, blkmem, set, fprow->size) );

   SCIPrationalSetReal(tmpval, SCIProwGetConstant(fprow));
   SCIP_CALL( SCIProwExactAddConstant(workrow, set, stat, lpexact, tmpval) );

   if( set->exact_cutmaxdenom > 0 )
   {
      SCIP_CALL( rowExactCreateFromRowLimitEncodingLength(fprow, workrow, set, blkmem, eventqueue, lpexact) );
      SCIProwRecalcNorms(fprow, set);
   }
   else
   {
      for( i = 0; i < SCIProwGetNNonz(fprow); i++ )
      {
         SCIP_COL* col;

         SCIPrationalSetReal(tmpval, rowvals[i]);
         col = SCIProwGetCols(fprow)[i];

         SCIP_CALL( SCIPvarAddToRowExact(SCIPcolGetVar(col), blkmem, set, stat, eventqueue, prob, lpexact, workrow, tmpval) );
         assert(SCIPrationalIsFpRepresentable(SCIProwExactGetVals(workrow)[SCIProwExactGetNNonz(workrow) -1]));
      }
   }

   SCIPrationalFreeBuffer(set->buffer, &tmplhs);
   SCIPrationalFreeBuffer(set->buffer, &tmpval);

   fprow->nlocks = oldnlocks; /*lint !e732*/

   return SCIP_OKAY;
}

/** populate data of two empty fp rows with data from exact row */
SCIP_RETCODE SCIProwExactGenerateFpRows(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< SCIP statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_ROWEXACT*        row,                /**< SCIP row */
   SCIP_ROW*             rowlhs,             /**< fp row-relaxation wrt lhs */
   SCIP_ROW*             rowrhs,             /**< fp row-relaxation wrt rhs */
   SCIP_Bool*            onerowrelax,        /**< is one row enough to represent the exact row */
   SCIP_Bool*            hasfprelax          /**< is it possible to generate relaxations at all for this row? */
   )
{
   SCIP_Real* valsrhsrelax;
   SCIP_Real* valslhsrelax;
   SCIP_Real rhsrelax;
   SCIP_Real lhsrelax;
   SCIP_VAR* var;
   SCIP_RATIONAL* ub;
   SCIP_RATIONAL* lb;
   SCIP_INTERVAL* rowexactvalsinterval;
   SCIP_Real lbreal = 0.0;
   SCIP_Real ubreal = 0.0;
   SCIP_ROUNDMODE roundmode;
   int i;
   int* sideindexpostprocess;
   int npostprocess;

   assert(row != NULL);
   assert(rowlhs != NULL);
   assert(rowrhs != NULL);

   SCIP_CALL( SCIPsetAllocBufferArray(set, &valsrhsrelax, row->len) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &valslhsrelax, row->len) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &sideindexpostprocess, row->len) );

   npostprocess = 0;
   *hasfprelax = TRUE;
   *onerowrelax = TRUE;

   rhsrelax = SCIPrationalRoundReal(row->rhs, SCIP_R_ROUND_UPWARDS);
   lhsrelax = SCIPrationalRoundReal(row->lhs, SCIP_R_ROUND_DOWNWARDS);
   roundmode = SCIPintervalGetRoundingMode();
   rowexactvalsinterval = row->valsinterval;

   /* we need to adapt the two fprows that are created, so that one is a relaxation wrt the lhs and the other wrt the rhs of the row */
   for( i = 0; i < row->len; i++ )
   {
      var = SCIPcolExactGetVar(row->cols[i]);
      ub = SCIPvarGetUbGlobalExact(var);
      lb = SCIPvarGetLbGlobalExact(var);

      /* coefficient is exactly representable as fp number */
      if( rowexactvalsinterval[i].inf == rowexactvalsinterval[i].sup )/*lint !e777*/
      {
         valslhsrelax[i] = rowexactvalsinterval[i].inf;
         valsrhsrelax[i] = rowexactvalsinterval[i].inf;
      }
      /* unbounded variable with non fp-representable coefficient: var would need to be split in pos/neg to be relaxable */
      else if( SCIPrationalIsInfinity(ub) && SCIPrationalIsNegInfinity(lb) )
      {
         *hasfprelax = FALSE;
         valslhsrelax[i] = rowexactvalsinterval[i].inf;
      }
      /* negative upper or positive lower bounds are good */
      else if( !SCIPrationalIsInfinity(ub) && SCIPrationalIsNegative(ub) )
      {
         *onerowrelax = FALSE;
         valslhsrelax[i] = rowexactvalsinterval[i].inf;
         valsrhsrelax[i] = rowexactvalsinterval[i].sup;
      }
      /* negative upper or positive lower bounds are good */
      else if( !SCIPrationalIsNegInfinity(lb) && SCIPrationalIsPositive(lb) )
      {
         *onerowrelax = FALSE;
         valslhsrelax[i] = rowexactvalsinterval[i].sup;
         valsrhsrelax[i] = rowexactvalsinterval[i].inf;
      }
      else if( !SCIPrationalIsInfinity(ub) )
      {
         *onerowrelax = FALSE;
         valslhsrelax[i] = rowexactvalsinterval[i].inf;
         valsrhsrelax[i] = rowexactvalsinterval[i].sup;
         sideindexpostprocess[npostprocess] = i;
         npostprocess++;
      }
      else
      {
         assert(!SCIPrationalIsInfinity(lb));
         *onerowrelax = FALSE;
         valslhsrelax[i] = rowexactvalsinterval[i].sup;
         valsrhsrelax[i] = rowexactvalsinterval[i].inf;
         sideindexpostprocess[npostprocess] = i;
         npostprocess++;
      }
   }

   SCIPintervalSetRoundingModeUpwards();

   /* change the sides where necessary (do not do it immediately to not change rounding mode too often) */
   for( i = 0; i < npostprocess; i++ )
   {
      int idx;
      idx = sideindexpostprocess[i];
      var = SCIPcolExactGetVar(row->cols[idx]);
      lbreal = SCIPvarGetLbGlobal(var);
      ubreal = SCIPvarGetUbGlobal(var);

      if( valslhsrelax[idx] == rowexactvalsinterval[idx].inf )/*lint !e777*/
         rhsrelax += ubreal >= 0 ? (rowexactvalsinterval[idx].sup - rowexactvalsinterval[idx].inf) * ubreal : 0;
      else
         rhsrelax -= lbreal <= 0 ? (rowexactvalsinterval[idx].sup - rowexactvalsinterval[idx].inf) * lbreal : 0;
   }

   SCIPintervalSetRoundingModeDownwards();
   for( i = 0; i < npostprocess; i++ )
   {
      int idx;
      idx = sideindexpostprocess[i];
      var = SCIPcolExactGetVar(row->cols[idx]);
      lbreal = SCIPvarGetLbGlobal(var);
      ubreal = SCIPvarGetUbGlobal(var);

      //  upper bound was used
      if( valslhsrelax[idx] == rowexactvalsinterval[idx].sup )/*lint !e777*/
         lhsrelax -= ubreal >= 0 ? (rowexactvalsinterval[i].sup - rowexactvalsinterval[i].inf) * ubreal : 0;
      else
         lhsrelax += lbreal <= 0 ? (rowexactvalsinterval[i].sup - rowexactvalsinterval[i].inf) * lbreal : 0;
   }

   SCIPintervalSetRoundingMode(roundmode);

   /* only create one row if possible, or if relaxation did not work at all */
   if( !(*hasfprelax) || *onerowrelax )
   {
      if( !SCIPsetIsInfinity(set, rhsrelax) )
      {
         SCIP_CALL( SCIProwChgRhs(rowlhs, blkmem, set, eventqueue, lpexact->fplp, rhsrelax) );
      }
      if( !SCIPsetIsInfinity(set, -lhsrelax) )
      {
         SCIP_CALL( SCIProwChgLhs(rowlhs, blkmem, set, eventqueue, lpexact->fplp, lhsrelax) );
      }

      for( i = 0; i < row->len; i++ )
      {
         SCIP_CALL( SCIPvarAddToRow(row->cols[i]->var, blkmem, set, stat,
            eventqueue, prob, lpexact->fplp, rowlhs, valslhsrelax[i]) );
      }

      /* we created the fprows directly from the exact row, so we should only have active variables inside it */
      assert(SCIProwGetConstant(rowlhs) == 0.0);
      SCIP_CALL( SCIProwChgConstant(rowlhs, blkmem, set, stat, eventqueue, lpexact->fplp, SCIPrationalRoundReal(row->constant, SCIP_R_ROUND_DOWNWARDS)) );

      SCIP_CALL( SCIProwRelease(&rowrhs, blkmem, set, lpexact->fplp) );
   }
   /* create two fp-rows for row relaxation */
   else
   {
      if( !SCIPsetIsInfinity(set, rhsrelax) )
      {
         SCIP_CALL( SCIProwChgRhs(rowlhs, blkmem, set, eventqueue, lpexact->fplp, rhsrelax) );
         SCIP_CALL( SCIProwChgRhs(rowrhs, blkmem, set, eventqueue, lpexact->fplp, rhsrelax) );
      }
      if( !SCIPsetIsInfinity(set, -lhsrelax) )
      {
         SCIP_CALL( SCIProwChgLhs(rowlhs, blkmem, set, eventqueue, lpexact->fplp, lhsrelax) );
         SCIP_CALL( SCIProwChgLhs(rowrhs, blkmem, set, eventqueue, lpexact->fplp, lhsrelax) );
      }

      for( i = 0; i < row->len; i++ )
      {
         SCIP_CALL( SCIPvarAddToRow(row->cols[i]->var, blkmem, set, stat,
            eventqueue, prob, lpexact->fplp, rowlhs, valslhsrelax[i]) );

         SCIP_CALL( SCIPvarAddToRow(row->cols[i]->var, blkmem, set, stat,
            eventqueue, prob, lpexact->fplp, rowrhs, valsrhsrelax[i]) );
      }

      /* we created the fprows directly from the exact row, so we should only have active variables inside it */
      assert(SCIProwGetConstant(rowlhs) == 0.0);
      assert(SCIProwGetConstant(rowrhs) == 0.0);
      SCIP_CALL( SCIProwChgConstant(rowlhs, blkmem, set, stat, eventqueue, lpexact->fplp, SCIPrationalRoundReal(row->constant, SCIP_R_ROUND_UPWARDS)) );
      SCIP_CALL( SCIProwChgConstant(rowrhs, blkmem, set, stat, eventqueue, lpexact->fplp, SCIPrationalRoundReal(row->constant, SCIP_R_ROUND_DOWNWARDS)) );
   }

   row->fprelaxable = *hasfprelax;

   SCIPsetFreeBufferArray(set, &sideindexpostprocess);
   SCIPsetFreeBufferArray(set, &valslhsrelax);
   SCIPsetFreeBufferArray(set, &valsrhsrelax);

   return SCIP_OKAY;
}

/** applies all cached changes to the LP solver */
SCIP_RETCODE SCIPlpExactFlush(
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   assert(lpexact != NULL);
   assert(blkmem != NULL);

   SCIPsetDebugMsg(set, "flushing exact LP changes: old (%d cols, %d rows), nchgcols=%d, nchgrows=%d, firstchgcol=%d, firstchgrow=%d, new (%d cols, %d rows), flushed=%u\n",
      lpexact->nlpicols, lpexact->nlpirows, lpexact->nchgcols, lpexact->nchgrows, lpexact->lpifirstchgcol, lpexact->lpifirstchgrow, lpexact->ncols, lpexact->nrows, lpexact->flushed);

   if( !lpexact->flushed )
   {
      lpexact->flushdeletedcols = FALSE;
      lpexact->flushaddedcols = FALSE;
      lpexact->flushdeletedrows = FALSE;
      lpexact->flushaddedrows = FALSE;

      SCIP_CALL( lpExactFlushDelCols(lpexact) );
      SCIP_CALL( lpExactFlushDelRows(lpexact, blkmem, set) );
      SCIP_CALL( lpExactFlushChgCols(lpexact, set) );
      SCIP_CALL( lpExactFlushChgRows(lpexact, set) );
      SCIP_CALL( lpExactFlushAddCols(lpexact, blkmem, set, eventqueue) );
      SCIP_CALL( lpExactFlushAddRows(lpexact, blkmem, set, eventqueue) );

      lpexact->flushed = TRUE;

      checkLinks(lpexact);
   }

   /* we can't retrieve the solution from Qsoptex after anything was changed, so we need to resolve the lp */
#ifdef SCIP_WITH_QSOPTEX
   lpexact->solved = FALSE;
   lpexact->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
#endif

   assert(lpexact->nlpicols == lpexact->ncols);
   assert(lpexact->lpifirstchgcol == lpexact->nlpicols);
   assert(lpexact->nlpirows == lpexact->nrows);
   assert(lpexact->lpifirstchgrow == lpexact->nlpirows);
   assert(lpexact->nchgcols == 0);
   assert(lpexact->nchgrows == 0);
#ifndef NDEBUG
   {
      int ncols;
      int nrows;

      SCIP_CALL( SCIPlpiExactGetNCols(lpexact->lpiexact, &ncols) );
      SCIP_CALL( SCIPlpiExactGetNRows(lpexact->lpiexact, &nrows) );
      assert(ncols == lpexact->ncols);
      assert(nrows == lpexact->nrows);
   }
#endif

   return SCIP_OKAY;
}

/** ensures all rows/columns are correctly updated, but changes are not yet communicated to the exact LP solver */
SCIP_RETCODE SCIPlpExactLink(
   SCIP_LPEXACT*         lpexact,                 /**< current exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   int c, r;
   SCIP_COLEXACT* col;
   SCIP_ROWEXACT* row;

   assert(lpexact != NULL);
   assert(blkmem != NULL);

   SCIPsetDebugMsg(set, "flushing exact LP changes: old (%d cols, %d rows), nchgcols=%d, nchgrows=%d, firstchgcol=%d, firstchgrow=%d, new (%d cols, %d rows), flushed=%u\n",
      lpexact->nlpicols, lpexact->nlpirows, lpexact->nchgcols, lpexact->nchgrows, lpexact->lpifirstchgcol, lpexact->lpifirstchgrow, lpexact->ncols, lpexact->nrows, lpexact->flushed);

   if( !lpexact->flushed )
   {
      lpexact->flushdeletedcols = FALSE;
      lpexact->flushaddedcols = FALSE;
      lpexact->flushdeletedrows = FALSE;
      lpexact->flushaddedrows = FALSE;

      /* we still flush added/deleted columns since this should only happen at the very start of the solve */
      SCIP_CALL( lpExactFlushDelCols(lpexact) );
      SCIP_CALL( lpExactFlushAddCols(lpexact, blkmem, set, eventqueue) );

      /* link new columns/rows */
      for( c = lpexact->nlpicols; c < lpexact->ncols; ++c )
      {
         col = lpexact->cols[c];
         assert(col != NULL);
         assert(col->var != NULL);
         assert(SCIPvarGetStatusExact(col->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetColExact(col->var) == col);
         assert(col->lppos == c);

         SCIPsetDebugMsg(set, "linking added column <%s>: ", SCIPvarGetName(col->var));
         SCIP_CALL( colExactLink(col, blkmem, set, eventqueue, lpexact) );
      }
      for( r = lpexact->nlpirows; r < lpexact->nrows; ++r )
      {
         row = lpexact->rows[r];
         assert(row != NULL);
         assert(row->lppos == r);

         SCIPsetDebugMsg(set, "linking added exact row <%s>: ", row->fprow->name);

         SCIP_CALL( rowExactLink(row, blkmem, set, eventqueue, lpexact) );
      }

      checkLinks(lpexact);
   }

   return SCIP_OKAY;
}

/*
 * lp methods
 */

/** creates the data needed for project and shift bounding method */
static
SCIP_RETCODE SCIPlpPsdataCreate(
   SCIP_LPEXACT*         lpexact,            /**< pointer to LP data object */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory buffers */
   )
{
   SCIP_PROJSHIFTDATA* projshiftdata;

   assert(lpexact != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &lpexact->projshiftdata) );

   projshiftdata = lpexact->projshiftdata;

   projshiftdata->lpiexact = NULL;
   projshiftdata->dvarmap = NULL;
   projshiftdata->ndvarmap = 0;
   projshiftdata->interiorpoint = NULL;
   projshiftdata->interiorray = NULL;
   projshiftdata->violation = NULL;
   projshiftdata->correction = NULL;
   projshiftdata->commonslack = NULL;
   projshiftdata->includedrows = NULL;
   projshiftdata->projshiftbasis = NULL;
#if defined SCIP_WITH_GMP && defined SCIP_WITH_EXACTSOLVE
   projshiftdata->rectfactor = (qsnum_factor_work*) NULL;
#endif

   projshiftdata->nextendedrows = 0;
   projshiftdata->projshiftbasisdim = 0;
   projshiftdata->violationsize = 0;

   projshiftdata->projshiftdatacon = FALSE;
   projshiftdata->projshiftdatafail = FALSE;
   projshiftdata->projshifthaspoint = FALSE;
   projshiftdata->projshifthasray = FALSE;
   projshiftdata->projshiftobjweight = FALSE;
   projshiftdata->scaleobj = FALSE;
   projshiftdata->projshiftuseintpoint = TRUE;

   return SCIP_OKAY;
}

/** frees the exact LPI in project-and-shift */
static
SCIP_RETCODE SCIPlpExactProjectShiftFreeLPIExact(
   SCIP_LPIEXACT**       lpiexact            /**< pointer to LPI object */
   )
{
   int nlpirows;
   int nlpicols;

   assert(lpiexact != NULL);
   assert(*lpiexact != NULL);

   SCIP_CALL( SCIPlpiExactGetNRows(*lpiexact, &nlpirows) );
   SCIP_CALL( SCIPlpiExactDelRows(*lpiexact, 0, nlpirows - 1) );

   SCIP_CALL( SCIPlpiExactGetNCols(*lpiexact, &nlpicols) );
   SCIP_CALL( SCIPlpiExactDelCols(*lpiexact, 0, nlpicols - 1) );

   SCIP_CALL( SCIPlpiExactClear(*lpiexact) );
   SCIP_CALL( SCIPlpiExactFree(lpiexact) );

   return SCIP_OKAY;
}

/** frees the data needed for project and shift bounding method */
static
SCIP_RETCODE SCIPlpExactProjectShiftFree(
   SCIP_LPEXACT*         lpexact,            /**< pointer to LP data object */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory buffers */
   )
{
   SCIP_PROJSHIFTDATA* projshiftdata;

   assert(lpexact != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   projshiftdata = lpexact->projshiftdata;

   if( projshiftdata->lpiexact != NULL )
   {
      SCIP_CALL( SCIPlpExactProjectShiftFreeLPIExact(&projshiftdata->lpiexact) );
   }
   assert(projshiftdata->lpiexact == NULL);

   BMSfreeBlockMemoryArrayNull(blkmem, &projshiftdata->dvarmap, projshiftdata->ndvarmap);

   if( projshiftdata->interiorpoint != NULL )
      SCIPrationalFreeBlockArray(blkmem, &projshiftdata->interiorpoint, projshiftdata->nextendedrows);
   if( projshiftdata->interiorray != NULL )
      SCIPrationalFreeBlockArray(blkmem, &projshiftdata->interiorray, projshiftdata->nextendedrows);
   if( projshiftdata->violation != NULL )
      SCIPrationalFreeBlockArray(blkmem, &projshiftdata->violation, projshiftdata->violationsize);
   if( projshiftdata->correction != NULL )
      SCIPrationalFreeBlockArray(blkmem, &projshiftdata->correction, projshiftdata->nextendedrows);
   if( projshiftdata->commonslack != NULL )
      SCIPrationalFreeBlock(blkmem, &projshiftdata->commonslack);

   BMSfreeBlockMemoryArrayNull(blkmem, &projshiftdata->includedrows, projshiftdata->nextendedrows);
   BMSfreeBlockMemoryArrayNull(blkmem, &projshiftdata->projshiftbasis, projshiftdata->nextendedrows);

#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_EXACTSOLVE)
   if( projshiftdata->rectfactor != NULL )
      RECTLUfreeFactorization(projshiftdata->rectfactor);
#endif

   assert(projshiftdata->interiorpoint == NULL);
   assert(projshiftdata->interiorray == NULL);
   assert(projshiftdata->includedrows == NULL);
   assert(projshiftdata->projshiftbasis == NULL);
   assert(projshiftdata->commonslack == NULL);

   BMSfreeBlockMemoryNull(blkmem, &lpexact->projshiftdata);

   return SCIP_OKAY;
}

/** returns whether the success rate of the Neumaier-Shcherbina safe bounding method is sufficiently high */
SCIP_Bool SCIPlpExactBoundShiftUseful(
   SCIP_LPEXACT*         lpexact              /**< pointer to LP data object */
   )
{
   assert(lpexact != NULL);

   return lpexact->boundshiftuseful;
}

/** returns whether it is possible to use project and shift bounding method */
SCIP_Bool SCIPlpExactProjectShiftPossible(
   SCIP_LPEXACT*         lpexact             /**< pointer to LP data object */
   )
{
   assert(lpexact != NULL);
   assert(lpexact->projshiftdata != NULL);

   return !(lpexact->projshiftdata->projshiftdatafail);
}

/** checks that lp and fplp are properly synced */
SCIP_Bool SCIPlpExactIsSynced(
   SCIP_LPEXACT*         lpexact,            /**< pointer to LP data object */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     msg                 /**< message handler for debug output */
   )
{
   assert(lpexact != NULL);
   assert(msg != NULL);

   return lpExactInSync(lpexact, set, msg);
}

/** creates empty LP data object */
SCIP_RETCODE SCIPlpExactCreate(
   SCIP_LPEXACT**        lpexact,            /**< pointer to LP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_LP*              fplp,               /**< the floating point LP */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name                /**< problem name */
   )
{
   assert(lpexact != NULL);
   assert(fplp != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(name != NULL);

   SCIP_ALLOC( BMSallocMemory(lpexact) );

   /* open LP Solver interface */
   SCIP_CALL( SCIPlpiExactCreate(&(*lpexact)->lpiexact, messagehdlr, name, SCIP_OBJSEN_MINIMIZE) );
   SCIP_CALL( SCIPlpPsdataCreate(*lpexact, set, blkmem) );

   (*lpexact)->fplp = fplp;
   fplp->lpexact = *lpexact;

   (*lpexact)->lpicols = NULL;
   (*lpexact)->lpirows = NULL;
   (*lpexact)->chgcols = NULL;
   (*lpexact)->chgrows = NULL;
   (*lpexact)->cols = NULL;
   (*lpexact)->rows = NULL;
   (*lpexact)->divechgsides = NULL;
   (*lpexact)->divechgsidetypes = NULL;
   (*lpexact)->divechgrows = NULL;
   (*lpexact)->divelpistate = NULL;
   (*lpexact)->storedsolvals = NULL;
   (*lpexact)->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   (*lpexact)->flushdeletedcols = FALSE;
   (*lpexact)->flushaddedcols = FALSE;
   (*lpexact)->flushdeletedrows = FALSE;
   (*lpexact)->flushaddedrows = FALSE;
   (*lpexact)->updateintegrality = TRUE;
   (*lpexact)->flushed = TRUE;
   (*lpexact)->solved = FALSE;
   (*lpexact)->primalfeasible = TRUE;
   (*lpexact)->primalchecked = TRUE;
   (*lpexact)->diving = FALSE;
   (*lpexact)->divelpwasprimfeas = TRUE;
   (*lpexact)->divelpwasprimchecked = TRUE;
   (*lpexact)->divelpwasdualfeas = TRUE;
   (*lpexact)->divelpwasdualchecked = TRUE;
   (*lpexact)->divingobjchg = FALSE;
   (*lpexact)->dualfeasible = TRUE;
   (*lpexact)->dualchecked = TRUE;
   (*lpexact)->solisbasic = FALSE;
   (*lpexact)->resolvelperror = FALSE;
   (*lpexact)->projshiftpossible = FALSE;
   (*lpexact)->boundshiftuseful = TRUE;
   (*lpexact)->forceexactsolve = FALSE;
   (*lpexact)->allowexactsolve = FALSE;
   (*lpexact)->forcesafebound = FALSE;
   (*lpexact)->wasforcedsafebound = FALSE;
   (*lpexact)->lpiscaling = set->lp_scaling;
   (*lpexact)->lpisolutionpolishing = (set->lp_solutionpolishing > 0);
   (*lpexact)->lpirefactorinterval = set->lp_refactorinterval;
   (*lpexact)->lpiitlim = INT_MAX;
   (*lpexact)->lpipricing = SCIP_PRICING_AUTO;
   (*lpexact)->lastlpalgo = SCIP_LPALGO_DUALSIMPLEX;
   (*lpexact)->lpitiming = (int) set->time_clocktype;
   (*lpexact)->lpirandomseed = set->random_randomseed;

   (*lpexact)->lpicolssize = 0;
   (*lpexact)->nlpicols = 0;
   (*lpexact)->lpirowssize = 0;
   (*lpexact)->nlpirows = 0;
   (*lpexact)->lpifirstchgcol = 0;
   (*lpexact)->lpifirstchgrow = 0;
   (*lpexact)->colssize = 0;
   (*lpexact)->ncols = 0;
   (*lpexact)->nloosevars = 0;
   (*lpexact)->rowssize = 0;
   (*lpexact)->nrows = 0;
   (*lpexact)->chgcolssize = 0;
   (*lpexact)->nchgcols = 0;
   (*lpexact)->chgrowssize = 0;
   (*lpexact)->nchgrows = 0;
   (*lpexact)->firstnewcol = 0;
   (*lpexact)->firstnewrow = 0;
   (*lpexact)->looseobjvalinf = 0;
   (*lpexact)->pseudoobjvalinf = 0;
   (*lpexact)->glbpseudoobjvalinf = 0;
   (*lpexact)->ndivingrows = 0;
   (*lpexact)->ndivechgsides = 0;
   (*lpexact)->nremovablerows = 0;
   (*lpexact)->lpiobjlim = SCIPlpiExactInfinity((*lpexact)->lpiexact);
   (*lpexact)->cutoffbound = SCIPsetInfinity(set);
   (*lpexact)->oldcutoffbound = SCIPsetInfinity(set);
   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &(*lpexact)->lpobjval) );
   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &(*lpexact)->pseudoobjval) );
   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &(*lpexact)->glbpseudoobjval) );
   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &(*lpexact)->looseobjval) );

   return SCIP_OKAY;
}

/** frees LP data object */
SCIP_RETCODE SCIPlpExactFree(
   SCIP_LPEXACT**        lpexact,            /**< pointer to exact LP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   if( !set->exact_enable )
      return SCIP_OKAY;

   assert(lpexact != NULL);
   assert(*lpexact != NULL);

   SCIP_CALL( SCIPlpExactProjectShiftFree(*lpexact, set, blkmem) );
   SCIP_CALL( SCIPlpExactClear(*lpexact, blkmem, set) );

   //freeDiveChgSideArrays(*lpexact);

   /* release LPI rows */
   for( i = 0; i < (*lpexact)->nlpirows; ++i )
   {
      SCIP_CALL( SCIProwExactRelease(&(*lpexact)->lpirows[i], blkmem, set, *lpexact) );
   }

   if( (*lpexact)->lpiexact != NULL )
   {
      SCIP_CALL( SCIPlpiExactFree(&(*lpexact)->lpiexact) );
   }

   SCIPrationalFreeBlock(blkmem, &(*lpexact)->lpobjval);
   SCIPrationalFreeBlock(blkmem, &(*lpexact)->pseudoobjval);
   SCIPrationalFreeBlock(blkmem, &(*lpexact)->glbpseudoobjval);
   SCIPrationalFreeBlock(blkmem, &(*lpexact)->looseobjval);

   if( (*lpexact)->storedsolvals != NULL )
   {
      SCIPrationalFreeBlock(blkmem, &(*lpexact)->storedsolvals->lpobjval);
      BMSfreeMemoryNull(&(*lpexact)->storedsolvals);
   }
   BMSfreeMemoryArrayNull(&(*lpexact)->lpicols);
   BMSfreeMemoryArrayNull(&(*lpexact)->lpirows);
   BMSfreeMemoryArrayNull(&(*lpexact)->chgcols);
   BMSfreeMemoryArrayNull(&(*lpexact)->chgrows);
   BMSfreeMemoryArrayNull(&(*lpexact)->cols);
   BMSfreeMemoryArrayNull(&(*lpexact)->rows);
   BMSfreeMemory(lpexact);

   return SCIP_OKAY;
}

/** adds a column to the LP and captures the variable */
SCIP_RETCODE SCIPlpExactAddCol(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COLEXACT*        col                 /**< LP column */
   )
{
   if( !set->exact_enable )
      return SCIP_OKAY;

   assert(lpexact != NULL);
   assert(!lpexact->fplp->diving);
   assert(col != NULL);
   assert(col->len == 0 || col->rows != NULL);
   assert(col->lppos == -1);
   assert(col->var != NULL);
   assert(SCIPvarGetStatusExact(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col->fpcol);
   assert(SCIPvarIsIntegral(col->var) == col->fpcol->integral);

   SCIPsetDebugMsg(set, "adding column <%s> to exact LP (%d rows, %d cols)\n", SCIPvarGetName(col->var), lpexact->nrows, lpexact->ncols);
#ifdef SCIP_DEBUG
      SCIPrationalDebugMessage("(obj: %q) [%q,%q]", col->obj, col->lb, col->ub);
      for( int i = 0; i < col->len; ++i )
         SCIPrationalDebugMessage(" %q<%s>", col->vals[i], col->rows[i]->fprow->name);
      SCIPsetDebugMsgPrint(set, "\n");
#endif

   SCIP_CALL( ensureColexsSize(lpexact, set, lpexact->ncols+1) );
   lpexact->cols[lpexact->ncols] = col;
   col->lppos = lpexact->ncols;
   lpexact->ncols++;

   /* mark the current LP unflushed */
   lpexact->flushed = FALSE;

   /* update column arrays of all linked rows */
   SCIP_CALL( colExactUpdateAddLP(col, set) );

   checkLinks(lpexact);

   return SCIP_OKAY;
}

/** adds a row to the LP and captures it */
SCIP_RETCODE SCIPlpExactAddRow(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROWEXACT*        rowexact            /**< LP row */
   )
{
   assert(lpexact != NULL);
   assert(rowexact != NULL);
   assert(rowexact->len == 0 || rowexact->cols != NULL);
   assert(rowexact->lppos == -1);
   assert(rowexact->fprow != NULL);

   SCIProwExactCapture(rowexact);

   SCIPsetDebugMsg(set, "adding row <%s> to LP (%d rows, %d cols)\n", rowexact->fprow->name, lpexact->nrows, lpexact->ncols);
#ifdef SCIP_DEBUG
   {
      int i;
      SCIPrationalDebugMessage("  %q <=", rowexact->lhs);
      for( i = 0; i < rowexact->len; ++i )
         SCIPrationalDebugMessage(" %q<%s>", rowexact->vals[i], SCIPvarGetName(rowexact->cols[i]->var));
      if( !SCIPrationalIsZero(rowexact->constant) )
         SCIPrationalDebugMessage(" %q", rowexact->constant);
      SCIPrationalDebugMessage(" <= %q\n", rowexact->rhs);
   }
#endif

   SCIP_CALL( ensureRowexsSize(lpexact, set, lpexact->nrows+1) );
   lpexact->rows[lpexact->nrows] = rowexact;
   rowexact->lppos = lpexact->nrows;
   lpexact->nrows++;

   /* mark the current LP unflushed */
   lpexact->flushed = FALSE;

   /* update row arrays of all linked columns */
   SCIP_CALL( rowExactUpdateAddLP(rowexact, set) );

   return SCIP_OKAY;
}

/** should the objective limit of the LP solver be disabled */
#define lpCutoffDisabled(set) (set->lp_disablecutoff == 1 || (set->nactivepricers > 0 && set->lp_disablecutoff == 2))

/** sets the upper objective limit of the exact LP solver */
SCIP_RETCODE SCIPlpExactSetCutoffbound(
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             cutoffbound         /**< new upper objective limit */
   )
{
   SCIP_RATIONAL* tmpobj;

   if( !set->exact_enable )
      return SCIP_OKAY;

   assert(lpexact != NULL);

   SCIPsetDebugMsg(set, "setting exact LP upper objective limit from %g to %g\n", lpexact->cutoffbound, cutoffbound);

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmpobj) );
   if( lpexact->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && lpexact->solved && lpexact->flushed )
      SCIPlpExactGetObjval(lpexact, set, tmpobj);

   /* if the cutoff bound is increased, and the LP was proved to exceed the old cutoff, it is no longer solved */
   if( lpexact->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT && cutoffbound > lpexact->cutoffbound )
   {
      /* mark the current solution invalid */
      lpexact->solved = FALSE;
      SCIPrationalSetInfinity(lpexact->lpobjval);
      lpexact->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }
   /* if the cutoff bound is decreased below the current optimal value, the LP now exceeds the objective limit;
    * if the objective limit in the LP solver was disabled, the solution status of the LP is not changed
    */
   else if( !lpCutoffDisabled(set) && lpexact->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && lpexact->solved && lpexact->flushed
            && SCIPrationalIsGEReal(tmpobj, cutoffbound) )
   {
      assert(lpexact->flushed);
      assert(lpexact->solved);
      lpexact->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
   }
   SCIPrationalFreeBuffer(set->buffer, &tmpobj);
   lpexact->cutoffbound = cutoffbound;

   return SCIP_OKAY;
}

/** maximal number of verblevel-high messages about numerical trouble in LP that will be printed
 * when this number is reached and display/verblevel is not full, then further messages are suppressed in this run
 */
#define MAXNUMTROUBLELPMSGS 10

/** prints message about numerical trouble
 *
 * If message has verblevel at most high and display/verblevel is not full,
 * then the message is not printed if already MAXNUMTROUBLELPMSGS messages
 * were printed before in the current run.
 */
static
void lpExactNumericalTroubleMessage(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VERBLEVEL        verblevel,          /**< verbosity level of message */
   const char*           formatstr,          /**< message format string */
   ...                                       /**< arguments to format string */
   )
{
   va_list ap;

   assert(verblevel > SCIP_VERBLEVEL_NONE);
   assert(verblevel <= SCIP_VERBLEVEL_FULL);
   assert(set->disp_verblevel <= SCIP_VERBLEVEL_FULL);

   if( set->disp_verblevel < SCIP_VERBLEVEL_FULL )
   {
      if( verblevel <= SCIP_VERBLEVEL_HIGH )
      {
         /* if already max number of messages about numerical trouble in LP on verblevel at most high, then skip message */
         if( stat->nnumtroublelpmsgs > MAXNUMTROUBLELPMSGS )
            return;

         /* increase count on messages with verblevel high */
         ++stat->nnumtroublelpmsgs ;
      }

      /* if messages wouldn't be printed, then return already */
      if( verblevel > set->disp_verblevel )
         return;
   }

   /* print common begin of message */
   SCIPmessagePrintInfo(messagehdlr,
      "(node %" SCIP_LONGINT_FORMAT ") numerical troubles in exact LP %" SCIP_LONGINT_FORMAT " -- ",
      stat->nnodes, stat->nexlp);

   /* print individual part of message */
   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintInfo(messagehdlr, NULL, formatstr, ap);
   va_end(ap);

   /* warn that further messages will be suppressed */
   if( set->disp_verblevel < SCIP_VERBLEVEL_FULL && verblevel <= SCIP_VERBLEVEL_HIGH && stat->nnumtroublelpmsgs > MAXNUMTROUBLELPMSGS )
   {
      SCIPmessagePrintInfo(messagehdlr, " -- further messages will be suppressed (use display/verblevel=5 to see all)");
   }

   /* print closing new-line */
   SCIPmessagePrintInfo(messagehdlr, "\n");
}

/** flushes the exact LP and solves it with the primal or dual simplex algorithm, depending on the current basis feasibility */
static
SCIP_RETCODE lpExactFlushAndSolve(
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   harditlim,          /**< maximal number of LP iterations to perform (hard limit for all LP calls), or -1 for no limit */
   SCIP_Bool             fromscratch,        /**< should the LP be solved from scratch without using current basis? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   int* cstat;
   int* rstat;
   SCIP_Bool solveagain;
   SCIP_Bool success;
   SCIP_RETCODE retcode;
   SCIP_Real lptimelimit;
   char algo;
   SCIP_LP* lp;
   SCIP_LPISTATE* lpistate;

   assert(lpexact != NULL);
   assert(lpexact->fplp != NULL);
   assert(set != NULL);
   assert(lperror != NULL);
   assert(set->exact_enable);

   SCIP_CALL( SCIPlpiExactSetIntpar(lpexact->lpiexact, SCIP_LPPAR_LPINFO, (int) set->exact_lpinfo) );
   algo = set->lp_initalgorithm;
   lp = lpexact->fplp;

   /* set up the exact lpi for the current node */
   SCIP_CALL( SCIPlpExactSyncLPs(lpexact, blkmem, set) );
   SCIP_CALL( SCIPlpExactFlush(lpexact, blkmem, set, eventqueue) );

   assert(SCIPlpExactIsSynced(lpexact, set, messagehdlr));

   /* check if a time limit is set, and set time limit for LP solver accordingly */
   lptimelimit = SCIPlpiExactInfinity(lpexact->lpiexact);
   if( set->istimelimitfinite )
      lptimelimit = set->limit_time - SCIPclockGetTime(stat->solvingtime);

   success = FALSE;
   if( lptimelimit > 0.0 )
      SCIP_CALL( lpExactSetRealpar(lpexact, SCIP_LPPAR_LPTILIM, lptimelimit, &success) );

   if( lptimelimit <= 0.0 || !success )
   {
      SCIPsetDebugMsg(set, "time limit of %f seconds could not be set\n", lptimelimit);
      *lperror = ((lptimelimit > 0.0) ? TRUE : FALSE);
      SCIPsetDebugMsg(set, "time limit exceeded before solving LP\n");
      lp->solved = TRUE;
      lpexact->lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
      lp->lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
      lp->lpobjval = -SCIPsetInfinity(set);
      return SCIP_OKAY;
   }

   /* get lpi state to check whether basis exists */
   SCIP_CALL( SCIPlpiGetState(lp->lpi, blkmem, &lpistate) );

   /* set the correct basis information for warmstart */
   if( !fromscratch && SCIPlpiHasStateBasis(lp->lpi, lpistate) )
   {
      SCIP_CALL( SCIPsetAllocBufferArray(set, &cstat, lpexact->nlpicols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &rstat, lpexact->nlpirows) );

      SCIP_CALL( SCIPlpiGetBase(lp->lpi, cstat, rstat) );
      SCIP_CALL( SCIPlpiExactSetBase(lpexact->lpiexact, cstat, rstat) );

      SCIPsetFreeBufferArray(set, &cstat);
      SCIPsetFreeBufferArray(set, &rstat);
   }
   else
   {
      SCIP_CALL( SCIPlpiExactSetIntpar(lpexact->lpiexact, SCIP_LPPAR_FROMSCRATCH, TRUE) );
   }

   SCIP_CALL( SCIPlpiFreeState(lp->lpi, blkmem, &lpistate) );

   /* solve with given settings (usually fast but imprecise) */
   if( SCIPsetIsInfinity(set, lpexact->cutoffbound) )
   {
      SCIP_CALL( lpExactSetObjlim(lpexact, set, lpexact->cutoffbound, &success) );
   }
   else
   {
      SCIP_CALL( lpExactSetObjlim(lpexact, set, lpexact->cutoffbound - SCIPrationalRoundReal(getFiniteLooseObjvalExact(lpexact, set, prob), SCIP_R_ROUND_DOWNWARDS), &success) );
   }
   SCIP_CALL( lpExactSetIterationLimit(lpexact, harditlim) );

   do {
      solveagain = FALSE;

      /* solve the lp exactly */
      if( algo != 's' )
      {
         SCIPerrorMessage("Lp-algorithm-type %d is not supported in exact solving mode \n", algo);
         SCIPABORT();
      }

      SCIPsetDebugMsg(set, "Calling SCIPlpiExactSolveDual()\n");
      retcode = SCIPlpiExactSolveDual(lpexact->lpiexact);

      if( retcode == SCIP_LPERROR )
      {
         *lperror = TRUE;
         lpexact->solved = FALSE;
         lpexact->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         SCIPdebugMessage("Error solving lp exactly. \n");
      }

      if( retcode != SCIP_LPERROR )
      {
         SCIP_CALL( SCIPlpiExactGetSolFeasibility(lpexact->lpiexact, &(lpexact->primalfeasible), &(lpexact->dualfeasible)) );
         lpexact->solisbasic = TRUE;
      }

      /* only one should return true */
      assert(!(SCIPlpiExactIsOptimal(lpexact->lpiexact) && SCIPlpiExactIsObjlimExc(lpexact->lpiexact) && SCIPlpiExactIsPrimalInfeasible(lpexact->lpiexact) &&
            SCIPlpiExactExistsPrimalRay(lpexact->lpiexact) && SCIPlpiExactIsIterlimExc(lpexact->lpiexact) && SCIPlpiExactIsTimelimExc(lpexact->lpiexact)));

      /* evaluate solution status */
      if( SCIPlpiExactIsOptimal(lpexact->lpiexact) )
      {
         assert(lpexact->primalfeasible && lpexact->dualfeasible);

         SCIP_CALL( SCIPlpiExactGetObjval(lpexact->lpiexact, lpexact->lpobjval) );
         SCIPdebugMessage("Exact lp solve terminated with optimal. Safe dual bound is %e, previous lp obj-val was %e \n",
               SCIPrationalRoundReal(lpexact->lpobjval, SCIP_R_ROUND_DOWNWARDS), lp->lpobjval);
         lpexact->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
         lp->validsollp = stat->lpcount;

         if( !SCIPsetIsInfinity(set, lpexact->lpiobjlim) && SCIPrationalIsGTReal(lpexact->lpobjval, lpexact->lpiobjlim) )
         {
            /* the solver may return the optimal value, even if this is greater or equal than the upper bound */
            SCIPrationalDebugMessage("optimal solution %q exceeds objective limit %.15g\n", lpexact->lpobjval, lp->lpiobjlim);
            lpexact->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
            SCIPrationalSetInfinity(lpexact->lpobjval);
         }
      }
      else if( SCIPlpiExactIsObjlimExc(lpexact->lpiexact) )
      {
         lpexact->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
         SCIPrationalSetInfinity(lpexact->lpobjval);
      }
      else if( SCIPlpiExactIsPrimalInfeasible(lpexact->lpiexact) )
      {
         SCIPrationalSetInfinity(lpexact->lpobjval);
         lpexact->lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;
      }
      else if( SCIPlpiExactIsPrimalUnbounded(lpexact->lpiexact) )
      {
         SCIPrationalSetNegInfinity(lpexact->lpobjval);
         lpexact->lpsolstat = SCIP_LPSOLSTAT_UNBOUNDEDRAY;
      }
      else if( SCIPlpiExactIsIterlimExc(lpexact->lpiexact) )
      {
         lpexact->lpsolstat = SCIP_LPSOLSTAT_ITERLIMIT;
         lp->lpsolstat = SCIP_LPSOLSTAT_ITERLIMIT;
      }
      else if( SCIPlpiExactIsTimelimExc(lpexact->lpiexact) )
      {
         lpexact->lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
         lp->lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
      }
      else
      {
         SCIPdebugMessage("(node %" SCIP_LONGINT_FORMAT ") error or unknown return status of %s in LP %" SCIP_LONGINT_FORMAT " (internal status: %d)\n",
            stat->nnodes, &algo, stat->nlps, SCIPlpiExactGetInternalStatus(lpexact->lpiexact));
         lpexact->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         lpexact->solved = FALSE;
         *lperror = TRUE;
         return SCIP_OKAY;
      }
   }
   while( solveagain == TRUE );

   lpexact->solved = TRUE;

   SCIPsetDebugMsg(set, "solving exact LP with %d returned solstat=%d (internal status: %d, primalfeasible=%u, dualfeasible=%u)\n",
      algo, lpexact->lpsolstat, SCIPlpiExactGetInternalStatus(lpexact->lpiexact),
      SCIPlpiExactIsPrimalFeasible(lpexact->lpiexact), SCIPlpiExactIsDualFeasible(lpexact->lpiexact));

   return SCIP_OKAY;
}

/** solves the LP with simplex algorithm, and copy the solution into the column's data */
SCIP_RETCODE SCIPlpExactSolveAndEval(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Longint          itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool             usefarkas           /**< are we aiming to prove infeasibility? */
   )
{
   SCIP_RETCODE retcode;
   SCIP_Bool overwritefplp;
   SCIP_Bool primalfeasible;
   SCIP_Bool dualfeasible;
   SCIP_Bool* primalfeaspointer;
   SCIP_Bool* dualfeaspointer;
   int harditlim;
   SCIP_Bool farkasvalid;
   SCIP_Bool fromscratch;
   int iterations;
   SCIP_Real previoustime;

   assert(lp != NULL);
   assert(lpexact != NULL);
   assert(prob != NULL);
   assert(prob->nvars >= lp->ncols);
   assert(lperror != NULL);

   retcode = SCIP_OKAY;
   *lperror = FALSE;

   /* to avoid complications, we just always overwrite the fp lp */
   overwritefplp = TRUE;

   /* compute the limit for the number of LP resolving iterations, if needed (i.e. if limitresolveiters == TRUE) */
   harditlim = (int) MIN(itlim, INT_MAX);

   if( usefarkas )
   {
      previoustime = SCIPclockGetTime(stat->provedinfeaslptime);
      SCIPclockStart(stat->provedinfeaslptime, set);
   }
   else
   {
      previoustime = SCIPclockGetTime(stat->provedfeaslptime);
      SCIPclockStart(stat->provedfeaslptime, set);
   }

   /* set initial LP solver settings */
   fromscratch = FALSE;
   primalfeasible = FALSE;
   dualfeasible = FALSE;

   /* solve the LP */
   SCIP_CALL( lpExactFlushAndSolve(lpexact, blkmem, set, messagehdlr, stat,
         prob, eventqueue, harditlim, fromscratch, lperror) );
   assert(!(*lperror) || !lpexact->solved);

   SCIP_CALL( SCIPlpExactGetIterations(lpexact, &iterations) );

   if( usefarkas )
      SCIPstatAdd(stat, set, niterationsexlpinf, iterations);
   else
      SCIPstatAdd(stat, set, niterationsexlp, iterations);

   /* if not already done, solve again from scratch */
   if( *lperror )
   {
      lpExactNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "solve exact lp again from scratch");
      *lperror = FALSE;
      SCIP_CALL( lpExactFlushAndSolve(lpexact, blkmem, set, messagehdlr, stat,
         prob, eventqueue, harditlim, TRUE, lperror) );
   }

   SCIP_CALL( SCIPlpExactGetIterations(lpexact, &iterations) );
   if( usefarkas )
      SCIPstatAdd(stat, set, niterationsexlpinf, iterations);
   else
      SCIPstatAdd(stat, set, niterationsexlp, iterations);

   /* check for error */
   if( *lperror )
   {
      retcode = SCIP_OKAY;
      lp->hasprovedbound = FALSE;
      goto TERMINATE;
   }

   /* evaluate solution status */
   switch( lpexact->lpsolstat )
   {
   case SCIP_LPSOLSTAT_OPTIMAL:
      /* get LP solution and possibly check the solution's feasibility again */
      if( set->lp_checkprimfeas )
      {
         primalfeaspointer = &primalfeasible;
         lp->primalchecked = TRUE;
      }
      else
      {
         /* believe in the primal feasibility of the LP solution */
         primalfeasible = TRUE;
         primalfeaspointer = NULL;
         lp->primalchecked = FALSE;
      }
      if( set->lp_checkdualfeas )
      {
         dualfeaspointer = &dualfeasible;
         lp->dualchecked = TRUE;
      }
      else
      {
         /* believe in the dual feasibility of the LP solution */
         dualfeasible = TRUE;
         dualfeaspointer = NULL;
         lp->dualchecked = FALSE;
      }

      overwritefplp = overwritefplp || (lpexact->lpsolstat != lp->lpsolstat);
      SCIP_CALL( SCIPlpExactGetSol(lpexact, set, stat, primalfeaspointer, dualfeaspointer, overwritefplp) );

      lpexact->primalfeasible = primalfeasible && lpexact->primalfeasible;
      lpexact->dualfeasible = dualfeasible && lpexact->dualfeasible;

      if( primalfeasible && dualfeasible )
      {
         lp->lpobjval = SCIPrationalRoundReal(lpexact->lpobjval, SCIP_R_ROUND_DOWNWARDS);
         lp->hasprovedbound = TRUE;
      }
      else
      {
         /* print common begin of message */
         SCIPmessagePrintInfo(messagehdlr, "(node %" SCIP_LONGINT_FORMAT ") numerical troubles in exact LP %" SCIP_LONGINT_FORMAT " -- ",stat->nnodes, stat->nlps);
         lp->solved = FALSE;
         lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         *lperror = TRUE;
      }
      break;

   case SCIP_LPSOLSTAT_INFEASIBLE:
      SCIPsetDebugMsg(set, " -> LP infeasible\n");
      if( SCIPlpiExactHasDualRay(lpexact->lpiexact) )
      {
         SCIP_CALL( SCIPlpExactGetDualfarkas(lpexact, set, stat, &farkasvalid, overwritefplp) );
         lp->solved = TRUE;
         lp->lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;
         lp->lpobjval = SCIPsetInfinity(set);
         lp->hasprovedbound = farkasvalid;
      }
      else
      {
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %" SCIP_LONGINT_FORMAT ") infeasibility of LP %" SCIP_LONGINT_FORMAT " could not be proven by dual ray\n", stat->nnodes, stat->nlps);
         lp->solved = FALSE;
         lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         lpexact->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         farkasvalid = FALSE;
         *lperror = TRUE;
      }

      /* if the LP solver does not provide a Farkas proof we don't want to resolve the LP */
      if( !farkasvalid && !(*lperror) )
      {
         /* the Farkas proof does not prove infeasibility (this can happen due to numerical problems) and nothing
            * helped forget about the LP at this node and mark it to be unsolved
            */
         SCIPmessagePrintInfo(messagehdlr, "(node %" SCIP_LONGINT_FORMAT ") numerical troubles in exakt LP %" SCIP_LONGINT_FORMAT " -- ",stat->nnodes, stat->nlps);
         lp->solved = FALSE;
         lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         *lperror = TRUE;
      }

      break;

   case SCIP_LPSOLSTAT_UNBOUNDEDRAY:
      SCIPsetDebugMsg(set, " -> LP unbounded\n");
      SCIPwarningMessage(set->scip, "Exact LP solver returned unbounded ray: handling not fully supported.\n");
      break;

   case SCIP_LPSOLSTAT_OBJLIMIT:
      assert(!lpCutoffDisabled(set));
      /* Some LP solvers, e.g. CPLEX With FASTMIP setting, do not apply the final pivot to reach the dual solution
       * exceeding the objective limit. In some cases like branch-and-price, however, we must make sure that a dual
       * feasible solution exists that exceeds the objective limit. Therefore, we have to continue solving it without
       * objective limit for at least one iteration. We first try to continue with FASTMIP for one additional simplex
       * iteration using the steepest edge pricing rule. If this does not fix the problem, we temporarily disable
       * FASTMIP and solve again. */
      {
         SCIP_RATIONAL* objval;

         SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &objval) );
         /* actually, SCIPsetIsGE(set, lp->lpobjval, lp->lpiuobjlim) should hold, but we are a bit less strict in
          * the assert by using !SCIPsetIsFeasNegative()
          */

         SCIP_CALL( SCIPlpiExactGetObjval(lpexact->lpiexact, objval) );

         /* do one additional simplex step if the computed dual solution doesn't exceed the objective limit */
         if( SCIPrationalIsLTReal(objval, lpexact->lpiobjlim) )
         {
            SCIP_Real tmpcutoff;
            char tmppricingchar;
            SCIP_LPSOLSTAT solstat;

            SCIPrationalDebugMessage("objval = %q < %f = lp->lpiobjlim, but status objlimit\n", objval, lp->lpiobjlim);

            /* temporarily disable cutoffbound, which also disables the objective limit */
            tmpcutoff = lpexact->cutoffbound;
            lpexact->cutoffbound = SCIPsetInfinity(set);

            /* set lp pricing strategy to steepest edge */
            SCIP_CALL( SCIPsetGetCharParam(set, "lp/pricing", &tmppricingchar) );
            SCIP_CALL( SCIPsetSetCharParam(set, messagehdlr, "lp/pricing", 's') );

            /* resolve LP with an iteration limit of 1 */
            SCIP_CALL( lpExactFlushAndSolve(lpexact, blkmem, set, messagehdlr, stat, prob, eventqueue, 1, FALSE, lperror) );

            /* reinstall old cutoff bound and lp pricing strategy */
            lpexact->cutoffbound = tmpcutoff;
            SCIP_CALL( SCIPsetSetCharParam(set, messagehdlr, "lp/pricing", tmppricingchar) );

            /* get objective value */
            SCIP_CALL( SCIPlpiExactGetObjval(lpexact->lpiexact, objval) );

            /* get solution status for the lp */
            solstat = lpexact->lpsolstat;
            assert(solstat != SCIP_LPSOLSTAT_OBJLIMIT);

            if( !(*lperror) && solstat != SCIP_LPSOLSTAT_ERROR && solstat != SCIP_LPSOLSTAT_NOTSOLVED )
            {
               SCIPrationalDebugMessage(" ---> new objval = %q (solstat: %d, 1 add. step)\n", objval, solstat);
            }

            /* check for lp errors */
            if( *lperror || solstat == SCIP_LPSOLSTAT_ERROR || solstat == SCIP_LPSOLSTAT_NOTSOLVED )
            {
               SCIPsetDebugMsg(set, "unresolved error while resolving LP in order to exceed the objlimit\n");
               lp->solved = FALSE;
               lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
               lpexact->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
               lp->hasprovedbound = FALSE;

               retcode = *lperror ? SCIP_OKAY : SCIP_LPERROR;
               SCIPrationalFreeBuffer(set->buffer, &objval);
               goto TERMINATE;
            }

            lpexact->solved = TRUE;
            lp->hasprovedbound = TRUE;

            /* optimal solution / objlimit with fastmip turned off / itlimit or timelimit, but objlimit exceeded */
            if( solstat == SCIP_LPSOLSTAT_OPTIMAL || solstat == SCIP_LPSOLSTAT_OBJLIMIT
               || ( (solstat == SCIP_LPSOLSTAT_ITERLIMIT || solstat == SCIP_LPSOLSTAT_TIMELIMIT)
                  &&  SCIPrationalIsGEReal(objval, lpexact->cutoffbound - SCIPrationalRoundReal(getFiniteLooseObjvalExact(lpexact, set, prob), SCIP_R_ROUND_DOWNWARDS)) ) )
            {
               /* get LP solution and possibly check the solution's feasibility again */
               if( set->lp_checkprimfeas )
               {
                  primalfeaspointer = &primalfeasible;
                  lp->primalchecked = TRUE;
               }
               else
               {
                  /* believe in the primal feasibility of the LP solution */
                  primalfeasible = TRUE;
                  primalfeaspointer = NULL;
                  lp->primalchecked = FALSE;
               }
               if( set->lp_checkdualfeas )
               {
                  dualfeaspointer = &dualfeasible;
                  lp->dualchecked = TRUE;
               }
               else
               {
                  /* believe in the dual feasibility of the LP solution */
                  dualfeasible = TRUE;
                  dualfeaspointer = NULL;
                  lp->dualchecked = FALSE;
               }

               SCIP_CALL( SCIPlpExactGetSol(lpexact, set, stat, primalfeaspointer, dualfeaspointer, TRUE) );

               /* if objective value is larger than the cutoff bound, set solution status to objective
                * limit reached and objective value to infinity, in case solstat = SCIP_LPSOLSTAT_OBJLIMIT,
                * this was already done in the lpSolve() method
                */
               if( SCIPrationalIsGEReal(objval, lp->cutoffbound - SCIPrationalRoundReal(getFiniteLooseObjvalExact(lpexact, set, prob), SCIP_R_ROUND_DOWNWARDS)) )
               {
                  lpexact->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
                  lp->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
                  lp->hasprovedbound = TRUE;
                  lp->lpobjval = SCIPsetInfinity(set);
               }

               /* LP solution is not feasible or objective limit was reached without the LP value really exceeding
                * the cutoffbound; mark the LP to be unsolved
                */
               if( !primalfeasible || !dualfeasible
                  || (solstat == SCIP_LPSOLSTAT_OBJLIMIT &&
                     !SCIPrationalIsGEReal(objval, lp->cutoffbound -  SCIPrationalRoundReal(getFiniteLooseObjvalExact(lpexact, set, prob), SCIP_R_ROUND_DOWNWARDS))) )
               {
                  SCIPmessagePrintInfo(messagehdlr, "(node %" SCIP_LONGINT_FORMAT ") numerical troubles exact in LP %" SCIP_LONGINT_FORMAT " \n ", stat->nnodes, stat->nlps);
                  lp->solved = FALSE;
                  lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
                  lp->hasprovedbound = FALSE;
                  *lperror = TRUE;
               }
            }
            /* infeasible solution */
            else if( solstat == SCIP_LPSOLSTAT_INFEASIBLE )
            {
               SCIPsetDebugMsg(set, " -> LPexact infeasible\n");

               if( SCIPlpiExactHasDualRay(lpexact->lpiexact) )
               {
                  SCIP_CALL( SCIPlpExactGetDualfarkas(lpexact, set, stat, &farkasvalid, TRUE) );
               }
               /* it might happen that we have no infeasibility proof for the current LP (e.g. if the LP was always solved
                * with the primal simplex due to numerical problems) - treat this case like an LP error
                */
               else
               {
                  SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                     "(node %" SCIP_LONGINT_FORMAT ") infeasibility of exact LP %" SCIP_LONGINT_FORMAT " could not be proven by dual ray\n", stat->nnodes, stat->nlps);
                  lp->solved = FALSE;
                  lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
                  lp->hasprovedbound = FALSE;
                  farkasvalid = FALSE;
                  *lperror = TRUE;
               }
               if( !farkasvalid )
               {
                  /* the Farkas proof does not prove infeasibility (this can happen due to numerical problems) and nothing
                   * helped forget about the LP at this node and mark it to be unsolved
                   */
                  SCIPmessagePrintInfo(messagehdlr, "(node %" SCIP_LONGINT_FORMAT ") numerical troubles exact in LP %" SCIP_LONGINT_FORMAT " \n ", stat->nnodes, stat->nlps);
                  lp->solved = FALSE;
                  lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
                  lp->hasprovedbound = FALSE;
                  *lperror = TRUE;
               }
            }
            /* unbounded solution */
            else if( solstat == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
            {
               SCIP_Bool rayfeasible;

               if( set->lp_checkprimfeas )
               {
                  /* get unbounded LP solution and check the solution's feasibility again */
                  SCIP_CALL( SCIPlpExactGetUnboundedSol(lpexact, set, stat, &primalfeasible, &rayfeasible) );

                  lp->primalchecked = TRUE;
               }
               else
               {
                  /* get unbounded LP solution believing in its feasibility */
                  SCIP_CALL( SCIPlpExactGetUnboundedSol(lpexact, set, stat, NULL, NULL) );

                  rayfeasible = TRUE;
                  primalfeasible = TRUE;
                  lp->primalchecked = FALSE;
               }

               SCIPsetDebugMsg(set, " -> exact LP has unbounded primal ray\n");

               if( !primalfeasible || !rayfeasible )
               {
                  /* unbounded solution is infeasible (this can happen due to numerical problems):
                   * forget about the LP at this node and mark it to be unsolved
                   *
                   * @todo: like in the default LP solving evaluation, solve without fastmip,
                   * with tighter feasibility tolerance and from scratch
                   */
                  SCIPmessagePrintInfo(messagehdlr, "(node %" SCIP_LONGINT_FORMAT ") numerical troubles exact in LP %" SCIP_LONGINT_FORMAT " \n ",
                     stat->nnodes, stat->nlps);
                  lp->solved = FALSE;
                  lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
                  lp->hasprovedbound = FALSE;
                  *lperror = TRUE;
               }
            }

            assert(lp->lpsolstat != SCIP_LPSOLSTAT_ITERLIMIT);
            assert(SCIPrationalIsGEReal(objval, lp->cutoffbound - SCIPrationalRoundReal(getFiniteLooseObjvalExact(lpexact, set, prob), SCIP_R_ROUND_DOWNWARDS))
               || lp->lpsolstat != SCIP_LPSOLSTAT_OBJLIMIT);
         }
         else
         {
            overwritefplp = overwritefplp || (lpexact->lpsolstat != lp->lpsolstat);
            SCIP_CALL( SCIPlpExactGetSol(lpexact, set, stat, NULL, NULL, overwritefplp) );
            lp->hasprovedbound = TRUE;
         }

         SCIPrationalFreeBuffer(set->buffer, &objval);
      }
      SCIPsetDebugMsg(set, " -> LP objective limit reached\n");
      break;

   case SCIP_LPSOLSTAT_ITERLIMIT:
      SCIPsetDebugMsg(set, " -> LP iteration limit exceeded\n");
      break;

   case SCIP_LPSOLSTAT_TIMELIMIT:
      SCIPsetDebugMsg(set, " -> LP time limit exceeded\n");

      /* make sure that we evaluate the time limit exactly in order to avoid erroneous warning */
      stat->nclockskipsleft = 0;
      if( !SCIPsolveIsStopped(set, stat, FALSE) )
      {
         SCIPmessagePrintWarning(messagehdlr, "LP solver reached time limit, but SCIP time limit is not exceeded yet; "
            "you might consider switching the clock type of SCIP\n");
         stat->status = SCIP_STATUS_TIMELIMIT;
      }

      /* set the status of the floating point lp also to timelimit to avoid using the uncorrected bound */
      lp->lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
      lp->solved = TRUE;
      break;

   case SCIP_LPSOLSTAT_ERROR:
   case SCIP_LPSOLSTAT_NOTSOLVED:
      SCIPerrorMessage("error in LP solver\n");
      retcode = SCIP_LPERROR;
      goto TERMINATE;

   default:
      SCIPerrorMessage("unknown LP solution status\n");
      retcode = SCIP_ERROR;
      goto TERMINATE;
   }

TERMINATE:

   /* stop timing and update number of calls and fails, and proved bound status */
   if( usefarkas )
   {
      SCIPclockStop(stat->provedinfeaslptime, set);
      stat->nexlpinf++;
      if( *lperror )
         stat->timefailexlpinf += SCIPclockGetTime(stat->provedinfeaslptime) - previoustime;
   }
   else
   {
      SCIPclockStop(stat->provedfeaslptime, set);
      stat->nexlp++;
      if( *lperror )
         stat->timefailexlp += SCIPclockGetTime(stat->provedfeaslptime) - previoustime;
   }

   return retcode;
}

/*
 * row methods
 */

/** increases usage counter of LP row */
void SCIProwExactCapture(
   SCIP_ROWEXACT*        row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->nuses >= 0);
   assert(row->nlocks <= (unsigned int)(row->nuses)); /*lint !e574*/

   SCIPdebugMessage("capture row <%s> with nuses=%d and nlocks=%u\n", row->fprow->name, row->nuses, row->nlocks);
   row->nuses++;
}

/** output column to file stream */
void SCIProwExactPrint(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int r;

   assert(row != NULL);
   assert(row->fprow != NULL);

   SCIPmessageFPrintInfo(messagehdlr, file, "%s: ", row->fprow->name);
   SCIPrationalMessage(messagehdlr, file, row->lhs);
   SCIPmessageFPrintInfo(messagehdlr, file, " <= ");

   /* print coefficients */
   if( row->len == 0 )
      SCIPmessageFPrintInfo(messagehdlr, file, "<empty>");
   for( r = 0; r < row->len; ++r )
   {
      assert(SCIPvarGetName(row->cols[r]->var) != NULL);
      assert(SCIPvarGetStatus(row->cols[r]->var) == SCIP_VARSTATUS_COLUMN);
      if( SCIPrationalIsPositive(row->vals[r]) )
         SCIPmessageFPrintInfo(messagehdlr, file, "+ ");

      SCIPrationalMessage(messagehdlr, file, row->vals[r]);
      SCIPmessageFPrintInfo(messagehdlr, file, "(%g)<%s> ", SCIPrationalGetReal(row->vals[r]), SCIPvarGetName(row->cols[r]->var));
   }

   /* print constant */
   if( !SCIPrationalIsZero(row->constant) )
   {
      if( SCIPrationalIsPositive(row->constant) )
         SCIPmessageFPrintInfo(messagehdlr, file, "+");
      SCIPrationalMessage(messagehdlr, file, row->constant);
   }

   SCIPmessageFPrintInfo(messagehdlr, file, "<= , ");
   SCIPrationalMessage(messagehdlr, file, row->rhs);
   SCIPmessageFPrintInfo(messagehdlr, file, "\n");
}

/** get the index of an exact row */
int SCIProwExactGetIndex(
   SCIP_ROWEXACT*        row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->index;
}

/** gets the length of a row */
int SCIProwExactGetNNonz(
   SCIP_ROWEXACT*        row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->len;
}

/** gets array with coefficients of nonzero entries */
SCIP_RATIONAL** SCIProwExactGetVals(
   SCIP_ROWEXACT*        row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->vals;
}

/** gets array of exact columns */
SCIP_COLEXACT** SCIProwExactGetCols(
   SCIP_ROWEXACT*        row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->cols;
}

/** returns TRUE iff row is member of current LP */
SCIP_Bool SCIProwExactIsInLP(
   SCIP_ROWEXACT*        row                 /**< LP row */
   )
{
   assert(row != NULL);

   return (row->lppos >= 0);
}

/** return TRUE iff row is modifiable */
SCIP_Bool SCIProwExactIsModifiable(
   SCIP_ROWEXACT*        row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->fprow != NULL);

   return row->fprow->modifiable;
}

/** returns true, if an exact row for this fprow was already created */
SCIP_Bool SCIProwHasExRow(
   SCIP_LPEXACT*         lpexact,            /**< exact lp data structure */
   SCIP_ROW*             row                 /**< SCIP row */
   )
{
   assert(row != NULL);
   assert(lpexact != NULL);

   return (NULL != row->rowexact);
}

/** returns fp row corresponding to exact row, if it exists. Otherwise returns NULL */
SCIP_ROW* SCIProwExactGetRow(
   SCIP_ROWEXACT*        row                 /**< SCIP row */
   )
{
   assert(row != NULL);

   return row->fprow;
}

/** returns rhs-relaxation part of exact row, if it exists. Otherwise returns NULL */
SCIP_ROW* SCIProwExactGetRowRhs(
   SCIP_ROWEXACT*        row                 /**< SCIP row */
   )
{
   assert(row != NULL);

   return row->fprowrhs;
}

/** true if row can be relaxed (possibly as two fp rows) */
SCIP_Bool SCIProwExactHasFpRelax(
   SCIP_ROWEXACT*        row                 /**< SCIP row */
   )
{
   assert(row != NULL);

   return row->fprelaxable;
}

/** returns exact col corresponding to fpcol, if it exists. Otherwise returns NULL */
SCIP_COLEXACT* SCIPcolGetColExact(
   SCIP_COL*             col                 /**< SCIP col */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(col->var->exactdata != NULL);

   return col->var->exactdata->colexact;
}

/** calculates the Farkas coefficient y^T A_i or reduced cost c - y^T A_i of a column i using the given dual Farkas vector y */
SCIP_RETCODE SCIPcolExactCalcFarkasRedcostCoef(
   SCIP_COLEXACT*        col,                /**< LP column */
   SCIP_SET*             set,                /**< SCIP settings pointer */
   SCIP_RATIONAL*        result,             /**< rational to store the result */
   SCIP_RATIONAL**       dual,               /**< dense dual vector, NULL to use internal row-values */
   SCIP_Bool             usefarkas           /**< should the farkas coefficient be computed ? */
   )
{
   SCIP_ROWEXACT* row;
   SCIP_RATIONAL* val;
   SCIP_RATIONAL* tmp;
   int i;

   assert(col != NULL);
   assert(SCIPvarGetStatusExact(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetColExact(col->var) == col);

   if( usefarkas )
      SCIPrationalSetFraction(result, 0LL, 1LL);
   else
      SCIPrationalSetRational(result, col->obj);

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );

   for( i = 0; i < col->nlprows; ++i )
   {
      row = col->rows[i];
      assert(row != NULL);
      assert(row->lppos >= 0);

      if( usefarkas )
         val = (dual == NULL) ? row->dualfarkas : dual[row->lppos];
      else
         val = (dual == NULL) ? row->dualsol : dual[row->lppos];

      assert(!SCIPrationalIsInfinity(val));

      SCIPrationalMult(tmp, col->vals[i], val);
      if( usefarkas )
         SCIPrationalAdd(result, result, tmp);
      else
         SCIPrationalDiff(result, result, tmp);
   }

   if( col->nunlinked > 0 )
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->lppos == -1 || col->linkpos[i] == -1);
         if( row->lppos >= 0 )
         {
            if( usefarkas )
               val = (dual == NULL) ? row->dualfarkas : dual[row->lppos];
            else
               val = (dual == NULL) ? row->dualsol : dual[row->lppos];

            SCIPrationalMult(tmp, col->vals[i], val);
            if( usefarkas )
               SCIPrationalAdd(result, result, tmp);
            else
               SCIPrationalDiff(result, result, tmp);
         }
      }
   }
#ifndef NDEBUG
   else
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->lppos == -1);
         assert(col->linkpos[i] >= 0);
         if( dual == NULL )
            assert((usefarkas && SCIPrationalIsZero(row->dualfarkas)) || SCIPrationalIsZero(row->dualsol));
      }
      assert(!SCIPrationalIsPositive(result) || !SCIPrationalIsInfinity(col->ub));
      assert(!SCIPrationalIsNegative(result) || !SCIPrationalIsNegInfinity(col->lb));
   }
#endif

   SCIPrationalFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** adds a previously non existing coefficient to an LP row */
SCIP_RETCODE SCIProwExactAddCoef(
   SCIP_ROWEXACT*        rowexact,           /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_COLEXACT*        colexact,           /**< LP column */
   SCIP_RATIONAL*        val                 /**< value of coefficient */
   )
{
   assert(rowexact != NULL);
   assert(colexact != NULL);
   assert(lpexact != NULL);

   assert(lpexact != NULL);
   assert(!lpexact->fplp->diving || rowexact->fprow->lppos == -1);

   SCIP_CALL( rowExactAddCoef(rowexact, blkmem, set, eventqueue, lpexact, colexact, val, -1) );

   checkLinks(lpexact);

   return SCIP_OKAY;
}

/** deletes coefficient from row */
SCIP_RETCODE SCIProwExactDelCoef(
   SCIP_ROWEXACT*        row,                /**< row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_COLEXACT*        col                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(row != NULL);
   assert(!row->delaysort);
   assert(lpexact != NULL);
   assert(!lpexact->fplp->diving || row->lppos == -1);
   assert(col != NULL);
   assert(col->var != NULL);

   /* search the position of the column in the row's col vector */
   pos = rowExactSearchCoef(row, col);
   if( pos == -1 )
   {
      SCIPerrorMessage("coefficient for column <%s> doesn't exist in row <%s>\n", SCIPvarGetName(col->var), row->fprow->name);
      return SCIP_INVALIDDATA;
   }
   assert(0 <= pos && pos < row->len);
   assert(row->cols[pos] == col);
   assert(row->cols_index[pos] == col->index);

   /* if column knows of the row, remove the row from the column's row vector */
   if( row->linkpos[pos] >= 0 )
   {
      assert(col->rows[row->linkpos[pos]] == row);
      assert(SCIPrationalIsEQ(col->vals[row->linkpos[pos]], row->vals[pos]));
      SCIP_CALL( colExactDelCoefPos(col, set, lpexact, row->linkpos[pos]) );
   }

   /* delete the column from the row's col vector */
   SCIP_CALL( rowExactDelCoefPos(row, set, lpexact, pos) );

   checkLinks(lpexact);

   return SCIP_OKAY;
}

/** changes or adds a coefficient to an LP row */
SCIP_RETCODE SCIProwExactChgCoef(
   SCIP_ROWEXACT*        row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_COLEXACT*        col,                /**< LP column */
   SCIP_RATIONAL*        val                 /**< value of coefficient */
   )
{
   int pos;

   assert(row != NULL);
   assert(!row->delaysort);
   assert(lpexact != NULL);
   assert(!lpexact->fplp->diving || row->lppos == -1);
   assert(col != NULL);

   /* search the position of the column in the row's col vector */
   pos = rowExactSearchCoef(row, col);

   /* check, if column already exists in the row's col vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( rowExactAddCoef(row, blkmem, set, eventqueue, lpexact, col, val, -1) );
   }
   else
   {
      /* modify already existing coefficient */
      assert(0 <= pos && pos < row->len);
      assert(row->cols[pos] == col);
      assert(row->cols_index[pos] == col->index);

      /* if column knows of the row, change the corresponding coefficient in the column */
      if( row->linkpos[pos] >= 0 )
      {
         assert(col->rows[row->linkpos[pos]] == row);
         assert(SCIPrationalIsEQ(col->vals[row->linkpos[pos]], row->vals[pos]));
         SCIP_CALL( colExactChgCoefPos(col, set, lpexact, row->linkpos[pos], val) );
      }

      /* change the coefficient in the row */
      SCIP_CALL( rowExactChgCoefPos(row, set, lpexact, pos, val) );
   }

   checkLinks(lpexact);

   return SCIP_OKAY;
}

/** increases value of an existing or non-existing coefficient in an LP row */
SCIP_RETCODE SCIProwExactIncCoef(
   SCIP_ROWEXACT*        row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_COLEXACT*        col,                /**< LP column */
   SCIP_RATIONAL*        incval              /**< value to add to the coefficient */
   )
{
   int pos;
   SCIP_RATIONAL* tmp;

   assert(row != NULL);
   assert(lpexact != NULL);
   assert(!lpexact->fplp->diving || row->lppos == -1);
   assert(col != NULL);

   if( SCIPrationalIsZero(incval) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );

   /* search the position of the column in the row's col vector */
   pos = rowExactSearchCoef(row, col);

   /* check, if column already exists in the row's col vector */
   if( pos == -1 )
   {
      /* coefficient doesn't exist, or sorting is delayed: add coefficient to the end of the row's arrays */
      SCIP_CALL( rowExactAddCoef(row, blkmem, set, eventqueue, lpexact, col, incval, -1) );
   }
   else
   {
      /* modify already existing coefficient */
      assert(0 <= pos && pos < row->len);
      assert(row->cols[pos] == col);
      assert(row->cols_index[pos] == col->index);

      SCIPrationalAdd(tmp, incval, row->vals[pos]);

      /* if column knows of the row, change the corresponding coefficient in the column */
      if( row->linkpos[pos] >= 0 )
      {
         assert(col->rows[row->linkpos[pos]] == row);
         assert(SCIPrationalIsEQ(col->vals[row->linkpos[pos]], row->vals[pos]));
         SCIP_CALL( colExactChgCoefPos(col, set, lpexact, row->linkpos[pos], tmp) );
      }

      /* change the coefficient in the row */
      SCIP_CALL( rowExactChgCoefPos(row, set, lpexact, pos, tmp) );
   }

   checkLinks(lpexact);

   /* invalid the activity */
   row->validactivitylp = -1;

   SCIPrationalFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** changes constant value of a row */
SCIP_RETCODE SCIProwExactChgConstant(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_RATIONAL*        constant            /**< new constant value */
   )
{
   assert(row != NULL);
   assert(SCIPrationalIsLE(row->lhs, row->rhs));
   assert(!SCIPrationalIsAbsInfinity(constant));
   assert(stat != NULL);
   assert(lpexact != NULL);
   assert(!lpexact->fplp->diving || row->fprow->lppos == -1);

   if( !SCIPrationalIsEQ(constant, row->constant) )
   {
      if( row->fprow->validpsactivitydomchg == stat->domchgcount )
      {
         assert(!SCIPrationalIsInfinity(row->pseudoactivity));
         SCIPrationalAdd(row->pseudoactivity, row->pseudoactivity, constant);
         SCIPrationalDiff(row->pseudoactivity, row->pseudoactivity, row->constant);
      }

      SCIPrationalSetRational(row->constant, constant);
      SCIPintervalSetRational(&row->constantreal, constant);
   }

   return SCIP_OKAY;
}

/** add constant value to a row */
SCIP_RETCODE SCIProwExactAddConstant(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_RATIONAL*        addval              /**< constant value to add to the row */
   )
{
   SCIP_RATIONAL* tmp;

   assert(row != NULL);
   assert(SCIPrationalIsLE(row->lhs, row->rhs));
   assert(!SCIPrationalIsAbsInfinity(addval));
   assert(stat != NULL);
   assert(lpexact != NULL);
   assert(!lpexact->fplp->diving || row->fprow->lppos == -1);

   if( !SCIPrationalIsZero(addval) )
   {
      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );
      SCIPrationalAdd(tmp, row->constant, addval);
      SCIP_CALL( SCIProwExactChgConstant(row, stat, lpexact, tmp) );

      SCIPrationalFreeBuffer(set->buffer, &tmp);
   }

   return SCIP_OKAY;
}

/** returns the feasibility of a row for the given solution */
SCIP_RETCODE SCIProwExactGetSolFeasibility(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_RATIONAL*        result              /**< result pointer */
   )
{
   SCIP_RATIONAL* temp1;
   SCIP_RATIONAL* temp2;

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &temp1) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &temp2) );

   assert(row != NULL);

   SCIP_CALL( SCIProwExactGetSolActivity(row, set, stat, sol, FALSE, result) );

   SCIPrationalDiff(temp1, row->rhs, result);
   SCIPrationalDiff(temp2, result, row->lhs);
   SCIPrationalMin(result, temp1, temp2);

   SCIPrationalFreeBuffer(set->buffer, &temp2);
   SCIPrationalFreeBuffer(set->buffer, &temp1);

   return SCIP_OKAY;
}

/** does activity computation with running error analysis for a row, return TRUE on success */
SCIP_Bool SCIProwExactGetSolActivityWithErrorbound(
   SCIP_ROWEXACT*        rowexact,           /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            activity,           /**< the approximate activity */
   SCIP_Real*            errorbound          /**< the error bound */
   )
{
   SCIP_ROW* row;
   SCIP_Real solval;
   SCIP_Real mu;
   SCIP_Real sum;
   int c;

   assert(rowexact->fprow != NULL);

   row = rowexact->fprow;

   if( row->len != rowexact->len )
      return FALSE;

   sum = 0.0;
   mu = 0.0;

   for( c = 0; c < row->len; c++ )
   {
      if( sol != NULL)
         solval = SCIPsolGetVal(sol, set, stat, SCIPcolGetVar(row->cols[c]));
      else
         solval = row->cols[c]->primsol;

      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         return FALSE;

      sum += row->vals[c] * solval;
      mu += REALABS(sum);
      /* the factor 3 + eps is needed to account for rounding errors in valsreal[v]/solval */
      mu += (3.0 + SCIP_REAL_UNITROUNDOFF) * REALABS(row->vals[c] * solval);
   }

   sum += row->constant;
   mu += (3.0 + SCIP_REAL_UNITROUNDOFF) * REALABS(row->constant);

   sum = MAX(sum, -SCIPsetInfinity(set)); /*lint !e666*/
   sum = MIN(sum, SCIPsetInfinity(set)); /*lint !e666*/

   *activity = sum;
   *errorbound = mu;

   return TRUE;
}

/** returns the activity of a row for a given solution */
SCIP_RETCODE SCIProwExactGetSolActivity(
   SCIP_ROWEXACT*        rowexact,           /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             useexact,           /**< should an exact solution be used */
   SCIP_RATIONAL*        result              /**< resulting activity */
   )
{
   SCIP_COLEXACT* colexact;
   SCIP_RATIONAL* solval;
   int i;

   assert(rowexact != NULL);

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &solval) );
   SCIPrationalSetRational(result, rowexact->constant);
   for( i = 0; i < rowexact->len; ++i )
   {
      colexact = rowexact->cols[i];

      assert(colexact != NULL);

      assert((i < rowexact->nlpcols) == (rowexact->linkpos[i] >= 0
         && colexact->lppos >= 0));

      if( useexact )
         SCIPsolGetValExact(solval, sol, set, stat, colexact->var);
      else
         SCIPrationalSetReal(solval, SCIPsolGetVal(sol, set, stat, colexact->var));

      if( SCIPrationalIsAbsInfinity(solval) ) /*lint !e777*/
      {
         if( SCIPrationalIsNegInfinity(rowexact->lhs) )
            SCIPrationalIsPositive(rowexact->vals[i]) ? SCIPrationalSetRational(solval, colexact->lb) : SCIPrationalSetRational(solval, colexact->ub);
         else if( SCIPrationalIsInfinity(rowexact->rhs) )
            SCIPrationalIsPositive(rowexact->vals[i]) ? SCIPrationalSetRational(solval, colexact->ub) : SCIPrationalSetRational(solval, colexact->lb);
         else
         {
            SCIPrationalAdd(solval, colexact->lb, colexact->ub);
            SCIPrationalMultReal(solval, solval, 0.5);
         }
      }

      SCIPrationalMult(solval, solval, rowexact->vals[i]);
      SCIPrationalAdd(result, result, solval);
   }

   SCIPrationalFreeBuffer(set->buffer, &solval);

   return SCIP_OKAY;
}

/** decreases usage counter of LP row, and frees memory if necessary */
SCIP_RETCODE SCIProwExactRelease(
   SCIP_ROWEXACT**       row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact             /**< current LP data */
   )
{
   assert(blkmem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->nuses >= 1);
   assert((*row)->nlocks < (unsigned int)((*row)->nuses)); /*lint !e574*/

   SCIPsetDebugMsg(set, "release exact row <%p> with nuses=%d and nlocks=%u\n",
      (void*) (*row), (*row)->nuses, (*row)->nlocks);
   (*row)->nuses--;
   if( (*row)->nuses == 0 )
   {
      SCIP_CALL( SCIProwExactFree(row, blkmem, set, lpexact) );
   }

   *row = NULL;

   return SCIP_OKAY;
}

/** frees an LP row */
SCIP_RETCODE SCIProwExactFree(
   SCIP_ROWEXACT**       row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact             /**< current LP data */
   )
{
   assert(blkmem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->nuses == 0);
   assert((*row)->lppos == -1);

   /* remove column indices from corresponding rows */
   SCIP_CALL( rowExactUnlink(*row, set, lpexact) );

   if( (*row)->storedsolvals != NULL )
   {
      SCIPrationalFreeBlock(blkmem, &(*row)->storedsolvals->activity);
      SCIPrationalFreeBlock(blkmem, &(*row)->storedsolvals->dualsol);
      BMSfreeBlockMemoryNull(blkmem, &(*row)->storedsolvals);
   }

   SCIPrationalFreeBlock(blkmem, &(*row)->constant);
   SCIPrationalFreeBlock(blkmem, &(*row)->lhs);
   SCIPrationalFreeBlock(blkmem, &(*row)->rhs);
   SCIPrationalFreeBlock(blkmem, &(*row)->flushedlhs);
   SCIPrationalFreeBlock(blkmem, &(*row)->flushedrhs);
   SCIPrationalFreeBlock(blkmem, &(*row)->objprod);
   SCIPrationalFreeBlock(blkmem, &(*row)->dualsol);
   SCIPrationalFreeBlock(blkmem, &(*row)->activity);
   SCIPrationalFreeBlock(blkmem, &(*row)->dualfarkas);
   SCIPrationalFreeBlock(blkmem, &(*row)->pseudoactivity);

   SCIPrationalFreeBlockArray(blkmem, &(*row)->vals, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->valsinterval, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->cols, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->cols_index, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->linkpos, (*row)->size);
   BMSfreeBlockMemory(blkmem, row);

   return SCIP_OKAY;
}

/** returns the feasibility of a row in the current LP solution: negative value means infeasibility */
SCIP_RETCODE SCIProwExactGetLPFeasibility(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_RATIONAL*        result              /**< rational pointer to store the result */
   )
{
   SCIP_RATIONAL* activity;
   SCIP_RATIONAL* actrhs;
   SCIP_RATIONAL* actlhs;

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &actrhs) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &actlhs) );
   assert(row != NULL);

   activity = SCIProwExactGetLPActivity(row, stat, lpexact);

   SCIPrationalDiff(actlhs, row->rhs, activity);
   SCIPrationalDiff(actrhs, activity, row->lhs);
   SCIPrationalMin(result, actrhs, actlhs);

   SCIPrationalFreeBuffer(set->buffer, &actlhs);
   SCIPrationalFreeBuffer(set->buffer, &actrhs);

   return SCIP_OKAY;
}

/** returns the pseudo feasibility of a row in the current pseudo solution: negative value means infeasibility */
SCIP_RETCODE SCIProwExactGetPseudoFeasibility(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_RATIONAL*        result              /**< rational pointer to store the result */
   )
{
   SCIP_RATIONAL* pseudoactivity;
   SCIP_RATIONAL* actrhs;
   SCIP_RATIONAL* actlhs;

   assert(row != NULL);

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &actrhs) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &actlhs) );

   pseudoactivity = SCIProwExactGetPseudoActivity(row, stat);

   SCIPrationalDiff(actlhs, row->rhs, pseudoactivity);
   SCIPrationalDiff(actrhs, pseudoactivity, row->lhs);
   SCIPrationalMin(result, actrhs, actlhs);

   SCIPrationalFreeBuffer(set->buffer, &actlhs);
   SCIPrationalFreeBuffer(set->buffer, &actrhs);

   return SCIP_OKAY;
}

/** returns the activity of a row in the current LP solution */
SCIP_RATIONAL* SCIProwExactGetLPActivity(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEXACT*         lpexact             /**< current LP data */
   )
{
   assert(row != NULL);
   assert(stat != NULL);
   assert(lpexact != NULL);
   assert(row->fprow->validactivitylp <= stat->lpcount);
   assert(lpexact->fplp->validsollp == stat->lpcount);

   if( row->fprow->validactivitylp != stat->lpcount )
      SCIProwExactRecalcLPActivity(row, stat);
   assert(row->fprow->validactivitylp == stat->lpcount);
   assert(row->fprow->activity < SCIP_INVALID);

   return row->activity;
}

/** returns the pseudo activity of a row in the current pseudo solution */
SCIP_RATIONAL* SCIProwExactGetPseudoActivity(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(row != NULL);
   assert(stat != NULL);
   assert(row->fprow->validpsactivitydomchg <= stat->domchgcount);

   /* check, if pseudo activity has to be calculated */
   if( row->fprow->validpsactivitydomchg != stat->domchgcount )
      SCIProwExactRecalcPseudoActivity(row, stat);
   assert(row->fprow->validpsactivitydomchg == stat->domchgcount);
   assert(row->fprow->pseudoactivity < SCIP_INVALID);

   return row->pseudoactivity;
}

/** sorts row entries such that LP columns precede non-LP columns and inside both parts lower column indices precede
 *  higher ones
 */
void SCIProwExactSort(
   SCIP_ROWEXACT*        row                 /**< row to be sorted */
   )
{
   assert(row != NULL);

   /* sort LP columns */
   rowExactSortLP(row);

   /* sort non-LP columns */
   rowExactSortNonLP(row);
}

/** sorts row, and merges equal column entries (resulting from lazy sorting and adding) into a single entry; removes
 *  zero entries from row the row must not be linked to the columns; otherwise, we would need to update the columns as
 *  well, which is too expensive
 */
static
void rowExactMerge(
   SCIP_ROWEXACT*        row,                /**< row to be sorted */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);
   assert(!row->delaysort);
   assert(row->nunlinked == row->len);
   assert(row->nlpcols == 0);

   SCIPsetDebugMsg(set, "merging row <%s>\n", row->fprow->name);

   /* do nothing on empty rows; if row is sorted, nothing has to be done */
   if( row->len > 0 && (!row->lpcolssorted || !row->nonlpcolssorted) )
   {
      SCIP_COLEXACT** cols;
      int* cols_index;
      SCIP_RATIONAL** vals;
      int s;
      int t;

      /* make sure, the row is sorted */
      SCIProwExactSort(row);
      assert(row->lpcolssorted);
      assert(row->nonlpcolssorted);

      /* merge equal columns, thereby recalculating whether the row's activity is always integral */
      cols = row->cols;
      cols_index = row->cols_index;
      vals = row->vals;
      assert(cols != NULL);
      assert(cols_index != NULL);
      assert(vals != NULL);

      t = 0;
      row->integral = TRUE;
      assert(!SCIPrationalIsZero(vals[0]));
      assert(row->linkpos[0] == -1);

      for( s = 1; s < row->len; ++s )
      {
         assert(!SCIPrationalIsZero(vals[s]));
         assert(row->linkpos[s] == -1);

         if( cols[s] == cols[t] )
         {
            /* merge entries with equal column */
            SCIPrationalAdd(vals[t], vals[t], vals[s]);
            SCIPintervalSetRational(&row->valsinterval[t], vals[t]);
         }
         else
         {
            /* go to the next entry, overwriting current entry if coefficient is zero */
            if( !SCIPrationalIsZero(vals[t]) )
            {
               row->integral = row->integral && SCIPcolIsIntegral(cols[t]->fpcol) && SCIPrationalIsIntegral(vals[t]);
               t++;
            }
            cols[t] = cols[s];
            cols_index[t] = cols_index[s];
            SCIPrationalSetRational(vals[t], vals[s]);
            SCIPintervalSetRational(&row->valsinterval[t], vals[t]);
         }
      }
      if( !SCIPrationalIsZero(vals[t]) )
      {
         row->integral = row->integral && SCIPcolIsIntegral(cols[t]->fpcol) && SCIPrationalIsIntegral(vals[t]);
         t++;
      }
      assert(s == row->len);
      assert(t <= row->len);

      row->len = t;
      row->nunlinked = t;
   }

#ifndef NDEBUG
   /* check for double entries */
   {
      int i;
      int j;

      for( i = 0; i < row->len; ++i )
      {
         assert(row->cols[i] != NULL);
         assert(row->cols[i]->index == row->cols_index[i]);
         for( j = i+1; j < row->len; ++j )
            assert(row->cols[i] != row->cols[j]);
      }
   }
#endif
}

/** enables delaying of row sorting */
void SCIProwExactDelaySort(
   SCIP_ROWEXACT*        rowexact            /**< LP rowexact */
   )
{
   assert(rowexact != NULL);
   assert(!rowexact->delaysort);

   rowexact->delaysort = TRUE;
}

/** disables delaying of row sorting, sorts row and merges coefficients with equal columns */
void SCIProwExactForceSort(
   SCIP_ROWEXACT*        rowexact,           /**< LP rowexact */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(rowexact != NULL);
   assert(rowexact->delaysort);

   rowexact->delaysort = FALSE;
   rowExactMerge(rowexact, set);
}

/** recalculates the current activity of a row */
void SCIProwExactRecalcLPActivity(
   SCIP_ROWEXACT*        rowexact,           /**< LP row */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_COLEXACT* colexact;
   SCIP_COL* col;
   SCIP_ROW* row;
   int c;

   assert(rowexact != NULL);

   row = rowexact->fprow;

   assert(row != NULL);
   assert(stat != NULL);

   SCIPrationalSetRational(rowexact->activity, rowexact->constant);
   for( c = 0; c < row->nlpcols; ++c )
   {
      colexact = rowexact->cols[c];
      col = row->cols[c];

      assert(col != NULL);
      assert(colexact != NULL);
      assert(!SCIPrationalIsInfinity(colexact->primsol));
      assert(col->lppos >= 0);
      assert(row->linkpos[c] >= 0);

      SCIPrationalAddProd(rowexact->activity, rowexact->vals[c], colexact->primsol);
   }

   if( row->nunlinked > 0 )
   {
      for( c = row->nlpcols; c < row->len; ++c )
      {
         col = row->cols[c];
         colexact = rowexact->cols[c];

         assert(col != NULL);
         assert(colexact != NULL);
         assert(col->lppos >= 0 || col->primsol == 0.0);
         assert(col->lppos == -1 || row->linkpos[c] == -1);
         if( col->lppos >= 0 )
            SCIPrationalAddProd(rowexact->activity, rowexact->vals[c], colexact->primsol);
      }
   }
#ifndef NDEBUG
   else
   {
      for( c = row->nlpcols; c < row->len; ++c )
      {
         col = row->cols[c];
         colexact = rowexact->cols[c];

         assert(col != NULL);
         assert(colexact != NULL);
         assert(SCIPrationalIsZero(colexact->primsol));
         assert(col->lppos == -1);
         assert(row->linkpos[c] >= 0);
      }
   }
#endif

   row->activity = SCIPrationalGetReal(rowexact->activity);
   row->validactivitylp = stat->lpcount;
}

/** calculates the current pseudo activity of a row */
void SCIProwExactRecalcPseudoActivity(
   SCIP_ROWEXACT*        rowexact,           /**< row data */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_COLEXACT* colexact;
   SCIP_ROW* row;

   int i;

   assert(rowexact != NULL);

   row = rowexact->fprow;

   assert(row != NULL);
   assert(stat != NULL);

   SCIPrationalSetRational(rowexact->pseudoactivity, rowexact->constant);
   for( i = 0; i < row->len; ++i )
   {
      colexact = rowexact->cols[i];

      assert(colexact->fpcol != NULL);
      assert(colexact != NULL);
      assert((i < row->nlpcols) == (row->linkpos[i] >= 0 && colexact->fpcol->lppos >= 0));
      assert(colexact->fpcol->var != NULL);
      assert(SCIPvarGetStatus(colexact->fpcol->var) == SCIP_VARSTATUS_COLUMN);

      SCIPrationalAddProd(rowexact->pseudoactivity, rowexact->vals[i], SCIPcolExactGetBestBound(colexact));
   }

   row->validpsactivitydomchg = stat->domchgcount;
   row->pseudoactivity = SCIPrationalGetReal(rowexact->pseudoactivity);
}

/** gets objective value of column */
SCIP_RATIONAL* SCIPcolExactGetObj(
   SCIP_COLEXACT*        col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->obj;
}

/** gets lower bound of column */
SCIP_RATIONAL* SCIPcolExactGetLb(
   SCIP_COLEXACT*        col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->lb;
}

/** gets upper bound of column */
SCIP_RATIONAL* SCIPcolExactGetUb(
   SCIP_COLEXACT*        col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->ub;
}

/** gets best bound of column with respect to the objective function */
SCIP_RATIONAL* SCIPcolExactGetBestBound(
   SCIP_COLEXACT*        col                 /**< LP column */
   )
{
   assert(col != NULL);

   if( SCIPrationalIsPositive(col->obj) || SCIPrationalIsZero(col->obj) )
      return col->lb;
   else
      return col->ub;
}

/** gets the primal LP solution of a column */
SCIP_RATIONAL* SCIPcolExactGetPrimsol(
   SCIP_COLEXACT*        col                 /**< LP column */
   )
{
   assert(col != NULL);

   if( col->fpcol->lppos >= 0 )
      return col->primsol;
   else
      return NULL;
}

/** gets variable this column represents */
SCIP_VAR* SCIPcolExactGetVar(
   SCIP_COLEXACT*        col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->var;
}

/** ensures, that column array of row can store at least num entries */
SCIP_RETCODE SCIProwExactEnsureSize(
   SCIP_ROWEXACT*        row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(row != NULL);
   assert(row->fprow != NULL);
   assert(row->len <= row->size);

   if( num > row->size )
   {
      int newsize;
      int i;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->cols, row->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->cols_index, row->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->vals, row->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->valsinterval, row->size, newsize) );
      for( i = row->size; i < newsize; ++i )
         SCIP_CALL( SCIPrationalCreateBlock(blkmem, &row->vals[i]) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->linkpos, row->size, newsize) );
      row->size = newsize;
   }
   assert(num <= row->size);

   return SCIP_OKAY;
}

/*
 * lp update methods
 */


/** compute the objective delta due the new objective coefficient */
static
SCIP_RETCODE getObjvalDeltaObjExact(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RATIONAL*        oldobj,             /**< old objective value of variable */
   SCIP_RATIONAL*        newobj,             /**< new objective value of variable */
   SCIP_RATIONAL*        lb,                 /**< lower bound of variable */
   SCIP_RATIONAL*        ub,                 /**< upper bound of variable */
   SCIP_RATIONAL*        deltaval,           /**< pointer to store the delta value */
   int*                  deltainf            /**< pointer to store the number of variables with infinite best bound */
   )
{
   SCIP_RATIONAL* tmp;
   assert(!SCIPrationalIsAbsInfinity(oldobj));
   assert(!SCIPrationalIsAbsInfinity(newobj));
   assert(!SCIPrationalIsInfinity(lb));
   assert(!SCIPrationalIsNegInfinity(ub));
   assert(!SCIPrationalIsEQ(oldobj, newobj));

   SCIPrationalSetReal(deltaval, 0.0);
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );
   (*deltainf) = 0;

   if( SCIPrationalIsPositive(oldobj) )
   {
      /* sign of objective did not change */
      if( SCIPrationalIsPositive(newobj) )
      {
         /* if the bound is finite, calculate the deltaval */
         if( !SCIPrationalIsNegInfinity(lb) )
         {
            SCIPrationalDiff(deltaval, newobj, oldobj);
            SCIPrationalMult(deltaval, deltaval, lb);
         }
      }
      /* sign of objective did change, so the best bound does change */
      else if( SCIPrationalIsNegative(newobj) )
      {
         if( SCIPrationalIsNegInfinity(lb) )
         {
            /* old best bound was infinite while new one is not */
            if( !SCIPrationalIsInfinity(ub) )
            {
               (*deltainf) = -1;
               SCIPrationalMult(deltaval, ub, newobj);
            }
         }
         else
         {
            /* new best bound is infinite while old one was not */
            if( SCIPrationalIsInfinity(ub) )
            {
               (*deltainf) = 1;
               SCIPrationalMult(deltaval, lb, oldobj);
               SCIPrationalNegate(deltaval, deltaval);
            }
            /* neither old nor new best bound is infinite, so just calculate the deltaval */
            else
            {
               SCIPrationalMult(tmp, lb, oldobj);
               SCIPrationalMult(deltaval, ub, newobj);

               SCIPrationalDiff(deltaval, deltaval, tmp);
            }
         }
      }
      /* new objective is 0.0 */
      else
      {
         if( SCIPrationalIsNegInfinity(lb) )
            (*deltainf) = -1;
         else
         {
            SCIPrationalMult(deltaval, lb, oldobj);
            SCIPrationalNegate(deltaval, deltaval);
         }
      }
   }
   else if( SCIPrationalIsNegative(oldobj) )
   {
      /* sign of objective did not change */
      if( SCIPrationalIsNegative(newobj) )
      {
         /* if the bound is finite, calculate the deltaval */
         if( !SCIPrationalIsInfinity(ub) )
         {
            SCIPrationalDiff(tmp, newobj, oldobj);
            SCIPrationalMult(deltaval, ub, tmp);
         }
      }
      /* sign of objective did change, so the best bound does change */
      else if( SCIPrationalIsPositive(newobj) )
      {
         if( SCIPrationalIsInfinity(ub) )
         {
            /* old best bound was infinite while new one is not */
            if( !SCIPrationalIsNegInfinity(lb) )
            {
               (*deltainf) = -1;
               SCIPrationalMult(deltaval, lb, newobj);
            }
         }
         else
         {
            /* new best bound is infinite while old one was not */
            if( SCIPrationalIsNegInfinity(lb) )
            {
               (*deltainf) = 1;
               SCIPrationalMult(deltaval, ub, oldobj);
               SCIPrationalNegate(deltaval, deltaval);
            }
            /* neither old nor new best bound is infinite, so just calculate the deltaval */
            else
            {
               SCIPrationalMult(tmp, ub, oldobj);
               SCIPrationalMult(deltaval, lb, newobj);
               SCIPrationalDiff(deltaval, deltaval, tmp);
            }
         }
      }
      /* new objective is 0.0 */
      else
      {
         if( SCIPrationalIsInfinity(ub) )
            (*deltainf) = -1;
         else
         {
            SCIPrationalMult(deltaval, ub, oldobj);
            SCIPrationalNegate(deltaval, deltaval);
         }
      }
   }
   /* old objective was 0.0 */
   else
   {
      if( SCIPrationalIsNegative(newobj) )
      {
         if( SCIPrationalIsInfinity(ub) )
            (*deltainf) = 1;
         else
            SCIPrationalMult(deltaval, ub, newobj);
      }
      else if( SCIPrationalIsPositive(newobj) )
      {
         if( SCIPrationalIsNegInfinity(lb) )
            (*deltainf) = 1;
         else
            SCIPrationalMult(deltaval, lb, newobj);
      }
   }

   SCIPrationalFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** returns the left hand side of the row */
SCIP_RATIONAL* SCIProwExactGetLhs(
   SCIP_ROWEXACT*        row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->lhs != NULL);

   return row->lhs;
}

/** returns the right hand side of the row */
SCIP_RATIONAL* SCIProwExactGetRhs(
   SCIP_ROWEXACT*        row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->rhs != NULL);

   return row->rhs;
}

/** returns the constant of the row */
SCIP_RATIONAL* SCIProwExactGetConstant(
   SCIP_ROWEXACT*        row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->constant != NULL);

   return row->constant;
}

/** compute the objective delta due the new lower bound */
static
void getObjvalDeltaLbExact(
   SCIP_RATIONAL*        obj,                /**< objective value of variable */
   SCIP_RATIONAL*        oldlb,              /**< old lower bound of variable */
   SCIP_RATIONAL*        newlb,              /**< new lower bound of variable */
   SCIP_RATIONAL*        deltaval,           /**< pointer to store the delta value */
   int*                  deltainf            /**< pointer to store the number of variables with infinite best bound */
   )
{
   assert(!SCIPrationalIsAbsInfinity(obj));
   assert(!SCIPrationalIsInfinity(oldlb));
   assert(!SCIPrationalIsNegInfinity(oldlb) || !SCIPrationalIsNegInfinity(newlb));
   assert(SCIPrationalIsPositive(obj)); /* we only need to update if the objective is positive */

   if( SCIPrationalIsNegInfinity(oldlb) )
   {
      if( !SCIPrationalIsInfinity(newlb) )
      {
         (*deltainf) = -1;
         SCIPrationalMult(deltaval, newlb, obj);
      }
      else
      {
         (*deltainf) = 0;
         SCIPrationalSetReal(deltaval, 0.0);
      }
   }
   else if( SCIPrationalIsAbsInfinity(newlb) )
   {
      (*deltainf) = 1;
      SCIPrationalMult(deltaval, oldlb, obj);
      SCIPrationalNegate(deltaval, deltaval);
   }
   else
   {
      (*deltainf) = 0;
      SCIPrationalDiff(deltaval, newlb, oldlb);
      SCIPrationalMult(deltaval, deltaval, obj);
   }
}

/** compute the objective delta due the new upper bound */
static
void getObjvalDeltaUbExact(
   SCIP_RATIONAL*        obj,                /**< objective value of variable */
   SCIP_RATIONAL*        oldub,              /**< old upper bound of variable */
   SCIP_RATIONAL*        newub,              /**< new upper bound of variable */
   SCIP_RATIONAL*        deltaval,           /**< pointer to store the delta value */
   int*                  deltainf            /**< pointer to store the number of variables with infinite best bound */
   )
{
   assert(!SCIPrationalIsAbsInfinity(obj));
   assert(!SCIPrationalIsNegInfinity(oldub));
   assert(!SCIPrationalIsInfinity(oldub) || !SCIPrationalIsInfinity(newub));
   assert(SCIPrationalIsNegative(obj)); /* we only need to update if the objective is negative */

   if( SCIPrationalIsInfinity(oldub) )
   {
      if( !SCIPrationalIsNegInfinity(newub) )
      {
         (*deltainf) = -1;
         SCIPrationalMult(deltaval, newub, obj);
      }
      else
      {
         (*deltainf) = 0;
         SCIPrationalSetReal(deltaval, 0.0);
      }
   }
   else if( SCIPrationalIsAbsInfinity(newub) )
   {
      (*deltainf) = 1;
      SCIPrationalMult(deltaval, oldub, obj);
      SCIPrationalNegate(deltaval, deltaval);
   }
   else
   {
      (*deltainf) = 0;
      SCIPrationalDiff(deltaval, newub, oldub);
      SCIPrationalMult(deltaval, deltaval, obj);
   }
}

/** updates current pseudo and loose objective values for a change in a variable's objective value or bounds */
static
void lpExactUpdateObjval(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_RATIONAL*        deltavalex,         /**< delta value in the objective function */
   int                   deltainf,           /**< delta value for the number of variables with infinite best bound */
   SCIP_Bool             local,              /**< should the local pseudo objective value be updated? */
   SCIP_Bool             loose,              /**< should the loose objective value be updated? */
   SCIP_Bool             global              /**< should the global pseudo objective value be updated? */
   )
{
   assert(lpexact != NULL);
   assert(lpexact->looseobjvalinf >= 0);
   assert(lpexact->pseudoobjvalinf >= 0);
   assert(lpexact->glbpseudoobjvalinf >= 0);

   /* update the pseudo objective value */
   if( local )
   {
      lpexact->pseudoobjvalinf += deltainf;

      SCIPrationalAdd(lpexact->pseudoobjval, lpexact->pseudoobjval, deltavalex);

      /* after changing a local bound on a LOOSE variable, we have to update the loose objective value, too */
      if( SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE )
         loose = TRUE;
   }
   /* update the loose objective value */
   if( loose )
   {
      lpexact->looseobjvalinf += deltainf;

      if( !SCIPrationalIsZero(deltavalex) )
         SCIPrationalAdd(lpexact->looseobjval, lpexact->looseobjval, deltavalex);
   }

   /* update the root pseudo objective values */
   if( global )
   {
      lpexact->glbpseudoobjvalinf += deltainf;

      SCIPrationalAdd(lpexact->glbpseudoobjval ,lpexact->glbpseudoobjval, deltavalex);
   }

   assert(lpexact->looseobjvalinf >= 0);
   assert(lpexact->pseudoobjvalinf >= 0);
   assert(lpexact->glbpseudoobjvalinf >= 0);
}

/** updates current pseudo and loose objective value for a change in a variable's objective value */
SCIP_RETCODE SCIPlpExactUpdateVarObj(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_RATIONAL*        oldobj,             /**< old objective value of variable */
   SCIP_RATIONAL*        newobj              /**< new objective value of variable */
   )
{
   assert(lpexact != NULL);
   assert(var != NULL);

   if( !SCIPrationalIsEQ(oldobj, newobj) )
   {
      SCIP_RATIONAL* deltaval;
      int deltainf = 0;

      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetProbindex(var) >= 0);
      /* the objective coefficient can only be changed during presolving, that implies that the global and local
       * domain of the variable are the same
       */
      assert(lpexact->fplp->probing || SCIPrationalIsEQ(SCIPvarGetLbGlobalExact(var), SCIPvarGetLbLocalExact(var)));
      assert(lpexact->fplp->probing || SCIPrationalIsEQ(SCIPvarGetUbGlobalExact(var), SCIPvarGetUbLocalExact(var)));

      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &deltaval) );

      /* compute the pseudo objective delta due the new objective coefficient */
      SCIP_CALL( getObjvalDeltaObjExact(set, oldobj, newobj, SCIPvarGetLbLocalExact(var),
            SCIPvarGetUbLocalExact(var), deltaval, &deltainf) );

      /* update the local pseudo objective value */
      lpExactUpdateObjval(lpexact, var, deltaval, deltainf, TRUE, FALSE, FALSE);

      /* compute the pseudo objective delta due the new objective coefficient */
      SCIP_CALL( getObjvalDeltaObjExact(set, oldobj, newobj, SCIPvarGetLbGlobalExact(var),
            SCIPvarGetUbGlobalExact(var), deltaval, &deltainf) );

      /* update the global pseudo objective value */
      lpExactUpdateObjval(lpexact, var, deltaval, deltainf, FALSE, FALSE, TRUE);

      SCIPrationalFreeBuffer(set->buffer, &deltaval);
   }

   return SCIP_OKAY;
}

/** updates current root pseudo objective value for a global change in a variable's lower bound */
SCIP_RETCODE SCIPlpExactUpdateVarLbGlobal(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_RATIONAL*        oldlb,              /**< old lower bound of variable */
   SCIP_RATIONAL*        newlb               /**< new lower bound of variable */
   )
{
   assert(lpexact != NULL);
   assert(set != NULL);
   assert(var != NULL);

   if( !SCIPrationalIsEQ(oldlb, newlb) && SCIPrationalIsPositive(SCIPvarGetObjExact(var)) )
   {
      SCIP_RATIONAL* deltaval;
      int deltainf;

      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &deltaval) );

      /* compute the pseudo objective delta due the new lower bound */
      getObjvalDeltaLbExact(SCIPvarGetObjExact(var), oldlb, newlb, deltaval, &deltainf);

      /* update the root pseudo objective values */
      lpExactUpdateObjval(lpexact, var, deltaval, deltainf, FALSE, FALSE, TRUE);

      SCIPrationalFreeBuffer(set->buffer, &deltaval);
   }

   return SCIP_OKAY;
}

/** updates current pseudo and loose objective value for a change in a variable's lower bound */
SCIP_RETCODE SCIPlpExactUpdateVarLb(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_RATIONAL*        oldlb,              /**< old lower bound of variable */
   SCIP_RATIONAL*        newlb               /**< new lower bound of variable */
   )
{
   assert(lpexact != NULL);
   assert(set != NULL);
   assert(var != NULL);

   if( !SCIPrationalIsEQ(oldlb, newlb) && SCIPrationalIsPositive(SCIPvarGetObjExact(var)) )
   {
      SCIP_RATIONAL* deltaval;
      int deltainf;

      assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetProbindex(var) >= 0);

      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &deltaval) );

      /* compute the pseudo objective delta due the new lower bound */
      getObjvalDeltaLbExact(SCIPvarGetObjExact(var), oldlb, newlb, deltaval, &deltainf);

      /* update the pseudo and loose objective values */
      lpExactUpdateObjval(lpexact, var, deltaval, deltainf, TRUE, FALSE, FALSE);

      SCIPrationalFreeBuffer(set->buffer, &deltaval);
   }

   return SCIP_OKAY;
}

/** updates current root pseudo objective value for a global change in a variable's upper bound */
SCIP_RETCODE SCIPlpExactUpdateVarUbGlobal(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_RATIONAL*        oldub,              /**< old upper bound of variable */
   SCIP_RATIONAL*        newub               /**< new upper bound of variable */
   )
{
   assert(lpexact != NULL);
   assert(set != NULL);
   assert(var != NULL);

   if( !SCIPrationalIsEQ(oldub, newub) && SCIPrationalIsNegative(SCIPvarGetObjExact(var)) )
   {
      SCIP_RATIONAL* deltaval;
      int deltainf;

      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &deltaval) );

      /* compute the pseudo objective delta due the new lower bound */
      getObjvalDeltaUbExact(SCIPvarGetObjExact(var), oldub, newub, deltaval, &deltainf);

      /* update the root pseudo objective values */
      lpExactUpdateObjval(lpexact, var, deltaval, deltainf, FALSE, FALSE, TRUE);

      SCIPrationalFreeBuffer(set->buffer, &deltaval);
   }

   return SCIP_OKAY;
}

/** updates current pseudo objective value for a change in a variable's upper bound */
SCIP_RETCODE SCIPlpExactUpdateVarUb(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_RATIONAL*        oldub,              /**< old upper bound of variable */
   SCIP_RATIONAL*        newub               /**< new upper bound of variable */
   )
{
   assert(lpexact != NULL);
   assert(set != NULL);
   assert(var != NULL);

   if( !SCIPrationalIsEQ(oldub, newub) && SCIPrationalIsNegative(SCIPvarGetObjExact(var)) )
   {
      SCIP_RATIONAL* deltaval;
      int deltainf;

      assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetProbindex(var) >= 0);

      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &deltaval) );

      /* compute the pseudo objective delta due the new lower bound */
      getObjvalDeltaUbExact(SCIPvarGetObjExact(var), oldub, newub, deltaval, &deltainf);

      /* update the pseudo and loose objective values */
      lpExactUpdateObjval(lpexact, var, deltaval, deltainf, TRUE, FALSE, FALSE);

      SCIPrationalFreeBuffer(set->buffer, &deltaval);
   }

   return SCIP_OKAY;
}

/** informs LP, that given variable was added to the problem */
SCIP_RETCODE SCIPlpExactUpdateAddVar(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that is now a LOOSE problem variable */
   )
{
   SCIP_RATIONAL* tmp;

   if( !set->exact_enable )
      return SCIP_OKAY;

   assert(lpexact != NULL);
   assert(set != NULL);
   assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );

   /* add the variable to the loose objective value sum */
   SCIP_CALL( SCIPlpExactUpdateVarObj(set, lpexact, var, tmp, SCIPvarGetObjExact(var)) );

   /* update the loose variables counter */
   if( SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE )
      lpexact->nloosevars++;

   SCIPrationalFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** informs LP, that given variable is to be deleted from the problem */
SCIP_RETCODE SCIPlpExactUpdateDelVar(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that will be deleted from the problem */
   )
{
   SCIP_RATIONAL* ratzero;

   assert(lpexact != NULL);
   assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &ratzero) );

   /* subtract the variable from the loose objective value sum */
   SCIP_CALL( SCIPlpExactUpdateVarObj(set, lpexact, var, SCIPvarGetObjExact(var), ratzero) );

   /* update the loose variables counter */
   if( SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE )
   {
      SCIPlpExactDecNLoosevars(lpexact);
   }

   SCIPrationalFreeBuffer(set->buffer, &ratzero);

   return SCIP_OKAY;
}

/** informs LP, that given formerly loose problem variable is now a column variable */
SCIP_RETCODE SCIPlpExactUpdateVarColumn(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   )
{
   SCIP_RATIONAL* tmp;
   SCIP_RATIONAL* obj;
   SCIP_RATIONAL* lb;
   SCIP_RATIONAL* ub;

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );

   assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(lpexact->looseobjvalinf >= 0);

   obj = SCIPvarGetObjExact(var);

   /* update loose objective value */
   if( SCIPrationalIsPositive(obj) )
   {
      lb = SCIPvarGetLbLocalExact(var);
      if( SCIPrationalIsNegInfinity(lb) )
         lpexact->looseobjvalinf--;
      else
      {
         SCIPrationalNegate(tmp, lb);
         SCIPrationalMult(tmp, tmp, obj);
         lpExactUpdateObjval(lpexact, var, tmp, 0, FALSE, TRUE, FALSE);
      }
   }
   else if( SCIPrationalIsNegative(obj) )
   {
      ub = SCIPvarGetUbLocalExact(var);
      if( SCIPrationalIsInfinity(ub) )
         lpexact->looseobjvalinf--;
      else
      {
         SCIPrationalNegate(tmp, ub);
         SCIPrationalMult(tmp, tmp, obj);
         lpExactUpdateObjval(lpexact, var, tmp, 0, FALSE, TRUE, FALSE);
      }
   }

   SCIPlpExactDecNLoosevars(lpexact);

   assert(lpexact->looseobjvalinf >= 0);

   SCIPrationalFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** informs LP, that given formerly column problem variable is now again a loose variable */
SCIP_RETCODE SCIPlpExactUpdateVarLoose(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   )
{
   SCIP_RATIONAL* tmp;
   SCIP_RATIONAL* obj;
   SCIP_RATIONAL* lb;
   SCIP_RATIONAL* ub;

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );

   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(lpexact->looseobjvalinf >= 0);

   obj = SCIPvarGetObjExact(var);

   /* update loose objective value corresponding to the addition of variable */
   if( SCIPrationalIsPositive(obj) )
   {
      lb = SCIPvarGetLbLocalExact(var);
      if( SCIPrationalIsNegInfinity(lb) )
         lpexact->looseobjvalinf++;
      else
      {
         SCIPrationalMult(tmp, lb, obj);
         lpExactUpdateObjval(lpexact, var, tmp, 0, FALSE, TRUE, FALSE);
      }
   }
   else if( SCIPrationalIsNegative(obj) )
   {
      ub = SCIPvarGetUbLocalExact(var);
      if( SCIPrationalIsInfinity(ub) )
         lpexact->looseobjvalinf++;
      else
      {
         SCIPrationalMult(tmp, ub, obj);
         lpExactUpdateObjval(lpexact, var, tmp, 0, FALSE, TRUE, FALSE);
      }
   }
   lpexact->nloosevars++;

   assert(lpexact->looseobjvalinf >= 0);

   SCIPrationalFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** decrease the number of loose variables by one */
void SCIPlpExactDecNLoosevars(
   SCIP_LPEXACT*         lpexact             /**< current LP data */
   )
{
   assert(lpexact != NULL);
   assert(lpexact->nloosevars > 0);

   lpexact->nloosevars--;

   /* get rid of numerical problems: set loose objective value explicitly to zero, if no loose variables remain */
   if( lpexact->nloosevars == 0 )
   {
      assert(lpexact->looseobjvalinf == 0);
      SCIPrationalSetReal(lpexact->looseobjval, 0.0);
   }
}

/** get the number of rows currently in the lp */
int SCIPlpExactGetNRows(
   SCIP_LPEXACT*         lpexact             /**< current LP data */
   )
{
   assert(lpexact != NULL);

   return lpexact->nrows;
}

#ifdef SCIP_DEBUG
static
SCIP_RETCODE lpexactComputeDualValidity(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RATIONAL**       dualsol,            /**< row dual multipliers */
   SCIP_RATIONAL**       redcost             /**< column reduced costs */
   )
{
   int r,c;
   SCIP_RATIONAL** obj;
   SCIP_RATIONAL* objval;

   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &obj, lpexact->ncols) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &objval) );

   for( c = 0; c < lpexact->nlpicols; c++ )
   {
      SCIPrationalSetRational(obj[c], lpexact->cols[c]->obj);
      SCIPrationalDiff(obj[c], obj[c], redcost[c]);

      if( SCIPrationalIsPositive(redcost[c]) )
         SCIPrationalDiffProd(objval, redcost[c], lpexact->cols[c]->lb);
      else if( SCIPrationalIsNegative(redcost[c]) )
         SCIPrationalAddProd(objval, redcost[c], lpexact->cols[c]->ub);
   }

   for( r = 0; r < lpexact->nlpirows; r++ )
   {
      SCIP_ROWEXACT* row = lpexact->lpirows[r];

      if( SCIPrationalIsPositive(dualsol[r]) )
         SCIPrationalDiffProd(objval, dualsol[r], row->lhs);
      else if( SCIPrationalIsNegative(dualsol[r]) )
         SCIPrationalAddProd(objval, dualsol[r], row->rhs);

      for( c = 0; c < row->len; c++ )
      {
         int idx = row->cols_index[c];
         SCIPrationalDiffProd(obj[idx], row->vals[c], dualsol[r]);
      }
   }

   for( c = 0; c < lpexact->ncols; c++ )
   {
      assert(SCIPrationalIsZero(obj[c]));
   }

   SCIPrationalFreeBuffer(set->buffer, &objval);
   SCIPrationalFreeBufferArray(set->buffer, &obj, lpexact->ncols);

   return SCIP_OKAY;
}
#endif

/** stores the LP solution in the columns and rows */
SCIP_RETCODE SCIPlpExactGetSol(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            dualfeasible,       /**< pointer to store whether the solution is dual feasible, or NULL */
   SCIP_Bool             overwritefplp       /**< should the floating point values be overwritten, e.g. if fp lp was infeasible */
   )
{
   SCIP_COLEXACT** lpicols;
   SCIP_ROWEXACT** lpirows;
   SCIP_RATIONAL** primsol;
   SCIP_RATIONAL** dualsol;
   SCIP_RATIONAL** activity;
   SCIP_RATIONAL** redcost;
   SCIP_RATIONAL* primalbound;
   SCIP_RATIONAL* dualbound;
   SCIP_RATIONAL* tmp;
   SCIP_Bool stillprimalfeasible;
   SCIP_Bool stilldualfeasible;
   int* cstat;
   int* rstat;
   SCIP_Longint lpcount;
   int nlpicols;
   int nlpirows;
   int c;
   int r;

   assert(lpexact != NULL);
   assert(lpexact->solved);
   assert(set != NULL);
   assert(stat != NULL);

   /* initialize return and feasibility flags; if primal oder dual feasibility shall not be checked, we set the
    * corresponding flag immediately to FALSE to skip all checks
    */
   if( primalfeasible == NULL )
      stillprimalfeasible = FALSE;
   else
   {
      *primalfeasible = TRUE;
      stillprimalfeasible = TRUE;
   }
   if( dualfeasible == NULL )
      stilldualfeasible = FALSE;
   else
   {
      *dualfeasible = TRUE;
      stilldualfeasible = TRUE;
   }

   SCIPsetDebugMsg(set, "getting new LP solution %" SCIP_LONGINT_FORMAT " for solstat %d\n",
      stat->lpcount, lpexact->lpsolstat);

   lpicols = lpexact->lpicols;
   lpirows = lpexact->lpirows;
   nlpicols = lpexact->nlpicols;
   nlpirows = lpexact->nlpirows;
   lpcount = stat->lpcount;

   /* get temporary memory */
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &primalbound) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &dualbound) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &primsol, nlpicols) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &dualsol, nlpirows) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &activity, nlpirows) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &redcost, nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cstat, nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rstat, nlpirows) );

   SCIP_CALL( SCIPlpiExactGetSol(lpexact->lpiexact, NULL, primsol, dualsol, activity, redcost) );

   /* avoid adding infinity to the bounding error */
   if( !SCIPrationalIsInfinity(lpexact->lpobjval) )
      stat->boundingerrorexlp += REALABS(lpexact->fplp->lpobjval - SCIPrationalRoundReal(lpexact->lpobjval, SCIP_R_ROUND_DOWNWARDS));
   if( overwritefplp )
   {
      lpexact->fplp->lpobjval = SCIPrationalRoundReal(lpexact->lpobjval, SCIP_R_ROUND_DOWNWARDS);
      lpexact->fplp->lpsolstat = lpexact->lpsolstat;
      lpexact->fplp->primalfeasible = lpexact->primalfeasible;
      lpexact->fplp->dualfeasible = lpexact->dualfeasible;
      lpexact->fplp->solved = lpexact->solved;
   }
   if( lpexact->solisbasic )
   {
      SCIP_CALL( SCIPlpiExactGetBase(lpexact->lpiexact, cstat, rstat) );
   }
   else
   {
      BMSclearMemoryArray(cstat, nlpicols);
      BMSclearMemoryArray(rstat, nlpirows);
   }

   SCIPrationalSetReal(primalbound, 0.0);
   SCIPrationalSetReal(dualbound, 0.0);

   SCIPdebug(SCIP_CALL( lpexactComputeDualValidity(lpexact, set, dualsol, redcost) ));

   /* copy primal solution and reduced costs into columns */
   for( c = 0; c < nlpicols; ++c )
   {
      assert( 0 <= cstat[c] && cstat[c] < 4 );
      SCIPrationalSetRational(lpicols[c]->primsol, primsol[c]);
      SCIPrationalSetRational(lpicols[c]->redcost, redcost[c]);
      lpicols[c]->basisstatus = (unsigned int) cstat[c];
      lpicols[c]->validredcostlp = lpcount;
      if( overwritefplp )
      {
         lpexact->fplp->lpicols[c]->primsol =  SCIPrationalGetReal(primsol[c]);
         lpexact->fplp->lpicols[c]->redcost =  SCIPrationalGetReal(redcost[c]);
         lpexact->fplp->lpicols[c]->basisstatus = (unsigned int) cstat[c];
         lpexact->fplp->lpicols[c]->validredcostlp = lpcount;
      }
      if( stillprimalfeasible )
      {
         stillprimalfeasible =
            (SCIPrationalIsNegInfinity(lpicols[c]->lb) || !SCIPrationalIsLT(lpicols[c]->primsol, lpicols[c]->lb))
            && (SCIPrationalIsInfinity(lpicols[c]->ub) || !SCIPrationalIsGT(lpicols[c]->primsol, lpicols[c]->ub));
         SCIPrationalAddProd(primalbound, lpicols[c]->primsol, lpicols[c]->obj);
      }

      /* complementary slackness means that if a variable is not at its lower or upper bound, its reduced costs
       * must be non-positive or non-negative, respectively; in particular, if a variable is strictly within its
       * bounds, its reduced cost must be zero
       */
      if( stilldualfeasible && (SCIPrationalIsNegInfinity(lpicols[c]->lb) || SCIPrationalIsGT(lpicols[c]->primsol, lpicols[c]->lb)) )
         stilldualfeasible = !SCIPrationalIsPositive(lpicols[c]->redcost);
      if( stilldualfeasible && (SCIPrationalIsInfinity(lpicols[c]->ub) || SCIPrationalIsLT(lpicols[c]->primsol, lpicols[c]->ub)) )
         stilldualfeasible = !SCIPrationalIsNegative(lpicols[c]->redcost);

      SCIPrationalDebugMessage("col <%s> [%q,%q]: primsol=%q, redcost=%q, pfeas=%u/%u(%u), dfeas=%d/%d(%u)\n",
         SCIPvarGetName(lpicols[c]->var), lpicols[c]->lb, lpicols[c]->ub, lpicols[c]->primsol, lpicols[c]->redcost,
         SCIPrationalIsGE(lpicols[c]->primsol, lpicols[c]->lb),
         SCIPrationalIsLE(lpicols[c]->primsol, lpicols[c]->ub),
         primalfeasible != NULL ? stillprimalfeasible : TRUE,
         !SCIPrationalIsGT(lpicols[c]->primsol, lpicols[c]->lb) || !SCIPrationalIsPositive(lpicols[c]->redcost),
         !SCIPrationalIsGT(lpicols[c]->primsol, lpicols[c]->ub) || !SCIPrationalIsNegative(lpicols[c]->redcost),
         dualfeasible != NULL ? stilldualfeasible : TRUE);

      /* we intentionally use an exact positive/negative check because ignoring small reduced cost values may lead to a
       * wrong bound value; if the corresponding bound is +/-infinity, we use zero reduced cost (if stilldualfeasible is
       * TRUE, we are in the case that the reduced cost is tiny with wrong sign)
       */
      if( stilldualfeasible )
      {
         if( SCIPrationalIsPositive(lpicols[c]->redcost) && !SCIPrationalIsNegInfinity(lpicols[c]->lb) )
         {
            SCIPrationalAddProd(dualbound, lpicols[c]->redcost, lpicols[c]->lb);
         }
         else if( SCIPrationalIsNegative(lpicols[c]->redcost) && !SCIPrationalIsInfinity(lpicols[c]->ub) )
         {
            SCIPrationalAddProd(dualbound, lpicols[c]->redcost, lpicols[c]->ub);
         }
      }
   }

   /* copy dual solution and activities into rows */
   for( r = 0; r < nlpirows; ++r )
   {
      assert( 0 <= rstat[r] && rstat[r] < 4 );
      SCIPrationalSetRational(lpirows[r]->dualsol, dualsol[r]);
      SCIPrationalAdd(lpirows[r]->activity, activity[r], lpirows[r]->constant);
      lpirows[r]->basisstatus = (unsigned int) rstat[r]; /*lint !e732*/
      lpirows[r]->validactivitylp = lpcount;
      if( overwritefplp )
      {
         SCIP_ROW* fprow;
         if( SCIProwIsInLP(lpirows[r]->fprow) )
            fprow = lpirows[r]->fprow;
         else
         {
            assert(SCIProwIsInLP(lpirows[r]->fprowrhs));
            fprow = lpirows[r]->fprowrhs;
         }
         fprow->dualsol = SCIPrationalGetReal(dualsol[r]);
         fprow->activity = SCIPrationalGetReal(lpirows[r]->activity);
         fprow->basisstatus = (unsigned int) rstat[r]; /*lint !e732*/
         fprow->validactivitylp = lpcount;
      }
      if( stillprimalfeasible )
      {
         stillprimalfeasible =
            (SCIPrationalIsNegInfinity(lpirows[r]->lhs) ||SCIPrationalIsGE(lpirows[r]->activity, lpirows[r]->lhs))
            && (SCIPrationalIsInfinity(lpirows[r]->rhs) || SCIPrationalIsLE(lpirows[r]->activity, lpirows[r]->rhs));
      }
      /* complementary slackness means that if the activity of a row is not at its left-hand or right-hand side,
       * its dual multiplier must be non-positive or non-negative, respectively; in particular, if the activity is
       * strictly within left-hand and right-hand side, its dual multiplier must be zero
       */
      if( stilldualfeasible &&
            (SCIPrationalIsNegInfinity(lpirows[r]->lhs) || SCIPrationalIsGT(lpirows[r]->activity, lpirows[r]->lhs)) )
         stilldualfeasible = !SCIPrationalIsPositive(lpirows[r]->dualsol);
      if( stilldualfeasible &&
            (SCIPrationalIsInfinity(lpirows[r]->rhs) || SCIPrationalIsLT(lpirows[r]->activity, lpirows[r]->rhs)) )
         stilldualfeasible = !SCIPrationalIsNegative(lpirows[r]->dualsol);

      SCIPrationalDebugMessage("<%s> [%q,%q] + %q: activity=%q, dualsol=%q, pfeas=%u/%u(%u), dfeas=%d/%d(%u)\n",
         lpirows[r]->fprow->name, lpirows[r]->lhs, lpirows[r]->rhs,
         lpirows[r]->constant, lpirows[r]->activity, lpirows[r]->dualsol,
         SCIPrationalIsGE(lpirows[r]->activity, lpirows[r]->lhs),
         SCIPrationalIsLE(lpirows[r]->activity, lpirows[r]->rhs),
         primalfeasible != NULL ? stillprimalfeasible : TRUE,
         !SCIPrationalIsGT(lpirows[r]->activity, lpirows[r]->lhs) || !SCIPrationalIsPositive(lpirows[r]->dualsol),
         !SCIPrationalIsLT(lpirows[r]->activity, lpirows[r]->rhs) || !SCIPrationalIsNegative(lpirows[r]->dualsol),
         dualfeasible != NULL ? stilldualfeasible : TRUE);

      /* we intentionally use an exact positive/negative check because ignoring small dual multipliers may lead to a
       * wrong bound value; if the corresponding side is +/-infinity, we use a zero dual multiplier (if
       * stilldualfeasible is TRUE, we are in the case that the dual multiplier is tiny with wrong sign)
       */
      if( stilldualfeasible )
      {
         if( SCIPrationalIsPositive(lpirows[r]->dualsol) && !SCIPrationalIsNegInfinity(lpirows[r]->lhs) )
         {
            SCIPrationalDiff(tmp, lpirows[r]->lhs, lpirows[r]->constant);
            SCIPrationalAddProd(dualbound, tmp, lpirows[r]->dualsol);
         }
         else if( SCIPrationalIsNegative(lpirows[r]->dualsol) && !SCIPrationalIsInfinity(lpirows[r]->rhs) )
         {
            SCIPrationalDiff(tmp, lpirows[r]->rhs, lpirows[r]->constant);
            SCIPrationalAddProd(dualbound, tmp, lpirows[r]->dualsol);
         }
      }
   }

   /* if the objective value returned by the LP solver is smaller than the internally computed primal bound, then we
    * declare the solution primal infeasible; we assume primalbound and lpexact->lpobjval to be equal if they are both +/-
    * infinity
    */
   /**@todo alternatively, if otherwise the LP solution is feasible, we could simply update the objective value */
   if( stillprimalfeasible && !(SCIPrationalIsInfinity(primalbound) && SCIPrationalIsInfinity(lpexact->lpobjval))
      && !(SCIPrationalIsNegInfinity(primalbound) && SCIPrationalIsNegInfinity(lpexact->lpobjval)) )
   {
      stillprimalfeasible = SCIPrationalIsLE(primalbound, lpexact->lpobjval);
      SCIPrationalDebugMessage(" primalbound=%q, lpbound=%q, pfeas=%u(%u)\n", primalbound, lpexact->lpobjval,
         SCIPrationalIsLE(primalbound, lpexact->lpobjval), primalfeasible != NULL ? stillprimalfeasible : TRUE);
   }

   /* if the objective value returned by the LP solver is smaller than the internally computed dual bound, we declare
    * the solution dual infeasible; we assume dualbound and lpexact->lpobjval to be equal if they are both +/- infinity
    */
   /**@todo alternatively, if otherwise the LP solution is feasible, we could simply update the objective value */
   if( stilldualfeasible && !(SCIPrationalIsInfinity(dualbound) && SCIPrationalIsInfinity(lpexact->lpobjval))
      && !(SCIPrationalIsNegInfinity(dualbound) && SCIPrationalIsNegInfinity(lpexact->lpobjval)) )
   {
      stilldualfeasible =  SCIPrationalIsGE(dualbound, lpexact->lpobjval);
      SCIPrationalDebugMessage(" dualbound=%q, lpbound=%q, dfeas=%u(%u)\n", dualbound, lpexact->lpobjval,
         SCIPrationalIsGE(dualbound, lpexact->lpobjval), dualfeasible != NULL ? stilldualfeasible : TRUE);
   }

   if( primalfeasible != NULL )
      *primalfeasible = stillprimalfeasible;
   if( dualfeasible != NULL )
      *dualfeasible = stilldualfeasible;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &rstat);
   SCIPsetFreeBufferArray(set, &cstat);
   SCIPrationalFreeBufferArray(set->buffer, &redcost, nlpicols);
   SCIPrationalFreeBufferArray(set->buffer, &activity, nlpirows);
   SCIPrationalFreeBufferArray(set->buffer, &dualsol, nlpirows);
   SCIPrationalFreeBufferArray(set->buffer, &primsol, nlpicols);
   SCIPrationalFreeBuffer(set->buffer, &tmp);
   SCIPrationalFreeBuffer(set->buffer, &dualbound);
   SCIPrationalFreeBuffer(set->buffer, &primalbound);

   return SCIP_OKAY;
}

/** stores LP solution with infinite objective value in the columns and rows */
SCIP_RETCODE SCIPlpExactGetUnboundedSol(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            rayfeasible         /**< pointer to store whether the primal ray is a feasible unboundedness proof, or NULL */
   )  /*lint --e{715}*/
{
   SCIPerrorMessage("Unbounded solution not implemented in exact solving mode.\n");
   return SCIP_ERROR;
}

/** returns primal ray proving the unboundedness of the current LP */
SCIP_RETCODE SCIPlpExactGetPrimalRay(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RATIONAL**       ray                 /**< array for storing primal ray values, they are stored w.r.t. the problem index of the variables,
                                              *   so the size of this array should be at least number of active variables
                                              *   (all entries have to be initialized to 0 before) */
   )
{
   SCIP_COLEXACT** lpicols;
   SCIP_RATIONAL** lpiray;
   SCIP_VAR* var;
   int nlpicols;
   int c;

   assert(lpexact != NULL);
   assert(set != NULL);
   assert(ray != NULL);
   assert(lpexact->flushed);
   assert(lpexact->solved);
   assert(lpexact->lpsolstat == SCIP_LPSOLSTAT_UNBOUNDEDRAY);
   assert(SCIPrationalIsNegInfinity(lpexact->lpobjval));

   /* check if the LP solver is able to provide a primal unbounded ray */
   if( !SCIPlpiExactHasPrimalRay(lpexact->lpiexact) )
   {
      SCIPerrorMessage("LP solver has no primal ray for unbounded LP\n");
      return SCIP_LPERROR;
   }

   /* get temporary memory */
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &lpiray, lpexact->nlpicols) );

   SCIPsetDebugMsg(set, "getting primal ray values\n");

   /* get primal unbounded ray */
   SCIP_CALL( SCIPlpiExactGetPrimalRay(lpexact->lpiexact, lpiray) );

   lpicols = lpexact->lpicols;
   nlpicols = lpexact->nlpicols;

   /* store the ray values of active problem variables */
   for( c = 0; c < nlpicols; c++ )
   {
      assert(lpicols[c] != NULL);

      var = lpicols[c]->var;
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) != -1);
      SCIPrationalSetRational(ray[SCIPvarGetProbindex(var)], lpiray[c]);
   }

   SCIPrationalFreeBufferArray(set->buffer, &lpiray, lpexact->nlpicols);

   return SCIP_OKAY;
}


/** stores the dual Farkas multipliers for infeasibility proof in rows. besides
 *
 *  @note The Farkas proof is checked for validity if lp/checkfarkas = TRUE and @p valid is not NULL.
 */
SCIP_RETCODE SCIPlpExactGetDualfarkas(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            valid,              /**< pointer to store whether the Farkas proof is valid or NULL */
   SCIP_Bool             overwritefplp       /**< should the floating point values be overwritten, e.g. if fp lp was infeasible */
   )
{
   SCIP_COLEXACT** lpicols;
   SCIP_ROWEXACT** lpirows;
   SCIP_RATIONAL** dualfarkas;
   SCIP_RATIONAL** farkascoefs;
   SCIP_RATIONAL* farkaslhs;
   SCIP_RATIONAL* maxactivity;
   SCIP_RATIONAL* tmp;
   SCIP_Bool checkfarkas;
   int nlpicols;
   int nlpirows;
   int c;
   int r;

   assert(lpexact != NULL);
   assert(lpexact->flushed);
   assert(lpexact->solved);
   assert(lpexact->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE);
   assert(set != NULL);
   assert(stat != NULL);

   if( valid != NULL )
      *valid = TRUE;

   farkascoefs = NULL;
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &maxactivity) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &farkaslhs) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );

   checkfarkas = (set->lp_checkfarkas && valid != NULL);

   /* get temporary memory */
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &dualfarkas, lpexact->nlpirows) );

   if( checkfarkas )
      SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &farkascoefs, lpexact->nlpicols) );

   /* get dual Farkas infeasibility proof */
   SCIP_CALL( SCIPlpiExactGetDualfarkas(lpexact->lpiexact, dualfarkas) );

   if( overwritefplp )
   {
      lpexact->fplp->lpobjval = SCIPsetInfinity(set);
      lpexact->fplp->lpsolstat = lpexact->lpsolstat;
   }

   lpicols = lpexact->lpicols;
   lpirows = lpexact->lpirows;
   nlpicols = lpexact->nlpicols;
   nlpirows = lpexact->nlpirows;

   /* store infeasibility proof in rows */
   SCIPsetDebugMsg(set, "LP is infeasible:\n");
   for( r = 0; r < nlpirows; ++r )
   {
      SCIPrationalDebugMessage(" row <%s>: dualfarkas=%q\n", lpirows[r]->fprow->name, dualfarkas[r]);
      SCIPrationalSetRational(lpirows[r]->dualfarkas, dualfarkas[r]);
      SCIPrationalSetInfinity(lpirows[r]->dualsol);
      SCIPrationalSetReal(lpirows[r]->activity, 0.0);
      lpirows[r]->validactivitylp = -1L;
      lpirows[r]->basisstatus = (unsigned int) SCIP_BASESTAT_BASIC;
      if( overwritefplp )
      {
         lpexact->fplp->lpirows[r]->dualfarkas = SCIPrationalGetReal(dualfarkas[r]);
         lpexact->fplp->lpirows[r]->dualsol = SCIPsetInfinity(set);
         lpexact->fplp->lpirows[r]->basisstatus = (unsigned int) SCIP_BASESTAT_BASIC;
         lpexact->fplp->lpirows[r]->validactivitylp = -1L;
      }

      if( checkfarkas )
      {
         assert(farkascoefs != NULL);

         /* the infeasibility proof would be invalid if
          *   (i)  dualfarkas[r] > 0 and lhs = -inf
          *   (ii) dualfarkas[r] < 0 and rhs = inf
          * however, due to numerics we accept slightly negative / positive values
          */
         if( (SCIPrationalIsPositive(dualfarkas[r]) && SCIPrationalIsNegInfinity(lpirows[r]->lhs))
            || (SCIPrationalIsNegative(dualfarkas[r]) && SCIPrationalIsInfinity(lpirows[r]->rhs)) )
         {
               SCIPrationalDebugMessage("farkas proof is invalid: row <%s>[lhs=%q,rhs=%q,c=%q] has multiplier %q\n",
               SCIProwGetName(lpirows[r]->fprow), lpirows[r]->lhs, lpirows[r]->rhs,
               lpirows[r]->constant, dualfarkas[r]);

               *valid = FALSE; /*lint --e{413,613}*/

            goto TERMINATE;
         }

         /* dual multipliers, for which the corresponding row side in infinite, are treated as zero if they are zero
          * within tolerances (see above) but slighty positive / negative
          */
         if( (SCIPrationalIsPositive(dualfarkas[r]) && SCIPrationalIsNegInfinity(lpirows[r]->lhs))
            || (SCIPrationalIsNegative(dualfarkas[r]) && SCIPrationalIsInfinity(lpirows[r]->rhs)) )
            continue;

         /* iterate over all columns and scale with dual solution */
         for( c = 0; c < lpirows[r]->len; c++ )
         {
            int pos = lpirows[r]->cols[c]->lppos;

            if( pos == -1 )
               continue;

            assert(pos >= 0 && pos < nlpicols);
            SCIPrationalAddProd(farkascoefs[pos], dualfarkas[r], lpirows[r]->vals[c]);
         }

         /* the row contributes with its left-hand side to the proof */
         if( SCIPrationalIsPositive(dualfarkas[r]) )
         {
            assert(!SCIPrationalIsNegInfinity(lpirows[r]->lhs));
            SCIPrationalDiff(tmp, lpirows[r]->lhs, lpirows[r]->constant);
            SCIPrationalAddProd(farkaslhs, tmp, dualfarkas[r]);
         }

         /* the row contributes with its right-hand side to the proof */
         else if( SCIPrationalIsNegative(dualfarkas[r]) )
         {
            assert(!SCIPrationalIsInfinity(lpirows[r]->rhs));
            SCIPrationalDiff(tmp, lpirows[r]->rhs, lpirows[r]->constant);
            SCIPrationalAddProd(farkaslhs, tmp, dualfarkas[r]);
         }
      }
   }

   /* set columns as invalid */
   for( c = 0; c < nlpicols; ++c )
   {
      SCIPrationalSetInfinity(lpicols[c]->primsol);
      SCIPrationalSetInfinity(lpicols[c]->redcost);
      lpicols[c]->validredcostlp = -1L;
      lpicols[c]->validfarkaslp = -1L;
      if( farkascoefs != NULL )
         SCIPrationalSetRational(lpicols[c]->farkascoef, farkascoefs[c]);

      if( overwritefplp )
      {
         lpexact->fplp->lpicols[c]->primsol =  SCIPsetInfinity(set);
         lpexact->fplp->lpicols[c]->redcost =  SCIPsetInfinity(set);
         lpexact->fplp->lpicols[c]->validredcostlp = -1L;
         lpexact->fplp->lpicols[c]->validfarkaslp = -1L;
      }

      if( checkfarkas )
      {
         assert(farkascoefs != NULL);
         assert(lpicols[c]->lppos == c);

         /* skip coefficients that are too close to zero */
         if( SCIPrationalIsZero(farkascoefs[c]) )
            continue;

         /* calculate the maximal activity */
         if( SCIPrationalIsPositive(farkascoefs[c]) )
         {
            SCIPrationalMult(tmp, farkascoefs[c], SCIPcolExactGetUb(lpicols[c]));
            SCIPrationalAdd(maxactivity, maxactivity, tmp);
         }
         else
         {
            SCIPrationalMult(tmp, farkascoefs[c], SCIPcolExactGetLb(lpicols[c]));
            SCIPrationalAdd(maxactivity, maxactivity, tmp);
         }
      }
   }

   /* check whether the farkasproof is valid
    * due to numerics, it might happen that the left-hand side of the aggregation is larger/smaller or equal than +/- infinity.
    * in that case, we declare the Farkas proof to be invalid.
    */
   if( checkfarkas && (SCIPrationalIsAbsInfinity(farkaslhs) || SCIPrationalIsGE(maxactivity, farkaslhs)) )
   {
      SCIPrationalDebugMessage("farkas proof is invalid: maxactivity=%q, lhs=%q\n", maxactivity, farkaslhs);

      *valid = FALSE; /*lint --e{413,613}*/
   }

  TERMINATE:
   /* free temporary memory */
   if( checkfarkas )
      SCIPrationalFreeBufferArray(set->buffer, &farkascoefs, nlpicols);

   SCIPrationalFreeBufferArray(set->buffer, &dualfarkas, nlpirows);
   SCIPrationalFreeBuffer(set->buffer, &tmp);
   SCIPrationalFreeBuffer(set->buffer, &farkaslhs);
   SCIPrationalFreeBuffer(set->buffer, &maxactivity);

   return SCIP_OKAY;
}

/** get number of iterations used in last LP solve */
SCIP_RETCODE SCIPlpExactGetIterations(
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   int*                  iterations          /**< pointer to store the iteration count */
   )
{
   assert(lpexact != NULL);

   SCIP_CALL( SCIPlpiExactGetIterations(lpexact->lpiexact, iterations) );

   return SCIP_OKAY;
}

/** gets objective value of current LP
 *
 *  @note This method returns the objective value of the current LP solution, which might be primal or dual infeasible
 *        if a limit was hit during solving. It must not be used as a dual bound if the LP solution status is
 *        SCIP_LPSOLSTAT_ITERLIMIT or SCIP_LPSOLSTAT_TIMELIMIT.
 */
void SCIPlpExactGetObjval(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RATIONAL*        res                 /**< result pointer to store rational */
   )
{
   assert(lpexact != NULL);
   assert(lpexact->fplp->hasprovedbound);
   assert((lpexact->nloosevars > 0) || (lpexact->looseobjvalinf == 0 && SCIPrationalIsZero(lpexact->looseobjval)));
   assert(set != NULL);

   if( lpexact->looseobjvalinf > 0 )
      SCIPrationalSetNegInfinity(res);
   else if( SCIPrationalIsAbsInfinity(lpexact->lpobjval) )
      SCIPrationalSetRational(res, lpexact->lpobjval);
   else
      SCIPrationalAdd(res, lpexact->lpobjval, lpexact->looseobjval);
}

/** gets the pseudo objective value for the current search node; that is all variables set to their best (w.r.t. the
 *  objective function) local bound
 */
void SCIPlpExactGetPseudoObjval(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RATIONAL*        res                 /**< result pointer to store rational */
   )
{
   assert(lpexact != NULL);
   assert(lpexact->pseudoobjvalinf >= 0);
   assert(set != NULL);

   if( lpexact->pseudoobjvalinf > 0 || set->nactivepricers > 0 )
      SCIPrationalSetNegInfinity(res);
   else
      SCIPrationalSetRational(res, lpexact->pseudoobjval);
}

/** removes all columns after the given number of cols from the LP */
SCIP_RETCODE SCIPlpExactShrinkCols(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   newncols            /**< new number of columns in the LP */
   )
{
   SCIP_COLEXACT* col;
   int c;

   assert(lpexact != NULL);

   SCIPsetDebugMsg(set, "shrinking LP from %d to %d columns\n", lpexact->ncols, newncols);
   assert(0 <= newncols);
   assert(newncols <= lpexact->ncols);

   if( newncols < lpexact->ncols )
   {
      assert(!lpexact->fplp->diving);

      for( c = lpexact->ncols-1; c >= newncols; --c )
      {
         col = lpexact->cols[c];
         assert(col != NULL);
         assert(col->len == 0 || col->rows != NULL);
         assert(col->var != NULL);
         assert(SCIPvarGetStatusExact(col->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetColExact(col->var) == lpexact->cols[c]);
         assert(col->lppos == c);

         /* mark column to be removed from the LP */
         col->lppos = -1;
         lpexact->ncols--;

         /* update column arrays of all linked rows */
         SCIP_CALL( colExactUpdateDelLP(col, set) );
      }

      assert(lpexact->ncols == newncols);
      lpexact->lpifirstchgcol = MIN(lpexact->lpifirstchgcol, newncols);

      lpexact->flushed = FALSE;
      checkLinks(lpexact);
   }

   return SCIP_OKAY;
}

/** removes and releases all rows after the given number of rows from the LP */
SCIP_RETCODE SCIPlpExactShrinkRows(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   newnrows            /**< new number of rows in the LP */
   )
{
   SCIP_ROWEXACT* row;
   int r;

   assert(lpexact != NULL);
   assert(0 <= newnrows && newnrows <= lpexact->nrows);

   SCIPsetDebugMsg(set, "shrinking exact LP from %d to %d rows\n", lpexact->nrows, newnrows);
   if( newnrows < lpexact->nrows )
   {
      for( r = lpexact->nrows-1; r >= newnrows; --r )
      {
         row = lpexact->rows[r];
         assert(row != NULL);
         assert(row->len == 0 || row->cols != NULL);
         assert(row->lppos == r);

         /* mark row to be removed from the LP */
         row->lppos = -1;
         row->lpdepth = -1;
         lpexact->nrows--;

         SCIP_CALL( rowExactUpdateDelLP(row, set) );

         SCIP_CALL( SCIProwExactRelease(&lpexact->rows[r], blkmem, set, lpexact) );
      }
      assert(lpexact->nrows == newnrows);
      lpexact->lpifirstchgrow = MIN(lpexact->lpifirstchgrow, newnrows);

      /* mark the current LP unflushed */
      lpexact->flushed = FALSE;

      checkLinks(lpexact);
   }

   return SCIP_OKAY;
}

/** resets the LP to the empty LP by removing all columns and rows from LP, releasing all rows, and flushing the
 *  changes to the LP solver
 */
SCIP_RETCODE SCIPlpExactReset(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   if( !set->exact_enable )
      return SCIP_OKAY;

   assert(stat != NULL);

   SCIP_CALL( SCIPlpExactClear(lpexact, blkmem, set) );
   SCIP_CALL( SCIPlpExactFlush(lpexact, blkmem, set, eventqueue) );

   /* mark the empty LP to be solved */
   lpexact->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   SCIPrationalSetReal(lpexact->lpobjval, 0.0);
   lpexact->solved = TRUE;
   lpexact->primalfeasible = TRUE;
   lpexact->primalchecked = TRUE;
   lpexact->dualfeasible = TRUE;
   lpexact->dualchecked = TRUE;
   lpexact->solisbasic = FALSE;
   lpexact->lastlpalgo = SCIP_LPALGO_DUALSIMPLEX;

   return SCIP_OKAY;
}

/** removes all columns and rows from LP, releases all rows */
SCIP_RETCODE SCIPlpExactClear(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(lpexact != NULL);
   assert(!lpexact->fplp->diving);

   SCIPsetDebugMsg(set, "clearing LP\n");
   SCIP_CALL( SCIPlpExactShrinkCols(lpexact, set, 0) );
   SCIP_CALL( SCIPlpExactShrinkRows(lpexact, blkmem, set, 0) );

   return SCIP_OKAY;
}

/** forces an exact lp to be solved in the next exact bound computation */
void SCIPlpExactForceExactSolve(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->exact_enable )
      return;

   assert(lpexact != NULL);

   lpexact->forceexactsolve = TRUE;
}

/** forces the next exact bound computation to be executed even in probing mode */
void SCIPlpExactForceSafeBound(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->exact_enable )
      return;

   assert(lpexact != NULL);

   lpexact->forcesafebound = TRUE;
}

/** allows an exact lp to be solved in the next exact bound computation */
void SCIPlpExactAllowExactSolve(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             allowexact          /**< TRUE if next safe bounding call should be allowed to be exact, FALSE otherwise */
   )
{
   assert(set != NULL);

   if( !set->exact_enable )
      return;

   assert(lpexact != NULL);

   lpexact->allowexactsolve = allowexact;
}

/** save current LP solution values stored in each column */
static
SCIP_RETCODE colExactStoreSolVals(
   SCIP_COLEXACT*        colexact,           /**< exact LP column */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_COLEXACTSOLVALS* storedsolvals;

   assert(colexact != NULL);
   assert(blkmem != NULL);

   /* allocate memory for storage */
   if( colexact->storedsolvals == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &colexact->storedsolvals) );

      storedsolvals = colexact->storedsolvals;

      /* store values */
      SCIP_CALL( SCIPrationalCopyBlock(blkmem, &(storedsolvals->primsol), colexact->primsol) );
      SCIP_CALL( SCIPrationalCopyBlock(blkmem, &(storedsolvals->redcost), colexact->redcost) );
      storedsolvals->basisstatus = colexact->basisstatus; /*lint !e641 !e732*/
   }
   else
   {
      storedsolvals = colexact->storedsolvals;

      /* store values */
      SCIPrationalSetRational(storedsolvals->primsol, colexact->primsol);
      SCIPrationalSetRational(storedsolvals->redcost, colexact->redcost);
      storedsolvals->basisstatus = colexact->basisstatus; /*lint !e641 !e732*/
   }

   return SCIP_OKAY;
}

/** restore LP solution values in column */
static
SCIP_RETCODE colExactRestoreSolVals(
   SCIP_COLEXACT*        colexact,           /**< exact LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Longint          validlp,            /**< number of lp for which restored values are valid */
   SCIP_Bool             freebuffer          /**< should buffer for LP solution values be freed? */
   )
{
   SCIP_COLEXACTSOLVALS* storedsolvals;

   assert(colexact != NULL);
   assert(blkmem != NULL);

   /* if stored values are available, restore them */
   storedsolvals = colexact->storedsolvals;
   if( storedsolvals != NULL )
   {
      SCIPrationalSetRational(colexact->primsol, storedsolvals->primsol);
      SCIPrationalSetRational(colexact->redcost, storedsolvals->redcost);
      colexact->validredcostlp = validlp;
      colexact->basisstatus = storedsolvals->basisstatus; /*lint !e641 !e732*/

      /* we do not save the farkas coefficient, since it can be recomputed; thus, we invalidate it here */
      colexact->validfarkaslp = -1;
   }
   /* if the column was created after performing the storage (possibly during probing), we treat it as implicitly zero;
    * we make sure to invalidate the reduced cost and farkas coefficient, which are not available
    */
   else
   {
      SCIPrationalSetReal(colexact->primsol, 0.0);
      colexact->validredcostlp = -1;
      colexact->validfarkaslp = -1;
      colexact->basisstatus = SCIP_BASESTAT_ZERO; /*lint !e641*/
   }

   /* free memory */
   if( freebuffer )
   {
      BMSfreeBlockMemoryNull(blkmem, &colexact->storedsolvals);
      assert(colexact->storedsolvals == NULL);
   }

   return SCIP_OKAY;
}

/** save current LP solution values stored in each column */
static
SCIP_RETCODE rowExactStoreSolVals(
   SCIP_ROWEXACT*        rowexact,           /**< exact LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             infeasible          /**< is the solution infeasible? */
   )
{
   SCIP_ROWEXACTSOLVALS* storedsolvals;

   assert(rowexact != NULL);
   assert(blkmem != NULL);

   /* allocate memory for storage */
   if( rowexact->storedsolvals == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &rowexact->storedsolvals) );

      storedsolvals = rowexact->storedsolvals;

      /* store values */
      if( infeasible )
      {
         SCIP_CALL( SCIPrationalCopyBlock(blkmem, &(storedsolvals->dualsol), rowexact->dualfarkas) );
         SCIP_CALL( SCIPrationalCreateBlock(blkmem, &(storedsolvals->activity)) );
         SCIPrationalSetInfinity(storedsolvals->activity);
         storedsolvals->basisstatus = SCIP_BASESTAT_BASIC;  /*lint !e641*/
      }
      else
      {
         SCIP_CALL( SCIPrationalCopyBlock(blkmem, &(storedsolvals->dualsol), rowexact->dualsol) );
         SCIP_CALL( SCIPrationalCopyBlock(blkmem, &(storedsolvals->activity), rowexact->activity) );
         storedsolvals->basisstatus = rowexact->basisstatus; /*lint !e641 !e732*/
      }
   }
   else
   {
      storedsolvals = rowexact->storedsolvals;

      /* store values */
      if( infeasible )
      {
         SCIPrationalSetRational(storedsolvals->dualsol, rowexact->dualfarkas);
         SCIPrationalSetInfinity(storedsolvals->activity);
         storedsolvals->basisstatus = SCIP_BASESTAT_BASIC;  /*lint !e641*/
      }
      else
      {
         SCIPrationalSetRational(storedsolvals->dualsol, rowexact->dualsol);
         SCIPrationalSetRational(storedsolvals->activity, rowexact->activity);
         storedsolvals->basisstatus = rowexact->basisstatus; /*lint !e641 !e732*/
      }
   }

   return SCIP_OKAY;
}

/** restore LP solution values in row */
static
SCIP_RETCODE rowExactRestoreSolVals(
   SCIP_ROWEXACT*        rowexact,           /**< exact LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Longint          validlp,            /**< number of lp for which restored values are valid */
   SCIP_Bool             freebuffer,         /**< should buffer for LP solution values be freed? */
   SCIP_Bool             infeasible          /**< is the solution infeasible? */
   )
{
   SCIP_ROWEXACTSOLVALS* storedsolvals;

   assert(rowexact != NULL);
   assert(blkmem != NULL);

   /* if stored values are available, restore them */
   storedsolvals = rowexact->storedsolvals;
   if( storedsolvals != NULL )
   {
      if( infeasible )
         SCIPrationalSetRational(rowexact->dualfarkas, storedsolvals->dualsol);
      else
         SCIPrationalSetRational(rowexact->dualsol, storedsolvals->dualsol);
      SCIPrationalSetRational(rowexact->activity, storedsolvals->activity);
      rowexact->validactivitylp = validlp;
      rowexact->basisstatus = storedsolvals->basisstatus; /*lint !e641 !e732*/
   }
   /* if the row was created after performing the storage (possibly during probing), we treat it as basic;
    * we make sure to invalidate the reduced cost and farkas coefficient, which are not available
    */
   else
   {
      SCIPrationalSetReal(rowexact->dualsol, 0.0);
      SCIPrationalSetReal(rowexact->dualfarkas, 0.0);
      SCIPrationalSetReal(rowexact->activity, SCIP_INVALID);
      rowexact->validactivitylp = -1;
      rowexact->basisstatus = SCIP_BASESTAT_BASIC; /*lint !e641*/
   }

   /* free memory */
   if( freebuffer )
   {
      BMSfreeBlockMemoryNull(blkmem, &rowexact->storedsolvals);
      assert(rowexact->storedsolvals == NULL);
   }

   return SCIP_OKAY;
}

/** save current LP values dependent on the solution */
static
SCIP_RETCODE lpExactStoreSolVals(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   SCIP_STAT*            stat,               /**< problem statistics */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_LPEXACTSOLVALS* storedsolvals;

   assert(lpexact != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);

   /* allocate memory for storage */
   if( lpexact->storedsolvals == NULL )
   {
      SCIP_ALLOC( BMSallocMemory(&lpexact->storedsolvals) );

      storedsolvals = lpexact->storedsolvals;

      /* store values */
      SCIP_CALL( SCIPrationalCopyBlock(blkmem, &(storedsolvals->lpobjval), lpexact->lpobjval));
      storedsolvals->lpsolstat = lpexact->lpsolstat;
      storedsolvals->primalfeasible = lpexact->primalfeasible;
      storedsolvals->primalchecked = lpexact->primalchecked;
      storedsolvals->dualfeasible = lpexact->dualfeasible;
      storedsolvals->dualchecked = lpexact->dualchecked;
      storedsolvals->solisbasic = lpexact->solisbasic;
      storedsolvals->lpissolved = lpexact->solved;
   }
   else
   {
      storedsolvals = lpexact->storedsolvals;

      /* store values */
      SCIPrationalSetRational(storedsolvals->lpobjval, lpexact->lpobjval);
      storedsolvals->lpsolstat = lpexact->lpsolstat;
      storedsolvals->primalfeasible = lpexact->primalfeasible;
      storedsolvals->primalchecked = lpexact->primalchecked;
      storedsolvals->dualfeasible = lpexact->dualfeasible;
      storedsolvals->dualchecked = lpexact->dualchecked;
      storedsolvals->solisbasic = lpexact->solisbasic;
      storedsolvals->lpissolved = lpexact->solved;
   }

   return SCIP_OKAY;
}

/** restore LP solution values in column */
static
SCIP_RETCODE lpExactRestoreSolVals(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_LPEXACTSOLVALS* storedsolvals;

   assert(lpexact != NULL);
   assert(blkmem != NULL);

   /* if stored values are available, restore them */
   storedsolvals = lpexact->storedsolvals;
   if( storedsolvals != NULL )
   {
#ifdef SCIP_WITH_QSOPTEX
      lpexact->solved = FALSE;
      lpexact->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
#else
      lpexact->solved = storedsolvals->lpissolved;
      lpexact->lpsolstat = storedsolvals->lpsolstat;
#endif
      SCIPrationalSetRational(lpexact->lpobjval, storedsolvals->lpobjval);
      lpexact->primalfeasible = storedsolvals->primalfeasible;
      lpexact->primalchecked = storedsolvals->primalchecked;
      lpexact->dualfeasible = storedsolvals->dualfeasible;
      lpexact->dualchecked = storedsolvals->dualchecked;
      lpexact->solisbasic = storedsolvals->solisbasic;

      /* solution values are stored only for LPs solved without error */
      assert(lpexact->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL ||
         lpexact->lpsolstat == SCIP_LPSOLSTAT_UNBOUNDEDRAY ||
         lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT ||
         lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_ITERLIMIT ||
         lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_TIMELIMIT ||
         lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE ||
         lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_NOTSOLVED);
   }
   /* no values available, mark LP as unsolved */
   else
   {
      lpexact->solved = FALSE;
      //lpexact->validsollp = -1;

      SCIPrationalSetInfinity(lpexact->lpobjval);
      lpexact->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
      lpexact->primalfeasible = FALSE;
      lpexact->primalchecked = FALSE;
      lpexact->dualfeasible = FALSE;
      lpexact->dualchecked = FALSE;
      lpexact->solisbasic = FALSE;
   }

   return SCIP_OKAY;
}

/** locks an unmodifiable row, which forbids further changes; has no effect on modifiable rows */
void SCIProwExactLock(
   SCIP_ROWEXACT*        row                 /**< exact LP row */
   )
{
   assert(row != NULL);

   /* check, if row is modifiable */
   if( !row->modifiable )
   {
      row->nlocks++;
   }
}

/** unlocks a lock of an unmodifiable row; a row with no sealed lock may be modified; has no effect on modifiable rows */
void SCIProwExactUnlock(
   SCIP_ROWEXACT*        row                 /**< exact LP row */
   )
{
   assert(row != NULL);

   /* check, if row is modifiable */
   if( !row->modifiable )
   {
      assert(row->nlocks > 0);
      row->nlocks--;
   }
}

/** ensures that chgrows array can store at least num entries */
static
SCIP_RETCODE ensureChgrowsSizeExact(
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lpexact->nchgrows <= lpexact->chgrowssize);

   if( num > lpexact->chgrowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpexact->chgrows, newsize) );
      lpexact->chgrowssize = newsize;
   }
   assert(num <= lpexact->chgrowssize);

   return SCIP_OKAY;
}


/** notifies exact LP row that its sides were changed */
static
SCIP_RETCODE rowExactSideChanged(
   SCIP_ROWEXACT*        rowexact,           /**< exact LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   SCIP_SIDETYPE         sidetype            /**< type of side: left or right hand side */
   )
{
   assert(rowexact != NULL);
   assert(lpexact != NULL);

   if( rowexact->lpipos >= 0 )
   {
      /* insert row in the chgrows list (if not already there) */
      if( !rowexact->lhschanged && !rowexact->rhschanged )
      {
         SCIP_CALL( ensureChgrowsSizeExact(lpexact, set, lpexact->nchgrows+1) );
         lpexact->chgrows[lpexact->nchgrows] = rowexact;
         lpexact->nchgrows++;
      }

      /* mark side change in the row */
      switch( sidetype )
      {
      case SCIP_SIDETYPE_LEFT:
         rowexact->lhschanged = TRUE;
         break;
      case SCIP_SIDETYPE_RIGHT:
         rowexact->rhschanged = TRUE;
         break;
      default:
         SCIPerrorMessage("unknown exact row side type\n");
         SCIPABORT();
         return SCIP_INVALIDDATA;  /*lint !e527*/
      }

      /* mark the current LP unflushed */
      lpexact->flushed = FALSE;

      assert(lpexact->nchgrows > 0);
   }

   return SCIP_OKAY;
}

/** changes left hand side of exact LP row */
SCIP_RETCODE SCIProwExactChgLhs(
   SCIP_ROWEXACT*        rowexact,           /**< exact LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   SCIP_RATIONAL*        lhs                 /**< new left hand side */
   )
{
   assert(rowexact != NULL);
   assert(lpexact != NULL);

   if( !SCIPrationalIsEQ(rowexact->lhs, lhs) )
   {
      SCIPrationalSetRational(rowexact->lhs, lhs);
      rowexact->lhsreal = SCIPrationalRoundReal(rowexact->lhs, SCIP_R_ROUND_DOWNWARDS);
      SCIP_CALL( rowExactSideChanged(rowexact, set, lpexact, SCIP_SIDETYPE_LEFT) );
   }

   return SCIP_OKAY;
}

/** changes right hand side of exact LP row */
SCIP_RETCODE SCIProwExactChgRhs(
   SCIP_ROWEXACT*        rowexact,           /**< exact LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   SCIP_RATIONAL*        rhs                 /**< new right hand side */
   )
{
   assert(rowexact != NULL);
   assert(lpexact != NULL);

   if( !SCIPrationalIsEQ(rowexact->rhs, rhs) )
   {
      SCIPrationalSetRational(rowexact->rhs, rhs);
      rowexact->rhsreal = SCIPrationalRoundReal(rowexact->rhs, SCIP_R_ROUND_UPWARDS);
      SCIP_CALL( rowExactSideChanged(rowexact, set, lpexact, SCIP_SIDETYPE_RIGHT) );
   }

   return SCIP_OKAY;
}

/** gets solution status of current exact LP */
SCIP_LPSOLSTAT SCIPlpExactGetSolstat(
   SCIP_LPEXACT*         lpexact             /**< current LP data */
   )
{
   assert(lpexact != NULL);

   return (lpexact->flushed ? lpexact->lpsolstat : SCIP_LPSOLSTAT_NOTSOLVED);
}

/** stores exact LP state (like basis information) into LP state object */
SCIP_RETCODE SCIPlpExactGetState(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
   assert(lpexact != NULL);
   assert(lpexact->flushed);
   assert(lpexact->solved);
   assert(blkmem != NULL);
   assert(lpistate != NULL);

   /* check whether there is no lp */
   if( lpexact->nlpicols == 0 && lpexact->nlpirows == 0 )
      *lpistate = NULL;
   else
   {
      SCIP_CALL( SCIPlpiExactGetState(lpexact->lpiexact, blkmem, lpistate) );
   }

   return SCIP_OKAY;
}

/** loads exact LP state (like basis information) into solver */
SCIP_RETCODE SCIPlpExactSetState(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPISTATE*        lpistate,           /**< LP state information (like basis information) */
   SCIP_Bool             wasprimfeas,        /**< primal feasibility when LP state information was stored */
   SCIP_Bool             wasprimchecked,     /**< true if the LP solution has passed the primal feasibility check */
   SCIP_Bool             wasdualfeas,        /**< dual feasibility when LP state information was stored */
   SCIP_Bool             wasdualchecked      /**< true if the LP solution has passed the dual feasibility check */
   )
{
   assert(lpexact != NULL);
   assert(blkmem != NULL);

   /* flush changes to the LP solver */
   SCIP_CALL( SCIPlpExactFlush(lpexact, blkmem, set, eventqueue) );
   assert(lpexact->flushed);

   if( lpexact->solved && lpexact->solisbasic )
      return SCIP_OKAY;

   /* set LPI state in the LP solver */
   if( lpistate == NULL )
      lpexact->solisbasic = FALSE;
   else
   {
      SCIP_CALL( SCIPlpiExactSetState(lpexact->lpiexact, blkmem, lpistate) );
      lpexact->solisbasic = SCIPlpiExactHasStateBasis(lpexact->lpiexact, lpistate);
   }
   /* @todo: setting feasibility to TRUE might be wrong because in probing mode, the state is even saved when the LP was
    *        flushed and solved, also, e.g., when we hit the iteration limit
    */
   lpexact->primalfeasible = wasprimfeas;
   lpexact->primalchecked = wasprimchecked;
   lpexact->dualfeasible = wasdualfeas;
   lpexact->dualchecked = wasdualchecked;

   return SCIP_OKAY;
}

/** frees exact LP state information */
SCIP_RETCODE SCIPlpExactFreeState(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
   assert(lpexact != NULL);

   if( *lpistate != NULL )
   {
      SCIP_CALL( SCIPlpiExactFreeState(lpexact->lpiexact, blkmem, lpistate) );
   }

   return SCIP_OKAY;
}

/** initiates exact LP diving */
SCIP_RETCODE SCIPlpExactStartDive(
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   int c;
   int r;

   assert(lpexact != NULL);
   assert(lpexact->flushed || !lpexact->solved);
   assert(lpexact->fplp->diving);
   assert(!lpexact->diving);
   assert(lpexact->divelpistate == NULL);
   assert(lpexact->divelpwasprimfeas);
   assert(lpexact->divelpwasdualfeas);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lpexact->ndivechgsides == 0);

   SCIPsetDebugMsg(set, "exact diving started (LP flushed: %u, LP solved: %u, solstat: %d)\n",
      lpexact->flushed, lpexact->solved, SCIPlpExactGetSolstat(lpexact));

#ifdef SCIP_MORE_DEBUG
   for( c = 0; c < lpexact->ncols; ++c )
   {
      assert(lpexact->cols[c] != NULL);
      assert(lpexact->cols[c]->var != NULL);
      assert(SCIPvarGetStatusExact(lpexact->cols[c]->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetColExact(lpexact->cols[c]->var) == lpexact->cols[c]);
      assert(SCIPrationalIsEQ(SCIPvarGetObjExact(lpexact->cols[c]->var), lpexact->cols[c]->obj));
      assert(SCIPrationalIsEQ(SCIPvarGetLbLocalExact(lpexact->cols[c]->var), lpexact->cols[c]->lb));
      assert(SCIPrationalIsEQ(SCIPvarGetUbLocalExact(lpexact->cols[c]->var), lpexact->cols[c]->ub));
   }
#endif

   /* save current LPI state (basis information) */
   SCIP_CALL( SCIPlpiExactGetState(lpexact->lpiexact, blkmem, &lpexact->divelpistate) );
   lpexact->divelpwasprimfeas = lpexact->primalfeasible;
   lpexact->divelpwasdualfeas = lpexact->dualfeasible;
   lpexact->divelpwasprimchecked = lpexact->primalchecked;
   lpexact->divelpwasdualchecked = lpexact->dualchecked;

   /* save current LP values dependent on the solution */
   SCIP_CALL( lpExactStoreSolVals(lpexact, stat, blkmem) );
   assert(lpexact->storedsolvals != NULL);
   if( !set->lp_resolverestore && lpexact->solved && lpexact->flushed )
   {
      SCIP_Bool store = TRUE;

      switch( lpexact->lpsolstat )
      {
      case SCIP_LPSOLSTAT_OPTIMAL:
         SCIP_CALL( SCIPlpExactGetSol(lpexact, set, stat, NULL, NULL, FALSE) );
         break;
      case SCIP_LPSOLSTAT_UNBOUNDEDRAY:
         SCIP_CALL( SCIPlpExactGetUnboundedSol(lpexact, set, stat, NULL, NULL) );
         break;
      case SCIP_LPSOLSTAT_OBJLIMIT:
      case SCIP_LPSOLSTAT_ITERLIMIT:
      case SCIP_LPSOLSTAT_TIMELIMIT:
         SCIP_CALL( SCIPlpExactGetSol(lpexact, set, stat, NULL, NULL, FALSE) );
         break;
      case SCIP_LPSOLSTAT_INFEASIBLE:
         SCIP_CALL( SCIPlpExactGetDualfarkas(lpexact, set, stat, NULL, FALSE) );
         break;
      case SCIP_LPSOLSTAT_NOTSOLVED:
      case SCIP_LPSOLSTAT_ERROR:
      default:
         store = FALSE;
      }

      if( store )
      {
         for( c = 0; c < lpexact->ncols; ++c )
         {
            SCIP_CALL( colExactStoreSolVals(lpexact->cols[c], blkmem) );
         }

         for( r = 0; r < lpexact->nrows; ++r )
         {
            SCIP_CALL( rowExactStoreSolVals(lpexact->rows[r], blkmem, lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE) );
         }
      }
   }

   /* store LPI iteration limit */
   SCIP_CALL( SCIPlpiExactGetIntpar(lpexact->lpiexact, SCIP_LPPAR_LPITLIM, &lpexact->divinglpiitlim) );

   /* remember the number of domain changes */
   lpexact->divenolddomchgs = stat->domchgcount;

   /* store current number of rows */
   lpexact->ndivingrows = lpexact->nrows;

   /* switch to diving mode */
   lpexact->diving = TRUE;

   return SCIP_OKAY;
}

/** quits exact LP diving and resets bounds and objective values of columns to the current node's values */
SCIP_RETCODE SCIPlpExactEndDive(
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR**            vars,               /**< array with all active variables */
   int                   nvars               /**< number of active variables */
   )
{
   SCIP_VAR* var;
   int v;

   assert(lpexact != NULL);
   assert(lpexact->diving);
   assert(blkmem != NULL);
   assert(nvars == 0 || vars != NULL);

   SCIPsetDebugMsg(set, "exact diving ended (LP flushed: %u, solstat: %d)\n", lpexact->flushed, SCIPlpExactGetSolstat(lpexact));

   /* reset all columns' objective values and bounds to its original values */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      {
         SCIP_CALL( SCIPcolExactChgObj(SCIPvarGetColExact(var), set, lpexact, SCIPvarGetObjExact(var)) );
         SCIP_CALL( SCIPcolExactChgLb(SCIPvarGetColExact(var), set, lpexact, SCIPvarGetLbLocalExact(var)) );
         SCIP_CALL( SCIPcolExactChgUb(SCIPvarGetColExact(var), set, lpexact, SCIPvarGetUbLocalExact(var)) );
      }
   }

   /* undo changes to left hand sides and right hand sides */
   while( lpexact->ndivechgsides > 0 )
   {
      SCIP_RATIONAL* oldside;
      SCIP_SIDETYPE sidetype;
      SCIP_ROWEXACT* row;

      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &oldside) );

      lpexact->ndivechgsides--;
      SCIPrationalSetRational(oldside, lpexact->divechgsides[lpexact->ndivechgsides]);
      sidetype = lpexact->divechgsidetypes[lpexact->ndivechgsides];
      row = lpexact->divechgrows[lpexact->ndivechgsides];

      if( sidetype == SCIP_SIDETYPE_LEFT )
      {
         SCIP_CALL( SCIProwExactChgLhs(row, set, lpexact, oldside) );
      }
      else
      {
         SCIP_CALL( SCIProwExactChgRhs(row, set, lpexact, oldside) );
      }

      SCIPrationalFreeBuffer(set->buffer, &oldside);
   }

   /* restore LPI iteration limit */
   SCIP_CALL( lpExactSetIterationLimit(lpexact, lpexact->divinglpiitlim) );

   /* reload LPI state saved at start of diving and free it afterwards; it may be NULL, in which case simply nothing
    * happens
    */
   SCIP_CALL( SCIPlpExactSetState(lpexact, blkmem, set, eventqueue, lpexact->divelpistate,
         lpexact->divelpwasprimfeas, lpexact->divelpwasprimchecked, lpexact->divelpwasdualfeas, lpexact->divelpwasdualchecked) );
   SCIP_CALL( SCIPlpExactFreeState(lpexact, blkmem, &lpexact->divelpistate) );
   lpexact->divelpwasprimfeas = TRUE;
   lpexact->divelpwasdualfeas = TRUE;
   lpexact->divelpwasprimchecked = TRUE;
   lpexact->divelpwasdualchecked = TRUE;
   assert(lpexact->divelpistate == NULL);

   /* switch to standard (non-diving) mode */
   lpexact->diving = FALSE;
   lpexact->divingobjchg = FALSE;

   assert(lpexact->storedsolvals != NULL);

   /* we can just always reload the buffered LP solution values at start of diving; this has the advantage that we
    * re-solve as above can lead to a different LP status
    */
   if( lpexact->storedsolvals->lpissolved )
   {
      int c;
      int r;

      /* restore LP solution values in lp data, columns and rows */
      if( lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL ||
          lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_UNBOUNDEDRAY ||
          lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT ||
          lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_ITERLIMIT ||
          lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_TIMELIMIT ||
          lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE )
      {
         SCIP_CALL( lpExactRestoreSolVals(lpexact, blkmem) );

         for( c = 0; c < lpexact->ncols; ++c )
         {
            SCIP_CALL( colExactRestoreSolVals(lpexact->cols[c], blkmem, stat->lpcount, set->lp_freesolvalbuffers) );
         }

         for( r = 0; r < lpexact->nrows; ++r )
         {
            SCIP_CALL( rowExactRestoreSolVals(lpexact->rows[r], blkmem, stat->lpcount, set->lp_freesolvalbuffers,
                  lpexact->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE) );
         }
      }
      else
      {
         SCIP_CALL( lpExactRestoreSolVals(lpexact, blkmem) );
      }
   }
   else
   {
      /* we still need to copy the exact lp objval back because the safe bounding result is saved there */
      SCIPrationalSetRational(lpexact->lpobjval, lpexact->storedsolvals->lpobjval);

      lpexact->solved = FALSE;
      lpexact->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

#ifdef SCIP_MORE_DEBUG
   {
      int c;
      for( c = 0; c < lpexact->ncols; ++c )
      {
         assert(lpexact->cols[c] != NULL);
         assert(lpexact->cols[c]->var != NULL);
         assert(SCIPvarGetStatusExact(lpexact->cols[c]->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetColExact(lpexact->cols[c]->var) == lpexact->cols[c]);
         assert(SCIPrationalIsEQ(SCIPvarGetObjExact(lpexact->cols[c]->var), lpexact->cols[c]->obj));
         assert(SCIPrationalIsEQ(SCIPvarGetLbLocalExact(lpexact->cols[c]->var), lpexact->cols[c]->lb));
         assert(SCIPrationalIsEQ(SCIPvarGetUbLocalExact(lpexact->cols[c]->var), lpexact->cols[c]->ub));
      }
   }
#endif

   return SCIP_OKAY;
}

/** returns whether the exact LP is in exact diving mode */
SCIP_Bool SCIPlpExactDiving(
   SCIP_LPEXACT*         lpexact             /**< current exact LP data */
   )
{
   if( lpexact == NULL )
      return FALSE;

   return lpexact->diving;
}

/** writes exact LP to a file */
SCIP_RETCODE SCIPlpExactWrite(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   const char*           fname               /**< file name */
   )
{
   assert(lpexact != NULL);
   assert(lpexact->flushed);
   assert(fname != NULL);

   SCIP_CALL( SCIPlpiExactWriteLP(lpexact->lpiexact, fname) );

   return SCIP_OKAY;
}

/** overwrites the dual values stored in the fp lp with exact values */
void SCIPlpExactOverwriteFpDualSol(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_Bool             dualfarkas          /**< TRUE if farkas proof, FALSE if dual sol? */
   )
{
   assert(lpexact != NULL);

   for( int c = 0; c < lpexact->ncols; ++c )
   {
      if( dualfarkas )
         lpexact->cols[c]->fpcol->farkascoef = SCIPrationalGetReal(lpexact->cols[c]->farkascoef);
      else
         lpexact->cols[c]->fpcol->redcost = SCIPrationalGetReal(lpexact->cols[c]->redcost);
   }

   for( int r = 0; r < lpexact->nrows; ++r )
   {
      if( dualfarkas )
         lpexact->rows[r]->fprow->dualfarkas = SCIPrationalGetReal(lpexact->rows[r]->dualfarkas);
      else
         lpexact->rows[r]->fprow->dualsol = SCIPrationalGetReal(lpexact->rows[r]->dualsol);
   }
}

/** synchronizes the exact LP with cuts from the floating-point LP */
SCIP_RETCODE SCIPlpExactSyncLPs(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_LP* fplp;
   SCIP_ROWEXACT* rowexact;
   int* rowdset;
   int nrowsfp;
   int nrowsex;
   int i;
   SCIP_Bool removecut;

   if( !set->exact_enable )
      return SCIP_OKAY;

   fplp = lpexact->fplp;

   assert(fplp != NULL);

   nrowsfp = SCIPlpGetNRows(fplp);
   nrowsex = SCIPlpExactGetNRows(lpexact);

   SCIP_CALL( SCIPsetAllocBufferArray(set, &rowdset, lpexact->nrows) );

   removecut = FALSE;

   /* this method should sync the fp-lp withe the exact lp */

   /* remove all rows from exact lp that are not in the floating point lp */
   for( i = 0; i < nrowsex; ++i )
   {
      SCIP_ROW* fprow = lpexact->rows[i]->fprow;

      if( fprow == NULL || !SCIProwIsInLP(fprow) || fprow->lppos != i || removecut )
      {
         removecut = TRUE;
         rowdset[i] = 1;
      }
      else
         rowdset[i] = 0;
   }

   SCIP_CALL( SCIPlpExactDelRowset(lpexact, blkmem, set, rowdset) );

   for( i = 0; i < nrowsfp; ++i )
   {
      rowexact = SCIProwGetRowExact(fplp->rows[i]);
      assert(rowexact != NULL);
      /* if the row is already in lp, do nothing */
      if( !SCIProwExactIsInLP(rowexact) )
      {
         /* add the exact row to the exact lp */
         SCIP_CALL( SCIPlpExactAddRow(lpexact, set, rowexact) );
      }
   }

   SCIPsetFreeBufferArray(set, &rowdset);

   return SCIP_OKAY;
}
