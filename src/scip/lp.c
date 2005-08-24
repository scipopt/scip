/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: lp.c,v 1.205 2005/08/24 17:26:48 bzfpfend Exp $"

/**@file   lp.c
 * @brief  LP management methods and datastructures
 * @author Tobias Achterberg
 * @author Kati Wolter
 *
 *  In LP management, we have to differ between the current LP and the SCIP_LP
 *  stored in the LP solver. All LP methods affect the current LP only. 
 *  Before solving the current LP with the LP solver or setting an LP state,
 *  the LP solvers data has to be updated to the current LP with a call to
 *  lpFlush().
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/intervalarith.h"
#include "scip/clock.h"
#include "scip/lpi.h"
#include "scip/misc.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/sol.h"


/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that chgcols array can store at least num entries */
static
SCIP_RETCODE ensureChgcolsSize(
   SCIP_LP*              lp,                 /**< current LP data */
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

/** ensures, that chgrows array can store at least num entries */
static
SCIP_RETCODE ensureChgrowsSize(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nchgrows <= lp->chgrowssize);
   
   if( num > lp->chgrowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->chgrows, newsize) );
      lp->chgrowssize = newsize;
   }
   assert(num <= lp->chgrowssize);

   return SCIP_OKAY;
}

/** ensures, that lpicols array can store at least num entries */
static
SCIP_RETCODE ensureLpicolsSize(
   SCIP_LP*              lp,                 /**< current LP data */
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
SCIP_RETCODE ensureLpirowsSize(
   SCIP_LP*              lp,                 /**< current LP data */
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

/** ensures, that cols array can store at least num entries */
static
SCIP_RETCODE ensureColsSize(
   SCIP_LP*              lp,                 /**< current LP data */
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

/** ensures, that rows array can store at least num entries */
static
SCIP_RETCODE ensureRowsSize(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nrows <= lp->rowssize);
   
   if( num > lp->rowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->rows, newsize) );
      lp->rowssize = newsize;
   }
   assert(num <= lp->rowssize);

   return SCIP_OKAY;
}

/** ensures, that row array of column can store at least num entries */
static
SCIP_RETCODE colEnsureSize(
   SCIP_COL*             col,                /**< LP column */
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

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->rows, col->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->vals, col->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->linkpos, col->size, newsize) );
      col->size = newsize;
   }
   assert(num <= col->size);

   return SCIP_OKAY;
}

/** ensures, that column array of row can store at least num entries */
SCIP_RETCODE SCIProwEnsureSize(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(row != NULL);
   assert(row->len <= row->size);
   
   if( num > row->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->cols, row->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->cols_index, row->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->vals, row->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->linkpos, row->size, newsize) );
      row->size = newsize;
   }
   assert(num <= row->size);

   return SCIP_OKAY;
}




/*
 * Sorting and searching rows and columns
 */

/** bubble sort part of rows in a column */
static
void colBSort(
   SCIP_COL*             col,                /**< LP column */
   int                   firstpos,           /**< first position to include in the sort */
   int                   lastpos             /**< last position to include in the sort */
   )
{
   SCIP_ROW** rows;
   SCIP_Real* vals;
   int* linkpos;
   SCIP_ROW* tmprow;
   SCIP_Real tmpval;
   int tmplinkpos;
   int tmpindex;
   int pos;
   int sortpos;

   assert(col != NULL);
   assert(0 <= firstpos && firstpos <= lastpos+1 && lastpos < col->len);

   /**@todo do a quick sort here, if many elements are unsorted (sorted-Bool -> sorted-Int?) */
   rows = col->rows;
   vals = col->vals;
   linkpos = col->linkpos;

   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && rows[pos]->index <= rows[pos+1]->index )
            pos++;
         if( pos >= lastpos )
            break;
         assert(rows[pos]->index > rows[pos+1]->index);
         tmprow = rows[pos];
         tmpval = vals[pos];
         tmplinkpos = linkpos[pos];
         tmpindex = tmprow->index;
         do
         {
            rows[pos] = rows[pos+1];
            vals[pos] = vals[pos+1];
            linkpos[pos] = linkpos[pos+1];
            pos++;
         }
         while( pos < lastpos && rows[pos+1]->index < tmpindex );
         rows[pos] = tmprow;
         vals[pos] = tmpval;
         linkpos[pos] = tmplinkpos;
         sortpos = pos;
         pos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      pos = lastpos;
      sortpos = lastpos;
      while( pos > firstpos )
      {
         while( pos > firstpos && rows[pos-1]->index <= rows[pos]->index )
            pos--;
         if( pos <= firstpos )
            break;
         assert(rows[pos-1]->index > rows[pos]->index);
         tmprow = rows[pos];
         tmpval = vals[pos];
         tmplinkpos = linkpos[pos];
         tmpindex = tmprow->index;
         do
         {
            rows[pos] = rows[pos-1];
            vals[pos] = vals[pos-1];
            linkpos[pos] = linkpos[pos-1];
            pos--;
         }
         while( pos > firstpos && rows[pos-1]->index > tmpindex );
         rows[pos] = tmprow;
         vals[pos] = tmpval;
         linkpos[pos] = tmplinkpos;
         sortpos = pos;
         pos--;
      }
      firstpos = sortpos+1;
   }
}

/** bubble sort part of columns in a row */
static
void rowBSort(
   SCIP_ROW*             row,                /**< LP row */
   int                   firstpos,           /**< first position to include in the sort */
   int                   lastpos             /**< last position to include in the sort */
   )
{
   SCIP_COL** cols;
   SCIP_Real* vals;
   int* idx;
   int* linkpos;
   SCIP_COL* tmpcol;
   SCIP_Real tmpval;
   int tmpidx;
   int tmplinkpos;
   int pos;
   int sortpos;

   assert(row != NULL);
   assert(0 <= firstpos && firstpos <= lastpos+1 && lastpos < row->len);

   /**@todo do a quick sort here, if many elements are unsorted (sorted-Bool -> sorted-Int?) */
   cols = row->cols;
   vals = row->vals;
   idx = row->cols_index;
   linkpos = row->linkpos;

#ifndef NDEBUG
   for( pos = 0; pos < row->len; ++pos )
      assert(idx[pos] == cols[pos]->index);
#endif

   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && idx[pos] <= idx[pos+1] )
            pos++;
         if( pos >= lastpos )
            break;
         assert(idx[pos] > idx[pos+1]);
         tmpcol = cols[pos];
         tmpidx = idx[pos];
         tmpval = vals[pos];
         tmplinkpos = linkpos[pos];
         do
         {
            cols[pos] = cols[pos+1];
            idx[pos] = idx[pos+1];
            vals[pos] = vals[pos+1];
            linkpos[pos] = linkpos[pos+1];
            pos++;
         }
         while( pos < lastpos && idx[pos+1] < tmpidx );
         cols[pos] = tmpcol;
         idx[pos] = tmpidx;
         vals[pos] = tmpval;
         linkpos[pos] = tmplinkpos;
         sortpos = pos;
         pos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      pos = lastpos;
      sortpos = lastpos;
      while( pos > firstpos )
      {
         while( pos > firstpos && idx[pos-1] <= idx[pos] )
            pos--;
         if( pos <= firstpos )
            break;
         assert(idx[pos-1] > idx[pos]);
         tmpcol = cols[pos];
         tmpidx = idx[pos];
         tmpval = vals[pos];
         tmplinkpos = linkpos[pos];
         do
         {
            cols[pos] = cols[pos-1];
            idx[pos] = idx[pos-1];
            vals[pos] = vals[pos-1];
            linkpos[pos] = linkpos[pos-1];
            pos--;
         }
         while( pos > firstpos && idx[pos-1] > tmpidx );
         cols[pos] = tmpcol;
         idx[pos] = tmpidx;
         vals[pos] = tmpval;
         linkpos[pos] = tmplinkpos;
         sortpos = pos;
         pos--;
      }
      firstpos = sortpos+1;
   }
}

/** sorts column entries of linked rows currently in the LP such that lower row indices precede higher ones */
static
void colSortLP(
   SCIP_COL*             col                 /**< column to be sorted */
   )
{
   int i;
   
   assert(col != NULL);

   /* check, if column is already sorted in the LP part */
   if( col->lprowssorted )
      return;

   /* sort coefficients */
   colBSort(col, 0, col->nlprows-1);

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
void colSortNonLP(
   SCIP_COL*             col                 /**< column to be sorted */
   )
{
   int i;

   assert(col != NULL);

   /* check, if column is already sorted in the non-LP part */
   if( col->nonlprowssorted )
      return;

   /* sort coefficients */
   colBSort(col, col->nlprows, col->len-1);

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
void rowSortLP(
   SCIP_ROW*             row                 /**< row to be sorted */
   )
{
   int i;

   assert(row != NULL);

   /* check, if row is already sorted in the LP part, or if the sorting should be delayed */
   if( row->lpcolssorted || row->delaysort )
      return;

   /* sort coefficients */
   rowBSort(row, 0, row->nlpcols-1);
   
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
void rowSortNonLP(
   SCIP_ROW*             row                 /**< row to be sorted */
   )
{
   int i;

   assert(row != NULL);

   /* check, if row is already sorted in the non-LP part, or if the sorting should be delayed */
   if( row->nonlpcolssorted || row->delaysort )
      return;

   /* sort coefficients */
   rowBSort(row, row->nlpcols, row->len-1);
   
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

/** searches coefficient in part of the column, returns position in col vector or -1 if not found */
static
int colSearchCoefPart(
   SCIP_COL*             col,                /**< column to be searched in */
   const SCIP_ROW*       row,                /**< coefficient to be searched for */
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
int colSearchCoef(
   SCIP_COL*             col,                /**< column to be searched in */
   const SCIP_ROW*       row                 /**< coefficient to be searched for */
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
      colSortLP(col);
      assert(col->lprowssorted);

      pos = colSearchCoefPart(col, row, 0, col->nlprows-1);
      if( pos >= 0 )
         return pos;
   }

   /* search in the non-LP/unlinked rows */
   if( row->lppos == -1 || col->nunlinked > 0 )
   {
      /* column has to be sorted, such that binary search works */
      colSortNonLP(col);
      assert(col->nonlprowssorted);

      pos = colSearchCoefPart(col, row, col->nlprows, col->len-1);
   }

   return pos;
}

/** searches coefficient in part of the row, returns position in col vector or -1 if not found */
static
int rowSearchCoefPart(
   SCIP_ROW*             row,                /**< row to be searched in */
   const SCIP_COL*       col,                /**< coefficient to be searched for */
   int                   minpos,             /**< first position of search range */
   int                   maxpos              /**< last position of search range */
   )
{
   int pos;
   int idx;
   int searchidx;

   assert(row != NULL);
   assert(col != NULL);

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
int rowSearchCoef(
   SCIP_ROW*             row,                /**< row to be searched in */
   const SCIP_COL*       col                 /**< coefficient to be searched for */
   )
{
   int pos;

   assert(row != NULL);
   assert(col != NULL);

   if( row->delaysort )
      return -1;

   pos = -1;

   /* search in the linked LP columns */
   if( col->lppos >= 0 )
   {
      /* row has to be sorted, such that binary search works */
      rowSortLP(row);
      assert(row->lpcolssorted);

      pos = rowSearchCoefPart(row, col, 0, row->nlpcols-1);
   }

   /* search in the non-LP/unlinked columns */
   if( pos == -1 && (col->lppos == -1 || row->nunlinked > 0) )
   {
      /* row has to be sorted, such that binary search works */
      rowSortNonLP(row);
      assert(row->nonlpcolssorted);

      pos = rowSearchCoefPart(row, col, row->nlpcols, row->len-1);
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

/** moves a coefficient in a column to a different place, and updates all corresponding data structures */
static
void colMoveCoef(
   SCIP_COL*             col,                /**< LP column */
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

   col->rows[newpos] = col->rows[oldpos];
   col->vals[newpos] = col->vals[oldpos];
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
void colSwapCoefs(
   SCIP_COL*             col,                /**< LP column */
   int                   pos1,               /**< position of first coefficient */
   int                   pos2                /**< position of second coefficient */
   )
{
   SCIP_ROW* tmprow;
   SCIP_Real tmpval;
   int tmplinkpos;
   
   assert(col != NULL);
   assert(0 <= pos1 && pos1 < col->len);
   assert(0 <= pos2 && pos2 < col->len);
   assert(col->rows[pos1] != NULL);

   if( pos1 == pos2 )
      return;

   /* swap coefficients */
   tmprow = col->rows[pos2];
   tmpval = col->vals[pos2];
   tmplinkpos = col->linkpos[pos2];

   col->rows[pos2] = col->rows[pos1];
   col->vals[pos2] = col->vals[pos1];
   col->linkpos[pos2] = col->linkpos[pos1];

   col->rows[pos1] = tmprow;
   col->vals[pos1] = tmpval;
   col->linkpos[pos1] = tmplinkpos;

   /* update link position in rowumns */
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
}

/** moves a coefficient in a row to a different place, and updates all corresponding data structures */
static
void rowMoveCoef(
   SCIP_ROW*             row,                /**< LP row */
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
   row->vals[newpos] = row->vals[oldpos];
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
void rowSwapCoefs(
   SCIP_ROW*             row,                /**< LP row */
   int                   pos1,               /**< position of first coefficient */
   int                   pos2                /**< position of second coefficient */
   )
{
   SCIP_COL* tmpcol;
   SCIP_Real tmpval;
   int tmpindex;
   int tmplinkpos;
   
   assert(row != NULL);
   assert(0 <= pos1 && pos1 < row->len);
   assert(0 <= pos2 && pos2 < row->len);
   assert(row->cols[pos1] != NULL);
   assert(row->cols[pos1]->index == row->cols_index[pos1]);

   if( pos1 == pos2 )
      return;

   /* swap coefficients */
   tmpcol = row->cols[pos2];
   tmpindex = row->cols_index[pos2];
   tmpval = row->vals[pos2];
   tmplinkpos = row->linkpos[pos2];

   row->cols[pos2] = row->cols[pos1];
   row->cols_index[pos2] = row->cols_index[pos1];
   row->vals[pos2] = row->vals[pos1];
   row->linkpos[pos2] = row->linkpos[pos1];

   row->cols[pos1] = tmpcol;
   row->cols_index[pos1] = tmpindex;
   row->vals[pos1] = tmpval;
   row->linkpos[pos1] = tmplinkpos;

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
}




#if 0

#ifdef NDEBUG
#define ASSERT(x) do { if( !(x) ) SCIPABORT(); } while( FALSE )
#else
#define ASSERT(x) assert(x)
#endif

static SCIP_Bool msgdisp = FALSE;

static
void checkLinks(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_COL* col;
   SCIP_ROW* row;
   int i;
   int j;

   ASSERT(lp != NULL);

   if( !msgdisp )
   {
      SCIPwarningMessage("LP LINK CHECKING ACTIVATED! THIS IS VERY SLOW!\n");
      msgdisp = TRUE;
   }

   for( i = 0; i < lp->ncols; ++i )
   {
      col = lp->cols[i];
      ASSERT(col != NULL);
      ASSERT(!lp->flushed || col->lppos >= 0 || col->primsol == 0.0);
      ASSERT(!lp->flushed || col->lppos >= 0 || col->farkascoef == 0.0);
      ASSERT(col->nlprows <= col->len);
      ASSERT(col->lppos == -1 || col->lppos >= lp->lpifirstchgcol || col->nunlinked == 0);

      for( j = 0; j < col->len; ++j )
      {
         row = col->rows[j];
         ASSERT(row != NULL);
         ASSERT(!lp->flushed || col->lppos == -1 || col->linkpos[j] >= 0);
         ASSERT(col->linkpos[j] == -1 || row->cols[col->linkpos[j]] == col);
         ASSERT(col->linkpos[j] == -1 || EPSEQ(row->vals[col->linkpos[j]], col->vals[j], 1e-6));
         ASSERT((j < col->nlprows) == (col->linkpos[j] >= 0 && row->lppos >= 0));
      }
   }

   for( i = 0; i < lp->nrows; ++i )
   {
      row = lp->rows[i];
      ASSERT(row != NULL);
      ASSERT(!lp->flushed || row->lppos >= 0 || row->dualsol == 0.0);
      ASSERT(!lp->flushed || row->lppos >= 0 || row->dualfarkas == 0.0);
      ASSERT(row->nlpcols <= row->len);
      ASSERT(row->lppos == -1 || row->lppos >= lp->lpifirstchgrow || row->nunlinked == 0);
      
      for( j = 0; j < row->len; ++j )
      {
         col = row->cols[j];
         ASSERT(col != NULL);
         ASSERT(!lp->flushed || row->lppos == -1 || row->linkpos[j] >= 0);
         ASSERT(row->linkpos[j] == -1 || col->rows[row->linkpos[j]] == row);
         ASSERT(row->linkpos[j] == -1 || EPSEQ(col->vals[row->linkpos[j]], row->vals[j], 1e-6));
         ASSERT((j < row->nlpcols) == (row->linkpos[j] >= 0 && col->lppos >= 0));
      }
   }
}

#undef ASSERT

#else
#define checkLinks(lp) /**/
#endif


/*
 * Changing announcements
 */

/** announces, that the given coefficient in the constraint matrix changed */
static
void coefChanged(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_COL*             col,                /**< LP col */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(row != NULL);
   assert(col != NULL);
   assert(lp != NULL);
   assert(!lp->diving);

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

   row->pseudoactivity = SCIP_INVALID;
   row->minactivity = SCIP_INVALID;
   row->maxactivity = SCIP_INVALID;
   row->validpsactivitydomchg = -1;
   row->validactivitybdsdomchg = -1;
}



/*
 * local column changing methods
 */

/* forward declaration for colAddCoef() */
static
SCIP_RETCODE rowAddCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             val,                /**< value of coefficient */
   int                   linkpos             /**< position of row in the column's row array, or -1 */
   );

/** adds a previously non existing coefficient to an LP column */
static
SCIP_RETCODE colAddCoef(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             val,                /**< value of coefficient */
   int                   linkpos             /**< position of column in the row's col array, or -1 */
   )
{
   int pos;

   assert(blkmem != NULL);
   assert(col != NULL);
   assert(col->nlprows <= col->len);
   assert(col->var != NULL);
   assert(row != NULL);
   assert(!SCIPsetIsZero(set, val));
   /*assert(colSearchCoef(col, row) == -1);*/ /* this assert would lead to slight differences in the solution process */

   SCIP_CALL( colEnsureSize(col, blkmem, set, col->len+1) );
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
         colMoveCoef(col, col->nlprows, pos);
         pos = col->nlprows;
      }
      col->nlprows++;
   }

   /* insert the row at the correct position and update the links */
   col->rows[pos] = row;
   col->vals[pos] = val;
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
         SCIP_CALL( rowAddCoef(row, blkmem, set, lp, col, val, pos) );
         if( row->lppos >= 0 )
            pos = col->nlprows-1;
         linkpos = col->linkpos[pos];

         assert(0 <= linkpos && linkpos < row->len);
         assert(row->cols[linkpos] == col);
         assert(col->rows[pos] == row);
         assert(col->rows[pos]->cols[col->linkpos[pos]] == col);
         assert(col->rows[pos]->linkpos[col->linkpos[pos]] == pos);
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
         rowSwapCoefs(row, linkpos, row->nlpcols-1);
      }
   }

   /* update the sorted flags */
   if( row->lppos >= 0 && linkpos >= 0 )
   {
      assert(col->nlprows >= 1);
      assert(col->rows[col->nlprows-1] == row);
      if( col->nlprows > 1 )
         col->lprowssorted = col->lprowssorted && (col->rows[col->nlprows-2]->index < row->index);
   }
   else
   {
      assert(col->len - col->nlprows >= 1);
      assert(col->rows[col->len-1] == row);
      if( col->len - col->nlprows > 1 )
         col->nonlprowssorted = col->nonlprowssorted && (col->rows[col->len-2]->index < row->index);
   }
   
   coefChanged(row, col, lp);

   SCIPdebugMessage("added coefficient %g * <%s> at position %d (%d/%d) to column <%s> (nunlinked=%d)\n",
      val, row->name, pos, col->nlprows, col->len, SCIPvarGetName(col->var), col->nunlinked);

   return SCIP_OKAY;
}

/** deletes coefficient at given position from column */
static
SCIP_RETCODE colDelCoefPos(
   SCIP_COL*             col,                /**< column to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   int                   pos                 /**< position in column vector to delete */
   )
{
   SCIP_ROW* row;

   assert(col != NULL);
   assert(col->var != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < col->len);
   assert(col->rows[pos] != NULL);
   assert(col->linkpos[pos] == -1 || col->rows[pos]->cols[col->linkpos[pos]] == col);
   assert((pos < col->nlprows) == (col->linkpos[pos] >= 0 && col->rows[pos]->lppos >= 0));

   row = col->rows[pos];
   assert((row->lppos >= 0) == (pos < col->nlprows));

   /*debugMessage("deleting coefficient %g * <%s> at position %d from column <%s>\n", 
     col->vals[pos], row->name, pos, SCIPvarGetName(col->var));*/

   if( col->linkpos[pos] == -1 )
      col->nunlinked--;

   /* if row is a linked LP row, move last linked LP coefficient to position of empty slot (deleted coefficient) */
   if( pos < col->nlprows )
   {
      colMoveCoef(col, col->nlprows-1, pos);
      col->nlprows--;
      pos = col->nlprows;
   }

   /* move last coefficient to position of empty slot */
   colMoveCoef(col, col->len-1, pos);
   col->len--;

   coefChanged(row, col, lp);

   return SCIP_OKAY;
}

/** changes a coefficient at given position of an LP column */
static
SCIP_RETCODE colChgCoefPos(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   int                   pos,                /**< position in column vector to change */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(0 <= pos && pos < col->len);
   assert(col->rows[pos] != NULL);
   assert(col->linkpos[pos] == -1 || col->rows[pos]->cols[col->linkpos[pos]] == col);

   /*debugMessage("changing coefficient %g * <%s> at position %d of column <%s> to %g\n", 
     col->vals[pos], col->rows[pos]->name, pos, SCIPvarGetName(col->var), val);*/

   if( SCIPsetIsZero(set, val) )
   {
      /* delete existing coefficient */
      SCIP_CALL( colDelCoefPos(col, set, lp, pos) );
   }
   else if( !SCIPsetIsEQ(set, col->vals[pos], val) )
   {
      /* change existing coefficient */
      col->vals[pos] = val;
      coefChanged(col->rows[pos], col, lp);
   }

   return SCIP_OKAY;
}




/*
 * local row changing methods
 */

/** update row norms after addition of coefficient */
static
void rowAddNorms(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COL*             col,                /**< column of added coefficient */
   SCIP_Real             val                 /**< value of added coefficient */
   )
{
   SCIP_Real absval;

   assert(row != NULL);
   assert(row->nummaxval >= 0);
   assert(row->numminval >= 0);
   assert(set != NULL);
   assert(col != NULL);

   absval = REALABS(val);
   assert(!SCIPsetIsZero(set, absval));

   /* update min/maxidx */
   row->minidx = MIN(row->minidx, col->index);
   row->maxidx = MAX(row->maxidx, col->index);

   /* update squared euclidean norm and sum norm */
   row->sqrnorm += SQR(absval);
   row->sumnorm += absval;

   /* update objective function scalar product */
   row->objprod += val * col->obj;

   /* update maximal and minimal non-zero value */
   if( row->nummaxval > 0 )
   {
      if( SCIPsetIsGT(set, absval, row->maxval) )
      {
         row->maxval = absval;
         row->nummaxval = 1;
      }
      else if( SCIPsetIsGE(set, absval, row->maxval) )
         row->nummaxval++;
   }
   if( row->numminval > 0 )
   {
      if( SCIPsetIsLT(set, absval, row->minval) )
      {
         row->minval = absval;
         row->numminval = 1;
      }
      else if( SCIPsetIsLE(set, absval, row->minval) )
         row->numminval++;
   }
}

/** update row norms after deletion of coefficient */
static
void rowDelNorms(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COL*             col,                /**< column of deleted coefficient */
   SCIP_Real             val,                /**< value of deleted coefficient */
   SCIP_Bool             updateindex         /**< should the minimal/maximal column index of row be updated? */
   )
{
   SCIP_Real absval;

   assert(row != NULL);
   assert(row->nummaxval >= 0);
   assert(row->numminval >= 0);
   assert(set != NULL);
   assert(col != NULL);

   absval = REALABS(val);
   assert(!SCIPsetIsZero(set, absval));
   assert(row->nummaxval == 0 || SCIPsetIsGE(set, row->maxval, absval));
   assert(row->numminval == 0 || SCIPsetIsLE(set, row->minval, absval));

   /* update min/maxidx validity */
   if( updateindex && (col->index == row->minidx || col->index == row->maxidx) )
      row->validminmaxidx = FALSE;

   /* update squared euclidean norm and sum norm */
   row->sqrnorm -= SQR(absval);
   row->sqrnorm = MAX(row->sqrnorm, 0.0);
   row->sumnorm -= absval;
   row->sumnorm = MAX(row->sumnorm, 0.0);

   /* update objective function scalar product */
   row->objprod -= val * col->obj;

   /* update maximal and minimal non-zero value */
   if( row->nummaxval > 0 )
   {
      if( SCIPsetIsGE(set, absval, row->maxval) )
         row->nummaxval--;
   }
   if( row->numminval > 0 )
   {
      if( SCIPsetIsLE(set, absval, row->minval) )
         row->numminval--;
   }
}

/** adds a previously non existing coefficient to an LP row */
static
SCIP_RETCODE rowAddCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             val,                /**< value of coefficient */
   int                   linkpos             /**< position of row in the column's row array, or -1 */
   )
{
   int pos;

   assert(row != NULL);
   assert(row->nlpcols <= row->len);
   assert(blkmem != NULL);
   assert(col != NULL);
   assert(col->var != NULL);
   assert(col->var_probindex == SCIPvarGetProbindex(col->var));
   assert(!SCIPsetIsZero(set, val));
   /*assert(rowSearchCoef(row, col) == -1);*/ /* this assert would lead to slight differences in the solution process */

   if( row->nlocks > 0 )
   {
      SCIPerrorMessage("cannot add a coefficient to the locked unmodifiable row <%s>\n", row->name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIProwEnsureSize(row, blkmem, set, row->len+1) );
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
         rowMoveCoef(row, row->nlpcols, pos);
         pos = row->nlpcols;
      }
      row->nlpcols++;
   }

   /* insert the column at the correct position and update the links */
   row->cols[pos] = col;
   row->cols_index[pos] = col->index;
   row->vals[pos] = val;
   row->linkpos[pos] = linkpos;
   row->integral = row->integral && SCIPcolIsIntegral(col) && SCIPsetIsIntegral(set, val);
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
         SCIP_CALL( colAddCoef(col, blkmem, set, lp, row, val, pos) );
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
         colSwapCoefs(col, linkpos, col->nlprows-1);
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
   
   rowAddNorms(row, set, col, val);

   coefChanged(row, col, lp);

   SCIPdebugMessage("added coefficient %g * <%s> at position %d (%d/%d) to row <%s> (nunlinked=%d)\n",
      val, SCIPvarGetName(col->var), pos, row->nlpcols, row->len, row->name, row->nunlinked);

   return SCIP_OKAY;
}

/** deletes coefficient at given position from row */
static
SCIP_RETCODE rowDelCoefPos(
   SCIP_ROW*             row,                /**< row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   int                   pos                 /**< position in row vector to delete */
   )
{
   SCIP_COL* col;
   SCIP_Real val;

   assert(row != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < row->len);
   assert(row->cols[pos] != NULL);
   assert(row->linkpos[pos] == -1 || row->cols[pos]->rows[row->linkpos[pos]] == row);
   assert((pos < row->nlpcols) == (row->linkpos[pos] >= 0 && row->cols[pos]->lppos >= 0));

   col = row->cols[pos];
   val = row->vals[pos];
   assert((pos < row->nlpcols) == (col->lppos >= 0 && row->linkpos[pos] >= 0));

   /*debugMessage("deleting coefficient %g * <%s> at position %d from row <%s>\n",
     val, SCIPvarGetName(col->var), pos, row->name);*/

   if( row->nlocks > 0 )
   {
      SCIPerrorMessage("cannot delete a coefficient from the locked unmodifiable row <%s>\n", row->name);
      return SCIP_INVALIDDATA;
   }

   if( row->linkpos[pos] == -1 )
      row->nunlinked--;
   
   /* if column is a linked LP column, move last linked LP coefficient to position of empty slot (deleted coefficient) */
   if( pos < row->nlpcols )
   {
      rowMoveCoef(row, row->nlpcols-1, pos);
      row->nlpcols--;
      pos = row->nlpcols;
   }

   /* move last coefficient to position of empty slot */
   rowMoveCoef(row, row->len-1, pos);
   row->len--;

   rowDelNorms(row, set, col, val, TRUE);

   coefChanged(row, col, lp);

   return SCIP_OKAY;
}

/** changes a coefficient at given position of an LP row */
static
SCIP_RETCODE rowChgCoefPos(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   int                   pos,                /**< position in row vector to change */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   assert(row != NULL);
   assert(0 <= pos && pos < row->len);
   assert(row->cols[pos] != NULL);
   assert(row->linkpos[pos] == -1 || row->cols[pos]->rows[row->linkpos[pos]] == row);

   /*debugMessage("changing coefficient %g * <%s> at position %d of row <%s> to %g\n", 
     row->vals[pos], SCIPvarGetName(row->cols[pos]->var), pos, row->name, val);*/

   if( row->nlocks > 0 )
   {
      SCIPerrorMessage("cannot change a coefficient of the locked unmodifiable row <%s>\n", row->name);
      return SCIP_INVALIDDATA;
   }

   if( SCIPsetIsZero(set, val) )
   {
      /* delete existing coefficient */
      SCIP_CALL( rowDelCoefPos(row, set, lp, pos) );
   }
   else if( !SCIPsetIsEQ(set, row->vals[pos], val) )
   {
      /* change existing coefficient */
      rowDelNorms(row, set, row->cols[pos], row->vals[pos], FALSE);
      row->vals[pos] = val;
      row->integral = row->integral && SCIPcolIsIntegral(row->cols[pos]) && SCIPsetIsIntegral(set, val);
      rowAddNorms(row, set, row->cols[pos], row->vals[pos]);
      coefChanged(row, row->cols[pos], lp);
   }

   return SCIP_OKAY;
}

/** notifies LP row, that its sides were changed */
static
SCIP_RETCODE rowSideChanged(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SIDETYPE         sidetype            /**< type of side: left or right hand side */
   )
{
   assert(row != NULL);
   assert(lp != NULL);

   if( row->lpipos >= 0 )
   {
      /* insert row in the chgrows list (if not already there) */
      if( !row->lhschanged && !row->rhschanged )
      {
         SCIP_CALL( ensureChgrowsSize(lp, set, lp->nchgrows+1) );
         lp->chgrows[lp->nchgrows] = row;
         lp->nchgrows++;
      }
      
      /* mark side change in the row */
      switch( sidetype )
      {
      case SCIP_SIDETYPE_LEFT:
         row->lhschanged = TRUE;
         break;
      case SCIP_SIDETYPE_RIGHT:
         row->rhschanged = TRUE;
         break;
      default:
         SCIPerrorMessage("unknown row side type\n");
         SCIPABORT();
      }

      /* mark the current LP unflushed */
      lp->flushed = FALSE;

      assert(lp->nchgrows > 0);
   }

   return SCIP_OKAY;
}




/*
 * double linked coefficient matrix methods 
 */

/** insert column coefficients in corresponding rows */
static
SCIP_RETCODE colLink(
   SCIP_COL*             col,                /**< column data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(col != NULL);
   assert(col->var != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( col->nunlinked > 0 )
   {
      SCIPdebugMessage("linking column <%s>\n", SCIPvarGetName(col->var));

      /* unlinked rows can only be in the non-LP/unlinked rows part of the rows array */
      for( i = col->nlprows; i < col->len; ++i )
      {
         assert(!SCIPsetIsZero(set, col->vals[i]));
         if( col->linkpos[i] == -1 )
         {
            /* this call might swap the current row with the first non-LP/not linked row, but this is of no harm */
            SCIP_CALL( rowAddCoef(col->rows[i], blkmem, set, lp, col, col->vals[i], i) );
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

/** removes column coefficients from corresponding rows */
static
SCIP_RETCODE colUnlink(
   SCIP_COL*             col,                /**< column data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(col != NULL);
   assert(col->var != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( col->nunlinked < col->len )
   {
      SCIPdebugMessage("unlinking column <%s>\n", SCIPvarGetName(col->var));
      for( i = 0; i < col->len; ++i )
      {
         if( col->linkpos[i] >= 0 )
         {
            assert(col->rows[i]->cols[col->linkpos[i]] == col);
            SCIP_CALL( rowDelCoefPos(col->rows[i], set, lp, col->linkpos[i]) );
            col->linkpos[i] = -1;
            col->nunlinked++;
         }
      }
   }
   assert(col->nunlinked == col->len);

   checkLinks(lp);

   return SCIP_OKAY;
}

/** insert row coefficients in corresponding columns */
static
SCIP_RETCODE rowLink(
   SCIP_ROW*             row,                /**< row data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(row != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( row->nunlinked > 0 )
   {
      SCIPdebugMessage("linking row <%s>\n", row->name);

      /* unlinked columns can only be in the non-LP/unlinked columns part of the cols array */
      for( i = row->nlpcols; i < row->len; ++i )
      {
         assert(!SCIPsetIsZero(set, row->vals[i]));
         if( row->linkpos[i] == -1 )
         {
            /* this call might swap the current column with the first non-LP/not linked column, but this is of no harm */
            SCIP_CALL( colAddCoef(row->cols[i], blkmem, set, lp, row, row->vals[i], i) );
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
SCIP_RETCODE rowUnlink(
   SCIP_ROW*             row,                /**< row data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(row != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( row->nunlinked < row->len )
   {
      SCIPdebugMessage("unlinking row <%s>\n", row->name);
      for( i = 0; i < row->len; ++i )
      {
         if( row->linkpos[i] >= 0 )
         {
            assert(row->cols[i]->rows[row->linkpos[i]] == row);
            SCIP_CALL( colDelCoefPos(row->cols[i], set, lp, row->linkpos[i]) );
            row->nunlinked++;
         }
      }
   }
   assert(row->nunlinked == row->len);

   checkLinks(lp);

   return SCIP_OKAY;
}




/*
 * local LP parameter methods
 */

/** sets parameter of type int in LP solver, ignoring unknown parameters */
static
SCIP_RETCODE lpSetIntpar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   int                   value,              /**< value to set parameter to */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   SCIP_RETCODE retcode;

   assert(lp != NULL);
   assert(success != NULL);

   retcode = SCIPlpiSetIntpar(lp->lpi, lpparam, value);

   /* check, if parameter is unknown */
   if( retcode == SCIP_PARAMETERUNKNOWN )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   *success = TRUE;

   return retcode;
}

/** sets parameter of type SCIP_Bool in LP solver, ignoring unknown parameters */
static
SCIP_RETCODE lpSetBoolpar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   SCIP_Bool             value,              /**< value to set parameter to */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   return lpSetIntpar(lp, lpparam, (int)value, success);
}

/** sets parameter of type SCIP_Real in LP solver, ignoring unknown parameters */
static
SCIP_RETCODE lpSetRealpar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   SCIP_Real             value,              /**< value to set parameter to */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   SCIP_RETCODE retcode;

   assert(lp != NULL);
   assert(success != NULL);

   retcode = SCIPlpiSetRealpar(lp->lpi, lpparam, value);

   /* check, if parameter is unknown */
   if( retcode == SCIP_PARAMETERUNKNOWN )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   *success = TRUE;

   return retcode;
}

#ifndef NDEBUG
/** checks, that parameter of type int in LP solver has the given value, ignoring unknown parameters */
static
SCIP_RETCODE lpCheckIntpar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   int                   value               /**< value parameter should have */
   )
{
   SCIP_RETCODE retcode;
   int lpivalue;

   assert(lp != NULL);

   retcode = SCIPlpiGetIntpar(lp->lpi, lpparam, &lpivalue);

   /* ignore unknown parameter error */
   if( retcode == SCIP_PARAMETERUNKNOWN )
      return SCIP_OKAY;

   /* check value */
   assert(lpivalue == value);

   return retcode;
}

/** checks, that parameter of type SCIP_Bool in LP solver has the given value, ignoring unknown parameters */
static
SCIP_RETCODE lpCheckBoolpar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   SCIP_Bool             value               /**< value parameter should have */
   )
{
   return lpCheckIntpar(lp, lpparam, (int)value);
}

/** checks, that parameter of type SCIP_Real in LP solver has the given value, ignoring unknown parameters */
static
SCIP_RETCODE lpCheckRealpar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   SCIP_Real             value               /**< value parameter should have */
   )
{
   SCIP_RETCODE retcode;
   SCIP_Real lpivalue;

   assert(lp != NULL);

   retcode = SCIPlpiGetRealpar(lp->lpi, lpparam, &lpivalue);

   /* ignore unknown parameter error */
   if( retcode == SCIP_PARAMETERUNKNOWN )
      return SCIP_OKAY;

   /* check value */
   if( lpivalue != value ) /*lint !e777*/
      return SCIP_LPERROR;

   return retcode;
}
#else
#define lpCheckIntpar(lp, lpparam, value) SCIP_OKAY
#define lpCheckBoolpar(lp, lpparam, value) SCIP_OKAY
#define lpCheckRealpar(lp, lpparam, value) SCIP_OKAY
#endif

/** sets the upper objective limit of the LP solver */
static
SCIP_RETCODE lpSetUobjlim(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             uobjlim             /**< new feasibility tolerance */
   )
{
   assert(lp != NULL);
   assert(set != NULL);

   /* if we want so solve exactly, we cannot rely on the LP solver's objective limit handling */
   if( set->misc_exactsolve )
      return SCIP_OKAY;

   SCIP_CALL( lpCheckRealpar(lp, SCIP_LPPAR_UOBJLIM, lp->lpiuobjlim) );

   if( uobjlim != lp->lpiuobjlim ) /*lint !e777*/
   {
      SCIP_Bool success;

      SCIP_CALL( lpSetRealpar(lp, SCIP_LPPAR_UOBJLIM, uobjlim, &success) );
      if( success )
      {
         /* mark the current solution invalid */
         lp->solved = FALSE;
         lp->primalfeasible = FALSE;
         lp->lpobjval = SCIP_INVALID;
         lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         lp->lpiuobjlim = uobjlim;
      }
   }

   return SCIP_OKAY;
}

/** sets the feasibility tolerance of the LP solver */
static
SCIP_RETCODE lpSetFeastol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             feastol,            /**< new feasibility tolerance */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(feastol >= 0.0);
   assert(success != NULL);

   SCIP_CALL( lpCheckRealpar(lp, SCIP_LPPAR_FEASTOL, lp->lpifeastol) );

   if( feastol != lp->lpifeastol ) /*lint !e777*/
   {
      SCIP_CALL( lpSetRealpar(lp, SCIP_LPPAR_FEASTOL, feastol, success) );
      if( *success )
      {
         if( lp->nrows > 0 && feastol < lp->lpifeastol )
         {
            /* mark the current solution invalid */
            lp->solved = FALSE;
            lp->primalfeasible = FALSE;
            lp->lpobjval = SCIP_INVALID;
            lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         }
         lp->lpifeastol = feastol;
      }
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the reduced costs feasibility tolerance of the LP solver */
static
SCIP_RETCODE lpSetDualfeastol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             dualfeastol,        /**< new reduced costs feasibility tolerance */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(dualfeastol >= 0.0);
   assert(success != NULL);

   SCIP_CALL( lpCheckRealpar(lp, SCIP_LPPAR_DUALFEASTOL, lp->lpidualfeastol) );

   if( dualfeastol != lp->lpidualfeastol ) /*lint !e777*/
   {
      SCIP_CALL( lpSetRealpar(lp, SCIP_LPPAR_DUALFEASTOL, dualfeastol, success) );
      if( *success )
      {
         if( lp->nrows > 0 && dualfeastol < lp->lpidualfeastol )
         {
            /* mark the current solution invalid */
            lp->solved = FALSE;
            lp->dualfeasible = FALSE;
            lp->lpobjval = SCIP_INVALID;
            lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         }
         lp->lpidualfeastol = dualfeastol;
      }
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the convergence tolerance used in barrier algorithm of the LP solver */
static
SCIP_RETCODE lpSetBarrierconvtol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             barrierconvtol,     /**< new convergence tolerance used in barrier algorithm */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(barrierconvtol >= 0.0);
   assert(success != NULL);

   SCIP_CALL( lpCheckRealpar(lp, SCIP_LPPAR_BARRIERCONVTOL, lp->lpibarrierconvtol) );

   if( barrierconvtol != lp->lpibarrierconvtol ) /*lint !e777*/
   {
      SCIP_CALL( lpSetRealpar(lp, SCIP_LPPAR_BARRIERCONVTOL, barrierconvtol, success) );
      if( *success )
      {
         if( lp->nrows > 0 && barrierconvtol < lp->lpibarrierconvtol
            && (lp->lastlpalgo == SCIP_LPALGO_BARRIER || lp->lastlpalgo == SCIP_LPALGO_BARRIERCROSSOVER) )
         {
            /* mark the current solution invalid */
            lp->solved = FALSE;
            lp->dualfeasible = FALSE;
            lp->lpobjval = SCIP_INVALID;
            lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         }
         lp->lpibarrierconvtol = barrierconvtol;
      }
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the FROMSCRATCH setting of the LP solver */
static
SCIP_RETCODE lpSetFromscratch(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             fromscratch,        /**< new FROMSCRATCH setting */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   SCIP_CALL( lpCheckBoolpar(lp, SCIP_LPPAR_FROMSCRATCH, lp->lpifromscratch) );

   if( fromscratch != lp->lpifromscratch )
   {
      SCIP_CALL( lpSetBoolpar(lp, SCIP_LPPAR_FROMSCRATCH, fromscratch, success) );
      if( *success )
         lp->lpifromscratch = fromscratch;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the FASTMIP setting of the LP solver */
static
SCIP_RETCODE lpSetFastmip(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             fastmip,            /**< new FASTMIP setting */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   SCIP_CALL( lpCheckBoolpar(lp, SCIP_LPPAR_FASTMIP, lp->lpifastmip) );

   if( fastmip != lp->lpifastmip )
   {
      SCIP_CALL( lpSetBoolpar(lp, SCIP_LPPAR_FASTMIP, fastmip, success) );
      if( *success )
         lp->lpifastmip = fastmip;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the SCALING setting of the LP solver */
static
SCIP_RETCODE lpSetScaling(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             scaling,            /**< new SCALING setting */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   SCIP_CALL( lpCheckBoolpar(lp, SCIP_LPPAR_SCALING, lp->lpiscaling) );

   if( scaling != lp->lpiscaling )
   {
      SCIP_CALL( lpSetBoolpar(lp, SCIP_LPPAR_SCALING, scaling, success) );
      if( *success )
         lp->lpiscaling = scaling;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the PRESOLVING setting of the LP solver */
static
SCIP_RETCODE lpSetPresolving(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             presolving,         /**< new PRESOLVING setting */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   SCIP_CALL( lpCheckBoolpar(lp, SCIP_LPPAR_PRESOLVING, lp->lpipresolving) );

   if( presolving != lp->lpipresolving )
   {
      SCIP_CALL( lpSetBoolpar(lp, SCIP_LPPAR_PRESOLVING, presolving, success) );
      if( *success )
         lp->lpipresolving = presolving;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the iteration limit of the LP solver */
static
SCIP_RETCODE lpSetIterationLimit(
   SCIP_LP*              lp,                 /**< current LP data */
   int                   itlim               /**< maximal number of LP iterations to perform, or -1 for no limit */
   )
{
   SCIP_Bool success;

   assert(lp != NULL);
   assert(itlim >= -1);

   if( itlim == -1 )
      itlim = INT_MAX;

   SCIP_CALL( lpCheckIntpar(lp, SCIP_LPPAR_LPITLIM, lp->lpiitlim) );

   if( itlim != lp->lpiitlim )
   {
      SCIP_CALL( lpSetIntpar(lp, SCIP_LPPAR_LPITLIM, itlim, &success) );
      if( success )
      {
         if( itlim > lp->lpiitlim )
         {
            /* mark the current solution invalid */
            lp->solved = FALSE;
            lp->lpobjval = SCIP_INVALID;
            lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         }
         lp->lpiitlim = itlim;
      }
   }

   return SCIP_OKAY;
}

/** sets the pricing strategy of the LP solver */
static
SCIP_RETCODE lpSetPricing(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_PRICING          pricing             /**< pricing strategy */
   )
{
   SCIP_Bool success;

   assert(lp != NULL);

   SCIP_CALL( lpCheckIntpar(lp, SCIP_LPPAR_PRICING, (int)lp->lpipricing) );

   if( pricing != lp->lpipricing )
   {
      SCIP_CALL( lpSetIntpar(lp, SCIP_LPPAR_PRICING, (int)pricing, &success) );
      if( success )
         lp->lpipricing = pricing;
   }

   return SCIP_OKAY;
}

/** sets the pricing strategy of the LP solver (given the character representation of the strategy) */
static
SCIP_RETCODE lpSetPricingChar(
   SCIP_LP*              lp,                 /**< current LP data */
   char                  pricingchar         /**< character representing the pricing strategy */
   )
{
   SCIP_PRICING pricing;

   switch( pricingchar )
   {
   case 'a':
      pricing = SCIP_PRICING_AUTO;
      break;
   case 'f':
      pricing = SCIP_PRICING_FULL;
      break;
   case 'p':
      pricing = SCIP_PRICING_PARTIAL;
      break;
   case 's':
      pricing = SCIP_PRICING_STEEP;
      break;
   case 'q':
      pricing = SCIP_PRICING_STEEPQSTART;
      break;
   case 'd':
      pricing = SCIP_PRICING_DEVEX;
      break;
   default:
      SCIPerrorMessage("invalid LP pricing parameter <%c>\n", pricingchar);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( lpSetPricing(lp, pricing) );

   return SCIP_OKAY;
}

/** sets the verbosity of the LP solver */
static
SCIP_RETCODE lpSetLPInfo(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             lpinfo              /**< should the LP solver display status messages? */
   )
{
   SCIP_Bool success;

   assert(lp != NULL);

   SCIP_CALL( lpCheckBoolpar(lp, SCIP_LPPAR_LPINFO, lp->lpilpinfo) );

   if( lpinfo != lp->lpilpinfo )
   {
      SCIP_CALL( lpSetBoolpar(lp, SCIP_LPPAR_LPINFO, lpinfo, &success) );
      if( success )
         lp->lpilpinfo = lpinfo;
   }

   return SCIP_OKAY;
}




/*
 * Column methods
 */

/** creates an LP column */
SCIP_RETCODE SCIPcolCreate(
   SCIP_COL**            col,                /**< pointer to column data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             var,                /**< variable, this column represents */
   int                   len,                /**< number of nonzeros in the column */
   SCIP_ROW**            rows,               /**< array with rows of column entries */
   SCIP_Real*            vals,               /**< array with coefficients of column entries */
   SCIP_Bool             removeable          /**< should the column be removed from the LP due to aging or cleanup? */
   )
{
   int i;

   assert(col != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(var != NULL);
   assert(len >= 0);
   assert(len == 0 || (rows != NULL && vals != NULL));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, col) );

   if( len > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*col)->rows, rows, len) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*col)->vals, vals, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*col)->linkpos, len) );

      for( i = 0; i < len; ++i )
      {
         assert(rows[i] != NULL);
         assert(!SCIPsetIsZero(set, vals[i]));
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
   (*col)->obj = SCIPvarGetObj(var);
   (*col)->lb = SCIPvarGetLbLocal(var);
   (*col)->ub = SCIPvarGetUbLocal(var);
   (*col)->flushedobj = 0.0;
   (*col)->flushedlb = 0.0;
   (*col)->flushedub = 0.0;
   (*col)->index = stat->ncolidx++;
   (*col)->size = len;
   (*col)->len = len;
   (*col)->nlprows = 0;
   (*col)->nunlinked = len;
   (*col)->lppos = -1;
   (*col)->lpipos = -1;
   (*col)->lpdepth = -1;
   (*col)->primsol = 0.0;
   (*col)->redcost = SCIP_INVALID;
   (*col)->farkascoef = SCIP_INVALID;
   (*col)->minprimsol = (*col)->ub;
   (*col)->maxprimsol = (*col)->lb;
   (*col)->sbdown = SCIP_INVALID;
   (*col)->sbup = SCIP_INVALID;
   (*col)->sbsolval  = SCIP_INVALID;
   (*col)->sblpobjval = SCIP_INVALID;
   (*col)->sbnode = -1;
   (*col)->validredcostlp = -1;
   (*col)->validfarkaslp = -1;
   (*col)->validsblp = -1;
   (*col)->sbitlim = -1;
   (*col)->nsbcalls = 0;
   (*col)->age = 0;
   (*col)->obsoletenode = -1;
   (*col)->var_probindex = SCIPvarGetProbindex(var);
   (*col)->basisstatus = SCIP_BASESTAT_ZERO; /*lint !e641*/
   (*col)->lprowssorted = TRUE;
   (*col)->nonlprowssorted = (len <= 1);
   (*col)->objchanged = FALSE;
   (*col)->lbchanged = FALSE;
   (*col)->ubchanged = FALSE;
   (*col)->coefchanged = FALSE;
   (*col)->integral = SCIPvarIsIntegral(var);
   (*col)->removeable = removeable;
   (*col)->sbdownvalid = FALSE;
   (*col)->sbupvalid = FALSE;

   return SCIP_OKAY;
}

/** frees an LP column */
SCIP_RETCODE SCIPcolFree(
   SCIP_COL**            col,                /**< pointer to LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(blkmem != NULL);
   assert(col != NULL);
   assert(*col != NULL);
   assert((*col)->var != NULL);
   assert(SCIPvarGetStatus((*col)->var) == SCIP_VARSTATUS_COLUMN);
   assert(&(*col)->var->data.col == col); /* SCIPcolFree() has to be called from SCIPvarFree() */
   assert((*col)->lppos == -1);
   assert((*col)->lpipos == -1);

   /* remove column indices from corresponding rows */
   SCIP_CALL( colUnlink(*col, blkmem, set, lp) );

   BMSfreeBlockMemoryArrayNull(blkmem, &(*col)->rows, (*col)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*col)->vals, (*col)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*col)->linkpos, (*col)->size);
   BMSfreeBlockMemory(blkmem, col);

   return SCIP_OKAY;
}

/** sorts column entries such that LP rows precede non-LP rows and inside both parts lower row indices precede higher ones
 */
void SCIPcolSort(
   SCIP_COL*             col                 /**< column to be sorted */
   )
{
   /* sort LP rows */
   colSortLP(col);

   /* sort non-LP rows */
   colSortNonLP(col);
}

/** adds a previously non existing coefficient to an LP column */
SCIP_RETCODE SCIPcolAddCoef(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);

   SCIP_CALL( colAddCoef(col, blkmem, set, lp, row, val, -1) );

   checkLinks(lp);

   return SCIP_OKAY;
}

/** deletes existing coefficient from column */
SCIP_RETCODE SCIPcolDelCoef(
   SCIP_COL*             col,                /**< column to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(col != NULL);
   assert(col->var != NULL);
   assert(lp != NULL);
   assert(!lp->diving);
   assert(row != NULL);

   /* search the position of the row in the column's row vector */
   pos = colSearchCoef(col, row);
   if( pos == -1 )
   {
      SCIPerrorMessage("coefficient for row <%s> doesn't exist in column <%s>\n", row->name, SCIPvarGetName(col->var));
      return SCIP_INVALIDDATA;
   }
   assert(0 <= pos && pos < col->len);
   assert(col->rows[pos] == row);

   /* if row knows of the column, remove the column from the row's col vector */
   if( col->linkpos[pos] >= 0 )
   {
      assert(row->cols[col->linkpos[pos]] == col);
      assert(row->cols_index[col->linkpos[pos]] == col->index);
      assert(SCIPsetIsEQ(set, row->vals[col->linkpos[pos]], col->vals[pos]));
      SCIP_CALL( rowDelCoefPos(row, set, lp, col->linkpos[pos]) );
   }

   /* delete the row from the column's row vector */
   SCIP_CALL( colDelCoefPos(col, set, lp, pos) );
   
   checkLinks(lp);

   return SCIP_OKAY;
}

/** changes or adds a coefficient to an LP column */
SCIP_RETCODE SCIPcolChgCoef(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   int pos;

   assert(col != NULL);
   assert(lp != NULL);
   assert(!lp->diving);
   assert(row != NULL);

   /* search the position of the row in the column's row vector */
   pos = colSearchCoef(col, row);

   /* check, if row already exists in the column's row vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( colAddCoef(col, blkmem, set, lp, row, val, -1) );
   }
   else
   {
      /* modifify already existing coefficient */
      assert(0 <= pos && pos < col->len);
      assert(col->rows[pos] == row);

      /* if row knows of the column, change the corresponding coefficient in the row */
      if( col->linkpos[pos] >= 0 )
      {
         assert(row->cols[col->linkpos[pos]] == col);
         assert(row->cols_index[col->linkpos[pos]] == col->index);
         assert(SCIPsetIsEQ(set, row->vals[col->linkpos[pos]], col->vals[pos]));
         SCIP_CALL( rowChgCoefPos(row, set, lp, col->linkpos[pos], val) );
      }

      /* change the coefficient in the column */
      SCIP_CALL( colChgCoefPos(col, set, lp, pos, val) );
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** increases value of an existing or nonexisting coefficient in an LP column */
SCIP_RETCODE SCIPcolIncCoef(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             incval              /**< value to add to the coefficient */
   )
{
   int pos;

   assert(col != NULL);
   assert(lp != NULL);
   assert(!lp->diving);
   assert(row != NULL);

   if( SCIPsetIsZero(set, incval) )
      return SCIP_OKAY;

   /* search the position of the row in the column's row vector */
   pos = colSearchCoef(col, row);

   /* check, if row already exists in the column's row vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( colAddCoef(col, blkmem, set, lp, row, incval, -1) );
   }
   else
   {
      /* modifify already existing coefficient */
      assert(0 <= pos && pos < col->len);
      assert(col->rows[pos] == row);

      /* if row knows of the column, change the corresponding coefficient in the row */
      if( col->linkpos[pos] >= 0 )
      {
         assert(row->cols[col->linkpos[pos]] == col);
         assert(row->cols_index[col->linkpos[pos]] == col->index);
         assert(SCIPsetIsEQ(set, row->vals[col->linkpos[pos]], col->vals[pos]));
         SCIP_CALL( rowChgCoefPos(row, set, lp, col->linkpos[pos], col->vals[pos] + incval) );
      }

      /* change the coefficient in the column */
      SCIP_CALL( colChgCoefPos(col, set, lp, pos, col->vals[pos] + incval) );
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** changes objective value of column */
SCIP_RETCODE SCIPcolChgObj(
   SCIP_COL*             col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newobj              /**< new objective value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(lp != NULL);
   
   SCIPdebugMessage("changing objective value of column <%s> from %f to %f\n", SCIPvarGetName(col->var), col->obj, newobj);

   if( col->lpipos >= 0 && !SCIPsetIsEQ(set, col->obj, newobj) )
   {
      /* insert column in the chgcols list (if not already there) */
      if( !col->objchanged && !col->lbchanged && !col->ubchanged )
      {
         SCIP_CALL( ensureChgcolsSize(lp, set, lp->nchgcols+1) );
         lp->chgcols[lp->nchgcols] = col;
         lp->nchgcols++;
      }
      
      /* mark objective value change in the column */
      col->objchanged = TRUE;
      
      /* mark the current LP unflushed */
      lp->flushed = FALSE;

      assert(lp->nchgcols > 0);
   }  

   /* store new objective function value */
   col->obj = newobj;

   return SCIP_OKAY;
}

/** changes lower bound of column */
SCIP_RETCODE SCIPcolChgLb(
   SCIP_COL*             col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newlb               /**< new lower bound value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(lp != NULL);
   
   SCIPdebugMessage("changing lower bound of column <%s> from %f to %f\n", SCIPvarGetName(col->var), col->lb, newlb);

   if( col->lpipos >= 0 && !SCIPsetIsEQ(set, col->lb, newlb) )
   {
      /* insert column in the chgcols list (if not already there) */
      if( !col->objchanged && !col->lbchanged && !col->ubchanged )
      {
         SCIP_CALL( ensureChgcolsSize(lp, set, lp->nchgcols+1) );
         lp->chgcols[lp->nchgcols] = col;
         lp->nchgcols++;
      }
      
      /* mark bound change in the column */
      col->lbchanged = TRUE;
      
      /* mark the current LP unflushed */
      lp->flushed = FALSE;

      assert(lp->nchgcols > 0);
   }  

   col->lb = newlb;

   return SCIP_OKAY;
}

/** changes upper bound of column */
SCIP_RETCODE SCIPcolChgUb(
   SCIP_COL*             col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newub               /**< new upper bound value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(lp != NULL);
   
   SCIPdebugMessage("changing upper bound of column <%s> from %f to %f\n", SCIPvarGetName(col->var), col->ub, newub);

   if( col->lpipos >= 0 && !SCIPsetIsEQ(set, col->ub, newub) )
   {
      /* insert column in the chgcols list (if not already there) */
      if( !col->objchanged && !col->lbchanged && !col->ubchanged )
      {
         SCIP_CALL( ensureChgcolsSize(lp, set, lp->nchgcols+1) );
         lp->chgcols[lp->nchgcols] = col;
         lp->nchgcols++;
      }
      
      /* mark bound change in the column */
      col->ubchanged = TRUE;
      
      /* mark the current LP unflushed */
      lp->flushed = FALSE;

      assert(lp->nchgcols > 0);
   }  

   col->ub = newub;

   return SCIP_OKAY;
}

/** calculates the reduced costs of a column using the given dual solution vector */
SCIP_Real SCIPcolCalcRedcost(
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real*            dualsol             /**< dual solution vector for current LP rows */
   )
{
   SCIP_ROW* row;
   SCIP_Real redcost;
   int i;
   
   assert(col != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(dualsol != NULL);

   redcost = col->obj;
   for( i = 0; i < col->nlprows; ++i )
   {
      row = col->rows[i];
      assert(row != NULL);
      assert(row->lppos >= 0);
      redcost -= col->vals[i] * dualsol[row->lppos];
   }

   if( col->nunlinked > 0 )
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->lppos == -1 || col->linkpos[i] == -1);
         if( row->lppos >= 0 )
            redcost -= col->vals[i] * dualsol[row->lppos];
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
      }
   }
#endif

   return redcost;
}

/** calculates the reduced costs of a column using the dual solution stored in the rows */
static
SCIP_Real colCalcInternalRedcost(
   SCIP_COL*             col                 /**< LP column */
   )
{
   SCIP_ROW* row;
   SCIP_Real redcost;
   int i;

   assert(col != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);

   redcost = col->obj;
   for( i = 0; i < col->nlprows; ++i )
   {
      row = col->rows[i];
      assert(row != NULL);
      assert(row->dualsol < SCIP_INVALID);
      assert(row->lppos >= 0);
      assert(col->linkpos[i] >= 0);
      redcost -= col->vals[i] * row->dualsol;
   }

   if( col->nunlinked > 0 )
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->lppos >= 0 || row->dualsol == 0.0);
         assert(row->lppos == -1 || col->linkpos[i] == -1);
         if( row->lppos >= 0 )
            redcost -= col->vals[i] * row->dualsol;
      }
   }
#ifndef NDEBUG
   else
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->dualsol == 0.0);
         assert(row->lppos == -1);
         assert(col->linkpos[i] >= 0);
      }
   }
#endif

   return redcost;
}

/** gets the reduced costs of a column in last LP or after recalculation */
SCIP_Real SCIPcolGetRedcost(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(col != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(col->validredcostlp <= stat->lpcount);
   assert(lp->validsollp == stat->lpcount);

   if( col->validredcostlp < stat->lpcount )
   {
      col->redcost = colCalcInternalRedcost(col);
      col->validredcostlp = stat->lpcount;
   }
   assert(col->validredcostlp == stat->lpcount);
   assert(col->redcost < SCIP_INVALID);

   return col->redcost;
}

/** gets the feasibility of (the dual row of) a column in last LP or after recalculation */
SCIP_Real SCIPcolGetFeasibility(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(col != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->validsollp == stat->lpcount);

   /* A column's reduced cost is defined as
    *   redcost  = obj - activity,  activity = y^T * col.   (activity = obj - redcost)
    * The activity is equal to the activity of the corresponding row in the dual LP.
    * The column's feasibility is the feasibility of the corresponding row in the dual LP.
    * The sides of the dual row depend on the bounds of the column:
    *  - lb == ub      :  dual row is a free row with infinite sides
    *  -  0 <= lb <  ub:         activity <= obj  =>  0 <= redcost
    *  - lb <   0 <  ub:  obj <= activity <= obj  =>  0 <= redcost <= 0
    *  - lb <  ub <=  0:  obj <= activity         =>       redcost <= 0
    */
   if( SCIPsetIsEQ(set, col->lb, col->ub) )
   {
      /* dual row is free */
      return SCIPsetInfinity(set);
   }
   else
   {
      SCIP_Real redcost;

      /* calculate reduced costs */
      redcost = SCIPcolGetRedcost(col, stat, lp);
   
      if( !SCIPsetIsNegative(set, col->lb) )
      {
         /* dual row is  activity <= obj  <=>  redcost >= 0 */
         return redcost;
      }
      else if( SCIPsetIsPositive(set, col->ub) )
      {
         /* dual row is  activity == obj  <=>  redcost == 0 */
         return -REALABS(redcost);
      }
      else
      {
         /* dual row is  activity >= obj  <=>  redcost <= 0 */
         return -redcost;
      }
   }
}

/** calculates the farkas coefficient y^T A_i of a column i using the given dual farkas vector y */
SCIP_Real SCIPcolCalcFarkasCoef(
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real*            dualfarkas          /**< dense dual farkas vector for current LP rows */
   )
{
   SCIP_ROW* row;
   SCIP_Real farkas;
   int i;

   assert(col != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(dualfarkas != NULL);

   farkas = 0.0;
   for( i = 0; i < col->nlprows; ++i )
   {
      row = col->rows[i];
      assert(row != NULL);
      assert(row->lppos >= 0);
      farkas += col->vals[i] * dualfarkas[row->lppos];
   }

   if( col->nunlinked > 0 )
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->lppos == -1 || col->linkpos[i] == -1);
         if( row->lppos >= 0 )
            farkas += col->vals[i] * dualfarkas[row->lppos];
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
      }
   }
#endif

   return farkas;
}

/** gets the farkas coefficient y^T A_i of a column i in last LP (which must be infeasible) */
static
SCIP_Real colCalcInternalFarkasCoef(
   SCIP_COL*             col                 /**< LP column */
   )
{
   SCIP_ROW* row;
   SCIP_Real farkas;
   int i;

   assert(col != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);

   farkas = 0.0;
   for( i = 0; i < col->nlprows; ++i )
   {
      row = col->rows[i];
      assert(row != NULL);
      assert(row->dualfarkas < SCIP_INVALID);
      assert(row->lppos >= 0);
      assert(col->linkpos[i] >= 0);
      farkas += col->vals[i] * row->dualfarkas;
   }

   if( col->nunlinked > 0 )
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->lppos >= 0 || row->dualfarkas == 0.0);
         assert(row->lppos == -1 || col->linkpos[i] == -1);
         if( row->lppos >= 0 )
            farkas += col->vals[i] * row->dualfarkas;
      }
   }
#ifndef NDEBUG
   else
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->dualfarkas == 0.0);
         assert(row->lppos == -1);
         assert(col->linkpos[i] >= 0);
      }
   }
#endif

   return farkas;
}

/** gets the farkas coefficient of a column in last LP (which must be infeasible) */
SCIP_Real SCIPcolGetFarkasCoef(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(col != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(col->validfarkaslp <= stat->lpcount);
   assert(lp->validfarkaslp == stat->lpcount);

   if( col->validfarkaslp < stat->lpcount )
   {
      col->farkascoef = colCalcInternalFarkasCoef(col);
      col->validfarkaslp = stat->lpcount;
   }
   assert(col->validfarkaslp == stat->lpcount);
   assert(col->farkascoef < SCIP_INVALID);

   return col->farkascoef;
}

/** gets the farkas value of a column in last LP (which must be infeasible), i.e. the farkas coefficient y^T A_i times
 *  the best bound for this coefficient, i.e. max{y^T A_i x_i | lb <= x_i <= ub}
 */
SCIP_Real SCIPcolGetFarkasValue(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Real farkascoef;

   assert(col != NULL);

   farkascoef = SCIPcolGetFarkasCoef(col, stat, lp);

   if( farkascoef > 0.0 )
      return col->ub * farkascoef;
   else
      return col->lb * farkascoef;
}

/** actually performs strong branching on the given variable */
static
SCIP_RETCODE colStrongbranch(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   SCIP_RETCODE retcode;
   SCIP_Real sbdown;
   SCIP_Real sbup;
   SCIP_Bool sbdownvalid;
   SCIP_Bool sbupvalid;
   int iter;

   assert(col != NULL);
   assert(SCIPcolIsIntegral(col));
   assert(SCIPvarIsIntegral(col->var));
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->solved);
   assert(lp->flushed);
   assert(lperror != NULL);

   SCIPdebugMessage("performing strong branching on variable <%s>(%g) with %d iterations\n", 
      SCIPvarGetName(col->var), col->primsol, itlim);

   /* start timing */
   SCIPclockStart(stat->strongbranchtime, set);
      
   /* call LPI strong branching */
   col->sbitlim = itlim;
   col->nsbcalls++;
   retcode = SCIPlpiStrongbranch(lp->lpi, col->lpipos, col->primsol, itlim,
      &sbdown, &sbup, &sbdownvalid, &sbupvalid, &iter);
   
   /* check return code for errors */
   if( retcode == SCIP_LPERROR )
   {
      *lperror = TRUE;
      col->sbdown = SCIP_INVALID;
      col->sbup = SCIP_INVALID;
      col->sbdownvalid = FALSE;
      col->sbupvalid = FALSE;
      col->validsblp = -1;
      col->sbsolval = SCIP_INVALID;
      col->sblpobjval = SCIP_INVALID;
      col->sbnode = -1;
   }
   else
   {
      *lperror = FALSE;
      SCIP_CALL( retcode );

      col->sbdown = MIN(sbdown + lp->looseobjval, lp->cutoffbound);
      col->sbup = MIN(sbup + lp->looseobjval, lp->cutoffbound);
      col->sbdownvalid = sbdownvalid;
      col->sbupvalid = sbupvalid;

      /* update strong branching statistics */
      if( iter == -1 )
      {
         /* calculate avergate iteration number */
         iter = stat->ndualresolvelps > 0 ? (int)(2*stat->ndualresolvelpiterations / stat->ndualresolvelps)
            : stat->nduallps > 0 ? (int)((stat->nduallpiterations / stat->nduallps) / 5)
            : stat->nprimalresolvelps > 0 ? (int)(2*stat->nprimalresolvelpiterations / stat->nprimalresolvelps)
            : stat->nprimallps > 0 ? (int)((stat->nprimallpiterations / stat->nprimallps) / 5)
            : 0;
         if( iter/2 >= itlim )
            iter = 2*itlim;
      }
      stat->nstrongbranchs++;
      stat->nsblpiterations += iter;
      if( stat->nnodes == 1 )
      {
         stat->nrootstrongbranchs++;
         stat->nrootsblpiterations += iter;
      }
   }

   /* stop timing */
   SCIPclockStop(stat->strongbranchtime, set);

   return SCIP_OKAY;
}

/** gets strong branching information on a column variable */
SCIP_RETCODE SCIPcolGetStrongbranch(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPcolIsIntegral(col));
   assert(SCIPvarIsIntegral(col->var));
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(col->primsol < SCIP_INVALID);
   assert(col->lpipos >= 0);
   assert(col->lppos >= 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL);
   assert(lp->validsollp == stat->lpcount);
   assert(col->lppos < lp->ncols);
   assert(lp->cols[col->lppos] == col);
   assert(itlim >= 1);
   assert(down != NULL);
   assert(up != NULL);
   assert(lperror != NULL);

   *lperror = FALSE;

   if( col->validsblp != stat->lpcount || itlim > col->sbitlim )
   {
      col->validsblp = stat->lpcount;
      col->sbsolval = col->primsol;
      col->sblpobjval = SCIPlpGetObjval(lp, set);
      col->sbnode = stat->nnodes;

      /* if a loose variables has an infinite best bound, the LP bound is -infinity and no gain can be achieved */
      if( lp->looseobjvalinf > 0 )
      {
         col->sbdown = -SCIPsetInfinity(set);
         col->sbup = -SCIPsetInfinity(set);
         col->sbdownvalid = FALSE;
         col->sbupvalid = FALSE;
      }
      else
      {
         /* perform the strong branching on the column */
         SCIP_CALL( colStrongbranch(col, set, stat, lp, itlim, lperror) );
      }
   }
   assert(*lperror || col->sbdown < SCIP_INVALID);
   assert(*lperror || col->sbup < SCIP_INVALID);

   *down = col->sbdown;
   *up = col->sbup;
   if( down != NULL )
      *downvalid = col->sbdownvalid;
   if( upvalid != NULL )
      *upvalid = col->sbupvalid;

   return SCIP_OKAY;
}

/** gets last strong branching information available for a column variable;
 *  returns values of SCIP_INVALID, if strong branching was not yet called on the given column;
 *  keep in mind, that the returned old values may have nothing to do with the current LP solution
 */
void SCIPcolGetStrongbranchLast(
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real*            down,               /**< stores dual bound after branching column down, or NULL */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up, or NULL */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Real*            solval,             /**< stores LP solution value of column at last strong branching call, or NULL */
   SCIP_Real*            lpobjval            /**< stores LP objective value at last strong branching call, or NULL */
   )
{
   assert(col != NULL);

   if( down != NULL )
      *down = col->sbdown;
   if( up != NULL )
      *up = col->sbup;
   if( downvalid != NULL )
      *downvalid = col->sbdownvalid;
   if( upvalid != NULL )
      *upvalid = col->sbupvalid;
   if( solval != NULL )
      *solval = col->sbsolval;
   if( lpobjval != NULL )
      *lpobjval = col->sblpobjval;
}

/** if strong branching was already applied on the column at the current node, returns the number of LPs solved after
 *  the LP where the strong branching on this column was applied;
 *  if strong branching was not yet applied on the column at the current node, returns INT_MAX
 */
int SCIPcolGetStrongbranchLPAge(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(col != NULL);
   assert(stat != NULL);

   return (col->sbnode != stat->nnodes ? INT_MAX : stat->lpcount - col->validsblp);
}

/** output column to file stream */
void SCIPcolPrint(
   SCIP_COL*             col,                /**< LP column */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int r;

   assert(col != NULL);
   assert(col->var != NULL);

   /* print bounds */
   SCIPmessageFPrintInfo(file, "(obj: %g) [%g,%g], ", col->obj, col->lb, col->ub);

   /* print coefficients */
   if( col->len == 0 )
      SCIPmessageFPrintInfo(file, "<empty>");
   for( r = 0; r < col->len; ++r )
   {
      assert(col->rows[r] != NULL);
      assert(col->rows[r]->name != NULL);
      SCIPmessageFPrintInfo(file, "%+g<%s> ", col->vals[r], col->rows[r]->name);
   }
   SCIPmessageFPrintInfo(file, "\n");
}




/*
 * Row methods
 */

/** calculates row norms and min/maxidx from scratch, and checks for sortation */
static
void rowCalcNorms(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(row != NULL);
   assert(set != NULL);

   row->sqrnorm = 0.0;
   row->sumnorm = 0.0;
   row->objprod = 0.0;
   row->maxval = 0.0;
   row->nummaxval = 1;
   row->minval = SCIPsetInfinity(set);
   row->numminval = 1;
   row->minidx = INT_MAX;
   row->maxidx = INT_MIN;
   row->validminmaxidx = TRUE;
   row->lpcolssorted = TRUE;
   row->nonlpcolssorted = TRUE;

   /* check, if row is sorted
    * calculate sqrnorm, sumnorm, maxval, minval, minidx, and maxidx
    */
   for( i = 0; i < row->nlpcols; ++i )
   {
      assert(row->cols[i] != NULL);
      assert(!SCIPsetIsZero(set, row->vals[i]));
      assert(row->cols[i]->lppos >= 0);
      assert(row->linkpos[i] >= 0);
      assert(row->cols[i]->index == row->cols_index[i]);

      rowAddNorms(row, set, row->cols[i], row->vals[i]);
      if( i > 0 )
      {
         assert(row->cols[i-1]->index == row->cols_index[i-1]);
         row->lpcolssorted = row->lpcolssorted && (row->cols_index[i-1] < row->cols_index[i]);
      }
   }
   for( i = row->nlpcols; i < row->len; ++i )
   {
      assert(row->cols[i] != NULL);
      assert(!SCIPsetIsZero(set, row->vals[i]));
      assert(row->cols[i]->lppos == -1 || row->linkpos[i] == -1);
      assert(row->cols[i]->index == row->cols_index[i]);

      rowAddNorms(row, set, row->cols[i], row->vals[i]);
      if( i > row->nlpcols )
      {
         assert(row->cols[i-1]->index == row->cols_index[i-1]);
         row->nonlpcolssorted = row->nonlpcolssorted && (row->cols_index[i-1] < row->cols_index[i]);
      }
   }
}

/** checks, whether the given scalar scales the given value to an integral number with error in the given bounds */
static
SCIP_Bool isIntegralScalar(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real*            intval              /**< pointer to store the scaled integral value, or NULL */
   )
{
   SCIP_Real sval;
   SCIP_Real downval;
   SCIP_Real upval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = EPSFLOOR(sval, 0.0);
   upval = EPSCEIL(sval, 0.0);
   
   if( SCIPrelDiff(sval, downval) <= maxdelta )
   {
      if( intval != NULL )
         *intval = downval;
      return TRUE;
   }
   else if( SCIPrelDiff(sval, upval) >= mindelta )
   {
      if( intval != NULL )
         *intval = upval;
      return TRUE;
   }
   
   return FALSE;
}

/** scales row with given factor, and rounds coefficients to integers if close enough;
 *  the constant is automatically moved to the sides;
 *  if the row's activity is proven to be integral, the sides are automatically rounded to the next integer
 */
static
SCIP_RETCODE rowScale(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             scaleval,           /**< value to scale row with */
   SCIP_Bool             integralcontvars,   /**< should the coefficients of the continuous variables also be made integral,
                                              *   if they are close to integral values? */
   SCIP_Real             minrounddelta,      /**< minimal relative difference of scaled coefficient s*c and integral i,
                                              *   upto which the integral is used instead of the scaled real coefficient */
   SCIP_Real             maxrounddelta       /**< maximal relative difference of scaled coefficient s*c and integral i
                                              *   upto which the integral is used instead of the scaled real coefficient */
   )
{
   SCIP_COL* col;
   SCIP_Real val;
   SCIP_Real newval;
   SCIP_Real intval;
   SCIP_Real mindelta;
   SCIP_Real maxdelta;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool mindeltainf;
   SCIP_Bool maxdeltainf;
   int c;

   assert(row != NULL);
   assert(row->len == 0 || row->cols != NULL);
   assert(row->len == 0 || row->vals != NULL);
   assert(SCIPsetIsPositive(set, scaleval));
   assert(-1.0 < minrounddelta && minrounddelta <= 0.0);
   assert(0.0 <= maxrounddelta && maxrounddelta < 1.0);

   SCIPdebugMessage("scale row <%s> with %g (tolerance=[%g,%g])\n", row->name, scaleval, minrounddelta, maxrounddelta);

   mindelta = 0.0;
   maxdelta = 0.0;
   mindeltainf = FALSE;
   maxdeltainf = FALSE;

   /* scale the row coefficients, thereby recalculating whether the row's activity is always integral;
    * if the row coefficients are rounded to the nearest integer value, calculate the maximal activity difference,
    * this rounding can lead to
    */
   row->integral = TRUE;
   for( c = 0; c < row->len; ++c )
   {
      col = row->cols[c];
      val = row->vals[c];
      assert(!SCIPsetIsZero(set, val));

      /* get local or global bounds for column, depending on the local or global feasibility of the row */
      if( row->local )
      {
         lb = col->lb;
         ub = col->ub;
      }
      else
      {
         lb = SCIPvarGetLbGlobal(col->var);
         ub = SCIPvarGetUbGlobal(col->var);
      }

      /* calculate scaled coefficient */
      newval = val * scaleval;
      if( (integralcontvars || SCIPcolIsIntegral(col) || SCIPsetIsIntegral(set, newval))
         && isIntegralScalar(val, scaleval, minrounddelta, maxrounddelta, &intval) )
      {
         if( intval < newval )
         {
            mindelta += (intval - newval)*ub;
            maxdelta += (intval - newval)*lb;
            mindeltainf = mindeltainf || SCIPsetIsInfinity(set, ub);
            maxdeltainf = maxdeltainf || SCIPsetIsInfinity(set, -lb);
         }
         else
         {
            mindelta += (intval - newval)*lb;
            maxdelta += (intval - newval)*ub;
            mindeltainf = mindeltainf || SCIPsetIsInfinity(set, -lb);
            maxdeltainf = maxdeltainf || SCIPsetIsInfinity(set, ub);
         }
         newval = intval;
      }

      if( !SCIPsetIsEQ(set, val, newval) )
      {
         /* if column knows of the row, change the corresponding coefficient in the column */
         if( row->linkpos[c] >= 0 )
         {
            assert(col->rows[row->linkpos[c]] == row);
            assert(SCIPsetIsEQ(set, col->vals[row->linkpos[c]], row->vals[c]));
            SCIP_CALL( colChgCoefPos(col, set, lp, row->linkpos[c], newval) );
         }
            
         /* change the coefficient in the row, and update the norms and integrality status */
         SCIP_CALL( rowChgCoefPos(row, set, lp, c, newval) );
      }
      else
         row->integral = row->integral && SCIPcolIsIntegral(col) && SCIPsetIsIntegral(set, val);
   }

   /* scale the row sides, and move the constant to the sides; relax the sides with accumulated delta in order
    * to not destroy feasibility due to rounding
    */
   if( !SCIPsetIsInfinity(set, -row->lhs) )
   {
      if( mindeltainf )
         newval = -SCIPsetInfinity(set);
      else
      {
         newval = (row->lhs - row->constant) * scaleval + mindelta;
         if( SCIPsetIsIntegral(set, newval) || (row->integral && !row->modifiable) )
            newval = SCIPsetCeil(set, newval);
      }
      SCIP_CALL( SCIProwChgLhs(row, set, lp, newval) );
   }
   if( !SCIPsetIsInfinity(set, row->rhs) )
   {
      if( maxdeltainf )
         newval = SCIPsetInfinity(set);
      else
      {
         newval = (row->rhs - row->constant) * scaleval + maxdelta;
         if( SCIPsetIsIntegral(set, newval) || (row->integral && !row->modifiable) )
            newval = SCIPsetFloor(set, newval);
      }
      SCIP_CALL( SCIProwChgRhs(row, set, lp, newval) );
   }

   /* clear the row constant */
   SCIP_CALL( SCIProwChgConstant(row, set, stat, lp, 0.0) );

   SCIPdebugMessage("scaled row <%s> (integral: %d)\n", row->name, row->integral);
   SCIPdebug(SCIProwPrint(row, NULL));

#ifdef SCIP_DEBUG
   /* check integrality status of row */
   for( c = 0; c < row->len && SCIPcolIsIntegral(row->cols[c]) && SCIPsetIsIntegral(set, row->vals[c]); ++c )
   {}
   assert(row->integral == (c == row->len));
#endif

   return SCIP_OKAY;
}

/** creates and captures an LP row */
SCIP_RETCODE SCIProwCreate(
   SCIP_ROW**            row,                /**< pointer to LP row data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name,               /**< name of row */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COL**            cols,               /**< array with columns of row entries */
   SCIP_Real*            vals,               /**< array with coefficients of row entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   assert(row != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(len >= 0);
   assert(len == 0 || (cols != NULL && vals != NULL));
   assert(lhs <= rhs);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, row) );

   (*row)->integral = TRUE;
   if( len > 0 )
   {
      SCIP_VAR* var;
      int i;

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*row)->cols, cols, len) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*row)->vals, vals, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*row)->cols_index, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*row)->linkpos, len) );

      for( i = 0; i < len; ++i )
      {
         assert(cols[i] != NULL);
         assert(!SCIPsetIsZero(set, vals[i]));

         var = cols[i]->var;
         (*row)->cols_index[i] = cols[i]->index;
         (*row)->linkpos[i] = -1;
         (*row)->integral = (*row)->integral && SCIPvarIsIntegral(var) && SCIPsetIsIntegral(set, vals[i]);
      }
   }
   else
   {
      (*row)->cols = NULL;
      (*row)->cols_index = NULL;
      (*row)->vals = NULL;
      (*row)->linkpos = NULL;
   }
   
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*row)->name, name, strlen(name)+1) );
   (*row)->constant = 0.0;
   (*row)->lhs = lhs;
   (*row)->rhs = rhs;
   (*row)->flushedlhs = -SCIPsetInfinity(set);
   (*row)->flushedrhs = SCIPsetInfinity(set);
   (*row)->sqrnorm = 0.0;
   (*row)->sumnorm = 0.0;
   (*row)->objprod = 0.0;
   (*row)->maxval = 0.0;
   (*row)->minval = SCIPsetInfinity(set);
   (*row)->dualsol = 0.0;
   (*row)->activity = SCIP_INVALID;
   (*row)->dualfarkas = 0.0;
   (*row)->pseudoactivity = SCIP_INVALID;
   (*row)->minactivity = SCIP_INVALID;
   (*row)->maxactivity = SCIP_INVALID;
   (*row)->index = stat->nrowidx++;
   (*row)->size = len;
   (*row)->len = len;
   (*row)->nlpcols = 0;
   (*row)->nunlinked = len;
   (*row)->nuses = 0;
   (*row)->lppos = -1;
   (*row)->lpipos = -1;
   (*row)->lpdepth = -1;
   (*row)->minidx = INT_MAX;
   (*row)->maxidx = INT_MIN;
   (*row)->nummaxval = 0;
   (*row)->numminval = 0;
   (*row)->validactivitylp = -1;
   (*row)->validpsactivitydomchg = -1;
   (*row)->validactivitybdsdomchg = -1;
   (*row)->age = 0;
   (*row)->obsoletenode = -1;
   (*row)->basisstatus = SCIP_BASESTAT_BASIC; /*lint !e641*/
   (*row)->lpcolssorted = TRUE;
   (*row)->nonlpcolssorted = (len <= 1);
   (*row)->delaysort = FALSE;
   (*row)->validminmaxidx = FALSE;
   (*row)->lhschanged = FALSE;
   (*row)->rhschanged = FALSE;
   (*row)->coefchanged = FALSE;
   (*row)->local = local;
   (*row)->modifiable = modifiable;
   (*row)->nlocks = 0;
   (*row)->removeable = removeable;

   /* calculate row norms and min/maxidx, and check if row is sorted */
   rowCalcNorms(*row, set);

   /* capture the row */
   SCIProwCapture(*row);

   return SCIP_OKAY;
}

/** frees an LP row */
SCIP_RETCODE SCIProwFree(
   SCIP_ROW**            row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(blkmem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->nuses == 0);
   assert((*row)->lppos == -1);

   /* remove column indices from corresponding rows */
   SCIP_CALL( rowUnlink(*row, blkmem, set, lp) );

   BMSfreeBlockMemoryArray(blkmem, &(*row)->name, strlen((*row)->name)+1);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->cols, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->cols_index, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->vals, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->linkpos, (*row)->size);
   BMSfreeBlockMemory(blkmem, row);

   return SCIP_OKAY;
}

/** increases usage counter of LP row */
void SCIProwCapture(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->nuses >= 0);
   assert(row->nlocks <= (unsigned int)(row->nuses));

   SCIPdebugMessage("capture row <%s> with nuses=%d and nlocks=%d\n", row->name, row->nuses, row->nlocks);
   row->nuses++;
}

/** decreases usage counter of LP row, and frees memory if necessary */
SCIP_RETCODE SCIProwRelease(
   SCIP_ROW**            row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(blkmem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->nuses >= 1);
   assert((*row)->nlocks < (unsigned int)((*row)->nuses));

   SCIPdebugMessage("release row <%s> with nuses=%d and nlocks=%d\n", (*row)->name, (*row)->nuses, (*row)->nlocks);
   (*row)->nuses--;
   if( (*row)->nuses == 0 )
   {
      SCIP_CALL( SCIProwFree(row, blkmem, set, lp) );
   }

   *row = NULL;

   return SCIP_OKAY;
}

/** locks an unmodifiable row, which forbids further changes; has no effect on modifiable rows */
void SCIProwLock(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   /* check, if row is modifiable */
   if( !row->modifiable )
   {
      SCIPdebugMessage("lock row <%s> with nuses=%d and nlocks=%d\n", row->name, row->nuses, row->nlocks);
      row->nlocks++;
   }
}

/** unlocks a lock of an unmodifiable row; a row with no sealed lock may be modified; has no effect on modifiable rows */
void SCIProwUnlock(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   /* check, if row is modifiable */
   if( !row->modifiable )
   {
      SCIPdebugMessage("unlock row <%s> with nuses=%d and nlocks=%d\n", row->name, row->nuses, row->nlocks);
      assert(row->nlocks > 0);
      row->nlocks--;
   }
}

/** adds a previously non existing coefficient to an LP row */
SCIP_RETCODE SCIProwAddCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);

   SCIP_CALL( rowAddCoef(row, blkmem, set, lp, col, val, -1) );

   checkLinks(lp);

   return SCIP_OKAY;
}

/** deletes coefficient from row */
SCIP_RETCODE SCIProwDelCoef(
   SCIP_ROW*             row,                /**< row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(row != NULL);
   assert(!row->delaysort);
   assert(lp != NULL);
   assert(!lp->diving);
   assert(col != NULL);
   assert(col->var != NULL);

   /* search the position of the column in the row's col vector */
   pos = rowSearchCoef(row, col);
   if( pos == -1 )
   {
      SCIPerrorMessage("coefficient for column <%s> doesn't exist in row <%s>\n", SCIPvarGetName(col->var), row->name);
      return SCIP_INVALIDDATA;
   }
   assert(0 <= pos && pos < row->len);
   assert(row->cols[pos] == col);
   assert(row->cols_index[pos] == col->index);

   /* if column knows of the row, remove the row from the column's row vector */
   if( row->linkpos[pos] >= 0 )
   {
      assert(col->rows[row->linkpos[pos]] == row);
      assert(SCIPsetIsEQ(set, col->vals[row->linkpos[pos]], row->vals[pos]));
      SCIP_CALL( colDelCoefPos(col, set, lp, row->linkpos[pos]) );
   }

   /* delete the column from the row's col vector */
   SCIP_CALL( rowDelCoefPos(row, set, lp, pos) );
   
   checkLinks(lp);

   return SCIP_OKAY;
}

/** changes or adds a coefficient to an LP row */
SCIP_RETCODE SCIProwChgCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   int pos;

   assert(row != NULL);
   assert(!row->delaysort);
   assert(lp != NULL);
   assert(!lp->diving);
   assert(col != NULL);

   /* search the position of the column in the row's col vector */
   pos = rowSearchCoef(row, col);

   /* check, if column already exists in the row's col vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( rowAddCoef(row, blkmem, set, lp, col, val, -1) );
   }
   else
   {
      /* modifify already existing coefficient */
      assert(0 <= pos && pos < row->len);
      assert(row->cols[pos] == col);
      assert(row->cols_index[pos] == col->index);

      /* if column knows of the row, change the corresponding coefficient in the column */
      if( row->linkpos[pos] >= 0 )
      {
         assert(col->rows[row->linkpos[pos]] == row);
         assert(SCIPsetIsEQ(set, col->vals[row->linkpos[pos]], row->vals[pos]));
         SCIP_CALL( colChgCoefPos(col, set, lp, row->linkpos[pos], val) );
      }

      /* change the coefficient in the row */
      SCIP_CALL( rowChgCoefPos(row, set, lp, pos, val) );
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** increases value of an existing or nonexisting coefficient in an LP row */
SCIP_RETCODE SCIProwIncCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             incval              /**< value to add to the coefficient */
   )
{
   int pos;

   assert(row != NULL);
   assert(lp != NULL);
   assert(!lp->diving);
   assert(col != NULL);

   if( SCIPsetIsZero(set, incval) )
      return SCIP_OKAY;

   /* search the position of the column in the row's col vector */
   pos = rowSearchCoef(row, col);

   /* check, if column already exists in the row's col vector */
   if( pos == -1 )
   {
      /* coefficient doesn't exist, or sorting is delayed: add coefficient to the end of the row's arrays */
      SCIP_CALL( rowAddCoef(row, blkmem, set, lp, col, incval, -1) );
   }
   else
   {
      /* modifify already existing coefficient */
      assert(0 <= pos && pos < row->len);
      assert(row->cols[pos] == col);
      assert(row->cols_index[pos] == col->index);

      /* if column knows of the row, change the corresponding coefficient in the column */
      if( row->linkpos[pos] >= 0 )
      {
         assert(col->rows[row->linkpos[pos]] == row);
         assert(SCIPsetIsEQ(set, col->vals[row->linkpos[pos]], row->vals[pos]));
         SCIP_CALL( colChgCoefPos(col, set, lp, row->linkpos[pos], row->vals[pos] + incval) );
      }

      /* change the coefficient in the row */
      SCIP_CALL( rowChgCoefPos(row, set, lp, pos, row->vals[pos] + incval) );
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** changes constant value of a row */
SCIP_RETCODE SCIProwChgConstant(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             constant            /**< new constant value */
   )
{
   assert(row != NULL);
   assert(row->lhs <= row->rhs);
   assert(!SCIPsetIsInfinity(set, REALABS(constant)));
   assert(stat != NULL);
   assert(lp != NULL);
   assert(!lp->diving);

   if( !SCIPsetIsEQ(set, constant, row->constant) )
   {
      if( row->validpsactivitydomchg == stat->domchgcount )
      {
         assert(row->pseudoactivity < SCIP_INVALID);
         row->pseudoactivity += constant - row->constant;
      }
      if( row->validactivitybdsdomchg == stat->domchgcount )
      {
         assert(row->minactivity < SCIP_INVALID);
         assert(row->maxactivity < SCIP_INVALID);
         row->minactivity += constant - row->constant;
         row->maxactivity += constant - row->constant;
      }

      if( !SCIPsetIsInfinity(set, -row->lhs) )
      {
         SCIP_CALL( rowSideChanged(row, set, lp, SCIP_SIDETYPE_LEFT) );
      }
      if( !SCIPsetIsInfinity(set, row->rhs) )
      {
         SCIP_CALL( rowSideChanged(row, set, lp, SCIP_SIDETYPE_RIGHT) );
      }

      row->constant = constant;
   }

   return SCIP_OKAY;
}

/** add constant value to a row */
SCIP_RETCODE SCIProwAddConstant(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             addval              /**< constant value to add to the row */
   )
{
   assert(row != NULL);
   assert(row->lhs <= row->rhs);
   assert(!SCIPsetIsInfinity(set, REALABS(addval)));
   assert(stat != NULL);
   assert(lp != NULL);
   assert(!lp->diving);

   if( !SCIPsetIsZero(set, addval) )
   {
      SCIP_CALL( SCIProwChgConstant(row, set, stat, lp, row->constant + addval) );
   }

   return SCIP_OKAY;
}

/** changes left hand side of LP row */
SCIP_RETCODE SCIProwChgLhs(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   assert(row != NULL);
   assert(lp != NULL);
   assert(!lp->diving);

   if( !SCIPsetIsEQ(set, row->lhs, lhs) )
   {
      row->lhs = lhs;
      SCIP_CALL( rowSideChanged(row, set, lp, SCIP_SIDETYPE_LEFT) );
   }

   return SCIP_OKAY;
}

/** changes right hand side of LP row */
SCIP_RETCODE SCIProwChgRhs(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   assert(row != NULL);
   assert(lp != NULL);
   assert(!lp->diving);

   if( !SCIPsetIsEQ(set, row->rhs, rhs) )
   {
      row->rhs = rhs;
      SCIP_CALL( rowSideChanged(row, set, lp, SCIP_SIDETYPE_RIGHT) );
   }

   return SCIP_OKAY;
}

/** additional scalars that are tried in integrality scaling */
static const SCIP_Real scalars[] = {3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0};
static const int nscalars = 9;

/** tries to find a value, such that all row coefficients, if scaled with this value become integral */
SCIP_RETCODE SCIProwCalcIntegralScalar(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal allowed scalar */
   SCIP_Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   SCIP_Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   SCIP_Bool*            success             /**< stores whether returned value is valid */
   )
{
   SCIP_COL* col;
   SCIP_Longint gcd;
   SCIP_Longint scm;
   SCIP_Longint nominator;
   SCIP_Longint denominator;
   SCIP_Real val;
   SCIP_Real absval;
   SCIP_Real minval;
   SCIP_Real scaleval;
   SCIP_Real twomultval;
   SCIP_Bool scalable;
   SCIP_Bool twomult;
   SCIP_Bool rational;
   int c;
   int s;

   /**@todo call misc.c:SCIPcalcIntegralScalar() instead - if usecontvars == FALSE, filter the integer variables first */
   assert(row != NULL);
   assert(row->len == 0 || row->cols != NULL);
   assert(row->len == 0 || row->cols_index != NULL);
   assert(row->len == 0 || row->vals != NULL);
   assert(maxdnom >= 1);
   assert(mindelta < 0.0);
   assert(maxdelta > 0.0);
   assert(success != NULL);

   SCIPdebugMessage("trying to find rational representation for row <%s> (contvars: %d)\n", SCIProwGetName(row), usecontvars);
   SCIPdebug( val = 0; ); /* avoid warning "val might be used uninitialized; see SCIPdebugMessage lastval=%g below */

   if( intscalar != NULL )
      *intscalar = SCIP_INVALID;
   *success = FALSE;

   /* get minimal absolute non-zero value */
   minval = SCIP_REAL_MAX;
   for( c = 0; c < row->len; ++c )
   {
      col = row->cols[c];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(col->var) == col);
      val = row->vals[c];
      assert(!SCIPsetIsZero(set, val));

      if( val < mindelta || val > maxdelta )
      {
         absval = REALABS(val);
         minval = MIN(minval, absval);
      }
   }
   if( minval == SCIP_REAL_MAX )
   {
      /* all coefficients are zero (inside tolerances) */
      if( intscalar != NULL )
         *intscalar = 1.0;
      *success = TRUE;
      SCIPdebugMessage(" -> all values are zero (inside tolerances)\n");

      return SCIP_OKAY;
   }
   assert(minval > MIN(-mindelta, maxdelta)); 
   assert(SCIPsetIsPositive(set, minval));
   assert(!SCIPsetIsInfinity(set, minval));

   /* try, if row coefficients can be made integral by multiplying them with the reciprocal of the smallest coefficient and
    * a power of 2
    */
   scalable = TRUE;
   scaleval = 1.0/minval;
   scalable = (scaleval <= maxscale);
   for( c = 0; c < row->len && scalable; ++c )
   {
      /* don't look at continuous variables, if we don't have to */
      if( !usecontvars && !SCIPcolIsIntegral(row->cols[c]) )
         continue;

      /* check, if the coefficient can be scaled with a simple scalar */
      val = row->vals[c];
      absval = REALABS(val);
      while( scaleval <= maxscale
         && (absval * scaleval < 0.5 || !isIntegralScalar(val, scaleval, mindelta, maxdelta, NULL)) )
      {
         for( s = 0; s < nscalars; ++s )
         {
            if( isIntegralScalar(val, scaleval * scalars[s], mindelta, maxdelta, NULL) )
            {
               scaleval *= scalars[s];
               break;
            }
         }
         if( s >= nscalars )
            scaleval *= 2.0;
      }
      scalable = (scaleval <= maxscale);
      SCIPdebugMessage(" -> val=%g, scaleval=%g, val*scaleval=%g, scalable=%d\n", 
         val, scaleval, val*scaleval, scalable);
   }
   if( scalable )
   {
      /* make row coefficients integral by dividing them by the smallest coefficient
       * (and multiplying them with a power of 2)
       */
      assert(scaleval <= maxscale);
      if( intscalar != NULL )
         *intscalar = scaleval;
      *success = TRUE;
      SCIPdebugMessage(" -> integrality can be achieved by scaling with %g (minval=%g)\n", scaleval, minval);

      return SCIP_OKAY;
   }

   /* try, if row coefficients can be made integral by multiplying them by a power of 2 */
   twomult = TRUE;
   twomultval = 1.0;
   twomult = (twomultval <= maxscale);
   for( c = 0; c < row->len && twomult; ++c )
   {
      /* don't look at continuous variables, if we don't have to */
      if( !usecontvars && !SCIPcolIsIntegral(row->cols[c]) )
         continue;

      /* check, if the coefficient can be scaled with a simple scalar */
      val = row->vals[c];
      absval = REALABS(val);
      while( twomultval <= maxscale
         && (absval * twomultval < 0.5 || !isIntegralScalar(val, twomultval, mindelta, maxdelta, NULL)) )
      {
         for( s = 0; s < nscalars; ++s )
         {
            if( isIntegralScalar(val, twomultval * scalars[s], mindelta, maxdelta, NULL) )
            {
               twomultval *= scalars[s];
               break;
            }
         }
         if( s >= nscalars )
            twomultval *= 2.0;
      }
      twomult = (twomultval <= maxscale);
      SCIPdebugMessage(" -> val=%g, twomult=%g, val*twomult=%g, twomultable=%d\n",
         val, twomultval, val*twomultval, twomult);
   }
   if( twomult )
   {
      /* make row coefficients integral by multiplying them with a power of 2 */
      assert(twomultval <= maxscale);
      if( intscalar != NULL )
         *intscalar = twomultval;
      *success = TRUE;
      SCIPdebugMessage(" -> integrality can be achieved by scaling with %g (power of 2)\n", twomultval);

      return SCIP_OKAY;
   }

   /* convert each coefficient into a rational number, calculate the greatest common divisor of the nominators
    * and the smallest common multiple of the denominators
    */
   gcd = 1;
   scm = 1;
   rational = TRUE;

   /* first coefficient (to initialize gcd) */
   for( c = 0; c < row->len && rational; ++c )
   {
      if( usecontvars || SCIPcolIsIntegral(row->cols[c]) )
      {
         val = row->vals[c];
         rational = SCIPrealToRational(val, mindelta, maxdelta, maxdnom, &nominator, &denominator);
         if( rational && nominator != 0 )
         {
            assert(denominator > 0);
            gcd = ABS(nominator);
            scm = denominator;
            rational = ((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
            SCIPdebugMessage(" -> first rational: val: %g == %lld/%lld, gcd=%lld, scm=%lld, rational=%d\n",
               val, nominator, denominator, gcd, scm, rational);
            break;
         }
      }
   }

   /* remaining coefficients */
   for( ++c; c < row->len && rational; ++c )
   {
      if( usecontvars || SCIPcolIsIntegral(row->cols[c]) )
      {
         val = row->vals[c];
         rational = SCIPrealToRational(val, mindelta, maxdelta, maxdnom, &nominator, &denominator);
         if( rational && nominator != 0 )
         {
            assert(denominator > 0);
            gcd = SCIPcalcGreComDiv(gcd, ABS(nominator));
            scm *= denominator / SCIPcalcGreComDiv(scm, denominator);
            rational = ((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
            SCIPdebugMessage(" -> next rational : val: %g == %lld/%lld, gcd=%lld, scm=%lld, rational=%d\n",
               val, nominator, denominator, gcd, scm, rational);
         }
      }
   }

   if( rational )
   {
      /* make row coefficients integral by multiplying them with the smallest common multiple of the denominators */
      assert((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
      if( intscalar != NULL )
         *intscalar = (SCIP_Real)scm/(SCIP_Real)gcd;
      *success = TRUE;
      SCIPdebugMessage(" -> integrality can be achieved by scaling with %g (rational:%lld/%lld)\n", 
         (SCIP_Real)scm/(SCIP_Real)gcd, scm, gcd);
   }
   else
   {
      assert(!(*success));
      SCIPdebugMessage(" -> rationalizing failed: gcd=%lld, scm=%lld, lastval=%g\n", gcd, scm, val); /*lint !e771*/
   }

   return SCIP_OKAY;
}

/** tries to scale row, s.t. all coefficients become integral */
SCIP_RETCODE SCIProwMakeIntegral(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal value to scale row with */
   SCIP_Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   SCIP_Bool*            success             /**< stores whether row could be made rational */
   )
{
   SCIP_Real intscalar;

   assert(success != NULL);

   /* calculate scalar to make coefficients integral */
   SCIP_CALL( SCIProwCalcIntegralScalar(row, set, mindelta, maxdelta, maxdnom, maxscale, usecontvars,
         &intscalar, success) );

   if( *success )
   {
      /* scale the row */
      SCIP_CALL( rowScale(row, set, stat, lp, intscalar, usecontvars, mindelta, maxdelta) );
   }

   return SCIP_OKAY;
}

/** sorts row entries such that LP columns precede non-LP columns and inside both parts lower column indices precede
 *  higher ones
 */
void SCIProwSort(
   SCIP_ROW*             row                 /**< row to be sorted */
   )
{
   assert(row != NULL);

   /* sort LP columns */
   rowSortLP(row);

   /* sort non-LP columns */
   rowSortNonLP(row);
}

/** sorts row, and merges equal column entries (resulting from lazy sorting and adding) into a single entry;
 *  removes zero entries from row
 *  the row must not be linked to the columns; otherwise, we would need to update the columns as well, which
 *  is too expensive
 */
static
void rowMerge(
   SCIP_ROW*             row,                /**< row to be sorted */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);
   assert(!row->delaysort);
   assert(row->nunlinked == row->len);
   assert(row->nlpcols == 0);

   SCIPdebugMessage("merging row <%s>\n", row->name);

   /* do nothing on empty rows; if row is sorted, nothing has to be done */
   if( row->len > 0 && (!row->lpcolssorted || !row->nonlpcolssorted) )
   {
      SCIP_COL** cols;
      int* cols_index;
      SCIP_Real* vals;
      int s;
      int t;
      
      /* make sure, the row is sorted */
      SCIProwSort(row);
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
      assert(!SCIPsetIsZero(set, vals[0]));
      assert(row->linkpos[0] == -1);

      for( s = 1; s < row->len; ++s )
      {
         assert(!SCIPsetIsZero(set, vals[s]));
         assert(row->linkpos[s] == -1);

         if( cols[s] == cols[t] )
         {
            /* merge entries with equal column */
            vals[t] += vals[s];
         }
         else
         {
            /* go to the next entry, overwriting current entry if coefficient is zero */
            if( !SCIPsetIsZero(set, vals[t]) )
            {
               row->integral = row->integral && SCIPcolIsIntegral(cols[t]) && SCIPsetIsIntegral(set, vals[t]);
               t++;
            }
            cols[t] = cols[s];
            cols_index[t] = cols_index[s];
            vals[t] = vals[s];
         }
      }
      if( !SCIPsetIsZero(set, vals[t]) )
      {
         row->integral = row->integral && SCIPcolIsIntegral(cols[t]) && SCIPsetIsIntegral(set, vals[t]);
         t++;
      }
      assert(s == row->len);
      assert(t <= row->len);

      row->len = t;
      row->nunlinked = t;

      /* if equal entries were merged, we have to recalculate the norms, since the squared euclidean norm is wrong */
      if( t < s )
         rowCalcNorms(row, set);
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
void SCIProwDelaySort(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(!row->delaysort);

   row->delaysort = TRUE;
}

/** disables delaying of row sorting, sorts row and merges coefficients with equal columns */
void SCIProwForceSort(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);
   assert(row->delaysort);

   row->delaysort = FALSE;
   rowMerge(row, set);
}

/** recalculates the current activity of a row */
void SCIProwRecalcLPActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_COL* col;
   int c;

   assert(row != NULL);
   assert(stat != NULL);

   row->activity = row->constant;
   for( c = 0; c < row->nlpcols; ++c )
   {
      col = row->cols[c];
      assert(col != NULL);
      assert(col->primsol < SCIP_INVALID);
      assert(col->lppos >= 0);
      assert(row->linkpos[c] >= 0);
      row->activity += row->vals[c] * col->primsol;
   }

   if( row->nunlinked > 0 )
   {
      for( c = row->nlpcols; c < row->len; ++c )
      {
         col = row->cols[c];
         assert(col != NULL);
         assert(col->lppos >= 0 || col->primsol == 0.0);
         assert(col->lppos == -1 || row->linkpos[c] == -1);
         if( col->lppos >= 0 )
            row->activity += row->vals[c] * col->primsol;
      }
   }
#ifndef NDEBUG
   else
   {
      for( c = row->nlpcols; c < row->len; ++c )
      {
         col = row->cols[c];
         assert(col != NULL);
         assert(col->primsol == 0.0);
         assert(col->lppos == -1);
         assert(row->linkpos[c] >= 0);
      }
   }
#endif

   row->validactivitylp = stat->lpcount;
}

/** returns the activity of a row in the current LP solution */
SCIP_Real SCIProwGetLPActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(row != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(row->validactivitylp <= stat->lpcount);
   assert(lp->validsollp == stat->lpcount);

   if( row->validactivitylp != stat->lpcount )
      SCIProwRecalcLPActivity(row, stat);
   assert(row->validactivitylp == stat->lpcount);
   assert(row->activity < SCIP_INVALID);

   return row->activity;
}

/** returns the feasibility of a row in the current LP solution: negative value means infeasibility */
SCIP_Real SCIProwGetLPFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Real activity;

   assert(row != NULL);

   activity = SCIProwGetLPActivity(row, stat, lp);

   return MIN(row->rhs - activity, activity - row->lhs);
}

/** calculates the current pseudo activity of a row */
void SCIProwRecalcPseudoActivity(
   SCIP_ROW*             row,                /**< row data */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_COL* col;
   int i;

   assert(row != NULL);
   assert(stat != NULL);

   row->pseudoactivity = row->constant;
   for( i = 0; i < row->len; ++i )
   {
      col = row->cols[i];
      assert(col != NULL);
      assert((i < row->nlpcols) == (row->linkpos[i] >= 0 && col->lppos >= 0));
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);

      row->pseudoactivity += SCIPcolGetBestBound(col) * row->vals[i];
   }
   row->validpsactivitydomchg = stat->domchgcount;
   assert(!row->integral || EPSISINT(row->pseudoactivity - row->constant, SCIP_DEFAULT_SUMEPSILON));
}

/** returns the pseudo activity of a row in the current pseudo solution */
SCIP_Real SCIProwGetPseudoActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(row != NULL);
   assert(stat != NULL);
   assert(row->validpsactivitydomchg <= stat->domchgcount);

   /* check, if activity bounds has to be calculated */
   if( row->validpsactivitydomchg != stat->domchgcount )
      SCIProwRecalcPseudoActivity(row, stat);
   assert(row->validpsactivitydomchg == stat->domchgcount);
   assert(row->pseudoactivity < SCIP_INVALID);

   return row->pseudoactivity;
}

/** returns the pseudo feasibility of a row in the current pseudo solution: negative value means infeasibility */
SCIP_Real SCIProwGetPseudoFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_Real pseudoactivity;

   assert(row != NULL);

   pseudoactivity = SCIProwGetPseudoActivity(row, stat);

   return MIN(row->rhs - pseudoactivity, pseudoactivity - row->lhs);
}

/** returns the activity of a row for a given solution */
SCIP_RETCODE SCIProwGetSolActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            solactivity         /**< pointer to store the row's activity for the solution */
   )
{
   SCIP_COL* col;
   SCIP_Real inf;
   int i;

   assert(row != NULL);
   assert(solactivity != NULL);

   *solactivity = row->constant;
   for( i = 0; i < row->len; ++i )
   {
      col = row->cols[i];
      assert(col != NULL);
      assert((i < row->nlpcols) == (row->linkpos[i] >= 0 && col->lppos >= 0));
      (*solactivity) += row->vals[i] * SCIPsolGetVal(sol, set, stat, col->var);
   }

   inf = SCIPsetInfinity(set);
   *solactivity = MAX(*solactivity, -inf);
   *solactivity = MIN(*solactivity, +inf);

   return SCIP_OKAY;
}

/** returns the feasibility of a row for the given solution */
SCIP_RETCODE SCIProwGetSolFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            solfeasibility      /**< pointer to store the row's feasibility for the solution */
   )
{
   SCIP_Real solactivity;

   assert(row != NULL);
   assert(solfeasibility != NULL);

   SCIP_CALL( SCIProwGetSolActivity(row, set, stat, sol, &solactivity) );

   *solfeasibility = MIN(row->rhs - solactivity, solactivity - row->lhs);

   return SCIP_OKAY;
}

/** calculates minimal and maximal activity of row w.r.t. the column's bounds */
static
void rowCalcActivityBounds(
   SCIP_ROW*             row,                /**< row data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   SCIP_COL* col;
   SCIP_Real val;
   SCIP_Bool mininfinite;
   SCIP_Bool maxinfinite;
   int i;
   
   assert(row != NULL);
   assert(!SCIPsetIsInfinity(set, REALABS(row->constant)));
   assert(stat != NULL);
   
   /* calculate activity bounds */
   mininfinite = FALSE;
   maxinfinite = FALSE;
   row->minactivity = row->constant;
   row->maxactivity = row->constant;
   for( i = 0; i < row->len && (!mininfinite || !maxinfinite); ++i )
   {
      col = row->cols[i];
      assert(col != NULL);
      assert((i < row->nlpcols) == (row->linkpos[i] >= 0 && col->lppos >= 0));
      val = row->vals[i];
      if( val >= 0.0 )
      {
         mininfinite = mininfinite || SCIPsetIsInfinity(set, -col->lb);
         maxinfinite = maxinfinite || SCIPsetIsInfinity(set, col->ub);
         if( !mininfinite )
            row->minactivity += val * col->lb;
         if( !maxinfinite )
            row->maxactivity += val * col->ub;
      }
      else
      {
         mininfinite = mininfinite || SCIPsetIsInfinity(set, col->ub);
         maxinfinite = maxinfinite || SCIPsetIsInfinity(set, -col->lb);
         if( !mininfinite )
            row->minactivity += val * col->ub;
         if( !maxinfinite )
            row->maxactivity += val * col->lb;
      }
   }

   if( mininfinite )
      row->minactivity = -SCIPsetInfinity(set);
   if( maxinfinite )
      row->maxactivity = SCIPsetInfinity(set);
   row->validactivitybdsdomchg = stat->domchgcount;

   assert(!row->integral || mininfinite || EPSISINT(row->minactivity - row->constant, SCIP_DEFAULT_SUMEPSILON));
   assert(!row->integral || maxinfinite || EPSISINT(row->maxactivity - row->constant, SCIP_DEFAULT_SUMEPSILON));
}

/** returns the minimal activity of a row w.r.t. the column's bounds */
SCIP_Real SCIProwGetMinActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(row != NULL);
   assert(stat != NULL);
   assert(row->validactivitybdsdomchg <= stat->domchgcount);

   /* check, if activity bounds has to be calculated */
   if( row->validactivitybdsdomchg != stat->domchgcount )
      rowCalcActivityBounds(row, set, stat);
   assert(row->validactivitybdsdomchg == stat->domchgcount);
   assert(row->minactivity < SCIP_INVALID);
   assert(row->maxactivity < SCIP_INVALID);

   return row->minactivity;
}
   
/** returns the maximal activity of a row w.r.t. the column's bounds */
SCIP_Real SCIProwGetMaxActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(row != NULL);
   assert(stat != NULL);
   assert(row->validactivitybdsdomchg <= stat->domchgcount);

   /* check, if activity bounds has to be calculated */
   if( row->validactivitybdsdomchg != stat->domchgcount )
      rowCalcActivityBounds(row, set, stat);
   assert(row->validactivitybdsdomchg == stat->domchgcount);
   assert(row->minactivity < SCIP_INVALID);
   assert(row->maxactivity < SCIP_INVALID);

   return row->maxactivity;
}
   
/** gets maximal absolute value of row vector coefficients */
SCIP_Real SCIProwGetMaxval(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);
   
   if( row->nummaxval == 0 )
      rowCalcNorms(row, set);
   assert(row->nummaxval > 0);
   assert(row->maxval >= 0.0);

   return row->maxval;
}

/** gets minimal absolute value of row vector's non-zero coefficients */
SCIP_Real SCIProwGetMinval(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);
   
   if( row->numminval == 0 )
      rowCalcNorms(row, set);
   assert(row->numminval >= 0);
   assert(row->minval >= 0.0);

   return row->minval;
}

/** returns row's efficacy with respect to the current LP solution: e = -feasibility/norm */
SCIP_Real SCIProwGetEfficacy(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Real norm;
   SCIP_Real feasibility;
   SCIP_Real eps;

   assert(set != NULL);

   switch( set->sepa_efficacynorm )
   {
   case 'e':
      norm = SCIProwGetNorm(row);
      break;
   case 'm':
      norm = SCIProwGetMaxval(row, set);
      break;
   case 's':
      norm = SCIProwGetSumNorm(row);
      break;
   case 'd':
      norm = (row->len == 0 ? 0.0 : 1.0);
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", set->sepa_efficacynorm);
      SCIPABORT();
      norm = 0.0;
   }

   eps = SCIPsetSumepsilon(set);
   norm = MAX(norm, eps);
   feasibility = SCIProwGetLPFeasibility(row, stat, lp);

   return -feasibility / norm;
}

/** returns whether the row's efficacy with respect to the current LP solution is greater than the minimal cut efficacy */
SCIP_Bool SCIProwIsEfficacious(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             root                /**< should the root's minimal cut efficacy be used? */
   )
{
   SCIP_Real efficacy;

   efficacy = SCIProwGetEfficacy(row, set, stat, lp);

   return SCIPsetIsEfficacious(set, root, efficacy);
}

/** returns the scalar product of the coefficient vectors of the two given rows */
SCIP_Real SCIProwGetScalarProduct(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2                /**< second LP row */
   )
{
   SCIP_Real scalarprod;
   int i1;
   int i2;

   assert(row1 != NULL);
   assert(row2 != NULL);

   /* make sure, the rows are sorted */
   SCIProwSort(row1);
   assert(row1->lpcolssorted);
   assert(row1->nonlpcolssorted);
   SCIProwSort(row2);
   assert(row2->lpcolssorted);
   assert(row2->nonlpcolssorted);

   /* calculate the scalar product */
   scalarprod = 0.0;
   i1 = 0;
   i2 = 0;
   while( i1 < row1->len && i2 < row2->len )
   {
      assert(row1->cols[i1]->index == row1->cols_index[i1]);
      assert(row2->cols[i2]->index == row2->cols_index[i2]);
      assert((row1->cols[i1] == row2->cols[i2]) == (row1->cols_index[i1] == row2->cols_index[i2]));
      if( row1->cols_index[i1] < row2->cols_index[i2] )
         i1++;
      else if( row1->cols_index[i1] > row2->cols_index[i2] )
         i2++;
      else
      {
         scalarprod += row1->vals[i1] * row2->vals[i2];
         i1++;
         i2++;
      }
   }

   return scalarprod;
}

/** returns the degree of parallelism between the hyperplanes defined by the two row vectors v, w:
 *  p = |v*w|/(|v|*|w|);
 *  the hyperplanes are parellel, iff p = 1, they are orthogonal, iff p = 0
 */
SCIP_Real SCIProwGetParallelism(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2                /**< second LP row */
   )
{
   SCIP_Real scalarprod;

   scalarprod = SCIProwGetScalarProduct(row1, row2);
   
   return (REALABS(scalarprod) / (SCIProwGetNorm(row1) * SCIProwGetNorm(row2)));
}

/** returns the degree of orthogonality between the hyperplanes defined by the two row vectors v, w:
 *  o = 1 - |v*w|/(|v|*|w|);
 *  the hyperplanes are orthogonal, iff p = 1, they are parallel, iff p = 0
 */
SCIP_Real SCIProwGetOrthogonality(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2                /**< second LP row */
   )
{
   return 1.0 - SCIProwGetParallelism(row1, row2);
}

/** gets parallelism of row with objective function: if the returned value is 1, the row is parellel to the objective
 *  function, if the value is 0, it is orthogonal to the objective function
 */
SCIP_Real SCIProwGetObjParallelism(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Real prod;
   SCIP_Real parallelism;

   assert(row != NULL);
   assert(lp != NULL);

   prod = row->sqrnorm * lp->objsqrnorm;

   parallelism = SCIPsetIsPositive(set, prod) ? REALABS(row->objprod) / SQRT(prod) : 0.0;
   assert(SCIPsetIsGE(set, parallelism, 0.0));
   assert(SCIPsetIsLE(set, parallelism, 1.0));

   return parallelism;
}

/** output row to file stream */
void SCIProwPrint(
   SCIP_ROW*             row,                /**< LP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int i;

   assert(row != NULL);

   /* print left hand side */
   SCIPmessageFPrintInfo(file, "%g <= ", row->lhs);

   /* print coefficients */
   if( row->len == 0 )
      SCIPmessageFPrintInfo(file, "0 ");
   for( i = 0; i < row->len; ++i )
   {
      assert(row->cols[i] != NULL);
      assert(row->cols[i]->var != NULL);
      assert(SCIPvarGetName(row->cols[i]->var) != NULL);
      assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
      SCIPmessageFPrintInfo(file, "%+g<%s> ", row->vals[i], SCIPvarGetName(row->cols[i]->var));
   }

   /* print constant */
   if( REALABS(row->constant) > SCIP_DEFAULT_EPSILON )
      SCIPmessageFPrintInfo(file, "%+g ", row->constant);

   /* print right hand side */
   SCIPmessageFPrintInfo(file, "<= %g\n", row->rhs);
}




/*
 * LP solver data update
 */

/** resets column data to represent a column not in the LP solver */
static
void markColDeleted(
   SCIP_COL*             col                 /**< column to be marked deleted */
   )
{
   assert(col != NULL);

   col->lpipos = -1;
   col->primsol = 0.0;
   col->redcost = SCIP_INVALID;
   col->farkascoef = SCIP_INVALID;
   col->sbdown = SCIP_INVALID;
   col->sbup = SCIP_INVALID;
   col->sbdownvalid = FALSE;
   col->sbupvalid = FALSE;
   col->validredcostlp = -1;
   col->validfarkaslp = -1;
   col->sbitlim = -1;
   col->basisstatus = SCIP_BASESTAT_ZERO; /*lint !e641*/
}

/** applies all cached column removals to the LP solver */
static
SCIP_RETCODE lpFlushDelCols(
   SCIP_LP*              lp                  /**< current LP data */
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

      assert(!lp->diving);
      SCIPdebugMessage("flushing col deletions: shrink LP from %d to %d colums\n", lp->nlpicols, lp->lpifirstchgcol);
      SCIP_CALL( SCIPlpiDelCols(lp->lpi, lp->lpifirstchgcol, lp->nlpicols-1) );
      for( i = lp->lpifirstchgcol; i < lp->nlpicols; ++i )
      {
         markColDeleted(lp->lpicols[i]);
      }
      lp->nlpicols = lp->lpifirstchgcol;
      lp->flushdeletedcols = TRUE;

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }
   assert(lp->nlpicols == lp->lpifirstchgcol);

   return SCIP_OKAY;
}

/** applies all cached column additions to the LP solver */
static
SCIP_RETCODE lpFlushAddCols(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   int* beg;
   int* ind;
   SCIP_Real* val;
   char** name;
   SCIP_COL* col;
   SCIP_Real infinity;
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
   assert(!lp->diving);
   assert(lp->ncols > lp->nlpicols);
   SCIP_CALL( ensureLpicolsSize(lp, set, lp->ncols) );

   /* get the solver's infinity value */
   infinity = SCIPlpiInfinity(lp->lpi);

   /* count the (maximal) number of added coefficients, calculate the number of added columns */
   naddcols = lp->ncols - lp->nlpicols;
   naddcoefs = 0;
   for( c = lp->nlpicols; c < lp->ncols; ++c )
      naddcoefs += lp->cols[c]->len;
   assert(naddcols > 0);

   /* get temporary memory for changes */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &obj, naddcols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lb, naddcols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ub, naddcols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &beg, naddcols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ind, naddcoefs) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &val, naddcoefs) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &name, naddcols) );
   
   /* fill temporary memory with column data */
   nnonz = 0;
   for( pos = 0, c = lp->nlpicols; c < lp->ncols; ++pos, ++c )
   {
      col = lp->cols[c];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(col->var) == col);
      assert(col->lppos == c);
      assert(nnonz + col->nlprows <= naddcoefs);

      SCIPdebugMessage("flushing added column <%s>: ", SCIPvarGetName(col->var));
      SCIPdebug( SCIPcolPrint(col, NULL) );

      /* Because the column becomes a member of the LP solver, it now can take values
       * different from zero. That means, we have to include the column in the corresponding
       * row vectors.
       */
      SCIP_CALL( colLink(col, blkmem, set, lp) );

      lp->lpicols[c] = col;
      col->lpipos = c;
      col->primsol = SCIP_INVALID;
      col->redcost = SCIP_INVALID;
      col->farkascoef = SCIP_INVALID;
      col->sbdown = SCIP_INVALID;
      col->sbup = SCIP_INVALID;
      col->sbdownvalid = FALSE;
      col->sbupvalid = FALSE;
      col->validredcostlp = -1;
      col->validfarkaslp = -1;
      col->sbitlim = -1;
      col->objchanged = FALSE;
      col->lbchanged = FALSE;
      col->ubchanged = FALSE;
      col->coefchanged = FALSE;
      obj[pos] = col->obj;
      if( SCIPsetIsInfinity(set, -col->lb) )
         lb[pos] = -infinity;
      else
         lb[pos] = col->lb;
      if( SCIPsetIsInfinity(set, col->ub) )
         ub[pos] = infinity;
      else
         ub[pos] = col->ub;
      beg[pos] = nnonz;
      name[pos] = (char*)SCIPvarGetName(col->var);

      col->flushedobj = obj[pos];
      col->flushedlb = lb[pos];
      col->flushedub = ub[pos];

      for( i = 0; i < col->nlprows; ++i )
      {
         assert(col->rows[i] != NULL);
         lpipos = col->rows[i]->lpipos;
         if( lpipos >= 0 )
         {
            assert(lpipos < lp->nrows);
            assert(nnonz < naddcoefs);
            ind[nnonz] = lpipos;
            val[nnonz] = col->vals[i];
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
   SCIPdebugMessage("flushing col additions: enlarge LP from %d to %d colums\n", lp->nlpicols, lp->ncols);
   SCIP_CALL( SCIPlpiAddCols(lp->lpi, naddcols, obj, lb, ub, name, nnonz, beg, ind, val) );
   lp->nlpicols = lp->ncols;
   lp->lpifirstchgcol = lp->nlpicols;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &name);
   SCIPsetFreeBufferArray(set, &val);
   SCIPsetFreeBufferArray(set, &ind);
   SCIPsetFreeBufferArray(set, &beg);
   SCIPsetFreeBufferArray(set, &ub);
   SCIPsetFreeBufferArray(set, &lb);
   SCIPsetFreeBufferArray(set, &obj);

   lp->flushaddedcols = TRUE;

   /* mark the LP unsolved */
   lp->solved = FALSE;
   lp->dualfeasible = FALSE;
   lp->lpobjval = SCIP_INVALID;
   lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

   return SCIP_OKAY;
}

/** resets row data to represent a row not in the LP solver */
static
void markRowDeleted(
   SCIP_ROW*             row                 /**< row to be marked deleted */
   )
{
   assert(row != NULL);

   row->lpipos = -1;
   row->dualsol = 0.0;
   row->activity = SCIP_INVALID;
   row->dualfarkas = 0.0;
   row->basisstatus = SCIP_BASESTAT_BASIC; /*lint !e641*/
   row->validactivitylp = -1;
}

/** applies all cached row removals to the LP solver */
static
SCIP_RETCODE lpFlushDelRows(
   SCIP_LP*              lp,                 /**< current LP data */
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

      assert(!lp->diving);
      SCIPdebugMessage("flushing row deletions: shrink LP from %d to %d rows\n", lp->nlpirows, lp->lpifirstchgrow);
      SCIP_CALL( SCIPlpiDelRows(lp->lpi, lp->lpifirstchgrow, lp->nlpirows-1) );
      for( i = lp->lpifirstchgrow; i < lp->nlpirows; ++i )
      {
         markRowDeleted(lp->lpirows[i]);
         SCIP_CALL( SCIProwRelease(&lp->lpirows[i], blkmem, set, lp) );
      }
      lp->nlpirows = lp->lpifirstchgrow;
      lp->flushdeletedrows = TRUE;

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }
   assert(lp->nlpirows == lp->lpifirstchgrow);

   return SCIP_OKAY;
}

/** applies all cached row additions and removals to the LP solver */
static
SCIP_RETCODE lpFlushAddRows(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_Real* lhs;
   SCIP_Real* rhs;
   int* beg;
   int* ind;
   SCIP_Real* val;
   char** name;
   SCIP_ROW* row;
   SCIP_Real infinity;
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
   assert(!lp->diving);
   assert(lp->nrows > lp->nlpirows);
   SCIP_CALL( ensureLpirowsSize(lp, set, lp->nrows) );

   /* get the solver's infinity value */
   infinity = SCIPlpiInfinity(lp->lpi);

   /* count the (maximal) number of added coefficients, calculate the number of added rows */
   naddrows = lp->nrows - lp->nlpirows;
   naddcoefs = 0;
   for( r = lp->nlpirows; r < lp->nrows; ++r )
      naddcoefs += lp->rows[r]->len;
   assert(naddrows > 0);

   /* get temporary memory for changes */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lhs, naddrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rhs, naddrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &beg, naddrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ind, naddcoefs) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &val, naddcoefs) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &name, naddrows) );
   
   /* fill temporary memory with row data */
   nnonz = 0;
   for( pos = 0, r = lp->nlpirows; r < lp->nrows; ++pos, ++r )
   {
      row = lp->rows[r];
      assert(row != NULL);
      assert(row->lppos == r);
      assert(nnonz + row->nlpcols <= naddcoefs);

      SCIPdebugMessage("flushing added row <%s>: ", row->name);
      SCIPdebug( SCIProwPrint(row, NULL) );

      /* Because the row becomes a member of the LP solver, its dual variable now can take values
       * different from zero. That means, we have to include the row in the corresponding
       * column vectors.
       */
      SCIP_CALL( rowLink(row, blkmem, set, lp) );

      SCIProwCapture(row);
      lp->lpirows[r] = row;
      row->lpipos = r;
      row->dualsol = SCIP_INVALID;
      row->activity = SCIP_INVALID;
      row->dualfarkas = SCIP_INVALID;
      row->validactivitylp = -1;
      row->lhschanged = FALSE;
      row->rhschanged = FALSE;
      row->coefchanged = FALSE;
      if( SCIPsetIsInfinity(set, -row->lhs) )
         lhs[pos] = -infinity;
      else
         lhs[pos] = row->lhs - row->constant;
      if( SCIPsetIsInfinity(set, row->rhs) )
         rhs[pos] = infinity;
      else
         rhs[pos] = row->rhs - row->constant;
      beg[pos] = nnonz;
      name[pos] = row->name;

      row->flushedlhs = lhs[pos];
      row->flushedrhs = rhs[pos];

      SCIPdebugMessage("flushing added row (SCIP_LPI): %+g <=", lhs[pos]);
      for( i = 0; i < row->nlpcols; ++i )
      {
         assert(row->cols[i] != NULL);
         lpipos = row->cols[i]->lpipos;
         if( lpipos >= 0 )
         {
            assert(lpipos < lp->ncols);
            assert(nnonz < naddcoefs);
            SCIPdebugPrintf(" %+gx%d(<%s>)", row->vals[i], lpipos+1, SCIPvarGetName(row->cols[i]->var));
            ind[nnonz] = lpipos;
            val[nnonz] = row->vals[i];
            nnonz++;
         }
      }
      SCIPdebugPrintf(" <= %+g\n", rhs[pos]);
#ifndef NDEBUG
      for( i = row->nlpcols; i < row->len; ++i )
      {
         assert(row->cols[i] != NULL);
         assert(row->cols[i]->lpipos == -1); /* because the column deletions are already performed */
      }
#endif
   }

   /* call LP interface */
   SCIPdebugMessage("flushing row additions: enlarge LP from %d to %d rows\n", lp->nlpirows, lp->nrows);
   SCIP_CALL( SCIPlpiAddRows(lp->lpi, naddrows, lhs, rhs, name, nnonz, beg, ind, val) );
   lp->nlpirows = lp->nrows;
   lp->lpifirstchgrow = lp->nlpirows;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &name);
   SCIPsetFreeBufferArray(set, &val);
   SCIPsetFreeBufferArray(set, &ind);
   SCIPsetFreeBufferArray(set, &beg);
   SCIPsetFreeBufferArray(set, &rhs);
   SCIPsetFreeBufferArray(set, &lhs);

   lp->flushaddedrows = TRUE;

   /* mark the LP unsolved */
   lp->solved = FALSE;
   lp->primalfeasible = FALSE;
   lp->lpobjval = SCIP_INVALID;
   lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   
   return SCIP_OKAY;
}

/** applies all cached column bound and objective changes to the LP */
static
SCIP_RETCODE lpFlushChgCols(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_COL* col;
   int* objind;
   int* bdind;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real infinity;
   int nobjchg;
   int nbdchg;
   int i;

   assert(lp != NULL);

   if( lp->nchgcols == 0 )
      return SCIP_OKAY;

   /* get the solver's infinity value */
   infinity = SCIPlpiInfinity(lp->lpi);

   /* get temporary memory for changes */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &objind, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &obj, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bdind, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lb, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ub, lp->ncols) );

   /* collect all cached bound and objective changes */
   nobjchg = 0;
   nbdchg = 0;
   for( i = 0; i < lp->nchgcols; ++i )
   {
      col = lp->chgcols[i];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(col->var) == col);

      if( col->lpipos >= 0 )
      {
#ifndef NDEBUG
         SCIP_Real lpiobj;
         SCIP_Real lpilb;
         SCIP_Real lpiub;

         SCIP_CALL( SCIPlpiGetObj(lp->lpi, col->lpipos, col->lpipos, &lpiobj) );
         SCIP_CALL( SCIPlpiGetBounds(lp->lpi, col->lpipos, col->lpipos, &lpilb, &lpiub) );
         assert(SCIPsetIsFeasEQ(set, lpiobj, col->flushedobj));
         assert(SCIPsetIsFeasEQ(set, lpilb, col->flushedlb));
         assert(SCIPsetIsFeasEQ(set, lpiub, col->flushedub));
#endif

         if( col->objchanged )
         {
            SCIP_Real newobj;

            newobj = col->obj;
            if( col->flushedobj != newobj ) /*lint !e777*/
            {
               assert(nobjchg < lp->ncols);
               objind[nobjchg] = col->lpipos;
               obj[nobjchg] = newobj;
               nobjchg++;
               col->flushedobj = newobj;
            }
            col->objchanged = FALSE;
         }
         if( col->lbchanged || col->ubchanged )
         {
            SCIP_Real newlb;
            SCIP_Real newub;

            newlb = (SCIPsetIsInfinity(set, -col->lb) ? -infinity : col->lb);
            newub = (SCIPsetIsInfinity(set, col->ub) ? infinity : col->ub);
            if( col->flushedlb != newlb || col->flushedub != newub ) /*lint !e777*/
            {
               assert(nbdchg < lp->ncols);
               bdind[nbdchg] = col->lpipos;
               lb[nbdchg] = newlb;
               ub[nbdchg] = newub;
               nbdchg++;
               col->flushedlb = newlb;
               col->flushedub = newub;
            }
            col->lbchanged = FALSE;
            col->ubchanged = FALSE;
         }
      }
   }

   /* change objective values in LP */
   if( nobjchg > 0 )
   {
      SCIPdebugMessage("flushing objective changes: change %d objective values of %d changed columns\n", nobjchg, lp->nchgcols);
      SCIP_CALL( SCIPlpiChgObj(lp->lpi, nobjchg, objind, obj) );

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   /* change bounds in LP */
   if( nbdchg > 0 )
   {
      SCIPdebugMessage("flushing bound changes: change %d bounds of %d changed columns\n", nbdchg, lp->nchgcols);
      SCIP_CALL( SCIPlpiChgBounds(lp->lpi, nbdchg, bdind, lb, ub) );

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   lp->nchgcols = 0;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &ub);
   SCIPsetFreeBufferArray(set, &lb);
   SCIPsetFreeBufferArray(set, &bdind);
   SCIPsetFreeBufferArray(set, &obj);
   SCIPsetFreeBufferArray(set, &objind);

   return SCIP_OKAY;
}

/** applies all cached row side changes to the LP */
static
SCIP_RETCODE lpFlushChgRows(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_ROW* row;
   int* ind;
   SCIP_Real* lhs;
   SCIP_Real* rhs;
   SCIP_Real infinity;
   int i;
   int nchg;

   assert(lp != NULL);

   if( lp->nchgrows == 0 )
      return SCIP_OKAY;

   assert(!lp->diving);

   /* get the solver's infinity value */
   infinity = SCIPlpiInfinity(lp->lpi);

   /* get temporary memory for changes */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ind, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lhs, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rhs, lp->nrows) );

   /* collect all cached left and right hand side changes */
   nchg = 0;
   for( i = 0; i < lp->nchgrows; ++i )
   {
      row = lp->chgrows[i];
      assert(row != NULL);

      if( row->lpipos >= 0 )
      {
#ifndef NDEBUG
         SCIP_Real lpilhs;
         SCIP_Real lpirhs;

         SCIP_CALL( SCIPlpiGetSides(lp->lpi, row->lpipos, row->lpipos, &lpilhs, &lpirhs) );
         assert(SCIPsetIsSumEQ(set, lpilhs, row->flushedlhs));
         assert(SCIPsetIsSumEQ(set, lpirhs, row->flushedrhs));
#endif
         if( row->lhschanged || row->rhschanged )
         {
            SCIP_Real newlhs;
            SCIP_Real newrhs;

            newlhs = (SCIPsetIsInfinity(set, -row->lhs) ? -infinity : row->lhs - row->constant);
            newrhs = (SCIPsetIsInfinity(set, row->rhs) ? infinity : row->rhs - row->constant);
            if( row->flushedlhs != newlhs || row->flushedrhs != newrhs ) /*lint !e777*/
            {
               assert(nchg < lp->nrows);
               ind[nchg] = row->lpipos;
               lhs[nchg] = newlhs;
               rhs[nchg] = newrhs;
               nchg++;
               row->flushedlhs = newlhs;
               row->flushedrhs = newrhs;
            }
            row->lhschanged = FALSE;
            row->rhschanged = FALSE;
         }
      }
   }

   /* change left and right hand sides in LP */
   if( nchg > 0 )
   {
      SCIPdebugMessage("flushing side changes: change %d sides of %d rows\n", nchg, lp->nchgrows);
      SCIP_CALL( SCIPlpiChgSides(lp->lpi, nchg, ind, lhs, rhs) );

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   lp->nchgrows = 0;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &rhs);
   SCIPsetFreeBufferArray(set, &lhs);
   SCIPsetFreeBufferArray(set, &ind);

   return SCIP_OKAY;
}

/** applies all cached changes to the LP solver */
SCIP_RETCODE SCIPlpFlush(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(blkmem != NULL);
   
   SCIPdebugMessage("flushing LP changes: old (%d cols, %d rows), nchgcols=%d, nchgrows=%d, firstchgcol=%d, firstchgrow=%d, new (%d cols, %d rows), flushed=%d\n",
      lp->nlpicols, lp->nlpirows, lp->nchgcols, lp->nchgrows, lp->lpifirstchgcol, lp->lpifirstchgrow, 
      lp->ncols, lp->nrows, lp->flushed);

   if( lp->flushed )
   {
      assert(lp->nlpicols == lp->ncols);
      assert(lp->lpifirstchgcol == lp->nlpicols);
      assert(lp->nlpirows == lp->nrows);
      assert(lp->lpifirstchgrow == lp->nlpirows);
      assert(lp->nchgcols == 0);

      return SCIP_OKAY;
   }

   lp->flushdeletedcols = FALSE;
   lp->flushaddedcols = FALSE;
   lp->flushdeletedrows = FALSE;
   lp->flushaddedrows = FALSE;

   SCIP_CALL( lpFlushDelCols(lp) );
   SCIP_CALL( lpFlushDelRows(lp, blkmem, set) );
   SCIP_CALL( lpFlushChgCols(lp, set) );
   SCIP_CALL( lpFlushChgRows(lp, set) );
   SCIP_CALL( lpFlushAddCols(lp, blkmem, set) );
   SCIP_CALL( lpFlushAddRows(lp, blkmem, set) );

   lp->flushed = TRUE;

   checkLinks(lp);

   return SCIP_OKAY;
}

/** marks the LP to be flushed, even if the LP thinks it is not flushed */
SCIP_RETCODE SCIPlpMarkFlushed(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(lp != NULL);

#ifndef NDEBUG
   /* check, if there are really no column or row deletions or coefficient changes left */
   while( lp->lpifirstchgcol < lp->nlpicols
      && lp->lpifirstchgcol < lp->ncols
      && lp->cols[lp->lpifirstchgcol]->lpipos == lp->lpifirstchgcol
      && !lp->cols[lp->lpifirstchgcol]->coefchanged )
   {
      assert(lp->cols[lp->lpifirstchgcol] == lp->lpicols[lp->lpifirstchgcol]);
      lp->lpifirstchgcol++;
   }
   assert(lp->nlpicols == lp->lpifirstchgcol);

   while( lp->lpifirstchgrow < lp->nlpirows
      && lp->lpifirstchgrow < lp->nrows
      && lp->rows[lp->lpifirstchgrow]->lpipos == lp->lpifirstchgrow
      && !lp->rows[lp->lpifirstchgrow]->coefchanged )
   {
      assert(lp->rows[lp->lpifirstchgrow] == lp->lpirows[lp->lpifirstchgrow]);
      lp->lpifirstchgrow++;
   }
   assert(lp->nlpirows == lp->lpifirstchgrow);
#endif

   lp->lpifirstchgcol = lp->nlpicols;
   lp->lpifirstchgrow = lp->nlpirows;

   /* check, if there are really no column or row additions left */
   assert(lp->ncols == lp->nlpicols);
   assert(lp->nrows == lp->nlpirows);

   /* mark the changed columns to be unchanged, and check, if this is really correct */
   for( i = 0; i < lp->nchgcols; ++i )
   {
      SCIP_COL* col;

      col = lp->chgcols[i];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(col->var) == col);

      if( col->lpipos >= 0 )
      {
#ifndef NDEBUG
         SCIP_Real lpiobj;
         SCIP_Real lpilb;
         SCIP_Real lpiub;

         SCIP_CALL( SCIPlpiGetObj(lp->lpi, col->lpipos, col->lpipos, &lpiobj) );
         SCIP_CALL( SCIPlpiGetBounds(lp->lpi, col->lpipos, col->lpipos, &lpilb, &lpiub) );
         assert(SCIPsetIsSumEQ(set, lpiobj, col->flushedobj));
         assert(SCIPsetIsSumEQ(set, lpilb, col->flushedlb));
         assert(SCIPsetIsSumEQ(set, lpiub, col->flushedub));
         assert(col->flushedobj == col->obj); /*lint !e777*/
         assert(col->flushedlb == (SCIPsetIsInfinity(set, -col->lb) ? -SCIPlpiInfinity(lp->lpi) : col->lb)); /*lint !e777*/
         assert(col->flushedub == (SCIPsetIsInfinity(set, col->ub) ? SCIPlpiInfinity(lp->lpi) : col->ub)); /*lint !e777*/
#endif
         col->objchanged = FALSE;
         col->lbchanged = FALSE;
         col->ubchanged = FALSE;
      }
   }
   lp->nchgcols = 0;

   /* mark the changed rows to be unchanged, and check, if this is really correct */
   for( i = 0; i < lp->nchgrows; ++i )
   {
      SCIP_ROW* row;

      row = lp->chgrows[i];
      assert(row != NULL);

      if( row->lpipos >= 0 )
      {
#ifndef NDEBUG
         SCIP_Real lpilhs;
         SCIP_Real lpirhs;

         SCIP_CALL( SCIPlpiGetSides(lp->lpi, row->lpipos, row->lpipos, &lpilhs, &lpirhs) );
         assert(SCIPsetIsSumEQ(set, lpilhs, row->flushedlhs));
         assert(SCIPsetIsSumEQ(set, lpirhs, row->flushedrhs));
         assert(row->flushedlhs ==
            (SCIPsetIsInfinity(set, -row->lhs) ? -SCIPlpiInfinity(lp->lpi) : row->lhs - row->constant)); /*lint !e777*/
         assert(row->flushedrhs ==
            (SCIPsetIsInfinity(set, row->rhs) ? SCIPlpiInfinity(lp->lpi) : row->rhs - row->constant)); /*lint !e777*/
#endif
         row->lhschanged = FALSE;
         row->rhschanged = FALSE;
      }
   }
   lp->nchgrows = 0;

   /* mark the LP to be flushed */
   lp->flushed = TRUE;

   checkLinks(lp);

   return SCIP_OKAY;
}




/*
 * LP methods
 */

/** updates link data after addition of column */
static
void colUpdateAddLP(
   SCIP_COL*             col                 /**< LP column */
   )
{
   SCIP_ROW* row;
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
         rowSwapCoefs(row, pos, row->nlpcols-1);
      }
   }
}

/** updates link data after addition of row */
static
void rowUpdateAddLP(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_COL* col;
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
         colSwapCoefs(col, pos, col->nlprows-1);
      }
   }
}

/** updates link data after removal of column */
static
void colUpdateDelLP(
   SCIP_COL*             col                 /**< LP column */
   )
{
   SCIP_ROW* row;
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
         rowSwapCoefs(row, pos, row->nlpcols);
      }
   }
}

/** updates link data after removal of row */
static
void rowUpdateDelLP(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_COL* col;
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
         colSwapCoefs(col, pos, col->nlprows);
      }
   }
}

/** creates empty LP data object */
SCIP_RETCODE SCIPlpCreate(
   SCIP_LP**             lp,                 /**< pointer to LP data object */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name                /**< problem name */
   )
{
   SCIP_Bool success;

   assert(lp != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(name != NULL);

   SCIP_ALLOC( BMSallocMemory(lp) );

   /* open LP Solver interface */
   SCIP_CALL( SCIPlpiCreate(&(*lp)->lpi, name, SCIP_OBJSEN_MINIMIZE) );

   (*lp)->lpicols = NULL;
   (*lp)->lpirows = NULL;
   (*lp)->chgcols = NULL;
   (*lp)->chgrows = NULL;
   (*lp)->cols = NULL;
   (*lp)->rows = NULL;
   (*lp)->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   (*lp)->lpobjval = 0.0;
   (*lp)->pseudoobjval = 0.0;
   (*lp)->pseudoobjvalinf = 0;
   (*lp)->looseobjval = 0.0;
   (*lp)->looseobjvalinf = 0;
   (*lp)->nloosevars = 0;
   (*lp)->cutoffbound = SCIPsetInfinity(set);
   (*lp)->objsqrnorm = 0.0;
   (*lp)->objsumnorm = 0.0;
   (*lp)->lpicolssize = 0;
   (*lp)->nlpicols = 0;
   (*lp)->lpirowssize = 0;
   (*lp)->nlpirows = 0;
   (*lp)->lpifirstchgcol = 0;
   (*lp)->lpifirstchgrow = 0;
   (*lp)->colssize = 0;
   (*lp)->ncols = 0;
   (*lp)->rowssize = 0;
   (*lp)->nrows = 0;
   (*lp)->chgcolssize = 0;
   (*lp)->nchgcols = 0;
   (*lp)->chgrowssize = 0;
   (*lp)->nchgrows = 0;
   (*lp)->firstnewcol = 0;
   (*lp)->firstnewrow = 0;
   (*lp)->nremoveablecols = 0;
   (*lp)->nremoveablerows = 0;
   (*lp)->validsollp = stat->lpcount; /* the initial (empty) SCIP_LP is solved with primal and dual solution of zero */
   (*lp)->validfarkaslp = -1;
   (*lp)->flushdeletedcols = FALSE;
   (*lp)->flushaddedcols = FALSE;
   (*lp)->flushdeletedrows = FALSE;
   (*lp)->flushaddedrows = FALSE;
   (*lp)->flushed = TRUE;
   (*lp)->solved = TRUE;
   (*lp)->primalfeasible = TRUE;
   (*lp)->dualfeasible = TRUE;
   (*lp)->solisbasic = FALSE;
   (*lp)->probing = FALSE;
   (*lp)->diving = FALSE;
   (*lp)->divingobjchg = FALSE;
   (*lp)->divelpistate = NULL;
   (*lp)->lpiuobjlim = SCIPsetInfinity(set);
   (*lp)->lpifeastol = SCIPsetFeastol(set);
   (*lp)->lpidualfeastol = SCIPsetDualfeastol(set);
   (*lp)->lpibarrierconvtol = SCIPsetBarrierconvtol(set);
   (*lp)->lpifromscratch = FALSE;
   (*lp)->lpifastmip = TRUE;
   (*lp)->lpiscaling = TRUE;
   (*lp)->lpipresolving = TRUE;
   (*lp)->lpilpinfo = FALSE;
   (*lp)->lpiitlim = INT_MAX;
   (*lp)->lpipricing = SCIP_PRICING_AUTO;
   (*lp)->lastlpalgo = SCIP_LPALGO_DUALSIMPLEX;

   /* set default parameters in LP solver */
   SCIP_CALL( lpSetRealpar(*lp, SCIP_LPPAR_UOBJLIM, (*lp)->lpiuobjlim, &success) );
   if( !success )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: upper objective limit cannot be set -- can lead to unnecessary simplex iterations\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetRealpar(*lp, SCIP_LPPAR_FEASTOL, (*lp)->lpifeastol, &success) );
   (*lp)->lpihasfeastol = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: primal feasibility tolerance cannot be set -- tolerance of SCIP and LP solver may differ\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetRealpar(*lp, SCIP_LPPAR_DUALFEASTOL, (*lp)->lpidualfeastol, &success) );
   (*lp)->lpihasdualfeastol = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: dual feasibility tolerance cannot be set -- tolerance of SCIP and LP solver may differ\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetRealpar(*lp, SCIP_LPPAR_BARRIERCONVTOL, (*lp)->lpibarrierconvtol, &success) );
   (*lp)->lpihasbarrierconvtol = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: barrier convergence tolerance cannot be set -- tolerance of SCIP and LP solver may differ\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetBoolpar(*lp, SCIP_LPPAR_FROMSCRATCH, (*lp)->lpifromscratch, &success) );
   SCIP_CALL( lpSetBoolpar(*lp, SCIP_LPPAR_FASTMIP, (*lp)->lpifastmip, &success) );
   (*lp)->lpihasfastmip = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: fastmip setting not available -- SCIP parameter has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetBoolpar(*lp, SCIP_LPPAR_SCALING, (*lp)->lpiscaling, &success) );
   (*lp)->lpihasscaling = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: scaling not available -- SCIP parameter has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetBoolpar(*lp, SCIP_LPPAR_PRESOLVING, (*lp)->lpipresolving, &success) );
   (*lp)->lpihaspresolving = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: presolving not available -- SCIP parameter has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetIntpar(*lp, SCIP_LPPAR_LPITLIM, (*lp)->lpiitlim, &success) );
   if( !success )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: iteration limit cannot be set -- can lead to unnecessary simplex iterations\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetIntpar(*lp, SCIP_LPPAR_PRICING, (int)(*lp)->lpipricing, &success) );
   if( !success )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: pricing strategy cannot be set -- SCIP parameter has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetBoolpar(*lp, SCIP_LPPAR_LPINFO, (*lp)->lpilpinfo, &success) );
   if( !success )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: lpinfo setting not available -- SCIP parameter has no effect\n",
         SCIPlpiGetSolverName());
   }

   return SCIP_OKAY;
}

/** frees LP data object */
SCIP_RETCODE SCIPlpFree(
   SCIP_LP**             lp,                 /**< pointer to LP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(lp != NULL);
   assert(*lp != NULL);
   
   SCIP_CALL( SCIPlpClear(*lp, blkmem, set) );

   /* release LPI rows */
   for( i = 0; i < (*lp)->nlpirows; ++i )
   {
      SCIP_CALL( SCIProwRelease(&(*lp)->lpirows[i], blkmem, set, *lp) );
   }

   if( (*lp)->lpi != NULL )
   {
      SCIP_CALL( SCIPlpiFree(&(*lp)->lpi) );
   }

   BMSfreeMemoryArrayNull(&(*lp)->lpicols);
   BMSfreeMemoryArrayNull(&(*lp)->lpirows);
   BMSfreeMemoryArrayNull(&(*lp)->chgcols);
   BMSfreeMemoryArrayNull(&(*lp)->cols);
   BMSfreeMemoryArrayNull(&(*lp)->rows);
   BMSfreeMemory(lp);

   return SCIP_OKAY;
}

/** resets the LP to the empty LP by removing all columns and rows from LP, releasing all rows, and flushing the
 *  changes to the LP solver
 */
SCIP_RETCODE SCIPlpReset(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(stat != NULL);

   SCIP_CALL( SCIPlpClear(lp, blkmem, set) );
   SCIP_CALL( SCIPlpFlush(lp, blkmem, set) );

   /* mark the empty LP to be solved */
   lp->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   lp->lpobjval = 0.0;
   lp->validsollp = stat->lpcount; /* the initial (empty) SCIP_LP is solved with primal and dual solution of zero */
   lp->validfarkaslp = -1;
   lp->solved = TRUE;
   lp->primalfeasible = TRUE;
   lp->dualfeasible = TRUE;
   lp->solisbasic = FALSE;
   lp->lastlpalgo = SCIP_LPALGO_DUALSIMPLEX;

   return SCIP_OKAY;
}

/** adds a column to the LP */
SCIP_RETCODE SCIPlpAddCol(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COL*             col,                /**< LP column */
   int                   depth               /**< depth in the tree where the column addition is performed */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);
   assert(col != NULL);
   assert(col->len == 0 || col->rows != NULL);
   assert(col->lppos == -1);
   assert(col->var != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(SCIPvarIsIntegral(col->var) == col->integral);
   
   SCIPdebugMessage("adding column <%s> to LP (%d rows, %d cols)\n", SCIPvarGetName(col->var), lp->nrows, lp->ncols);
#ifdef SCIP_DEBUG
   {
      int i;
      SCIPdebugPrintf("  (obj: %g) [%g,%g]", col->obj, col->lb, col->ub);
      for( i = 0; i < col->len; ++i )
         SCIPdebugPrintf(" %+g<%s>", col->vals[i], col->rows[i]->name);
      SCIPdebugPrintf("\n");
   }
#endif

   SCIP_CALL( ensureColsSize(lp, set, lp->ncols+1) );
   lp->cols[lp->ncols] = col;
   col->lppos = lp->ncols;
   col->lpdepth = depth;
   col->age = 0;
   lp->ncols++;
   if( col->removeable )
      lp->nremoveablecols++;

   /* mark the current LP unflushed */
   lp->flushed = FALSE;

   /* update column arrays of all linked rows */
   colUpdateAddLP(col);

   checkLinks(lp);

   return SCIP_OKAY;
}

/** adds a row to the LP and captures it */
SCIP_RETCODE SCIPlpAddRow(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row,                /**< LP row */
   int                   depth               /**< depth in the tree where the row addition is performed */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);
   assert(row != NULL);
   assert(row->len == 0 || row->cols != NULL);
   assert(row->lppos == -1);

   SCIProwCapture(row);
   SCIProwLock(row);

   SCIPdebugMessage("adding row <%s> to LP (%d rows, %d cols)\n", row->name, lp->nrows, lp->ncols);
#ifdef SCIP_DEBUG
   {
      int i;
      SCIPdebugPrintf("  %g <=", row->lhs);
      for( i = 0; i < row->len; ++i )
         SCIPdebugPrintf(" %+g<%s>", row->vals[i], SCIPvarGetName(row->cols[i]->var));
      if( !SCIPsetIsZero(set, row->constant) )
         SCIPdebugPrintf(" %+g", row->constant);
      SCIPdebugPrintf(" <= %g\n", row->rhs);
   }
#endif

   SCIP_CALL( ensureRowsSize(lp, set, lp->nrows+1) );
   lp->rows[lp->nrows] = row;
   row->lppos = lp->nrows;
   row->lpdepth = depth;
   row->age = 0;
   lp->nrows++;
   if( row->removeable )
      lp->nremoveablerows++;

   /* mark the current LP unflushed */
   lp->flushed = FALSE;

   /* update row arrays of all linked columns */
   rowUpdateAddLP(row);

   checkLinks(lp);

   return SCIP_OKAY;
}

/** removes all columns after the given number of cols from the LP */
SCIP_RETCODE SCIPlpShrinkCols(
   SCIP_LP*              lp,                 /**< LP data */
   int                   newncols            /**< new number of columns in the LP */
   )
{
   SCIP_COL* col;
   int c;

   assert(lp != NULL);
   SCIPdebugMessage("shrinking LP from %d to %d columns\n", lp->ncols, newncols);
   assert(0 <= newncols);
   assert(newncols <= lp->ncols);

   if( newncols < lp->ncols )
   {
      assert(!lp->diving);

      for( c = lp->ncols-1; c >= newncols; --c )
      {
         col = lp->cols[c];
         assert(col != NULL);
         assert(col->len == 0 || col->rows != NULL);
         assert(col->var != NULL);
         assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetCol(col->var) == lp->cols[c]);
         assert(col->lppos == c);
         
         /* mark column to be removed from the LP */
         col->lppos = -1;
         col->lpdepth = -1;
         lp->ncols--;

         /* count removeable columns */
         if( col->removeable )
            lp->nremoveablecols--;

         /* update column arrays of all linked rows */
         colUpdateDelLP(col);
      }
      assert(lp->ncols == newncols);
      lp->lpifirstchgcol = MIN(lp->lpifirstchgcol, newncols);

      /* mark the current LP unflushed */
      lp->flushed = FALSE;

      checkLinks(lp);
   }
   assert(lp->nremoveablecols <= lp->ncols);

   return SCIP_OKAY;
}

/** removes and releases all rows after the given number of rows from the LP */
SCIP_RETCODE SCIPlpShrinkRows(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   newnrows            /**< new number of rows in the LP */
   )
{
   SCIP_ROW* row;
   int r;

   assert(lp != NULL);
   assert(0 <= newnrows && newnrows <= lp->nrows);

   SCIPdebugMessage("shrinking LP from %d to %d rows\n", lp->nrows, newnrows);
   if( newnrows < lp->nrows )
   {
      assert(!lp->diving);

      for( r = lp->nrows-1; r >= newnrows; --r )
      {
         row = lp->rows[r];
         assert(row != NULL);
         assert(row->len == 0 || row->cols != NULL);
         assert(row->lppos == r);

         /* mark row to be removed from the LP */
         row->lppos = -1;
         row->lpdepth = -1;
         lp->nrows--;

         /* count removeable rows */
         if( row->removeable )
            lp->nremoveablerows--;

         /* update row arrays of all linked columns */
         rowUpdateDelLP(row);

         SCIProwUnlock(lp->rows[r]);
         SCIP_CALL( SCIProwRelease(&lp->rows[r], blkmem, set, lp) );
      }
      assert(lp->nrows == newnrows);
      lp->lpifirstchgrow = MIN(lp->lpifirstchgrow, newnrows);

      /* mark the current LP unflushed */
      lp->flushed = FALSE;

      checkLinks(lp);
   }
   assert(lp->nremoveablerows <= lp->nrows);

   return SCIP_OKAY;
}

/** removes all columns and rows from LP, releases all rows */
SCIP_RETCODE SCIPlpClear(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);

   SCIPdebugMessage("clearing LP\n");
   SCIP_CALL( SCIPlpShrinkCols(lp, 0) );
   SCIP_CALL( SCIPlpShrinkRows(lp, blkmem, set, 0) );

   return SCIP_OKAY;
}

/** remembers number of columns and rows to track the newly added ones */
void SCIPlpMarkSize(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);

   lp->firstnewcol = lp->ncols;
   lp->firstnewrow = lp->nrows;
}

/** gets all indices of basic columns and rows: index i >= 0 corresponds to column i, index i < 0 to row -i-1 */
SCIP_RETCODE SCIPlpGetBasisInd(
   SCIP_LP*              lp,                 /**< LP data */
   int*                  basisind            /**< pointer to store the basis indices */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->solisbasic);
   assert(basisind != NULL);

   SCIP_CALL( SCIPlpiGetBasisInd(lp->lpi, basisind) );

   return SCIP_OKAY;
}

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpGetBase(
   SCIP_LP*              lp,                 /**< LP data */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->solisbasic);

   SCIP_CALL( SCIPlpiGetBase(lp->lpi, cstat, rstat) );

   return SCIP_OKAY;
}

/** gets a row from the inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpGetBInvRow(
   SCIP_LP*              lp,                 /**< LP data */
   int                   r,                  /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->solisbasic);
   assert(0 <= r && r < lp->nrows);  /* the basis matrix is nrows x nrows */
   assert(coef != NULL);

   SCIP_CALL( SCIPlpiGetBInvRow(lp->lpi, r, coef) );

   return SCIP_OKAY;
}

/** gets a row from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A) */
SCIP_RETCODE SCIPlpGetBInvARow(
   SCIP_LP*              lp,                 /**< LP data */
   int                   r,                  /**< row number */
   SCIP_Real*            binvrow,            /**< row in B^-1 from prior call to SCIPlpGetBInvRow(), or NULL */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->solisbasic);
   assert(0 <= r && r < lp->nrows);  /* the basis matrix is nrows x nrows */
   assert(coef != NULL);

   SCIP_CALL( SCIPlpiGetBInvARow(lp->lpi, r, binvrow, coef) );

   return SCIP_OKAY;
}

/** calculates a weighted sum of all LP rows; for negative weights, the left and right hand side of the corresponding
 *  LP row are swapped in the summation
 */
SCIP_RETCODE SCIPlpSumRows(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_REALARRAY*       sumcoef,            /**< array to store sum coefficients indexed by variables' probindex */
   SCIP_Real*            sumlhs,             /**< pointer to store the left hand side of the row summation */
   SCIP_Real*            sumrhs              /**< pointer to store the right hand side of the row summation */
   )
{
   SCIP_ROW* row;
   int r;
   int i;
   int idx;
   SCIP_Bool lhsinfinite;
   SCIP_Bool rhsinfinite;

   assert(lp != NULL);
   assert(prob != NULL);
   assert(weights != NULL);
   assert(sumcoef != NULL);
   assert(sumlhs != NULL);
   assert(sumrhs != NULL);

   /**@todo test, if a column based summation is faster */

   SCIP_CALL( SCIPrealarrayClear(sumcoef) );
   SCIP_CALL( SCIPrealarrayExtend(sumcoef, set, 0, prob->nvars-1) );
   *sumlhs = 0.0;
   *sumrhs = 0.0;
   lhsinfinite = FALSE;
   rhsinfinite = FALSE;
   for( r = 0; r < lp->nrows; ++r )
   {
      if( !SCIPsetIsZero(set, weights[r]) )
      {
         row = lp->rows[r];
         assert(row != NULL);
         assert(row->len == 0 || row->cols != NULL);
         assert(row->len == 0 || row->cols_index != NULL);
         assert(row->len == 0 || row->vals != NULL);

         /* add the row coefficients to the sum */
         for( i = 0; i < row->len; ++i )
         {
            assert(row->cols[i] != NULL);
            assert(row->cols[i]->var != NULL);
            assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
            assert(SCIPvarGetCol(row->cols[i]->var) == row->cols[i]);
            assert(SCIPvarGetProbindex(row->cols[i]->var) == row->cols[i]->var_probindex);
            idx = row->cols[i]->var_probindex;
            assert(0 <= idx && idx < prob->nvars);
            SCIP_CALL( SCIPrealarrayIncVal(sumcoef, set, idx, weights[r] * row->vals[i]) );
         }
         
         /* add the row sides to the sum, depending on the sign of the weight */
         if( weights[r] > 0.0 )
         {
            lhsinfinite = lhsinfinite || SCIPsetIsInfinity(set, -row->lhs);
            if( !lhsinfinite )
               (*sumlhs) += weights[r] * (row->lhs - row->constant);
            rhsinfinite = rhsinfinite || SCIPsetIsInfinity(set, row->rhs);
            if( !rhsinfinite )
               (*sumrhs) += weights[r] * (row->rhs - row->constant);
         }
         else
         {
            lhsinfinite = lhsinfinite || SCIPsetIsInfinity(set, row->rhs);
            if( !lhsinfinite )
               (*sumlhs) += weights[r] * (row->rhs - row->constant);
            rhsinfinite = rhsinfinite || SCIPsetIsInfinity(set, -row->lhs);
            if( !rhsinfinite )
               (*sumrhs) += weights[r] * (row->lhs - row->constant);
         }
      }
      else
         weights[r] = 0.0;
   }

   *sumlhs = -SCIPsetInfinity(set);
   *sumrhs = SCIPsetInfinity(set);

   return SCIP_OKAY;
}

/** builds a weighted sum of rows, and decides whether to use the left or right hand side of the rows in summation */
static
void sumMIRRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Bool             allowlocal,         /**< should local rows be included, resulting in a locally valid summation? */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size prob->nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   SCIP_Bool*            emptyrow,           /**< pointer to store whether the returned row is empty */
   SCIP_Bool*            localrowsused       /**< pointer to store whether local rows were used in summation */
   )
{
   SCIP_ROW* row;
   SCIP_Real weight;
   SCIP_Real absweight;
   SCIP_Real maxweight;
   int idx;
   int r;
   int i;

   assert(prob != NULL);
   assert(lp != NULL);
   assert(weights != NULL);
   assert(SCIPsetIsPositive(set, scale));
   assert(maxweightrange >= 1.0);
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(slacksign != NULL);
   assert(emptyrow != NULL);
   assert(localrowsused != NULL);

   /* search the maximal absolute weight */
   maxweight = 0.0;
   for( r = 0; r < lp->nrows; ++r )
   {
      weight = scale * weights[r];
      absweight = REALABS(weight);
      maxweight = MAX(maxweight, absweight);
   }

   /* calculate the row summation */
   BMSclearMemoryArray(mircoef, prob->nvars);
   *mirrhs = 0.0;
   *emptyrow = TRUE;
   *localrowsused = FALSE;
   for( r = 0; r < lp->nrows; ++r )
   {
      row = lp->rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* modifiable rows cannot be part of a MIR row summation;
       * local rows are only included, if the allowlocal flag is set;
       * close to zero weights or weights outside the maximal range are ignored
       */
      weight = scale * weights[r];
      absweight = REALABS(weight);
      if( !row->modifiable && (allowlocal || !row->local)
         && absweight * maxweightrange >= maxweight && !SCIPsetIsSumZero(set, weight) )
      {
         *emptyrow = FALSE;
         *localrowsused = *localrowsused || row->local;

         /* Decide, if we want to use the left or the right hand side of the row in the summation.
          * If possible, use the side that leads to a positive slack value in the summation.
          */
         if( SCIPsetIsInfinity(set, row->rhs) || (!SCIPsetIsInfinity(set, -row->lhs) && weight < 0.0) )
         {
            slacksign[r] = -1;
            (*mirrhs) += weight * (row->lhs - row->constant);
         }
         else
         {
            slacksign[r] = +1;
            (*mirrhs) += weight * (row->rhs - row->constant);
         }

         /* add the row coefficients to the sum */
         for( i = 0; i < row->len; ++i )
         {
            assert(row->cols[i] != NULL);
            assert(row->cols[i]->var != NULL);
            assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
            assert(SCIPvarGetCol(row->cols[i]->var) == row->cols[i]);
            assert(SCIPvarGetProbindex(row->cols[i]->var) == row->cols[i]->var_probindex);
            idx = row->cols[i]->var_probindex;
            assert(0 <= idx && idx < prob->nvars);
            mircoef[idx] += weight * row->vals[i];
         }

         SCIPdebugMessage("MIR: %d: row <%s>, lhs = %g, rhs = %g, scale = %g, weight = %g, slacksign = %d -> rhs = %g\n",
            r, SCIProwGetName(row), row->lhs - row->constant, row->rhs - row->constant, 
            scale, weights[r], slacksign[r], *mirrhs);
         SCIPdebug(SCIProwPrint(row, NULL));
      }
      else
      {
         slacksign[r] = 0;
         weights[r] = 0.0;
      }
   }
}

/** removes all nearly-zero coefficients from MIR row and relaxes the right hand side correspondingly in order to
 *  prevent numerical rounding errors
 */
static
void cleanupMIRRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   SCIP_Bool             cutislocal          /**< is the cut only valid locally? */
   )
{
   SCIP_Real bd;
   SCIP_Bool rhsinf;
   int v;

   assert(prob != NULL);
   assert(mircoef != NULL);
   assert(mirrhs != NULL);

   rhsinf = SCIPsetIsInfinity(set, *mirrhs);
   for( v = 0; v < prob->nvars && !rhsinf; ++v )
   {
      if( mircoef[v] != 0.0 && SCIPsetIsSumZero(set, mircoef[v]) )
      {
         SCIPdebugMessage("coefficient of <%s> in transformed MIR row is too small: %.12f\n",
            SCIPvarGetName(prob->vars[v]), mircoef[v]);

         if( mircoef[v] > 0.0 )
         {
            bd = cutislocal ? SCIPvarGetLbLocal(prob->vars[v]) : SCIPvarGetLbGlobal(prob->vars[v]);
            rhsinf = rhsinf || SCIPsetIsInfinity(set, -bd);
         }
         else
         {
            bd = cutislocal ? SCIPvarGetUbLocal(prob->vars[v]) : SCIPvarGetUbGlobal(prob->vars[v]);
            rhsinf = rhsinf || SCIPsetIsInfinity(set, bd);
         }
         *mirrhs -= bd * mircoef[v];
         mircoef[v] = 0.0;
      }
   }
   if( rhsinf )
      *mirrhs = SCIPsetInfinity(set);
}

/** returns LP solution value and index of variable lower bound that is closest to variable's current LP solution value */
static
void getClosestVlb(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_Real*            closestvlb,         /**< pointer to store the value of the closest variable lower bound */
   int*                  closestvlbidx       /**< pointer to store the index of the closest variable lower bound */
   )
{
   SCIP_VAR** vlbvars;
   SCIP_Real* vlbcoefs;
   SCIP_Real* vlbconsts;
   int nvlbs;
   int i;

   assert(closestvlb != NULL);
   assert(closestvlbidx != NULL);

   *closestvlbidx = -1;
   *closestvlb = SCIP_REAL_MIN;

   nvlbs = SCIPvarGetNVlbs(var);
   if( nvlbs > 0 )
   {
      vlbvars = SCIPvarGetVlbVars(var);
      vlbcoefs = SCIPvarGetVlbCoefs(var);
      vlbconsts = SCIPvarGetVlbConstants(var);
      
      for( i = 0; i < nvlbs; i++ )
      {
         SCIP_Real vlbsol;
         
         vlbsol = vlbcoefs[i] * SCIPvarGetLPSol(vlbvars[i]) + vlbconsts[i];
         if( vlbsol > *closestvlb )
         {
            *closestvlb = vlbsol;
            *closestvlbidx = i;
         }
      }
   }
}

/** returns LP solution value and index of variable upper bound that is closest to variable's current LP solution value */
static
void getClosestVub(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_Real*            closestvub,         /**< pointer to store the value of the closest variable upper bound */
   int*                  closestvubidx       /**< pointer to store the index of the closest variable upper bound */
   )
{
   SCIP_VAR** vubvars;
   SCIP_Real* vubcoefs;
   SCIP_Real* vubconsts;
   int nvubs;
   int i;

   assert(closestvub != NULL);
   assert(closestvubidx != NULL);

   *closestvubidx = -1;
   *closestvub = SCIP_REAL_MAX;

   nvubs = SCIPvarGetNVubs(var);
   if( nvubs > 0 )
   {
      vubvars = SCIPvarGetVubVars(var);
      vubcoefs = SCIPvarGetVubCoefs(var);
      vubconsts = SCIPvarGetVubConstants(var);
      
      for( i = 0; i < nvubs; i++ )
      {
         SCIP_Real vubsol;
         
         vubsol = vubcoefs[i] * SCIPvarGetLPSol(vubvars[i]) + vubconsts[i];
         if( vubsol < *closestvub )
         {
            *closestvub = vubsol;
            *closestvubidx = i;
         }
      }
   }
}

/** Transform equation  a*x == b, lb <= x <= ub  into standard form
 *    a'*x' == b, 0 <= x' <= ub'.
 *  
 *  Transform variables (lb or ub):
 *    x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
 *    x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
 *  and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
 *
 *  Transform variables (vlb or vub):
 *    x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
 *    x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
 *  move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
 *    a_{zl_j} := a_{zl_j} + a_j * bl_j, or
 *    a_{zu_j} := a_{zu_j} + a_j * bu_j
 */
static
void transformMIRRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Bool*            freevariable,       /**< stores whether a free variable was found in MIR row -> invalid summation */
   SCIP_Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{
   SCIP_VAR* var;
   SCIP_Real varsol;
   SCIP_Real bestlb;
   SCIP_Real bestub;
   int bestlbtype;
   int bestubtype;
   int v;

   assert(prob != NULL);
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(varsign != NULL);
   assert(boundtype != NULL);
   assert(freevariable != NULL);
   assert(localbdsused != NULL);

   *freevariable = FALSE;
   *localbdsused = FALSE;
   
   /* substitute continuous variables with best standard or variable bound (lb, ub, vlb or vub),
    * substitute integral variables with best standard bound (lb, ub);
    * start with continuous variables, because using variable bounds can affect the untransformed integral
    * variables, and these changes have to be incorporated in the transformation of the integral variables
    */
   for( v = prob->nvars-1; v >= 0; --v )
   {
      var = prob->vars[v];
      assert(v == SCIPvarGetProbindex(var));

      /* ignore variables that don't exist in the MIR row */
      if( SCIPsetIsZero(set, mircoef[v]) )
      {
         varsign[v] = +1;
         boundtype[v] = -2;
         continue;
      }

      /* find closest lower bound in standard lower bound (and variable lower bounds for continuous variables) */
      bestlb = SCIPvarGetLbGlobal(var);
      bestlbtype = -1;
      if( allowlocal )
      {
         SCIP_Real loclb;

         loclb = SCIPvarGetLbLocal(var);
         if( SCIPsetIsGT(set, loclb, bestlb) )
         {
            bestlb = loclb;
            bestlbtype = -2;
         }
      }
      if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_Real bestvlb;
         int bestvlbidx;

         getClosestVlb(var, &bestvlb, &bestvlbidx);
         if( bestvlbidx >= 0
            && (bestvlb > bestlb || (bestlbtype < 0 && SCIPsetIsGE(set, bestvlb, bestlb))) )
         {
            bestlb = bestvlb;
            bestlbtype = bestvlbidx;
         }
      }

      /* find closest upper bound in standard upper bound (and variable upper bounds for continuous variables) */
      bestub = SCIPvarGetUbGlobal(var);
      bestubtype = -1;
      if( allowlocal )
      {
         SCIP_Real locub;

         locub = SCIPvarGetUbLocal(var);
         if( SCIPsetIsLT(set, locub, bestub) )
         {
            bestub = locub;
            bestubtype = -2;
         }
      }
      if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_Real bestvub;
         int bestvubidx;

         getClosestVub(var, &bestvub, &bestvubidx);
         if( bestvubidx >= 0
            && (bestvub < bestub || (bestubtype < 0 && SCIPsetIsLE(set, bestvub, bestub))) )
         {
            bestub = bestvub;
            bestubtype = bestvubidx;
         }
      }

      /* check, if variable is free variable */
      if( SCIPsetIsInfinity(set, -bestlb) && SCIPsetIsInfinity(set, bestub) )
      {
         /* we found a free variable in the row with non-zero coefficient
          *  -> MIR row can't be transformed in standard form
          */
         *freevariable = TRUE;
         return;
      }

      /* select transformation bound */
      varsol = SCIPvarGetLPSol(var);
      if( varsol <= (1.0 - boundswitch) * bestlb + boundswitch * bestub )
      {
         /* use lower bound as transformation bound: x'_j := x_j - lb_j */
         boundtype[v] = bestlbtype;
         varsign[v] = +1;
    
         /* standard (bestlbtype < 0) or variable (bestlbtype >= 0) lower bound? */
         if( bestlbtype < 0 )
         {
            (*mirrhs) -= mircoef[v] * bestlb;
            *localbdsused = *localbdsused || (bestlbtype == -2);
         }
         else
         {
            SCIP_VAR** vlbvars = SCIPvarGetVlbVars(var);
            SCIP_Real* vlbcoefs = SCIPvarGetVlbCoefs(var);
            SCIP_Real* vlbconsts = SCIPvarGetVlbConstants(var);
            int zidx;

            assert(0 <= bestlbtype && bestlbtype < SCIPvarGetNVlbs(var));
            zidx = SCIPvarGetProbindex(vlbvars[bestlbtype]);
            assert(0 <= zidx && zidx < v);
               
            (*mirrhs) -= mircoef[v] * vlbconsts[bestlbtype];
            mircoef[zidx] += mircoef[v] * vlbcoefs[bestlbtype];
         }
      }
      else
      {
         /* use upper bound as transformation bound: x'_j := ub_j - x_j */
         boundtype[v] = bestubtype;
         varsign[v] = -1;
         
         /* standard (bestubtype < 0) or variable (bestubtype >= 0) upper bound? */
         if( bestubtype < 0 )
         {
            (*mirrhs) -= mircoef[v] * bestub;
            *localbdsused = *localbdsused || (bestubtype == -2);
         }
         else
         {
            SCIP_VAR** vubvars = SCIPvarGetVubVars(var);
            SCIP_Real* vubcoefs = SCIPvarGetVubCoefs(var);
            SCIP_Real* vubconsts = SCIPvarGetVubConstants(var);
            int zidx;

            assert(0 <= bestubtype && bestubtype < SCIPvarGetNVubs(var));
            zidx = SCIPvarGetProbindex(vubvars[bestubtype]);
            assert(zidx >= 0);
               
            (*mirrhs) -= mircoef[v] * vubconsts[bestubtype];
            mircoef[zidx] += mircoef[v] * vubcoefs[bestubtype];
         }
      }

      SCIPdebugMessage("MIR var <%s>: varsign=%d, boundtype=%d, mircoef=%g, lb=%g, ub=%g, vlb=%g<%s>%+g, vub=%g<%s>%+g -> rhs=%g\n", 
         SCIPvarGetName(var), varsign[v], boundtype[v], mircoef[v], bestlb, bestub,
         bestlbtype >= 0 ? SCIPvarGetVlbCoefs(var)[bestlbtype] : 0.0,
         bestlbtype >= 0 ? SCIPvarGetName(SCIPvarGetVlbVars(var)[bestlbtype]) : "-",
         bestlbtype >= 0 ? SCIPvarGetVlbConstants(var)[bestlbtype] : bestlb,
         bestubtype >= 0 ? SCIPvarGetVubCoefs(var)[bestubtype] : 0.0,
         bestubtype >= 0 ? SCIPvarGetName(SCIPvarGetVubVars(var)[bestubtype]) : "-",
         bestubtype >= 0 ? SCIPvarGetVubConstants(var)[bestubtype] : bestub,
         *mirrhs);
   }
}

/** Calculate fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j) , and derive MIR cut
 *    a~*x' <= down(b)
 *  integers :  a~_j = down(a'_j)                      , if f_j <= f_0
 *              a~_j = down(a'_j) + (f_j - f0)/(1 - f0), if f_j >  f_0
 *  continuous: a~_j = 0                               , if a'_j >= 0
 *              a~_j = a'_j/(1 - f0)                   , if a'_j <  0
 *
 *  Transform inequality back to a°*x <= rhs:
 *
 *  (lb or ub):
 *    x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a°_j :=  a~_j,   if lb was used in transformation
 *    x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a°_j := -a~_j,   if ub was used in transformation
 *  and move the constant terms
 *    -a~_j * lb_j == -a°_j * lb_j, or
 *     a~_j * ub_j == -a°_j * ub_j
 *  to the rhs.
 *
 *  (vlb or vub):
 *    x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a°_j :=  a~_j,   (vlb)
 *    x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a°_j := -a~_j,   (vub)
 *  move the constant terms
 *    -a~_j * dl_j == -a°_j * dl_j, or
 *     a~_j * du_j == -a°_j * du_j
 *  to the rhs, and update the VB variable coefficients:
 *    a°_{zl_j} := a°_{zl_j} - a~_j * bl_j == a°_{zl_j} - a°_j * bl_j, or
 *    a°_{zu_j} := a°_{zu_j} + a~_j * bu_j == a°_{zu_j} - a°_j * bu_j
 */
static
void roundMIRRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable (vlb/vub_idx or -1 for lb/ub)*/
   SCIP_Real             f0,                 /**< fracional value of rhs */
   SCIP_Real*            cutactivity         /**< pointer to store the activity of the resulting cut */
   )
{
   SCIP_VAR* var;
   SCIP_Real onedivoneminusf0;
   SCIP_Real aj;
   SCIP_Real downaj;
   SCIP_Real cutaj;
   SCIP_Real fj;
   int nintvars;
   int v;

   assert(prob != NULL);
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(varsign != NULL);
   assert(0.0 < f0 && f0 < 1.0);
   assert(cutactivity != NULL);

   *cutactivity = 0.0;
   onedivoneminusf0 = 1.0 / (1.0 - f0);

   nintvars = prob->nvars - prob->ncontvars;

   /* integer variables */
   for( v = 0; v < nintvars; ++v )
   {
      var = prob->vars[v];
      assert(var != NULL);
      assert(SCIPvarIsIntegral(var));
      assert(SCIPvarGetProbindex(var) == v);
      assert(boundtype[v] == -1 || boundtype[v] == -2);
      assert(varsign[v] == +1 || varsign[v] == -1);

      /* calculate the coefficient in the retransformed cut */
      aj = varsign[v] * mircoef[v]; /* a'_j */
      downaj = SCIPsetFloor(set, aj);
      fj = aj - downaj;
      
      if( SCIPsetIsSumLE(set, fj, f0) )
         cutaj = varsign[v] * downaj; /* a°_j */
      else
         cutaj = varsign[v] * (downaj + (fj - f0) * onedivoneminusf0); /* a°_j */

      if( SCIPsetIsZero(set, cutaj) )
         mircoef[v] = 0.0;
      else
      {
         mircoef[v] = cutaj;
         (*cutactivity) += cutaj * SCIPvarGetLPSol(var);

         /* move the constant term  -a~_j * lb_j == -a°_j * lb_j , or  a~_j * ub_j == -a°_j * ub_j  to the rhs */
         if( varsign[v] == +1 )
         {
            /* lower bound was used */
            if( boundtype[v] == -1 )
            {
               assert(!SCIPsetIsInfinity(set, -SCIPvarGetLbGlobal(var)));
               (*mirrhs) += cutaj * SCIPvarGetLbGlobal(var);
            }
            else
            {
               assert(!SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)));
               (*mirrhs) += cutaj * SCIPvarGetLbLocal(var);
            }
         }
         else
         {
            /* upper bound was used */
            if( boundtype[v] == -1 )
            {
               assert(!SCIPsetIsInfinity(set, SCIPvarGetUbGlobal(var)));
               (*mirrhs) += cutaj * SCIPvarGetUbGlobal(var);
            }
            else
            {
               assert(!SCIPsetIsInfinity(set, SCIPvarGetUbLocal(var)));
               (*mirrhs) += cutaj * SCIPvarGetUbLocal(var);
            }
         }
      }
   }

   /* continuous variables */
   for( v = nintvars; v < prob->nvars; ++v )
   {
      var = prob->vars[v];
      assert(var != NULL);
      assert(!SCIPvarIsIntegral(var));
      assert(SCIPvarGetProbindex(var) == v);
      assert(varsign[v] == +1 || varsign[v] == -1);

      /* calculate the coefficient in the retransformed cut */
      aj = varsign[v] * mircoef[v]; /* a'_j */
      if( aj >= 0.0 )
      {
         mircoef[v] = 0.0;
         continue;
      }
      cutaj = varsign[v] * aj * onedivoneminusf0; /* a°_j */
      if( SCIPsetIsZero(set, cutaj) )
      {
         mircoef[v] = 0.0;
         continue;
      }
      mircoef[v] = cutaj;
      (*cutactivity) += cutaj * SCIPvarGetLPSol(var);
         
      /* check for variable bound use */
      if( boundtype[v] < 0 )
      {
         /* standard bound */

         /* move the constant term  -a~_j * lb_j == -a°_j * lb_j , or  a~_j * ub_j == -a°_j * ub_j  to the rhs */
         if( varsign[v] == +1 )
         {
            /* lower bound was used */
            if( boundtype[v] == -1 )
            {
               assert(!SCIPsetIsInfinity(set, -SCIPvarGetLbGlobal(var)));
               (*mirrhs) += cutaj * SCIPvarGetLbGlobal(var);
            }
            else
            {
               assert(!SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)));
               (*mirrhs) += cutaj * SCIPvarGetLbLocal(var);
            }
         }
         else
         {
            /* upper bound was used */
            if( boundtype[v] == -1 )
            {
               assert(!SCIPsetIsInfinity(set, SCIPvarGetUbGlobal(var)));
               (*mirrhs) += cutaj * SCIPvarGetUbGlobal(var);
            }
            else
            {
               assert(!SCIPsetIsInfinity(set, SCIPvarGetUbLocal(var)));
               (*mirrhs) += cutaj * SCIPvarGetUbLocal(var);
            }
         }
      }
      else
      {
         SCIP_VAR** vbz;
         SCIP_Real* vbb;
         SCIP_Real* vbd;
         int vbidx;
         int zidx;

         /* variable bound */
         vbidx = boundtype[v];

         /* change mirrhs and cutaj of integer variable z_j of variable bound */
         if( varsign[v] == +1 )
         {
            /* variable lower bound was used */
            assert(0 <= vbidx && vbidx < SCIPvarGetNVlbs(var));
            vbz = SCIPvarGetVlbVars(var);
            vbb = SCIPvarGetVlbCoefs(var);
            vbd = SCIPvarGetVlbConstants(var);
         }
         else
         {
            /* variable upper bound was used */
            assert(0 <= vbidx && vbidx < SCIPvarGetNVubs(var));
            vbz = SCIPvarGetVubVars(var);
            vbb = SCIPvarGetVubCoefs(var);
            vbd = SCIPvarGetVubConstants(var);
         }
         zidx =  SCIPvarGetProbindex(vbz[vbidx]);
         assert(0 <= zidx && zidx < v);

         (*mirrhs) += cutaj * vbd[vbidx];
         mircoef[zidx] -= cutaj * vbb[vbidx];
         (*cutactivity) -= cutaj * vbb[vbidx] * SCIPvarGetLPSol(vbz[vbidx]);
      }
   }
}

/** substitute aggregated slack variables:
 *
 *  The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
 *  variable only appears in its own row:
 *     a'_r = scale * weight[r] * slacksign[r].
 *
 *  Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
 *    integers :  a°_r = a~_r = down(a'_r)                      , if f_r <= f0
 *                a°_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
 *    continuous: a°_r = a~_r = 0                               , if a'_r >= 0
 *                a°_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
 *
 *  Substitute a°_r * s_r by adding a°_r times the slack's definition to the cut.
 */
static
void substituteMIRRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   SCIP_Real             f0,                 /**< fracional value of rhs */
   SCIP_Real*            cutactivity         /**< pointer to update the activity of the resulting cut */
   )
{
   SCIP_ROW* row;
   SCIP_Real onedivoneminusf0;
   SCIP_Real ar;
   SCIP_Real downar;
   SCIP_Real cutar;
   SCIP_Real fr;
   SCIP_Real mul;
   int idx;
   int r;
   int i;

   assert(lp != NULL);
   assert(weights != NULL);
   assert(SCIPsetIsPositive(set, scale));
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(slacksign != NULL);
   assert(0.0 < f0 && f0 < 1.0);
   assert(cutactivity != NULL);

   onedivoneminusf0 = 1.0 / (1.0 - f0);
   for( r = 0; r < lp->nrows; ++r )
   {
      /* unused slacks can be ignored */
      if( slacksign[r] == 0 )
         continue;

      assert(slacksign[r] == -1 || slacksign[r] == +1);
      assert(!SCIPsetIsZero(set, weights[r]));

      row = lp->rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* get the slack's coefficient a'_r in the aggregated row */
      ar = slacksign[r] * scale * weights[r];

      /* calculate slack variable's coefficient a°_r in the cut */
      if( row->integral
         && ((slacksign[r] == +1 && SCIPsetIsFeasIntegral(set, row->rhs - row->constant))
            || (slacksign[r] == -1 && SCIPsetIsFeasIntegral(set, row->lhs - row->constant))) )
      {
         /* slack variable is always integral:
          *    a°_r = a~_r = down(a'_r)                      , if f_r <= f0
          *    a°_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
          */
         downar = SCIPsetFloor(set, ar);
         fr = ar - downar;
         if( SCIPsetIsLE(set, fr, f0) )
            cutar = downar;
         else
            cutar = downar + (fr - f0) * onedivoneminusf0;
      }
      else
      {
         /* slack variable is continuous:
          *    a°_r = a~_r = 0                               , if a'_r >= 0
          *    a°_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
          */
         if( ar >= 0.0 )
            continue; /* slack can be ignored, because its coefficient is reduced to 0.0 */
         else
            cutar = ar * onedivoneminusf0;
      }

      /* if the coefficient was reduced to zero, ignore the slack variable */
      if( SCIPsetIsZero(set, cutar) )
         continue;

      /* depending on the slack's sign, we have
       *   a*x + c + s == rhs  =>  s == - a*x - c + rhs,  or  a*x + c - s == lhs  =>  s == a*x + c - lhs
       * substitute a°_r * s_r by adding a°_r times the slack's definition to the cut.
       */
      mul = -slacksign[r] * cutar;

      /* add the slack's definition multiplied with a°_j to the cut */
      for( i = 0; i < row->len; ++i )
      {
         assert(row->cols[i] != NULL);
         assert(row->cols[i]->var != NULL);
         assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetCol(row->cols[i]->var) == row->cols[i]);
         assert(SCIPvarGetProbindex(row->cols[i]->var) == row->cols[i]->var_probindex);
         idx = row->cols[i]->var_probindex;
         mircoef[idx] += mul * row->vals[i];
      }

      /* update the activity: we have to add  mul * a*x^  to the cut's activity (row activity = a*x^ + c) */
      (*cutactivity) += mul * (SCIProwGetLPActivity(row, stat, lp) - row->constant);

      /* move slack's constant to the right hand side */
      if( slacksign[r] == +1 )
      {
         /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a°_r * (rhs - c) to the right hand side */
         assert(!SCIPsetIsInfinity(set, row->rhs));
         (*mirrhs) -= cutar * (row->rhs - row->constant);
      }
      else
      {
         /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a°_r * (c - lhs) to the right hand side */
         assert(!SCIPsetIsInfinity(set, -row->lhs));
         (*mirrhs) -= cutar * (row->constant - row->lhs);
      }
   }

   /* set rhs to zero, if it's very close to */
   if( SCIPsetIsZero(set, *mirrhs) )
      *mirrhs = 0.0;
}

#ifdef SCIP_DEBUG
static
void printMIR(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real*            mircoef,            /**< MIR coefficients */
   SCIP_Real             mirrhs              /**< right hand side of the MIR row */
   )
{
   int i;

   assert(prob != NULL);

   SCIPdebugPrintf("MIR:");
   for( i = 0; i < prob->nvars; ++i )
   {
      if( mircoef[i] != 0.0 )
         SCIPdebugPrintf(" %+g<%s>", mircoef[i], SCIPvarGetName(prob->vars[i]));
   }
   SCIPdebugPrintf(" <= %g\n", mirrhs);
}
#endif

/* calculates a MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a MIR cut.
 */
SCIP_RETCODE SCIPlpCalcMIR(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   SCIP_Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid MIR cut */
   SCIP_Bool*            cutislocal          /**< pointer to store whether the returned cut is only valid locally */
   )
{
   int* slacksign;
   int* varsign;
   int* boundtype;
   SCIP_Real rhs;
   SCIP_Real downrhs;
   SCIP_Real f0;
   SCIP_Bool emptyrow;
   SCIP_Bool freevariable;
   SCIP_Bool localrowsused;
   SCIP_Bool localbdsused;

   assert(lp != NULL);
   assert(lp->solved);
   assert(prob != NULL);
   assert(weights != NULL);
   assert(SCIPsetIsPositive(set, scale));
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(cutactivity != NULL);
   assert(success != NULL);
   assert(cutislocal != NULL);

   SCIPdebugMessage("calculating MIR cut (scale: %g)\n", scale);

   /**@todo test, if a column based summation is faster */

   *success = FALSE;
   *cutislocal = FALSE;

   /* allocate temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &slacksign, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varsign, prob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &boundtype, prob->nvars) );

   /* calculate the row summation */
   sumMIRRow(set, prob, lp, weights, scale, allowlocal, 
      maxweightrange, mircoef, &rhs, slacksign, &emptyrow, &localrowsused);
   assert(allowlocal || !localrowsused);
   *cutislocal = *cutislocal || localrowsused;
   if( emptyrow )
      goto TERMINATE;
   SCIPdebug(printMIR(prob, mircoef, rhs));
   
   /* remove all nearly-zero coefficients from MIR row and relaxes the right hand side correspondingly in order to
    *  prevent numerical rounding errors
    */
   cleanupMIRRow(set, prob, mircoef, &rhs, *cutislocal);
   SCIPdebug(printMIR(prob, mircoef, rhs));

   /* Transform equation  a*x == b, lb <= x <= ub  into standard form
    *   a'*x' == b, 0 <= x' <= ub'.
    * 
    * Transform variables (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
    * and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
    *
    * Transform variables (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
    * move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
    *   a_{zl_j} := a_{zl_j} + a_j * bl_j, or
    *   a_{zu_j} := a_{zu_j} + a_j * bu_j
    */
   transformMIRRow(set, prob, boundswitch, usevbds, allowlocal, 
      mircoef, &rhs, varsign, boundtype, &freevariable, &localbdsused);
   assert(allowlocal || !localbdsused);
   *cutislocal = *cutislocal || localbdsused;
   if( freevariable )
      goto TERMINATE;
   SCIPdebug(printMIR(prob, mircoef, rhs));

   /* Calculate fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j) , and derive MIR cut
    *   a~*x' <= down(b)
    * integers :  a~_j = down(a'_j)                      , if f_j <= f_0
    *             a~_j = down(a'_j) + (f_j - f0)/(1 - f0), if f_j >  f_0
    * continuous: a~_j = 0                               , if a'_j >= 0
    *             a~_j = a'_j/(1 - f0)                   , if a'_j <  0
    *
    * Transform inequality back to a°*x <= rhs:
    *
    * (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a°_j :=  a~_j,   if lb was used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a°_j := -a~_j,   if ub was used in transformation
    * and move the constant terms
    *   -a~_j * lb_j == -a°_j * lb_j, or
    *    a~_j * ub_j == -a°_j * ub_j
    * to the rhs.
    *
    * (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a°_j :=  a~_j,   (vlb)
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a°_j := -a~_j,   (vub)
    * move the constant terms
    *   -a~_j * dl_j == -a°_j * dl_j, or
    *    a~_j * du_j == -a°_j * du_j
    * to the rhs, and update the VB variable coefficients:
    *   a°_{zl_j} := a°_{zl_j} - a~_j * bl_j == a°_{zl_j} - a°_j * bl_j, or
    *   a°_{zu_j} := a°_{zu_j} + a~_j * bu_j == a°_{zu_j} - a°_j * bu_j
    */
   downrhs = SCIPsetFloor(set, rhs);
   f0 = rhs - downrhs;
   if( f0 < minfrac )
      goto TERMINATE;

   *mirrhs = downrhs;
   roundMIRRow(set, prob, mircoef, mirrhs, varsign, boundtype, f0, cutactivity);
   SCIPdebug(printMIR(prob, mircoef, *mirrhs));

   /* substitute aggregated slack variables:
    *
    * The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
    * variable only appears in its own row:
    *    a'_r = scale * weight[r] * slacksign[r].
    *
    * Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
    *   integers :  a°_r = a~_r = down(a'_r)                      , if f_r <= f0
    *               a°_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
    *   continuous: a°_r = a~_r = 0                               , if a'_r >= 0
    *               a°_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
    *
    * Substitute a°_r * s_r by adding a°_r times the slack's definition to the cut.
    */
   substituteMIRRow(set, stat, lp, weights, scale, mircoef, mirrhs, slacksign, f0, cutactivity);
   SCIPdebug(printMIR(prob, mircoef, *mirrhs));

   *success = TRUE;

#ifndef NDEBUG
   {
      SCIP_Real act;
      int i;

      act = 0.0;
      for( i = 0; i < prob->nvars; ++i )
         act += mircoef[i] * SCIPvarGetLPSol(prob->vars[i]);
      assert(SCIPsetIsFeasEQ(set, act/10.0, *cutactivity/10.0));
   }
#endif

 TERMINATE:
   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &boundtype);
   SCIPsetFreeBufferArray(set, &varsign);
   SCIPsetFreeBufferArray(set, &slacksign);

   return SCIP_OKAY;
}

/** builds a weighted sum of rows, and decides whether to use the left or right hand side of the rows in summation */
static
void sumStrongCGRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Bool             allowlocal,         /**< should local rows be included, resulting in a locally valid summation? */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size prob->nvars */
   SCIP_Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   SCIP_Bool*            emptyrow,           /**< pointer to store whether the returned row is empty */
   SCIP_Bool*            localrowsused       /**< pointer to store whether local rows were used in summation */
   )
{
   SCIP_ROW* row;
   SCIP_Real weight;
   SCIP_Real absweight;
   SCIP_Real maxweight;
   int idx;
   int r;
   int i;

   assert(prob != NULL);
   assert(lp != NULL);
   assert(weights != NULL);
   assert(SCIPsetIsPositive(set, scale));
   assert(maxweightrange >= 1.0);
   assert(strongcgcoef != NULL);
   assert(strongcgrhs != NULL);
   assert(slacksign != NULL);
   assert(emptyrow != NULL);
   assert(localrowsused != NULL);

   /* search the maximal absolute weight */
   maxweight = 0.0;
   for( r = 0; r < lp->nrows; ++r )
   {
      weight = scale * weights[r];
      absweight = ABS(weight);
      maxweight = MAX(maxweight, absweight);
   }

   /* calculate the row summation */
   BMSclearMemoryArray(strongcgcoef, prob->nvars);
   *strongcgrhs = 0.0;
   *emptyrow = TRUE;
   *localrowsused = FALSE;
   for( r = 0; r < lp->nrows; ++r )
   {
      row = lp->rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* modifiable rows cannot be part of a strong CG row summation;
       * local rows are only included, if the allowlocal flag is set;
       * close to zero weights or weights outside the maximal range are ignored
       */
      weight = scale * weights[r];
      absweight = ABS(weight);
      if( !row->modifiable && (allowlocal || !row->local)
         && absweight * maxweightrange >= maxweight && !SCIPsetIsSumZero(set, weight) )
      {
         *emptyrow = FALSE;
         *localrowsused = *localrowsused || row->local;

         if( row->integral )
         { 
            /* Row is integral:
             * Decide, if we want to use the left or the right hand side of the row in the summation.
             * If possible, use the side that leads to a positive slack value in the summation. 
             */
            if( SCIPsetIsInfinity(set, row->rhs) || (!SCIPsetIsInfinity(set, -row->lhs) && weight < 0.0) )
            {
               slacksign[r] = -1;
               (*strongcgrhs) += weight * SCIPsetFeasCeil(set, row->lhs - row->constant); /* row is integral: round left hand side up */
            }
            else
            {
               slacksign[r] = +1;
               (*strongcgrhs) += weight * SCIPsetFeasFloor(set, row->rhs - row->constant); /* row is integral: round right hand side down */
            }
         }
         else
         {
            /* Row is NOT integral:
             * Decide, if we have to use the left or the right hand side of the row in the summation,
             * in order to get a positive slack variable in the summation.
             * If not possible, ignore row in summation.
             */
            if( weight > 0.0 && !SCIPsetIsInfinity(set, row->rhs) )
            {
               slacksign[r] = +1;
               (*strongcgrhs) += weight * (row->rhs - row->constant);
            }
            else if( weight < 0.0 && !SCIPsetIsInfinity(set, -row->lhs) )
            {
               slacksign[r] = -1;
               (*strongcgrhs) += weight * (row->lhs - row->constant);
            }
            else
            {
               slacksign[r] = 0;
               weights[r] = 0.0;
               continue;
            }
         }

         /* add the row coefficients to the sum */
         for( i = 0; i < row->len; ++i )
         {
            assert(row->cols[i] != NULL);
            assert(row->cols[i]->var != NULL);
            assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
            assert(SCIPvarGetCol(row->cols[i]->var) == row->cols[i]);
            assert(SCIPvarGetProbindex(row->cols[i]->var) == row->cols[i]->var_probindex);
            idx = row->cols[i]->var_probindex;
            assert(0 <= idx && idx < prob->nvars);
            strongcgcoef[idx] += weight * row->vals[i];
         }

         SCIPdebugMessage("strong CG: %d: row <%s>, lhs = %g, rhs = %g, scale = %g, weight = %g, slacksign = %d -> rhs = %g\n",
            r, SCIProwGetName(row), row->lhs - row->constant, row->rhs - row->constant, 
            scale, weights[r], slacksign[r], *strongcgrhs);
         SCIPdebug(SCIProwPrint(row, NULL));
      }
      else
      {
         slacksign[r] = 0;
         weights[r] = 0.0;
      }
   }
}

/** Transform equation  a*x == b, lb <= x <= ub  into standard form
 *    a'*x' == b, 0 <= x' <= ub'.
 *  
 *  Transform variables (lb or ub):
 *    x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
 *    x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
 *  and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
 *
 *  Transform variables (vlb or vub):
 *    x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
 *    x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
 *  move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
 *    a_{zl_j} := a_{zl_j} + a_j * bl_j, or
 *    a_{zu_j} := a_{zu_j} + a_j * bu_j
 */
static
void transformStrongCGRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size nvars */
   SCIP_Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Bool*            freevariable,       /**< stores whether a free variable was found in strong CG row -> invalid summation */
   SCIP_Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{
   SCIP_VAR* var;
   SCIP_Real varsol;
   SCIP_Real bestlb;
   SCIP_Real bestub;
   int bestlbtype;
   int bestubtype;
   int v;

   assert(prob != NULL);
   assert(strongcgcoef != NULL);
   assert(strongcgrhs != NULL);
   assert(varsign != NULL);
   assert(boundtype != NULL);
   assert(freevariable != NULL);
   assert(localbdsused != NULL);

   *freevariable = FALSE;
   *localbdsused = FALSE;
   
   /* substitute continuous variables with best standard or variable bound (lb, ub, vlb or vub),
    * substitute integral variables with best standard bound (lb, ub);
    * start with continuous variables, because using variable bounds can affect the untransformed integral
    * variables, and these changes have to be incorporated in the transformation of the integral variables
    */
   for( v = prob->nvars-1; v >= 0; --v )
   {
      var = prob->vars[v];
      assert(v == SCIPvarGetProbindex(var));

      /* ignore variables that don't exist in the strong CG row */
      if( SCIPsetIsZero(set, strongcgcoef[v]) )
      {
         varsign[v] = +1;
         boundtype[v] = -2;
         continue;
      }

      /* find closest lower bound in standard lower bound (and variable lower bounds for continuous variables) */
      bestlb = SCIPvarGetLbGlobal(var);
      bestlbtype = -1;
      if( allowlocal )
      {
         SCIP_Real loclb;

         loclb = SCIPvarGetLbLocal(var);
         if( SCIPsetIsGT(set, loclb, bestlb) )
         {
            bestlb = loclb;
            bestlbtype = -2;
         }
      }
      if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_Real bestvlb;
         int bestvlbidx;

         getClosestVlb(var, &bestvlb, &bestvlbidx);
         if( bestvlbidx >= 0
            && (bestvlb > bestlb || (bestlbtype < 0 && SCIPsetIsGE(set, bestvlb, bestlb))) )
         {
            bestlb = bestvlb;
            bestlbtype = bestvlbidx;
         }
      }

      /* find closest upper bound in standard upper bound (and variable upper bounds for continuous variables) */
      bestub = SCIPvarGetUbGlobal(var);
      bestubtype = -1;
      if( allowlocal )
      {
         SCIP_Real locub;

         locub = SCIPvarGetUbLocal(var);
         if( SCIPsetIsLT(set, locub, bestub) )
         {
            bestub = locub;
            bestubtype = -2;
         }
      }
      if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_Real bestvub;
         int bestvubidx;

         getClosestVub(var, &bestvub, &bestvubidx);
         if( bestvubidx >= 0
            && (bestvub < bestub || (bestubtype < 0 && SCIPsetIsLE(set, bestvub, bestub))) )
         {
            bestub = bestvub;
            bestubtype = bestvubidx;
         }
      }

      /* check, if variable is free variable 
       * (for continuous variable use bound, so that coefficient will be nonnegative ) */
      if( ( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS 
            && SCIPsetIsInfinity(set, -bestlb) && SCIPsetIsInfinity(set, bestub) )
         || ( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS 
            && ( ( strongcgcoef[v] > 0.0 && SCIPsetIsInfinity(set, -bestlb) )
               || ( strongcgcoef[v] < 0.0 && SCIPsetIsInfinity(set, bestub) ) ) ) )
      {
         /* we found a free variable in the row with non-zero coefficient
          *  -> strong CG row can't be transformed in standard form
          */
         *freevariable = TRUE;
         return;
      }

      /* select transformation bound 
       * (for continuous variable use bound, so that coefficient will be nonnegative ) */
      varsol = SCIPvarGetLPSol(var);

      if( ( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS && strongcgcoef[v] > 0.0 ) 
         || ( SCIPvarGetType(var)!= SCIP_VARTYPE_CONTINUOUS && 
            varsol <= (1.0 - boundswitch) * bestlb + boundswitch * bestub ) )
      {
         /* use lower bound as transformation bound: x'_j := x_j - lb_j */
         boundtype[v] = bestlbtype;
         varsign[v] = +1;
    
         /* standard (bestlbtype < 0) or variable (bestlbtype >= 0) lower bound? */
         if( bestlbtype < 0 )
         {
            (*strongcgrhs) -= strongcgcoef[v] * bestlb;
            *localbdsused = *localbdsused || (bestlbtype == -2);
         }
         else
         {
            SCIP_VAR** vlbvars = SCIPvarGetVlbVars(var);
            SCIP_Real* vlbcoefs = SCIPvarGetVlbCoefs(var);
            SCIP_Real* vlbconsts = SCIPvarGetVlbConstants(var);
            int zidx;

            assert(0 <= bestlbtype && bestlbtype < SCIPvarGetNVlbs(var));
            zidx = SCIPvarGetProbindex(vlbvars[bestlbtype]);
            assert(0 <= zidx && zidx < v);
               
            (*strongcgrhs) -= strongcgcoef[v] * vlbconsts[bestlbtype];
            strongcgcoef[zidx] += strongcgcoef[v] * vlbcoefs[bestlbtype];
         }
      }
      else
      {
         /* use upper bound as transformation bound: x'_j := ub_j - x_j */
         boundtype[v] = bestubtype;
         varsign[v] = -1;
         
         /* standard (bestubtype < 0) or variable (bestubtype >= 0) upper bound? */
         if( bestubtype < 0 )
         {
            (*strongcgrhs) -= strongcgcoef[v] * bestub;
            *localbdsused = *localbdsused || (bestubtype == -2);
         }
         else
         {
            SCIP_VAR** vubvars = SCIPvarGetVubVars(var);
            SCIP_Real* vubcoefs = SCIPvarGetVubCoefs(var);
            SCIP_Real* vubconsts = SCIPvarGetVubConstants(var);
            int zidx;

            assert(0 <= bestubtype && bestubtype < SCIPvarGetNVubs(var));
            zidx = SCIPvarGetProbindex(vubvars[bestubtype]);
            assert(zidx >= 0);
               
            (*strongcgrhs) -= strongcgcoef[v] * vubconsts[bestubtype];
            strongcgcoef[zidx] += strongcgcoef[v] * vubcoefs[bestubtype];
         }
      }

      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         assert(strongcgcoef[v]*varsign[v] > 0.0);

      SCIPdebugMessage("strong CG var <%s>: varsign=%d, boundtype=%d, strongcgcoef=%g, lb=%g, ub=%g -> rhs=%g\n", 
         SCIPvarGetName(var), varsign[v], boundtype[v], strongcgcoef[v], bestlb, bestub, *strongcgrhs);
   }
}

/** Calculate 
 *   fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j) and 
 *   integer k >= 1 with 1/(k + 1) <= f_0 < 1/k and 
 *   (=> k = up(1/f_0) + 1)
 *   integer 1 <= p_j <= k with f_0 + ((p_j - 1) * (1 - f_0)/k) < f_j <= f_0 + (p_j * (1 - f_0)/k) 
 *   (=> p_j = up( k*(f_j - f_0)/(1 - f_0) ))
 * and derive strong CG cut 
 *   a~*x' <= down(b)
 * integers :  a~_j = down(a'_j)                , if f_j <= f_0
 *             a~_j = down(a'_j) + p_j/(k + 1)  , if f_j >  f_0  
 * continuous: a~_j = 0                         , if a'_j >= 0
 *             no strong CG cut found          , if a'_j <  0
 *
 * Transform inequality back to a°*x <= rhs:
 *
 *  (lb or ub):
 *    x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a°_j :=  a~_j,   if lb was used in transformation
 *    x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a°_j := -a~_j,   if ub was used in transformation
 *  and move the constant terms
 *    -a~_j * lb_j == -a°_j * lb_j, or
 *     a~_j * ub_j == -a°_j * ub_j
 *  to the rhs.
 *
 *  (vlb or vub):
 *    x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a°_j :=  a~_j,   (vlb)
 *    x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a°_j := -a~_j,   (vub)
 *  move the constant terms
 *    -a~_j * dl_j == -a°_j * dl_j, or
 *     a~_j * du_j == -a°_j * du_j
 *  to the rhs, and update the VB variable coefficients:
 *    a°_{zl_j} := a°_{zl_j} - a~_j * bl_j == a°_{zl_j} - a°_j * bl_j, or
 *    a°_{zu_j} := a°_{zu_j} + a~_j * bu_j == a°_{zu_j} - a°_j * bu_j
 */
static
void roundStrongCGRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size nvars */
   SCIP_Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable (vlb/vub_idx or -1 for lb/ub)*/
   SCIP_Real             f0,                 /**< fracional value of rhs */
   SCIP_Real             k,                  /**< factor to strengthen strongcg cut */
   SCIP_Real*            cutactivity         /**< pointer to store the activity of the resulting cut */
   )
{
   SCIP_VAR* var;
   SCIP_Real onedivoneminusf0;
   SCIP_Real pj;
   SCIP_Real aj;
   SCIP_Real downaj;
   SCIP_Real cutaj;
   SCIP_Real fj;
   int nintvars;
   int v;

   assert(prob != NULL);
   assert(strongcgcoef != NULL);
   assert(strongcgrhs != NULL);
   assert(varsign != NULL);
   assert(0.0 < f0 && f0 < 1.0);
   assert(cutactivity != NULL);
 
   *cutactivity = 0.0;
   onedivoneminusf0 = 1.0 / (1.0 - f0);
   nintvars = prob->nvars - prob->ncontvars;

   /* integer variables */
   for( v = 0; v < nintvars; ++v )
   {
      var = prob->vars[v];
      assert(var != NULL);
      assert(SCIPvarIsIntegral(var));
      assert(SCIPvarGetProbindex(var) == v);
      assert(boundtype[v] == -1 || boundtype[v] == -2);
      assert(varsign[v] == +1 || varsign[v] == -1);
      
      /* calculate the coefficient in the retransformed cut */
      aj = varsign[v] * strongcgcoef[v]; /* a'_j */
      downaj = SCIPsetFloor(set, aj);
      fj = aj - downaj;

      if( SCIPsetIsSumLE(set, fj, f0) )
         cutaj = varsign[v] * downaj; /* a°_j */
      else
      {
         pj = SCIPsetCeil(set, k * (fj - f0) * onedivoneminusf0);
         assert(pj >= 0); /* should be >= 1, but due to rounding bias can be 0 if fj almost equal to f0 */ 
         assert(pj <= k);
         cutaj = varsign[v] * (downaj + pj / (k + 1)); /* a°_j */
      }
      if( SCIPsetIsZero(set, cutaj) )
         strongcgcoef[v] = 0.0;
      else
      {
         strongcgcoef[v] = cutaj;
         (*cutactivity) += cutaj * SCIPvarGetLPSol(var);

         /* move the constant term  -a~_j * lb_j == -a°_j * lb_j , or  a~_j * ub_j == -a°_j * ub_j  to the rhs */
         if( varsign[v] == +1 )
         {
            /* lower bound was used */
            if( boundtype[v] == -1 )
            {
               assert(!SCIPsetIsInfinity(set, -SCIPvarGetLbGlobal(var)));
               (*strongcgrhs) += cutaj * SCIPvarGetLbGlobal(var);
            }
            else
            {
               assert(!SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)));
               (*strongcgrhs) += cutaj * SCIPvarGetLbLocal(var);
            }
         }
         else
         {
            /* upper bound was used */
            if( boundtype[v] == -1 )
            {
               assert(!SCIPsetIsInfinity(set, SCIPvarGetUbGlobal(var)));
               (*strongcgrhs) += cutaj * SCIPvarGetUbGlobal(var);
            }
            else
            {
               assert(!SCIPsetIsInfinity(set, SCIPvarGetUbLocal(var)));
               (*strongcgrhs) += cutaj * SCIPvarGetUbLocal(var);
            }
         }
      }
   }

   /* continuous variables */
   for( v = nintvars; v < prob->nvars; ++v )
   {
      var = prob->vars[v];
      assert(var != NULL);
      assert(!SCIPvarIsIntegral(var));
      assert(SCIPvarGetProbindex(var) == v);
      assert(varsign[v] == +1 || varsign[v] == -1);

      /* calculate the coefficient in the retransformed cut */
      aj = varsign[v] * strongcgcoef[v]; /* a'_j */

      assert( aj >= 0.0 );
      strongcgcoef[v] = 0.0;
   }
}

/** substitute aggregated slack variables:
 *
 *  The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
 *  variable only appears in its own row:
 *     a'_r = scale * weight[r] * slacksign[r].
 *
 *  Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
 *    integers:   a°_r = a~_r = down(a'_r)                  , if f_r <= f0
 *                a°_r = a~_r = down(a'_r) + p_r/(k + 1)    , if f_r >  f0
 *    continuous: a°_r = a~_r = 0                           , if a'_r >= 0
 *                no strong CG cut found                    , if a'_r <  0
 *
 *  Substitute a°_r * s_r by adding a°_r times the slack's definition to the cut.
 */
static
void substituteStrongCGRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size nvars */
   SCIP_Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   SCIP_Real             f0,                 /**< fracional value of rhs */
   SCIP_Real             k,                  /**< factor to strengthen strongcg cut */
   SCIP_Real*            cutactivity         /**< pointer to update the activity of the resulting cut */
   )
{
   SCIP_ROW* row;
   SCIP_Real onedivoneminusf0;
   SCIP_Real pr;
   SCIP_Real ar;
   SCIP_Real downar;
   SCIP_Real cutar;
   SCIP_Real fr;
   SCIP_Real mul;
   int idx;
   int r;
   int i;

   assert(lp != NULL);
   assert(weights != NULL);
   assert(SCIPsetIsPositive(set, scale));
   assert(strongcgcoef != NULL);
   assert(strongcgrhs != NULL);
   assert(slacksign != NULL);
   assert(0.0 < f0 && f0 < 1.0);
   assert(cutactivity != NULL);
 
   onedivoneminusf0 = 1.0 / (1.0 - f0);
  
   for( r = 0; r < lp->nrows; ++r )
   {
      /* unused slacks can be ignored */
      if( slacksign[r] == 0 )
         continue;

      assert(slacksign[r] == -1 || slacksign[r] == +1);
      assert(!SCIPsetIsZero(set, weights[r]));

      row = lp->rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* get the slack's coefficient a'_r in the aggregated row */
      ar = slacksign[r] * scale * weights[r];

      /* calculate slack variable's coefficient a°_r in the cut */
      if( row->integral )
      {
         /* slack variable is always integral: */
         downar = SCIPsetFloor(set, ar);
         fr = ar - downar;

         if( SCIPsetIsLE(set, fr, f0) )
            cutar = downar;
         else
         {
            pr = SCIPsetCeil(set, k * (fr - f0) * onedivoneminusf0);
            assert(pr >= 0); /* should be >= 1, but due to rounding bias can be 0 if fr almost equal to f0 */ 
            assert(pr <= k);
            cutar = downar + pr / (k + 1);
         }
      }
      else
      {
         /* slack variable is continuous: */
         assert(ar >= 0.0);
         continue; /* slack can be ignored, because its coefficient is reduced to 0.0 */
      }

      /* if the coefficient was reduced to zero, ignore the slack variable */
      if( SCIPsetIsZero(set, cutar) )
         continue;

      /* depending on the slack's sign, we have
       *   a*x + c + s == rhs  =>  s == - a*x - c + rhs,  or  a*x + c - s == lhs  =>  s == a*x + c - lhs
       * substitute a°_r * s_r by adding a°_r times the slack's definition to the cut.
       */
      mul = -slacksign[r] * cutar;

      /* add the slack's definition multiplied with a°_j to the cut */
      for( i = 0; i < row->len; ++i )
      {
         assert(row->cols[i] != NULL);
         assert(row->cols[i]->var != NULL);
         assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetCol(row->cols[i]->var) == row->cols[i]);
         assert(SCIPvarGetProbindex(row->cols[i]->var) == row->cols[i]->var_probindex);
         idx = row->cols[i]->var_probindex;
         strongcgcoef[idx] += mul * row->vals[i];
      }

      /* update the activity: we have to add  mul * a*x^  to the cut's activity (row activity = a*x^ + c) */
      *cutactivity += mul * (SCIProwGetLPActivity(row, stat, lp) - row->constant);

      /* move slack's constant to the right hand side */
      if( slacksign[r] == +1 )
      {
 	 SCIP_Real rhs;

         /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a°_r * (rhs - c) to the right hand side */
         assert(!SCIPsetIsInfinity(set, row->rhs));
	 rhs = row->rhs - row->constant;
	 if( row->integral )
	 {
	    /* the right hand side was implicitly rounded down in row aggregation */
	    rhs = SCIPsetFeasFloor(set, rhs);
	 }
         *strongcgrhs -= cutar * rhs;
	 
      }
      else
      {
 	 SCIP_Real lhs;

         /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a°_r * (c - lhs) to the right hand side */
         assert(!SCIPsetIsInfinity(set, -row->lhs));
	 lhs = row->lhs - row->constant;
         if( row->integral )
	 {
	    /* the left hand side was implicitly rounded up in row aggregation */
	    lhs = SCIPsetFeasCeil(set, lhs);
	 }
	 *strongcgrhs += cutar * lhs;
      }
   }

   /* set rhs to zero, if it's very close to */
   if( SCIPsetIsZero(set, *strongcgrhs) )
      *strongcgrhs = 0.0;
}

/* calculates a strong CG cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a strong CG cut.
 */
SCIP_RETCODE SCIPlpCalcStrongCG(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce strong CG cut for */
   SCIP_Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size nvars */
   SCIP_Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   SCIP_Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid strong CG cut */
   SCIP_Bool*            cutislocal          /**< pointer to store whether the returned cut is only valid locally */
   )
{
   int* slacksign;
   int* varsign;
   int* boundtype;
   SCIP_Real rhs;
   SCIP_Real downrhs;
   SCIP_Real f0;
   SCIP_Real k;
   SCIP_Bool emptyrow;
   SCIP_Bool freevariable;
   SCIP_Bool localrowsused;
   SCIP_Bool localbdsused;

   assert(lp != NULL);
   assert(lp->solved);
   assert(prob != NULL);
   assert(weights != NULL);
   assert(SCIPsetIsPositive(set, scale));
   assert(strongcgcoef != NULL);
   assert(strongcgrhs != NULL);
   assert(cutactivity != NULL);
   assert(success != NULL);
   assert(cutislocal != NULL);

   SCIPdebugMessage("calculating strong CG cut (scale: %g)\n", scale);

   /**@todo test, if a column based summation is faster */

   *success = FALSE;
   *cutislocal = FALSE;

   /* allocate temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &slacksign, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varsign, prob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &boundtype, prob->nvars) );

   /* calculate the row summation */
   sumStrongCGRow(set, prob, lp, weights, scale, allowlocal, 
      maxweightrange, strongcgcoef, &rhs, slacksign, &emptyrow, &localrowsused);
   assert(allowlocal || !localrowsused);
   *cutislocal = *cutislocal || localrowsused;
   if( emptyrow )
      goto TERMINATE;
   SCIPdebug(printMIR(prob, strongcgcoef, rhs));
   
   /* remove all nearly-zero coefficients from strong CG row and relaxes the right hand side correspondingly in order to
    *  prevent numerical rounding errors
    */
   cleanupMIRRow(set, prob, strongcgcoef, &rhs, *cutislocal);
   SCIPdebug(printMIR(prob, strongcgcoef, rhs));

   /* Transform equation  a*x == b, lb <= x <= ub  into standard form
    *   a'*x' == b, 0 <= x' <= ub'.
    * 
    * Transform variables (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
    * and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
    *
    * Transform variables (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
    * move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
    *   a_{zl_j} := a_{zl_j} + a_j * bl_j, or
    *   a_{zu_j} := a_{zu_j} + a_j * bu_j
    */
   transformStrongCGRow(set, prob, boundswitch, usevbds, allowlocal, 
      strongcgcoef, &rhs, varsign, boundtype, &freevariable, &localbdsused);
   assert(allowlocal || !localbdsused);
   *cutislocal = *cutislocal || localbdsused;
   if( freevariable )
      goto TERMINATE;
   SCIPdebug(printMIR(prob, strongcgcoef, rhs));

   /* Calculate 
    *   fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j)  
    *   integer k >= 1 with 1/(k + 1) <= f_0 < 1/k  
    *   (=> k = up(1/f_0) + 1)
    *   integer 1 <= p_j <= k with f_0 + ((p_j - 1) * (1 - f_0)/k) < f_j <= f_0 + (p_j * (1 - f_0)/k) 
    *   (=> p_j = up( (f_j - f_0)/((1 - f_0)/k) ))
    * and derive strong CG cut 
    *   a~*x' <= (k+1) * down(b)
    * integers :  a~_j = down(a'_j)                , if f_j <= f_0
    *             a~_j = down(a'_j) + p_j/(k + 1)  , if f_j >  f_0  
    * continuous: a~_j = 0                         , if a'_j >= 0
    *             no strong CG cut found          , if a'_j <  0 
    *
    * Transform inequality back to a°*x <= rhs:
    *
    * (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a°_j :=  a~_j,   if lb was used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a°_j := -a~_j,   if ub was used in transformation
    * and move the constant terms
    *   -a~_j * lb_j == -a°_j * lb_j, or
    *    a~_j * ub_j == -a°_j * ub_j
    * to the rhs.
    *
    * (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a°_j :=  a~_j,   (vlb)
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a°_j := -a~_j,   (vub)
    * move the constant terms
    *   -a~_j * dl_j == -a°_j * dl_j, or
    *    a~_j * du_j == -a°_j * du_j
    * to the rhs, and update the VB variable coefficients:
    *   a°_{zl_j} := a°_{zl_j} - a~_j * bl_j == a°_{zl_j} - a°_j * bl_j, or
    *   a°_{zu_j} := a°_{zu_j} + a~_j * bu_j == a°_{zu_j} - a°_j * bu_j
    */
   downrhs = SCIPsetFloor(set, rhs);
   f0 = rhs - downrhs;
   if( f0 < minfrac )
      goto TERMINATE;
   k = SCIPsetCeil(set, 1.0 / f0) - 1;

   *strongcgrhs = downrhs;
   roundStrongCGRow(set, prob, strongcgcoef, strongcgrhs, varsign, boundtype, f0, k, cutactivity);
   SCIPdebug(printMIR(prob, strongcgcoef, *strongcgrhs));

   /* substitute aggregated slack variables:
    *
    * The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
    * variable only appears in its own row:
    *    a'_r = scale * weight[r] * slacksign[r].
    *
    * Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
    *   integers :  a°_r = a~_r = (k + 1) * down(a'_r)        , if f_r <= f0
    *               a°_r = a~_r = (k + 1) * down(a'_r) + p_r  , if f_r >  f0
    *   continuous: a°_r = a~_r = 0                           , if a'_r >= 0
    *               a°_r = a~_r = a'_r/(1 - f0)               , if a'_r <  0
    *
    * Substitute a°_r * s_r by adding a°_r times the slack's definition to the cut.
    */
   substituteStrongCGRow(set, stat, lp, weights, scale, strongcgcoef, strongcgrhs, slacksign, f0, k, cutactivity);
   SCIPdebug(printMIR(prob, strongcgcoef, *strongcgrhs));

   *success = TRUE;
   
#ifndef NDEBUG
   {
      SCIP_Real act;
      int i;

      act = 0.0;
      for( i = 0; i < prob->nvars; ++i )
         act += strongcgcoef[i] * SCIPvarGetLPSol(prob->vars[i]);
      assert(SCIPsetIsFeasEQ(set, act, *cutactivity));
   }
#endif

 TERMINATE:
   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &boundtype);
   SCIPsetFreeBufferArray(set, &varsign);
   SCIPsetFreeBufferArray(set, &slacksign);

   return SCIP_OKAY;
}

/** stores LP state (like basis information) into LP state object */
SCIP_RETCODE SCIPlpGetState(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(blkmem != NULL);
   assert(lpistate != NULL);

   SCIP_CALL( SCIPlpiGetState(lp->lpi, blkmem, lpistate) );

   return SCIP_OKAY;
}

/** loads LP state (like basis information) into solver */
SCIP_RETCODE SCIPlpSetState(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{
   assert(lp != NULL);
   assert(blkmem != NULL);

   /* flush changes to the LP solver */
   SCIP_CALL( SCIPlpFlush(lp, blkmem, set) );
   assert(lp->flushed);

   /* set LPI state in the LP solver */
   SCIP_CALL( SCIPlpiSetState(lp->lpi, blkmem, lpistate) );
   lp->primalfeasible = TRUE;
   lp->dualfeasible = TRUE;
   lp->solisbasic = SCIPlpiHasStateBasis(lp->lpi, lpistate);

   return SCIP_OKAY;
}

/** frees LP state information */
SCIP_RETCODE SCIPlpFreeState(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
   assert(lp != NULL);

   SCIP_CALL( SCIPlpiFreeState(lp->lpi, blkmem, lpistate) );

   return SCIP_OKAY;
}

/** sets the upper objective limit of the LP solver */
SCIP_RETCODE SCIPlpSetCutoffbound(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             cutoffbound         /**< new upper objective limit */
   )
{
   assert(lp != NULL);

   SCIPdebugMessage("setting LP upper objective limit from %g to %g\n", lp->cutoffbound, cutoffbound);
   
   /* if the objective function was changed in diving, the cutoff bound has no meaning (it will be set correctly
    * in SCIPendDive())
    */
   if( SCIPlpDivingObjChanged(lp) )
   {
      assert(SCIPsetIsInfinity(set, lp->cutoffbound));
      return SCIP_OKAY;
   }

   /* if the cutoff bound is increased, and the LP was proved to exceed the old cutoff, it is no longer solved;
    * if the cutoff bound is decreased below the current optimal value, the LP now exceeds the objective limit
    */
   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT && cutoffbound > lp->cutoffbound )
   {
      /* mark the current solution invalid */
      lp->solved = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }
   else if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && SCIPlpGetObjval(lp, set) >= cutoffbound )
   {
      assert(lp->flushed);
      assert(lp->solved);
      lp->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
   }

   lp->cutoffbound = cutoffbound;

   return SCIP_OKAY;
}

/** returns the name of the given LP algorithm */
static
const char* lpalgoName(
   SCIP_LPALGO           lpalgo              /**< LP algorithm */
)
{
   switch( lpalgo )
   {
   case SCIP_LPALGO_PRIMALSIMPLEX:
      return "primal simplex";
   case SCIP_LPALGO_DUALSIMPLEX:
      return "dual simplex";
   case SCIP_LPALGO_BARRIER:
      return "barrier";
   case SCIP_LPALGO_BARRIERCROSSOVER:
      return "barrier/crossover";
   default:
      SCIPerrorMessage("invalid LP algorithm\n");
      SCIPABORT();
      return "invalid";
   }
}

/** calls LPI to perform primal simplex, measures time and counts iterations, gets basis feasibility status */
static
SCIP_RETCODE lpPrimalSimplex(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             resolve,            /**< is this a resolving call (starting with feasible basis)? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   SCIP_RETCODE retcode;
   int iterations;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);

   SCIPdebugMessage("solving primal LP %d (LP %d, %d cols, %d rows)\n", 
      stat->nprimallps+1, stat->nlps+1, lp->ncols, lp->nrows);

   *lperror = FALSE;

#if 0
   if( stat->nnodes == 1 && !lp->diving && !lp->probing )
   {
      char fname[SCIP_MAXSTRLEN];
      sprintf(fname, "lp%lld_%d.mps", stat->nnodes, stat->lpcount);
      SCIP_CALL( SCIPlpWrite(lp, fname) );
      SCIPmessagePrintInfo("wrote LP to file <%s> (primal simplex, uobjlim=%g, feastol=%g/%g, fromscratch=%d, fastmip=%d, scaling=%d, presolving=%d)\n", 
         fname, lp->lpiuobjlim, lp->lpifeastol, lp->lpidualfeastol,
         lp->lpifromscratch, lp->lpifastmip, lp->lpiscaling, lp->lpipresolving);
   }
#endif

   /* start timing */
   if( lp->diving || lp->probing )
      SCIPclockStart(stat->divinglptime, set);
   else
      SCIPclockStart(stat->primallptime, set);

   /* call primal simplex */
   retcode = SCIPlpiSolvePrimal(lp->lpi);
   if( retcode == SCIP_LPERROR )
   {
      *lperror = TRUE;
      SCIPdebugMessage("(node %lld) primal simplex solving error in LP %d\n", stat->nnodes, stat->nlps);
   }
   else
   {
      SCIP_CALL( retcode );
   }
   lp->lastlpalgo = SCIP_LPALGO_PRIMALSIMPLEX;
   lp->solisbasic = TRUE;

   /* stop timing */
   if( lp->diving || lp->probing )
      SCIPclockStop(stat->divinglptime, set);
   else
      SCIPclockStop(stat->primallptime, set);

   /* count number of iterations */
   stat->lpcount++;
   SCIP_CALL( SCIPlpGetIterations(lp, &iterations) );
   if( iterations > 0 ) /* don't count the resolves after removing unused columns/rows */
   {
      stat->nlps++;
      stat->nlpiterations += iterations;
      if( resolve && !lp->lpifromscratch && stat->nlps > 1 )
      {
         stat->nprimalresolvelps++;
         stat->nprimalresolvelpiterations += iterations;
      }
      if( lp->diving || lp->probing )
      {
         stat->lastdivenode = stat->nnodes;
         stat->ndivinglps++;
         stat->ndivinglpiterations += iterations;
      }
      else
      {
         stat->nprimallps++;
         stat->nprimallpiterations += iterations;
      }
   }
   else if( keepsol && !(*lperror) )
   {
      /* the solution didn't change: if the solution was valid before resolve, it is still valid */
      if( lp->validsollp == stat->lpcount-1 )
         lp->validsollp = stat->lpcount;
      if( lp->validfarkaslp == stat->lpcount-1 )
         lp->validfarkaslp = stat->lpcount;
   }

   SCIPdebugMessage("solved primal LP %d in %d iterations\n", stat->lpcount, iterations);

   return SCIP_OKAY;
}

/** calls LPI to perform dual simplex, measures time and counts iterations */
static
SCIP_RETCODE lpDualSimplex(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             resolve,            /**< is this a resolving call (starting with feasible basis)? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   SCIP_RETCODE retcode;
   int iterations;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);

   SCIPdebugMessage("solving dual LP %d (LP %d, %d cols, %d rows)\n", 
      stat->nduallps+1, stat->nlps+1, lp->ncols, lp->nrows);

   *lperror = FALSE;

#if 0
   if( stat->nnodes == 1 && !lp->diving && !lp->probing )
   {
      char fname[SCIP_MAXSTRLEN];
      sprintf(fname, "lp%lld_%d.mps", stat->nnodes, stat->lpcount);
      SCIP_CALL( SCIPlpWrite(lp, fname) );
      SCIPmessagePrintInfo("wrote LP to file <%s> (dual simplex, uobjlim=%g, feastol=%g/%g, fromscratch=%d, fastmip=%d, scaling=%d, presolving=%d)\n", 
         fname, lp->lpiuobjlim, lp->lpifeastol, lp->lpidualfeastol, 
         lp->lpifromscratch, lp->lpifastmip, lp->lpiscaling, lp->lpipresolving);
   }
#endif

   /* start timing */
   if( lp->diving || lp->probing )
      SCIPclockStart(stat->divinglptime, set);
   else
      SCIPclockStart(stat->duallptime, set);

   /* call dual simplex */
   retcode = SCIPlpiSolveDual(lp->lpi);
   if( retcode == SCIP_LPERROR )
   {
      *lperror = TRUE;
      SCIPdebugMessage("(node %lld) dual simplex solving error in LP %d\n", stat->nnodes, stat->nlps);
   }
   else
   {
      SCIP_CALL( retcode );
   }
   lp->lastlpalgo = SCIP_LPALGO_DUALSIMPLEX;
   lp->solisbasic = TRUE;

   /* stop timing */
   if( lp->diving || lp->probing )
      SCIPclockStop(stat->divinglptime, set);
   else
      SCIPclockStop(stat->duallptime, set);

   /* count number of iterations */
   stat->lpcount++;
   SCIP_CALL( SCIPlpGetIterations(lp, &iterations) );
   if( iterations > 0 ) /* don't count the resolves after removing unused columns/rows */
   {
      stat->nlps++;
      stat->nlpiterations += iterations;
      if( resolve && !lp->lpifromscratch && stat->nlps > 1  )
      {
         stat->ndualresolvelps++;
         stat->ndualresolvelpiterations += iterations;
      }
      if( lp->diving || lp->probing )
      {
         stat->lastdivenode = stat->nnodes;
         stat->ndivinglps++;
         stat->ndivinglpiterations += iterations;
      }
      else
      {
         stat->nduallps++;
         stat->nduallpiterations += iterations;
      }
   }
   else if( keepsol && !(*lperror) )
   {
      /* the solution didn't change: if the solution was valid before resolve, it is still valid */
      if( lp->validsollp == stat->lpcount-1 )
         lp->validsollp = stat->lpcount;
      if( lp->validfarkaslp == stat->lpcount-1 )
         lp->validfarkaslp = stat->lpcount;
   }

   SCIPdebugMessage("solved dual LP %d in %d iterations\n", stat->lpcount, iterations);

   return SCIP_OKAY;
}

/** calls LPI to perform barrier, measures time and counts iterations, gets basis feasibility status */
static
SCIP_RETCODE lpBarrier(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             crossover,          /**< should crossover be performed? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   SCIP_RETCODE retcode;
   int iterations;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);

   SCIPdebugMessage("solving barrier%s LP %d (LP %d, %d cols, %d rows)\n", 
      crossover ? "/crossover" : "", stat->nbarrierlps+1, stat->nlps+1, lp->ncols, lp->nrows);

   *lperror = FALSE;

#if 0
   if( stat->nnodes == 1 && !lp->diving && !lp->probing )
   {
      char fname[SCIP_MAXSTRLEN];
      sprintf(fname, "lp%lld_%d.mps", stat->nnodes, stat->lpcount);
      SCIP_CALL( SCIPlpWrite(lp, fname) );
      SCIPmessagePrintInfo("wrote LP to file <%s> (barrier, uobjlim=%g, feastol=%g/%g, convtol=%g, fromscratch=%d, fastmip=%d, scaling=%d, presolving=%d)\n", 
         fname, lp->lpiuobjlim, lp->lpifeastol, lp->lpidualfeastol, lp->lpibarrierconvtol,
         lp->lpifromscratch, lp->lpifastmip, lp->lpiscaling, lp->lpipresolving);
   }
#endif

   /* start timing */
   if( lp->diving || lp->probing )
      SCIPclockStart(stat->divinglptime, set);
   else
      SCIPclockStart(stat->barrierlptime, set);

   /* call barrier algorithm */
   retcode = SCIPlpiSolveBarrier(lp->lpi, crossover);
   if( retcode == SCIP_LPERROR )
   {
      *lperror = TRUE;
      SCIPdebugMessage("(node %lld) barrier solving error in LP %d\n", stat->nnodes, stat->nlps);
   }
   else
   {
      SCIP_CALL( retcode );
   }
   lp->lastlpalgo = (crossover ? SCIP_LPALGO_BARRIERCROSSOVER : SCIP_LPALGO_BARRIER);
   lp->solisbasic = crossover;

   /* stop timing */
   if( lp->diving || lp->probing )
      SCIPclockStop(stat->divinglptime, set);
   else
      SCIPclockStop(stat->barrierlptime, set);

   /* count number of iterations */
   stat->lpcount++;
   SCIP_CALL( SCIPlpGetIterations(lp, &iterations) );
   if( iterations > 0 ) /* don't count the resolves after removing unused columns/rows */
   {
      stat->nlps++;
      stat->nlpiterations += iterations;
      if( lp->diving || lp->probing )
      {
         stat->lastdivenode = stat->nnodes;
         stat->ndivinglps++;
         stat->ndivinglpiterations += iterations;
      }
      else
      {
         stat->nbarrierlps++;
         stat->nbarrierlpiterations += iterations;
      }
   }
   else if( keepsol && !(*lperror) )
   {
      /* the solution didn't change: if the solution was valid before resolve, it is still valid */
      if( lp->validsollp == stat->lpcount-1 )
         lp->validsollp = stat->lpcount;
      if( lp->validfarkaslp == stat->lpcount-1 )
         lp->validfarkaslp = stat->lpcount;
   }

   SCIPdebugMessage("solved barrier LP %d in %d iterations\n", stat->lpcount, iterations);

   return SCIP_OKAY;
}

/** solves the LP with the given algorithm */
static
SCIP_RETCODE lpAlgorithm(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPALGO           lpalgo,             /**< LP algorithm that should be applied */
   SCIP_Bool             resolve,            /**< is this a resolving call (starting with feasible basis)? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lperror != NULL);

   SCIPdebugMessage("calling LP algorithm <%s>\n", lpalgoName(lpalgo));

   /* call appropriate LP algorithm */
   switch( lpalgo )
   {
   case SCIP_LPALGO_PRIMALSIMPLEX:
      SCIP_CALL( lpPrimalSimplex(lp, set, stat, resolve, keepsol, lperror) );
      break;

   case SCIP_LPALGO_DUALSIMPLEX:
      SCIP_CALL( lpDualSimplex(lp, set, stat, resolve, keepsol, lperror) );
      break;

   case SCIP_LPALGO_BARRIER:
      SCIP_CALL( lpBarrier(lp, set, stat, FALSE, keepsol, lperror) );
      break;

   case SCIP_LPALGO_BARRIERCROSSOVER:
      SCIP_CALL( lpBarrier(lp, set, stat, TRUE, keepsol, lperror) );
      break;

   default:
      SCIPerrorMessage("invalid LP algorithm\n");
      return SCIP_INVALIDDATA;
   }

   if( !(*lperror) )
   {
      /* check for primal and dual feasibility */
      SCIP_CALL( SCIPlpiGetSolFeasibility(lp->lpi, &lp->primalfeasible, &lp->dualfeasible) );
      
      SCIPdebugMessage("LP feasibility: primalfeasible=%d, dualfeasible=%d\n", lp->primalfeasible, lp->dualfeasible);
   }

   return SCIP_OKAY;
}

#define FEASTOLTIGHTFAC 0.001
/** solves the LP with the given LP algorithm, and tries to resolve numerical problems */
static
SCIP_RETCODE lpSolveStable(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPALGO           lpalgo,             /**< LP algorithm that should be applied */
   SCIP_Bool             resolve,            /**< is this a resolving call (starting with feasible basis)? */
   SCIP_Bool             fastmip,            /**< should the FASTMIP setting of the LP solver be activated? */
   SCIP_Bool             tightfeastol,       /**< should a tighter feasibility tolerance be used? */
   SCIP_Bool             fromscratch,        /**< should the LP be solved from scratch without using current basis? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   SCIP_Bool success;
   SCIP_Bool success2;
   SCIP_Bool simplex;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->looseobjvalinf == 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);

   *lperror = FALSE;

   /* check, whether we solve with a simplex algorithm */
   simplex = (lpalgo == SCIP_LPALGO_PRIMALSIMPLEX || lpalgo == SCIP_LPALGO_DUALSIMPLEX);

   /* solve with given settings (usually fast but unprecise) */
   SCIP_CALL( lpSetUobjlim(lp, set, lp->cutoffbound - lp->looseobjval) );
   SCIP_CALL( lpSetFeastol(lp, tightfeastol ? FEASTOLTIGHTFAC * SCIPsetFeastol(set) : SCIPsetFeastol(set), &success) );
   SCIP_CALL( lpSetDualfeastol(lp, tightfeastol ? FEASTOLTIGHTFAC * SCIPsetDualfeastol(set) : SCIPsetDualfeastol(set), 
         &success) );
   SCIP_CALL( lpSetBarrierconvtol(lp, tightfeastol ? FEASTOLTIGHTFAC * SCIPsetBarrierconvtol(set)
         : SCIPsetBarrierconvtol(set), &success) );
   SCIP_CALL( lpSetFromscratch(lp, fromscratch, &success) );
   SCIP_CALL( lpSetFastmip(lp, fastmip, &success) );
   SCIP_CALL( lpSetScaling(lp, set->lp_scaling, &success) );
   SCIP_CALL( lpSetPresolving(lp, set->lp_presolving, &success) );
   SCIP_CALL( lpSetPricingChar(lp, set->lp_pricing) );
   SCIP_CALL( lpSetLPInfo(lp, set->disp_lpinfo) );
   SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, lperror) );
   resolve = FALSE; /* only the first solve should be counted as resolving call */

   /* check for stability */
   if( !(*lperror) && SCIPlpiIsStable(lp->lpi) )
      return SCIP_OKAY;
   else if( !set->lp_checkstability )
   {
      SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
      if( success )
      {
         SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %lld) numerical troubles in LP %d -- ignoring instability of %s\n",
            stat->nnodes, stat->nlps, lpalgoName(lpalgo));
         return SCIP_OKAY;
      }
   }

   /* if FASTMIP is turned on, solve again without FASTMIP */
   if( fastmip && simplex )
   {
      SCIP_CALL( lpSetFastmip(lp, FALSE, &success) );
      if( success )
      {
         SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %lld) numerical troubles in LP %d -- solve again without FASTMIP with %s\n", 
            stat->nnodes, stat->nlps, lpalgoName(lpalgo));
         SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, lperror) );
         
         /* check for stability */
         if( !(*lperror) && SCIPlpiIsStable(lp->lpi) )
            return SCIP_OKAY;
         else if( !set->lp_checkstability )
         {
            SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
            if( success )
            {
               SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) numerical troubles in LP %d -- ignoring instability of %s\n",
                  stat->nnodes, stat->nlps, lpalgoName(lpalgo));
               return SCIP_OKAY;
            }
         }
      }
   }

   /* solve again with opposite scaling setting */
   SCIP_CALL( lpSetScaling(lp, !set->lp_scaling, &success) );
   if( success )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles in LP %d -- solve again with %s %s scaling\n", 
         stat->nnodes, stat->nlps, lpalgoName(lpalgo), !set->lp_scaling ? "with" : "without");
      SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, lperror) );
   
      /* check for stability */
      if( !(*lperror) && SCIPlpiIsStable(lp->lpi) )
         return SCIP_OKAY;
      else if( !set->lp_checkstability )
      {
         SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
         if( success )
         {
            SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "(node %lld) numerical troubles in LP %d -- ignoring instability of %s\n",
               stat->nnodes, stat->nlps, lpalgoName(lpalgo));
            return SCIP_OKAY;
         }
      }

      /* reset scaling */
      SCIP_CALL( lpSetScaling(lp, set->lp_scaling, &success) );
      assert(success);
   }
      
   /* solve again with opposite presolving setting */
   SCIP_CALL( lpSetPresolving(lp, !set->lp_presolving, &success) );
   if( success )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles in LP %d -- solve again with %s %s presolving\n", 
         stat->nnodes, stat->nlps, lpalgoName(lpalgo), !set->lp_presolving ? "with" : "without");
      SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, lperror) );
   
      /* check for stability */
      if( !(*lperror) && SCIPlpiIsStable(lp->lpi) )
         return SCIP_OKAY;
      else if( !set->lp_checkstability )
      {
         SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
         if( success )
         {
            SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "(node %lld) numerical troubles in LP %d -- ignoring instability of %s\n",
               stat->nnodes, stat->nlps, lpalgoName(lpalgo));
            return SCIP_OKAY;
         }
      }

      /* reset presolving */
      SCIP_CALL( lpSetPresolving(lp, set->lp_presolving, &success) );
      assert(success);
   }
      
   /* solve again with a tighter feasibility tolerance */
   if( !tightfeastol )
   {
      SCIP_CALL( lpSetFeastol(lp, FEASTOLTIGHTFAC * SCIPsetFeastol(set), &success) );
      SCIP_CALL( lpSetDualfeastol(lp, FEASTOLTIGHTFAC * SCIPsetDualfeastol(set), &success2) );
      SCIP_CALL( lpSetBarrierconvtol(lp, FEASTOLTIGHTFAC * SCIPsetBarrierconvtol(set), &success2) );
      if( success || success2 )
      {
         SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %lld) numerical troubles in LP %d -- solve again with %s with tighter feasibility tolerance\n", 
            stat->nnodes, stat->nlps, lpalgoName(lpalgo));
         SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, lperror) );
      
         /* check for stability */
         if( !(*lperror) && SCIPlpiIsStable(lp->lpi) )
            return SCIP_OKAY;
         else if( !set->lp_checkstability )
         {
            SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
            if( success )
            {
               SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) numerical troubles in LP %d -- ignoring instability of %s\n",
                  stat->nnodes, stat->nlps, lpalgoName(lpalgo));
               return SCIP_OKAY;
            }
         }

         /* reset feasibility tolerance */
         SCIP_CALL( lpSetFeastol(lp, SCIPsetFeastol(set), &success) );
         SCIP_CALL( lpSetDualfeastol(lp, SCIPsetDualfeastol(set), &success) );
         SCIP_CALL( lpSetBarrierconvtol(lp, SCIPsetBarrierconvtol(set), &success) );
      }
   }

   /* if not already done, solve again from scratch */
   if( !fromscratch && simplex )
   {
      SCIP_CALL( lpSetFromscratch(lp, TRUE, &success) );
      if( success )
      {
         SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s\n", 
            stat->nnodes, stat->nlps, lpalgoName(lpalgo));
         SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, lperror) );
         
         /* check for stability */
         if( !(*lperror) && SCIPlpiIsStable(lp->lpi) )
            return SCIP_OKAY;
         else if( !set->lp_checkstability )
         {
            SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
            if( success )
            {
               SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) numerical troubles in LP %d -- ignoring instability of %s\n",
                  stat->nnodes, stat->nlps, lpalgoName(lpalgo));
               return SCIP_OKAY;
            }
         }
      }
   }

   /* solve again, use other simplex this time */
   if( simplex )
   {
      lpalgo = (lpalgo == SCIP_LPALGO_PRIMALSIMPLEX ? SCIP_LPALGO_DUALSIMPLEX : SCIP_LPALGO_PRIMALSIMPLEX);
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s\n", 
         stat->nnodes, stat->nlps, lpalgoName(lpalgo));
      SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, lperror) );

      /* check for stability */
      if( !(*lperror) && SCIPlpiIsStable(lp->lpi) )
         return SCIP_OKAY;
      else if( !set->lp_checkstability )
      {
         SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
         if( success )
         {
            SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "(node %lld) numerical troubles in LP %d -- ignoring instability of %s\n",
               stat->nnodes, stat->nlps, lpalgoName(lpalgo));
            return SCIP_OKAY;
         }
      }

      /* solve again with opposite scaling and other simplex */
      SCIP_CALL( lpSetScaling(lp, !set->lp_scaling, &success) );
      if( success )
      {
         SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s %s scaling\n", 
            stat->nnodes, stat->nlps, lpalgoName(lpalgo), !set->lp_scaling ? "with" : "without");
         SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, lperror) );
         
         /* check for stability */
         if( !(*lperror) && SCIPlpiIsStable(lp->lpi) )
            return SCIP_OKAY;
         else if( !set->lp_checkstability )
         {
            SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
            if( success )
            {
               SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) numerical troubles in LP %d -- ignoring instability of %s\n",
                  stat->nnodes, stat->nlps, lpalgoName(lpalgo));
               return SCIP_OKAY;
            }
         }
         
         /* reset scaling */
         SCIP_CALL( lpSetScaling(lp, set->lp_scaling, &success) );
         assert(success);
      }

      /* solve again with opposite presolving and other simplex */
      SCIP_CALL( lpSetPresolving(lp, !set->lp_presolving, &success) );
      if( success )
      {
         SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s %s presolving\n", 
            stat->nnodes, stat->nlps, lpalgoName(lpalgo), !set->lp_presolving ? "with" : "without");
         SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, lperror) );
         
         /* check for stability */
         if( !(*lperror) && SCIPlpiIsStable(lp->lpi) )
            return SCIP_OKAY;
         else if( !set->lp_checkstability )
         {
            SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
            if( success )
            {
               SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) numerical troubles in LP %d -- ignoring instability of %s\n",
                  stat->nnodes, stat->nlps, lpalgoName(lpalgo));
               return SCIP_OKAY;
            }
         }
         
         /* reset presolving */
         SCIP_CALL( lpSetPresolving(lp, set->lp_presolving, &success) );
         assert(success);
      }
      
      /* solve again with tighter feasibility tolerance, use other simplex this time */
      if( !tightfeastol )
      {
         SCIP_CALL( lpSetFeastol(lp, FEASTOLTIGHTFAC * SCIPsetFeastol(set), &success) );
         SCIP_CALL( lpSetDualfeastol(lp, FEASTOLTIGHTFAC * SCIPsetDualfeastol(set), &success2) );
         SCIP_CALL( lpSetBarrierconvtol(lp, FEASTOLTIGHTFAC * SCIPsetBarrierconvtol(set), &success2) );
         if( success || success2 )
         {
            SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s with tighter feasibility tolerance\n", 
               stat->nnodes, stat->nlps, lpalgoName(lpalgo));
            SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, lperror) );
         
            /* check for stability */
            if( !(*lperror) && SCIPlpiIsStable(lp->lpi) )
               return SCIP_OKAY;
            else if( !set->lp_checkstability )
            {
               SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
               if( success )
               {
                  SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                     "(node %lld) numerical troubles in LP %d -- ignoring instability of %s\n",
                     stat->nnodes, stat->nlps, lpalgoName(lpalgo));
                  return SCIP_OKAY;
               }
            }
         
            /* reset feasibility tolerance */
            SCIP_CALL( lpSetFeastol(lp, SCIPsetFeastol(set), &success) );
            SCIP_CALL( lpSetDualfeastol(lp, SCIPsetDualfeastol(set), &success) );
            SCIP_CALL( lpSetBarrierconvtol(lp, SCIPsetBarrierconvtol(set), &success) );
         }
      }
   }

   /* nothing worked -- store the instable LP to a file and exit with an LPERROR */
   SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "(node %lld) unresolved numerical troubles in LP %d\n", 
      stat->nnodes, stat->nlps);
   *lperror = TRUE;

   return SCIP_OKAY;
}

/** solves the LP with the given algorithm and evaluates return status */
static
SCIP_RETCODE lpSolve(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPALGO           lpalgo,             /**< LP algorithm that should be applied */
   SCIP_Bool             needprimalray,      /**< if the LP is unbounded, do we need a primal ray? */
   SCIP_Bool             needdualray,        /**< if the LP is infeasible, do we need a dual ray? */
   SCIP_Bool             resolve,            /**< is this a resolving call (starting with feasible basis)? */
   SCIP_Bool             fastmip,            /**< should the FASTMIP setting of the LP solver be activated? */
   SCIP_Bool             tightfeastol,       /**< should a tighter feasibility tolerance be used? */
   SCIP_Bool             fromscratch,        /**< should the LP be solved from scratch without using current basis? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   SCIP_Bool solvedprimal;
   SCIP_Bool solveddual;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);

   checkLinks(lp);

   solvedprimal = FALSE;
   solveddual = FALSE;

 SOLVEAGAIN:
   /* call simplex */
   SCIP_CALL( lpSolveStable(lp, set, stat, lpalgo, resolve, fastmip, tightfeastol, fromscratch, keepsol, lperror) );
   resolve = FALSE; /* only the first solve should be counted as resolving call */
   solvedprimal = solvedprimal || (lpalgo == SCIP_LPALGO_PRIMALSIMPLEX);
   solveddual = solveddual || (lpalgo == SCIP_LPALGO_DUALSIMPLEX);
      
   /* check, if an error occured */
   if( *lperror )
   {
      SCIPdebugMessage("unresolved error while solving LP with %s\n", lpalgoName(lp->lastlpalgo));
      lp->solved = FALSE;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
      return SCIP_OKAY;
   }

   /* evaluate solution status */
   if( SCIPlpiIsOptimal(lp->lpi) )
   {
      assert(lp->primalfeasible);
      assert(lp->dualfeasible);
      lp->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
      SCIP_CALL( SCIPlpiGetObjval(lp->lpi, &lp->lpobjval) );
      if( SCIPsetIsRelGE(set, lp->lpobjval, lp->lpiuobjlim) )
      {
         /* the solver may return the optimal value, even if this is greater or equal than the upper bound */
         SCIPdebugMessage("optimal solution %g exceeds objective limit %g\n", lp->lpobjval, lp->lpiuobjlim);
         lp->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
         lp->lpobjval = SCIPsetInfinity(set);
      }
   }
   else if( SCIPlpiIsObjlimExc(lp->lpi) )
   {
#if 0 /* SOPLEX may return with objective limit reached in any case, because it doesn't distinct btw. primal and dual */
      if( lp->lastlpalgo == SCIP_LPALGO_PRIMALSIMPLEX )
      {
         SCIPerrorMessage("objective limit exceeded in primal simplex - this should not happen, because no lower limit exists\n");
         lp->lpsolstat = SCIP_LPSOLSTAT_ERROR;
         lp->lpobjval = -SCIPsetInfinity(set);
         return SCIP_LPERROR;
      }
#endif
      lp->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
      lp->lpobjval = SCIPsetInfinity(set);
   }
   else if( SCIPlpiIsPrimalInfeasible(lp->lpi) )
   {
      if( needdualray && !SCIPlpiHasDualRay(lp->lpi) && !solveddual )
      {
         assert(lpalgo != SCIP_LPALGO_DUALSIMPLEX);
         lpalgo = SCIP_LPALGO_DUALSIMPLEX;
         goto SOLVEAGAIN;
      }
      lp->lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;
      lp->lpobjval = SCIPsetInfinity(set);
   }
   else if( SCIPlpiExistsPrimalRay(lp->lpi) )
   {
      if( needprimalray && !SCIPlpiHasPrimalRay(lp->lpi) && !solvedprimal )
      {
         assert(lpalgo != SCIP_LPALGO_PRIMALSIMPLEX);
         lpalgo = SCIP_LPALGO_PRIMALSIMPLEX;
         goto SOLVEAGAIN;
      }
      lp->lpsolstat = SCIP_LPSOLSTAT_UNBOUNDEDRAY;
      lp->lpobjval = -SCIPsetInfinity(set);
   }
   else if( SCIPlpiIsIterlimExc(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_ITERLIMIT;
      lp->lpobjval = -SCIPsetInfinity(set);
   }
   else if( SCIPlpiIsTimelimExc(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
      lp->lpobjval = -SCIPsetInfinity(set);
   }
   else if( !solveddual )
   {
      assert(lpalgo != SCIP_LPALGO_DUALSIMPLEX);
      lpalgo = SCIP_LPALGO_DUALSIMPLEX;
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) solution status of LP %d couldn't be proved (internal status:%d) -- solve again with %s\n", 
         stat->nnodes, stat->nlps, lpalgoName(lpalgo));
      goto SOLVEAGAIN;
   }
   else if( !solvedprimal )
   {
      assert(lpalgo != SCIP_LPALGO_PRIMALSIMPLEX);
      lpalgo = SCIP_LPALGO_PRIMALSIMPLEX;
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) solution status of LP %d couldn't be proved (internal status:%d) -- solve again with %s\n", 
         stat->nnodes, stat->nlps, lpalgoName(lpalgo));
      goto SOLVEAGAIN;
   }
   else
   {
      SCIPerrorMessage("(node %lld) error or unknown return status of %s in LP %d (internal status: %d)\n", 
         stat->nnodes, lpalgoName(lp->lastlpalgo), stat->nlps, SCIPlpiGetInternalStatus(lp->lpi));
      lp->lpsolstat = SCIP_LPSOLSTAT_ERROR;
      return SCIP_LPERROR;
   }

   lp->solved = TRUE;

   SCIPdebugMessage("solving LP with %s returned solstat=%d (internal status: %d, pfeas=%d, dfeas=%d)\n",
      lpalgoName(lp->lastlpalgo), lp->lpsolstat, SCIPlpiGetInternalStatus(lp->lpi),
      SCIPlpiIsPrimalFeasible(lp->lpi), SCIPlpiIsDualFeasible(lp->lpi));

   return SCIP_OKAY;
}

/** flushes the LP and solves it with the primal or dual simplex algorithm, depending on the current basis feasibility */
static
SCIP_RETCODE lpFlushAndSolve(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             needprimalray,      /**< if the LP is unbounded, do we need a primal ray? */
   SCIP_Bool             needdualray,        /**< if the LP is infeasible, do we need a dual ray? */
   SCIP_Bool             fastmip,            /**< should the FASTMIP setting of the LP solver be activated? */
   SCIP_Bool             tightfeastol,       /**< should a tighter feasibility tolerance be used? */
   SCIP_Bool             fromscratch,        /**< should the LP be solved from scratch without using current basis? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   SCIP_Bool resolve;
   char algo;

   assert(lp != NULL);
   assert(set != NULL);
   assert(lperror != NULL);

   /* flush changes to the LP solver */
   SCIP_CALL( SCIPlpFlush(lp, blkmem, set) );
   fastmip = fastmip && !lp->flushaddedcols && !lp->flushdeletedcols; /* turn off FASTMIP if columns were changed */
   
   /* select LP algorithm to apply */
   resolve = lp->solisbasic && (lp->dualfeasible || lp->primalfeasible);
   algo = resolve ? set->lp_resolvealgorithm : set->lp_initalgorithm;
   switch( algo )
   {
   case 's':
      /* select simplex method */
      if( lp->dualfeasible || !lp->primalfeasible )
      {
         SCIPdebugMessage("solving dual LP\n");
         SCIP_CALL( lpSolve(lp, set, stat, SCIP_LPALGO_DUALSIMPLEX, needprimalray, needdualray, resolve, fastmip,
               tightfeastol, fromscratch, keepsol, lperror) );
      }
      else
      {
         SCIPdebugMessage("solving primal LP\n");
         SCIP_CALL( lpSolve(lp, set, stat, SCIP_LPALGO_PRIMALSIMPLEX, needprimalray, needdualray, resolve, fastmip,
               tightfeastol, fromscratch, keepsol, lperror) );
      }
      break;

   case 'b':
      SCIPdebugMessage("solving barrier LP\n");
      SCIP_CALL( lpSolve(lp, set, stat, SCIP_LPALGO_BARRIER, needprimalray, needdualray, resolve, fastmip,
            tightfeastol, fromscratch, keepsol, lperror) );
      break;

   case 'c':
      SCIPdebugMessage("solving barrier LP with crossover\n");
      SCIP_CALL( lpSolve(lp, set, stat, SCIP_LPALGO_BARRIERCROSSOVER, needprimalray, needdualray, resolve, fastmip,
            tightfeastol, fromscratch, keepsol, lperror) );
      break;

   default:
      SCIPerrorMessage("invalid parameter setting <%c> for LP algorithm\n", algo);
      return SCIP_PARAMETERWRONGVAL;
   }
   assert(!(*lperror) || !lp->solved);

   return SCIP_OKAY;
}

/** solves the LP with simplex algorithm, and copy the solution into the column's data */
SCIP_RETCODE SCIPlpSolveAndEval(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   int                   itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool             aging,              /**< should aging and removal of obsolete cols/rows be applied? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   SCIP_Bool needprimalray;
   SCIP_Bool needdualray;

   assert(lp != NULL);
   assert(prob != NULL);
   assert(prob->nvars >= lp->ncols);
   assert(lperror != NULL);

   SCIPdebugMessage("solving LP: %d rows, %d cols, primalfeasible=%d, dualfeasible=%d, solved=%d, diving=%d, probing=%d, cutoff=%g\n", 
      lp->nrows, lp->ncols, lp->primalfeasible, lp->dualfeasible, lp->solved, lp->diving, lp->probing, lp->cutoffbound);

   *lperror = FALSE;

   /* check whether we need a proof of unboundness or infeasibility by a primal or dual ray */
   needprimalray = TRUE;
   needdualray = (!SCIPprobAllColsInLP(prob, set, lp) || set->misc_exactsolve || set->conf_uselp);

   /* flush changes to the LP solver */
   SCIP_CALL( SCIPlpFlush(lp, blkmem, set) );

   /* set iteration limit for this solve */
   SCIP_CALL( lpSetIterationLimit(lp, itlim) );

   if( !lp->solved )
   {
      SCIP_Bool primalfeasible;
      SCIP_Bool dualfeasible;
      SCIP_Bool fastmip;
      SCIP_Bool tightfeastol;
      SCIP_Bool fromscratch;
      int oldnlps;

      /* set initial LP solver settings */
      fastmip = set->lp_fastmip && lp->lpihasfastmip && !lp->flushaddedcols && !lp->flushdeletedcols && stat->nnodes > 1;
      tightfeastol = FALSE;
      fromscratch = FALSE;

   SOLVEAGAIN:
      /* solve the LP */
      oldnlps = stat->nlps;
      SCIP_CALL( lpFlushAndSolve(lp, blkmem, set, stat, needprimalray, needdualray, fastmip, tightfeastol, fromscratch,
            keepsol, lperror) );
      SCIPdebugMessage("lpFlushAndSolve() returned solstat %d (error=%d)\n", SCIPlpGetSolstat(lp), *lperror);
      assert(!(*lperror) || !lp->solved);

      /* check for error */
      if( *lperror )
         return SCIP_OKAY;

      /* evaluate solution status */
      switch( SCIPlpGetSolstat(lp) )
      {
      case SCIP_LPSOLSTAT_OPTIMAL:
         if( set->lp_checkfeas )
         {
            /* get LP solution and check the solution's feasibility again */
            SCIP_CALL( SCIPlpGetSol(lp, set, stat, &primalfeasible, &dualfeasible) );
         }
         else
         {
            /* get LP solution believing in the feasibility of the LP solution */
            SCIP_CALL( SCIPlpGetSol(lp, set, stat, NULL, NULL) );
            primalfeasible = TRUE;
            dualfeasible = TRUE;
         }
         if( primalfeasible && dualfeasible && aging && !lp->diving && !lp->probing && stat->nlps > oldnlps )
         {
            /* update ages and remove obsolete columns and rows from LP */
            SCIP_CALL( SCIPlpUpdateAges(lp, stat) );
            SCIP_CALL( SCIPlpRemoveNewObsoletes(lp, blkmem, set, stat) );
            
            if( !lp->solved )
            {
               /* resolve LP after removing obsolete columns and rows */
               SCIPdebugMessage("removed obsoletes - resolve LP again: %d rows, %d cols\n", lp->nrows, lp->ncols);
               aging = FALSE; /* to prevent infinite loops */
               goto SOLVEAGAIN;
            }
         }
         if( !primalfeasible || !dualfeasible )
         {
            SCIP_Bool simplex = (lp->lastlpalgo == SCIP_LPALGO_PRIMALSIMPLEX || lp->lastlpalgo == SCIP_LPALGO_DUALSIMPLEX);

            if( fastmip && simplex )
            {
               /* solution is infeasible (this can happen due to numerical problems): solve again without FASTMIP */
               SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) solution of LP %d not optimal (pfeas=%d, dfeas=%d) -- solving again without FASTMIP\n",
                  stat->nnodes, stat->nlps, primalfeasible, dualfeasible);
               fastmip = FALSE;
               goto SOLVEAGAIN;
            }
            else if( !tightfeastol )
            {
               /* solution is infeasible (this can happen due to numerical problems): solve again with tighter feasibility
                * tolerance
                */
               SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) solution of LP %d not optimal (pfeas=%d, dfeas=%d) -- solving again with tighter feasibility tolerance\n",
                  stat->nnodes, stat->nlps, primalfeasible, dualfeasible);
               tightfeastol = TRUE;
               goto SOLVEAGAIN;
            }
            else if( !fromscratch && simplex )
            {
               /* solution is infeasible (this can happen due to numerical problems): solve again from scratch */
               SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) solution of LP %d not optimal (pfeas=%d, dfeas=%d) -- solving again from scratch\n",
                  stat->nnodes, stat->nlps, primalfeasible, dualfeasible);
               fromscratch = TRUE;
               goto SOLVEAGAIN;
            }
            else
            {
               SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) unresolved numerical troubles in LP %d\n", stat->nnodes, stat->nlps);
               lp->solved = FALSE;
               lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
               *lperror = TRUE;
            }
         }
         SCIPdebugMessage(" -> LP objective value: %g + %g = %g (solstat=%d, cutoff=%g)\n",
            lp->lpobjval, lp->looseobjval, lp->lpobjval + lp->looseobjval, lp->lpsolstat, lp->cutoffbound);
         break;

      case SCIP_LPSOLSTAT_INFEASIBLE:
	 if( !SCIPprobAllColsInLP(prob, set, lp) || set->misc_exactsolve )
         {
            SCIP_CALL( SCIPlpGetDualfarkas(lp, set, stat) );
         }
         SCIPdebugMessage(" -> LP infeasible\n");
         break;

      case SCIP_LPSOLSTAT_UNBOUNDEDRAY:
         SCIP_CALL( SCIPlpGetUnboundedSol(lp, set, stat) );
         SCIPdebugMessage(" -> LP has unbounded primal ray\n");
         break;

      case SCIP_LPSOLSTAT_OBJLIMIT:
         if( !SCIPprobAllColsInLP(prob, set, lp) || set->misc_exactsolve )
         {
            SCIP_CALL( SCIPlpGetSol(lp, set, stat, NULL, NULL) );
         }
         SCIPdebugMessage(" -> LP objective limit reached\n");
         break;

      case SCIP_LPSOLSTAT_ITERLIMIT:
         break;

      case SCIP_LPSOLSTAT_TIMELIMIT:
         /**@todo time limit exceeded processing */
         SCIPerrorMessage("LP time limit exceeded -- case not implemented yet\n");
         return SCIP_ERROR;

      case SCIP_LPSOLSTAT_ERROR:
      case SCIP_LPSOLSTAT_NOTSOLVED:
         SCIPerrorMessage("error in LP solver\n");
         return SCIP_LPERROR;

      default:
         SCIPerrorMessage("unknown LP solution status\n");
         return SCIP_ERROR;
      }
   }
   assert(!(*lperror) || !lp->solved);

   return SCIP_OKAY;
}

/** gets solution status of current LP */
SCIP_LPSOLSTAT SCIPlpGetSolstat(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->solved || lp->lpsolstat == SCIP_LPSOLSTAT_NOTSOLVED);

   return (lp->flushed ? lp->lpsolstat : SCIP_LPSOLSTAT_NOTSOLVED);
}

/** gets objective value of current LP */
SCIP_Real SCIPlpGetObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert((lp->nloosevars > 0) || (lp->looseobjvalinf == 0 && lp->looseobjval == 0.0));
   assert(set != NULL);

   if( !lp->flushed )
      return SCIP_INVALID;
   else if( SCIPsetIsInfinity(set, lp->lpobjval) )
      return lp->lpobjval;
   else if( lp->looseobjvalinf > 0 )
      return -SCIPsetInfinity(set);
   else
      return lp->lpobjval + lp->looseobjval;
}

/** gets part of objective value of current LP that results from COLUMN variables only */
SCIP_Real SCIPlpGetColumnObjval(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->solved);

   return (lp->flushed ? lp->lpobjval : SCIP_INVALID);
}

/** gets part of objective value of current LP that results from LOOSE variables only */
SCIP_Real SCIPlpGetLooseObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert((lp->nloosevars > 0) || (lp->looseobjvalinf == 0 && lp->looseobjval == 0.0));
   assert(set != NULL);

   if( !lp->flushed )
      return SCIP_INVALID;
   else if( lp->looseobjvalinf > 0 )
      return -SCIPsetInfinity(set);
   else
      return lp->looseobjval;
}

/** gets current pseudo objective value */
SCIP_Real SCIPlpGetPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(lp->pseudoobjvalinf >= 0);
   assert(set != NULL);

   if( lp->pseudoobjvalinf > 0 )
      return -SCIPsetInfinity(set);
   else
      return lp->pseudoobjval;
}

/** gets pseudo objective value, if a bound of the given variable would be modified in the given way */
SCIP_Real SCIPlpGetModifiedPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   SCIP_Real pseudoobjval;
   int pseudoobjvalinf;
   
   pseudoobjval = lp->pseudoobjval;
   pseudoobjvalinf = lp->pseudoobjvalinf;
   if( boundtype == SCIPvarGetBestBoundType(var) )
   {
      if( SCIPsetIsInfinity(set, REALABS(oldbound)) )
         pseudoobjvalinf--;
      else
         pseudoobjval -= oldbound * SCIPvarGetObj(var);
      assert(pseudoobjvalinf >= 0);
      if( SCIPsetIsInfinity(set, REALABS(newbound)) )
         pseudoobjvalinf++;
      else
         pseudoobjval += newbound * SCIPvarGetObj(var);
   }
   assert(pseudoobjvalinf >= 0);

   if( pseudoobjvalinf > 0 )
      return -SCIPsetInfinity(set);
   else
      return pseudoobjval;
}

/** gets pseudo objective value, if a bound of the given variable would be modified in the given way;
 *  perform calculations with interval arithmetic to get an exact lower bound
 */
SCIP_Real SCIPlpGetModifiedProvedPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   SCIP_Real pseudoobjval;
   int pseudoobjvalinf;
   
   pseudoobjval = lp->pseudoobjval;
   pseudoobjvalinf = lp->pseudoobjvalinf;
   if( boundtype == SCIPvarGetBestBoundType(var) )
   {
      SCIP_INTERVAL obj;
      SCIP_INTERVAL bd;
      SCIP_INTERVAL prod;
      SCIP_INTERVAL psval;

      SCIPintervalSet(&psval, pseudoobjval);
      SCIPintervalSet(&obj, SCIPvarGetObj(var));

      if( SCIPsetIsInfinity(set, REALABS(oldbound)) )
         pseudoobjvalinf--;
      else
      {
         SCIPintervalSet(&bd, oldbound);
         SCIPintervalMul(&prod, bd, obj);
         SCIPintervalSub(&psval, psval, prod);
      }
      assert(pseudoobjvalinf >= 0);
      if( SCIPsetIsInfinity(set, REALABS(newbound)) )
         pseudoobjvalinf++;
      else
      {
         SCIPintervalSet(&bd, newbound);
         SCIPintervalMul(&prod, bd, obj);
         SCIPintervalAdd(&psval, psval, prod);
      }

      pseudoobjval = SCIPintervalGetInf(psval);
   }
   assert(pseudoobjvalinf >= 0);

   if( pseudoobjvalinf > 0 )
      return -SCIPsetInfinity(set);
   else
      return pseudoobjval;
}

/** updates current pseudo and loose objective values for a change in a variable's objective value or bounds */
static
SCIP_RETCODE lpUpdateVar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldobj,             /**< old objective value of variable */
   SCIP_Real             oldlb,              /**< old objective value of variable */
   SCIP_Real             oldub,              /**< old objective value of variable */
   SCIP_Real             newobj,             /**< new objective value of variable */
   SCIP_Real             newlb,              /**< new objective value of variable */
   SCIP_Real             newub               /**< new objective value of variable */
   )
{
   SCIP_Real deltaval;
   int deltainf;

   assert(lp != NULL);
   assert(lp->pseudoobjvalinf >= 0);
   assert(lp->looseobjvalinf >= 0);
   assert(!SCIPsetIsInfinity(set, REALABS(oldobj)));
   assert(!SCIPsetIsInfinity(set, oldlb));
   assert(!SCIPsetIsInfinity(set, -oldub));
   assert(!SCIPsetIsInfinity(set, REALABS(newobj)));
   assert(!SCIPsetIsInfinity(set, newlb));
   assert(!SCIPsetIsInfinity(set, -newub));
   assert(var != NULL);

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE && SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      SCIPerrorMessage("LP was informed of an objective change of a non-mutable variable\n");
      return SCIP_INVALIDDATA;
   }

   assert(SCIPvarGetProbindex(var) >= 0);

   deltaval = 0.0;
   deltainf = 0;

   /* subtract old pseudo objective value */
   if( SCIPsetIsPositive(set, oldobj) )
   {
      if( SCIPsetIsInfinity(set, -oldlb) )
         deltainf--;
      else
         deltaval -= oldlb * oldobj;
   }
   else if( SCIPsetIsNegative(set, oldobj) )
   {
      if( SCIPsetIsInfinity(set, oldub) )
         deltainf--;
      else
         deltaval -= oldub * oldobj;
   }

   /* add new pseudo objective value */
   if( SCIPsetIsPositive(set, newobj) )
   {
      if( SCIPsetIsInfinity(set, -newlb) )
         deltainf++;
      else
         deltaval += newlb * newobj;
   }
   else if( SCIPsetIsNegative(set, newobj) )
   {
      if( SCIPsetIsInfinity(set, newub) )
         deltainf++;
      else
         deltaval += newub * newobj;
   }

   /* update the pseudo and loose objective values */
   lp->pseudoobjval += deltaval;
   lp->pseudoobjvalinf += deltainf;
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
   {
      lp->looseobjval += deltaval;
      lp->looseobjvalinf += deltainf;
   }

   assert(lp->pseudoobjvalinf >= 0);
   assert(lp->looseobjvalinf >= 0);

   /* update squared euclidean norm and sum norm of objective function vector */
   lp->objsqrnorm += SQR(newobj) - SQR(oldobj);
   lp->objsqrnorm = MAX(lp->objsqrnorm, 0.0);
   lp->objsumnorm += REALABS(newobj) - REALABS(oldobj);
   lp->objsumnorm = MAX(lp->objsumnorm, 0.0);

   return SCIP_OKAY;
}

/** updates current pseudo and loose objective values for a change in a variable's objective value or bounds;
 *  pseudo objective value is calculated with interval arithmetics to get a proved lower bound
 */
static
SCIP_RETCODE lpUpdateVarProved(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldobj,             /**< old objective value of variable */
   SCIP_Real             oldlb,              /**< old objective value of variable */
   SCIP_Real             oldub,              /**< old objective value of variable */
   SCIP_Real             newobj,             /**< new objective value of variable */
   SCIP_Real             newlb,              /**< new objective value of variable */
   SCIP_Real             newub               /**< new objective value of variable */
   )
{
   SCIP_INTERVAL deltaval;
   SCIP_INTERVAL bd;
   SCIP_INTERVAL obj;
   SCIP_INTERVAL prod;
   SCIP_INTERVAL psval;
   int deltainf;

   assert(lp != NULL);
   assert(lp->pseudoobjvalinf >= 0);
   assert(lp->looseobjvalinf >= 0);
   assert(!SCIPsetIsInfinity(set, REALABS(oldobj)));
   assert(!SCIPsetIsInfinity(set, oldlb));
   assert(!SCIPsetIsInfinity(set, -oldub));
   assert(!SCIPsetIsInfinity(set, REALABS(newobj)));
   assert(!SCIPsetIsInfinity(set, newlb));
   assert(!SCIPsetIsInfinity(set, -newub));
   assert(var != NULL);

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE && SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      SCIPerrorMessage("LP was informed of an objective change of a non-mutable variable\n");
      return SCIP_INVALIDDATA;
   }

   assert(SCIPvarGetProbindex(var) >= 0);

   SCIPintervalSet(&deltaval, 0.0);
   deltainf = 0;

   /* subtract old pseudo objective value */
   if( oldobj > 0.0 )
   {
      if( SCIPsetIsInfinity(set, -oldlb) )
         deltainf--;
      else
      {
         SCIPintervalSet(&bd, oldlb);
         SCIPintervalSet(&obj, oldobj);
         SCIPintervalMul(&prod, bd, obj);
         SCIPintervalSub(&deltaval, deltaval, prod);  /* deltaval -= oldlb * oldobj; */
      }
   }
   else if( oldobj < 0.0 )
   {
      if( SCIPsetIsInfinity(set, oldub) )
         deltainf--;
      else
      {
         SCIPintervalSet(&bd, oldub);
         SCIPintervalSet(&obj, oldobj);
         SCIPintervalMul(&prod, bd, obj);
         SCIPintervalSub(&deltaval, deltaval, prod);  /* deltaval -= oldub * oldobj; */
      }
   }

   /* add new pseudo objective value */
   if( newobj > 0.0 )
   {
      if( SCIPsetIsInfinity(set, -newlb) )
         deltainf++;
      else
      {
         SCIPintervalSet(&bd, newlb);
         SCIPintervalSet(&obj, newobj);
         SCIPintervalMul(&prod, bd, obj);
         SCIPintervalAdd(&deltaval, deltaval, prod);  /* deltaval += newlb * newobj; */
      }
   }
   else if( newobj < 0.0 )
   {
      if( SCIPsetIsInfinity(set, newub) )
         deltainf++;
      else
      {
         SCIPintervalSet(&bd, newub);
         SCIPintervalSet(&obj, newobj);
         SCIPintervalMul(&prod, bd, obj);
         SCIPintervalAdd(&deltaval, deltaval, prod);  /* deltaval += newub * newobj; */
      }
   }

   /* update the pseudo and loose objective values */
   SCIPintervalSet(&psval, lp->pseudoobjval);
   SCIPintervalAdd(&psval, psval, deltaval);
   lp->pseudoobjval = SCIPintervalGetInf(psval);
   lp->pseudoobjvalinf += deltainf;
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
   {
      SCIPintervalSet(&psval, lp->looseobjval);
      SCIPintervalAdd(&psval, psval, deltaval);
      lp->looseobjval = SCIPintervalGetInf(psval);
      lp->looseobjvalinf += deltainf;
   }

   assert(lp->pseudoobjvalinf >= 0);
   assert(lp->looseobjvalinf >= 0);

   return SCIP_OKAY;
}

/** updates current pseudo and loose objective value for a change in a variable's objective value */
SCIP_RETCODE SCIPlpUpdateVarObj(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldobj,             /**< old objective value of variable */
   SCIP_Real             newobj              /**< new objective value of variable */
   )
{
   assert(set != NULL);
   assert(var != NULL);

   if( set->misc_exactsolve )
   {
      if( oldobj != newobj ) /*lint !e777*/
      {
         SCIP_CALL( lpUpdateVarProved(lp, set, var, oldobj, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
                        newobj, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) );
      }
   }
   else
   {
      if( !SCIPsetIsEQ(set, oldobj, newobj) )
      {
         SCIP_CALL( lpUpdateVar(lp, set, var, oldobj, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
                        newobj, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) );
      }
   }
   
   return SCIP_OKAY;
}

/** updates current pseudo and loose objective value for a change in a variable's lower bound */
SCIP_RETCODE SCIPlpUpdateVarLb(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb               /**< new lower bound of variable */
   )
{
   assert(set != NULL);
   assert(var != NULL);

   if( set->misc_exactsolve )
   {
      if( oldlb != newlb && SCIPvarGetObj(var) > 0.0 ) /*lint !e777*/
      {
         SCIP_CALL( lpUpdateVarProved(lp, set, var, SCIPvarGetObj(var), oldlb, SCIPvarGetUbLocal(var), 
                        SCIPvarGetObj(var), newlb, SCIPvarGetUbLocal(var)) );
      }
   }
   else
   {
      if( !SCIPsetIsEQ(set, oldlb, newlb) && SCIPsetIsPositive(set, SCIPvarGetObj(var)) )
      {
         SCIP_CALL( lpUpdateVar(lp, set, var, SCIPvarGetObj(var), oldlb, SCIPvarGetUbLocal(var), 
                        SCIPvarGetObj(var), newlb, SCIPvarGetUbLocal(var)) );
      }
   }
   
   return SCIP_OKAY;
}

/** updates current pseudo objective value for a change in a variable's upper bound */
SCIP_RETCODE SCIPlpUpdateVarUb(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldub,              /**< old upper bound of variable */
   SCIP_Real             newub               /**< new upper bound of variable */
   )
{
   assert(set != NULL);
   assert(var != NULL);

   if( set->misc_exactsolve )
   {
      if( oldub != newub && SCIPvarGetObj(var) < 0.0 ) /*lint !e777*/
      {
         SCIP_CALL( lpUpdateVarProved(lp, set, var, SCIPvarGetObj(var), SCIPvarGetLbLocal(var), oldub, 
                        SCIPvarGetObj(var), SCIPvarGetLbLocal(var), newub) );
      }
   }
   else
   {
      if( !SCIPsetIsEQ(set, oldub, newub) && SCIPsetIsNegative(set, SCIPvarGetObj(var)) )
      {
         SCIP_CALL( lpUpdateVar(lp, set, var, SCIPvarGetObj(var), SCIPvarGetLbLocal(var), oldub, 
                        SCIPvarGetObj(var), SCIPvarGetLbLocal(var), newub) );
      }
   }

   return SCIP_OKAY;
}

/** informs LP, that given variable was added to the problem */
SCIP_RETCODE SCIPlpUpdateAddVar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that is now a LOOSE problem variable */
   )
{
   assert(lp != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);

   /* add the variable to the loose objective value sum */
   SCIP_CALL( SCIPlpUpdateVarObj(lp, set, var, 0.0, SCIPvarGetObj(var)) );

   /* update the loose variables counter */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      lp->nloosevars++;

   return SCIP_OKAY;
}

/** informs LP, that given variable is to be deleted from the problem */
SCIP_RETCODE SCIPlpUpdateDelVar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that will be deleted from the problem */
   )
{
   assert(lp != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);

   /* subtract the variable from the loose objective value sum */
   SCIP_CALL( SCIPlpUpdateVarObj(lp, set, var, SCIPvarGetObj(var), 0.0) );

   /* update the loose variables counter */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      lp->nloosevars--;

   return SCIP_OKAY;
}

/** informs LP, that given formerly loose problem variable is now a column variable */
static
SCIP_RETCODE lpUpdateVarColumn(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   )
{
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(lp != NULL);
   assert(lp->nloosevars > 0);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);

   obj = SCIPvarGetObj(var);

   /* update loose objective value corresponding to the deletion of variable */
   if( SCIPsetIsPositive(set, obj) )
   {
      lb = SCIPvarGetLbLocal(var);
      if( SCIPsetIsInfinity(set, -lb) )
         lp->looseobjvalinf--;
      else
         lp->looseobjval -= lb * obj;
   }
   else if( SCIPsetIsNegative(set, obj) )
   {
      ub = SCIPvarGetUbLocal(var);
      if( SCIPsetIsInfinity(set, ub) )
         lp->looseobjvalinf--;
      else
         lp->looseobjval -= ub * obj;
   }
   lp->nloosevars--;

   /* get rid of numerical problems: set loose objective value explicitly to zero, if no loose variables remain */
   if( lp->nloosevars == 0 )
   {
      assert(lp->looseobjvalinf == 0);
      lp->looseobjval = 0.0;
   }

   return SCIP_OKAY;
}

/** informs LP, that given formerly loose problem variable is now a column variable
 *  pseudo objective value is calculated with interval arithmetics to get a proved lower bound
 */
static
SCIP_RETCODE lpUpdateVarColumnProved(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   )
{
   SCIP_INTERVAL bd;
   SCIP_INTERVAL ob;
   SCIP_INTERVAL prod;
   SCIP_INTERVAL loose;
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(lp != NULL);
   assert(lp->nloosevars > 0);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);

   obj = SCIPvarGetObj(var);

   SCIPintervalSet(&loose, lp->looseobjval);

   /* update loose objective value corresponding to the deletion of variable */
   if( obj > 0.0 )
   {
      lb = SCIPvarGetLbLocal(var);
      if( SCIPsetIsInfinity(set, -lb) )
         lp->looseobjvalinf--;
      else
      {
         SCIPintervalSet(&bd, lb);
         SCIPintervalSet(&ob, obj);
         SCIPintervalMul(&prod, bd, ob);
         SCIPintervalSub(&loose, loose, prod);  /* lp->looseobjval -= lb * obj; */
      }
   }
   else if( SCIPsetIsNegative(set, obj) )
   {
      ub = SCIPvarGetUbLocal(var);
      if( SCIPsetIsInfinity(set, ub) )
         lp->looseobjvalinf--;
      else
      {
         SCIPintervalSet(&bd, ub);
         SCIPintervalSet(&ob, obj);
         SCIPintervalMul(&prod, bd, ob);
         SCIPintervalSub(&loose, loose, prod);  /* lp->looseobjval -= ub * obj; */
      }
   }
   lp->nloosevars--;

   /* get rid of numerical problems: set loose objective value explicitly to zero, if no loose variables remain */
   if( lp->nloosevars == 0 )
   {
      assert(lp->looseobjvalinf == 0);
      lp->looseobjval = 0.0;
   }
   else
      lp->looseobjval = SCIPintervalGetInf(loose);

   return SCIP_OKAY;
}

/** informs LP, that given formerly loose problem variable is now a column variable */
SCIP_RETCODE SCIPlpUpdateVarColumn(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   )
{
   assert(set != NULL);

   if( set->misc_exactsolve )
   {
      SCIP_CALL( lpUpdateVarColumnProved(lp, set, var) );
   }
   else
   {
      SCIP_CALL( lpUpdateVarColumn(lp, set, var) );
   }

   return SCIP_OKAY;
}

/** informs LP, that given formerly column problem variable is now again a loose variable */
static
SCIP_RETCODE lpUpdateVarLoose(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   )
{
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(lp != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetProbindex(var) >= 0);

   obj = SCIPvarGetObj(var);

   /* update loose objective value corresponding to the addition of variable */
   if( SCIPsetIsPositive(set, obj) )
   {
      lb = SCIPvarGetLbLocal(var);
      if( SCIPsetIsInfinity(set, -lb) )
         lp->looseobjvalinf++;
      else
         lp->looseobjval += lb * obj;
   }
   else if( SCIPsetIsNegative(set, obj) )
   {
      ub = SCIPvarGetUbLocal(var);
      if( SCIPsetIsInfinity(set, ub) )
         lp->looseobjvalinf++;
      else
         lp->looseobjval += ub * obj;
   }
   lp->nloosevars++;

   return SCIP_OKAY;
}

/** informs LP, that given formerly column problem variable is now again a loose variable
 *  pseudo objective value is calculated with interval arithmetics to get a proved lower bound
 */
static
SCIP_RETCODE lpUpdateVarLooseProved(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   )
{
   SCIP_INTERVAL bd;
   SCIP_INTERVAL ob;
   SCIP_INTERVAL prod;
   SCIP_INTERVAL loose;
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(lp != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetProbindex(var) >= 0);

   obj = SCIPvarGetObj(var);

   SCIPintervalSet(&loose, lp->looseobjval);

   /* update loose objective value corresponding to the deletion of variable */
   if( obj > 0.0 )
   {
      lb = SCIPvarGetLbLocal(var);
      if( SCIPsetIsInfinity(set, -lb) )
         lp->looseobjvalinf++;
      else
      {
         SCIPintervalSet(&bd, lb);
         SCIPintervalSet(&ob, obj);
         SCIPintervalMul(&prod, bd, ob);
         SCIPintervalAdd(&loose, loose, prod);  /* lp->looseobjval += lb * obj; */
      }
   }
   else if( SCIPsetIsNegative(set, obj) )
   {
      ub = SCIPvarGetUbLocal(var);
      if( SCIPsetIsInfinity(set, ub) )
         lp->looseobjvalinf++;
      else
      {
         SCIPintervalSet(&bd, ub);
         SCIPintervalSet(&ob, obj);
         SCIPintervalMul(&prod, bd, ob);
         SCIPintervalAdd(&loose, loose, prod);  /* lp->looseobjval += ub * obj; */
      }
   }
   lp->nloosevars++;

   lp->looseobjval = SCIPintervalGetInf(loose);

   return SCIP_OKAY;
}

/** informs LP, that given formerly column problem variable is now again a loose variable */
SCIP_RETCODE SCIPlpUpdateVarLoose(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   )
{
   assert(set != NULL);

   if( set->misc_exactsolve )
   {
      SCIP_CALL( lpUpdateVarLooseProved(lp, set, var) );
   }
   else
   {
      SCIP_CALL( lpUpdateVarLoose(lp, set, var) );
   }

   return SCIP_OKAY;
}

/** stores the LP solution in the columns and rows */
SCIP_RETCODE SCIPlpGetSol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            dualfeasible        /**< pointer to store whether the solution is dual feasible, or NULL */
   )
{
   SCIP_COL** lpicols;
   SCIP_ROW** lpirows;
   SCIP_Real* primsol;
   SCIP_Real* dualsol;
   SCIP_Real* activity;
   SCIP_Real* redcost;
   int* cstat;
   int* rstat;
   int nlpicols;
   int nlpirows;
   int lpcount;
   int c;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp->validsollp <= stat->lpcount);

   if( primalfeasible != NULL )
      *primalfeasible = TRUE;
   if( dualfeasible != NULL )
      *dualfeasible = TRUE;

   /* check if the values are already calculated */
   if( lp->validsollp == stat->lpcount )
      return SCIP_OKAY;
   lp->validsollp = stat->lpcount;

   SCIPdebugMessage("getting new LP solution %d\n", stat->lpcount);

   lpicols = lp->lpicols;
   lpirows = lp->lpirows;
   nlpicols = lp->nlpicols;
   nlpirows = lp->nlpirows;
   lpcount = stat->lpcount;

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &primsol, nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualsol, nlpirows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &activity, nlpirows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &redcost, nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cstat, nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rstat, nlpirows) );
   
   SCIP_CALL( SCIPlpiGetSol(lp->lpi, NULL, primsol, dualsol, activity, redcost) );
   if( lp->solisbasic )
   {
      SCIP_CALL( SCIPlpiGetBase(lp->lpi, cstat, rstat) );
   }
   else
   {
      BMSclearMemoryArray(cstat, nlpicols);
      BMSclearMemoryArray(rstat, nlpirows);
   }

   /* copy primal solution and reduced costs into columns */
   for( c = 0; c < nlpicols; ++c )
   {
      lpicols[c]->primsol = primsol[c];
      lpicols[c]->minprimsol = MIN(lpicols[c]->minprimsol, primsol[c]);
      lpicols[c]->maxprimsol = MAX(lpicols[c]->maxprimsol, primsol[c]);
      lpicols[c]->redcost = redcost[c];
      lpicols[c]->basisstatus = cstat[c]; /*lint !e732*/
      lpicols[c]->validredcostlp = lpcount;
      if( primalfeasible != NULL )
         *primalfeasible = *primalfeasible
            && !SCIPsetIsFeasNegative(set, lpicols[c]->primsol - lpicols[c]->lb)
            && !SCIPsetIsFeasPositive(set, lpicols[c]->primsol - lpicols[c]->ub);
      if( dualfeasible != NULL )
      {
         if( SCIPsetIsFeasGT(set, lpicols[c]->primsol, lpicols[c]->lb) )
            *dualfeasible = *dualfeasible && !SCIPsetIsFeasPositive(set, lpicols[c]->redcost);
         if( SCIPsetIsFeasLT(set, lpicols[c]->primsol, lpicols[c]->ub) )
            *dualfeasible = *dualfeasible && !SCIPsetIsFeasNegative(set, lpicols[c]->redcost);
      }
      SCIPdebugMessage(" col <%s> [%g,%g]: primsol=%.9f, redcost=%.9f, pfeas=%d/%d(%d), dfeas=%d(%d)\n",
         SCIPvarGetName(lpicols[c]->var), lpicols[c]->lb, lpicols[c]->ub, lpicols[c]->primsol, lpicols[c]->redcost,
         SCIPsetIsFeasGE(set, lpicols[c]->primsol, lpicols[c]->lb),
         SCIPsetIsFeasLE(set, lpicols[c]->primsol, lpicols[c]->ub), 
         primalfeasible != NULL ? *primalfeasible : TRUE,
         !SCIPsetIsFeasNegative(set, lpicols[c]->redcost),
         dualfeasible != NULL ? *dualfeasible : TRUE);
   }
   
   /* copy dual solution and activities into rows */
   for( r = 0; r < nlpirows; ++r )
   {
      lpirows[r]->dualsol = dualsol[r];
      lpirows[r]->activity = activity[r] + lpirows[r]->constant;
      lpirows[r]->basisstatus = rstat[r]; /*lint !e732*/
      lpirows[r]->validactivitylp = lpcount;
      if( primalfeasible != NULL )
         *primalfeasible = *primalfeasible
            && SCIPsetIsFeasGE(set, lpirows[r]->activity, lpirows[r]->lhs)
            && SCIPsetIsFeasLE(set, lpirows[r]->activity, lpirows[r]->rhs);
      if( dualfeasible != NULL )
      {
         if( SCIPsetIsInfinity(set, -lpirows[r]->lhs) )
            *dualfeasible = *dualfeasible && !SCIPsetIsFeasPositive(set, lpirows[r]->dualsol);
         if( SCIPsetIsInfinity(set, lpirows[r]->rhs) )
            *dualfeasible = *dualfeasible && !SCIPsetIsFeasNegative(set, lpirows[r]->dualsol);
      }
      SCIPdebugMessage(" row <%s> [%g,%g]: dualsol=%.9f, activity=%.9f, pfeas=%d/%d(%d), dfeas=%d/%d(%d)\n", 
         lpirows[r]->name, lpirows[r]->lhs, lpirows[r]->rhs, lpirows[r]->dualsol, lpirows[r]->activity,
         SCIPsetIsFeasGE(set, lpirows[r]->activity, lpirows[r]->lhs),
         SCIPsetIsFeasLE(set, lpirows[r]->activity, lpirows[r]->rhs),
         primalfeasible != NULL ? *primalfeasible : TRUE,
         SCIPsetIsInfinity(set, -lpirows[r]->lhs) ? !SCIPsetIsFeasPositive(set, lpirows[r]->dualsol) : TRUE,
         SCIPsetIsInfinity(set, lpirows[r]->rhs) ? !SCIPsetIsFeasNegative(set, lpirows[r]->dualsol) : TRUE,
         dualfeasible != NULL ? *dualfeasible : TRUE);
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &rstat);
   SCIPsetFreeBufferArray(set, &cstat);
   SCIPsetFreeBufferArray(set, &redcost);
   SCIPsetFreeBufferArray(set, &activity);
   SCIPsetFreeBufferArray(set, &dualsol);
   SCIPsetFreeBufferArray(set, &primsol);

   return SCIP_OKAY;
}

/** stores LP solution with infinite objective value in the columns and rows */
SCIP_RETCODE SCIPlpGetUnboundedSol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_COL** lpicols;
   SCIP_ROW** lpirows;
   SCIP_Real* primsol;
   SCIP_Real* activity;
   SCIP_Real* ray;
   SCIP_Real rayobjval;
   SCIP_Real rayscale;
   int nlpicols;
   int nlpirows;
   int lpcount;
   int c;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_UNBOUNDEDRAY);
   assert(SCIPsetIsInfinity(set, -lp->lpobjval));
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp->validsollp <= stat->lpcount);

   /* check if the values are already calculated */
   if( lp->validsollp == stat->lpcount )
      return SCIP_OKAY;
   lp->validsollp = stat->lpcount;

   SCIPdebugMessage("getting new unbounded LP solution %d\n", stat->lpcount);

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &primsol, lp->nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &activity, lp->nlpirows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ray, lp->nlpicols) );

   /* get primal feasible point */
   SCIP_CALL( SCIPlpiGetSol(lp->lpi, NULL, primsol, NULL, activity, NULL) );

   /* get primal unbounded ray */
   SCIP_CALL( SCIPlpiGetPrimalRay(lp->lpi, ray) );
   
   lpicols = lp->lpicols;
   lpirows = lp->lpirows;
   nlpicols = lp->nlpicols;
   nlpirows = lp->nlpirows;
   lpcount = stat->lpcount;

   /* calculate the objective value decrease of the ray */
   rayobjval = 0.0;
   for( c = 0; c < nlpicols; ++c )
   {
      assert(lpicols[c] != NULL);
      assert(lpicols[c]->var != NULL);
      rayobjval += ray[c] * lpicols[c]->obj;
   }
   assert(SCIPsetIsNegative(set, rayobjval));

   /* scale the ray, such that the resulting point has infinite objective value */
   rayscale = -2*SCIPsetInfinity(set)/rayobjval;

   /* calculate the unbounded point: x' = x + rayscale * ray */
   SCIPdebugMessage("unbounded LP solution: rayobjval=%f, rayscale=%f\n", rayobjval, rayscale);

   for( c = 0; c < nlpicols; ++c )
   {
      lpicols[c]->primsol = primsol[c] + rayscale * ray[c];
      lpicols[c]->redcost = SCIP_INVALID;
      lpicols[c]->validredcostlp = -1;
      /*debugMessage(" col <%s>: basesol=%f, ray=%f, unbdsol=%f\n", 
        SCIPvarGetName(lpicols[c]->var), primsol[c], ray[c], lpicols[c]->primsol);*/
   }

   for( r = 0; r < nlpirows; ++r )
   {
      lpirows[r]->dualsol = SCIP_INVALID;
      lpirows[r]->activity = activity[r] + lpirows[r]->constant;
      lpirows[r]->validactivitylp = lpcount;
      /*debugMessage(" row <%s>: activity=%f\n", lpirows[r]->name, lpirows[r]->activity);*/
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &ray);
   SCIPsetFreeBufferArray(set, &activity);
   SCIPsetFreeBufferArray(set, &primsol);

   return SCIP_OKAY;
}

/** stores the dual farkas multipliers for infeasibility proof in rows */
SCIP_RETCODE SCIPlpGetDualfarkas(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_ROW** lpirows;
   SCIP_Real* dualfarkas;
   int nlpirows;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp->validfarkaslp <= stat->lpcount);

   /* check if the values are already calculated */
   if( lp->validfarkaslp == stat->lpcount )
      return SCIP_OKAY;
   lp->validfarkaslp = stat->lpcount;

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualfarkas, lp->nlpirows) );

   /* get dual farkas infeasibility proof */
   SCIP_CALL( SCIPlpiGetDualfarkas(lp->lpi, dualfarkas) );

   lpirows = lp->lpirows;
   nlpirows = lp->nlpirows;

   /* store infeasibility proof in rows */
   SCIPdebugMessage("LP is infeasible:\n");
   for( r = 0; r < nlpirows; ++r )
   {
      SCIPdebugMessage(" row <%s>: dualfarkas=%f\n", lpirows[r]->name, dualfarkas[r]);
      lpirows[r]->dualfarkas = dualfarkas[r];
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &dualfarkas);

   return SCIP_OKAY;
}

/** get number of iterations used in last LP solve */
SCIP_RETCODE SCIPlpGetIterations(
   SCIP_LP*              lp,                 /**< current LP data */
   int*                  iterations          /**< pointer to store the iteration count */
   )
{
   assert(lp != NULL);

   SCIP_CALL( SCIPlpiGetIterations(lp->lpi, iterations) );

   return SCIP_OKAY;
}

/** increases age of columns with solution value 0.0 and basic rows with activity not at its bounds,
 *  resets age of non-zero columns and sharp rows
 */
SCIP_RETCODE SCIPlpUpdateAges(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_COL** lpicols;
   SCIP_ROW** lpirows;
   int nlpicols;
   int nlpirows;
   int c;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->nlpicols == lp->ncols);
   assert(lp->nlpirows == lp->nrows);
   assert(stat != NULL);
   assert(lp->validsollp == stat->lpcount);

   SCIPdebugMessage("updating LP ages\n");

   lpicols = lp->lpicols;
   lpirows = lp->lpirows;
   nlpicols = lp->nlpicols;
   nlpirows = lp->nlpirows;

   for( c = 0; c < nlpicols; ++c )
   {
      assert(lpicols[c] == lp->cols[c]);
      if( lpicols[c]->primsol == 0.0 )  /* non-basic columns to remove are exactly at 0.0 */
         lpicols[c]->age++;
      else
         lpicols[c]->age = 0;
      /*debugMessage(" -> col <%s>: primsol=%f, age=%d\n", 
        SCIPvarGetName(lpicols[c]->var), lpicols[c]->primsol, lpicols[c]->age);*/
   }

   for( r = 0; r < nlpirows; ++r )
   {
      assert(lpirows[r] == lp->rows[r]);
      if( lpirows[r]->dualsol == 0.0 )  /* basic rows to remove are exactly at 0.0 */
         lpirows[r]->age++;
      else
         lpirows[r]->age = 0;
      /*debugMessage(" -> row <%s>: activity=%f, age=%d\n", lpirows[r]->name, lpirows[r]->activity, lpirows[r]->age);*/
   }

   return SCIP_OKAY;
}

/* deletes the marked columns from the LP and the LP interface */
static
SCIP_RETCODE lpDelColset(
   SCIP_LP*              lp,                 /**< current LP data */
   int*                  coldstat            /**< deletion status of columns:  1 if column should be deleted, 0 if not */
   )
{
   SCIP_COL* col;
   int ncols;
   int c;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->ncols == lp->nlpicols);
   assert(!lp->diving);
   assert(coldstat != NULL);

   ncols = lp->ncols;

   /* delete columns in LP solver */
   SCIP_CALL( SCIPlpiDelColset(lp->lpi, coldstat) );

   /* update LP data respectively */
   for( c = 0; c < ncols; ++c )
   {
      col = lp->cols[c];
      assert(col == lp->lpicols[c]);
      assert(coldstat[c] <= c);
      col->lppos = coldstat[c];
      if( coldstat[c] == -1 )
      {
         assert(col != NULL);
         assert(col->removeable);

         /* mark column to be deleted from the LPI and update column arrays of all linked rows */
         markColDeleted(col);
         colUpdateDelLP(col);
         col->lpdepth = -1;

         lp->cols[c] = NULL;
         lp->lpicols[c] = NULL;
         lp->ncols--;
         lp->nremoveablecols--;
         lp->nlpicols--;
      }
      else if( coldstat[c] < c )
      {
         assert(lp->cols[coldstat[c]] == NULL);
         assert(lp->lpicols[coldstat[c]] == NULL);
         lp->cols[coldstat[c]] = col;
         lp->lpicols[coldstat[c]] = col;
         lp->cols[coldstat[c]]->lppos = coldstat[c];
         lp->cols[coldstat[c]]->lpipos = coldstat[c];
         lp->cols[c] = NULL;
         lp->lpicols[c] = NULL;
      }
   }

   /* mark LP to be unsolved */
   if( lp->ncols < ncols )
   {
      assert(lp->ncols == lp->nlpicols);
      assert(lp->nchgcols == 0);
      assert(lp->flushed);

      lp->lpifirstchgcol = lp->nlpicols;

      /* mark the current solution invalid */
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/* deletes the marked rows from the LP and the LP interface */
static
SCIP_RETCODE lpDelRowset(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   int*                  rowdstat            /**< deletion status of rows:  1 if row should be deleted, 0 if not */
   )
{
   SCIP_ROW* row;
   int nrows;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->nrows == lp->nlpirows);
   assert(!lp->diving);
   assert(rowdstat != NULL);

   nrows = lp->nrows;

   /* delete rows in LP solver */
   SCIP_CALL( SCIPlpiDelRowset(lp->lpi, rowdstat) );

   /* update LP data respectively */
   for( r = 0; r < nrows; ++r )
   {
      row = lp->rows[r];
      assert(row == lp->lpirows[r]);
      assert(rowdstat[r] <= r);
      row->lppos = rowdstat[r];
      if( rowdstat[r] == -1 )
      {
         assert(row != NULL);
         assert(row->removeable);

         /* mark row to be deleted from the LPI and update row arrays of all linked columns */
         markRowDeleted(row);
         rowUpdateDelLP(row);
         row->lpdepth = -1;

         SCIP_CALL( SCIProwRelease(&lp->lpirows[r], blkmem, set, lp) );
         SCIProwUnlock(lp->rows[r]);
         SCIP_CALL( SCIProwRelease(&lp->rows[r], blkmem, set, lp) );
         assert(lp->lpirows[r] == NULL);
         assert(lp->rows[r] == NULL);
         lp->nrows--;
         lp->nremoveablerows--;
         lp->nlpirows--;
      }
      else if( rowdstat[r] < r )
      {
         assert(lp->rows[rowdstat[r]] == NULL);
         assert(lp->lpirows[rowdstat[r]] == NULL);
         lp->rows[rowdstat[r]] = row;
         lp->lpirows[rowdstat[r]] = row;
         lp->rows[rowdstat[r]]->lppos = rowdstat[r];
         lp->rows[rowdstat[r]]->lpipos = rowdstat[r];
         lp->rows[r] = NULL;
         lp->lpirows[r] = NULL;
      }
   }

   /* mark LP to be unsolved */
   if( lp->nrows < nrows )
   {
      assert(lp->nrows == lp->nlpirows);
      assert(lp->nchgrows == 0);
      assert(lp->flushed);

      lp->lpifirstchgrow = lp->nlpirows;

      /* mark the current solution invalid */
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** removes all non-basic columns, that are too old, beginning with the given firstcol */
static
SCIP_RETCODE lpRemoveObsoleteCols(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   int                   firstcol            /**< first column to check for clean up */
   )
{
   SCIP_COL** cols;
   SCIP_COL** lpicols;
   int* coldstat;
   int ncols;
   int ndelcols;
   int c;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->ncols == lp->nlpicols);
   assert(lp->nremoveablecols <= lp->ncols);
   assert(!lp->diving);
   assert(set != NULL);
   assert(stat != NULL);

   if( lp->nremoveablecols == 0 || set->lp_colagelimit == -1 || !lp->solisbasic )
      return SCIP_OKAY;

   ncols = lp->ncols;
   cols = lp->cols;
   lpicols = lp->lpicols;

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &coldstat, ncols) );

   /* mark obsolete columns to be deleted */
   ndelcols = 0;
   BMSclearMemoryArray(coldstat, ncols);
   for( c = firstcol; c < ncols; ++c )
   {
      assert(cols[c] == lpicols[c]);
      assert(cols[c]->lppos == c);
      assert(cols[c]->lpipos == c);
      if( cols[c]->removeable
         && cols[c]->obsoletenode != stat->nnodes /* don't remove column a second time from same node (avoid cycling) */
         && cols[c]->age > set->lp_colagelimit
         && (SCIP_BASESTAT)cols[c]->basisstatus != SCIP_BASESTAT_BASIC
         && SCIPsetIsZero(set, SCIPcolGetBestBound(cols[c])) ) /* bestbd != 0 -> column would be priced in next time */
      {
         assert(cols[c]->primsol == 0.0);
         coldstat[c] = 1;
         ndelcols++;
         cols[c]->obsoletenode = stat->nnodes;
         SCIPdebugMessage("removing obsolete col <%s>: primsol=%f, bounds=[%g,%g]\n", 
            SCIPvarGetName(cols[c]->var), cols[c]->primsol, cols[c]->lb, cols[c]->ub);
      }
   }

   SCIPdebugMessage("removing %d/%d obsolete columns from LP\n", ndelcols, ncols);

   /* delete the marked columns in the LP solver interface, update the LP respectively */
   if( ndelcols > 0 )
   {
      SCIP_CALL( lpDelColset(lp, coldstat) );
   }
   assert(lp->ncols == ncols - ndelcols);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &coldstat);
      
   return SCIP_OKAY;
}

/** removes all basic rows, that are too old, beginning with the given firstrow */
static
SCIP_RETCODE lpRemoveObsoleteRows(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   int                   firstrow            /**< first row to check for clean up */
   )
{
   SCIP_ROW** rows;
   SCIP_ROW** lpirows;
   int* rowdstat;
   int nrows;
   int ndelrows;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->nrows == lp->nlpirows);
   assert(lp->nremoveablerows <= lp->nrows);
   assert(!lp->diving);
   assert(set != NULL);
   assert(stat != NULL);

   if( lp->nremoveablerows == 0 || set->lp_rowagelimit == -1 || !lp->solisbasic )
      return SCIP_OKAY;

   nrows = lp->nrows;
   rows = lp->rows;
   lpirows = lp->lpirows;

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rowdstat, nrows) );

   /* mark obsolete rows to be deleted */
   ndelrows = 0;
   BMSclearMemoryArray(rowdstat, nrows);
   for( r = firstrow; r < nrows; ++r )
   {
      assert(rows[r] == lpirows[r]);
      assert(rows[r]->lppos == r);
      assert(rows[r]->lpipos == r);
      if( rows[r]->removeable
         && rows[r]->obsoletenode != stat->nnodes  /* don't remove row a second time from same node (avoid cycling) */
         && rows[r]->age > set->lp_rowagelimit
         && (SCIP_BASESTAT)rows[r]->basisstatus == SCIP_BASESTAT_BASIC )
      {
         rowdstat[r] = 1;
         ndelrows++;
         rows[r]->obsoletenode = stat->nnodes;
         SCIPdebugMessage("removing obsolete row <%s>: activity=%f, sides=[%g,%g]\n", 
            rows[r]->name, rows[r]->activity, rows[r]->lhs, rows[r]->rhs);
      }
   }

   SCIPdebugMessage("removing %d/%d obsolete rows from LP\n", ndelrows, nrows);

   /* delete the marked rows in the LP solver interface, update the LP respectively */
   if( ndelrows > 0 )
   {
      SCIP_CALL( lpDelRowset(lp, blkmem, set, rowdstat) );
   }
   assert(lp->nrows == nrows - ndelrows);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &rowdstat);
      
   return SCIP_OKAY;
}

/** removes all non-basic columns and basic rows in the part of the LP created at the current node, that are too old */
SCIP_RETCODE SCIPlpRemoveNewObsoletes(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert(!lp->diving);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(set != NULL);

   SCIPdebugMessage("removing obsolete columns starting with %d/%d, obsolete rows starting with %d/%d\n",
      lp->firstnewcol, lp->ncols, lp->firstnewrow, lp->nrows);

   if( lp->firstnewcol < lp->ncols )
   {
      SCIP_CALL( lpRemoveObsoleteCols(lp, set, stat, lp->firstnewcol) );
   }
   if( lp->firstnewrow < lp->nrows )
   {
      SCIP_CALL( lpRemoveObsoleteRows(lp, blkmem, set, stat, lp->firstnewrow) );
   }

   return SCIP_OKAY;
}

/** removes all non-basic columns and basic rows in whole LP, that are too old */
SCIP_RETCODE SCIPlpRemoveAllObsoletes(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert(!lp->diving);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(set != NULL);

   SCIPdebugMessage("removing all obsolete columns and rows\n");

   if( 0 < lp->ncols )
   {
      SCIP_CALL( lpRemoveObsoleteCols(lp, set, stat, 0) );
   }
   if( 0 < lp->nrows )
   {
      SCIP_CALL( lpRemoveObsoleteRows(lp, blkmem, set, stat, 0) );
   }

   return SCIP_OKAY;
}

/** removes all non-basic columns at 0.0 beginning with the given firstcol */
static
SCIP_RETCODE lpCleanupCols(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   int                   firstcol            /**< first column to check for clean up */
   )
{
   SCIP_COL** cols;
   SCIP_COL** lpicols;
   int* coldstat;
   int ncols;
   int ndelcols;
   int c;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->ncols == lp->nlpicols);
   assert(!lp->diving);
   assert(stat != NULL);
   assert(lp->validsollp == stat->lpcount);
   assert(0 <= firstcol && firstcol < lp->ncols);

   if( lp->nremoveablecols == 0 || !lp->solisbasic )
      return SCIP_OKAY;

   ncols = lp->ncols;
   cols = lp->cols;
   lpicols = lp->lpicols;

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &coldstat, ncols) );

   /* mark unused columns to be deleted */
   ndelcols = 0;
   BMSclearMemoryArray(coldstat, ncols);
   for( c = firstcol; c < ncols; ++c )
   {
      assert(cols[c] == lpicols[c]);
      assert(cols[c]->lppos == c);
      assert(cols[c]->lpipos == c);
      if( lpicols[c]->removeable
         && (SCIP_BASESTAT)lpicols[c]->basisstatus != SCIP_BASESTAT_BASIC
         && lpicols[c]->primsol == 0.0 /* non-basic columns to remove are exactly at 0.0 */
         && SCIPsetIsZero(set, SCIPcolGetBestBound(cols[c])) ) /* bestbd != 0 -> column would be priced in next time */
      {
         coldstat[c] = 1;
         ndelcols++;
      }
   }

   SCIPdebugMessage("removing %d/%d unused columns from LP\n", ndelcols, ncols);

   /* delete the marked columns in the LP solver interface, update the LP respectively */
   if( ndelcols > 0 )
   {
      SCIP_CALL( lpDelColset(lp, coldstat) );
   }
   assert(lp->ncols == ncols - ndelcols);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &coldstat);
      
   return SCIP_OKAY;
}

/** removes all basic rows beginning with the given firstrow */
static
SCIP_RETCODE lpCleanupRows(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   int                   firstrow            /**< first row to check for clean up */
   )
{
   SCIP_ROW** rows;
   SCIP_ROW** lpirows;
   int* rowdstat;
   int nrows;
   int ndelrows;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->ncols == lp->nlpicols);
   assert(lp->nrows == lp->nlpirows);
   assert(!lp->diving);
   assert(stat != NULL);
   assert(lp->validsollp == stat->lpcount);
   assert(0 <= firstrow && firstrow < lp->nrows);

   if( lp->nremoveablerows == 0 || !lp->solisbasic  )
      return SCIP_OKAY;

   nrows = lp->nrows;
   rows = lp->rows;
   lpirows = lp->lpirows;

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rowdstat, nrows) );

   /* mark unused rows to be deleted */
   ndelrows = 0;
   BMSclearMemoryArray(rowdstat, nrows);
   for( r = firstrow; r < nrows; ++r )
   {
      assert(rows[r] == lpirows[r]);
      assert(rows[r]->lppos == r);
      assert(rows[r]->lpipos == r);
      if( lpirows[r]->removeable && (SCIP_BASESTAT)lpirows[r]->basisstatus == SCIP_BASESTAT_BASIC )
      {
         rowdstat[r] = 1;
         ndelrows++;
      }
   }

   SCIPdebugMessage("removing %d/%d unused rows from LP\n", ndelrows, nrows);

   /* delete the marked rows in the LP solver interface, update the LP respectively */
   if( ndelrows > 0 )
   {
      SCIP_CALL( lpDelRowset(lp, blkmem, set, rowdstat) );
   }
   assert(lp->nrows == nrows - ndelrows);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &rowdstat);

   return SCIP_OKAY;
}

/** removes all non-basic columns at 0.0 and basic rows in the part of the LP created at the current node */
SCIP_RETCODE SCIPlpCleanupNew(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             root                /**< are we at the root node? */
   )
{
   SCIP_Bool cleanupcols;
   SCIP_Bool cleanuprows;

   assert(lp != NULL);
   assert(lp->solved);
   assert(!lp->diving);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(set != NULL);

   /* check, if we want to clean up the columns and rows */
   cleanupcols = (root ? set->lp_cleanupcolsroot : set->lp_cleanupcols);
   cleanuprows = (root ? set->lp_cleanuprowsroot : set->lp_cleanuprows);

   SCIPdebugMessage("removing unused columns starting with %d/%d (%d), unused rows starting with %d/%d (%d), LP algo: %d, basic sol: %d\n",
      lp->firstnewcol, lp->ncols, cleanupcols, lp->firstnewrow, lp->nrows, cleanuprows, lp->lastlpalgo, lp->solisbasic);

   if( cleanupcols && lp->firstnewcol < lp->ncols )
   {
      SCIP_CALL( lpCleanupCols(lp, set, stat, lp->firstnewcol) );
   }
   if( cleanuprows && lp->firstnewrow < lp->nrows )
   {
      SCIP_CALL( lpCleanupRows(lp, blkmem, set, stat, lp->firstnewrow) );
   }

   return SCIP_OKAY;
}

/** removes all non-basic columns at 0.0 and basic rows in the whole LP */
SCIP_RETCODE SCIPlpCleanupAll(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             root                /**< are we at the root node? */
   )
{
   SCIP_Bool cleanupcols;
   SCIP_Bool cleanuprows;

   assert(lp != NULL);
   assert(lp->solved);
   assert(!lp->diving);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(set != NULL);

   /* check, if we want to clean up the columns and rows */
   cleanupcols = (root ? set->lp_cleanupcolsroot : set->lp_cleanupcols);
   cleanuprows = (root ? set->lp_cleanuprowsroot : set->lp_cleanuprows);

   SCIPdebugMessage("removing all unused columns (%d) and rows (%d), LP algo: %d, basic sol: %d\n", 
      cleanupcols, cleanuprows, lp->lastlpalgo, lp->solisbasic);

   if( cleanupcols && 0 < lp->ncols )
   {
      SCIP_CALL( lpCleanupCols(lp, set, stat, 0) );
   }
   if( cleanuprows && 0 < lp->nrows )
   {
      SCIP_CALL( lpCleanupRows(lp, blkmem, set, stat, 0) );
   }

   return SCIP_OKAY;
}

/** initiates LP diving */
SCIP_RETCODE SCIPlpStartDive(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert(lp->flushed);
   assert(!lp->diving);
   assert(!lp->probing);
   assert(lp->divelpistate == NULL);
   assert(set != NULL);

   SCIPdebugMessage("diving started (LP flushed: %d, LP solved: %d, solstat: %d)\n",
      lp->flushed, lp->solved, SCIPlpGetSolstat(lp));

#ifndef NDEBUG
   {
      int c;
      for( c = 0; c < lp->ncols; ++c )
      {
         assert(lp->cols[c] != NULL);
         assert(lp->cols[c]->var != NULL);
         assert(SCIPvarGetStatus(lp->cols[c]->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetCol(lp->cols[c]->var) == lp->cols[c]);
         assert(SCIPsetIsFeasEQ(set, SCIPvarGetObj(lp->cols[c]->var), lp->cols[c]->obj));
         assert(SCIPsetIsFeasEQ(set, SCIPvarGetLbLocal(lp->cols[c]->var), lp->cols[c]->lb));
         assert(SCIPsetIsFeasEQ(set, SCIPvarGetUbLocal(lp->cols[c]->var), lp->cols[c]->ub));
      }
   }
#endif

   /* save current LPI state (basis information) */
   SCIP_CALL( SCIPlpiGetState(lp->lpi, blkmem, &lp->divelpistate) );

   /* switch to diving mode */
   lp->diving = TRUE;

   return SCIP_OKAY;
}

/** quits LP diving and resets bounds and objective values of columns to the current node's values */
SCIP_RETCODE SCIPlpEndDive(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_VAR**            vars,               /**< array with all active variables */
   int                   nvars               /**< number of active variables */
   )
{
   SCIP_VAR* var;
   SCIP_Bool lperror;
   int v;

   assert(lp != NULL);
   assert(lp->diving);
   assert(lp->divelpistate != NULL);
   assert(nvars == 0 || vars != NULL);

   SCIPdebugMessage("diving ended (LP flushed: %d, solstat: %d)\n",
      lp->flushed, SCIPlpGetSolstat(lp));
   
   /* reset all columns' objective values and bounds to its original values */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      {
         SCIP_CALL( SCIPcolChgObj(SCIPvarGetCol(var), set, lp, SCIPvarGetObj(var)) );
         SCIP_CALL( SCIPcolChgLb(SCIPvarGetCol(var), set, lp, SCIPvarGetLbLocal(var)) );
         SCIP_CALL( SCIPcolChgUb(SCIPvarGetCol(var), set, lp, SCIPvarGetUbLocal(var)) );
      }
   }

   /* reload LPI state saved at start of diving, free LPI state afterwards */
   SCIP_CALL( SCIPlpSetState(lp, blkmem, set, lp->divelpistate) );
   SCIP_CALL( SCIPlpFreeState(lp, blkmem, &lp->divelpistate) );
   assert(lp->divelpistate == NULL);

   /**@todo Get rid of resolving after diving: use separate data fields in columns to store all diving
    *       information (bounds, obj, solution values) and create calls SCIPvarGetDiveSol() etc.
    *       to access this diving LP information (not applicable on LOOSE variables).
    *       Just declare the LP to be solved at this point (remember the LP solution status beforehand).
    */
   /* resolve LP to reset solution */
   SCIP_CALL( SCIPlpSolveAndEval(lp, blkmem, set, stat, prob, -1, FALSE, FALSE, &lperror) );
   if( lperror )
   {
      SCIPmessagePrintVerbInfo(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) unresolved numerical troubles while resolving LP %d after diving\n", stat->nnodes, stat->nlps);
   }

   /* switch to standard (non-diving) mode */
   lp->diving = FALSE;
   lp->divingobjchg = FALSE;

#ifndef NDEBUG
   {
      int c;
      for( c = 0; c < lp->ncols; ++c )
      {
         assert(lp->cols[c] != NULL);
         assert(lp->cols[c]->var != NULL);
         assert(SCIPvarGetStatus(lp->cols[c]->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetCol(lp->cols[c]->var) == lp->cols[c]);
         assert(SCIPsetIsEQ(set, SCIPvarGetObj(lp->cols[c]->var), lp->cols[c]->obj));
         assert(SCIPsetIsEQ(set, SCIPvarGetLbLocal(lp->cols[c]->var), lp->cols[c]->lb));
         assert(SCIPsetIsEQ(set, SCIPvarGetUbLocal(lp->cols[c]->var), lp->cols[c]->ub));
      }
   }
#endif

   return SCIP_OKAY;
}

/** informs the LP that probing mode was initiated */
SCIP_RETCODE SCIPlpStartProbing(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(!lp->probing);

   lp->probing = TRUE;

   return SCIP_OKAY;
}

/** informs the LP that probing mode was finished */
SCIP_RETCODE SCIPlpEndProbing(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->probing);

   lp->probing = FALSE;

   return SCIP_OKAY;
}

/** calculates y*b + min{(c - y*A)*x | lb <= x <= ub} for given vectors y and c;
 *  the vector b is defined with b[i] = lhs[i] if y[i] >= 0, b[i] = rhs[i] if y[i] < 0
 *  Calculating this value in interval arithmetics gives a proved lower LP bound for the following reason (assuming,
 *  we have only left hand sides):
 *           min{cx       |  b <=  Ax, lb <= x <= ub}
 *   >=      min{cx       | yb <= yAx, lb <= x <= ub}   (restriction in minimum is relaxed)
 *   == yb + min{cx - yb  | yb <= yAx, lb <= x <= ub}   (added yb - yb == 0)
 *   >= yb + min{cx - yAx | yb <= yAx, lb <= x <= ub}   (because yAx >= yb inside minimum)
 *   >= yb + min{cx - yAx |            lb <= x <= ub}   (restriction in minimum is relaxed)
 */
static
SCIP_RETCODE provedBound(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             usefarkas,          /**< use y = dual farkas and c = 0 instead of y = dual solution and c = obj? */
   SCIP_Real*            bound               /**< result of interval arithmetic minimization */
   )
{
   SCIP_INTERVAL* yinter;
   SCIP_INTERVAL b;
   SCIP_INTERVAL ytb;
   SCIP_INTERVAL prod;
   SCIP_INTERVAL diff;
   SCIP_INTERVAL x;
   SCIP_INTERVAL minprod;
   SCIP_INTERVAL a;
   SCIP_ROW* row;
   SCIP_COL* col;
   SCIP_Real y;
   SCIP_Real c;
   int i;
   int j;

   assert(lp != NULL);
   assert(lp->solved);
   assert(set != NULL);
   assert(bound != NULL);

   /* allocate buffer for storing y in interval arithmetic */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &yinter, lp->nrows) );

   /* create y vector in interval arithmetic, setting near zeros to zero; calculate y^Tb */
   SCIPintervalSet(&ytb, 0.0);
   for( j = 0; j < lp->nrows; ++j )
   {
      row = lp->rows[j];
      assert(row != NULL);

      y = (usefarkas ? row->dualfarkas : row->dualsol);
         
      if( SCIPsetIsFeasPositive(set, y) )
      {
         SCIPintervalSet(&yinter[j], y);
         SCIPintervalSet(&b, row->lhs - row->constant);
      }
      else if( SCIPsetIsFeasNegative(set, y) )
      {
         SCIPintervalSet(&yinter[j], y);
         SCIPintervalSet(&b, row->rhs - row->constant);
      }
      else
      {
         SCIPintervalSet(&yinter[j], 0.0);
         SCIPintervalSet(&b, 0.0);
      }
      
      SCIPintervalMul(&prod, yinter[j], b);
      SCIPintervalAdd(&ytb, ytb, prod);
   }

   /* calculate min{(c^T - y^TA)x} */
   SCIPintervalSet(&minprod, 0.0);
   for( j = 0; j < lp->ncols; ++j )
   {
      col = lp->cols[j];
      assert(col != NULL);
      assert(col->nunlinked == 0);

      SCIPintervalSetBounds(&x, SCIPcolGetLb(col), SCIPcolGetUb(col));

      c = usefarkas ? 0.0 : col->obj;
      SCIPintervalSet(&diff, c);

      for( i = 0; i < col->nlprows; ++i )
      {
         assert(col->rows[i] != NULL);
         assert(col->rows[i]->lppos >= 0);
         assert(col->linkpos[i] >= 0);
         SCIPintervalSet(&a, col->vals[i]);
         SCIPintervalMul(&prod, yinter[col->rows[i]->lppos], a);
         SCIPintervalSub(&diff, diff, prod);
      }

#ifndef NDEBUG
      for( i = col->nlprows; i < col->len; ++i )
      {
         assert(col->rows[i] != NULL);
         assert(col->rows[i]->lppos == -1);
         assert(col->rows[i]->dualsol == 0.0);
         assert(col->rows[i]->dualfarkas == 0.0);
         assert(col->linkpos[i] >= 0);
      }
#endif

      SCIPintervalSetBounds(&x, col->lb, col->ub);
      SCIPintervalMul(&diff, diff, x);
      SCIPintervalAdd(&minprod, minprod, diff);
   }

   /* add y^Tb */
   SCIPintervalAdd(&minprod, minprod, ytb);

   /* free buffer for storing y in interval arithmetic */
   SCIPsetFreeBufferArray(set, &yinter);

   *bound = SCIPintervalGetInf(minprod);

   return SCIP_OKAY;
}

/** gets proven lower (dual) bound of last LP solution */
SCIP_RETCODE SCIPlpGetProvedLowerbound(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            bound               /**< pointer to store proven dual bound */
   )
{
   SCIP_CALL( provedBound(lp, set, FALSE, bound) );

   SCIPdebugMessage("proved lower bound of LP: %g\n", *bound);

   return SCIP_OKAY;
}

/** gets proven dual bound of last LP solution */
SCIP_RETCODE SCIPlpIsInfeasibilityProved(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool*            proved              /**< pointer to store whether infeasibility is proven */
   )
{
   SCIP_Real bound;

   assert(proved != NULL);

   SCIP_CALL( provedBound(lp, set, TRUE, &bound) );

   *proved = (bound > 0.0);

   SCIPdebugMessage("proved farkas value of LP: %g -> infeasibility %sproved\n", bound, *proved ? "" : "not ");

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpWrite(
   SCIP_LP*              lp,                 /**< current LP data */
   const char*           fname               /**< file name */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(fname != NULL);

   SCIP_CALL( SCIPlpiWriteLP(lp->lpi, fname) );

#if 0
   {
      FILE* f;
      int c;

      SCIPmessagePrintInfo("storing integrality conditions in <%s>\n", fname);
      f = fopen(fname, "a");
      SCIPmessageFPrintInfo(f, "General\n");
      for( c = 0; c < lp->ncols; ++c )
      {
         assert(lp->cols[c] == lp->lpicols[c]);
         if( SCIPvarGetType(SCIPcolGetVar(lp->cols[c])) != SCIP_VARTYPE_CONTINUOUS )
            SCIPmessageFPrintInfo(f, "%s\n", SCIPvarGetName(SCIPcolGetVar(lp->cols[c])));
      }
      SCIPmessageFPrintInfo(f, "End\n");
      fclose(f);
   }
#endif

   return SCIP_OKAY;
}




/*
 * simple functions implemented as defines
 */

/* In SCIPdebug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPcolGetObj
#undef SCIPcolGetLb
#undef SCIPcolGetUb
#undef SCIPcolGetBestBound
#undef SCIPcolGetPrimsol
#undef SCIPcolGetBasisStatus
#undef SCIPcolGetVar
#undef SCIPcolIsIntegral
#undef SCIPcolIsRemoveable
#undef SCIPcolGetLPPos
#undef SCIPcolGetLPDepth
#undef SCIPcolIsInLP
#undef SCIPcolGetNNonz
#undef SCIPcolGetNLPNonz
#undef SCIPcolGetRows
#undef SCIPcolGetVals
#undef SCIPcolGetStrongbranchNode
#undef SCIPcolGetNStrongbranchs
#undef SCIPboundtypeOpposite
#undef SCIProwGetNNonz
#undef SCIProwGetNLPNonz
#undef SCIProwGetCols
#undef SCIProwGetVals
#undef SCIProwGetConstant
#undef SCIProwGetNorm
#undef SCIProwGetSumNorm
#undef SCIProwGetLhs
#undef SCIProwGetRhs
#undef SCIProwGetDualsol
#undef SCIProwGetDualfarkas
#undef SCIProwGetBasisStatus
#undef SCIProwGetName
#undef SCIProwGetIndex
#undef SCIProwIsIntegral
#undef SCIProwIsLocal
#undef SCIProwIsModifiable
#undef SCIProwIsRemoveable
#undef SCIProwGetLPPos
#undef SCIProwGetLPDepth
#undef SCIProwIsInLP
#undef SCIPlpGetCols
#undef SCIPlpGetNCols
#undef SCIPlpGetRows
#undef SCIPlpGetNRows
#undef SCIPlpGetNewcols
#undef SCIPlpGetNNewcols
#undef SCIPlpGetNewrows
#undef SCIPlpGetNNewrows
#undef SCIPlpGetObjNorm
#undef SCIPlpGetLPI
#undef SCIPlpIsSolved
#undef SCIPlpIsSolBasic
#undef SCIPlpDiving
#undef SCIPlpDivingObjChanged
#undef SCIPlpMarkDivingObjChanged

/** gets objective value of column */
SCIP_Real SCIPcolGetObj(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->obj;
}

/** gets lower bound of column */
SCIP_Real SCIPcolGetLb(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->lb;
}

/** gets upper bound of column */
SCIP_Real SCIPcolGetUb(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->ub;
}

/** gets best bound of column with respect to the objective function */
SCIP_Real SCIPcolGetBestBound(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   if( col->obj >= 0.0 )
      return col->lb;
   else
      return col->ub;
}

/** gets the primal LP solution of a column */
SCIP_Real SCIPcolGetPrimsol(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   if( col->lppos >= 0 )
      return col->primsol;
   else
      return 0.0;
}

/** gets the basis status of a column in the LP solution; only valid for LPs with status SCIP_LPSOLSTAT_OPTIMAL
 *  and with SCIPisLPSolBasic(scip) == TRUE; returns SCIP_BASESTAT_ZERO for columns not in the current SCIP_LP
 */
SCIP_BASESTAT SCIPcolGetBasisStatus(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert(col->lppos >= 0 || (SCIP_BASESTAT)col->basisstatus == SCIP_BASESTAT_ZERO);

   return (SCIP_BASESTAT)col->basisstatus;
}

/** gets variable this column represents */
SCIP_VAR* SCIPcolGetVar(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->var;
}

/** returns whether the associated variable is of integral type (binary, integer, implicit integer) */
SCIP_Bool SCIPcolIsIntegral(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert(SCIPvarIsIntegral(col->var) == col->integral);

   return col->integral;
}

/** returns TRUE iff column is removeable from the LP (due to aging or cleanup) */
SCIP_Bool SCIPcolIsRemoveable(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->removeable;
}

/** gets position of column in current LP, or -1 if it is not in LP */
int SCIPcolGetLPPos(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert((col->lppos == -1) == (col->lpdepth == -1));

   return col->lppos;
}

/** gets depth in the tree where the column entered the LP, or -1 if it is not in LP */
int SCIPcolGetLPDepth(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert((col->lppos == -1) == (col->lpdepth == -1));

   return col->lpdepth;
}

/** returns TRUE iff column is member of current LP */
SCIP_Bool SCIPcolIsInLP(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert((col->lppos == -1) == (col->lpdepth == -1));

   return (col->lppos >= 0);
}

/** get number of nonzero entries in column vector */
int SCIPcolGetNNonz(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->len;
}

/** get number of nonzero entries in column vector, that correspond to rows currently in the SCIP_LP;
 *  Warning! This method is only applicable on columns, that are completely linked to their rows (e.g. a column
 *  that is in the current LP and the LP was solved, or a column that was in a solved LP and didn't change afterwards
 */
int SCIPcolGetNLPNonz(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert(col->nunlinked == 0);

   return col->nlprows;
}

/** gets array with rows of nonzero entries */
SCIP_ROW** SCIPcolGetRows(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->rows;
}

/** gets array with coefficients of nonzero entries */
SCIP_Real* SCIPcolGetVals(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->vals;
}

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given column, or -1 if strong branching was never applied to the column in current run
 */
SCIP_Longint SCIPcolGetStrongbranchNode(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->sbnode;
}

/** gets number of times, strong branching was applied in current run on the given column */
int SCIPcolGetNStrongbranchs(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->nsbcalls;
}

/** gets opposite bound type of given bound type */
SCIP_BOUNDTYPE SCIPboundtypeOpposite(
   SCIP_BOUNDTYPE        boundtype           /**< type of bound (lower or upper) */
   )
{
   assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

   return (boundtype == SCIP_BOUNDTYPE_LOWER ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
}

/** get number of nonzero entries in row vector */
int SCIProwGetNNonz(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->len;
}

/** get number of nonzero entries in row vector, that correspond to columns currently in the SCIP_LP;
 *  Warning! This method is only applicable on rows, that are completely linked to their columns (e.g. a row
 *  that is in the current LP and the LP was solved, or a row that was in a solved LP and didn't change afterwards
 */
int SCIProwGetNLPNonz(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->nunlinked == 0);

   return row->nlpcols;
}

/** gets array with columns of nonzero entries */
SCIP_COL** SCIProwGetCols(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->cols;
}

/** gets array with coefficients of nonzero entries */
SCIP_Real* SCIProwGetVals(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->vals;
}

/** gets constant shift of row */
SCIP_Real SCIProwGetConstant(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->constant;
}

/** gets euclidean norm of row vector */
SCIP_Real SCIProwGetNorm(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return sqrt(row->sqrnorm);
}

/** gets sum norm of row vector (sum of absolute values of coefficients) */
SCIP_Real SCIProwGetSumNorm(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->sumnorm;
}

/** returns the left hand side of the row */
SCIP_Real SCIProwGetLhs(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->lhs;
}

/** returns the right hand side of the row */
SCIP_Real SCIProwGetRhs(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->rhs;
}

/** gets the dual LP solution of a row */
SCIP_Real SCIProwGetDualsol(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   if( row->lppos >= 0 )
      return row->dualsol;
   else
      return 0.0;
}

/** gets the dual farkas coefficient of a row in an infeasible LP */
SCIP_Real SCIProwGetDualfarkas(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   if( row->lppos >= 0 )
      return row->dualfarkas;
   else
      return 0.0;
}

/** gets the basis status of a row in the LP solution; only valid for LPs with status SCIP_LPSOLSTAT_OPTIMAL
 *  and with SCIPisLPSolBasic(scip) == TRUE; returns SCIP_BASESTAT_BASIC for rows not in the current SCIP_LP
 */
SCIP_BASESTAT SCIProwGetBasisStatus(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->lppos >= 0 || (SCIP_BASESTAT)row->basisstatus == SCIP_BASESTAT_BASIC);

   return (SCIP_BASESTAT)row->basisstatus;
}

/** returns the name of the row */
const char* SCIProwGetName(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->name;
}

/** gets unique index of row */
int SCIProwGetIndex(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->index;
}

/** returns TRUE iff the activity of the row (without the row's constant) is always integral in a feasible solution */
SCIP_Bool SCIProwIsIntegral(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->integral;
}

/** returns TRUE iff row is only valid locally */
SCIP_Bool SCIProwIsLocal(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->local;
}

/** returns TRUE iff row is modifiable during node processing (subject to column generation) */
SCIP_Bool SCIProwIsModifiable(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->modifiable;
}

/** returns TRUE iff row is removeable from the LP (due to aging or cleanup) */
SCIP_Bool SCIProwIsRemoveable(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->removeable;
}

/** gets position of row in current LP, or -1 if it is not in LP */
int SCIProwGetLPPos(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert((row->lppos == -1) == (row->lpdepth == -1));

   return row->lppos;
}

/** gets depth in the tree where the row entered the LP, or -1 if it is not in LP */
int SCIProwGetLPDepth(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert((row->lppos == -1) == (row->lpdepth == -1));

   return row->lpdepth;
}

/** returns TRUE iff row is member of current LP */
SCIP_Bool SCIProwIsInLP(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert((row->lppos == -1) == (row->lpdepth == -1));

   return (row->lppos >= 0);
}

/** gets array with columns of the LP */
SCIP_COL** SCIPlpGetCols(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->cols;
}

/** gets current number of columns in LP */
int SCIPlpGetNCols(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->ncols;
}

/** gets array with rows of the LP */
SCIP_ROW** SCIPlpGetRows(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->rows;
}

/** gets current number of rows in LP */
int SCIPlpGetNRows(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->nrows;
}

/** gets array with newly added columns after the last mark */
SCIP_COL** SCIPlpGetNewcols(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewcol && lp->firstnewcol <= lp->ncols);

   return &(lp->cols[lp->firstnewcol]);
}

/** gets number of newly added columns after the last mark */
int SCIPlpGetNNewcols(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewcol && lp->firstnewcol <= lp->ncols);

   return lp->ncols - lp->firstnewcol;
}

/** gets array with newly added rows after the last mark */
SCIP_ROW** SCIPlpGetNewrows(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewrow && lp->firstnewrow <= lp->nrows);

   return &(lp->rows[lp->firstnewrow]);
}

/** gets number of newly added rows after the last mark */
int SCIPlpGetNNewrows(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewrow && lp->firstnewrow <= lp->nrows);

   return lp->nrows - lp->firstnewrow;
}

/** gets euclidean norm of objective function vector of column variables */
SCIP_Real SCIPlpGetObjNorm(
   SCIP_LP*              lp                  /**< LP data */
   )
{
   assert(lp != NULL);

   return SQRT(lp->objsqrnorm);
}

/** gets the LP solver interface */
SCIP_LPI* SCIPlpGetLPI(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->lpi;
}

/** returns whether the current LP is flushed and solved */
SCIP_Bool SCIPlpIsSolved(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->flushed && lp->solved;
}

/** returns whether the current LP solution is a basic solution */
SCIP_Bool SCIPlpIsSolBasic(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->solisbasic;
}

/** returns whether the LP is in diving mode */
SCIP_Bool SCIPlpDiving(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->diving;
}

/** returns whether the LP is in diving mode and the objective value of at least one column was changed */
SCIP_Bool SCIPlpDivingObjChanged(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->divingobjchg;
}

/** marks the diving LP to have a changed objective function */
void SCIPlpMarkDivingObjChanged(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->diving);

   lp->divingobjchg = TRUE;
}
