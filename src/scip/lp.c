/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: lp.c,v 1.161 2004/11/03 13:26:41 bzfwolte Exp $"

/**@file   lp.c
 * @brief  LP management methods and datastructures
 * @author Tobias Achterberg
 * @author Kati Wolter
 *
 *  In LP management, we have to differ between the current LP and the LP
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

#include "def.h"
#include "message.h"
#include "set.h"
#include "stat.h"
#include "intervalarith.h"
#include "clock.h"
#include "lpi.h"
#include "misc.h"
#include "lp.h"
#include "var.h"
#include "prob.h"
#include "sol.h"


/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that chgcols array can store at least num entries */
static
RETCODE ensureChgcolsSize(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nchgcols <= lp->chgcolssize);
   
   if( num > lp->chgcolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&lp->chgcols, newsize) );
      lp->chgcolssize = newsize;
   }
   assert(num <= lp->chgcolssize);

   return SCIP_OKAY;
}

/** ensures, that chgrows array can store at least num entries */
static
RETCODE ensureChgrowsSize(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nchgrows <= lp->chgrowssize);
   
   if( num > lp->chgrowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&lp->chgrows, newsize) );
      lp->chgrowssize = newsize;
   }
   assert(num <= lp->chgrowssize);

   return SCIP_OKAY;
}

/** ensures, that lpicols array can store at least num entries */
static
RETCODE ensureLpicolsSize(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nlpicols <= lp->lpicolssize);
   
   if( num > lp->lpicolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&lp->lpicols, newsize) );
      lp->lpicolssize = newsize;
   }
   assert(num <= lp->lpicolssize);

   return SCIP_OKAY;
}

/** ensures, that lpirows array can store at least num entries */
static
RETCODE ensureLpirowsSize(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nlpirows <= lp->lpirowssize);
   
   if( num > lp->lpirowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&lp->lpirows, newsize) );
      lp->lpirowssize = newsize;
   }
   assert(num <= lp->lpirowssize);

   return SCIP_OKAY;
}

/** ensures, that cols array can store at least num entries */
static
RETCODE ensureColsSize(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->ncols <= lp->colssize);
   
   if( num > lp->colssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&lp->cols, newsize) );
      lp->colssize = newsize;
   }
   assert(num <= lp->colssize);

   return SCIP_OKAY;
}

/** ensures, that rows array can store at least num entries */
static
RETCODE ensureRowsSize(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nrows <= lp->rowssize);
   
   if( num > lp->rowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&lp->rows, newsize) );
      lp->rowssize = newsize;
   }
   assert(num <= lp->rowssize);

   return SCIP_OKAY;
}

/** ensures, that row array of column can store at least num entries */
static
RETCODE colEnsureSize(
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(col != NULL);
   assert(col->len <= col->size);
   
   if( num > col->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &col->rows, col->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &col->vals, col->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &col->linkpos, col->size, newsize) );
      col->size = newsize;
   }
   assert(num <= col->size);

   return SCIP_OKAY;
}

/** ensures, that column array of row can store at least num entries */
RETCODE SCIProwEnsureSize(
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(row != NULL);
   assert(row->len <= row->size);
   
   if( num > row->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &row->cols, row->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &row->cols_index, row->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &row->vals, row->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &row->linkpos, row->size, newsize) );
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
   COL*             col,                /**< LP column */
   int              firstpos,           /**< first position to include in the sort */
   int              lastpos             /**< last position to include in the sort */
   )
{
   ROW** rows;
   Real* vals;
   int* linkpos;
   ROW* tmprow;
   Real tmpval;
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
   ROW*             row,                /**< LP row */
   int              firstpos,           /**< first position to include in the sort */
   int              lastpos             /**< last position to include in the sort */
   )
{
   COL** cols;
   Real* vals;
   int* index;
   int* linkpos;
   COL* tmpcol;
   Real tmpval;
   int tmpindex;
   int tmplinkpos;
   int pos;
   int sortpos;

   assert(row != NULL);
   assert(0 <= firstpos && firstpos <= lastpos+1 && lastpos < row->len);

   /**@todo do a quick sort here, if many elements are unsorted (sorted-Bool -> sorted-Int?) */
   cols = row->cols;
   vals = row->vals;
   index = row->cols_index;
   linkpos = row->linkpos;

#ifndef NDEBUG
   for( pos = 0; pos < row->len; ++pos )
      assert(index[pos] == cols[pos]->index);
#endif

   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && index[pos] <= index[pos+1] )
            pos++;
         if( pos >= lastpos )
            break;
         assert(index[pos] > index[pos+1]);
         tmpcol = cols[pos];
         tmpindex = index[pos];
         tmpval = vals[pos];
         tmplinkpos = linkpos[pos];
         do
         {
            cols[pos] = cols[pos+1];
            index[pos] = index[pos+1];
            vals[pos] = vals[pos+1];
            linkpos[pos] = linkpos[pos+1];
            pos++;
         }
         while( pos < lastpos && index[pos+1] < tmpindex );
         cols[pos] = tmpcol;
         index[pos] = tmpindex;
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
         while( pos > firstpos && index[pos-1] <= index[pos] )
            pos--;
         if( pos <= firstpos )
            break;
         assert(index[pos-1] > index[pos]);
         tmpcol = cols[pos];
         tmpindex = index[pos];
         tmpval = vals[pos];
         tmplinkpos = linkpos[pos];
         do
         {
            cols[pos] = cols[pos-1];
            index[pos] = index[pos-1];
            vals[pos] = vals[pos-1];
            linkpos[pos] = linkpos[pos-1];
            pos--;
         }
         while( pos > firstpos && index[pos-1] > tmpindex );
         cols[pos] = tmpcol;
         index[pos] = tmpindex;
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
   COL*             col                 /**< column to be sorted */
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
   COL*             col                 /**< column to be sorted */
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
   ROW*             row                 /**< row to be sorted */
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
   ROW*             row                 /**< row to be sorted */
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
   COL*             col,                /**< column to be searched in */
   const ROW*       row,                /**< coefficient to be searched for */
   int              minpos,             /**< first position of search range */
   int              maxpos              /**< last position of search range */
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
   COL*             col,                /**< column to be searched in */
   const ROW*       row                 /**< coefficient to be searched for */
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
   ROW*             row,                /**< row to be searched in */
   const COL*       col,                /**< coefficient to be searched for */
   int              minpos,             /**< first position of search range */
   int              maxpos              /**< last position of search range */
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
   ROW*             row,                /**< row to be searched in */
   const COL*       col                 /**< coefficient to be searched for */
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
   COL*             col,                /**< LP column */
   int              oldpos,             /**< old position of coefficient */
   int              newpos              /**< new position of coefficient */
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
   COL*             col,                /**< LP column */
   int              pos1,               /**< position of first coefficient */
   int              pos2                /**< position of second coefficient */
   )
{
   ROW* tmprow;
   Real tmpval;
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
   ROW*             row,                /**< LP row */
   int              oldpos,             /**< old position of coefficient */
   int              newpos              /**< new position of coefficient */
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
   ROW*             row,                /**< LP row */
   int              pos1,               /**< position of first coefficient */
   int              pos2                /**< position of second coefficient */
   )
{
   COL* tmpcol;
   Real tmpval;
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
#define ASSERT(x) do { if( !(x) ) abort(); } while( FALSE )
#else
#define ASSERT(x) assert(x)
#endif

static Bool msgdisp = FALSE;

static
void checkLinks(
   LP*              lp                  /**< current LP data */
   )
{
   COL* col;
   ROW* row;
   int i;
   int j;

   ASSERT(lp != NULL);

   if( !msgdisp )
   {
      warningMessage("LP LINK CHECKING ACTIVATED! THIS IS VERY SLOW!\n");
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
   ROW*             row,                /**< LP row */
   COL*             col,                /**< LP col */
   LP*              lp                  /**< current LP data */
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
RETCODE rowAddCoef(
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   COL*             col,                /**< LP column */
   Real             val,                /**< value of coefficient */
   int              linkpos             /**< position of row in the column's row array, or -1 */
   );

/** adds a previously non existing coefficient to an LP column */
static
RETCODE colAddCoef(
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   ROW*             row,                /**< LP row */
   Real             val,                /**< value of coefficient */
   int              linkpos             /**< position of column in the row's col array, or -1 */
   )
{
   int pos;

   assert(memhdr != NULL);
   assert(col != NULL);
   assert(col->nlprows <= col->len);
   assert(col->var != NULL);
   assert(row != NULL);
   assert(!SCIPsetIsZero(set, val));
   /*assert(colSearchCoef(col, row) == -1);*/ /* this assert would lead to slight differences in the solution process */

   CHECK_OKAY( colEnsureSize(col, memhdr, set, col->len+1) );
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
         CHECK_OKAY( rowAddCoef(row, memhdr, set, lp, col, val, pos) );
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

   debugMessage("added coefficient %g * <%s> at position %d (%d/%d) to column <%s> (nunlinked=%d)\n",
      val, row->name, pos, col->nlprows, col->len, SCIPvarGetName(col->var), col->nunlinked);

   return SCIP_OKAY;
}

/** deletes coefficient at given position from column */
static
RETCODE colDelCoefPos(
   COL*             col,                /**< column to be changed */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   int              pos                 /**< position in column vector to delete */
   )
{
   ROW* row;

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
RETCODE colChgCoefPos(
   COL*             col,                /**< LP column */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   int              pos,                /**< position in column vector to change */
   Real             val                 /**< value of coefficient */
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
      CHECK_OKAY( colDelCoefPos(col, set, lp, pos) );
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
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   COL*             col,                /**< column of added coefficient */
   Real             val                 /**< value of added coefficient */
   )
{
   Real absval;

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
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   COL*             col,                /**< column of deleted coefficient */
   Real             val,                /**< value of deleted coefficient */
   Bool             updateindex         /**< should the minimal/maximal column index of row be updated? */
   )
{
   Real absval;

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
RETCODE rowAddCoef(
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   COL*             col,                /**< LP column */
   Real             val,                /**< value of coefficient */
   int              linkpos             /**< position of row in the column's row array, or -1 */
   )
{
   int pos;

   assert(row != NULL);
   assert(row->nlpcols <= row->len);
   assert(memhdr != NULL);
   assert(col != NULL);
   assert(col->var != NULL);
   assert(col->var_probindex == SCIPvarGetProbindex(col->var));
   assert(!SCIPsetIsZero(set, val));
   /*assert(rowSearchCoef(row, col) == -1);*/ /* this assert would lead to slight differences in the solution process */

   if( row->nlocks > 0 )
   {
      errorMessage("cannot add a coefficient to the locked unmodifiable row <%s>\n", row->name);
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIProwEnsureSize(row, memhdr, set, row->len+1) );
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
         CHECK_OKAY( colAddCoef(col, memhdr, set, lp, row, val, pos) );
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

   debugMessage("added coefficient %g * <%s> at position %d (%d/%d) to row <%s> (nunlinked=%d)\n",
      val, SCIPvarGetName(col->var), pos, row->nlpcols, row->len, row->name, row->nunlinked);

   return SCIP_OKAY;
}

/** deletes coefficient at given position from row */
static
RETCODE rowDelCoefPos(
   ROW*             row,                /**< row to be changed */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   int              pos                 /**< position in row vector to delete */
   )
{
   COL* col;
   Real val;

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
      errorMessage("cannot delete a coefficient from the locked unmodifiable row <%s>\n", row->name);
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
RETCODE rowChgCoefPos(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   int              pos,                /**< position in row vector to change */
   Real             val                 /**< value of coefficient */
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
      errorMessage("cannot change a coefficient of the locked unmodifiable row <%s>\n", row->name);
      return SCIP_INVALIDDATA;
   }

   if( SCIPsetIsZero(set, val) )
   {
      /* delete existing coefficient */
      CHECK_OKAY( rowDelCoefPos(row, set, lp, pos) );
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
RETCODE rowSideChanged(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   SIDETYPE         sidetype            /**< type of side: left or right hand side */
   )
{
   assert(row != NULL);
   assert(lp != NULL);

   if( row->lpipos >= 0 )
   {
      /* insert row in the chgrows list (if not already there) */
      if( !row->lhschanged && !row->rhschanged )
      {
         CHECK_OKAY( ensureChgrowsSize(lp, set, lp->nchgrows+1) );
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
         errorMessage("unknown row side type\n");
         abort();
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
RETCODE colLink(
   COL*             col,                /**< column data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(col != NULL);
   assert(col->var != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( col->nunlinked > 0 )
   {
      debugMessage("linking column <%s>\n", SCIPvarGetName(col->var));

      /* unlinked rows can only be in the non-LP/unlinked rows part of the rows array */
      for( i = col->nlprows; i < col->len; ++i )
      {
         assert(!SCIPsetIsZero(set, col->vals[i]));
         if( col->linkpos[i] == -1 )
         {
            /* this call might swap the current row with the first non-LP/not linked row, but this is of no harm */
            CHECK_OKAY( rowAddCoef(col->rows[i], memhdr, set, lp, col, col->vals[i], i) );
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
RETCODE colUnlink(
   COL*             col,                /**< column data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(col != NULL);
   assert(col->var != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( col->nunlinked < col->len )
   {
      debugMessage("unlinking column <%s>\n", SCIPvarGetName(col->var));
      for( i = 0; i < col->len; ++i )
      {
         if( col->linkpos[i] >= 0 )
         {
            assert(col->rows[i]->cols[col->linkpos[i]] == col);
            CHECK_OKAY( rowDelCoefPos(col->rows[i], set, lp, col->linkpos[i]) );
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
RETCODE rowLink(
   ROW*             row,                /**< row data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(row != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( row->nunlinked > 0 )
   {
      debugMessage("linking row <%s>\n", row->name);

      /* unlinked columns can only be in the non-LP/unlinked columns part of the cols array */
      for( i = row->nlpcols; i < row->len; ++i )
      {
         assert(!SCIPsetIsZero(set, row->vals[i]));
         if( row->linkpos[i] == -1 )
         {
            /* this call might swap the current column with the first non-LP/not linked column, but this is of no harm */
            CHECK_OKAY( colAddCoef(row->cols[i], memhdr, set, lp, row, row->vals[i], i) );
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
RETCODE rowUnlink(
   ROW*             row,                /**< row data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(row != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( row->nunlinked < row->len )
   {
      debugMessage("unlinking row <%s>\n", row->name);
      for( i = 0; i < row->len; ++i )
      {
         if( row->linkpos[i] >= 0 )
         {
            assert(row->cols[i]->rows[row->linkpos[i]] == row);
            CHECK_OKAY( colDelCoefPos(row->cols[i], set, lp, row->linkpos[i]) );
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
RETCODE lpSetIntpar(
   LP*              lp,                 /**< current LP data */
   LPPARAM          lpparam,            /**< LP parameter */
   int              value,              /**< value to set parameter to */
   Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   RETCODE retcode;

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

/** sets parameter of type Real in LP solver, ignoring unknown parameters */
static
RETCODE lpSetRealpar(
   LP*              lp,                 /**< current LP data */
   LPPARAM          lpparam,            /**< LP parameter */
   Real             value,              /**< value to set parameter to */
   Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   RETCODE retcode;

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
RETCODE lpCheckIntpar(
   LP*              lp,                 /**< current LP data */
   LPPARAM          lpparam,            /**< LP parameter */
   int              value               /**< value parameter should have */
   )
{
   RETCODE retcode;
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

/** checks, that parameter of type Real in LP solver has the given value, ignoring unknown parameters */
static
RETCODE lpCheckRealpar(
   LP*              lp,                 /**< current LP data */
   LPPARAM          lpparam,            /**< LP parameter */
   Real             value               /**< value parameter should have */
   )
{
   RETCODE retcode;
   Real lpivalue;

   assert(lp != NULL);

   retcode = SCIPlpiGetRealpar(lp->lpi, lpparam, &lpivalue);

   /* ignore unknown parameter error */
   if( retcode == SCIP_PARAMETERUNKNOWN )
      return SCIP_OKAY;

   /* check value */
   if( lpivalue != value )
      return SCIP_LPERROR;

   return retcode;
}
#else
#define lpCheckIntpar(lp, lpparam, value) SCIP_OKAY
#define lpCheckRealpar(lp, lpparam, value) SCIP_OKAY
#endif

/** sets the upper objective limit of the LP solver */
static
RETCODE lpSetUobjlim(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   Real             uobjlim             /**< new feasibility tolerance */
   )
{
   assert(lp != NULL);
   assert(set != NULL);

   /* if we want so solve exactly, we cannot rely on the LP solver's objective limit handling */
   if( set->misc_exactsolve )
      return SCIP_OKAY;

   CHECK_OKAY( lpCheckRealpar(lp, SCIP_LPPAR_UOBJLIM, lp->lpiuobjlim) );

   if( uobjlim != lp->lpiuobjlim )
   {
      Bool success;

      CHECK_OKAY( lpSetRealpar(lp, SCIP_LPPAR_UOBJLIM, uobjlim, &success) );
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
RETCODE lpSetFeastol(
   LP*              lp,                 /**< current LP data */
   Real             feastol,            /**< new feasibility tolerance */
   Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(feastol >= 0.0);
   assert(success != NULL);

   CHECK_OKAY( lpCheckRealpar(lp, SCIP_LPPAR_FEASTOL, lp->lpifeastol) );

   if( feastol != lp->lpifeastol )
   {
      CHECK_OKAY( lpSetRealpar(lp, SCIP_LPPAR_FEASTOL, feastol, success) );
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
RETCODE lpSetDualFeastol(
   LP*              lp,                 /**< current LP data */
   Real             dualfeastol,        /**< new reduced costs feasibility tolerance */
   Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(dualfeastol >= 0.0);
   assert(success != NULL);

   CHECK_OKAY( lpCheckRealpar(lp, SCIP_LPPAR_DUALFEASTOL, lp->lpidualfeastol) );

   if( dualfeastol != lp->lpidualfeastol )
   {
      CHECK_OKAY( lpSetRealpar(lp, SCIP_LPPAR_DUALFEASTOL, dualfeastol, success) );
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

/** sets the FROMSCRATCH setting of the LP solver */
static
RETCODE lpSetFromscratch(
   LP*              lp,                 /**< current LP data */
   Bool             fromscratch,        /**< new FROMSCRATCH setting */
   Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   CHECK_OKAY( lpCheckIntpar(lp, SCIP_LPPAR_FROMSCRATCH, lp->lpifromscratch) );

   if( fromscratch != lp->lpifromscratch )
   {
      CHECK_OKAY( lpSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, fromscratch, success) );
      if( *success )
         lp->lpifromscratch = fromscratch;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the FASTMIP setting of the LP solver */
static
RETCODE lpSetFastmip(
   LP*              lp,                 /**< current LP data */
   Bool             fastmip,            /**< new FASTMIP setting */
   Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   CHECK_OKAY( lpCheckIntpar(lp, SCIP_LPPAR_FASTMIP, lp->lpifastmip) );

   if( fastmip != lp->lpifastmip )
   {
      CHECK_OKAY( lpSetIntpar(lp, SCIP_LPPAR_FASTMIP, fastmip, success) );
      if( *success )
         lp->lpifastmip = fastmip;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the SCALING setting of the LP solver */
static
RETCODE lpSetScaling(
   LP*              lp,                 /**< current LP data */
   Bool             scaling,            /**< new SCALING setting */
   Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   CHECK_OKAY( lpCheckIntpar(lp, SCIP_LPPAR_SCALING, lp->lpiscaling) );

   if( scaling != lp->lpiscaling )
   {
      CHECK_OKAY( lpSetIntpar(lp, SCIP_LPPAR_SCALING, scaling, success) );
      if( *success )
         lp->lpiscaling = scaling;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the PRESOLVING setting of the LP solver */
static
RETCODE lpSetPresolving(
   LP*              lp,                 /**< current LP data */
   Bool             presolving,         /**< new PRESOLVING setting */
   Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   CHECK_OKAY( lpCheckIntpar(lp, SCIP_LPPAR_PRESOLVING, lp->lpipresolving) );

   if( presolving != lp->lpipresolving )
   {
      CHECK_OKAY( lpSetIntpar(lp, SCIP_LPPAR_PRESOLVING, presolving, success) );
      if( *success )
         lp->lpipresolving = presolving;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the iteration limit of the LP solver */
static
RETCODE lpSetIterationLimit(
   LP*              lp,                 /**< current LP data */
   int              itlim               /**< maximal number of LP iterations to perform, or -1 for no limit */
   )
{
   Bool success;

   assert(lp != NULL);
   assert(itlim >= -1);

   if( itlim == -1 )
      itlim = INT_MAX;

   CHECK_OKAY( lpCheckIntpar(lp, SCIP_LPPAR_LPITLIM, lp->lpiitlim) );

   if( itlim != lp->lpiitlim )
   {
      CHECK_OKAY( lpSetIntpar(lp, SCIP_LPPAR_LPITLIM, itlim, &success) );
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

/** sets the verbosity of the LP solver */
static
RETCODE lpSetLPInfo(
   LP*              lp,                 /**< current LP data */
   Bool             lpinfo              /**< should the LP solver display status messages? */
   )
{
   Bool success;

   assert(lp != NULL);

   CHECK_OKAY( lpCheckIntpar(lp, SCIP_LPPAR_LPINFO, lp->lpilpinfo) );

   if( lpinfo != lp->lpilpinfo )
   {
      CHECK_OKAY( lpSetIntpar(lp, SCIP_LPPAR_LPINFO, lpinfo, &success) );
      if( success )
         lp->lpilpinfo = lpinfo;
   }

   return SCIP_OKAY;
}




/*
 * Column methods
 */

/** creates an LP column */
RETCODE SCIPcolCreate(
   COL**            col,                /**< pointer to column data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR*             var,                /**< variable, this column represents */
   int              len,                /**< number of nonzeros in the column */
   ROW**            rows,               /**< array with rows of column entries */
   Real*            vals,               /**< array with coefficients of column entries */
   Bool             removeable          /**< should the column be removed from the LP due to aging or cleanup? */
   )
{
   int i;

   assert(col != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(var != NULL);
   assert(len >= 0);
   assert(len == 0 || (rows != NULL && vals != NULL));

   ALLOC_OKAY( allocBlockMemory(memhdr, col) );

   if( len > 0 )
   {
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*col)->rows, rows, len) );
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*col)->vals, vals, len) );
      ALLOC_OKAY( allocBlockMemoryArray(memhdr, &(*col)->linkpos, len) );

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
   (*col)->strongbranchdown = SCIP_INVALID;
   (*col)->strongbranchup = SCIP_INVALID;
   (*col)->strongbranchsolval  = SCIP_INVALID;
   (*col)->strongbranchlpobjval = SCIP_INVALID;
   (*col)->strongbranchnode = -1;
   (*col)->validredcostlp = -1;
   (*col)->validfarkaslp = -1;
   (*col)->validstrongbranchlp = -1;
   (*col)->strongbranchitlim = -1;
   (*col)->age = 0;
   (*col)->obsoletenode = -1;
   (*col)->var_probindex = SCIPvarGetProbindex(var);
   (*col)->lprowssorted = TRUE;
   (*col)->nonlprowssorted = (len <= 1);
   (*col)->objchanged = FALSE;
   (*col)->lbchanged = FALSE;
   (*col)->ubchanged = FALSE;
   (*col)->coefchanged = FALSE;
   (*col)->integral = SCIPvarIsIntegral(var);
   (*col)->removeable = removeable;

   return SCIP_OKAY;
}

/** frees an LP column */
RETCODE SCIPcolFree(
   COL**            col,                /**< pointer to LP column */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   assert(memhdr != NULL);
   assert(col != NULL);
   assert(*col != NULL);
   assert((*col)->var != NULL);
   assert(SCIPvarGetStatus((*col)->var) == SCIP_VARSTATUS_COLUMN);
   assert(&(*col)->var->data.col == col); /* SCIPcolFree() has to be called from SCIPvarFree() */
   assert((*col)->lppos == -1);
   assert((*col)->lpipos == -1);

   /* remove column indices from corresponding rows */
   CHECK_OKAY( colUnlink(*col, memhdr, set, lp) );

   freeBlockMemoryArrayNull(memhdr, &(*col)->rows, (*col)->size);
   freeBlockMemoryArrayNull(memhdr, &(*col)->vals, (*col)->size);
   freeBlockMemoryArrayNull(memhdr, &(*col)->linkpos, (*col)->size);
   freeBlockMemory(memhdr, col);

   return SCIP_OKAY;
}

/** sorts column entries such that LP rows precede non-LP rows and inside both parts lower row indices precede higher ones
 */
void SCIPcolSort(
   COL*             col                 /**< column to be sorted */
   )
{
   /* sort LP rows */
   colSortLP(col);

   /* sort non-LP rows */
   colSortNonLP(col);
}

/** adds a previously non existing coefficient to an LP column */
RETCODE SCIPcolAddCoef(
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);

   CHECK_OKAY( colAddCoef(col, memhdr, set, lp, row, val, -1) );

   checkLinks(lp);

   return SCIP_OKAY;
}

/** deletes existing coefficient from column */
RETCODE SCIPcolDelCoef(
   COL*             col,                /**< column to be changed */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   ROW*             row                 /**< coefficient to be deleted */
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
      errorMessage("coefficient for row <%s> doesn't exist in column <%s>\n", row->name, SCIPvarGetName(col->var));
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
      CHECK_OKAY( rowDelCoefPos(row, set, lp, col->linkpos[pos]) );
   }

   /* delete the row from the column's row vector */
   CHECK_OKAY( colDelCoefPos(col, set, lp, pos) );
   
   checkLinks(lp);

   return SCIP_OKAY;
}

/** changes or adds a coefficient to an LP column */
RETCODE SCIPcolChgCoef(
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
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
      CHECK_OKAY( colAddCoef(col, memhdr, set, lp, row, val, -1) );
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
         CHECK_OKAY( rowChgCoefPos(row, set, lp, col->linkpos[pos], val) );
      }

      /* change the coefficient in the column */
      CHECK_OKAY( colChgCoefPos(col, set, lp, pos, val) );
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** increases value of an existing or nonexisting coefficient in an LP column */
RETCODE SCIPcolIncCoef(
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   ROW*             row,                /**< LP row */
   Real             incval              /**< value to add to the coefficient */
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
      CHECK_OKAY( colAddCoef(col, memhdr, set, lp, row, incval, -1) );
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
         CHECK_OKAY( rowChgCoefPos(row, set, lp, col->linkpos[pos], col->vals[pos] + incval) );
      }

      /* change the coefficient in the column */
      CHECK_OKAY( colChgCoefPos(col, set, lp, pos, col->vals[pos] + incval) );
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** changes objective value of column */
RETCODE SCIPcolChgObj(
   COL*             col,                /**< LP column to change */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             newobj              /**< new objective value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(lp != NULL);
   
   debugMessage("changing objective value of column <%s> from %f to %f\n", SCIPvarGetName(col->var), col->obj, newobj);

   if( col->lpipos >= 0 && !SCIPsetIsEQ(set, col->obj, newobj) )
   {
      /* insert column in the chgcols list (if not already there) */
      if( !col->objchanged && !col->lbchanged && !col->ubchanged )
      {
         CHECK_OKAY( ensureChgcolsSize(lp, set, lp->nchgcols+1) );
         lp->chgcols[lp->nchgcols] = col;
         lp->nchgcols++;
      }
      
      /* mark objective value change in the column */
      col->objchanged = TRUE;
      
      /* mark the current LP unflushed */
      lp->flushed = FALSE;

      assert(lp->nchgcols > 0);
   }  

   /* update squared euclidean norm and sum norm of objective function vector */
   lp->objsqrnorm += SQR(newobj) - SQR(col->obj);
   lp->objsumnorm += REALABS(newobj) - REALABS(col->obj);

   /* store new objective function value */
   col->obj = newobj;

   return SCIP_OKAY;
}

/** changes lower bound of column */
RETCODE SCIPcolChgLb(
   COL*             col,                /**< LP column to change */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             newlb               /**< new lower bound value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(lp != NULL);
   
   debugMessage("changing lower bound of column <%s> from %f to %f\n", SCIPvarGetName(col->var), col->lb, newlb);

   if( col->lpipos >= 0 && !SCIPsetIsEQ(set, col->lb, newlb) )
   {
      /* insert column in the chgcols list (if not already there) */
      if( !col->objchanged && !col->lbchanged && !col->ubchanged )
      {
         CHECK_OKAY( ensureChgcolsSize(lp, set, lp->nchgcols+1) );
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
RETCODE SCIPcolChgUb(
   COL*             col,                /**< LP column to change */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             newub               /**< new upper bound value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(lp != NULL);
   
   debugMessage("changing upper bound of column <%s> from %f to %f\n", SCIPvarGetName(col->var), col->ub, newub);

   if( col->lpipos >= 0 && !SCIPsetIsEQ(set, col->ub, newub) )
   {
      /* insert column in the chgcols list (if not already there) */
      if( !col->objchanged && !col->lbchanged && !col->ubchanged )
      {
         CHECK_OKAY( ensureChgcolsSize(lp, set, lp->nchgcols+1) );
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
Real SCIPcolCalcRedcost(
   COL*             col,                /**< LP column */
   Real*            dualsol             /**< dual solution vector for current LP rows */
   )
{
   ROW* row;
   Real redcost;
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
Real colCalcInternalRedcost(
   COL*             col                 /**< LP column */
   )
{
   ROW* row;
   Real redcost;
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
Real SCIPcolGetRedcost(
   COL*             col,                /**< LP column */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< current LP data */
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
Real SCIPcolGetFeasibility(
   COL*             col,                /**< LP column */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< current LP data */
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
      Real redcost;

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
Real SCIPcolCalcFarkasCoef(
   COL*             col,                /**< LP column */
   Real*            dualfarkas          /**< dense dual farkas vector for current LP rows */
   )
{
   ROW* row;
   Real farkas;
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
Real colCalcInternalFarkasCoef(
   COL*             col                 /**< LP column */
   )
{
   ROW* row;
   Real farkas;
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
Real SCIPcolGetFarkasCoef(
   COL*             col,                /**< LP column */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< current LP data */
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
Real SCIPcolGetFarkasValue(
   COL*             col,                /**< LP column */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< current LP data */
   )
{
   Real farkascoef;

   assert(col != NULL);

   farkascoef = SCIPcolGetFarkasCoef(col, stat, lp);

   if( farkascoef > 0.0 )
      return col->ub * farkascoef;
   else
      return col->lb * farkascoef;
}

/** actually performs strong branching on the given variable */
static
RETCODE colStrongbranch(
   COL*             col,                /**< LP column */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   int              itlim,              /**< iteration limit for strong branchings */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   RETCODE retcode;
   Real strongbranchdown;
   Real strongbranchup;
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

   debugMessage("performing strong branching on variable <%s>(%g) with %d iterations\n", 
      SCIPvarGetName(col->var), col->primsol, itlim);

   /* start timing */
   SCIPclockStart(stat->strongbranchtime, set);
      
   /* call LPI strong branching */
   col->strongbranchitlim = itlim;
   retcode = SCIPlpiStrongbranch(lp->lpi, col->lpipos, col->primsol, itlim,
      &strongbranchdown, &strongbranchup, &iter);
   
   /* check return code for errors */
   if( retcode == SCIP_LPERROR )
   {
      *lperror = TRUE;
      col->strongbranchdown = SCIP_INVALID;
      col->strongbranchup = SCIP_INVALID;
      col->validstrongbranchlp = -1;
      col->strongbranchsolval = SCIP_INVALID;
      col->strongbranchlpobjval = SCIP_INVALID;
      col->strongbranchnode = -1;
   }
   else
   {
      *lperror = FALSE;
      CHECK_OKAY( retcode );

      col->strongbranchdown = MIN(strongbranchdown + lp->looseobjval, lp->cutoffbound);
      col->strongbranchup = MIN(strongbranchup + lp->looseobjval, lp->cutoffbound);
            
      /* update strong branching statistics */
      if( iter == -1 )
      {
         /* calculate avergate iteration number */
         iter = stat->nresolvelps > 0 ? (int)(2*stat->nresolvelpiterations / stat->nresolvelps)
            : stat->nlps > 0 ? (int)((stat->nlpiterations / stat->nlps) / 5) : 0;
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
RETCODE SCIPcolGetStrongbranch(
   COL*             col,                /**< LP column */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   LP*              lp,                 /**< LP data */
   int              itlim,              /**< iteration limit for strong branchings */
   Real*            down,               /**< stores dual bound after branching column down */
   Real*            up,                 /**< stores dual bound after branching column up */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
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

   if( col->validstrongbranchlp != stat->lpcount || itlim > col->strongbranchitlim )
   {
      col->validstrongbranchlp = stat->lpcount;
      col->strongbranchsolval = col->primsol;
      col->strongbranchlpobjval = SCIPlpGetObjval(lp, set);
      col->strongbranchnode = stat->nnodes;

      /* if a loose variables has an infinite best bound, the LP bound is -infinity and no gain can be achieved */
      if( lp->looseobjvalinf > 0 )
      {
         col->strongbranchdown = -SCIPsetInfinity(set);
         col->strongbranchup = -SCIPsetInfinity(set);
      }
      else
      {
         /* perform the strong branching on the column */
         CHECK_OKAY( colStrongbranch(col, set, stat, lp, itlim, lperror) );
      }
   }
   assert(*lperror || col->strongbranchdown < SCIP_INVALID);
   assert(*lperror || col->strongbranchup < SCIP_INVALID);

   *down = col->strongbranchdown;
   *up = col->strongbranchup;

   return SCIP_OKAY;
}

/** gets last strong branching information available for a column variable;
 *  returns values of SCIP_INVALID, if strong branching was not yet called on the given column;
 *  keep in mind, that the returned old values may have nothing to do with the current LP solution
 */
void SCIPcolGetStrongbranchLast(
   COL*             col,                /**< LP column */
   Real*            down,               /**< stores dual bound after branching column down, or NULL */
   Real*            up,                 /**< stores dual bound after branching column up, or NULL */
   Real*            solval,             /**< stores LP solution value of column at last strong branching call, or NULL */
   Real*            lpobjval            /**< stores LP objective value at last strong branching call, or NULL */
   )
{
   assert(col != NULL);

   if( down != NULL )
      *down = col->strongbranchdown;
   if( up != NULL )
      *up = col->strongbranchup;
   if( solval != NULL )
      *solval = col->strongbranchsolval;
   if( lpobjval != NULL )
      *lpobjval = col->strongbranchlpobjval;
}

/** output column to file stream */
void SCIPcolPrint(
   COL*             col,                /**< LP column */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int r;

   assert(col != NULL);
   assert(col->var != NULL);

   if( file == NULL )
      file = stdout;

   /* print bounds */
   fprintf(file, "(obj: %g) [%g,%g], ", col->obj, col->lb, col->ub);

   /* print coefficients */
   if( col->len == 0 )
      fprintf(file, "<empty>");
   for( r = 0; r < col->len; ++r )
   {
      assert(col->rows[r] != NULL);
      assert(col->rows[r]->name != NULL);
      fprintf(file, "%+g<%s> ", col->vals[r], col->rows[r]->name);
   }
   fprintf(file, "\n");
}


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets objective value of column */
Real SCIPcolGetObj(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->obj;
}

/** gets lower bound of column */
Real SCIPcolGetLb(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->lb;
}

/** gets upper bound of column */
Real SCIPcolGetUb(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->ub;
}

/** gets best bound of column with respect to the objective function */
Real SCIPcolGetBestBound(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   if( col->obj >= 0.0 )
      return col->lb;
   else
      return col->ub;
}

/** gets the primal LP solution of a column */
Real SCIPcolGetPrimsol(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   if( col->lppos >= 0 )
      return col->primsol;
   else
      return 0.0;
}

/** gets variable this column represents */
VAR* SCIPcolGetVar(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->var;
}

/** returns whether the associated variable is of integral type (binary, integer, implicit integer) */
Bool SCIPcolIsIntegral(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert(SCIPvarIsIntegral(col->var) == col->integral);

   return col->integral;
}

/** returns TRUE iff column is removeable from the LP (due to aging or cleanup) */
Bool SCIPcolIsRemoveable(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->removeable;
}

/** gets position of column in current LP, or -1 if it is not in LP */
int SCIPcolGetLPPos(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert((col->lppos == -1) == (col->lpdepth == -1));

   return col->lppos;
}

/** gets depth in the tree where the column entered the LP, or -1 if it is not in LP */
int SCIPcolGetLPDepth(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert((col->lppos == -1) == (col->lpdepth == -1));

   return col->lpdepth;
}

/** returns TRUE iff column is member of current LP */
Bool SCIPcolIsInLP(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert((col->lppos == -1) == (col->lpdepth == -1));

   return (col->lppos >= 0);
}

/** get number of nonzero entries in column vector */
int SCIPcolGetNNonz(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->len;
}

/** get number of nonzero entries in column vector, that correspond to rows currently in the LP;
 *  Warning! This method is only applicable on columns, that are completely linked to their rows (e.g. a column
 *  that is in the current LP and the LP was solved, or a column that was in a solved LP and didn't change afterwards
 */
int SCIPcolGetNLPNonz(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert(col->nunlinked == 0);

   return col->nlprows;
}

/** gets array with rows of nonzero entries */
ROW** SCIPcolGetRows(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->rows;
}

/** gets array with coefficients of nonzero entries */
Real* SCIPcolGetVals(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->vals;
}

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given column, or -1 if strong branching was never applied to the column in current run
 */
Longint SCIPcolGetStrongbranchNode(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->strongbranchnode;
}

/** gets opposite bound type of given bound type */
BOUNDTYPE SCIPboundtypeOpposite(
   BOUNDTYPE        boundtype           /**< type of bound (lower or upper) */
   )
{
   assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

   return (boundtype == SCIP_BOUNDTYPE_LOWER ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
}

#endif




/*
 * Row methods
 */

/** calculates row norms and min/maxidx from scratch, and checks for sortation */
static
void rowCalcNorms(
   ROW*             row,                /**< LP row */
   SET*             set                 /**< global SCIP settings */
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
Bool isIntegralScalar(
   Real             val,                /**< value that should be scaled to an integral value */
   Real             scalar,             /**< scalar that should be tried */
   Real             mindelta,           /**< minimal allowed difference s*c - i of scaled coefficient s*c and integral i */
   Real             maxdelta            /**< maximal allowed difference s*c - i of scaled coefficient s*c and integral i */
   )
{
   Real sval;
   Real ival;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   ival = EPSFLOOR(sval, -mindelta);
   
   return (sval - ival <= maxdelta);
}

/** scales row with given factor, and rounds coefficients to integers if close enough;
 *  the constant is automatically moved to the sides;
 *  if the row's activity is proven to be integral, the sides are automatically rounded to the next integer
 */
static
RETCODE rowScale(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   Real             scaleval,           /**< value to scale row with */
   Bool             integralcontvars,   /**< should the coefficients of the continuous variables also be made integral,
                                         *   if they are close to integral values? */
   Real             minrounddelta,      /**< minimal allowed difference s*c - i of scaled coefficient s*c and integral i,
                                         *   upto which the integral is used instead of the scaled real coefficient */
   Real             maxrounddelta       /**< maximal allowed difference s*c - i of scaled coefficient s*c and integral i
                                         *   upto which the integral is used instead of the scaled real coefficient */
   )
{
   COL* col;
   Real val;
   Real newval;
   Real intval;
   Real mindelta;
   Real maxdelta;
   Real absscaleval;
   Real lb;
   Real ub;
   Bool mindeltainf;
   Bool maxdeltainf;
   int pos;
   int c;

   assert(row != NULL);
   assert(row->len == 0 || row->cols != NULL);
   assert(row->len == 0 || row->vals != NULL);
   assert(SCIPsetIsPositive(set, scaleval));
   assert(-1.0 < minrounddelta && minrounddelta <= 0.0);
   assert(0.0 <= maxrounddelta && maxrounddelta < 1.0);

   debugMessage("scale row <%s> with %g (tolerance=[%g,%g])\n", row->name, scaleval, minrounddelta, maxrounddelta);

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
      if( ((integralcontvars || SCIPcolIsIntegral(col))
            && isIntegralScalar(val, scaleval, minrounddelta, maxrounddelta))
         || SCIPsetIsIntegral(set, newval) )
      {
         intval = EPSFLOOR(newval, -minrounddelta);
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
            CHECK_OKAY( colChgCoefPos(col, set, lp, row->linkpos[c], newval) );
         }
            
         /* change the coefficient in the row, and update the norms and integrality status */
         CHECK_OKAY( rowChgCoefPos(row, set, lp, c, newval) );
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
      CHECK_OKAY( SCIProwChgLhs(row, set, lp, newval) );
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
      CHECK_OKAY( SCIProwChgRhs(row, set, lp, newval) );
   }

   /* clear the row constant */
   CHECK_OKAY( SCIProwChgConstant(row, set, stat, lp, 0.0) );

   debug(SCIProwPrint(row, NULL));

   return SCIP_OKAY;
}

/** creates and captures an LP row */
RETCODE SCIProwCreate(
   ROW**            row,                /**< pointer to LP row data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            cols,               /**< array with columns of row entries */
   Real*            vals,               /**< array with coefficients of row entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             local,              /**< is row only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   assert(row != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);
   assert(len >= 0);
   assert(len == 0 || (cols != NULL && vals != NULL));
   assert(lhs <= rhs);

   ALLOC_OKAY( allocBlockMemory(memhdr, row) );

   (*row)->integral = TRUE;
   if( len > 0 )
   {
      VAR* var;
      int i;

      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*row)->cols, cols, len) );
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*row)->vals, vals, len) );
      ALLOC_OKAY( allocBlockMemoryArray(memhdr, &(*row)->cols_index, len) );
      ALLOC_OKAY( allocBlockMemoryArray(memhdr, &(*row)->linkpos, len) );

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
   
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*row)->name, name, strlen(name)+1) );
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
RETCODE SCIProwFree(
   ROW**            row,                /**< pointer to LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   assert(memhdr != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->nuses == 0);
   assert((*row)->lppos == -1);

   /* remove column indices from corresponding rows */
   CHECK_OKAY( rowUnlink(*row, memhdr, set, lp) );

   freeBlockMemoryArray(memhdr, &(*row)->name, strlen((*row)->name)+1);
   freeBlockMemoryArrayNull(memhdr, &(*row)->cols, (*row)->size);
   freeBlockMemoryArrayNull(memhdr, &(*row)->cols_index, (*row)->size);
   freeBlockMemoryArrayNull(memhdr, &(*row)->vals, (*row)->size);
   freeBlockMemoryArrayNull(memhdr, &(*row)->linkpos, (*row)->size);
   freeBlockMemory(memhdr, row);

   return SCIP_OKAY;
}

/** increases usage counter of LP row */
void SCIProwCapture(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->nuses >= 0);
   assert(row->nlocks <= (unsigned int)(row->nuses));

   debugMessage("capture row <%s> with nuses=%d and nlocks=%d\n", row->name, row->nuses, row->nlocks);
   row->nuses++;
}

/** decreases usage counter of LP row, and frees memory if necessary */
RETCODE SCIProwRelease(
   ROW**            row,                /**< pointer to LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   assert(memhdr != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->nuses >= 1);
   assert((*row)->nlocks < (unsigned int)((*row)->nuses));

   debugMessage("release row <%s> with nuses=%d and nlocks=%d\n", (*row)->name, (*row)->nuses, (*row)->nlocks);
   (*row)->nuses--;
   if( (*row)->nuses == 0 )
   {
      CHECK_OKAY( SCIProwFree(row, memhdr, set, lp) );
   }

   *row = NULL;

   return SCIP_OKAY;
}

/** locks an unmodifiable row, which forbids further changes; has no effect on modifiable rows */
void SCIProwLock(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   /* check, if row is modifiable */
   if( !row->modifiable )
   {
      debugMessage("lock row <%s> with nuses=%d and nlocks=%d\n", row->name, row->nuses, row->nlocks);
      row->nlocks++;
   }
}

/** unlocks a lock of an unmodifiable row; a row with no sealed lock may be modified; has no effect on modifiable rows */
void SCIProwUnlock(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   /* check, if row is modifiable */
   if( !row->modifiable )
   {
      debugMessage("unlock row <%s> with nuses=%d and nlocks=%d\n", row->name, row->nuses, row->nlocks);
      assert(row->nlocks > 0);
      row->nlocks--;
   }
}

/** adds a previously non existing coefficient to an LP row */
RETCODE SCIProwAddCoef(
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   COL*             col,                /**< LP column */
   Real             val                 /**< value of coefficient */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);

   CHECK_OKAY( rowAddCoef(row, memhdr, set, lp, col, val, -1) );

   checkLinks(lp);

   return SCIP_OKAY;
}

/** deletes coefficient from row */
RETCODE SCIProwDelCoef(
   ROW*             row,                /**< row to be changed */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   COL*             col                 /**< coefficient to be deleted */
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
      errorMessage("coefficient for column <%s> doesn't exist in row <%s>\n", SCIPvarGetName(col->var), row->name);
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
      CHECK_OKAY( colDelCoefPos(col, set, lp, row->linkpos[pos]) );
   }

   /* delete the column from the row's col vector */
   CHECK_OKAY( rowDelCoefPos(row, set, lp, pos) );
   
   checkLinks(lp);

   return SCIP_OKAY;
}

/** changes or adds a coefficient to an LP row */
RETCODE SCIProwChgCoef(
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   COL*             col,                /**< LP column */
   Real             val                 /**< value of coefficient */
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
      CHECK_OKAY( rowAddCoef(row, memhdr, set, lp, col, val, -1) );
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
         CHECK_OKAY( colChgCoefPos(col, set, lp, row->linkpos[pos], val) );
      }

      /* change the coefficient in the row */
      CHECK_OKAY( rowChgCoefPos(row, set, lp, pos, val) );
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** increases value of an existing or nonexisting coefficient in an LP row */
RETCODE SCIProwIncCoef(
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   COL*             col,                /**< LP column */
   Real             incval              /**< value to add to the coefficient */
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
      CHECK_OKAY( rowAddCoef(row, memhdr, set, lp, col, incval, -1) );
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
         CHECK_OKAY( colChgCoefPos(col, set, lp, row->linkpos[pos], row->vals[pos] + incval) );
      }

      /* change the coefficient in the row */
      CHECK_OKAY( rowChgCoefPos(row, set, lp, pos, row->vals[pos] + incval) );
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** changes constant value of a row */
RETCODE SCIProwChgConstant(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   Real             constant            /**< new constant value */
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
         CHECK_OKAY( rowSideChanged(row, set, lp, SCIP_SIDETYPE_LEFT) );
      }
      if( !SCIPsetIsInfinity(set, row->rhs) )
      {
         CHECK_OKAY( rowSideChanged(row, set, lp, SCIP_SIDETYPE_RIGHT) );
      }

      row->constant = constant;
   }

   return SCIP_OKAY;
}

/** add constant value to a row */
RETCODE SCIProwAddConstant(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   Real             addval              /**< constant value to add to the row */
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
      CHECK_OKAY( SCIProwChgConstant(row, set, stat, lp, row->constant + addval) );
   }

   return SCIP_OKAY;
}

/** changes left hand side of LP row */
RETCODE SCIProwChgLhs(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             lhs                 /**< new left hand side */
   )
{
   assert(row != NULL);
   assert(lp != NULL);
   assert(!lp->diving);

   if( !SCIPsetIsEQ(set, row->lhs, lhs) )
   {
      row->lhs = lhs;
      CHECK_OKAY( rowSideChanged(row, set, lp, SCIP_SIDETYPE_LEFT) );
   }

   return SCIP_OKAY;
}

/** changes right hand side of LP row */
RETCODE SCIProwChgRhs(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             rhs                 /**< new right hand side */
   )
{
   assert(row != NULL);
   assert(lp != NULL);
   assert(!lp->diving);

   if( !SCIPsetIsEQ(set, row->rhs, rhs) )
   {
      row->rhs = rhs;
      CHECK_OKAY( rowSideChanged(row, set, lp, SCIP_SIDETYPE_RIGHT) );
   }

   return SCIP_OKAY;
}

/** additional scalars that are tried in integrality scaling */
static const Real scalars[] = {3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0};
static const int nscalars = 9;

/** tries to find a value, such that all row coefficients, if scaled with this value become integral */
RETCODE SCIProwCalcIntegralScalar(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   Real             mindelta,           /**< minimal allowed difference s*c - i of scaled coefficient s*c and integral i */
   Real             maxdelta,           /**< maximal allowed difference s*c - i of scaled coefficient s*c and integral i */
   Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   Real             maxscale,           /**< maximal allowed scalar */
   Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   Bool*            success             /**< stores whether returned value is valid */
   )
{
   COL* col;
   Real val;
   Real minval;
   Real usedtol;
   Bool fractional;
   int c;

   assert(row != NULL);
   assert(row->len == 0 || row->cols != NULL);
   assert(row->len == 0 || row->cols_index != NULL);
   assert(row->len == 0 || row->vals != NULL);
   assert(maxdnom >= 1);
   assert(mindelta < 0.0);
   assert(maxdelta > 0.0);
   assert(success != NULL);

   debugMessage("trying to find rational representation for row <%s> (contvars: %d)\n", SCIProwGetName(row), usecontvars);

   if( intscalar != NULL )
      *intscalar = SCIP_INVALID;
   *success = FALSE;

   /* nothing to do, if row is empty */
   if( row->len == 0 )
   {
      if( intscalar != NULL )
         *intscalar = 1.0;
      *success = TRUE;
      return SCIP_OKAY;
   }

   /* check, if there are fractional coefficients in the row that have to be made integral, 
    * and get minimal value that has to be made integral
    */
   fractional = FALSE;
   minval = SCIPsetInfinity(set);
   for( c = 0; c < row->len; ++c )
   {
      col = row->cols[c];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(col->var) == col);
      val = row->vals[c];
      assert(!SCIPsetIsZero(set, val));

      if( (usecontvars || SCIPcolIsIntegral(col)) && !SCIPsetIsIntegral(set, val) )
      {
         Real absval;

         fractional = TRUE;
         absval = REALABS(val);
         minval = MIN(minval, absval);
      }
   }

   /* if fractional coefficients exist, try to find a rational representation */
   if( fractional )
   {
      Bool scalable;
      Bool twomult;
      Real absval;
      Real scaleval;
      Real twomultval;
      int s;

      /* try, if row coefficients can be made integral by 
       *  - multiplying them with the reciprocal of the smallest coefficient and a power of 2
       *  - by multiplying them by a power of 2
       */
      assert(SCIPsetIsPositive(set, minval));
      assert(!SCIPsetIsInfinity(set, minval));
      scalable = TRUE;
      scaleval = 1.0/minval;
      twomult = TRUE;
      twomultval = 1.0;
      for( c = 0; c < row->len && (scalable || twomult); ++c )
      {
         /* don't look at continuous variables, if we don't have to */
         if( !usecontvars && !SCIPcolIsIntegral(row->cols[c]) )
            continue;

         /* check, if the coefficient can be scaled with a simple scalar */
         val = row->vals[c];
         absval = REALABS(val);
         if( scalable )
         {
            while( scaleval <= maxscale
               && (absval * scaleval < 0.5 || !isIntegralScalar(val, scaleval, mindelta, maxdelta)) )
            {
               for( s = 0; s < nscalars; ++s )
               {
                  if( isIntegralScalar(val, scaleval * scalars[s], mindelta, maxdelta) )
                  {
                     scaleval *= scalars[s];
                     break;
                  }
               }
               if( s >= nscalars )
                  scaleval *= 2.0;
            }
            scalable = (scaleval <= maxscale);
            debugMessage(" -> val=%g, scaleval=%g, val*scaleval=%g, scalable=%d\n", 
               val, scaleval, val*scaleval, scalable);
         }
         if( twomult )
         {
            while( twomultval <= maxscale
               && (absval * twomultval < 0.5 || !isIntegralScalar(val, twomultval, mindelta, maxdelta)) )
            {
               for( s = 0; s < nscalars; ++s )
               {
                  if( isIntegralScalar(val, twomultval * scalars[s], mindelta, maxdelta) )
                  {
                     twomultval *= scalars[s];
                     break;
                  }
               }
               if( s >= nscalars )
                  twomultval *= 2.0;
            }
            twomult = (twomultval <= maxscale);
            debugMessage(" -> val=%g, twomult=%g, val*twomult=%g, twomultable=%d\n",
               val, twomultval, val*twomultval, twomult);
         }
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
         debugMessage(" -> integrality can be achieved by scaling with %g (minval=%g)\n", scaleval, minval);
      }
      else if( twomult )
      {
         /* make row coefficients integral by multiplying them with a power of 2 */
         assert(twomultval <= maxscale);
         if( intscalar != NULL )
            *intscalar = twomultval;
         *success = TRUE;
         debugMessage(" -> integrality can be achieved by scaling with %g (power of 2)\n", twomultval);
      }
      else
      {
         Longint gcd;
         Longint scm;
         Longint nominator;
         Longint denominator;
         Bool rational;
         
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
                  rational = ((Real)scm/(Real)gcd <= maxscale);
                  debugMessage(" -> first rational: val: %g == %lld/%lld, gcd=%lld, scm=%lld, rational=%d\n",
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
                  rational = ((Real)scm/(Real)gcd <= maxscale);
                  debugMessage(" -> next rational : val: %g == %lld/%lld, gcd=%lld, scm=%lld, rational=%d\n",
                     val, nominator, denominator, gcd, scm, rational);
               }
            }
         }

         if( rational )
         {
            /* make row coefficients integral by multiplying them with the smallest common multiple of the denominators */
            assert((Real)scm/(Real)gcd <= maxscale);
            if( intscalar != NULL )
               *intscalar = (Real)scm/(Real)gcd;
            *success = TRUE;
            debugMessage(" -> integrality can be achieved by scaling with %g (rational:%lld/%lld)\n", 
               (Real)scm/(Real)gcd, scm, gcd);
         }
         else
         {
            assert(!(*success));
            debugMessage(" -> rationalizing failed: gcd=%lld, scm=%lld, lastval=%g\n", gcd, scm, val);
         }
      }
   }
   else
   {
      /* all coefficients are already integral */
      if( intscalar != NULL )
         *intscalar = 1.0;
      *success = TRUE;
      debugMessage(" -> row is already integral\n");
   }

   return SCIP_OKAY;
}

/** tries to scale row, s.t. all coefficients become integral */
RETCODE SCIProwMakeIntegral(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   Real             mindelta,           /**< minimal allowed difference s*c - i of scaled coefficient s*c and integral i */
   Real             maxdelta,           /**< maximal allowed difference s*c - i of scaled coefficient s*c and integral i */
   Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   Real             maxscale,           /**< maximal value to scale row with */
   Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   Bool*            success             /**< stores whether row could be made rational */
   )
{
   Real intscalar;

   assert(success != NULL);

   /* calculate scalar to make coefficients integral */
   CHECK_OKAY( SCIProwCalcIntegralScalar(row, set, mindelta, maxdelta, maxdnom, maxscale, usecontvars,
         &intscalar, success) );

   if( *success )
   {
      /* scale the row */
      CHECK_OKAY( rowScale(row, set, stat, lp, intscalar, usecontvars, mindelta, maxdelta) );
   }

   return SCIP_OKAY;
}

/** sorts row entries such that LP columns precede non-LP columns and inside both parts lower column indices precede
 *  higher ones
 */
void SCIProwSort(
   ROW*             row                 /**< row to be sorted */
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
   ROW*             row,                /**< row to be sorted */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);
   assert(!row->delaysort);
   assert(row->nunlinked == row->len);
   assert(row->nlpcols == 0);

   debugMessage("merging row <%s>\n", row->name);

   /* do nothing on empty rows; if row is sorted, nothing has to be done */
   if( row->len > 0 && (!row->lpcolssorted || !row->nonlpcolssorted) )
   {
      COL** cols;
      int* cols_index;
      Real* vals;
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
void SCIProwDelaySort(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(!row->delaysort);

   row->delaysort = TRUE;
}

/** disables delaying of row sorting, sorts row and merges coefficients with equal columns */
void SCIProwForceSort(
   ROW*             row,                /**< LP row */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);
   assert(row->delaysort);

   row->delaysort = FALSE;
   rowMerge(row, set);
}

/** recalculates the current activity of a row */
void SCIProwRecalcLPActivity(
   ROW*             row,                /**< LP row */
   STAT*            stat                /**< problem statistics */
   )
{
   COL* col;
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
Real SCIProwGetLPActivity(
   ROW*             row,                /**< LP row */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< current LP data */
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
Real SCIProwGetLPFeasibility(
   ROW*             row,                /**< LP row */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< current LP data */
   )
{
   Real activity;

   assert(row != NULL);

   activity = SCIProwGetLPActivity(row, stat, lp);

   return MIN(row->rhs - activity, activity - row->lhs);
}

/** calculates the current pseudo activity of a row */
void SCIProwRecalcPseudoActivity(
   ROW*             row,                /**< row data */
   STAT*            stat                /**< problem statistics */
   )
{
   COL* col;
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
Real SCIProwGetPseudoActivity(
   ROW*             row,                /**< LP row */
   STAT*            stat                /**< problem statistics */
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
Real SCIProwGetPseudoFeasibility(
   ROW*             row,                /**< LP row */
   STAT*            stat                /**< problem statistics */
   )
{
   Real pseudoactivity;

   assert(row != NULL);

   pseudoactivity = SCIProwGetPseudoActivity(row, stat);

   return MIN(row->rhs - pseudoactivity, pseudoactivity - row->lhs);
}

/** returns the activity of a row for a given solution */
RETCODE SCIProwGetSolActivity(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   SOL*             sol,                /**< primal CIP solution */
   Real*            solactivity         /**< pointer to store the row's activity for the solution */
   )
{
   COL* col;
   int i;

   assert(row != NULL);
   assert(solactivity != NULL);

   *solactivity = row->constant;
   for( i = 0; i < row->len; ++i )
   {
      col = row->cols[i];
      assert(col != NULL);
      assert((i < row->nlpcols) == (row->linkpos[i] >= 0 && col->lppos >= 0));
      (*solactivity) += row->vals[i] * SCIPsolGetVal(sol, stat, col->var);
   }

   *solactivity = MAX(*solactivity, -SCIPsetInfinity(set));
   *solactivity = MIN(*solactivity, +SCIPsetInfinity(set));

   return SCIP_OKAY;
}

/** returns the feasibility of a row for the given solution */
RETCODE SCIProwGetSolFeasibility(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   SOL*             sol,                /**< primal CIP solution */
   Real*            solfeasibility      /**< pointer to store the row's feasibility for the solution */
   )
{
   Real solactivity;

   assert(row != NULL);
   assert(solfeasibility != NULL);

   CHECK_OKAY( SCIProwGetSolActivity(row, set, stat, sol, &solactivity) );

   *solfeasibility = MIN(row->rhs - solactivity, solactivity - row->lhs);

   return SCIP_OKAY;
}

/** calculates minimal and maximal activity of row w.r.t. the column's bounds */
static
void rowCalcActivityBounds(
   ROW*             row,                /**< row data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics data */
   )
{
   COL* col;
   Real val;
   Bool mininfinite;
   Bool maxinfinite;
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
Real SCIProwGetMinActivity(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics data */
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
Real SCIProwGetMaxActivity(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics data */
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
Real SCIProwGetMaxval(
   ROW*             row,                /**< LP row */
   SET*             set                 /**< global SCIP settings */
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
Real SCIProwGetMinval(
   ROW*             row,                /**< LP row */
   SET*             set                 /**< global SCIP settings */
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
Real SCIProwGetEfficacy(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp                  /**< current LP data */
   )
{
   Real norm;
   Real feasibility;

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
      errorMessage("invalid efficacy norm parameter '%c'\n", set->sepa_efficacynorm);
      abort();
   }

   norm = MAX(norm, SCIPsetSumepsilon(set));
   feasibility = SCIProwGetLPFeasibility(row, stat, lp);

   return -feasibility / norm;
}

/** returns whether the row's efficacy with respect to the current LP solution is greater than the minimal cut efficacy */
Bool SCIProwIsEfficacious(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< current LP data */
   Bool             root                /**< should the root's minimal cut efficacy be used? */
   )
{
   Real efficacy;

   efficacy = SCIProwGetEfficacy(row, set, stat, lp);

   return SCIPsetIsEfficacious(set, root, efficacy);
}

/** returns the scalar product of the coefficient vectors of the two given rows */
Real SCIProwGetScalarProduct(
   ROW*             row1,               /**< first LP row */
   ROW*             row2                /**< second LP row */
   )
{
   Real scalarprod;
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
Real SCIProwGetParallelism(
   ROW*             row1,               /**< first LP row */
   ROW*             row2                /**< second LP row */
   )
{
   Real scalarprod;

   scalarprod = SCIProwGetScalarProduct(row1, row2);
   
   return (REALABS(scalarprod) / (SCIProwGetNorm(row1) * SCIProwGetNorm(row2)));
}

/** returns the degree of orthogonality between the hyperplanes defined by the two row vectors v, w:
 *  o = 1 - |v*w|/(|v|*|w|);
 *  the hyperplanes are orthogonal, iff p = 1, they are parallel, iff p = 0
 */
Real SCIProwGetOrthogonality(
   ROW*             row1,               /**< first LP row */
   ROW*             row2                /**< second LP row */
   )
{
   return 1.0 - SCIProwGetParallelism(row1, row2);
}

/** gets parallelism of row with objective function: if the returned value is 1, the row is parellel to the objective
 *  function, if the value is 0, it is orthogonal to the objective function
 */
Real SCIProwGetObjParallelism(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   Real prod;
   Real parallelism;

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
   ROW*             row,                /**< LP row */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int i;

   assert(row != NULL);

   if( file == NULL )
      file = stdout;

   /* print left hand side */
   fprintf(file, "%g <= ", row->lhs);

   /* print coefficients */
   if( row->len == 0 )
      fprintf(file, "0 ");
   for( i = 0; i < row->len; ++i )
   {
      assert(row->cols[i] != NULL);
      assert(row->cols[i]->var != NULL);
      assert(SCIPvarGetName(row->cols[i]->var) != NULL);
      assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
      fprintf(file, "%+g<%s> ", row->vals[i], SCIPvarGetName(row->cols[i]->var));
   }

   /* print constant */
   if( REALABS(row->constant) > SCIP_DEFAULT_EPSILON )
      fprintf(file, "%+g ", row->constant);

   /* print right hand side */
   fprintf(file, "<= %g\n", row->rhs);
}


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** get number of nonzero entries in row vector */
int SCIProwGetNNonz(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->len;
}

/** get number of nonzero entries in row vector, that correspond to columns currently in the LP;
 *  Warning! This method is only applicable on rows, that are completely linked to their columns (e.g. a row
 *  that is in the current LP and the LP was solved, or a row that was in a solved LP and didn't change afterwards
 */
int SCIProwGetNLPNonz(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->nunlinked == 0);

   return row->nlpcols;
}

/** gets array with columns of nonzero entries */
COL** SCIProwGetCols(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->cols;
}

/** gets array with coefficients of nonzero entries */
Real* SCIProwGetVals(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->vals;
}

/** gets constant shift of row */
Real SCIProwGetConstant(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->constant;
}

/** gets euclidean norm of row vector */
Real SCIProwGetNorm(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return sqrt(row->sqrnorm);
}

/** gets sum norm of row vector (sum of absolute values of coefficients) */
Real SCIProwGetSumNorm(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->sumnorm;
}

/** returns the left hand side of the row */
Real SCIProwGetLhs(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->lhs;
}

/** returns the right hand side of the row */
Real SCIProwGetRhs(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->rhs;
}

/** gets the dual LP solution of a row */
Real SCIProwGetDualsol(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   if( row->lppos >= 0 )
      return row->dualsol;
   else
      return 0.0;
}

/** returns the name of the row */
const char* SCIProwGetName(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->name;
}

/** gets unique index of row */
int SCIProwGetIndex(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->index;
}

/** returns TRUE iff the activity of the row (without the row's constant) is always integral in a feasible solution */
Bool SCIProwIsIntegral(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->integral;
}

/** returns TRUE iff row is only valid locally */
Bool SCIProwIsLocal(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->local;
}

/** returns TRUE iff row is modifiable during node processing (subject to column generation) */
Bool SCIProwIsModifiable(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->modifiable;
}

/** returns TRUE iff row is removeable from the LP (due to aging or cleanup) */
Bool SCIProwIsRemoveable(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->removeable;
}

/** gets position of row in current LP, or -1 if it is not in LP */
int SCIProwGetLPPos(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert((row->lppos == -1) == (row->lpdepth == -1));

   return row->lppos;
}

/** gets depth in the tree where the row entered the LP, or -1 if it is not in LP */
int SCIProwGetLPDepth(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert((row->lppos == -1) == (row->lpdepth == -1));

   return row->lpdepth;
}

/** returns TRUE iff row is member of current LP */
Bool SCIProwIsInLP(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert((row->lppos == -1) == (row->lpdepth == -1));

   return (row->lppos >= 0);
}

#endif




/*
 * LP solver data update
 */

/** resets column data to represent a column not in the LP solver */
static
void markColDeleted(
   COL*             col                 /**< column to be marked deleted */
   )
{
   assert(col != NULL);

   col->lpipos = -1;
   col->primsol = 0.0;
   col->redcost = SCIP_INVALID;
   col->farkascoef = SCIP_INVALID;
   col->strongbranchdown = SCIP_INVALID;
   col->strongbranchup = SCIP_INVALID;
   col->validredcostlp = -1;
   col->validfarkaslp = -1;
   col->strongbranchitlim = -1;
}

/** applies all cached column removals to the LP solver */
static
RETCODE lpFlushDelCols(
   LP*              lp                  /**< current LP data */
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
      debugMessage("flushing col deletions: shrink LP from %d to %d colums\n", lp->nlpicols, lp->lpifirstchgcol);
      CHECK_OKAY( SCIPlpiDelCols(lp->lpi, lp->lpifirstchgcol, lp->nlpicols-1) );
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
RETCODE lpFlushAddCols(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   Real* obj;
   Real* lb;
   Real* ub;
   int* beg;
   int* ind;
   Real* val;
   char** name;
   COL* col;
   Real infinity;
   int c;
   int pos;
   int nnonz;
   int naddcols;
   int naddcoefs;
   int i;
   int lpipos;

   assert(lp != NULL);
   assert(lp->lpifirstchgcol == lp->nlpicols);
   assert(memhdr != NULL);
   assert(set != NULL);

   /* if there are no columns to add, we are ready */
   if( lp->ncols == lp->nlpicols )
      return SCIP_OKAY;

   /* add the additional columns */
   assert(!lp->diving);
   assert(lp->ncols > lp->nlpicols);
   CHECK_OKAY( ensureLpicolsSize(lp, set, lp->ncols) );

   /* get the solver's infinity value */
   infinity = SCIPlpiInfinity(lp->lpi);

   /* count the (maximal) number of added coefficients, calculate the number of added columns */
   naddcols = lp->ncols - lp->nlpicols;
   naddcoefs = 0;
   for( c = lp->nlpicols; c < lp->ncols; ++c )
      naddcoefs += lp->cols[c]->len;
   assert(naddcols > 0);

   /* get temporary memory for changes */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &obj, naddcols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &lb, naddcols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &ub, naddcols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &beg, naddcols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &ind, naddcoefs) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &val, naddcoefs) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &name, naddcols) );
   
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

      debugMessage("flushing added column <%s>: ", SCIPvarGetName(col->var));
      debug( SCIPcolPrint(col, NULL) );

      /* Because the column becomes a member of the LP solver, it now can take values
       * different from zero. That means, we have to include the column in the corresponding
       * row vectors.
       */
      CHECK_OKAY( colLink(col, memhdr, set, lp) );

      lp->lpicols[c] = col;
      col->lpipos = c;
      col->primsol = SCIP_INVALID;
      col->redcost = SCIP_INVALID;
      col->farkascoef = SCIP_INVALID;
      col->strongbranchdown = SCIP_INVALID;
      col->strongbranchup = SCIP_INVALID;
      col->validredcostlp = -1;
      col->validfarkaslp = -1;
      col->strongbranchitlim = -1;
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
   debugMessage("flushing col additions: enlarge LP from %d to %d colums\n", lp->nlpicols, lp->ncols);
   CHECK_OKAY( SCIPlpiAddCols(lp->lpi, naddcols, obj, lb, ub, name, nnonz, beg, ind, val) );
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
   ROW*             row                 /**< row to be marked deleted */
   )
{
   assert(row != NULL);

   row->lpipos = -1;
   row->dualsol = 0.0;
   row->activity = SCIP_INVALID;
   row->dualfarkas = 0.0;
   row->validactivitylp = -1;
}

/** applies all cached row removals to the LP solver */
static
RETCODE lpFlushDelRows(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
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
      debugMessage("flushing row deletions: shrink LP from %d to %d rows\n", lp->nlpirows, lp->lpifirstchgrow);
      CHECK_OKAY( SCIPlpiDelRows(lp->lpi, lp->lpifirstchgrow, lp->nlpirows-1) );
      for( i = lp->lpifirstchgrow; i < lp->nlpirows; ++i )
      {
         markRowDeleted(lp->lpirows[i]);
         CHECK_OKAY( SCIProwRelease(&lp->lpirows[i], memhdr, set, lp) );
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
RETCODE lpFlushAddRows(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   Real* lhs;
   Real* rhs;
   int* beg;
   int* ind;
   Real* val;
   char** name;
   ROW* row;
   Real infinity;
   int r;
   int pos;
   int nnonz;
   int naddrows;
   int naddcoefs;
   int i;
   int lpipos;

   assert(lp != NULL);
   assert(lp->lpifirstchgrow == lp->nlpirows);
   assert(memhdr != NULL);
      
   /* if there are no rows to add, we are ready */
   if( lp->nrows == lp->nlpirows )
      return SCIP_OKAY;

   /* add the additional rows */
   assert(!lp->diving);
   assert(lp->nrows > lp->nlpirows);
   CHECK_OKAY( ensureLpirowsSize(lp, set, lp->nrows) );

   /* get the solver's infinity value */
   infinity = SCIPlpiInfinity(lp->lpi);

   /* count the (maximal) number of added coefficients, calculate the number of added rows */
   naddrows = lp->nrows - lp->nlpirows;
   naddcoefs = 0;
   for( r = lp->nlpirows; r < lp->nrows; ++r )
      naddcoefs += lp->rows[r]->len;
   assert(naddrows > 0);

   /* get temporary memory for changes */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &lhs, naddrows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &rhs, naddrows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &beg, naddrows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &ind, naddcoefs) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &val, naddcoefs) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &name, naddrows) );
   
   /* fill temporary memory with row data */
   nnonz = 0;
   for( pos = 0, r = lp->nlpirows; r < lp->nrows; ++pos, ++r )
   {
      row = lp->rows[r];
      assert(row != NULL);
      assert(row->lppos == r);
      assert(nnonz + row->nlpcols <= naddcoefs);

      debugMessage("flushing added row <%s>: ", row->name);
      debug( SCIProwPrint(row, NULL) );

      /* Because the row becomes a member of the LP solver, its dual variable now can take values
       * different from zero. That means, we have to include the row in the corresponding
       * column vectors.
       */
      CHECK_OKAY( rowLink(row, memhdr, set, lp) );

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

      debugMessage("flushing added row (LPI): %+g <=", lhs[pos]);
      for( i = 0; i < row->nlpcols; ++i )
      {
         assert(row->cols[i] != NULL);
         lpipos = row->cols[i]->lpipos;
         if( lpipos >= 0 )
         {
            assert(lpipos < lp->ncols);
            assert(nnonz < naddcoefs);
            debug( printf(" %+gx%d(<%s>)", row->vals[i], lpipos+1, SCIPvarGetName(row->cols[i]->var)) );
            ind[nnonz] = lpipos;
            val[nnonz] = row->vals[i];
            nnonz++;
         }
      }
      debug( printf(" <= %+g\n", rhs[pos]) );
#ifndef NDEBUG
      for( i = row->nlpcols; i < row->len; ++i )
      {
         assert(row->cols[i] != NULL);
         assert(row->cols[i]->lpipos == -1); /* because the column deletions are already performed */
      }
#endif
   }

   /* call LP interface */
   debugMessage("flushing row additions: enlarge LP from %d to %d rows\n", lp->nlpirows, lp->nrows);
   CHECK_OKAY( SCIPlpiAddRows(lp->lpi, naddrows, lhs, rhs, name, nnonz, beg, ind, val) );
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
RETCODE lpFlushChgCols(
   LP*              lp,                 /**< current LP data */
   SET*             set                 /**< global SCIP settings */
   )
{
   COL* col;
   int* objind;
   int* bdind;
   Real* obj;
   Real* lb;
   Real* ub;
   Real infinity;
   int nobjchg;
   int nbdchg;
   int i;

   assert(lp != NULL);

   if( lp->nchgcols == 0 )
      return SCIP_OKAY;

   /* get the solver's infinity value */
   infinity = SCIPlpiInfinity(lp->lpi);

   /* get temporary memory for changes */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &objind, lp->ncols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &obj, lp->ncols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &bdind, lp->ncols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &lb, lp->ncols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &ub, lp->ncols) );

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
         Real lpiobj;
         Real lpilb;
         Real lpiub;

         CHECK_OKAY( SCIPlpiGetObj(lp->lpi, col->lpipos, col->lpipos, &lpiobj) );
         CHECK_OKAY( SCIPlpiGetBounds(lp->lpi, col->lpipos, col->lpipos, &lpilb, &lpiub) );
         assert(SCIPsetIsFeasEQ(set, lpiobj, col->flushedobj));
         assert(SCIPsetIsFeasEQ(set, lpilb, col->flushedlb));
         assert(SCIPsetIsFeasEQ(set, lpiub, col->flushedub));
#endif
         if( col->objchanged )
         {
            Real newobj;

            newobj = col->obj;
            if( col->flushedobj != newobj )
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
            Real newlb;
            Real newub;

            newlb = (SCIPsetIsInfinity(set, -col->lb) ? -infinity : col->lb);
            newub = (SCIPsetIsInfinity(set, col->ub) ? infinity : col->ub);
            if( col->flushedlb != newlb || col->flushedub != newub )
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
      debugMessage("flushing bound changes: change %d objective values of %d changed columns\n", nobjchg, lp->nchgcols);
      CHECK_OKAY( SCIPlpiChgObj(lp->lpi, nobjchg, objind, obj) );

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   /* change bounds in LP */
   if( nbdchg > 0 )
   {
      debugMessage("flushing bound changes: change %d bounds of %d changed columns\n", nbdchg, lp->nchgcols);
      CHECK_OKAY( SCIPlpiChgBounds(lp->lpi, nbdchg, bdind, lb, ub) );

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
RETCODE lpFlushChgRows(
   LP*              lp,                 /**< current LP data */
   SET*             set                 /**< global SCIP settings */
   )
{
   ROW* row;
   int* ind;
   Real* lhs;
   Real* rhs;
   Real infinity;
   int i;
   int nchg;

   assert(lp != NULL);

   if( lp->nchgrows == 0 )
      return SCIP_OKAY;

   assert(!lp->diving);

   /* get the solver's infinity value */
   infinity = SCIPlpiInfinity(lp->lpi);

   /* get temporary memory for changes */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &ind, lp->nrows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &lhs, lp->nrows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &rhs, lp->nrows) );

   /* collect all cached left and right hand side changes */
   nchg = 0;
   for( i = 0; i < lp->nchgrows; ++i )
   {
      row = lp->chgrows[i];
      assert(row != NULL);

      if( row->lpipos >= 0 )
      {
#ifndef NDEBUG
         Real lpilhs;
         Real lpirhs;

         CHECK_OKAY( SCIPlpiGetSides(lp->lpi, row->lpipos, row->lpipos, &lpilhs, &lpirhs) );
         assert(SCIPsetIsSumEQ(set, lpilhs, row->flushedlhs));
         assert(SCIPsetIsSumEQ(set, lpirhs, row->flushedrhs));
#endif
         if( row->lhschanged || row->rhschanged )
         {
            Real newlhs;
            Real newrhs;

            newlhs = (SCIPsetIsInfinity(set, -row->lhs) ? -infinity : row->lhs - row->constant);
            newrhs = (SCIPsetIsInfinity(set, row->rhs) ? infinity : row->rhs - row->constant);
            if( row->flushedlhs != newlhs || row->flushedrhs != newrhs )
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
      debugMessage("flushing side changes: change %d sides of %d rows\n", nchg, lp->nchgrows);
      CHECK_OKAY( SCIPlpiChgSides(lp->lpi, nchg, ind, lhs, rhs) );

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
RETCODE SCIPlpFlush(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(memhdr != NULL);
   
   debugMessage("flushing LP changes: old (%d cols, %d rows), chgcol=%d, chgrow=%d, new (%d cols, %d rows), flushed=%d\n",
      lp->nlpicols, lp->nlpirows, lp->lpifirstchgcol, lp->lpifirstchgrow, lp->ncols, lp->nrows, lp->flushed);

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

   CHECK_OKAY( lpFlushDelCols(lp) );
   CHECK_OKAY( lpFlushDelRows(lp, memhdr, set) );
   CHECK_OKAY( lpFlushChgCols(lp, set) );
   CHECK_OKAY( lpFlushChgRows(lp, set) );
   CHECK_OKAY( lpFlushAddCols(lp, memhdr, set) );
   CHECK_OKAY( lpFlushAddRows(lp, memhdr, set) );

   lp->flushed = TRUE;

   checkLinks(lp);

   return SCIP_OKAY;
}

/** marks the LP to be flushed, even if the LP thinks it is not flushed */
RETCODE SCIPlpMarkFlushed(
   LP*              lp,                 /**< current LP data */
   SET*             set                 /**< global SCIP settings */
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
      COL* col;

      col = lp->chgcols[i];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(col->var) == col);

      if( col->lpipos >= 0 )
      {
#ifndef NDEBUG
         Real lpiobj;
         Real lpilb;
         Real lpiub;

         CHECK_OKAY( SCIPlpiGetObj(lp->lpi, col->lpipos, col->lpipos, &lpiobj) );
         CHECK_OKAY( SCIPlpiGetBounds(lp->lpi, col->lpipos, col->lpipos, &lpilb, &lpiub) );
         assert(SCIPsetIsSumEQ(set, lpiobj, col->flushedobj));
         assert(SCIPsetIsSumEQ(set, lpilb, col->flushedlb));
         assert(SCIPsetIsSumEQ(set, lpiub, col->flushedub));
         assert(col->flushedobj == col->obj);
         assert(col->flushedlb == (SCIPsetIsInfinity(set, -col->lb) ? -SCIPlpiInfinity(lp->lpi) : col->lb));
         assert(col->flushedub == (SCIPsetIsInfinity(set, col->ub) ? SCIPlpiInfinity(lp->lpi) : col->ub));
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
      ROW* row;

      row = lp->chgrows[i];
      assert(row != NULL);

      if( row->lpipos >= 0 )
      {
#ifndef NDEBUG
         Real lpilhs;
         Real lpirhs;

         CHECK_OKAY( SCIPlpiGetSides(lp->lpi, row->lpipos, row->lpipos, &lpilhs, &lpirhs) );
         assert(SCIPsetIsSumEQ(set, lpilhs, row->flushedlhs));
         assert(SCIPsetIsSumEQ(set, lpirhs, row->flushedrhs));
         assert(row->flushedlhs ==
            (SCIPsetIsInfinity(set, -row->lhs) ? -SCIPlpiInfinity(lp->lpi) : row->lhs - row->constant));
         assert(row->flushedrhs ==
            (SCIPsetIsInfinity(set, row->rhs) ? SCIPlpiInfinity(lp->lpi) : row->rhs - row->constant));
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
   COL*             col                 /**< LP column */
   )
{
   ROW* row;
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
   ROW*             row                 /**< LP row */
   )
{
   COL* col;
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
   COL*             col                 /**< LP column */
   )
{
   ROW* row;
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
   ROW*             row                 /**< LP row */
   )
{
   COL* col;
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
RETCODE SCIPlpCreate(
   LP**             lp,                 /**< pointer to LP data object */
   SET*             set,                /**< global SCIP settings */
   const char*      name                /**< problem name */
   )
{
   Bool success;

   assert(lp != NULL);
   assert(set != NULL);
   assert(name != NULL);

   ALLOC_OKAY( allocMemory(lp) );

   /* open LP Solver interface */
   CHECK_OKAY( SCIPlpiCreate(&(*lp)->lpi, name, SCIP_OBJSEN_MINIMIZE) );

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
   (*lp)->validsollp = 0; /* the initial (empty) LP is solved with primal and dual solution of zero */
   (*lp)->validfarkaslp = -1;
   (*lp)->flushdeletedcols = FALSE;
   (*lp)->flushaddedcols = FALSE;
   (*lp)->flushdeletedrows = FALSE;
   (*lp)->flushaddedrows = FALSE;
   (*lp)->flushed = TRUE;
   (*lp)->solved = TRUE;
   (*lp)->primalfeasible = TRUE;
   (*lp)->dualfeasible = TRUE;
   (*lp)->diving = FALSE;
   (*lp)->divingobjchg = FALSE;
   (*lp)->divelpistate = NULL;
   (*lp)->lpiuobjlim = SCIPsetInfinity(set);
   (*lp)->lpifeastol = SCIPsetFeastol(set);
   (*lp)->lpidualfeastol = SCIPsetDualfeastol(set);
   (*lp)->lpifromscratch = FALSE;
   (*lp)->lpifastmip = TRUE;
   (*lp)->lpiscaling = TRUE;
   (*lp)->lpipresolving = TRUE;
   (*lp)->lpilpinfo = FALSE;
   (*lp)->lpiitlim = INT_MAX;
   (*lp)->lastwasprimal = FALSE;

   /* set default parameters in LP solver */
   CHECK_OKAY( lpSetRealpar(*lp, SCIP_LPPAR_UOBJLIM, (*lp)->lpiuobjlim, &success) );
   if( !success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver: upper objective limit cannot be set -- can lead to unnecessary simplex iterations\n");
   }
   CHECK_OKAY( lpSetRealpar(*lp, SCIP_LPPAR_FEASTOL, (*lp)->lpifeastol, &success) );
   if( !success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver: primal feasibility tolerance cannot be set -- tolerance of SCIP and LP solver may differ\n");
   }
   CHECK_OKAY( lpSetRealpar(*lp, SCIP_LPPAR_DUALFEASTOL, (*lp)->lpidualfeastol, &success) );
   if( !success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver: dual feasibility tolerance cannot be set -- tolerance of SCIP and LP solver may differ\n");
   }
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_FROMSCRATCH, (*lp)->lpifromscratch, &success) );
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_FASTMIP, (*lp)->lpifastmip, &success) );
   if( !success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver: fastmip setting not available -- SCIP parameter has no effect\n");
   }
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_SCALING, (*lp)->lpiscaling, &success) );
   if( !success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver: scaling not available -- SCIP parameter has no effect\n");
   }
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_PRESOLVING, (*lp)->lpipresolving, &success) );
   if( !success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver: presolving not available -- SCIP parameter has no effect\n");
   }
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_LPITLIM, (*lp)->lpiitlim, &success) );
   if( !success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver: iteration limit cannot be set -- can lead to unnecessary simplex iterations\n");
   }
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_LPINFO, (*lp)->lpilpinfo, &success) );
   if( !success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver: lpinfo setting not available -- SCIP parameter has no effect\n");
   }
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_PRICING, SCIP_PRICING_AUTO, &success) ); /*lint !e641*/

   return SCIP_OKAY;
}

/** frees LP data object */
RETCODE SCIPlpFree(
   LP**             lp,                 /**< pointer to LP data object */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(lp != NULL);
   assert(*lp != NULL);
   
   CHECK_OKAY( SCIPlpClear(*lp, memhdr, set) );

   /* release LPI rows */
   for( i = 0; i < (*lp)->nlpirows; ++i )
   {
      CHECK_OKAY( SCIProwRelease(&(*lp)->lpirows[i], memhdr, set, *lp) );
   }

   if( (*lp)->lpi != NULL )
   {
      CHECK_OKAY( SCIPlpiFree(&(*lp)->lpi) );
   }

   freeMemoryArrayNull(&(*lp)->lpicols);
   freeMemoryArrayNull(&(*lp)->lpirows);
   freeMemoryArrayNull(&(*lp)->chgcols);
   freeMemoryArrayNull(&(*lp)->cols);
   freeMemoryArrayNull(&(*lp)->rows);
   freeMemory(lp);

   return SCIP_OKAY;
}

/** adds a column to the LP */
RETCODE SCIPlpAddCol(
   LP*              lp,                 /**< LP data */
   SET*             set,                /**< global SCIP settings */
   COL*             col,                /**< LP column */
   int              depth               /**< depth in the tree where the column addition is performed */
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
   
   debugMessage("adding column <%s> to LP (%d rows, %d cols)\n", SCIPvarGetName(col->var), lp->nrows, lp->ncols);
#ifdef DEBUG
   {
      int i;
      printf("  (obj: %g) [%g,%g]", col->obj, col->lb, col->ub);
      for( i = 0; i < col->len; ++i )
         printf(" %+g<%s>", col->vals[i], col->rows[i]->name);
      printf("\n");
   }
#endif

   CHECK_OKAY( ensureColsSize(lp, set, lp->ncols+1) );
   lp->cols[lp->ncols] = col;
   col->lppos = lp->ncols;
   col->lpdepth = depth;
   col->age = 0;
   lp->ncols++;
   if( col->removeable )
      lp->nremoveablecols++;

   /* mark the current LP unflushed */
   lp->flushed = FALSE;

   /* update squared euclidean norm and sum norm of objective function vector */
   lp->objsqrnorm += SQR(col->obj);
   lp->objsumnorm += REALABS(col->obj);

   /* update column arrays of all linked rows */
   colUpdateAddLP(col);

   checkLinks(lp);

   return SCIP_OKAY;
}

/** adds a row to the LP and captures it */
RETCODE SCIPlpAddRow(
   LP*              lp,                 /**< LP data */
   SET*             set,                /**< global SCIP settings */
   ROW*             row,                /**< LP row */
   int              depth               /**< depth in the tree where the row addition is performed */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);
   assert(row != NULL);
   assert(row->len == 0 || row->cols != NULL);
   assert(row->lppos == -1);

   SCIProwCapture(row);
   SCIProwLock(row);

   debugMessage("adding row <%s> to LP (%d rows, %d cols)\n", row->name, lp->nrows, lp->ncols);
#ifdef DEBUG
   {
      int i;
      printf("  %g <=", row->lhs);
      for( i = 0; i < row->len; ++i )
         printf(" %+g<%s>", row->vals[i], SCIPvarGetName(row->cols[i]->var));
      if( !SCIPsetIsZero(set, row->constant) )
         printf(" %+g", row->constant);
      printf(" <= %g\n", row->rhs);
   }
#endif

   CHECK_OKAY( ensureRowsSize(lp, set, lp->nrows+1) );
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
RETCODE SCIPlpShrinkCols(
   LP*              lp,                 /**< LP data */
   int              newncols            /**< new number of columns in the LP */
   )
{
   COL* col;
   int c;

   assert(lp != NULL);
   debugMessage("shrinking LP from %d to %d columns\n", lp->ncols, newncols);
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

         /* update squared euclidean norm and sum norm of objective function vector */
         lp->objsqrnorm -= SQR(col->obj);
         lp->objsqrnorm = MAX(lp->objsqrnorm, 0.0);
         lp->objsumnorm -= REALABS(col->obj);
         lp->objsumnorm = MAX(lp->objsumnorm, 0.0);

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
RETCODE SCIPlpShrinkRows(
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              newnrows            /**< new number of rows in the LP */
   )
{
   ROW* row;
   int r;

   assert(lp != NULL);
   assert(0 <= newnrows && newnrows <= lp->nrows);

   debugMessage("shrinking LP from %d to %d rows\n", lp->nrows, newnrows);
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
         CHECK_OKAY( SCIProwRelease(&lp->rows[r], memhdr, set, lp) );
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
RETCODE SCIPlpClear(
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);

   debugMessage("clearing LP\n");
   CHECK_OKAY( SCIPlpShrinkCols(lp, 0) );
   CHECK_OKAY( SCIPlpShrinkRows(lp, memhdr, set, 0) );

   return SCIP_OKAY;
}

/** remembers number of columns and rows to track the newly added ones */
void SCIPlpMarkSize(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);

   lp->firstnewcol = lp->ncols;
   lp->firstnewrow = lp->nrows;
}

/** gets all indices of basic columns and rows: index i >= 0 corresponds to column i, index i < 0 to row -i-1 */
RETCODE SCIPlpGetBasisInd(
   LP*              lp,                 /**< LP data */
   int*             basisind            /**< pointer to store the basis indices */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(basisind != NULL);

   CHECK_OKAY( SCIPlpiGetBasisInd(lp->lpi, basisind) );

   return SCIP_OKAY;
}

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
RETCODE SCIPlpGetBase(
   LP*              lp,                 /**< LP data */
   int*             cstat,              /**< array to store column basis status, or NULL */
   int*             rstat               /**< array to store row basis status, or NULL */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   CHECK_OKAY( SCIPlpiGetBase(lp->lpi, cstat, rstat) );

   return SCIP_OKAY;
}

/** gets a row from the inverse basis matrix B^-1 */
RETCODE SCIPlpGetBInvRow(
   LP*              lp,                 /**< LP data */
   int              r,                  /**< row number */
   Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(0 <= r && r < lp->nrows);  /* the basis matrix is nrows x nrows */
   assert(coef != NULL);

   CHECK_OKAY( SCIPlpiGetBInvRow(lp->lpi, r, coef) );

   return SCIP_OKAY;
}

/** gets a row from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A) */
RETCODE SCIPlpGetBInvARow(
   LP*              lp,                 /**< LP data */
   int              r,                  /**< row number */
   Real*            binvrow,            /**< row in B^-1 from prior call to SCIPlpGetBInvRow(), or NULL */
   Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(0 <= r && r < lp->nrows);  /* the basis matrix is nrows x nrows */
   assert(coef != NULL);

   CHECK_OKAY( SCIPlpiGetBInvARow(lp->lpi, r, binvrow, coef) );

   return SCIP_OKAY;
}

/** calculates a weighted sum of all LP rows; for negative weights, the left and right hand side of the corresponding
 *  LP row are swapped in the summation
 */
RETCODE SCIPlpSumRows(
   LP*              lp,                 /**< LP data */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   Real*            weights,            /**< row weights in row summation */
   REALARRAY*       sumcoef,            /**< array to store sum coefficients indexed by variables' probindex */
   Real*            sumlhs,             /**< pointer to store the left hand side of the row summation */
   Real*            sumrhs              /**< pointer to store the right hand side of the row summation */
   )
{
   ROW* row;
   int r;
   int i;
   int idx;
   Bool lhsinfinite;
   Bool rhsinfinite;

   assert(lp != NULL);
   assert(prob != NULL);
   assert(weights != NULL);
   assert(sumcoef != NULL);
   assert(sumlhs != NULL);
   assert(sumrhs != NULL);

   /**@todo test, if a column based summation is faster */

   CHECK_OKAY( SCIPrealarrayClear(sumcoef) );
   CHECK_OKAY( SCIPrealarrayExtend(sumcoef, set, 0, prob->nvars-1) );
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
            CHECK_OKAY( SCIPrealarrayIncVal(sumcoef, set, idx, weights[r] * row->vals[i]) );
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
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< LP data */
   Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   Real             scale,              /**< additional scaling factor multiplied to all rows */
   Bool             allowlocal,         /**< should local rows be included, resulting in a locally valid summation? */
   Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size prob->nvars */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*             slacksign,          /**< stores the sign of the row's slack variable in summation */
   Bool*            emptyrow,           /**< pointer to store whether the returned row is empty */
   Bool*            localrowsused       /**< pointer to store whether local rows were used in summation */
   )
{
   ROW* row;
   Real weight;
   Real absweight;
   Real maxweight;
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
   clearMemoryArray(mircoef, prob->nvars);
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

         debugMessage("MIR: %d: row <%s>, lhs = %g, rhs = %g, scale = %g, weight = %g, slacksign = %d -> rhs = %g\n",
            r, SCIProwGetName(row), row->lhs - row->constant, row->rhs - row->constant, 
            scale, weights[r], slacksign[r], *mirrhs);
         debug(SCIProwPrint(row, NULL));
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
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   Bool             cutislocal          /**< is the cut only valid locally? */
   )
{
   Real bd;
   Bool rhsinf;
   int v;

   assert(prob != NULL);
   assert(mircoef != NULL);
   assert(mirrhs != NULL);

   rhsinf = SCIPsetIsInfinity(set, *mirrhs);
   for( v = 0; v < prob->nvars && !rhsinf; ++v )
   {
      if( mircoef[v] != 0.0 && SCIPsetIsSumZero(set, mircoef[v]) )
      {
         debugMessage("coefficient of <%s> in transformed MIR row is too small: %.12f\n",
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
   VAR*             var,                /**< active problem variable */
   Real*            closestvlb,         /**< pointer to store the value of the closest variable lower bound */
   int*             closestvlbidx       /**< pointer to store the index of the closest variable lower bound */
   )
{
   VAR** vlbvars;
   Real* vlbcoefs;
   Real* vlbconsts;
   int nvlbs;
   int i;

   assert(closestvlb != NULL);
   assert(closestvlbidx != NULL);

   *closestvlbidx = -1;
   *closestvlb = REAL_MIN;

   nvlbs = SCIPvarGetNVlbs(var);
   if( nvlbs > 0 )
   {
      vlbvars = SCIPvarGetVlbVars(var);
      vlbcoefs = SCIPvarGetVlbCoefs(var);
      vlbconsts = SCIPvarGetVlbConstants(var);
      
      for( i = 0; i < nvlbs; i++ )
      {
         Real vlbsol;
         
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
   VAR*             var,                /**< active problem variable */
   Real*            closestvub,         /**< pointer to store the value of the closest variable upper bound */
   int*             closestvubidx       /**< pointer to store the index of the closest variable upper bound */
   )
{
   VAR** vubvars;
   Real* vubcoefs;
   Real* vubconsts;
   int nvubs;
   int i;

   assert(closestvub != NULL);
   assert(closestvubidx != NULL);

   *closestvubidx = -1;
   *closestvub = REAL_MAX;

   nvubs = SCIPvarGetNVubs(var);
   if( nvubs > 0 )
   {
      vubvars = SCIPvarGetVubVars(var);
      vubcoefs = SCIPvarGetVubCoefs(var);
      vubconsts = SCIPvarGetVubConstants(var);
      
      for( i = 0; i < nvubs; i++ )
      {
         Real vubsol;
         
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
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*             varsign,            /**< stores the sign of the transformed variable in summation */
   int*             boundtype,          /**< stores the bound used for transformed variable:
                                         *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   Bool*            freevariable,       /**< stores whether a free variable was found in MIR row -> invalid summation */
   Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{
   VAR* var;
   Real varsol;
   Real bestlb;
   Real bestub;
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
         Real loclb;

         loclb = SCIPvarGetLbLocal(var);
         if( SCIPsetIsGT(set, loclb, bestlb) )
         {
            bestlb = loclb;
            bestlbtype = -2;
         }
      }
      if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         Real bestvlb;
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
         Real locub;

         locub = SCIPvarGetUbLocal(var);
         if( SCIPsetIsLT(set, locub, bestub) )
         {
            bestub = locub;
            bestubtype = -2;
         }
      }
      if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         Real bestvub;
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
            VAR** vlbvars = SCIPvarGetVlbVars(var);
            Real* vlbcoefs = SCIPvarGetVlbCoefs(var);
            Real* vlbconsts = SCIPvarGetVlbConstants(var);
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
            VAR** vubvars = SCIPvarGetVubVars(var);
            Real* vubcoefs = SCIPvarGetVubCoefs(var);
            Real* vubconsts = SCIPvarGetVubConstants(var);
            int zidx;

            assert(0 <= bestubtype && bestubtype < SCIPvarGetNVubs(var));
            zidx = SCIPvarGetProbindex(vubvars[bestubtype]);
            assert(zidx >= 0);
               
            (*mirrhs) -= mircoef[v] * vubconsts[bestubtype];
            mircoef[zidx] += mircoef[v] * vubcoefs[bestubtype];
         }
      }

      debugMessage("MIR var <%s>: varsign=%d, boundtype=%d, mircoef=%g, lb=%g, ub=%g -> rhs=%g\n", 
         SCIPvarGetName(var), varsign[v], boundtype[v], mircoef[v], bestlb, bestub, *mirrhs);
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
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*             varsign,            /**< stores the sign of the transformed variable in summation */
   int*             boundtype,          /**< stores the bound used for transformed variable (vlb/vub_idx or -1 for lb/ub)*/
   Real             f0,                 /**< fracional value of rhs */
   Real*            cutactivity         /**< pointer to store the activity of the resulting cut */
   )
{
   VAR* var;
   Real onedivoneminusf0;
   Real aj;
   Real downaj;
   Real cutaj;
   Real fj;
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
         VAR** vbz;
         Real* vbb;
         Real* vbd;
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
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< LP data */
   Real*            weights,            /**< row weights in row summation */
   Real             scale,              /**< additional scaling factor multiplied to all rows */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*             slacksign,          /**< stores the sign of the row's slack variable in summation */
   Real             f0,                 /**< fracional value of rhs */
   Real*            cutactivity         /**< pointer to update the activity of the resulting cut */
   )
{
   ROW* row;
   Real onedivoneminusf0;
   Real ar;
   Real downar;
   Real cutar;
   Real fr;
   Real mul;
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
         && ((slacksign[r] == +1 && SCIPsetIsIntegral(set, row->rhs - row->constant))
            || (slacksign[r] == -1 && SCIPsetIsIntegral(set, row->lhs - row->constant))) )
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

#ifdef DEBUG
static
void printMIR(
   PROB*            prob,               /**< problem data */
   Real*            mircoef,            /**< MIR coefficients */
   Real             mirrhs              /**< right hand side of the MIR row */
   )
{
   int i;

   assert(prob != NULL);

   printf("MIR:");
   for( i = 0; i < prob->nvars; ++i )
   {
      if( mircoef[i] != 0.0 )
         printf(" %+g<%s>", mircoef[i], SCIPvarGetName(prob->vars[i]));
   }
   printf(" <= %g\n", mirrhs);
}
#endif

/* calculates a MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a MIR cut.
 */
RETCODE SCIPlpCalcMIR(
   LP*              lp,                 /**< LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   Real             scale,              /**< additional scaling factor multiplied to all rows */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   Bool*            success,            /**< pointer to store whether the returned coefficients are a valid MIR cut */
   Bool*            cutislocal          /**< pointer to store whether the returned cut is only valid locally */
   )
{
   int* slacksign;
   int* varsign;
   int* boundtype;
   Real rhs;
   Real downrhs;
   Real f0;
   Bool emptyrow;
   Bool freevariable;
   Bool localrowsused;
   Bool localbdsused;

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

   debugMessage("calculating MIR cut (scale: %g)\n", scale);

   /**@todo test, if a column based summation is faster */

   *success = FALSE;
   *cutislocal = FALSE;

   /* allocate temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &slacksign, lp->nrows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &varsign, prob->nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &boundtype, prob->nvars) );

   /* calculate the row summation */
   sumMIRRow(set, stat, prob, lp, weights, scale, allowlocal, 
      maxweightrange, mircoef, &rhs, slacksign, &emptyrow, &localrowsused);
   assert(allowlocal || !localrowsused);
   *cutislocal = *cutislocal || localrowsused;
   if( emptyrow )
      goto TERMINATE;
   debug(printMIR(prob, mircoef, rhs));
   
   /* remove all nearly-zero coefficients from MIR row and relaxes the right hand side correspondingly in order to
    *  prevent numerical rounding errors
    */
   cleanupMIRRow(set, prob, mircoef, &rhs, *cutislocal);
   debug(printMIR(prob, mircoef, rhs));

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
   debug(printMIR(prob, mircoef, rhs));

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
   debug(printMIR(prob, mircoef, *mirrhs));

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
   debug(printMIR(prob, mircoef, *mirrhs));

   *success = TRUE;

#ifndef NDEBUG
   {
      Real act;
      int i;

      act = 0.0;
      for( i = 0; i < prob->nvars; ++i )
         act += mircoef[i] * SCIPvarGetLPSol(prob->vars[i]);
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

/** builds a weighted sum of rows, and decides whether to use the left or right hand side of the rows in summation */
static
void sumStrongCGRow(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< LP data */
   Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   Real             scale,              /**< additional scaling factor multiplied to all rows */
   Bool             allowlocal,         /**< should local rows be included, resulting in a locally valid summation? */
   Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size prob->nvars */
   Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   int*             slacksign,          /**< stores the sign of the row's slack variable in summation */
   Bool*            emptyrow,           /**< pointer to store whether the returned row is empty */
   Bool*            localrowsused       /**< pointer to store whether local rows were used in summation */
   )
{
   ROW* row;
   Real weight;
   Real absweight;
   Real maxweight;
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
   clearMemoryArray(strongcgcoef, prob->nvars);
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
               (*strongcgrhs) += weight * (row->lhs - row->constant);
            }
            else
            {
               slacksign[r] = +1;
               (*strongcgrhs) += weight * (row->rhs - row->constant);
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

         debugMessage("strong CG: %d: row <%s>, lhs = %g, rhs = %g, scale = %g, weight = %g, slacksign = %d -> rhs = %g\n",
            r, SCIProwGetName(row), row->lhs - row->constant, row->rhs - row->constant, 
            scale, weights[r], slacksign[r], *strongcgrhs);
         debug(SCIProwPrint(row, NULL));
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
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size nvars */
   Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   int*             varsign,            /**< stores the sign of the transformed variable in summation */
   int*             boundtype,          /**< stores the bound used for transformed variable:
                                         *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   Bool*            freevariable,       /**< stores whether a free variable was found in strong CG row -> invalid summation */
   Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{
   VAR* var;
   Real varsol;
   Real bestlb;
   Real bestub;
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
         Real loclb;

         loclb = SCIPvarGetLbLocal(var);
         if( SCIPsetIsGT(set, loclb, bestlb) )
         {
            bestlb = loclb;
            bestlbtype = -2;
         }
      }
      if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         Real bestvlb;
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
         Real locub;

         locub = SCIPvarGetUbLocal(var);
         if( SCIPsetIsLT(set, locub, bestub) )
         {
            bestub = locub;
            bestubtype = -2;
         }
      }
      if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         Real bestvub;
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
            VAR** vlbvars = SCIPvarGetVlbVars(var);
            Real* vlbcoefs = SCIPvarGetVlbCoefs(var);
            Real* vlbconsts = SCIPvarGetVlbConstants(var);
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
            VAR** vubvars = SCIPvarGetVubVars(var);
            Real* vubcoefs = SCIPvarGetVubCoefs(var);
            Real* vubconsts = SCIPvarGetVubConstants(var);
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

      debugMessage("strong CG var <%s>: varsign=%d, boundtype=%d, strongcgcoef=%g, lb=%g, ub=%g -> rhs=%g\n", 
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
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   Real*            strongcgcoef,            /**< array to store strong CG coefficients: must be of size nvars */
   Real*            strongcgrhs,             /**< pointer to store the right hand side of the strong CG row */
   int*             varsign,            /**< stores the sign of the transformed variable in summation */
   int*             boundtype,          /**< stores the bound used for transformed variable (vlb/vub_idx or -1 for lb/ub)*/
   Real             f0,                 /**< fracional value of rhs */
   Real             k,                  /**< factor to strengthen strongcg cut */
   Real*            cutactivity         /**< pointer to store the activity of the resulting cut */
   )
{
   VAR* var;
   Real onedivoneminusf0;
   Real pj;
   Real aj;
   Real downaj;
   Real cutaj;
   Real fj;
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
 *                no strong CG cut found                   , if a'_r <  0
 *
 *  Substitute a°_r * s_r by adding a°_r times the slack's definition to the cut.
 */
static
void substituteStrongCGRow(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< LP data */
   Real*            weights,            /**< row weights in row summation */
   Real             scale,              /**< additional scaling factor multiplied to all rows */
   Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size nvars */
   Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   int*             slacksign,          /**< stores the sign of the row's slack variable in summation */
   Real             f0,                 /**< fracional value of rhs */
   Real             k,                  /**< factor to strengthen strongcg cut */
   Real*            cutactivity         /**< pointer to update the activity of the resulting cut */
   )
{
   ROW* row;
   Real onedivoneminusf0;
   Real pr;
   Real ar;
   Real downar;
   Real cutar;
   Real fr;
   Real mul;
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
      if( row->integral
         && ((slacksign[r] == +1 && SCIPsetIsIntegral(set, row->rhs - row->constant))
            || (slacksign[r] == -1 && SCIPsetIsIntegral(set, row->lhs - row->constant))) )
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
         assert( ar >= 0.0 );
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
      (*cutactivity) += mul * (SCIProwGetLPActivity(row, stat, lp) - row->constant);

      /* move slack's constant to the right hand side */
      if( slacksign[r] == +1 )
      {
         /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a°_r * (rhs - c) to the right hand side */
         assert(!SCIPsetIsInfinity(set, row->rhs));
         (*strongcgrhs) -= cutar * (row->rhs - row->constant);
      }
      else
      {
         /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a°_r * (c - lhs) to the right hand side */
         assert(!SCIPsetIsInfinity(set, -row->lhs));
         (*strongcgrhs) -= cutar * (row->constant - row->lhs);
      }
   }

   /* set rhs to zero, if it's very close to */
   if( SCIPsetIsZero(set, *strongcgrhs) )
      *strongcgrhs = 0.0;
}

/* calculates a strong CG cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a strong CG cut.
 */
RETCODE SCIPlpCalcStrongCG(
   LP*              lp,                 /**< LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   Real             minfrac,            /**< minimal fractionality of rhs to produce strong CG cut for */
   Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   Real             scale,              /**< additional scaling factor multiplied to all rows */
   Real*            strongcgcoef,       /**< array to store strong CG coefficients: must be of size nvars */
   Real*            strongcgrhs,        /**< pointer to store the right hand side of the strong CG row */
   Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   Bool*            success,            /**< pointer to store whether the returned coefficients are a valid strong CG cut */
   Bool*            cutislocal          /**< pointer to store whether the returned cut is only valid locally */
   )
{
   int* slacksign;
   int* varsign;
   int* boundtype;
   Real rhs;
   Real downrhs;
   Real f0;
   Real k;
   Bool emptyrow;
   Bool freevariable;
   Bool localrowsused;
   Bool localbdsused;

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

   debugMessage("calculating strong CG cut (scale: %g)\n", scale);

   /**@todo test, if a column based summation is faster */

   *success = FALSE;
   *cutislocal = FALSE;

   /* allocate temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &slacksign, lp->nrows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &varsign, prob->nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &boundtype, prob->nvars) );

   /* calculate the row summation */
   sumStrongCGRow(set, stat, prob, lp, weights, scale, allowlocal, 
      maxweightrange, strongcgcoef, &rhs, slacksign, &emptyrow, &localrowsused);
   assert(allowlocal || !localrowsused);
   *cutislocal = *cutislocal || localrowsused;
   if( emptyrow )
      goto TERMINATE;
   debug(printMIR(prob, strongcgcoef, rhs));
   
   /* remove all nearly-zero coefficients from strong CG row and relaxes the right hand side correspondingly in order to
    *  prevent numerical rounding errors
    */
   cleanupMIRRow(set, prob, strongcgcoef, &rhs, *cutislocal);
   debug(printMIR(prob, strongcgcoef, rhs));

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
   debug(printMIR(prob, strongcgcoef, rhs));

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
   debug(printMIR(prob, strongcgcoef, *strongcgrhs));

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
   debug(printMIR(prob, strongcgcoef, *strongcgrhs));

   *success = TRUE;
   
#ifndef NDEBUG
   {
      Real act;
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
RETCODE SCIPlpGetState(
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(memhdr != NULL);
   assert(lpistate != NULL);

   CHECK_OKAY( SCIPlpiGetState(lp->lpi, memhdr, lpistate) );

   return SCIP_OKAY;
}

/** loads LP state (like basis information) into solver */
RETCODE SCIPlpSetState(
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{
   assert(lp != NULL);
   assert(memhdr != NULL);
   assert(lpistate != NULL);

   /* flush changes to the LP solver */
   CHECK_OKAY( SCIPlpFlush(lp, memhdr, set) );
   assert(lp->flushed);

   /* set LPI state in the LP solver */
   CHECK_OKAY( SCIPlpiSetState(lp->lpi, memhdr, lpistate) );
   lp->primalfeasible = TRUE;
   lp->dualfeasible = TRUE;

   return SCIP_OKAY;
}

/** frees LP state information */
RETCODE SCIPlpFreeState(
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
   assert(lp != NULL);

   CHECK_OKAY( SCIPlpiFreeState(lp->lpi, memhdr, lpistate) );

   return SCIP_OKAY;
}

/** sets the upper objective limit of the LP solver */
RETCODE SCIPlpSetCutoffbound(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   Real             cutoffbound         /**< new upper objective limit */
   )
{
   assert(lp != NULL);

   debugMessage("setting LP upper objective limit from %g to %g\n", lp->cutoffbound, cutoffbound);
   
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

/** calls LPI to perform primal simplex, measures time and counts iterations, gets basis feasibility status */
static
RETCODE lpPrimalSimplex(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   )
{
   int iterations;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);

   debugMessage("solving primal LP %d (LP %d, %d cols, %d rows)\n", 
      stat->nprimallps+1, stat->nlps+1, lp->ncols, lp->nrows);

#if 0 /*???????????????????????*/
   if( stat->nnodes == 1 && !lp->diving  )
   {
      char fname[MAXSTRLEN];
      sprintf(fname, "lp%lld_%d.lp", stat->nnodes, stat->lpcount);
      CHECK_OKAY( SCIPlpWrite(lp, fname) );
      printf("wrote LP to file <%s> (primal simplex, uobjlim=%g, feastol=%g/%g, fromscratch=%d, fastmip=%d, scaling=%d, presolving=%d)\n", 
         fname, lp->lpiuobjlim, lp->lpifeastol, lp->lpidualfeastol,
         lp->lpifromscratch, lp->lpifastmip, lp->lpiscaling, lp->lpipresolving);
   }
#endif

   /* start timing */
   if( lp->diving )
      SCIPclockStart(stat->divinglptime, set);
   else
      SCIPclockStart(stat->primallptime, set);

   /* call primal simplex */
   CHECK_OKAY( SCIPlpiSolvePrimal(lp->lpi) );
   lp->lastwasprimal = TRUE;

   /* stop timing */
   if( lp->diving )
      SCIPclockStop(stat->divinglptime, set);
   else
      SCIPclockStop(stat->primallptime, set);

   /* count number of iterations */
   stat->lpcount++;
   CHECK_OKAY( SCIPlpGetIterations(lp, &iterations) );
   if( iterations > 0 ) /* don't count the resolves after removing unused columns/rows */
   {
      stat->nlps++;
      stat->nlpiterations += iterations;
      if( !lp->lpifromscratch && stat->nlps > 1 )
      {
         stat->nresolvelps++;
         stat->nresolvelpiterations += iterations;
      }
      if( lp->diving )
      {
         stat->ndivinglps++;
         stat->ndivinglpiterations += iterations;
      }
      else
      {
         stat->nprimallps++;
         stat->nprimallpiterations += iterations;
      }
   }

   debugMessage("solved primal LP %d in %d iterations\n", stat->lpcount, iterations);

   return SCIP_OKAY;
}

/** calls LPI to perform dual simplex, measures time and counts iterations */
static
RETCODE lpDualSimplex(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   )
{
   int iterations;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);

   debugMessage("solving dual LP %d (LP %d, %d cols, %d rows)\n", 
      stat->nduallps+1, stat->nlps+1, lp->ncols, lp->nrows);

#if 0 /*???????????????????????*/
   if( stat->nnodes == 1 && !lp->diving )
   {
      char fname[MAXSTRLEN];
      sprintf(fname, "lp%lld_%d.lp", stat->nnodes, stat->lpcount);
      CHECK_OKAY( SCIPlpWrite(lp, fname) );
      printf("wrote LP to file <%s> (dual simplex, uobjlim=%g, feastol=%g/%g, fromscratch=%d, fastmip=%d, scaling=%d, presolving=%d)\n", 
         fname, lp->lpiuobjlim, lp->lpifeastol, lp->lpidualfeastol, 
         lp->lpifromscratch, lp->lpifastmip, lp->lpiscaling, lp->lpipresolving);
   }
#endif

   /* start timing */
   if( lp->diving )
      SCIPclockStart(stat->divinglptime, set);
   else
      SCIPclockStart(stat->duallptime, set);

   /* call dual simplex */
   CHECK_OKAY( SCIPlpiSolveDual(lp->lpi) );
   lp->lastwasprimal = FALSE;

   /* stop timing */
   if( lp->diving )
      SCIPclockStop(stat->divinglptime, set);
   else
      SCIPclockStop(stat->duallptime, set);

   /* count number of iterations */
   stat->lpcount++;
   CHECK_OKAY( SCIPlpGetIterations(lp, &iterations) );
   if( iterations > 0 ) /* don't count the resolves after removing unused columns/rows */
   {
      stat->nlps++;
      stat->nlpiterations += iterations;
      if( !lp->lpifromscratch && stat->nlps > 1  )
      {
         stat->nresolvelps++;
         stat->nresolvelpiterations += iterations;
      }
      if( lp->diving )
      {
         stat->ndivinglps++;
         stat->ndivinglpiterations += iterations;
      }
      else
      {
         stat->nduallps++;
         stat->nduallpiterations += iterations;
      }
   }

   debugMessage("solved dual LP %d in %d iterations\n", stat->lpcount, iterations);

   return SCIP_OKAY;
}

/** solves the LP with the primal or dual simplex algorithm */
static
RETCODE lpSimplex(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   Bool             useprimal           /**< should the primal simplex be used? */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);

   /* call appropriate simplex */
   if( useprimal )
   {
      CHECK_OKAY( lpPrimalSimplex(lp, set, stat) );
   }
   else
   {
      CHECK_OKAY( lpDualSimplex(lp, set, stat) );
   }
   
   /* check for primal and dual feasibility */
   CHECK_OKAY( SCIPlpiGetBasisFeasibility(lp->lpi, &lp->primalfeasible, &lp->dualfeasible) );

   debugMessage("LP feasibility: primalfeasible=%d, dualfeasible=%d\n", lp->primalfeasible, lp->dualfeasible);

   return SCIP_OKAY;
}

/** solves the LP with the simplex algorithm, and tries to resolve numerical problems */
static
RETCODE lpSolveStable(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   Bool             fastmip,            /**< should the FASTMIP setting of the LP solver be activated? */
   Bool             fromscratch,        /**< should the LP be solved from scratch without using current basis? */
   Bool             useprimal,          /**< should the primal simplex be used? */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   Bool success;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->looseobjvalinf == 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);

   *lperror = FALSE;

   /* solve with given settings (usually fast but unprecise) */
   CHECK_OKAY( lpSetUobjlim(lp, set, lp->cutoffbound - lp->looseobjval) );
   CHECK_OKAY( lpSetFeastol(lp, SCIPsetFeastol(set), &success) );
   CHECK_OKAY( lpSetDualFeastol(lp, SCIPsetDualfeastol(set), &success) );
   CHECK_OKAY( lpSetFromscratch(lp, fromscratch, &success) );
   CHECK_OKAY( lpSetFastmip(lp, fastmip, &success) );
   CHECK_OKAY( lpSetScaling(lp, set->lp_scaling, &success) );
   CHECK_OKAY( lpSetPresolving(lp, set->lp_presolving, &success) );
   CHECK_OKAY( lpSetLPInfo(lp, set->disp_lpinfo) );
   CHECK_OKAY( lpSimplex(lp, set, stat, useprimal) );

   /* check for stability */
   if( SCIPlpiIsStable(lp->lpi) )
      return SCIP_OKAY;

   /* if FASTMIP is turned on, solve again without FASTMIP */
   if( fastmip )
   {
      CHECK_OKAY( lpSetFastmip(lp, FALSE, &success) );
      if( success )
      {
         infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %lld) numerical troubles in LP %d -- solve again without FASTMIP with %s simplex\n", 
            stat->nnodes, stat->nlps, useprimal ? "primal" : "dual");
         CHECK_OKAY( lpSimplex(lp, set, stat, useprimal) );
         
         /* check for stability */
         if( SCIPlpiIsStable(lp->lpi) )
            return SCIP_OKAY;
      }
   }

   /* solve again with opposite scaling setting */
   CHECK_OKAY( lpSetScaling(lp, !set->lp_scaling, &success) );
   if( success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles in LP %d -- solve again with %s simplex %s scaling\n", 
         stat->nnodes, stat->nlps, useprimal ? "primal" : "dual", !set->lp_scaling ? "with" : "without");
      CHECK_OKAY( lpSimplex(lp, set, stat, useprimal) );
   
      /* check for stability */
      if( SCIPlpiIsStable(lp->lpi) )
         return SCIP_OKAY;

      /* reset scaling */
      CHECK_OKAY( lpSetScaling(lp, set->lp_scaling, &success) );
      assert(success);
   }
      
   /* solve again with opposite presolving setting */
   CHECK_OKAY( lpSetPresolving(lp, !set->lp_presolving, &success) );
   if( success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles in LP %d -- solve again with %s simplex %s presolving\n", 
         stat->nnodes, stat->nlps, useprimal ? "primal" : "dual", !set->lp_presolving ? "with" : "without");
      CHECK_OKAY( lpSimplex(lp, set, stat, useprimal) );
   
      /* check for stability */
      if( SCIPlpiIsStable(lp->lpi) )
         return SCIP_OKAY;

      /* reset presolving */
      CHECK_OKAY( lpSetPresolving(lp, set->lp_presolving, &success) );
      assert(success);
   }
      
   /* solve again with a tighter feasibility tolerance */
   CHECK_OKAY( lpSetFeastol(lp, 0.001*SCIPsetFeastol(set), &success) );
   if( success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles in LP %d -- solve again with tighter feasibility tolerance with %s simplex\n", 
         stat->nnodes, stat->nlps, useprimal ? "primal" : "dual");
      CHECK_OKAY( lpSimplex(lp, set, stat, useprimal) );
      
      /* check for stability */
      if( SCIPlpiIsStable(lp->lpi) )
         return SCIP_OKAY;

      /* reset feasibility tolerance */
      CHECK_OKAY( lpSetFeastol(lp, SCIPsetFeastol(set), &success) );
      assert(success);
   }

   /* if not already done, solve again from scratch */
   if( !fromscratch )
   {
      CHECK_OKAY( lpSetFromscratch(lp, TRUE, &success) );
      if( success )
      {
         infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s simplex\n", 
            stat->nnodes, stat->nlps, useprimal ? "primal" : "dual");
         CHECK_OKAY( lpSimplex(lp, set, stat, useprimal) );
         
         /* check for stability */
         if( SCIPlpiIsStable(lp->lpi) )
            return SCIP_OKAY;
      }
   }

   /* solve again, use other simplex this time */
   infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
      "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s simplex\n", 
      stat->nnodes, stat->nlps, !useprimal ? "primal" : "dual");
   CHECK_OKAY( lpSimplex(lp, set, stat, !useprimal) );

   /* check for stability */
   if( SCIPlpiIsStable(lp->lpi) )
      return SCIP_OKAY;

   /* solve again with opposite scaling and other simplex */
   CHECK_OKAY( lpSetScaling(lp, !set->lp_scaling, &success) );
   if( success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s simplex %s scaling\n", 
         stat->nnodes, stat->nlps, !useprimal ? "primal" : "dual", !set->lp_scaling ? "with" : "without");
      CHECK_OKAY( lpSimplex(lp, set, stat, !useprimal) );
      
      /* check for stability */
      if( SCIPlpiIsStable(lp->lpi) )
         return SCIP_OKAY;

      /* reset scaling */
      CHECK_OKAY( lpSetScaling(lp, set->lp_scaling, &success) );
      assert(success);
   }

   /* solve again with opposite presolving and other simplex */
   CHECK_OKAY( lpSetPresolving(lp, !set->lp_presolving, &success) );
   if( success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s simplex %s presolving\n", 
         stat->nnodes, stat->nlps, !useprimal ? "primal" : "dual", !set->lp_presolving ? "with" : "without");
      CHECK_OKAY( lpSimplex(lp, set, stat, !useprimal) );
      
      /* check for stability */
      if( SCIPlpiIsStable(lp->lpi) )
         return SCIP_OKAY;

      /* reset presolving */
      CHECK_OKAY( lpSetPresolving(lp, set->lp_presolving, &success) );
      assert(success);
   }

   /* solve again with tighter feasibility tolerance, use other simplex this time */
   CHECK_OKAY( lpSetFeastol(lp, 0.001*SCIPsetFeastol(set), &success) );
   if( success )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles in LP %d -- solve again from scratch with tighter feasibility tolerance with %s simplex\n", 
         stat->nnodes, stat->nlps, !useprimal ? "primal" : "dual");
      CHECK_OKAY( lpSimplex(lp, set, stat, !useprimal) );

      /* check for stability */
      if( SCIPlpiIsStable(lp->lpi) )
         return SCIP_OKAY;

      /* reset feasibility tolerance */
      CHECK_OKAY( lpSetFeastol(lp, SCIPsetFeastol(set), &success) );
      assert(success);
   }

   /* nothing worked -- store the instable LP to a file and exit with an LPERROR */
   infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "(node %lld) unresolved numerical troubles in LP %d\n", 
      stat->nnodes, stat->nlps);
   *lperror = TRUE;

   return SCIP_OKAY;
}

/** solves the LP with the primal or dual simplex algorithm and evaluates return status */
static
RETCODE lpSolve(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   Bool             fastmip,            /**< should the FASTMIP setting of the LP solver be activated? */
   Bool             fromscratch,        /**< should the LP be solved from scratch without using current basis? */
   Bool             useprimal,          /**< should the primal simplex be used? */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   Bool solvedagain;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);

   checkLinks(lp);

   solvedagain = FALSE;

 SOLVEAGAIN:
   /* call simplex */
   CHECK_OKAY( lpSolveStable(lp, set, stat, fastmip, fromscratch, useprimal, lperror) );

   /* check, if an error occured */
   if( *lperror )
   {
      debugMessage("unresolved error while solving %s LP\n", lp->lastwasprimal ? "primal" : "dual");
      return SCIP_OKAY;
   }

   /* evaluate solution status */
   if( SCIPlpiIsOptimal(lp->lpi) )
   {
      assert(lp->primalfeasible);
      assert(lp->dualfeasible);
      lp->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
      CHECK_OKAY( SCIPlpiGetObjval(lp->lpi, &lp->lpobjval) );
      if( SCIPsetIsRelGE(set, lp->lpobjval, lp->lpiuobjlim) )
      {
         /* the solver may return the optimal value, even if this is greater or equal than the upper bound */
         debugMessage("optimal solution %g exceeds objective limit %g\n", lp->lpobjval, lp->lpiuobjlim);
         lp->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
         lp->lpobjval = SCIPsetInfinity(set);
      }
   }
   else if( SCIPlpiIsObjlimExc(lp->lpi) )
   {
#if 0 /* SOPLEX may return with objective limit reached in any case, because it doesn't distinct btw. primal and dual */
      if( lp->lastwasprimal )
      {
         errorMessage("objective limit exceeded in primal simplex - this should not happen, because no lower limit exists\n");
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
      lp->lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;
      lp->lpobjval = SCIPsetInfinity(set);
   }
   else if( SCIPlpiHasPrimalRay(lp->lpi) )
   {
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
   else if( !solvedagain )
   {
      useprimal = !useprimal;
      solvedagain = TRUE;
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) solution status of LP %d couldn't be proved (internal status:%d) -- solve again with %s simplex\n", 
         stat->nnodes, stat->nlps, useprimal ? "primal" : "dual");
      goto SOLVEAGAIN;
   }
   else
   {
      errorMessage("(node %lld) error or unknown return status of %s simplex in LP %d (internal status: %d)\n", 
         stat->nnodes, lp->lastwasprimal ? "primal" : "dual", stat->nlps, SCIPlpiGetInternalStatus(lp->lpi));
      lp->lpsolstat = SCIP_LPSOLSTAT_ERROR;
      return SCIP_LPERROR;
   }

   lp->solved = TRUE;

   debugMessage("solving %s LP returned solstat=%d (internal status: %d, pfeas=%d, dfeas=%d)\n",
      lp->lastwasprimal ? "primal" : "dual", lp->lpsolstat, SCIPlpiGetInternalStatus(lp->lpi),
      SCIPlpiIsPrimalFeasible(lp->lpi), SCIPlpiIsDualFeasible(lp->lpi));

   return SCIP_OKAY;
}

/** solves the LP with the primal or dual simplex algorithm, depending on the current basis feasibility */
RETCODE SCIPlpSolve(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   Bool             fastmip,            /**< should the FASTMIP setting of the LP solver be activated? */
   Bool             fromscratch,        /**< should the LP be solved from scratch without using current basis? */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   assert(lp != NULL);
   assert(lperror != NULL);

   /* flush changes to the LP solver */
   CHECK_OKAY( SCIPlpFlush(lp, memhdr, set) );
   fastmip = fastmip && !lp->flushaddedcols && !lp->flushdeletedcols; /* turn off FASTMIP if columns were changed */
   
   /* select simplex method */
   if( lp->dualfeasible || !lp->primalfeasible )
   {
      debugMessage("solving dual LP\n");
      CHECK_OKAY( lpSolve(lp, set, stat, fastmip, fromscratch, FALSE, lperror) );
   }
   else
   {
      debugMessage("solving primal LP\n");
      CHECK_OKAY( lpSolve(lp, set, stat, fastmip, fromscratch, TRUE, lperror) );
   }

   return SCIP_OKAY;
}

/** solves the LP with simplex algorithm, and copy the solution into the column's data */
RETCODE SCIPlpSolveAndEval(
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   int              itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   Bool             aging,              /**< should aging and removal of obsolete cols/rows be applied? */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   assert(lp != NULL);
   assert(prob != NULL);
   assert(prob->nvars >= lp->ncols);

   debugMessage("solving LP: %d rows, %d cols, primalfeasible=%d, dualfeasible=%d, solved=%d, diving=%d, cutoff=%g\n", 
      lp->nrows, lp->ncols, lp->primalfeasible, lp->dualfeasible, lp->solved, lp->diving, lp->cutoffbound);

   /* flush changes to the LP solver */
   CHECK_OKAY( SCIPlpFlush(lp, memhdr, set) );

   /* set iteration limit for this solve */
   CHECK_OKAY( lpSetIterationLimit(lp, itlim) );

   if( !lp->solved )
   {
      Bool primalfeasible;
      Bool dualfeasible;
      Bool fastmip;
      Bool fromscratch;

      /* set initial LP solver settings */
      fastmip = set->lp_fastmip && !lp->flushaddedcols && !lp->flushdeletedcols && stat->nnodes > 1;
      fromscratch = FALSE;

   SOLVEAGAIN:
      /* solve the LP */
      CHECK_OKAY( SCIPlpSolve(lp, memhdr, set, stat, fastmip, fromscratch, lperror) );
      debugMessage("SCIPlpSolve() returned solstat %d (error=%d)\n", SCIPlpGetSolstat(lp), *lperror);

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
            CHECK_OKAY( SCIPlpGetSol(lp, memhdr, set, stat, &primalfeasible, &dualfeasible) );
         }
         else
         {
            /* get LP solution believing in the feasibility of the LP solution */
            CHECK_OKAY( SCIPlpGetSol(lp, memhdr, set, stat, NULL, NULL) );
            primalfeasible = TRUE;
            dualfeasible = TRUE;
         }
         if( primalfeasible && dualfeasible && aging && !lp->diving )
         {
            /* update ages and remove obsolete columns and rows from LP */
            CHECK_OKAY( SCIPlpUpdateAges(lp, set, stat) );
            CHECK_OKAY( SCIPlpRemoveNewObsoletes(lp, memhdr, set, stat) );
            
            if( !lp->solved )
            {
               /* resolve LP after removing obsolete columns and rows */
               debugMessage("removed obsoletes - resolve LP again: %d rows, %d cols\n", lp->nrows, lp->ncols);
               aging = FALSE; /* to prevent infinite loops */
               goto SOLVEAGAIN;
            }
         }
         if( !primalfeasible || !dualfeasible )
         {
            if( fastmip )
            {
               /* solution is infeasible (this can happen due to numerical problems): solve again without FASTMIP */
               infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) solution of LP %d not optimal (pfeas=%d, dfeas=%d) -- solving again without FASTMIP\n",
                  stat->nnodes, stat->nlps, primalfeasible, dualfeasible);
               fastmip = FALSE;
               goto SOLVEAGAIN;
            }
            else if( !fromscratch )
            {
               /* solution is infeasible (this can happen due to numerical problems): solve again from scratch */
               infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) solution of LP %d not optimal (pfeas=%d, dfeas=%d) -- solving again from scratch\n",
                  stat->nnodes, stat->nlps, primalfeasible, dualfeasible);
               fromscratch = TRUE;
               goto SOLVEAGAIN;
            }
            else
            {
               warningMessage("(node %lld) unresolved numerical troubles in LP %d\n", stat->nnodes, stat->nlps);
            }
         }
         debugMessage(" -> LP objective value: %g + %g = %g (solstat=%d, cutoff=%g)\n",
            lp->lpobjval, lp->looseobjval, lp->lpobjval + lp->looseobjval, lp->lpsolstat, lp->cutoffbound);
         break;

      case SCIP_LPSOLSTAT_INFEASIBLE:
         if( !SCIPprobAllColsInLP(prob, set, lp) || set->misc_exactsolve )
         {
            CHECK_OKAY( SCIPlpGetDualfarkas(lp, memhdr, set, stat) );
         }
         debugMessage(" -> LP infeasible\n");
         break;

      case SCIP_LPSOLSTAT_UNBOUNDEDRAY:
         CHECK_OKAY( SCIPlpGetUnboundedSol(lp, memhdr, set, stat) );
         debugMessage(" -> LP has unbounded primal ray\n");
         break;

      case SCIP_LPSOLSTAT_OBJLIMIT:
         if( !SCIPprobAllColsInLP(prob, set, lp) || set->misc_exactsolve )
         {
            CHECK_OKAY( SCIPlpGetSol(lp, memhdr, set, stat, NULL, NULL) );
         }
         debugMessage(" -> LP objective limit reached\n");
         break;

      case SCIP_LPSOLSTAT_ITERLIMIT:
         break;

      case SCIP_LPSOLSTAT_TIMELIMIT:
         /**@todo time limit exceeded processing */
         errorMessage("LP time limit exceeded -- case not implemented yet\n");
         return SCIP_ERROR;

      case SCIP_LPSOLSTAT_ERROR:
      case SCIP_LPSOLSTAT_NOTSOLVED:
         errorMessage("error in LP solver\n");
         return SCIP_LPERROR;

      default:
         errorMessage("unknown LP solution status\n");
         return SCIP_ERROR;
      }
   }

   return SCIP_OKAY;
}

/** gets solution status of current LP */
LPSOLSTAT SCIPlpGetSolstat(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->solved || lp->lpsolstat == SCIP_LPSOLSTAT_NOTSOLVED);

   return (lp->flushed ? lp->lpsolstat : SCIP_LPSOLSTAT_NOTSOLVED);
}

/** gets objective value of current LP */
Real SCIPlpGetObjval(
   LP*              lp,                 /**< current LP data */
   SET*             set                 /**< global SCIP settings */
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
Real SCIPlpGetColumnObjval(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->solved);

   return (lp->flushed ? lp->lpobjval : SCIP_INVALID);
}

/** gets part of objective value of current LP that results from LOOSE variables only */
Real SCIPlpGetLooseObjval(
   LP*              lp,                 /**< current LP data */
   SET*             set                 /**< global SCIP settings */
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
Real SCIPlpGetPseudoObjval(
   LP*              lp,                 /**< current LP data */
   SET*             set                 /**< global SCIP settings */
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
Real SCIPlpGetModifiedPseudoObjval(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable */
   Real             oldbound,           /**< old value for bound */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   Real pseudoobjval;
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
Real SCIPlpGetModifiedProvedPseudoObjval(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable */
   Real             oldbound,           /**< old value for bound */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   Real pseudoobjval;
   int pseudoobjvalinf;
   
   pseudoobjval = lp->pseudoobjval;
   pseudoobjvalinf = lp->pseudoobjvalinf;
   if( boundtype == SCIPvarGetBestBoundType(var) )
   {
      INTERVAL obj;
      INTERVAL bd;
      INTERVAL prod;
      INTERVAL psval;

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
RETCODE lpUpdateVar(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   Real             oldobj,             /**< old objective value of variable */
   Real             oldlb,              /**< old objective value of variable */
   Real             oldub,              /**< old objective value of variable */
   Real             newobj,             /**< new objective value of variable */
   Real             newlb,              /**< new objective value of variable */
   Real             newub               /**< new objective value of variable */
   )
{
   Real deltaval;
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
      errorMessage("LP was informed of an objective change of a non-mutable variable\n");
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

   return SCIP_OKAY;
}

/** updates current pseudo and loose objective values for a change in a variable's objective value or bounds;
 *  pseudo objective value is calculated with interval arithmetics to get a proved lower bound
 */
static
RETCODE lpUpdateVarProved(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   Real             oldobj,             /**< old objective value of variable */
   Real             oldlb,              /**< old objective value of variable */
   Real             oldub,              /**< old objective value of variable */
   Real             newobj,             /**< new objective value of variable */
   Real             newlb,              /**< new objective value of variable */
   Real             newub               /**< new objective value of variable */
   )
{
   INTERVAL deltaval;
   INTERVAL bd;
   INTERVAL obj;
   INTERVAL prod;
   INTERVAL psval;
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
      errorMessage("LP was informed of an objective change of a non-mutable variable\n");
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
RETCODE SCIPlpUpdateVarObj(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   Real             oldobj,             /**< old objective value of variable */
   Real             newobj              /**< new objective value of variable */
   )
{
   assert(set != NULL);
   assert(var != NULL);

   if( set->misc_exactsolve )
   {
      if( oldobj != newobj )
      {
         CHECK_OKAY( lpUpdateVarProved(lp, set, var, oldobj, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
                        newobj, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) );
      }
   }
   else
   {
      if( !SCIPsetIsEQ(set, oldobj, newobj) )
      {
         CHECK_OKAY( lpUpdateVar(lp, set, var, oldobj, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
                        newobj, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) );
      }
   }
   
   return SCIP_OKAY;
}

/** updates current pseudo and loose objective value for a change in a variable's lower bound */
RETCODE SCIPlpUpdateVarLb(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   Real             oldlb,              /**< old lower bound of variable */
   Real             newlb               /**< new lower bound of variable */
   )
{
   assert(set != NULL);
   assert(var != NULL);

   if( set->misc_exactsolve )
   {
      if( oldlb != newlb && SCIPvarGetObj(var) > 0.0 )
      {
         CHECK_OKAY( lpUpdateVarProved(lp, set, var, SCIPvarGetObj(var), oldlb, SCIPvarGetUbLocal(var), 
                        SCIPvarGetObj(var), newlb, SCIPvarGetUbLocal(var)) );
      }
   }
   else
   {
      if( !SCIPsetIsEQ(set, oldlb, newlb) && SCIPsetIsPositive(set, SCIPvarGetObj(var)) )
      {
         CHECK_OKAY( lpUpdateVar(lp, set, var, SCIPvarGetObj(var), oldlb, SCIPvarGetUbLocal(var), 
                        SCIPvarGetObj(var), newlb, SCIPvarGetUbLocal(var)) );
      }
   }
   
   return SCIP_OKAY;
}

/** updates current pseudo objective value for a change in a variable's upper bound */
RETCODE SCIPlpUpdateVarUb(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   Real             oldub,              /**< old upper bound of variable */
   Real             newub               /**< new upper bound of variable */
   )
{
   assert(set != NULL);
   assert(var != NULL);

   if( set->misc_exactsolve )
   {
      if( oldub != newub && SCIPvarGetObj(var) < 0.0 )
      {
         CHECK_OKAY( lpUpdateVarProved(lp, set, var, SCIPvarGetObj(var), SCIPvarGetLbLocal(var), oldub, 
                        SCIPvarGetObj(var), SCIPvarGetLbLocal(var), newub) );
      }
   }
   else
   {
      if( !SCIPsetIsEQ(set, oldub, newub) && SCIPsetIsNegative(set, SCIPvarGetObj(var)) )
      {
         CHECK_OKAY( lpUpdateVar(lp, set, var, SCIPvarGetObj(var), SCIPvarGetLbLocal(var), oldub, 
                        SCIPvarGetObj(var), SCIPvarGetLbLocal(var), newub) );
      }
   }

   return SCIP_OKAY;
}

/** informs LP, that given variable was added to the problem */
RETCODE SCIPlpUpdateAddVar(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< variable that is now a LOOSE problem variable */
   )
{
   assert(lp != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);

   /* add the variable to the loose objective value sum */
   CHECK_OKAY( SCIPlpUpdateVarObj(lp, set, var, 0.0, SCIPvarGetObj(var)) );

   /* update the loose variables counter */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      lp->nloosevars++;

   return SCIP_OKAY;
}

/** informs LP, that given formerly loose problem variable is now a column variable */
static
RETCODE lpUpdateVarColumn(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   )
{
   Real obj;
   Real lb;
   Real ub;

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
RETCODE lpUpdateVarColumnProved(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   )
{
   INTERVAL bd;
   INTERVAL ob;
   INTERVAL prod;
   INTERVAL loose;
   Real obj;
   Real lb;
   Real ub;

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
RETCODE SCIPlpUpdateVarColumn(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   )
{
   assert(set != NULL);

   if( set->misc_exactsolve )
   {
      CHECK_OKAY( lpUpdateVarColumnProved(lp, set, var) );
   }
   else
   {
      CHECK_OKAY( lpUpdateVarColumn(lp, set, var) );
   }

   return SCIP_OKAY;
}

/** informs LP, that given formerly column problem variable is now again a loose variable */
static
RETCODE lpUpdateVarLoose(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   )
{
   Real obj;
   Real lb;
   Real ub;

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
RETCODE lpUpdateVarLooseProved(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   )
{
   INTERVAL bd;
   INTERVAL ob;
   INTERVAL prod;
   INTERVAL loose;
   Real obj;
   Real lb;
   Real ub;

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
RETCODE SCIPlpUpdateVarLoose(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   )
{
   assert(set != NULL);

   if( set->misc_exactsolve )
   {
      CHECK_OKAY( lpUpdateVarLooseProved(lp, set, var) );
   }
   else
   {
      CHECK_OKAY( lpUpdateVarLoose(lp, set, var) );
   }

   return SCIP_OKAY;
}

/** stores the LP solution in the columns and rows */
RETCODE SCIPlpGetSol(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   Bool*            dualfeasible        /**< pointer to store whether the solution is dual feasible, or NULL */
   )
{
   COL** lpicols;
   ROW** lpirows;
   Real* primsol;
   Real* dualsol;
   Real* activity;
   Real* redcost;
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

   debugMessage("getting new LP solution %d\n", stat->lpcount);

   /* get temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &primsol, lp->nlpicols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &dualsol, lp->nlpirows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &activity, lp->nlpirows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &redcost, lp->nlpicols) );
   
   CHECK_OKAY( SCIPlpiGetSol(lp->lpi, NULL, primsol, dualsol, activity, redcost) );

   lpicols = lp->lpicols;
   lpirows = lp->lpirows;
   nlpicols = lp->nlpicols;
   nlpirows = lp->nlpirows;
   lpcount = stat->lpcount;

   /* copy primal solution and reduced costs into columns */
   for( c = 0; c < nlpicols; ++c )
   {
      lpicols[c]->primsol = primsol[c];
      if( primsol[c] < lpicols[c]->minprimsol )
         lpicols[c]->minprimsol = primsol[c];
      if( primsol[c] > lpicols[c]->maxprimsol )
         lpicols[c]->maxprimsol = primsol[c];
      lpicols[c]->redcost = redcost[c];
      lpicols[c]->validredcostlp = lpcount;
      if( primalfeasible != NULL )
         *primalfeasible = *primalfeasible
            && SCIPsetIsFeasGE(set, lpicols[c]->primsol, lpicols[c]->lb)
            && SCIPsetIsFeasLE(set, lpicols[c]->primsol, lpicols[c]->ub);
      if( dualfeasible != NULL )
      {
         if( SCIPsetIsGT(set, lpicols[c]->primsol, lpicols[c]->lb) )
            *dualfeasible = *dualfeasible && !SCIPsetIsFeasPositive(set, lpicols[c]->redcost);
         if( SCIPsetIsLT(set, lpicols[c]->primsol, lpicols[c]->ub) )
            *dualfeasible = *dualfeasible && !SCIPsetIsFeasNegative(set, lpicols[c]->redcost);
      }
      debugMessage(" col <%s> [%g,%g]: primsol=%.9f, redcost=%.9f, pfeas=%d/%d, dfeas=%d\n",
         SCIPvarGetName(lpicols[c]->var), lpicols[c]->lb, lpicols[c]->ub, lpicols[c]->primsol, lpicols[c]->redcost,
         SCIPsetIsFeasGE(set, lpicols[c]->primsol, lpicols[c]->lb),
         SCIPsetIsFeasLE(set, lpicols[c]->primsol, lpicols[c]->ub),
         !SCIPsetIsFeasNegative(set, lpicols[c]->redcost));
   }

   /* copy dual solution and activities into rows */
   for( r = 0; r < nlpirows; ++r )
   {
      lpirows[r]->dualsol = dualsol[r];
      lpirows[r]->activity = activity[r] + lpirows[r]->constant;
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
      debugMessage(" row <%s> [%g,%g]: dualsol=%.9f, activity=%.9f, pfeas=%d/%d, dfeas=%d/%d\n", 
         lpirows[r]->name, lpirows[r]->lhs, lpirows[r]->rhs, lpirows[r]->dualsol, lpirows[r]->activity,
         SCIPsetIsFeasGE(set, lpirows[r]->activity, lpirows[r]->lhs),
         SCIPsetIsFeasLE(set, lpirows[r]->activity, lpirows[r]->rhs),
         SCIPsetIsInfinity(set, -lpirows[r]->lhs) ? !SCIPsetIsFeasPositive(set, lpirows[r]->dualsol) : TRUE,
         SCIPsetIsInfinity(set, lpirows[r]->rhs) ? !SCIPsetIsFeasNegative(set, lpirows[r]->dualsol) : TRUE);
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &redcost);
   SCIPsetFreeBufferArray(set, &activity);
   SCIPsetFreeBufferArray(set, &dualsol);
   SCIPsetFreeBufferArray(set, &primsol);

   return SCIP_OKAY;
}

/** stores LP solution with infinite objective value in the columns and rows */
RETCODE SCIPlpGetUnboundedSol(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   )
{
   COL** lpicols;
   ROW** lpirows;
   Real* primsol;
   Real* activity;
   Real* ray;
   Real rayobjval;
   Real rayscale;
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

   debugMessage("getting new unbounded LP solution %d\n", stat->lpcount);

   /* get temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &primsol, lp->nlpicols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &activity, lp->nlpirows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &ray, lp->nlpicols) );

   /* get primal feasible point */
   CHECK_OKAY( SCIPlpiGetSol(lp->lpi, NULL, primsol, NULL, activity, NULL) );

   /* get primal unbounded ray */
   CHECK_OKAY( SCIPlpiGetPrimalRay(lp->lpi, ray) );
   
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
   debugMessage("unbounded LP solution: rayobjval=%f, rayscale=%f\n", rayobjval, rayscale);

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
RETCODE SCIPlpGetDualfarkas(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   )
{
   ROW** lpirows;
   Real* dualfarkas;
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
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &dualfarkas, lp->nlpirows) );

   /* get dual farkas infeasibility proof */
   CHECK_OKAY( SCIPlpiGetDualfarkas(lp->lpi, dualfarkas) );

   lpirows = lp->lpirows;
   nlpirows = lp->nlpirows;

   /* store infeasibility proof in rows */
   debugMessage("LP is infeasible:\n");
   for( r = 0; r < nlpirows; ++r )
   {
      /*debugMessage(" row <%s>: dualfarkas=%f\n", lpirows[r]->name, dualfarkas[r]);*/
      lpirows[r]->dualfarkas = dualfarkas[r];
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &dualfarkas);

   return SCIP_OKAY;
}

/** get number of iterations used in last LP solve */
RETCODE SCIPlpGetIterations(
   LP*              lp,                 /**< current LP data */
   int*             iterations          /**< pointer to store the iteration count */
   )
{
   assert(lp != NULL);

   CHECK_OKAY( SCIPlpiGetIterations(lp->lpi, iterations) );

   return SCIP_OKAY;
}

/** increases age of columns with solution value 0.0 and rows with activity not at its bounds,
 *  resets age of non-zero columns and sharp rows
 */
RETCODE SCIPlpUpdateAges(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   )
{
   COL** lpicols;
   ROW** lpirows;
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

   debugMessage("updating LP ages\n");

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
RETCODE lpDelColset(
   LP*              lp,                 /**< current LP data */
   int*             coldstat            /**< deletion status of columns:  1 if column should be deleted, 0 if not */
   )
{
   COL* col;
   int ncols;
   int c;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->ncols == lp->nlpicols);
   assert(!lp->diving);
   assert(coldstat != NULL);

   ncols = lp->ncols;

   /* delete columns in LP solver */
   CHECK_OKAY( SCIPlpiDelColset(lp->lpi, coldstat) );

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
RETCODE lpDelRowset(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   int*             rowdstat            /**< deletion status of rows:  1 if row should be deleted, 0 if not */
   )
{
   ROW* row;
   int nrows;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->nrows == lp->nlpirows);
   assert(!lp->diving);
   assert(rowdstat != NULL);

   nrows = lp->nrows;

   /* delete rows in LP solver */
   CHECK_OKAY( SCIPlpiDelRowset(lp->lpi, rowdstat) );

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

         CHECK_OKAY( SCIProwRelease(&lp->lpirows[r], memhdr, set, lp) );
         SCIProwUnlock(lp->rows[r]);
         CHECK_OKAY( SCIProwRelease(&lp->rows[r], memhdr, set, lp) );
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

/** removes all columns, that are too old, beginning with the given firstcol */
static
RETCODE lpRemoveObsoleteCols(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   int              firstcol            /**< first column to check for clean up */
   )
{
   COL** cols;
   COL** lpicols;
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

   if( lp->nremoveablecols == 0 || set->lp_colagelimit == -1 )
      return SCIP_OKAY;

   ncols = lp->ncols;
   cols = lp->cols;
   lpicols = lp->lpicols;

   /* get temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &coldstat, ncols) );

   /* mark obsolete columns to be deleted */
   ndelcols = 0;
   clearMemoryArray(coldstat, ncols);
   for( c = firstcol; c < ncols; ++c )
   {
      assert(cols[c] == lpicols[c]);
      assert(cols[c]->lppos == c);
      assert(cols[c]->lpipos == c);
      if( cols[c]->removeable
         && cols[c]->obsoletenode != stat->nnodes /* don't remove column a second time from same node (avoid cycling) */
         && cols[c]->age > set->lp_colagelimit
         && SCIPsetIsZero(set, SCIPcolGetBestBound(cols[c])) ) /* bestbd != 0 -> column would be priced in next time */
      {
         coldstat[c] = 1;
         ndelcols++;
         cols[c]->obsoletenode = stat->nnodes;
         debugMessage("removing obsolete col <%s>: primsol=%f, bounds=[%g,%g]\n", 
            SCIPvarGetName(cols[c]->var), cols[c]->primsol, cols[c]->lb, cols[c]->ub);
      }
   }

   debugMessage("removing %d/%d obsolete columns from LP\n", ndelcols, ncols);

   /* delete the marked columns in the LP solver interface, update the LP respectively */
   if( ndelcols > 0 )
   {
      CHECK_OKAY( lpDelColset(lp, coldstat) );
   }
   assert(lp->ncols == ncols - ndelcols);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &coldstat);
      
   return SCIP_OKAY;
}

/** removes all rows, that are too old, beginning with the given firstrow */
static
RETCODE lpRemoveObsoleteRows(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   int              firstrow            /**< first row to check for clean up */
   )
{
   ROW** rows;
   ROW** lpirows;
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

   if( lp->nremoveablerows == 0 || set->lp_rowagelimit == -1 )
      return SCIP_OKAY;

   nrows = lp->nrows;
   rows = lp->rows;
   lpirows = lp->lpirows;

   /* get temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &rowdstat, nrows) );

   /* mark obsolete rows to be deleted */
   ndelrows = 0;
   clearMemoryArray(rowdstat, nrows);
   for( r = firstrow; r < nrows; ++r )
   {
      assert(rows[r] == lpirows[r]);
      assert(rows[r]->lppos == r);
      assert(rows[r]->lpipos == r);
      if( rows[r]->removeable
         && rows[r]->obsoletenode != stat->nnodes  /* don't remove row a second time from same node (avoid cycling) */
         && rows[r]->age > set->lp_rowagelimit )
      {
         rowdstat[r] = 1;
         ndelrows++;
         rows[r]->obsoletenode = stat->nnodes;
         debugMessage("removing obsolete row <%s>: activity=%f, sides=[%g,%g]\n", 
            rows[r]->name, rows[r]->activity, rows[r]->lhs, rows[r]->rhs);
      }
   }

   debugMessage("removing %d/%d obsolete rows from LP\n", ndelrows, nrows);

   /* delete the marked rows in the LP solver interface, update the LP respectively */
   if( ndelrows > 0 )
   {
      CHECK_OKAY( lpDelRowset(lp, memhdr, set, rowdstat) );
   }
   assert(lp->nrows == nrows - ndelrows);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &rowdstat);
      
   return SCIP_OKAY;
}

/** removes all columns and rows in the part of the LP created at the current node, that are too old */
RETCODE SCIPlpRemoveNewObsoletes(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);
   assert(set != NULL);

   debugMessage("removing obsolete columns starting with %d/%d, obsolete rows starting with %d/%d\n",
      lp->firstnewcol, lp->ncols, lp->firstnewrow, lp->nrows);

   if( lp->firstnewcol < lp->ncols )
   {
      CHECK_OKAY( lpRemoveObsoleteCols(lp, set, stat, lp->firstnewcol) );
   }
   if( lp->firstnewrow < lp->nrows )
   {
      CHECK_OKAY( lpRemoveObsoleteRows(lp, memhdr, set, stat, lp->firstnewrow) );
   }

   return SCIP_OKAY;
}

/** removes all columns and rows in whole LP, that are too old */
RETCODE SCIPlpRemoveAllObsoletes(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);
   assert(set != NULL);

   debugMessage("removing all obsolete columns and rows\n");

   if( 0 < lp->ncols )
   {
      CHECK_OKAY( lpRemoveObsoleteCols(lp, set, stat, 0) );
   }
   if( 0 < lp->nrows )
   {
      CHECK_OKAY( lpRemoveObsoleteRows(lp, memhdr, set, stat, 0) );
   }

   return SCIP_OKAY;
}

/** removes all columns at 0.0 beginning with the given firstcol */
static
RETCODE lpCleanupCols(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   int              firstcol            /**< first column to check for clean up */
   )
{
   COL** cols;
   COL** lpicols;
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

   if( lp->nremoveablecols == 0 )
      return SCIP_OKAY;

   ncols = lp->ncols;
   cols = lp->cols;
   lpicols = lp->lpicols;

   /* get temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &coldstat, ncols) );

   /* mark unused columns to be deleted */
   ndelcols = 0;
   clearMemoryArray(coldstat, ncols);
   for( c = firstcol; c < ncols; ++c )
   {
      assert(cols[c] == lpicols[c]);
      assert(cols[c]->lppos == c);
      assert(cols[c]->lpipos == c);
      if( lpicols[c]->removeable
         && lpicols[c]->primsol == 0.0 /* non-basic columns to remove are exactly at 0.0 */
         && SCIPsetIsZero(set, SCIPcolGetBestBound(cols[c])) ) /* bestbd != 0 -> column would be priced in next time */
      {
         coldstat[c] = 1;
         ndelcols++;
      }
   }

   debugMessage("removing %d/%d unused columns from LP\n", ndelcols, ncols);

   /* delete the marked columns in the LP solver interface, update the LP respectively */
   if( ndelcols > 0 )
   {
      CHECK_OKAY( lpDelColset(lp, coldstat) );
   }
   assert(lp->ncols == ncols - ndelcols);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &coldstat);
      
   return SCIP_OKAY;
}

/** removes all rows not at one of their bounds beginning with the given firstrow */
static
RETCODE lpCleanupRows(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   int              firstrow            /**< first row to check for clean up */
   )
{
   ROW** rows;
   ROW** lpirows;
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

   if( lp->nremoveablerows == 0 )
      return SCIP_OKAY;

   nrows = lp->nrows;
   rows = lp->rows;
   lpirows = lp->lpirows;

   /* get temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &rowdstat, nrows) );

   /* mark unused rows to be deleted */
   ndelrows = 0;
   clearMemoryArray(rowdstat, nrows);
   for( r = firstrow; r < nrows; ++r )
   {
      assert(rows[r] == lpirows[r]);
      assert(rows[r]->lppos == r);
      assert(rows[r]->lpipos == r);
      if( lpirows[r]->removeable && lpirows[r]->dualsol == 0.0 ) /* basic rows to remove are exactly at 0.0 */
      {
         rowdstat[r] = 1;
         ndelrows++;
      }
   }

   debugMessage("removing %d/%d unused rows from LP\n", ndelrows, nrows);

   /* delete the marked rows in the LP solver interface, update the LP respectively */
   if( ndelrows > 0 )
   {
      CHECK_OKAY( lpDelRowset(lp, memhdr, set, rowdstat) );
   }
   assert(lp->nrows == nrows - ndelrows);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &rowdstat);

   return SCIP_OKAY;
}

/** removes all columns at 0.0 and rows not at their bound in the part of the LP created at the current node */
RETCODE SCIPlpCleanupNew(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert(!lp->diving);
   assert(set != NULL);

   debugMessage("removing unused columns starting with %d/%d (%d), unused rows starting with %d/%d (%d)\n",
      lp->firstnewcol, lp->ncols, set->lp_cleanupcols, lp->firstnewrow, lp->nrows, set->lp_cleanuprows);

   if( set->lp_cleanupcols && lp->firstnewcol < lp->ncols )
   {
      CHECK_OKAY( lpCleanupCols(lp, set, stat, lp->firstnewcol) );
   }
   if( set->lp_cleanuprows && lp->firstnewrow < lp->nrows )
   {
      CHECK_OKAY( lpCleanupRows(lp, memhdr, set, stat, lp->firstnewrow) );
   }

   return SCIP_OKAY;
}

/** removes all columns at 0.0 and rows not at their bound in the whole LP */
RETCODE SCIPlpCleanupAll(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert(!lp->diving);
   assert(set != NULL);

   debugMessage("removing all unused columns and rows\n");

   if( /*set->lp_cleanupcols &&*/ 0 < lp->ncols )
   {
      CHECK_OKAY( lpCleanupCols(lp, set, stat, 0) );
   }
   if( /*set->lp_cleanuprows &&*/ 0 < lp->nrows )
   {
      CHECK_OKAY( lpCleanupRows(lp, memhdr, set, stat, 0) );
   }

   return SCIP_OKAY;
}

/** initiates LP diving */
RETCODE SCIPlpStartDive(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert(lp->flushed);
   assert(!lp->diving);
   assert(lp->divelpistate == NULL);
   assert(set != NULL);

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
   CHECK_OKAY( SCIPlpiGetState(lp->lpi, memhdr, &lp->divelpistate) );

   /* switch to diving mode */
   lp->diving = TRUE;

   return SCIP_OKAY;
}

/** quits LP diving and resets bounds and objective values of columns to the current node's values */
RETCODE SCIPlpEndDive(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   VAR**            vars,               /**< array with all active variables */
   int              nvars               /**< number of active variables */
   )
{
   VAR* var;
   Bool lperror;
   int v;

   assert(lp != NULL);
   assert(lp->diving);
   assert(lp->divelpistate != NULL);
   assert(nvars == 0 || vars != NULL);

   /* reset all columns' objective values and bounds to its original values */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      {
         CHECK_OKAY( SCIPcolChgObj(SCIPvarGetCol(var), set, lp, SCIPvarGetObj(var)) );
         CHECK_OKAY( SCIPcolChgLb(SCIPvarGetCol(var), set, lp, SCIPvarGetLbLocal(var)) );
         CHECK_OKAY( SCIPcolChgUb(SCIPvarGetCol(var), set, lp, SCIPvarGetUbLocal(var)) );
      }
   }

   /* reload LPI state saved at start of diving, free LPI state afterwards */
   CHECK_OKAY( SCIPlpSetState(lp, memhdr, set, lp->divelpistate) );
   CHECK_OKAY( SCIPlpFreeState(lp, memhdr, &lp->divelpistate) );
   assert(lp->divelpistate == NULL);

   /**@todo Get rid of resolving after diving: use separate data fields in columns to store all diving
    *       information (bounds, obj, solution values) and create calls SCIPvarGetDiveSol() etc.
    *       to access this diving LP information (not applicable on LOOSE variables).
    *       Just declare the LP to be solved at this point (remember the LP solution status beforehand).
    */
   /* resolve LP to reset solution */
   CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, -1, FALSE, &lperror) );
   if( lperror )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) unresolved numerical troubles while resolving LP %d after diving\n", stat->nnodes, stat->nlps);
   }

   /* switch to standard (non-diving) mode and remember the diving node */
   lp->diving = FALSE;
   lp->divingobjchg = FALSE;
   stat->lastdivenode = stat->nnodes;

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
RETCODE provedBound(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   Bool             usefarkas,          /**< use y = dual farkas and c = 0 instead of y = dual solution and c = obj? */
   Real*            bound               /**< result of interval arithmetic minimization */
   )
{
   INTERVAL* yinter;
   INTERVAL b;
   INTERVAL ytb;
   INTERVAL prod;
   INTERVAL diff;
   INTERVAL x;
   INTERVAL minprod;
   INTERVAL a;
   ROW* row;
   COL* col;
   Real y;
   Real c;
   int i;
   int j;

   assert(lp != NULL);
   assert(lp->solved);
   assert(set != NULL);
   assert(bound != NULL);

   /* allocate buffer for storing y in interval arithmetic */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &yinter, lp->nrows) );

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
RETCODE SCIPlpGetProvedLowerbound(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   Real*            bound               /**< pointer to store proven dual bound */
   )
{
   CHECK_OKAY( provedBound(lp, set, FALSE, bound) );

   debugMessage("proved lower bound of LP: %g\n", *bound);

   return SCIP_OKAY;
}

/** gets proven dual bound of last LP solution */
RETCODE SCIPlpIsInfeasibilityProved(
   LP*              lp,                 /**< current LP data */
   SET*             set,                /**< global SCIP settings */
   Bool*            proved              /**< pointer to store whether infeasibility is proven */
   )
{
   Real bound;

   assert(proved != NULL);

   CHECK_OKAY( provedBound(lp, set, TRUE, &bound) );

   *proved = (bound > 0.0);

   debugMessage("proved farkas value of LP: %g -> infeasibility %sproved\n", bound, *proved ? "" : "not ");

   return SCIP_OKAY;
}

/** writes LP to a file */
RETCODE SCIPlpWrite(
   LP*              lp,                 /**< current LP data */
   const char*      fname               /**< file name */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(fname != NULL);

   CHECK_OKAY( SCIPlpiWriteLP(lp->lpi, fname) );

#if 0 /*???????????????????????????*/
   {
      FILE* f;
      int c;

      printf("storing integrality conditions in <%s>\n", fname);
      f = fopen(fname, "a");
      fprintf(f, "General\n");
      for( c = 0; c < lp->ncols; ++c )
      {
         assert(lp->cols[c] == lp->lpicols[c]);
         if( SCIPvarGetType(SCIPcolGetVar(lp->cols[c])) != SCIP_VARTYPE_CONTINUOUS )
            fprintf(f, "%s\n", SCIPvarGetName(SCIPcolGetVar(lp->cols[c])));
      }
      fprintf(f, "End\n");
      fclose(f);
   }
#endif

   return SCIP_OKAY;
}


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets array with columns of the LP */
COL** SCIPlpGetCols(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->cols;
}

/** gets current number of columns in LP */
int SCIPlpGetNCols(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->ncols;
}

/** gets array with rows of the LP */
ROW** SCIPlpGetRows(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->rows;
}

/** gets current number of rows in LP */
int SCIPlpGetNRows(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->nrows;
}

/** gets array with newly added columns after the last mark */
COL** SCIPlpGetNewcols(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewcol && lp->firstnewcol <= lp->ncols);

   return &(lp->cols[lp->firstnewcol]);
}

/** gets number of newly added columns after the last mark */
int SCIPlpGetNNewcols(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewcol && lp->firstnewcol <= lp->ncols);

   return lp->ncols - lp->firstnewcol;
}

/** gets array with newly added rows after the last mark */
ROW** SCIPlpGetNewrows(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewrow && lp->firstnewrow <= lp->nrows);

   return &(lp->rows[lp->firstnewrow]);
}

/** gets number of newly added rows after the last mark */
int SCIPlpGetNNewrows(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewrow && lp->firstnewrow <= lp->nrows);

   return lp->nrows - lp->firstnewrow;
}

/** gets the LP solver interface */
LPI* SCIPlpGetLPI(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->lpi;
}

/** returns whether the LP is in diving mode */
Bool SCIPlpDiving(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->diving;
}

/** returns whether the LP is in diving mode and the objective value of at least one column was changed */
Bool SCIPlpDivingObjChanged(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->divingobjchg;
}

/** marks the diving LP to have a changed objective function */
void SCIPlpMarkDivingObjChanged(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->diving);

   lp->divingobjchg = TRUE;
}

#endif
