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
#pragma ident "@(#) $Id: lp.c,v 1.119 2004/05/14 13:43:54 bzfpfend Exp $"

/**@file   lp.c
 * @brief  LP management methods and datastructures
 * @author Tobias Achterberg
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
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &row->cols_probindex, row->size, newsize) );
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
   int* probindex;
   int* linkpos;
   COL* tmpcol;
   Real tmpval;
   int tmpprobindex;
   int tmplinkpos;
   int tmpindex;
   int pos;
   int sortpos;

   assert(row != NULL);
   assert(0 <= firstpos && firstpos <= lastpos+1 && lastpos < row->len);

   /**@todo do a quick sort here, if many elements are unsorted (sorted-Bool -> sorted-Int?) */
   cols = row->cols;
   vals = row->vals;
   probindex = row->cols_probindex;
   linkpos = row->linkpos;

   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && cols[pos]->index <= cols[pos+1]->index )
            pos++;
         if( pos >= lastpos )
            break;
         assert(cols[pos]->index > cols[pos+1]->index);
         tmpcol = cols[pos];
         tmpprobindex = probindex[pos];
         tmpval = vals[pos];
         tmplinkpos = linkpos[pos];
         tmpindex = tmpcol->index;
         do
         {
            cols[pos] = cols[pos+1];
            probindex[pos] = probindex[pos+1];
            vals[pos] = vals[pos+1];
            linkpos[pos] = linkpos[pos+1];
            pos++;
         }
         while( pos < lastpos && cols[pos+1]->index < tmpindex );
         cols[pos] = tmpcol;
         probindex[pos] = tmpprobindex;
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
         while( pos > firstpos && cols[pos-1]->index <= cols[pos]->index )
            pos--;
         if( pos <= firstpos )
            break;
         assert(cols[pos-1]->index > cols[pos]->index);
         tmpcol = cols[pos];
         tmpprobindex = probindex[pos];
         tmpval = vals[pos];
         tmplinkpos = linkpos[pos];
         tmpindex = tmpcol->index;
         do
         {
            cols[pos] = cols[pos-1];
            probindex[pos] = probindex[pos-1];
            vals[pos] = vals[pos-1];
            linkpos[pos] = linkpos[pos-1];
            pos--;
         }
         while( pos > firstpos && cols[pos-1]->index > tmpindex );
         cols[pos] = tmpcol;
         probindex[pos] = tmpprobindex;
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
int colSearchCoeffPart(
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

      pos = colSearchCoeffPart(col, row, 0, col->nlprows-1);
      if( pos >= 0 )
         return pos;
   }

   /* search in the non-LP/unlinked rows */
   if( row->lppos == -1 || col->nunlinked > 0 )
   {
      /* column has to be sorted, such that binary search works */
      colSortNonLP(col);
      assert(col->nonlprowssorted);

      pos = colSearchCoeffPart(col, row, col->nlprows, col->len-1);
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
      idx = row->cols[pos]->index;
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
void colSwapCoeffs(
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
   row->cols_probindex[newpos] = row->cols_probindex[oldpos];
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
void rowSwapCoeffs(
   ROW*             row,                /**< LP row */
   int              pos1,               /**< position of first coefficient */
   int              pos2                /**< position of second coefficient */
   )
{
   COL* tmpcol;
   Real tmpval;
   int tmpprobindex;
   int tmplinkpos;
   
   assert(row != NULL);
   assert(0 <= pos1 && pos1 < row->len);
   assert(0 <= pos2 && pos2 < row->len);
   assert(row->cols[pos1] != NULL);

   if( pos1 == pos2 )
      return;

   /* swap coefficients */
   tmpcol = row->cols[pos2];
   tmpprobindex = row->cols_probindex[pos2];
   tmpval = row->vals[pos2];
   tmplinkpos = row->linkpos[pos2];

   row->cols[pos2] = row->cols[pos1];
   row->cols_probindex[pos2] = row->cols_probindex[pos1];
   row->vals[pos2] = row->vals[pos1];
   row->linkpos[pos2] = row->linkpos[pos1];

   row->cols[pos1] = tmpcol;
   row->cols_probindex[pos1] = tmpprobindex;
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

   assert(lp != NULL);

   if( !msgdisp )
   {
      warningMessage("LP LINK CHECKING ACTIVATED! THIS IS VERY SLOW!\n");
      msgdisp = TRUE;
   }

   for( i = 0; i < lp->ncols; ++i )
   {
      col = lp->cols[i];
      assert(col != NULL);
      assert(col->lppos >= 0 || col->primsol == 0.0);
      assert(col->lppos >= 0 || col->farkas == 0.0);
      assert(col->nlprows <= col->len);

      for( j = 0; j < col->len; ++j )
      {
         row = col->rows[j];
         assert(row != NULL);
         assert(!lp->flushed || col->lppos == -1 || col->linkpos[j] >= 0);
         assert(col->linkpos[j] == -1 || row->cols[col->linkpos[j]] == col);
         assert(col->linkpos[j] == -1 || EPSEQ(row->vals[col->linkpos[j]], col->vals[j], 1e-6));
         assert((j < col->nlprows) == (col->linkpos[j] >= 0 && row->lppos >= 0));
      }
   }

   for( i = 0; i < lp->nrows; ++i )
   {
      row = lp->rows[i];
      assert(row != NULL);
      assert(row->lppos >= 0 || row->dualsol == 0.0);
      assert(row->lppos >= 0 || row->dualfarkas == 0.0);
      assert(row->nlpcols <= row->len);
      
      for( j = 0; j < row->len; ++j )
      {
         col = row->cols[j];
         assert(col != NULL);
         assert(!lp->flushed || row->lppos == -1 || row->linkpos[j] >= 0);
         assert(row->linkpos[j] == -1 || col->rows[row->linkpos[j]] == row);
         assert(row->linkpos[j] == -1 || EPSEQ(col->vals[row->linkpos[j]], row->vals[j], 1e-6));
         assert((j < row->nlpcols) == (row->linkpos[j] >= 0 && col->lppos >= 0));
      }
   }
}
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
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->primalfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   row->pseudoactivity = SCIP_INVALID;
   row->minactivity = SCIP_INVALID;
   row->maxactivity = SCIP_INVALID;
   row->validpsactivitybdchg = -1;
   row->validactivitybdsbdchg = -1;
}




/*
 * local column changing methods
 */

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
      col->nunlinked++;
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
         rowSwapCoeffs(row, linkpos, row->nlpcols-1);
      }
   }

   /* update the sorted flags */
   if( row->lppos >= 0 && linkpos >= 0 )
   {
      if( col->nlprows > 1 )
         col->lprowssorted = col->lprowssorted && (col->rows[col->nlprows-2]->index < row->index);
   }
   else
   {
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
RETCODE colDelCoeffPos(
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
RETCODE colChgCoeffPos(
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   int              pos,                /**< position in column vector to change */
   Real             val                 /**< value of coefficient */
   )
{
   assert(memhdr != NULL);
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
      CHECK_OKAY( colDelCoeffPos(col, set, lp, pos) );
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

/** update row norms after addition of new coefficient */
static
void rowAddNorms(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   int              colidx,             /**< column index of new coefficient, or -1 */
   Real             val                 /**< value of new coefficient */
   )
{
   Real absval;

   assert(row != NULL);
   assert(row->nummaxval >= 0);
   assert(row->numminval >= 0);
   assert(set != NULL);
   assert(colidx >= -1);

   absval = ABS(val);
   assert(!SCIPsetIsZero(set, absval));

   /* update min/maxidx */
   if( colidx >= 0 )
   {
      row->minidx = MIN(row->minidx, colidx);
      row->maxidx = MAX(row->maxidx, colidx);
   }

   /* update squared euclidean norm */
   row->sqrnorm += SQR(absval);

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
   int              colidx,             /**< column index of deleted coefficient, or -1 */
   Real             val                 /**< value of deleted coefficient */
   )
{
   Real absval;

   assert(row != NULL);
   assert(row->nummaxval >= 0);
   assert(row->numminval >= 0);
   assert(set != NULL);
   assert(colidx >= -1);

   absval = ABS(val);
   assert(!SCIPsetIsZero(set, absval));
   assert(row->nummaxval == 0 || SCIPsetIsGE(set, row->maxval, absval));
   assert(row->numminval == 0 || SCIPsetIsLE(set, row->minval, absval));

   /* update min/maxidx validity */
   if( colidx >= 0 )
   {
      if( colidx == row->minidx || colidx == row->maxidx )
         row->validminmaxidx = FALSE;
   }

   /* update squared euclidean norm */
   row->sqrnorm -= SQR(absval);
   row->sqrnorm = MAX(row->sqrnorm, 0.0);

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
   row->cols_probindex[pos] = col->var_probindex;
   row->vals[pos] = val;
   row->linkpos[pos] = linkpos;
   row->integral = row->integral && SCIPcolIsIntegral(col) && SCIPsetIsIntegral(set, val);
   if( linkpos == -1 )
      row->nunlinked++;
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
         colSwapCoeffs(col, linkpos, col->nlprows-1);
      }
   }

   /* update the sorted flags */
   if( col->lppos >= 0 && linkpos >= 0 )
   {
      if( row->nlpcols > 1 )
         row->lpcolssorted = row->lpcolssorted && (row->cols[row->nlpcols-2]->index < col->index);
   }
   else
   {
      if( row->len - row->nlpcols > 1 )
         row->nonlpcolssorted = row->nonlpcolssorted && (row->cols[row->len-2]->index < col->index);
   }
   
   rowAddNorms(row, set, col->index, val);

   coefChanged(row, col, lp);

   debugMessage("added coefficient %g * <%s> at position %d (%d/%d) to row <%s> (nunlinked=%d)\n",
      val, SCIPvarGetName(col->var), pos, row->nlpcols, row->len, row->name, row->nunlinked);

   return SCIP_OKAY;
}

/** deletes coefficient at given position from row */
static
RETCODE rowDelCoeffPos(
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

   rowDelNorms(row, set, col->index, val);

   coefChanged(row, col, lp);

   return SCIP_OKAY;
}

/** changes a coefficient at given position of an LP row */
static
RETCODE rowChgCoeffPos(
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   int              pos,                /**< position in row vector to change */
   Real             val                 /**< value of coefficient */
   )
{
   assert(memhdr != NULL);
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
      CHECK_OKAY( rowDelCoeffPos(row, set, lp, pos) );
   }
   else if( !SCIPsetIsEQ(set, row->vals[pos], val) )
   {
      /* change existing coefficient */
      rowDelNorms(row, set, -1, row->vals[pos]);
      row->vals[pos] = val;
      row->integral = row->integral && SCIPsetIsIntegral(set, val);
      rowAddNorms(row, set, -1, row->vals[pos]);
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
         errorMessage("Unknown row side type\n");
         abort();
      }
      
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

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
            CHECK_OKAY( rowDelCoeffPos(col->rows[i], set, lp, col->linkpos[i]) );
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
            CHECK_OKAY( colDelCoeffPos(row->cols[i], set, lp, row->linkpos[i]) );
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
   int              value               /**< value to set parameter to */
   )
{
   RETCODE retcode;

   assert(lp != NULL);

   retcode = SCIPlpiSetIntpar(lp->lpi, lpparam, value);

   /* ignore unknown parameter error */
   if( retcode == SCIP_PARAMETERUNKNOWN )
      return SCIP_OKAY;

   return retcode;
}

/** sets parameter of type Real in LP solver, ignoring unknown parameters */
static
RETCODE lpSetRealpar(
   LP*              lp,                 /**< current LP data */
   LPPARAM          lpparam,            /**< LP parameter */
   Real             value               /**< value to set parameter to */
   )
{
   RETCODE retcode;

   assert(lp != NULL);

   retcode = SCIPlpiSetRealpar(lp->lpi, lpparam, value);

   /* ignore unknown parameter error */
   if( retcode == SCIP_PARAMETERUNKNOWN )
      return SCIP_OKAY;

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
   if( set->exactsolve )
      return SCIP_OKAY;

   CHECK_OKAY( lpCheckRealpar(lp, SCIP_LPPAR_UOBJLIM, lp->lpiuobjlim) );

   if( uobjlim != lp->lpiuobjlim )
   {
      CHECK_OKAY( lpSetRealpar(lp, SCIP_LPPAR_UOBJLIM, uobjlim) );
      lp->solved = FALSE;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
      lp->primalfeasible = FALSE;
      lp->lpiuobjlim = uobjlim;
   }

   return SCIP_OKAY;
}

/** sets the feasibility tolerance of the LP solver */
static
RETCODE lpSetFeastol(
   LP*              lp,                 /**< current LP data */
   Real             feastol             /**< new feasibility tolerance */
   )
{
   assert(lp != NULL);
   assert(feastol >= 0.0);
   
   CHECK_OKAY( lpCheckRealpar(lp, SCIP_LPPAR_FEASTOL, lp->lpifeastol) );

   if( feastol != lp->lpifeastol )
   {
      CHECK_OKAY( lpSetRealpar(lp, SCIP_LPPAR_FEASTOL, feastol) );
      if( lp->nrows > 0 && feastol < lp->lpifeastol )
      {
         lp->solved = FALSE;
         lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         lp->primalfeasible = FALSE;
      }
      lp->lpifeastol = feastol;
   }

   return SCIP_OKAY;
}

/** sets the reduced costs feasibility tolerance of the LP solver */
static
RETCODE lpSetDualFeastol(
   LP*              lp,                 /**< current LP data */
   Real             dualfeastol         /**< new reduced costs feasibility tolerance */
   )
{
   assert(lp != NULL);
   assert(dualfeastol >= 0.0);

   CHECK_OKAY( lpCheckRealpar(lp, SCIP_LPPAR_DUALFEASTOL, lp->lpidualfeastol) );

   if( dualfeastol != lp->lpidualfeastol )
   {
      CHECK_OKAY( lpSetRealpar(lp, SCIP_LPPAR_DUALFEASTOL, dualfeastol) );
      if( lp->nrows > 0 && dualfeastol < lp->lpidualfeastol )
      {
         lp->solved = FALSE;
         lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         lp->primalfeasible = FALSE;
      }
      lp->lpidualfeastol = dualfeastol;
   }

   return SCIP_OKAY;
}

/** sets the FROMSCRATCH setting of the LP solver */
static
RETCODE lpSetFromscratch(
   LP*              lp,                 /**< current LP data */
   Bool             fromscratch         /**< new FROMSCRATCH setting */
   )
{
   assert(lp != NULL);

   CHECK_OKAY( lpCheckIntpar(lp, SCIP_LPPAR_FROMSCRATCH, lp->lpifromscratch) );

   if( fromscratch != lp->lpifromscratch )
   {
      CHECK_OKAY( lpSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, fromscratch) );
      lp->lpifromscratch = fromscratch;
   }
   
   return SCIP_OKAY;
}

/** sets the FASTMIP setting of the LP solver */
static
RETCODE lpSetFastmip(
   LP*              lp,                 /**< current LP data */
   Bool             fastmip             /**< new FASTMIP setting */
   )
{
   assert(lp != NULL);

   CHECK_OKAY( lpCheckIntpar(lp, SCIP_LPPAR_FASTMIP, lp->lpifastmip) );

   if( fastmip != lp->lpifastmip )
   {
      CHECK_OKAY( lpSetIntpar(lp, SCIP_LPPAR_FASTMIP, fastmip) );
      lp->lpifastmip = fastmip;
   }
   
   return SCIP_OKAY;
}

/** sets the SCALING setting of the LP solver */
static
RETCODE lpSetScaling(
   LP*              lp,                 /**< current LP data */
   Bool             scaling             /**< new SCALING setting */
   )
{
   assert(lp != NULL);

   CHECK_OKAY( lpCheckIntpar(lp, SCIP_LPPAR_SCALING, lp->lpiscaling) );

   if( scaling != lp->lpiscaling )
   {
      CHECK_OKAY( lpSetIntpar(lp, SCIP_LPPAR_SCALING, scaling) );
      lp->lpiscaling = scaling;
   }
   
   return SCIP_OKAY;
}

/** sets the iteration limit of the LP solver */
static
RETCODE lpSetIterationLimit(
   LP*              lp,                 /**< current LP data */
   int              itlim               /**< maximal number of LP iterations to perform, or -1 for no limit */
   )
{
   assert(lp != NULL);
   assert(itlim >= -1);

   if( itlim == -1 )
      itlim = INT_MAX;

   CHECK_OKAY( lpCheckIntpar(lp, SCIP_LPPAR_LPITLIM, lp->lpiitlim) );

   if( itlim != lp->lpiitlim )
   {
      if( itlim > lp->lpiitlim )
         lp->solved = FALSE;

      CHECK_OKAY( lpSetIntpar(lp, SCIP_LPPAR_LPITLIM, itlim) );
      lp->lpiitlim = itlim;
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
   (*col)->index = stat->ncolidx++;
   (*col)->size = len;
   (*col)->len = len;
   (*col)->nlprows = 0;
   (*col)->nunlinked = len;
   (*col)->lppos = -1;
   (*col)->lpipos = -1;
   (*col)->primsol = 0.0;
   (*col)->redcost = SCIP_INVALID;
   (*col)->farkas = SCIP_INVALID;
   (*col)->minprimsol = (*col)->ub;
   (*col)->maxprimsol = (*col)->lb;
   (*col)->strongbranchdown = SCIP_INVALID;
   (*col)->strongbranchup = SCIP_INVALID;
   (*col)->strongbranchsolval  = SCIP_INVALID;
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
      assert(row->cols_probindex[col->linkpos[pos]] == col->var_probindex);
      assert(SCIPsetIsEQ(set, row->vals[col->linkpos[pos]], col->vals[pos]));
      CHECK_OKAY( rowDelCoeffPos(row, set, lp, col->linkpos[pos]) );
   }

   /* delete the row from the column's row vector */
   CHECK_OKAY( colDelCoeffPos(col, set, lp, pos) );
   
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
         assert(row->cols_probindex[col->linkpos[pos]] == col->var_probindex);
         assert(SCIPsetIsEQ(set, row->vals[col->linkpos[pos]], col->vals[pos]));
         CHECK_OKAY( rowChgCoeffPos(row, memhdr, set, lp, col->linkpos[pos], val) );
      }

      /* change the coefficient in the column */
      CHECK_OKAY( colChgCoeffPos(col, memhdr, set, lp, pos, val) );
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
         assert(row->cols_probindex[col->linkpos[pos]] == col->var_probindex);
         assert(SCIPsetIsEQ(set, row->vals[col->linkpos[pos]], col->vals[pos]));
         CHECK_OKAY( rowChgCoeffPos(row, memhdr, set, lp, col->linkpos[pos], col->vals[pos] + incval) );
      }

      /* change the coefficient in the column */
      CHECK_OKAY( colChgCoeffPos(col, memhdr, set, lp, pos, col->vals[pos] + incval) );
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
      
      /* invalidate LP solution */
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

      assert(lp->nchgcols > 0);
   }  

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
      
      /* invalidate LP solution */
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

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
      
      /* invalidate LP solution */
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

      assert(lp->nchgcols > 0);
   }  

   col->ub = newub;

   return SCIP_OKAY;
}

/** calculates the reduced costs of a column */
static
void colCalcRedcost(
   COL*             col,                /**< LP column */
   STAT*            stat                /**< problem statistics */
   )
{
   ROW* row;
   int r;

   assert(col != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(stat != NULL);

   col->redcost = col->obj;
   for( r = 0; r < col->nlprows; ++r )
   {
      row = col->rows[r];
      assert(row != NULL);
      assert(row->dualsol < SCIP_INVALID);
      assert(row->lppos >= 0);
      assert(col->linkpos[r] >= 0);
      col->redcost -= col->vals[r] * row->dualsol;
   }

   if( col->nunlinked > 0 )
   {
      for( r = col->nlprows; r < col->len; ++r )
      {
         row = col->rows[r];
         assert(row != NULL);
         assert(row->lppos >= 0 || row->dualsol == 0.0);
         assert(row->lppos == -1 || col->linkpos[r] == -1);
         if( row->lppos >= 0 )
            col->redcost -= col->vals[r] * row->dualsol;
      }
   }
#ifndef NDEBUG
   else
   {
      for( r = col->nlprows; r < col->len; ++r )
      {
         row = col->rows[r];
         assert(row != NULL);
         assert(row->dualsol == 0.0);
         assert(row->lppos == -1);
         assert(col->linkpos[r] >= 0);
      }
   }
#endif

   col->validredcostlp = stat->lpcount;
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
      colCalcRedcost(col, stat);
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
      return set->infinity;
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
         return -ABS(redcost);
      }
      else
      {
         /* dual row is  activity >= obj  <=>  redcost <= 0 */
         return -redcost;
      }
   }
}

/** calculates the farkas value of a column */
static
void colCalcFarkas(
   COL*             col,                /**< LP column */
   STAT*            stat                /**< problem statistics */
   )
{
   ROW* row;
   int r;

   assert(col != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(stat != NULL);

   col->farkas = 0.0;
   for( r = 0; r < col->nlprows; ++r )
   {
      row = col->rows[r];
      assert(row != NULL);
      assert(row->dualfarkas < SCIP_INVALID);
      assert(row->lppos >= 0);
      assert(col->linkpos[r] >= 0);
      col->farkas += col->vals[r] * row->dualfarkas;
   }

   if( col->nunlinked > 0 )
   {
      for( r = col->nlprows; r < col->len; ++r )
      {
         row = col->rows[r];
         assert(row != NULL);
         assert(row->lppos >= 0 || row->dualfarkas == 0.0);
         assert(row->lppos == -1 || col->linkpos[r] == -1);
         if( row->lppos >= 0 )
            col->farkas += col->vals[r] * row->dualfarkas;
      }
   }
#ifndef NDEBUG
   else
   {
      for( r = col->nlprows; r < col->len; ++r )
      {
         row = col->rows[r];
         assert(row != NULL);
         assert(row->dualfarkas == 0.0);
         assert(row->lppos == -1);
         assert(col->linkpos[r] >= 0);
      }
   }
#endif

   if( col->farkas > 0.0 )
      col->farkas *= col->ub;
   else
      col->farkas *= col->lb;

   col->validfarkaslp = stat->lpcount;
}

/** gets the farkas value of a column in last LP (which must be infeasible) */
Real SCIPcolGetFarkas(
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
      colCalcFarkas(col, stat);
   assert(col->validfarkaslp == stat->lpcount);
   assert(col->farkas < SCIP_INVALID);

   return col->farkas;
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
   assert(stat != NULL);
   assert(lp != NULL);
   assert(lperror != NULL);

   debugMessage("performing strong branching on variable <%s>(%g) with %d iterations\n", 
      SCIPvarGetName(col->var), col->primsol, itlim);

   /* start timing */
   SCIPclockStart(stat->strongbranchtime, set);
      
   /* call LPI strong branching */
   stat->nstrongbranchs++;
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
      col->strongbranchnode = -1;
   }
   else
   {
      CHECK_OKAY( retcode );
      col->strongbranchdown = MIN(strongbranchdown + lp->looseobjval, lp->cutoffbound);
      col->strongbranchup = MIN(strongbranchup + lp->looseobjval, lp->cutoffbound);
            
      /* update strong branching statistics */
      if( iter == -1 )
      {
         /* calculate avergate iteration number */
         iter = stat->nlps > 0 ? (int)(2*stat->nlpiterations / stat->nlps) : 0;
         if( iter/2 >= itlim )
            iter = 2*itlim;
      }
      stat->nsblpiterations += iter;
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
      col->strongbranchnode = stat->nnodes;

      /* if a loose variables has an infinite best bound, the LP bound is -infinity and no gain can be achieved */
      if( lp->looseobjvalinf > 0 )
      {
         col->strongbranchdown = -set->infinity;
         col->strongbranchup = -set->infinity;
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
   Real*            solval              /**< stores LP solution value of column at last strong branching call, or NULL */
   )
{
   assert(col != NULL);

   if( down != NULL )
      *down = col->strongbranchdown;
   if( up != NULL )
      *up = col->strongbranchup;
   if( solval != NULL )
      *solval = col->strongbranchsolval;
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
   fprintf(file, "[%f,%f], ", col->lb, col->ub);

   /* print coefficients */
   if( col->len == 0 )
      fprintf(file, "<empty>");
   for( r = 0; r < col->len; ++r )
   {
      assert(col->rows[r] != NULL);
      assert(col->rows[r]->name != NULL);
      fprintf(file, "%+f%s ", col->vals[r], col->rows[r]->name);
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

   return col->lppos;
}

/** returns TRUE iff column is member of current LP */
Bool SCIPcolIsInLP(
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

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
   row->maxval = 0.0;
   row->nummaxval = 1;
   row->minval = set->infinity;
   row->numminval = 1;
   row->minidx = INT_MAX;
   row->maxidx = INT_MIN;
   row->validminmaxidx = TRUE;
   row->lpcolssorted = TRUE;
   row->nonlpcolssorted = TRUE;

   /* check, if row is sorted
    * calculate sqrnorm, maxval, minval, minidx, and maxidx
    */
   for( i = 0; i < row->nlpcols; ++i )
   {
      assert(row->cols[i] != NULL);
      assert(!SCIPsetIsZero(set, row->vals[i]));
      assert(row->cols[i]->lppos >= 0);
      assert(row->linkpos[i] >= 0);

      rowAddNorms(row, set, row->cols[i]->index, row->vals[i]);
      if( i > 0 )
         row->lpcolssorted = row->lpcolssorted && (row->cols[i-1]->index < row->cols[i]->index);
   }
   for( i = row->nlpcols; i < row->len; ++i )
   {
      assert(row->cols[i] != NULL);
      assert(!SCIPsetIsZero(set, row->vals[i]));
      assert(row->cols[i]->lppos == -1 || row->linkpos[i] == -1);

      rowAddNorms(row, set, row->cols[i]->index, row->vals[i]);
      if( i > row->nlpcols )
         row->nonlpcolssorted = row->nonlpcolssorted && (row->cols[i-1]->index < row->cols[i]->index);
   }
}

/** scales row with given factor, and rounds coefficients to integers if close enough;
 *  the constant is automatically moved to the sides
 */
static
RETCODE rowScale(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   Real             scaleval,           /**< value to scale row with */
   Real             roundtol            /**< rounding tolerance, upto which values are rounded to next integer */
   )
{
   COL* col;
   Real val;
   Real newval;
   int pos;
   int c;

   assert(row != NULL);
   assert(row->len == 0 || row->cols != NULL);
   assert(row->len == 0 || row->vals != NULL);
   assert(SCIPsetIsPositive(set, scaleval));
   assert(!SCIPsetIsNegative(set, roundtol));

   debugMessage("scale row <%s> with %g (tolerance=%g)\n", row->name, scaleval, roundtol);

   /* scale the row coefficients, thereby recalculating whether the row's activity is always integral */
   row->integral = TRUE;
   for( c = 0; c < row->len; ++c )
   {
      col = row->cols[c];
      val = row->vals[c];
      assert(!SCIPsetIsZero(set, val));

      newval = val * scaleval;
      if( EPSISINT(newval, roundtol) )
         newval = EPSFLOOR(newval, roundtol);
      assert(!SCIPsetIsZero(set, newval));

      row->vals[c] = newval;
      row->integral = row->integral && SCIPcolIsIntegral(col) && SCIPsetIsIntegral(set, newval);

      /* update the norms of the row */
      rowDelNorms(row, set, -1, val);
      rowAddNorms(row, set, -1, newval);

      /* update the value in the corresponding column vector, if already linked */
      pos = row->linkpos[c];
      if( pos >= 0 )
      {
         assert(col->rows != NULL);
         assert(col->vals != NULL);
         assert(col->rows[pos] == row);
         assert(SCIPsetIsEQ(set, col->vals[pos], val));
         col->vals[pos] = newval;
      }

      /* mark the coefficient changed */
      coefChanged(row, col, lp);
   }

   /* scale the row sides, and move the constant to the sides */
   if( !SCIPsetIsInfinity(set, -row->lhs) )
   {
      newval = (row->lhs - row->constant) * scaleval;
      if( EPSISINT(newval, roundtol) )
         newval = EPSFLOOR(newval, roundtol);
      CHECK_OKAY( SCIProwChgLhs(row, set, lp, newval) );
   }
   if( !SCIPsetIsInfinity(set, row->rhs) )
   {
      newval = (row->rhs - row->constant) * scaleval;
      if( EPSISINT(newval, roundtol) )
         newval = EPSFLOOR(newval, roundtol);
      CHECK_OKAY( SCIProwChgRhs(row, set, lp, newval) );
   }

   /* clear the row constant */
   CHECK_OKAY( SCIProwChgConstant(row, set, stat, lp, 0.0) );

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
      ALLOC_OKAY( allocBlockMemoryArray(memhdr, &(*row)->cols_probindex, len) );
      ALLOC_OKAY( allocBlockMemoryArray(memhdr, &(*row)->linkpos, len) );

      for( i = 0; i < len; ++i )
      {
         assert(cols[i] != NULL);
         assert(!SCIPsetIsZero(set, vals[i]));

         var = cols[i]->var;
         assert(cols[i]->var_probindex == SCIPvarGetProbindex(var));

         (*row)->cols_probindex[i] = cols[i]->var_probindex;
         (*row)->linkpos[i] = -1;
         (*row)->integral = (*row)->integral && SCIPvarIsIntegral(var) && SCIPsetIsIntegral(set, vals[i]);
      }
   }
   else
   {
      (*row)->cols = NULL;
      (*row)->cols_probindex = NULL;
      (*row)->vals = NULL;
      (*row)->linkpos = NULL;
   }
   
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*row)->name, name, strlen(name)+1) );
   (*row)->constant = 0.0;
   (*row)->lhs = lhs;
   (*row)->rhs = rhs;
   (*row)->sqrnorm = 0.0;
   (*row)->maxval = 0.0;
   (*row)->minval = set->infinity;
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
   (*row)->minidx = INT_MAX;
   (*row)->maxidx = INT_MIN;
   (*row)->nummaxval = 0;
   (*row)->numminval = 0;
   (*row)->validactivitylp = -1;
   (*row)->validpsactivitybdchg = -1;
   (*row)->validactivitybdsbdchg = -1;
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
   freeBlockMemoryArrayNull(memhdr, &(*row)->cols_probindex, (*row)->size);
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

/** locks an unmodifiable row, which forbids further changes */
RETCODE SCIProwLock(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   debugMessage("lock row <%s> with nuses=%d and nlocks=%d\n", row->name, row->nuses, row->nlocks);

   /* check, if row is modifiable */
   if( row->modifiable )
   {
      errorMessage("cannot lock the modifiable row <%s>\n", row->name);
      return SCIP_INVALIDDATA;
   }
   
   row->nlocks++;

   return SCIP_OKAY;
}

/** unlocks a lock of a row; a row with no sealed lock may be modified */
RETCODE SCIProwUnlock(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   debugMessage("unlock row <%s> with nuses=%d and nlocks=%d\n", row->name, row->nuses, row->nlocks);

   /* check, if row is modifiable */
   if( row->modifiable )
   {
      errorMessage("cannot unlock the modifiable row <%s>\n", row->name);
      return SCIP_INVALIDDATA;
   }
   
   /* check, if row is locked */
   if( row->nlocks == 0 )
   {
      errorMessage("row <%s> has no sealed lock\n", row->name);
      return SCIP_INVALIDDATA;
   }

   row->nlocks--;

   return SCIP_OKAY;
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
   assert(row->cols_probindex[pos] == col->var_probindex);

   /* if column knows of the row, remove the row from the column's row vector */
   if( row->linkpos[pos] >= 0 )
   {
      assert(col->rows[row->linkpos[pos]] == row);
      assert(SCIPsetIsEQ(set, col->vals[row->linkpos[pos]], row->vals[pos]));
      CHECK_OKAY( colDelCoeffPos(col, set, lp, row->linkpos[pos]) );
   }

   /* delete the column from the row's col vector */
   CHECK_OKAY( rowDelCoeffPos(row, set, lp, pos) );
   
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
      assert(row->cols_probindex[pos] == col->var_probindex);

      /* if column knows of the row, change the corresponding coefficient in the column */
      if( row->linkpos[pos] >= 0 )
      {
         assert(col->rows[row->linkpos[pos]] == row);
         assert(SCIPsetIsEQ(set, col->vals[row->linkpos[pos]], row->vals[pos]));
         CHECK_OKAY( colChgCoeffPos(col, memhdr, set, lp, row->linkpos[pos], val) );
      }

      /* change the coefficient in the row */
      CHECK_OKAY( rowChgCoeffPos(row, memhdr, set, lp, pos, val) );
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
      assert(row->cols_probindex[pos] == col->var_probindex);

      /* if column knows of the row, change the corresponding coefficient in the column */
      if( row->linkpos[pos] >= 0 )
      {
         assert(col->rows[row->linkpos[pos]] == row);
         assert(SCIPsetIsEQ(set, col->vals[row->linkpos[pos]], row->vals[pos]));
         CHECK_OKAY( colChgCoeffPos(col, memhdr, set, lp, row->linkpos[pos], row->vals[pos] + incval) );
      }

      /* change the coefficient in the row */
      CHECK_OKAY( rowChgCoeffPos(row, memhdr, set, lp, pos, row->vals[pos] + incval) );
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
   assert(!SCIPsetIsInfinity(set, ABS(constant)));
   assert(stat != NULL);
   assert(lp != NULL);
   assert(!lp->diving);

   if( !SCIPsetIsEQ(set, constant, row->constant) )
   {
      if( row->validpsactivitybdchg == stat->nboundchanges )
      {
         assert(row->pseudoactivity < SCIP_INVALID);
         row->pseudoactivity += constant - row->constant;
      }
      if( row->validactivitybdsbdchg == stat->nboundchanges )
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
   assert(!SCIPsetIsInfinity(set, ABS(addval)));
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



#define DIVTOL      (1e+06*set->epsilon)
#define TWOMULTTOL  (1e+03*set->epsilon)
#define RATIONALTOL (1e+02*set->epsilon)
/** tries to find a rational representation of the row and multiplies coefficients with common denominator */
RETCODE SCIProwMakeRational(
   ROW*             row,                /**< LP row */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   Real             maxscale,           /**< maximal value to scale row with */
   Bool*            success             /**< stores whether row could be made rational */
   )
{
   COL* col;
   Real val;
   Real minval;
   Real maxval;
   Real usedtol;
   Bool contvars;
   Bool fractional;
   int c;

   assert(row != NULL);
   assert(row->len == 0 || row->cols != NULL);
   assert(row->len == 0 || row->cols_probindex != NULL);
   assert(row->len == 0 || row->vals != NULL);
   assert(maxdnom >= 1);
   assert(success != NULL);

   *success = FALSE;

   /* nothing to do, if row is empty */
   if( row->len == 0 )
   {
      *success = TRUE;
      return SCIP_OKAY;
   }

   /* get minimal and maximal non-zero coefficient of row */
   minval = SCIProwGetMinval(row, set);
   maxval = SCIProwGetMaxval(row, set);
   assert(SCIPsetIsPositive(set, minval));
   assert(SCIPsetIsPositive(set, maxval));

   /* check, if there are fractional coefficients and continuous variables in the row */
   contvars = FALSE;
   fractional = FALSE;
   for( c = 0; c < row->len; ++c )
   {
      col = row->cols[c];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(col->var) == col);
      val = row->vals[c];
      assert(!SCIPsetIsZero(set, val));
      
      contvars = contvars || !SCIPcolIsIntegral(col);
      fractional = fractional || !SCIPsetIsIntegral(set, val);
   }

   /* if fractional coefficients exist, try to find a rational representation */
   if( fractional )
   {
      Bool scalable;
      Bool twomult;
      Real scaleval;
      Real twomultval;

      /* try, if row coefficients can be made integral by 
       *  - multiplying them with the reciprocal of the smallest coefficient and a power of 2
       *  - by multiplying them by a power of 2
       */
      scalable = TRUE;
      scaleval = 1.0/minval;
      twomult = TRUE;
      twomultval = 1.0;
      for( c = 0; c < row->len && (scalable || twomult); ++c )
      {
         val = row->vals[c];
         if( scalable )
         {
            while( scaleval <= maxscale && !EPSISINT(val * scaleval, DIVTOL) )
               scaleval *= 2.0;
            scalable = (scaleval <= maxscale);
         }
         if( twomult )
         {
            while( twomultval <= maxscale && !EPSISINT(val * twomultval, TWOMULTTOL) )
               twomultval *= 2.0;
            twomult = (twomultval <= maxscale);
         }
      }

      if( scalable )
      {
         /* make row coefficients integral by dividing them by the smallest coefficient */
         assert(scaleval <= maxscale);
         CHECK_OKAY( rowScale(row, set, stat, lp, scaleval, DIVTOL) );
         *success = TRUE;
      }
      else if( twomult )
      {
         /* make row coefficients integral by multiplying them with a power of 2 */
         assert(twomultval <= maxscale);
         CHECK_OKAY( rowScale(row, set, stat, lp, twomultval, TWOMULTTOL) );
         *success = TRUE;
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
         if( row->len > 0 )
         {
            /* first coefficient (to initialize gcd) */
            val = row->vals[0];
            rational = SCIPrealToRational(val, RATIONALTOL, maxdnom, &nominator, &denominator);
            if( rational )
            {
               assert(denominator > 0);
               gcd = (nominator == 0 ? 1 : ABS(nominator));
               scm = denominator;
               rational = ((Real)scm/(Real)gcd <= maxscale);
            }
            /* remaining coefficients */
            for( c = 1; c < row->len && rational; ++c )
            {
               val = row->vals[c];
               rational = SCIPrealToRational(val, RATIONALTOL, maxdnom, &nominator, &denominator);
               if( rational )
               {
                  assert(denominator > 0);
                  if( nominator != 0 )
                     gcd = SCIPcalcGreComDiv(gcd, ABS(nominator));
                  scm *= denominator / SCIPcalcGreComDiv(scm, denominator);
                  rational = ((Real)scm/(Real)gcd <= maxscale);
               }
            }
         }

#if 0 /*??????????????????????*/
         /* convert the sides into a rational number, calculate the greatest common divisor of the nominators
          * and the smallest common multiple of the denominators
          */
         if( rational && !SCIPsetIsInfinity(set, -row->lhs) )
         {
            rational = SCIPrealToRational(row->lhs - row->constant, RATIONALTOL, maxdnom, &nominator, &denominator);
            if( rational && nominator != 0  )
            {
               assert(denominator > 0);
               if( nominator != 0 )
                  gcd = SCIPcalcGreComDiv(gcd, ABS(nominator));
               scm *= denominator / SCIPcalcGreComDiv(scm, denominator);
               rational = ((Real)scm/(Real)gcd <= maxscale);
            }
         }
         if( rational && !SCIPsetIsInfinity(set, row->rhs) )
         {
            rational = SCIPrealToRational(row->rhs - row->constant, RATIONALTOL, maxdnom, &nominator, &denominator);
            if( rational && nominator != 0 )
            {
               assert(denominator > 0);
               if( nominator != 0 )
                  gcd = SCIPcalcGreComDiv(gcd, ABS(nominator));
               scm *= denominator / SCIPcalcGreComDiv(scm, denominator);
               rational = ((Real)scm/(Real)gcd <= maxscale);
            }
         }
#endif

         if( rational )
         {
            /* make row coefficients integral by multiplying them with the smallest common multiple of the denominators */
            assert((Real)scm/(Real)gcd <= maxscale);
            CHECK_OKAY( rowScale(row, set, stat, lp, (Real)scm/(Real)gcd, RATIONALTOL) );
            *success = TRUE;
         }
      }
   }
   else
   {
      /* all coefficients are integral: we have nothing to do except moving the constant to the sides */
      CHECK_OKAY( SCIProwChgLhs(row, set, lp, row->lhs - row->constant) );
      CHECK_OKAY( SCIProwChgRhs(row, set, lp, row->rhs - row->constant) );
      CHECK_OKAY( SCIProwChgConstant(row, set, stat, lp, 0.0) );
      *success = TRUE;
   }

   /* clean up the row sides */
   if( *success )
   {
      assert(row->constant == 0.0); /* in rowScale(), the constant should be moved to the sides */
      if( !contvars )
      {
         /* no continuous variables exist in the row, all coefficients of the new row are integral -> round sides */
         if( !SCIPsetIsInfinity(set, -row->lhs) )
            row->lhs = SCIPsetCeil(set, row->lhs);
         if( !SCIPsetIsInfinity(set, row->rhs) )
            row->rhs = SCIPsetFloor(set, row->rhs);
      }
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
      int* cols_probindex;
      Real* vals;
      int s;
      int t;
      
      /* make sure, the row is sorted */
      SCIProwSort(row);
      assert(row->lpcolssorted);
      assert(row->nonlpcolssorted);
      
      /* merge equal columns, thereby recalculating whether the row's activity is always integral */
      cols = row->cols;
      cols_probindex = row->cols_probindex;
      vals = row->vals;
      assert(cols != NULL);
      assert(cols_probindex != NULL);
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
               row->integral = row->integral && SCIPvarIsIntegral(cols[t]->var) && SCIPsetIsIntegral(set, vals[t]);
               t++;
            }
            cols[t] = cols[s];
            cols_probindex[t] = cols_probindex[s];
            vals[t] = vals[s];
         }
      }
      if( !SCIPsetIsZero(set, vals[t]) )
      {
         row->integral = row->integral && SCIPvarIsIntegral(cols[t]->var) && SCIPsetIsIntegral(set, vals[t]);
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
static
void rowCalcLPActivity(
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
      rowCalcLPActivity(row, stat);
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
static
void rowCalcPseudoActivity(
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
   row->validpsactivitybdchg = stat->nboundchanges;
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
   assert(row->validpsactivitybdchg <= stat->nboundchanges);

   /* check, if activity bounds has to be calculated */
   if( row->validpsactivitybdchg != stat->nboundchanges )
      rowCalcPseudoActivity(row, stat);
   assert(row->validpsactivitybdchg == stat->nboundchanges);
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

   *solactivity = MAX(*solactivity, -set->infinity);
   *solactivity = MIN(*solactivity, +set->infinity);

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
   assert(!SCIPsetIsInfinity(set, ABS(row->constant)));
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
         mininfinite |= SCIPsetIsInfinity(set, -col->lb);
         maxinfinite |= SCIPsetIsInfinity(set, col->ub);
         if( !mininfinite )
            row->minactivity += val * col->lb;
         if( !maxinfinite )
            row->maxactivity += val * col->ub;
      }
      else
      {
         mininfinite |= SCIPsetIsInfinity(set, col->ub);
         maxinfinite |= SCIPsetIsInfinity(set, -col->lb);
         if( !mininfinite )
            row->minactivity += val * col->ub;
         if( !maxinfinite )
            row->maxactivity += val * col->lb;
      }
   }

   if( mininfinite )
      row->minactivity = -set->infinity;
   if( maxinfinite )
      row->maxactivity = set->infinity;
   row->validactivitybdsbdchg = stat->nboundchanges;

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
   assert(row->validactivitybdsbdchg <= stat->nboundchanges);

   /* check, if activity bounds has to be calculated */
   if( row->validactivitybdsbdchg != stat->nboundchanges )
      rowCalcActivityBounds(row, set, stat);
   assert(row->validactivitybdsbdchg == stat->nboundchanges);
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
   assert(row->validactivitybdsbdchg <= stat->nboundchanges);

   /* check, if activity bounds has to be calculated */
   if( row->validactivitybdsbdchg != stat->nboundchanges )
      rowCalcActivityBounds(row, set, stat);
   assert(row->validactivitybdsbdchg == stat->nboundchanges);
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
      fprintf(file, "%+g%s ", row->vals[i], SCIPvarGetName(row->cols[i]->var));
   }

   /* print constant */
   if( ABS(row->constant) > SCIP_DEFAULT_EPSILON )
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

/** get euclidean norm of row vector */
Real SCIProwGetNorm(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return sqrt(row->sqrnorm);
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

   return row->lppos;
}

/** returns TRUE iff row is member of current LP */
Bool SCIProwIsInLP(
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

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
   col->farkas = SCIP_INVALID;
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

      /*debugMessage("flushing added column <%s>:", SCIPvarGetName(col->var));*/
      /*debug( SCIPcolPrint(col, NULL) );*/

      /* Because the column becomes a member of the LP solver, it now can take values
       * different from zero. That means, we have to include the column in the corresponding
       * row vectors.
       */
      CHECK_OKAY( colLink(col, memhdr, set, lp) );

      lp->lpicols[c] = col;
      col->lpipos = c;
      col->primsol = SCIP_INVALID;
      col->redcost = SCIP_INVALID;
      col->farkas = SCIP_INVALID;
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
   LP*              lp                  /**< current LP data */
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
      }
      lp->nlpirows = lp->lpifirstchgrow;
      lp->flushdeletedrows = TRUE;
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

      debugMessage("flushing added row:");
      debug( SCIProwPrint(row, NULL) );

      /* Because the row becomes a member of the LP solver, its dual variable now can take values
       * different from zero. That means, we have to include the row in the corresponding
       * column vectors.
       */
      CHECK_OKAY( rowLink(row, memhdr, set, lp) );

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
   
   return SCIP_OKAY;
}

/** applies all cached column bound and objective changes to the LP */
static
RETCODE lpFlushChgCols(
   LP*              lp,                 /**< current LP data */
   MEMHDR*          memhdr,             /**< block memory */
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
   assert(memhdr != NULL);

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
         if( col->objchanged )
         {
            assert(nobjchg < lp->ncols);
            objind[nobjchg] = col->lpipos;
            obj[nobjchg] = col->obj;
            nobjchg++;
            col->objchanged = FALSE;
         }
         if( col->lbchanged || col->ubchanged )
         {
            assert(nbdchg < lp->ncols);
            bdind[nbdchg] = col->lpipos;
            if( SCIPsetIsInfinity(set, -col->lb) )
               lb[nbdchg] = -infinity;
            else
               lb[nbdchg] = col->lb;
            if( SCIPsetIsInfinity(set, col->ub) )
               ub[nbdchg] = infinity;
            else
               ub[nbdchg] = col->ub;
            nbdchg++;
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
   }

   /* change bounds in LP */
   if( nbdchg > 0 )
   {
      debugMessage("flushing bound changes: change %d bounds of %d changed columns\n", nbdchg, lp->nchgcols);
      CHECK_OKAY( SCIPlpiChgBounds(lp->lpi, nbdchg, bdind, lb, ub) );
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
   MEMHDR*          memhdr,             /**< block memory */
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
   assert(memhdr != NULL);

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
         if( row->lhschanged || row->rhschanged )
         {
            assert(nchg < lp->nrows);
            ind[nchg] = row->lpipos;
            if( SCIPsetIsInfinity(set, -row->lhs) )
               lhs[nchg] = -infinity;
            else
               lhs[nchg] = row->lhs - row->constant;
            if( SCIPsetIsInfinity(set, row->rhs) )
               rhs[nchg] = infinity;
            else
               rhs[nchg] = row->rhs - row->constant;
            nchg++;
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
   
   assert(!lp->solved);

   lp->flushdeletedcols = FALSE;
   lp->flushaddedcols = FALSE;
   lp->flushdeletedrows = FALSE;
   lp->flushaddedrows = FALSE;

   CHECK_OKAY( lpFlushDelCols(lp) );
   CHECK_OKAY( lpFlushDelRows(lp) );
   CHECK_OKAY( lpFlushChgCols(lp, memhdr, set) );
   CHECK_OKAY( lpFlushChgRows(lp, memhdr, set) );
   CHECK_OKAY( lpFlushAddCols(lp, memhdr, set) );
   CHECK_OKAY( lpFlushAddRows(lp, memhdr, set) );

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
         rowSwapCoeffs(row, pos, row->nlpcols-1);
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
         colSwapCoeffs(col, pos, col->nlprows-1);
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
         rowSwapCoeffs(row, pos, row->nlpcols);
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
         assert(0 <= pos && pos < col->nlprows);

         col->nlprows--;
         colSwapCoeffs(col, pos, col->nlprows);
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
   assert(lp != NULL);
   assert(set != NULL);
   assert(name != NULL);

   ALLOC_OKAY( allocMemory(lp) );

   /* open LP Solver interface */
   CHECK_OKAY( SCIPlpiCreate(&(*lp)->lpi, name) );

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
   (*lp)->cutoffbound = set->infinity;
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
   (*lp)->lpiuobjlim = set->infinity;
   (*lp)->lpifeastol = set->feastol;
   (*lp)->lpidualfeastol = set->dualfeastol;
   (*lp)->lpifromscratch = FALSE;
   (*lp)->lpifastmip = TRUE;
   (*lp)->lpiscaling = TRUE;
   (*lp)->lpiitlim = INT_MAX;
   (*lp)->lastwasprimal = FALSE;

   /* set objective sense */
   CHECK_OKAY( SCIPlpiChgObjsen((*lp)->lpi, SCIP_OBJSEN_MINIMIZE) );

   /* set default parameters in LP solver */
   CHECK_OKAY( lpSetRealpar(*lp, SCIP_LPPAR_UOBJLIM, (*lp)->lpiuobjlim) );
   CHECK_OKAY( lpSetRealpar(*lp, SCIP_LPPAR_FEASTOL, (*lp)->lpifeastol) );
   CHECK_OKAY( lpSetRealpar(*lp, SCIP_LPPAR_DUALFEASTOL, (*lp)->lpidualfeastol) );
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_FROMSCRATCH, (*lp)->lpifromscratch) );
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_FASTMIP, (*lp)->lpifastmip) );
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_SCALING, (*lp)->lpiscaling) );
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_LPITLIM, (*lp)->lpiitlim) );
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_PRICING, SCIP_PRICING_AUTO) ); /*lint !e641*/
   CHECK_OKAY( lpSetIntpar(*lp, SCIP_LPPAR_LPINFO, FALSE) );

   return SCIP_OKAY;
}

/** frees LP data object */
RETCODE SCIPlpFree(
   LP**             lp,                 /**< pointer to LP data object */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(*lp != NULL);
   
   CHECK_OKAY( SCIPlpClear(*lp, memhdr, set) );

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
   COL*             col                 /**< LP column */
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
   CHECK_OKAY( ensureColsSize(lp, set, lp->ncols+1) );
   lp->cols[lp->ncols] = col;
   col->lppos = lp->ncols;
   col->age = 0;
   lp->ncols++;
   if( col->removeable )
      lp->nremoveablecols++;
   lp->flushed = FALSE;
   lp->solved = FALSE;
   lp->dualfeasible = FALSE;
   lp->lpobjval = SCIP_INVALID;
   lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

   /* update column arrays of all linked rows */
   colUpdateAddLP(col);

   return SCIP_OKAY;
}

/** adds a row to the LP and captures it */
RETCODE SCIPlpAddRow(
   LP*              lp,                 /**< LP data */
   SET*             set,                /**< global SCIP settings */
   ROW*             row                 /**< LP row */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);
   assert(row != NULL);
   assert(row->len == 0 || row->cols != NULL);
   assert(row->lppos == -1);

   SCIProwCapture(row);

   debugMessage("adding row <%s> to LP (%d rows, %d cols)\n", row->name, lp->nrows, lp->ncols);
   CHECK_OKAY( ensureRowsSize(lp, set, lp->nrows+1) );
   lp->rows[lp->nrows] = row;
   row->lppos = lp->nrows;
   row->age = 0;
   lp->nrows++;
   if( row->removeable )
      lp->nremoveablerows++;
   lp->flushed = FALSE;
   lp->solved = FALSE;
   lp->primalfeasible = FALSE;
   lp->lpobjval = SCIP_INVALID;
   lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

   /* update row arrays of all linked columns */
   rowUpdateAddLP(row);

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
         lp->ncols--;

         /* count removeable columns */
         if( col->removeable )
            lp->nremoveablecols--;

         /* update column arrays of all linked rows */
         colUpdateDelLP(col);
      }
      assert(lp->ncols == newncols);
      lp->lpifirstchgcol = MIN(lp->lpifirstchgcol, newncols);
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
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
         lp->nrows--;

         /* count removeable rows */
         if( row->removeable )
            lp->nremoveablerows--;

         /* update row arrays of all linked columns */
         rowUpdateDelLP(row);

         CHECK_OKAY( SCIProwRelease(&lp->rows[r], memhdr, set, lp) );
      }
      assert(lp->nrows == newnrows);
      lp->lpifirstchgrow = MIN(lp->lpifirstchgrow, newnrows);
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
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
   int              nvars,              /**< number of active variables in the problem */
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
   assert(weights != NULL);
   assert(sumcoef != NULL);
   assert(sumlhs != NULL);
   assert(sumrhs != NULL);

   /**@todo test, if a column based summation is faster */

   CHECK_OKAY( SCIPrealarrayClear(sumcoef) );
   CHECK_OKAY( SCIPrealarrayExtend(sumcoef, set, 0, nvars-1) );
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
         assert(row->len == 0 || row->cols_probindex != NULL);
         assert(row->len == 0 || row->vals != NULL);

         /* add the row coefficients to the sum */
         for( i = 0; i < row->len; ++i )
         {
            assert(row->cols[i] != NULL);
            assert(row->cols[i]->var != NULL);
            assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
            assert(SCIPvarGetCol(row->cols[i]->var) == row->cols[i]);
            assert(SCIPvarGetProbindex(row->cols[i]->var) == row->cols[i]->var_probindex);
            assert(SCIPvarGetProbindex(row->cols[i]->var) == row->cols_probindex[i]);
            idx = row->cols_probindex[i];
            assert(0 <= idx && idx < nvars);
            CHECK_OKAY( SCIPrealarrayIncVal(sumcoef, set, idx, weights[r] * row->vals[i]) );
         }
         
         /* add the row sides to the sum, depending on the sign of the weight */
         if( weights[r] > 0.0 )
         {
            lhsinfinite |= SCIPsetIsInfinity(set, -row->lhs);
            if( !lhsinfinite )
               (*sumlhs) += weights[r] * (row->lhs - row->constant);
            rhsinfinite |= SCIPsetIsInfinity(set, row->rhs);
            if( !rhsinfinite )
               (*sumrhs) += weights[r] * (row->rhs - row->constant);
         }
         else
         {
            lhsinfinite |= SCIPsetIsInfinity(set, row->rhs);
            if( !lhsinfinite )
               (*sumlhs) += weights[r] * (row->rhs - row->constant);
            rhsinfinite |= SCIPsetIsInfinity(set, -row->lhs);
            if( !rhsinfinite )
               (*sumrhs) += weights[r] * (row->lhs - row->constant);
         }
      }
      else
         weights[r] = 0.0;
   }

   *sumlhs = -set->infinity;
   *sumrhs = set->infinity;

   return SCIP_OKAY;
}

/** builds a weighted sum of rows, and decides whether to use the left or right hand side of the rows in summation */
static
void sumMIRRow(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< LP data */
   int              nvars,              /**< number of active variables in the problem */
   Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   Real             scale,              /**< additional scaling factor multiplied to all rows */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*             slacksign,          /**< stores the sign of the row's slack variable in summation */
   Bool*            emptyrow            /**< pointer to store whether the returned row is empty */
   )
{
   ROW* row;
   Real rowactivity;
   int idx;
   int r;
   int i;

   assert(lp != NULL);
   assert(weights != NULL);
   assert(SCIPsetIsPositive(set, scale));
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(slacksign != NULL);
   assert(emptyrow != NULL);

   clearMemoryArray(mircoef, nvars);
   *mirrhs = 0.0;
   *emptyrow = TRUE;
   for( r = 0; r < lp->nrows; ++r )
   {
      row = lp->rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_probindex != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* modifiable rows cannot be part of a MIR row summation; close to zero weights are ignored */
      if( !row->modifiable && !SCIPsetIsZero(set, weights[r]) )
      {
         /* Decide, if we want to use the left or the right hand side of the row in the summation.
          * If the current row activity is closer to the left hand side, we use the  lhs <= a*x  part of the row,
          * and treat it implicitly as  a*x - s == lhs. Otherwise, we use the  a*x <= rhs  part of the row,
          * and treat it implicitly as  a*x + s == rhs. We have to remember, which sign the implicit slack variable
          * has.
          */
         *emptyrow = FALSE;
         rowactivity = SCIProwGetLPActivity(row, stat, lp);
         assert(SCIPsetIsFeasGE(set, rowactivity, row->lhs));
         assert(SCIPsetIsFeasLE(set, rowactivity, row->rhs));
            
         if( rowactivity < (row->lhs + row->rhs)/2.0 )
         {
            slacksign[r] = -1;
            (*mirrhs) += scale * weights[r] * (row->lhs - row->constant);
         }
         else
         {
            slacksign[r] = +1;
            (*mirrhs) += scale * weights[r] * (row->rhs - row->constant);
         }

         /* add the row coefficients to the sum */
         for( i = 0; i < row->len; ++i )
         {
            assert(row->cols[i] != NULL);
            assert(row->cols[i]->var != NULL);
            assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
            assert(SCIPvarGetCol(row->cols[i]->var) == row->cols[i]);
            assert(SCIPvarGetProbindex(row->cols[i]->var) == row->cols[i]->var_probindex);
            assert(SCIPvarGetProbindex(row->cols[i]->var) == row->cols_probindex[i]);
            idx = row->cols_probindex[i];
            assert(0 <= idx && idx < nvars);
            mircoef[idx] += scale * weights[r] * row->vals[i];
         }
      }
      else
      {
         slacksign[r] = 0;
         weights[r] = 0.0;
      }
   }
}

#define BOUNDSWITCH 0.9999
/** Transform equation  a*x == b, lb <= x <= ub  into standard form
 *    a'*x' == b, 0 <= x' <= ub'.
 *  Transform variables:
 *    x'_j := x_j - lb_j,       x_j == x'_j + lb_j,       a'_j =  a_j,      if x^_j is closer to lb
 *    x'_j := ub_j - x_j,       x_j == ub_j - x'_j,       a'_j = -a_j,      if x^_j is closer to ub
 *  and move the constant terms "a_j * lb_j" and "a_j * ub_j" to the rhs.
 */
static
void transformMIRRow(
   SET*             set,                /**< global SCIP settings */
   int              nvars,              /**< number of active variables in the problem */
   VAR**            vars,               /**< active variables in the problem */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*             varsign,            /**< stores the sign of the transformed variable in summation */
   Bool*            freevariable        /**< stores whether a free variable was found in MIR row -> invalid summation */
   )
{
   VAR* var;
   Real lb;
   Real ub;
   int idx;
   int v;

   assert(vars != NULL);
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(varsign != NULL);
   assert(freevariable != NULL);

   *freevariable = FALSE;

   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      idx = SCIPvarGetProbindex(var);
      assert(0 <= idx && idx < nvars);

      if( SCIPsetIsZero(set, mircoef[idx]) )
      {
         varsign[idx] = +1;
         continue;
      }

      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         if( SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_LOWER )
         {
            varsign[idx] = +1;
            (*mirrhs) -= mircoef[idx] * lb;
         }
         else
         {
            varsign[idx] = -1;
            (*mirrhs) -= mircoef[idx] * ub;
         }
      }
      else if( !SCIPsetIsInfinity(set, -lb) && !SCIPsetIsInfinity(set, ub) )
      {
         if( SCIPvarGetCol(var)->primsol <= (1.0-BOUNDSWITCH)*lb + BOUNDSWITCH*ub )
         {
            varsign[idx] = +1;
            (*mirrhs) -= mircoef[idx] * lb;
         }
         else
         {
            varsign[idx] = -1;
            (*mirrhs) -= mircoef[idx] * ub;
         }
      }
      else if( !SCIPsetIsInfinity(set, -lb) )
      {
         varsign[idx] = +1;
         (*mirrhs) -= mircoef[idx] * lb;
      }
      else if( !SCIPsetIsInfinity(set, ub) )
      {
         varsign[idx] = -1;
         (*mirrhs) -= mircoef[idx] * ub;
      }
      else
      {
         /* we found a free variable in the row with non-zero coefficient
          *  -> the MIR row cannot be transformed in standard form
          */
         *freevariable = TRUE;
         return;
      }            
   }
}

/** Calculate fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j) , and derive MIR cut
 *    a~*x' <= down(b)
 *  integers :  a~_j = down(a'_j)                      , if f_j <= f0
 *              a~_j = down(a'_j) + (f_j - f0)/(1 - f0), if f_j >  f0
 *  continuous: a~_j = 0                               , if a'_j >= 0
 *              a~_j = a'_j/(1 - f0)                   , if a'_j <  0
 *
 *  Transform inequality back to a°*x <= down(b):
 *    x'_j := x_j - lb_j,       x_j == x'_j + lb_j,       a'_j =  a_j,      if x^_j is closer to lb
 *    x'_j := ub_j - x_j,       x_j == ub_j - x'_j,       a'_j = -a_j,      if x^_j is closer to ub
 *    a°_j :=  a~_j, if x^_j is closer to lb
 *    a°_j := -a~_j, if x^_j is closer to ub
 *  and move the constant terms
 *    -a~_j * lb_j == -a°_j * lb_j, or
 *     a~_j * ub_j == -a°_j * ub_j
 *  to the rhs.
 */
static
void roundMIRRow(
   SET*             set,                /**< global SCIP settings */
   int              nvars,              /**< number of active variables in the problem */
   VAR**            vars,               /**< active variables in the problem */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*             varsign,            /**< stores the sign of the transformed variable in summation */
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
   int idx;
   int v;

   assert(vars != NULL);
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(varsign != NULL);
   assert(0.0 < f0 && f0 < 1.0);
   assert(cutactivity != NULL);

   *cutactivity = 0.0;
   onedivoneminusf0 = 1.0 / (1.0 - f0);

   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      idx = SCIPvarGetProbindex(var);
      assert(0 <= idx && idx < nvars);

      /* calculate the coefficient in the retransformed cut */
      aj = varsign[idx] * mircoef[idx];
      if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      {
         /* integer variable */
         downaj = SCIPsetFloor(set, aj);
         fj = aj - downaj;
         if( SCIPsetIsSumLE(set, fj, f0) )
            cutaj = varsign[idx] * downaj;
         else
            cutaj = varsign[idx] * (downaj + (fj - f0) * onedivoneminusf0);
      }
      else
      {
         /* continuous variable */
         if( SCIPsetIsSumGE(set, aj, 0.0) )
            cutaj = 0.0;
         else
            cutaj = varsign[idx] * aj * onedivoneminusf0;
      }

      if( SCIPsetIsZero(set, cutaj) )
         mircoef[idx] = 0.0;
      else
      {
         mircoef[idx] = cutaj;
         (*cutactivity) += cutaj * SCIPvarGetLPSol(var);

         /* move the constant term  -a~_j * lb_j == -a°_j * lb_j , or  a~_j * ub_j == -a°_j * ub_j  to the rhs */
         if( varsign[idx] == +1 )
         {
            assert(!SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)));
            (*mirrhs) += cutaj * SCIPvarGetLbLocal(var);
         }
         else
         {
            assert(!SCIPsetIsInfinity(set, SCIPvarGetUbLocal(var)));
            (*mirrhs) += cutaj * SCIPvarGetUbLocal(var);
         }
      }
   }
}

/** substitute aggregated slack variables:
 *
 *  The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
 *  variable only appears in its own row:
 *     a'_r = weight[r] * slacksign[r].
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
      /* unused rows can be ignored */
      if( slacksign[r] == 0 )
         continue;

      assert(!SCIPsetIsZero(set, weights[r]));

      row = lp->rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_probindex != NULL);
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
         if( !SCIPsetIsNegative(set, ar) )
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
         assert(SCIPvarGetProbindex(row->cols[i]->var) == row->cols_probindex[i]);
         idx = row->cols_probindex[i];
         mircoef[idx] += mul * row->vals[i];
      }

      /* update the activity: we have to add  mul * a*x^  to the cut's activity (row activity = a*x^ + c) */
      (*cutactivity) += mul * (SCIProwGetLPActivity(row, stat, lp) - row->constant);

      /* move slack's constant to the right hand side */
      if( slacksign[r] == +1 )
      {
         /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a°_r * (rhs - c) to the right hand side */
         (*mirrhs) -= cutar * (row->rhs - row->constant);
      }
      else
      {
         /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a°_r * (c - lhs) to the right hand side */
         (*mirrhs) -= cutar * (row->constant - row->lhs);
      }
   }

   /* set rhs to zero, if it's very close to */
   if( SCIPsetIsZero(set, *mirrhs) )
      *mirrhs = 0.0;
}

/* calculates a MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a MIR cut.
 */
RETCODE SCIPlpCalcMIR(
   LP*              lp,                 /**< LP data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   int              nvars,              /**< number of active variables in the problem */
   VAR**            vars,               /**< active variables in the problem */
   Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   Real             scale,              /**< additional scaling factor multiplied to all rows */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size nvars */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   Bool*            success             /**< pointer to store whether the returned coefficients are a valid MIR cut */
   )
{
   int* slacksign;
   int* varsign;
   Real rhs;
   Real downrhs;
   Real f0;
   Bool emptyrow;
   Bool freevariable;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(vars != NULL);
   assert(weights != NULL);
   assert(SCIPsetIsPositive(set, scale));
   assert(mircoef != NULL);
   assert(mirrhs != NULL);
   assert(cutactivity != NULL);
   assert(success != NULL);

   /**@todo test, if a column based summation is faster */

   *success = FALSE;

   /* allocate temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &slacksign, lp->nrows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &varsign, nvars) );

   /* calculate the row summation */
   sumMIRRow(set, stat, lp, nvars, weights, scale, mircoef, &rhs, slacksign, &emptyrow);
   if( emptyrow )
      goto TERMINATE;

   /* Transform equation  a*x == b, lb <= x <= ub  into standard form
    *   a*x' == b, 0 <= x' <= ub'.
    * Transform variables:
    *   x'_j := x_j - lb_j,       x_j == x'_j + lb_j,       if x^_j is closer to lb
    *   x'_j := ub_j - x_j,       x_j == ub_j - x'_j,       if x^_j is closer to ub
    * and move the constant terms "a_j * lb_j" and "a_j * ub_j" to the rhs.
    */
   transformMIRRow(set, nvars, vars, mircoef, &rhs, varsign, &freevariable);
   if( freevariable )
      goto TERMINATE;

   /* Calculate fractionalities  f0 := b - down(b), f_j := a_j - down(a_j) , and derive MIR cut
    *   a~*x' <= down(b)
    * integers:   a~_j = down(a_j)                      , if f_j <= f0
    *             a~_j = down(a_j) + (f_j - f0)/(1 - f0), if f_j >  f0
    * continuous: a~_j = 0                              , if a_j >= 0
    *             a~_j = a_j/(1 - f0)                   , if a_j <  0
    * Keep in mind, that the varsign has to be implicitly incorporated into a~_j.
    * Transform inequality back to a°*x <= down(b):
    *   x'_j := x_j - lb_j,       x_j == x'_j + lb_j,       if x^_j is closer to lb
    *   x'_j := ub_j - x_j,       x_j == ub_j - x'_j,       if x^_j is closer to ub
    *   a°_j :=  a~_j, if x^_j is closer to lb
    *   a°_j := -a~_j, if x^_j is closer to ub
    * and move the constant terms
    *   -a~_j * lb_j == -a°_j * lb_j, or
    *    a~_j * ub_j == -a°_j * ub_j
    * to the rhs.
    */
   downrhs = SCIPsetFloor(set, rhs);
   f0 = rhs - downrhs;
   if( f0 < minfrac )
      goto TERMINATE;

   *mirrhs = downrhs;
   roundMIRRow(set, nvars, vars, mircoef, mirrhs, varsign, f0, cutactivity);

   /* substitute negatively aggregated slack variables:
    * - if row was aggregated with a positive factor (weight * slacksign), the a_j for the continuous
    *   slack variable is a_j > 0, which leads to a°_j = 0, so we can ignore the slack variable in
    *   the resulting cut
    * - if row was aggregated with a negative factor (weight * slacksign), the a_j for the continuous
    *   slack variable is a_j < 0, which leads to a°_j = a_j/(1 - f0), so we have to subtract 
    *   a°_j times the row to the cut to eliminate the slack variable
    */
   substituteMIRRow(set, stat, lp, weights, scale, mircoef, mirrhs, slacksign, f0, cutactivity);

   *success = TRUE;

 TERMINATE:
   /* free temporary memory */
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

   CHECK_OKAY( SCIPlpFlush(lp, memhdr, set) );

   CHECK_OKAY( SCIPlpiSetState(lp->lpi, memhdr, lpistate) );
   lp->primalfeasible = TRUE;
   lp->dualfeasible = TRUE;

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
   if( lp->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT && cutoffbound > lp->cutoffbound )
      lp->solved = FALSE;
   else if( lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && SCIPlpGetObjval(lp, set) >= cutoffbound )
      lp->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;

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
      printf("wrote LP to file <%s> (primal simplex, uobjlim=%g, feastol=%g/%g, fromscratch=%d, fastmip=%d, scaling=%d)\n", 
         fname, lp->lpiuobjlim, lp->lpifeastol, lp->lpidualfeastol, lp->lpifromscratch, lp->lpifastmip, lp->lpiscaling);
   }
#endif

   /* start timing */
   SCIPclockStart(stat->primallptime, set);

   /* call primal simplex */
   CHECK_OKAY( SCIPlpiSolvePrimal(lp->lpi) );
   lp->lastwasprimal = TRUE;

   /* stop timing */
   SCIPclockStop(stat->primallptime, set);

   /* count number of iterations */
   stat->lpcount++;
   CHECK_OKAY( SCIPlpGetIterations(lp, &iterations) );
   if( iterations > 0 ) /* don't count the resolves after removing unused columns/rows */
   {
      stat->nlps++;
      stat->nprimallps++;
      stat->nlpiterations += iterations;
      stat->nprimallpiterations += iterations;
      if( lp->diving )
      {
         stat->ndivinglps++;
         stat->ndivinglpiterations += iterations;
      }
   }

   debugMessage("solved primal LP in %d iterations\n", iterations);

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
      printf("wrote LP to file <%s> (dual simplex, uobjlim=%g, feastol=%g/%g, fromscratch=%d, fastmip=%d, scaling=%d)\n", 
         fname, lp->lpiuobjlim, lp->lpifeastol, lp->lpidualfeastol, lp->lpifromscratch, lp->lpifastmip, lp->lpiscaling);
   }
#endif

   /* start timing */
   SCIPclockStart(stat->duallptime, set);

   /* call dual simplex */
   CHECK_OKAY( SCIPlpiSolveDual(lp->lpi) );
   lp->lastwasprimal = FALSE;

   /* stop timing */
   SCIPclockStop(stat->duallptime, set);

   /* count number of iterations */
   stat->lpcount++;
   CHECK_OKAY( SCIPlpGetIterations(lp, &iterations) );
   if( iterations > 0 ) /* don't count the resolves after removing unused columns/rows */
   {
      stat->nlps++;
      stat->nduallps++;
      stat->nlpiterations += iterations;
      stat->nduallpiterations += iterations;
      if( lp->diving )
      {
         stat->ndivinglps++;
         stat->ndivinglpiterations += iterations;
      }
   }

   debugMessage("solved dual LP in %d iterations\n", iterations);

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
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->looseobjvalinf == 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);

   *lperror = FALSE;

   /* solve with given settings (usually fast but unprecise) */
   CHECK_OKAY( lpSetUobjlim(lp, set, lp->cutoffbound - lp->looseobjval) );
   CHECK_OKAY( lpSetFeastol(lp, set->feastol) );
   CHECK_OKAY( lpSetDualFeastol(lp, set->dualfeastol) );
   CHECK_OKAY( lpSetFromscratch(lp, fromscratch) );
   CHECK_OKAY( lpSetFastmip(lp, fastmip) );
   CHECK_OKAY( lpSetScaling(lp, set->scaling) );
   CHECK_OKAY( lpSimplex(lp, set, stat, useprimal) );

   /* check for stability */
   if( SCIPlpiIsStable(lp->lpi) )
      return SCIP_OKAY;

   /* if FASTMIP is turned on, solve again without FASTMIP */
   if( fastmip )
   {
      infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles in LP %d -- solve again without FASTMIP with %s simplex\n", 
         stat->nnodes, stat->nlps, useprimal ? "primal" : "dual");
      CHECK_OKAY( lpSetFastmip(lp, FALSE) );
      CHECK_OKAY( lpSimplex(lp, set, stat, useprimal) );

      /* check for stability */
      if( SCIPlpiIsStable(lp->lpi) )
         return SCIP_OKAY;
   }

   /* if not already done, solve again from scratch */
   if( !fromscratch )
   {
      infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
         "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s simplex\n", 
         stat->nnodes, stat->nlps, useprimal ? "primal" : "dual");
      CHECK_OKAY( lpSetFromscratch(lp, TRUE) );
      CHECK_OKAY( lpSimplex(lp, set, stat, useprimal) );
      
      /* check for stability */
      if( SCIPlpiIsStable(lp->lpi) )
         return SCIP_OKAY;
   }

   /* solve again with a tighter feasibility tolerance */
   infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
      "(node %lld) numerical troubles in LP %d -- solve again with tighter feasibility tolerance with %s simplex\n", 
      stat->nnodes, stat->nlps, useprimal ? "primal" : "dual");
   CHECK_OKAY( lpSetFeastol(lp, 0.001*set->feastol) );
   CHECK_OKAY( lpSimplex(lp, set, stat, useprimal) );

   /* check for stability */
   if( SCIPlpiIsStable(lp->lpi) )
      return SCIP_OKAY;

   /* solve again, use other simplex this time */
   infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
      "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s simplex\n", 
      stat->nnodes, stat->nlps, !useprimal ? "primal" : "dual");
   CHECK_OKAY( lpSetFeastol(lp, set->feastol) );
   CHECK_OKAY( lpSimplex(lp, set, stat, !useprimal) );

   /* check for stability */
   if( SCIPlpiIsStable(lp->lpi) )
      return SCIP_OKAY;

   /* solve again with tighter feasibility tolerance, use other simplex this time */
   infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
      "(node %lld) numerical troubles in LP %d -- solve again with tighter feasibility tolerance with %s simplex\n", 
      stat->nnodes, stat->nlps, !useprimal ? "primal" : "dual");
   CHECK_OKAY( lpSetFeastol(lp, 0.001*set->feastol) );
   CHECK_OKAY( lpSimplex(lp, set, stat, !useprimal) );

   /* check for stability */
   if( SCIPlpiIsStable(lp->lpi) )
      return SCIP_OKAY;

   /* solve again with opposite scaling setting */
   infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
      "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s simplex %s scaling\n", 
      stat->nnodes, stat->nlps, useprimal ? "primal" : "dual", !set->scaling ? "with" : "without");
   CHECK_OKAY( lpSetScaling(lp, !set->scaling) );
   CHECK_OKAY( lpSimplex(lp, set, stat, useprimal) );
   
   /* check for stability */
   if( SCIPlpiIsStable(lp->lpi) )
      return SCIP_OKAY;

   /* solve again with opposite scaling, use other simplex this time */
   infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
      "(node %lld) numerical troubles in LP %d -- solve again from scratch with %s simplex %s scaling\n", 
      stat->nnodes, stat->nlps, !useprimal ? "primal" : "dual", !set->scaling ? "with" : "without");
   CHECK_OKAY( lpSimplex(lp, set, stat, !useprimal) );

   /* check for stability */
   if( !SCIPlpiIsStable(lp->lpi) )
   {
      /* nothing worked -- store the instable LP to a file and exit with an LPERROR */
      infoMessage(set->verblevel, SCIP_VERBLEVEL_HIGH, "(node %lld) unresolved numerical troubles in LP %d\n", 
         stat->nnodes, stat->nlps);
      *lperror = TRUE;
   }

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
   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);

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
         lp->lpobjval = set->infinity;
      }
   }
   else if( SCIPlpiIsObjlimExc(lp->lpi) )
   {
#if 0 /* SOPLEX may return with objective limit reached in any case, because it doesn't distinct btw. primal and dual */
      if( lp->lastwasprimal )
      {
         errorMessage("Objective limit exceeded in primal simplex - this should not happen, because no lower limit exists\n");
         lp->lpsolstat = SCIP_LPSOLSTAT_ERROR;
         lp->lpobjval = -set->infinity;
         return SCIP_LPERROR;
      }
#endif
      lp->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
      lp->lpobjval = set->infinity;
   }
   else if( SCIPlpiIsPrimalInfeasible(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;
      lp->lpobjval = set->infinity;
   }
   else if( SCIPlpiIsPrimalUnbounded(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_UNBOUNDED;
      lp->lpobjval = -set->infinity;
   }
   else if( SCIPlpiIsIterlimExc(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_ITERLIMIT;
      lp->lpobjval = -set->infinity;
   }
   else if( SCIPlpiIsTimelimExc(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
      lp->lpobjval = -set->infinity;
   }
   else
   {
      errorMessage("Unknown return status of %s simplex (internal status: %d)\n", 
         lp->lastwasprimal ? "primal" : "dual", SCIPlpiGetInternalStatus(lp->lpi));
      lp->lpsolstat = SCIP_LPSOLSTAT_ERROR;
      return SCIP_LPERROR;
   }

   lp->solved = TRUE;

   debugMessage("solving %s LP returned solstat=%d\n", lp->lastwasprimal ? "primal" : "dual", lp->lpsolstat);

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

   CHECK_OKAY( lpSetIterationLimit(lp, itlim) );

   if( !lp->solved )
   {
      Bool primalfeasible;
      Bool dualfeasible;
      Bool fastmip;
      Bool fromscratch;

      /* flush changes to the LP solver */
      CHECK_OKAY( SCIPlpFlush(lp, memhdr, set) );

      /* set initial LP solver settings */
      fastmip = set->fastmip && !lp->flushaddedcols && !lp->flushdeletedcols;
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
         if( set->checklpfeas )
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
               infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %lld) solution of LP %d not optimal (pfeas=%d, dfeas=%d) -- solving again without FASTMIP\n",
                  stat->nnodes, stat->nlps, primalfeasible, dualfeasible);
               fastmip = FALSE;
               goto SOLVEAGAIN;
            }
            else if( !fromscratch )
            {
               /* solution is infeasible (this can happen due to numerical problems): solve again from scratch */
               infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
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
         if( !SCIPprobAllColsInLP(prob, set, lp) || set->exactsolve )
         {
            CHECK_OKAY( SCIPlpGetDualfarkas(lp, memhdr, set, stat) );
         }
         debugMessage(" -> LP infeasible\n");
         break;

      case SCIP_LPSOLSTAT_UNBOUNDED:
         CHECK_OKAY( SCIPlpGetUnboundedSol(lp, memhdr, set, stat) );
         debugMessage(" -> LP unbounded\n");
         break;

      case SCIP_LPSOLSTAT_OBJLIMIT:
         if( !SCIPprobAllColsInLP(prob, set, lp) || set->exactsolve )
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
         errorMessage("Error in LP solver\n");
         return SCIP_LPERROR;

      default:
         errorMessage("Unknown LP solution status\n");
         return SCIP_ERROR;
      }
   }

   return SCIP_OKAY;
}

/** gets solution status of last solve call */
LPSOLSTAT SCIPlpGetSolstat(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->solved || lp->lpsolstat == SCIP_LPSOLSTAT_NOTSOLVED);

   return lp->lpsolstat;
}

/** gets objective value of last solution */
Real SCIPlpGetObjval(
   LP*              lp,                 /**< current LP data */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert((lp->nloosevars > 0) || (lp->looseobjvalinf == 0 && lp->looseobjval == 0.0));
   assert(set != NULL);

   if( SCIPsetIsInfinity(set, lp->lpobjval) )
      return lp->lpobjval;
   else if( lp->looseobjvalinf > 0 )
      return -set->infinity;
   else
      return lp->lpobjval + lp->looseobjval;
}

/** gets part of objective value of last solution that results from COLUMN variables only */
Real SCIPlpGetColumnObjval(
   LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->solved);

   return lp->lpobjval;
}

/** gets part of objective value of last solution that results from LOOSE variables only */
Real SCIPlpGetLooseObjval(
   LP*              lp,                 /**< current LP data */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert((lp->nloosevars > 0) || (lp->looseobjvalinf == 0 && lp->looseobjval == 0.0));
   assert(set != NULL);

   if( lp->looseobjvalinf > 0 )
      return -set->infinity;
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
      return -set->infinity;
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
      if( SCIPsetIsInfinity(set, ABS(oldbound)) )
         pseudoobjvalinf--;
      else
         pseudoobjval -= oldbound * SCIPvarGetObj(var);
      assert(pseudoobjvalinf >= 0);
      if( SCIPsetIsInfinity(set, ABS(newbound)) )
         pseudoobjvalinf++;
      else
         pseudoobjval += newbound * SCIPvarGetObj(var);
   }
   assert(pseudoobjvalinf >= 0);

   if( pseudoobjvalinf > 0 )
      return -set->infinity;
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

      if( SCIPsetIsInfinity(set, ABS(oldbound)) )
         pseudoobjvalinf--;
      else
      {
         SCIPintervalSet(&bd, oldbound);
         SCIPintervalMul(&prod, bd, obj);
         SCIPintervalSub(&psval, psval, prod);
      }
      assert(pseudoobjvalinf >= 0);
      if( SCIPsetIsInfinity(set, ABS(newbound)) )
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
      return -set->infinity;
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
   assert(!SCIPsetIsInfinity(set, ABS(oldobj)));
   assert(!SCIPsetIsInfinity(set, oldlb));
   assert(!SCIPsetIsInfinity(set, -oldub));
   assert(!SCIPsetIsInfinity(set, ABS(newobj)));
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
   assert(!SCIPsetIsInfinity(set, ABS(oldobj)));
   assert(!SCIPsetIsInfinity(set, oldlb));
   assert(!SCIPsetIsInfinity(set, -oldub));
   assert(!SCIPsetIsInfinity(set, ABS(newobj)));
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

   if( set->exactsolve )
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

   if( set->exactsolve )
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

   if( set->exactsolve )
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

   if( set->exactsolve )
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

   if( set->exactsolve )
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

   debugMessage("getting new LP solution\n");

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
      /*debugMessage(" col <%s> [%g,%g]: primsol=%.9f, redcost=%.9f, pfeas=%d/%d, dfeas=%d\n",
        SCIPvarGetName(lpicols[c]->var), lpicols[c]->lb, lpicols[c]->ub, lpicols[c]->primsol, lpicols[c]->redcost,
        SCIPsetIsFeasGE(set, lpicols[c]->primsol, lpicols[c]->lb),
        SCIPsetIsFeasLE(set, lpicols[c]->primsol, lpicols[c]->ub),
        !SCIPsetIsFeasNegative(set, lpicols[c]->redcost));*/
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
      /*debugMessage(" row <%s> [%g,%g]: dualsol=%.9f, activity=%.9f, pfeas=%d/%d, dfeas=%d/%d\n", 
        lpirows[r]->name, lpirows[r]->lhs, lpirows[r]->rhs, lpirows[r]->dualsol, lpirows[r]->activity,
        SCIPsetIsFeasGE(set, lpirows[r]->activity, lpirows[r]->lhs),
        SCIPsetIsFeasLE(set, lpirows[r]->activity, lpirows[r]->rhs),
        SCIPsetIsInfinity(set, -lpirows[r]->lhs) ? !SCIPsetIsFeasPositive(set, lpirows[r]->dualsol) : TRUE,
        SCIPsetIsInfinity(set, lpirows[r]->rhs) ? !SCIPsetIsFeasNegative(set, lpirows[r]->dualsol) : TRUE);
      */
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
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_UNBOUNDED);
   assert(SCIPsetIsInfinity(set, -lp->lpobjval));
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp->validsollp <= stat->lpcount);

   /* check if the values are already calculated */
   if( lp->validsollp == stat->lpcount )
      return SCIP_OKAY;
   lp->validsollp = stat->lpcount;

   debugMessage("getting new unbounded LP solution\n");

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
   rayscale = -2*set->infinity/rayobjval;

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
   assert(iterations != NULL);

   CHECK_OKAY( SCIPlpiGetIntpar(lp->lpi, SCIP_LPPAR_LPITER, iterations) );

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
      if( SCIPsetIsGT(set, lpirows[r]->activity, lpirows[r]->lhs)
         && SCIPsetIsLT(set, lpirows[r]->activity, lpirows[r]->rhs) )
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
      assert(lp->flushed == TRUE);
      lp->lpifirstchgcol = lp->nlpicols;
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

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

         CHECK_OKAY( SCIProwRelease(&lp->rows[r], memhdr, set, lp) );
         assert(lp->rows[r] == NULL);
         lp->lpirows[r] = NULL;
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
      assert(lp->flushed == TRUE);
      lp->lpifirstchgrow = lp->nlpirows;
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

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

   if( lp->nremoveablecols == 0 )
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
         && cols[c]->age > set->colagelimit
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

   if( lp->nremoveablerows == 0 )
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
         && rows[r]->age > set->rowagelimit )
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
      if( lpirows[r]->removeable
         && SCIPsetIsGT(set, lpirows[r]->activity, lpirows[r]->lhs)
         && SCIPsetIsLT(set, lpirows[r]->activity, lpirows[r]->rhs) )
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
      lp->firstnewcol, lp->ncols, set->cleanupcols, lp->firstnewrow, lp->nrows, set->cleanuprows);

   if( set->cleanupcols && lp->firstnewcol < lp->ncols )
   {
      CHECK_OKAY( lpCleanupCols(lp, set, stat, lp->firstnewcol) );
   }
   if( set->cleanuprows && lp->firstnewrow < lp->nrows )
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

   if( /*set->cleanupcols &&*/ 0 < lp->ncols )
   {
      CHECK_OKAY( lpCleanupCols(lp, set, stat, 0) );
   }
   if( /*set->cleanuprows &&*/ 0 < lp->nrows )
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
   assert(!lp->diving);
   assert(lp->divelpistate == NULL);

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
   CHECK_OKAY( SCIPlpiSetState(lp->lpi, memhdr, lp->divelpistate) );
   CHECK_OKAY( SCIPlpiFreeState(lp->lpi, memhdr, &lp->divelpistate) );
   assert(lp->divelpistate == NULL);

   /* resolve LP to reset solution */
   CHECK_OKAY( SCIPlpSolveAndEval(lp, memhdr, set, stat, prob, -1, FALSE, &lperror) );
   if( lperror )
   {
      infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL,
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
   assert(lp->flushed);
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
         SCIPintervalSet(&b, row->lhs);
      }
      else if( SCIPsetIsFeasNegative(set, y) )
      {
         SCIPintervalSet(&yinter[j], y);
         SCIPintervalSet(&b, row->rhs);
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

#endif
