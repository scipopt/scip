/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lp.c
 * @brief  LP management and variable's domains datastructures and methods
 * @author Tobias Achterberg
 */

/** The main datastructures for storing an LP are the rows and the columns.
 *  A row can live on its own (if it was created by a separator), or as LP
 *  relaxation of a constraint. Thus, it has a numuses-counter, and is
 *  deleted, if not needed any more.
 *  A column cannot live on its own. It is always connected to a problem
 *  variable. Because pricing is always problem specific, it cannot create
 *  LP columns without introducing new variables. Thus, each column is
 *  connected to exactly one variable, and is deleted, if the variable
 *  is deleted.
 *
 *  In LP management, we have to differ between the actual LP and the LP
 *  stored in the LP solver. All LP methods affect the actual LP only. 
 *  Before solving the actual LP with the LP solver or setting an LP state,
 *  the LP solvers data has to be updated to the actual LP with a call to
 *  lpFlush().
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "sort.h"
#include "lp.h"


/** list of columns */
struct ColList
{
   COL*             col;                /**< pointer to this column */
   COLLIST*         next;               /**< pointer to next collist entry */
};

/** list of rows */
struct RowList
{
   ROW*             row;                /**< pointer to this row */
   ROWLIST*         next;               /**< pointer to next rowlist entry */
};


/*
 * memory growing methods for dynamically allocated arrays
 */

static
RETCODE ensureChgbdsSize(               /**< ensures, that chgbds array can store at least num entries */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nchgbds <= lp->chgbdssize);
   
   if( num > lp->chgbdssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(lp->chgbds, newsize) );
      lp->chgbdssize = newsize;
   }
   assert(num <= lp->chgbdssize);

   return SCIP_OKAY;
}

static
RETCODE ensureLpicolsSize(              /**< ensures, that lpicols array can store at least num entries */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nlpicols <= lp->lpicolssize);
   
   if( num > lp->lpicolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(lp->lpicols, newsize) );
      lp->lpicolssize = newsize;
   }
   assert(num <= lp->lpicolssize);

   return SCIP_OKAY;
}

static
RETCODE ensureLpirowsSize(              /**< ensures, that lpirows array can store at least num entries */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nlpirows <= lp->lpirowssize);
   
   if( num > lp->lpirowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(lp->lpirows, newsize) );
      lp->lpirowssize = newsize;
   }
   assert(num <= lp->lpirowssize);

   return SCIP_OKAY;
}

static
RETCODE ensureColsSize(                 /**< ensures, that cols array can store at least num entries */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->ncols <= lp->colssize);
   
   if( num > lp->colssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(lp->cols, newsize) );
      lp->colssize = newsize;
   }
   assert(num <= lp->colssize);

   return SCIP_OKAY;
}

static
RETCODE ensureRowsSize(                 /**< ensures, that rows array can store at least num entries */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nrows <= lp->rowssize);
   
   if( num > lp->rowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(lp->rows, newsize) );
      lp->rowssize = newsize;
   }
   assert(num <= lp->rowssize);

   return SCIP_OKAY;
}



/*
 * double linked coefficient matrix methods 
 */

static
RETCODE colLink(                        /**< insert column coefficients in corresponding rows */
   COL*             col,                /**< column data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   int i;

   assert(col != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( !col->linked )
   {
      for( i = 0; i < col->len; ++i )
      {
         assert(!SCIPsetIsZero(set, col->val[i]));
         
         CHECK_OKAY( SCIProwAddCoeff(col->row[i], memhdr, set, lp, col, col->val[i]) );
      }
      col->linked = TRUE;
   }

   return SCIP_OKAY;
}

static
void colUnlink(                         /**< removes column coefficients from corresponding rows */
   COL*             col,                /**< column data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   int i;

   assert(col != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( col->linked )
   {
      for( i = 0; i < col->len; ++i )
         SCIProwDelCoeff(col->row[i], set, lp, col);
      col->linked = FALSE;
   }
}

static
RETCODE rowLink(                        /**< insert row coefficients in corresponding columns */
   ROW*             row,                /**< row data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   int i;

   assert(row != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( !row->linked )
   {
      for( i = 0; i < row->len; ++i )
      {
         assert(!SCIPsetIsZero(set, row->val[i]));
         
         CHECK_OKAY( SCIPcolAddCoeff(row->col[i], memhdr, set, lp, row, row->val[i]) );
      }
      row->linked = TRUE;
   }

   return SCIP_OKAY;
}

static
void rowUnlink(                         /**< removes row coefficients from corresponding columns */
   ROW*             row,                /**< row data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   int i;

   assert(row != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( row->linked )
   {
      for( i = 0; i < row->len; ++i )
         SCIPcolDelCoeff(row->col[i], set, lp, row);
      row->linked = FALSE;
   }
}



/*
 * LP solver data update
 */

static
void coefChanged(                       /**< announces, that the given coefficient in the constraint matrix changed */
   ROW*             row,                /**< LP row */
   COL*             col,                /**< LP col */
   LP*              lp                  /**< actual LP data */
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
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->objval = SCIP_INVALID;
   }
}
   
static
RETCODE lpFlushDelCols(                 /**< applies all cached column removals to the LP solver */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   assert(lp->lpifirstchgcol <= lp->nlpicols);
      
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

      CHECK_OKAY( SCIPlpiDelCols(lp->lpi, lp->lpifirstchgcol, lp->nlpicols-1) );
      for( i = lp->lpifirstchgcol; i < lp->nlpicols; ++i )
      {
         assert(!lp->lpicols[i]->inLP);
         lp->lpicols[i]->lpipos = -1;
         lp->lpicols[i]->primsol = 0.0;
         lp->lpicols[i]->redcost = SCIP_INVALID;
      }
      lp->nlpicols = lp->lpifirstchgcol;
   }
   assert(lp->nlpicols == lp->lpifirstchgcol);

   return SCIP_OKAY;
}

static
RETCODE lpFlushAddCols(                 /**< applies all cached column additions to the LP solver */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr              /**< block memory */
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
   const VAR* var;
   int c;
   int nnonz;
   int naddcols;
   int naddcoefs;
   int i;
   int lpipos;
   int tmpsizecols;
   int tmpsizecoefs;

   assert(lp != NULL);
   assert(lp->lpifirstchgcol == lp->nlpicols);
   assert(memhdr != NULL);
   assert(set != NULL);

   /* if there are no columns to add, we are ready */
   if( lp->ncols == lp->nlpicols )
      return SCIP_OKAY;

   /* add the additional columns */
   assert(lp->ncols > lp->nlpicols);
   CHECK_OKAY( ensureLpicolsSize(lp, set, lp->ncols) );

   /* count the (maximal) number of added coefficients, calculate the number of added columns */
   naddcols = lp->ncols - lp->nlpicols;
   naddcoefs = 0;
   for( c = lp->nlpicols; c < lp->ncols; ++c )
      naddcoefs += lp->cols[c]->len;
   assert(naddcols > 0);

   /* get temporary memory for changes */
   tmpsizecols = SCIPsetCalcMemGrowSize(set, naddcols);     /* use standard sizes to reuse memory more often */
   tmpsizecoefs = SCIPsetCalcMemGrowSize(set, naddcoefs);
   assert(tmpsizecols >= naddcols);
   assert(tmpsizecoefs >= naddcoefs);
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, obj, tmpsizecols) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, lb, tmpsizecols) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, ub, tmpsizecols) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, beg, tmpsizecols) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, ind, tmpsizecoefs) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, val, tmpsizecoefs) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, name, tmpsizecols) );
   
   /* fill temporary memory with column data */
   nnonz = 0;
   for( c = lp->nlpicols; c < lp->ncols; ++c )
   {
      col = lp->cols[c];
      assert(col != NULL);
      assert(col->inLP);
      assert(nnonz + col->len <= naddcoefs);

      /* Because the column becomes a member of the LP solver, it now can take values
       * different from zero. That means, we have to include the column in the corresponding
       * row vectors.
       */
      CHECK_OKAY( colLink(col, memhdr, set, lp) );

      var = col->var;
      assert(var != NULL);
      assert(var->varstatus == SCIP_VARSTATUS_COLUMN);
      assert(var->data.col == col);

      lp->lpicols[c] = col;
      col->lpipos = c;
      col->primsol = SCIP_INVALID;
      col->redcost = SCIP_INVALID;
      col->lbchanged = FALSE;
      col->ubchanged = FALSE;
      col->coefchanged = FALSE;
      obj[c] = var->obj;
      lb[c] = var->dom.lb;
      ub[c] = var->dom.ub;
      beg[c] = nnonz;
      name[c] = var->name;

      for( i = 0; i < col->len; ++i )
      {
         lpipos = col->row[i]->lpipos;
         if( lpipos >= 0 )
         {
            assert(lpipos < lp->nrows);
            ind[nnonz] = lpipos;
            val[nnonz] = col->val[i];
            nnonz++;
         }
      }
   }

   /* call LP interface */
   CHECK_OKAY( SCIPlpiAddCols(lp->lpi, naddcols, nnonz, obj, lb, ub, beg, ind, val, name) );
   lp->nlpicols = lp->ncols;
   lp->lpifirstchgcol = lp->nlpicols;

   /* free temporary memory */
   freeBlockMemoryArray(memhdr, obj, tmpsizecols);
   freeBlockMemoryArray(memhdr, lb, tmpsizecols);
   freeBlockMemoryArray(memhdr, ub, tmpsizecols);
   freeBlockMemoryArray(memhdr, beg, tmpsizecols);
   freeBlockMemoryArray(memhdr, ind, tmpsizecoefs);
   freeBlockMemoryArray(memhdr, val, tmpsizecoefs);
   freeBlockMemoryArray(memhdr, name, tmpsizecols);

   return SCIP_OKAY;
}

static
RETCODE lpFlushDelRows(                 /**< applies all cached row removals to the LP solver */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   assert(lp->lpifirstchgrow <= lp->nlpirows);
      
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

      CHECK_OKAY( SCIPlpiDelRows(lp->lpi, lp->lpifirstchgrow, lp->nlpirows-1) );
      for( i = lp->lpifirstchgrow; i < lp->nlpirows; ++i )
      {
         assert(!lp->lpirows[i]->inLP);
         lp->lpirows[i]->lpipos = -1;
         lp->lpirows[i]->dualsol = 0.0;
         lp->lpirows[i]->slack = SCIP_INVALID;
      }
      lp->nlpirows = lp->lpifirstchgrow;
   }
   assert(lp->nlpirows == lp->lpifirstchgrow);

   return SCIP_OKAY;
}

static
RETCODE lpFlushAddRows(                 /**< applies all cached row additions and removals to the LP solver */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   Real* rhs;
   char* sen;
   int* beg;
   int* ind;
   Real* val;
   char** name;
   ROW* row;
   int r;
   int nnonz;
   int naddrows;
   int naddcoefs;
   int i;
   int lpipos;
   int tmpsizerows;
   int tmpsizecoefs;

   assert(lp != NULL);
   assert(lp->lpifirstchgrow == lp->nlpirows);
   assert(memhdr != NULL);
      
   /* if there are no rows to add, we are ready */
   if( lp->nrows == lp->nlpirows )
      return SCIP_OKAY;

   /* add the additional rows */
   assert(lp->nrows > lp->nlpirows);
   CHECK_OKAY( ensureLpirowsSize(lp, set, lp->nrows) );

   /* count the (maximal) number of added coefficients, calculate the number of added rows */
   naddrows = lp->nrows - lp->nlpirows;
   naddcoefs = 0;
   for( r = lp->nlpirows; r < lp->nrows; ++r )
      naddcoefs += lp->rows[r]->len;
   assert(naddrows > 0);

   /* get temporary memory for changes */
   tmpsizerows = SCIPsetCalcMemGrowSize(set, naddrows);     /* use standard sizes to reuse memory more often */
   tmpsizecoefs = SCIPsetCalcMemGrowSize(set, naddcoefs);
   assert(tmpsizerows >= naddrows);
   assert(tmpsizecoefs >= naddcoefs);
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, rhs, tmpsizerows) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, sen, tmpsizerows) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, beg, tmpsizerows) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, ind, tmpsizecoefs) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, val, tmpsizecoefs) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, name, tmpsizerows) );
   
   /* fill temporary memory with row data */
   nnonz = 0;
   for( r = lp->nlpirows; r < lp->nrows; ++r )
   {
      row = lp->rows[r];
      assert(row != NULL);
      assert(row->inLP);
      assert(nnonz + row->len <= naddcoefs);

      /* Because the row becomes a member of the LP solver, its dual variable now can take values
       * different from zero. That means, we have to include the row in the corresponding
       * column vectors.
       */
      CHECK_OKAY( rowLink(row, memhdr, set, lp) );

      lp->lpirows[r] = row;
      row->lpipos = r;
      row->dualsol = SCIP_INVALID;
      row->slack = SCIP_INVALID;
      row->coefchanged = FALSE;
      rhs[r] = row->rhs;
      sen[r] = row->equality ? 'E' : 'L';
      beg[r] = nnonz;
      name[r] = row->name;

      for( i = 0; i < row->len; ++i )
      {
         lpipos = row->col[i]->lpipos;
         if( lpipos >= 0 )
         {
            assert(lpipos < lp->ncols);
            ind[nnonz] = lpipos;
            val[nnonz] = row->val[i];
            nnonz++;
         }
      }
   }

   /* call LP interface */
   CHECK_OKAY( SCIPlpiAddRows(lp->lpi, naddrows, nnonz, rhs, sen, beg, ind, val, name) );
   lp->nlpirows = lp->nrows;
   lp->lpifirstchgrow = lp->nlpirows;

   /* free temporary memory */
   freeBlockMemoryArray(memhdr, rhs, tmpsizerows);
   freeBlockMemoryArray(memhdr, sen, tmpsizerows);
   freeBlockMemoryArray(memhdr, beg, tmpsizerows);
   freeBlockMemoryArray(memhdr, ind, tmpsizecoefs);
   freeBlockMemoryArray(memhdr, val, tmpsizecoefs);
   freeBlockMemoryArray(memhdr, name, tmpsizerows);
   
   return SCIP_OKAY;
}

static
RETCODE lpFlushChgbds(                  /**< applies all cached bound changes to the LP */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   COL* col;
   const VAR* var;
   int* ind;
   char* lu;
   Real* bd;
   int i;
   int nchg;
   int tmpsizecols;

   assert(lp != NULL);
   assert(memhdr != NULL);

   if( lp->nchgbds == 0 )
      return SCIP_OKAY;

   /* get temporary memory for changes */
   tmpsizecols = SCIPsetCalcMemGrowSize(set, 2*lp->ncols);     /* use standard sizes to reuse memory more often */
   assert(tmpsizecols >= 2*lp->ncols);
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, ind, tmpsizecols) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, lu, tmpsizecols) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, bd, tmpsizecols) );

   /* collect all cached bound changes */
   nchg = 0;
   for( i = 0; i < lp->nchgbds; ++i )
   {
      col = lp->chgbds[i];
      assert(col != NULL);

      var = col->var;
      assert(var != NULL);
      assert(var->varstatus == SCIP_VARSTATUS_COLUMN);
      assert(var->data.col == col);

      if( col->lpipos >= 0 )
      {
         if( col->lbchanged )
         {
            assert(nchg < 2*lp->ncols);
            ind[nchg] = col->lpipos;
            lu[nchg] = 'L';
            bd[nchg] = var->dom.lb;
            nchg++;
         }
         if( col->ubchanged )
         {
            assert(nchg < 2*lp->ncols);
            ind[nchg] = col->lpipos;
            lu[nchg] = 'U';
            bd[nchg] = var->dom.ub;
            nchg++;
         }
      }
   }

   /* change bounds in LP */
   if( nchg > 0 )
   {
      CHECK_OKAY( SCIPlpiChgBd(lp->lpi, nchg, ind, lu, bd) );
   }

   lp->nchgbds = 0;

   /* free temporary memory */
   freeBlockMemoryArray(memhdr, ind, tmpsizecols);
   freeBlockMemoryArray(memhdr, lu, tmpsizecols);
   freeBlockMemoryArray(memhdr, bd, tmpsizecols);

   return SCIP_OKAY;
}

static
RETCODE lpFlush(                        /**< applies all cached changes to the LP solver */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(lp != NULL);
   assert(memhdr != NULL);
   
   if( lp->flushed )
   {
      assert(lp->nlpicols == lp->ncols);
      assert(lp->lpifirstchgcol == lp->nlpicols);
      assert(lp->nlpirows == lp->nrows);
      assert(lp->lpifirstchgrow == lp->nlpirows);
      assert(lp->nchgbds == 0);

      return SCIP_OKAY;
   }

   CHECK_OKAY( lpFlushDelCols(lp) );
   CHECK_OKAY( lpFlushDelRows(lp) );
   CHECK_OKAY( lpFlushChgbds(lp, set, memhdr) );
   CHECK_OKAY( lpFlushAddCols(lp, set, memhdr) );
   CHECK_OKAY( lpFlushAddRows(lp, set, memhdr) );

   lp->flushed = TRUE;

   return SCIP_OKAY;
}






static
void colAddSign(                        /**< update column sign after addition of new coefficient */
   COL*             col,                /**< LP column */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value of new coefficient */
   )
{   
   assert(col != NULL);
   assert(col->numpos >= 0 && col->numneg >= 0);
   assert(!SCIPsetIsZero(set, val));

   if( SCIPsetIsPos(set, val) )
      col->numpos++;
   else
   {
      assert(SCIPsetIsNeg(set, val));
      col->numneg++;
   }
}

static
void colDelSign(                        /**< update column sign after deletion of coefficient */
   COL*             col,                /**< LP column */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value of deleted coefficient */
   )
{
   assert(col != NULL);
   assert(!SCIPsetIsZero(set, val));

   if( SCIPsetIsPos(set, val) )
      col->numpos--;
   else
   {
      assert(SCIPsetIsNeg(set, val));
      col->numneg--;
   }

   assert(col->numpos >= 0 && col->numneg >= 0);
}

static
void rowAddNorms(                       /**< update row norms after addition of new coefficient */
   ROW*             row,                /**< LP row */
   const SET*       set,                /**< global SCIP settings */
   int              colidx,             /**< column index of new coefficient */
   Real             val                 /**< value of new coefficient */
   )
{
   Real absval;

   assert(row != NULL);
   assert(row->nummaxval >= 0);
   assert(set != NULL);

   absval = ABS(val);
   assert(!SCIPsetIsZero(set, absval));

   /* update min/maxidx */
   row->minidx = MIN(row->minidx, colidx);
   row->maxidx = MAX(row->maxidx, colidx);

   /* update squared euclidean norm */
   row->sqrnorm += SQR(absval);

   /* update maximum norm */
   if( SCIPsetIsG(set, absval, row->maxval) )
   {
      row->maxval = absval;
      row->nummaxval = 1;
   }
   else if( SCIPsetIsGE(set, absval, row->maxval) )
   {
      assert(row->nummaxval >= 1);
      row->nummaxval++;
   }
}

static
void rowCalcNorms(                      /**< calculates row norms and min/maxidx from scratch, and checks for sortation */
   ROW*             row,                /**< LP row */
   const SET*       set                 /**< global SCIP settings */
   )
{
   int i;
   int idx;

   assert(row != NULL);
   assert(set != NULL);

   row->sqrnorm = 0.0;
   row->maxval = 0.0;
   row->nummaxval = 0;
   row->minidx = INT_MAX;
   row->maxidx = INT_MIN;
   row->validminmaxidx = TRUE;
   row->sorted = TRUE;

   /* check, if row is sorted
    * calculate sqrnorm, maxval, minidx, and maxidx
    */
   for( i = 0; i < row->len; ++i )
   {
      assert(!SCIPsetIsZero(set, row->val[i]));
      idx = row->col[i]->index;
      rowAddNorms(row, set, idx, row->val[i]);
      row->sorted &= (i == 0 || row->col[i-1]->index < idx);
   }
}

static
void rowDelNorms(                       /**< update row norms after deletion of coefficient */
   ROW*             row,                /**< LP row */
   const SET*       set,                /**< global SCIP settings */
   int              colidx,             /**< column index of deleted coefficient */
   Real             val                 /**< value of deleted coefficient */
   )
{
   Real absval;

   assert(row != NULL);
   assert(row->nummaxval > 0);
   assert(set != NULL);

   absval = ABS(val);
   assert(!SCIPsetIsZero(set, absval));
   assert(SCIPsetIsGE(set, row->maxval, absval));

   /* update min/maxidx validity */
   if( colidx == row->minidx || colidx == row->maxidx )
      row->validminmaxidx = FALSE;

   /* update squared euclidean norm */
   row->sqrnorm -= SQR(absval);
   assert(SCIPsetIsGE(set, row->sqrnorm, 0.0));

   /* update maximum norm */
   if( SCIPsetIsGE(set, absval, row->maxval) )
   {
      row->nummaxval--;
      if( row->nummaxval == 0 )
         rowCalcNorms(row, set);
   }
}

static
DECL_SORTPTRCOMP(cmpCol)
{
   return ((COL*)elem1)->index - ((COL*)elem2)->index;
}

static
DECL_SORTPTRCOMP(cmpRow)
{
   return ((ROW*)elem1)->index - ((ROW*)elem2)->index;
}

void SCIPcolSort(                       /**< sorts column entries by row index */
   COL* col                             /**< column to be sorted */
   )
{
   if( !col->sorted )
   {
      SCIPbsortPtrDbl((void**)(col->row), col->val, col->len, &cmpCol);
      col->sorted = TRUE;
   }
}

void SCIProwSort(                       /**< sorts row entries by column index */
   ROW* row                             /**< row to be sorted */
   )
{
   if( !row->sorted )
   {
      SCIPbsortPtrDbl((void**)(row->col), row->val, row->len, &cmpRow);
      row->sorted = TRUE;
   }
}

static
int colSearchCoeff(                     /**< searches coefficient in column, returns position in col vector or -1 */
   COL* col,                            /**< column to be searched in */
   const ROW* row                       /**< coefficient to be searched for */
   )
{
   int actpos;
   int minpos;
   int maxpos;
   int actidx;
   int searchidx;

   assert(col != NULL);
   assert(row != NULL);

   /* row has to be sorted, such that binary search works */
   if( !col->sorted )
      SCIPcolSort(col);
   assert(col->sorted);

   /* binary search */
   searchidx = row->index;
   minpos = 0;
   maxpos = col->len-1;
   actpos = 0;
   while(minpos <= maxpos)
   {
      assert(0 <= actpos && actpos < col->len);
      actpos = (minpos + maxpos)/2;
      actidx = col->row[actpos]->index;
      if( searchidx == actidx )
         return actpos;
      else if( searchidx < actidx )
         maxpos = actpos-1;
      else
         minpos = actpos+1;
   }

   return -1;
}

static
int rowSearchCoeff(                     /**< searches coefficient in row, returns position in row vector or -1 */
   ROW* row,                            /**< row to be searched in */
   const COL* col                       /**< coefficient to be searched for */
   )
{
   int actpos;
   int minpos;
   int maxpos;
   int actidx;
   int searchidx;

   assert(row != NULL);
   assert(col != NULL);

   /* row has to be sorted, such that binary search works */
   if( !row->sorted )
      SCIProwSort(row);
   assert(row->sorted);

   /* binary search */
   searchidx = col->index;
   minpos = 0;
   maxpos = row->len-1;
   actpos = 0;
   while(minpos <= maxpos)
   {
      assert(0 <= actpos && actpos < row->len);
      actpos = (minpos + maxpos)/2;
      actidx = row->col[actpos]->index;
      if( searchidx == actidx )
         return actpos;
      else if( searchidx < actidx )
         maxpos = actpos-1;
      else
         minpos = actpos+1;
   }

   return -1;
}

static
RETCODE ensureColSize(                  /**< ensures, that row array of column can store at least num additional entries */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   COL*             col,                /**< LP column */
   int              num                 /**< minimum number of additional entries to store */
   )
{
   assert(col->len <= col->size);
   
   if( num > col->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, col->row, col->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, col->val, col->size, newsize) );
      col->size = newsize;
   }
   assert(num <= col->size);

   return SCIP_OKAY;
}

static
RETCODE ensureRowSize(                  /**< ensures, that column array of row can store at least num additional entries */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row,                /**< LP row */
   int              num                 /**< minimum number of additional entries to store */
   )
{
   assert(row->len <= row->size);
   
   if( num > row->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, row->col, row->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, row->val, row->size, newsize) );
      row->size = newsize;
   }
   assert(num <= row->size);

   return SCIP_OKAY;
}

RETCODE SCIPcolAddCoeff(                /**< adds a previously non existing coefficient to an LP column */
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   )
{
   assert(memhdr != NULL);
   assert(col != NULL);
   assert(row != NULL);
   assert(!SCIPsetIsZero(set, val));
   assert(colSearchCoeff(col, row) == -1);

   if( col->len > 0 )
      col->sorted &= (col->row[col->len-1]->index < row->index);

   CHECK_OKAY( ensureColSize(memhdr, set, col, col->len+1) );
   assert(col->row != NULL);
   assert(col->val != NULL);

   col->row[col->len] = row;
   col->val[col->len] = val;
   col->len++;

   colAddSign(col, set, val);
   coefChanged(row, col, lp);

   return SCIP_OKAY;
}

RETCODE SCIProwAddCoeff(                /**< adds a previously non existing coefficient to an LP row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
   Real             val                 /**< value of coefficient */
   )
{
   assert(memhdr != NULL);
   assert(row != NULL);
   assert(col != NULL);
   assert(!SCIPsetIsZero(set, val));
   assert(rowSearchCoeff(row, col) == -1);

   if( row->len > 0 )
      row->sorted &= (row->col[row->len-1]->index < col->index);

   CHECK_OKAY( ensureRowSize(memhdr, set, row, row->len+1) );
   assert(row->col != NULL);
   assert(row->val != NULL);

   row->col[row->len] = col;
   row->val[row->len] = val;
   row->len++;

   rowAddNorms(row, set, col->index, val);
   coefChanged(row, col, lp);

   return SCIP_OKAY;
}

static
void colDelCoeffPos(                    /**< deletes coefficient at given position from column */
   COL*             col,                /**< column to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   int              pos                 /**< position in column vector to delete */
   )
{
   ROW* row;
   Real val;

   assert(col != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < col->len);

   row = col->row[pos];
   val = col->val[pos];

   if( pos < col->len-1 )
   {
      /* move last coefficient to position of deleted coefficient */
      col->row[pos] = col->row[col->len-1];
      col->val[pos] = col->val[col->len-1];
      col->sorted = FALSE;
   }
   col->len--;

   colDelSign(col, set, val);
   coefChanged(row, col, lp);
}

void SCIPcolDelCoeff(                   /**< deletes existing coefficient from column */
   COL*             col,                /**< column to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(col != NULL);
   assert(col->len > 0);
   assert(row != NULL);

   pos = colSearchCoeff(col, row);
   assert(pos >= 0);
   assert(col->row[pos] == row);

   colDelCoeffPos(col, set, lp, pos);
}

static
void rowDelCoeffPos(                    /**< deletes coefficient at given position from row */
   ROW*             row,                /**< row to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   int              pos                 /**< position in row vector to delete */
   )
{
   COL* col;
   Real val;

   assert(row != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < row->len);

   col = row->col[pos];
   val = row->val[pos];

   if( pos < row->len-1 )
   {
      /* move last coefficient to position of deleted coefficient */
      row->col[pos] = row->col[row->len-1];
      row->val[pos] = row->val[row->len-1];
      row->sorted = FALSE;
   }
   row->len--;
   
   rowDelNorms(row, set, col->index, val);
   coefChanged(row, col, lp);
}

void SCIProwDelCoeff(                   /**< deletes coefficient from row */
   ROW*             row,                /**< row to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(row != NULL);
   assert(row->len > 0);
   assert(col != NULL);

   pos = rowSearchCoeff(row, col);
   assert(pos >= 0);
   assert(row->col[pos] == col);

   rowDelCoeffPos(row, set, lp, pos);
}

RETCODE SCIPcolChgCoeff(                /**< changes or adds a coefficient to an LP column */
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   )
{
   int pos;
   Bool isZero;

   assert(memhdr != NULL);
   assert(col != NULL);
   assert(row != NULL);

   pos = colSearchCoeff(col, row);
   isZero = SCIPsetIsZero(set, val);

   if( pos >= 0 )
   {
      if( isZero )
      {
         /* delete existing coefficient */
         colDelCoeffPos(col, set, lp, pos);
      }
      else if( !SCIPsetIsEQ(set, col->val[pos], val) )
      {
         /* change existing coefficient */
         col->val[pos] = val;
         coefChanged(row, col, lp);
      }
   }
   else if( !isZero )
   {
      /* add non existing coefficient */
      CHECK_OKAY( SCIPcolAddCoeff(col, memhdr, set, lp, row, val) );
   }

   return SCIP_OKAY;
}

RETCODE SCIProwChgCoeff(                /**< changes or adds a coefficient to an LP row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
   Real             val                 /**< value of coefficient */
   )
{
   int pos;
   Bool isZero;

   assert(memhdr != NULL);
   assert(row != NULL);
   assert(col != NULL);

   pos = rowSearchCoeff(row, col);
   isZero = SCIPsetIsZero(set, val);

   if( pos >= 0 )
   {
      if( isZero )
      {
         /* delete existing coefficient */
         rowDelCoeffPos(row, set, lp, pos);
      }
      else if( !SCIPsetIsEQ(set, row->val[pos], val) )
      {
         /* change existing coefficient */
         row->val[pos] = val;
         coefChanged(row, col, lp);
      }
   }
   else if( !isZero )
   {
      /* add non existing coefficient */
      CHECK_OKAY( SCIProwAddCoeff(row, memhdr, set, lp, col, val) );
   }

   return SCIP_OKAY;
}

RETCODE SCIPcolCreate(                  /**< creates an LP column */
   COL**            col,                /**< pointer to column data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   const VAR*       var,                /**< variable, this column represents */
   int              len,                /**< number of nonzeros in the column */
   ROW**            row,                /**< array with rows of column entries */
   Real*            val                 /**< array with coefficients of column entries */
   )
{
   int i;
   int idx;

   assert(col != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(stat != NULL);
   assert(var != NULL);
   assert(len >= 0);
   assert(len == 0 || (row != NULL && val != NULL));

   ALLOC_OKAY( allocBlockMemory(memhdr, *col) );

   if( len > 0 )
   {
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*col)->row, row, len) );
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*col)->val, val, len) );
   }
   else
   {
      (*col)->row = NULL;
      (*col)->val = NULL;
   }

   (*col)->var = var;
   (*col)->index = stat->numcolidx++;
   (*col)->size = len;
   (*col)->len = len;
   (*col)->lpipos = -1;
   (*col)->primsol = 0.0;
   (*col)->redcost = SCIP_INVALID;
   (*col)->numpos = 0;
   (*col)->numneg = 0;
   (*col)->sorted = TRUE;
   (*col)->lbchanged = FALSE;
   (*col)->ubchanged = FALSE;
   (*col)->coefchanged = FALSE;
   (*col)->linked = FALSE;
   (*col)->inLP = FALSE;

   /* check, if column is sorted
    * update number of positive/negative entries
    */
   for( i = 0; i < len; ++i )
   {
      assert(!SCIPsetIsZero(set, (*col)->val[i]));
      (*col)->sorted &= (i == 0 || (*col)->row[i-1]->index < (*col)->row[i]->index);
      colAddSign(*col, set, (*col)->val[i]);
   }

   return SCIP_OKAY;
}

void SCIPcolFree(                       /**< frees an LP column */
   COL**            col,                /**< pointer to LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   int i;

   assert(memhdr != NULL);
   assert(col != NULL);
   assert(*col != NULL);
   assert((*col)->var != NULL);
   assert((*col)->var->varstatus == SCIP_VARSTATUS_COLUMN);
   assert(&(*col)->var->data.col == col); /* SCIPcolFree() has to be called from SCIPvarFree() */

   /* remove column indices from corresponding rows */
   colUnlink(*col, memhdr, set, lp);

   freeBlockMemoryArrayNull(memhdr, (*col)->row, (*col)->size);
   freeBlockMemoryArrayNull(memhdr, (*col)->val, (*col)->size);
   freeBlockMemory(memhdr, *col);
}

void SCIPcolCalcRedcost(                /**< calculates the reduced costs of a column */
   COL*             col                 /**< LP column */
   )
{
   ROW* row;
   int r;

   assert(col != NULL);
   assert(col->var->varstatus == SCIP_VARSTATUS_COLUMN);
   assert(col->var->data.col == col);

   col->redcost = col->var->obj;
   for( r = 0; r < col->len; ++r )
   {
      row = col->row[r];
      assert(row->dualsol < SCIP_INVALID);
      col->redcost -= col->val[r] * row->dualsol;
   }
}

RETCODE SCIProwCreate(                  /**< creates an LP row */
   ROW**            row,                /**< pointer to LP row data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             rhs,                /**< right hand side of row */
   Real             epsilon,            /**< maximal normed violation of row */
   Bool             equality            /**< is row an equality? otherwise, it is a lower or equal inequality */
   )
{
   int i;
   int idx;

   assert(row != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);
   assert(len >= 0);
   assert(len == 0 || (col != NULL && val != NULL));
   assert(epsilon >= 0.0);

   ALLOC_OKAY( allocBlockMemory(memhdr, *row) );

   if( len > 0 )
   {
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*row)->col, col, len) );
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*row)->val, val, len) );
   }
   else
   {
      (*row)->col = NULL;
      (*row)->val = NULL;
   }
   
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*row)->name, name, strlen(name)+1) );
   (*row)->rhs = rhs;
   (*row)->epsilon = epsilon;
   (*row)->index = stat->numrowidx++;
   (*row)->size = len;
   (*row)->len = len;
   (*row)->numuses = 0;
   (*row)->lpipos = -1;
   (*row)->dualsol = 0.0;
   (*row)->slack = SCIP_INVALID;
   (*row)->coefchanged = FALSE;
   (*row)->equality = equality;
   (*row)->linked = FALSE;
   (*row)->inLP = FALSE;

   /* calculate row norms and min/maxidx, and check if row is sorted */
   rowCalcNorms(*row, set);

   return SCIP_OKAY;
}

void SCIProwFree(                       /**< frees an LP row */
   ROW**            row,                /**< pointer to LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   int i;

   assert(memhdr != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->numuses == 0);
   
   /* remove column indices from corresponding rows */
   rowUnlink(*row, memhdr, set, lp);

   freeBlockMemoryArrayNull(memhdr, (*row)->col, (*row)->size);
   freeBlockMemoryArrayNull(memhdr, (*row)->val, (*row)->size);
   freeBlockMemory(memhdr, *row);
}

void SCIProwCapture(                    /**< increases usage counter of LP row */
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->numuses >= 0);

   row->numuses++;
}

void SCIProwRelease(                    /**< decreases usage counter of LP row, and frees memory if necessary */
   ROW**            row,                /**< pointer to LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(memhdr != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->numuses >= 1);

   (*row)->numuses--;
   if( (*row)->numuses == 0 )
      SCIProwFree(row, memhdr, set, lp);
}

RETCODE SCIPcolBoundChanged(            /**< notifies LP column, that its bounds were changed */
   COL*             col,                /**< LP column */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   assert(col != NULL);
   assert(lp != NULL);

   /* insert column in the chgbds list (if not already there) */
   if( !col->lbchanged && !col->ubchanged )
   {
      CHECK_OKAY( ensureChgbdsSize(lp, set, lp->nchgbds+1) );
      lp->chgbds[lp->nchgbds] = col;
      lp->nchgbds++;
   }
   
   /* mark bound change in the column */
   switch( boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      col->lbchanged = TRUE;
      break;
   case SCIP_BOUNDTYPE_UPPER:
      col->ubchanged = TRUE;
      break;
   default:
      errorMessage("Unknown bound type");
      abort();
   }

   return SCIP_OKAY;
}


void SCIPlpMarkSize(                    /**< remembers number of columns and rows to track the newly added ones */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   lp->firstnewcol = lp->ncols;
   lp->firstnewrow = lp->nrows;
}

COL** SCIPlpGetNewcols(                 /**< get array with newly added columns after the last mark */
   const LP*        lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewcol && lp->firstnewcol <= lp->ncols);

   return &(lp->cols[lp->firstnewcol]);
}

int SCIPlpGetNumNewcols(                /**< get number of newly added columns after the last mark */
   const LP*        lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewcol && lp->firstnewcol <= lp->ncols);

   return lp->ncols - lp->firstnewcol;
}

ROW** SCIPlpGetNewrows(                 /**< get array with newly added rows after the last mark */
   const LP*        lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewrow && lp->firstnewrow <= lp->nrows);

   return &(lp->rows[lp->firstnewrow]);
}

int SCIPlpGetNumNewrows(                /**< get number of newly added rows after the last mark */
   const LP*        lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewrow && lp->firstnewrow <= lp->nrows);

   return lp->nrows - lp->firstnewrow;
}

RETCODE SCIPlpGetState(                 /**< stores LP state (like basis information) into LP state object */
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   LPSTATE**        lpstate             /**< pointer to LP state information (like basis information) */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(memhdr != NULL);
   assert(lpstate != NULL);

   return SCIPlpiGetState(lp->lpi, memhdr, lpstate);
}

RETCODE SCIPlpSetState(                 /**< loads LP state (like basis information) into solver */
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LPSTATE*         lpstate             /**< LP state information (like basis information) */
   )
{
   assert(lp != NULL);
   assert(memhdr != NULL);
   assert(lpstate != NULL);

   lpFlush(lp, set, memhdr);

   return SCIPlpiSetState(lp->lpi, memhdr, lpstate);
}

RETCODE SCIPlpAddCol(                   /**< adds a column to the LP */
   LP*              lp,                 /**< LP data */
   const SET*       set,                /**< global SCIP settings */
   COL*             col                 /**< LP column */
   )
{
   assert(lp != NULL);
   assert(col != NULL);
   assert(!col->inLP);

   CHECK_OKAY( ensureColsSize(lp, set, lp->ncols+1) );
   lp->cols[lp->ncols] = col;
   lp->ncols++;
   lp->flushed = FALSE;
   lp->solved = FALSE;
   lp->objval = SCIP_INVALID;
   col->inLP = TRUE;

   return SCIP_OKAY;
}

RETCODE SCIPlpAddRow(                   /**< adds a row to the LP */
   LP*              lp,                 /**< LP data */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row                 /**< LP row */
   )
{
   assert(lp != NULL);
   assert(row != NULL);
   assert(!row->inLP);

   CHECK_OKAY( ensureRowsSize(lp, set, lp->nrows+1) );
   lp->rows[lp->nrows] = row;
   lp->nrows++;
   lp->flushed = FALSE;
   lp->solved = FALSE;
   lp->objval = SCIP_INVALID;
   row->inLP = TRUE;

   return SCIP_OKAY;
}

RETCODE SCIPlpShrinkCols(               /**< removes all columns after the given number of columns from the LP */
   LP*              lp,                 /**< LP data */
   int              newncols            /**< new number of columns in the LP */
   )
{
   int c;

   assert(lp != NULL);
   assert(0 <= newncols && newncols <= lp->ncols);

   if( newncols < lp->ncols )
   {
      for( c = newncols; c < lp->ncols; ++c )
      {
         assert(lp->cols[c]->inLP);
         lp->cols[c]->inLP = FALSE;
      }
      lp->ncols = newncols;
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->objval = SCIP_INVALID;
   }

   return SCIP_OKAY;
}

RETCODE SCIPlpShrinkRows(               /**< removes all rows after the given number of rows from the LP */
   LP*              lp,                 /**< LP data */
   int              newnrows            /**< new number of rows in the LP */
   )
{
   int r;

   assert(lp != NULL);
   assert(0 <= newnrows && newnrows <= lp->nrows);

   if( newnrows < lp->nrows )
   {
      for( r = newnrows; r < lp->nrows; ++r )
      {
         assert(lp->rows[r]->inLP);
         lp->rows[r]->inLP = FALSE;
      }
      lp->nrows = newnrows;
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->objval = SCIP_INVALID;
   }

   return SCIP_OKAY;
}

RETCODE SCIPlpClear(                    /**< removes all columns and rows from LP */
   LP*              lp                  /**< LP data */
   )
{
   assert(lp != NULL);

   CHECK_OKAY( SCIPlpShrinkCols(lp, 0) );
   CHECK_OKAY( SCIPlpShrinkRows(lp, 0) );

   return SCIP_OKAY;
}

RETCODE SCIPlpCreate(                   /**< creates empty LP data object */
   LP**             lp,                 /**< pointer to LP data object */
   const char*      name                /**< problem name */
   )
{
   ALLOC_OKAY( allocMemory(*lp) );
   
   CHECK_OKAY( SCIPlpiOpen(&((*lp)->lpi), name) );

   (*lp)->lpicols = NULL;
   (*lp)->lpirows = NULL;
   (*lp)->chgbds = NULL;
   (*lp)->cols = NULL;
   (*lp)->rows = NULL;
   (*lp)->objoffset = 0.0;
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
   (*lp)->chgbdssize = 0;
   (*lp)->nchgbds = 0;
   (*lp)->firstnewcol = 0;
   (*lp)->firstnewrow = 0;
   (*lp)->flushed = FALSE;
   (*lp)->solved = FALSE;
   (*lp)->objval = SCIP_INVALID;

   return SCIP_OKAY;
}

RETCODE SCIPlpFree(                     /**< frees LP data object */
   LP**             lp                  /**< pointer to LP data object */
   )
{
   assert(lp != NULL);
   assert(*lp != NULL);

   if( (*lp)->lpi != NULL )
   {
      CHECK_OKAY( SCIPlpiClose(&((*lp)->lpi)) );
   }
   freeMemoryArrayNull((*lp)->lpicols);
   freeMemoryArrayNull((*lp)->lpirows);
   freeMemoryArrayNull((*lp)->chgbds);
   freeMemoryArrayNull((*lp)->cols);
   freeMemoryArrayNull((*lp)->rows);
   freeMemory(*lp);

   return SCIP_OKAY;
}

RETCODE SCIPlpSolvePrimal(              /**< solves the LP with the primal simplex algorithm */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   LPSOLSTAT*       lpsolstat           /**< pointer to store the LP solution status */
   )
{
   assert(lp != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);

   /* flush changes to the LP solver */
   CHECK_OKAY( lpFlush(lp, set, memhdr) );

   /* call primal simplex */
   CHECK_OKAY( SCIPlpiSolvePrimal(lp->lpi) );

   /* evaluate solution status */
   if( SCIPlpiIsOptimal(lp->lpi) )
      *lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   else if( SCIPlpiIsPrimalInfeasible(lp->lpi) )
      *lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;
   else if( SCIPlpiIsPrimalUnbounded(lp->lpi) )
      *lpsolstat = SCIP_LPSOLSTAT_UNBOUNDED;
   else if( SCIPlpiIsIterlimExc(lp->lpi) )
      *lpsolstat = SCIP_LPSOLSTAT_ITERLIMIT;
   else if( SCIPlpiIsTimelimExc(lp->lpi) )
      *lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
   else if( SCIPlpiIsObjlimExc(lp->lpi) )
   {
      errorMessage("Objective limit exceeded in primal simplex - this should not happen");
      *lpsolstat = SCIP_LPSOLSTAT_ERROR;
      return SCIP_LPERROR;
   }
   else if( SCIPlpiIsError(lp->lpi) )
   {
      errorMessage("Error in primal simplex");
      *lpsolstat = SCIP_LPSOLSTAT_ERROR;
      return SCIP_LPERROR;
   }
   else
   {
      errorMessage("Unknown return status of primal simplex");
      *lpsolstat = SCIP_LPSOLSTAT_ERROR;
      return SCIP_LPERROR;
   }

   lp->solved = TRUE;

   return SCIP_OKAY;
}

RETCODE SCIPlpSolveDual(                /**< solves the LP with the dual simplex algorithm */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   LPSOLSTAT*       lpsolstat           /**< pointer to store the LP solution status */
   )
{
   assert(lp != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);

   /* flush changes to the LP solver */
   CHECK_OKAY( lpFlush(lp, set, memhdr) );

   /* call primal simplex */
   CHECK_OKAY( SCIPlpiSolveDual(lp->lpi) );

   /* evaluate solution status */
   if( SCIPlpiIsOptimal(lp->lpi) )
      *lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   else if( SCIPlpiIsObjlimExc(lp->lpi) )
      *lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;
   else if( SCIPlpiIsPrimalInfeasible(lp->lpi) )
      *lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;
   else if( SCIPlpiIsPrimalUnbounded(lp->lpi) )
      *lpsolstat = SCIP_LPSOLSTAT_UNBOUNDED;
   else if( SCIPlpiIsIterlimExc(lp->lpi) )
      *lpsolstat = SCIP_LPSOLSTAT_ITERLIMIT;
   else if( SCIPlpiIsTimelimExc(lp->lpi) )
      *lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
   else if( SCIPlpiIsError(lp->lpi) )
   {
      errorMessage("Error in dual simplex");
      *lpsolstat = SCIP_LPSOLSTAT_ERROR;
      return SCIP_LPERROR;
   }
   else
   {
      errorMessage("Unknown return status of dual simplex");
      *lpsolstat = SCIP_LPSOLSTAT_ERROR;
      return SCIP_LPERROR;
   }

   lp->solved = TRUE;

   return SCIP_OKAY;
}

RETCODE SCIPlpGetSol(                   /**< stores the LP solution in the columns and rows */
   LP*              lp,                 /**< actual LP data */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr              /**< block memory buffers */
   )
{
   Real* primsol;
   Real* dualsol;
   Real* slack;
   Real* redcost;
   int colsize;
   int rowsize;
   int c;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(set != NULL);
   assert(memhdr != NULL);

   colsize = SCIPsetCalcMemGrowSize(set, lp->nlpicols);
   rowsize = SCIPsetCalcMemGrowSize(set, lp->nlpirows);

   ALLOC_OKAY( allocBlockMemoryArray(memhdr, primsol, colsize) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, dualsol, rowsize) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, slack, rowsize) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, redcost, colsize) );

   CHECK_OKAY( SCIPlpiGetSol(lp->lpi, &lp->objval, primsol, dualsol, slack, redcost) );

   debugMessage("LP solution: obj=%f\n", lp->objval);

   for( c = 0; c < lp->nlpicols; ++c )
   {
      lp->lpicols[c]->primsol = primsol[c];
      lp->lpicols[c]->redcost = redcost[c];
   }

   for( r = 0; r < lp->nlpirows; ++r )
   {
      lp->lpirows[r]->dualsol = dualsol[r];
      lp->lpirows[r]->slack = slack[r];
   }
   
   return SCIP_OKAY;
}
