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
RETCODE ensureChgcolsSize(              /**< ensures, that chgcols array can store at least num entries */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nchgcols <= lp->chgcolssize);
   
   if( num > lp->chgcolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(lp->chgcols, newsize) );
      lp->chgcolssize = newsize;
   }
   assert(num <= lp->chgcolssize);

   return SCIP_OKAY;
}

static
RETCODE ensureChgrowsSize(               /**< ensures, that chgrows array can store at least num entries */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nchgrows <= lp->chgrowssize);
   
   if( num > lp->chgrowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(lp->chgrows, newsize) );
      lp->chgrowssize = newsize;
   }
   assert(num <= lp->chgrowssize);

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
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, col->linked, col->size, newsize) );
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
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, row->linked, row->size, newsize) );
      row->size = newsize;
   }
   assert(num <= row->size);

   return SCIP_OKAY;
}



/*
 * double linked coefficient matrix methods 
 */

/* forward references, because these methods call each other cyclic */

static
RETCODE rowAddCoeffLinked(              /**< adds a previously non existing coefficient to an LP row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
   Real             val,                /**< value of coefficient */
   Bool             linkupdate          /**< is this an update of the link? */
   );

static
void rowDelCoeffLinked(                 /**< deletes coefficient from row */
   ROW*             row,                /**< row to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< coefficient to be deleted */
   Bool             linkupdate          /**< is this an update of the link? */
   );

static
RETCODE colAddCoeffLinked(              /**< adds a previously non existing coefficient to an LP column */
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< LP row */
   Real             val,                /**< value of coefficient */
   Bool             linkupdate          /**< is this an update of the link? */
   );

static
void colDelCoeffLinked(                 /**< deletes existing coefficient from column */
   COL*             col,                /**< column to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< coefficient to be deleted */
   Bool             linkupdate          /**< is this an update of the link? */
   );

static
RETCODE rowChgCoeffLinked(              /**< changes or adds a coefficient to an LP row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
   Real             val,                /**< value of coefficient */
   Bool             linkupdate          /**< is this an update of the link? */
   );

static
RETCODE rowIncCoeffLinked(              /**< increases value of an existing or nonexisting coefficient in an LP row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
   Real             incval,             /**< value to add to the coefficient */
   Bool             linkupdate          /**< is this an update of the link? */
   );


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

   if( !col->alllinked )
   {
      debugMessage("linking column <%s>\n", col->var->name);
      for( i = 0; i < col->len; ++i )
      {
         assert(!SCIPsetIsZero(set, col->val[i]));
         if( !col->linked[i] )
         {
            CHECK_OKAY( rowAddCoeffLinked(col->row[i], memhdr, set, lp, col, col->val[i], TRUE) );
            col->linked[i] = TRUE;
         }
      }
      col->alllinked = TRUE;
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

   debugMessage("unlinking column <%s>\n", col->var->name);
   for( i = 0; i < col->len; ++i )
   {
      if( col->linked[i] )
         rowDelCoeffLinked(col->row[i], set, lp, col, TRUE);
      col->linked[i] = FALSE;
   }
   col->alllinked = FALSE;
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

   if( !row->alllinked )
   {
      debugMessage("linking row <%s>\n", row->name);
      for( i = 0; i < row->len; ++i )
      {
         assert(!SCIPsetIsZero(set, row->val[i]));
         if( !row->linked[i] )
         {
            CHECK_OKAY( colAddCoeffLinked(row->col[i], memhdr, set, lp, row, row->val[i], TRUE) );
            row->linked[i] = TRUE;
         }
      }
      row->alllinked = TRUE;
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

   debugMessage("unlinking row <%s>\n", row->name);
   for( i = 0; i < row->len; ++i )
   {
      if( row->linked[i] )
         colDelCoeffLinked(row->col[i], set, lp, row, TRUE);
      row->linked[i] = FALSE;
   }
   row->alllinked = FALSE;
}



/*
 * Changing announcements
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
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }
}

   


/*
 * Column methods
 */

static
DECL_SORTPTRCOMP(cmpCol)
{
   return ((COL*)elem1)->index - ((COL*)elem2)->index;
}

static
int colSearchCoeff(                     /**< searches coefficient in column, returns position in col vector or -1 */
   COL*             col,                /**< column to be searched in */
   const ROW*       row                 /**< coefficient to be searched for */
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
RETCODE colAddCoeffLinked(              /**< adds a previously non existing coefficient to an LP column */
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< LP row */
   Real             val,                /**< value of coefficient */
   Bool             linkupdate          /**< is this an update of the link? */
   )
{
   assert(memhdr != NULL);
   assert(col != NULL);
   assert(row != NULL);
   assert(!SCIPsetIsZero(set, val));
   assert(colSearchCoeff(col, row) == -1);

   debugMessage("adding coefficient %g * <%s> to column <%s>, linkupdate=%d\n", 
      val, row->name, col->var->name, linkupdate);

   if( col->len > 0 )
      col->sorted &= (col->row[col->len-1]->index < row->index);

   CHECK_OKAY( ensureColSize(memhdr, set, col, col->len+1) );
   assert(col->row != NULL);
   assert(col->val != NULL);
   assert(col->linked != NULL);

   col->row[col->len] = row;
   col->val[col->len] = val;
   col->linked[col->len] = linkupdate;
   col->alllinked &= linkupdate;
   col->len++;

   colAddSign(col, set, val);

   if( !linkupdate )
      coefChanged(row, col, lp);
   
   return SCIP_OKAY;
}

static
void colDelCoeffPos(                    /**< deletes coefficient at given position from column */
   COL*             col,                /**< column to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   int              pos,                /**< position in column vector to delete */
   Bool             linkupdate          /**< is this an update of the link? */
   )
{
   ROW* row;
   Real val;

   assert(col != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < col->len);

   row = col->row[pos];
   val = col->val[pos];

   debugMessage("deleting coefficient %g * <%s> from column <%s>, linkupdate=%d\n", 
      val, row->name, col->var->name, linkupdate);

   /* update the corresponding row vector, if necessary */
   if( !linkupdate && col->linked[pos] )
   {
      rowDelCoeffLinked(col->row[pos], set, lp, col, TRUE);
   }

   if( pos < col->len-1 )
   {
      /* move last coefficient to position of deleted coefficient */
      col->row[pos] = col->row[col->len-1];
      col->val[pos] = col->val[col->len-1];
      col->linked[pos] = col->linked[col->len-1];
      col->sorted = FALSE;
   }
   col->len--;

   colDelSign(col, set, val);
   
   if( !linkupdate )
      coefChanged(row, col, lp);
}

static
void colDelCoeffLinked(                 /**< deletes existing coefficient from column */
   COL*             col,                /**< column to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< coefficient to be deleted */
   Bool             linkupdate          /**< is this an update of the link? */
   )
{
   int pos;

   assert(col != NULL);
   assert(col->len > 0);
   assert(row != NULL);

   pos = colSearchCoeff(col, row);
   assert(pos >= 0);
   assert(col->row[pos] == row);

   colDelCoeffPos(col, set, lp, pos, linkupdate);
}

static
RETCODE colChgCoeffLinked(              /**< changes or adds a coefficient to an LP column */
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< LP row */
   Real             val,                /**< value of coefficient */
   Bool             linkupdate          /**< is this an update of the link? */
   )
{
   int pos;
   Bool isZero;

   assert(memhdr != NULL);
   assert(col != NULL);
   assert(row != NULL);

   debugMessage("changing coefficient <%s> of column <%s> to %g, linkupdate=%d\n", 
      row->name, col->var->name, val, linkupdate);

   pos = colSearchCoeff(col, row);
   isZero = SCIPsetIsZero(set, val);

   if( pos >= 0 )
   {
      if( !linkupdate && col->linked[pos] )
      {
         CHECK_OKAY( rowChgCoeffLinked(col->row[pos], memhdr, set, lp, col, val, TRUE) );
      }

      if( isZero )
      {
         /* delete existing coefficient */
         colDelCoeffPos(col, set, lp, pos, linkupdate);
      }
      else if( !SCIPsetIsEQ(set, col->val[pos], val) )
      {
         /* change existing coefficient */
         col->val[pos] = val;
         if( !linkupdate )
            coefChanged(row, col, lp);
      }
   }
   else if( !isZero )
   {
      /* add non existing coefficient */
      CHECK_OKAY( colAddCoeffLinked(col, memhdr, set, lp, row, val, linkupdate) );
   }

   return SCIP_OKAY;
}

static
RETCODE colIncCoeffLinked(              /**< increases value of an existing or nonexisting coefficient in an LP column */
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< LP row */
   Real             incval,             /**< value to add to the coefficient */
   Bool             linkupdate          /**< is this an update of the link? */
   )
{
   int pos;

   assert(memhdr != NULL);
   assert(col != NULL);
   assert(row != NULL);

   debugMessage("increasing coefficient <%s> of column <%s> with %g, linkupdate=%d\n", 
      row->name, col->var->name, incval, linkupdate);

   if( SCIPsetIsZero(set, incval) )
      return SCIP_OKAY;

   pos = colSearchCoeff(col, row);

   if( pos >= 0 )
   {
      if( !linkupdate && col->linked[pos] )
      {
         CHECK_OKAY( rowIncCoeffLinked(col->row[pos], memhdr, set, lp, col, incval, TRUE) );
      }

      if( SCIPsetIsZero(set, col->val[pos] + incval) )
      {
         /* delete existing coefficient */
         colDelCoeffPos(col, set, lp, pos, linkupdate);
      }
      else
      {
         /* increase existing coefficient */
         col->val[pos] += incval;
         if( !linkupdate )
            coefChanged(row, col, lp);
      }
   }
   else
   {
      /* add non existing coefficient */
      CHECK_OKAY( colAddCoeffLinked(col, memhdr, set, lp, row, incval, linkupdate) );
   }

   return SCIP_OKAY;
}


RETCODE SCIPcolCreate(                  /**< creates an LP column */
   COL**            col,                /**< pointer to column data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   VAR*             var,                /**< variable, this column represents */
   int              len,                /**< number of nonzeros in the column */
   ROW**            row,                /**< array with rows of column entries */
   Real*            val                 /**< array with coefficients of column entries */
   )
{
   int i;

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
      int i;

      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*col)->row, row, len) );
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*col)->val, val, len) );
      ALLOC_OKAY( allocBlockMemoryArray(memhdr, (*col)->linked, len) );
      for( i = 0; i < len; ++i )
         (*row)->linked[i] = FALSE;
   }
   else
   {
      (*col)->row = NULL;
      (*col)->val = NULL;
      (*col)->linked = NULL;
   }

   (*col)->var = var;
   (*col)->index = stat->ncolidx++;
   (*col)->size = len;
   (*col)->len = len;
   (*col)->lpipos = -1;
   (*col)->primsol = 0.0;
   (*col)->redcost = SCIP_INVALID;
   (*col)->farkas = SCIP_INVALID;
   (*col)->numpos = 0;
   (*col)->numneg = 0;
   (*col)->sorted = TRUE;
   (*col)->lbchanged = FALSE;
   (*col)->ubchanged = FALSE;
   (*col)->coefchanged = FALSE;
   (*col)->alllinked = (len == 0);
   (*col)->inlp = FALSE;

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
   assert(memhdr != NULL);
   assert(col != NULL);
   assert(*col != NULL);
   assert((*col)->var != NULL);
   assert((*col)->var->varstatus == SCIP_VARSTATUS_COLUMN);
   assert(&(*col)->var->data.col == col); /* SCIPcolFree() has to be called from SCIPvarFree() */
   assert(!(*col)->inlp);

   /* remove column indices from corresponding rows */
   colUnlink(*col, memhdr, set, lp);

   freeBlockMemoryArrayNull(memhdr, (*col)->row, (*col)->size);
   freeBlockMemoryArrayNull(memhdr, (*col)->val, (*col)->size);
   freeBlockMemoryArrayNull(memhdr, (*col)->linked, (*col)->size);
   freeBlockMemory(memhdr, *col);
}

void SCIPcolSort(                       /**< sorts column entries by row index */
   COL*             col                 /**< column to be sorted */
   )
{
   if( !col->sorted )
   {
      SCIPbsortPtrDbl((void**)(col->row), col->val, col->len, &cmpCol);
      col->sorted = TRUE;
   }
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
   return colAddCoeffLinked(col, memhdr, set, lp, row, val, FALSE);
}

void SCIPcolDelCoeff(                   /**< deletes existing coefficient from column */
   COL*             col,                /**< column to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row                 /**< coefficient to be deleted */
   )
{
   colDelCoeffLinked(col, set, lp, row, FALSE);
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
   return colChgCoeffLinked(col, memhdr, set, lp, row, val, FALSE);
}

RETCODE SCIPcolIncCoeff(                /**< increases value of an existing or nonexisting coefficient in an LP column */
   COL*             col,                /**< LP column */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< LP row */
   Real             incval              /**< value to add to the coefficient */
   )
{
   return colIncCoeffLinked(col, memhdr, set, lp, row, incval, FALSE);
}

RETCODE SCIPcolBoundChanged(            /**< notifies LP, that the bounds of a column were changed */
   COL*             col,                /**< LP column that changed */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   assert(col != NULL);
   assert(lp != NULL);
   
   if( col->lpipos >= 0 )
   {
      /* insert column in the chgcols list (if not already there) */
      if( !col->lbchanged && !col->ubchanged )
      {
         CHECK_OKAY( ensureChgcolsSize(lp, set, lp->nchgcols+1) );
         lp->chgcols[lp->nchgcols] = col;
         lp->nchgcols++;
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
         return SCIP_INVALIDDATA;
      }
      
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->objval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

      assert(lp->nchgcols > 0);
   }  

   return SCIP_OKAY;
}

static
void colCalcRedcost(                    /**< calculates the reduced costs of a column */
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

Real SCIPcolGetRedcost(                 /**< gets the reduced costs of a column in last LP or after recalculation */
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   if( col->redcost >= SCIP_INVALID )
      colCalcRedcost(col);
   assert(col->redcost < SCIP_INVALID);

   return col->redcost;
}

Real SCIPcolGetFeasibility(             /**< gets the feasibility of a column in last LP or after recalculation */
   COL*             col                 /**< LP column */
   )
{
   Real redcost;

   assert(col != NULL);
   assert(col->var != NULL);

   redcost = SCIPcolGetRedcost(col);

   if( col->var->dom.lb < 0 )
      return -ABS(redcost);
   else
      return redcost;
}

static
void colCalcFarkas(                     /**< calculates the farkas value of a column */
   COL*             col                 /**< LP column */
   )
{
   ROW* row;
   int r;

   assert(col != NULL);
   assert(col->var->varstatus == SCIP_VARSTATUS_COLUMN);
   assert(col->var->data.col == col);

   col->farkas = 0.0;
   for( r = 0; r < col->len; ++r )
   {
      row = col->row[r];
      assert(row->dualfarkas < SCIP_INVALID);
      col->farkas += col->val[r] * row->dualfarkas;
   }
   if( col->farkas > 0.0 )
      col->farkas *= col->var->dom.ub;
   else
      col->farkas *= col->var->dom.lb;
}

Real SCIPcolGetFarkas(                  /**< gets the farkas value of a column in last LP (which must be infeasible) */
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   if( col->farkas >= SCIP_INVALID )
      colCalcFarkas(col);
   assert(col->farkas < SCIP_INVALID);
   
   return col->farkas;
}

Bool SCIPcolIsInLP(                     /**< returns TRUE iff column is member of actual LP */
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->inlp;
}

void SCIPcolPrint(                      /**< output column to file stream */
   COL*             col,                /**< LP column */
   const SET*       set,                /**< global SCIP settings */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int r;

   assert(col != NULL);
   assert(col->var != NULL);

   if( file == NULL )
      file = stdout;

   /* print bounds */
   fprintf(file, "[%f,%f], ", col->var->dom.lb, col->var->dom.ub);

   /* print coefficients */
   if( col->len == 0 )
      fprintf(file, "<empty>");
   for( r = 0; r < col->len; ++r )
   {
      assert(col->row[r] != NULL);
      assert(col->row[r]->name != NULL);
      fprintf(file, "%+f%s ", col->val[r], col->row[r]->name);
   }
   fprintf(file, "\n");
}




/*
 * Row methods
 */

static
DECL_SORTPTRCOMP(cmpRow)
{
   return ((ROW*)elem1)->index - ((ROW*)elem2)->index;
}

static
int rowSearchCoeff(                     /**< searches coefficient in row, returns position in row vector or -1 */
   ROW*             row,                /**< row to be searched in */
   const COL*       col                 /**< coefficient to be searched for */
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
RETCODE rowAddCoeffLinked(              /**< adds a previously non existing coefficient to an LP row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
   Real             val,                /**< value of coefficient */
   Bool             linkupdate          /**< is this an update of the link? */
   )
{
   assert(row != NULL);
   assert(memhdr != NULL);
   assert(col != NULL);
   assert(!SCIPsetIsZero(set, val));
   assert(rowSearchCoeff(row, col) == -1);

   debugMessage("adding coefficient %g * <%s> to row <%s>, linkupdate=%d\n", val, col->var->name, row->name, linkupdate);

   if( row->len > 0 )
      row->sorted &= (row->col[row->len-1]->index < col->index);

   CHECK_OKAY( ensureRowSize(memhdr, set, row, row->len+1) );
   assert(row->col != NULL);
   assert(row->val != NULL);

   row->col[row->len] = col;
   row->val[row->len] = val;
   row->linked[row->len] = linkupdate;
   row->alllinked &= linkupdate;
   row->len++;

   rowAddNorms(row, set, col->index, val);

   if( !linkupdate )
      coefChanged(row, col, lp);

   return SCIP_OKAY;
}

static
void rowDelCoeffPos(                    /**< deletes coefficient at given position from row */
   ROW*             row,                /**< row to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   int              pos,                /**< position in row vector to delete */
   Bool             linkupdate          /**< is this an update of the link? */
   )
{
   COL* col;
   Real val;

   assert(row != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < row->len);

   col = row->col[pos];
   val = row->val[pos];

   debugMessage("deleting coefficient %g * <%s> from row <%s>, linkupdate=%d\n", 
      val, col->var->name, row->name, linkupdate);

   /* update the corresponding column vector, if necessary */
   if( !linkupdate && row->linked[pos] )
   {
      colDelCoeffLinked(row->col[pos], set, lp, row, TRUE);
   }

   if( pos < row->len-1 )
   {
      /* move last coefficient to position of deleted coefficient */
      row->col[pos] = row->col[row->len-1];
      row->val[pos] = row->val[row->len-1];
      row->linked[pos] = row->linked[row->len-1];
      row->sorted = FALSE;
   }
   row->len--;
   
   rowDelNorms(row, set, col->index, val);

   if( !linkupdate )
      coefChanged(row, col, lp);
}

static
void rowDelCoeffLinked(                 /**< deletes coefficient from row */
   ROW*             row,                /**< row to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< coefficient to be deleted */
   Bool             linkupdate          /**< is this an update of the link? */
   )
{
   int pos;

   assert(row != NULL);
   assert(row->len > 0);
   assert(col != NULL);

   pos = rowSearchCoeff(row, col);
   assert(pos >= 0);
   assert(row->col[pos] == col);

   rowDelCoeffPos(row, set, lp, pos, linkupdate);
}

static
RETCODE rowChgCoeffLinked(              /**< changes or adds a coefficient to an LP row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
   Real             val,                /**< value of coefficient */
   Bool             linkupdate          /**< is this an update of the link? */
   )
{
   int pos;
   Bool isZero;

   assert(memhdr != NULL);
   assert(row != NULL);
   assert(col != NULL);

   debugMessage("changing coefficient <%s> of row <%s> to %g, linkupdate=%d\n", 
      col->var->name, row->name, val, linkupdate);

   pos = rowSearchCoeff(row, col);
   isZero = SCIPsetIsZero(set, val);

   if( pos >= 0 )
   {
      if( !linkupdate && row->linked[pos] )
      {
         CHECK_OKAY( colChgCoeffLinked(row->col[pos], memhdr, set, lp, row, val, TRUE) );
      }

      if( isZero )
      {
         /* delete existing coefficient */
         rowDelCoeffPos(row, set, lp, pos, linkupdate);
      }
      else if( !SCIPsetIsEQ(set, row->val[pos], val) )
      {
         /* change existing coefficient */
         row->sqrnorm -= SQR(row->val[pos]);
         row->val[pos] = val;
         row->sqrnorm += SQR(row->val[pos]);
         if( !linkupdate )
            coefChanged(row, col, lp);
      }
   }
   else if( !isZero )
   {
      /* add non existing coefficient */
      CHECK_OKAY( rowAddCoeffLinked(row, memhdr, set, lp, col, val, linkupdate) );
   }

   return SCIP_OKAY;
}

static
RETCODE rowIncCoeffLinked(              /**< increases value of an existing or nonexisting coefficient in an LP row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
   Real             incval,             /**< value to add to the coefficient */
   Bool             linkupdate          /**< is this an update of the link? */
   )
{
   int pos;

   assert(memhdr != NULL);
   assert(row != NULL);
   assert(col != NULL);

   debugMessage("increasing coefficient <%s> of row <%s> with %g, linkupdate=%d\n", 
      col->var->name, row->name, incval, linkupdate);

   if( SCIPsetIsZero(set, incval) )
      return SCIP_OKAY;

   pos = rowSearchCoeff(row, col);

   if( pos >= 0 )
   {
      if( !linkupdate && row->linked[pos] )
      {
         CHECK_OKAY( colIncCoeffLinked(row->col[pos], memhdr, set, lp, row, incval, TRUE) );
      }

      if( SCIPsetIsZero(set, row->val[pos] + incval) )
      {
         /* delete existing coefficient */
         rowDelCoeffPos(row, set, lp, pos, linkupdate);
      }
      else
      {
         /* change existing coefficient */
         row->sqrnorm -= SQR(row->val[pos]);
         row->val[pos] += incval;
         row->sqrnorm += SQR(row->val[pos]);
         if( !linkupdate )
            coefChanged(row, col, lp);
      }
   }
   else
   {
      /* add non existing coefficient */
      CHECK_OKAY( rowAddCoeffLinked(row, memhdr, set, lp, col, incval, linkupdate) );
   }

   return SCIP_OKAY;
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
   Real             lhs,                /**< left hand side of row */
   Real             rhs                 /**< right hand side of row */
   )
{
   assert(row != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);
   assert(len >= 0);
   assert(len == 0 || (col != NULL && val != NULL));
   assert(lhs <= rhs);

   ALLOC_OKAY( allocBlockMemory(memhdr, *row) );

   if( len > 0 )
   {
      int i;

      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*row)->col, col, len) );
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*row)->val, val, len) );
      ALLOC_OKAY( allocBlockMemoryArray(memhdr, (*row)->linked, len) );
      for( i = 0; i < len; ++i )
         (*row)->linked[i] = FALSE;
   }
   else
   {
      (*row)->col = NULL;
      (*row)->val = NULL;
      (*row)->linked = NULL;
   }
   
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*row)->name, name, strlen(name)+1) );
   (*row)->lhs = lhs;
   (*row)->rhs = rhs;
   (*row)->sqrnorm = 0.0;
   (*row)->maxval = 0.0;
   (*row)->dualsol = 0.0;
   (*row)->slack = SCIP_INVALID;
   (*row)->dualfarkas = SCIP_INVALID;
   (*row)->index = stat->nrowidx++;
   (*row)->size = len;
   (*row)->len = len;
   (*row)->numuses = 0;
   (*row)->lpipos = -1;
   (*row)->minidx = INT_MAX;
   (*row)->maxidx = INT_MIN;
   (*row)->nummaxval = 0;
   (*row)->sorted = FALSE;
   (*row)->validminmaxidx = FALSE;
   (*row)->lhschanged = FALSE;
   (*row)->rhschanged = FALSE;
   (*row)->coefchanged = FALSE;
   (*row)->alllinked = (len == 0);
   (*row)->inlp = FALSE;

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
   assert(memhdr != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->numuses == 0);
   assert(!(*row)->inlp);

   /* remove column indices from corresponding rows */
   rowUnlink(*row, memhdr, set, lp);

   freeBlockMemoryArray(memhdr, (*row)->name, strlen((*row)->name)+1);
   freeBlockMemoryArrayNull(memhdr, (*row)->col, (*row)->size);
   freeBlockMemoryArrayNull(memhdr, (*row)->val, (*row)->size);
   freeBlockMemoryArrayNull(memhdr, (*row)->linked, (*row)->size);
   freeBlockMemory(memhdr, *row);
}

void SCIProwCapture(                    /**< increases usage counter of LP row */
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->numuses >= 0);

   debugMessage("capture row <%s> with numuses=%d\n", row->name, row->numuses);
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

   debugMessage("release row <%s> with numuses=%d\n", (*row)->name, (*row)->numuses);
   (*row)->numuses--;
   if( (*row)->numuses == 0 )
      SCIProwFree(row, memhdr, set, lp);

   *row = NULL;
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
   return rowAddCoeffLinked(row, memhdr, set, lp, col, val, FALSE);
}

void SCIProwDelCoeff(                   /**< deletes coefficient from row */
   ROW*             row,                /**< row to be changed */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col                 /**< coefficient to be deleted */
   )
{
   rowDelCoeffLinked(row, set, lp, col, FALSE);
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
   return rowChgCoeffLinked(row, memhdr, set, lp, col, val, FALSE);
}

RETCODE SCIProwIncCoeff(                /**< increases value of an existing or nonexisting coefficient in an LP row */
   ROW*             row,                /**< LP row */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   COL*             col,                /**< LP column */
   Real             incval              /**< value to add to the coefficient */
   )
{
   return rowIncCoeffLinked(row, memhdr, set, lp, col, incval, FALSE);
}

RETCODE SCIProwAddConst(                /**< add constant value to a row, i.e. subtract value from lhs and rhs */
   ROW*             row,                /**< LP row */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             constant            /**< constant value to add to the row */
   )
{
   assert(row != NULL);
   assert(row->lhs <= row->rhs);
   assert(!SCIPsetIsInfinity(set, ABS(row->rhs)));
   assert(!SCIPsetIsInfinity(set, ABS(constant)));

   if( !SCIPsetIsInfinity(set, -row->lhs) )
   {
      row->lhs -= constant;
      CHECK_OKAY( SCIProwSideChanged(row, set, lp, SCIP_SIDETYPE_LEFT) );
   }

   row->rhs -= constant;
   CHECK_OKAY( SCIProwSideChanged(row, set, lp, SCIP_SIDETYPE_RIGHT) );

   return SCIP_OKAY;
}

void SCIProwSort(                       /**< sorts row entries by column index */
   ROW*             row                 /**< row to be sorted */
   )
{
   if( !row->sorted )
   {
      SCIPbsortPtrDbl((void**)(row->col), row->val, row->len, &cmpRow);
      row->sorted = TRUE;
   }
}

void SCIProwCalcSlack(                  /**< recalculates the actual slack of a row */
   ROW*             row                 /**< LP row */
   )
{
   COL* col;
   int c;

   assert(row != NULL);

   row->slack = row->rhs;
   for( c = 0; c < row->len; ++c )
   {
      col = row->col[c];
      assert(col->primsol < SCIP_INVALID);
      row->slack -= row->val[c] * col->primsol;
   }
}

Real SCIProwGetSlack(                   /**< returns the slack of a row in the last LP solution or after recalculation */
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   if( row->slack >= SCIP_INVALID )
      SCIProwCalcSlack(row);
   assert(row->slack < SCIP_INVALID);

   return row->slack;
}

Real SCIProwGetFeasibility(             /**< returns the feasibility of a row in the last solution or after recalc */
   ROW*             row                 /**< LP row */
   )
{
   Real slack;

   assert(row != NULL);

   slack = SCIProwGetSlack(row);

   return MIN(slack, row->rhs - row->lhs - slack);
}

int SCIProwGetNNonz(                    /**< get number of nonzero entries in row vector */
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->len;
}

Real SCIProwGetNorm(                    /**< get euclidean norm of row vector */
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return sqrt(row->sqrnorm);
}

RETCODE SCIProwSideChanged(             /**< notifies LP row, that its sides were changed */
   ROW*             row,                /**< LP row */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
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
         errorMessage("Unknown row side type");
         abort();
      }
      
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->objval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

      assert(lp->nchgrows > 0);
   }

   return SCIP_OKAY;
}

const char* SCIProwGetName(             /**< returns the name of the row */
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->name;
}

Bool SCIProwIsInLP(                     /**< returns TRUE iff row is member of actual LP */
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->inlp;
}

void SCIProwPrint(                      /**< output row to file stream */
   ROW*             row,                /**< LP row */
   const SET*       set,                /**< global SCIP settings */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int c;

   assert(row != NULL);

   if( file == NULL )
      file = stdout;

   /* print left hand side */
   if( !SCIPsetIsInfinity(set, -row->lhs) )
      fprintf(file, "%+f <= ", row->lhs);

   /* print coefficients */
   if( row->len == 0 )
      fprintf(file, "0 ");
   for( c = 0; c < row->len; ++c )
   {
      assert(row->col[c] != NULL);
      assert(row->col[c]->var != NULL);
      assert(row->col[c]->var->name != NULL);
      assert(row->col[c]->var->varstatus == SCIP_VARSTATUS_COLUMN);
      fprintf(file, "%+f%s ", row->val[c], row->col[c]->var->name);
   }

   /* print right hand side */
   fprintf(file, "<= %+f\n", row->rhs);
}



/*
 * LP solver data update
 */

static
RETCODE lpFlushDelCols(                 /**< applies all cached column removals to the LP solver */
   LP*              lp                  /**< actual LP data */
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

      debugMessage("flushing col deletions: shrink LP from %d to %d colums\n", lp->nlpicols, lp->lpifirstchgcol);
      CHECK_OKAY( SCIPlpiDelCols(lp->lpi, lp->lpifirstchgcol, lp->nlpicols-1) );
      for( i = lp->lpifirstchgcol; i < lp->nlpicols; ++i )
      {
         assert(!lp->lpicols[i]->inlp);
         lp->lpicols[i]->lpipos = -1;
         lp->lpicols[i]->primsol = 0.0;
         lp->lpicols[i]->redcost = SCIP_INVALID;
         lp->lpicols[i]->farkas = SCIP_INVALID;
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
   int pos;
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
   for( pos = 0, c = lp->nlpicols; c < lp->ncols; ++pos, ++c )
   {
      col = lp->cols[c];
      assert(col != NULL);
      assert(col->inlp);
      assert(nnonz + col->len <= naddcoefs);

      debugMessage("flushing added column:");
      debug( SCIPcolPrint(col, set, NULL) );

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
      col->farkas = SCIP_INVALID;
      col->lbchanged = FALSE;
      col->ubchanged = FALSE;
      col->coefchanged = FALSE;
      obj[pos] = var->obj;
      lb[pos] = var->dom.lb;
      ub[pos] = var->dom.ub;
      beg[pos] = nnonz;
      name[pos] = var->name;

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
   debugMessage("flushing col additions: enlarge LP from %d to %d colums\n", lp->nlpicols, lp->ncols);
   CHECK_OKAY( SCIPlpiAddCols(lp->lpi, naddcols, nnonz, obj, lb, ub, beg, ind, val, name, set->infinity) );
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

      debugMessage("flushing row deletions: shrink LP from %d to %d rows\n", lp->nlpirows, lp->lpifirstchgrow);
      CHECK_OKAY( SCIPlpiDelRows(lp->lpi, lp->lpifirstchgrow, lp->nlpirows-1) );
      for( i = lp->lpifirstchgrow; i < lp->nlpirows; ++i )
      {
         assert(!lp->lpirows[i]->inlp);
         lp->lpirows[i]->lpipos = -1;
         lp->lpirows[i]->dualsol = 0.0;
         lp->lpirows[i]->slack = SCIP_INVALID;
         lp->lpirows[i]->dualfarkas = SCIP_INVALID;
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
   Real* lhs;
   Real* rhs;
   int* beg;
   int* ind;
   Real* val;
   char** name;
   ROW* row;
   int r;
   int pos;
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
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, lhs, tmpsizerows) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, rhs, tmpsizerows) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, beg, tmpsizerows) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, ind, tmpsizecoefs) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, val, tmpsizecoefs) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, name, tmpsizerows) );
   
   /* fill temporary memory with row data */
   nnonz = 0;
   for( pos = 0, r = lp->nlpirows; r < lp->nrows; ++pos, ++r )
   {
      row = lp->rows[r];
      assert(row != NULL);
      assert(row->inlp);
      assert(nnonz + row->len <= naddcoefs);

      debugMessage("flushing added row:");
      debug( SCIProwPrint(row, set, NULL) );

      /* Because the row becomes a member of the LP solver, its dual variable now can take values
       * different from zero. That means, we have to include the row in the corresponding
       * column vectors.
       */
      CHECK_OKAY( rowLink(row, memhdr, set, lp) );

      lp->lpirows[r] = row;
      row->lpipos = r;
      row->dualsol = SCIP_INVALID;
      row->slack = SCIP_INVALID;
      row->dualfarkas = SCIP_INVALID;
      row->lhschanged = FALSE;
      row->rhschanged = FALSE;
      row->coefchanged = FALSE;
      lhs[pos] = row->lhs;
      rhs[pos] = row->rhs;
      beg[pos] = nnonz;
      name[pos] = row->name;

      debugMessage("flushing added row (LPI): %+g <=", lhs[pos]);
      for( i = 0; i < row->len; ++i )
      {
         lpipos = row->col[i]->lpipos;
         debug( printf(" %+gx%d(<%s>)", row->val[i], lpipos+1, row->col[i]->var->name) );
         if( lpipos >= 0 )
         {
            assert(lpipos < lp->ncols);
            ind[nnonz] = lpipos;
            val[nnonz] = row->val[i];
            nnonz++;
         }
      }
      debug( printf(" <= %+g\n", rhs[pos]) );
   }

   /* call LP interface */
   debugMessage("flushing row additions: enlarge LP from %d to %d rows\n", lp->nlpirows, lp->nrows);
   CHECK_OKAY( SCIPlpiAddRows(lp->lpi, naddrows, nnonz, lhs, rhs, beg, ind, val, name, set->infinity) );
   lp->nlpirows = lp->nrows;
   lp->lpifirstchgrow = lp->nlpirows;

   /* free temporary memory */
   freeBlockMemoryArray(memhdr, lhs, tmpsizerows);
   freeBlockMemoryArray(memhdr, rhs, tmpsizerows);
   freeBlockMemoryArray(memhdr, beg, tmpsizerows);
   freeBlockMemoryArray(memhdr, ind, tmpsizecoefs);
   freeBlockMemoryArray(memhdr, val, tmpsizecoefs);
   freeBlockMemoryArray(memhdr, name, tmpsizerows);
   
   return SCIP_OKAY;
}

static
RETCODE lpFlushChgcols(                 /**< applies all cached column changes to the LP */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   COL* col;
   const VAR* var;
   int* ind;
   Real* lb;
   Real* ub;
   int i;
   int nchg;
   int tmpsizecols;

   assert(lp != NULL);
   assert(memhdr != NULL);

   if( lp->nchgcols == 0 )
      return SCIP_OKAY;

   /* get temporary memory for changes */
   tmpsizecols = SCIPsetCalcMemGrowSize(set, lp->ncols);   /* use standard sizes to reuse memory more often */
   assert(tmpsizecols >= lp->ncols);
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, ind, tmpsizecols) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, lb, tmpsizecols) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, ub, tmpsizecols) );

   /* collect all cached bound changes */
   nchg = 0;
   for( i = 0; i < lp->nchgcols; ++i )
   {
      col = lp->chgcols[i];
      assert(col != NULL);

      var = col->var;
      assert(var != NULL);
      assert(var->varstatus == SCIP_VARSTATUS_COLUMN);
      assert(var->data.col == col);

      if( col->lpipos >= 0 )
      {
         if( col->lbchanged || col->ubchanged )
         {
            assert(nchg < lp->ncols);
            ind[nchg] = col->lpipos;
            lb[nchg] = var->dom.lb;
            ub[nchg] = var->dom.ub;
            nchg++;
            col->lbchanged = FALSE;
            col->ubchanged = FALSE;
         }
      }
   }

   /* change sides in LP */
   if( nchg > 0 )
   {
      debugMessage("flushing bound changes: change %d bounds of %d columns\n", nchg, lp->nchgcols);
      CHECK_OKAY( SCIPlpiChgBounds(lp->lpi, nchg, ind, lb, ub, set->infinity) );
   }

   lp->nchgcols = 0;

   /* free temporary memory */
   freeBlockMemoryArray(memhdr, ind, tmpsizecols);
   freeBlockMemoryArray(memhdr, lb, tmpsizecols);
   freeBlockMemoryArray(memhdr, ub, tmpsizecols);

   return SCIP_OKAY;
}

static
RETCODE lpFlushChgrows(                 /**< applies all cached row changes to the LP */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   ROW* row;
   int* ind;
   Real* lhs;
   Real* rhs;
   int i;
   int nchg;
   int tmpsizerows;

   assert(lp != NULL);
   assert(memhdr != NULL);

   if( lp->nchgrows == 0 )
      return SCIP_OKAY;

   /* get temporary memory for changes */
   tmpsizerows = SCIPsetCalcMemGrowSize(set, lp->nrows);   /* use standard sizes to reuse memory more often */
   assert(tmpsizerows >= lp->nrows);
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, ind, tmpsizerows) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, lhs, tmpsizerows) );
   ALLOC_OKAY( allocBlockMemoryArray(memhdr, rhs, tmpsizerows) );

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
            lhs[nchg] = row->lhs;
            rhs[nchg] = row->rhs;
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
      CHECK_OKAY( SCIPlpiChgSides(lp->lpi, nchg, ind, lhs, rhs, set->infinity) );
   }

   lp->nchgrows = 0;

   /* free temporary memory */
   freeBlockMemoryArray(memhdr, ind, tmpsizerows);
   freeBlockMemoryArray(memhdr, lhs, tmpsizerows);
   freeBlockMemoryArray(memhdr, rhs, tmpsizerows);

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
   
   debugMessage("flushing LP changes: old (%d cols, %d rows), chgcol=%d, chgrow=%d, new (%d cols, %d rows)\n",
      lp->nlpicols, lp->nlpirows, lp->lpifirstchgcol, lp->lpifirstchgrow, lp->ncols, lp->nrows);
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

   CHECK_OKAY( lpFlushDelCols(lp) );
   CHECK_OKAY( lpFlushDelRows(lp) );
   CHECK_OKAY( lpFlushChgcols(lp, set, memhdr) );
   CHECK_OKAY( lpFlushChgrows(lp, set, memhdr) );
   CHECK_OKAY( lpFlushAddCols(lp, set, memhdr) );
   CHECK_OKAY( lpFlushAddRows(lp, set, memhdr) );

   lp->flushed = TRUE;

   return SCIP_OKAY;
}



/*
 * LP methods
 */

RETCODE SCIPlpCreate(                   /**< creates empty LP data object */
   LP**             lp,                 /**< pointer to LP data object */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< problem name */
   )
{
   assert(lp != NULL);
   assert(set != NULL);
   assert(name != NULL);

   ALLOC_OKAY( allocMemory(*lp) );

   /* open LP Solver interface */
   CHECK_OKAY( SCIPlpiCreate(&((*lp)->lpi), name) );

   (*lp)->lpicols = NULL;
   (*lp)->lpirows = NULL;
   (*lp)->chgcols = NULL;
   (*lp)->chgrows = NULL;
   (*lp)->cols = NULL;
   (*lp)->rows = NULL;
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
   (*lp)->flushed = TRUE;
   (*lp)->solved = TRUE;
   (*lp)->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   (*lp)->objval = 0.0;

   /* set default parameters in LP solver */
   CHECK_OKAY( SCIPlpSetFeastol(*lp, set->feastol) );

   return SCIP_OKAY;
}

RETCODE SCIPlpFree(                     /**< frees LP data object */
   LP**             lp,                 /**< pointer to LP data object */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(*lp != NULL);
   
   CHECK_OKAY( SCIPlpClear(*lp, memhdr, set) );

   if( (*lp)->lpi != NULL )
   {
      CHECK_OKAY( SCIPlpiFree(&((*lp)->lpi)) );
   }

   freeMemoryArrayNull((*lp)->lpicols);
   freeMemoryArrayNull((*lp)->lpirows);
   freeMemoryArrayNull((*lp)->chgcols);
   freeMemoryArrayNull((*lp)->cols);
   freeMemoryArrayNull((*lp)->rows);
   freeMemory(*lp);

   return SCIP_OKAY;
}

RETCODE SCIPlpAddCol(                   /**< adds a column to the LP */
   LP*              lp,                 /**< LP data */
   const SET*       set,                /**< global SCIP settings */
   COL*             col                 /**< LP column */
   )
{
   assert(lp != NULL);
   assert(col != NULL);
   assert(!col->inlp);
   assert(col->var != NULL);
   assert(col->var->varstatus == SCIP_VARSTATUS_COLUMN);
   assert(col->var->data.col == col);
   
   debugMessage("adding column <%s> to LP\n", col->var->name);
   CHECK_OKAY( ensureColsSize(lp, set, lp->ncols+1) );
   lp->cols[lp->ncols] = col;
   lp->ncols++;
   lp->flushed = FALSE;
   lp->solved = FALSE;
   lp->objval = SCIP_INVALID;
   lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   col->inlp = TRUE;

   return SCIP_OKAY;
}

RETCODE SCIPlpAddRow(                   /**< adds a row to the LP and captures it */
   LP*              lp,                 /**< LP data */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row                 /**< LP row */
   )
{
   assert(lp != NULL);
   assert(row != NULL);
   assert(!row->inlp);

   SCIProwCapture(row);

   debugMessage("adding row <%s> to LP\n", row->name);
   CHECK_OKAY( ensureRowsSize(lp, set, lp->nrows+1) );
   lp->rows[lp->nrows] = row;
   lp->nrows++;
   lp->flushed = FALSE;
   lp->solved = FALSE;
   lp->objval = SCIP_INVALID;
   lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   row->inlp = TRUE;

   return SCIP_OKAY;
}

RETCODE SCIPlpShrinkCols(               /**< removes all columns after the given number of cols from the LP */
   LP*              lp,                 /**< LP data */
   int              newncols            /**< new number of columns in the LP */
   )
{
   int c;

   assert(lp != NULL);
   assert(0 <= newncols && newncols <= lp->ncols);

   debugMessage("shrinking LP from %d to %d columns\n", lp->ncols, newncols);
   if( newncols < lp->ncols )
   {
      for( c = newncols; c < lp->ncols; ++c )
      {
         assert(lp->cols[c]->var != NULL);
         assert(lp->cols[c]->var->varstatus == SCIP_VARSTATUS_COLUMN);
         assert(lp->cols[c]->var->data.col == lp->cols[c]);
         assert(lp->cols[c]->inlp);
         
         lp->cols[c]->inlp = FALSE;
      }
      lp->ncols = newncols;
      lp->lpifirstchgcol = MIN(lp->lpifirstchgcol, newncols);
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->objval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   return SCIP_OKAY;
}

RETCODE SCIPlpShrinkRows(               /**< removes and releases all rows after the given number of rows from the LP */
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   int              newnrows            /**< new number of rows in the LP */
   )
{
   int r;

   assert(lp != NULL);
   assert(0 <= newnrows && newnrows <= lp->nrows);

   debugMessage("shrinking LP from %d to %d rows\n", lp->nrows, newnrows);
   if( newnrows < lp->nrows )
   {
      for( r = newnrows; r < lp->nrows; ++r )
      {
         assert(lp->rows[r]->inlp);
         lp->rows[r]->inlp = FALSE;
         SCIProwRelease(&lp->rows[r], memhdr, set, lp);
      }
      lp->nrows = newnrows;
      lp->lpifirstchgrow = MIN(lp->lpifirstchgrow, newnrows);
      lp->flushed = FALSE;
      lp->solved = FALSE;
      lp->objval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   return SCIP_OKAY;
}

RETCODE SCIPlpClear(                    /**< removes all columns and rows from LP, releases all rows */
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);

   debugMessage("clearing LP\n");
   CHECK_OKAY( SCIPlpShrinkCols(lp, 0) );
   CHECK_OKAY( SCIPlpShrinkRows(lp, memhdr, set, 0) );

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
   LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(memhdr != NULL);
   assert(lpistate != NULL);

   return SCIPlpiGetState(lp->lpi, memhdr, lpistate);
}

RETCODE SCIPlpSetState(                 /**< loads LP state (like basis information) into solver */
   LP*              lp,                 /**< LP data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{
   assert(lp != NULL);
   assert(memhdr != NULL);
   assert(lpistate != NULL);

   lpFlush(lp, set, memhdr);

   return SCIPlpiSetState(lp->lpi, memhdr, lpistate);
}

RETCODE SCIPlpSetFeastol(               /**< sets the feasibility tolerance of the LP solver */
   LP*              lp,                 /**< actual LP data */
   Real             feastol             /**< new feasibility tolerance */
   )
{
   assert(lp != NULL);
   assert(feastol >= 0.0);

   CHECK_OKAY( SCIPlpiSetRealpar(lp->lpi, SCIP_LPPAR_FEASTOL, feastol) );
   if( lp->nrows > 0 )
   {
      lp->solved = FALSE;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   return SCIP_OKAY;
}

RETCODE SCIPlpSetUpperbound(            /**< sets the upper objective limit of the LP solver */
   LP*              lp,                 /**< actual LP data */
   Real             upperbound          /**< new upper objective limit */
   )
{
   assert(lp != NULL);

   debugMessage("setting LP upper objective limit to %g\n", upperbound);
   CHECK_OKAY( SCIPlpiSetRealpar(lp->lpi, SCIP_LPPAR_UOBJLIM, upperbound) );

   return SCIP_OKAY;
}

RETCODE SCIPlpSolvePrimal(              /**< solves the LP with the primal simplex algorithm */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat                /**< problem statistics */
   )
{
   assert(lp != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(stat != NULL);

   debugMessage("solving primal LP %d (LP %d)\n", stat->nprimallp+1, stat->nlp+1);

   /* flush changes to the LP solver */
   CHECK_OKAY( lpFlush(lp, set, memhdr) );

   /* call primal simplex */
   CHECK_OKAY( SCIPlpiSolvePrimal(lp->lpi) );

   /* evaluate solution status */
   if( SCIPlpiIsOptimal(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
      CHECK_OKAY( SCIPlpiGetObjval(lp->lpi, &lp->objval) );
   }
   else if( SCIPlpiIsPrimalInfeasible(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;
      lp->objval = set->infinity;
   }
   else if( SCIPlpiIsPrimalUnbounded(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_UNBOUNDED;
      lp->objval = -set->infinity;
   }
   else if( SCIPlpiIsIterlimExc(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_ITERLIMIT;
      lp->objval = -set->infinity;
   }
   else if( SCIPlpiIsTimelimExc(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
      lp->objval = -set->infinity;
   }
   else if( SCIPlpiIsObjlimExc(lp->lpi) )
   {
      errorMessage("Objective limit exceeded in primal simplex - this should not happen");
      lp->lpsolstat = SCIP_LPSOLSTAT_ERROR;
      lp->objval = -set->infinity;
      return SCIP_LPERROR;
   }
   else
   {
      errorMessage("Unknown return status of primal simplex");
      lp->lpsolstat = SCIP_LPSOLSTAT_ERROR;
      return SCIP_LPERROR;
   }

   lp->solved = TRUE;

   stat->nlp++;
   stat->nprimallp++;

   debugMessage("solving primal LP returned solstat=%d\n", lp->lpsolstat);

   return SCIP_OKAY;
}

RETCODE SCIPlpSolveDual(                /**< solves the LP with the dual simplex algorithm */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat                /**< problem statistics */
   )
{
   assert(lp != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(stat != NULL);

   debugMessage("solving dual LP %d (LP %d)\n", stat->nduallp+1, stat->nlp+1);

   /* flush changes to the LP solver */
   CHECK_OKAY( lpFlush(lp, set, memhdr) );

   /* call primal simplex */
   CHECK_OKAY( SCIPlpiSolveDual(lp->lpi) );

   /* evaluate solution status */
   if( SCIPlpiIsOptimal(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
      CHECK_OKAY( SCIPlpiGetObjval(lp->lpi, &lp->objval) );
   }
   else if( SCIPlpiIsObjlimExc(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
      lp->objval = set->infinity;
   }
   else if( SCIPlpiIsPrimalInfeasible(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;
      lp->objval = set->infinity;
   }
   else if( SCIPlpiIsPrimalUnbounded(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_UNBOUNDED;
      lp->objval = -set->infinity;
   }
   else if( SCIPlpiIsIterlimExc(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_ITERLIMIT;
      CHECK_OKAY( SCIPlpiGetObjval(lp->lpi, &lp->objval) );
   }
   else if( SCIPlpiIsTimelimExc(lp->lpi) )
   {
      lp->lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
      CHECK_OKAY( SCIPlpiGetObjval(lp->lpi, &lp->objval) );
   }
   else
   {
      errorMessage("Unknown return status of dual simplex");
      lp->lpsolstat = SCIP_LPSOLSTAT_ERROR;
      lp->objval = -set->infinity;
      return SCIP_LPERROR;
   }

   lp->solved = TRUE;

   stat->nlp++;
   stat->nduallp++;

   debugMessage("solving dual LP returned solstat=%d\n", lp->lpsolstat);

   return SCIP_OKAY;
}

LPSOLSTAT SCIPlpGetSolstat(             /**< gets solution status of last solve call */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   assert(lp->solved || lp->lpsolstat == SCIP_LPSOLSTAT_NOTSOLVED);

   return lp->lpsolstat;
}

Real SCIPlpGetObjval(                   /**< gets objective value of last solution */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   assert(lp->solved);

   return lp->objval;
}


RETCODE SCIPlpGetSol(                   /**< stores the LP solution in the columns and rows */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
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
      debugMessage(" col <%s>: primsol=%f, redcost=%f\n", lp->lpicols[c]->var->name, primsol[c], redcost[c]);
      lp->lpicols[c]->primsol = primsol[c];
      lp->lpicols[c]->redcost = redcost[c];
   }

   for( r = 0; r < lp->nlpirows; ++r )
   {
      debugMessage(" row <%s>: dualsol=%f, slack=%f\n", lp->lpirows[r]->name, dualsol[r], slack[r]);
      lp->lpirows[r]->dualsol = dualsol[r];
      lp->lpirows[r]->slack = slack[r];
   }

   freeBlockMemoryArray(memhdr, primsol, colsize);
   freeBlockMemoryArray(memhdr, dualsol, rowsize);
   freeBlockMemoryArray(memhdr, slack, rowsize);
   freeBlockMemoryArray(memhdr, redcost, colsize);
   
   return SCIP_OKAY;
}

RETCODE SCIPlpGetDualfarkas(            /**< stores the dual farkas multipliers for infeasibility proof in rows */
   LP*              lp,                 /**< actual LP data */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr              /**< block memory buffers */
   )
{
   Real* dualfarkas;
   int rowsize;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(set != NULL);
   assert(memhdr != NULL);

   rowsize = SCIPsetCalcMemGrowSize(set, lp->nlpirows);

   ALLOC_OKAY( allocBlockMemoryArray(memhdr, dualfarkas, rowsize) );

   CHECK_OKAY( SCIPlpiGetDualfarkas(lp->lpi, dualfarkas) );

   debugMessage("LP is infeasible:\n");

   for( r = 0; r < lp->nlpirows; ++r )
   {
      debugMessage(" row <%s>: dualfarkas=%f\n", lp->lpirows[r]->name, dualfarkas[r]);
      lp->lpirows[r]->dualfarkas = dualfarkas[r];
   }

   freeBlockMemoryArray(memhdr, dualfarkas, rowsize);
   
   return SCIP_OKAY;
}
