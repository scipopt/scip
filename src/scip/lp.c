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

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <math.h>
#include <limits.h>

#include "stat.h"
#include "lp.h"


struct Hole                             /**< hole in a domain of an integer variable */
{
   int              first;              /**< first value of hole */
   int              last;               /**< last value of hole */
};

struct Holelist                         /**< list of domain holes */
{
   HOLE             hole;               /**< this hole */
   HOLELIST*        next;               /**< next hole in list */
};

struct Vardom                           /**< domain of a variable */
{
   HOLELIST*        holelist;           /**< list of holes (only for the integer variables) */
   double           lb;                 /**< lower bounds of variables */
   double           ub;                 /**< upper bounds of variables */
};

struct BoundChg                         /**< change in one bound of a variable */
{
   COL*             col;                /**< column to change the bounds for */
   double           newbound;           /**< new value for bound */
   double           oldbound;           /**< old value for bound */
   BOUNDTYPE        boundtype;          /**< type of bound: lower or upper bound */
};

struct HoleChg                          /**< change in a hole list */
{
   HOLELIST**       ptr;                /**< changed list pointer */
   HOLELIST*        newlist;            /**< new value of list pointer */
   HOLELIST*        oldlist;            /**< old value of list pointer */
};

struct VardomChg                        /**< tracks changes of the variable's domains */
{
   BOUNDCHG*        boundchg;           /**< array with changes in bounds of variables */
   HOLECHG*         holechg;            /**< array with changes in hole lists */
   int              nboundchg;          /**< number of bound changes */
   int              nholechg;           /**< number of hole list changes */
};

struct Col                              /**< variable of the problem and corresponding LP column */
{
   ROW**            row;                /**< rows of column entries */
   double*          val;                /**< coefficients of column entries */
   VARDOM           vardom;             /**< domain of variable */
   double           obj;                /**< objective function value of variable */
   int              index;              /**< consecutively numbered variable identifier */
   int              size;               /**< size of the row- and val-arrays */
   int              len;                /**< number of nonzeros in column */
   int              numuses;            /**< number of times, this column is referenced */
   int              lppos;              /**< column position number in LP, or -1 if not in LP */
   unsigned int     vartype:2;          /**< type of variable: binary, integer, implicit integer, continous */
   unsigned int     sorted:1;           /**< TRUE iff row indices are sorted in increasing order */
   unsigned int     validsign:1;        /**< TRUE iff positive- and negative-flag is valid */
   unsigned int     positive:1;         /**< TRUE iff column has any positive entries */
   unsigned int     negative:1;         /**< TRUE iff column has any negative entries */
};

struct ColList                          /**< list of columns */
{
   COL*             col;                /**< pointer to this column */
   COLLIST*         next;               /**< pointer to next collist entry */
};

struct Row                              /**< row of the LP */
{
   COL**            col;                /**< columns of row entries */
   double*          val;                /**< coefficients of row entries */
   double           rhs;                /**< right hand side of row */
   double           epsilon;            /**< maximal normed violation of row */
   double           eucnorm;            /**< euclidean norm of row vector */
   double           maxval;             /**< maximal absolute value of row vector */
   int              index;              /**< consecutively numbered row identifier */
   int              size;               /**< size of the col- and val-arrays */
   int              len;                /**< number of nonzeros in row */
   int              numuses;            /**< number of times, this row is referenced */
   int              lppos;              /**< row position number in LP, or -1 if not in LP */
   int              minidx;             /**< minimal column index of row entries */
   int              maxidx;             /**< maximal column index of row entries */
   unsigned int     equality:1;         /**< TRUE iff row is an equality, FALSE iff row is a less or equal inequality */
   unsigned int     sorted:1;           /**< TRUE iff column indices are sorted in increasing order */
   unsigned int     validminmaxidx:1;   /**< TRUE iff minimal and maximal column index is valid */
};

struct RowList                          /**< list of rows */
{
   ROW*             row;                /**< pointer to this row */
   ROWLIST*         next;               /**< pointer to next rowlist entry */
};


static
RETCODE ensureColSize(                  /**< ensures, that row array of column can store at least num entries */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   COL*             col,                /**< LP column */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(col->len <= col->size);
   
   if( col->len + num > col->size )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, col->len + num);
      ALLOC_OKAY( reallocBlockMemoryArray(mem->lpmem, col->row, col->size, newsize) );
      col->size = newsize;
   }
   assert(col->len + num <= col->size);

   return SCIP_OKAY;
}

static
RETCODE ensureRowSize(                  /**< ensures, that column array of row can store at least num entries */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row,                /**< LP row */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(set->memGrowAdd >= 1);
   assert(row->len <= row->size);
   
   if( row->len + num > row->size )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, row->len + num);
      ALLOC_OKAY( reallocBlockMemoryArray(mem->lpmem, row->col, row->size, newsize) );
      row->size = newsize;
   }
   assert(row->len + num <= row->size);

   return SCIP_OKAY;
}

static
RETCODE addColCoeff(                    /**< adds a coefficient to an LP column */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   COL*             col,                /**< LP column */
   ROW*             row,                /**< LP row */
   double           val                 /**< value of coefficient */
   )
{
   assert(mem != NULL);
   assert(col != NULL);
   assert(row != NULL);
   assert(ABS(val) > set->epsZero);

   CHECK_OKAY( ensureColSize(mem, set, col, col->len+1) );
   col->row[col->len] = row;
   col->val[col->len] = val;
   
   if( col->validsign )
   {
      col->positive |= (val > 0.0);
      col->negative |= (val < 0.0);
   }

   col->len++;

   return SCIP_OKAY;
}

static
RETCODE addRowCoeff(                    /**< adds a coefficient to an LP row */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row,                /**< LP row */
   COL*             col,                /**< LP column */
   double           val                 /**< value of coefficient */
   )
{
   assert(mem != NULL);
   assert(row != NULL);
   assert(col != NULL);
   assert(ABS(val) > set->epsZero);

   CHECK_OKAY( ensureRowSize(mem, set, row, row->len+1) );
   row->col[row->len] = col;
   row->val[row->len] = val;
   
   if( row->validminmaxidx )
   {
      row->minidx = MIN(row->minidx, col->index);
      row->maxidx = MAX(row->maxidx, col->index);
   }
   row->sorted &= (row->len == 0 || row->col[row->len-1]->index < col->index);

   row->len++;

   return SCIP_OKAY;
}

COL* SCIPcolCreate(                     /**< creates an LP column */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   int              len,                /**< number of nonzeros in the column */
   ROW**            row,                /**< array with rows of column entries */
   double*          val,                /**< array with coefficients of column entries */
   double           lb,                 /**< lower bound of variable */
   double           ub,                 /**< upper bound of variable */
   double           obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   )
{
   COL* col;
   int i;
   int idx;

   assert(mem != NULL);
   assert(stat != NULL);
   assert(len >= 0);
   assert(row != NULL);
   assert(val != NULL);

   ALLOC_NULL( allocBlockMemory(mem->lpmem, col) );
   ALLOC_NULL( duplicateBlockMemoryArray(mem->lpmem, col->row, row, len) );
   ALLOC_NULL( duplicateBlockMemoryArray(mem->lpmem, col->val, val, len) );
   col->obj = obj;
   col->index = stat->numcolidx++;
   col->size = len;
   col->len = len;
   col->numuses = 0;
   col->lppos = -1;
   col->vartype = vartype;
   col->sorted = TRUE;
   col->validsign = TRUE;
   col->positive = FALSE;
   col->negative = FALSE;

   /* check, if column is sorted
    * calculate sign of column
    * insert coefficients in corresponding rows */
   for( i = 0; i < len; ++i )
   {
      assert(ABS(col->val[i]) > set->epsZero);
      col->sorted &= (i == 0 || col->row[i-1]->index < col->row[i]->index);
      col->positive |= (col->val[i] > 0.0);
      col->negative |= (col->val[i] < 0.0);
      addRowCoeff(mem, set, col->row[i], col, col->val[i]);
   }

   return col;
}

void SCIPcolFree(                       /**< frees an LP column */
   COL**            col,                /**< pointer to LP column */
   MEM*             mem                 /**< block memory buffers */
   )
{
   int i;

   assert(mem != NULL);
   assert(col != NULL);
   assert(*col != NULL);
   assert((*col)->numuses == 0);
   
   /* delete coefficients in corresponding rows */
   for( i = 0; i < (*col)->len; ++i )
      deleteRowCoeff((*col)->row[i], *col);

   freeBlockMemoryArray(mem->lpmem, (*col)->row, (*col)->size);
   freeBlockMemoryArray(mem->lpmem, (*col)->val, (*col)->size);
   freeBlockMemory(mem->lpmem, *col);
}

void SCIPcolCapture(                    /**< increases usage counter of LP column */
   COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert(col->numuses >= 0);

   col->numuses++;
}

void SCIPcolRelease(                    /**< decreases usage counter of LP column, and frees memory if necessary */
   COL**            col,                /**< pointer to LP column */
   MEM*             mem                 /**< block memory buffers */
   )
{
   assert(mem != NULL);
   assert(col != NULL);
   assert(*col != NULL);
   assert((*col)->numuses >= 1);

   (*col)->numuses--;
   if( (*col)->numuses == 0 )
      SCIPcolFree(col, mem);
}

ROW* SCIProwCreate(                     /**< creates an LP row */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   double*          val,                /**< array with coefficients of row entries */
   Bool             equality,           /**< TRUE iff row is an equality, FALSE iff row is a less or equal inequality */
   double           rhs,                /**< right hand side of row */
   double           epsilon             /**< maximal normed violation of row */
   )
{
   ROW* row;
   int i;
   int idx;

   assert(mem != NULL);
   assert(stat != NULL);
   assert(len >= 0);
   assert(col != NULL);
   assert(val != NULL);
   assert(epsilon >= 0.0);

   ALLOC_NULL( allocBlockMemory(mem->lpmem, row) );
   ALLOC_NULL( duplicateBlockMemoryArray(mem->lpmem, row->col, col, len) );
   ALLOC_NULL( duplicateBlockMemoryArray(mem->lpmem, row->val, val, len) );
   row->rhs = rhs;
   row->epsilon = epsilon;
   row->eucnorm = 0.0;
   row->maxval = 0.0;
   row->index = stat->numrowidx++;
   row->size = len;
   row->len = len;
   row->numuses = 0;
   row->lppos = -1;
   row->minidx = INT_MAX;
   row->maxidx = INT_MIN;
   row->equality = equality;
   row->sorted = TRUE;
   row->validminmaxidx = TRUE;

   /* check, if row is sorted
    * calculate eucnorm, maxval, minidx, and maxidx
    * insert coefficients in corresponding columns */
   for( i = 0; i < len; ++i )
   {
      assert(ABS(row->val[i]) > set->epsZero);
      row->eucnorm += SQR(row->val[i]);
      row->maxval = MAX(row->maxval, ABS(row->val[i]));
      idx = row->col[i]->index;
      row->minidx = MIN(row->minidx, idx);
      row->maxidx = MAX(row->maxidx, idx);
      row->sorted &= (i == 0 || row->col[i-1]->index < idx);
      addColCoeff(mem, set, row->col[i], row, row->val[i]);
   }
   row->eucnorm = SQRT(row->eucnorm);

   return row;
}

void SCIProwFree(                       /**< frees an LP row */
   ROW**            row,                /**< pointer to LP row */
   MEM*             mem                 /**< block memory buffers */
   )
{
   int i;

   assert(mem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->numuses == 0);
   
   /* delete coefficients in corresponding columns */
   for( i = 0; i < (*row)->len; ++i )
      deleteColCoeff((*row)->col[i], *row);

   freeBlockMemoryArray(mem->lpmem, (*row)->col, (*row)->size);
   freeBlockMemoryArray(mem->lpmem, (*row)->val, (*row)->size);
   freeBlockMemory(mem->lpmem, *row);
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
   MEM*             mem                 /**< block memory buffers */
   )
{
   assert(mem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->numuses >= 1);

   (*row)->numuses--;
   if( (*row)->numuses == 0 )
      SCIProwFree(row, mem);
}

