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
 * @brief  LP management datastructures and methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <math.h>

#include "lp.h"


struct Row
{
   int*             ind;                /**< column indices of row entries */
   double*          val;                /**< coefficients of row entries */
   double           rhs;                /**< right hand side of row */
   double           epsilon;            /**< maximal normed violation of row */
   double           eucnorm;            /**< euclidean norm of row vector */
   double           maxval;             /**< maximal absolute value of row vector */
   int              size;               /**< size of the ind- and val-arrays */
   int              len;                /**< number of nonzeros in row */
   int              numuses;            /**< number of times, this row is referenced */
   unsigned int     equality:1;         /**< TRUE iff row is an equality, FALSE iff row is a less or equal inequality */
   unsigned int     inLP:1;             /**< TRUE iff row is included in actual LP */
   unsigned int     sorted:1;           /**< TRUE iff column indices are sorted in increasing order */
};

struct RowList
{
   ROW*             row;                /**< pointer to this row */
   ROWLIST*         next;               /**< pointer to next rowlist entry */
};


ROW* SCIPcreateRow(                     /**< creates an LP row */
   MEM*             mem,                /**< block memory buffers */
   int              len,                /**< number of nonzeros in the row */
   int*             ind,                /**< column indices of row entries */
   double*          val,                /**< coefficients of row entries */
   Bool             equality,           /**< TRUE iff row is an equality, FALSE iff row is a less or equal inequality */
   double           rhs,                /**< right hand side of row */
   double           epsilon             /**< maximal normed violation of row */
   )
{
   ROW* row;
   int i;

   assert(mem != NULL);
   assert(len >= 0);
   assert(ind != NULL);
   assert(val != NULL);
   assert(epsilon >= 0.0);

   ALLOC_NULL( allocBlockMemory(mem->lpmem, row) );
   ALLOC_NULL( duplicateBlockMemoryArray(mem->lpmem, row->ind, ind, len) );
   ALLOC_NULL( duplicateBlockMemoryArray(mem->lpmem, row->val, val, len) );
   row->rhs = rhs;
   row->epsilon = epsilon;
   row->size = len;
   row->len = len;
   row->numuses = 0;
   row->equality = equality;
   row->inLP = FALSE;

   row->eucnorm = 0.0;
   row->maxval = 0.0;
   row->sorted = TRUE;
   for( i = 0; i < len; ++i )
   {
      row->eucnorm += SQR(row->val[i]);
      row->maxval = MAX(row->maxval, ABS(row->val[i]));
      row->sorted &= (i == 0 || row->ind[i-1] < row->ind[i]);
   }
   row->eucnorm = SQRT(row->eucnorm);

   return row;
}

void SCIPfreeRow(                       /**< frees an LP row */
   MEM*             mem,                /**< block memory buffers */
   ROW**            row                 /**< pointer to LP row */
   )
{
   assert(mem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->numuses == 0);
   
   freeBlockMemoryArray(mem->lpmem, (*row)->ind, (*row)->size);
   freeBlockMemoryArray(mem->lpmem, (*row)->val, (*row)->size);
   freeBlockMemory(mem->lpmem, *row);
}

void SCIPcaptureRow(                    /**< increases usage counter of LP row */
   ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->numuses >= 0);

   row->numuses++;
}

void SCIPreleaseRow(                    /**< decreases usage counter of LP row, and frees memory if necessary */
   MEM*             mem,                /**< block memory buffers */
   ROW**            row                 /**< pointer to LP row */
   )
{
   assert(mem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->numuses >= 1);

   (*row)->numuses--;
   if( (*row)->numuses == 0 )
      SCIPfreeRow(mem, row);
}

