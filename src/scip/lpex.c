


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lp.c
 * @brief  LP management methods and data structures
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gerald Gamrath
 *
 *  In LP management, we have to differ between the current LP and the SCIP_LP
 *  stored in the LP solver. All LP methods affect the current LP only.
 *  Before solving the current LP with the LP solver or setting an LP state,
 *  the LP solvers data has to be updated to the current LP with a call to
 *  lpFlush().
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "lpi/lpi.h"
#include "lpi/lpiex.h"
#include "scip/clock.h"
#include "scip/cons.h"
#include "scip/event.h"
#include "scip/intervalarith.h"
#include "scip/lp.h"
#include "scip/lpex.h"
#include "scip/misc.h"
#include "scip/prob.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_var.h"
#include "scip/pub_varex.h"
#include "scip/rational.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/solex.h"
#include "scip/solve.h"
#include "scip/stat.h"
#include "scip/struct_event.h"
#include "scip/struct_lpex.h"
#include "scip/struct_prob.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_var.h"
#include "scip/var.h"
#include "scip/varex.h"
#include <string.h>

static
SCIP_DECL_HASHKEYVAL(hashValExRow)
{
   /*lint --e{715}*/
   SCIP_ROW* row;
   int i;
   SCIP_Real scale;
   SCIP_SET* set;
   uint64_t hash;

   set = (SCIP_SET*) userptr;
   row = (SCIP_ROW*)key;
   assert(row != NULL);
   assert(row->len > 0);

   scale = 1.0 / SCIProwGetMaxval(row, set);
   if( SCIPsetIsInfinity(set, row->rhs) )
      scale = -scale;

   hash = (uint64_t) (long) row->len;

   for( i = 0; i < row->len; ++i )
   {
      SCIP_Real val = scale * row->vals[i];

      hash += SCIPhashTwo(SCIPrealHashCode(val), row->cols_index[i]);
   }

   return hash;
}

static
SCIP_DECL_HASHGETKEY(hashKeyExRow)
{
   /*lint --e{715}*/
   SCIP_ROWEX* row;

   row = (SCIP_ROWEX*)elem;
   assert(row != NULL);
   assert(row->fprow != NULL);

   /* the key of a cut is the row */
   return row->fprow;
}

static
SCIP_DECL_HASHKEYEQ(hashEqExRow)
{
   /*lint --e{715}*/
   SCIP_ROW* row1;
   SCIP_ROW* row2;

   row1 = (SCIP_ROW*)key1;
   row2 = (SCIP_ROW*)key2;
   assert(row1 != NULL);
   assert(row2 != NULL);

   /* return true if the row is the same */
   if( row1 == row2 )
      return TRUE;
   else
      return FALSE;
}

/** comparison method for sorting rows by non-decreasing index */
SCIP_DECL_SORTPTRCOMP(SCIProwexComp)
{
   assert(elem1 != NULL);
   assert(elem2 != NULL);

   assert(((SCIP_ROWEX*)elem1)->fprow != NULL);
   assert(((SCIP_ROWEX*)elem2)->fprow != NULL);

   if( ((SCIP_ROWEX*)elem1)->index < ((SCIP_ROWEX*)elem2)->index )
      return -1;
   else if( ((SCIP_ROWEX*)elem1)->index > ((SCIP_ROWEX*)elem2)->index )
      return +1;
   else
   {
      assert(SCIProwexGetIndex(((SCIP_ROWEX*)elem1))
         == SCIProwexGetIndex(((SCIP_ROWEX*)elem2)));
      return 0;
   }
}

#ifndef NDEBUG

/** checks if the exact column and its fpcol are consistent */
static
SCIP_Bool colexInSync(
   SCIP_COLEX*           colex,              /**< exact column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     msg                 /**< message handler */
   )
{
   int i;
   SCIP_COL* fpcol;

   assert(colex != NULL);

   fpcol = colex->fpcol;
   assert(fpcol != NULL);

   assert(colex->len == fpcol->len);
   assert(colex->size == fpcol->size);
   assert(colex->var == fpcol->var);

   assert(RisEqualReal(colex->obj, fpcol->obj));
   assert(RisEqualReal(colex->lb, fpcol->lb)
         || (RisNegInfinity(colex->lb) && SCIPsetIsInfinity(set, -fpcol->lb)));
   assert(RisEqualReal(colex->ub, fpcol->ub)
         || (RisInfinity(colex->ub) && SCIPsetIsInfinity(set, fpcol->ub)));

   for( i = 0; i < colex->len; ++i )
   {
      assert(RisEqualReal(colex->vals[i], fpcol->vals[i]));
   }

   return TRUE;
}

/** checks if the exact row and its fprow are consistent */
static
SCIP_Bool rowexInSync(
   SCIP_ROWEX*           rowex,              /**< exact row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     msg                 /**< message handler */
   )
{
   int i;
   SCIP_ROW* fprow;
   SCIP_Bool synced = TRUE;

   assert(rowex != NULL);

   fprow = rowex->fprow;

   assert(fprow != NULL);

   assert(rowex->len == fprow->len);

   synced = RisEqualReal(rowex->lhs, fprow->lhs)
            || (RisNegInfinity(rowex->lhs) && SCIPsetIsInfinity(set, -fprow->lhs));
   synced = synced && (RisEqualReal(rowex->rhs, fprow->rhs)
            || (RisInfinity(rowex->rhs) && SCIPsetIsInfinity(set, fprow->rhs)));
   synced = synced && (RisEqualReal(rowex->constant, fprow->constant) && rowex->origin == fprow->origin);
   if( !synced )
   {
      SCIPdebug(SCIProwPrint(rowex->fprow, msg, NULL));
      SCIPdebug(SCIProwexPrint(rowex, msg, NULL));
      SCIPABORT();
   }

   for( i = 0; i < rowex->len; ++i )
   {
      if( !RisEqualReal(rowex->vals[i], fprow->vals[i]) )
      {
            SCIProwPrint(rowex->fprow, msg, NULL);
            SCIProwexPrint(rowex, msg, NULL);
            SCIPABORT();
      }
   }

   for( i = 0; i < rowex->len; i++ )
   {
      if( !colexInSync(rowex->cols[i], set, msg) )
      {
         SCIPcolexPrint(rowex->cols[i], msg, NULL);
         SCIPcolPrint(fprow->cols[i], msg, NULL);
         SCIPABORT();
      }
   }

   return TRUE;
}

/** checks if the exact lp and lp are consistent */
static
SCIP_Bool lpexInSync(
   SCIP_LPEX*            lpex,               /**< exact lp */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     msg                 /**< message handler */
   )
{
   int i;
   SCIP_LP* fplp;
   assert(lpex != NULL);

   fplp = lpex->fplp;
   assert(fplp != NULL);

   assert(lpex->nrows == fplp->nrows);
   for( i = 0; i < lpex->nrows; i++)
   {
      assert(rowexInSync(lpex->rows[i], set, msg));
   }

   assert(lpex->ncols == fplp->ncols);
   for( i = 0; i < lpex->ncols; i++)
   {
      assert(colexInSync(lpex->cols[i], set, msg));
   }

   return TRUE;
}
#endif

/** ensures, that rows array can store at least num entries */
static
SCIP_RETCODE ensureRowexsSize(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lpex->nrows <= lpex->rowssize);

   if( num > lpex->rowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpex->rows, newsize) );
      lpex->rowssize = newsize;
   }
   assert(num <= lpex->rowssize);

   return SCIP_OKAY;
}

/** sorts column entries of linked rows currently in the LP such that lower row indices precede higher ones */
static
void colexSortLP(
   SCIP_COLEX*           col                 /**< column to be sorted */
   )
{
   int i;

   assert(col != NULL);

   /* check, if column is already sorted in the LP part */
   if( col->lprowssorted )
      return;

   /* sort coefficients */
   SCIPsortPtrPtrInt((void**)col->rows, (void**)col->vals, col->linkpos, SCIProwexComp, col->nlprows );

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
void colexSortNonLP(
   SCIP_COLEX*           col                   /**< column to be sorted */
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
                     &(col->linkpos[col->nlprows]), SCIProwexComp,
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
void rowexSortLP(
   SCIP_ROWEX*           row                 /**< row to be sorted */
   )
{
   int i;
   int pos;

   assert(row != NULL);

   /* check, if row is already sorted in the LP part, or if the sorting should be delayed */
   if( row->lpcolssorted || row->delaysort )
      return;

   /* sort coefficients */

   /** todo exip: add correct sorting methods so this works in one go */
   SCIPsortIntIntPtrPtr(row->cols_index, row->linkpos, (void**)row->cols,
                        (void**)row->vals, row->nlpcols);

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
void rowexSortNonLP(
   SCIP_ROWEX*           row                 /**< row to be sorted */
   )
{
   int i;
   assert(row != NULL);

   /* check, if row is already sorted in the non-LP part, or if the sorting should be delayed */
   if( row->nonlpcolssorted || row->delaysort )
      return;

/** todo exip: add correct sorting methods so this works in one go */
   /* sort coefficients */
   SCIPsortIntIntPtrPtr(&(row->cols_index[row->nlpcols]), &(row->linkpos[row->nlpcols]),
                        (void**)(&(row->cols[row->nlpcols])), (void**)&(row->vals[row->nlpcols]),
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

/** sorts row entries such that LP columns precede non-LP columns and inside both parts lower column indices precede
 *  higher ones
 */
void SCIProwexSort(
   SCIP_ROWEX*           row                 /**< row to be sorted */
   )
{
   assert(row != NULL);

   /* sort LP columns */
   rowexSortLP(row);

   /* sort non-LP columns */
   rowexSortNonLP(row);

#ifdef SCIP_MORE_DEBUG
   /* check the sorting */
   {
      int c;
      if( !row->delaysort )
      {
         for( c = 1; c < row->nlpcols; ++c )
            assert(row->cols[c]->index >= row->cols[c-1]->index);
         for( c = row->nlpcols + 1; c < row->len; ++c )
            assert(row->cols[c]->index >= row->cols[c-1]->index);
      }
   }
#endif
}

/** searches coefficient in part of the column, returns position in col vector or -1 if not found */
static
int colexSearchCoefPart(
   SCIP_COLEX*           col,                /**< column to be searched in */
   const SCIP_ROWEX*     row,                /**< coefficient to be searched for */
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
int colexSearchCoef(
   SCIP_COLEX*           col,                /**< column to be searched in */
   const SCIP_ROWEX*     row                 /**< coefficient to be searched for */
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
      colexSortLP(col);
      assert(col->lprowssorted);

      pos = colexSearchCoefPart(col, row, 0, col->nlprows-1);
      if( pos >= 0 )
         return pos;
   }

   /* search in the non-LP/unlinked rows */
   if( row->lppos == -1 || col->nunlinked > 0 )
   {
      /* column has to be sorted, such that binary search works */
      colexSortNonLP(col);
      assert(col->nonlprowssorted);

      pos = colexSearchCoefPart(col, row, col->nlprows, col->len-1);
   }

   return pos;
}

/** searches coefficient in part of the row, returns position in col vector or -1 if not found */
static
int rowexSearchCoefPart(
   SCIP_ROWEX*           row,                /**< row to be searched in */
   const SCIP_COLEX*     col,                /**< coefficient to be searched for */
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
int rowexSearchCoef(
   SCIP_ROWEX*           row,                /**< row to be searched in */
   const SCIP_COLEX*     col                 /**< coefficient to be searched for */
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
      rowexSortLP(row);
      assert(row->lpcolssorted);

      pos = rowexSearchCoefPart(row, col, 0, row->nlpcols-1);
   }

   /* search in the non-LP/unlinked columns */
   if( pos == -1 && (col->lppos == -1 || row->nunlinked > 0) )
   {
      /* row has to be sorted, such that binary search works */
      rowexSortNonLP(row);
      assert(row->nonlpcolssorted);

      pos = rowexSearchCoefPart(row, col, row->nlpcols, row->len-1);
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
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_COLEX*           col,                /**< LP col */
   SCIP_LPEX*            lp                  /**< current LP data */
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

   RsetString(row->maxactivity, "inf");
   RsetString(row->minactivity, "inf");
   RsetString(row->pseudoactivity, "inf");

}

/*
 * local column changing methods
 */

/** moves a coefficient in a column to a different place, and updates all corresponding data structures */
static
void colexMoveCoef(
   SCIP_COLEX*           col,                /**< LP column */
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

   Rset(col->vals[newpos], col->vals[oldpos]);
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
void colexSwapCoefs(
   SCIP_COLEX*           col,                /**< LP column */
   BMS_BUFMEM*           buffer,             /**< buffer for temp real */
   int                   pos1,               /**< position of first coefficient */
   int                   pos2                /**< position of second coefficient */
   )
{
   SCIP_ROWEX* tmprow;
   SCIP_Rational* tmpval;
   int tmplinkpos;

   assert(col != NULL);
   assert(0 <= pos1 && pos1 < col->len);
   assert(0 <= pos2 && pos2 < col->len);
   assert(col->rows[pos1] != NULL);

   if( pos1 == pos2 )
      return;

   tmpval = RcreateTemp(buffer);
   /* swap coefficients */
   tmprow = col->rows[pos2];
   Rset(tmpval, col->vals[pos2]);
   tmplinkpos = col->linkpos[pos2];

   col->rows[pos2] = col->rows[pos1];
   Rset(col->vals[pos2], col->vals[pos1]);
   col->linkpos[pos2] = col->linkpos[pos1];

   col->rows[pos1] = tmprow;
   Rset(col->vals[pos1], tmpval);
   col->linkpos[pos1] = tmplinkpos;

   RdeleteTemp(buffer, &tmpval);

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
}

/** moves a coefficient in a row to a different place, and updates all corresponding data structures */
static
void rowexMoveCoef(
   SCIP_ROWEX*           row,                /**< LP row */
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
   Rset(row->vals[newpos], row->vals[oldpos]);
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
void rowexSwapCoefs(
   SCIP_ROWEX*           row,                /**< LP row */
   BMS_BUFMEM*           buffer,             /**< buffer for temp real */
   int                   pos1,               /**< position of first coefficient */
   int                   pos2                /**< position of second coefficient */
   )
{
   SCIP_COLEX* tmpcol;
   SCIP_Rational* tmpval;
   int tmpindex;
   int tmplinkpos;

   assert(row != NULL);
   assert(0 <= pos1 && pos1 < row->len);
   assert(0 <= pos2 && pos2 < row->len);
   assert(row->cols[pos1] != NULL);
   assert(row->cols[pos1]->index == row->cols_index[pos1]);

   if( pos1 == pos2 )
      return;

   tmpval = RcreateTemp(buffer);
   /* swap coefficients */
   tmpcol = row->cols[pos2];
   tmpindex = row->cols_index[pos2];
   Rset(tmpval, row->vals[pos2]);
   tmplinkpos = row->linkpos[pos2];

   row->cols[pos2] = row->cols[pos1];
   row->cols_index[pos2] = row->cols_index[pos1];
   Rset(row->vals[pos2], row->vals[pos1]);
   row->linkpos[pos2] = row->linkpos[pos1];

   row->cols[pos1] = tmpcol;
   row->cols_index[pos1] = tmpindex;
   Rset(row->vals[pos1], tmpval);
   row->linkpos[pos1] = tmplinkpos;

   RdeleteTemp(buffer, &tmpval);

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

/** ensures, that row array of column can store at least num entries */
static
SCIP_RETCODE colexEnsureSize(
   SCIP_COLEX*           col,                /**< LP column */
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

      /* realloc fpcol */
      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->rows, col->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->vals, col->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->linkpos, col->size, newsize) );

      /* realloc colex */
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->vals, col->size, newsize) );

      col->size = newsize;
   }
   assert(num <= col->size);

   return SCIP_OKAY;
}

/** ensures, that cols array can store at least num entries */
static
SCIP_RETCODE ensureColexsSize(
   SCIP_LPEX*            lp,                 /**< current LP data */
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

/* forward declaration for colAddCoef() */
static
SCIP_RETCODE rowexAddCoef(
   SCIP_ROWEX*           row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_COLEX*           col,                /**< LP column */
   SCIP_Rational*        val,                /**< value of coefficient */
   int                   linkpos             /**< position of row in the column's row array, or -1 */
   );

/** adds a previously non existing coefficient to an LP column */
static
SCIP_RETCODE colexAddCoef(
   SCIP_COLEX*           col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_Rational*        val,                /**< value of coefficient */
   int                   linkpos             /**< position of column in the row's col array, or -1 */
   )
{
   int pos;

   assert(blkmem != NULL);
   assert(col != NULL);
   assert(col->nlprows <= col->len);
   assert(col->var != NULL);
   assert(!RisZero(val));
   /*assert(colSearchCoef(col, row) == -1);*/ /* this assert would lead to slight differences in the solution process */

   SCIP_CALL( colexEnsureSize(col, blkmem, set, col->len+1) );
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
         colexMoveCoef(col, col->nlprows, pos);
         pos = col->nlprows;
      }
      col->nlprows++;
   }

   /* insert the row at the correct position and update the links */
   col->rows[pos] = row;

   if( col->vals[pos] != NULL )
      Rset(col->vals[pos], val);
   else
      col->vals[pos] = Rcopy(blkmem, val);

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
         SCIP_CALL( rowexAddCoef(row, blkmem, set, eventqueue, lp, col, val, pos) );
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
         rowexSwapCoefs(row, set->buffer, linkpos, row->nlpcols-1);

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

   SCIPsetDebugMsg(set, "added coefficient %g * <%s> at position %d (%d/%d) to column <%s> (nunlinked=%d)\n",
      RgetRealApprox(val), row->fprow->name, pos, col->nlprows, col->len,
      SCIPvarGetName(col->var), col->nunlinked);

   return SCIP_OKAY;
}

/** deletes coefficient at given position from column */
static
SCIP_RETCODE colexDelCoefPos(
   SCIP_COLEX*           col,                /**< column to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEX*            lpex,               /**< current LP data */
   int                   pos                 /**< position in column vector to delete */
   )
{
   SCIP_ROWEX* row;

   assert(lpex != NULL);
   assert(col != NULL);
   assert(col->var != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < col->len);
   assert(col->rows[pos] != NULL);
   assert(col->linkpos[pos] == -1 || col->rows[pos]->cols[col->linkpos[pos]] == col);
   assert((pos < col->nlprows) == (col->linkpos[pos] >= 0 && col->rows[pos]->lppos >= 0));

   row = col->rows[pos];
   assert((row->lppos >= 0) == (pos < col->nlprows));

   /*SCIPsetDebugMsg(set, "deleting coefficient %g * <%s> at position %d from column <%s>\n",
     col->vals[pos], row->name, pos, SCIPvarGetName(col->var));*/

   if( col->linkpos[pos] == -1 )
      col->nunlinked--;

   /* if row is a linked LP row, move last linked LP coefficient to position of empty slot (deleted coefficient) */
   if( pos < col->nlprows )
   {
      colexMoveCoef(col, col->nlprows-1, pos);
      col->nlprows--;
      pos = col->nlprows;
   }

   /* move last coefficient to position of empty slot */
   colexMoveCoef(col, col->len-1, pos);
   col->len--;

   coefChangedExact(row, col, lpex);

   return SCIP_OKAY;
}

/** changes a coefficient at given position of an LP column */
static
SCIP_RETCODE colexChgCoefPos(
   SCIP_COLEX*           col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEX*            lp,                 /**< current LP data */
   int                   pos,                /**< position in column vector to change */
   SCIP_Rational*        val                 /**< value of coefficient */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(0 <= pos && pos < col->len);
   assert(col->rows[pos] != NULL);
   assert(col->linkpos[pos] == -1 || col->rows[pos]->cols[col->linkpos[pos]] == col);

   /*debugMsg(scip, "changing coefficient %g * <%s> at position %d of column <%s> to %g\n",
     col->vals[pos], col->rows[pos]->name, pos, SCIPvarGetName(col->var), val);*/

   if( RisZero(val) )
   {
      /* delete existing coefficient */
      SCIP_CALL( colexDelCoefPos(col, set, lp, pos) );
   }
   else if( !RisEqual(col->vals[pos], val) )
   {
      /* change existing coefficient */
      Rset(col->vals[pos], val);
      coefChangedExact(col->rows[pos], col, lp);
   }

   return SCIP_OKAY;
}


/* forward declaration for colAddCoef() */
static
SCIP_RETCODE rowexAddCoef(
   SCIP_ROWEX*           row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_COLEX*           col,                /**< LP column */
   SCIP_Rational*        val,                /**< value of coefficient */
   int                   linkpos             /**< position of row in the column's row array, or -1 */
   )
{
   int pos;

   assert(row != NULL);
   assert(row->nlpcols <= row->len);
   assert(blkmem != NULL);
   assert(col != NULL);
   assert(col->var != NULL);
   assert(!RisZero(val));

   if( row->nlocks > 0 )
   {
      SCIPerrorMessage("cannot add a coefficient to the locked unmodifiable row <%s>\n", row->name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIProwexEnsureSize(row, blkmem, set, row->len+1) );
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
         rowexMoveCoef(row, row->nlpcols, pos);
         pos = row->nlpcols;
      }
      row->nlpcols++;
   }

   /* insert the column at the correct position and update the links */
   row->cols[pos] = col;
   row->cols_index[pos] = col->index;
   if( row->vals[pos] == NULL )
      row->vals[pos] = Rcopy(blkmem, val);
   else
      Rset(row->vals[pos], val);

   row->linkpos[pos] = linkpos;
   row->integral = row->integral && SCIPcolIsIntegral(col->fpcol) && RisIntegral(val);
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
         SCIP_CALL( colexAddCoef(col, blkmem, set, eventqueue, lp, row, val, pos) );
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
         colexSwapCoefs(col, set->buffer, linkpos, col->nlprows-1);

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

   /* update row norm */
   //rowAddNorms(row, set, col, val, TRUE);

   coefChangedExact(row, col, lp);

   SCIPsetDebugMsg(set, "added coefficient %g * <%s> at position %d (%d/%d) to row <%s> (nunlinked=%d)\n",
      RgetRealApprox(val), SCIPvarGetName(col->var), pos, row->nlpcols, row->len, row->name, row->nunlinked);

   /* issue row coefficient changed event */
   //SCIP_CALL( rowEventCoefChanged(row, blkmem, set, eventqueue, col, 0.0, val) );

   return SCIP_OKAY;
}

/** deletes coefficient at given position from row */
static
SCIP_RETCODE rowexDelCoefPos(
   SCIP_ROWEX*           row,                /**< row to be changed */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   int                   pos                 /**< position in row vector to delete */
   )
{
   SCIP_COLEX* col;

   assert(row != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < row->len);
   assert(row->cols[pos] != NULL);
   assert((pos < row->nlpcols) == (row->linkpos[pos] >= 0 && row->cols[pos]->lppos >= 0));

   col = row->cols[pos];

   assert((pos < row->nlpcols) == (col->lppos >= 0 && row->linkpos[pos] >= 0));

   /*SCIPsetDebugMsg(set, "deleting coefficient %g * <%s> at position %d from row <%s>\n",
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
      rowexMoveCoef(row, row->nlpcols-1, pos);
      assert(!row->lpcolssorted);
      row->nlpcols--;
      pos = row->nlpcols;
   }

   /* move last coefficient to position of empty slot */
   rowexMoveCoef(row, row->len-1, pos);
   row->len--;

   /* update norms */
   /* todo: exip do we need this? */
   //rowexDelNorms(row, set, col, val, FALSE, TRUE, TRUE);

   coefChangedExact(row, col, lp);

   /* issue row coefficient changed event */
   /*  todo exip: do we need this? */
   //   SCIP_CALL( rowexEventCoefChanged(rowex, blkmem, set, eventqueue, colex, valex, NULL) );

   return SCIP_OKAY;
}

/** changes a coefficient at given position of an LP row */
static
SCIP_RETCODE rowexChgCoefPos(
   SCIP_ROWEX*           row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   int                   pos,                /**< position in row vector to change */
   SCIP_Rational*        val                 /**< value of coefficient */
   )
{
   SCIP_COLEX* col;

   assert(row != NULL);
   assert(lp != NULL);

   col = row->cols[pos];

   assert(col != NULL);
   assert(0 <= pos && pos < row->len);

   /*SCIPsetDebugMsg(set, "changing coefficient %g * <%s> at position %d of row <%s> to %g\n",
     row->vals[pos], SCIPvarGetName(row->cols[pos]->var), pos, row->name, val);*/

   if( row->nlocks > 0 )
   {
      SCIPerrorMessage("cannot change a coefficient of the locked unmodifiable row <%s>\n", row->name);
      return SCIP_INVALIDDATA;
   }

   assert(col != NULL);

   if( RisZero(val) )
   {
      /* delete existing coefficient */
      SCIP_CALL( rowexDelCoefPos(row, blkmem, set, eventqueue, lp, pos) );
   }
   else if( !RisEqual(row->vals[pos], val) )
   {
      /* change existing coefficient */
      //rowDelNorms(row, set, col, row->vals[pos], FALSE, FALSE, TRUE);
      Rset(row->vals[pos], val);
      row->integral = row->integral && SCIPcolIsIntegral(col->fpcol) && RisIntegral(val);
      //rowAddNorms(row, set, col, row->vals[pos], TRUE);
      coefChangedExact(row, col, lp);

      /* issue row coefficient changed event */
      //SCIP_CALL( rowEventCoefChanged(row, blkmem, set, eventqueue, col, oldval, val) );
   }

   return SCIP_OKAY;
}

/*
 * Column methods
 */

/** creates an LP column */
extern
SCIP_RETCODE SCIPcolexCreate(
   SCIP_COLEX**          col,                /**< pointer to column data */
   SCIP_COL*             fpcol,              /**< the corresponding fp col */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             var,                /**< variable, this column represents */
   int                   len,                /**< number of nonzeros in the column */
   SCIP_ROWEX**          rows,               /**< array with rows of column entries */
   SCIP_Rational**       vals,               /**< array with coefficients of column entries */
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
      SCIP_ALLOC( RcopyArray(blkmem, &(*col)->vals, vals, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*col)->linkpos, len) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*col)->rows, rows, len) );

      for( i = 0; i < len; ++i )
      {
         assert(!RisZero(vals[i]));
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
   (*col)->obj = Rcopy(blkmem, SCIPvarGetObjExact(var));
   (*col)->lb = Rcopy(blkmem, SCIPvarGetLbLocalExact(var));
   (*col)->ub = Rcopy(blkmem, SCIPvarGetUbLocalExact(var));
   //(*col)->unchangedobj = Rcopy(blkmem, SCIPvarGetUnchangedObjExact(var));
   //(*col)->lazylb = Rcopy(blkmem, SCIPvarGetLbLazyExact(var));
   //(*col)->lazyub = Rcopy(blkmem, SCIPvarGetUbLazyExact(var));
   (*col)->index = (*col)->fpcol->index;
   (*col)->flushedobj = Rcreate(blkmem);
   (*col)->flushedlb = Rcreate(blkmem);
   (*col)->flushedub = Rcreate(blkmem);
   (*col)->primsol = Rcreate(blkmem);
   (*col)->redcost = RcreateString(blkmem, "inf");
   (*col)->farkascoef = RcreateString(blkmem, "inf");
   (*col)->minprimsol = Rcopy(blkmem, (*col)->ub);
   (*col)->maxprimsol = Rcopy(blkmem, (*col)->lb);
   (*col)->removable = FALSE;
   //(*col)->sbdown = Rcreate(blkmem);
   //(*col)->sbup = Rcreate(blkmem);
   //(*col)->sbsolval = Rcreate(blkmem);
   //(*col)->sblpobjval = Rcreate(blkmem);
   (*col)->integral = SCIPvarIsIntegral(var);


   (*col)->size = len;
   (*col)->len = len;
   (*col)->nlprows = 0;
   (*col)->nunlinked = len;
   (*col)->lppos = -1;
   (*col)->lpipos = -1;
   

   return SCIP_OKAY;
}

/** frees an LP column */
extern
SCIP_RETCODE SCIPcolexFree(
   SCIP_COLEX**          col,                /**< pointer to LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp                  /**< current LP data */
   )
{
   assert(blkmem != NULL);
   assert(col != NULL);
   assert(*col != NULL);
   assert((*col)->fpcol != NULL);

   if( (*col)->size > 0 )
      RdeleteArray(blkmem, &(*col)->vals, (*col)->size);
   else
      assert((*col)->vals == NULL);

   Rdelete(blkmem, &(*col)->obj);
   Rdelete(blkmem, &(*col)->lb);
   Rdelete(blkmem, &(*col)->ub);
   //Rdelete(blkmem, &(*col)->unchangedobj);
   Rdelete(blkmem, &(*col)->flushedobj);
   Rdelete(blkmem, &(*col)->flushedlb);
   Rdelete(blkmem, &(*col)->flushedub);
   Rdelete(blkmem, &(*col)->primsol);
   Rdelete(blkmem, &(*col)->redcost);
   Rdelete(blkmem, &(*col)->farkascoef);
   Rdelete(blkmem, &(*col)->minprimsol);
   Rdelete(blkmem, &(*col)->maxprimsol);

   BMSfreeBlockMemory(blkmem, col);

   return SCIP_OKAY;
}

/** output column to file stream */
extern
void SCIPcolexPrint(
   SCIP_COLEX*           col,                /**< LP column */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int r;
   char buf[SCIP_MAXSTRLEN];

   assert(col != NULL);
   assert(col->fpcol != NULL);
   assert(col->fpcol->var != NULL);

   RtoString(col->obj, buf);
   SCIPmessageFPrintInfo(messagehdlr, file, "(obj: %s) ", buf);
   RtoString(col->lb, buf);
   SCIPmessageFPrintInfo(messagehdlr, file, "[%s, ", buf);
   RtoString(col->ub, buf);
   SCIPmessageFPrintInfo(messagehdlr, file, ",%s], ", buf);

   /* print coefficients */
   if( col->len == 0 )
      SCIPmessageFPrintInfo(messagehdlr, file, "<empty>");
   for( r = 0; r < col->len; ++r )
   {
      assert(col->rows[r] != NULL);
      assert(col->rows[r]->fprow->name != NULL);

      RtoString(col->vals[r], buf);
      if( RisPositive(col->vals[r]) )
         SCIPmessageFPrintInfo(messagehdlr, file, "+%s<%s> ", buf, col->rows[r]->fprow->name);
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "%s<%s> ", buf, col->rows[r]->fprow->name);
   }
   SCIPmessageFPrintInfo(messagehdlr, file, "\n");
}

/** adds a previously non existing coefficient to an LP column */
extern
SCIP_RETCODE SCIPcolexAddCoef(
   SCIP_COLEX*           col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_Rational*        val                 /**< value of coefficient */
   )
{
   assert(lp != NULL);

   SCIP_CALL( colexAddCoef(col, blkmem, set, eventqueue, lp, row, val, -1) );

   return SCIP_OKAY;
}

/** deletes coefficient from column */
extern
SCIP_RETCODE SCIPcolexDelCoef(
   SCIP_COLEX*           col,                /**< column to be changed */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_ROWEX*           row                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(col != NULL);
   assert(col->var != NULL);
   assert(lp != NULL);
   assert(row != NULL);

   /* search the position of the row in the column's row vector */
   pos = colexSearchCoef(col, row);
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
      assert(RisEqual(row->vals[col->linkpos[pos]], col->vals[pos]));
      SCIP_CALL( rowexDelCoefPos(row, blkmem, set, eventqueue, lp, col->linkpos[pos]) );
   }

   /* delete the row from the column's row vector */
   SCIP_CALL( colexDelCoefPos(col, set, lp, pos) );

   //checkLinks(lp);

   return SCIP_OKAY;
}

/** changes or adds a coefficient to an LP column */
extern
SCIP_RETCODE SCIPcolexChgCoef(
   SCIP_COLEX*           col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_Rational*        val                 /**< value of coefficient */
   )
{
   int pos;

   assert(col != NULL);
   assert(lp != NULL);
   assert(!lp->fplp->diving);
   assert(row != NULL);

   /* search the position of the row in the column's row vector */
   pos = colexSearchCoef(col, row);

   /* check, if row already exists in the column's row vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( colexAddCoef(col, blkmem, set, eventqueue, lp, row, val, -1) );
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
         assert(RisEqual(row->vals[col->linkpos[pos]], col->vals[pos]));
         SCIP_CALL( rowexChgCoefPos(row, blkmem, set, eventqueue, lp, col->linkpos[pos], val) );
      }

      /* change the coefficient in the column */
      SCIP_CALL( colexChgCoefPos(col, set, lp, pos, val) );
   }

   //   checkLinks(lp);

   return SCIP_OKAY;
}

/** increases value of an existing or nonexisting coefficient in an LP column */
extern
SCIP_RETCODE SCIPcolexIncCoef(
   SCIP_COLEX*           col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_Rational*        incval              /**< value to add to the coefficient */
   )
{
   int pos;

   assert(col != NULL);
   assert(lp != NULL);
   assert(!lp->fplp->diving);
   assert(row != NULL);

   if( RisZero(incval) )
      return SCIP_OKAY;

   /* search the position of the row in the column's row vector */
   pos = colexSearchCoef(col, row);

   /* check, if row already exists in the column's row vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( colexAddCoef(col, blkmem, set, eventqueue, lp, row, incval, -1) );
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
         assert(RisEqual(row->vals[col->linkpos[pos]], col->vals[pos]));

         Radd(incval, incval, col->vals[pos]);
         SCIP_CALL( rowexChgCoefPos(row, blkmem, set, eventqueue, lp, col->linkpos[pos], incval) );
      }

      /* change the coefficient in the column */
      SCIP_CALL( colexChgCoefPos(col, set, lp, pos, incval) );
   }

   //   checkLinks(lp);

   return SCIP_OKAY;
}

#if 0
/** changes objective value of column */
extern
SCIP_RETCODE SCIPcolexChgObj(
   SCIP_COLEX*           col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_Rational*        newobj              /**< new objective value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(lp != NULL);

   SCIPsetDebugMsg(set, "changing objective value of column <%s> from %f to %f\n", SCIPvarGetName(col->var), col->obj, newobj);

   /* only add actual changes */
   if( !SCIPsetIsEQ(set, col->obj, newobj) )
   {
      /* only variables with a real position in the LPI can be inserted */
      if( col->lpipos >= 0 )
      {
         /* insert column in the chgcols list (if not already there) */
         SCIP_CALL( insertColChgcols(col, set, lp) );

         /* mark objective value change in the column */
         col->objchanged = TRUE;

         assert(lp->nchgcols > 0);
      }
      /* in any case, when the sign of the objective (and thereby the best bound) changes, the variable has to enter the
       * LP and the LP has to be flushed
       */
      else if( (col->obj < 0.0 && newobj >= 0.0 && SCIPsetIsZero(set, col->ub))
         || (col->obj >= 0.0 && newobj < 0.0 && SCIPsetIsZero(set, col->lb)) )
      {
         /* mark the LP unflushed */
         lp->flushed = FALSE;
      }
   }

   /* store new objective function value */
   col->obj = newobj;

   /* update original objective value, as long as we are not in diving or probing and changed objective values */
   if( !lp->divingobjchg )
   {

      SCIP_Real oldobj = col->unchangedobj;

      assert(SCIPsetIsEQ(set, newobj, SCIPvarGetUnchangedObj(col->var)));
      col->unchangedobj = newobj;

      /* update the objective function vector norms */
      lpUpdateObjNorms(lp, set, oldobj, newobj);
   }

   return SCIP_OKAY;
}

/** changes lower bound of column */
extern
SCIP_RETCODE SCIPcolexChgLb(
   SCIP_COLEX*           col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_Rational*        newlb               /**< new lower bound value */
   );

/** changes upper bound of column */
extern
SCIP_RETCODE SCIPcolexChgUb(
   SCIP_COLEX*           col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_Rational*        newub               /**< new upper bound value */
   );
#endif

/** creates and captures an LP row */
SCIP_RETCODE SCIProwCreateExact(
   SCIP_ROWEX**          row,                /**< pointer to LP row data */
   SCIP_ROW*             fprow,              /**< corresponding fp row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEX*            lp,                 /**< current LP data */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COLEX**          cols,               /**< array with columns of row entries */
   SCIP_Rational**       vals,               /**< array with coefficients of row entries */
   SCIP_Rational*        lhs,                /**< left hand side of row */
   SCIP_Rational*        rhs,                /**< right hand side of row */
   SCIP_ROWORIGINTYPE    origintype,         /**< type of origin of row */
   void*                 origin              /**< pointer to constraint handler or separator who created the row (NULL if unkown) */
   )
{
   assert(row != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(len >= 0);
   assert(len == 0 || (cols != NULL && vals != NULL));
   assert(SCIProwGetNNonz(fprow) == len || len == 0);
   /* note, that the assert tries to avoid numerical troubles in the LP solver.
    * in case, for example, lhs > rhs but they are equal with tolerances, one could pass lhs=rhs=lhs+rhs/2 to
    * SCIProwCreate() (see cons_linear.c: detectRedundantConstraints())
    */
   assert(RisLE(lhs, rhs));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, row) );

   (*row)->integral = TRUE;
   (*row)->fprow = fprow;

   if( len > 0 )
   {
      SCIP_VAR* var;
      int i;

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*row)->cols, cols, len) );
      //SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*row)->vals, vals, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*row)->cols_index, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*row)->linkpos, len) );
      SCIP_ALLOC( RcopyArray(blkmem, &(*row)->vals, vals, len) );

      for( i = 0; i < len; ++i )
      {
         assert(cols[i] != NULL);
         assert(!RisZero(vals[i]));

         var = cols[i]->var;
         (*row)->cols_index[i] = cols[i]->index;
         (*row)->linkpos[i] = -1;

         if( RisIntegral((*row)->vals[i]) )
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
      (*row)->linkpos = NULL;
      (*row)->cols_index = NULL;
   }

   (*row)->lhs = Rcopy(blkmem, lhs);
   (*row)->rhs = Rcopy(blkmem, rhs);
   (*row)->flushedlhs = RcreateString(blkmem, "-inf");
   (*row)->flushedrhs = RcreateString(blkmem, "inf");
   (*row)->objprod = RcreateString(blkmem, "0");
   (*row)->maxval = RcreateString(blkmem, "0");
   (*row)->minval = RcreateString(blkmem, "inf");
   (*row)->dualsol = RcreateString(blkmem, "0");
   (*row)->activity = RcreateString(blkmem, "inf");
   (*row)->dualfarkas = RcreateString(blkmem, "0");
   (*row)->pseudoactivity = RcreateString(blkmem, "inf");
   (*row)->minactivity = RcreateString(blkmem, "inf");
   (*row)->maxactivity = RcreateString(blkmem, "inf");
   (*row)->constant = RcreateString(blkmem, "0");

   (*row)->origin = origin;
   (*row)->eventfilter = NULL;
   (*row)->index = stat->nrowidx;
   SCIPstatIncrement(stat, set, nrowidx);
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
   (*row)->numintcols = -1;
   (*row)->validactivitylp = -1;
   (*row)->age = 0;
   (*row)->rank = 0;
   (*row)->delaysort = FALSE;
   (*row)->lpcolssorted = TRUE;
   (*row)->nonlpcolssorted = (len <= 1);
   (*row)->delaysort = FALSE;
   (*row)->nlocks = 0;
   (*row)->removable = FALSE;

   //SCIPhashtableInsert(lp->exrowhash, (void*) (*row));

   /* capture the row */
   SCIProwexCapture(*row);

   return SCIP_OKAY;
} /*lint !e715*/

/*
 * lp methods
 */

/** updates link data after addition of row */
static
void rowexUpdateAddLP(
   SCIP_ROWEX*           row,                /**< LP row */
   BMS_BUFMEM*           buffer              /**< buffer for rationals */
   )
{
   SCIP_COLEX* col;
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
         colexSwapCoefs(col, buffer, pos, col->nlprows - 1);

         /* if no swap was necessary, mark lprows to be unsorted */
         if( pos == col->nlprows-1 )
            col->lprowssorted = FALSE;
      }
   }
}

/** updates link data after addition of column */
static
void colexUpdateAddLP(
   SCIP_COLEX*           col,                /**< LP column */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_ROWEX* row;
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
         rowexSwapCoefs(row, set->buffer, pos, row->nlpcols-1);
         assert(row->cols[row->nlpcols-1] == col);

         /* if no swap was necessary, mark lpcols to be unsorted */
         if( pos == row->nlpcols-1 )
            row->lpcolssorted = FALSE;

         /* update norms */
         //rowAddNorms(row, set, col, row->vals[row->nlpcols-1], FALSE);
      }
   }
}

/** checks that lp and fplp are properly synced */
extern
SCIP_Bool SCIPlpexIsSynced(
   SCIP_LPEX*            lp,                 /**< pointer to LP data object */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     msg                 /**< message handler */
   )
{
   assert(lp != NULL);
   assert(msg != NULL);

   return lpexInSync(lp, set, msg);
}

/** creates empty LP data object */
extern
SCIP_RETCODE SCIPlpexCreate(
   SCIP_LPEX**           lp,                 /**< pointer to LP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_LP*              fplp,               /**< the floating point LP */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name                /**< problem name */
   )
{
   SCIP_Bool success;

   assert(lp != NULL);
   assert(fplp != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(name != NULL);

   SCIP_ALLOC( BMSallocMemory(lp) );

   /* open LP Solver interface */
   SCIP_CALL( SCIPlpiexCreate(&(*lp)->lpiex, messagehdlr, name, SCIP_OBJSEN_MINIMIZE) );

   (*lp)->fplp = fplp;
   fplp->lpex = *lp;

   (*lp)->lpicols = NULL;
   (*lp)->lpirows = NULL;
   (*lp)->chgcols = NULL;
   (*lp)->chgrows = NULL;
   (*lp)->cols = NULL;
   (*lp)->soldirection = NULL;
   (*lp)->lazycols = NULL;
   (*lp)->rows = NULL;
   (*lp)->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   (*lp)->validsoldirsol = NULL;
   (*lp)->flushdeletedcols = FALSE;
   (*lp)->flushaddedcols = FALSE;
   (*lp)->flushdeletedrows = FALSE;
   (*lp)->flushaddedrows = FALSE;
   (*lp)->updateintegrality = TRUE;
   (*lp)->flushed = TRUE;
   (*lp)->solved = TRUE;
   (*lp)->primalfeasible = TRUE;
   (*lp)->primalchecked = TRUE;
   (*lp)->dualfeasible = TRUE;
   (*lp)->dualchecked = TRUE;
   (*lp)->solisbasic = FALSE;
   (*lp)->rootlpisrelax = TRUE;
   (*lp)->isrelax = TRUE;
   (*lp)->installing = FALSE;
   (*lp)->resolvelperror = FALSE;
   (*lp)->adjustlpval = FALSE;
   (*lp)->lpiscaling = set->lp_scaling;
   (*lp)->lpisolutionpolishing = (set->lp_solutionpolishing > 0);
   (*lp)->lpirefactorinterval = set->lp_refactorinterval;
   (*lp)->lpiitlim = INT_MAX;
   (*lp)->lpipricing = SCIP_PRICING_AUTO;
   (*lp)->lastlpalgo = SCIP_LPALGO_DUALSIMPLEX;
   (*lp)->lpitiming = (int) set->time_clocktype;
   (*lp)->lpirandomseed = set->random_randomseed;
   (*lp)->storedsolvals = NULL;

   (*lp)->lpicolssize = 0;
   (*lp)->nlpicols = 0;
   (*lp)->lpirowssize = 0;
   (*lp)->nlpirows = 0;
   (*lp)->lpifirstchgcol = 0;
   (*lp)->lpifirstchgrow = 0;
   (*lp)->colssize = 0;
   (*lp)->soldirectionsize = 0;
   (*lp)->ncols = 0;
   //(*lp)->lazycolssize = 0;
   //(*lp)->nlazycols = 0;
   (*lp)->rowssize = 0;
   (*lp)->nrows = 0;
   (*lp)->chgcolssize = 0;
   (*lp)->nchgcols = 0;
   (*lp)->chgrowssize = 0;
   (*lp)->nchgrows = 0;
   (*lp)->firstnewcol = 0;
   (*lp)->firstnewrow = 0;
   (*lp)->nremovablecols = 0;
   (*lp)->nremovablerows = 0;
   (*lp)->looseobjvalinf = 0;
   (*lp)->pseudoobjvalinf = 0;
   (*lp)->glbpseudoobjvalinf = 0;
   (*lp)->lpobjval = Rcreate(blkmem);
   (*lp)->pseudoobjval = Rcreate(blkmem);
   (*lp)->glbpseudoobjval = Rcreate(blkmem);
   (*lp)->looseobjval = Rcreate(blkmem);
   (*lp)->relglbpseudoobjval = Rcreate(blkmem);
   (*lp)->rellooseobjval = Rcreate(blkmem);
   (*lp)->relpseudoobjval = Rcreate(blkmem);
   (*lp)->rootlpobjval = Rcreate(blkmem);
   (*lp)->rootlooseobjval = Rcreate(blkmem);
   (*lp)->cutoffbound = Rcreate(blkmem);
   (*lp)->lpiobjlim = Rcreate(blkmem);

   //(*lp)->validsollp = stat->lpcount; /* the initial (empty) SCIP_LP is solved with primal and dual solution of zero */
   //(*lp)->validfarkaslp = -1;
   //(*lp)->validsoldirlp = -1;


   SCIPhashtableCreate(&(*lp)->exrowhash, blkmem, (set->misc_usesmalltables ? SCIP_HASHSIZE_CUTPOOLS_SMALL : SCIP_HASHSIZE_CUTPOOLS),
                     hashKeyExRow, hashEqExRow, hashValExRow, (void*) set );

   /** todo: exip: set the right defaults in lp solver */
   return SCIP_OKAY;
}

/** frees LP data object */
SCIP_RETCODE SCIPlpexFree(
   SCIP_LPEX**           lp,                 /**< pointer to LP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   )
{
   int i;

   assert(lp != NULL);
   assert(*lp != NULL);

   SCIP_CALL( SCIPlpexClear(*lp, blkmem, set, eventqueue, eventfilter) );

   SCIPhashtableFree(&(*lp)->exrowhash);
   //freeDiveChgSideArrays(*lp);

   /* release LPI rows */
   for( i = 0; i < (*lp)->nlpirows; ++i )
   {
      SCIP_CALL( SCIProwexRelease(&(*lp)->lpirows[i], blkmem, set, *lp) );
   }

   if( (*lp)->lpiex != NULL )
   {
      SCIP_CALL( SCIPlpiexFree(&(*lp)->lpiex) );
   }

   Rdelete(blkmem, &(*lp)->lpobjval);
   Rdelete(blkmem, &(*lp)->pseudoobjval);
   Rdelete(blkmem, &(*lp)->glbpseudoobjval);
   Rdelete(blkmem, &(*lp)->looseobjval);
   Rdelete(blkmem, &(*lp)->rellooseobjval);
   Rdelete(blkmem, &(*lp)->relglbpseudoobjval);
   Rdelete(blkmem, &(*lp)->relpseudoobjval);
   Rdelete(blkmem, &(*lp)->rootlpobjval);
   Rdelete(blkmem, &(*lp)->rootlooseobjval);
   Rdelete(blkmem, &(*lp)->cutoffbound);
   Rdelete(blkmem, &(*lp)->lpiobjlim);

   BMSfreeMemoryNull(&(*lp)->storedsolvals);
   BMSfreeMemoryArrayNull(&(*lp)->lpicols);
   BMSfreeMemoryArrayNull(&(*lp)->lpirows);
   BMSfreeMemoryArrayNull(&(*lp)->chgcols);
   BMSfreeMemoryArrayNull(&(*lp)->chgrows);
   BMSfreeMemoryArrayNull(&(*lp)->lazycols);
   BMSfreeMemoryArrayNull(&(*lp)->cols);
   BMSfreeMemoryArrayNull(&(*lp)->rows);
   BMSfreeMemoryArrayNull(&(*lp)->soldirection);
   BMSfreeMemory(lp);

   return SCIP_OKAY;
}




/** adds a column to the LP and captures the variable */
SCIP_RETCODE SCIPlpexAddCol(
   SCIP_LPEX*            lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COLEX*           col,                /**< LP column */
   int                   depth               /**< depth in the tree where the column addition is performed */
   )
{
   assert(lp != NULL);
   assert(!lp->fplp->diving);
   assert(col != NULL);
   assert(col->len == 0 || col->rows != NULL);
   assert(col->lppos == -1);
   assert(col->var != NULL);
   assert(SCIPvarGetStatusExact(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col->fpcol);
   assert(SCIPvarIsIntegral(col->var) == col->integral);

   SCIPsetDebugMsg(set, "adding column <%s> to exact LP (%d rows, %d cols)\n", SCIPvarGetName(col->var), lp->nrows, lp->ncols);
#ifdef SCIP_DEBUG
   /* {
      int i;
      SCIPsetDebugMsgPrint(set, "  (obj: %g) [%g,%g]", col->obj, col->lb, col->ub);
      for( i = 0; i < col->len; ++i )
         SCIPsetDebugMsgPrint(set, " %+g<%s>", col->vals[i], col->rows[i]->name);
      SCIPsetDebugMsgPrint(set, "\n");
   } */
#endif

   SCIP_CALL( ensureColexsSize(lp, set, lp->ncols+1) );
   lp->cols[lp->ncols] = col;
   col->lppos = lp->ncols;
   lp->ncols++;
   if( col->removable )
      lp->nremovablecols++;

   /** todo: exip: do we need this? if( !RisNegInfinity(col->lazylb) || !RisInfinity(col->lazyub) )
   {
      SCIP_CALL( ensureLazycolesSize(lp, set, lp->nlazycols+1) );
      lp->lazycols[lp->nlazycols] = col;
      lp->nlazycols++;
   } */

   /* mark the current LP unflushed */
   lp->flushed = FALSE;

   /* update column arrays of all linked rows */
   colexUpdateAddLP(col, set);

   /* update the objective function vector norms */
   //lpUpdateObjNorms(lp, set, 0.0, col->unchangedobj);

   //checkLinks(lp);

   return SCIP_OKAY;
}

/** adds a row to the LP and captures it */
SCIP_RETCODE SCIPlpexAddRow(
   SCIP_LPEX*            lpex,               /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_ROWEX*           rowex,              /**< LP row */
   int                   depth               /**< depth in the tree where the row addition is performed */
   )
{
   assert(lpex != NULL);
   assert(rowex != NULL);
   assert(rowex->len == 0 || rowex->cols != NULL);
   assert(rowex->lppos == -1);
   assert(rowex->fprow != NULL);

   /** todo: exip do we need capture and locks on exact rows? */
   SCIProwexCapture(rowex);
   //SCIProwLock(row);

   SCIPsetDebugMsg(set, "adding row <%s> to LP (%d rows, %d cols)\n", rowex->fprow->name, lpex->nrows, lpex->ncols);
#ifdef SCIP_DEBUG
   {
      int i;
      SCIPsetDebugMsgPrint(set, "  %g <=", rowex->fprow->lhs);
      for( i = 0; i < row->len; ++i )
         SCIPsetDebugMsgPrint(set, " %+g<%s>", rowex->fprow->vals[i], SCIPvarGetName(rowex->fprow->cols[i]->var));
      if( !SCIPsetIsZero(set, row->constant) )
         SCIPsetDebugMsgPrint(set, " %+g", rowex->fprow->constant);
      SCIPsetDebugMsgPrint(set, " <= %g\n", row->fprow->rhs);
   }
#endif

   SCIP_CALL( ensureRowexsSize(lpex, set, lpex->nrows+1) );
   lpex->rows[lpex->nrows] = rowex;
   rowex->lppos = lpex->nrows;
   lpex->nrows++;

   /* mark the current LP unflushed */
   lpex->flushed = FALSE;

   /* update row arrays of all linked columns */
   rowexUpdateAddLP(rowex, set->buffer);

   //rowCalcNorms(rowex, set);

   /* check, if row addition to LP events are tracked
    * if so, issue ROWADDEDLP event
    */
   //if( (eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWADDEDLP) != 0) )
   //{
   //   SCIP_EVENT* event;
   //   SCIP_CALL( SCIPeventCreateRowAddedLP(&event, blkmem, rowex) );
   //   SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
   //}
   return SCIP_OKAY;
}

/*
 * row mehods
 */

/** increases usage counter of LP row */
void SCIProwexCapture(
   SCIP_ROWEX*           row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->nuses >= 0);
   assert(row->nlocks <= (unsigned int)(row->nuses)); /*lint !e574*/

   SCIPdebugMessage("capture row <%s> with nuses=%d and nlocks=%u\n", row->name, row->nuses, row->nlocks);
   row->nuses++;
}

/** output column to file stream */
extern
void SCIProwexPrint(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int r;
   char buf[SCIP_MAXSTRLEN];

   assert(row != NULL);
   assert(row->fprow != NULL);

   SCIPmessageFPrintInfo(messagehdlr, file, "%s: ", row->fprow->name);
   RtoString(row->lhs, buf);
   SCIPmessageFPrintInfo(messagehdlr, file, "%s <= ", buf);

   /* print coefficients */
   if( row->len == 0 )
      SCIPmessageFPrintInfo(messagehdlr, file, "<empty>");
   for( r = 0; r < row->len; ++r )
   {
      RtoString(row->vals[r], buf);
      assert(SCIPvarGetName(row->cols[r]->var) != NULL);
      assert(SCIPvarGetStatus(row->cols[r]->var) == SCIP_VARSTATUS_COLUMN);
      if( RisPositive(row->vals[r]) )
         SCIPmessageFPrintInfo(messagehdlr, file, "+%s<%s> ", buf, SCIPvarGetName(row->cols[r]->var));
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "%s<%s> ", buf, SCIPvarGetName(row->cols[r]->var));
   }

   RtoString(row->rhs, buf);
   SCIPmessageFPrintInfo(messagehdlr, file, "<= %s, ", buf);
   SCIPmessageFPrintInfo(messagehdlr, file, "\n");
}

/** get the index of an exact row */
int SCIProwexGetIndex(
   SCIP_ROWEX*           row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->index;
}

/** returns TRUE iff row is member of current LP */
SCIP_Bool SCIProwexIsInLP(
   SCIP_ROWEX*           row                 /**< LP row */
   )
{
   assert(row != NULL);

   return (row->lppos >= 0);
}

/** returns true, if an exact row for this fprow was already created */
SCIP_Bool SCIProwHasExRow(
   SCIP_LPEX*            lpex,               /**< exact lp data structure */
   SCIP_ROW*             row                 /**< SCIP row */
   )
{
   assert(row != NULL);
   assert(lpex != NULL);
   assert(lpex->exrowhash != NULL);

   return (NULL != SCIPhashtableRetrieve(lpex->exrowhash, row));
}

/** returns exact row corresponding to fprow, if it exists. Otherwise returns NULL */
SCIP_ROWEX* SCIProwGetExRow(
   SCIP_LPEX*            lpex,               /**< exact lp data structure */
   SCIP_ROW*             row                 /**< SCIP row */
   )
{
   assert(row != NULL);
   assert(lpex != NULL);
   assert(lpex->exrowhash != NULL);

   return SCIPhashtableRetrieve(lpex->exrowhash, row);
}

/** adds a previously non existing coefficient to an LP row */
SCIP_RETCODE SCIProwexAddCoef(
   SCIP_ROWEX*           rowex,              /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_COLEX*           colex,              /**< LP column */
   SCIP_Rational*        val                 /**< value of coefficient */
   )
{
   SCIP_ROW* row;
   SCIP_COL* col;

   assert(rowex != NULL);
   assert(colex != NULL);
   assert(lp != NULL);

   row = rowex->fprow;
   col = colex->fpcol;

   assert(lp != NULL);
   assert(!lp->fplp->diving || row->lppos == -1);

   SCIP_CALL( rowexAddCoef(rowex, blkmem, set, eventqueue, lp, colex, val, -1) );

   //checkLinks(lp);

   return SCIP_OKAY;
}

/** deletes coefficient from row */
SCIP_RETCODE SCIProwexDelCoef(
   SCIP_ROWEX*           row,                /**< row to be changed */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_COLEX*           col                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(row != NULL);
   assert(!row->delaysort);
   assert(lp != NULL);
   assert(!lp->fplp->diving || row->lppos == -1);
   assert(col != NULL);
   assert(col->var != NULL);

   /* search the position of the column in the row's col vector */
   pos = rowexSearchCoef(row, col);
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
      assert(RisEqual(col->vals[row->linkpos[pos]], row->vals[pos]));
      SCIP_CALL( colexDelCoefPos(col, set, lp, row->linkpos[pos]) );
   }

   /* delete the column from the row's col vector */
   SCIP_CALL( rowexDelCoefPos(row, blkmem, set, eventqueue, lp, pos) );

   //checkLinks(lp);

   return SCIP_OKAY;
}

/** changes or adds a coefficient to an LP row */
SCIP_RETCODE SCIProwexChgCoef(
   SCIP_ROWEX*           row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_COLEX*           col,                /**< LP column */
   SCIP_Rational*        val                 /**< value of coefficient */
   )
{
   int pos;

   assert(row != NULL);
   assert(!row->delaysort);
   assert(lp != NULL);
   assert(!lp->fplp->diving || row->lppos == -1);
   assert(col != NULL);

   /* search the position of the column in the row's col vector */
   pos = rowexSearchCoef(row, col);

   /* check, if column already exists in the row's col vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( rowexAddCoef(row, blkmem, set, eventqueue, lp, col, val, -1) );
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
         assert(RisEqual(col->vals[row->linkpos[pos]], row->vals[pos]));
         SCIP_CALL( colexChgCoefPos(col, set, lp, row->linkpos[pos], val) );
      }

      /* change the coefficient in the row */
      SCIP_CALL( rowexChgCoefPos(row, blkmem, set, eventqueue, lp, pos, val) );
   }

   //checkLinks(lp);

   return SCIP_OKAY;
}

/** increases value of an existing or non-existing coefficient in an LP row */
SCIP_RETCODE SCIProwexIncCoef(
   SCIP_ROWEX*           row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_COLEX*           col,                /**< LP column */
   SCIP_Rational*        incval              /**< value to add to the coefficient */
   )
{
   int pos;

   assert(row != NULL);
   assert(lp != NULL);
   assert(!lp->fplp->diving || row->lppos == -1);
   assert(col != NULL);

   if( RisZero(incval) )
      return SCIP_OKAY;

   /* search the position of the column in the row's col vector */
   pos = rowexSearchCoef(row, col);

   /* check, if column already exists in the row's col vector */
   if( pos == -1 )
   {
      /* coefficient doesn't exist, or sorting is delayed: add coefficient to the end of the row's arrays */
      SCIP_CALL( rowexAddCoef(row, blkmem, set, eventqueue, lp, col, incval, -1) );
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
         assert(RisEqual(col->vals[row->linkpos[pos]], row->vals[pos]));
         Radd(incval, incval, row->vals[pos]);
         SCIP_CALL( colexChgCoefPos(col, set, lp, row->linkpos[pos], incval) );
      }

      /* change the coefficient in the row */
      SCIP_CALL( rowexChgCoefPos(row, blkmem, set, eventqueue, lp, pos, incval) );
   }

   //checkLinks(lp);

   /* invalid the activity */
   row->validactivitylp = -1;

   return SCIP_OKAY;
}

/** changes constant value of a row */
SCIP_RETCODE SCIProwexChgConstant(
   SCIP_ROWEX*           row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_Rational*        constant            /**< new constant value */
   )
{
   assert(row != NULL);
   assert(row->lhs <= row->rhs);
   assert(!RisAbsInfinity(constant));
   assert(stat != NULL);
   assert(lp != NULL);
   assert(!lp->fplp->diving || row->fprow->lppos == -1);

   if( !RisEqual(constant, row->constant) )
   {
      if( row->fprow->validpsactivitydomchg == stat->domchgcount )
      {
         assert(!RisInfinity(row->pseudoactivity));
         Radd(row->pseudoactivity, row->pseudoactivity, constant);
         Rdiff(row->pseudoactivity, row->pseudoactivity, row->constant);
      }
      if( row->fprow->validactivitybdsdomchg == stat->domchgcount )
      {
         assert(!RisInfinity(row->minactivity));
         assert(!RisInfinity(row->maxactivity));
         Radd(row->minactivity, row->minactivity, constant);
         Radd(row->maxactivity, row->maxactivity, constant);
         Rdiff(row->minactivity, row->minactivity, row->constant);
         Rdiff(row->maxactivity, row->maxactivity, row->constant);
      }
      /* todo exip do we need the changed event here? */
   }

   return SCIP_OKAY;
}

/** add constant value to a row */
SCIP_RETCODE SCIProwexAddConstant(
   SCIP_ROWEX*           row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_Rational*        addval              /**< constant value to add to the row */
   )
{
   SCIP_Rational* tmp;

   assert(row != NULL);
   assert(row->lhs <= row->rhs);
   assert(!RisAbsInfinity(addval));
   assert(stat != NULL);
   assert(lp != NULL);
   assert(!lp->fplp->diving || row->fprow->lppos == -1);

   if( !RisZero(addval) )
   {
      tmp = RcreateTemp(set->buffer);
      Radd(tmp, row->constant, addval);
      SCIP_CALL( SCIProwexChgConstant(row, blkmem, set, stat, eventqueue, lp, tmp) );

      RdeleteTemp(set->buffer, &tmp);
   }

   return SCIP_OKAY;
}

/** returns the feasibility of a row for the given solution */
void SCIProwexGetSolFeasibility(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Rational*        result              /**< result pointer */
   )
{
   SCIP_Real activity;
   SCIP_Rational* temp1 = RcreateTemp(set->buffer);
   SCIP_Rational* temp2 = RcreateTemp(set->buffer);

   assert(row != NULL);

   SCIProwexGetSolActivity(row, set, stat, sol, result);

   Rdiff(temp1, row->rhs, result);
   Rdiff(temp2, result, row->lhs);
   Rmin(result, temp1, temp2);

   RdeleteTemp(set->buffer, &temp1);
   RdeleteTemp(set->buffer, &temp2);
}

/** returns the activity of a row for a given solution */
void SCIProwexGetSolActivity(
   SCIP_ROWEX*           rowex,              /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Rational*        result              /**< result pointer */
   )
{
   /** todo: exip: rational solution might be necessary */
   SCIP_COLEX* colex;
   SCIP_Real inf;
   SCIP_Rational* solval = RcreateTemp(set->buffer);
   int i;

   assert(rowex != NULL);

   Rset(result, rowex->constant);
   for( i = 0; i < rowex->len; ++i )
   {
      colex = rowex->cols[i];

      assert(colex != NULL);

      assert((i < rowex->nlpcols) == (rowex->linkpos[i] >= 0
         && colex->lppos >= 0));

      RsetReal(solval, SCIPsolGetVal(sol, set, stat, colex->var));

      if( RisAbsInfinity(solval) ) /*lint !e777*/
      {
         if( RisNegInfinity(rowex->lhs) )
            RisPositive(rowex->vals[i]) ? Rset(solval, colex->lb) : Rset(solval, colex->ub);
         else if( RisInfinity(rowex->rhs) )
            RisPositive(rowex->vals[i]) ? Rset(solval, colex->ub) : Rset(solval, colex->lb);
         else
            Radd(solval, colex->lb, colex->ub);
            RmultReal(solval, solval, 0.5);
      }

      Rmult(solval, solval, rowex->vals[i]);
      Radd(result, result, solval);
   }

   RdeleteTemp(set->buffer, &solval);
}


/** decreases usage counter of LP row, and frees memory if necessary */
SCIP_RETCODE SCIProwexRelease(
   SCIP_ROWEX**          row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEX*            lp                  /**< current LP data */
   )
{
   assert(blkmem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row) != NULL);
   assert((*row)->nuses >= 1);
   assert((*row)->nlocks < (unsigned int)((*row)->nuses)); /*lint !e574*/

   SCIPsetDebugMsg(set, "release row <%s> with nuses=%d and nlocks=%u\n",
      (*row)->fprow->name, (*row)->nuses, (*row)->nlocks);
   (*row)->nuses--;
   if( (*row)->nuses == 0 )
   {
      SCIP_CALL( SCIProwexFree(row, blkmem, set, lp) );
   }

   *row = NULL;

   return SCIP_OKAY;
}

/** frees an LP row */
SCIP_RETCODE SCIProwexFree(
   SCIP_ROWEX**          row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEX*            lp                  /**< current LP data */
   )
{
   assert(blkmem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->nuses == 0);
   assert((*row)->lppos == -1);
   //assert((*row)->eventfilter != NULL);

   Rdelete(blkmem, &(*row)->constant);
   Rdelete(blkmem, &(*row)->lhs);
   Rdelete(blkmem, &(*row)->rhs);
   Rdelete(blkmem, &(*row)->flushedlhs);
   Rdelete(blkmem, &(*row)->flushedrhs);
   Rdelete(blkmem, &(*row)->objprod);
   Rdelete(blkmem, &(*row)->maxval);
   Rdelete(blkmem, &(*row)->minval);
   Rdelete(blkmem, &(*row)->dualsol);
   Rdelete(blkmem, &(*row)->activity);
   Rdelete(blkmem, &(*row)->dualfarkas);
   Rdelete(blkmem, &(*row)->pseudoactivity);
   Rdelete(blkmem, &(*row)->minactivity);
   Rdelete(blkmem, &(*row)->maxactivity);
   
   RdeleteArray(blkmem, &(*row)->vals, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->cols, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->cols_index, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->linkpos, (*row)->size);
   BMSfreeBlockMemory(blkmem, row);

   return SCIP_OKAY;
}

/** returns the feasibility of a row in the current LP solution: negative value means infeasibility */
void SCIProwexGetLPFeasibility(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_Rational*        result              /**< rational pointer to store the result */
   )
{
   SCIP_Rational* activity;
   SCIP_Rational* actrhs = RcreateTemp(set->buffer);
   SCIP_Rational* actlhs = RcreateTemp(set->buffer);

   assert(row != NULL);

   activity = SCIProwexGetLPActivity(row, set, stat, lp);

   Rdiff(actlhs, row->rhs, activity);
   Rdiff(actrhs, activity, row->lhs);
   Rmin(result, actrhs, actlhs);

   RdeleteTemp(set->buffer, &actrhs);
   RdeleteTemp(set->buffer, &actlhs);
}

/** returns the pseudo feasibility of a row in the current pseudo solution: negative value means infeasibility */
void SCIProwexGetPseudoFeasibility(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Rational*        result              /**< rational pointer to store the result */
   )
{
   SCIP_Rational* pseudoactivity;
   SCIP_Rational* actrhs = RcreateTemp(set->buffer);
   SCIP_Rational* actlhs = RcreateTemp(set->buffer);

   assert(row != NULL);

   pseudoactivity = SCIProwexGetPseudoActivity(row, set, stat);

   Rdiff(actlhs, row->rhs, pseudoactivity);
   Rdiff(actrhs, pseudoactivity, row->lhs);
   Rmin(result, actrhs, actlhs);

   RdeleteTemp(set->buffer, &actrhs);
   RdeleteTemp(set->buffer, &actlhs);
}

/** returns the activity of a row in the current LP solution */
SCIP_Rational* SCIProwexGetLPActivity(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEX*            lp                  /**< current LP data */
   )
{
   SCIP_Real inf;

   assert(row != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(row->fprow->validactivitylp <= stat->lpcount);
   assert(lp->fplp->validsollp == stat->lpcount);

   if( row->fprow->validactivitylp != stat->lpcount )
      SCIProwexRecalcLPActivity(row, set, stat);
   assert(row->fprow->validactivitylp == stat->lpcount);
   assert(row->fprow->activity < SCIP_INVALID);

   return row->activity;
}

/** returns the pseudo activity of a row in the current pseudo solution */
SCIP_Rational* SCIProwexGetPseudoActivity(
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(row != NULL);
   assert(stat != NULL);
   assert(row->fprow->validpsactivitydomchg <= stat->domchgcount);

   /* check, if pseudo activity has to be calculated */
   if( row->fprow->validpsactivitydomchg != stat->domchgcount )
      SCIProwexRecalcPseudoActivity(row, set, stat);
   assert(row->fprow->validpsactivitydomchg == stat->domchgcount);
   assert(row->fprow->pseudoactivity < SCIP_INVALID);

   return row->pseudoactivity;
}

/** recalculates the current activity of a row */
void SCIProwexRecalcLPActivity(
   SCIP_ROWEX*           rowex,              /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_COLEX* colex;
   SCIP_COL* col;
   SCIP_ROW* row;
   SCIP_Rational*  tmp = RcreateTemp(set->buffer);
   int c;

   assert(rowex != NULL);

   row = rowex->fprow;

   assert(row != NULL);
   assert(stat != NULL);


   Rset(rowex->activity, rowex->constant);
   for( c = 0; c < row->nlpcols; ++c )
   {
      colex = rowex->cols[c];
      col = row->cols[c];

      assert(col != NULL);
      assert(colex != NULL);
      assert(!RisInfinity(colex->primsol));
      assert(col->lppos >= 0);
      assert(row->linkpos[c] >= 0);
      Rmult(tmp, rowex->vals[c], colex->primsol);
      Radd(rowex->activity, rowex->activity, tmp);
   }

   if( row->nunlinked > 0 )
   {
      for( c = row->nlpcols; c < row->len; ++c )
      {
         col = row->cols[c];
         colex = rowex->cols[c];

         assert(col != NULL);
         assert(colex != NULL);
         assert(col->lppos >= 0 || col->primsol == 0.0);
         assert(col->lppos == -1 || row->linkpos[c] == -1);
         if( col->lppos >= 0 )
         {
            Rmult(tmp, rowex->vals[c], colex->primsol);
            Radd(rowex->activity, rowex->activity, tmp);
         }
      }
   }
#ifndef NDEBUG
   else
   {
      for( c = row->nlpcols; c < row->len; ++c )
      {
         col = row->cols[c];
         colex = rowex->cols[c];

         assert(col != NULL);
         assert(colex != NULL);
         assert(RisZero(colex->primsol));
         assert(col->lppos == -1);
         assert(row->linkpos[c] >= 0);
      }
   }
#endif

   row->activity = RgetRealApprox(rowex->activity);
   row->validactivitylp = stat->lpcount;

   RdeleteTemp(set->buffer, &tmp);
}

 /** calculates the current pseudo activity of a row */
void SCIProwexRecalcPseudoActivity(
   SCIP_ROWEX*           rowex,              /**< row data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_COLEX* colex;
   SCIP_COL* col;
   SCIP_ROW* row;

   SCIP_Rational* tmp = RcreateTemp(set->buffer);
   int i;

   assert(rowex != NULL);

   row = rowex->fprow;

   assert(row != NULL);
   assert(stat != NULL);

   Rset(rowex->pseudoactivity, rowex->constant);
   for( i = 0; i < row->len; ++i )
   {
      col = row->cols[i];
      colex = rowex->cols[i];

      assert(col != NULL);
      assert(colex != NULL);
      assert((i < row->nlpcols) == (row->linkpos[i] >= 0
         && col->lppos >= 0));
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);

      Rmult(tmp, rowex->vals[i], SCIPcolexGetBestBound(colex));
      Radd(rowex->pseudoactivity, rowex->pseudoactivity, tmp);
   }
   row->validpsactivitydomchg = stat->domchgcount;
   row->pseudoactivity = RgetRealApprox(rowex->pseudoactivity);

   RdeleteTemp(set->buffer, &tmp);
}

/** gets objective value of column */
SCIP_Rational* SCIPcolexGetObj(
   SCIP_COLEX*           col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->obj;
}

/** gets lower bound of column */
SCIP_Rational* SCIPcolexGetLb(
   SCIP_COLEX*           col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->lb;
}

/** gets upper bound of column */
SCIP_Rational* SCIPcolexGetUb(
   SCIP_COLEX*           col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->ub;
}

/** gets best bound of column with respect to the objective function */
SCIP_Rational* SCIPcolexGetBestBound(
   SCIP_COLEX*           col                 /**< LP column */
   )
{
   assert(col != NULL);

   if( RisPositive(col->obj) || RisZero(col->obj) )
      return col->lb;
   else
      return col->ub;
}

/** gets the primal LP solution of a column */
SCIP_Rational* SCIPcolexGetPrimsol(
   SCIP_COLEX*           col                 /**< LP column */
   )
{
   assert(col != NULL);

   if( col->fpcol->lppos >= 0 )
      return col->primsol;
   else
      return NULL;
}

/** gets the minimal LP solution value, this column ever assumed */
SCIP_Rational* SCIPcolexGetMinPrimsol(
   SCIP_COLEX*           col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->minprimsol;
}

/** gets the maximal LP solution value, this column ever assumed */
SCIP_Rational* SCIPcolexGetMaxPrimsol(
   SCIP_COLEX*           col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->maxprimsol;
}

/** ensures, that column array of row can store at least num entries */
SCIP_RETCODE SCIProwexEnsureSize(
   SCIP_ROWEX*           row,                /**< LP row */
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
      for( i = row->size; i < newsize; ++i )
         row->vals[i] = Rcreate(blkmem);
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
void getObjvalDeltaObjExact(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Rational*        oldobj,             /**< old objective value of variable */
   SCIP_Rational*        newobj,             /**< new objective value of variable */
   SCIP_Rational*        lb,                 /**< lower bound of variable */
   SCIP_Rational*        ub,                 /**< upper bound of variable */
   SCIP_Rational*        deltaval,           /**< pointer to store the delta value */
   int*                  deltainf            /**< pointer to store the number of variables with infinite best bound */
   )
{
   SCIP_Rational* tmp;
   assert(!RisAbsInfinity(oldobj));
   assert(!RisAbsInfinity(newobj));
   assert(!RisInfinity(lb));
   assert(!RisNegInfinity(ub));
   assert(!RisEqual(oldobj, newobj));

   RsetReal(deltaval, 0);
   tmp = RcreateTemp(set->buffer);
   (*deltainf) = 0;

   if( RisPositive(oldobj) )
   {
      /* sign of objective did not change */
      if( RisPositive(newobj) )
      {
         /* if the bound is finite, calculate the deltaval */
         if( !RisNegInfinity(lb) )
         {
            Rdiff(deltaval, newobj, oldobj);
            Rmult(deltaval, lb, lb);
         }
      }
      /* sign of objective did change, so the best bound does change */
      else if( RisNegative(newobj) )
      {
         if( RisNegInfinity(lb) )
         {
            /* old best bound was infinite while new one is not */
            if( !RisInfinity(ub) )
            {
               (*deltainf) = -1;
               Rmult(deltaval, ub, newobj);
            }
         }
         else
         {
            /* new best bound is infinite while old one was not */
            if( RisInfinity(ub) )
            {
               (*deltainf) = 1;
               Rmult(deltaval, lb, oldobj);
               Rneg(deltaval, deltaval);
            }
            /* neither old nor new best bound is infinite, so just calculate the deltaval */
            else
            {
               Rmult(tmp, lb, oldobj);
               Rmult(deltaval, ub, newobj);

               Rdiff(deltaval, deltaval, tmp);
            }
         }
      }
      /* new objective is 0.0 */
      else
      {
         if( RisNegInfinity(lb) )
            (*deltainf) = -1;
         else
         {
            Rmult(deltaval, lb, oldobj);
            Rneg(deltaval, deltaval);
         }
      }
   }
   else if( RisNegative(oldobj) )
   {
      /* sign of objective did not change */
      if( RisNegative(newobj) )
      {
         /* if the bound is finite, calculate the deltaval */
         if( !RisInfinity(ub) )
         {
            Rdiff(tmp, newobj, oldobj);
            Rmult(deltaval, ub, tmp);
         }
      }
      /* sign of objective did change, so the best bound does change */
      else if( RisPositive(newobj) )
      {
         if( RisInfinity(ub) )
         {
            /* old best bound was infinite while new one is not */
            if( !RisNegInfinity(lb) )
            {
               (*deltainf) = -1;
               Rmult(deltaval, lb, newobj);
            }
         }
         else
         {
            /* new best bound is infinite while old one was not */
            if( RisNegInfinity(lb) )
            {
               (*deltainf) = 1;
               Rmult(deltaval, ub, oldobj);
               Rneg(deltaval, deltaval);
            }
            /* neither old nor new best bound is infinite, so just calculate the deltaval */
            else
            {
               Rmult(tmp, ub, oldobj);
               Rmult(deltaval, lb, newobj);
               Rdiff(deltaval, deltaval, tmp);
            }
         }
      }
      /* new objective is 0.0 */
      else
      {
         if( RisInfinity(ub) )
            (*deltainf) = -1;
         else
         {
            Rmult(deltaval, ub, oldobj);
            Rneg(deltaval, deltaval);
         }
      }
   }
   /* old objective was 0.0 */
   else
   {
      if( RisNegative(newobj) )
      {
         if( RisInfinity(ub) )
            (*deltainf) = 1;
         else
            Rmult(deltaval, ub, newobj);
      }
      else if( RisPositive(newobj) )
      {
         if( RisNegInfinity(lb) )
            (*deltainf) = 1;
         else
            Rmult(deltaval, lb, newobj);
      }
   }

   RdeleteTemp(set->buffer, &tmp);
}

/** returns the left hand side of the row */
SCIP_Rational* SCIProwexGetLhs(
   SCIP_ROWEX*           row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->lhs != NULL);

   return row->lhs;
}

/** returns the right hand side of the row */
SCIP_Rational* SCIProwexGetRhs(
   SCIP_ROWEX*           row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->rhs != NULL);

   return row->rhs;
}

/** compute the objective delta due the new lower bound */
static
void getObjvalDeltaLbExact(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Rational*        obj,                /**< objective value of variable */
   SCIP_Rational*        oldlb,              /**< old lower bound of variable */
   SCIP_Rational*        newlb,              /**< new lower bound of variable */
   SCIP_Rational*        deltaval,           /**< pointer to store the delta value */
   int*                  deltainf            /**< pointer to store the number of variables with infinite best bound */
   )
{
   assert(!RisAbsInfinity(obj));
   assert(!RisInfinity(oldlb));
   assert(!RisNegInfinity(oldlb) || !RisNegInfinity(newlb));
   assert(RisPositive(obj)); /* we only need to update if the objective is positive */

   if( RisNegInfinity(oldlb) )
   {
      if( !RisInfinity(newlb) )
      {
         (*deltainf) = -1;
         Rmult(deltaval, newlb, obj);
      }
      else
      {
         (*deltainf) = 0;
         RsetReal(deltaval, 0.0);
      }
   }
   else if( RisAbsInfinity(newlb) )
   {
      (*deltainf) = 1;
      Rmult(deltaval, oldlb, obj);
      Rneg(deltaval, deltaval);
   }
   else
   {
      (*deltainf) = 0;
      Rdiff(deltaval, newlb, oldlb);
      Rmult(deltaval, deltaval, obj);
   }
}

/** compute the objective delta due the new upper bound */
static
void getObjvalDeltaUbExact(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Rational*        obj,                /**< objective value of variable */
   SCIP_Rational*        oldub,              /**< old upper bound of variable */
   SCIP_Rational*        newub,              /**< new upper bound of variable */
   SCIP_Rational*        deltaval,           /**< pointer to store the delta value */
   int*                  deltainf            /**< pointer to store the number of variables with infinite best bound */
   )
{
   assert(!RisAbsInfinity(obj));
   assert(!RisNegInfinity(oldub));
   assert(!RisInfinity(oldub) || !RisInfinity(newub));
   assert(RisNegative(obj)); /* we only need to update if the objective is negative */

   if( RisInfinity(oldub) )
   {
      if( !RisNegInfinity(newub) )
      {
         (*deltainf) = -1;
         Rmult(deltaval, newub, obj);
      }
      else
      {
         (*deltainf) = 0;
         RsetReal(deltaval, 0.0);
      }
   }
   else if( RisAbsInfinity(newub) )
   {
      (*deltainf) = 1;
      Rmult(deltaval, oldub, obj);
      Rneg(deltaval, deltaval);
   }
   else
   {
      (*deltainf) = 0;
      Rdiff(deltaval, newub, oldub);
      Rmult(deltaval, deltaval, obj);
   }
}

/** updates current pseudo and loose objective values for a change in a variable's objective value or bounds */
static
void lpexUpdateObjval(
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        deltavalex,         /**< delta value in the objective function */
   int                   deltainf,           /**< delta value for the number of variables with infinite best bound */
   SCIP_Bool             local,              /**< should the local pseudo objective value be updated? */
   SCIP_Bool             loose,              /**< should the loose objective value be updated? */
   SCIP_Bool             global              /**< should the global pseudo objective value be updated? */
   )
{
   assert(lp != NULL);
   assert(lp->looseobjvalinf >= 0);
   assert(lp->pseudoobjvalinf >= 0);
   assert(lp->glbpseudoobjvalinf >= 0);

   /* update the pseudo objective value */
   if( local )
   {
      lp->pseudoobjvalinf += deltainf;

      Radd(lp->pseudoobjval, lp->pseudoobjval, deltavalex);

      /* after changing a local bound on a LOOSE variable, we have to update the loose objective value, too */
      if( SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE )
         loose = TRUE;
   }
   /* update the loose objective value */
   if( loose )
   {
      lp->looseobjvalinf += deltainf;

      if( !RisZero(deltavalex) )
         Radd(lp->looseobjval, lp->looseobjval, deltavalex);
   }

   /* update the root pseudo objective values */
   if( global )
   {
      lp->glbpseudoobjvalinf += deltainf;

      Radd(lp->glbpseudoobjval ,lp->glbpseudoobjval, deltavalex);
   }

   assert(lp->looseobjvalinf >= 0);
   assert(lp->pseudoobjvalinf >= 0);
   assert(lp->glbpseudoobjvalinf >= 0);
}

/** updates current pseudo and loose objective value for a change in a variable's objective value */
SCIP_RETCODE SCIPlpexUpdateVarObj(
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldobj,             /**< old objective value of variable */
   SCIP_Rational*        newobj              /**< new objective value of variable */
   )
{
   assert(lp != NULL);
   assert(set != NULL);
   assert(var != NULL);

   if( !RisEqual(oldobj, newobj) )
   {
      SCIP_Rational* deltaval = RcreateTemp(set->buffer);
      int deltainf;

      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetProbindex(var) >= 0);
      /* the objective coefficient can only be changed during presolving, that implies that the global and local
       * domain of the variable are the same
       */
      assert(lp->fplp->probing || RisEqual(SCIPvarGetLbGlobalExact(var), SCIPvarGetLbLocalExact(var)));
      assert(lp->fplp->probing || RisEqual(SCIPvarGetUbGlobalExact(var), SCIPvarGetUbLocalExact(var)));

      /* compute the pseudo objective delta due the new objective coefficient */
      getObjvalDeltaObjExact(set, oldobj, newobj, SCIPvarGetLbLocalExact(var),
          SCIPvarGetUbLocalExact(var), deltaval, &deltainf);

      /* update the local pseudo objective value */
      lpexUpdateObjval(lp, set, var, deltaval, deltainf, TRUE, FALSE, FALSE);

      /* compute the pseudo objective delta due the new objective coefficient */
      getObjvalDeltaObjExact(set, oldobj, newobj, SCIPvarGetLbGlobalExact(var),
          SCIPvarGetUbGlobalExact(var), deltaval, &deltainf);

      /* update the global pseudo objective value */
      lpexUpdateObjval(lp, set, var, deltaval, deltainf, FALSE, FALSE, TRUE);

      RdeleteTemp(set->buffer, &deltaval);
   }

   return SCIP_OKAY;
}

/** updates current root pseudo objective value for a global change in a variable's lower bound */
SCIP_RETCODE SCIPlpexUpdateVarLbGlobal(
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldlb,              /**< old lower bound of variable */
   SCIP_Rational*        newlb               /**< new lower bound of variable */
   )
{
   assert(lp != NULL);
   assert(set != NULL);
   assert(var != NULL);

   if( !RisEqual(oldlb, newlb) && RisPositive(SCIPvarGetObjExact(var)) )
   {
      SCIP_Rational* deltaval = RcreateTemp(set->buffer);
      int deltainf;

      /* compute the pseudo objective delta due the new lower bound */
      getObjvalDeltaLbExact(set, SCIPvarGetObjExact(var), oldlb, newlb, deltaval, &deltainf);

      /* update the root pseudo objective values */
      lpexUpdateObjval(lp, set, var, deltaval, deltainf, FALSE, FALSE, TRUE);

      RdeleteTemp(set->buffer, &deltaval);
   }

   return SCIP_OKAY;
}

/** updates current pseudo and loose objective value for a change in a variable's lower bound */
SCIP_RETCODE SCIPlpexUpdateVarLb(
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldlb,              /**< old lower bound of variable */
   SCIP_Rational*        newlb               /**< new lower bound of variable */
   )
{
   assert(lp != NULL);
   assert(set != NULL);
   assert(var != NULL);

   if( !RisEqual(oldlb, newlb) && RisPositive(SCIPvarGetObjExact(var)) )
   {
      SCIP_Rational* deltaval = RcreateTemp(set->buffer);
      int deltainf;

      assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetProbindex(var) >= 0);

      /* compute the pseudo objective delta due the new lower bound */
      getObjvalDeltaLbExact(set, SCIPvarGetObjExact(var), oldlb, newlb, deltaval, &deltainf);

      /* update the pseudo and loose objective values */
      lpexUpdateObjval(lp, set, var, deltaval, deltainf, TRUE, FALSE, FALSE);

      RdeleteTemp(set->buffer, &deltaval);
   }
}

/** updates current root pseudo objective value for a global change in a variable's upper bound */
SCIP_RETCODE SCIPlpexUpdateVarUbGlobal(
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldub,              /**< old upper bound of variable */
   SCIP_Rational*        newub               /**< new upper bound of variable */
   )
{
   assert(lp != NULL);
   assert(set != NULL);
   assert(var != NULL);

   if( !RisEqual(oldub, newub) && RisNegative(SCIPvarGetObjExact(var)) )
   {
      SCIP_Rational* deltaval = RcreateTemp(set->buffer);
      int deltainf;

      /* compute the pseudo objective delta due the new lower bound */
      getObjvalDeltaUbExact(set, SCIPvarGetObjExact(var), oldub, newub, deltaval, &deltainf);

      /* update the root pseudo objective values */
      lpexUpdateObjval(lp, set, var, deltaval, deltainf, FALSE, FALSE, TRUE);

      RdeleteTemp(set->buffer, &deltaval);
   }

   return SCIP_OKAY;
}

/** updates current pseudo objective value for a change in a variable's upper bound */
SCIP_RETCODE SCIPlpexUpdateVarUb(
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldub,              /**< old upper bound of variable */
   SCIP_Rational*        newub               /**< new upper bound of variable */
   )
{
   assert(lp != NULL);
   assert(set != NULL);
   assert(var != NULL);

   if( !RisEqual(oldub, newub) && RisNegative(SCIPvarGetObjExact(var)) )
   {
      SCIP_Rational* deltaval = RcreateTemp(set->buffer);
      int deltainf;

      assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetProbindex(var) >= 0);

      /* compute the pseudo objective delta due the new lower bound */
      getObjvalDeltaUbExact(set, SCIPvarGetObjExact(var), oldub, newub, deltaval, &deltainf);

      /* update the pseudo and loose objective values */
      lpexUpdateObjval(lp, set, var, deltaval, deltainf, TRUE, FALSE, FALSE);

      RdeleteTemp(set->buffer, &deltaval);
   }
}

/** informs LP, that given variable was added to the problem */
SCIP_RETCODE SCIPlpexUpdateAddVar(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that is now a LOOSE problem variable */
   )
{
   SCIP_Rational* tmp;

   if( !set->misc_exactsolve )
      return SCIP_OKAY;

   assert(lpex != NULL);
   assert(set != NULL);
   assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);

   tmp = RcreateTemp(set->buffer);

   /* add the variable to the loose objective value sum */
   SCIP_CALL( SCIPlpexUpdateVarObj(lpex, set, var, tmp, SCIPvarGetObjExact(var)) );

   /* update the loose variables counter */
   if( SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE )
      lpex->nloosevars++;

   RdeleteTemp(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** informs LP, that given variable is to be deleted from the problem */
SCIP_RETCODE SCIPlpexUpdateDelVar(
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that will be deleted from the problem */
   )
{
   assert(lp != NULL);
   assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);

   /* subtract the variable from the loose objective value sum */
   SCIP_CALL( SCIPlpexUpdateVarObj(lp, set, var, SCIPvarGetObjExact(var), NULL) );

   /* update the loose variables counter */
   if( SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE )
   {
      SCIPlpexDecNLoosevars(lp);
   }

   return SCIP_OKAY;
}

/** informs LP, that given formerly loose problem variable is now a column variable */
SCIP_RETCODE SCIPlpexUpdateVarColumn(
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   )
{
   SCIP_Rational* tmp;
   SCIP_Rational* obj;
   SCIP_Rational* lb;
   SCIP_Rational* ub;

   tmp = RcreateTemp(set->buffer);

   assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(lp->looseobjvalinf >= 0);

   obj = SCIPvarGetObjExact(var);

   /* update loose objective value */
   if( RisPositive(obj) )
   {
      lb = SCIPvarGetLbLocalExact(var);
      if( RisNegInfinity(lb) )
         lp->looseobjvalinf--;
      else
      {
         Rneg(tmp, lb);
         Rmult(tmp, tmp, obj);
         lpexUpdateObjval(lp, set, var, tmp, 0, FALSE, TRUE, FALSE);
      }
   }
   else if( RisNegative(obj) )
   {
      ub = SCIPvarGetUbLocalExact(var);
      if( RisInfinity(ub) )
         lp->looseobjvalinf--;
      else
      {
         Rneg(tmp, ub);
         Rmult(tmp, tmp, obj);
         lpexUpdateObjval(lp, set, var, tmp, 0, FALSE, TRUE, FALSE);
      }
   }

   assert(lp->looseobjvalinf >= 0);

   RdeleteTemp(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** informs LP, that given formerly column problem variable is now again a loose variable */
SCIP_RETCODE SCIPlpexUpdateVarLoose(
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   )
{
   SCIP_Rational* tmp;
   SCIP_Rational* obj;
   SCIP_Rational* lb;
   SCIP_Rational* ub;

   tmp = RcreateTemp(set->buffer);

   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(lp->looseobjvalinf >= 0);

   obj = SCIPvarGetObjExact(var);

   /* update loose objective value corresponding to the addition of variable */
   if( RisPositive(obj) )
   {
      lb = SCIPvarGetLbLocalExact(var);
      if( RisNegInfinity(lb) )
         lp->looseobjvalinf++;
      else
      {
         Rmult(tmp, lb, obj);
         lpexUpdateObjval(lp, set, var, tmp, 0, FALSE, TRUE, FALSE);
      }
   }
   else if( RisNegative(obj) )
   {
      ub = SCIPvarGetUbLocalExact(var);
      if( RisInfinity(ub) )
         lp->looseobjvalinf++;
      else
      {
         Rmult(tmp, ub, obj);
         lpexUpdateObjval(lp, set, var, tmp, 0, FALSE, TRUE, FALSE);
      }
   }
   lp->nloosevars++;

   assert(lp->looseobjvalinf >= 0);

   RdeleteTemp(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** decrease the number of loose variables by one */
void SCIPlpexDecNLoosevars(
   SCIP_LPEX*            lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->nloosevars > 0);

   lp->nloosevars--;

   /* get rid of numerical problems: set loose objective value explicitly to zero, if no loose variables remain */
   if( lp->nloosevars == 0 )
   {
      assert(lp->looseobjvalinf == 0);
      RsetReal(lp->looseobjval, 0.0);
   }
}

/** get the number of rows currently in the lp */
SCIP_RETCODE SCIPlexGetNRows(
   SCIP_LPEX*            lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->nrows;
}

/** stores the LP solution in the columns and rows */
SCIP_RETCODE SCIPlpexGetSol(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            dualfeasible        /**< pointer to store whether the solution is dual feasible, or NULL */
   )
{
   return SCIP_ERROR;
}

/** stores LP solution with infinite objective value in the columns and rows */
SCIP_RETCODE SCIPlpexGetUnboundedSol(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            rayfeasible         /**< pointer to store whether the primal ray is a feasible unboundedness proof, or NULL */
   );

/** returns primal ray proving the unboundedness of the current LP */
SCIP_RETCODE SCIPlpexGetPrimalRay(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Rational**       ray                 /**< array for storing primal ray values, they are stored w.r.t. the problem index of the variables,
                                              *   so the size of this array should be at least number of active variables
                                              *   (all entries have to be initialized to 0 before) */
   );

/** stores the dual Farkas multipliers for infeasibility proof in rows. besides, the proof is checked for validity if
 *  lp/checkfarkas = TRUE.
 *
 *  @note the check will not be performed if @p valid is NULL.
 */
SCIP_RETCODE SCIPlpexGetDualfarkas(
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            valid               /**< pointer to store whether the Farkas proof is valid  or NULL */
   );

/** removes all columns after the given number of cols from the LP */
SCIP_RETCODE SCIPlpexShrinkCols(
   SCIP_LPEX*            lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   newncols            /**< new number of columns in the LP */
   )
{
   SCIP_COLEX* col;
   int c;

   assert(lp != NULL);

   SCIPsetDebugMsg(set, "shrinking LP from %d to %d columns\n", lp->ncols, newncols);
   assert(0 <= newncols);
   assert(newncols <= lp->ncols);

   if( newncols < lp->ncols )
   {
      assert(!lp->fplp->diving);

      for( c = lp->ncols-1; c >= newncols; --c )
      {
         col = lp->cols[c];
         assert(col != NULL);
         assert(col->len == 0 || col->rows != NULL);
         assert(col->var != NULL);
         assert(SCIPvarGetStatusExact(col->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetColExact(col->var) == lp->cols[c]);
         assert(col->lppos == c);

         /* mark column to be removed from the LP */
         col->lppos = -1;
         lp->ncols--;

         /* count removable columns */
         if( col->removable )
            lp->nremovablecols--;

         /* update column arrays of all linked rows */
         /* colUpdateDelLP(col, set); */
      }

   }

   return SCIP_OKAY;
}

/** removes and releases all rows after the given number of rows from the LP */
SCIP_RETCODE SCIPlpexShrinkRows(
   SCIP_LPEX*            lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   int                   newnrows            /**< new number of rows in the LP */
   )
{
   SCIP_ROWEX* row;
   int r;

   assert(lp != NULL);
   assert(0 <= newnrows && newnrows <= lp->nrows);

   SCIPsetDebugMsg(set, "shrinking exact LP from %d to %d rows\n", lp->nrows, newnrows);
   if( newnrows < lp->nrows )
   {
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

         /* count removable rows */
         if( row->removable )
            lp->nremovablerows--;

         /* rowexUpdateDelLP(row);

         SCIProwexUnlocK(row); */
         SCIP_CALL( SCIProwexRelease(&lp->rows[r], blkmem, set, lp) );
      }
   }

   return SCIP_OKAY;
}

/** resets the LP to the empty LP by removing all columns and rows from LP, releasing all rows, and flushing the
 *  changes to the LP solver
 */
SCIP_RETCODE SCIPlpexReset(
   SCIP_LPEX*            lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   )
{
   if( !set->misc_exactsolve )
      return SCIP_OKAY;
      
   assert(stat != NULL);

   SCIP_CALL( SCIPlpexClear(lp, blkmem, set, eventqueue, eventfilter) );
   //SCIP_CALL( SCIPlpexFlush(lp, blkmem, set, eventqueue) );

   /* mark the empty LP to be solved */
   lp->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   RsetReal(lp->lpobjval, 0.0);
   //lp->validsollp = stat->lpcount; /* the initial (empty) SCIP_LP is solved with primal and dual solution of zero */
   //lp->validfarkaslp = -1;
   //lp->validsoldirlp = -1;
   lp->validsoldirsol = NULL;
   lp->solved = TRUE;
   lp->primalfeasible = TRUE;
   lp->primalchecked = TRUE;
   lp->dualfeasible = TRUE;
   lp->dualchecked = TRUE;
   lp->solisbasic = FALSE;
   lp->lastlpalgo = SCIP_LPALGO_DUALSIMPLEX;

   return SCIP_OKAY;
}

/** removes all columns and rows from LP, releases all rows */
SCIP_RETCODE SCIPlpexClear(
   SCIP_LPEX*            lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   )
{
   assert(lp != NULL);
   assert(!lp->fplp->diving);

   SCIPsetDebugMsg(set, "clearing LP\n");
   SCIP_CALL( SCIPlpexShrinkCols(lp, set, 0) );
   SCIP_CALL( SCIPlpexShrinkRows(lp, blkmem, set, eventqueue, eventfilter, 0) );

   return SCIP_OKAY;
}
