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

#include "sort.h"
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

struct BoundChg                         /**< change in one bound of a variable */
{
   COL*             col;                /**< column to change the bounds for */
   Real             newbound;           /**< new value for bound */
   Real             oldbound;           /**< old value for bound */
   BOUNDTYPE        boundtype;          /**< type of bound: lower or upper bound */
};

struct HoleChg                          /**< change in a hole list */
{
   HOLELIST**       ptr;                /**< changed list pointer */
   HOLELIST*        newlist;            /**< new value of list pointer */
   HOLELIST*        oldlist;            /**< old value of list pointer */
};

struct DomChg                           /**< tracks changes of the variable's domains (fixed sized arrays) */
{
   BOUNDCHG*        boundchg;           /**< array with changes in bounds of variables */
   HOLECHG*         holechg;            /**< array with changes in hole lists */
   int              nboundchg;          /**< number of bound changes */
   int              nholechg;           /**< number of hole list changes */
};

struct DomChgDyn                        /**< tracks changes of the variable's domains (dynamically sized arrays) */
{
   DOMCHG           domchg;             /**< domain changes */
   int              boundchgsize;       /**< size of bound change array */
   int              holechgsize;        /**< size of hole change array */
};

struct ColList                          /**< list of columns */
{
   COL*             col;                /**< pointer to this column */
   COLLIST*         next;               /**< pointer to next collist entry */
};

struct RowList                          /**< list of rows */
{
   ROW*             row;                /**< pointer to this row */
   ROWLIST*         next;               /**< pointer to next rowlist entry */
};


/*
 * memory growing methods for dynamically allocated arrays
 */

static
RETCODE ensureBoundchgSize(             /**< ensures, that boundchg array can store at least num entries */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(domchgdyn->domchg.nboundchg <= domchgdyn->boundchgsize);
   
   if( num > domchgdyn->boundchgsize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(mem->dommem, domchgdyn->domchg.boundchg, domchgdyn->boundchgsize, newsize) );
      domchgdyn->boundchgsize = newsize;
   }
   assert(num <= domchgdyn->boundchgsize);

   return SCIP_OKAY;
}

static
RETCODE ensureHolechgSize(              /**< ensures, that holechg array can store at least num additional entries */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of additional entries to store */
   )
{
   assert(domchgdyn->domchg.nholechg <= domchgdyn->holechgsize);
   
   if( num > domchgdyn->holechgsize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(mem->dommem, domchgdyn->domchg.holechg, domchgdyn->holechgsize, newsize) );
      domchgdyn->holechgsize = newsize;
   }
   assert(num <= domchgdyn->holechgsize);

   return SCIP_OKAY;
}

static
RETCODE ensureChgbdsSize(               /**< ensures, that chgbds array can store at least num entries */
   LP*              lp,                 /**< actual LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nchgbds <= lp->chgbdssize);
   
   if( num > lp->chgbdssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(mem->lpmem, lp->chgbds, lp->chgbdssize, newsize) );
      lp->chgbdssize = newsize;
   }
   assert(num <= lp->chgbdssize);

   return SCIP_OKAY;
}

static
RETCODE ensureAddColsSize(              /**< ensures, that addcols array can store at least num entries */
   LP*              lp,                 /**< actual LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->naddcols <= lp->addcolssize);
   
   if( num > lp->addcolssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(mem->lpmem, lp->addcols, lp->addcolssize, newsize) );
      lp->addcolssize = newsize;
   }
   assert(num <= lp->addcolssize);

   return SCIP_OKAY;
}

static
RETCODE ensureAddRowsSize(              /**< ensures, that addrows array can store at least num entries */
   LP*              lp,                 /**< actual LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(lp->naddrows <= lp->addrowssize);
   
   if( num > lp->addrowssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(mem->lpmem, lp->addrows, lp->addrowssize, newsize) );
      lp->addrowssize = newsize;
   }
   assert(num <= lp->addrowssize);

   return SCIP_OKAY;
}


/*
 * domain changes
 */

DOMCHGDYN* SCIPdomchgdynCreate(         /**< creates a dynamically sized domain change data structure */
   MEM*             mem                 /**< block memory buffers */   
   )
{
   DOMCHGDYN* domchgdyn;

   assert(mem != NULL);

   ALLOC_NULL( allocBlockMemory(mem->dommem, domchgdyn) );
   domchgdyn->domchg.boundchg = NULL;
   domchgdyn->domchg.holechg = NULL;
   domchgdyn->domchg.nboundchg = 0;
   domchgdyn->domchg.nholechg = 0;
   domchgdyn->boundchgsize = 0;
   domchgdyn->holechgsize = 0;

   return domchgdyn;
};

void SCIPdomchgdynFree(                 /**< frees a dynamically sized domain change data structure */
   DOMCHGDYN**      domchgdyn,          /**< pointer to dynamically sized domain change data structure */
   MEM*             mem                 /**< block memory buffers */   
   )
{
   assert(domchgdyn != NULL);
   assert(*domchgdyn != NULL);
   assert(mem != NULL);

   freeBlockMemoryArrayNull(mem->dommem, (*domchgdyn)->domchg.boundchg, (*domchgdyn)->boundchgsize);
   freeBlockMemoryArrayNull(mem->dommem, (*domchgdyn)->domchg.holechg, (*domchgdyn)->holechgsize);
   freeBlockMemory(mem->dommem, *domchgdyn);
}

RETCODE SCIPdomchgdynAddBoundchg(       /**< adds bound change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEM*             mem,                /**< block memory buffers */   
   const SET*       set,                /**< global SCIP settings */
   COL*             col,                /**< column to change the bounds for */
   Real             newbound,           /**< new value for bound */
   Real             oldbound,           /**< old value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   assert(domchgdyn != NULL);
   assert(mem != NULL);
   assert(col != NULL);

   CHECK_OKAY( ensureBoundchgSize(domchgdyn, mem, set, domchgdyn->boundchgsize+1) );
   domchgdyn->domchg.boundchg[domchgdyn->domchg.nboundchg].col = col;
   domchgdyn->domchg.boundchg[domchgdyn->domchg.nboundchg].newbound = newbound;
   domchgdyn->domchg.boundchg[domchgdyn->domchg.nboundchg].oldbound = oldbound;
   domchgdyn->domchg.boundchg[domchgdyn->domchg.nboundchg].boundtype = boundtype;

   return SCIP_OKAY;
}

RETCODE SCIPdomchgdynAddHolechg(        /**< adds hole change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEM*             mem,                /**< block memory buffers */   
   const SET*       set,                /**< global SCIP settings */
   HOLELIST**       ptr,                /**< changed list pointer */
   HOLELIST*        newlist,            /**< new value of list pointer */
   HOLELIST*        oldlist             /**< old value of list pointer */
   )
{
   assert(domchgdyn != NULL);
   assert(mem != NULL);
   assert(ptr != NULL);

   CHECK_OKAY( ensureHolechgSize(domchgdyn, mem, set, domchgdyn->holechgsize+1) );
   domchgdyn->domchg.holechg[domchgdyn->domchg.nholechg].ptr = ptr;
   domchgdyn->domchg.holechg[domchgdyn->domchg.nholechg].newlist = newlist;
   domchgdyn->domchg.holechg[domchgdyn->domchg.nholechg].oldlist = oldlist;

   return SCIP_OKAY;
}

DOMCHG* SCIPdomchgCreate(               /**< creates domain change data (fixed size) from dynamically sized data */
   MEM*             mem,                /**< block memory buffers */   
   const DOMCHGDYN* domchgdyn           /**< dynamically sized domain change data structure */
   )
{
   DOMCHG* domchg;

   assert(mem != NULL);
   assert(domchgdyn != NULL);

   ALLOC_NULL( allocBlockMemory(mem->dommem, domchg) );

   if( domchgdyn->domchg.nboundchg > 0 )
   {
      ALLOC_NULL( duplicateBlockMemoryArray(mem->dommem, domchg->boundchg, domchgdyn->domchg.boundchg,
                     domchgdyn->domchg.nboundchg) );
   }
   else
      domchg->boundchg = NULL;

   if( domchgdyn->domchg.nholechg > 0 )
   {
      ALLOC_NULL( duplicateBlockMemoryArray(mem->dommem, domchg->holechg, domchgdyn->domchg.holechg,
                     domchgdyn->domchg.nholechg) );
   }
   else
      domchg->holechg = NULL;

   domchg->nboundchg = domchgdyn->domchg.nboundchg;
   domchg->nholechg = domchgdyn->domchg.nholechg;

   return domchg;
}

void SCIPdomchgFree(                    /**< frees domain change data */
   DOMCHG**         domchg,             /**< pointer to domain change */
   MEM*             mem                 /**< block memory buffers */
   )
{
   assert(domchg != NULL);
   assert(*domchg != NULL);
   assert(mem != NULL);

   freeBlockMemoryArrayNull(mem->dommem, (*domchg)->boundchg, (*domchg)->nboundchg);
   freeBlockMemoryArrayNull(mem->dommem, (*domchg)->holechg, (*domchg)->nholechg);
   freeBlockMemory(mem->dommem, *domchg);
}

RETCODE SCIPlpApplyDomchg(              /**< applies domain change */
   LP*              lp,                 /**< actual LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   const DOMCHG*    domchg              /**< domain change to apply */
   )
{
   COL* col;
   int i;

   assert(lp != NULL);
   assert(domchg != NULL);

   /* apply bound changes */
   for( i = 0; i < domchg->nboundchg; ++i )
   {
      col = domchg->boundchg[i].col;
      /* insert col in the chgbds list (if not already there) */
      if( !col->lbchanged && !col->ubchanged )
      {
         CHECK_OKAY( ensureChgbdsSize(lp, mem, set, lp->nchgbds+1) );
         lp->chgbds[lp->nchgbds] = col;
         lp->nchgbds++;
      }

      /* apply bound change to the LP data */
      switch( domchg->boundchg[i].boundtype )
      {
      case SCIP_BOUNDTYPE_LOWER:
         col->dom.lb = domchg->boundchg[i].newbound;
         col->lbchanged = TRUE;
         break;
      case SCIP_BOUNDTYPE_UPPER:
         col->dom.ub = domchg->boundchg[i].newbound;
         col->ubchanged = TRUE;
         break;
      default:
         errorMessage("Unknown bound type");
         abort();
      }
   }

   /* apply holelist changes */
   for( i = 0; i < domchg->nholechg; ++i )
      *(domchg->holechg[i].ptr) = domchg->holechg[i].newlist;

   return SCIP_OKAY;
}
   
RETCODE SCIPlpUndoDomchg(               /**< undoes domain change */
   LP*              lp,                 /**< actual LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   const DOMCHG*    domchg              /**< domain change to remove */
   )
{
   COL* col;
   int i;

   assert(lp != NULL);
   assert(domchg != NULL);

   /* undo bound changes */
   for( i = domchg->nboundchg-1; i >= 0; --i )
   {
      col = domchg->boundchg[i].col;
      /* insert col in the chgbds list (if not already there) */
      if( !col->lbchanged && !col->ubchanged )
      {
         CHECK_OKAY( ensureChgbdsSize(lp, mem, set, lp->nchgbds+1) );
         lp->chgbds[lp->nchgbds] = col;
         lp->nchgbds++;
      }

      /* apply bound change to the LP data */
      switch( domchg->boundchg[i].boundtype )
      {
      case SCIP_BOUNDTYPE_LOWER:
         col->dom.lb = domchg->boundchg[i].oldbound;
         col->lbchanged = TRUE;
         break;
      case SCIP_BOUNDTYPE_UPPER:
         col->dom.ub = domchg->boundchg[i].oldbound;
         col->ubchanged = TRUE;
         break;
      default:
         errorMessage("Unknown bound type");
         abort();
      }
   }

   /* undo holelist changes */
   for( i = domchg->nholechg-1; i >= 0; --i )
      *(domchg->holechg[i].ptr) = domchg->holechg[i].oldlist;

   return SCIP_OKAY;
}
   
static
RETCODE lpFlushDomchg(                  /**< applies all cached domain changes to the LP */
   LP*              lp,                 /**< actual LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   )
{
   COL* col;
   int* ind;
   char* lu;
   Real* bd;
   int i;
   int nchg;

   assert(lp != NULL);
   assert(mem != NULL);

   if( lp->nchgbds == 0 )
      return SCIP_OKAY;

   /* get temporary memory for bound changes */
   ALLOC_OKAY( ind = SCIPmemGetIntbuf(mem, set, 2*lp->ncols) );
   ALLOC_OKAY( lu = SCIPmemGetCharbuf(mem, set, 2*lp->ncols) );
   ALLOC_OKAY( bd = SCIPmemGetRealbuf(mem, set, 2*lp->ncols) );

   /* collect all cached bound changes */
   nchg = 0;
   for( i = 0; i < lp->nchgbds; ++i )
   {
      col = lp->chgbds[i];
      assert(col != NULL);

      if( col->lppos >= 0 )
      {
         if( col->lbchanged )
         {
            assert(nchg < 2*lp->ncols);
            ind[nchg] = col->lppos;
            lu[nchg] = 'L';
            bd[nchg] = col->dom.lb;
            nchg++;
         }
         if( col->ubchanged )
         {
            assert(nchg < 2*lp->ncols);
            ind[nchg] = col->lppos;
            lu[nchg] = 'U';
            bd[nchg] = col->dom.ub;
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

   return SCIP_OKAY;
}

static
RETCODE lpFlushAddcols(                 /**< applies all cached column additions to the LP */
   LP*              lp,                 /**< actual LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   )
{
   Real* realbuf;
   int* intbuf;
   Real* obj;
   Real* lb;
   Real* ub;
   int* beg;
   int* ind;
   Real* val;
   char** name;
   COL* col;
   int c;
   int nnonz;
   int i;
   int lppos;

   assert(lp != NULL);
   assert(mem != NULL);

   if( lp->naddcols == 0 )
      return SCIP_OKAY;

   /* get temporary memory for bound changes */
   ALLOC_OKAY( realbuf = SCIPmemGetRealbuf(mem, set, 3*lp->naddcols + lp->addcolscoefs) );
   obj = &(realbuf[0*lp->naddcols]);
   lb = &(realbuf[1*lp->naddcols]);
   ub = &(realbuf[2*lp->naddcols]);
   val = &(realbuf[3*lp->naddcols]);
   ALLOC_OKAY( intbuf = SCIPmemGetIntbuf(mem, set, lp->naddcols + lp->addcolscoefs) );
   beg = &(intbuf[0*lp->naddcols]);
   ind = &(intbuf[1*lp->naddcols]);
   ALLOC_OKAY( name = (char**)(SCIPmemGetPtrbuf(mem, set, lp->naddcols)) );
   
   /* fill temporary memory with column data */
   nnonz = 0;
   for( c = 0; c < lp->naddcols; ++c )
   {
      col = lp->addcols[c];
      obj[c] = col->obj;
      lb[c] = col->dom.lb;
      ub[c] = col->dom.ub;
      col->lbchanged = FALSE;
      col->ubchanged = FALSE;
      beg[c] = nnonz;
      name[c] = col->name;

      for( i = 0; i < col->len; ++i )
      {
         lppos = col->row[i]->lppos;
         if( lppos >= 0 )
         {
            assert(lppos < lp->nrows);
            ind[nnonz] = lppos;
            val[nnonz] = col->val[i];
            nnonz++;
         }
      }
   }

   /* call LP interface */
   CHECK_OKAY( SCIPlpiAddCols(lp->lpi, lp->naddcols, nnonz, obj, lb, ub, beg, ind, val, name) );

   lp->naddcols = 0;

   return SCIP_OKAY;
}

static
RETCODE lpFlushAddrows(                 /**< applies all cached row additions to the LP */
   LP*              lp,                 /**< actual LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   )
{
   Real* realbuf;
   int* intbuf;
   Real* rhs;
   char* sen;
   int* beg;
   int* ind;
   Real* val;
   char** name;
   ROW* row;
   int r;
   int nnonz;
   int i;
   int lppos;

   assert(lp != NULL);
   assert(mem != NULL);

   if( lp->naddrows == 0 )
      return SCIP_OKAY;

   /* get temporary memory for bound changes */
   ALLOC_OKAY( realbuf = SCIPmemGetRealbuf(mem, set, lp->naddrows + lp->addrowscoefs) );
   rhs = &(realbuf[0*lp->naddrows]);
   val = &(realbuf[1*lp->naddrows]);
   ALLOC_OKAY( intbuf = SCIPmemGetIntbuf(mem, set, lp->naddrows + lp->addrowscoefs) );
   beg = &(intbuf[0*lp->naddrows]);
   ind = &(intbuf[1*lp->naddrows]);
   ALLOC_OKAY( sen = SCIPmemGetCharbuf(mem, set, lp->naddrows) );
   ALLOC_OKAY( name = (char**)(SCIPmemGetPtrbuf(mem, set, lp->naddrows)) );
   
   /* fill temporary memory with row data */
   nnonz = 0;
   for( r = 0; r < lp->naddrows; ++r )
   {
      row = lp->addrows[r];
      rhs[r] = row->rhs;
      sen[r] = row->equality ? 'E' : 'L';
      beg[c] = nnonz;
      name[c] = row->name;

      for( i = 0; i < row->len; ++i )
      {
         lppos = row->col[i]->lppos;
         if( lppos >= 0 )
         {
            assert(lppos < lp->ncols);
            ind[nnonz] = lppos;
            val[nnonz] = row->val[i];
            nnonz++;
         }
      }
   }

   /* call LP interface */
   CHECK_OKAY( SCIPlpiAddRows(lp->lpi, lp->naddrows, nnonz, rhs, sen, beg, ind, val, name) );

   lp->naddrows = 0;

   return SCIP_OKAY;
}

RETCODE SCIPlpFlush(                    /**< applies all cached changes to the LP */
   LP*              lp,                 /**< actual LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(mem != NULL);
   
   if( lp->flushed )
      return SCIP_OKAY;

   CHECK_OKAY( lpFlushDomchg(lp, mem, set) );
   CHECK_OKAY( lpFlushAddcols(lp, mem, set) );
   CHECK_OKAY( lpFlushAddrows(lp, mem, set) );

   lp->flushed = TRUE;

   return SCIP_OKAY;
}


static
void addColSign(                        /**< update column sign after addition of new coefficient */
   COL*             col,                /**< LP column */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value of new coefficient */
   )
{   
   assert(col != NULL);
   assert(col->numpos >= 0 && col->numneg >= 0);
   assert(!SCIPisZero(set, val));

   if( SCIPisPos(set, val) )
      col->numpos++;
   else
   {
      assert(SCIPisNeg(set, val));
      col->numneg++;
   }
}

static
void delColSign(                        /**< update column sign after deletion of coefficient */
   COL*             col,                /**< LP column */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value of deleted coefficient */
   )
{
   assert(col != NULL);
   assert(!SCIPisZero(set, val));

   if( SCIPisPos(set, val) )
      col->numpos--;
   else
   {
      assert(SCIPisNeg(set, val));
      col->numneg--;
   }

   assert(col->numpos >= 0 && col->numneg >= 0);
}

static
void addRowNorms(                       /**< update row norms after addition of new coefficient */
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
   assert(!SCIPisZero(set, absval));

   /* update min/maxidx */
   row->minidx = MIN(row->minidx, colidx);
   row->maxidx = MAX(row->maxidx, colidx);

   /* update squared euclidean norm */
   row->sqrnorm += SQR(absval);

   /* update maximum norm */
   if( SCIPisG(set, absval, row->maxval) )
   {
      row->maxval = absval;
      row->nummaxval = 1;
   }
   else if( SCIPisGE(set, absval, row->maxval) )
   {
      assert(row->nummaxval >= 1);
      row->nummaxval++;
   }
}

static
void calcRowNorms(                      /**< calculates row norms and min/maxidx from scratch, and checks for sortation */
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
      assert(!SCIPisZero(set, row->val[i]));
      idx = row->col[i]->index;
      addRowNorms(row, set, idx, row->val[i]);
      row->sorted &= (i == 0 || row->col[i-1]->index < idx);
   }
}

static
void delRowNorms(                       /**< update row norms after deletion of coefficient */
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
   assert(!SCIPisZero(set, absval));
   assert(SCIPisGE(set, row->maxval, absval));

   /* update min/maxidx validity */
   if( colidx == row->minidx || colidx == row->maxidx )
      row->validminmaxidx = FALSE;

   /* update squared euclidean norm */
   row->sqrnorm -= SQR(absval);
   assert(SCIPisGE(set, row->sqrnorm, 0.0));

   /* update maximum norm */
   if( SCIPisGE(set, absval, row->maxval) )
   {
      row->nummaxval--;
      if( row->nummaxval == 0 )
         calcRowNorms(row, set);
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
int searchColCoeff(                     /**< searches existing coefficient in column, returns position in col vector */
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
   SCIPcolSort(col);
   assert(col->sorted);

   /* binary search */
   searchidx = row->index;
   minpos = 0;
   maxpos = col->len-1;
   actpos = 0;
   while(minpos < maxpos)
   {
      actpos = (minpos + maxpos)/2;
      actidx = col->row[actpos]->index;
      if( searchidx == actidx )
         break;
      else if( searchidx < actidx )
         maxpos = actpos-1;
      else
         minpos = actpos+1;
      assert(minpos <= maxpos);
   }
   assert(minpos <= actpos && actpos <= maxpos);
   assert(col->row[actpos]->index == searchidx);

   return actpos;
}

static
int searchRowCoeff(                     /**< searches existing coefficient in row, returns position in row vector */
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
   SCIProwSort(row);
   assert(row->sorted);

   /* binary search */
   searchidx = col->index;
   minpos = 0;
   maxpos = row->len-1;
   actpos = 0;
   while(minpos < maxpos)
   {
      actpos = (minpos + maxpos)/2;
      actidx = row->col[actpos]->index;
      if( searchidx == actidx )
         break;
      else if( searchidx < actidx )
         maxpos = actpos-1;
      else
         minpos = actpos+1;
      assert(minpos <= maxpos);
   }
   assert(minpos <= actpos && actpos <= maxpos);
   assert(0 <= actpos && actpos < row->len);
   assert(row->col[actpos]->index == searchidx);

   return actpos;
}

static
RETCODE ensureColSize(                  /**< ensures, that row array of column can store at least num additional entries */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   COL*             col,                /**< LP column */
   int              num                 /**< minimum number of additional entries to store */
   )
{
   assert(col->len <= col->size);
   
   if( num > col->size )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(mem->lpmem, col->row, col->size, newsize) );
      col->size = newsize;
   }
   assert(num <= col->size);

   return SCIP_OKAY;
}

static
RETCODE ensureRowSize(                  /**< ensures, that column array of row can store at least num additional entries */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row,                /**< LP row */
   int              num                 /**< minimum number of additional entries to store */
   )
{
   assert(row->len <= row->size);
   
   if( num > row->size )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(mem->lpmem, row->col, row->size, newsize) );
      row->size = newsize;
   }
   assert(num <= row->size);

   return SCIP_OKAY;
}

static
RETCODE addColCoeff(                    /**< adds a coefficient to an LP column */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   COL*             col,                /**< LP column */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   )
{
   assert(mem != NULL);
   assert(col != NULL);
   assert(row != NULL);
   assert(!SCIPisZero(set, val));

   CHECK_OKAY( ensureColSize(mem, set, col, col->len+1) );
   col->row[col->len] = row;
   col->val[col->len] = val;
   col->len++;

   addColSign(col, set, val);

   return SCIP_OKAY;
}

static
RETCODE addRowCoeff(                    /**< adds a coefficient to an LP row */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row,                /**< LP row */
   COL*             col,                /**< LP column */
   Real             val                 /**< value of coefficient */
   )
{
   assert(mem != NULL);
   assert(row != NULL);
   assert(col != NULL);
   assert(!SCIPisZero(set, val));

   CHECK_OKAY( ensureRowSize(mem, set, row, row->len+1) );
   row->col[row->len] = col;
   row->val[row->len] = val;
   row->len++;

   addRowNorms(row, set, col->index, val);

   row->sorted &= (row->len == 0 || row->col[row->len-1]->index < col->index);

   return SCIP_OKAY;
}

static
void deleteColCoeff(                    /**< deletes coefficient from column */
   COL*             col,                /**< column to be changed */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row                 /**< coefficient to be deleted */
   )
{
   int pos;
   Real val;

   assert(col != NULL);
   assert(col->len > 0);
   assert(row != NULL);

   pos = searchColCoeff(col, row);
   assert(col->row[pos] == row);
   val = col->val[pos];

   if( pos < col->len-1 )
   {
      /* move last coefficient to position of deleted coefficient */
      col->row[pos] = col->row[col->len-1];
      col->val[pos] = col->val[col->len-1];
      col->sorted = FALSE;
   }
   col->len--;

   delColSign(col, set, val);
}

static
void deleteRowCoeff(                    /**< deletes coefficient from row */
   ROW*             row,                /**< row to be changed */
   const SET*       set,                /**< global SCIP settings */
   COL*             col                 /**< coefficient to be deleted */
   )
{
   int pos;
   Real val;

   assert(row != NULL);
   assert(row->len > 0);
   assert(col != NULL);

   pos = searchRowCoeff(row, col);
   assert(row->col[pos] == col);
   val = row->val[pos];

   if( pos < row->len-1 )
   {
      /* move last coefficient to position of deleted coefficient */
      row->col[pos] = row->col[row->len-1];
      row->val[pos] = row->val[row->len-1];
      row->sorted = FALSE;
   }
   row->len--;
   
   delRowNorms(row, set, col->index, val);
}

COL* SCIPcolCreate(                     /**< creates an LP column */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   char*            name,               /**< name of column */
   int              len,                /**< number of nonzeros in the column */
   ROW**            row,                /**< array with rows of column entries */
   Real*            val,                /**< array with coefficients of column entries */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   )
{
   COL* col;
   int i;
   int idx;

   assert(mem != NULL);
   assert(stat != NULL);
   assert(len >= 0);
   assert(len == 0 || (row != NULL && val != NULL));

   ALLOC_NULL( allocBlockMemory(mem->lpmem, col) );

   if( len > 0 )
   {
      ALLOC_NULL( duplicateBlockMemoryArray(mem->lpmem, col->row, row, len) );
      ALLOC_NULL( duplicateBlockMemoryArray(mem->lpmem, col->val, val, len) );
   }
   else
   {
      col->row = NULL;
      col->val = NULL;
   }

   col->name = name;
   col->dom.holelist = NULL;
   col->dom.lb = lb;
   col->dom.ub = ub;
   col->problb = lb;
   col->probub = ub;
   col->obj = obj;
   col->index = stat->numcolidx++;
   col->size = len;
   col->len = len;
   col->numuses = 0;
   col->lppos = -1;
   col->numpos = 0;
   col->numneg = 0;
   col->vartype = vartype;
   col->sorted = TRUE;
   col->lbchanged = FALSE;
   col->ubchanged = FALSE;

   /* check, if column is sorted
    * update number of positive/negative entries
    * insert coefficients in corresponding rows
    */
   for( i = 0; i < len; ++i )
   {
      assert(!SCIPisZero(set, col->val[i]));
      col->sorted &= (i == 0 || col->row[i-1]->index < col->row[i]->index);
      addColSign(col, set, col->val[i]);
      addRowCoeff(mem, set, col->row[i], col, col->val[i]);
   }

   return col;
}

void SCIPcolFree(                       /**< frees an LP column */
   COL**            col,                /**< pointer to LP column */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   )
{
   int i;

   assert(mem != NULL);
   assert(col != NULL);
   assert(*col != NULL);
   assert((*col)->numuses == 0);
   
   /* delete coefficients in corresponding rows */
   for( i = 0; i < (*col)->len; ++i )
      deleteRowCoeff((*col)->row[i], set, *col);

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
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(mem != NULL);
   assert(col != NULL);
   assert(*col != NULL);
   assert((*col)->numuses >= 1);

   (*col)->numuses--;
   if( (*col)->numuses == 0 )
      SCIPcolFree(col, mem, set);
}

ROW* SCIProwCreate(                     /**< creates an LP row */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   char*            name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Bool             equality,           /**< TRUE iff row is an equality, FALSE iff row is a less or equal inequality */
   Real             rhs,                /**< right hand side of row */
   Real             epsilon             /**< maximal normed violation of row */
   )
{
   ROW* row;
   int i;
   int idx;

   assert(mem != NULL);
   assert(stat != NULL);
   assert(len >= 0);
   assert(len == 0 || (col != NULL && val != NULL));
   assert(epsilon >= 0.0);

   ALLOC_NULL( allocBlockMemory(mem->lpmem, row) );

   if( len > 0 )
   {
      ALLOC_NULL( duplicateBlockMemoryArray(mem->lpmem, row->col, col, len) );
      ALLOC_NULL( duplicateBlockMemoryArray(mem->lpmem, row->val, val, len) );
   }
   else
   {
      row->col = NULL;
      row->val = NULL;
   }
   
   row->name = name;
   row->rhs = rhs;
   row->epsilon = epsilon;
   row->index = stat->numrowidx++;
   row->size = len;
   row->len = len;
   row->numuses = 0;
   row->lppos = -1;
   row->equality = equality;

   /* calculate row norms and min/maxidx, and check if row is sorted */
   calcRowNorms(row, set);

   /* add coefficients to columns */
   for( i = 0; i < len; ++i )
      addColCoeff(mem, set, row->col[i], row, row->val[i]);
   
   return row;
}

void SCIProwFree(                       /**< frees an LP row */
   ROW**            row,                /**< pointer to LP row */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   )
{
   int i;

   assert(mem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->numuses == 0);
   
   /* delete coefficients in corresponding columns */
   for( i = 0; i < (*row)->len; ++i )
      deleteColCoeff((*row)->col[i], set, *row);

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
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(mem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->numuses >= 1);

   (*row)->numuses--;
   if( (*row)->numuses == 0 )
      SCIProwFree(row, mem, set);
}

COL** SCIPlpGetNewcols(                 /**< get array with newly added columns */
   const LP*        lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewcol && lp->firstnewcol <= lp->ncols);

   return &(lp->cols[lp->firstnewcol]);
}

int SCIPlpGetNumNewcols(                /**< get number of newly added columns */
   const LP*        lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewcol && lp->firstnewcol <= lp->ncols);

   return lp->ncols - lp->firstnewcol;
}

ROW** SCIPlpGetNewrows(                 /**< get array with newly added rows */
   const LP*        lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewrow && lp->firstnewrow <= lp->nrows);

   return &(lp->rows[lp->firstnewrow]);
}

int SCIPlpGetNumNewrows(                /**< get number of newly added rows */
   const LP*        lp                  /**< actual LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewrow && lp->firstnewrow <= lp->nrows);

   return lp->nrows - lp->firstnewrow;
}

RETCODE SCIPlpGetState(                 /**< stores LP state (like basis information) into LP state object */
   LP*              lp,                 /**< LP data */
   MEM*             mem,                /**< block memory buffers */
   LPSTATE**        lpstate             /**< pointer to LP state information (like basis information) */
   )
{
   assert(lp != NULL);
   assert(mem != NULL);
   assert(lpstate != NULL);

   return SCIPlpiGetState(lp->lpi, mem, lpstate);
}

RETCODE SCIPlpSetState(                 /**< loads LP state (like basis information) into solver */
   LP*              lp,                 /**< LP data */
   MEM*             mem,                /**< block memory buffers */
   LPSTATE*         lpstate             /**< LP state information (like basis information) */
   )
{
   assert(lp != NULL);
   assert(mem != NULL);
   assert(lpstate != NULL);

   return SCIPlpiSetState(lp->lpi, mem, lpstate);
}

RETCODE SCIPlpAddCol(                   /**< adds a column to the LP */
   LP*              lp,                 /**< LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   COL*             col                 /**< LP column */
   )
{
   assert(lp != NULL);
   assert(col != NULL);
   assert(col->lppos == -1);
   
   CHECK_OKAY( ensureAddColsSize(lp, mem, set, lp->naddcols+1) );
   lp->addcols[lp->naddcols] = col;
   lp->naddcols++;
   lp->addcolscoefs += col->len;
   lp->flushed = FALSE;

   return SCIP_OKAY;
}

RETCODE SCIPlpAddRow(                   /**< adds a row to the LP */
   LP*              lp,                 /**< LP data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   ROW*             row                 /**< LP row */
   )
{
   assert(lp != NULL);
   assert(row != NULL);
   assert(row->lppos == -1);
   
   CHECK_OKAY( ensureAddRowsSize(lp, mem, set, lp->naddrows+1) );
   lp->addrows[lp->naddrows] = row;
   lp->naddrows++;
   lp->addrowscoefs += row->len;
   lp->flushed = FALSE;

   return SCIP_OKAY;
}
