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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cutpool.c,v 1.42 2005/02/14 13:35:42 bzfpfend Exp $"

/**@file   cutpool.c
 * @brief  methods for storing cuts in a cut pool
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/misc.h"
#include "scip/lp.h"
#include "scip/sepastore.h"
#include "scip/cutpool.h"

#include "scip/struct_cutpool.h"



/*
 * Hash functions
 */

/** gets the hash key of a cut */
static
DECL_HASHGETKEY(hashGetKeyCut)
{  /*lint --e{715}*/
   CUT* cut;

   cut = (CUT*)elem;
   assert(cut != NULL);
   assert(cut->row != NULL);

   /* the key of a cut is the row */
   return cut->row;
}

/** returns TRUE iff both cuts are identical */
static
DECL_HASHKEYEQ(hashKeyEqCut)
{  /*lint --e{715}*/
   /* Warning: The comparison of real values is made against default epsilon.
    *          This is ugly, but we have no settings at hand.
    */

   int i;

   ROW* row1;
   ROW* row2;

   row1 = (ROW*)key1;
   row2 = (ROW*)key2;
   assert(row1 != NULL);
   assert(row2 != NULL);

   /* sort the column indices of a row */
   SCIProwSort(row1);
   SCIProwSort(row2);
   assert(row1->lpcolssorted);
   assert(row1->nonlpcolssorted);
   assert(row1->validminmaxidx);
   assert(row2->lpcolssorted);
   assert(row2->nonlpcolssorted);
   assert(row2->validminmaxidx);

   /* compare the trivial characteristics of the rows */
   if( row1->len != row2->len
      || row1->minidx != row2->minidx
      || row1->maxidx != row2->maxidx
      || row1->nummaxval != row2->nummaxval
      || REALABS(row1->lhs - row2->lhs) > SCIP_DEFAULT_EPSILON
      || REALABS(row1->rhs - row2->rhs) > SCIP_DEFAULT_EPSILON
      || REALABS(row1->sqrnorm - row2->sqrnorm) > SCIP_DEFAULT_SUMEPSILON
      || REALABS(row1->sumnorm - row2->sumnorm) > SCIP_DEFAULT_SUMEPSILON
      || REALABS(row1->maxval - row2->maxval) > SCIP_DEFAULT_EPSILON
       )
      return FALSE;

   /* compare the columns of the rows */
   for( i = 0; i < row1->len; ++i )
   {
      if( row1->cols[i] != row2->cols[i] )
         return FALSE;
   }

   /* compare the coefficients of the rows */
   for( i = 0; i < row1->len; ++i )
   {
      if( REALABS(row1->vals[i] - row2->vals[i]) > SCIP_DEFAULT_EPSILON )
         return FALSE;
   }

   return TRUE;
}

static
DECL_HASHKEYVAL(hashKeyValCut)
{  /*lint --e{715}*/
   ROW* row;
   unsigned int keyval;

   row = (ROW*)key;
   assert(row != NULL);

   keyval = (row->nummaxval << 29) + (row->len << 22) + (row->minidx << 11) + row->maxidx; /*lint !e701*/

   return keyval;
}



/*
 * dynamic memory arrays
 */

/** resizes cuts array to be able to store at least num entries */
static
RETCODE cutpoolEnsureCutsMem(
   CUTPOOL*         cutpool,            /**< cut pool */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(cutpool != NULL);
   assert(set != NULL);

   if( num > cutpool->cutssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&cutpool->cuts, newsize) );
      cutpool->cutssize = newsize;
   }
   assert(num <= cutpool->cutssize);

   return SCIP_OKAY;
}



/*
 * Cut methods
 */

/** creates a cut and captures the row */
static
RETCODE cutCreate(
   CUT**            cut,                /**< pointer to store the cut */
   BLKMEM*          blkmem,             /**< block memory */
   ROW*             row                 /**< row this cut represents */
   )
{
   assert(cut != NULL);
   assert(blkmem != NULL);
   assert(row != NULL);

   /* allocate cut memory */
   ALLOC_OKAY( allocBlockMemory(blkmem, cut) );
   (*cut)->row = row;
   (*cut)->age = 0;
   (*cut)->processedlp = -1;
   (*cut)->pos = -1;

   /* capture row */
   SCIProwCapture(row);

   return SCIP_OKAY;
}

/** frees a cut and releases the row */
static
RETCODE cutFree(
   CUT**            cut,                /**< pointer to store the cut */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   assert(cut != NULL);
   assert(*cut != NULL);
   assert((*cut)->row != NULL);
   assert(blkmem != NULL);
   
   /* release row */
   CHECK_OKAY( SCIProwRelease(&(*cut)->row, blkmem, set, lp) );

   /* free cut memory */
   freeBlockMemory(blkmem, cut);

   return SCIP_OKAY;
}

/** returns whether the cut's age exceeds the age limit */
static
Bool cutIsAged(
   CUT*             cut,                /**< cut to check */
   int              agelimit            /**< maximum age a cut can reach before it is deleted from the pool, or -1 */
   )
{
   assert(cut != NULL);

   return (agelimit >= 0 && cut->age > agelimit);
}




/*
 * Cutpool methods
 */

/** creates cut pool */
RETCODE SCIPcutpoolCreate(
   CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   int              agelimit            /**< maximum age a cut can reach before it is deleted from the pool */
   )
{
   assert(cutpool != NULL);
   assert(agelimit >= -1);

   ALLOC_OKAY( allocMemory(cutpool) );

   CHECK_OKAY( SCIPclockCreate(&(*cutpool)->clock, SCIP_CLOCKTYPE_DEFAULT) );

   CHECK_OKAY( SCIPhashtableCreate(&(*cutpool)->hashtable, blkmem, SCIP_HASHSIZE_CUTPOOLS,
                  hashGetKeyCut, hashKeyEqCut, hashKeyValCut) );

   (*cutpool)->cuts = NULL;
   (*cutpool)->cutssize = 0;
   (*cutpool)->ncuts = 0;
   (*cutpool)->agelimit = agelimit;
   (*cutpool)->processedlp = -1;
   (*cutpool)->firstunprocessed = 0;
   (*cutpool)->maxncuts = 0;
   (*cutpool)->ncalls = 0;
   (*cutpool)->ncutsfound = 0;

   return SCIP_OKAY;
}

/** frees cut pool */
RETCODE SCIPcutpoolFree(
   CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   assert(cutpool != NULL);
   assert(*cutpool != NULL);

   /* remove all cuts from the pool */
   CHECK_OKAY( SCIPcutpoolClear(*cutpool, blkmem, set, lp) );

   /* free clock */
   SCIPclockFree(&(*cutpool)->clock);

   /* free hash table */
   SCIPhashtableFree(&(*cutpool)->hashtable);
   
   freeMemoryArrayNull(&(*cutpool)->cuts);
   freeMemory(cutpool);

   return SCIP_OKAY;
}

/** removes all rows from the cut pool */
RETCODE SCIPcutpoolClear(
   CUTPOOL*         cutpool,            /**< cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(cutpool != NULL);

   /* free cuts */
   for( i = 0; i < cutpool->ncuts; ++i )
   {
      SCIProwUnlock(cutpool->cuts[i]->row);
      CHECK_OKAY( cutFree(&cutpool->cuts[i], blkmem, set, lp) );
   }
   cutpool->ncuts = 0;

   return SCIP_OKAY;
}

/** if not already existing, adds row to cut pool and captures it */
RETCODE SCIPcutpoolAddRow(
   CUTPOOL*         cutpool,            /**< cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   ROW*             row                 /**< cutting plane to add */
   )
{
   assert(cutpool != NULL);
   assert(row != NULL);

   /* check in hash table, if cut already exists in the pool */
   if( SCIPhashtableRetrieve(cutpool->hashtable, (void*)row) == NULL )
   {
      CHECK_OKAY( SCIPcutpoolAddNewRow(cutpool, blkmem, set, row) );
   }

   return SCIP_OKAY;
}

/** adds row to cut pool and captures it; doesn't check for multiple cuts */
RETCODE SCIPcutpoolAddNewRow(
   CUTPOOL*         cutpool,            /**< cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   ROW*             row                 /**< cutting plane to add */
   )
{
   CUT* cut;

   assert(cutpool != NULL);
   assert(row != NULL);

   /* check, if row is modifiable or local */
   if( SCIProwIsModifiable(row) )
   {
      errorMessage("cannot store modifiable row <%s> in a cut pool\n", SCIProwGetName(row));
      return SCIP_INVALIDDATA;
   }
   if( SCIProwIsLocal(row) )
   {
      errorMessage("cannot store locally valid row <%s> in a cut pool\n", SCIProwGetName(row));
      return SCIP_INVALIDDATA;
   }

   /* create the cut */
   CHECK_OKAY( cutCreate(&cut, blkmem, row) );
   cut->pos = cutpool->ncuts;

   /* add cut to the pool */
   CHECK_OKAY( cutpoolEnsureCutsMem(cutpool, set, cutpool->ncuts+1) );
   cutpool->cuts[cutpool->ncuts] = cut;
   cutpool->ncuts++;
   cutpool->maxncuts = MAX(cutpool->maxncuts, cutpool->ncuts);

   /* insert cut in the hash table */
   CHECK_OKAY( SCIPhashtableInsert(cutpool->hashtable, (void*)cut) );

   /* lock the row */
   SCIProwLock(row);

   return SCIP_OKAY;
}

/** removes the cut from the cut pool */
static
RETCODE cutpoolDelCut(
   CUTPOOL*         cutpool,            /**< cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< current LP data */
   CUT*             cut                 /**< cut to remove */
   )
{
   int pos;

   assert(cutpool != NULL);
   assert(cutpool->firstunprocessed <= cutpool->ncuts);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(cutpool->processedlp <= stat->lpcount);
   assert(cut != NULL);
   assert(cut->row != NULL);

   pos = cut->pos;
   assert(0 <= pos && pos < cutpool->ncuts);
   assert(cutpool->cuts[pos] == cut);

   /* unlock the row */
   SCIProwUnlock(cut->row);

   /* remove the cut from the hash table */
   CHECK_OKAY( SCIPhashtableRemove(cutpool->hashtable, (void*)cut) );

   /* free the cut */
   CHECK_OKAY( cutFree(&cutpool->cuts[pos], blkmem, set, lp) );
   
   /* move the last cut of the pool to the free position */
   if( pos < cutpool->ncuts-1 )
   {
      cutpool->cuts[pos] = cutpool->cuts[cutpool->ncuts-1];
      cutpool->cuts[pos]->pos = pos;
      assert(cutpool->cuts[pos]->processedlp <= stat->lpcount);
      if( cutpool->cuts[pos]->processedlp < stat->lpcount )
         cutpool->firstunprocessed = MIN(cutpool->firstunprocessed, pos);
   }
   else
      cutpool->firstunprocessed = MIN(cutpool->firstunprocessed, cutpool->ncuts-1);

   cutpool->ncuts--;

   return SCIP_OKAY;
}

/** removes the LP row from the cut pool */
RETCODE SCIPcutpoolDelRow(
   CUTPOOL*         cutpool,            /**< cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< current LP data */
   ROW*             row                 /**< row to remove */
   )
{
   CUT* cut;

   assert(cutpool != NULL);
   assert(row != NULL);

   /* find the cut in hash table */
   cut = (CUT*)SCIPhashtableRetrieve(cutpool->hashtable, (void*)row);
   if( cut == NULL )
   {
      errorMessage("row <%s> is not existing in cutpool %p\n", SCIProwGetName(row), cutpool);
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( cutpoolDelCut(cutpool, blkmem, set, stat, lp, cut) );

   return SCIP_OKAY;
}


/** separates cuts of the cut pool */
RETCODE SCIPcutpoolSeparate(
   CUTPOOL*         cutpool,            /**< cut pool */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< current LP data */
   SEPASTORE*       sepastore,          /**< separation storage */
   Bool             root,               /**< are we at the root node? */
   RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   CUT* cut;
   Bool found;
   int oldncutsstored;
   int c;

   assert(cutpool != NULL);
   assert(stat != NULL);
   assert(cutpool->processedlp <= stat->lpcount);
   assert(cutpool->firstunprocessed <= cutpool->ncuts);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* don't separate cut pool in the root node */
   if( root )
      return SCIP_OKAY;

   if( cutpool->processedlp < stat->lpcount )
      cutpool->firstunprocessed = 0;
   if( cutpool->firstunprocessed == cutpool->ncuts )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;
   cutpool->ncalls++;
   found = FALSE;

   debugMessage("separating cut pool %p with %d cuts, beginning with cut %d\n",
      cutpool, cutpool->ncuts, cutpool->firstunprocessed);

   /* start timing */
   SCIPclockStart(cutpool->clock, set);

   /* remember the current total number of found cuts */
   oldncutsstored = SCIPsepastoreGetNCutsStored(sepastore);

   /* process all unprocessed cuts in the pool */
   for( c = cutpool->firstunprocessed; c < cutpool->ncuts; ++c )
   {
      cut = cutpool->cuts[c];
      assert(cut != NULL);
      assert(cut->processedlp <= stat->lpcount);
      assert(cut->pos == c);

      if( cut->processedlp < stat->lpcount )
      {
         ROW* row;

         cut->processedlp = stat->lpcount;
         row = cut->row;

         if( !SCIProwIsInLP(row) )
         {         
            if( SCIProwIsEfficacious(row, set, stat, lp, root) )
            {
               /* insert cut in separation storage */
               debugMessage(" -> separated cut <%s> from the cut pool (feasibility: %g)\n",
                  SCIProwGetName(row), SCIProwGetLPFeasibility(row, stat, lp));
               CHECK_OKAY( SCIPsepastoreAddCut(sepastore, blkmem, set, stat, lp, row, FALSE, root) );
               found = TRUE;
            }
            else
            {
               cut->age++;
               if( cutIsAged(cut, cutpool->agelimit) )
               {
                  CHECK_OKAY( cutpoolDelCut(cutpool, blkmem, set, stat, lp, cut) );
               }
            }
         }
      }
   }
   
   cutpool->processedlp = stat->lpcount;
   cutpool->firstunprocessed = cutpool->ncuts;

   /* update the number of found cuts */
   cutpool->ncutsfound += SCIPsepastoreGetNCutsStored(sepastore) - oldncutsstored; /*lint !e776*/

   /* stop timing */
   SCIPclockStop(cutpool->clock, set);

   if( found )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/** get number of cuts in the cut pool */
int SCIPcutpoolGetNCuts(
   CUTPOOL*         cutpool             /**< cut pool */
   )
{
   assert(cutpool != NULL);

   return cutpool->ncuts;
}

/** get maximum number of cuts that were stored in the cut pool at the same time */
int SCIPcutpoolGetMaxNCuts(
   CUTPOOL*         cutpool             /**< cut pool */
   )
{
   assert(cutpool != NULL);

   return cutpool->maxncuts;
}

/** gets time in seconds used for separating cuts from the pool */
Real SCIPcutpoolGetTime(
   CUTPOOL*         cutpool             /**< cut pool */
   )
{
   assert(cutpool != NULL);

   return SCIPclockGetTime(cutpool->clock);
}

/** get number of times, the cut pool was separated */
Longint SCIPcutpoolGetNCalls(
   CUTPOOL*         cutpool             /**< cut pool */
   )
{
   assert(cutpool != NULL);

   return cutpool->ncalls;
}

/** get total number of cuts that were separated from the cut pool */
Longint SCIPcutpoolGetNCutsFound(
   CUTPOOL*         cutpool             /**< cut pool */
   )
{
   assert(cutpool != NULL);

   return cutpool->ncutsfound;
}

