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

/**@file   sepa.c
 * @brief  methods and datastructures for separating cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "set.h"
#include "mem.h"
#include "prob.h"
#include "lp.h"
#include "stat.h"
#include "var.h"
#include "scip.h"
#include "cons_linear.h"
#include "sepa.h"


/** storage for separated cuts */
struct Sepa
{
   ROW**            cuts;               /**< array with separated cuts with violated slacks sorted by score */
   Real*            score;              /**< score for each separated cut (e.g. |slack|/(eucnorm * #nonzeros)) */
   Bool*            pool;               /**< should the cut be used in the global cut pool? Cut must be global valid! */
   int              cutssize;           /**< size of cuts, score, and pool arrays */
   int              ncuts;              /**< number of separated cuts (max. is set->maxsepacuts) */
};


/*
 * dynamic memory arrays
 */

static
RETCODE sepaEnsureCutsMem(              /**< resizes cuts and score arrays to be able to store at least num entries */
   SEPA*            sepa,               /**< separation storage */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   if( num > sepa->cutssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(sepa->cuts, newsize) );
      ALLOC_OKAY( reallocMemoryArray(sepa->score, newsize) );
      ALLOC_OKAY( reallocMemoryArray(sepa->pool, newsize) );
      sepa->cutssize = newsize;
   }
   assert(num <= sepa->cutssize);

   return SCIP_OKAY;
}




RETCODE SCIPsepaCreate(                 /**< creates separation storage */
   SEPA**           sepa                /**< pointer to store separation storage */
   )
{
   assert(sepa != NULL);
   
   ALLOC_OKAY( allocMemory(*sepa) );
   
   (*sepa)->cuts = NULL;
   (*sepa)->score = NULL;
   (*sepa)->pool = NULL;
   (*sepa)->cutssize = 0;
   (*sepa)->ncuts = 0;

   return SCIP_OKAY;
}

RETCODE SCIPsepaFree(                   /**< frees separation storage */
   SEPA**           sepa                /**< pointer to store separation storage */
   )
{
   assert(sepa != NULL);

   freeMemoryArrayNull((*sepa)->cuts);
   freeMemoryArrayNull((*sepa)->score);
   freeMemoryArrayNull((*sepa)->pool);
   freeMemory(*sepa);

   return SCIP_OKAY;
}

RETCODE SCIPsepaAddCut(                 /**< adds cut to separation storage and captures it */
   SEPA*            sepa,               /**< separation storage */
   const SET*       set,                /**< global SCIP settings */
   ROW*             cut,                /**< separated cut */
   Real             score,              /**< separation score of cut (the larger, the better the cut) */
   Bool             pool                /**< should the cut be used in the global cut pool? Cut must be global valid! */
   )
{
   int c;
   int i;

   assert(sepa != NULL);
   assert(set != NULL);
   assert(cut != NULL);

   /* get enough memory to store the cut */
   CHECK_OKAY( sepaEnsureCutsMem(sepa, set, sepa->ncuts+1) );
   assert(sepa->ncuts <= sepa->cutssize);

   /* capture cut */
   SCIProwCapture(cut);

   /* check, if cut belongs to the best "maxsepacuts" separation cuts */
   if( sepa->ncuts < set->maxsepacuts || score > sepa->score[set->maxsepacuts-1] )
   {
      /* search the correct position of the cut in the cuts array */
      for( c = 0; c < sepa->ncuts && score <= sepa->score[c]; ++c )
      {
      }
      assert(c < set->maxsepacuts);

      /* if the array consists of at least "maxsepacuts" cuts, move the worst of the best "maxsepacuts" cuts
       * to the end of the array
       */
      if( sepa->ncuts >= set->maxsepacuts )
      {
         sepa->cuts[sepa->ncuts] = sepa->cuts[set->maxsepacuts-1];
         sepa->score[sepa->ncuts] = sepa->score[set->maxsepacuts-1];
         sepa->pool[sepa->ncuts] = sepa->pool[set->maxsepacuts-1];
      }

      /* insert cut in the sorted arrays */
      for( i = MIN(sepa->ncuts, set->maxsepacuts-1); i > c; --i )
      {
         sepa->cuts[i] = sepa->cuts[i-1];
         sepa->score[i] = sepa->score[i-1];
         sepa->pool[i] = sepa->pool[i-1];
      }
      sepa->cuts[c] = cut;
      sepa->score[c] = score;
      sepa->pool[c] = pool;
   }
   else
   {
      /* cut is not under the best "maxsepacuts" cuts: put it to the end of the arrays */
      sepa->cuts[sepa->ncuts] = cut;
      sepa->score[sepa->ncuts] = score;
      sepa->pool[sepa->ncuts] = pool;
   }
    
   sepa->ncuts++;

   return SCIP_OKAY;
}

RETCODE SCIPsepaApplyCuts(              /**< adds cuts to the LP and clears separation storage */
   SEPA*            sepa,               /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< LP data */
   )
{
   int c;

   assert(sepa != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   for( c = 0; c < MIN(sepa->ncuts, set->maxsepacuts); ++c )
   {
      debugMessage("apply cut: ");
      debug( SCIProwPrint(sepa->cuts[c], set, NULL) );

      /* add cut to the LP and capture it */
      CHECK_OKAY( SCIPlpAddRow(lp, set, sepa->cuts[c]) );

      if( sepa->pool[c] )
      {
         CONS* cons;

         /* add the cut to the global cut pool and capture it: create a global linear constraint */
         debugMessage("  -> put it as global linear constraint in the cut pool\n");
         CHECK_OKAY( SCIPcreateConsLPRow(set->scip, &cons, sepa->cuts[c], FALSE) );
         CHECK_OKAY( SCIPtreeAddGlobalCons(tree, memhdr, set, cons) );
      }
   }

   /* release all cuts of the pricing storage and clear the storage */
   for( c = 0; c < sepa->ncuts; ++c )
      SCIProwRelease(&sepa->cuts[c], memhdr, set, lp);
   sepa->ncuts = 0;

   return SCIP_OKAY;
}

int SCIPsepaGetNCuts(                   /**< get number of cuts in the separation storage */
   SEPA*            sepa                /**< separation storage */
   )
{
   assert(sepa != NULL);

   return sepa->ncuts;
}
