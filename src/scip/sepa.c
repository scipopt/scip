/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
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

#include "sepa.h"
#include "prob.h"
#include "stat.h"
#include "var.h"
#include "scip.h"
#include "cons_linear.h"


/** storage for separated cuts */
struct Sepa
{
   ROW**            cuts;               /**< array with separated cuts sorted by score */
   Real*            score;              /**< score for each separated cut (e.g. violation/(eucnorm * #nonzeros)) */
   int              cutssize;           /**< size of cuts and score arrays */
   int              ncuts;              /**< number of separated cuts (max. is set->maxsepacuts) */
};


/*
 * dynamic memory arrays
 */

/** resizes cuts and score arrays to be able to store at least num entries */
static
RETCODE sepaEnsureCutsMem(
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
      ALLOC_OKAY( reallocMemoryArray(&sepa->cuts, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&sepa->score, newsize) );
      sepa->cutssize = newsize;
   }
   assert(num <= sepa->cutssize);

   return SCIP_OKAY;
}




/** creates separation storage */
RETCODE SCIPsepaCreate(
   SEPA**           sepa                /**< pointer to store separation storage */
   )
{
   assert(sepa != NULL);
   
   ALLOC_OKAY( allocMemory(sepa) );
   
   (*sepa)->cuts = NULL;
   (*sepa)->score = NULL;
   (*sepa)->cutssize = 0;
   (*sepa)->ncuts = 0;

   return SCIP_OKAY;
}

/** frees separation storage */
RETCODE SCIPsepaFree(
   SEPA**           sepa                /**< pointer to store separation storage */
   )
{
   assert(sepa != NULL);

   freeMemoryArrayNull(&(*sepa)->cuts);
   freeMemoryArrayNull(&(*sepa)->score);
   freeMemory(sepa);

   return SCIP_OKAY;
}

/** adds cut to separation storage and captures it */
RETCODE SCIPsepaAddCut(
   SEPA*            sepa,               /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   ROW*             cut,                /**< separated cut */
   Real             score,              /**< separation score of cut (the larger, the better the cut) */
   Bool             root                /**< are we at the root node? */
   )
{
   int maxsepacuts;
   int c;
   int i;

   assert(sepa != NULL);
   assert(set != NULL);
   assert(cut != NULL);

   maxsepacuts = SCIPsetGetMaxsepacuts(set, root);
   assert(sepa->ncuts <= maxsepacuts);
   if( maxsepacuts == 0 )
      return SCIP_OKAY;

   /* get enough memory to store the cut */
   CHECK_OKAY( sepaEnsureCutsMem(sepa, set, sepa->ncuts+1) );
   assert(sepa->ncuts <= sepa->cutssize);

   /* check, if cut belongs to the best "maxsepacuts" separation cuts */
   if( sepa->ncuts < maxsepacuts || score > sepa->score[maxsepacuts-1] )
   {
      debugMessage("adding cut to separation storage of size %d/%d\n", sepa->ncuts, maxsepacuts);

      /* capture the cut */
      SCIProwCapture(cut);

      /* search the correct position of the cut in the cuts array */
      for( c = 0; c < sepa->ncuts && score <= sepa->score[c]; ++c )
      {
      }
      assert(c <= sepa->ncuts);
      assert(c < maxsepacuts);

      /* if the array consists of "maxsepacuts" cuts, release the worst cut */
      if( sepa->ncuts == maxsepacuts )
      {
         CHECK_OKAY( SCIProwRelease(&sepa->cuts[sepa->ncuts-1], memhdr, set, lp) );
         sepa->ncuts--;
      }
      assert(sepa->ncuts < maxsepacuts);

      /* insert cut in the sorted arrays */
      for( i = sepa->ncuts; i > c; --i )
      {
         sepa->cuts[i] = sepa->cuts[i-1];
         sepa->score[i] = sepa->score[i-1];
      }
      sepa->cuts[c] = cut;
      sepa->score[c] = score;
      sepa->ncuts++;
   }

   return SCIP_OKAY;
}

/** adds cuts to the LP and clears separation storage */
RETCODE SCIPsepaApplyCuts(
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

   debugMessage("applying %d cuts\n", sepa->ncuts);

   for( c = 0; c < sepa->ncuts; ++c )
   {
      debugMessage("apply cut: ");
      debug( SCIProwPrint(sepa->cuts[c], NULL) );

      /* add cut to the LP and capture it */
      CHECK_OKAY( SCIPlpAddRow(lp, set, sepa->cuts[c]) );

      /* release the row */
      CHECK_OKAY( SCIProwRelease(&sepa->cuts[c], memhdr, set, lp) );
   }

   /* clear the separation storage */
   sepa->ncuts = 0;

   return SCIP_OKAY;
}

/** clears the separation storage without adding the cuts to the LP */
RETCODE SCIPsepaClearCuts(
   SEPA*            sepa,               /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< LP data */
   )
{
   int c;

   assert(sepa != NULL);

   debugMessage("clearing %d cuts\n", sepa->ncuts);

   for( c = 0; c < sepa->ncuts; ++c )
   {
      /* release the row */
      CHECK_OKAY( SCIProwRelease(&sepa->cuts[c], memhdr, set, lp) );
   }

   /* clear the separation storage */
   sepa->ncuts = 0;

   return SCIP_OKAY;
}

/** get number of cuts in the separation storage */
int SCIPsepaGetNCuts(
   SEPA*            sepa                /**< separation storage */
   )
{
   assert(sepa != NULL);

   return sepa->ncuts;
}
