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

/**@file   sepastore.c
 * @brief  methods and datastructures for storing separated cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "sepastore.h"
#include "prob.h"
#include "stat.h"
#include "var.h"
#include "scip.h"
#include "cons_linear.h"


/** storage for separated cuts */
struct Sepastore
{
   ROW**            cuts;               /**< array with separated cuts sorted by score */
   Real*            score;              /**< score for each separated cut (e.g. violation/(eucnorm * #nonzeros)) */
   int              cutssize;           /**< size of cuts and score arrays */
   int              ncuts;              /**< number of separated cuts (max. is set->maxsepacuts) */
   int              ncutsfound;         /**< total number of cuts found so far */
};


/*
 * dynamic memory arrays
 */

/** resizes cuts and score arrays to be able to store at least num entries */
static
RETCODE sepastoreEnsureCutsMem(
   SEPASTORE*       sepastore,          /**< separation storage */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(sepastore != NULL);
   assert(set != NULL);

   if( num > sepastore->cutssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&sepastore->cuts, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&sepastore->score, newsize) );
      sepastore->cutssize = newsize;
   }
   assert(num <= sepastore->cutssize);

   return SCIP_OKAY;
}




/** creates separation storage */
RETCODE SCIPsepastoreCreate(
   SEPASTORE**           sepastore                /**< pointer to store separation storage */
   )
{
   assert(sepastore != NULL);
   
   ALLOC_OKAY( allocMemory(sepastore) );
   
   (*sepastore)->cuts = NULL;
   (*sepastore)->score = NULL;
   (*sepastore)->cutssize = 0;
   (*sepastore)->ncuts = 0;
   (*sepastore)->ncutsfound = 0;

   return SCIP_OKAY;
}

/** frees separation storage */
RETCODE SCIPsepastoreFree(
   SEPASTORE**           sepastore                /**< pointer to store separation storage */
   )
{
   assert(sepastore != NULL);

   freeMemoryArrayNull(&(*sepastore)->cuts);
   freeMemoryArrayNull(&(*sepastore)->score);
   freeMemory(sepastore);

   return SCIP_OKAY;
}

/** adds cut to separation storage and captures it */
RETCODE SCIPsepastoreAddCut(
   SEPASTORE*       sepastore,          /**< separation storage */
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

   assert(sepastore != NULL);
   assert(set != NULL);
   assert(cut != NULL);

   /* update statistics of total number of found cuts */
   sepastore->ncutsfound++;

   /* get maximum of separated cuts at this node */
   maxsepacuts = SCIPsetGetMaxsepacuts(set, root);
   assert(sepastore->ncuts <= maxsepacuts);
   if( maxsepacuts == 0 )
      return SCIP_OKAY;

   /* get enough memory to store the cut */
   CHECK_OKAY( sepastoreEnsureCutsMem(sepastore, set, sepastore->ncuts+1) );
   assert(sepastore->ncuts <= sepastore->cutssize);

   /* check, if cut belongs to the best "maxsepacuts" separation cuts */
   if( sepastore->ncuts < maxsepacuts || score > sepastore->score[maxsepacuts-1] )
   {
      debugMessage("adding cut to separation storage of size %d/%d\n", sepastore->ncuts, maxsepacuts);

      /* capture the cut */
      SCIProwCapture(cut);

      /* search the correct position of the cut in the cuts array */
      for( c = 0; c < sepastore->ncuts && score <= sepastore->score[c]; ++c )
      {
      }
      assert(c <= sepastore->ncuts);
      assert(c < maxsepacuts);

      /* if the array consists of "maxsepacuts" cuts, release the worst cut */
      if( sepastore->ncuts == maxsepacuts )
      {
         CHECK_OKAY( SCIProwRelease(&sepastore->cuts[sepastore->ncuts-1], memhdr, set, lp) );
         sepastore->ncuts--;
      }
      assert(sepastore->ncuts < maxsepacuts);

      /* insert cut in the sorted arrays */
      for( i = sepastore->ncuts; i > c; --i )
      {
         sepastore->cuts[i] = sepastore->cuts[i-1];
         sepastore->score[i] = sepastore->score[i-1];
      }
      sepastore->cuts[c] = cut;
      sepastore->score[c] = score;
      sepastore->ncuts++;
   }

   return SCIP_OKAY;
}

/** adds cuts to the LP and clears separation storage */
RETCODE SCIPsepastoreApplyCuts(
   SEPASTORE*       sepastore,          /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< LP data */
   )
{
   int c;

   assert(sepastore != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   debugMessage("applying %d cuts\n", sepastore->ncuts);

   for( c = 0; c < sepastore->ncuts; ++c )
   {
      debugMessage("apply cut: ");
      debug( SCIProwPrint(sepastore->cuts[c], NULL) );

      /* add cut to the LP and capture it */
      CHECK_OKAY( SCIPlpAddRow(lp, set, sepastore->cuts[c]) );

      /* release the row */
      CHECK_OKAY( SCIProwRelease(&sepastore->cuts[c], memhdr, set, lp) );
   }

   /* clear the separation storage */
   sepastore->ncuts = 0;

   return SCIP_OKAY;
}

/** clears the separation storage without adding the cuts to the LP */
RETCODE SCIPsepastoreClearCuts(
   SEPASTORE*       sepastore,          /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< LP data */
   )
{
   int c;

   assert(sepastore != NULL);

   debugMessage("clearing %d cuts\n", sepastore->ncuts);

   for( c = 0; c < sepastore->ncuts; ++c )
   {
      /* release the row */
      CHECK_OKAY( SCIProwRelease(&sepastore->cuts[c], memhdr, set, lp) );
   }

   /* clear the separation storage */
   sepastore->ncuts = 0;

   return SCIP_OKAY;
}

/** get number of cuts in the separation storage */
int SCIPsepastoreGetNCuts(
   SEPASTORE*            sepastore                /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncuts;
}

/** get total number of cuts found so far */
int SCIPsepastoreGetNCutsFound(
   SEPASTORE*            sepastore                /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsfound;
}
