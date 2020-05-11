/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepastoreexact.c
 * @brief  internal methods for storing separated exact cuts
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/tree.h"
#include "scip/reopt.h"
#include "scip/sepastoreexact.h"
#include "scip/event.h"
#include "scip/sepa.h"
#include "scip/cons.h"
#include "scip/debug.h"
#include "scip/scip.h"
#include "scip/cuts.h"
#include "scip/struct_sepastore.h"
#include "scip/misc.h"
#include "scip/lpexact.h"
#include "scip/rational.h"
#include "scip/pub_lpexact.h"

/** resizes cuts and score arrays to be able to store at least num entries */
static
SCIP_RETCODE sepastoreExactEnsureCutsMem(
   SCIP_SEPASTOREEXACT*  sepastoreexact,     /**< separation storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(sepastoreexact != NULL);
   assert(set != NULL);

   if( num > sepastoreexact->cutssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&sepastoreexact->cuts, newsize) );
      sepastoreexact->cutssize = newsize;
   }
   assert(num <= sepastoreexact->cutssize);

   return SCIP_OKAY;
}

/** creates separation storage */
SCIP_RETCODE SCIPsepastoreExactCreate(
   SCIP_SEPASTOREEXACT** sepastoreexact,     /**< pointer to store separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   if( !set->misc_exactsolve )
      return SCIP_OKAY;

   assert(sepastoreexact != NULL);

   SCIP_ALLOC( BMSallocMemory(sepastoreexact) );

   (*sepastoreexact)->cuts = NULL;
   (*sepastoreexact)->cutssize = 0;
   (*sepastoreexact)->ncuts = 0;
   (*sepastoreexact)->ncutsfound = 0;
   (*sepastoreexact)->ncutsfoundround = 0;
   (*sepastoreexact)->ncutsapplied = 0;

   (*sepastoreexact)->initiallp = FALSE;

   return SCIP_OKAY;
}

/** frees separation storage */
SCIP_RETCODE SCIPsepastoreExactFree(
   SCIP_SEPASTOREEXACT** sepastoreexact,     /**< pointer to store separation storage */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(sepastoreexact != NULL);
   assert(*sepastoreexact != NULL);
   assert((*sepastoreexact)->ncuts == 0);

   BMSfreeMemoryArrayNull(&(*sepastoreexact)->cuts);
   BMSfreeMemory(sepastoreexact);

   return SCIP_OKAY;
}

/** informs separation storage that the setup of the initial LP starts now */
void SCIPsepastoreExactStartInitialLP(
   SCIP_SEPASTOREEXACT*  sepastoreexact         /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);
   assert(!sepastoreexact->initiallp);
   assert(sepastoreexact->ncuts == 0);

   sepastoreexact->initiallp = TRUE;
}

/** informs separation storage that the setup of the initial LP is now finished */
void SCIPsepastoreExactEndInitialLP(
   SCIP_SEPASTOREEXACT*  sepastoreexact         /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);
   assert(sepastoreexact->initiallp);
   assert(sepastoreexact->ncuts == 0);

   sepastoreexact->initiallp = FALSE;
}

/** adds cut to separation storage and captures it */
SCIP_RETCODE SCIPsepastoreexAddCut(
   SCIP_SEPASTOREEXACT*  sepastoreexact,     /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LPEXACT*         lp,                 /**< LP data */
   SCIP_ROWEXACT*        cut,                /**< separated cut */
   SCIP_Bool*            infeasible          /**< pointer to store whether the cut is infeasible */
   )
{
   SCIP_Bool redundant;
   int pos;

   assert(sepastoreexact != NULL);
   assert(set != NULL);
   assert(cut != NULL);
   assert(!RatIsNegInfinity(SCIProwExactGetLhs(cut)) || !RatIsInfinity(SCIProwExactGetRhs(cut)));
   assert(eventqueue != NULL);
   assert(eventfilter != NULL);

   /* debug: check cut for feasibility */
   SCIP_CALL( SCIPdebugCheckRow(set, cut) ); /*lint !e506 !e774*/

   /* update statistics of total number of found cuts */
   if( !sepastoreexact->initiallp )
   {
      sepastoreexact->ncutsfound++;
      sepastoreexact->ncutsfoundround++;
   }

   /* get enough memory to store the cut */
   SCIP_CALL( sepastoreExactEnsureCutsMem(sepastoreexact, set, sepastoreexact->ncuts+1) );
   assert(sepastoreexact->ncuts < sepastoreexact->cutssize);

   SCIPsetDebugMsg(set, "adding cut <%s> to exact separation storage of size %d (forcecut=%u, len=%d)\n",
      SCIProwGetName(cut->fprow), sepastoreexact->ncuts, SCIProwGetNNonz(cut->fprow));
   /*SCIP_CALL( SCIPprintRow(set->scip, cut, NULL) );*/

   /* capture the cut */
   SCIProwExactCapture(cut);

   pos = sepastoreexact->ncuts;

   sepastoreexact->cuts[pos] = cut;
   sepastoreexact->ncuts++;

   /* If the duals need to be collected, then the infeasible flag is set to FALSE. This ensures that the LP is solved */
   if( set->lp_alwaysgetduals && sepastoreexact->initiallp )
      (*infeasible) = FALSE;

   return SCIP_OKAY;
}

/** adds cuts to the LP and clears separation storage */
SCIP_RETCODE SCIPsepastoreExactSyncLPs(
   SCIP_SEPASTOREEXACT*  sepastoreexact,     /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   )
{
   SCIP_LP* fplp;
   SCIP_ROW** fprows;
   SCIP_ROWEXACT* rowexact;
   SCIP_CONS* origcons;
   int nrowsfp;
   int nrowsex;
   int nreleases;
   int nadded;
   int i;

   if( !set->misc_exactsolve )
      return SCIP_OKAY;

   fplp = lpexact->fplp;
   nreleases = 0;
   nadded = 0;

   assert(fplp != NULL);

   fprows = SCIPlpGetRows(fplp);
   nrowsfp = SCIPlpGetNRows(fplp);
   nrowsex = SCIPlexGetNRows(lpexact);

   assert(fprows != NULL);

   /* this method should sync the fp-lp withe the exact lp */

   /* remove all rows from exact lp that are not in the floating point lp */
   for( i = nrowsex - 1; i >= 0; --i )
   {
      SCIP_ROW* fprow =lpexact->rows[i]->fprow;
      assert(fprow != NULL);

      if( !SCIProwIsInLP(fprow) )
      {
         nreleases++;
         assert(i == nrowsex - nreleases);
      }
   }
   SCIPlpExactshrinkRows(lpexact, blkmem, set, eventqueue, eventfilter, lpexact->nrows - nreleases);

   for( i = 0; i < nrowsfp; ++i )
   {
      rowexact = SCIProwGetExRow(lpexact, fplp->rows[i]);
      if( rowexact != NULL )
      {
         /* if the row is already in lp, do nothing */
         if( !SCIProwExactIsInLP(rowexact) )
         {
            /* add the exact row to the exact lp */
            SCIP_CALL( SCIPlpExactAddRow(lpexact, blkmem, set, eventqueue,
                eventfilter, rowexact, 0) );
         }
      }
      else
      {
         SCIPerrorMessage("exact cut has not been created \n");
         return SCIP_OKAY;
      }
   }

   //assert(SCIPlpExactIsSynced(lpexact, set, SCIPgetMessagehdlr(set->scip)));

   SCIP_CALL( SCIPsepastoreExactClearCuts(sepastoreexact, blkmem, set, eventqueue, eventfilter, lpexact) );

   return SCIP_OKAY;
}

/** clears the separation storage without adding the cuts to the LP */
SCIP_RETCODE SCIPsepastoreExactClearCuts(
   SCIP_SEPASTOREEXACT*  sepastoreexact,     /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LPEXACT*         lp                  /**< LP data */
   )
{
   int c;

   if( !set->misc_exactsolve )
      return SCIP_OKAY;

   assert(sepastoreexact != NULL);

   SCIPsetDebugMsg(set, "clearing %d cuts\n", sepastoreexact->ncuts);

   /* release cuts */
   for( c = 0; c < sepastoreexact->ncuts; ++c )
   {
      SCIP_CALL( SCIProwExactRelease(&sepastoreexact->cuts[c], blkmem, set, lp) );
   }

   /* reset counters */
   sepastoreexact->ncuts = 0;
   sepastoreexact->ncutsfoundround = 0;

   /* if we have just finished the initial LP construction, free the (potentially large) cuts array */
   if( sepastoreexact->initiallp )
   {
      BMSfreeMemoryArrayNull(&sepastoreexact->cuts);
      sepastoreexact->cutssize = 0;
   }

   return SCIP_OKAY;
}


/** get cuts in the separation storage */
SCIP_ROWEXACT** SCIPsepastoreExactGetCuts(
   SCIP_SEPASTOREEXACT*  sepastoreexact         /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);

   return sepastoreexact->cuts;
}

/** get number of cuts in the separation storage */
int SCIPsepastoreExactGetNCuts(
   SCIP_SEPASTOREEXACT*  sepastoreexact         /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);

   return sepastoreexact->ncuts;
}

/** get total number of cuts found so far */
int SCIPsepastoreExactGetNCutsFound(
   SCIP_SEPASTOREEXACT*  sepastoreexact         /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);

   return sepastoreexact->ncutsfound;
}

/** get number of cuts found so far in current separation round */
int SCIPsepastoreExactGetNCutsFoundRound(
   SCIP_SEPASTOREEXACT*  sepastoreexact         /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);

   return sepastoreexact->ncutsfoundround;
}

/** get total number of cuts applied to the LPs */
int SCIPsepastoreExactGetNCutsApplied(
   SCIP_SEPASTOREEXACT*  sepastoreexact         /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);

   return sepastoreexact->ncutsapplied;
}
