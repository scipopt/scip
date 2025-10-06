/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepastoreexact.c
 * @brief  internal methods for storing separated exact cuts
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>

#include "scip/def.h"
#include "scip/cons.h"
#include "scip/cuts.h"
#include "scip/debug.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/lpexact.h"
#include "scip/misc.h"
#include "scip/rational.h"
#include "scip/reopt.h"
#include "scip/sepa.h"
#include "scip/sepastoreexact.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/struct_lpexact.h"
#include "scip/struct_sepastore.h"
#include "scip/tree.h"
#include "scip/var.h"

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
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   if( !set->exact_enable )
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
   SCIP_SEPASTOREEXACT** sepastoreexact      /**< pointer to store separation storage */
   )
{
   assert(sepastoreexact != NULL);
   assert(*sepastoreexact != NULL);
   assert((*sepastoreexact)->ncuts == 0);

   BMSfreeMemoryArrayNull(&(*sepastoreexact)->cuts);
   BMSfreeMemory(sepastoreexact);

   return SCIP_OKAY;
}

#ifdef SCIP_DISABLED_CODE
/** informs separation storage that the setup of the initial LP starts now */
void SCIPsepastoreExactStartInitialLP(
   SCIP_SEPASTOREEXACT*  sepastoreexact      /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);
   assert(!sepastoreexact->initiallp);
   assert(sepastoreexact->ncuts == 0);

   sepastoreexact->initiallp = TRUE;
}

/** informs separation storage that the setup of the initial LP is now finished */
void SCIPsepastoreExactEndInitialLP(
   SCIP_SEPASTOREEXACT*  sepastoreexact      /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);
   assert(sepastoreexact->initiallp);
   assert(sepastoreexact->ncuts == 0);

   sepastoreexact->initiallp = FALSE;
}
#endif

/** adds cut to separation storage and captures it */
SCIP_RETCODE SCIPsepastoreExactAddCut(
   SCIP_SEPASTOREEXACT*  sepastoreexact,     /**< separation storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_ROWEXACT*        cut                 /**< separated cut */
   )
{
   int pos;

   assert(sepastoreexact != NULL);
   assert(set != NULL);
   assert(cut != NULL);
   assert(!SCIPrationalIsNegInfinity(SCIProwExactGetLhs(cut)) || !SCIPrationalIsInfinity(SCIProwExactGetRhs(cut)));
   assert(eventqueue != NULL);

   /* debug: check cut for feasibility; note that this just checks for fp feasibility and could be extended */
   SCIP_CALL( SCIPdebugCheckRow(set, cut->fprow) ); /*lint !e506 !e774*/

   /* update statistics of total number of found cuts */
   if( !sepastoreexact->initiallp )
   {
      sepastoreexact->ncutsfound++;
      sepastoreexact->ncutsfoundround++;
   }

   /* get enough memory to store the cut */
   SCIP_CALL( sepastoreExactEnsureCutsMem(sepastoreexact, set, sepastoreexact->ncuts+1) );
   assert(sepastoreexact->ncuts < sepastoreexact->cutssize);

   SCIPsetDebugMsg(set, "adding cut <%s> to exact separation storage of size %d (len=%d)\n",
      SCIProwGetName(cut->fprow), sepastoreexact->ncuts, SCIProwGetNNonz(cut->fprow));
   /*SCIP_CALL( SCIPprintRow(set->scip, cut, NULL) );*/

   /* capture the cut */
   SCIProwExactCapture(cut);

   pos = sepastoreexact->ncuts;

   sepastoreexact->cuts[pos] = cut;
   sepastoreexact->ncuts++;

   return SCIP_OKAY;
}

/** clears the separation storage without adding the cuts to the LP */
SCIP_RETCODE SCIPsepastoreExactClearCuts(
   SCIP_SEPASTOREEXACT*  sepastoreexact,     /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lp                  /**< LP data */
   )
{
   int c;

   if( !set->exact_enable )
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
   SCIP_SEPASTOREEXACT*  sepastoreexact      /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);

   return sepastoreexact->cuts;
}

/** get number of cuts in the separation storage */
int SCIPsepastoreExactGetNCuts(
   SCIP_SEPASTOREEXACT*  sepastoreexact      /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);

   return sepastoreexact->ncuts;
}

/** get total number of cuts found so far */
int SCIPsepastoreExactGetNCutsFound(
   SCIP_SEPASTOREEXACT*  sepastoreexact      /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);

   return sepastoreexact->ncutsfound;
}

/** get number of cuts found so far in current separation round */
int SCIPsepastoreExactGetNCutsFoundRound(
   SCIP_SEPASTOREEXACT*  sepastoreexact      /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);

   return sepastoreexact->ncutsfoundround;
}

/** get total number of cuts applied to the LPs */
int SCIPsepastoreExactGetNCutsApplied(
   SCIP_SEPASTOREEXACT*  sepastoreexact      /**< separation storage */
   )
{
   assert(sepastoreexact != NULL);

   return sepastoreexact->ncutsapplied;
}
