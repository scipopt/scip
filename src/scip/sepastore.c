/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepastore.c
 * @brief  methods for storing separated cuts
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/tree.h"
#include "scip/sepastore.h"
#include "scip/event.h"
#include "scip/sepa.h"
#include "scip/cons.h"
#include "scip/debug.h"

#include "scip/struct_sepastore.h"



/*
 * dynamic memory arrays
 */

/** resizes cuts and score arrays to be able to store at least num entries */
static
SCIP_RETCODE sepastoreEnsureCutsMem(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(sepastore != NULL);
   assert(set != NULL);

   if( num > sepastore->cutssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&sepastore->cuts, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&sepastore->efficacies, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&sepastore->objparallelisms, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&sepastore->orthogonalities, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&sepastore->scores, newsize) );
      sepastore->cutssize = newsize;
   }
   assert(num <= sepastore->cutssize);

   return SCIP_OKAY;
}




/** creates separation storage */
SCIP_RETCODE SCIPsepastoreCreate(
   SCIP_SEPASTORE**      sepastore           /**< pointer to store separation storage */
   )
{
   assert(sepastore != NULL);

   SCIP_ALLOC( BMSallocMemory(sepastore) );

   (*sepastore)->cuts = NULL;
   (*sepastore)->efficacies = NULL;
   (*sepastore)->objparallelisms = NULL;
   (*sepastore)->orthogonalities = NULL;
   (*sepastore)->scores = NULL;
   (*sepastore)->cutssize = 0;
   (*sepastore)->ncuts = 0;
   (*sepastore)->nforcedcuts = 0;
   (*sepastore)->ncutsfound = 0;
   (*sepastore)->ncutsfoundround = 0;
   (*sepastore)->ncutsapplied = 0;
   (*sepastore)->initiallp = FALSE;
   (*sepastore)->forcecuts = FALSE;

   return SCIP_OKAY;
}

/** frees separation storage */
SCIP_RETCODE SCIPsepastoreFree(
   SCIP_SEPASTORE**      sepastore           /**< pointer to store separation storage */
   )
{
   assert(sepastore != NULL);
   assert(*sepastore != NULL);
   assert((*sepastore)->ncuts == 0);

   BMSfreeMemoryArrayNull(&(*sepastore)->cuts);
   BMSfreeMemoryArrayNull(&(*sepastore)->efficacies);
   BMSfreeMemoryArrayNull(&(*sepastore)->objparallelisms);
   BMSfreeMemoryArrayNull(&(*sepastore)->orthogonalities);
   BMSfreeMemoryArrayNull(&(*sepastore)->scores);
   BMSfreeMemory(sepastore);

   return SCIP_OKAY;
}

/** informs separation storage, that the setup of the initial LP starts now */
void SCIPsepastoreStartInitialLP(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);
   assert(!sepastore->initiallp);
   assert(sepastore->ncuts == 0);

   sepastore->initiallp = TRUE;
}

/** informs separation storage, that the setup of the initial LP is now finished */
void SCIPsepastoreEndInitialLP(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);
   assert(sepastore->initiallp);
   assert(sepastore->ncuts == 0);

   sepastore->initiallp = FALSE;
}

/** informs separation storage, that the following cuts should be used in any case */
void SCIPsepastoreStartForceCuts(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);
   assert(!sepastore->forcecuts);

   sepastore->forcecuts = TRUE;
}

/** informs separation storage, that the following cuts should no longer be used in any case */
void SCIPsepastoreEndForceCuts(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);
   assert(sepastore->forcecuts);

   sepastore->forcecuts = FALSE;
}

/** checks cut for redundancy due to activity bounds */
static
SCIP_Bool sepastoreIsCutRedundant(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_ROW*             cut                 /**< separated cut */
   )
{
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(sepastore != NULL);
   assert(cut != NULL);

   /* modifiable cuts cannot be declared redundant, since we don't know all coefficients */
   if( SCIProwIsModifiable(cut) )
      return FALSE;

   /* check for activity redundancy */
   lhs = SCIProwGetLhs(cut);
   rhs = SCIProwGetRhs(cut);
   minactivity = SCIProwGetMinActivity(cut, set, stat);
   maxactivity = SCIProwGetMaxActivity(cut, set, stat);
   if( SCIPsetIsLE(set, lhs, minactivity) && SCIPsetIsLE(set, maxactivity, rhs) )
   {
      SCIPdebugMessage("ignoring activity redundant cut <%s> (sides=[%g,%g], act=[%g,%g]\n",
         SCIProwGetName(cut), lhs, rhs, minactivity, maxactivity);
      /*SCIPdebug(SCIProwPrint(cut, set->scip->messagehdlr, NULL));*/
      return TRUE;
   }

   return FALSE;
}

/** adds cut stored as LP row to separation storage and captures it;
 *  if the cut should be forced to be used, an infinite score has to be used
 */
static
SCIP_RETCODE sepastoreAddCut(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SOL*             sol,                /**< primal solution that was separated, or NULL for LP solution */
   SCIP_ROW*             cut,                /**< separated cut */
   SCIP_Bool             forcecut,           /**< should the cut be forced to enter the LP? */
   SCIP_Bool             root                /**< are we at the root node? */
   )
{
   SCIP_Real cutefficacy;
   SCIP_Real cutobjparallelism;
   SCIP_Real cutscore;
   int pos;

   assert(sepastore != NULL);
   assert(sepastore->nforcedcuts <= sepastore->ncuts);
   assert(set != NULL);
   assert(cut != NULL);
   assert(sol != NULL || !SCIProwIsInLP(cut));
   assert(!SCIPsetIsInfinity(set, -SCIProwGetLhs(cut)) || !SCIPsetIsInfinity(set, SCIProwGetRhs(cut)));
   assert(eventqueue != NULL);
   assert(eventfilter != NULL);

   /* in the root node, every local cut is a global cut, and global cuts are nicer in many ways...*/
   if( root && SCIProwIsLocal(cut) )
   {
      SCIPdebugMessage("change local flag of cut <%s> to FALSE due to addition in root node\n", SCIProwGetName(cut));

      SCIP_CALL( SCIProwChgLocal(cut, FALSE) );

      assert(!SCIProwIsLocal(cut));
   }

   /* check cut for redundancy
    * in each separation round, make sure that at least one (even redundant) cut enters the LP to avoid cycling
    */
   if( !forcecut && sepastore->ncuts > 0 && sepastoreIsCutRedundant(sepastore, set, stat, cut) )
      return SCIP_OKAY;

   /* if only one cut is currently present in the cut store, it could be redundant; in this case, it can now be removed
    * again, because now a non redundant cut enters the store
    */
   if( sepastore->ncuts == 1 && sepastoreIsCutRedundant(sepastore, set, stat, sepastore->cuts[0]) )
   {
      /* check, if the row deletions from separation storage events are tracked
       * if so, issue ROWDELETEDSEPA event
       */
      if( eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWDELETEDSEPA) != 0 )
      {
         SCIP_EVENT* event;

         SCIP_CALL( SCIPeventCreateRowDeletedSepa(&event, blkmem, sepastore->cuts[0]) );
         SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
      }
      
      SCIP_CALL( SCIProwRelease(&sepastore->cuts[0], blkmem, set, lp) );
      sepastore->ncuts = 0;
      sepastore->nforcedcuts = 0;
   }

   /* a cut is forced to enter the LP if
    *  - we construct the initial LP, or
    *  - it has infinite score factor, or
    *  - it is a bound change
    * if it is a non-forced cut and no cuts should be added, abort
    */
   forcecut = forcecut || sepastore->initiallp || sepastore->forcecuts
      || (!SCIProwIsModifiable(cut) && SCIProwGetNNonz(cut) == 1);
   if( !forcecut && SCIPsetGetSepaMaxcuts(set, root) == 0 )
      return SCIP_OKAY;

   /* get enough memory to store the cut */
   SCIP_CALL( sepastoreEnsureCutsMem(sepastore, set, sepastore->ncuts+1) );
   assert(sepastore->ncuts < sepastore->cutssize);

   if( forcecut )
   {
      cutefficacy = SCIPsetInfinity(set);
      cutscore = SCIPsetInfinity(set);
      cutobjparallelism = 1.0;
   }
   else
   {
      /* initialize values to invalid (will be initialized during cut filtering) */
      cutefficacy = SCIP_INVALID;
      cutscore = SCIP_INVALID;

      /* initialize parallelism to objective (constant throughout filtering) */
      if( set->sepa_objparalfac > 0.0 )
         cutobjparallelism = SCIProwGetObjParallelism(cut, set, lp);
      else
         cutobjparallelism = 0.0; /* no need to calculate it */
   }

   SCIPdebugMessage("adding cut <%s> to separation storage of size %d (forcecut=%u, len=%d)\n",
      SCIProwGetName(cut), sepastore->ncuts, forcecut, SCIProwGetNNonz(cut));
   /*SCIPdebug(SCIProwPrint(cut, set->scip->messagehdlr, NULL));*/

   /* capture the cut */
   SCIProwCapture(cut);

   /* add cut to arrays */
   if( forcecut )
   {
      /* make room at the beginning of the array for forced cut */
      pos = sepastore->nforcedcuts;
      sepastore->cuts[sepastore->ncuts] = sepastore->cuts[pos];
      sepastore->efficacies[sepastore->ncuts] = sepastore->efficacies[pos];
      sepastore->objparallelisms[sepastore->ncuts] = sepastore->objparallelisms[pos];
      sepastore->orthogonalities[sepastore->ncuts] = sepastore->orthogonalities[pos];
      sepastore->scores[sepastore->ncuts] = sepastore->scores[pos];
      sepastore->nforcedcuts++;
   }
   else
      pos = sepastore->ncuts;

   sepastore->cuts[pos] = cut;
   sepastore->efficacies[pos] = cutefficacy;
   sepastore->objparallelisms[pos] = cutobjparallelism;
   sepastore->orthogonalities[pos] = 1.0;
   sepastore->scores[pos] = cutscore;
   sepastore->ncuts++;

   /* check, if the row addition to separation storage events are tracked
    * if so, issue ROWADDEDSEPA event
    */
   if( eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWADDEDSEPA) != 0 )
   {
      SCIP_EVENT* event;

      SCIP_CALL( SCIPeventCreateRowAddedSepa(&event, blkmem, cut) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
   }

   return SCIP_OKAY;
}

/** removes a non-forced cut from the separation storage */
static
SCIP_RETCODE sepastoreDelCut(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< LP data */
   int                   pos                 /**< position of cut to delete */
   )
{
   assert(sepastore != NULL);
   assert(sepastore->cuts != NULL);
   assert(sepastore->nforcedcuts <= pos && pos < sepastore->ncuts);

   /* check, if the row deletions from separation storage events are tracked
    * if so, issue ROWDELETEDSEPA event
    */
   if( eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWDELETEDSEPA) != 0 )
   {
      SCIP_EVENT* event;

      SCIP_CALL( SCIPeventCreateRowDeletedSepa(&event, blkmem, sepastore->cuts[pos]) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
   }

   /* release the row */
   SCIP_CALL( SCIProwRelease(&sepastore->cuts[pos], blkmem, set, lp) );

   /* move last cut to the empty position */
   sepastore->cuts[pos] = sepastore->cuts[sepastore->ncuts-1];
   sepastore->efficacies[pos] = sepastore->efficacies[sepastore->ncuts-1];
   sepastore->objparallelisms[pos] = sepastore->objparallelisms[sepastore->ncuts-1];
   sepastore->orthogonalities[pos] = sepastore->orthogonalities[sepastore->ncuts-1];
   sepastore->scores[pos] = sepastore->scores[sepastore->ncuts-1];
   sepastore->ncuts--;

   return SCIP_OKAY;
}

/** adds cut to separation storage and captures it;
 *  if the cut should be forced to enter the LP, an infinite score has to be used
 */
SCIP_RETCODE SCIPsepastoreAddCut(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SOL*             sol,                /**< primal solution that was separated, or NULL for LP solution */
   SCIP_ROW*             cut,                /**< separated cut */
   SCIP_Bool             forcecut,           /**< should the cut be forced to enter the LP? */
   SCIP_Bool             root                /**< are we at the root node? */
   )
{
   assert(sepastore != NULL);
   assert(cut != NULL);
   assert(!SCIProwIsInLP(cut));
   assert(!SCIPsetIsInfinity(set, -SCIProwGetLhs(cut)) || !SCIPsetIsInfinity(set, SCIProwGetRhs(cut)));

   /* debug: check cut for feasibility */
   SCIP_CALL( SCIPdebugCheckRow(set, cut) ); /*lint !e506 !e774*/

   /* update statistics of total number of found cuts */
   if( !sepastore->initiallp )
   {
      sepastore->ncutsfound++;
      sepastore->ncutsfoundround++;
   }

   /* add LP row cut to separation storage */
   SCIP_CALL( sepastoreAddCut(sepastore, blkmem, set, stat, eventqueue, eventfilter, lp, sol, cut, forcecut, root) );

   return SCIP_OKAY;
}

/** applies a lower bound change */
static
SCIP_RETCODE sepastoreApplyLb(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             bound,              /**< new lower bound of variable */
   SCIP_Bool             local,              /**< is it a local bound change? (otherwise global) */
   SCIP_Bool*            cutoff              /**< pointer to store TRUE, if an infeasibility has been detected */
   )
{
   assert(sepastore != NULL);
   assert(cutoff != NULL);

   if( local )
   {
      /* apply the local bound change or detect a cutoff */
      if( SCIPsetIsGT(set, bound, SCIPvarGetLbLocal(var)) )
      {
         SCIPdebugMessage(" -> applying bound change: <%s>: [%g,%g] -> [%g,%g]\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), bound, SCIPvarGetUbLocal(var));

         if( SCIPsetIsFeasLE(set, bound, SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(tree), blkmem, set, stat, prob, tree, lp, branchcand,
                  eventqueue, var, bound, SCIP_BOUNDTYPE_LOWER, FALSE) );
         }
         else
            *cutoff = TRUE;

         /* count the bound change as applied cut */
         if( !sepastore->initiallp )
            sepastore->ncutsapplied++;
      }
   }
   else
   {
      /* apply the global bound change or detect a global cutoff which means we can cutoff the root node */
      if( SCIPsetIsGT(set, bound, SCIPvarGetLbGlobal(var)) )
      {
         SCIPdebugMessage(" -> applying global bound change: <%s>: [%g,%g] -> [%g,%g]\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), bound, SCIPvarGetUbGlobal(var));

         if( SCIPsetIsFeasLE(set, bound, SCIPvarGetUbGlobal(var)) )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, prob, tree, lp, branchcand,
                  eventqueue, var, bound, SCIP_BOUNDTYPE_LOWER, FALSE) );
         }
         else
         {
            /* we are done with solving since a global bound change is infeasible */
            SCIPnodeCutoff(SCIPtreeGetRootNode(tree), set, stat, tree);
            *cutoff = TRUE;
         }

         /* count the bound change as applied cut */
         if( !sepastore->initiallp )
            sepastore->ncutsapplied++;
      }
   }

   return SCIP_OKAY;
}

/** applies an upper bound change */
static
SCIP_RETCODE sepastoreApplyUb(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             bound,              /**< new upper bound of variable */
   SCIP_Bool             local,              /**< is it a local bound change? (otherwise global) */
   SCIP_Bool*            cutoff              /**< pointer to store TRUE, if an infeasibility has been detected */
   )
{
   assert(sepastore != NULL);
   assert(cutoff != NULL);

   if( local )
   {
      /* apply the local bound change or detect a cutoff */
      if( SCIPsetIsLT(set, bound, SCIPvarGetUbLocal(var)) )
      {
         SCIPdebugMessage(" -> applying bound change: <%s>: [%g,%g] -> [%g,%g]\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var), bound);

         if( SCIPsetIsFeasGE(set, bound, SCIPvarGetLbLocal(var)) )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(tree), blkmem, set, stat, prob, tree, lp, branchcand,
                  eventqueue, var, bound, SCIP_BOUNDTYPE_UPPER, FALSE) );
         }
         else
            *cutoff = TRUE;

         /* count the bound change as applied cut */
         if( !sepastore->initiallp )
            sepastore->ncutsapplied++;
      }
   }
   else
   {
      /* apply the global bound change or detect a global cutoff which means we can cutoff the root node */
      if( SCIPsetIsLT(set, bound, SCIPvarGetUbGlobal(var)) )
      {
         SCIPdebugMessage(" -> applying global bound change: <%s>: [%g,%g] -> [%g,%g]\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), SCIPvarGetLbGlobal(var), bound);

         if( SCIPsetIsFeasGE(set, bound, SCIPvarGetLbGlobal(var)) )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, prob, tree, lp, branchcand,
                  eventqueue, var, bound, SCIP_BOUNDTYPE_UPPER, FALSE) );
         }
         else
         {
            /* we are done with solving since a global bound change is infeasible */
            SCIPnodeCutoff(SCIPtreeGetRootNode(tree), set, stat, tree);
            *cutoff = TRUE;
         }

         /* count the bound change as applied cut */
         if( !sepastore->initiallp )
            sepastore->ncutsapplied++;
      }
   }

   return SCIP_OKAY;
}

/** applies a cut that is a bound change directly as bound change instead of adding it as row to the LP */
static
SCIP_RETCODE sepastoreApplyBdchg(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_ROW*             cut,                /**< cut with a single variable */
   SCIP_Bool*            cutoff              /**< pointer to store whether an empty domain was created */
   )
{
   SCIP_COL** cols;
   SCIP_Real* vals;
   SCIP_VAR* var;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool local;

   assert(sepastore != NULL);
   assert(!SCIProwIsModifiable(cut));
   assert(SCIProwGetNNonz(cut) == 1);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* get the single variable and its coefficient of the cut */
   cols = SCIProwGetCols(cut);
   assert(cols != NULL);

   var = SCIPcolGetVar(cols[0]);
   vals = SCIProwGetVals(cut);
   assert(vals != NULL);
   assert(!SCIPsetIsZero(set, vals[0]));

   local = SCIProwIsLocal(cut);

   /* if the coefficient is nearly zero, we better ignore this cut for numerical reasons */
   if( SCIPsetIsFeasZero(set, vals[0]) )
      return SCIP_OKAY;

   /* get the left hand side of the cut and convert it to a bound */
   lhs = SCIProwGetLhs(cut);
   if( !SCIPsetIsInfinity(set, -lhs) )
   {
      lhs -= SCIProwGetConstant(cut);
      if( vals[0] > 0.0 )
      {
         /* coefficient is positive -> lhs corresponds to lower bound */
         SCIP_CALL( sepastoreApplyLb(sepastore, blkmem, set, stat, prob, tree, lp, branchcand, eventqueue,
               var, lhs/vals[0], local, cutoff) );
      }
      else
      {
         /* coefficient is negative -> lhs corresponds to upper bound */
         SCIP_CALL( sepastoreApplyUb(sepastore, blkmem, set, stat, prob, tree, lp, branchcand, eventqueue,
               var, lhs/vals[0], local, cutoff) );
      }
   }

   /* get the right hand side of the cut and convert it to a bound */
   rhs = SCIProwGetRhs(cut);
   if( !SCIPsetIsInfinity(set, rhs) )
   {
      rhs -= SCIProwGetConstant(cut);
      if( vals[0] > 0.0 )
      {
         /* coefficient is positive -> rhs corresponds to upper bound */
         SCIP_CALL( sepastoreApplyUb(sepastore, blkmem, set, stat, prob, tree, lp, branchcand, eventqueue,
               var, rhs/vals[0], local, cutoff) );
      }
      else
      {
         /* coefficient is negative -> rhs corresponds to lower bound */
         SCIP_CALL( sepastoreApplyLb(sepastore, blkmem, set, stat, prob, tree, lp, branchcand, eventqueue,
               var, rhs/vals[0], local, cutoff) );
      }
   }

   return SCIP_OKAY;
}

/** updates the orthogonalities and scores of the non-forced cuts after the given cut was added to the LP */
static
SCIP_RETCODE sepastoreUpdateOrthogonalities(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_ROW*             cut,                /**< cut that was applied */
   SCIP_Real             mincutorthogonality /**< minimal orthogonality of cuts to apply to LP */
   )
{
   int pos;

   assert(sepastore != NULL);

   pos = sepastore->nforcedcuts;
   while( pos < sepastore->ncuts )
   {
      SCIP_Real thisortho;
      
      /* update orthogonality */
      thisortho = SCIProwGetOrthogonality(cut, sepastore->cuts[pos], set->sepa_orthofunc);
      if( thisortho < sepastore->orthogonalities[pos] )
      {
         if( thisortho < mincutorthogonality )
         {
            /* cut is too parallel: release the row and delete the cut */
            SCIPdebugMessage("    -> deleting parallel cut <%s> after adding <%s> (pos=%d, len=%d, orthogonality=%g, score=%g)\n",
               SCIProwGetName(sepastore->cuts[pos]), SCIProwGetName(cut), pos, SCIProwGetNNonz(cut), thisortho, sepastore->scores[pos]);
            SCIP_CALL( sepastoreDelCut(sepastore, blkmem, set, eventqueue, eventfilter, lp, pos) );
            continue;
         }
         else
         {
            /* recalculate score */
            sepastore->orthogonalities[pos] = thisortho;
            assert( sepastore->objparallelisms[pos] != SCIP_INVALID ); /*lint !e777*/
            assert( sepastore->scores[pos] != SCIP_INVALID ); /*lint !e777*/
            assert( sepastore->efficacies[pos] != SCIP_INVALID ); /*lint !e777*/
            sepastore->scores[pos] = sepastore->efficacies[pos]
               + set->sepa_objparalfac * sepastore->objparallelisms[pos]
               + set->sepa_orthofac * thisortho;
         }
      }
      pos++;
   }

   return SCIP_OKAY;
}

/** applies the given cut to the LP and updates the orthogonalities and scores of remaining cuts */
static
SCIP_RETCODE sepastoreApplyCut(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_ROW*             cut,                /**< cut to apply to the LP */
   SCIP_Real             mincutorthogonality,/**< minimal orthogonality of cuts to apply to LP */
   int                   depth,              /**< depth of current node */
   int*                  ncutsapplied        /**< pointer to count the number of applied cuts */
   )
{
   assert(sepastore != NULL);
   assert(ncutsapplied != NULL);

   /* a row could have been added twice to the separation store; add it only once! */
   if( !SCIProwIsInLP(cut) )
   {
      /* add cut to the LP and capture it */
      SCIP_CALL( SCIPlpAddRow(lp, blkmem, set, eventqueue, eventfilter, cut, depth) );

      /* update statistics -> only if we are not in the initial lp (cuts are only counted if added during run) */
      if( !sepastore->initiallp )
      {
         sepastore->ncutsapplied++;

         /* increase count of applied cuts for origins of row */
         switch ( cut->origintype )
         {
         case SCIP_ROWORIGINTYPE_CONS:
            assert( cut->origin != NULL );
            SCIPconshdlrIncNAppliedCuts((SCIP_CONSHDLR*) cut->origin);
            break;
         case SCIP_ROWORIGINTYPE_SEPA:
            assert( cut->origin != NULL );
            SCIPsepaIncNAppliedCuts((SCIP_SEPA*) cut->origin);
            break;
         case SCIP_ROWORIGINTYPE_UNSPEC:
            /* do nothing - cannot update statistics */
            break;
         default:
            SCIPerrorMessage("unkown type of row origin.\n");
            return SCIP_INVALIDDATA;
         }
      }

      /* update the orthogonalities */
      SCIP_CALL( sepastoreUpdateOrthogonalities(sepastore, blkmem, set, eventqueue, eventfilter, lp, cut, mincutorthogonality) );
      (*ncutsapplied)++;
   }

   return SCIP_OKAY;
}

/** returns the position of the best non-forced cut in the cuts array */
static
int sepastoreGetBestCut(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   SCIP_Real bestscore;
   int bestpos;
   int pos;

   assert(sepastore != NULL);

   bestscore = SCIP_REAL_MIN;
   bestpos = -1;
   for( pos = sepastore->nforcedcuts; pos < sepastore->ncuts; pos++ )
   {
      /* check if cut is current best cut */
      assert( sepastore->scores[pos] != SCIP_INVALID ); /*lint !e777*/
      if( sepastore->scores[pos] > bestscore )
      {
         bestscore = sepastore->scores[pos];
         bestpos = pos;
      }
   }

   return bestpos;
}

/** computes score for current LP solution and initialized orthogonalities */
static
SCIP_RETCODE computeScore(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             handlepool,         /**< whether the efficacy of cuts in the pool should be reduced  */
   int                   pos,                /**< position of cut to handle */
   SCIP_EFFICIACYCHOICE  efficiacychoice     /**< type of solution to base efficiacy computation on */
   )
{
   SCIP_ROW* cut;
   SCIP_Real cutefficacy;
   SCIP_Real cutscore;

   cut = sepastore->cuts[pos];

   /* calculate cut's efficacy */
   switch ( efficiacychoice )
   {
   case SCIP_EFFICIACYCHOICE_LP:
      cutefficacy = SCIProwGetLPEfficacy(cut, set, stat, lp);
      break;
   case SCIP_EFFICIACYCHOICE_RELAX:
      cutefficacy = SCIProwGetRelaxEfficacy(cut, set, stat);
      break;
   case SCIP_EFFICIACYCHOICE_NLP:
      cutefficacy = SCIProwGetNLPEfficacy(cut, set, stat);
      break;
   default:
      SCIPerrorMessage("Invalid efficiacy choice.\n");
      return SCIP_INVALIDCALL;
   }

   /* If a cut is not member of the cut pool, we slightly decrease its score to prefer identical
    * cuts which are in the cut pool.  This is because the conversion of cuts into linear
    * constraints after a restart looks at the cut pool and cannot find tight non-pool cuts.
    */
   if( handlepool && !SCIProwIsInGlobalCutpool(cut) )
      cutefficacy *= 0.9999;

   /* calculate resulting score */
   assert( sepastore->objparallelisms[pos] != SCIP_INVALID ); /*lint !e777*/
   cutscore = cutefficacy + set->sepa_objparalfac * sepastore->objparallelisms[pos] + set->sepa_orthofac * 1.0;
   assert( !SCIPsetIsInfinity(set, cutscore) );

   sepastore->efficacies[pos] = cutefficacy;
   sepastore->scores[pos] = cutscore;

   /* make sure that the orthogonalities are initialized to 1.0 */
   sepastore->orthogonalities[pos] = 1.0;

   return SCIP_OKAY;
}

/** adds cuts to the LP and clears separation storage */
SCIP_RETCODE SCIPsepastoreApplyCuts(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_EFFICIACYCHOICE  efficiacychoice,    /**< type of solution to base efficiacy computation on */
   SCIP_Bool*            cutoff              /**< pointer to store whether an empty domain was created */
   )
{
   SCIP_NODE* node;
   SCIP_Real mincutorthogonality;
   int depth;
   int maxsepacuts;
   int ncutsapplied;
   int pos;

   assert(sepastore != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   SCIPdebugMessage("applying %d cuts\n", sepastore->ncuts);

   node = SCIPtreeGetCurrentNode(tree);
   assert(node != NULL);

   /* get maximal number of cuts to add to the LP */
   maxsepacuts = SCIPsetGetSepaMaxcuts(set, root);
   ncutsapplied = 0;

   /* get depth of current node */
   depth = SCIPnodeGetDepth(node);

   /* calculate minimal cut orthogonality */
   mincutorthogonality = (root ? set->sepa_minorthoroot : set->sepa_minortho);
   mincutorthogonality = MAX(mincutorthogonality, set->num_epsilon);

   /* Compute scores for all non-forced cuts and initialize orthogonalities - make sure all cuts are initialized again for the current LP solution */
   for( pos = sepastore->nforcedcuts; pos < sepastore->ncuts; pos++ )
   {
      SCIP_CALL( computeScore(sepastore, set, stat, lp, TRUE, pos, efficiacychoice) );
   }

   /* apply all forced cuts */
   for( pos = 0; pos < sepastore->nforcedcuts && !(*cutoff); pos++ )
   {
      SCIP_ROW* cut;

      cut = sepastore->cuts[pos];
      assert(SCIPsetIsInfinity(set, sepastore->scores[pos]));

      /* if the cut is a bound change (i.e. a row with only one variable), add it as bound change instead of LP row */
      if( !SCIProwIsModifiable(cut) && SCIProwGetNNonz(cut) == 1 )
      {
         SCIPdebugMessage(" -> applying forced cut <%s> as boundchange\n", SCIProwGetName(cut));
         SCIP_CALL( sepastoreApplyBdchg(sepastore, blkmem, set, stat, prob, tree, lp, branchcand, eventqueue, cut, cutoff) );
      }
      else
      {
         /* add cut to the LP and update orthogonalities */
         SCIPdebugMessage(" -> applying forced cut <%s>\n", SCIProwGetName(cut));
         /*SCIPdebug( SCIProwPrint(cut, set->scip->messagehdlr, NULL));*/
         SCIP_CALL( sepastoreApplyCut(sepastore, blkmem, set, eventqueue, eventfilter, lp, cut, mincutorthogonality, depth, &ncutsapplied) );
      }
   }

   /* apply non-forced cuts */
   while( ncutsapplied < maxsepacuts && sepastore->ncuts > sepastore->nforcedcuts && !(*cutoff) )
   {
      SCIP_ROW* cut;
      int bestpos;
      
      /* get best non-forced cut */
      bestpos = sepastoreGetBestCut(sepastore);
      assert(sepastore->nforcedcuts <= bestpos && bestpos < sepastore->ncuts);
      assert(sepastore->scores[bestpos] != SCIP_INVALID ); /*lint !e777*/
      assert(sepastore->efficacies[bestpos] != SCIP_INVALID ); /*lint !e777*/
      cut = sepastore->cuts[bestpos];
      assert(SCIProwIsModifiable(cut) || SCIProwGetNNonz(cut) != 1); /* bound changes are forced cuts */
      assert(!SCIPsetIsInfinity(set, sepastore->scores[bestpos]));
      
      SCIPdebugMessage(" -> applying cut <%s> (pos=%d/%d, len=%d, efficacy=%g, objparallelism=%g, orthogonality=%g, score=%g)\n",
         SCIProwGetName(cut), bestpos, sepastore->ncuts, SCIProwGetNNonz(cut), sepastore->efficacies[bestpos], sepastore->objparallelisms[bestpos],
         sepastore->orthogonalities[bestpos], sepastore->scores[bestpos]);
      /*SCIPdebug(SCIProwPrint(cut, set->scip->messagehdlr, NULL));*/

      /* capture cut such that it is not destroyed in sepastoreDelCut() */
      SCIProwCapture(cut);

      /* release the row and delete the cut (also issuing ROWDELETEDSEPA event) */
      SCIP_CALL( sepastoreDelCut(sepastore, blkmem, set, eventqueue, eventfilter, lp, bestpos) );

      /* Do not add (non-forced) non-violated cuts.
       * Note: do not take SCIPsetIsEfficacious(), because constraint handlers often add cuts w.r.t. SCIPsetIsFeasPositive().
       */
      if( SCIPsetIsFeasPositive(set, sepastore->efficacies[bestpos]) )
      {
         /* add cut to the LP and update orthogonalities */
         SCIP_CALL( sepastoreApplyCut(sepastore, blkmem, set, eventqueue, eventfilter, lp, cut, mincutorthogonality, depth, &ncutsapplied) );
      }

      /* release cut */
      SCIP_CALL( SCIProwRelease(&cut, blkmem, set, lp) );
   }

   /* clear the separation storage and reset statistics for separation round */
   SCIP_CALL( SCIPsepastoreClearCuts(sepastore, blkmem, set, eventqueue, eventfilter, lp) );

   return SCIP_OKAY;
}

/** clears the separation storage without adding the cuts to the LP */
SCIP_RETCODE SCIPsepastoreClearCuts(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp                  /**< LP data */
   )
{
   int c;

   assert(sepastore != NULL);

   SCIPdebugMessage("clearing %d cuts\n", sepastore->ncuts);

   /* release cuts */
   for( c = 0; c < sepastore->ncuts; ++c )
   {
      /* check, if the row deletions from separation storage events are tracked
       * if so, issue ROWDELETEDSEPA event
       */
      if( eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWDELETEDSEPA) != 0 )
      {
         SCIP_EVENT* event;

         SCIP_CALL( SCIPeventCreateRowDeletedSepa(&event, blkmem, sepastore->cuts[c]) );
         SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
      }
      
      SCIP_CALL( SCIProwRelease(&sepastore->cuts[c], blkmem, set, lp) );
   }

   /* reset counters */
   sepastore->ncuts = 0;
   sepastore->nforcedcuts = 0;
   sepastore->ncutsfoundround = 0;

   /* if we have just finished the initial LP construction, free the (potentially large) cuts array */
   if( sepastore->initiallp )
   {
      BMSfreeMemoryArrayNull(&sepastore->cuts);
      sepastore->cutssize = 0;
   }

   return SCIP_OKAY;
}

/** removes cuts that are inefficacious w.r.t. the current LP solution from separation storage without adding the cuts to the LP */
SCIP_RETCODE SCIPsepastoreRemoveInefficaciousCuts(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_EFFICIACYCHOICE  efficiacychoice     /**< type of solution to base efficiacy computation on */
   )
{
   int cnt;
   int c;

   assert( sepastore != NULL );

   /* check non-forced cuts only */
   cnt = 0;
   c = sepastore->nforcedcuts;
   while( c < sepastore->ncuts )
   {
      assert( sepastore->efficacies[c] == SCIP_INVALID ); /*lint !e777*/
      SCIP_CALL( computeScore(sepastore, set, stat, lp, FALSE, c, efficiacychoice) );
      if( !SCIPsetIsEfficacious(set, root, sepastore->efficacies[c]) )
      {
         SCIP_CALL( sepastoreDelCut(sepastore, blkmem, set, eventqueue, eventfilter, lp, c) );
         ++cnt;
      }
      else
         ++c;
   }
   SCIPdebugMessage("removed %d non-efficacious cuts\n", cnt);

   return SCIP_OKAY;
}

/** get cuts in the separation storage */
SCIP_ROW** SCIPsepastoreGetCuts(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->cuts;
}

/** get number of cuts in the separation storage */
int SCIPsepastoreGetNCuts(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncuts;
}

/** get total number of cuts found so far */
int SCIPsepastoreGetNCutsFound(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsfound;
}

/** get number of cuts found so far in current separation round */
int SCIPsepastoreGetNCutsFoundRound(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsfoundround;
}

/** get total number of cuts applied to the LPs */
int SCIPsepastoreGetNCutsApplied(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsapplied;
}
