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

/**@file   bounding_exact.h
 * @brief  safe exact rational bounding methods
 * @author Leon Eifler
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef __SCIP_BOUNDING_EXACT_C__
#define __SCIP_BOUNDING_EXACT_C__

#include <stdio.h>
#include <assert.h>

#include "scip/bounding_exact.h"
#include "scip/struct_set.h"
#include "scip/pub_message.h"
#include "scip/stat.h"
#include "scip/set.h"
#include "scip/rational.h"
#include "scip/clock.h"
#include "lpi/lpi.h"
#include "scip/lp.h"
#include "lpi/lpiex.h"
#include "scip/scip_prob.h"
#include "scip/scip.h"
#include "scip/pub_varex.h"
#include "scip/sepastoreex.h"
#include "scip/struct_scip.h"


static
char chooseBoundingMethod(
   SCIP_LPEX*            lpex,
   SCIP_SET*             set,
   SCIP_Bool             infeaslp,           /**< will dual bound method be applied to safely verify infeasible LP? */
   SCIP_PROB*            prob
   )
{
   char dualboundmethod;

   assert(!lpex->fplp->hasprovedbound);

   if( infeaslp )
   {
      /* check if neumair-scher is possible */
      if( SCIPlpexBSpossible(lpex) )
         dualboundmethod = 'n';
      /* check if project and shift is possible */
      else if( SCIPlpexPSpossible(lpex) )
         dualboundmethod = 'p';
      /* otherwise solve exactly */
      else
         dualboundmethod = 'e';
   }
   else
   {
      /* decide whether we want to interleave with exact LP call/basis verification
       * - given freq
       * or
       * - Neumair Shcherbina bound only nearly able to cutoff node
       */
      if( (lpex->interleavedbfreq > 0 && !SCIPsetIsInfinity(set, SCIPlpGetCutoffbound(lpex->fplp)) && SCIPgetDepth(set->scip) > 0
            && SCIPgetDepth(set->scip) % (lpex->interleavedbfreq) == 0)
         || (lpex->interleavedbfreq == 0 && SCIPsetIsGE(set, SCIPlpGetObjval(lpex->fplp, set, prob), SCIPlpGetCutoffbound(lpex->fplp))
            && SCIPlpGetObjval(lpex->fplp, set, prob) < SCIPlpGetCutoffbound(lpex->fplp)) )
      {
         dualboundmethod = 'e';
      }
      else
      {
         /* check if neumair-scher is possible */
         if( SCIPlpexBSpossible(lpex) )
            dualboundmethod = 'n';
         /* check if project and shift is possible */
         else if( SCIPlpexPSpossible(lpex) )
            dualboundmethod = 'p';
         /* otherwise solve exactly */
         else
            dualboundmethod = 'e';
      }
   }
   return dualboundmethod;
}

/** calculates y*b + min{(c - y*A)*x | lb <= x <= ub} for given vectors y and c;
 *  the vector b is defined with b[i] = lhs[i] if y[i] >= 0, b[i] = rhs[i] if y[i] < 0
 *  Calculating this value in interval arithmetics gives a proved lower LP bound for the following reason (assuming,
 *  we have only left hand sides):
 *           min{cx       |  b <=  Ax, lb <= x <= ub}
 *   >=      min{cx       | yb <= yAx, lb <= x <= ub}   (restriction in minimum is relaxed)
 *   == yb + min{cx - yb  | yb <= yAx, lb <= x <= ub}   (added yb - yb == 0)
 *   >= yb + min{cx - yAx | yb <= yAx, lb <= x <= ub}   (because yAx >= yb inside minimum)
 *   >= yb + min{cx - yAx |            lb <= x <= ub}   (restriction in minimum is relaxed)
 */
static
SCIP_RETCODE boundShift(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< Exact LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Bool             usefarkas,
   SCIP_Real*            safebound
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_INTERVAL* rhsinter;
   SCIP_INTERVAL* constantinter;
   SCIP_INTERVAL* xinter;
   SCIP_INTERVAL* ainter;
   SCIP_INTERVAL* atyinter;
   SCIP_INTERVAL* cinter;
   SCIP_INTERVAL ytb;
   SCIP_INTERVAL minprod;
   SCIP_ROW* row;
   SCIP_COL* col;
   SCIP_Real* y;
   SCIP_Real* ycol;
   SCIP_Real c;
   int i;
   int j;

   assert(lpex != NULL);
   assert(lp != NULL);
   assert(lp->solved);
   assert(set != NULL);
   assert(safebound != NULL);
   /* start timing */
   if ( usefarkas )
      SCIPclockStart(stat->provedinfeasbstime, set);
   else
      SCIPclockStart(stat->provedfeasbstime, set);

   /* allocate temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &y, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rhsinter, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &constantinter, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ycol, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ainter, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &atyinter, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cinter, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &xinter, lp->ncols) );

   SCIPdebugMessage("calling proved bound for %s LP\n", usefarkas ? "infeasible" : "feasible");

   /* reset proved bound status */
   lp->hasprovedbound = FALSE;

   /* calculate y^Tb */
   SCIPintervalSet(&ytb, 0.0);
   SCIPdebugMessage("ytb intervall computation with vectors:\n");

   /* create y, rhs and constant vector in interval arithmetic */
   for( j = 0; j < lp->nrows; ++j )
   {
      row = lp->rows[j];
      assert(row != NULL);

      /* create y vector in interval arithmetic, setting near zeros to zero */
      y[j] = (usefarkas ? row->dualfarkas : row->dualsol);

      if( SCIPlpiIsInfinity(lp->lpi, y[j]) )
	      y[j] = SCIPsetInfinity(set);

      if( SCIPlpiIsInfinity(lp->lpi, -y[j]) )
	      y[j] = -SCIPsetInfinity(set);

      /** @todo exiptodo: dual bounding improvement
       *  - should we also set nonzero values of y to zero if corresponding lhs/rhs is not finite (to improve dual bound)?
       *  - do such situations come up?
       */
      /* create rhs and constant vectors in interval arithmetic */
      if( SCIPsetIsFeasPositive(set, y[j]) )
      {
         SCIPintervalSet(&rhsinter[j], row->lhs);
         SCIPintervalSet(&constantinter[j], -1.0 * row->constant);
      }
      else if( SCIPsetIsFeasNegative(set, y[j]) )
      {
         SCIPintervalSet(&rhsinter[j], row->rhs);
         SCIPintervalSet(&constantinter[j], -1.0 * row->constant);
      }
      else
      {
         y[j] = 0.0;
         SCIPintervalSet(&rhsinter[j], 0.0);
         SCIPintervalSet(&constantinter[j], 0.0);
      }

      SCIPdebugMessage("   j=%d: b=[%g,%g] (lhs=%g, rhs=%g, const=%g, y=%g)\n", j, rhsinter[j].inf, rhsinter[j].sup, row->lhs,
            row->rhs, row->constant, y[j]);
   }
   /* substract constant from rhs in interval arithmetic and calculate y^Tb */
   SCIPintervalAddVectors(SCIPsetInfinity(set), rhsinter, lp->nrows, rhsinter, constantinter);
   SCIPintervalScalprodScalars(SCIPsetInfinity(set), &ytb, lp->nrows, rhsinter, y);

   SCIPdebugMessage("   resulting ytb=[%g,%g]\n", SCIPintervalGetInf(ytb), SCIPintervalGetSup(ytb));

   /* calculate min{(c^T - y^TA)x} */

   /* compute infimums of -A^Ty */
   roundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeDownwards();
   for( j = 0; j < lp->ncols; ++j )
   {
      col = lp->cols[j];
      assert(col != NULL);
      assert(col->nunlinked == 0);

      /* create -a.j vector in interval arithmetic and corresponding y vector and compute infimum of vector -a.j^Ty */
      for( i = 0; i < col->nlprows; ++i )
      {
         SCIP_INTERVAL val;
         SCIP_ROWEX* rowex;

         assert(col->rows[i] != NULL);
         assert(col->rows[i]->lppos >= 0);
         assert(col->linkpos[i] >= 0);

         rowex = SCIProwGetExRow(lpex, col->rows[i]);

         val = rowex->valsinterval[col->linkpos[i]];
         assert(val.inf <= col->vals[i] && col->vals[i] <= val.sup);

         SCIPintervalSetBounds(&ainter[i], -val.sup, -val.inf);
         ycol[i] = y[col->rows[i]->lppos];
      }
      atyinter[j].inf = 0.0;
      SCIPintervalScalprodScalarsInf(SCIPsetInfinity(set), &atyinter[j], col->nlprows, ainter, ycol);

#ifndef NDEBUG
      for( i = col->nlprows; i < col->len; ++i )
      {
         assert(col->rows[i] != NULL);
         assert(col->rows[i]->lppos == -1);
         assert(col->rows[i]->dualsol == 0.0);
         assert(col->rows[i]->dualfarkas == 0.0);
         assert(col->linkpos[i] >= 0);
      }
#endif
   }

   /* compute supremums of -A^Ty */
   SCIPintervalSetRoundingModeUpwards();
   for( j = 0; j < lp->ncols; ++j )
   {
      col = lp->cols[j];
      assert(col != NULL);
      assert(col->nunlinked == 0);

      /* create -a.j vector in interval arithmetic and corresponding y vector and compute supremums of vector -a.j^Ty */
      for( i = 0; i < col->nlprows; ++i )
      {
         SCIP_INTERVAL val;
         SCIP_ROWEX* rowex;

         assert(col->rows[i] != NULL);
         assert(col->rows[i]->lppos >= 0);
         assert(col->linkpos[i] >= 0);

         rowex = SCIProwGetExRow(lpex, col->rows[i]);

         val = rowex->valsinterval[col->linkpos[i]];

         assert(val.inf <= col->vals[i] && col->vals[i] <= val.sup);

         SCIPintervalSetBounds(&ainter[i], -val.sup, -val.inf);
         ycol[i] = y[col->rows[i]->lppos];
      }
      atyinter[j].sup = 0.0;
      SCIPintervalScalprodScalarsSup(SCIPsetInfinity(set), &atyinter[j], col->nlprows, ainter, ycol);

#ifndef NDEBUG
      for( i = col->nlprows; i < col->len; ++i )
      {
         assert(col->rows[i] != NULL);
         assert(col->rows[i]->lppos == -1);
         assert(col->rows[i]->dualsol == 0.0);
         assert(col->rows[i]->dualfarkas == 0.0);
         assert(col->linkpos[i] >= 0);
      }
#endif
   }
   SCIPintervalSetRoundingMode(roundmode);

   /* create c vector and x vector in interval arithmetic and compute min{(c^T - y^TA)x} */
   for( j = 0; j < lp->ncols; ++j )
   {
      //assert(!SCIPsetIsInfinity(set, -SCIPcolGetLb(col)));
      //assert(!SCIPsetIsInfinity(set, SCIPcolGetUb(col)));
      col = lp->cols[j];
      assert(col != NULL);
      assert(col->nunlinked == 0);

      if( usefarkas )
         SCIPintervalSet(&cinter[j], 0);
      else
      {
         if( RisFpRepresentable(SCIPvarGetObjExact(SCIPcolGetVar(col))) )
            SCIPintervalSet(&cinter[j], col->obj);
         else
         {
            SCIPintervalSetRational(&cinter[j], SCIPvarGetObjExact(SCIPcolGetVar(col)));
         }
      }
      /** @todo exip: get exact column bounds ? */
      SCIPintervalSetBounds(&xinter[j], SCIPcolGetLb(col), SCIPcolGetUb(col));
      if( (SCIPsetIsInfinity(set, -SCIPcolGetLb(col)) || SCIPsetIsInfinity(set, SCIPcolGetUb(col))) )
      {
         SCIPmessagePrintWarning(messagehdlr, "warning: trying bound shift with unbounded column variable. Column %d, lb: %e, ub %e \n",
               SCIPcolGetIndex(col), SCIPcolGetLb(col) ,SCIPcolGetUb(col) );
         SCIPmessagePrintWarning(messagehdlr, "Multiplied with interval: min %e,  max %e \n",
               atyinter[j].inf + cinter[j].inf, atyinter[j].sup + cinter[j].sup);
      }
   }
   SCIPintervalAddVectors(SCIPsetInfinity(set), atyinter, lp->ncols, atyinter, cinter);
   SCIPintervalScalprod(SCIPsetInfinity(set), &minprod, lp->ncols, atyinter, xinter);

   /* add y^Tb */
   SCIPintervalAdd(SCIPsetInfinity(set), &minprod, minprod, ytb);

   /* free buffer for storing y in interval arithmetic */
   SCIPsetFreeBufferArray(set, &xinter);
   SCIPsetFreeBufferArray(set, &cinter);
   SCIPsetFreeBufferArray(set, &atyinter);
   SCIPsetFreeBufferArray(set, &ainter);
   SCIPsetFreeBufferArray(set, &ycol);
   SCIPsetFreeBufferArray(set, &constantinter);
   SCIPsetFreeBufferArray(set, &rhsinter);
   SCIPsetFreeBufferArray(set, &y);

   *safebound = SCIPintervalGetInf(minprod);
   //printf("safebound computed: %e, previous fp-bound: %e, difference %e \n", *safebound, lp->lpobjval, *safebound - lp->lpobjval);

   /* stop timing and update number of calls and fails, and proved bound status */
   if ( usefarkas )
   {
      SCIPclockStop(stat->provedinfeasbstime, set);
      stat->nboundshiftinf++;
      if( *safebound <= 0.0 )
      {
         stat->nfailboundshiftinf++;
         assert(!lp->hasprovedbound);
      }
      else
         lp->hasprovedbound = TRUE;

   }
   else
   {
      SCIPclockStop(stat->provedfeasbstime, set);
      stat->nboundshift++;
      if( !SCIPsetIsInfinity(set, -1.0 * (*safebound)) )
      {
#ifdef WITH_EXACTSOLVE
         /* SCIP_CONS** conss;

         conss = SCIPgetConss(set->scip);
         assert(conss != NULL);
         assert(SCIPgetNConss(set->scip) == 1); */

         //SCIP_CALL( SCIPcomputeDualboundQuality(set->scip, conss[0], *safebound) );
#endif
         lp->hasprovedbound = TRUE;
      }
      else
      {
         stat->nfailboundshift++;
         assert(!lp->hasprovedbound);
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE projectShiftInterval(
   void
)
{
   return SCIP_OKAY;
}

static
SCIP_RETCODE projectShiftRational(
   void
)
{
   return SCIP_OKAY;
}

static
SCIP_RETCODE basisVerification(
   void
)
{
   return SCIP_OKAY;
}

static
SCIP_RETCODE solveLpExact(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< Exact LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Longint          itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool             usefarkas,
   SCIP_Real*            safebound
   )
{
   int* cstat;
   int* rstat;
   SCIP_LPALGO lpalgo = SCIP_LPALGO_DUALSIMPLEX;
   SCIP_RETCODE retcode;
   int niterations = 0;

   assert(lp != NULL);
   assert(lpex != NULL);
   assert(set != NULL);
   assert(set->misc_exactsolve);

   if ( usefarkas )
      SCIPclockStart(stat->provedinfeaslptime, set);
   else
      SCIPclockStart(stat->provedfeaslptime, set);

   /* set up the exact lpi for the current node */
   SCIP_CALL( SCIPsepastoreexApplyCuts(set->scip->sepastoreex, blkmem, set, stat, lpex, eventqueue, eventfilter) );
   SCIP_CALL( SCIPlpexFlush(lp->lpex, blkmem, set, eventqueue) );

   assert(SCIPlpexIsSynced(lpex, set, messagehdlr));

   SCIP_CALL( SCIPsetAllocBufferArray(set, &cstat, lp->nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rstat, lp->nlpirows) );

   /* set the correct basis information for warmstart */
   SCIP_CALL( SCIPlpiGetBase(lp->lpi, cstat, rstat) );
   SCIP_CALL( SCIPlpiexSetBase(lpex->lpiex, cstat, rstat) );

   /* solve the lp exactly */
   switch(lpalgo)
   {
      case SCIP_LPALGO_PRIMALSIMPLEX:
         retcode = SCIPlpiexSolvePrimal(lpex->lpiex);
         break;
      case SCIP_LPALGO_DUALSIMPLEX:
         retcode = SCIPlpiexSolveDual(lpex->lpiex);
         break;
      default:
         SCIPerrorMessage("Lp-algorithm-type %d is not supported in exact solving mode \n", lpalgo);
         SCIPABORT();
   }
   if( retcode == SCIP_LPERROR )
   {
      *lperror = TRUE;
      SCIPdebugMessage("Error solving lp exactly in node %"SCIP_LONGINT_FORMAT" \n", SCIPnodeGetNumber(SCIPgetCurrentNode(set->scip)));
   }

   SCIPlpiexGetIterations(lpex->lpiex, &niterations);
   if( usefarkas )
      stat->niterationsexlpinf += niterations;
   else
      stat->niterationsexlp += niterations;

   if( SCIPlpiexIsOptimal(lpex->lpiex) )
   {
      if( usefarkas )
         SCIPdebugMessage("Failed to prove infeasibility, exact lp is solved optimally \n");

      /* evaluate solution status and set safe bound correctly */
      SCIP_CALL( SCIPlpiexGetSol(lpex->lpiex, lpex->lpobjval, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPlpiexGetObjval(lpex->lpiex, lpex->lpobjval) );
      SCIPdebugMessage("Exact lp solve terminated with optimal. Safe dual bound is %e, previous lp obj-val was %e \n", 
            RgetRealRelax(lpex->lpobjval, SCIP_ROUND_DOWNWARDS), lp->lpobjval);
      lp->lpobjval = RgetRealRelax(lpex->lpobjval, SCIP_ROUND_DOWNWARDS);
      lp->hasprovedbound = TRUE;
   }

   /* stop timing and update number of calls and fails, and proved bound status */
   if ( usefarkas )
   {
      if( !SCIPlpiexIsPrimalInfeasible(lpex->lpiex) )
         stat->nfailexlpinf++;
      else
         lp->hasprovedbound = TRUE;

      SCIPclockStop(stat->provedinfeaslptime, set);
      stat->nexlpinf++;
   }
   else
   {
      if( !SCIPlpiexIsOptimal(lpex->lpiex) )
         stat->nfailexlp++;
   
      SCIPclockStop(stat->provedfeaslptime, set);
      stat->nexlp++;
      if( !SCIPsetIsInfinity(set, -1.0 * (lp->lpobjval)) )
      {
#ifdef WITH_EXACTSOLVE
         /* SCIP_CONS** conss;

         conss = SCIPgetConss(set->scip);
         assert(conss != NULL);
         assert(SCIPgetNConss(set->scip) == 1); */

         //SCIP_CALL( SCIPcomputeDualboundQuality(set->scip, conss[0], *safebound) );
#endif
         lp->hasprovedbound = TRUE;
      }
      else
      {
         assert(!lp->hasprovedbound);
      }
   }

   SCIPsetFreeBufferArray(set, &rstat);
   SCIPsetFreeBufferArray(set, &cstat);

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPcomputeSafeBound(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< Exact LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Longint          itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool             dualfarkas,
   SCIP_Real*            safebound
   )
{
   char dualboundmethod;

   /* If we are not in exact solving mode, just return */
   if( !set->misc_exactsolve )
      return SCIP_OKAY;

   assert(set->misc_exactsolve);

   /* choose which bounding method to use. only needed if automatic is enabled. */
   if( set->misc_dbmethod == 'a' )
      dualboundmethod = chooseBoundingMethod(lpex, set, dualfarkas, prob);
   else
      dualboundmethod = set->misc_dbmethod;

   switch(dualboundmethod)
   {
      /* neumaier and scherbina */
      case 'n':
         SCIP_CALL( boundShift(lp, lpex, set, messagehdlr, blkmem, stat, eventqueue, eventfilter,
                        prob, dualfarkas, safebound) );
         break;
      /* basis verification */
      case 'v':
         SCIPerrorMessage("bounding method %c not implemented yet \n", set->misc_dbmethod);
         SCIPABORT();
         break;
      /* repair lp basis */
      case 'r':
         SCIPerrorMessage("bounding method %c not implemented yet \n", set->misc_dbmethod);
         SCIPABORT();
         break;
      /* project and scale */
      case 'p':
         SCIPerrorMessage("bounding method %c not implemented yet \n", set->misc_dbmethod);
         SCIPABORT();
         break;
      /* exact LP */
      case 'e':
         SCIP_CALL( solveLpExact(lp, lpex, set, messagehdlr, blkmem, stat, eventqueue, eventfilter,
                        prob, itlim, lperror, dualfarkas, safebound) );
         break;
      /* interval neumaier and scherbina */
      case 'i':
         SCIPerrorMessage("bounding method %c not implemented yet \n", set->misc_dbmethod);
         SCIPABORT();
         break;
      /* exact neumaier and scherbina */
      case 'x':
         SCIPerrorMessage("bounding method %c not implemented yet \n", set->misc_dbmethod);
         SCIPABORT();
         break;
      default:
         SCIPerrorMessage("bounding method %c not implemented yet \n", set->misc_dbmethod);
         SCIPABORT();
         break;
   }

   /* choose which bounding method should be calles and return a safe objective bound */
   return SCIP_OKAY;
}

#endif