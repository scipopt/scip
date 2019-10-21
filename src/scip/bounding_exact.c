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
#include <time.h>
 
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
#include "scip/prob.h"
#include "scip/scip.h"
#include "scip/solex.h"
#include "scip/sol.h"
#include "scip/pub_varex.h"
#include "scip/primal.h"
#include "scip/sepastoreex.h"
#include "scip/struct_scip.h"
#include "rectlu/rectlu.h"

#define PSBIGM                      100
#define PSWARMSTARTAUXPROB         TRUE
#define PSPOSTPROCESSDUALSOL       TRUE

static
SCIP_Bool fpLPisIntFeasible(
   SCIP_LP*              lp,
   SCIP_SET*             set,
   SCIP_STAT*            stat
   )
{
   SCIP_Real primsol;
   SCIP_Bool feasible;
   int nfracs;
   int i;

   assert( SCIPlpIsSolved(lp));

   feasible = TRUE;
   for( i = 0; i < lp->ncols && feasible; i++ )
   {
      SCIP_COL* col;
      col = lp->cols[i];
      if( SCIPvarGetType(SCIPcolGetVar(col)) == SCIP_VARTYPE_CONTINUOUS )
         continue;
      primsol = SCIPcolGetPrimsol(col);
      if( !SCIPsetIsIntegral(set, primsol) )
         feasible = FALSE;
   }

   return feasible;
}

/** evaluates the result of the exact LP */
static
SCIP_RETCODE evaluateLPEX(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< Exact LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_RESULT*          result              /**< pointer to store the result of the lp enforcement call */
   )
{
   /** @todo exiptodo:
    *  - in a similar function for LP in lp.c the case SCIPlpiIsPrimalUnbounded() is not explicitely handled. why? if case
    *    should be added, include it here as well.
    */
   /* evaluate solution status */
   if( SCIPlpiexIsOptimal(lpex->lpiex) )
   {
      SCIPdebugMessage("   exact LP solved to optimality\n");

#ifndef NDEBUG
      {
         SCIP_Bool primalfeasible;
         SCIP_Bool dualfeasible;

         SCIP_CALL( SCIPlpiexGetSolFeasibility(lpex->lpiex, &primalfeasible, &dualfeasible) );
         assert(primalfeasible);
         assert(dualfeasible);
      }
#endif
      /* check whether exact LP solution is feasible for the MIP; if it is feasible the solution is stored
       * and the current node is cut off otherwise a branching is created
       */
      assert(*result == SCIP_CUTOFF || *result == SCIP_BRANCHED || *result == SCIP_SOLVELP );
   }
   else if( SCIPlpiexIsObjlimExc(lpex->lpiex) )
   {
      /* CERT: for now we don't trust the exact LP solver on this, so we do not need to print anything here */
      SCIPdebugMessage("   exact LP exceeds upper objective limit\n");
      *result = SCIP_CUTOFF;
   }
   else if( SCIPlpiexIsPrimalInfeasible(lpex->lpiex) )
   {
      SCIPdebugMessage("   exact LP is primal infeasible\n");
      *result = SCIP_CUTOFF;
   }
   else if( SCIPlpiexExistsPrimalRay(lpex->lpiex) )
   {
      /* CERT: for now we don't need to handle the unbounded case */
      SCIPerrorMessage("exact LP has primal ray: case not handled yet\n");
      return SCIP_ERROR;
   }
   else if( SCIPlpiexIsIterlimExc(lpex->lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds iteration limit: case not handled yet\n");
      return SCIP_ERROR;
   }
   else if( SCIPlpiexIsTimelimExc(lpex->lpiex) )
   {
      /* CERT: for now we should not activate the buggy time limit support of QSopt_ex */
      SCIPdebugMessage("   exact LP exceeds time limit\n");
      *result = SCIP_INFEASIBLE;
   }
   else
   {
      SCIPerrorMessage(" error or unknown return status in current exact LP (internal status: %d)\n",
         SCIPlpiexGetInternalStatus(lpex->lpiex));
      return SCIP_LPERROR;
   }

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
   SCIP_Bool             usefarkas           /**< are we aiming to prove infeasibility? */
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
   SCIP_CALL( SCIPsepastoreexSyncLPs(set->scip->sepastoreex, blkmem, set, stat, lpex, eventqueue, eventfilter) );
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

   lpex->solved = TRUE;

   SCIPlpiexGetIterations(lpex->lpiex, &niterations);
   if( usefarkas )
      stat->niterationsexlpinf += niterations;
   else
      stat->niterationsexlp += niterations;

   if( SCIPlpiexIsOptimal(lpex->lpiex) )
   {
      /* evaluate solution status and set safe bound correctly */
      SCIP_CALL( SCIPlpexGetSol(lpex, set, stat, NULL, NULL) );
      SCIP_CALL( SCIPlpiexGetObjval(lpex->lpiex, lpex->lpobjval) );
      SCIPdebugMessage("Exact lp solve terminated with optimal. Safe dual bound is %e, previous lp obj-val was %e \n",
            RatRoundReal(lpex->lpobjval, SCIP_ROUND_DOWNWARDS), lp->lpobjval);
      lp->lpobjval = RatRoundReal(lpex->lpobjval, SCIP_ROUND_DOWNWARDS);
      lp->hasprovedbound = TRUE;
      lpex->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   }
   else if( SCIPlpiexIsPrimalUnbounded(lpex->lpiex) )
   {
      /** @todo exip: where to save the ray? */
      SCIP_CALL( SCIPlpexGetPrimalRay(lpex, set, NULL) );
      lp->hasprovedbound = TRUE;
      lpex->lpsolstat = SCIP_LPSOLSTAT_UNBOUNDEDRAY;
   }
   else if( SCIPlpiexIsPrimalInfeasible(lpex->lpiex) )
   {
      SCIP_Bool valid;
      lpex->lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;

      SCIP_CALL( SCIPlpexGetDualfarkas(lpex, set, stat, &valid) );
      if( valid )
      {
         lp->hasprovedbound = TRUE;
         lp->lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;
      }
   }
   else
   {
      lp->hasprovedbound = FALSE;
      SCIPdebugMessage("Exact lp solve failed. Terminated with status %d \n", SCIPlpiexGetInternalStatus(lpex->lpiex));
      if( usefarkas )
         stat->nfailexlpinf++;
      else
         stat->nfailexlp++;
   }
   /* stop timing and update number of calls and fails, and proved bound status */
   if ( usefarkas )
   {
      SCIPclockStop(stat->provedinfeaslptime, set);
      stat->nexlpinf++;
   }
   else
   {
      SCIPclockStop(stat->provedfeaslptime, set);
      stat->nexlp++;
   }

   SCIPsetFreeBufferArray(set, &rstat);
   SCIPsetFreeBufferArray(set, &cstat);

   return SCIP_OKAY;
}

/** helper method, compute number of nonzeros in lp */
static
int getNNonz(
   SCIP_LPEX* lpex
   )
{
   int ret = 0;
   int i = 0;

   assert(lpex != NULL);

   for( i = 0; i < lpex->nrows; ++i )
      ret+= lpex->rows[i]->len;

   return ret;
}

/** allocate memory for ps interior point/ray computation */
static
SCIP_RETCODE allocIntMem(
   SCIP_SET*             set,                /**< scip settings */
   SCIP_Rational***      psobj,              /**< obj function */
   SCIP_Rational***      pslb,               /**< lower bounds */
   SCIP_Rational***      psub,               /**< upper bounds */
   SCIP_Rational***      pslhs,              /**< left hand sides */
   SCIP_Rational***      psrhs,              /**< right hand sides */
   SCIP_Rational***      psval,              /**< matrix values */
   SCIP_Rational***      sol,                /**< solution array */
   SCIP_Rational**       objval,             /**< objective values */
   int**                 psbeg,              /**< beg of each row in val array */
   int**                 pslen,              /**< len of each row */
   int**                 psind,              /**< index of val-entries */
   char***               colnames,           /**< column names */
   int                   psncols,            /**< number of columsns*/
   int                   psnrows,            /**< number of rows */
   int                   psnnonz             /**< number of non-zeros */
   )
{
   int i;
   assert(set != NULL);

   /* allocate memory for aux problem */
   SCIP_CALL( RatCreateBufferArray(set->buffer, psobj, psncols) );
   SCIP_CALL( RatCreateBufferArray(set->buffer, pslb, psncols) );
   SCIP_CALL( RatCreateBufferArray(set->buffer, psub, psncols) );
   SCIP_CALL( RatCreateBufferArray(set->buffer, pslhs, psnrows) );
   SCIP_CALL( RatCreateBufferArray(set->buffer, psrhs, psnrows) );
   SCIP_CALL( RatCreateBufferArray(set->buffer, psval, psnnonz) );
   SCIP_CALL( RatCreateBufferArray(set->buffer, sol, psncols) );
   SCIP_CALL( RatCreateBuffer(set->buffer, objval) );

   SCIP_CALL( SCIPsetAllocBufferArray(set, psbeg, psnrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, pslen, psnrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, psind, psnnonz) );

   SCIP_CALL( SCIPsetAllocBufferArray(set, colnames, psncols) );
   for( i = 0; i < psncols; i++ )
   {
      SCIP_CALL( SCIPsetAllocBufferArray(set, &((*colnames)[i]), SCIP_MAXSTRLEN ) );
      (void) SCIPsnprintf( (*colnames)[i] , SCIP_MAXSTRLEN, "var%d", i);
   }

   return SCIP_OKAY;
}

/** subroutine of constructPSdata(); chooses which columns of the matrix are designated as set S, used for projections */
static
SCIP_RETCODE psChooseS(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< exact LP data */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_STAT*            stat,               /**< statistics pointer */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem              /**< block memory struct */
   )
{
   int i;
   int nrows;
   int ncols;
   int nextendedrows;

   /* solution information for exact root LP */
   SCIP_Rational** rootactivity;
   SCIP_Rational** rootprimal;
   SCIP_Bool lperror;
   SCIP_PSDATA* psdata;

   nrows = lpex->nrows;
   ncols = lpex->ncols;
   psdata = lpex->psdata;

   assert(psdata != NULL);

   nextendedrows = psdata->nextendedrows;

   /* build includedrows vector based on psdualcolselection, this determines the matrix D */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &psdata->includedrows, nextendedrows) );
   for( i = 0; i < nextendedrows; i++ )
      psdata->includedrows[i] = 0;
   if( psdata->psdualcolselection == PS_DUALCOSTSEL_NO || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE )
   {
      /* determine which dual variables to included in the problem
       * (ones with finite dual objective coef. in [lhs',-rhs',lb',-ub'])
       */
      for( i = 0; i < nrows; i++ )
      {
         if( !RatIsNegInfinity(lpex->rows[i]->lhs) )
            psdata->includedrows[i] = 1;
         if( !RatIsInfinity(lpex->rows[i]->rhs) )
            psdata->includedrows[nrows + i] = 1;
      }
      for( i = 0; i < ncols; i++ )
      {
         if( !RatIsNegInfinity(lpex->cols[i]->lb) )
            psdata->includedrows[2*nrows + i] = 1;
         if( !RatIsInfinity(lpex->cols[i]->ub) )
            psdata->includedrows[2*nrows + ncols + i] = 1;
      }
   }
   else if( psdata->psdualcolselection == PS_DUALCOSTSEL_ACTIVE_EXLP )
   {
      /* determone which dual variables to include in the problem (in this case we choose dual variables whose primal
       * constraints are active at the solution of the exact LP at the root node)
       */

      solveLpExact(lp, lpex, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, prob, 100, &lperror, FALSE);

      SCIP_CALL( RatCreateBufferArray(set->buffer, &rootprimal, ncols) );
      SCIP_CALL( RatCreateBufferArray(set->buffer, &rootactivity, nrows) );

      /* get the primal solution and activity */
      SCIP_CALL( SCIPlpiexGetSol(lpex->lpiex, NULL, rootprimal, NULL, rootactivity, NULL) );

      /* determine which dual variables to include in the problem
       * (primal constraints active at optimal solution found at root node)
       */
      for( i = 0; i < nrows; i++ )
      {
         if( RatIsEqual(rootactivity[i], lpex->rows[i]->lhs) )
            psdata->includedrows[i] = 1;
         if( RatIsEqual(rootactivity[i], lpex->rows[i]->rhs) )
            psdata->includedrows[nrows + i] = 1;
      }
      for( i = 0; i < ncols; i++ )
      {
         if( RatIsEqual(rootprimal[i], lpex->cols[i]->lb) )
            psdata->includedrows[2*nrows + i] = 1;
         if( RatIsEqual(rootprimal[i], lpex->cols[i]->ub) )
            psdata->includedrows[2*nrows + ncols + i] = 1;
      }

      RatFreeBufferArray(set->buffer, &rootactivity, nrows);
      RatFreeBufferArray(set->buffer, &rootprimal, ncols);
   }
   else if( psdata->psdualcolselection == PS_DUALCOSTSEL_ACTIVE_FPLP )
   {
      /* determine which dual variables to include in the problem (in this case we choose dual variables whose primal
       * constraints are active at the solution of the exact LP at the root node)
       */

      assert(lp->nrows == nrows);
      for( i = 0; i < nrows; i++ )
      {
         if( SCIPsetIsFeasEQ(set, SCIProwGetLPActivity(lp->rows[i], set, stat, lp), SCIProwGetLhs(lp->rows[i])) )
            psdata->includedrows[i] = 1;
         if( SCIPsetIsFeasEQ(set, SCIProwGetLPActivity(lp->rows[i], set, stat, lp), SCIProwGetRhs(lp->rows[i])) )
            psdata->includedrows[nrows + i] = 1;
      }

      assert(ncols == lp->ncols);
      for( i = 0; i < ncols; i++ )
      {
         if( SCIPsetIsFeasEQ(set, SCIPcolGetPrimsol(lp->cols[i]), SCIPcolGetLb(lp->cols[i])) )
            psdata->includedrows[2*nrows + i] = 1;
         if( SCIPsetIsFeasEQ(set, SCIPcolGetPrimsol(lp->cols[i]), SCIPcolGetUb(lp->cols[i])) )
            psdata->includedrows[2*nrows + ncols + i] = 1;
      }
   }
   else
   {
      SCIPerrorMessage("Invald value for parameter psdualcolselection \n");
   }
   return SCIP_OKAY;
}

/** subroutine of constructPSdata(); computes the LU factorization used by the project-and-shift method */
static
SCIP_RETCODE psFactorizeD(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< exact LP data */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             findintpoint        /**< TRUE, if we search int point, FALSE if we search for ray */
   )
{

   int i;
   int j;
   int rval;
   int pos;
   int nrows;
   int ncols;
   int nextendedrows;
   int nnonz;
   mpq_t* projvalgmp;

   /* sparse representation of the matrix used for the LU factorization */
   int* projbeg;
   int* projlen;
   int* projind;
   SCIP_Rational** projval;
   SCIP_PSDATA* psdata;

   psdata = lpex->psdata;

   nrows = lpex->nrows;
   ncols = lpex->ncols;
   nextendedrows = psdata->nextendedrows;
   nnonz = getNNonz(lpex);

   /* allocate memory for the projection factorization */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &projbeg, nextendedrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &projlen, nextendedrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &projind, 2*nnonz + 2*ncols) );
   for( i = 0; i < 2*nnonz + 2*ncols; i++)
      projind[i]=0;
   SCIP_CALL( RatCreateBufferArray(set->buffer, &projval, 2*nnonz + 2*ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &projvalgmp, 2*nnonz + 2*ncols) );

   /* allocate memory for the basis mapping */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &psdata->psbasis, nextendedrows) );

   /* use includedrows to construct psbasis, a description/mapping for D it has length npsbasis and psbasis[i] tells what
    * column (out of the original nextendecons) the ith column in D is
    */
   pos = 0;
   for( i = 0; i < nextendedrows; i++ )
   {
      if( psdata->includedrows[i] )
      {
         psdata->psbasis[pos] = i;
         pos++;
      }
   }
   psdata->npsbasis = pos;

   /* build the sparse representation of D that will be passed to the RECTLU code for factorization */
   pos = 0;
   for( i = 0; i < nextendedrows; i++ )
   {
      /* A part (lhs constraints) */
      if(i < nrows)
      {
         projlen[i] = lpex->rows[i]->len;
         projbeg[i] = pos;
         for(j = 0; j < projlen[i]; j++)
         {
            projind[ projbeg[i] + j ] = lpex->rows[i]->cols_index[j];
            RatSet( projval[ projbeg[i] + j], lpex->rows[i]->vals[j]);
         }
         pos += projlen[i];
      }
      /* -A part (rhs constraints) */
      else if(i < 2 * nrows)
      {
         projlen[i] = lpex->rows[i - nrows]->len;
         projbeg[i] = pos;
         for(j = 0; j < projlen[i]; j++)
         {
            projind[ projbeg[i] + j ] = lpex->rows[i - nrows]->cols_index[j];
            RatNegate( projval[ projbeg[i] + j], lpex->rows[i - nrows]->vals[j]);
         }
         pos += projlen[i];
      }
      /* I part (lb constraints) */
      else if (i < 2*nrows + ncols)
      {
         projbeg[i] = pos;
         projlen[i] = 1;
         projind[pos] = i - 2*nrows;
         RatSetInt(projval[pos], 1, 1);
         pos ++;
      }
      /* -I part (ub constraints) */
      else
      {
         projbeg[i] = pos;
         projlen[i] = 1;
         projind[pos] = i - (2*nrows + ncols);
         RatSetInt(projval[pos], -1, 1);
         pos ++;
      }
   }

#ifdef PS_OUT
   printf("factoring matrix: ncols=%d, npsbasis=%d\n", ncols, psdata->npsbasis);
   for( i = 0; i < nextendedrows; i++ )
      printf("   j=%d:\t projbeg=<%d>,\t projlen=<%d>\n", i, projbeg[i], projlen[i]);

   for( i = 0; i < 2*nnonz + 2*ncols; i++ )
   {
      printf("   i=%d:\t projind=<%d>,\t projval=<", i, projind[i]);
      RatPrint(projval[i]);
      printf(">\n");
   }
#endif

   /* factorize projection matrix D
    * - psbasis stores a mapping to tell us what D is, i.e. the dual columns corresponding to dual valuse that have a
    *   strictly positive value in the relative interior point
    * - D is equal to a subset of [A',-A',I,-I] and is given to the factor code in sparse column representation
    */
   RatSetGMPArray(projvalgmp, projval, 2 * nnonz + 2 * ncols);

   rval = RECTLUbuildFactorization(&psdata->rectfactor, ncols, psdata->npsbasis,
      psdata->psbasis, projvalgmp, projind, projbeg, projlen);

   /* if rval != 0 then RECTLUbuildFactorization has failed. In this case the project-and-shift method will not work and
    * we will return failure
    */
   if( rval )
   {
      psdata->psdatafail = TRUE;
      SCIPdebugMessage("factorization of matrix for project-and-shift method failed.\n");
   }

#ifdef PS_OUT
   printf("   matrix factorization complete: %s\n", rval ? "failed" : "correct termination");
#endif

   RatClearGMPArray(projvalgmp, 2 * nnonz+ 2 * ncols);
   SCIPsetFreeBufferArray(set, &projvalgmp);
   RatFreeBufferArray(set->buffer, &projval, 2*nnonz + 2*ncols);
   SCIPsetFreeBufferArray(set, &projind);
   SCIPsetFreeBufferArray(set, &projlen);
   SCIPsetFreeBufferArray(set, &projbeg);

   return SCIP_OKAY;
}

/** prints error related to the current lpiex status, if there is one */
static
SCIP_RETCODE printlpiexerr(
   SCIP_LPIEX*           lpiex              /**< lpiex interface */
   )
{
   if( SCIPlpiexIsOptimal(lpiex) )
   {
      return SCIP_OKAY;
   }
   else if( SCIPlpiexIsObjlimExc(lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds objlimit: case not handled yet\n");
   }
   else if( SCIPlpiexIsPrimalInfeasible(lpiex) )
   {
      SCIPerrorMessage(" Exact LP infeas.\n");
   }
   else if( SCIPlpiexExistsPrimalRay(lpiex) )
   {
      SCIPerrorMessage("exact LP has primal ray: case not handled yet\n");
   }
   else if( SCIPlpiexIsIterlimExc(lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds iteration limit: case not handled yet\n");
   }
   else if( SCIPlpiexIsTimelimExc(lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds time limit: case not handled yet\n");
   }
   else
   {
      SCIPerrorMessage("lpiex not solved, or other error\n");
   }
   return SCIP_OKAY;
}

/** setup the data for ps in optimal version */
static
SCIP_RETCODE setupPsOpt(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< exact LP data */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Rational**       psobj,              /**< obj function */
   SCIP_Rational**       psub,               /**< upper bounds */
   SCIP_Rational**       pslb,               /**< lower bounds */
   SCIP_Rational**       pslhs,              /**< left hand sides */
   SCIP_Rational**       psrhs,              /**< right hand sides */
   SCIP_Rational**       psval,              /**< matrix values */
   int*                  pslen,              /**< lengths of rows */
   int*                  psind,              /**< indices of vals */
   int*                  psbeg,              /**< beginning of rows in vals */
   int*                  dvarincidence,      /**< indicates cols with bounds */
   int*                  dvarmap,            /**< mapping between red-sized prob and original one */
   SCIP_Rational*        alpha,              /**< scale obj */
   SCIP_Rational*        beta,               /**< scale obj */
   SCIP_Rational*        tmp,                /**< helper rational */
   int                   psnrows,            /**< numer of rows */
   int                   psnnonz,            /**< numer of non-zeros */
   int                   psncols,            /**< numer of columns */
   int                   ndvarmap,           /**< numer of cols with bounds in extended prob */
   int                   nrows,              /**< numer of rows */
   int                   ncols,              /**< numer of cols */
   SCIP_Bool             findintpoint        /**< TRUE, if we search int point, FALSE if we search for ray */
   )
{
   SCIP_PSDATA* psdata;
   int pos;
   int i;
   int j;
   int indx;

   psdata = lpex->psdata;

   /* set up the objective */
   pos = 0;
   for( i = 0; i < nrows; i++ )
   {
      if( dvarincidence[i] )
      {
         RatSet(psobj[pos], lpex->rows[i]->lhs);
         pos++;
      }
   }
   for( i = 0; i < nrows; i++ )
   {
      if( dvarincidence[nrows + i] )
      {
         RatNegate(psobj[pos], lpex->rows[i]->rhs);
         pos++;
      }
   }
   for( i = 0; i < ncols; i++ )
   {
      if( dvarincidence[2*nrows + i] )
      {
         RatSet(psobj[pos], lpex->cols[i]->lb);
         pos++;
      }
   }
   for( i = 0; i < ncols; i++ )
   {
      if( dvarincidence[2*nrows + ncols + i])
      {
         RatNegate(psobj[pos], lpex->cols[i]->ub);
         pos++;
      }
   }
   assert(pos == ndvarmap);

   /* set alpha and beta. */
   RatSetReal(alpha, psdata->psobjweight);
   RatSetInt(beta, 1, 1);

   if( RatIsPositive(alpha) )
   {
      RatDiff(beta, beta, alpha);

      /*  beta = (1-alpha)*|OBJ|   Where OBJ = optimal objective value of root LP, if |OBJ|<1 use 1 instead */
      if( REALABS(SCIPlpGetObjval(lp, set, prob)) > 1 )
      {
         RatSetReal(tmp, REALABS(SCIPlpGetObjval(lp, set, prob)));
         RatMult(beta, beta, tmp);
      }
      /* divide through by alpha and round beta to be a power of 2 */
      RatDiv(beta, beta, alpha);
      RatSetInt(alpha, 1, 1);
      RatSetReal(beta, pow(2, (int) (log(RatApproxReal(beta))/log(2))));
   }

   /* set objective to normalized value */
   for( i = 0; i < ndvarmap; i ++)
      RatMult(psobj[i], psobj[i], alpha);
   RatSet(psobj[ndvarmap], beta);

   /* set variable bounds */
   for( i = 0; i < ndvarmap; i++ )
   {
      RatSetString(psub[i], "inf");
      RatSetInt(pslb[i], 0, 1);
   }
   RatSetInt(psub[ndvarmap], PSBIGM, 1);
   RatSetInt(pslb[ndvarmap], 0 ,1);

   /* set up constraint bounds */
   for( i = 0; i < ncols; i++ )
   {
      RatSet(pslhs[i], lpex->cols[i]->obj);
      RatSet(psrhs[i], lpex->cols[i]->obj);
   }
   for( i = 0; i < psdata->npsbasis; i++ )
   {
      RatSetInt(pslhs[ncols + i], 0, 1);
      RatSetString(psrhs[ncols + i], "inf");
   }

   /* set up constraint matrix: this involves transposing the constraint matrix */

   /* count the length of each constraint */
   for( i = 0; i < psnrows; i++ )
      pslen[i] = 0;
   for( i = 0; i < ndvarmap; i++ )
   {
      indx = dvarmap[i];
      if( indx < 2*nrows )
      {
         if( indx >= nrows )
            indx -= nrows;
         for( j = 0; j < lpex->rows[indx]->len; j++ )
         {
            pslen[lpex->rows[indx]->cols_index[j]]++;
         }
      }
      else
      {
         if ( indx < 2*nrows + ncols )
            indx -= 2 * nrows;
         else
            indx -= (2*nrows + ncols);
         pslen[indx]++;
      }
   }
   for( i = 0; i < psdata->npsbasis; i++ )
   {
      pslen[ncols + i] = 2;
   }
   /* set up the beg array */
   pos = 0;
   for( i = 0; i < psnrows; i++ )
   {
      psbeg[i] = pos;
      pos += pslen[i];
   }
   assert(pos == psnnonz);

   /* reset the length array and build it up as entries are added one by one by scanning through matrix. */
   for( i = 0; i < ncols; i++ )
      pslen[i] = 0;
   for( i = 0; i < ndvarmap; i++ )
   {
      indx = dvarmap[i];
      if( indx < 2*nrows )
      {
         if( indx >= nrows )
            indx -= nrows;
         for(j = 0; j < lpex->rows[indx]->len; j++)
         {
            pos = psbeg[lpex->rows[indx]->cols_index[j]] + pslen[lpex->rows[indx]->cols_index[j]];
            psind[pos] = i;
            if(dvarmap[i]<nrows)
               RatSet(psval[pos], lpex->rows[indx]->vals[j]);
            else
               RatNegate(psval[pos], lpex->rows[indx]->vals[j]);
            pslen[lpex->rows[indx]->cols_index[j]]++;
         }
      }
      else
      {
         if ( indx < 2*nrows + ncols )
            indx -= 2 * nrows;
         else
            indx -= (2*nrows + ncols);
         pos = psbeg[indx] + pslen[indx];
         psind[pos] = i;
         if( dvarmap[i] < 2*nrows + ncols)
            RatSetInt(psval[pos], 1, 1);
         else
            RatSetInt(psval[pos], -1, 1);
         pslen[indx]++;
      }
   }
   /* set up the last npsbasis rows */
   pos = ncols;
   for( i = 0; i < ndvarmap; i++ )
   {
      indx = dvarmap[i];
      if( psdata->includedrows[indx] )
      {
         psind[psbeg[pos]] = i;
         RatSetInt(psval[psbeg[pos]], 1, 1);
         psind[psbeg[pos] + 1] = psncols - 1;
         RatSetInt(psval[psbeg[pos] + 1], -1, 1);
         pos++;
      }
   }
   assert( pos == psnrows);

   if( !findintpoint )
   {
      /* in this case we want to find an interior ray instead of an interior point
      * the problem will be modified to the following problem:
      * max:  [OBJ, 0]*[y,d]'
      * s.t.: [0] <= [ A~ |  0]   [y] <= [  0   ]
      *       [0] <= [ I* | -1] * [d] <= [\infty] <-- only for dual vars from includecons
      * bounds:     0 <= y <= \infty
      *             1 <= d <= \infty
      * y is a vector of length (ndvarmap) and d is a single variable
      * and A~ is the submatrix of [A',-A',I,-I] using columns in dvarmap
      * and OBJ is the subvector of [lhs,-rhs,lb,-ub] using columns in dvarmap
      *
      * the parts that change are the objective function, the RHS/LHS of the first constraint set
      * and the lower bound for d
      */

      RatSetInt(psobj[ndvarmap], 0, 1);

      /* update the rhs/lhs */
      for( i = 0; i < ncols; i++ )
      {
         RatSetInt(pslhs[i], 0, 1);
         RatSetInt(psrhs[i], 0, 1);
      }

      /* update bounds on d */
      RatSetString(psub[ndvarmap], "inf");
      RatSetInt(pslb[ndvarmap], 1 ,1);
   }

   return SCIP_OKAY;
}

/** setup the data for ps in arbitrary-point version */
static
SCIP_RETCODE setupPsArb(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< exact LP data */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Rational**       psobj,              /**< obj function */
   SCIP_Rational**       psub,               /**< upper bounds */
   SCIP_Rational**       pslb,               /**< lower bounds */
   SCIP_Rational**       pslhs,              /**< left hand sides */
   SCIP_Rational**       psrhs,              /**< right hand sides */
   SCIP_Rational**       psval,              /**< matrix values */
   int*                  pslen,              /**< lengths of rows */
   int*                  psind,              /**< indices of vals */
   int*                  psbeg,              /**< beginning of rows in vals */
   int*                  dvarincidence,      /**< indicates cols with bounds */
   int*                  dvarmap,            /**< mapping between red-sized prob and original one */
   SCIP_Rational*        alpha,              /**< scale obj */
   SCIP_Rational*        beta,               /**< scale obj */
   SCIP_Rational*        tmp,                /**< helper rational */
   int                   psnrows,            /**< numer of rows */
   int                   psnnonz,            /**< numer of non-zeros */
   int                   psncols,            /**< numer of columns */
   int                   ndvarmap,           /**< numer of cols with bounds in extended prob */
   int                   nrows,              /**< numer of rows */
   int                   ncols,              /**< numer of cols */
   int                   nobjnz,
   SCIP_Bool             findintpoint        /**< TRUE, if we search int point, FALSE if we search for ray */
   )
{
   SCIP_PSDATA* psdata;
   int pos;
   int i;
   int j;
   int indx;
   int nextendedrows;

   psdata = lpex->psdata;
   nextendedrows = psdata->nextendedrows;

   /* set objective */
      for( i = 0; i < ncols + ndvarmap; i++ )
         RatSetInt( psobj[i ], 0, 1);
      RatSetInt( psobj[ncols + ndvarmap],-1,1);
      for( i = ncols + ndvarmap + 1; i < psncols; i++ )
         RatSetInt( psobj[i], 1, 1);

      /* set variable bounds */
      for( i = 0; i < psncols; i++ )
         RatSetString(psub[i], "inf");

      for( i = 0; i < psncols; i++ )
      {
         if( i < ncols)
            RatSetString( pslb[i], "-inf" );
         else
            RatSetInt( pslb[i], 0, 1);
      }

      /* set up constraint bounds */
      for( i = 0; i < psnrows; i++ )
      {
         if( i < ndvarmap )
         {
            RatSetInt(pslhs[i], 0, 1);
            RatSetInt(psrhs[i], 0, 1);
         }
         else if( i == psnrows - 1 )
         {
            RatSetInt(pslhs[i], 0, 1);
            RatSetString(psrhs[i], "inf");
         }
         else
         {
            RatSetInt(pslhs[i], 1, 1);
            RatSetString(psrhs[i], "inf");
         }
      }

      /* set up constraint matrix */

      /* set up the first ndvarmap rows*/
      pos = 0;
      for( i = 0; i < ndvarmap; i++ )
      {
         indx = dvarmap[i];
         /* current row comes from lhs/rhs constraints of original problem */
         if( indx < 2*nrows )
         {
            if(indx >= nrows)
               indx -= nrows;
            pslen[i] = lpex->rows[indx]->len + 1;
            psbeg[i] = pos;

            /* set A,-A part of row */
            for(j = 0; j < pslen[i] - 1; j++)
            {
               psind[ psbeg[i] + j ] = lpex->rows[indx]->cols_index[j];
               /* if it is an LHS constraint */
               if( dvarmap[i] < nrows )
                  RatSet( psval[ psbeg[i] + j], lpex->rows[indx]->vals[j]);
               /* otherwise it is the RHS of the constraint */
               else
                  RatNegate( psval[ psbeg[i] + j], lpex->rows[indx]->vals[j]);
            }
            /* set I part of row */
            psind[ psbeg[i] + pslen[i] - 1] = ncols + i ;
            RatSetInt( psval[psbeg[i] + pslen[i] - 1 ], -1, 1);

            /* update pos */
            pos += lpex->rows[indx]->len + 1;
         }
         /* current row comes from lower bound constraints of original problem */
         else if ( indx < 2*nrows + ncols )
         {
            indx -= 2 * nrows;
            psbeg[i] = pos;
            pslen[i] = 2;
            psind[pos] = indx;
            psind[pos + 1] = ncols + i;
            RatSetInt(psval[pos], 1, 1);
            RatSetInt(psval[pos+1], -1, 1);
            pos += 2;
         }
         /* current row comes from upper bound constraints of original problem */
         else
         {
            indx -= (2*nrows + ncols);
            psbeg[i] = pos;
            pslen[i] = 2;
            psind[pos] = indx;
            psind[pos + 1] = ncols + i;
            RatSetInt(psval[pos], -1, 1);
            RatSetInt(psval[pos+1], -1, 1);
            pos += 2;
         }
      }

      /* set up the next ndvarmap rows */
      for( i = 0; i < ndvarmap; i++ )
      {
         psbeg[ndvarmap + i] = pos;
         pslen[ndvarmap + i] = 2;
         psind[pos] = ncols + i;
         psind[pos + 1] = ncols + ndvarmap + 1 + i;
         RatSetInt(psval[pos], 1,1);
         RatSetInt(psval[pos+1],1,1);
         pos += 2;
      }

      /* last row */
      psbeg[ psnrows - 1 ] = pos;
      pslen[ psnrows - 1 ] = nobjnz + 1; /* this = objective length + 1 */
      j = psbeg[ 2 * nextendedrows ];
      for( i = 0; i < ncols; i++ )
      {
         if( !RatIsZero(lpex->cols[i]->obj) )
         {
            RatNegate(psval[pos], lpex->cols[i]->obj);
            psind[pos] = i;
            pos++;
         }
      }
      RatSetInt(psval[pos], -1,1);
      psind[pos] = ncols + ndvarmap;
      pos++;
      assert(pos == psnnonz);


   return SCIP_OKAY;
}

/** setup the data for ps in arb-dual version */
static
SCIP_RETCODE setupPsArbDual(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< exact LP data */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Rational**       psobj,              /**< obj function */
   SCIP_Rational**       psub,               /**< upper bounds */
   SCIP_Rational**       pslb,               /**< lower bounds */
   SCIP_Rational**       pslhs,              /**< left hand sides */
   SCIP_Rational**       psrhs,              /**< right hand sides */
   SCIP_Rational**       psval,              /**< matrix values */
   int*                  pslen,              /**< lengths of rows */
   int*                  psind,              /**< indices of vals */
   int*                  psbeg,              /**< beginning of rows in vals */
   int*                  dvarincidence,      /**< indicates cols with bounds */
   int*                  dvarmap,            /**< mapping between red-sized prob and original one */
   SCIP_Rational*        alpha,              /**< scale obj */
   SCIP_Rational*        beta,               /**< scale obj */
   SCIP_Rational*        tmp,                /**< helper rational */
   int                   psnrows,            /**< numer of rows */
   int                   psnnonz,            /**< numer of non-zeros */
   int                   psncols,            /**< numer of columns */
   int                   ndvarmap,           /**< numer of cols with bounds in extended prob */
   int                   nrows,              /**< numer of rows */
   int                   ncols,              /**< numer of cols */
   SCIP_Bool             findintpoint        /**< TRUE, if we search int point, FALSE if we search for ray */
   )
{
   SCIP_PSDATA* psdata;
   int pos;
   int i;
   int j;
   int indx;

   psdata = lpex->psdata;

   /* set up the objective*/
   for( i = 0; i < ndvarmap + 1; i++ )
   {
      RatSetInt(psobj[i], 0, 1);
   }
   for( i = ndvarmap + 1; i < psncols; i++ )
   {
      RatSetInt(psobj[i], -1, 1);
   }

   /* set variable bounds */
   for( i = 0; i < ndvarmap; i++ )
   {
      RatSetString(psub[i], "inf");
      RatSetInt(pslb[i], 0, 1);
   }
   RatSetString(psub[ndvarmap], "inf");
   RatSetInt(pslb[ndvarmap], 1 ,1);
   for( i = ndvarmap + 1; i < psncols; i++ )
   {
      RatSetInt(psub[i], 1, 1);
      RatSetInt(pslb[i], 0, 1);
   }

   /* set up constraint bounds */
   for( i = 0; i < ncols; i++ )
   {
      RatSetInt(pslhs[i], 0, 1);
      RatSetInt(psrhs[i], 0, 1);
   }
   for( i = 0; i < psdata->npsbasis; i++ )
   {
      RatSetInt(pslhs[ncols + i], 0, 1);
      RatSetString(psrhs[ncols + i], "inf");
   }

   /* set up constraint matrix: this involves transposing the constraint matrix */
   SCIPdebugMessage("setting up constraint matrix\n");

   /* count the length of each constraint */
   for( i = 0; i < psnrows; i++ )
      pslen[i] = 0;
   for( i = 0; i < ndvarmap; i++ )
   {
      indx = dvarmap[i];
      if( indx < 2*nrows )
      {
         if( indx >= nrows )
            indx -= nrows;
         for( j = 0 ; j <lpex->rows[indx]->len; j++)
         {
            pslen[lpex->rows[indx]->cols_index[j]]++;
         }
      }
      else
      {
         if ( indx < 2*nrows + ncols )
            indx -= 2 * nrows;
         else
            indx -= (2*nrows + ncols);
         pslen[indx]++;
      }
   }
   for( i = 0; i < psdata->npsbasis; i++ )
   {
      pslen[ncols + i] = 2;
   }

   /* add another element to the first nvar rows for the c vector */
   for( i = 0; i < ncols; i++ )
   {
      pslen[i]++;
   }

   /* set up the beg array */
   pos = 0;
   for( i = 0; i < psnrows; i++ )
   {
      psbeg[i] = pos;
      pos += pslen[i];
   }
   assert(pos == psnnonz);

   /* reset the length array and build it up as entries are added one by one by scanning through matrix. */
   for( i = 0; i < ncols; i++ )
      pslen[i] = 0;
   for( i = 0; i < ndvarmap; i++ )
   {
      indx = dvarmap[i];
      if( indx < 2*nrows )
      {
         if( indx >= nrows )
            indx -= nrows;
         for(j = 0; j < lpex->rows[indx]->len; j++)
         {
            pos = psbeg[lpex->rows[indx]->cols_index[j]] + pslen[lpex->rows[indx]->cols_index[j]];
            psind[pos] = i;
            if(dvarmap[i]<nrows)
               RatSet(psval[pos], lpex->rows[indx]->vals[j]);
            else
               RatNegate(psval[pos], lpex->rows[indx]->vals[j]);
            pslen[lpex->rows[indx]->cols_index[j]]++;
         }
      }
      else
      {
         if ( indx < 2*nrows + ncols )
            indx -= 2 * nrows;
         else
            indx -= (2*nrows + ncols);
         pos = psbeg[indx] + pslen[indx];
         psind[pos] = i;
         if( dvarmap[i] < 2*nrows + ncols)
            RatSetInt(psval[pos], 1, 1);
         else
            RatSetInt(psval[pos], -1, 1);
         pslen[indx]++;
      }
   }
   for( i = 0; i < ncols; i++ )
   {
      RatNegate(psval[psbeg[i] + pslen[i]], lpex->cols[i]->obj);
      psind[psbeg[i] + pslen[i]] = ndvarmap;
      pslen[i]++;
   }

   /* set up the last npsbasis rows */
   pos = ncols;
   for( i = 0; i < ndvarmap; i++ )
   {
      indx = dvarmap[i];
      if( psdata->includedrows[indx] )
      {
         psind[psbeg[pos]] = i;
         RatSetInt(psval[psbeg[pos]], 1, 1);
         psind[psbeg[pos] + 1] = psncols - psnrows + pos;
         RatSetInt(psval[psbeg[pos] + 1], -1, 1);
         pos++;
      }
   }
   assert( pos == psnrows);

   return SCIP_OKAY;
}

/** setup the data for ps in two-stage version */
static
SCIP_RETCODE setupPsTwoStage(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< exact LP data */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Rational**       psobj,              /**< obj function */
   SCIP_Rational**       psub,               /**< upper bounds */
   SCIP_Rational**       pslb,               /**< lower bounds */
   SCIP_Rational**       pslhs,              /**< left hand sides */
   SCIP_Rational**       psrhs,              /**< right hand sides */
   SCIP_Rational**       psval,              /**< matrix values */
   int*                  pslen,              /**< lengths of rows */
   int*                  psind,              /**< indices of vals */
   int*                  psbeg,              /**< beginning of rows in vals */
   int*                  dvarincidence,      /**< indicates cols with bounds */
   int*                  dvarmap,            /**< mapping between red-sized prob and original one */
   SCIP_Rational*        alpha,              /**< scale obj */
   SCIP_Rational*        beta,               /**< scale obj */
   SCIP_Rational*        tmp,                /**< helper rational */
   int                   psnrows,            /**< numer of rows */
   int                   psnnonz,            /**< numer of non-zeros */
   int                   psncols,            /**< numer of columns */
   int                   ndvarmap,           /**< numer of cols with bounds in extended prob */
   int                   nrows,              /**< numer of rows */
   int                   ncols,              /**< numer of cols */
   SCIP_Bool             findintpoint        /**< TRUE, if we search int point, FALSE if we search for ray */
   )
{
   SCIP_PSDATA* psdata;
   int pos;
   int i;
   int j;
   int indx;

   psdata = lpex->psdata;

   /* the representation of the problem will be:
       * max:              [0,1]*[y|d]'
       * s.t.: [c] <= [ A~ |  0]   [y] <= [  c   ]
       *       [0] <= [ I* | -1] * [d] <= [\infty] <-- only for dual vars from includecons
       * bounds:     0 <= y <= \infty
       *             0 <= d <= M
       * y is a vector of length (ndvarmap) and d is a single variable and A~ is the submatrix of [A',-A',I,-I] using
       * columns in dvarmap and OBJ is the subvector of [lhs,-rhs,lb,-ub] using columns in dvarmap
       */


      /* solve the aux problem in two stages, first maximize the interiorness of the point, and then in the second phase
       * problem, move the interiorness to the constraint bounds and optimize over the objective.
       */
   for( i = 0; i < ndvarmap; i ++)
      RatSetInt(psobj[i], 0, 1);
   RatSetInt(psobj[ndvarmap], 1, 1);

   /* set variable bounds */
   for( i = 0; i < ndvarmap; i++ )
   {
      RatSetString(psub[i], "inf");
      RatSetInt(pslb[i], 0, 1);
   }
   RatSetInt(psub[ndvarmap], PSBIGM, 1);
   RatSetInt(pslb[ndvarmap], 0 ,1);

   /* set up constraint bounds */
   for( i = 0; i < ncols; i++ )
   {
      RatSet(pslhs[i], lpex->cols[i]->obj);
      RatSet(psrhs[i], lpex->cols[i]->obj);
   }
   for( i = 0; i < psdata->npsbasis; i++ )
   {
      RatSetInt(pslhs[ncols + i], 0, 1);
      RatSetString(psrhs[ncols + i], "inf");
   }

   /* set up constraint matrix: this involves transposing the constraint matrix */

   /* count the length of each constraint */
   for( i = 0; i < psnrows; i++ )
      pslen[i] = 0;
   for( i = 0; i < ndvarmap; i++ )
   {
      indx = dvarmap[i];
      if( indx < 2*nrows )
      {
         if( indx >= nrows )
            indx -= nrows;
         for(j = 0; j < lpex->rows[indx]->len; j++)
         {
            pslen[lpex->rows[indx]->cols_index[j]]++;
         }
      }
      else
      {
         if ( indx < 2*nrows + ncols )
            indx -= 2 * nrows;
         else
            indx -= (2*nrows + ncols);
         pslen[indx]++;
      }
   }
   for( i = 0; i < psdata->npsbasis; i++ )
   {
      pslen[ncols + i] = 2;
   }

   /* set up the beg array */
   pos = 0;
   for( i = 0; i < psnrows; i++ )
   {
      psbeg[i] = pos;
      pos += pslen[i];
   }
   assert(pos == psnnonz);

   /* reset the length array and build it up as entries are added one by one by scanning through matrix */
   for( i = 0; i < ncols; i++ )
      pslen[i] = 0;
   for( i = 0; i < ndvarmap; i++ )
   {
      indx = dvarmap[i];
      if( indx < 2*nrows )
      {
         if( indx >= nrows )
            indx -= nrows;
         for(j = 0; j < lpex->rows[indx]->len; j++)
         {
            pos = psbeg[lpex->rows[indx]->cols_index[j]] + pslen[lpex->rows[indx]->cols_index[j]];
            psind[pos] = i;
            if(dvarmap[i]<nrows)
               RatSet(psval[pos], lpex->rows[indx]->vals[j]);
            else
               RatNegate(psval[pos], lpex->rows[indx]->vals[j]);
            pslen[lpex->rows[indx]->cols_index[j]]++;
         }
      }
      else
      {
         if ( indx < 2*nrows + ncols )
            indx -= 2 * nrows;
         else
            indx -= (2*nrows + ncols);
         pos = psbeg[indx] + pslen[indx];
         psind[pos] = i;
         if( dvarmap[i] < 2*nrows + ncols)
            RatSetInt(psval[pos], 1, 1);
         else
            RatSetInt(psval[pos], -1, 1);
         pslen[indx]++;
      }
   }

   /* set up the last npsbasis rows */
   pos = ncols;
   for( i = 0; i < ndvarmap; i++ )
   {
      indx = dvarmap[i];
      if( psdata->includedrows[indx] )
      {
         psind[psbeg[pos]] = i;
         RatSetInt(psval[psbeg[pos]], 1, 1);
         psind[psbeg[pos] + 1] = psncols - 1;
         RatSetInt(psval[psbeg[pos] + 1], -1, 1);
         pos++;
      }
   }
   assert( pos == psnrows);

   return SCIP_OKAY;
}

/** subroutine to compute number of nonzeros in lp-matrix */
static
int computePsNnonz(
   SCIP_LPEX*            lpex,               /**< the exact lp */
   int*                  dvarincidence       /**< the columns with existing bounds */
   )
{
   int ret = 0;
   int i;
   int nrows = lpex->nrows;
   int ncols = lpex->ncols;

   for( i = 0; i < nrows; i++ )
   {
      if( dvarincidence[i] )
         ret += lpex->rows[i]->len;
      if( dvarincidence[nrows + i] )
         ret += lpex->rows[i]->len;
   }
   for( i = 0; i < ncols; i++ )
   {
      if( dvarincidence[2*nrows + i] )
         ret++;
      if( dvarincidence[2*nrows + ncols + i] )
         ret++;
   }

   return ret;
}

/** subroutine of constructPSdata(); computes S-interior point or ray which is used to do the shifting step */
static
SCIP_RETCODE psComputeSintPointRay(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< exact LP data */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             findintpoint        /**< TRUE, if we search int point, FALSE if we search for ray */
   )
{
   int i;
   int j;
   int pos;
   int nrows;
   int ncols;
   int nextendedrows; /* number of extended constraints, # of cols in [A',-A',I,-I] */
   int indx;
   int psncols;
   int psnrows;
   int psnnonz;
   int nobjnz;
   SCIP_Rational* tmp;   
   SCIP_Rational* alpha; 
   SCIP_Rational* beta;  

   /* lpiex and data used for the aux. problem */
   SCIP_LPIEX* pslpiex;
   SCIP_LPISTATE* lpistate;
   SCIP_Rational** psobj = NULL;
   SCIP_Rational** pslb = NULL;
   SCIP_Rational** psub = NULL;
   SCIP_PSDATA* psdata;
   SCIP_ROWEX** lprows;
   SCIP_COLEX** lpcols;

   SCIP_Rational** pslhs = NULL;
   SCIP_Rational** psrhs = NULL;

   int* psbeg;
   int* pslen;
   int* psind;
   SCIP_Rational** psval = NULL;
   SCIP_Rational** sol = NULL; /** either primal or dualsol */
   char ** colnames = NULL;
   SCIP_Rational* objval;

   /* mapping between variables used in the aux. problem and the original problem */
   int ndvarmap;
   int* dvarmap;
   int* dvarincidence;

   psdata = lpex->psdata;
   nrows = lpex->nrows;
   ncols = lpex->ncols;

   assert(psdata != NULL);
   nextendedrows = psdata->nextendedrows;

   psncols = 0;
   psnrows = 0;
   psnnonz = 0;
   lprows = lpex->rows;
   lpcols = lpex->cols;

   /* set up dvarmap - mapping between variables and original problem
    * - use the rows that are used for aux. problem
    * - dvarmap[i] is the index in the original problem of the i^th constraint in the reduced size problem
    *   (reduced from nextendedrows to ndvarmap)
    * - dvarincidence gives the incidence vector of variables used in aux problem
    */
   SCIP_CALL( RatCreateBuffer(set->buffer, &tmp) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &alpha) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &beta) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dvarmap, nextendedrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dvarincidence, nextendedrows) );
   {
      /* if the aux. lp is not reduced then expand the selection for dvarmap to include all dual vars with finite cost */
      for( i = 0; i < nextendedrows; i++ )
         dvarincidence[i] = 0;
      for( i = 0; i < nrows; i++ )
      {
         if( !RatIsNegInfinity(lprows[i]->lhs) )
            dvarincidence[i] = 1;
         if( !RatIsInfinity(lprows[i]->rhs) )
            dvarincidence[nrows + i] = 1;
      }
      for( i = 0; i < ncols; i++ )
      {
         if( !RatIsNegInfinity(lpcols[i]->lb) )
            dvarincidence[2*nrows + i] = 1;
         if( !RatIsInfinity(lpcols[i]->ub) )
            dvarincidence[2*nrows + ncols + i] = 1;
      }
   }
   pos = 0;
   for( i = 0; i < nextendedrows; i++ )
   {
      if(dvarincidence[i])
      {
         dvarmap[pos] = i;
         pos++;
      }
   }
   ndvarmap = pos;

   /* if we are finding an interior ray, always use the optimized selection */
   if( psdata->psintpointselection == PS_INTPOINTSEL_OPT || !findintpoint )
   {
      /* in this case we will find an optimized interior point for which we will try to push it interior and
       * optimize over its objective value.  To do this we will solve the following problem
       * max \alpha * [lhs,-rhs,lb,ub] * y + \beta d
       *              s.t. [A,-A,I,-I] * y        = c
       *                                 y_i - d >= 0 for each i \in S
       *                                     y   >= 0
       *                                  M >= d >= 0
       * M is a bound on how interior we will let the point be, S is the set of dual columns chosen earlier
       * which could have nonzero values for the S-interior point.  The parameter psreduceauxlp=TRUE then we
       * exclude all dual variables y_i that are not in S to not be present in this problem.
       *
       * After solving this y will be the S-interior point and d will be the common slack.
       * Here we actually construct the dual in row representation so it can be solved directly.
       */

      psncols =  ndvarmap + 1;
      psnrows = ncols + psdata->npsbasis;
      psnnonz = computePsNnonz(lpex, dvarincidence);
      psnnonz += 2*psdata->npsbasis;

      SCIP_CALL( allocIntMem(set, &psobj, &pslb, &psub, &pslhs, &psrhs, &psval, &sol,
         &objval, &psbeg, &pslen, &psind, &colnames, psncols, psnrows, psnnonz) );

      /* the representation of the problem will be:
       * max:  [\alpha*OBJ, \beta]*[y,d]'
       * s.t.: [c] <= [ A~ |  0]   [y] <= [  c   ]
       *       [0] <= [ I* | -1] * [d] <= [\infty] <-- only for dual vars from includecons
       * bounds:     0 <= y <= \infty
       *             0 <= d <= M
       * y is a vector of length (ndvarmap) and d is a single variable
       * and A~ is the submatrix of [A',-A',I,-I] using columns in dvarmap
       * and OBJ is the subvector of [lhs,-rhs,lb,-ub] using columns in dvarmap
       *
       * beta is set equal to the param psobjweight and alpha is set equal to
       * alpha := (1-beta)/||OBJ||
       */

      SCIP_CALL( setupPsOpt(lp, lpex, set, prob, psobj, psub, pslb, pslhs, psrhs, psval,              /**< matrix values */
         pslen, psind, psbeg, dvarincidence, dvarmap, alpha, beta, tmp, psnrows, psnnonz,            /**< numer of non-zeros */
         psncols, ndvarmap, nrows, ncols, findintpoint) );

      /* build aux LP using the exact LP interface */
      SCIP_CALL( SCIPlpiexCreate(&pslpiex, NULL, "pslpiex", SCIP_OBJSEN_MAXIMIZE) );

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiexAddCols(pslpiex, psncols, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiexAddRows(pslpiex, psnrows, pslhs, psrhs, NULL, psnnonz, psbeg, psind, psval) );

      /* solve the LP */
      SCIP_CALL( SCIPlpiexSolveDual(pslpiex) );

      /* recover the optimal solution and set interior point and slack in constraint handler data */
      if( SCIPlpiexIsOptimal(pslpiex) )
      {
         SCIPdebugMessage("   exact LP solved to optimality\n");
         /* get optimal dual solution */
         SCIP_CALL( SCIPlpiexGetSol(pslpiex, objval, sol, NULL, NULL, NULL) );

         RatSet(psdata->commonslack, sol[psncols - 1]);
         if( RatIsZero(psdata->commonslack) )
         {
            /* if commonslack == 0, point/ray is not interior */
            SCIPdebugMessage("   --> project-and-shift failed to find interior point/ray\n");
         }
         else
         {
            /* assign interior point solution to constraint handler data */
            for( i = 0; i < ndvarmap; i++ )
            {
               if( findintpoint )
                  RatSet( psdata->interiorpt[dvarmap[i]], sol[i]);
               else
                  RatSet( psdata->interiorray[dvarmap[i]], sol[i]);

            }
            if( findintpoint )
               psdata->pshaspoint = TRUE;
            else
               psdata->pshasray = TRUE;
         }
      }
      else
      {
         SCIP_CALL( printlpiexerr( pslpiex ) );
      }
   }
   /* use 'a'rbitrary interior point */
   else if( psdata->psintpointselection == PS_INTPOINTSEL_ARB )
   {
      SCIPdebugMessage("building aux. problem with arbitrary interior point\n");

      /* the aux problem here that we want to solve can be written in the following way.
       * First let A# be the submatrix of [A',-A',I,-I] defined by dvarmap.  Then we want to solve:
       *
       * max   \sum \delta_i
       * s.t.:  A# * y - c*\lambda = 0
       *             y_i >= \delta_i for each i in S
       *               y_i >= 0
       *           1 >= \delta_i >=0
       *                \lambda >= 1
       *
       * solving this problem determines an interior point to the dual problem (which is y/\lambda)
       * it maximizes the number of components which are interior using the \delta_i's.
       *
       * However, instead of solving it in this form, we construct and solve the dual of this problem
       * which can be written like this:
       *
       *   min      [ 0 | 0 |-1 | 1 ] * [x,y,z,w]'
       *   s.t 0 <= [A#'|-I | 0 | 0 ]              <= 0
       *       1 <= [ 0 | I | 0 | I ] * [x,y,z,w]' <= \infty
       *       0 <= [-c'| 0 |-1 | 0 ]              <= \infty
       *             x free, y,z,w >= 0
       *
       * This problem is solved and the dual multipliers for the first set of rows give us the values of y and the
       * next block of rows tell us which components were nonzero (\delta_i) and the final row tells us what the
       * scale factor \lambda of the c in the original problem was.
       */

      /* figure out the dimensions of the aux problem */
      psncols = ncols + 2 * ndvarmap + 1;
      psnrows = 2 * ndvarmap + 1;
      nobjnz = 0;

      /* count the number of nonzeros of the aux problem: psnnonz*/
      for( i = 0; i < ncols; i++ )
      {
         if( !RatIsZero(lpex->cols[i]->obj) )
         {
            nobjnz++;
         }
      }

      psnnonz = computePsNnonz(lpex, dvarincidence);
      psnnonz += nobjnz + 1 + 3 * ndvarmap;

      /* allocate memory for aux problem */
      SCIP_CALL( allocIntMem(set, &psobj, &pslb, &psub, &pslhs, &psrhs, &psval, &sol,
         &objval, &psbeg, &pslen, &psind, &colnames, psncols, psnrows, psnnonz) );

      SCIP_CALL( setupPsArb(lp, lpex, set, prob, psobj, psub, pslb, pslhs, psrhs, psval,              /**< matrix values */
         pslen, psind, psbeg, dvarincidence, dvarmap, alpha, beta, tmp, psnrows, psnnonz,            /**< numer of non-zeros */
         psncols, ndvarmap, nrows, ncols, nobjnz, findintpoint) );

      SCIPdebugMessage("building LPIEX for aux. problem\n");

      /* build aux LP using the exact LP interface */
      SCIP_CALL( SCIPlpiexCreate(&pslpiex, NULL, "lpiexps", SCIP_OBJSEN_MINIMIZE) );

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiexAddCols(pslpiex, psncols, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiexAddRows(pslpiex, psnrows, pslhs, psrhs,              /**< right hand sides */
            NULL, psnnonz, psbeg, psind, psval) );

      /* solve the LP */
      SCIPdebugMessage("solving aux. problem\n");
      SCIP_CALL( SCIPlpiexSolveDual(pslpiex) );

      if( SCIPlpiexIsOptimal(pslpiex) )
      {
         SCIPdebugMessage("   exact LP solved to optimality\n");

         /* compute 1/lambda (lambda is the dual variable corresponding to the last row in the aux LP) */
         SCIP_CALL( SCIPlpiexGetSol(pslpiex, objval, NULL, sol, NULL, NULL) );
         if( !RatIsZero(sol[psnrows - 1]) )
            RatInvert(psdata->commonslack, sol[psnrows - 1]);
         else
         {
            RatSetInt(psdata->commonslack, 0, 1);
         }
         if( RatIsZero(psdata->commonslack) )
         {
            SCIPdebugMessage("   --> project-and-shift did not find S-interior point/ray\n");
         }

         /* interior point is y/lambda */
         psdata->pshaspoint = TRUE;
         for( i = 0; i < ndvarmap; i++ )
         {
            if( psdata->includedrows[dvarmap[i]] && RatIsZero(sol[i]) )
            {
               SCIPdebugMessage("   --> project-and-shift did not find S-interior point/ray\n");
               psdata->pshaspoint = FALSE;
               i = ndvarmap;
            }
            else
               RatDiv( psdata->interiorpt[dvarmap[i]], sol[i], sol[psnrows - 1]);
         }
      }
      else
      {
         SCIP_CALL( printlpiexerr( pslpiex ) );
      }
   }
   /* use 'A'rbitrary interior point in transposed form*/
   else if( psdata->psintpointselection == PS_INTPOINTSEL_ARBDUAL )
   {
      SCIPdebugMessage("building new version of arbitrary interior point aux. problem\n");

      /* the aux problem here that we want to solve can be written in the following way.
       * First let A# be the submatrix of [A',-A',I,-I] defined by dvarmap.  Then we want to solve:
       *
       * max   \sum \delta_i
       * s.t.:  A# * y - c*\lambda = 0
       *             y_i >= \delta_i for each i in S
       *               y_i >= 0
       *           1 >= \delta_i >=0
       *                \lambda >= 1
       *
       *  the representation of the problem will be:
       * min:         [  0 | 0 | -1 ] * [y,z,w]'
       * s.t.: [0] <= [ A~ | -c|  0 ]   [y] <= [  0   ]
       *       [0] <= [ I* | 0 | -I*] * [z] <= [\infty] <-- only for dual vars from includecons
       *                                [w]
       * bounds:     0 <= y <= \infty
       *             1 <= z <= \infty
       *             0 <= w <= 1
       * y is a vector of length (ndvarmap) and d is a single variable
       * and A~ is the submatrix of [A',-A',I,-I] using columns in dvarmap
       */
      psncols =  ndvarmap + 1 + psdata->npsbasis;
      psnrows = ncols + psdata->npsbasis;
      psnnonz = computePsNnonz(lpex, dvarincidence);
      psnnonz += 2*psdata->npsbasis + ncols;

      SCIP_CALL( allocIntMem(set, &psobj, &pslb, &psub, &pslhs, &psrhs, &psval, &sol,
         &objval, &psbeg, &pslen, &psind, &colnames, psncols, psnrows, psnnonz) );

      SCIP_CALL( setupPsArbDual(lp, lpex, set, prob, psobj, psub, pslb, pslhs, psrhs, psval,              /**< matrix values */
         pslen, psind, psbeg, dvarincidence, dvarmap, alpha, beta, tmp, psnrows, psnnonz,            /**< numer of non-zeros */
         psncols, ndvarmap, nrows, ncols, findintpoint) );

      SCIPdebugMessage("building LPIEX for aux. problem\n");

      /* build aux LP using the exact LP interface */
      SCIP_CALL( SCIPlpiexCreate(&pslpiex, NULL, "pslpiex", SCIP_OBJSEN_MINIMIZE) );

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiexAddCols(pslpiex, psncols, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiexAddRows(pslpiex, psnrows, pslhs, psrhs,              /**< right hand sides */
            NULL, psnnonz, psbeg, psind, psval) );

      /* solve the LP */
      SCIPdebugMessage("solving aux. problem\n");
      SCIP_CALL( SCIPlpiexSolveDual(pslpiex) );

      if( SCIPlpiexIsOptimal(pslpiex) )
      {
         SCIPdebugMessage("   exact LP solved to optimality\n");

         /* compute 1/lambda (lambda is the dual variable corresponding to the last row in the aux LP) */
         SCIP_CALL( SCIPlpiexGetSol(pslpiex, objval, sol, NULL, NULL, NULL) );
         if( !RatIsZero(sol[ndvarmap]) )
            RatInvert(psdata->commonslack, sol[ndvarmap]);
         else
         {
            RatSetInt(psdata->commonslack, 0, 1);
         }
         if( RatIsZero(psdata->commonslack) == 0 )
         {
            SCIPdebugMessage("   --> interior point not found\n");
         }

         /* interior point is y/lambda */
         psdata->pshaspoint = TRUE;
         for( i = 0; i < ndvarmap; i++ )
         {
            if( psdata->includedrows[dvarmap[i]] && RatIsZero(sol[i]) )
            {
               psdata->pshaspoint = FALSE;
               SCIPdebugMessage("   --> interior point not found\n");
               i = ndvarmap;
            }
            else
               RatDiv( psdata->interiorpt[dvarmap[i]], sol[i], sol[ndvarmap]);
         }

#ifdef PS_OUT
         printf("constraints all satisfied by slack=");
         RatPrint(psdata->commonslack);
         printf("\n");

         printf("objective value of aux problem=");
         RatPrint(objval);
         printf("\n");

         printf("relative interior solution:\n");
         for( i = 0; i <  psdata->nextendedrows; i++ )
         {
            printf("   i=%d: ", i);
            RatPrint(psdata->interiorpt[i]);
            printf("\n");
         }
#endif
      }
      else
      {
         printlpiexerr( pslpiex );
         psdata->psdatafail = TRUE;
      }
   }
   else if( psdata->psintpointselection == PS_INTPOINTSEL_TWOSTAGE )
   {
      /* in this case we will find an optimized interior point for which we will try to push it interior and
       * optimize over its objective value.  To do this we will solve the following two problems
       * max                                   d
       *              s.t. [A,-A,I,-I] * y        = c
       *                                 y_i - d >= 0 for each i \in S
       *                                     y   >= 0
       *                                  M >= d >= 0
       * max          [lhs,-rhs,lb,ub] * y
       *              s.t. [A,-A,I,-I] * y        = c
       *                                 y_i - d >= 0 for each i \in S
       *                                     y   >= 0
       *                                       d >= d* <-- where d* is the optimal solution from the first problem
       * M is a bound on how interior we will let the point be, S is the set of dual columns chosen earlier
       * which could have nonzero values for the S-interior point.  The parameter psreduceauxlp=TRUE then we
       * exclude all dual variables y_i that are not in S to not be present in this problem.
       *
       * After solving this y will be the S-interior point and d will be the common slack.
       * Here we actually construct the dual in row representation so it can be solved directly.
       */
      psncols =  ndvarmap + 1;
      psnrows = ncols + psdata->npsbasis;
      psnnonz = computePsNnonz(lpex, dvarincidence);
      psnnonz += 2*psdata->npsbasis;

      SCIP_CALL( allocIntMem(set, &psobj, &pslb, &psub, &pslhs, &psrhs, &psval, &sol,
         &objval, &psbeg, &pslen, &psind, &colnames, psncols, psnrows, psnnonz) );

      SCIP_CALL( setupPsTwoStage(lp, lpex, set, prob, psobj, psub, pslb, pslhs, psrhs, psval,              /**< matrix values */
         pslen, psind, psbeg, dvarincidence, dvarmap, alpha, beta, tmp, psnrows, psnnonz,            /**< numer of non-zeros */
         psncols, ndvarmap, nrows, ncols, findintpoint) );

      /* build aux LP using the exact LP interface */
      if( pslpiex != NULL )
      {
         SCIP_CALL( SCIPlpiexFree(&pslpiex) );
      }
      pslpiex = NULL;
      SCIP_CALL( SCIPlpiexCreate(&pslpiex, NULL, "pslpiex", SCIP_OBJSEN_MAXIMIZE) );

#ifdef PS_OUT
      /* activate extra output from exact LP solver */
      SCIPlpiexSetIntpar(pslpiex, SCIP_LPPAR_LPINFO, TRUE);
#endif

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiexAddCols(pslpiex, psncols, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiexAddRows(pslpiex, psnrows, pslhs, psrhs,              /**< right hand sides */
            NULL, psnnonz, psbeg, psind, psval) );

      /* solve the LP */
      SCIP_CALL( SCIPlpiexSolveDual(pslpiex) );

      /* get state and solution of lpiex that was just solved */
      SCIP_CALL( SCIPlpiexGetState(pslpiex, blkmem, &lpistate) );
      SCIP_CALL( SCIPlpiexGetSol(pslpiex, objval, NULL, NULL, NULL, NULL) );

      /* now reset the objective value to be the original objective */
      pos = 0;
      for( i = 0; i < nrows; i++ )
      {
         if( dvarincidence[i] )
         {
            RatSet(psobj[pos], lpex->rows[i]->lhs);
            pos++;
         }
      }
      for( i = 0; i < nrows; i++ )
      {
         if( dvarincidence[nrows + i] )
         {
            RatNegate(psobj[pos], lpex->rows[i]->rhs);
            pos++;
         }
      }
      for( i = 0; i < ncols; i++ )
      {
         if( dvarincidence[2*nrows + i] )
         {
            RatSet(psobj[pos], lpex->cols[i]->lb);
            pos++;
         }
      }
      for( i = 0; i < ncols; i++ )
      {
         if( dvarincidence[2*nrows + ncols + i])
         {
            RatNegate(psobj[pos], lpex->cols[i]->ub);
            pos++;
         }
      }
      assert(pos == ndvarmap);
      RatSetInt(psobj[ndvarmap], 0, 1);

      /* set the lower bound on the interiorness based on the objective value */
      RatSet(pslb[ndvarmap], objval);

      /* reuse the psind array to pass indices to updated the bounds and objective */
      for( i = 0; i < psncols; i++ )
         psind[i] = i;
      SCIP_CALL( SCIPlpiexChgBounds(pslpiex, psncols, psind, pslb, NULL) );
      SCIP_CALL( SCIPlpiexChgObj(pslpiex, psncols, psind, psobj) );

      /* reload state and solve new LP */
      SCIP_CALL( SCIPlpiexSetState(pslpiex, blkmem, lpistate) );

      /* reoptimizing using primal simplex seems MUCH faster here, warm start basis is primal feasible */
      SCIP_CALL( SCIPlpiexSolvePrimal(pslpiex) );
      SCIP_CALL( SCIPlpiexFreeState(pslpiex, blkmem, &lpistate) );

        /* recover the optimal solution and set interior point and slack in constraint handler data */
      if( SCIPlpiexIsOptimal(pslpiex) )
      {
         SCIPdebugMessage("   exact LP solved to optimality\n");

         /* get optimal dual solution */
         SCIP_CALL( SCIPlpiexGetSol(pslpiex, objval, sol, NULL, NULL, NULL) );

         /* assign interior point solution to constraint handler data */
         RatSet(psdata->commonslack, sol[psncols - 1]);
         if( RatIsZero(psdata->commonslack) )
         {
            SCIPdebugMessage("   --> interior point not found \n");
         }
         else
         {
            for( i = 0; i < ndvarmap; i++ )
            {
               RatSet( psdata->interiorpt[dvarmap[i]], sol[i]);
            }
            psdata->pshaspoint = TRUE;
         }
      }
      else
      {
         SCIP_CALL( printlpiexerr(pslpiex) );
      }
   }
   else
   {
      SCIPerrorMessage("invalid parameter setting <%d> for selection method to compute interior point\n",
         psdata->psintpointselection);
      return SCIP_PARAMETERWRONGVAL;
   }
   
   for( i = 0; i < ndvarmap; i++ )
   {
      if( psdata->pshaspoint )
         RatCanonicalize(psdata->interiorpt[i]);
      if( psdata->pshasray )
         RatCanonicalize(psdata->interiorray[i]);
   }

   /* free memory */
   if( pslpiex != NULL )
   {
      int nlpirows, nlpicols;
      SCIP_CALL( SCIPlpiexGetNRows(pslpiex, &nlpirows) );
      SCIP_CALL( SCIPlpiexDelRows(pslpiex, 0, nlpirows - 1) );

      SCIP_CALL( SCIPlpiexGetNCols(pslpiex, &nlpicols) );
      SCIP_CALL( SCIPlpiexDelCols(pslpiex, 0, nlpicols - 1) );

      SCIP_CALL( SCIPlpiexClear(pslpiex) );
      SCIP_CALL( SCIPlpiexFree(&pslpiex) );
   }
   assert(pslpiex == NULL);
   for( i = psncols - 1; i >= 0; i--)
      SCIPsetFreeBufferArray(set, &colnames[i] );
   SCIPsetFreeBufferArray(set, &colnames);

   SCIPsetFreeBufferArray(set, &psind);
   SCIPsetFreeBufferArray(set, &pslen);
   SCIPsetFreeBufferArray(set, &psbeg);

   RatFreeBuffer(set->buffer, &objval);
   RatFreeBufferArray(set->buffer, &sol, psncols);
   RatFreeBufferArray(set->buffer, &psval, psnnonz);
   RatFreeBufferArray(set->buffer, &psrhs, psnrows);
   RatFreeBufferArray(set->buffer, &pslhs, psnrows);
   RatFreeBufferArray(set->buffer, &psub, psncols);
   RatFreeBufferArray(set->buffer, &pslb, psncols);
   RatFreeBufferArray(set->buffer, &psobj, psncols);

   SCIPsetFreeBufferArray(set, &dvarincidence);
   SCIPsetFreeBufferArray(set, &dvarmap);

   RatFreeBuffer(set->buffer, &beta);
   RatFreeBuffer(set->buffer, &alpha);
   RatFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** constructs datas used to compute dual bounds by the project-and-shift method */
static
SCIP_RETCODE constructPSData(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< exact LP data */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_STAT*            stat,               /**< statistics pointer */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem
   )
{
   int i;
   SCIP_PSDATA* psdata;

   assert(lpex != NULL);
   assert(lpex->psdata != NULL);

   psdata = lpex->psdata;

   /* consider the primal problem as
    * min c'x
    * lhs <= Ax <= rhs
    *  lb <=  x <= ub
    *
    * and the dual of the form
    *     [  A',  -A',  I,  -I] y =  c
    *                           y >= 0
    *
    * A subset S of the dual columns are chosen to give a submatrix D of [A',-A',I,-I], which is then LU factorized using
    * rectlu code then an S-interior point is found (a dual solution that is strictly positive for each column in S).
    * this data is then reused throughout the tree where the LU factorization can be used to correct feasibility of
    * the equality constraints of the dual, and a convex combination with the S-interior point can correct any
    * infeasibility coming from negative variables.
    */

   /* if the ps data was already constructed, exit */
   if( psdata->psdatacon )
      return SCIP_OKAY;

   SCIPclockStart(stat->provedfeaspstime, set);
   /* now mark that this function has been called */
   psdata->psdatacon = TRUE;

   SCIPdebugMessage("calling constructPSdata()\n");
   SCIP_CALL( RatCreateBlock(blkmem, &psdata->commonslack) );

   /* process the bound changes */
   SCIP_CALL( SCIPsepastoreexSyncLPs(set->scip->sepastoreex, blkmem, set, stat, lpex, eventqueue, eventfilter) );
   SCIP_CALL( SCIPlpexFlush(lp->lpex, blkmem, set, eventqueue) );

   assert(lpex->nrows > 0);

   psdata->nextendedrows = 2*lpex->nrows + 2*lpex->ncols;

   /* call function to select the set S */
   SCIP_CALL( psChooseS(lp, lpex, set, stat, messagehdlr, eventqueue, eventfilter, prob, blkmem) );

   /* compute LU factorization of D == A|_S */
   SCIP_CALL( psFactorizeD(lp, lpex, set, prob, blkmem, psdata->psuseintpoint) );

   /* if no fail in LU factorization, compute S-interior point and/or ray */
   if( !psdata->psdatafail )
   {
      if( TRUE || !psdata->psuseintpoint )
      {
         /* try to compute the S-interior ray if we want to use it for bounding or infeasibility */
         SCIP_CALL( RatCreateBlockArray(blkmem, &psdata->interiorray, psdata->nextendedrows) );
         SCIP_CALL( psComputeSintPointRay(lp, lpex, set, prob, blkmem, FALSE) );
      }
      if( psdata->psuseintpoint || !psdata->pshasray )
      {
         /* now, compute S-interior point if we need it OR if the ray construction failed */
         SCIP_CALL( RatCreateBlockArray(blkmem, &psdata->interiorpt, psdata->nextendedrows) );
         SCIP_CALL( psComputeSintPointRay(lp, lpex, set, prob, blkmem, TRUE) );
      }
   }

   /* if construction of both point and ray has failed, mark psdatafail as true. */
   if( !psdata->pshaspoint && !psdata->pshasray )
   {
      psdata->psdatafail = TRUE;
   }

   SCIP_CALL( RatCreateBlockArray(blkmem, &psdata->violation, lpex->ncols) );
   SCIP_CALL( RatCreateBlockArray(blkmem, &psdata->correction, psdata->nextendedrows) );
   psdata->violationsize = lpex->ncols;

   SCIPclockStop(stat->provedfeaspstime, set);
   SCIPdebugMessage("exiting constructPSdata()\n");

   return SCIP_OKAY;
}

/** computes safe dual bound by project-and-shift method or corrects dual ray for infeasibility proof */
static
SCIP_RETCODE getPSdualTwo(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< exact LP data */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_STAT*            stat,               /**< statistics pointer */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             usefarkas           /**< to we aim to prove infeasibility ? */
   )
{
   SCIP_COL** cols;
   SCIP_Rational** approxdual;
   SCIP_Rational** costvect;                 /**< dual cost vector for expanded problem*/
   SCIP_Rational** violation;
   SCIP_Rational** correction;
   SCIP_Bool useinteriorpoint;
   SCIP_Rational* tmp;
   SCIP_Rational* tmp2;
   SCIP_Rational* lambda1;
   SCIP_Rational* lambda2;
   SCIP_Rational* maxv;
   SCIP_Rational* dualbound;
   SCIP_PSDATA* psdata;
   mpq_t* violationgmp = NULL;
   mpq_t* correctiongmp = NULL;
   int i;
   int j;
   int rval;
   int nextendedrows;
   int nrows;
   int nrowsps; /* number of rows used for ps. this can be lower than nrows, if e.g. cuts were added */
   int ncols;
   int currentrow;
   SCIP_Bool isfeas;

   /* project-and-shift method:
    * 1. projection step (to ensure that equalities are satisfied):
    *   - compute error in equalities: r=c-Ay^
    *   - backsolve system of equations to find correction of error: z with Dz=r
    *   - add corretion to approximate dual solution: bold(y)=y^+[z 0]
    * 2. shifing step (to ensure that inequalities are satisfied):
    *   - take convex combination of projected approximate point bold(y) with interior point y*
    * 3. compute dual objective value of feasible dual solution and set bound
    */

   /* constructpsdata should always be called first */
   if( usefarkas )
   {
      stat->nprojshift++;
      SCIPclockStart(stat->provedinfeaspstime, set);
   }
   else
   {
      stat->nprojshiftinf++;
      SCIPclockStart(stat->provedfeaspstime, set);
   }
   psdata = lpex->psdata;

   assert(psdata->psdatacon);

   /* if data has not been constructed, or it failed, then exit */
   if( (psdata->psdatafail && !usefarkas) || (usefarkas && !psdata->pshasray) )
   {
      lp->hasprovedbound = FALSE;
      return SCIP_OKAY;
   }

   lp->hasprovedbound = TRUE;

   SCIPdebugMessage("calling getPSdualbound()\n");

   /* decide if we should use ray or point to compute bound */
   if( !usefarkas && psdata->psuseintpoint && psdata->pshaspoint )
   {
      /* if we are supposed to use the interior point, and it exists use it */
      useinteriorpoint = TRUE;
   }
   else
   {
      /* in this case, since psdatafail != TRUE, pshasray should be true -- use it */
      assert( psdata->pshasray );
      useinteriorpoint = FALSE;
   }

   SCIP_CALL( RatCreateBuffer(set->buffer, &tmp) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &tmp2) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &lambda1) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &lambda2) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &maxv) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &dualbound) );

   /* flush exact lp */
   /* set up the exact lpi for the current node */
   SCIP_CALL( SCIPsepastoreexSyncLPs(set->scip->sepastoreex, blkmem, set, stat, lpex, eventqueue, eventfilter) );
   SCIP_CALL( SCIPlpexFlush(lp->lpex, blkmem, set, eventqueue) );

   nextendedrows = psdata->nextendedrows;
   nrows = lpex->nrows;
   ncols = lpex->ncols;
   nrowsps = psdata->nextendedrows/2 - ncols;

   /* allocate memory for approximate dual solution, dual cost vector, violation and correction */
   SCIP_CALL( RatCreateBufferArray(set->buffer, &approxdual, nextendedrows) );
   SCIP_CALL( RatCreateBufferArray(set->buffer, &costvect, nextendedrows) );
   SCIP_CALL( RatCreateBufferArray(set->buffer, &violation, ncols) );
   SCIP_CALL( RatCreateBufferArray(set->buffer, &correction, nextendedrows) );

   /* recover the objective coefs and approximate solution value of dual solution;
    * dual vars of lhs constraints (including -inf) and rhs constraints (including +inf),
    * dual vars of lb constraint (including -inf) and ub constraints (including +inf)
    */
   for( i = 0; i < nrowsps; i++ )
   {
      /* in case we want to prove infeasibility it might be that we were not able to compute a dual solution
       * with bound exceeding objective value; in this case dual solution is set to SCIP_INVALID
       */
      if( (!usefarkas && SCIProwGetDualsol(lp->rows[i]) == SCIP_INVALID) || (usefarkas && SCIProwGetDualfarkas(lp->rows[i]) == SCIP_INVALID) )
      {
         SCIPdebugMessage("  no valid unbounded approx dual sol given\n");
         lp->hasprovedbound = FALSE;
         if( usefarkas )
            stat->nfailprojshiftinf++;
         else
            stat->nfailprojshift++;

         goto TERMINATE;
      }

      if( !usefarkas )
         RatSetReal(tmp, SCIProwGetDualsol(lp->rows[i]));
      else
         RatSetReal(tmp, SCIProwGetDualfarkas(lp->rows[i]));

      /* lhs constraint */
      if( RatIsPositive(tmp) )
         RatSet(approxdual[i], tmp);
      /* rhs constraint */
      else
         RatNegate(approxdual[i+nrowsps], tmp);

      RatSet(costvect[i], lpex->rows[i]->lhs);
      RatNegate(costvect[i+nrowsps], lpex->rows[i]->rhs);
   }

   cols = SCIPlpGetCols(lp);
   for( i = 0; i < ncols; i++ )
   {
      if( !usefarkas )
         RatSetReal(tmp, SCIPcolGetRedcost(cols[i], stat, lp));
      else
         RatSetReal(tmp, -SCIPcolGetFarkasCoef(cols[i], stat, lp));
      /* lb constraint */
      if( RatIsPositive(tmp) )
         RatSet(approxdual[i + 2*nrowsps], tmp);
      /* ub constraint */
      else
         RatNegate(approxdual[i + 2*nrowsps + ncols], tmp);

      RatSet(costvect[i + 2*nrowsps], lpex->cols[i]->lb);
      RatNegate(costvect[i + 2*nrowsps + ncols], lpex->cols[i]->ub);
   }

#ifdef PS_OUT
   printf("approximate dual solution:\n");

   RatSetInt(dualbound, 0, 1);
   for( i = 0; i < nextendedrows; i++ )
   {
      printf("   i=%d: ", i);
      RatPrint(approxdual[i]);
      printf(" * ");
      RatPrint(costvect[i]);
      printf("\n");
      if( RatIsAbsInfinity(costvect[i]) )
       assert(RatIsZero(approxdual[i]));
      else
      {
         RatAddProd(dualbound, approxdual[i], costvect[i]);
      }
   }

   printf("   objective value=%.20f (", RgetRealApprox(dualbound));
   RatPrint(dualbound);
   printf(")\n");
#endif

   /* first, ensure nonnegativity of dual solution and fix artificial dual variables to zero */
   for( i = 0; i < nextendedrows; i++ )
   {
      if( RatIsNegative(approxdual[i]) )
         RatSetInt(approxdual[i], 0, 1);
      else if( RatIsNegInfinity(costvect[i]) )
      {
         if( i < nextendedrows )
            assert(!psdata->includedrows[i]);

         RatSetInt(approxdual[i], 0, 1);
      }
   }

   /* calculate violation of equality constraints r=c-A^ty */
   for( i = 0; i < ncols; i++ )
   {
      if( !usefarkas )
         RatSet(violation[i], lpex->cols[i]->obj);
      else
         RatSetInt(violation[i], 0, 1);
   }

   /* A^ty for y corresponding to primal constraints */
   for( i = 0; i < nrowsps; i++ )
   {
      for( j = 0; j < lpex->rows[i]->len; j++)
      {
         currentrow = lpex->rows[i]->cols_index[j];
         /* avoid expensive null-sum */
         if( RatIsEqual(approxdual[i], approxdual[i+nrowsps]) )
            continue;

         RatDiffProd(violation[currentrow], approxdual[i], lpex->rows[i]->vals[j]);
         RatAddProd(violation[currentrow], approxdual[i+nrowsps], lpex->rows[i]->vals[j]);
      }
   }

   /* add violations coming from non-initial rows */
   for( i = nrowsps; i < nrows; ++i )
   {
      for( j = 0; j < lpex->rows[i]->len; j++)
      {
         if( !usefarkas )
            RatSetReal(tmp, SCIProwGetDualsol(lp->rows[i]));
         else
            RatSetReal(tmp, SCIProwGetDualfarkas(lp->rows[i]));

         currentrow = lpex->rows[i]->cols_index[j];
         RatDiffProd(violation[currentrow], tmp, lpex->rows[i]->vals[j]);
      }
   }

   /* A^ty for y corresponding to bound constraints */
   for( i = 0; i < ncols; i++ )
   {
      /* avoid expensive null-sum */
      if( RatIsEqual(approxdual[i+2*nrowsps], approxdual[i+2*nrowsps+ncols]) )
            continue;
      RatDiff(violation[i], violation[i], approxdual[i+2*nrowsps]);
      RatAdd(violation[i], violation[i], approxdual[i+2*nrowsps+ncols]);
   }


   /* project solution */

#ifdef PS_OUT
   printf("violation of solution:\n");
   for( i = 0; i < ncols; i++ )
   {
      printf("   i=%d: ", i);
      RatPrint(violation[i]);
      printf("\n");
   }
#endif

   /* if there is no violation of the constraints, then skip the projection */
   isfeas = TRUE;
   for( i = 0; i < ncols && isfeas; i++ )
   {
      if( !RatIsZero(violation[i]) )
         isfeas = FALSE;
   }

   /* isfeas is equal to one only if approximate dual solution is already feasible for the dual */
   if( !isfeas )
   {
      /* compute [z] with Dz=r (D depends on psdata->psdualcolselection) */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &violationgmp, ncols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &correctiongmp, nextendedrows) );

      RatSetGMPArray(violationgmp, violation, ncols);
      RatSetGMPArray(correctiongmp, correction, nextendedrows);
      rval = RECTLUsolveSystem( psdata->rectfactor, ncols, nextendedrows, violationgmp, correctiongmp);
      if( rval )
      {
         lp->hasprovedbound = FALSE;
         if( usefarkas )
            stat->nfailprojshiftinf++;
         else
            stat->nfailprojshift++;

         goto TERMINATE;
      }
      RatSetArrayGMP(correction, correctiongmp, nextendedrows);

#ifdef PS_OUT
      printf("correction of solution:\n");
      for( i = 0; i < psdata->npsbasis; i++ )
      {
         printf("   i=%d: ", i);
         RatPrint(correction[i]);
         printf(", position=%d\n", psdata->psbasis[i]);
      }
#endif

      /* projection step: compute bold(y)=y^+[z 0];
       * correct only components corresponding to D (npsbasis=#of columns in D
       */
      for( i = 0; i < psdata->npsbasis; i++ )
      {
         RatAdd(approxdual[psdata->psbasis[i]], approxdual[psdata->psbasis[i]], correction[i]);
      }
   }

#ifdef PS_OUT
   printf("updated dual solution:\n");
   for( i = 0; i < nextendedrows; i++ )
   {
      printf("   i=%d: ", i);
      RatPrint(approxdual[i]);
      printf("\n");
   }
#endif

   if( useinteriorpoint )
   {
      assert(!usefarkas);
      /* shifting step (scale solution with interior point to be dual feasible):
       * y' = lambda1 bold(y) + lambda2 y*, where
       *   lambda1 = ( slack of int point)/ (slack of int point + max violation) = d/m+d
       *   lambda2 = 1 - lambda1
       */

      /* compute lambda1 componentwise (set lambda1 = 1 and lower it if necessary) */
      RatSetInt(lambda1, 1, 1);
      for( i = 0; i < nextendedrows; i++ )
      {
         if( RatIsNegative(approxdual[i]) )
         {
            RatSet(tmp2, psdata->interiorpt[i]);
            RatDiff(tmp, psdata->interiorpt[i], approxdual[i]);
            RatDiv(tmp2, tmp2, tmp);
            if( RatIsGT(lambda1, tmp2) )
               RatSet(lambda1, tmp2);
         }
      }
      RatSetInt(lambda2, 1, 1);
      RatDiff(lambda2, lambda2, lambda1);
   }
   else
   {
      /* in this case we are using an interior ray that can be added freely to the solution */
      RatSetInt(lambda1, 1, 1);

      /* compute lambda values: compute lambda1 componentwise (set lambda1 = 1 and lower it if necessary) */
      RatSetInt(lambda1, 1, 1);
      for( i = 0; i < nextendedrows; i++ )
      {
         if( RatIsNegative(approxdual[i]) && psdata->includedrows[i] )
         {
            RatDiv(tmp, approxdual[i], psdata->interiorray[i]);
            RatNegate(tmp, tmp);
            if( RatIsLT(lambda2, tmp) ) /* changed to < */
               RatSet(lambda2, tmp);
         }
      }
   }

   /* perform shift */
   if( !RatIsZero(lambda2) )
   {
      for( i = 0; i < nextendedrows; i++ )
      {
         RatMult(approxdual[i], approxdual[i], lambda1);
      }
      for( i =0; i < nextendedrows; i++ )
      {
         RatAddProd(approxdual[i], useinteriorpoint ? psdata->interiorpt[i] : psdata->interiorray[i], lambda2);
      }
   }

   /* postprocess dual solution to reduce values when both sides of constraint used;
    * if y(lhs) and y(rhs) are both nonzero shift them such that one becomes zero
    * this will tighten the solution and improve the objective value, there is no way this can hurt
    */
      /* y(lhs) and y(rhs) corresponding to primal constraints */
   for( i = 0; i < nrowsps; i++ )
   {
      /* find the min value of y(lhs) and y(rhs) */
      RatMIN(tmp, approxdual[i], approxdual[i+nrowsps]);

      /* shift if both are nonzero */
      if( RatIsPositive(tmp) )
      {
         RatDiff(approxdual[i], approxdual[i], tmp);
         RatDiff(approxdual[i+nrowsps], approxdual[i+nrowsps], tmp);
      }
   }
   /* y(lhs) and y(rhs) corresponding to bound constraints */
   for( i = 0; i < ncols; i++ )
   {
      /* find the min value of y(lhs) and y(rhs) */
      RatMIN(tmp, approxdual[i + 2*nrowsps], approxdual[i + 2*nrowsps + ncols]);

      /* shift if both are nonzero */
      if( RatIsPositive(tmp) )
      {
         RatDiff(approxdual[i+2*nrowsps], approxdual[i+2*nrowsps], tmp);
         RatDiff(approxdual[i+2*nrowsps + ncols], approxdual[i+2*nrowsps + ncols], tmp);
      }
   }

#ifdef PS_OUT
   printf("projected and shifted dual solution (should be an exact dual feasible solution)\n");
   for( i = 0; i < nextendedrows; i++ )
   {
      printf("   i=%d: ", i);
      RatPrint(approxdual[i]);
      printf("\n");
   }
#endif


#ifndef NDEBUG
   SCIPdebugMessage("debug test: verifying feasibility of dual solution:\n");

   /* calculate violation of equality constraints: subtract Ax to get violation b-Ax, subtract A(approxdual) */
   rval = 0;
   for( i = 0; i < ncols; i++ )
   {
      if( !usefarkas )
        RatSet(violation[i], lpex->cols[i]->obj);
      else
         RatSetInt(violation[i], 0, 1);
   }
   for( i = 0; i < nrowsps; i++ )
   {
      for( j = 0; j < lpex->rows[i]->len; j++ )
      {
         currentrow = lpex->rows[i]->cols_index[j];
         /* avoid expensive null-sum */
         if( RatIsEqual(approxdual[i], approxdual[i+nrowsps]) )
            continue;
         RatDiffProd(violation[currentrow], approxdual[i], lpex->rows[i]->vals[j]);
         RatAddProd(violation[currentrow], approxdual[i+nrowsps], lpex->rows[i]->vals[j]);
      }
   }
   for( i = 0; i < ncols; i++ )
   {
      /* avoid expensive null-sum */
      if( RatIsEqual(approxdual[i+2*nrowsps], approxdual[i+2*nrowsps+ncols]) )
         continue;
      RatDiff(violation[i], violation[i], approxdual[i+2*nrowsps]);
      RatAdd(violation[i], violation[i], approxdual[i+2*nrowsps+ncols]);
   }
   /* add violations coming from non-initial rows */
   for( i = nrowsps; i < nrows; ++i )
   {
      for( j = 0; j < lpex->rows[i]->len; j++)
      {
         if( !usefarkas )
            RatSetReal(tmp, SCIProwGetDualsol(lp->rows[i]));
         else
            RatSetReal(tmp, SCIProwGetDualfarkas(lp->rows[i]));

         currentrow = lpex->rows[i]->cols_index[j];
         RatDiffProd(violation[currentrow], tmp, lpex->rows[i]->vals[j]);
      }
   }

   for( i = 0; i < ncols; i++ )
   {
      if( !RatIsZero(violation[i]) )
      {
         SCIPdebugMessage("   dual solution incorrect, violates equalties\n");
         i = ncols;
         rval = 1;
      }
   }
   for( i = 0; i < nextendedrows; i++ )
   {
      if( RatIsNegative(approxdual[i]) )
      {
         SCIPdebugMessage("   dual solution incorrect, negative components\n");
         i = nextendedrows;
         rval = 1;
      }
   }
   if( !rval )
      SCIPdebugMessage("   dual solution verified\n");
   assert(!rval);
#endif

   RatSetInt(dualbound, 0, 1);
   for( i = 0; i < nextendedrows; i++ )
   {
      RatAddProd(dualbound, approxdual[i], costvect[i]);
   }

   if( !usefarkas )
   {
      lp->lpobjval = RatRoundReal(dualbound, SCIP_ROUND_DOWNWARDS);
      lp->hasprovedbound = TRUE;
   }
   else
   {
      /* if the objective value of the corrected ray is positive we can prune node, otherwise not */
      if( RatIsPositive(dualbound) )
         lp->hasprovedbound = TRUE;
      else
      {
         stat->nfailprojshiftinf++;
         lp->hasprovedbound = FALSE;
      }
   }

#ifdef PS_OUT
   printf("   common slack=%.20f (", RgetRealApprox(psdata->commonslack));
   RatPrint(psdata->commonslack);
   printf(")\n");

   printf("   max violation=%.20f (", RgetRealApprox(maxv));
   RatPrint(maxv);
   printf(")\n");

   printf("   lambda (use of interior point)=%.20f (", RgetRealApprox(lambda2));
   RatPrint(lambda2);
   printf(")\n");

   printf("   dual objective value=%.20f (", RgetRealApprox(dualbound));
   RatPrint(dualbound);
   printf(")\n");
#endif

 TERMINATE:
   /* free memory */
   if( correctiongmp != NULL )
   {
      RatClearGMPArray(correctiongmp, nextendedrows);
      SCIPsetFreeBufferArray(set, &correctiongmp);
   }
   if( violationgmp != NULL )
   {
      RatClearGMPArray(violationgmp, ncols);
      SCIPsetFreeBufferArray(set, &violationgmp);
   }

   RatFreeBufferArray(set->buffer, &correction, nextendedrows);
   RatFreeBufferArray(set->buffer, &violation, ncols);
   RatFreeBufferArray(set->buffer, &costvect, nextendedrows);
   RatFreeBufferArray(set->buffer, &approxdual, nextendedrows);

   RatFreeBuffer(set->buffer, &dualbound);
   RatFreeBuffer(set->buffer, &maxv);
   RatFreeBuffer(set->buffer, &lambda2);
   RatFreeBuffer(set->buffer, &lambda1);
   RatFreeBuffer(set->buffer, &tmp2);
   RatFreeBuffer(set->buffer, &tmp);

   if( usefarkas )
      SCIPclockStop(stat->provedinfeaspstime, set);
   else
      SCIPclockStop(stat->provedfeaspstime, set);

   return SCIP_OKAY;
}

/** computes safe dual bound by project-and-shift method or corrects dual ray for infeasibility proof */
static
SCIP_RETCODE getPSdual(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< exact LP data */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_STAT*            stat,               /**< statistics pointer */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             usefarkas           /**< to we aim to prove infeasibility ? */
   )
{
   SCIP_COL** cols;
   SCIP_Rational** approxdual;
   SCIP_Rational** violation;
   SCIP_Rational** correction;
   SCIP_Bool useinteriorpoint;
   SCIP_Rational* tmp;
   SCIP_Rational* tmp2;
   SCIP_Rational* lambda1;
   SCIP_Rational* lambda2;
   SCIP_Rational* maxv;
   SCIP_Rational* dualbound;
   SCIP_PSDATA* psdata;
   mpq_t* violationgmp = NULL;
   mpq_t* correctiongmp = NULL;
   SCIP_Bool* isupper;
   int i;
   int j;
   int rval;
   int nextendedrows;
   int nrows;
   int nrowsps; /* number of rows used for ps. this can be lower than nrows, if e.g. cuts were added */
   int ncols;
   int currentrow;
   int shift;
   int startt, endt, setupt, violt, rectlut, projt, shiftt, cleanupt;
   SCIP_Bool isfeas;

   /* project-and-shift method:
    * 1. projection step (to ensure that equalities are satisfied):
    *   - compute error in equalities: r=c-Ay^
    *   - backsolve system of equations to find correction of error: z with Dz=r
    *   - add corretion to approximate dual solution: bold(y)=y^+[z 0]
    * 2. shifing step (to ensure that inequalities are satisfied):
    *   - take convex combination of projected approximate point bold(y) with interior point y*
    * 3. compute dual objective value of feasible dual solution and set bound
    */

   /* if data has not been constructed, or it failed, then exit */
   psdata = lpex->psdata;

   assert(psdata->psdatacon);
   if( (psdata->psdatafail && !usefarkas) || (usefarkas && !psdata->pshasray) )
   {
      lp->hasprovedbound = FALSE;
      return SCIP_OKAY;
   }

   /* constructpsdata should always be called first */
   if( usefarkas )
   {
      stat->nprojshift++;
      SCIPclockStart(stat->provedinfeaspstime, set);
   }
   else
   {
      stat->nprojshiftinf++;
      SCIPclockStart(stat->provedfeaspstime, set);
   }

   startt = clock();

   lp->hasprovedbound = TRUE;

   SCIPdebugMessage("calling getPSdual()\n");

   /* decide if we should use ray or point to compute bound */
   if( !usefarkas && psdata->psuseintpoint && psdata->pshaspoint )
   {
      /* if we are supposed to use the interior point, and it exists use it */
      useinteriorpoint = TRUE;
   }
   else
   {
      /* in this case, since psdatafail != TRUE, pshasray should be true -- use it */
      assert( psdata->pshasray );
      useinteriorpoint = FALSE;
   }

   SCIP_CALL( RatCreateBuffer(set->buffer, &tmp) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &tmp2) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &lambda1) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &lambda2) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &maxv) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &dualbound) );

   /* flush exact lp */
   /* set up the exact lpi for the current node */
   SCIP_CALL( SCIPsepastoreexSyncLPs(set->scip->sepastoreex, blkmem, set, stat, lpex, eventqueue, eventfilter) );
   SCIP_CALL( SCIPlpexFlush(lp->lpex, blkmem, set, eventqueue) );

   nextendedrows = psdata->nextendedrows;
   nrows = lpex->nrows;
   ncols = lpex->ncols;
   nrowsps = psdata->nextendedrows/2 - ncols;
   shift = nrows - nrowsps;

   assert(ncols == psdata->violationsize);

   /* allocate memory for approximate dual solution, dual cost vector, violation and correction */
   SCIPsetAllocBufferArray(set, &approxdual, nrows + ncols);
   violation = psdata->violation;
   correction = psdata->correction;

   /** @todo: exip this could be removed */
   for( i = 0; i < nrows + ncols; ++i )
   {
      if( i < nrows )
         approxdual[i] = usefarkas ? lpex->rows[i]->dualfarkas : lpex->rows[i]->dualsol;
      else
         approxdual[i] = usefarkas ? lpex->cols[i - nrows]->redcost : lpex->cols[i - nrows]->farkascoef;

      RatSetInt(approxdual[i], 0, 1);
      if( i < ncols )
         RatSetInt(violation[i], 0, 1);
   }
   for( i = 0; i < nextendedrows; ++i )
   {
      RatSetInt(correction[i], 0, 1);
   }

   SCIP_CALL( SCIPsetAllocBufferArray(set, &isupper, nrows + ncols) );

   setupt = clock();

   /* recover the objective coefs and approximate solution value of dual solution;
    * dual vars of lhs constraints (including -inf) and rhs constraints (including +inf),
    * dual vars of lb constraint (including -inf) and ub constraints (including +inf)
    */
   for( i = 0; i < nrows; i++ )
   {
      /* in case we want to prove infeasibility it might be that we were not able to compute a dual solution
       * with bound exceeding objective value; in this case dual solution is set to SCIP_INVALID
       */
      if( (!usefarkas && SCIProwGetDualsol(lp->rows[i]) == SCIP_INVALID) || (usefarkas && SCIProwGetDualfarkas(lp->rows[i]) == SCIP_INVALID) )
      {
         SCIPdebugMessage("  no valid unbounded approx dual sol given\n");
         lp->hasprovedbound = FALSE;
         if( usefarkas )
            stat->nfailprojshiftinf++;
         else
            stat->nfailprojshift++;

         goto TERMINATE;
      }

      if( !usefarkas )
         RatSetReal(approxdual[i], SCIProwGetDualsol(lp->rows[i]));
      else
         RatSetReal(approxdual[i], SCIProwGetDualfarkas(lp->rows[i]));

      if( RatIsPositive(approxdual[i]) )
         isupper[i] = FALSE;
      else
         isupper[i] = TRUE;
   }

   cols = SCIPlpGetCols(lp);
   for( i = 0; i < ncols; i++ )
   {
      if( !usefarkas )
         RatSetReal(approxdual[i+nrows], SCIPcolGetRedcost(cols[i], stat, lp));
      else
         RatSetReal(approxdual[i+nrows], -SCIPcolGetFarkasCoef(cols[i], stat, lp));

      if( RatIsPositive(approxdual[i+nrows]) )
         isupper[i+nrows] = FALSE;
      else
         isupper[i+nrows] = TRUE;
   }

   /* first, fix artificial dual variables to zero */
   for( i = 0; i < nrows + ncols; i++ )
   {
      SCIP_Rational* val;

       if( !isupper[i] )
         val = i < nrows ? lpex->rows[i]->lhs : lpex->cols[i - nrows]->lb;
      else
         val = i < nrows ? lpex->rows[i]->rhs : lpex->cols[i - nrows]->ub;

      if( RatIsAbsInfinity(val) )
         RatSetInt(approxdual[i], 0, 1);
   }

#ifdef PS_OUT
   printf("approximate dual solution:\n");

   RatSetInt(dualbound, 0, 1);
   for( i = 0; i < nrows + ncols; i++ )
   {
      SCIP_Rational* val;

      if( !isupper[i] )
         val = i < nrows ? lpex->rows[i]->lhs : lpex->cols[i - nrows]->lb;
      else
         val = i < nrows ? lpex->rows[i]->rhs : lpex->cols[i - nrows]->ub;

      printf("   i=%d: ", i);
      RatPrint(approxdual[i]);
      printf(" * ");
      RatPrint(val)
      printf("\n");
      if( RatIsAbsInfinity(val) )
         assert(RatIsZero(approxdual[i]));
      else
      {
         RatMult(tmp, approxdual[i], val);
         RatAdd(dualbound, dualbound, tmp);
      }
   }

   printf("   objective value=%.20f (", RgetRealApprox(dualbound));
   RatPrint(dualbound);
   printf(")\n");
#endif

   /* calculate violation of equality constraints r=c-A^ty */
   for( i = 0; i < ncols; i++ )
   {
      /* instead of setting and then subtracting the A^ty corresponding to bound constraints, we can do this directly */
      if( !usefarkas )
      {
         /* set to obj - bound-redcost */
         RatDiff(violation[i], lpex->cols[i]->obj, approxdual[i + nrows]);
      }
      else
      {
         /* set to 0 - bound-redcost */
         RatNegate(violation[i], approxdual[i+nrows]);
      }
   }

   /* A^ty for y corresponding to primal constraints */
   for( i = 0; i < nrows; i++ )
   {
      for( j = 0; j < lpex->rows[i]->len; j++)
      {
         currentrow = lpex->rows[i]->cols_index[j];
         RatMult(tmp, approxdual[i], lpex->rows[i]->vals[j]);
         RatDiff(violation[currentrow], violation[currentrow], tmp);
      }
   }

   /* project solution */

#ifdef PS_OUT
   printf("violation of solution:\n");
   for( i = 0; i < ncols; i++ )
   {
      printf("   i=%d: ", i);
      RatPrint(violation[i]);
      printf("\n");
   }
#endif

   /* if there is no violation of the constraints, then skip the projection */
   isfeas = TRUE;
   for( i = 0; i < ncols && isfeas; i++ )
   {
      if( !RatIsZero(violation[i]) )
         isfeas = FALSE;
   }

   violt = clock();

   /* isfeas is equal to one only if approximate dual solution is already feasible for the dual */
   if( !isfeas )
   {
      /* compute [z] with Dz=r (D depends on psdata->psdualcolselection) */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &violationgmp, ncols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &correctiongmp, nextendedrows) );
      RatSetGMPArray(violationgmp, violation, ncols);

      for ( i = 0; i < nextendedrows; i++)
      {
         mpq_init(correctiongmp[i]);
      }

      rval = RECTLUsolveSystem(psdata->rectfactor, ncols, nextendedrows, violationgmp, correctiongmp);
      if( rval )
      {
         lp->hasprovedbound = FALSE;
         if( usefarkas )
            stat->nfailprojshiftinf++;
         else
            stat->nfailprojshift++;

         goto TERMINATE;
      }
      RatSetArrayGMP(correction, correctiongmp, nextendedrows);

#ifdef PS_OUT
      printf("correction of solution:\n");
      for( i = 0; i < psdata->npsbasis; i++ )
      {
         printf("   i=%d: ", i);
         RatPrint(correction[i]);
         printf(", position=%d\n", psdata->psbasis[i]);
      }
#endif

      rectlut = clock();

      /* projection step: compute bold(y)=y^+[z 0];
       * save the corrected components in the correction vector
       */
      for( i = 0; i < psdata->npsbasis; i++ )
      {
         // map is the point in the extended space -> transform it back to the original space
         int map = psdata->psbasis[i];
         if( map < nrowsps )
         {
            if( !isupper[map] )
            {
               RatAdd(correction[i], correction[i], approxdual[map]);
               RatSetInt(approxdual[map], 0, 1);
            }
         }
         else if( map < 2 * nrowsps )
         {
            if( isupper[map - nrowsps] )
            {
               RatDiff(correction[i], correction[i], approxdual[map - nrowsps]);
               RatSetInt(approxdual[map - nrowsps], 0, 1);
            }
         }
         else if( map < 2 * nrowsps + ncols )
         {
            if( !isupper[map - nrowsps + shift] )
            {
               RatAdd(correction[i], correction[i], approxdual[map - nrowsps + shift]);
               RatSetInt(approxdual[map - nrowsps + shift], 0, 1);
            }
         }
         else
         {
            if( isupper[map - nrowsps - ncols  + shift] )
            {
               RatDiff(correction[i], correction[i], approxdual[map - nrowsps - ncols + shift]);
               RatSetInt(approxdual[map - nrowsps - ncols + shift], 0, 1);
            }
         }
      }
      projt = clock();

#ifdef PS_OUT
      printf("updated dual solution:\n");
      for( i = 0; i < psdata->npsbasis; i++ )
      {
         printf("   i=%d: ", i);
         RatPrint(correction[i]);
         printf(", position=%d\n", psdata->psbasis[i]);
      }
#endif

      if( useinteriorpoint )
      {
         assert(!usefarkas);
         /* shifting step (scale solution with interior point to be dual feasible):
         * y' = lambda1 bold(y) + lambda2 y*, where
         *   lambda1 = ( slack of int point)/ (slack of int point + max violation) = d/m+d
         *   lambda2 = 1 - lambda1
         */

         /* compute lambda1 componentwise (set lambda1 = 1 and lower it if necessary) */
         RatSetInt(lambda1, 1, 1);
         for( i = 0; i < psdata->npsbasis; i++ )
         {
            if( RatIsNegative(correction[i]) )
            {
               int map = psdata->psbasis[i];

               RatSet(tmp2, psdata->interiorpt[map]);
               RatDiff(tmp, psdata->interiorpt[map], correction[i]);
               RatDiv(tmp2, tmp2, tmp);
               RatMIN(lambda1, lambda1, tmp2);
            }
         }
         RatSetInt(lambda2, 1, 1);
         RatDiff(lambda2, lambda2, lambda1);
      }
      else
      {
         /* in this case we are using an interior ray that can be added freely to the solution */
         /* compute lambda values: compute lambda1 componentwise (set lambda1 = 1 and lower it if necessary) */
         RatSetInt(lambda1, 1, 1);
         //RatSetInt(lambda2, 0, 1);
         for( i = 0; i < psdata->npsbasis; i++ )
         {
            int map = psdata->psbasis[i];
            if( RatIsNegative(correction[i]) && psdata->includedrows[map] )
            {
               RatDiv(tmp, correction[i], psdata->interiorray[map]);
               RatNegate(tmp, tmp);
               RatMAX(lambda2, lambda2, tmp);
            }
         }
      }

#ifdef PS_OUT
      printf("transformed projected dual solution:\n");

      RatSetInt(dualbound, 0, 1);
      for( i = 0; i < nrows + ncols; i++ )
      {
         SCIP_Rational* val;

         printf("   i=%d: ", i);
         RatPrint(approxdual[i]);
         printf("\n");
      }

      printf("   lambda1: ");
      RatPrint(lambda1);
      printf(")\n");
#endif

      /* tranfsorm correction back to approxdual */
      for( i = 0; i < psdata->npsbasis; i++ )
      {
         int map = psdata->psbasis[i];
         if( map < nrowsps )
            RatAdd(approxdual[map], approxdual[map], correction[i]);
         else if( map < 2 * nrowsps )
            RatDiff(approxdual[map - nrowsps], approxdual[map - nrowsps], correction[i]);
         else if ( map < 2 * nrowsps + ncols )
            RatAdd(approxdual[map - nrowsps + shift], approxdual[map - nrowsps + shift], correction[i]);
         else
            RatDiff(approxdual[map - nrowsps - ncols + shift], approxdual[map - nrowsps - ncols + shift], correction[i]);
      }

#ifdef PS_OUT
      printf("transformed projected dual solution:\n");

      RatSetInt(dualbound, 0, 1);
      for( i = 0; i < nrows + ncols; i++ )
      {
         SCIP_Rational* val;

         printf("   i=%d: ", i);
         RatPrint(approxdual[i]);
         printf("\n");
      }

      printf("   lambda1: ");
      RatPrint(lambda1);
      printf(")\n");
#endif

      /* perform shift */
      if( !RatIsZero(lambda2) )
      {
         for( i = 0; i < nrows + ncols; i++ )
         {
            if( i < nrows && i >= nrowsps )
               continue;
            RatMult(approxdual[i], approxdual[i], lambda1);
         }
         for( i = 0; i < nrows + ncols; i++ )
         {
            /* todo @exip: refactor this somehow. explanation: when the number of lp-rows increased 
            * the number of rows in the ps-data does not. so we have [1,...,nrows, ... extrarows ..., 1, ... ncols]
            * so if we map to the column part in the extended space, we have to subtract the difference */
            int map;
            if( i < nrows && i >= nrowsps )
               continue;
            map = (i < nrowsps) ? i + nrowsps : i + nrowsps + ncols - shift;
            RatMult(tmp, useinteriorpoint ? psdata->interiorpt[map] : psdata->interiorray[map], lambda2);
            RatDiff(approxdual[i], approxdual[i], tmp);
            map = (i < nrowsps) ? i : i + nrowsps - shift;
            RatMult(tmp, useinteriorpoint ? psdata->interiorpt[map] : psdata->interiorray[map], lambda2);
            RatAdd(approxdual[i], approxdual[i], tmp);
         }
      }
      shiftt = clock();

#ifdef PS_OUT
      printf("projected and shifted dual solution (should be an exact dual feasible solution)\n");
      for( i = 0; i < nrows+ncols; i++ )
      {
         printf("   i=%d: ", i);
         RatPrint(approxdual[i]);
         printf("\n");
      }
#endif
   }


#ifndef NDEBUG
   SCIPdebugMessage("debug test: verifying feasibility of dual solution:\n");

   /* calculate violation of equality constraints: subtract Ax to get violation b-Ax, subtract A(approxdual) */
   rval = 0;
   for( i = 0; i < ncols; i++ )
   {
      if( !usefarkas )
        RatSet(violation[i], lpex->cols[i]->obj);
      else
         RatSetInt(violation[i], 0, 1);
   }
   for( i = 0; i < nrows; i++ )
   {
      for( j = 0; j < lpex->rows[i]->len; j++ )
      {
         currentrow = lpex->rows[i]->cols_index[j];
         RatMult(tmp, approxdual[i], lpex->rows[i]->vals[j]);
         RatDiff(violation[currentrow], violation[currentrow], tmp);
      }
   }
   for( i = 0; i < ncols; i++ )
   {
         RatDiff(violation[i], violation[i], approxdual[i + nrows]);
   }

   for( i = 0; i < ncols && rval == 0; i++ )
   {
      if( !RatIsZero(violation[i]) )
      {
         SCIPdebugMessage("   dual solution incorrect, violates equalties\n");
         rval = 1;
      }
   }
   if( !rval )
      SCIPdebugMessage("   dual solution verified\n");
   assert(!rval);
#endif

   RatSetInt(dualbound, 0, 1);
   for( i = 0; i < nrows + ncols; i++ )
   {
      SCIP_Rational* val;
      if( RatIsPositive(approxdual[i]) )
         val = i < nrows ? lpex->rows[i]->lhs : lpex->cols[i - nrows]->lb;
      else
         val = i < nrows ? lpex->rows[i]->rhs : lpex->cols[i - nrows]->ub;

      RatMult(tmp, approxdual[i], val);
      RatAdd(dualbound, dualbound, tmp);
   }

   if( !usefarkas )
   {
      RatSet(lpex->lpobjval, dualbound);
      lp->lpobjval = RatRoundReal(dualbound, SCIP_ROUND_DOWNWARDS);
      lp->hasprovedbound = TRUE;
   }
   else
   {
      /* if the objective value of the corrected ray is positive we can prune node, otherwise not */
      if( RatIsPositive(dualbound) )
      {
         RatSetString(lpex->lpobjval, "inf");
         lp->lpobjval = SCIPsetInfinity(set);
         lp->hasprovedbound = TRUE;
      }
      else
      {
         stat->nfailprojshiftinf++;
         lp->hasprovedbound = FALSE;
      }
   }

#ifdef PS_OUT
   printf("   common slack=%.20f (", RgetRealApprox(psdata->commonslack));
   RatPrint(psdata->commonslack);
   printf(")\n");

   printf("   max violation=%.20f (", RgetRealApprox(maxv));
   RatPrint(maxv);
   printf(")\n");

   printf("   lambda (use of interior point)=%.20f (", RgetRealApprox(lambda2));
   RatPrint(lambda2);
   printf(")\n");

   printf("   dual objective value=%.20f (", RgetRealApprox(dualbound));
   RatPrint(dualbound);
   printf(")\n");
#endif

 TERMINATE:
   /* free memory */
   if( correctiongmp != NULL )
   {
      RatClearGMPArray(correctiongmp, nextendedrows);
      SCIPsetFreeBufferArray(set, &correctiongmp);
   }
   if( violationgmp != NULL )
   {
      RatClearGMPArray(violationgmp, ncols);
      SCIPsetFreeBufferArray(set, &violationgmp);
   }

   SCIPsetFreeBufferArray(set, &isupper);
   SCIPsetFreeBufferArray(set, &approxdual);

   RatFreeBuffer(set->buffer, &dualbound);
   RatFreeBuffer(set->buffer, &maxv);
   RatFreeBuffer(set->buffer, &lambda2);
   RatFreeBuffer(set->buffer, &lambda1);
   RatFreeBuffer(set->buffer, &tmp2);
   RatFreeBuffer(set->buffer, &tmp);

   endt = clock();
   // startt, endt, setupt, violt, rectlut, projt, shiftt;
   //printf("Time used: %e \n", ((double) (endt - startt)) / CLOCKS_PER_SEC);
   //printf("Setup time    : %e, percentage %e \n", ((double) (setupt - startt)) / CLOCKS_PER_SEC, ((double) (setupt - startt)/(endt - startt)) );
   //printf("Violation time: %e, percentage %e \n", ((double) (violt - setupt)) / CLOCKS_PER_SEC, ((double) (violt - setupt)/(endt - startt))   );
   //printf("Rectlu time   : %e, percentage %e \n", ((double) (rectlut - violt)) / CLOCKS_PER_SEC, ((double) (rectlut - violt)/(endt - startt)) );
   //printf("Proj time     : %e, percentage %e \n", ((double) (projt - rectlut)) / CLOCKS_PER_SEC, ((double) (projt - rectlut)/(endt - startt)) );
   //printf("Shifting time : %e, percentage %e \n", ((double) (shiftt - projt)) / CLOCKS_PER_SEC, ((double) (shiftt - projt)/(endt - startt))   );
   //printf("Cleanup time  : %e, percentage %e \n", ((double) (endt - shiftt)) / CLOCKS_PER_SEC, ((double) (endt - shiftt)/(endt - startt))     );

   if( usefarkas )
      SCIPclockStop(stat->provedinfeaspstime, set);
   else
      SCIPclockStop(stat->provedfeaspstime, set);

   return SCIP_OKAY;
}

static
char chooseBoundingMethod(
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

   assert(!lpex->fplp->hasprovedbound);

   if( dualfarkas )
   {
      /* check if neumair-scher is possible */
      if( SCIPlpexBSpossible(lpex) && !SCIPisCertificateActive(set->scip) )
         dualboundmethod = 'n';
      /* check if project and shift is possible */
      else if( SCIPlpexPSpossible(lpex) )
      {
         SCIP_CALL( constructPSData(lp, lpex, set, stat, messagehdlr, eventqueue, eventfilter,
                        prob, blkmem) );
         if( !lpex->psdata->psdatafail )
            dualboundmethod = 'p';
         else
            dualboundmethod = 'e';
      }
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
      if( (lpex->interleavedbfreq > 0 && SCIPsetIsInfinity(set, SCIPlpGetCutoffbound(lpex->fplp)) && SCIPgetDepth(set->scip) > 0
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
         {
            SCIP_CALL( constructPSData(lp, lpex, set, stat, messagehdlr, eventqueue, eventfilter,
                           prob, blkmem) );
            if( !lpex->psdata->psdatafail )
               dualboundmethod = 'p';
            else
               dualboundmethod = 'e';
         }
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
         if( RatIsFpRepresentable(SCIPvarGetObjExact(SCIPcolGetVar(col))) )
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
   SCIPdebugMessage("safebound computed: %e, previous fp-bound: %e.17, difference %e.17 \n", *safebound, lp->lpobjval, *safebound - lp->lpobjval);

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
         lp->lpobjval = *safebound;
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

SCIP_RETCODE SCIPlpexComputeSafeBound(
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
   assert(!lp->hasprovedbound);

   /* if the lp was solved to optimality and there are no fractional variables we solve exactly to generate a feasible solution */
   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && fpLPisIntFeasible(lp, set, stat) )
      lpex->forceexactsolve = TRUE;

   /* choose which bounding method to use. only needed if automatic is enabled. */
   if( lpex->forceexactsolve )
   {
      dualboundmethod = 'e';
      lpex->forceexactsolve = FALSE;
   }
   else if( set->misc_dbmethod == 'a' )
   {
      dualboundmethod = chooseBoundingMethod(lp, lpex, set, messagehdlr, blkmem, stat, eventqueue, eventfilter,
                              prob, itlim, lperror, dualfarkas, safebound);
   }
   else
      dualboundmethod = set->misc_dbmethod;

   SCIPdebugMessage("Computing safe bound for lp with status: %d using bounding method %c \n", SCIPlpGetSolstat(lp), dualboundmethod);

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
      /* project and shift */
      case 'p':
         SCIP_CALL( constructPSData(lp, lpex, set, stat, messagehdlr, eventqueue, eventfilter,
                        prob, blkmem) );
         SCIP_CALL( getPSdual(lp, lpex, set, stat, messagehdlr, eventqueue, eventfilter,
                        prob, blkmem, dualfarkas) );
         break;
      /* exact LP */
      case 'e':
         SCIP_CALL( solveLpExact(lp, lpex, set, messagehdlr, blkmem, stat, eventqueue, eventfilter,
                        prob, itlim, lperror, dualfarkas) );
         break;
      default:
         SCIPerrorMessage("bounding method %c not implemented yet \n", set->misc_dbmethod);
         SCIPABORT();
         break;
   }

   if( !lp->hasprovedbound )
   {
      SCIP_CALL( solveLpExact(lp, lpex, set, messagehdlr, blkmem, stat, eventqueue, eventfilter,
                        prob, itlim, lperror, dualfarkas) );
   }

   /* choose which bounding method should be calles and return a safe objective bound */
   return SCIP_OKAY;
}

#endif