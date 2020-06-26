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
#include "scip/lpexact.h"
#include "lpi/lpiexact.h"
#include "scip/scip_prob.h"
#include "scip/prob.h"
#include "scip/scip.h"
#include "scip/sol.h"
#include "scip/primal.h"
#include "scip/sepastoreexact.h"
#include "scip/struct_scip.h"
#include "rectlu/rectlu.h"

#define PSBIGM                      100
#define PSWARMSTARTAUXPROB         TRUE
#define PSPOSTPROCESSDUALSOL       TRUE


#ifdef SCIP_WITH_BOOST
static
SCIP_Bool fpLPisIntFeasible(
   SCIP_LP*              lp,
   SCIP_SET*             set,
   SCIP_STAT*            stat
   )
{
   SCIP_Real primsol;
   SCIP_Bool feasible;
   SCIP_Real frac;
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

      // we can't use SCIPsetIsIntegral since we have to chech here in the same way as in SCIPcalcBranchCands
      primsol = SCIPcolGetPrimsol(col);
      frac = SCIPsetFeasFrac(set, primsol);
      if( !SCIPsetIsFeasFracIntegral(set, frac) )
      {
         feasible = FALSE;
         break;
      }
   }

   return feasible;
}

/** evaluates the result of the exact LP */
static
SCIP_RETCODE evaluateLPExact(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEXACT*         lpexact,            /**< Exact LP data */
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
   /** @todo exip:
    *  - in a similar function for LP in lp.c the case SCIPlpiIsPrimalUnbounded() is not explicitely handled. why? if case
    *    should be added, include it here as well.
    */
   /* evaluate solution status */
   if( SCIPlpiExactIsOptimal(lpexact->lpiexact) )
   {
      SCIPdebugMessage("   exact LP solved to optimality\n");

#ifndef NDEBUG
      {
         SCIP_Bool primalfeasible;
         SCIP_Bool dualfeasible;

         SCIP_CALL( SCIPlpiExactGetSolFeasibility(lpexact->lpiexact, &primalfeasible, &dualfeasible) );
         assert(primalfeasible);
         assert(dualfeasible);
      }
#endif
      /* check whether exact LP solution is feasible for the MIP; if it is feasible the solution is stored
       * and the current node is cut off otherwise a branching is created
       */
      assert(*result == SCIP_CUTOFF || *result == SCIP_BRANCHED || *result == SCIP_SOLVELP );
   }
   else if( SCIPlpiExactIsObjlimExc(lpexact->lpiexact) )
   {
      /* CERT: for now we don't trust the exact LP solver on this, so we do not need to print anything here */
      SCIPdebugMessage("   exact LP exceeds upper objective limit\n");
      *result = SCIP_CUTOFF;
   }
   else if( SCIPlpiExactIsPrimalInfeasible(lpexact->lpiexact) )
   {
      SCIPdebugMessage("   exact LP is primal infeasible\n");
      *result = SCIP_CUTOFF;
   }
   else if( SCIPlpiExactExistsPrimalRay(lpexact->lpiexact) )
   {
      /* CERT: for now we don't need to handle the unbounded case */
      SCIPerrorMessage("exact LP has primal ray: case not handled yet\n");
      return SCIP_ERROR;
   }
   else if( SCIPlpiExactIsIterlimExc(lpexact->lpiexact) )
   {
      SCIPerrorMessage("exact LP exceeds iteration limit: case not handled yet\n");
      return SCIP_ERROR;
   }
   else if( SCIPlpiExactIsTimelimExc(lpexact->lpiexact) )
   {
      /* CERT: for now we should not activate the buggy time limit support of QSopt_ex */
      SCIPdebugMessage("   exact LP exceeds time limit\n");
      *result = SCIP_INFEASIBLE;
   }
   else
   {
      SCIPerrorMessage(" error or unknown return status in current exact LP (internal status: %d)\n",
         SCIPlpiExactGetInternalStatus(lpexact->lpiexact));
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** helper method, compute number of nonzeros in lp */
static
int getNNonz(
   SCIP_LPEXACT* lpexact
   )
{
   int ret = 0;
   int i = 0;

   assert(lpexact != NULL);

   for( i = 0; i < lpexact->nrows; ++i )
      ret+= lpexact->rows[i]->len;

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

#ifdef SCIP_WITH_GMP
/** subroutine of constructPSdata(); chooses which columns of the matrix are designated as set S, used for projections */
static
SCIP_RETCODE psChooseS(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
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
   SCIP_PROJSHIFTDATA* projshiftdata;
   SCIP_Bool dualfeasible;
   SCIP_Bool primalfeasible;

   nrows = lpexact->nrows;
   ncols = lpexact->ncols;
   projshiftdata = lpexact->projshiftdata;

   assert(projshiftdata != NULL);

   nextendedrows = projshiftdata->nextendedrows;

   /* build includedrows vector based on psfpdualcolwiseselection, this determines the matrix D */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &projshiftdata->includedrows, nextendedrows) );
   for( i = 0; i < nextendedrows; i++ )
      projshiftdata->includedrows[i] = 0;
   if( projshiftdata->psdualcolselection == PS_DUALCOSTSEL_NO || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE )
   {
      /* determine which dual variables to included in the problem
       * (ones with finite dual objective coef. in [lhs',-rhs',lb',-ub'])
       */
      for( i = 0; i < nrows; i++ )
      {
         if( !RatIsNegInfinity(lpexact->rows[i]->lhs) )
            projshiftdata->includedrows[i] = 1;
         if( !RatIsInfinity(lpexact->rows[i]->rhs) )
            projshiftdata->includedrows[nrows + i] = 1;
      }
      for( i = 0; i < ncols; i++ )
      {
         if( !RatIsNegInfinity(lpexact->cols[i]->lb) )
            projshiftdata->includedrows[2*nrows + i] = 1;
         if( !RatIsInfinity(lpexact->cols[i]->ub) )
            projshiftdata->includedrows[2*nrows + ncols + i] = 1;
      }
   }
   else if( projshiftdata->psdualcolselection == PS_DUALCOSTSEL_ACTIVE_EXLP )
   {
      /* determone which dual variables to include in the problem (in this case we choose dual variables whose primal
       * constraints are active at the solution of the exact LP at the root node)
       */

      SCIP_CALL( SCIPlpExactSolveAndEval(lpexact, lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, prob, 100,
               &lperror, FALSE) );

      SCIP_CALL( RatCreateBufferArray(set->buffer, &rootprimal, ncols) );
      SCIP_CALL( RatCreateBufferArray(set->buffer, &rootactivity, nrows) );

      /* get the primal solution and activity */
      SCIP_CALL( SCIPlpiExactGetSol(lpexact->lpiexact, NULL, rootprimal, NULL, rootactivity, NULL) );

      /* determine which dual variables to include in the problem
       * (primal constraints active at optimal solution found at root node)
       */
      for( i = 0; i < nrows; i++ )
      {
         if( RatIsEqual(rootactivity[i], lpexact->rows[i]->lhs) )
            projshiftdata->includedrows[i] = 1;
         if( RatIsEqual(rootactivity[i], lpexact->rows[i]->rhs) )
            projshiftdata->includedrows[nrows + i] = 1;
      }
      for( i = 0; i < ncols; i++ )
      {
         if( RatIsEqual(rootprimal[i], lpexact->cols[i]->lb) )
            projshiftdata->includedrows[2*nrows + i] = 1;
         if( RatIsEqual(rootprimal[i], lpexact->cols[i]->ub) )
            projshiftdata->includedrows[2*nrows + ncols + i] = 1;
      }

      RatFreeBufferArray(set->buffer, &rootactivity, nrows);
      RatFreeBufferArray(set->buffer, &rootprimal, ncols);
   }
   else if( projshiftdata->psdualcolselection == PS_DUALCOSTSEL_ACTIVE_FPLP )
   {
      /* determine which dual variables to include in the problem (in this case we choose dual variables whose primal
       * constraints are active at the solution of the exact LP at the root node)
       */

      assert(lp->nrows == nrows);
      for( i = 0; i < nrows; i++ )
      {
         if( SCIPsetIsFeasEQ(set, SCIProwGetLPActivity(lp->rows[i], set, stat, lp), SCIProwGetLhs(lp->rows[i])) )
            projshiftdata->includedrows[i] = 1;
         if( SCIPsetIsFeasEQ(set, SCIProwGetLPActivity(lp->rows[i], set, stat, lp), SCIProwGetRhs(lp->rows[i])) )
            projshiftdata->includedrows[nrows + i] = 1;
      }

      assert(ncols == lp->ncols);
      for( i = 0; i < ncols; i++ )
      {
         if( SCIPsetIsFeasEQ(set, SCIPcolGetPrimsol(lp->cols[i]), SCIPcolGetLb(lp->cols[i])) )
            projshiftdata->includedrows[2*nrows + i] = 1;
         if( SCIPsetIsFeasEQ(set, SCIPcolGetPrimsol(lp->cols[i]), SCIPcolGetUb(lp->cols[i])) )
            projshiftdata->includedrows[2*nrows + ncols + i] = 1;
      }
   }
   else
   {
      SCIPerrorMessage("Invald value for parameter psfpdualcolwiseselection \n");
   }
   return SCIP_OKAY;
}

/** subroutine of constructPSdata(); computes the LU factorization used by the project-and-shift method */
static
SCIP_RETCODE psFactorizeD(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
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
   SCIP_PROJSHIFTDATA* projshiftdata;

   projshiftdata = lpexact->projshiftdata;

   nrows = lpexact->nrows;
   ncols = lpexact->ncols;
   nextendedrows = projshiftdata->nextendedrows;
   nnonz = getNNonz(lpexact);

   /* allocate memory for the projection factorization */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &projbeg, nextendedrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &projlen, nextendedrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &projind, 2*nnonz + 2*ncols) );
   for( i = 0; i < 2*nnonz + 2*ncols; i++)
      projind[i]=0;
   SCIP_CALL( RatCreateBufferArray(set->buffer, &projval, 2*nnonz + 2*ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &projvalgmp, 2*nnonz + 2*ncols) );

   /* allocate memory for the basis mapping */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &projshiftdata->projshiftbasis, nextendedrows) );

   /* use includedrows to construct projshiftbasis, a description/mapping for D it has length projshiftbasisdim and projshiftbasis[i] tells what
    * column (out of the original nextendecons) the ith column in D is
    */
   pos = 0;
   for( i = 0; i < nextendedrows; i++ )
   {
      if( projshiftdata->includedrows[i] )
      {
         projshiftdata->projshiftbasis[pos] = i;
         pos++;
      }
   }
   projshiftdata->projshiftbasisdim = pos;

   /* build the sparse representation of D that will be passed to the RECTLU code for factorization */
   pos = 0;
   for( i = 0; i < nextendedrows; i++ )
   {
      /* A part (lhs constraints) */
      if(i < nrows)
      {
         projlen[i] = lpexact->rows[i]->len;
         projbeg[i] = pos;
         for(j = 0; j < projlen[i]; j++)
         {
            projind[ projbeg[i] + j ] = lpexact->rows[i]->cols_index[j];
            RatSet( projval[ projbeg[i] + j], lpexact->rows[i]->vals[j]);
         }
         pos += projlen[i];
      }
      /* -A part (rhs constraints) */
      else if(i < 2 * nrows)
      {
         projlen[i] = lpexact->rows[i - nrows]->len;
         projbeg[i] = pos;
         for(j = 0; j < projlen[i]; j++)
         {
            projind[ projbeg[i] + j ] = lpexact->rows[i - nrows]->cols_index[j];
            RatNegate( projval[ projbeg[i] + j], lpexact->rows[i - nrows]->vals[j]);
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
   printf("factoring matrix: ncols=%d, projshiftbasisdim=%d\n", ncols, projshiftdata->projshiftbasisdim);
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
    * - projshiftbasis stores a mapping to tell us what D is, i.e. the dual columns corresponding to dual valuse that have a
    *   strictly positive value in the relative interior point
    * - D is equal to a subset of [A',-A',I,-I] and is given to the factor code in sparse column representation
    */
   RatSetGMPArray(projvalgmp, projval, 2 * nnonz + 2 * ncols);

   rval = RECTLUbuildFactorization(&projshiftdata->rectfactor, ncols, projshiftdata->projshiftbasisdim,
      projshiftdata->projshiftbasis, projvalgmp, projind, projbeg, projlen);

   /* if rval != 0 then RECTLUbuildFactorization has failed. In this case the project-and-shift method will not work and
    * we will return failure
    */
   if( rval )
   {
      projshiftdata->projshiftdatafail = TRUE;
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

/** prints error related to the current lpiexact status, if there is one */
static
SCIP_RETCODE printlpiexacterr(
   SCIP_LPIEXACT*        lpiexact            /**< lpiexact interface */
   )
{
   if( SCIPlpiExactIsOptimal(lpiexact) )
   {
      return SCIP_OKAY;
   }
   else if( SCIPlpiExactIsObjlimExc(lpiexact) )
   {
      SCIPerrorMessage("exact LP exceeds objlimit: case not handled yet\n");
   }
   else if( SCIPlpiExactIsPrimalInfeasible(lpiexact) )
   {
      SCIPerrorMessage(" Exact LP infeas.\n");
   }
   else if( SCIPlpiExactExistsPrimalRay(lpiexact) )
   {
      SCIPerrorMessage("exact LP has primal ray: case not handled yet\n");
   }
   else if( SCIPlpiExactIsIterlimExc(lpiexact) )
   {
      SCIPerrorMessage("exact LP exceeds iteration limit: case not handled yet\n");
   }
   else if( SCIPlpiExactIsTimelimExc(lpiexact) )
   {
      SCIPerrorMessage("exact LP exceeds time limit: case not handled yet\n");
   }
   else
   {
      SCIPerrorMessage("lpiexact not solved, or other error\n");
   }
   return SCIP_OKAY;
}

/** setup the data for ps in optimal version */
static
SCIP_RETCODE setupPsOpt(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
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
   SCIP_PROJSHIFTDATA* projshiftdata;
   int pos;
   int i;
   int j;
   int indx;

   projshiftdata = lpexact->projshiftdata;

   /* set up the objective */
   pos = 0;
   for( i = 0; i < nrows; i++ )
   {
      if( dvarincidence[i] )
      {
         RatSet(psobj[pos], lpexact->rows[i]->lhs);
         pos++;
      }
   }
   for( i = 0; i < nrows; i++ )
   {
      if( dvarincidence[nrows + i] )
      {
         RatNegate(psobj[pos], lpexact->rows[i]->rhs);
         pos++;
      }
   }
   for( i = 0; i < ncols; i++ )
   {
      if( dvarincidence[2*nrows + i] )
      {
         RatSet(psobj[pos], lpexact->cols[i]->lb);
         pos++;
      }
   }
   for( i = 0; i < ncols; i++ )
   {
      if( dvarincidence[2*nrows + ncols + i])
      {
         RatNegate(psobj[pos], lpexact->cols[i]->ub);
         pos++;
      }
   }
   assert(pos == ndvarmap);

   /* set alpha and beta. */
   RatSetReal(alpha, projshiftdata->projshiftobjweight);
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
      RatSet(pslhs[i], lpexact->cols[i]->obj);
      RatSet(psrhs[i], lpexact->cols[i]->obj);
   }
   for( i = 0; i < projshiftdata->projshiftbasisdim; i++ )
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
         for( j = 0; j < lpexact->rows[indx]->len; j++ )
         {
            pslen[lpexact->rows[indx]->cols_index[j]]++;
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
   for( i = 0; i < projshiftdata->projshiftbasisdim; i++ )
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
         for(j = 0; j < lpexact->rows[indx]->len; j++)
         {
            pos = psbeg[lpexact->rows[indx]->cols_index[j]] + pslen[lpexact->rows[indx]->cols_index[j]];
            psind[pos] = i;
            if(dvarmap[i]<nrows)
               RatSet(psval[pos], lpexact->rows[indx]->vals[j]);
            else
               RatNegate(psval[pos], lpexact->rows[indx]->vals[j]);
            pslen[lpexact->rows[indx]->cols_index[j]]++;
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
   /* set up the last projshiftbasisdim rows */
   pos = ncols;
   for( i = 0; i < ndvarmap; i++ )
   {
      indx = dvarmap[i];
      if( projshiftdata->includedrows[indx] )
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
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
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
   SCIP_PROJSHIFTDATA* projshiftdata;
   int pos;
   int i;
   int j;
   int indx;
   int nextendedrows;

   projshiftdata = lpexact->projshiftdata;
   nextendedrows = projshiftdata->nextendedrows;

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
         pslen[i] = lpexact->rows[indx]->len + 1;
         psbeg[i] = pos;

         /* set A,-A part of row */
         for(j = 0; j < pslen[i] - 1; j++)
         {
            psind[ psbeg[i] + j ] = lpexact->rows[indx]->cols_index[j];
            /* if it is an LHS constraint */
            if( dvarmap[i] < nrows )
               RatSet( psval[ psbeg[i] + j], lpexact->rows[indx]->vals[j]);
            /* otherwise it is the RHS of the constraint */
            else
               RatNegate( psval[ psbeg[i] + j], lpexact->rows[indx]->vals[j]);
         }
         /* set I part of row */
         psind[ psbeg[i] + pslen[i] - 1] = ncols + i ;
         RatSetInt( psval[psbeg[i] + pslen[i] - 1 ], -1, 1);

         /* update pos */
         pos += lpexact->rows[indx]->len + 1;
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
      if( !RatIsZero(lpexact->cols[i]->obj) )
      {
         RatNegate(psval[pos], lpexact->cols[i]->obj);
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
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
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
   SCIP_PROJSHIFTDATA* projshiftdata;
   int pos;
   int i;
   int j;
   int indx;

   projshiftdata = lpexact->projshiftdata;

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
   for( i = 0; i < projshiftdata->projshiftbasisdim; i++ )
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
         for( j = 0 ; j <lpexact->rows[indx]->len; j++)
         {
            pslen[lpexact->rows[indx]->cols_index[j]]++;
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
   for( i = 0; i < projshiftdata->projshiftbasisdim; i++ )
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
         for(j = 0; j < lpexact->rows[indx]->len; j++)
         {
            pos = psbeg[lpexact->rows[indx]->cols_index[j]] + pslen[lpexact->rows[indx]->cols_index[j]];
            psind[pos] = i;
            if(dvarmap[i]<nrows)
               RatSet(psval[pos], lpexact->rows[indx]->vals[j]);
            else
               RatNegate(psval[pos], lpexact->rows[indx]->vals[j]);
            pslen[lpexact->rows[indx]->cols_index[j]]++;
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
      RatNegate(psval[psbeg[i] + pslen[i]], lpexact->cols[i]->obj);
      psind[psbeg[i] + pslen[i]] = ndvarmap;
      pslen[i]++;
   }

   /* set up the last projshiftbasisdim rows */
   pos = ncols;
   for( i = 0; i < ndvarmap; i++ )
   {
      indx = dvarmap[i];
      if( projshiftdata->includedrows[indx] )
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
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
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
   SCIP_PROJSHIFTDATA* projshiftdata;
   int pos;
   int i;
   int j;
   int indx;

   projshiftdata = lpexact->projshiftdata;

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
      RatSet(pslhs[i], lpexact->cols[i]->obj);
      RatSet(psrhs[i], lpexact->cols[i]->obj);
   }
   for( i = 0; i < projshiftdata->projshiftbasisdim; i++ )
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
         for(j = 0; j < lpexact->rows[indx]->len; j++)
         {
            pslen[lpexact->rows[indx]->cols_index[j]]++;
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
   for( i = 0; i < projshiftdata->projshiftbasisdim; i++ )
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
         for(j = 0; j < lpexact->rows[indx]->len; j++)
         {
            pos = psbeg[lpexact->rows[indx]->cols_index[j]] + pslen[lpexact->rows[indx]->cols_index[j]];
            psind[pos] = i;
            if(dvarmap[i]<nrows)
               RatSet(psval[pos], lpexact->rows[indx]->vals[j]);
            else
               RatNegate(psval[pos], lpexact->rows[indx]->vals[j]);
            pslen[lpexact->rows[indx]->cols_index[j]]++;
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

   /* set up the last projshiftbasisdim rows */
   pos = ncols;
   for( i = 0; i < ndvarmap; i++ )
   {
      indx = dvarmap[i];
      if( projshiftdata->includedrows[indx] )
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
   SCIP_LPEXACT*         lpexact,            /**< the exact lp */
   int*                  dvarincidence       /**< the columns with existing bounds */
   )
{
   int ret = 0;
   int i;
   int nrows = lpexact->nrows;
   int ncols = lpexact->ncols;

   for( i = 0; i < nrows; i++ )
   {
      if( dvarincidence[i] )
         ret += lpexact->rows[i]->len;
      if( dvarincidence[nrows + i] )
         ret += lpexact->rows[i]->len;
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
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
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

   /* lpiexact and data used for the aux. problem */
   SCIP_LPIEXACT* pslpiexact;
   SCIP_LPISTATE* lpistate;
   SCIP_Rational** psobj = NULL;
   SCIP_Rational** pslb = NULL;
   SCIP_Rational** psub = NULL;
   SCIP_PROJSHIFTDATA* projshiftdata;
   SCIP_ROWEXACT** lprows;
   SCIP_COLEXACT** lpcols;

   SCIP_Rational** pslhs = NULL;
   SCIP_Rational** psrhs = NULL;

   int* psbeg;
   int* pslen;
   int* psind;
   SCIP_Rational** psval = NULL;
   SCIP_Rational** sol = NULL; /* either primal or dualsol */
   char ** colnames = NULL;
   SCIP_Rational* objval;

   /* mapping between variables used in the aux. problem and the original problem */
   int ndvarmap;
   int* dvarmap;
   int* dvarincidence;

   projshiftdata = lpexact->projshiftdata;
   nrows = lpexact->nrows;
   ncols = lpexact->ncols;

   assert(projshiftdata != NULL);
   nextendedrows = projshiftdata->nextendedrows;

   psncols = 0;
   psnrows = 0;
   psnnonz = 0;
   lprows = lpexact->rows;
   lpcols = lpexact->cols;

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
   if( projshiftdata->psintpointselection == PS_INTPOINTSEL_OPT || !findintpoint )
   {
      /* in this case we will find an optimized interior point for which we will try to push it interior and
       * optimize over its objective value.  To do this we will solve the following problem
       * max \alpha * [lhs,-rhs,lb,ub] * y + \beta d
       *              s.t. [A,-A,I,-I] * y        = c
       *                                 y_i - d >= 0 for each i \in S
       *                                     y   >= 0
       *                                  M >= d >= 0
       * M is a bound on how interior we will let the point be, S is the set of dual columns chosen earlier
       * which could have nonzero values for the S-interior point.
       *
       * After solving this y will be the S-interior point and d will be the common slack.
       * Here we actually construct the dual in row representation so it can be solved directly.
       */

      psncols =  ndvarmap + 1;
      psnrows = ncols + projshiftdata->projshiftbasisdim;
      psnnonz = computePsNnonz(lpexact, dvarincidence);
      psnnonz += 2*projshiftdata->projshiftbasisdim;

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
       * beta is set equal to the param projshiftobjweight and alpha is set equal to
       * alpha := (1-beta)/||OBJ||
       */

      SCIP_CALL( setupPsOpt(lp, lpexact, set, prob, psobj, psub, pslb, pslhs, psrhs, psval,
         pslen, psind, psbeg, dvarincidence, dvarmap, alpha, beta, tmp, psnrows, psnnonz,
         psncols, ndvarmap, nrows, ncols, findintpoint) );

      /* build aux LP using the exact LP interface */
      SCIP_CALL( SCIPlpiExactCreate(&pslpiexact, NULL, "pslpiexact", SCIP_OBJSEN_MAXIMIZE) );

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiExactAddCols(pslpiexact, psncols, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiExactAddRows(pslpiexact, psnrows, pslhs, psrhs, NULL, psnnonz, psbeg, psind, psval) );

      /* solve the LP */
      SCIP_CALL( SCIPlpiExactSolveDual(pslpiexact) );

      /* recover the optimal solution and set interior point and slack in constraint handler data */
      if( SCIPlpiExactIsOptimal(pslpiexact) )
      {
         SCIPdebugMessage("   exact LP solved to optimality\n");
         /* get optimal dual solution */
         SCIP_CALL( SCIPlpiExactGetSol(pslpiexact, objval, sol, NULL, NULL, NULL) );

         RatSet(projshiftdata->commonslack, sol[psncols - 1]);
         if( RatIsZero(projshiftdata->commonslack) )
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
                  RatSet( projshiftdata->interiorpoint[dvarmap[i]], sol[i]);
               else
                  RatSet( projshiftdata->interiorray[dvarmap[i]], sol[i]);
            }
            if( findintpoint )
               projshiftdata->projshifthaspoint = TRUE;
            else
               projshiftdata->projshifthasray = TRUE;
         }
      }
      else
      {
         SCIP_CALL( printlpiexacterr( pslpiexact ) );
      }
   }
   /* use 'a'rbitrary interior point */
   else if( projshiftdata->psintpointselection == PS_INTPOINTSEL_ARB )
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
         if( !RatIsZero(lpexact->cols[i]->obj) )
         {
            nobjnz++;
         }
      }

      psnnonz = computePsNnonz(lpexact, dvarincidence);
      psnnonz += nobjnz + 1 + 3 * ndvarmap;

      /* allocate memory for aux problem */
      SCIP_CALL( allocIntMem(set, &psobj, &pslb, &psub, &pslhs, &psrhs, &psval, &sol,
         &objval, &psbeg, &pslen, &psind, &colnames, psncols, psnrows, psnnonz) );

      SCIP_CALL( setupPsArb(lp, lpexact, set, prob, psobj, psub, pslb, pslhs, psrhs, psval,
         pslen, psind, psbeg, dvarincidence, dvarmap, alpha, beta, tmp, psnrows, psnnonz,
         psncols, ndvarmap, nrows, ncols, nobjnz, findintpoint) );

      SCIPdebugMessage("building LPIEXACT for aux. problem\n");

      /* build aux LP using the exact LP interface */
      SCIP_CALL( SCIPlpiExactCreate(&pslpiexact, NULL, "lpiexactps", SCIP_OBJSEN_MINIMIZE) );

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiExactAddCols(pslpiexact, psncols, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiExactAddRows(pslpiexact, psnrows, pslhs, psrhs,
            NULL, psnnonz, psbeg, psind, psval) );

      /* solve the LP */
      SCIPdebugMessage("solving aux. problem\n");
      SCIP_CALL( SCIPlpiExactSolveDual(pslpiexact) );

      if( SCIPlpiExactIsOptimal(pslpiexact) )
      {
         SCIPdebugMessage("   exact LP solved to optimality\n");

         /* compute 1/lambda (lambda is the dual variable corresponding to the last row in the aux LP) */
         SCIP_CALL( SCIPlpiExactGetSol(pslpiexact, objval, NULL, sol, NULL, NULL) );
         if( !RatIsZero(sol[psnrows - 1]) )
            RatInvert(projshiftdata->commonslack, sol[psnrows - 1]);
         else
         {
            RatSetInt(projshiftdata->commonslack, 0, 1);
         }
         if( RatIsZero(projshiftdata->commonslack) )
         {
            SCIPdebugMessage("   --> project-and-shift did not find S-interior point/ray\n");
         }

         /* interior point is y/lambda */
         projshiftdata->projshifthaspoint = TRUE;
         for( i = 0; i < ndvarmap; i++ )
         {
            if( projshiftdata->includedrows[dvarmap[i]] && RatIsZero(sol[i]) )
            {
               SCIPdebugMessage("   --> project-and-shift did not find S-interior point/ray\n");
               projshiftdata->projshifthaspoint = FALSE;
               i = ndvarmap;
            }
            else
               RatDiv( projshiftdata->interiorpoint[dvarmap[i]], sol[i], sol[psnrows - 1]);
         }
      }
      else
      {
         SCIP_CALL( printlpiexacterr( pslpiexact ) );
      }
   }
   /* use 'A'rbitrary interior point in transposed form*/
   else if( projshiftdata->psintpointselection == PS_INTPOINTSEL_ARBDUAL )
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
      psncols =  ndvarmap + 1 + projshiftdata->projshiftbasisdim;
      psnrows = ncols + projshiftdata->projshiftbasisdim;
      psnnonz = computePsNnonz(lpexact, dvarincidence);
      psnnonz += 2*projshiftdata->projshiftbasisdim + ncols;

      SCIP_CALL( allocIntMem(set, &psobj, &pslb, &psub, &pslhs, &psrhs, &psval, &sol,
         &objval, &psbeg, &pslen, &psind, &colnames, psncols, psnrows, psnnonz) );

      SCIP_CALL( setupPsArbDual(lp, lpexact, set, prob, psobj, psub, pslb, pslhs, psrhs, psval,
         pslen, psind, psbeg, dvarincidence, dvarmap, alpha, beta, tmp, psnrows, psnnonz,
         psncols, ndvarmap, nrows, ncols, findintpoint) );

      SCIPdebugMessage("building LPIEXACT for aux. problem\n");

      /* build aux LP using the exact LP interface */
      SCIP_CALL( SCIPlpiExactCreate(&pslpiexact, NULL, "pslpiexact", SCIP_OBJSEN_MINIMIZE) );

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiExactAddCols(pslpiexact, psncols, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiExactAddRows(pslpiexact, psnrows, pslhs, psrhs,
            NULL, psnnonz, psbeg, psind, psval) );

      /* solve the LP */
      SCIPdebugMessage("solving aux. problem\n");
      SCIP_CALL( SCIPlpiExactSolveDual(pslpiexact) );

      if( SCIPlpiExactIsOptimal(pslpiexact) )
      {
         SCIPdebugMessage("   exact LP solved to optimality\n");

         /* compute 1/lambda (lambda is the dual variable corresponding to the last row in the aux LP) */
         SCIP_CALL( SCIPlpiExactGetSol(pslpiexact, objval, sol, NULL, NULL, NULL) );
         if( !RatIsZero(sol[ndvarmap]) )
            RatInvert(projshiftdata->commonslack, sol[ndvarmap]);
         else
         {
            RatSetInt(projshiftdata->commonslack, 0, 1);
         }
         if( RatIsZero(projshiftdata->commonslack) == 0 )
         {
            SCIPdebugMessage("   --> interior point not found\n");
         }

         /* interior point is y/lambda */
         projshiftdata->projshifthaspoint = TRUE;
         for( i = 0; i < ndvarmap; i++ )
         {
            if( projshiftdata->includedrows[dvarmap[i]] && RatIsZero(sol[i]) )
            {
               projshiftdata->projshifthaspoint = FALSE;
               SCIPdebugMessage("   --> interior point not found\n");
               i = ndvarmap;
            }
            else
               RatDiv( projshiftdata->interiorpoint[dvarmap[i]], sol[i], sol[ndvarmap]);
         }

#ifdef PS_OUT
         printf("constraints all satisfied by slack=");
         RatPrint(projshiftdata->commonslack);
         printf("\n");

         printf("objective value of aux problem=");
         RatPrint(objval);
         printf("\n");

         printf("relative interior solution:\n");
         for( i = 0; i <  projshiftdata->nextendedrows; i++ )
         {
            printf("   i=%d: ", i);
            RatPrint(projshiftdata->interiorpoint[i]);
            printf("\n");
         }
#endif
      }
      else
      {
         printlpiexacterr( pslpiexact );
         projshiftdata->projshiftdatafail = TRUE;
      }
   }
   else if( projshiftdata->psintpointselection == PS_INTPOINTSEL_TWOSTAGE )
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
       * which could have nonzero values for the S-interior point.
       *
       * After solving this y will be the S-interior point and d will be the common slack.
       * Here we actually construct the dual in row representation so it can be solved directly.
       */
      psncols =  ndvarmap + 1;
      psnrows = ncols + projshiftdata->projshiftbasisdim;
      psnnonz = computePsNnonz(lpexact, dvarincidence);
      psnnonz += 2*projshiftdata->projshiftbasisdim;

      SCIP_CALL( allocIntMem(set, &psobj, &pslb, &psub, &pslhs, &psrhs, &psval, &sol,
         &objval, &psbeg, &pslen, &psind, &colnames, psncols, psnrows, psnnonz) );

      SCIP_CALL( setupPsTwoStage(lp, lpexact, set, prob, psobj, psub, pslb, pslhs, psrhs, psval,
         pslen, psind, psbeg, dvarincidence, dvarmap, alpha, beta, tmp, psnrows, psnnonz,
         psncols, ndvarmap, nrows, ncols, findintpoint) );

      /* build aux LP using the exact LP interface */
      if( pslpiexact != NULL )
      {
         SCIP_CALL( SCIPlpiExactFree(&pslpiexact) );
      }
      pslpiexact = NULL;
      SCIP_CALL( SCIPlpiExactCreate(&pslpiexact, NULL, "pslpiexact", SCIP_OBJSEN_MAXIMIZE) );

#ifdef PS_OUT
      /* activate extra output from exact LP solver */
      SCIPlpiExactSetIntpar(pslpiexact, SCIP_LPPAR_LPINFO, TRUE);
#endif

      /* add all columns to the exact LP */
      SCIP_CALL( SCIPlpiExactAddCols(pslpiexact, psncols, psobj, pslb, psub, colnames, 0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiExactAddRows(pslpiexact, psnrows, pslhs, psrhs,
            NULL, psnnonz, psbeg, psind, psval) );

      /* solve the LP */
      SCIP_CALL( SCIPlpiExactSolveDual(pslpiexact) );

      /* get state and solution of lpiexact that was just solved */
      SCIP_CALL( SCIPlpiExactGetState(pslpiexact, blkmem, &lpistate) );
      SCIP_CALL( SCIPlpiExactGetSol(pslpiexact, objval, NULL, NULL, NULL, NULL) );

      /* now reset the objective value to be the original objective */
      pos = 0;
      for( i = 0; i < nrows; i++ )
      {
         if( dvarincidence[i] )
         {
            RatSet(psobj[pos], lpexact->rows[i]->lhs);
            pos++;
         }
      }
      for( i = 0; i < nrows; i++ )
      {
         if( dvarincidence[nrows + i] )
         {
            RatNegate(psobj[pos], lpexact->rows[i]->rhs);
            pos++;
         }
      }
      for( i = 0; i < ncols; i++ )
      {
         if( dvarincidence[2*nrows + i] )
         {
            RatSet(psobj[pos], lpexact->cols[i]->lb);
            pos++;
         }
      }
      for( i = 0; i < ncols; i++ )
      {
         if( dvarincidence[2*nrows + ncols + i])
         {
            RatNegate(psobj[pos], lpexact->cols[i]->ub);
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
      SCIP_CALL( SCIPlpiExactChgBounds(pslpiexact, psncols, psind, pslb, NULL) );
      SCIP_CALL( SCIPlpiExactChgObj(pslpiexact, psncols, psind, psobj) );

      /* reload state and solve new LP */
      SCIP_CALL( SCIPlpiExactSetState(pslpiexact, blkmem, lpistate) );

      /* reoptimizing using primal simplex seems MUCH faster here, warm start basis is primal feasible */
      SCIP_CALL( SCIPlpiExactSolvePrimal(pslpiexact) );
      SCIP_CALL( SCIPlpiExactFreeState(pslpiexact, blkmem, &lpistate) );

        /* recover the optimal solution and set interior point and slack in constraint handler data */
      if( SCIPlpiExactIsOptimal(pslpiexact) )
      {
         SCIPdebugMessage("   exact LP solved to optimality\n");

         /* get optimal dual solution */
         SCIP_CALL( SCIPlpiExactGetSol(pslpiexact, objval, sol, NULL, NULL, NULL) );

         /* assign interior point solution to constraint handler data */
         RatSet(projshiftdata->commonslack, sol[psncols - 1]);
         if( RatIsZero(projshiftdata->commonslack) )
         {
            SCIPdebugMessage("   --> interior point not found \n");
         }
         else
         {
            for( i = 0; i < ndvarmap; i++ )
            {
               RatSet( projshiftdata->interiorpoint[dvarmap[i]], sol[i]);
            }
            projshiftdata->projshifthaspoint = TRUE;
         }
      }
      else
      {
         SCIP_CALL( printlpiexacterr(pslpiexact) );
      }
   }
   else
   {
      SCIPerrorMessage("invalid parameter setting <%d> for selection method to compute interior point\n",
         projshiftdata->psintpointselection);
      return SCIP_PARAMETERWRONGVAL;
   }

   for( i = 0; i < ndvarmap; i++ )
   {
      if( projshiftdata->projshifthaspoint )
         RatCanonicalize(projshiftdata->interiorpoint[i]);
      if( projshiftdata->projshifthasray )
         RatCanonicalize(projshiftdata->interiorray[i]);
   }

   /* free memory */
   if( pslpiexact != NULL )
   {
      int nlpirows, nlpicols;
      SCIP_CALL( SCIPlpiExactGetNRows(pslpiexact, &nlpirows) );
      SCIP_CALL( SCIPlpiExactDelRows(pslpiexact, 0, nlpirows - 1) );

      SCIP_CALL( SCIPlpiExactGetNCols(pslpiexact, &nlpicols) );
      SCIP_CALL( SCIPlpiExactDelCols(pslpiexact, 0, nlpicols - 1) );

      SCIP_CALL( SCIPlpiExactClear(pslpiexact) );
      SCIP_CALL( SCIPlpiExactFree(&pslpiexact) );
   }
   assert(pslpiexact == NULL);
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
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
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
   SCIP_PROJSHIFTDATA* projshiftdata;

   assert(lpexact != NULL);
   assert(lpexact->projshiftdata != NULL);

   projshiftdata = lpexact->projshiftdata;

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
   if( projshiftdata->projshiftdatacon )
      return SCIP_OKAY;

   SCIPclockStart(stat->provedfeaspstime, set);
   /* now mark that this function has been called */
   projshiftdata->projshiftdatacon = TRUE;

   SCIPdebugMessage("calling constructPSdata()\n");
   SCIP_CALL( RatCreateBlock(blkmem, &projshiftdata->commonslack) );

   /* process the bound changes */
   SCIP_CALL( SCIPsepastoreExactSyncLPs(set->scip->sepastoreexact, blkmem, set, stat, lpexact, eventqueue) );
   SCIP_CALL( SCIPlpExactFlush(lp->lpexact, blkmem, set, eventqueue) );

   assert(lpexact->nrows > 0);

   projshiftdata->nextendedrows = 2*lpexact->nrows + 2*lpexact->ncols;

   /* call function to select the set S */
   SCIP_CALL( psChooseS(lp, lpexact, set, stat, messagehdlr, eventqueue, eventfilter, prob, blkmem) );

   /* compute LU factorization of D == A|_S */
   SCIP_CALL( psFactorizeD(lp, lpexact, set, prob, blkmem, projshiftdata->projshiftuseintpoint) );

   /* if no fail in LU factorization, compute S-interior point and/or ray */
   if( !projshiftdata->projshiftdatafail )
   {
      if( TRUE || !projshiftdata->projshiftuseintpoint )
      {
         /* try to compute the S-interior ray if we want to use it for bounding or infeasibility */
         SCIP_CALL( RatCreateBlockArray(blkmem, &projshiftdata->interiorray, projshiftdata->nextendedrows) );
         SCIP_CALL( psComputeSintPointRay(lp, lpexact, set, prob, blkmem, FALSE) );
      }
      if( projshiftdata->projshiftuseintpoint || !projshiftdata->projshifthasray )
      {
         /* now, compute S-interior point if we need it OR if the ray construction failed */
         SCIP_CALL( RatCreateBlockArray(blkmem, &projshiftdata->interiorpoint, projshiftdata->nextendedrows) );
         SCIP_CALL( psComputeSintPointRay(lp, lpexact, set, prob, blkmem, TRUE) );
      }
   }

   /* if construction of both point and ray has failed, mark projshiftdatafail as true. */
   if( !projshiftdata->projshifthaspoint && !projshiftdata->projshifthasray )
   {
      projshiftdata->projshiftdatafail = TRUE;
   }

   SCIP_CALL( RatCreateBlockArray(blkmem, &projshiftdata->violation, lpexact->ncols) );
   SCIP_CALL( RatCreateBlockArray(blkmem, &projshiftdata->correction, projshiftdata->nextendedrows) );
   projshiftdata->violationsize = lpexact->ncols;

   SCIPclockStop(stat->provedfeaspstime, set);
   SCIPdebugMessage("exiting constructPSdata()\n");

   return SCIP_OKAY;
}

/** computes safe dual bound by project-and-shift method or corrects dual ray for infeasibility proof */
static
SCIP_RETCODE getPSdual(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_STAT*            stat,               /**< statistics pointer */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             usefarkas,          /**< to we aim to prove infeasibility ? */
   SCIP_Real*            safebound           /**< store the calculated safebound here */
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
   SCIP_PROJSHIFTDATA* projshiftdata;
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
   projshiftdata = lpexact->projshiftdata;

   assert(projshiftdata->projshiftdatacon);
   if( (projshiftdata->projshiftdatafail && !usefarkas) || (usefarkas && !projshiftdata->projshifthasray) )
   {
      lp->hasprovedbound = FALSE;
      return SCIP_OKAY;
   }

   /* constructprojshiftdata should always be called first */
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
   if( !usefarkas && projshiftdata->projshiftuseintpoint && projshiftdata->projshifthaspoint )
   {
      /* if we are supposed to use the interior point, and it exists use it */
      useinteriorpoint = TRUE;
   }
   else
   {
      /* in this case, since projshiftdatafail != TRUE, projshifthasray should be true -- use it */
      assert( projshiftdata->projshifthasray );
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
   SCIP_CALL( SCIPsepastoreExactSyncLPs(set->scip->sepastoreexact, blkmem, set, stat, lpexact, eventqueue) );
   SCIP_CALL( SCIPlpExactFlush(lp->lpexact, blkmem, set, eventqueue) );

   nextendedrows = projshiftdata->nextendedrows;
   nrows = lpexact->nrows;
   ncols = lpexact->ncols;
   nrowsps = projshiftdata->nextendedrows/2 - ncols;
   shift = nrows - nrowsps;

   assert(ncols == projshiftdata->violationsize);

   /* allocate memory for approximate dual solution, dual cost vector, violation and correction */
   SCIPsetAllocBufferArray(set, &approxdual, nrows + ncols);
   violation = projshiftdata->violation;
   correction = projshiftdata->correction;

   /** @todo: exip this could be removed */
   for( i = 0; i < nrows + ncols; ++i )
   {
      if( i < nrows )
         approxdual[i] = usefarkas ? lpexact->rows[i]->dualfarkas : lpexact->rows[i]->dualsol;
      else
         approxdual[i] = usefarkas ? lpexact->cols[i - nrows]->redcost : lpexact->cols[i - nrows]->farkascoef;

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
         val = i < nrows ? lpexact->rows[i]->lhs : lpexact->cols[i - nrows]->lb;
      else
         val = i < nrows ? lpexact->rows[i]->rhs : lpexact->cols[i - nrows]->ub;

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
         val = i < nrows ? lpexact->rows[i]->lhs : lpexact->cols[i - nrows]->lb;
      else
         val = i < nrows ? lpexact->rows[i]->rhs : lpexact->cols[i - nrows]->ub;

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
         RatDiff(violation[i], lpexact->cols[i]->obj, approxdual[i + nrows]);
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
      for( j = 0; j < lpexact->rows[i]->len; j++)
      {
         currentrow = lpexact->rows[i]->cols_index[j];
         RatMult(tmp, approxdual[i], lpexact->rows[i]->vals[j]);
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
      /* compute [z] with Dz=r (D depends on projshiftdata->psfpdualcolwiseselection) */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &violationgmp, ncols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &correctiongmp, nextendedrows) );
      RatSetGMPArray(violationgmp, violation, ncols);

      for ( i = 0; i < nextendedrows; i++)
      {
         mpq_init(correctiongmp[i]);
      }

      rval = RECTLUsolveSystem(projshiftdata->rectfactor, ncols, nextendedrows, violationgmp, correctiongmp);
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
      for( i = 0; i < projshiftdata->projshiftbasisdim; i++ )
      {
         printf("   i=%d: ", i);
         RatPrint(correction[i]);
         printf(", position=%d\n", projshiftdata->projshiftbasis[i]);
      }
#endif

      rectlut = clock();

      /* projection step: compute bold(y)=y^+[z 0];
       * save the corrected components in the correction vector
       */
      for( i = 0; i < projshiftdata->projshiftbasisdim; i++ )
      {
         // map is the point in the extended space -> transform it back to the original space
         int map = projshiftdata->projshiftbasis[i];
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
      for( i = 0; i < projshiftdata->projshiftbasisdim; i++ )
      {
         printf("   i=%d: ", i);
         RatPrint(correction[i]);
         printf(", position=%d\n", projshiftdata->projshiftbasis[i]);
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
         for( i = 0; i < projshiftdata->projshiftbasisdim; i++ )
         {
            if( RatIsNegative(correction[i]) )
            {
               int map = projshiftdata->projshiftbasis[i];

               RatSet(tmp2, projshiftdata->interiorpoint[map]);
               RatDiff(tmp, projshiftdata->interiorpoint[map], correction[i]);
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
         for( i = 0; i < projshiftdata->projshiftbasisdim; i++ )
         {
            int map = projshiftdata->projshiftbasis[i];
            if( RatIsNegative(correction[i]) && projshiftdata->includedrows[map] )
            {
               RatDiv(tmp, correction[i], projshiftdata->interiorray[map]);
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
      for( i = 0; i < projshiftdata->projshiftbasisdim; i++ )
      {
         int map = projshiftdata->projshiftbasis[i];
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
            RatMult(tmp, useinteriorpoint ? projshiftdata->interiorpoint[map] : projshiftdata->interiorray[map], lambda2);
            RatDiff(approxdual[i], approxdual[i], tmp);
            map = (i < nrowsps) ? i : i + nrowsps - shift;
            RatMult(tmp, useinteriorpoint ? projshiftdata->interiorpoint[map] : projshiftdata->interiorray[map], lambda2);
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
        RatSet(violation[i], lpexact->cols[i]->obj);
      else
         RatSetInt(violation[i], 0, 1);
   }
   for( i = 0; i < nrows; i++ )
   {
      for( j = 0; j < lpexact->rows[i]->len; j++ )
      {
         currentrow = lpexact->rows[i]->cols_index[j];
         RatMult(tmp, approxdual[i], lpexact->rows[i]->vals[j]);
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
         val = i < nrows ? lpexact->rows[i]->lhs : lpexact->cols[i - nrows]->lb;
      else
         val = i < nrows ? lpexact->rows[i]->rhs : lpexact->cols[i - nrows]->ub;

      RatMult(tmp, approxdual[i], val);
      RatAdd(dualbound, dualbound, tmp);
   }

   if( !usefarkas )
   {
      RatSet(lpexact->lpobjval, dualbound);
      lp->lpobjval = RatRoundReal(dualbound, SCIP_ROUND_DOWNWARDS);
      lp->hasprovedbound = TRUE;
   }
   else
   {
      /* if the objective value of the corrected ray is positive we can prune node, otherwise not */
      if( RatIsPositive(dualbound) )
      {
         RatSetString(lpexact->lpobjval, "inf");
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
   printf("   common slack=%.20f (", RgetRealApprox(projshiftdata->commonslack));
   RatPrint(projshiftdata->commonslack);
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
#endif

static
char chooseBoundingMethod(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEXACT*         lpexact,            /**< Exact LP data */
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

   assert(!lpexact->fplp->hasprovedbound);

   if( dualfarkas )
   {
      /* check if neumair-scher is possible */
      if( SCIPlpExactBSpossible(lpexact) )
         dualboundmethod = 'n';
      /* check if project and shift is possible */
#ifdef SCIP_WITH_GMP
      else if( SCIPlpExactPSpossible(lpexact) )
      {
         SCIP_CALL( constructPSData(lp, lpexact, set, stat, messagehdlr, eventqueue, eventfilter,
                        prob, blkmem) );
         if( !lpexact->projshiftdata->projshiftdatafail )
            dualboundmethod = 'p';
         else
            dualboundmethod = 'e';
      }
#endif
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
      if( (lpexact->interleavedbfreq > 0 && SCIPsetIsInfinity(set, SCIPlpGetCutoffbound(lpexact->fplp)) && SCIPgetDepth(set->scip) > 0
            && SCIPgetDepth(set->scip) % (lpexact->interleavedbfreq) == 0)
         || (lpexact->interleavedbfreq == 0 && SCIPsetIsGE(set, SCIPlpGetObjval(lpexact->fplp, set, prob), SCIPlpGetCutoffbound(lpexact->fplp))
            && SCIPlpGetObjval(lpexact->fplp, set, prob) < SCIPlpGetCutoffbound(lpexact->fplp)) )
      {
         dualboundmethod = 'e';
      }
      else
      {
         /* check if neumair-scher is possible */
         if( SCIPlpExactBSpossible(lpexact) )
            dualboundmethod = 'n';
#ifdef SCIP_WITH_GMP
         /* check if project and shift is possible */
         else if( SCIPlpExactPSpossible(lpexact) )
         {
            SCIP_CALL( constructPSData(lp, lpexact, set, stat, messagehdlr, eventqueue, eventfilter,
                           prob, blkmem) );
            if( !lpexact->projshiftdata->projshiftdatafail )
               dualboundmethod = 'p';
            else
               dualboundmethod = 'e';
         }
#endif
         /* otherwise solve exactly */
         else
            dualboundmethod = 'e';
      }
   }
   return dualboundmethod;
}

/** calculates a valid dual bound/farkas proof if all variables have lower and upper bounds 
 * Let (y,z) be the dual variables, y corresponding to primal rows, z to variable bounds. 
 * An exactly feasible dual solution is computed with y' = max{0,y}, z'=max{0,(c-A^Ty')}.
 * The bound is then computed as b^Ty'+s^Tz', with b being the lhs/rhs and s the lb/ub depending on the
 * sign of the dual value.
 * To avoid rational computations everything is done in interval arithmetic.
 */
static
SCIP_RETCODE boundShift(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEXACT*         lpexact,            /**< Exact LP data */
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
   SCIP_INTERVAL* rhslhsrow; // the rhs or lhs of row depending on sign of dual solution
   SCIP_INTERVAL* ublbcol; // the ub or lb of col depending on sign of dual solution
   SCIP_INTERVAL* constantinter;
   SCIP_INTERVAL* lpcolvals; // values in a column of the constraint matrix
   SCIP_INTERVAL* productcoldualval; // scalar product of lp column with dual solution vector
   SCIP_INTERVAL* obj; // objective (or 0 in farkas proof)
   SCIP_INTERVAL productsidedualval; //scalar product of sides with dual solution vector
   SCIP_INTERVAL safeboundinterval;
   SCIP_ROW* row;
   SCIP_COL* col;
   SCIP_Real* fpdual;
   SCIP_Real* fpdualcolwise;
   SCIP_Real c;
   int i;
   int j;

   assert(lpexact != NULL);
   assert(lp != NULL);
   assert(lp->solved);
   assert(set != NULL);
   assert(safebound != NULL);

   if( !SCIPlpExactBSpossible(lpexact) )
      return SCIP_OKAY;

   /* start timing */
   if ( usefarkas )
      SCIPclockStart(stat->provedinfeasbstime, set);
   else
      SCIPclockStart(stat->provedfeasbstime, set);

   /* allocate temporarfpdual memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &fpdual, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rhslhsrow, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &constantinter, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &fpdualcolwise, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lpcolvals, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &productcoldualval, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &obj, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ublbcol, lp->ncols) );

   SCIPdebugMessage("calling proved bound for %s LP\n", usefarkas ? "infeasible" : "feasible");

   /** @todo exip: actually we only need to link the rows and cols in the exact lp. So possible performance improvement if we don't
    * flush it to the lpiexact */
   SCIP_CALL( SCIPsepastoreExactSyncLPs(set->scip->sepastoreexact, blkmem, set, stat, lpexact, eventqueue) );
   SCIP_CALL( SCIPlpExactFlush(lpexact, blkmem, set, eventqueue) );

   /* reset proved bound status */
   lp->hasprovedbound = FALSE;

   /* calculate y^Tb */
   SCIPintervalSet(&productsidedualval, 0.0);
   SCIPdebugMessage("productsidedualval intervall computation with vectors:\n");

   /* create dual vector, sides and constant vector in interval arithmetic */
   for( j = 0; j < lp->nrows; ++j )
   {
      row = lp->rows[j];
      assert(row != NULL);

      /* create dual vector in interval arithmetic, setting near zeros to zero */
      fpdual[j] = (usefarkas ? row->dualfarkas : row->dualsol);

      if( SCIPlpiIsInfinity(lp->lpi, fpdual[j]) )
	      fpdual[j] = SCIPsetInfinity(set);

      if( SCIPlpiIsInfinity(lp->lpi, -fpdual[j]) )
	      fpdual[j] = -SCIPsetInfinity(set);

      /** @todo exip: dual bounding improvement
       *  - should we also set nonzero values of y to zero if corresponding lhs/rhs is not finite (to improve dual bound)?
       *  - do such situations come up?
       */
      /* create sides and constant vectors in interval arithmetic */
      if( SCIPsetIsFeasPositive(set, fpdual[j]) )
      {
         SCIPintervalSet(&rhslhsrow[j], row->lhs);
         SCIPintervalSet(&constantinter[j], -1.0 * row->constant);
      }
      else if( SCIPsetIsFeasNegative(set, fpdual[j]) )
      {
         SCIPintervalSet(&rhslhsrow[j], row->rhs);
         SCIPintervalSet(&constantinter[j], -1.0 * row->constant);
      }
      else
      {
         fpdual[j] = 0.0;
         SCIPintervalSet(&rhslhsrow[j], 0.0);
         SCIPintervalSet(&constantinter[j], 0.0);
      }

      SCIPdebugMessage("   j=%d: b=[%g,%g] (lhs=%g, rhs=%g, const=%g, fpdual=%g)\n", j, rhslhsrow[j].inf, rhslhsrow[j].sup, row->lhs,
            row->rhs, row->constant, fpdual[j]);
   }
   /* substract constant from sides in interval arithmetic and calculate fpdual * side */
   SCIPintervalAddVectors(SCIPsetInfinity(set), rhslhsrow, lp->nrows, rhslhsrow, constantinter);
   SCIPintervalScalprodScalars(SCIPsetInfinity(set), &productsidedualval, lp->nrows, rhslhsrow, fpdual);

   SCIPdebugMessage("   resulting scalar product=[%g,%g]\n", SCIPintervalGetInf(productsidedualval), SCIPintervalGetSup(productsidedualval));

   /* calculate min{(obj - dual^TMatrix)redcost} */

   /* compute infimums of -dual^TMatrix */
   roundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeDownwards();
   for( j = 0; j < lp->ncols; ++j )
   {
      SCIP_COLEXACT* colexact;

      col = lp->cols[j];
      assert(col != NULL);
      assert(col->nunlinked == 0);

      colexact = SCIPcolGetColExact(col);

      assert(colexact != NULL);

      /* create -Matrix.j vector in interval arithmetic and corresponding dual vector and compute infimum of vector -Matrix.j^Tdual */
      for( i = 0; i < colexact->nlprows; ++i )
      {
         SCIP_INTERVAL val;
         SCIP_ROWEXACT* rowexact;

         assert(colexact->rows[i] != NULL);
         assert(colexact->rows[i]->lppos >= 0);
         assert(colexact->linkpos[i] >= 0);

         rowexact = colexact->rows[i];

         val = rowexact->valsinterval[colexact->linkpos[i]];
         assert(RatIsGEReal(colexact->vals[i], val.inf) && RatIsLEReal(colexact->vals[i], val.sup));

         SCIPintervalSetBounds(&lpcolvals[i], -val.sup, -val.inf);
         fpdualcolwise[i] = fpdual[colexact->rows[i]->lppos];
      }
      productcoldualval[j].inf = 0.0;
      SCIPintervalScalprodScalarsInf(SCIPsetInfinity(set), &productcoldualval[j], col->nlprows, lpcolvals, fpdualcolwise);

#ifndef NDEBUG
      for( i = colexact->nlprows; i < col->len; ++i )
      {
         assert(colexact->rows[i] != NULL);
         assert(colexact->rows[i]->lppos == -1);
         assert(colexact->linkpos[i] >= 0);
      }
#endif
   }

   /* compute supremums of -dual^TMatrix */
   SCIPintervalSetRoundingModeUpwards();
   for( j = 0; j < lp->ncols; ++j )
   {
      SCIP_COLEXACT* colexact;

      col = lp->cols[j];
      assert(col != NULL);
      assert(col->nunlinked == 0);

      /* colexact can be longer than col (epsilon entries) so we have to work with colexact */
      colexact = SCIPcolGetColExact(col);

      assert(colexact->nlprows >= col->nlprows);

      /* create -Matrix.j vector in interval arithmetic and corresponding dual vector and compute supremums of vector -a.j^Ty */
      for( i = 0; i < colexact->nlprows; ++i )
      {
         SCIP_INTERVAL val;
         SCIP_ROWEXACT* rowexact;

         assert(colexact->rows[i] != NULL);
         assert(colexact->rows[i]->lppos >= 0);
         assert(colexact->linkpos[i] >= 0);

         rowexact = colexact->rows[i];

         val = rowexact->valsinterval[colexact->linkpos[i]];

         assert(RatIsGEReal(colexact->vals[i], val.inf) && RatIsLEReal(colexact->vals[i], val.sup));

         SCIPintervalSetBounds(&lpcolvals[i], -val.sup, -val.inf);
         fpdualcolwise[i] = fpdual[colexact->rows[i]->lppos];
      }
      productcoldualval[j].sup = 0.0;
      SCIPintervalScalprodScalarsSup(SCIPsetInfinity(set), &productcoldualval[j], col->nlprows, lpcolvals, fpdualcolwise);

#ifndef NDEBUG
      for( i = col->nlprows; i < col->len; ++i )
      {
         assert(col->rows[i] != NULL);
         assert(col->rows[i]->lppos == -1);
         assert(col->linkpos[i] >= 0);
      }
#endif
   }
   SCIPintervalSetRoundingMode(roundmode);

   /* create objective vector and lb/ub vector in interval arithmetic and compute min{(obj^T - dual^TMatrix)lb/ub} */
   for( j = 0; j < lp->ncols; ++j )
   {
      //assert(!SCIPsetIsInfinity(set, -SCIPcolGetLb(col)));
      //assert(!SCIPsetIsInfinity(set, SCIPcolGetUb(col)));
      col = lp->cols[j];
      assert(col != NULL);
      assert(col->nunlinked == 0);

      if( usefarkas )
         SCIPintervalSet(&obj[j], 0);
      else
      {
         if( RatIsFpRepresentable(SCIPvarGetObjExact(SCIPcolGetVar(col))) )
            SCIPintervalSet(&obj[j], col->obj);
         else
         {
            SCIPintervalSetRational(&obj[j], SCIPvarGetObjExact(SCIPcolGetVar(col)));
         }
      }

      assert(SCIPcolGetLb(col) <= RatRoundReal(SCIPvarGetLbLocalExact(col->var), SCIP_ROUND_DOWNWARDS));
      assert(SCIPcolGetUb(col) >= RatRoundReal(SCIPvarGetUbLocalExact(col->var), SCIP_ROUND_UPWARDS));
      SCIPintervalSetBounds(&ublbcol[j], SCIPcolGetLb(col), SCIPcolGetUb(col));

      /* opt out if there are infinity bounds and a non-infinte value */
      if( (SCIPsetIsInfinity(set, -SCIPcolGetLb(col)) || SCIPsetIsInfinity(set, SCIPcolGetUb(col))) )
      {
         if( productcoldualval[j].inf + obj[j].inf != 0 || productcoldualval[j].sup + obj[j].sup != 0 )
         {
            SCIPdebugMessage("trying bound shift with unbounded column variable. Column %d, lb: %e, ub %e \n",
               SCIPcolGetIndex(col), SCIPcolGetLb(col) ,SCIPcolGetUb(col) );
            SCIPdebugMessage("Multiplied with interval: min %e,  max %e \n",
               productcoldualval[j].inf + obj[j].inf, productcoldualval[j].sup + obj[j].sup);

            lp->hasprovedbound = FALSE;
            if( usefarkas )
            {
               stat->nboundshiftinf++;
               stat->nfailboundshiftinf++;
            }
            else
            {
               stat->nboundshift++;
               stat->nfailboundshift++;
            }
            goto CLEANUP;
         }
      }
   }
   SCIPintervalAddVectors(SCIPsetInfinity(set), productcoldualval, lp->ncols, productcoldualval, obj);
   SCIPintervalScalprod(SCIPsetInfinity(set), &safeboundinterval, lp->ncols, productcoldualval, ublbcol);

   /* add dualsol * rhs/lhs (or farkas * rhs/lhs) */
   SCIPintervalAdd(SCIPsetInfinity(set), &safeboundinterval, safeboundinterval, productsidedualval);

   *safebound = SCIPintervalGetInf(safeboundinterval);
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

   /* if certificate is active, save the corrected dual solution into the lpexact data */
   if( SCIPisCertificateActive(set->scip) && lp->hasprovedbound )
   {
      SCIP_INTERVAL tmp, tmp2;
      SCIP_Real cand1, cand2;
      SCIP_Real value;
      /* set up the exact lpi for the current node */
      SCIP_CALL( SCIPsepastoreExactSyncLPs(set->scip->sepastoreexact, blkmem, set, stat, lpexact, eventqueue) );
      SCIP_CALL( SCIPlpExactFlush(lp->lpexact, blkmem, set, eventqueue) );
      for( j = 0; j < lpexact->nrows; j++ )
      {
         if( usefarkas )
            RatSetReal(lpexact->rows[j]->dualfarkas, fpdual[j]);
         else
            RatSetReal(lpexact->rows[j]->dualsol, fpdual[j]);
      }
      for( j = 0; j < lpexact->ncols; j++ )
      {
         cand1 = productcoldualval[j].inf;
         cand2 = productcoldualval[j].sup;
         SCIPintervalMulScalar(SCIPsetInfinity(set), &tmp, ublbcol[j], cand1);
         SCIPintervalMulScalar(SCIPsetInfinity(set), &tmp2, ublbcol[j], cand2);
         if( ((tmp.inf) < (tmp2.inf)) == RatIsPositive(lpexact->cols[j]->obj))
            value = cand1;
         else
            value = cand2;

         if( usefarkas )
            RatSetReal(lpexact->cols[j]->farkascoef, value);
         else
            RatSetReal(lpexact->cols[j]->redcost, value);
      }
   }

   RatSetReal(lpexact->lpobjval, SCIPintervalGetInf(safeboundinterval));

CLEANUP:

   /* if the fail percentage is higher than 20 % we do not want to waste time trying bound shift again and again */
   if( stat->nboundshift + stat->nboundshiftinf > 10
      && (1.0 * stat->nfailboundshift + stat->nfailboundshiftinf) / stat->nboundshift + stat->nboundshiftinf > 0.2 )
   {
      lpexact->boundshiftviable = FALSE;
   }
   /* free buffer for storing y in interval arithmetic */
   SCIPsetFreeBufferArray(set, &ublbcol);
   SCIPsetFreeBufferArray(set, &obj);
   SCIPsetFreeBufferArray(set, &productcoldualval);
   SCIPsetFreeBufferArray(set, &lpcolvals);
   SCIPsetFreeBufferArray(set, &fpdualcolwise);
   SCIPsetFreeBufferArray(set, &constantinter);
   SCIPsetFreeBufferArray(set, &rhslhsrow);
   SCIPsetFreeBufferArray(set, &fpdual);

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
#endif

SCIP_RETCODE SCIPlpExactComputeSafeBound(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEXACT*         lpexact,            /**< Exact LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Longint          itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool             dualfarkas,         /**< should infeasiblity be proven? */
   SCIP_Real*            safebound,          /**< store the calculated safebound here */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            dualfeasible        /**< pointer to store whether the solution is dual feasible, or NULL */
   )
{
   char dualboundmethod;

   /* If we are not in exact solving mode, just return */
   if( !set->misc_exactsolve || lp->diving || lp->probing || lp->strongbranchprobing )
      return SCIP_OKAY;

#ifdef SCIP_WITH_BOOST
   assert(set->misc_exactsolve);
   assert(!lp->hasprovedbound);

   /* if the lp was solved to optimality and there are no fractional variables we solve exactly to generate a feasible solution */
   if( lpexact->forceexactsolve || (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && fpLPisIntFeasible(lp, set, stat)) )
      lpexact->forceexactsolve = TRUE;

   /* if the lp is fp-objlim then we solve exactly to hopefully achieve the same (exact) result
      @todo exip maybe it makes more sense to just disable objlimit alltogether? */
   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT )
      lpexact->forceexactsolve = TRUE;

   /* choose which bounding method to use. only needed if automatic is enabled. */
   if( lpexact->forceexactsolve )
   {
      dualboundmethod = 'e';
      lpexact->forceexactsolve = FALSE;
   }
   else if( set->misc_dbmethod == 'a' )
   {
      dualboundmethod = chooseBoundingMethod(lp, lpexact, set, messagehdlr, blkmem, stat, eventqueue, eventfilter,
                              prob, itlim, lperror, dualfarkas, safebound);
   }
   else
      dualboundmethod = set->misc_dbmethod;

   SCIPdebugMessage("Computing safe bound for lp with status: %d using bounding method %c \n", SCIPlpGetSolstat(lp), dualboundmethod);

   switch(dualboundmethod)
   {
      /* neumaier and scherbina */
      case 'n':
         SCIP_CALL( boundShift(lp, lpexact, set, messagehdlr, blkmem, stat, eventqueue, eventfilter,
                        prob, dualfarkas, safebound) );
         if( !lp->hasprovedbound )
         {
            SCIP_CALL( constructPSData(lp, lpexact, set, stat, messagehdlr, eventqueue, eventfilter,
                        prob, blkmem) );
            SCIP_CALL( getPSdual(lp, lpexact, set, stat, messagehdlr, eventqueue, eventfilter,
                        prob, blkmem, dualfarkas, safebound) );
         }
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
#ifdef SCIP_WITH_GMP
      /* project and shift */
      case 'p':
         SCIP_CALL( constructPSData(lp, lpexact, set, stat, messagehdlr, eventqueue, eventfilter,
                        prob, blkmem) );
         SCIP_CALL( getPSdual(lp, lpexact, set, stat, messagehdlr, eventqueue, eventfilter,
                        prob, blkmem, dualfarkas, safebound) );
         break;
#endif
      /* exact LP */
      case 'e':
         SCIP_CALL( SCIPlpExactSolveAndEval(lpexact, lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, prob, set->lp_iterlim,
               lperror, dualfarkas) );
         *primalfeasible = lpexact->primalfeasible;
         *dualfeasible = lpexact->dualfeasible;

         break;
      default:
         SCIPerrorMessage("bounding method %c not implemented yet \n", set->misc_dbmethod);
         SCIPABORT();
         break;
   }

   if( !lp->hasprovedbound && !(*lperror) )
   {
      SCIP_CALL( SCIPlpExactSolveAndEval(lpexact, lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, prob, set->lp_iterlim,
               lperror, dualfarkas) );
      *primalfeasible = lpexact->primalfeasible;
      *dualfeasible = lpexact->dualfeasible;
   }
#endif

   /* choose which bounding method should be calles and return a safe objective bound */
   return SCIP_OKAY;
}

#endif