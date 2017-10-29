/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file   presol_sparsify.c
 * @brief  cancel non-zeros of the constraint matrix
 * @author Dieter Weninger
 * @author Robert Lion Gottwald
 *
 * This presolver attempts to cancel non-zero entries of the constraint
 * matrix by adding equalities to other constraints.
 * It supports two cases:
 * a) the variable index set from non-zero entries of the other constraint
 *    is a superset of the variable index set of the non-zero entries of the
 *    equality
 * b) same as case a) except for one variable index
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/pub_matrix.h"
#include "scip/cons_linear.h"
#include "scip/presol_sparsify.h"


#define PRESOL_NAME            "sparsify"
#define PRESOL_DESC            "eliminate non-zero coefficients"

#define PRESOL_PRIORITY            -24000    /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               -1    /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_MAX_CONT_FILLIN         0    /**< default value for the maximal fillin for continuous variables */
#define DEFAULT_MAX_BIN_FILLIN          0    /**< default value for the maximal fillin for binary variables */
#define DEFAULT_MAX_INT_FILLIN          0    /**< default value for the maximal fillin for integer variables */
#define DEFAULT_MAXNONZEROS            70    /**< maximal support of one equality to be used for cancelling (-1: no limit) */
#define DEFAULT_MAXCONSIDEREDNONZEROS  70    /**< maximal number of considered non-zeros within one row (-1: no limit) */
#define DEFAULT_ROWSORT               'd'    /**< order in which to process inequalities ('n'o sorting, 'i'ncreasing nonzeros, 'd'ecreasing nonzeros) */
#define DEFAULT_MAXRETRIEVEFAC      100.0    /**< limit on the number of useless vs. useful hashtable retrieves as a multiple of the number of constraints */

#define MAXSCALE                   1000.0    /**< maximal allowed scale for cancelling non-zeros */


/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   int                   nfailures;          /**< number of calls to presolver without success */
   int                   nwaitingcalls;      /**< number of presolver calls until next real execution */
   int                   maxcontfillin;      /**< maximal fillin for continuous variables */
   int                   maxintfillin;       /**< maximal fillin for integer variables*/
   int                   maxbinfillin;       /**< maximal fillin for binary variables */
   int                   maxnonzeros;        /**< maximal support of one equality to be used for cancelling (-1: no limit) */
   int                   maxconsiderednonzeros;/**< maximal number of considered non-zeros within one row (-1: no limit) */
   SCIP_Real             maxretrievefac;     /**< limit on the number of useless vs. useful hashtable retrieves as a multiple of the number of constraints */
   char                  rowsort;            /**< order in which to process inequalities ('n'o sorting, 'i'ncreasing nonzeros, 'd'ecreasing nonzeros) */
};


/*
 * Local methods
 */

struct RowVarPair
{
   int rowindex;
   int varindex1;
   int varindex2;
   SCIP_Real varcoef1;
   SCIP_Real varcoef2;
};

typedef struct RowVarPair ROWVARPAIR;

/** returns TRUE iff both keys are equal */
static
SCIP_DECL_HASHKEYEQ(varPairsEqual)
{  /*lint --e{715}*/
   ROWVARPAIR* varpair1;
   ROWVARPAIR* varpair2;
   SCIP_Real ratio1;
   SCIP_Real ratio2;

   varpair1 = (ROWVARPAIR*) key1;
   varpair2 = (ROWVARPAIR*) key2;

   if( varpair1->varindex1 != varpair2->varindex1 )
      return FALSE;

   if( varpair1->varindex2 != varpair2->varindex2 )
      return FALSE;

   ratio1 = varpair1->varcoef2 / varpair1->varcoef1;
   ratio2 = varpair2->varcoef2 / varpair2->varcoef1;

   if( !EPSEQ(ratio1, ratio2, SCIP_DEFAULT_EPSILON) )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(varPairHashval)
{  /*lint --e{715}*/
   ROWVARPAIR* varpair;

   varpair = (ROWVARPAIR*) key;

   return SCIPhashTwo(SCIPcombineTwoInt(varpair->varindex1, varpair->varindex2),
                      SCIPrealHashCode(varpair->varcoef2 / varpair->varcoef1));
}

/** non-zero cancellation of rows */
static
SCIP_RETCODE cancelRow(
   SCIP*                 scip,
   SCIP_MATRIX*          matrix,
   SCIP_HASHTABLE*       pairtable,
   int                   rowidx,
   int                   maxcontfillin,
   int                   maxintfillin,
   int                   maxbinfillin,
   int                   maxconsiderednonzeros,
   SCIP_Longint*         nuseless,
   int*                  nchgcoefs,
   int*                  ncanceled,
   int*                  nfillin
   )
{
   int* cancelrowinds;
   SCIP_Real* cancelrowvals;
   SCIP_Real cancellhs;
   SCIP_Real cancelrhs;
   int* tmpinds;
   int* locks;
   SCIP_Real* tmpvals;
   int cancelrowlen;
   int* rowidxptr;
   SCIP_Real* rowvalptr;
   int nchgcoef;
   int nretrieves;
   int bestnfillin;
   SCIP_Bool rowiseq;

   rowiseq = SCIPisEQ(scip, SCIPmatrixGetRowLhs(matrix, rowidx), SCIPmatrixGetRowRhs(matrix, rowidx));

   cancelrowlen = SCIPmatrixGetRowNNonzs(matrix, rowidx);
   rowidxptr = SCIPmatrixGetRowIdxPtr(matrix, rowidx);
   rowvalptr = SCIPmatrixGetRowValPtr(matrix, rowidx);

   SCIP_CALL( SCIPduplicateBufferArray(scip, &cancelrowinds, rowidxptr, cancelrowlen) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &cancelrowvals, rowvalptr, cancelrowlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpinds, cancelrowlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvals, cancelrowlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &locks, cancelrowlen) );

   cancellhs = SCIPmatrixGetRowLhs(matrix, rowidx);
   cancelrhs = SCIPmatrixGetRowRhs(matrix, rowidx);

   nchgcoef = 0;
   nretrieves = 0;
   while( TRUE ) /*lint !e716 */
   {
      SCIP_Real bestscale;
      int bestcand;
      int bestncancel;
      SCIP_Bool bestdeccond;
      int i;
      int j;
      ROWVARPAIR rowvarpair;
      int maxlen;

      bestscale = 1.0;
      bestcand = -1;
      bestncancel = 0;
      bestnfillin = 0;
      bestdeccond = FALSE;

      for( i = 0; i < cancelrowlen; ++i )
      {
         tmpinds[i] = i;
         locks[i] = SCIPmatrixGetColNDownlocks(matrix, cancelrowinds[i]) + SCIPmatrixGetColNUplocks(matrix, cancelrowinds[i]);
      }

      SCIPsortIntInt(locks, tmpinds, cancelrowlen);

      maxlen = cancelrowlen;
      if( maxconsiderednonzeros >= 0 )
         maxlen = MIN(cancelrowlen, maxconsiderednonzeros);

      for( i = 0; i < maxlen; ++i )
      {
         for( j = i + 1; j < maxlen; ++j )
         {
            int a,b;
            int ncancel;
            int ncontfillin;
            int nintfillin;
            int nbinfillin;
            int ntotfillin;
            int eqrowlen;
            ROWVARPAIR* eqrowvarpair;
            SCIP_Real* eqrowvals;
            int* eqrowinds;
            SCIP_Real eqrowside;
            SCIP_Real scale;
            int i1,i2;
            SCIP_Bool deccond;

            i1 = tmpinds[i];
            i2 = tmpinds[j];

            assert(cancelrowinds[i] < cancelrowinds[j]);

            if( cancelrowinds[i1] < cancelrowinds[i2] )
            {
               rowvarpair.varindex1 = cancelrowinds[i1];
               rowvarpair.varindex2 = cancelrowinds[i2];
               rowvarpair.varcoef1 = cancelrowvals[i1];
               rowvarpair.varcoef2 = cancelrowvals[i2];
            }
            else
            {
               rowvarpair.varindex1 = cancelrowinds[i2];
               rowvarpair.varindex2 = cancelrowinds[i1];
               rowvarpair.varcoef1 = cancelrowvals[i2];
               rowvarpair.varcoef2 = cancelrowvals[i1];
            }

            eqrowvarpair = (ROWVARPAIR*)SCIPhashtableRetrieve(pairtable, (void*) &rowvarpair);
            nretrieves++;

            if( eqrowvarpair == NULL || eqrowvarpair->rowindex == rowidx )
               continue;

            /* if the row we want to cancel is an equality, we will only use equalities
             * for canceling with less non-zeros and if the number of non-zeros is equal we use the
             * rowindex as tie-breaker to avoid cyclic non-zero cancellation
             */
            eqrowlen = SCIPmatrixGetRowNNonzs(matrix, eqrowvarpair->rowindex);
            if( rowiseq && (cancelrowlen < eqrowlen || (cancelrowlen == eqrowlen && rowidx < eqrowvarpair->rowindex)) )
               continue;

            eqrowvals = SCIPmatrixGetRowValPtr(matrix, eqrowvarpair->rowindex);
            eqrowinds = SCIPmatrixGetRowIdxPtr(matrix, eqrowvarpair->rowindex);
            eqrowside = SCIPmatrixGetRowRhs(matrix, eqrowvarpair->rowindex);

            scale = -rowvarpair.varcoef1 / eqrowvarpair->varcoef1;

            if( scale > MAXSCALE )
               continue;

            a = 0;
            b = 0;
            ncancel = 0;
            ncontfillin = 0;
            nintfillin = 0;
            nbinfillin = 0;
            while( a < cancelrowlen && b < eqrowlen )
            {
               if( cancelrowinds[a] == eqrowinds[b] )
               {
                  if( SCIPisZero(scip, cancelrowvals[a] + scale * eqrowvals[b]) )
                     ++ncancel;

                  ++a;
                  ++b;
               }
               else if( cancelrowinds[a] < eqrowinds[b] )
               {
                  ++a;
               }
               else
               {
                  SCIP_VAR* var = SCIPmatrixGetVar(matrix, eqrowinds[b]);
                  ++b;
                  if( SCIPvarIsIntegral(var) )
                  {
                     if( ++nintfillin > maxintfillin )
                        break;
                     if( SCIPvarIsBinary(var) )
                     {
                        if( ++nbinfillin > maxbinfillin )
                           break;
                     }
                  }
                  else
                  {
                     if( ++ncontfillin > maxcontfillin )
                        break;
                  }
               }
            }

            if( ncontfillin > maxcontfillin || nbinfillin > maxbinfillin || nintfillin > maxintfillin )
               continue;

            while( b < eqrowlen )
            {
               SCIP_VAR* var = SCIPmatrixGetVar(matrix, eqrowinds[b]);
               ++b;
               if( SCIPvarIsBinary(var) )
               {
                  if( ++nbinfillin > maxbinfillin )
                     break;
               }
               else if( SCIPvarIsIntegral(var) )
               {
                  if( ++nintfillin > maxintfillin )
                     break;
               }
               else
               {
                  if( ++ncontfillin > maxcontfillin )
                     break;
               }
            }

            if( ncontfillin > maxcontfillin || nbinfillin > maxbinfillin || nintfillin > maxintfillin )
               continue;

            /* secure that integer lhs stays integer after scaling */
            if( !SCIPisInfinity(scip, -cancellhs) )
               if( SCIPisIntegral(scip, cancellhs) && !SCIPisIntegral(scip, cancellhs + scale * eqrowside) )
                  continue;

            /* secure that integer rhs stays integer after scaling */
            if( !SCIPisInfinity(scip, cancelrhs) )
               if( SCIPisIntegral(scip, cancelrhs) && !SCIPisIntegral(scip, cancelrhs + scale * eqrowside) )
                  continue;

            ntotfillin = nbinfillin + nintfillin + ncontfillin;

            /* the case if all non-zeros of the equation for cancelling are canceled in the
             * cancel row is a necessary condition for generating independent components
             */
            deccond = (ntotfillin == 0 && ncancel == eqrowlen);

            if( (deccond && !bestdeccond) ||
                ((ncancel - ntotfillin) >  (bestncancel - bestnfillin) ||
                ((ncancel - ntotfillin) == (bestncancel - bestnfillin) && ntotfillin < bestnfillin) ) )
            {
               bestncancel = ncancel;
               bestnfillin = ntotfillin;
               bestcand = eqrowvarpair->rowindex;
               bestscale = scale;
               bestdeccond = deccond;
            }
         }
      }

      if( bestcand != -1 )
      {
         int a;
         int b;
         SCIP_Real* eqrowvals;
         int* eqrowinds;
         int eqrowlen;
         int tmprowlen;
         SCIP_Real eqrhs;

         eqrowvals = SCIPmatrixGetRowValPtr(matrix, bestcand);
         eqrowinds = SCIPmatrixGetRowIdxPtr(matrix, bestcand);
         eqrowlen = SCIPmatrixGetRowNNonzs(matrix, bestcand);
         eqrhs = SCIPmatrixGetRowRhs(matrix, bestcand);

         a = 0;
         b = 0;
         tmprowlen = 0;

         if( !SCIPisZero(scip, eqrhs) )
         {
            if( !SCIPisInfinity(scip, -cancellhs) )
               cancellhs += bestscale * eqrhs;
            if( !SCIPisInfinity(scip, cancelrhs) )
               cancelrhs += bestscale * eqrhs;
         }

         while( a < cancelrowlen && b < eqrowlen )
         {
            if( cancelrowinds[a] == eqrowinds[b] )
            {
               SCIP_Real val = cancelrowvals[a] + bestscale * eqrowvals[b];

               if( !SCIPisZero(scip, val) )
               {
                  tmpinds[tmprowlen] = cancelrowinds[a];
                  tmpvals[tmprowlen] = val;
                  ++tmprowlen;
               }
               ++nchgcoef;

               ++a;
               ++b;
            }
            else if( cancelrowinds[a] < eqrowinds[b] )
            {
               tmpinds[tmprowlen] = cancelrowinds[a];
               tmpvals[tmprowlen] = cancelrowvals[a];
               ++tmprowlen;
               ++a;
            }
            else
            {
               tmpinds[tmprowlen] = eqrowinds[b];
               tmpvals[tmprowlen] = eqrowvals[b] * bestscale;
               ++nchgcoef;
               ++tmprowlen;
               ++b;
            }
         }

         while( a < cancelrowlen )
         {
            tmpinds[tmprowlen] = cancelrowinds[a];
            tmpvals[tmprowlen] = cancelrowvals[a];
            ++tmprowlen;
            ++a;
         }

         while( b < eqrowlen )
         {
            tmpinds[tmprowlen] = eqrowinds[b];
            tmpvals[tmprowlen] = eqrowvals[b] * bestscale;
            ++nchgcoef;
            ++tmprowlen;
            ++b;
         }

         SCIPswapPointers((void**) &tmpinds, (void**) &cancelrowinds);
         SCIPswapPointers((void**) &tmpvals, (void**) &cancelrowvals);
         cancelrowlen = tmprowlen;
      }
      else
         break;
   }

   if( nchgcoef != 0 )
   {
      SCIP_CONS* cons;
      SCIP_VAR** consvars;

      int i;

      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, cancelrowlen) );

      for( i = 0; i < cancelrowlen; ++i )
         consvars[i] = SCIPmatrixGetVar(matrix, cancelrowinds[i]);

      /* create sparsified constraint and add it to scip */
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, SCIPmatrixGetRowName(matrix, rowidx), cancelrowlen, consvars, cancelrowvals,
                                      cancellhs, cancelrhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPdelCons(scip, SCIPmatrixGetCons(matrix, rowidx)) );
      SCIP_CALL( SCIPaddCons(scip, cons) );

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "########\n");
      SCIPdebugMsg(scip, "old:\n");
      SCIPmatrixPrintRow(scip, matrix, rowidx);
      SCIPdebugMsg(scip, "new:\n");
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIPdebugMsg(scip, "########\n");
#endif

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* update counters */
      *nchgcoefs += nchgcoef;
      *ncanceled += SCIPmatrixGetRowNNonzs(matrix, rowidx) - cancelrowlen;
      *nfillin += bestnfillin;

      /* if successful, decrease the useless hashtable retrieves counter; the rationale here is that we want to keep
       * going if, after many useless calls that almost exceeded the budget, we finally reach a useful section; but we
       * don't allow a negative build-up for the case that the useful section is all at the beginning and we just want
       * to quit quickly afterwards
       */
      *nuseless -= nretrieves;
      *nuseless = MAX(*nuseless, 0);

      SCIPfreeBufferArray(scip, &consvars);
   }
   else
   {
      /* if not successful, increase useless hashtable retrieves counter */
      *nuseless += nretrieves;
   }

   SCIPfreeBufferArray(scip, &locks);
   SCIPfreeBufferArray(scip, &tmpvals);
   SCIPfreeBufferArray(scip, &tmpinds);
   SCIPfreeBufferArray(scip, &cancelrowvals);
   SCIPfreeBufferArray(scip, &cancelrowinds);

   return SCIP_OKAY;
}

/** updates failure counter after one execution */
static
void updateFailureStatistic(
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_Bool             success             /**< was this execution successful? */
   )
{
   assert(presoldata != NULL);

   if( success )
   {
      presoldata->nfailures = 0;
      presoldata->nwaitingcalls = 0;
   }
   else
   {
      presoldata->nfailures++;
      presoldata->nwaitingcalls = 2*presoldata->nfailures;
   }
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecSparsify)
{  /*lint --e{715}*/
   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;
   int nrows;
   int r;
   int i;
   int j;
   int numcancel;
   int nfillin;
   int* locks;
   int* perm;
   int* rowidxsorted;
   int* rowsparsity;
   SCIP_HASHTABLE* pairtable;
   ROWVARPAIR* varpairs;
   int nvarpairs;
   int varpairssize;
   SCIP_PRESOLDATA* presoldata;
   SCIP_Longint maxuseless;
   SCIP_Longint nuseless;
   int noldchgcoefs;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   presoldata = SCIPpresolGetData(presol);

   if( presoldata->nwaitingcalls > 0 )
   {
      presoldata->nwaitingcalls--;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   /* we only work on pure MIPs currently */
   if( initialized && complete )
   {
      nrows = SCIPmatrixGetNRows(matrix);

      /* sort rows by column indices */
      for( i = 0; i < nrows; i++ )
      {
         int* rowpnt = SCIPmatrixGetRowIdxPtr(matrix, i);
         SCIP_Real* valpnt = SCIPmatrixGetRowValPtr(matrix, i);
         SCIPsortIntReal(rowpnt, valpnt, SCIPmatrixGetRowNNonzs(matrix, i));
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &locks, SCIPmatrixGetNColumns(matrix)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, SCIPmatrixGetNColumns(matrix)) );

      /* loop over all rows and create var pairs */
      numcancel = 0;
      nfillin = 0;
      varpairssize = 0;
      nvarpairs = 0;
      varpairs = NULL;
      SCIP_CALL( SCIPhashtableCreate(&pairtable, SCIPblkmem(scip), 1, SCIPhashGetKeyStandard, varPairsEqual, varPairHashval, NULL) );

      /* collect equalities and their number of non-zeros */
      for( r = 0; r < nrows; r++ )
      {
         int nnonz;

         nnonz = SCIPmatrixGetRowNNonzs(matrix, r);

         /* consider equalities with support at most maxnonzeros; skip singleton equalities, because these are faster
          * processed by trivial presolving
          */
         if( nnonz >= 2 && (presoldata->maxnonzeros < 0 || nnonz <= presoldata->maxnonzeros)
            && SCIPisEQ(scip, SCIPmatrixGetRowRhs(matrix, r), SCIPmatrixGetRowLhs(matrix, r)) )
         {
            int* rowinds;
            SCIP_Real* rowvals;
            int npairs;

            rowinds = SCIPmatrixGetRowIdxPtr(matrix, r);
            rowvals = SCIPmatrixGetRowValPtr(matrix, r);

            for( i = 0; i < nnonz; ++i )
            {
               perm[i] = i;
               locks[i] = SCIPmatrixGetColNDownlocks(matrix, rowinds[i]) + SCIPmatrixGetColNUplocks(matrix, rowinds[i]);
            }

            SCIPsortIntInt(locks, perm, nnonz);

            if( presoldata->maxconsiderednonzeros >= 0 )
               nnonz = MIN(nnonz, presoldata->maxconsiderednonzeros);

            npairs = (nnonz * (nnonz - 1)) / 2;

            if( nvarpairs + npairs > varpairssize )
            {
               int newsize = SCIPcalcMemGrowSize(scip, nvarpairs + npairs);
               SCIP_CALL( SCIPreallocBufferArray(scip, &varpairs, newsize) );
               varpairssize = newsize;
            }

            for( i = 0; i < nnonz; ++i )
            {
               for( j = i + 1; j < nnonz; ++j )
               {
                  int i1,i2;

                  assert(nvarpairs < varpairssize);
                  assert(varpairs != NULL);

                  i1 = perm[i];
                  i2 = perm[j];
                  varpairs[nvarpairs].rowindex = r;

                  if( rowinds[i1] < rowinds[i2])
                  {
                     varpairs[nvarpairs].varindex1 = rowinds[i1];
                     varpairs[nvarpairs].varindex2 = rowinds[i2];
                     varpairs[nvarpairs].varcoef1 = rowvals[i1];
                     varpairs[nvarpairs].varcoef2 = rowvals[i2];
                  }
                  else
                  {
                     varpairs[nvarpairs].varindex1 = rowinds[i2];
                     varpairs[nvarpairs].varindex2 = rowinds[i1];
                     varpairs[nvarpairs].varcoef1 = rowvals[i2];
                     varpairs[nvarpairs].varcoef2 = rowvals[i1];
                  }
                  ++nvarpairs;
               }
            }
         }
      }

      /* insert varpairs into hash table */
      for( r = 0; r < nvarpairs; ++r )
      {
         ROWVARPAIR* othervarpair;

         assert(varpairs != NULL);

         othervarpair = (ROWVARPAIR*)SCIPhashtableRetrieve(pairtable, (void*) &varpairs[r]);

         if( othervarpair != NULL && SCIPmatrixGetRowNNonzs(matrix, othervarpair->rowindex) <= SCIPmatrixGetRowNNonzs(matrix, varpairs[r].rowindex) )
            continue;

         SCIP_CALL( SCIPhashtableInsert(pairtable, (void*) &varpairs[r]) );
      }

      /* sort rows according to parameter value */
      if( presoldata->rowsort == 'i' || presoldata->rowsort == 'd' )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &rowidxsorted, nrows) );
         SCIP_CALL( SCIPallocBufferArray(scip, &rowsparsity, nrows) );
         for( r = 0; r < nrows; ++r )
            rowidxsorted[r] = r;
         if( presoldata->rowsort == 'i' )
         {
            for( r = 0; r < nrows; ++r )
               rowsparsity[r] = SCIPmatrixGetRowNNonzs(matrix, r);
         }
         else if( presoldata->rowsort == 'd' )
         {
            for( r = 0; r < nrows; ++r )
               rowsparsity[r] = -SCIPmatrixGetRowNNonzs(matrix, r);
         }
         SCIPsortIntInt(rowsparsity, rowidxsorted, nrows);
      }
      else
      {
         assert(presoldata->rowsort == 'n');
         rowidxsorted = NULL;
         rowsparsity = NULL;
      }

      /* loop over the rows and cancel non-zeros until maximum number of retrieves is reached */
      maxuseless = (SCIP_Longint)(presoldata->maxretrievefac * (SCIP_Real)nrows);
      nuseless = 0;
      noldchgcoefs = *nchgcoefs;
      for( r = 0; r < nrows && nuseless <= maxuseless; r++ )
      {
         int rowidx;

         rowidx = rowidxsorted != NULL ? rowidxsorted[r] : r;
         SCIP_CALL( cancelRow(scip, matrix, pairtable, rowidx, \
               presoldata->maxcontfillin, presoldata->maxintfillin, presoldata->maxbinfillin, \
               presoldata->maxconsiderednonzeros, \
               &nuseless, nchgcoefs, &numcancel, &nfillin) );
      }

      SCIPfreeBufferArrayNull(scip, &rowsparsity);
      SCIPfreeBufferArrayNull(scip, &rowidxsorted);

      SCIPhashtableFree(&pairtable);
      SCIPfreeBufferArrayNull(scip, &varpairs);

      SCIPfreeBufferArray(scip, &perm);
      SCIPfreeBufferArray(scip, &locks);

      /* update result */
      if( numcancel > 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "   (%.1fs) sparsify %s: %d/%d (%.1f%%) nonzeros canceled, %d changed coefficients, %d added nonzeros\n",
            SCIPgetSolvingTime(scip), (nuseless > maxuseless ? "aborted" : "finished"), numcancel, SCIPmatrixGetNNonzs(matrix),
            100.0*(SCIP_Real)numcancel/(SCIP_Real)SCIPmatrixGetNNonzs(matrix), *nchgcoefs - noldchgcoefs, nfillin);
         *result = SCIP_SUCCESS;
      }
      updateFailureStatistic(presoldata, numcancel > 0);
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeSparsify)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** initialization method of presolver (called after problem was transformed) */
static
SCIP_DECL_PRESOLINIT(presolInitSparsify)
{
   SCIP_PRESOLDATA* presoldata;

   /* we set the nfailures counter in the init (and not in the initpre) callback, because we want it to persist across
    * restarts
    */
   presoldata = SCIPpresolGetData(presol);
   presoldata->nfailures = 0;
   presoldata->nwaitingcalls = 0;

   return SCIP_OKAY;
}

/** creates the sparsify presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolSparsify(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create sparsify presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecSparsify, presoldata) );

   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeSparsify) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitSparsify) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/sparsify/maxcontfillin",
         "maximal fillin for continuous variables",
         &presoldata->maxcontfillin, FALSE, DEFAULT_MAX_CONT_FILLIN, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/sparsify/maxbinfillin",
         "maximal fillin for binary variables",
         &presoldata->maxbinfillin, FALSE, DEFAULT_MAX_BIN_FILLIN, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/sparsify/maxintfillin",
         "maximal fillin for integer variables",
         &presoldata->maxintfillin, FALSE, DEFAULT_MAX_INT_FILLIN, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/sparsify/maxnonzeros",
         "maximal support of one equality to be used for cancelling (-1: no limit)",
         &presoldata->maxnonzeros, TRUE, DEFAULT_MAXNONZEROS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/sparsify/maxconsiderednonzeros",
         "maximal number of considered non-zeros within one row (-1: no limit)",
         &presoldata->maxconsiderednonzeros, TRUE, DEFAULT_MAXCONSIDEREDNONZEROS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip,
         "presolving/sparsify/rowsort",
         "order in which to process inequalities ('n'o sorting, 'i'ncreasing nonzeros, 'd'ecreasing nonzeros)",
         &presoldata->rowsort, TRUE, DEFAULT_ROWSORT, "nid", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/sparsify/maxretrievefac",
         "limit on the number of useless vs. useful hashtable retrieves as a multiple of the number of constraints",
         &presoldata->maxretrievefac, TRUE, DEFAULT_MAXRETRIEVEFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
