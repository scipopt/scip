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

#define MIN_EQS_NONZEROS                3    /**< minimal number of non-zeros of the equalities */
#define CHECK_INTERVAL               1000    /**< number of observed row pairs after which the success is verified */
#define THRESHOLD        ((SCIP_Real)0.02)   /**< minimal ratio of non-zero cancellation per observed row pairs */
#define DEFAULT_MAX_CONT_FILLIN         0    /**< default value for the maximal fillin for continuous variables */
#define DEFAULT_MAX_BIN_FILLIN          0    /**< default value for the maximal fillin for binary variables */
#define DEFAULT_MAX_INT_FILLIN          0    /**< default value for the maximal fillin for integer variables */
#define DEFAULT_FULL_SEARCH             0    /**< default value for full search */


#define MAXCONSIDEREDNONZEROS 70
/*
 * Data structures
  */

/** control parameters */
struct SCIP_PresolData
{
   int                   maxcontfillin;      /**< maximal fillin for continuous variables */
   int                   maxintfillin;       /**< maximal fillin for integer variables*/
   int                   maxbinfillin;       /**< maximal fillin for binary variables */
   SCIP_Bool             fullsearch;         /**< flag indicating that full sparsification is required */
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

static
SCIP_DECL_HASHKEYEQ(varPairsEqual)
{
   ROWVARPAIR* varpair1;
   ROWVARPAIR* varpair2;
   SCIP_Real scale;

   varpair1 = (ROWVARPAIR*) key1;
   varpair2 = (ROWVARPAIR*) key2;

   if( varpair1->varindex1 != varpair2->varindex1 )
      return FALSE;

   if( varpair1->varindex2 != varpair2->varindex2 )
      return FALSE;

   scale = varpair1->varcoef1 / varpair2->varcoef1;

   if( !EPSEQ(varpair1->varcoef2, scale * varpair2->varcoef2, SCIP_DEFAULT_EPSILON) )
      return FALSE;

   return TRUE;
}

static
SCIP_DECL_HASHKEYVAL(varPairHashval)
{
   ROWVARPAIR* varpair;

   varpair = (ROWVARPAIR*) key;

   return SCIPhashTwo(SCIPcombineTwoInt(varpair->varindex1, varpair->varindex2),
                      SCIPrealHashCode(varpair->varcoef2 / varpair->varcoef1));
}

static
SCIP_RETCODE cancelRow(
   SCIP*                 scip,
   SCIP_MATRIX*          matrix,
   SCIP_HASHTABLE*       pairtable,
   int                   rowidx,
   unsigned int          maxcontfillin,
   unsigned int          maxintfillin,
   unsigned int          maxbinfillin,
   int*                  nchgcoefs,
   int*                  ncanceled
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
   int nchgcoef;
   SCIP_Bool rowiseq;

   rowiseq = SCIPisEQ(scip, SCIPmatrixGetRowLhs(matrix, rowidx), SCIPmatrixGetRowRhs(matrix, rowidx));

   cancelrowlen = SCIPmatrixGetRowNNonzs(matrix, rowidx);

   SCIPduplicateBufferArray(scip, &cancelrowinds, SCIPmatrixGetRowIdxPtr(matrix, rowidx), cancelrowlen);
   SCIPduplicateBufferArray(scip, &cancelrowvals, SCIPmatrixGetRowValPtr(matrix, rowidx), cancelrowlen);
   SCIPallocBufferArray(scip, &tmpinds, cancelrowlen);
   SCIPallocBufferArray(scip, &tmpvals, cancelrowlen);
   SCIPallocBufferArray(scip, &locks, cancelrowlen);

   cancellhs = SCIPmatrixGetRowLhs(matrix, rowidx);
   cancelrhs = SCIPmatrixGetRowRhs(matrix, rowidx);

   nchgcoef = 0;
   while(TRUE)
   {
      SCIP_Real bestscale;
      int bestcand = -1;
      int bestncancel = 0;
      int bestnfillin = 0;
      int i;
      int j;
      ROWVARPAIR rowvarpair;
      int maxlen;

      for( i = 0; i < cancelrowlen; ++i )
      {
         tmpinds[i] = i;
         locks[i] = SCIPmatrixGetColNDownlocks(matrix, cancelrowinds[i]) + SCIPmatrixGetColNUplocks(matrix, cancelrowinds[i]);
      }

      SCIPsortIntInt(locks, tmpinds, cancelrowlen);

      maxlen = MIN(cancelrowlen, MAXCONSIDEREDNONZEROS);

      for( i = 0; i < maxlen; ++i )
      {
         for( j = i + 1; j < maxlen; ++j )
         {
            int a,b;
            int ncancel;
            unsigned int ncontfillin;
            unsigned int nintfillin;
            unsigned int nbinfillin;
            int nfillin;
            int eqrowlen;
            ROWVARPAIR* eqrowvarpair;
            SCIP_Real* eqrowvals;
            int* eqrowinds;
            SCIP_Real scale;
            int i1,i2;

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

            eqrowvarpair = SCIPhashtableRetrieve(pairtable, (void*) &rowvarpair);

            /* if the row we want to cancel is an equality, we will only use equalities
             * with less non-zeros and if the number of non-zeros is equal we use the
             * rowindex as tie-breaker to avoid cyclic non-zero cancellation
             */
            if( eqrowvarpair == NULL || eqrowvarpair->rowindex == rowidx )
               continue;

            eqrowlen = SCIPmatrixGetRowNNonzs(matrix, eqrowvarpair->rowindex);
            if( rowiseq && (cancelrowlen < eqrowlen || (cancelrowlen == eqrowlen && rowidx < eqrowvarpair->rowindex)) )
               continue;

            eqrowvals = SCIPmatrixGetRowValPtr(matrix, eqrowvarpair->rowindex);
            eqrowinds = SCIPmatrixGetRowIdxPtr(matrix, eqrowvarpair->rowindex);

            scale = -rowvarpair.varcoef1 / eqrowvarpair->varcoef1;

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

            nfillin = nbinfillin + nintfillin + ncontfillin;

            if( ((ncancel - nfillin) > (bestncancel - bestnfillin) ||
               ((ncancel - nfillin) == (bestncancel - bestnfillin) && nfillin < bestnfillin) ) )
            {
               bestncancel = ncancel;
               bestnfillin = nfillin;
               bestcand = eqrowvarpair->rowindex;
               bestscale = scale;
            }
         }
      }

      if( bestcand != -1 )
      {
         int a,b;
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

      SCIPallocBufferArray(scip, &consvars, cancelrowlen);

      for( i = 0; i < cancelrowlen; ++i )
         consvars[i] = SCIPmatrixGetVar(matrix, cancelrowinds[i]);

      /* create sparsified constraint and add it to scip */
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, SCIPmatrixGetRowName(matrix, rowidx), cancelrowlen, consvars, cancelrowvals,
                                      cancellhs, cancelrhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPdelCons(scip, SCIPmatrixGetCons(matrix, rowidx)) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      #ifdef SCIP_MORE_DEBUG
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIPdebugMsg(scip, "########\n");
      #endif

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      *nchgcoefs += nchgcoef;
      *ncanceled += SCIPmatrixGetRowNNonzs(matrix, rowidx) - cancelrowlen;

      SCIPfreeBufferArray(scip, &consvars);
   }

   SCIPfreeBufferArray(scip, &locks);
   SCIPfreeBufferArray(scip, &tmpvals);
   SCIPfreeBufferArray(scip, &tmpinds);
   SCIPfreeBufferArray(scip, &cancelrowvals);
   SCIPfreeBufferArray(scip, &cancelrowinds);

   return SCIP_OKAY;
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
   int* locks;
   int* perm;
   SCIP_HASHTABLE* pairtable;
   ROWVARPAIR* varpairs;
   int nvarpairs;
   int varpairssize;
   SCIP_PRESOLDATA* presoldata;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   presoldata = SCIPpresolGetData(presol);

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

      SCIPallocBufferArray(scip, &locks, SCIPmatrixGetNColumns(matrix));
      SCIPallocBufferArray(scip, &perm, SCIPmatrixGetNColumns(matrix));

      /* loop over all rows and create var pairs */
      numcancel = 0;
      varpairssize = 0;
      nvarpairs = 0;
      varpairs = NULL;
      SCIP_CALL( SCIPhashtableCreate(&pairtable, SCIPblkmem(scip), 1, SCIPhashGetKeyStandard, varPairsEqual, varPairHashval, NULL) );
      /* collect equalities and their number of non-zeros */
      for( r = 0; r < nrows; r++ )
      {
         if( SCIPisEQ(scip, SCIPmatrixGetRowRhs(matrix,r), SCIPmatrixGetRowLhs(matrix,r)) )
         {
            int* rowinds;
            SCIP_Real* rowvals;
            int nnonz;
            int npairs;

            nnonz = SCIPmatrixGetRowNNonzs(matrix, r);
            rowinds = SCIPmatrixGetRowIdxPtr(matrix, r);
            rowvals = SCIPmatrixGetRowValPtr(matrix, r);

            for( i = 0; i < nnonz; ++i )
            {
               perm[i] = i;
               locks[i] = SCIPmatrixGetColNDownlocks(matrix, rowinds[i]) + SCIPmatrixGetColNUplocks(matrix, rowinds[i]);
            }

            SCIPsortIntInt(locks, perm, nnonz);

            nnonz = MIN(nnonz, MAXCONSIDEREDNONZEROS);

            npairs = (nnonz * (nnonz - 1)) / 2;

            if( nvarpairs + npairs > varpairssize )
            {
               int newsize = SCIPcalcMemGrowSize(scip, nvarpairs + npairs);
               SCIPreallocBufferArray(scip, &varpairs, newsize);
               varpairssize = newsize;
            }

            for( i = 0; i < nnonz; ++i )
            {
               for( j = i + 1; j < nnonz; ++j )
               {
                  int i1,i2;
                  assert(nvarpairs < varpairssize);

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

         othervarpair = SCIPhashtableRetrieve(pairtable, (void*) &varpairs[r]);

         if( othervarpair != NULL )
         {
            int thisvarpairlocks;
            int othervarpairlocks;

            thisvarpairlocks = SCIPmatrixGetColNDownlocks(matrix, varpairs[r].varindex1) +
                               SCIPmatrixGetColNDownlocks(matrix, varpairs[r].varindex2) +
                               SCIPmatrixGetColNUplocks(matrix, varpairs[r].varindex1) +
                               SCIPmatrixGetColNUplocks(matrix, varpairs[r].varindex2);

            othervarpairlocks = SCIPmatrixGetColNDownlocks(matrix, othervarpair->varindex1) +
                                SCIPmatrixGetColNDownlocks(matrix, othervarpair->varindex2) +
                                SCIPmatrixGetColNUplocks(matrix, othervarpair->varindex1) +
                                SCIPmatrixGetColNUplocks(matrix, othervarpair->varindex2);

            /* only override old var pair if this one has more locks or if it has the same number of locks but its row is sparser */
            if( othervarpairlocks < thisvarpairlocks || (othervarpairlocks == thisvarpairlocks &&
               SCIPmatrixGetRowNNonzs(matrix, othervarpair->rowindex) <= SCIPmatrixGetRowNNonzs(matrix, varpairs[r].rowindex)) )
               continue;
         }

         SCIP_CALL( SCIPhashtableInsert(pairtable, (void*) &varpairs[r]) );
      }

      /* loop over the rows and cancel non-zeros */
      for( r = 0; r < nrows; r++ )
      {
         SCIP_CALL( cancelRow(scip, matrix, pairtable, r, \
                              (unsigned int)presoldata->maxcontfillin, (unsigned int)presoldata->maxintfillin, (unsigned int)presoldata->maxbinfillin, \
                              nchgcoefs, &numcancel) );
      }

      SCIPhashtableFree(&pairtable);
      SCIPfreeBufferArrayNull(scip, &varpairs);

      SCIPfreeBufferArray(scip, &perm);
      SCIPfreeBufferArray(scip, &locks);

      /* update result */
      if( numcancel > 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "cancelled %.1f%%(%i) non-zeros\n", 100 * numcancel / (SCIP_Real)SCIPmatrixGetNNonzs(matrix), numcancel);
         *result = SCIP_SUCCESS;
      }
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

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/sparsify/fullsearch",
         "require full search for sparsification",
         &presoldata->fullsearch, TRUE, DEFAULT_FULL_SEARCH, NULL, NULL) );

   return SCIP_OKAY;
}
