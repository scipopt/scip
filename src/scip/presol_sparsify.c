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
#define SCIP_DEBUG
/**@file   presol_sparsify.c
 * @brief  cancel non-zeros of the constraint matrix
 * @author Dieter Weninger
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
#include "scip/sgtrie.h"
#include "presol_sparsify.h"

#define PRESOL_NAME            "sparsify"
#define PRESOL_DESC            "eliminate non-zero coefficients"
#define PRESOL_PRIORITY            -24000    /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               -1    /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define MIN_EQS_NONZEROS                3    /**< minimal number of non-zeros of the equalities */
#define CHECK_INTERVAL               1000    /**< number of observed row pairs after which the success is verified */
#define THRESHOLD        ((SCIP_Real)0.02)   /**< minimal ratio of non-zero cancellation per observed row pairs */
#define DEFAULT_MAX_SUPERSET_MISSES     1    /**< default value for the maximal number of superset misses */
#define DEFAULT_FULL_SEARCH             0    /**< default value for full search */

/*
 * Data structures
  */

/** control parameters */
struct SCIP_PresolData
{
   int                   maxsupersetmisses;  /**< maximal number of superset misses */
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
   int*                  nchgcoefs,
   int*                  ncanceled
   )
{
   int* cancelrowinds;
   SCIP_Real* cancelrowvals;
   SCIP_Real cancellhs;
   SCIP_Real cancelrhs;
   int* tmpinds;
   SCIP_Real* tmpvals;
   int cancelrowlen;
   int nchgcoef;

   cancelrowlen = SCIPmatrixGetRowNNonzs(matrix, rowidx);

   SCIPduplicateBufferArray(scip, &cancelrowinds, SCIPmatrixGetRowIdxPtr(matrix, rowidx), cancelrowlen);
   SCIPduplicateBufferArray(scip, &cancelrowvals, SCIPmatrixGetRowValPtr(matrix, rowidx), cancelrowlen);
   SCIPallocBufferArray(scip, &tmpinds, cancelrowlen);
   SCIPallocBufferArray(scip, &tmpvals, cancelrowlen);

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

      for( i = 0; i < cancelrowlen; ++i )
      {
         for( j = i + 1; j < cancelrowlen; ++j )
         {
            int a,b;
            int ncancel;
            int nfillin;
            int eqrowlen;
            ROWVARPAIR* eqrowvarpair;
            SCIP_Real* eqrowvals;
            int* eqrowinds;
            SCIP_Real scale;

            assert(cancelrowinds[i] < cancelrowinds[j]);

            rowvarpair.varindex1 = cancelrowinds[i];
            rowvarpair.varindex2 = cancelrowinds[j];
            rowvarpair.varcoef1 = cancelrowvals[i];
            rowvarpair.varcoef2 = cancelrowvals[j];

            eqrowvarpair = SCIPhashtableRetrieve(pairtable, (void*) &rowvarpair);

            if( eqrowvarpair == NULL || eqrowvarpair->rowindex == rowidx )
               continue;

            eqrowvals = SCIPmatrixGetRowValPtr(matrix, eqrowvarpair->rowindex);
            eqrowinds = SCIPmatrixGetRowIdxPtr(matrix, eqrowvarpair->rowindex);
            eqrowlen = SCIPmatrixGetRowNNonzs(matrix, eqrowvarpair->rowindex);

            scale = -rowvarpair.varcoef1 / eqrowvarpair->varcoef1;

            a = 0;
            b = 0;
            ncancel = 0;
            nfillin = 0;
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
                  ++nfillin;
                  ++b;
               }
            }

            nfillin += (eqrowlen - b);

            if( (ncancel - nfillin) > (bestncancel - bestnfillin) ||
               ((ncancel - nfillin) == (bestncancel - bestnfillin) && nfillin < bestnfillin) )
            {
               bestncancel = ncancel;
               bestnfillin = nfillin;
               bestcand = eqrowvarpair->rowindex;
               bestscale = scale;
            }

            ++i;
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
   SCIP_HASHTABLE* pairtable;
   ROWVARPAIR* varpairs;
   int nvarpairs;
   int varpairssize;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

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
                  assert(nvarpairs < varpairssize);
                  varpairs[nvarpairs].rowindex = r;
                  varpairs[nvarpairs].varindex1 = rowinds[i];
                  varpairs[nvarpairs].varindex2 = rowinds[j];
                  varpairs[nvarpairs].varcoef1 = rowvals[i];
                  varpairs[nvarpairs].varcoef2 = rowvals[j];
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
            if( SCIPmatrixGetRowNNonzs(matrix, othervarpair->rowindex) <= SCIPmatrixGetRowNNonzs(matrix, varpairs[r].rowindex) )
               continue;
         }

         SCIP_CALL( SCIPhashtableInsert(pairtable, (void*) &varpairs[r]) );
      }

      /* loop over the rows and cancel non-zeros */
      for( r = 0; r < nrows; r++ )
      {
         SCIP_CALL( cancelRow(scip, matrix, pairtable, r, nchgcoefs, &numcancel) );
      }

      SCIPhashtableFree(&pairtable);
      SCIPfreeBufferArrayNull(scip, &varpairs);

      {
         SCIP_Real percentagenzcancelled = 100 * numcancel / (SCIP_Real)SCIPmatrixGetNNonzs(matrix);
         if( percentagenzcancelled >= 0.1 )
            SCIPinfoMessage(scip, NULL, "cancelled %.1f%% non-zeros\n", 100 * numcancel / (SCIP_Real)SCIPmatrixGetNNonzs(matrix));
      }

      /* update result */
      if( numcancel > 0 )
         *result = SCIP_SUCCESS;
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
         "presolving/sparsify/maxsupersetmisses",
         "maximal number of superset misses",
         &presoldata->maxsupersetmisses, TRUE, DEFAULT_MAX_SUPERSET_MISSES, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/sparsify/fullsearch",
         "require full search for sparsification",
         &presoldata->fullsearch, TRUE, DEFAULT_FULL_SEARCH, NULL, NULL) );

   return SCIP_OKAY;
}
