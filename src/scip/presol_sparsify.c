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
#define UPSHOT_INTERVAL             50000    /**< number of observed equalities after which the success is verified */
#define HIT_THRESHOLD    ((SCIP_Real)0.0001) /**< minimal ratio of non-zero cancellation per observed equalities */
#define DEFAULT_MAX_SUPERSET_MISSES     1    /**< default value for the maximal number of superset misses */

/*
 * Data structures
  */

/** control parameters */
struct SCIP_PresolData
{
   int                   maxsupersetmisses;  /**< maximal number of superset misses */
   SCIP_Bool             fullsearch;         /**< flag indicating that full sparsification is required */
};


typedef struct MatrixRow
{
   SCIP_Real*            val;                /**< coefficients of row */
   int*                  ind;                /**< column indexes of row */
   int                   cnt;                /**< number of column entries (non-zeros) */
   int                   idx;                /**< row index */
} MATRIXROW;


/*
 * Local methods
 */

/** signature function */
static
SCIP_DECL_GETSIGNATURE(matrixrowGetSignature)
{
   MATRIXROW* row;
   int k;
   uint64_t signature = 0;

   row = (MATRIXROW*) a;

   /* update the signature with each column index in the row */
   for( k = 0; k < row->cnt; ++k )
      UPDATE_SIGNATURE(signature, row->ind[k]);

   return signature;
}


/* setcmp callback for the matrixrow */
static
SCIP_DECL_SETCMP(cmpMatrixrow)
{
   MATRIXROW* x;
   MATRIXROW* y;
   int i,j;

   x = (MATRIXROW*) a;
   y = (MATRIXROW*) b;

   i = 0;
   j = 0;

   switch(op)
   {
   case SCIP_SGTRIE_SETEQ:
      /* x and y should be equal */
      if( x->cnt != y->cnt )
         return FALSE;

      while( i < x->cnt && j < y->cnt )
      {
         if( x->ind[i++] != y->ind[j++] )
            return FALSE;
      }

      return TRUE;
   case SCIP_SGTRIE_SUBSET:
      /* x should be subset of y */
      if( x->cnt > y->cnt )
         return FALSE;

      while( i < x->cnt && j < y->cnt )
      {
         if( x->ind[i] < y->ind[j] )
            return FALSE;
         else if( x->ind[i] > y->ind[j] )
            ++j;
         else
         {
            ++i;
            ++j;
         }
      }

      if( i < x->cnt )
         return FALSE;

      return TRUE;
   case SCIP_SGTRIE_SUBSETPLUSONE:
      {
         /* x should be subset of y except one entry */
         int distance;

         if( x->cnt > y->cnt + 1 )
            return FALSE;

         distance = 0;

         while( i < x->cnt && j < y->cnt )
         {
            if( x->ind[i] < y->ind[j] )
            {
               ++distance;
               if( distance > 1 )
                  return FALSE;
               ++i;
            }
            else if( x->ind[i] > y->ind[j] )
               ++j;
            else
            {
               ++i;
               ++j;
            }
         }

         while( i < x->cnt )
         {
            ++distance;
            if( distance > 1 )
               return FALSE;
            ++i;
         }

         return TRUE;
      }
   }

   assert(FALSE);
   return FALSE;
}


/** add the scaled equality to the constraint for sparsification */
static
SCIP_RETCODE sparsifyCons(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   otherrow,           /**< index of the other constraint */
   int                   equality,           /**< index of the equality */
   SCIP_Real             scalefac,           /**< scalefactor for sparsification */
   SCIP_VAR**            consvars,           /**< helper array for constraint variables */
   SCIP_Real*            consvals,           /**< helper array for constraint coefficients */
   SCIP_Bool*            nonzerosotherrow,   /**< mask of non-zero entries of the other constraint */
   int*                  scatterotherrow,    /**< scatter array holding cleaning information of the other constraint */
   int                   numscatterotherrow, /**< number of entries in scatter array of the other constraint */
   SCIP_Real*            coefsotherrow,      /**< coefficients of the other constraint */
   SCIP_Bool*            nonzerosequality,   /**< mask of non-zero entries of the equality */
   int*                  scatterequality,    /**< scatter array holding cleaning information of the equality */
   int                   numscatterequality, /**< number of entries in scatter array of the equality */
   SCIP_Real*            coefsequality,      /**< coefficients of the equality */
   int*                  numcancel,          /**< number of deleted non-zeros */
   int*                  numchangedcoefs     /**< number of changed coefficients */
   )
{
   int ncols;
   int i;
   SCIP_CONS* cons;
   char consname[SCIP_MAXSTRLEN];
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real eqside;
   int nconsvars;
   int numaddednonzeros;
   SCIP_VAR* var;
   int ncancel;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= otherrow && otherrow < SCIPmatrixGetNRows(matrix));
   assert(0 <= equality && equality < SCIPmatrixGetNRows(matrix));
   assert(!SCIPisZero(scip, scalefac));
   assert(consvars != NULL);
   assert(consvals != NULL);
   assert(nonzerosotherrow != NULL);
   assert(scatterotherrow != NULL);
   assert(coefsotherrow != NULL);
   assert(nonzerosequality != NULL);
   assert(scatterequality != NULL);
   assert(coefsequality != NULL);
   assert(numcancel != NULL);
   assert(numchangedcoefs != NULL);

   /* get other constraint */
   cons = SCIPmatrixGetCons(matrix, otherrow);
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s", SCIPmatrixGetRowName(matrix, otherrow));
   lhs = SCIPmatrixGetRowLhs(matrix, otherrow);
   rhs = SCIPmatrixGetRowRhs(matrix, otherrow);

   /* get side of equality */
   assert(SCIPisEQ(scip, SCIPmatrixGetRowLhs(matrix, equality), SCIPmatrixGetRowRhs(matrix, equality)));
   eqside = SCIPmatrixGetRowLhs(matrix, equality);

   nconsvars = 0;
   numaddednonzeros = 0;
   ncols = SCIPmatrixGetNColumns(matrix);
   ncancel = 0;

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMsg(scip, "########\n");
   SCIPmatrixPrintRow(scip, matrix, otherrow);
   SCIPmatrixPrintRow(scip, matrix, equality);
   SCIPdebugMsg(scip, "### scalefac: %.2f\n",scalefac);
#endif

   for( i = 0; i < ncols; i++ )
   {
      var = SCIPmatrixGetVar(matrix, i);

      if( nonzerosotherrow[i] && nonzerosequality[i] )
      {
         /* calculate new coefficient of both non-zero entries */
         consvals[nconsvars] = coefsotherrow[i] + scalefac * coefsequality[i];

         /* count complete number of modified coefficients in the matrix */
         if( !SCIPisEQ(scip, consvals[nconsvars], coefsotherrow[i]) )
            (*numchangedcoefs)++;

         if( !SCIPisZero(scip,consvals[nconsvars]) )
         {
            /* update coefficient */
            consvars[nconsvars] = var;
            nconsvars++;
         }
         else
         {
            /* non-zero cancelling happend */
            ncancel++;
         }
      }
      else if( !nonzerosotherrow[i] && nonzerosequality[i] )
      {
         /* fill-in occur */
         numaddednonzeros++;
         consvals[nconsvars] = scalefac * coefsequality[i];
         consvars[nconsvars] = var;
         nconsvars++;
         (*numchangedcoefs)++;
      }
      else if( nonzerosotherrow[i] && !nonzerosequality[i] )
      {
         /* copy */
         consvals[nconsvars] = coefsotherrow[i];
         consvars[nconsvars] = var;
         nconsvars++;
      }
   }

   /* secure that non-zero cancellation has happened and not fill-in */
   assert(ncancel > numaddednonzeros);
   (*numcancel) += ncancel;

   /* determine new lhs and rhs */
   if( !SCIPisZero(scip, eqside) )
   {
      if( !SCIPisInfinity(scip, -lhs) )
         lhs += scalefac * eqside;
      if( !SCIPisInfinity(scip, rhs) )
         rhs += scalefac * eqside;
   }

   /* create sparsified constraint and add it to scip */
   SCIPdelCons(scip, cons);
   SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, nconsvars, consvars, consvals,
         lhs, rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

#ifdef SCIP_MORE_DEBUG
   SCIPdebugPrintCons(scip, cons, NULL);
   SCIPdebugMsg(scip, "########\n");
#endif

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}


/** determine the minimal number of variables indexes being superset misses of the equality */
static
int getMinNumSupersetMisses(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   otherrow,           /**< index of the other constraint */
   int                   equality,           /**< index of the equality */
   SCIP_Bool*            nonzerosotherrow,   /**< mask of non-zero entries of the other constraint */
   int*                  scatterotherrow,    /**< scatter array holding cleaning information of the other constraint */
   int*                  numscatterotherrow, /**< number of entries in scatter array of the other constraint */
   SCIP_Real*            coefsotherrow,      /**< coefficients of the other constraint */
   SCIP_Bool*            nonzerosequality,   /**< mask of non-zero entries of the equality */
   int*                  scatterequality,    /**< scatter array holding cleaning information of the equality */
   int*                  numscatterequality, /**< number of entries in scatter array of the equality */
   SCIP_Real*            coefsequality,      /**< coefficients of the equality */
   int*                  numoverlap,         /**< number of overlapping variable indexes */
   int*                  overlapidxs,        /**< overlapping indexes */
   SCIP_Real*            ratios              /**< coefficient ratios */
   )
{
   int nsupersetmisses;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   int ncols;
   int col;
   int i;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= otherrow && otherrow < SCIPmatrixGetNRows(matrix));
   assert(0 <= equality && equality < SCIPmatrixGetNRows(matrix));
   assert(nonzerosotherrow != NULL);
   assert(scatterotherrow != NULL);
   assert(numscatterotherrow != NULL);
   assert(coefsotherrow != NULL);
   assert(nonzerosequality != NULL);
   assert(scatterequality != NULL);
   assert(numscatterequality != NULL);
   assert(coefsequality != NULL);
   assert(numoverlap != NULL);
   assert(overlapidxs != NULL);
   assert(ratios != NULL);

   ncols = SCIPmatrixGetNColumns(matrix);
   nsupersetmisses = 0;
   *numoverlap = 0;
   *numscatterotherrow = 0;
   *numscatterequality = 0;

   /* read data of other row and scatter it */
   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, otherrow);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, otherrow);
   valpnt = SCIPmatrixGetRowValPtr(matrix, otherrow);
   for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
   {
      col = *rowpnt;
      nonzerosotherrow[col] = TRUE;
      coefsotherrow[col] = *valpnt;
      scatterotherrow[*numscatterotherrow] = col;
      (*numscatterotherrow)++;
   }

   /* read data of equality and scatter it */
   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, equality);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, equality);
   valpnt = SCIPmatrixGetRowValPtr(matrix, equality);
   for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
   {
      col = *rowpnt;
      nonzerosequality[col] = TRUE;
      coefsequality[col] = *valpnt;
      scatterequality[*numscatterequality] = col;
      (*numscatterequality)++;
   }

   /* determine overlap and superset misses */
   for( i = 0; i < ncols; i++ )
   {
      if( nonzerosotherrow[i] && nonzerosequality[i] )
      {
         overlapidxs[*numoverlap] = i;
         ratios[*numoverlap] = coefsotherrow[i] / coefsequality[i];
         (*numoverlap)++;
      }
      else if( !nonzerosotherrow[i] && nonzerosequality[i] )
      {
         nsupersetmisses++;

         /* we stop counting the misses if more than one superset miss is present */
         if( nsupersetmisses > 1 )
            break;
      }
   }

   return nsupersetmisses;
}


/** determine scalefac for non-zero cancellation */
static
SCIP_Real getScalefac(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   otherrow,           /**< index of the other constraint */
   int                   equality,           /**< index of the equality */
   int                   numoverlap,         /**< number of overlapping coefficients */
   int*                  overlapidxs,        /**< indexes of overlapping columns */
   SCIP_Real*            ratios,             /**< ratios of overlapping columns */
   int*                  nonzcancellations   /**< number of possible non-zeros cancellations */
   )
{
   SCIP_VAR* var;
   SCIP_Real bestbndscalefac;
   SCIP_Real bestallscalefac;
   int bestbndcancelled;
   int bestallcancelled;
   SCIP_Real scalefac;
   int cntbndvars;
   int cntallvars;
   SCIP_Real currentratio;
   int i;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= otherrow && otherrow < SCIPmatrixGetNRows(matrix));
   assert(0 <= equality && equality < SCIPmatrixGetNRows(matrix));
   assert(numoverlap > 0);
   assert(overlapidxs != NULL);
   assert(ratios != NULL);
   assert(nonzcancellations != NULL);

   *nonzcancellations = 0;

   /* sort ratios */
   SCIPsortRealInt(ratios, overlapidxs, numoverlap);

   cntbndvars = 0;
   cntallvars = 0;
   currentratio = ratios[0];
   bestbndscalefac = 0.0;
   bestallscalefac = 0.0;
   bestbndcancelled = 0;
   bestallcancelled = 0;

   /* parse the sorted ratios for an appropriate scale factor */
   for( i = 0; i < numoverlap; i++ )
   {
      var = SCIPmatrixGetVar(matrix, overlapidxs[i]);
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      if( !SCIPisEQ(scip, currentratio, ratios[i]) )
      {
         cntbndvars = 0;
         cntallvars = 0;
         currentratio = ratios[i];
      }

      if( SCIPisInfinity(scip, ub) || SCIPisInfinity(scip,-lb) )
      {
         cntbndvars++;

         if( bestbndcancelled < cntbndvars )
         {
            bestbndscalefac = -ratios[i];
            bestbndcancelled = cntbndvars;
         }
      }

      cntallvars++;

      if( bestallcancelled < cntallvars )
      {
         bestallscalefac = -ratios[i];
         bestallcancelled = cntallvars;
      }
   }

   /* we prefer cancelling of coefficients with infinite variable bounds */
   if( bestbndcancelled > 1 )
   {
      assert(!SCIPisZero(scip, bestbndscalefac));
      scalefac = bestbndscalefac;
      *nonzcancellations = bestbndcancelled;
   }
   else
   {
      assert(!SCIPisZero(scip, bestallscalefac));
      scalefac = bestallscalefac;
      *nonzcancellations = bestallcancelled;
   }

   return scalefac;
}


/** clean non-zeros entries in scatter array */
static
void rescattering(
   SCIP_Bool*            nonzerosrow,        /**< non-zero indexes */
   int*                  scatterrow,         /**< cleaning information */
   int*                  numscatterrow       /**< number of scatter entries */
   )
{
   int i;
   for( i = 0; i < *numscatterrow; i++ )
      nonzerosrow[scatterrow[i]] = FALSE;

   *numscatterrow = 0;
}


/** fill matrix data into sets for sgtrie */
static
void getMatrixRowSets(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   MATRIXROW*            msets,              /**< matrix row sets */
   int*                  mindensity,         /**< minimal number of non-zeros in constraint */
   int*                  maxdensity          /**< maximal number of non-zeros in constraint */
   )
{
   int nrows;
   int r;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(msets != NULL);
   assert(mindensity != NULL);
   assert(maxdensity != NULL);

   nrows = SCIPmatrixGetNRows(matrix);
   *mindensity = nrows + 1;
   *maxdensity = -1;

   for( r = 0; r < nrows; r++ )
   {
      msets[r].val = SCIPmatrixGetRowValPtr(matrix, r);
      msets[r].ind = SCIPmatrixGetRowIdxPtr(matrix, r);
      msets[r].cnt = SCIPmatrixGetRowNNonzs(matrix, r);
      msets[r].idx = r;

      if( SCIPmatrixGetRowNNonzs(matrix, r) < *mindensity )
         *mindensity = SCIPmatrixGetRowNNonzs(matrix, r);

      if( *maxdensity < SCIPmatrixGetRowNNonzs(matrix, r) )
         *maxdensity = SCIPmatrixGetRowNNonzs(matrix, r);
   }
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecSparsify)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   /* we only work on pure MIPs currently */
   if( initialized && complete )
   {
      int ncols;
      int nrows;
      int r;
      int i;
      int numcancel;
      int numchangedcoefs;
      SCIP_VAR** consvars;
      SCIP_Real* consvals;
      int* eqnonzs;
      int* eqidxs;
      int numeqs;
      SCIP_Bool* rowsparsified;
      SCIP_Bool eqprocessed;
      int nmisses;
      SCIP_Real* ratios;
      SCIP_Bool* nonzerosotherrow;
      int* scatterotherrow;
      int numscatterotherrow;
      SCIP_Real* coefsotherrow;
      SCIP_Bool* nonzerosequality;
      int* scatterequality;
      int numscatterequality;
      SCIP_Real* coefsequality;
      int numoverlap;
      int* overlapidxs;
      SCIP_Real scalefac;
      int minnonzcancellations;
      int numobservedequalities;
      MATRIXROW* msets;
      SCIP_SGTRIE* sgtrie;
      MATRIXROW** matches;
      int nmatches;
      int rowidx;
      int missfails;
      int scalefails;
      int mindensity;
      int maxdensity;

      numcancel = 0;
      numchangedcoefs = 0;
      numobservedequalities = 0;
      missfails = 0;
      scalefails = 0;
      mindensity = 0;
      maxdensity = 0;

      nrows = SCIPmatrixGetNRows(matrix);
      ncols = SCIPmatrixGetNColumns(matrix);

      SCIP_CALL( SCIPallocBufferArray(scip, &eqnonzs, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &eqidxs, nrows) );

      SCIP_CALL( SCIPallocBufferArray(scip, &rowsparsified, nrows) );
      BMSclearMemoryArray(rowsparsified, nrows);

      SCIP_CALL( SCIPallocBufferArray(scip, &ratios, ncols) );

      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consvals, ncols) );

      SCIP_CALL( SCIPallocBufferArray(scip, &coefsotherrow, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &scatterotherrow, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nonzerosotherrow, ncols) );
      BMSclearMemoryArray(nonzerosotherrow, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &coefsequality, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &scatterequality, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nonzerosequality, ncols) );
      BMSclearMemoryArray(nonzerosequality, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &overlapidxs, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &msets, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &matches, nrows));

      SCIP_CALL( SCIPsgtrieCreate(&sgtrie, SCIPblkmem(scip), SCIPbuffer(scip), matrixrowGetSignature, cmpMatrixrow) );

      getMatrixRowSets(scip, matrix, msets, &mindensity, &maxdensity);

      /* insert rows into signature trie */
      for( r = 0; r < nrows; r++ )
      {
         SCIP_RETCODE rc;
         rc = SCIPsgtrieInsert(sgtrie, &msets[r]);
#ifdef SCIP_DEBUG
         if( rc == SCIP_KEYALREADYEXISTING )
         {
            MATRIXROW* tmpset;
            tmpset = SCIPsgtrieFindEq(sgtrie, &msets[r]);
            SCIPdebugMsg(scip, "### Warning: key already present in sgtrie!\n");
            SCIPmatrixPrintRow(scip, matrix, msets[r].idx);
            SCIPmatrixPrintRow(scip, matrix, tmpset->idx);
         }
#endif
      }

      /* collect equalities and their number of non-zeros */
      numeqs = 0;
      for( r = 0; r < nrows; r++ )
      {
         /* we do not consider equalities with too less non-zeros if full search is not required */
         if( !presoldata->fullsearch && SCIPmatrixGetRowNNonzs(matrix, r) < MIN_EQS_NONZEROS )
            continue;

         if( SCIPisEQ(scip, SCIPmatrixGetRowRhs(matrix,r), SCIPmatrixGetRowLhs(matrix,r)) )
         {
            eqnonzs[numeqs] = SCIPmatrixGetRowNNonzs(matrix, r);
            eqidxs[numeqs] = r;
            numeqs++;
         }
      }

      /* sort the equalities by their number of non-zeros */
      SCIPsortIntInt(eqnonzs, eqidxs, numeqs);

      /* loop over the equalities by increasing number of non-zeros */
      for( r = 0; r < numeqs; r++ )
      {
         eqprocessed = FALSE;
         numobservedequalities++;

         /* use signature trie to retrieve matching row candidates */
         if( presoldata->maxsupersetmisses == 1 )
            SCIP_CALL( SCIPsgtrieFindSupersetsPlusOne(sgtrie, (void*)&msets[eqidxs[r]], (void**)matches, &nmatches) );
         else
            SCIP_CALL( SCIPsgtrieFindSupersets(sgtrie, (void*)&msets[eqidxs[r]], (void**)matches, &nmatches) );

         /* loop over all signature trie matches */
         for( i = 0; i < nmatches; i++ )
         {
            rowidx = matches[i]->idx;

            /* we treat only different constraint pairs and rows which are not already sparsified */
            if( rowidx != eqidxs[r] && !rowsparsified[rowidx] && !rowsparsified[eqidxs[r]] )
            {
               /* calculate exact number of non-zero misses */
               nmisses = getMinNumSupersetMisses(scip, matrix, rowidx, eqidxs[r],
                  nonzerosotherrow, scatterotherrow, &numscatterotherrow, coefsotherrow,
                  nonzerosequality, scatterequality, &numscatterequality, coefsequality,
                  &numoverlap, overlapidxs, ratios);

               if( nmisses <= 1 )
               {
                  /* determine scalefac which removes minimal one non-zero */
                  scalefac = getScalefac(scip, matrix, rowidx, eqidxs[r], numoverlap, overlapidxs, ratios, &minnonzcancellations);
                  assert(minnonzcancellations > 0);

                  /* process row only if no fill-in occur */
                  if( minnonzcancellations > nmisses )
                  {
                     SCIP_CALL( sparsifyCons(scip, matrix, rowidx, eqidxs[r], scalefac, consvars, consvals,
                           nonzerosotherrow, scatterotherrow, numscatterotherrow, coefsotherrow,
                           nonzerosequality, scatterequality, numscatterequality, coefsequality,
                           &numcancel, &numchangedcoefs) );

                     rowsparsified[rowidx] = TRUE;

                     /* if full search is not required, treat next equality */
                     if( !presoldata->fullsearch )
                        eqprocessed = TRUE;
                  }
                  else
                     scalefails++;
               }
               else
                  missfails++;

               /* clear scatter arrays */
               rescattering(nonzerosotherrow, scatterotherrow, &numscatterotherrow);
               rescattering(nonzerosequality, scatterequality, &numscatterequality);
            }

            if( eqprocessed )
               break;
         }

         if( !presoldata->fullsearch )
         {
            /* check success of this presolver */
            if( !(numobservedequalities % UPSHOT_INTERVAL) )
               if( (SCIP_Real)numcancel/(SCIP_Real)numobservedequalities < HIT_THRESHOLD )
                  break;
         }
      }

      /* free local memory */
      SCIPsgtrieFree(&sgtrie);
      SCIPfreeBufferArray(scip, &matches);
      SCIPfreeBufferArray(scip, &msets);
      SCIPfreeBufferArray(scip, &overlapidxs);
      SCIPfreeBufferArray(scip, &nonzerosequality);
      SCIPfreeBufferArray(scip, &scatterequality);
      SCIPfreeBufferArray(scip, &coefsequality);
      SCIPfreeBufferArray(scip, &nonzerosotherrow);
      SCIPfreeBufferArray(scip, &scatterotherrow);
      SCIPfreeBufferArray(scip, &coefsotherrow);
      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
      SCIPfreeBufferArray(scip, &ratios);
      SCIPfreeBufferArray(scip, &rowsparsified);
      SCIPfreeBufferArray(scip, &eqidxs);
      SCIPfreeBufferArray(scip, &eqnonzs);

      /* update complete number of modified coefficients */
      (*nchgcoefs) += numchangedcoefs;

      if( numcancel > 0 )
         *result = SCIP_SUCCESS;

#ifdef SCIP_DEBUG
      /* print out number of changed coefficients and non-zero cancellations */
      SCIPdebugMsg(scip, "### nrows=%d, eqs=%d, mfails=%d, sfails=%d, mindensity=%d, maxdensity=%d, chgcoefs=%d, nzcancel=%d\n",
         nrows, numeqs, missfails, scalefails, mindensity, maxdensity, numchangedcoefs, numcancel);
#endif
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopySparsify)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolSparsify(scip) );

   return SCIP_OKAY;
}

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

   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopySparsify) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeSparsify) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/sparsify/maxsupersetmisses",
         "maximal number of superset misses",
         &presoldata->maxsupersetmisses, TRUE, DEFAULT_MAX_SUPERSET_MISSES, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/sparsify/fullsearch",
         "require full search for sparsification",
         &presoldata->fullsearch, TRUE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
