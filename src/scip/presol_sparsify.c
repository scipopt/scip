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
SCIP_DECL_ISSETEQ(matrixRowEq)
{
   MATRIXROW* x;
   MATRIXROW* y;
   int i,j,dist;

   x = (MATRIXROW*) a;
   y = (MATRIXROW*) b;

   if( x->cnt > y->cnt )
      SCIPswapPointers((void**) &x, (void**) &y);

   i = 0;
   j = 0;
   dist = 0;

   /* x and y should be equal (with maxdist error) */
   if( x->cnt + maxdist < y->cnt )
      return FALSE;

   while( i < x->cnt && j < y->cnt )
   {
      if( x->ind[i++] != y->ind[j++] )
      {
         ++dist;
         if( dist > maxdist )
            return FALSE;
      }
   }

   return TRUE;
}

/* setcmp callback for the matrixrow */
static
SCIP_DECL_ISSETEQ(matrixRowSubset)
{
   MATRIXROW* x;
   MATRIXROW* y;
   int i,j,dist;

   x = (MATRIXROW*) a;
   y = (MATRIXROW*) b;

   i = 0;
   j = 0;

   /* x should be subset of y except one entry */
   if( x->cnt > y->cnt + maxdist )
      return FALSE;

   dist = 0;

   while( i < x->cnt && j < y->cnt )
   {
      if( x->ind[i] < y->ind[j] )
      {
         ++dist;
         if( dist > maxdist )
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
      ++dist;
      if( dist > maxdist )
         return FALSE;
      ++i;
   }

   return TRUE;
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
   uint8_t*              nonzerosotherrow,   /**< mask of non-zero entries of the other constraint */
   int*                  scatterotherrow,    /**< scatter array holding cleaning information of the other constraint */
   int                   numscatterotherrow, /**< number of entries in scatter array of the other constraint */
   SCIP_Real*            coefsotherrow,      /**< coefficients of the other constraint */
   uint8_t*              nonzerosequality,   /**< mask of non-zero entries of the equality */
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
   uint8_t*              nonzerosotherrow,   /**< mask of non-zero entries of the other constraint */
   int*                  scatterotherrow,    /**< scatter array holding cleaning information of the other constraint */
   int*                  numscatterotherrow, /**< number of entries in scatter array of the other constraint */
   SCIP_Real*            coefsotherrow,      /**< coefficients of the other constraint */
   uint8_t*              nonzerosequality,   /**< mask of non-zero entries of the equality */
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

   /* we prefer cancelling of coefficients of variables with infinite variable bounds
    * if there are not at least twice as much other cancellings */
   if( bestbndcancelled > 1 && (2 * bestallcancelled) >= bestbndcancelled )
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
   uint8_t*              nonzerosrow,        /**< non-zero indexes */
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
void prepareMatrixRows(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   MATRIXROW*            msets,              /**< matrix row sets */
   int*                  mindensity,         /**< minimal number of non-zeros in constraint */
   int*                  maxdensity          /**< maximal number of non-zeros in constraint */
   )
{
   int nrows;
   int i;
   int* rowpnt;
   SCIP_Real* valpnt;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(msets != NULL);
   assert(mindensity != NULL);
   assert(maxdensity != NULL);

   nrows = SCIPmatrixGetNRows(matrix);
   *mindensity = nrows + 1;
   *maxdensity = -1;

   /* sort rows of row major format by variable index */
   for( i = 0; i < nrows; i++ )
   {
      rowpnt = SCIPmatrixGetRowIdxPtr(matrix, i);
      valpnt = SCIPmatrixGetRowValPtr(matrix, i);
      SCIPsortIntReal(rowpnt, valpnt, SCIPmatrixGetRowNNonzs(matrix, i));
   }

   /* fill rows into structure for sgtrie */
   for( i = 0; i < nrows; i++ )
   {
      msets[i].val = SCIPmatrixGetRowValPtr(matrix, i);
      msets[i].ind = SCIPmatrixGetRowIdxPtr(matrix, i);
      msets[i].cnt = SCIPmatrixGetRowNNonzs(matrix, i);
      msets[i].idx = i;

      if( SCIPmatrixGetRowNNonzs(matrix, i) < *mindensity )
         *mindensity = SCIPmatrixGetRowNNonzs(matrix, i);

      if( *maxdensity < SCIPmatrixGetRowNNonzs(matrix, i) )
         *maxdensity = SCIPmatrixGetRowNNonzs(matrix, i);
   }
}

#if 0
/** count the number of arcs of the digraph represented by the constraint matrix */
static
int getNumberDigraphArcs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix              /**< matrix containing the constraints */
   )
{
   int nrows;
   int r;
   int narcs;

   narcs = 0;

   nrows = SCIPmatrixGetNRows(matrix);
   for( r = 0; r < nrows; r++ )
      narcs += SCIPmatrixGetRowNNonzs(matrix, r);

   narcs -= nrows;

   return narcs;
}
#endif

/** loop over constraints and fill directed graph */
static
SCIP_RETCODE fillDigraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int*                  firstvaridxperrow   /**< array to store for each row the index of the first variable */
   )
{
   int nrows;
   int r;
   int* rowpnt;
   int* rowend;
   int lastcolidx;
   int colidx;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(digraph != NULL);
   assert(firstvaridxperrow != NULL);

   nrows = SCIPmatrixGetNRows(matrix);

   for( r = 0; r < nrows; r++ )
   {
      lastcolidx = -1;

      rowpnt = SCIPmatrixGetRowIdxPtr(matrix, r);
      rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, r);
      assert(SCIPmatrixGetRowNNonzs(matrix, r) > 0);

      for( ; (rowpnt < rowend); rowpnt++ )
      {
         colidx = *rowpnt;
         assert(colidx >= 0);

         if( lastcolidx == -1 )
            firstvaridxperrow[r] = colidx;
         else
            SCIP_CALL( SCIPdigraphAddArc(digraph, lastcolidx, colidx, NULL) );

         lastcolidx = colidx;
      }
   }

   return SCIP_OKAY;
}

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
      int c;
      int numcancel;
      int numcancelcurrent;
      int numcancellast;
      int numchangedcoefs;
      SCIP_VAR** consvars;
      SCIP_Real* consvals;
      int* eqnonzs;
      int* eqidxs;
      int numeqs;
      uint8_t* rowsparsified;
      int nmisses;
      SCIP_Real* ratios;
      uint8_t* nonzerosotherrow;
      int* scatterotherrow;
      int numscatterotherrow;
      SCIP_Real* coefsotherrow;
      uint8_t* nonzerosequality;
      int* scatterequality;
      int numscatterequality;
      SCIP_Real* coefsequality;
      int numoverlap;
      int* overlapidxs;
      SCIP_Real scalefac;
      int minnonzcancellations;
      int numobservedrowpairs;
      MATRIXROW* msets;
      SCIP_SGTRIE* sgtrie;
      MATRIXROW** matches;
      int nmatches;
      int rowidx;
      int missfails;
      int scalefails;
      int mindensity;
      int maxdensity;
      int multiplekeyspresent;
      SCIP_Bool stop;
      SCIP_DIGRAPH* digraph;
      int* firstvaridxpercons;
      int* varlocks;
      int* components;
      int ncomponents;
      int* rowtocomp;
      int* rowindexes;
      int* rowscomp;
      int nrowscomp;
      int compnumber;
      int rowcnt;

      numcancel = 0;
      numcancellast = 0;
      numcancelcurrent = 0;
      numchangedcoefs = 0;
      numobservedrowpairs = 0;
      missfails = 0;
      scalefails = 0;
      mindensity = 0;
      maxdensity = 0;
      multiplekeyspresent = 0;
      stop = FALSE;

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

      SCIP_CALL( SCIPallocBufferArray(scip, &varlocks, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &firstvaridxpercons, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &components, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowtocomp, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowindexes, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowscomp, nrows) );

      /* prepare matrix rows for insertion to signature trie(s) */
      prepareMatrixRows(scip, matrix, msets, &mindensity, &maxdensity);

      /* create and fill digraph */
      SCIP_CALL( SCIPdigraphCreate(&digraph, SCIPblkmem(scip), ncols) );
      for( i = 0; i < ncols; i++ )
      {
         varlocks[i] = 4 * (SCIPmatrixGetColNUplocks(matrix, i) + SCIPmatrixGetColNDownlocks(matrix, i));
      }
      SCIP_CALL( SCIPdigraphSetSizes(digraph, varlocks) );
      SCIP_CALL( fillDigraph(scip, matrix, digraph, firstvaridxpercons) );

      /* compute independent components */
      SCIP_CALL( SCIPdigraphComputeUndirectedComponents(digraph, 1, components, &ncomponents) );
      assert(ncomponents >= 1);

      /* assign row to component */
      for( r = 0; r < nrows; r++ )
      {
         rowtocomp[r] = components[firstvaridxpercons[r]];
         rowindexes[r] = r;
      }

      /* sort row indexes by the component number */
      SCIPsortIntInt(rowtocomp, rowindexes, nrows);

      /* loop over all independent components */
      rowcnt = 0;
      int nhashcancel;
      nhashcancel = 0;
      for( c = 0; c < ncomponents; c++ )
      {
         compnumber = rowtocomp[rowcnt];
         nrowscomp = 0;
         while( rowtocomp[rowcnt] == compnumber && rowcnt < nrows )
         {
            rowscomp[nrowscomp] = rowindexes[rowcnt];
            rowcnt++;
            nrowscomp++;
         }

         if( !presoldata->fullsearch )
         {
            /* do no sparsification for large components */
            if( nrowscomp >= 500000 )
               continue;
         }

         /* create signature trie for one component */
         SCIP_CALL( SCIPsgtrieCreate(&sgtrie, SCIPblkmem(scip), SCIPbuffer(scip),
         matrixrowGetSignature, matrixRowEq /*NULL*/, matrixRowSubset /*NULL*/) );

#if 0
         /* insert rows of one component into signature trie */
         for( r = 0; r < nrowscomp; r++ )
         {
            SCIP_RETCODE rc;
            rc = SCIPsgtrieInsert(sgtrie, &msets[rowscomp[r]]);

            if( rc == SCIP_KEYALREADYEXISTING )
            {
#ifdef SCIP_MORE_DEBUG
               MATRIXROW* tmpset;
               tmpset = SCIPsgtrieFindEq(sgtrie, &msets[rowscomp[r]]);
               SCIPdebugMsg(scip, "### Warning: key already present in signature trie!\n");
               SCIPmatrixPrintRow(scip, matrix, msets[rowscomp[r]].idx);
               SCIPmatrixPrintRow(scip, matrix, tmpset->idx);
#endif
               multiplekeyspresent = 1;
            }
         }
#endif

         {
            SCIP_HASHTABLE* pairtable;
            ROWVARPAIR* varpairs;
            int nvarpairs;
            int varpairssize;
            varpairssize = 0;
            nvarpairs = 0;
            varpairs = NULL;
            SCIP_CALL( SCIPhashtableCreate(&pairtable, SCIPblkmem(scip), 1, SCIPhashGetKeyStandard, varPairsEqual, varPairHashval, NULL) );
         /* collect equalities and their number of non-zeros */
         numeqs = 0;
         for( r = 0; r < nrowscomp; r++ )
         {
            if( SCIPisEQ(scip, SCIPmatrixGetRowRhs(matrix,rowscomp[r]), SCIPmatrixGetRowLhs(matrix,rowscomp[r])) )
            {
               int* rowinds;
               SCIP_Real* rowvals;
               int i,j;
               int nnonz;
               int npairs;

               nnonz = SCIPmatrixGetRowNNonzs(matrix, rowscomp[r]);
               rowinds = SCIPmatrixGetRowIdxPtr(matrix, rowscomp[r]);
               rowvals = SCIPmatrixGetRowValPtr(matrix, rowscomp[r]);

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
                     varpairs[nvarpairs].rowindex = rowscomp[r];
                     varpairs[nvarpairs].varindex1 = rowinds[i];
                     varpairs[nvarpairs].varindex2 = rowinds[j];
                     varpairs[nvarpairs].varcoef1 = rowvals[i];
                     varpairs[nvarpairs].varcoef2 = rowvals[j];
                     ++nvarpairs;
                  }
               }
            }
         }

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

         /* sort the equalities by their number of non-zeros */
         /* TODO: check if we actually need this sorting */
//          SCIPsortIntInt(eqnonzs, eqidxs, numeqs);

         /* loop over the rows of the component */

         for( r = 0; r < nrowscomp; r++ )
         {

            SCIP_CALL( cancelRow(scip, matrix, pairtable, rowscomp[r], nchgcoefs, &nhashcancel) );
            continue;

            /* use signature trie to retrieve matching equality row candidates */
#if 0
            SCIP_CALL( SCIPsgtrieFind(sgtrie, (void*)&msets[rowscomp[r]], SCIP_SGTRIE_SUBSET, presoldata->maxsupersetmisses, (void**)matches, &nmatches) );
#else
            nmatches = 0;
#endif
            /* loop over all signature trie matches until the row is sparsified */
            for( i = 0; i < nmatches && !rowsparsified[rowscomp[r]]; i++ )
            {
               rowidx = matches[i]->idx;
               numobservedrowpairs++;

               /* we treat only different constraint pairs and rows which are not already sparsified */
               if( rowidx != rowscomp[r] && !rowsparsified[rowidx] )//&& !rowsparsified[eqidxs[r]] )
               {
                  /* calculate exact number of non-zero misses */
                  nmisses = getMinNumSupersetMisses(scip, matrix, rowscomp[r], rowidx,
                     nonzerosotherrow, scatterotherrow, &numscatterotherrow, coefsotherrow,
                     nonzerosequality, scatterequality, &numscatterequality, coefsequality,
                     &numoverlap, overlapidxs, ratios);

                  if( nmisses <= presoldata->maxsupersetmisses )
                  {
                     /* determine scalefac which removes minimal one non-zero */
                     scalefac = getScalefac(scip, matrix, rowscomp[r], rowidx, numoverlap, overlapidxs, ratios, &minnonzcancellations);
                     assert(minnonzcancellations > 0);

                     /* process row only if no fill-in occur */
                     if( minnonzcancellations > nmisses )
                     {
                        SCIP_CALL( sparsifyCons(scip, matrix, rowscomp[r], rowidx, scalefac, consvars, consvals,
                              nonzerosotherrow, scatterotherrow, numscatterotherrow, coefsotherrow,
                              nonzerosequality, scatterequality, numscatterequality, coefsequality,
                              &numcancel, &numchangedcoefs) );

                        rowsparsified[rowscomp[r]] = TRUE;
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

               if( !presoldata->fullsearch )
               {
                  /* check success of this presolver, stop otherwise */
                  if( !(numobservedrowpairs % CHECK_INTERVAL) )
                  {
                     numcancellast = numcancel - numcancelcurrent;
                     numcancelcurrent = numcancel;
#ifdef SCIP_MORE_DEBUG
                     SCIPdebugMsg(scip, "### cancellast=%d, sucratio=%.6f\n",numcancellast,(SCIP_Real)numcancellast/(SCIP_Real)CHECK_INTERVAL);
#endif
                     if( ((SCIP_Real)numcancellast/(SCIP_Real)CHECK_INTERVAL) < THRESHOLD )
                     {
                        stop = TRUE;
                        break;
                     }

                     numcancellast = 0;
                  }
               }
            }

            if( stop )
               break;
         }


         SCIPhashtableFree(&pairtable);
         SCIPfreeBufferArrayNull(scip, &varpairs);
         }
         /* free signature trie */
         SCIPsgtrieFree(&sgtrie);

         if( stop )
            break;
      }

      {
         SCIP_Real percentagenzcancelled = 100 * nhashcancel / (SCIP_Real)SCIPmatrixGetNNonzs(matrix);
         if( percentagenzcancelled >= 0.1 )
            SCIPinfoMessage(scip, NULL, "cancelled %.1f%% non-zeros\n", 100 * nhashcancel / (SCIP_Real)SCIPmatrixGetNNonzs(matrix));
      }

      /* free local memory */
      SCIPdigraphFree(&digraph);
      SCIPfreeBufferArray(scip, &rowscomp);
      SCIPfreeBufferArray(scip, &rowindexes);
      SCIPfreeBufferArray(scip, &rowtocomp);
      SCIPfreeBufferArray(scip, &components);
      SCIPfreeBufferArray(scip, &firstvaridxpercons);
      SCIPfreeBufferArray(scip, &varlocks);
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
//       SCIPdebugMsg(scip, "### nrows=%d, ntries=%d, mfails=%d, sfails=%d, minrowdensity=%d, maxrowdensity=%d, sgmultkeys=%d, sucratio=%.4f, chgcoefs=%d, nzcancel=%d\n",
//          nrows, ncomponents, missfails, scalefails, mindensity, maxdensity, multiplekeyspresent, (SCIP_Real)numcancel/(SCIP_Real)numobservedrowpairs, numchangedcoefs, numcancel);
#endif
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
